/**
 *    @author A. Muennich, CERN (original TauFinder).
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP (adaptation).
 *
 *  Procedure:
 *    - (Optional:) Discard background particles (p_t_cut_, max_cos_theta_).
 *    - Identify seeds (min_p_t_seed_).
 *    - Define tau and isolation cones centered at the  seed
 *       (search_cone_angle_, isolation_cone_angle_).
 *    - Merge in the case of close by seeds.
 *    - Collect input variables for e.g. a BDT:
 *       + invariant tau mass
 *       + isolation energy
 *       + leading track D0
 *       + leading track Z0
 *       + number of particle tracks (charged particles)
 *       + number of particle signatures (all particles)
 *
 *  Good values for the first parameters are found by a grid search.
 *  E.g. a BDT output can be used as discriminating variable build from the
 *  input variables.
 */
// -- C++ STL headers.

// -- ROOT headers.

// -- LCIO headers.
#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Track.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCRelationImpl.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "UTIL/LCRelationNavigator.h"

// -- Marlin headers.
#include "marlin/Exceptions.h"

// -- Header for this processor and other project-specific headers.
#include "ref_util.h"
#include "tau_cones_bdt_input_processor.h"

// -- Using-declarations and global constants.
// Only in .cc files, never in .h header files!
using RP = EVENT::ReconstructedParticle;
using RPImpl = IMPL::ReconstructedParticleImpl;
using MCP = EVENT::MCParticle;
using Tlv = ROOT::Math::XYZTVector;

// This line allows to register your processor in marlin when calling
// "Marlin steering_file.xml".
TauConesBDTInputProcessor aTauConesBDTInputProcessor;
// ----------------------------------------------------------------------------

TauConesBDTInputProcessor::TauConesBDTInputProcessor() :
    marlin::Processor("TauConesBDTInputProcessor") {
  // _description comes from marlin::Processor.
  _description = "Find tau remnants with the tau/isolation cone approach. ";

  registerInputCollection(
    LCIO::RECONSTRUCTEDPARTICLE,
    "PfoCollection",
    "The Pandora PFO collection name.",
    pfo_collection_name,
    std::string("PandoraPFOs"));

  registerProcessorParameter(
    "SearchConeAngle",
    "Opening angle of the search cone for the tau jet in rad.",
    cd.search_cone_angle,
    0.1f);

  registerProcessorParameter(
    "IsolationConeAngle",
    "Outer isolation cone around the search cone of the tau jet in rad.",
    cd.isolation_cone_angle,
    0.2f);

  registerProcessorParameter(
    "PtCut",
    "Cut on Pt to suppress background.",
    cd.p_t_cut,
    0.2f);

  registerProcessorParameter(
    "MaxCosTheta",
    "Cut on the cosine of theta to suppress background.",
    cd.max_cos_theta,
    0.99f);

  registerProcessorParameter(
    "MinPtSeed",
    "Minimum tranverse momentum of tau seed.",
    cd.min_p_t_seed,
    5.0f);

  registerProcessorParameter(
    "OutputRootFile",
    "Name of the output root file.",
    bdt_inputs_name,
    std::string("tau_cones_bdt_input"));
}
// ----------------------------------------------------------------------------

void TauConesBDTInputProcessor::init() {
  printParameters(); // Method from marlin::Processor.
  // This is also the place where a root file would be opened
  // and histograms would be defined.
  TString fnn(bdt_inputs_name.c_str());
  fnn += ".root";
  bdt_inputs = new TFile(fnn, "update");
  tau_parameters = new TNtuple("tau_parameters","tau_parameters",
      "n_charged:n_remnants:isolation_energy:n_in_isolation:m_invariant:"
      "is_tau");
}
// ----------------------------------------------------------------------------

void TauConesBDTInputProcessor::processEvent(EVENT::LCEvent* event) {
  streamlog_out(DEBUG) << "Processing event no " << event->getEventNumber()
    << " - run " << event->getEventNumber() << std::endl;
  EVENT::LCCollection* pfo_collection = nullptr;
  try {
    pfo_collection = event->getCollection(pfo_collection_name);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "RP collection " << pfo_collection_name
      << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }
  // Relation collection needed to classify the sample.
  UTIL::LCRelationNavigator* relation_navigator = nullptr;
  EVENT::LCCollection* relation_collection = nullptr;
  try {
    relation_collection = event->getCollection("RecoMCTruthLink");
    relation_navigator = new UTIL::LCRelationNavigator(relation_collection);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "The relation collection " << "RecoMCTruthLink"
      << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }

 // Create dedicated vectors for charged and neutral RPs. Discard particles
  // outside of the acceptance volume.
  EventVector rpv = EventVector(pfo_collection, cd);

 // Fill the vector of tau candidates (each candidate being a vector of RPs).
  FindAllTaus(rpv, cd.min_p_t_seed, cd.search_cone_angle);
  MergeCloseByTaus(rpv, cd.search_cone_angle);

 // Fill the bdt parameters.
  for (TauCandidate tau : rpv.taus) {
    TauParameters tau_bdt{};
   // Mass, tracks.
    tau_bdt.n_charged = tau.n_charged;
    tau_bdt.n_remnants = tau.parts.size();
    tau_bdt.m_invariant = tau.tlv.M();
   // Isolation.
    for (RP* particle : rpv.all) {
      Tlv particle_tlv = ref_util::getTlv(particle);
      double angle_to_tau = acos(ROOT::Math::VectorUtil::CosTheta(
          tau.tlv, particle_tlv));
      if (angle_to_tau > cd.search_cone_angle
      && angle_to_tau < cd.isolation_cone_angle) {
        ++tau_bdt.n_in_isolation;
        tau_bdt.isolation_energy += particle_tlv.E();
      }
    }
    bool is_tau = APartFromTau(tau, relation_navigator);
  tau_parameters->Fill(tau_bdt.n_charged, tau_bdt.n_remnants,
      tau_bdt.isolation_energy, tau_bdt.n_in_isolation, tau_bdt.m_invariant,
      is_tau);
  }
}
// ----------------------------------------------------------------------------

void TauConesBDTInputProcessor::end() {
  bdt_inputs->cd();
  tau_parameters->Write();
  bdt_inputs->Write(0);
  bdt_inputs->Close();
}
// ----------------------------------------------------------------------------

// Take as a tau any candidate with at least half its energy stemming from a
// Monte Carlo tau right out of the simulation (e.g. not H->tautau).
bool TauConesBDTInputProcessor::APartFromTau(TauCandidate tau_candidate,
    UTIL::LCRelationNavigator* relation_navigator) {
  int pdg = ref_util::getPrimaryPdg(tau_candidate.seed, relation_navigator);
  if (pdg == 15 || pdg == -15) {
    return true;
  }
  return false;
}