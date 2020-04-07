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
#include <cassert>

// -- ROOT headers.

// -- LCIO headers.
#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Track.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCRelationImpl.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include <UTIL/LCIterator.h>
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

  registerInputCollection(
	LCIO::MCPARTICLE,
    "MCParticleCollection",
    "This MCParticle collection should link to the specified ."
      "ReconstructedParticle collection.",
    mc_collection_name,
    std::string("MCParticlesSkimmed"));

  registerInputCollection(
	LCIO::LCRELATION,
	"RelationCollection",
    "Collection of relations that link the Monte Carlo particles with the "
    "reconstructed particles.",
	relation_collection_name,
    std::string("RecoMCTruthLink"));

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
  quality_control = new TNtuple("quality_control","quality_control",
      "tau_p_in_event:tau_m_in_event:taus_per_event:not_tau_in_event:"
      "total_in_event:charge_mixed_up");
}
// ----------------------------------------------------------------------------

void TauConesBDTInputProcessor::processRunHeader(EVENT::LCRunHeader* run) {
  streamlog_out(MESSAGE) << "Starting run no " << run->getRunNumber()
    << std::endl;
}
// ----------------------------------------------------------------------------

void TauConesBDTInputProcessor::processEvent(EVENT::LCEvent* event) {
  streamlog_out(DEBUG) << "Processing event no " << event->getEventNumber()
    << std::endl;
  EVENT::LCCollection* pfo_collection = nullptr;
  try {
    pfo_collection = event->getCollection(pfo_collection_name);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "RP collection " << pfo_collection_name
      << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }

 // Create dedicated vectors for charged and neutral RPs. Discard particles
  // outside of the acceptance volume.
  EventVector rpv = EventVector(pfo_collection, cd);

 // Fill the vector of tau candidates (each candidate being a vector of RPs).
  ////FindAllTaus(rpv, cd.min_p_t_seed, cd.search_cone_angle);
  ////MergeCloseByTaus(rpv, cd.search_cone_angle);

 // Fill the bdt parameters.
  int n_tau_p = 0, n_tau_m = 0, n_tau = 0, n_not_tau = 0, n_charge_mixed_up = 0;
  std::vector<RP*> tau_seeds_pm = TauPMSeed(event,
      mc_collection_name, relation_collection_name);
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
    bool is_tau = false;
    if (tau.seed == tau_seeds_pm[0] || tau.seed == tau_seeds_pm[1]) {
      is_tau = true;
      ++n_tau;
      int charge = tau.rp_impl->getCharge();
      if ((charge > 0 && tau.seed != tau_seeds_pm[0])
      || (charge < 0 && tau.seed != tau_seeds_pm[1])) {
        ++n_charge_mixed_up;
      } else {
        if (charge > 0) ++n_tau_p;
        if (charge < 0) ++n_tau_m;
      }
    } else {
      ++n_not_tau;
    }
    tau_parameters->Fill(tau_bdt.n_charged, tau_bdt.n_remnants,
        tau_bdt.isolation_energy, tau_bdt.n_in_isolation, tau_bdt.m_invariant,
        is_tau);
  }
  quality_control->Fill(n_tau_p, n_tau_m, n_tau, n_not_tau,
      n_tau + n_not_tau, n_charge_mixed_up);

  // Checks. Might be too strong: It is imaginable that the tau seed is
  // disguised by a nearby higher-pT background particle.
  if ((tau_seeds_pm[0] && tau_seeds_pm[1]) && (n_tau != 2)) {
    std::cout << "Two seeds are provided by TauPMSeed, but only "
      << n_tau << " result in a is_tau flag." << std::endl;
  }
  if ((tau_seeds_pm[0] || tau_seeds_pm[1]) && (n_tau <= 0)) {
    std::cout << "At least one seed is provided by TauPMSeed, but only "
      << n_tau << " result in a is_tau flag." << std::endl;
  }

  TrackCheck(event);
}
// ----------------------------------------------------------------------------

void TauConesBDTInputProcessor::end() {
  bdt_inputs->cd();
  tau_parameters->Write();
  quality_control->Write();
  bdt_inputs->Write(0);
  bdt_inputs->Close();
}
// ----------------------------------------------------------------------------

// Identify the should-be MCTruth seeds for the taus from the recoiling Z decay.
// Takes into account the parton mixing.
// A significant fraction of the MC taus will not be able to find a seed.
// We called ref_util::printFamilyTree to validate that this (mainly) is due to
// taus where the large majority of the energy is carried by neutral particles
// (oftentimes especially neutrinos), or (seldomly) a seed with almost no pT.
std::vector<RP*> TauConesBDTInputProcessor::TauPMSeed(EVENT::LCEvent* event,
    std::string mc_col_name, std::string relation_col_name) {
  // Relation collection needed to classify the sample.
  UTIL::LCRelationNavigator* relation_navigator = nullptr;
  EVENT::LCCollection* relation_collection = nullptr;
  try {
    relation_collection = event->getCollection(relation_col_name);
    relation_navigator = new UTIL::LCRelationNavigator(relation_collection);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "The relation collection "
      << relation_col_name << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }
  // MonteCarlo collection that the ReconstructedParticles are related to.
  EVENT::LCCollection* mc_collection = nullptr;
  try {
    mc_collection = event->getCollection(mc_col_name);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "MC collection " << mc_col_name
      << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }
  ref_util::printFamilyTree(mc_collection);
 // Find the three (starting) MCParticles of interest: 25, 15, -15.
  MCP *tau_p = nullptr, *tau_m = nullptr, *higgs = nullptr;
  for(UTIL::LCIterator<MCP> it(mc_collection); MCP* mcp = it.next();) {
    if (mcp->getParents().size() != 0) continue;
    if (mcp->getPDG() == 15) {
      if (tau_m) streamlog_out(ERROR) << "There was already a starting tau_m "
        << "before this one in the same event. How?" << std::endl;
      tau_m = mcp;
    } else if (mcp->getPDG() == -15) {
      if (tau_p) streamlog_out(ERROR) << "There was already a starting tau_p "
        << "before this one in the same event. How?" << std::endl;
      tau_p = mcp;
    } else if (mcp->getPDG() == 25) {
      if (higgs) streamlog_out(ERROR) << "There was already a starting higgs "
        << "before this one in the same event. How?" << std::endl;
      higgs = mcp;
    }  // Other starting MCParticles (e+/e-, quarks, photons) don't bother us.
  }

 // Events without a Higgs are not signal for us. Thus, any taus in them should
 // be treated as background.
  if (not higgs) return std::vector<RP*>{nullptr, nullptr};
  if (not tau_p || not tau_m) return std::vector<RP*>{nullptr, nullptr};


 // Get the RP seeds for tau_p / tau_m.
  // There can be parton level mixing between the tau_p and the tau_m (~1/3 of
  // the events). Account for this by following the decay chain up to
  // the last occurrance of a tau_p / tau_m.
  while (tau_p->getDaughters().size()) {
    MCP* tau_d = tau_p->getDaughters()[0];
    int d_pdg = tau_d->getPDG();
    if (d_pdg == -15) {
      tau_p = tau_d;
    } else if (d_pdg == 94) {
      std::vector<MCP*> pm_mcs = tau_d->getDaughters();
      if (pm_mcs.size() != 2 || abs(pm_mcs[0]->getPDG()) != 15
      || abs(pm_mcs[1]->getPDG()) != 15) {
        streamlog_out(ERROR) << "We would expect to get a tau+/- pair after the"
          <<" parton mixing. What happend here?" << std::endl;
      } else if (pm_mcs[0]->getPDG() == -15) {
        tau_p = pm_mcs[0];
        tau_m = pm_mcs[1];
      } else {
        tau_p = pm_mcs[1];
        tau_m = pm_mcs[0];
      }
    } else if (d_pdg == -16) {
      // Land here if we start going down a proper tau decay chain.
      break;
    } else {
      streamlog_out(ERROR) << "Unexpected tau daughter of PDG " << d_pdg << "."
        << std::endl;
    }
  }

  std::vector<RP*> tau_seeds_pm{};
  for (MCP* tau_mc : std::vector<MCP*>{tau_p, tau_m}) {
    std::vector<MCP*> remnants{tau_mc};
    RP* seed = nullptr;
    double seed_pt = 0;
    while (remnants.size()) {
      MCP* rem = remnants.back();
      remnants.pop_back();
      for (auto tau_rem : relation_navigator->getRelatedFromObjects(rem)) {
        RP* seed_candidate = static_cast<RP*>(tau_rem);
        if (seed_candidate->getCharge() == 0) continue;
        double pt = sqrt(pow(seed_candidate->getMomentum()[0], 2)
                       + pow(seed_candidate->getMomentum()[1], 2));
        if (pt > seed_pt) {
          // Still "veto" a seed if it has no chance to be found (does not pass
          // our criteria for seeds / RPs in general).
          if (pt < fmax(cd.min_p_t_seed, cd.p_t_cut)) continue;
          double pz = fabs(seed_candidate->getMomentum()[2]);
          if (pz / (pt*pt+pz*pz) > 0.99) continue;

          seed = seed_candidate;
        }
      }
      for (MCP* d : rem->getDaughters()) {
        remnants.push_back(d);
      }
    }
    tau_seeds_pm.push_back(seed);
  }
  if (not tau_seeds_pm[0] || not tau_seeds_pm[1]) {
  streamlog_out(MESSAGE) << "No seed RP for tau_p or tau_m could be "
    << "identified with the 'TauPMSeed' algorithm." << std::endl;
  ////ref_util::printFamilyTree(mc_collection);
  }
  if (not tau_seeds_pm[0] && not tau_seeds_pm[1]) {
  streamlog_out(MESSAGE) << "This is the fault of both of them." << std::endl;
  }

  return tau_seeds_pm;

}
// ----------------------------------------------------------------------------

// Use the tracks directly.
void TauConesBDTInputProcessor::TrackCheck(EVENT::LCEvent* event) {
  std::string track_collection_name = "MarlinTrkTracks";
  EVENT::LCCollection* track_collection = nullptr;
  try {
    track_collection = event->getCollection(track_collection_name);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "RP collection " << track_collection_name
      << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }
  for (int e = 0; e < track_collection->getNumberOfElements(); ++e) {
    EVENT::Track* track = static_cast<EVENT::Track*>(
      track_collection->getElementAt(e));
  }

}