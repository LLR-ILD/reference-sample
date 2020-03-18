/**
 *    @author A. Muennich, CERN (original TauFinder).
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP (adaptation).
*/
// -- C++ STL headers.
#include <cassert>

// -- ROOT headers.

// -- LCIO headers.
#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCRelationImpl.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "UTIL/LCRelationNavigator.h"

// -- Marlin headers.
#include "marlin/Exceptions.h"

// -- Header for this processor and other project-specific headers.
#include "ref_util.h"
#include "tau_cones_processor.h"

// -- Using-declarations and global constants.
// Only in .cc files, never in .h header files!
using ColVec = IMPL::LCCollectionVec;
using RP = EVENT::ReconstructedParticle;
using RPImpl = IMPL::ReconstructedParticleImpl;
using Tlv = ROOT::Math::XYZTVector;
const int kPrintInfoOfEventNumber = 44;
const bool kAssertCollectionNumber = true;

// This line allows to register your processor in marlin when calling
// "Marlin steering_file.xml".
TauConesProcessor aTauConesProcessor;
// ----------------------------------------------------------------------------

TauConesProcessor::TauConesProcessor() :
    marlin::Processor("TauConesProcessor") {
  // _description comes from marlin::Processor.
  _description = "Find tau remnants with the tau/isolation cone approach. ";

  registerInputCollection(
    LCIO::RECONSTRUCTEDPARTICLE,
    "PfoCollection",
    "The Pandora PFO collection name.",
    pfo_collection_name_,
    std::string("PandoraPFOs"));

  registerOutputCollection(
    LCIO::RECONSTRUCTEDPARTICLE,
	"TauCollection",
    "Name of the new PFO collection of tau candidates.",
    tau_collection_name_,
    std::string("TauPFOs"));

  registerOutputCollection(
    LCIO::RECONSTRUCTEDPARTICLE,
	"TauRestCollection",
	"Name of the new PFO collection containing all those particles that are "
      "not part of a tau candidate.",
    rest_collection_name_,
    std::string("TauRestPFOs"));

  registerOutputCollection(
    LCIO::LCRELATION,
	"TauRelationCollection",
	"Name of the relation collection linking the tau candidates to their "
    "counterparts in the ReconstructedParticle collection." ,
	tau_relation_collection_name_,
	std::string("TauRelation"));

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
    "IsolationEnergy",
    "Maximal allowed energy within the isolation cone region.",
    tau_cut.isolation_energy,
    2.0f);

  registerProcessorParameter(
    "MaxIsoParticles",
    "Upper limit on the number of particles in the isolation cone (without tau "
    "cone).",
    tau_cut.n_in_isolation,
    100);

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
    "MaxInvariantMass",
    "Upper limit on the invariant mass of the tau candidate.",
    tau_cut.m_invariant,
    2.0f);

  registerProcessorParameter(
    "MaxChargedRemnants",
    "Upper limit on the number of charged tracks that can belong to a tau "
    "candidate.",
    tau_cut.n_charged,
    4);

  registerProcessorParameter(
    "MaxRemnants",
    "Upper limit on the number of remanants of a tau candidate.",
    tau_cut.n_remnants,
    10);

  registerProcessorParameter(
    "OutputRootFile",
    "Name of the output root file.",
    out_root_filename_,
    std::string("tau_cones"));
}
// ----------------------------------------------------------------------------

void TauConesProcessor::init() {
  printParameters();
  root_out_ = new TFile((out_root_filename_+".root").c_str(), "update");
  fail_reason_tuple_ = new TNtuple("fail_reason_tuple_","fail_reason_tuple_",
      "m_invariant_too_high:m_invariant_negative:wrong_track_number"
        ":not_isolated:tried_to_merge:taus_identified"
        ":background_suppressed");
  n_events_total = 0;
}
// ----------------------------------------------------------------------------

void TauConesProcessor::processEvent(EVENT::LCEvent* event) {
  streamlog_out(DEBUG) << "Processing event no " << event->getEventNumber()
    << " - run " << event->getEventNumber() << std::endl;
  ++n_events_total;
  CandidateCounts event_count;

 // Prepare the input and output collections.
  EVENT::LCCollection* pfo_collection = nullptr;
  try {
    pfo_collection = event->getCollection(pfo_collection_name_);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "RP collection " << pfo_collection_name_
      << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }
  ColVec* tau_collection = new ColVec(LCIO::RECONSTRUCTEDPARTICLE);
  ColVec* rest_collection = new ColVec(LCIO::RECONSTRUCTEDPARTICLE);
  rest_collection->setSubset(1); // Is subset of pfo_collection.
  ColVec* tau_relation_collection = new ColVec(LCIO::LCRELATION);
  tau_relation_collection->parameters().setValue(
      std::string("FromType"), LCIO::RECONSTRUCTEDPARTICLE);
  tau_relation_collection->parameters().setValue(
      std::string("ToType"), LCIO::RECONSTRUCTEDPARTICLE);

 // Create dedicated vectors for charged and neutral RPs. Discard particles
  // outside of the acceptance volume.
  EventVector rpv = EventVector(pfo_collection, cd);
  event_count.n_background_suppressed = rpv.n_background_suppressed;

 // Fill the vector of tau candidates (each candidate being a vector of RPs).
  FindAllTaus(rpv, cd.min_p_t_seed, cd.search_cone_angle,
      kPrintInfoOfEventNumber == n_events_total); // Verbose output for one ev.
  for (RP* particle: rpv.all) rest_collection->addElement(particle);

 // Refine the candidate selection.
  MassTracksRejection(rpv, event_count, rest_collection);

  int n_merged = MergeCloseByTaus(rpv, cd.search_cone_angle,
      n_events_total== kPrintInfoOfEventNumber);
  event_count.n_tried_to_merge += n_merged;

  int add_reject = MassTracksRejection(rpv, event_count, rest_collection);
  if (n_events_total== kPrintInfoOfEventNumber and add_reject > 0) {
    streamlog_out(MESSAGE) << "Some merged object does not qualify as a "
      "tau anymore." << std::endl;
    }

  Isolation(rpv, event_count, rest_collection);

  // What was not rejected up to now will be considered part of a tau. Fill
  // the two tau related collections.
  for (TauCandidate tau_struct : rpv.taus) {
    // These are finally our tau candidates.
    RPImpl* tau = tau_struct.rp_impl;
    ++event_count.n_taus_identified;
    tau_collection->addElement(tau);
    for (RP* tau_remnant: tau->getParticles()) {
	  IMPL::LCRelationImpl* relation = new LCRelationImpl(tau, tau_remnant);
	  tau_relation_collection->addElement(relation);
    }
    if (n_events_total== kPrintInfoOfEventNumber) {
      Tlv tau_tlv = ref_util::getTlv(tau);
      streamlog_out(MESSAGE) << "Identified a tau: E=" << tau_tlv.E()
        << ", charged tracks: "
        << tau_struct.n_charged
        << ", particles: " << tau->getParticles().size() << "." << std::endl;
    }
  }

  // Wrapping up. Check that all particles found their place in the collections.
  if (kAssertCollectionNumber) {
    streamlog_out(DEBUG) << "Collection member count: "
      << pfo_collection->getNumberOfElements() << " ?= "
      << rest_collection->getNumberOfElements()
          + tau_relation_collection->getNumberOfElements()
          + event_count.n_background_suppressed << std::endl << "  = "
      << rest_collection->getNumberOfElements() << " + "
      << tau_relation_collection->getNumberOfElements() << " + "
      << event_count.n_background_suppressed << std::endl;
    assert(pfo_collection->getNumberOfElements()
        == rest_collection->getNumberOfElements()
          + tau_relation_collection->getNumberOfElements()
          + event_count.n_background_suppressed);
    std::vector<RP*> relation;
    UTIL::LCRelationNavigator* relation_navigator
        = new UTIL::LCRelationNavigator(tau_relation_collection);
    int n_relations_from{0};
    for (int i = 0; i < tau_collection->getNumberOfElements(); ++i) {
      RP* tau = static_cast<RP*>(tau_collection->getElementAt(i));
      n_relations_from += relation_navigator->getRelatedToObjects(tau).size();
    }
    assert(n_relations_from == tau_relation_collection->getNumberOfElements());
  }
  fail_reason_tuple_->Fill(
      event_count.n_m_invariant_too_high,
      event_count.n_m_invariant_negative,
      event_count.n_wrong_track_number,
      event_count.n_not_isolated,
      event_count.n_tried_to_merge,
      event_count.n_taus_identified,
      event_count.n_background_suppressed);
  total_count += event_count;
  event->addCollection(tau_collection, tau_collection_name_);
  event->addCollection(rest_collection, rest_collection_name_);
  event->addCollection(tau_relation_collection, tau_relation_collection_name_);
}
// ----------------------------------------------------------------------------

void TauConesProcessor::end() {
  // Fill the root file.
  root_out_->cd();
  fail_reason_tuple_->Write();
  root_out_->Write(0);
  root_out_->Close();
  // Print the tuple's information.
  streamlog_out(MESSAGE) << "end().  " << name() << " processed "
    << n_events_total << " events." << std::endl
    << "Number of found tau candidates:  " << total_count.n_taus_identified << std::endl
    << "Reasons for failure:  " <<std::endl
    << "  High invariant mass:  " << total_count.n_m_invariant_too_high << std::endl
    << "  Negative inv. mass :  " << total_count.n_m_invariant_negative << std::endl
    << "  No/too many tracks :  " << total_count.n_wrong_track_number << std::endl
    << "  No isolation       :  " << total_count.n_not_isolated << std::endl
    << "  Tried to merge     :  " << total_count.n_tried_to_merge << std::endl
    << "  Particles that were not considered due to background suppression:  "
      << total_count.n_background_suppressed << std::endl;

 ////std::cout
 ////<< tau_cut.m_invariant  << std::endl
 ////<< tau_cut.isolation_energy  << std::endl
 ////<< tau_cut.n_charged << std::endl
 ////<< tau_cut.n_in_isolation << std::endl
 ////<< tau_cut.n_remnants << std::endl;
}
// ----------------------------------------------------------------------------

// Reject candidates based on their mass and their number of (charged/neutral)
// signatures. Returns the number of rejected particles. The rejection region
// is updated in the event_count struct.
int TauConesProcessor::MassTracksRejection(EventVector &rpv,
    CandidateCounts &event_count, ColVec* rest_collection) {
  int n_rejected = 0;
  std::vector<TauCandidate> tc_to_remove;
  for (auto iter = rpv.taus.begin(); iter != rpv.taus.end(); iter++) {
    TauCandidate tau_c = *iter;
    if ((tau_c.tlv.M() > tau_cut.m_invariant)
    || (tau_c.tlv.M() < 0)
    || (tau_c.rp_impl->getParticles().size() > size_t(tau_cut.n_remnants))
    || (tau_c.n_charged > tau_cut.n_charged)
    || (tau_c.n_charged == 0)
    ) {
      // Reject these events. Add the rejection reason, move the RPs into the
      // rest collection and maybe print some information.
      ++n_rejected;
      if (tau_c.tlv.M() > tau_cut.m_invariant) {
        ++event_count.n_m_invariant_too_high;
      } else if (tau_c.tlv.M() < 0) {
        ++event_count.n_m_invariant_negative;
      } else {
        ++event_count.n_wrong_track_number;
      }
      if (n_events_total== kPrintInfoOfEventNumber) {
        streamlog_out(MESSAGE) << "A tau candidate failed:"
          << "  inv.Mass=" << tau_c.tlv.M() << ",  Pt=" << tau_c.tlv.Pt()
          << ",  E="<<tau_c.tlv.E() << std::endl
          << "  Charged tracks: " << tau_c.n_charged
          << ",  particles: " << tau_c.rp_impl->getParticles().size()
          << ",  phi=" << tau_c.tlv.Phi() << ",  theta=" << tau_c.tlv.Theta()
          << std::endl << " (Note: The candidate will not be considered any "
            "further. This includes not being considered for the tau numbering "
            "in the merger information." << std::endl;
	  }
      for (RP* part: tau_c.parts) {
        rest_collection->addElement(part);
      }
      tc_to_remove.push_back(tau_c);
    }
  }
  for (auto tc : tc_to_remove) rpv.taus.remove(tc);
  return n_rejected;
}
// ----------------------------------------------------------------------------

// A (final) basic cut-based candidate rejection based on Isolation (cone) and
// number of remnants.
int TauConesProcessor::Isolation(EventVector &rpv,
    CandidateCounts &event_count, ColVec* rest_collection) {
  int n_rejected = 0;
  std::vector<TauCandidate> tc_to_remove;
  for (TauCandidate tau : rpv.taus) {
    float isolation_energy = 0.f;
    int n_isolation_particles = 0;
    for (RP* particle : rpv.all) {
      Tlv particle_tlv = ref_util::getTlv(particle);
      double angle_to_tau = acos(ROOT::Math::VectorUtil::CosTheta(
          tau.tlv, particle_tlv));
      if (angle_to_tau > cd.search_cone_angle
      && angle_to_tau < cd.isolation_cone_angle) {
        ++n_isolation_particles;
        isolation_energy += particle_tlv.E();
      }
    }
    if (isolation_energy > tau_cut.isolation_energy) {
      ++event_count.n_not_isolated;
      ++n_rejected;
      if (n_events_total== kPrintInfoOfEventNumber) {
        streamlog_out(MESSAGE) << "Isolation failed:  tau E=" << tau.tlv.E()
          << ", isolation E=" << isolation_energy << " from "
          << n_isolation_particles << " particles." << std::endl;
      }
    } else {
      continue;  // Exactly the non-rejected particles arrive in this case.
    }
    for (RP* part: tau.parts) {
      rest_collection->addElement(part);
    }
    tc_to_remove.push_back(tau);
  }
  for (auto tc : tc_to_remove) rpv.taus.remove(tc);
  return n_rejected;
}