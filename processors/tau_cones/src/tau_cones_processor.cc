/**
 *    @author A. Muennich, CERN (original TauFinder).
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP (adaptation).
*/
// -- C++ STL headers.
#include <cassert>

// -- ROOT headers.

// -- LCIO headers.
#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
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
using RP = EVENT::ReconstructedParticle;
using RPImpl = IMPL::ReconstructedParticleImpl;
using MCP = EVENT::MCParticle;
using Tlv = ROOT::Math::XYZTVector;
const int kMaxChargedTracks = 4;
const int kMaxOverallTracks = 10; // TODO: maybe those two as processor parameters?
const int kPrintInfoOfEventNumber = 5;
const bool kAddUpCollectionNumbers = true;

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
    search_cone_angle_,
    0.05f);

  registerProcessorParameter(
    "IsolationConeAngle",
    "Outer isolation cone around the search cone of the tau jet in rad.",
    isolation_cone_angle_,
    0.07f);

  registerProcessorParameter(
    "IsolationEnergy",
    "Maximal allowed energy within the isolation cone region.",
    max_isolation_energy_,
    5.0f);

  registerProcessorParameter(
    "PtCut",
    "Cut on Pt to suppress background.",
    p_t_cut_,
    0.2f);

  registerProcessorParameter(
    "MaxCosTheta",
    "Cut on the cosine of theta to suppress background.",
    max_cos_theta_,
    0.99f);

  registerProcessorParameter(
    "MinPtSeed",
    "Minimum tranverse momentum of tau seed.",
    min_p_t_seed_,
    5.0f);

  registerProcessorParameter(
    "MaxInvariantMass",
    "Upper limit on the invariant mass of the tau candidate.",
    max_m_invariant_,
    2.0f);
}

// ----------------------------------------------------------------------------
void TauConesProcessor::init() {
  // Usually a good idea to print parameters.
  printParameters(); // method from marlin::Processor.
  // This is also the place where a root file would be opened
  // and histograms would be defined.
  TString fnn(out_root_filename_.c_str());
  fnn += ".root";
  root_out_ = new TFile(fnn, "recreate");
  fail_reason_tuple_ = new TNtuple("fail_reason_tuple_","fail_reason_tuple_",
      ":m_invariant_too_high:m_invariant_negative: wrong_track_number"
        ":not_isolated:tried_to_merge:taus_identified");
  n_m_invariant_too_high_ = 0;
  n_m_invariant_negative_ = 0;
  n_wrong_track_number_ = 0;
  n_not_isolated_ = 0;
  n_tried_to_merge_ = 0;
  n_taus_identified_ = 0;
  n_events_total_ = 0;
  n_runs_ = 0;
}

// ----------------------------------------------------------------------------
void TauConesProcessor::processRunHeader(EVENT::LCRunHeader* run) {
  streamlog_out(MESSAGE) << "Starting run no " << run->getRunNumber()
    << std::endl;
  ++n_runs_;
}

// ----------------------------------------------------------------------------
void TauConesProcessor::processEvent(EVENT::LCEvent* event) {
  streamlog_out(DEBUG) << "Processing event no " << event->getEventNumber()
    << " - run " << event->getEventNumber() << std::endl;
  ++n_events_total_;
  // Prepare the input and output collections.
  EVENT::LCCollection* pfo_collection = nullptr;
  try {
    pfo_collection = event->getCollection(pfo_collection_name_);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "RP collection" << pfo_collection_name_
      << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }
  IMPL::LCCollectionVec* tau_collection
      = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  IMPL::LCCollectionVec* rest_collection
      = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  rest_collection->setSubset(1); // Is subset of pfo_collection.
  IMPL::LCCollectionVec* tau_relation_collection
      = new IMPL::LCCollectionVec(LCIO::LCRELATION);
  tau_relation_collection->parameters().setValue(
      std::string("FromType"), LCIO::RECONSTRUCTEDPARTICLE);
  tau_relation_collection->parameters().setValue(
      std::string("ToType"), LCIO::RECONSTRUCTEDPARTICLE);

  // Create dedicated vectors for charged and neutral RPs. Discard particles
  // outside of the acceptance volume.
  std::vector<RP*> all_rps;
  std::vector<RP*> charged_rps;
  std::vector<RP*> neutral_rps;
  for (int e = 0; e < pfo_collection->getNumberOfElements(); ++e) {
    RP* particle = static_cast<RP*>(pfo_collection->getElementAt(e));
    Tlv tlv = ref_util::getTlv(particle);
	if ((tlv.Pt() < p_t_cut_) && (cos(tlv.Theta()) > max_cos_theta_)) {
        all_rps.push_back(particle);
        int charge = particle->getCharge();
        if (charge == 0) {
          neutral_rps.push_back(particle);
        } else if (fabs(charge) == 1) {
          charged_rps.push_back(particle);
        } else {
          streamlog_out(ERROR) << "The reconstructed particle has charge="
            << charge << ". How that? Only charges 0, +-1 should appear."
            << std::endl;
          throw marlin::StopProcessingException(this);
        }
    }
  }
  // The algorithm depends on the order of the charged RPs (the possible seeds).
  // To have a reliable result, order them by energy.
  std::sort(charged_rps.begin(), charged_rps.end(), ref_util::rpEnergySort);

  // Fill the vector of tau candidates (each candidate being a vector of RPs).
  std::vector<std::vector<RP*>> taus;
  bool found_a_tau = true;
  while(charged_rps.size() && found_a_tau) {
    found_a_tau = FindTau(charged_rps, neutral_rps, taus);
  }
  for (RP* particle: charged_rps) rest_collection->addElement(particle);
  for (RP* particle: neutral_rps) rest_collection->addElement(particle);

  // Combine the remants of each tau back together into one particle.
  // While at it, we will reject some candidates based on track number, mass.
  std::vector<RPImpl*> newly_formed_taus;
  IntVec n_charged_tracks_per_tau;
  for (std::vector<RP*> tau_parts: taus) {
    RPImpl* tau
        = new RPImpl();
    int n_charged_tracks{0};
    int charge{0};
    Tlv tau_tlv = Tlv(0, 0, 0, 0);
    for (RP* part: tau_parts) {
      Tlv part_tlv = ref_util::getTlv(part);
      tau_tlv += part_tlv;
      int partial_charge = part->getCharge();
      charge += partial_charge;
      if (partial_charge) {
        ++n_charged_tracks;
      }
      tau->addParticle(part);
    }
    if ((tau_tlv.M() > max_m_invariant_) || (tau_tlv.M() < 0)
    || (n_charged_tracks > kMaxChargedTracks) || (n_charged_tracks == 0)) {
      // Reject these events. Add the rejection reason, move the RPs into the
      // rest collection and maybe print some information.
      if (tau_tlv.M() > max_m_invariant_) {
          ++n_m_invariant_too_high_;
      } else if (tau_tlv.M() < 0) {
        ++n_m_invariant_negative_;
      } else {
        ++n_wrong_track_number_;
      }
      for (RP* part: tau_parts) {
        rest_collection->addElement(part);
      }
      if (n_events_total_== kPrintInfoOfEventNumber) {
        streamlog_out(DEBUG) << "A tau candidate failed:"
          << "  inv.Mass=" << tau_tlv.M() << ",  Pt=" << tau_tlv.Pt()
          << ",  E="<<tau_tlv.E() << ",  charged tracks: "<< n_charged_tracks
          << ",  particles: " << tau->getParticles().size()
          << ",  phi=" << tau_tlv.Phi() << ",  theta=" << tau_tlv.Theta()
          << std::endl;
	  }
    } else {
      // These RPs are accepted as true taus.
      int tau_pdg = 15;
      if (charge < 0) tau_pdg = -15;
      // TODO: Maybe replace the 2 lines above by int tau_pdg = charge*15. For this we would first have to investigate which charge values are actually taken by our constructed taus.
      tau->setType(tau_pdg);
      tau->setEnergy(tau_tlv.E());
      double momentum[3];
      tau_tlv.GetCoordinates(momentum);
      tau->setMomentum(momentum);
      newly_formed_taus.push_back(tau);
      n_charged_tracks_per_tau.push_back(n_charged_tracks);
    }
  }
  // It turns out removing the merger parts right away complicates things
  // (keeping track of the positions in the track-number-counting vectors as
  // well as the RP vector itself). It is therefore performed in an additional
  // step later on.
  std::vector<std::vector<RPImpl*>::iterator> reject_as_taus;

  // Sometimes the FindTau algorithm splits a tau. Improve the performance by
  // merging taus that are very close together.
  if (newly_formed_taus.size() > 1) {
    for (std::vector<RPImpl*>::iterator iter_tau1 = newly_formed_taus.begin();
    // Two conditions: did not reach end yet and the tau was not already merged.
    iter_tau1 < newly_formed_taus.end() && (std::find(reject_as_taus.begin(),
        reject_as_taus.end(), iter_tau1) != reject_as_taus.end());
    ++iter_tau1) {
      for (std::vector<RPImpl*>::iterator iter_tau2 = iter_tau1 + 1;
      iter_tau2 < newly_formed_taus.end(); ++iter_tau2) {
        RPImpl* tau1 = *iter_tau1;
        Tlv tau1_tlv = ref_util::getTlv(tau1);
        RPImpl* tau2 = *iter_tau2;
        Tlv tau2_tlv = ref_util::getTlv(tau2);
        double angle_between_taus = acos(ROOT::Math::VectorUtil::CosTheta(
        tau1_tlv, tau2_tlv));
        // The angular criterion for merging: The reconstructed tracks of both
        // taus are closer together than the search_cone_angle_.
        if (angle_between_taus > search_cone_angle_) continue;
        // Now let's merge.
        ++n_tried_to_merge_;
        tau1_tlv += tau2_tlv;
        tau1->setEnergy(tau1_tlv.E());
        double momentum[3];
        tau1_tlv.GetCoordinates(momentum);
        tau1->setMomentum(momentum);
		tau1->setCharge(tau1->getCharge()+tau2->getCharge());
        reject_as_taus.push_back(iter_tau2);
        n_charged_tracks_per_tau[iter_tau1-newly_formed_taus.begin()]
            += n_charged_tracks_per_tau[iter_tau2-newly_formed_taus.begin()];
        if (n_events_total_== kPrintInfoOfEventNumber) {
          int n_charged_tracks_tau2 = n_charged_tracks_per_tau[
              iter_tau2-newly_formed_taus.begin()];
          int n_charged_tracks_tau1_only = n_charged_tracks_per_tau[
              iter_tau1-newly_formed_taus.begin()] - n_charged_tracks_tau2;
          int n_neutral_particles_tau1_only = tau1->getParticles().size()
              - tau2->getParticles().size();
          streamlog_out(DEBUG) << "Tau Merging: "
            << iter_tau1 - newly_formed_taus.begin() << " with "
            << iter_tau2 - newly_formed_taus.begin() << "." << std::endl;
          streamlog_out(DEBUG) << "Angle was: " << angle_between_taus
            << ". Combined kinematics: " << "E: " << tau1_tlv.E()
            << ", theta: " << tau1_tlv.Theta() << ", phi: " << tau1_tlv.Phi()
            << ", charged tracks: " << n_charged_tracks_tau1_only << "+"
            << n_charged_tracks_tau2 << ", particles: "
            << tau1->getParticles().size() << "+" << tau2->getParticles().size()
            << std::endl;
        }
        // TODO: Simlar to merging conditions above. Consider writing a merge-function.
        int charged_tracks_in_merged = n_charged_tracks_per_tau[
            iter_tau1-newly_formed_taus.begin()];
        if ((tau1_tlv.M() > max_m_invariant_) || (tau1_tlv.M() < 0)
        ||  charged_tracks_in_merged > kMaxChargedTracks) {
          // After merging, it is not a proper tau candidate any more.
          reject_as_taus.push_back(iter_tau1);
          if (tau1_tlv.M() > max_m_invariant_) {
              ++n_m_invariant_too_high_;
          } else if (tau1_tlv.M() < 0) {
            ++n_m_invariant_negative_;
          } else {
            ++n_wrong_track_number_;
          }
        } else {  // The tau conditions are satisfied for the merged tau.
          for (RP* tau2_part: tau2->getParticles()) {
            tau1->addParticle(tau2_part);
            ++n_taus_identified_;
          }
        }
      }
    }
  }

  // A (final) candidate rejection based on Isolation (cone) and number of
  // remnants.
  for (std::vector<RPImpl*>::iterator iter_tau = newly_formed_taus.begin();
  // Two conditions: did not reach end yet and the tau was not already merged.
  iter_tau < newly_formed_taus.end() && (std::find(reject_as_taus.begin(),
      reject_as_taus.end(), iter_tau) != reject_as_taus.end());
  ++iter_tau) {
    RPImpl* tau = *iter_tau;
    Tlv tau_tlv = ref_util::getTlv(tau);
    float isolation_energy{0.f};
    int n_isolation_particles{0};
    int n_charged_tracks = n_charged_tracks_per_tau[
        iter_tau-newly_formed_taus.begin()];
    // Make sure the number of particles is realistic for a tau decay.
    if ((tau->getParticles().size() > kMaxOverallTracks)
    ||  (n_charged_tracks > kMaxChargedTracks)) {
      ++n_wrong_track_number_;
      reject_as_taus.push_back(iter_tau);
      continue; // Do not have to check for isolation as we already rejected.
    }
    // Now check for isolation.
    for (RP* particle : all_rps) {
        Tlv particle_tlv = ref_util::getTlv(particle);
        double angle_to_tau = acos(ROOT::Math::VectorUtil::CosTheta(
            tau_tlv, particle_tlv));
        if (angle_to_tau > search_cone_angle_
        && angle_to_tau < isolation_cone_angle_) {
          ++n_isolation_particles;
          isolation_energy += particle_tlv.E();
        }
    }
    if (isolation_energy > max_isolation_energy_) {
      ++n_not_isolated_;
      reject_as_taus.push_back(iter_tau);
      if (n_events_total_== kPrintInfoOfEventNumber) {
        streamlog_out(DEBUG) << "Isolation failed:  tau E=" << tau_tlv.E()
          << ", isolation E=" << isolation_energy << " from "
          << n_isolation_particles << "particles." << std::endl;
      }
    }
  }

  // Put the  tau candidate particles in the right collections, depending on
  // wether they ended up being rejected or not.
  for (std::vector<RPImpl*>::iterator iter_tau = newly_formed_taus.begin();
  iter_tau < newly_formed_taus.end(); ++iter_tau) {
    RPImpl* tau = *iter_tau;
    if (std::find(reject_as_taus.begin(), reject_as_taus.end(), iter_tau)
    == reject_as_taus.end()) {
      // These are finally our tau candidates.
      tau_collection->addElement(tau);
      for (RP* tau_remnant: tau->getParticles()) {
	    IMPL::LCRelationImpl* relation = new LCRelationImpl(tau, tau_remnant);
	    tau_relation_collection->addElement(relation);
      }
      if (n_events_total_== kPrintInfoOfEventNumber) {
        Tlv tau_tlv = ref_util::getTlv(tau);
        streamlog_out(DEBUG) << "Identified a tau: E=" << tau_tlv.E()
          << ", charged tracks: "
          << n_charged_tracks_per_tau[iter_tau-newly_formed_taus.begin()]
          << ", particles: " << tau->getParticles().size() << std::endl;
      }
    } else {
      for (RP* not_tau_part: tau->getParticles()) {
        rest_collection->addElement(not_tau_part);
      }
    }
  }

  // Wrapping up. Check that all particles found their place in the collections.
  if (kAddUpCollectionNumbers) {
    assert(pfo_collection->getNumberOfElements()
        == rest_collection->getNumberOfElements()
          + tau_relation_collection->getNumberOfElements());
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
  event->addCollection(tau_collection, tau_collection_name_);
  event->addCollection(rest_collection, rest_collection_name_);
  event->addCollection(tau_relation_collection, tau_relation_collection_name_);
}

// ----------------------------------------------------------------------------
void TauConesProcessor::end() {
  // Fill the root file.
  fail_reason_tuple_->Fill(n_m_invariant_too_high_, n_m_invariant_negative_,
      n_wrong_track_number_, n_not_isolated_, n_tried_to_merge_,
      n_taus_identified_);
  root_out_->cd();
  fail_reason_tuple_->Write();
  root_out_->Write(0);
  root_out_->Close();
  // Print the tuple's information.
  streamlog_out(MESSAGE) << "end()  " << name() << " processed "
    << n_events_total_ << " events in " << n_runs_ << " runs." << std::endl;
  streamlog_out(MESSAGE) << "Number of found tau candidates:  "
    << n_taus_identified_ << std::endl;
  streamlog_out(MESSAGE) << "Reasons for failure:  " <<std::endl;
  streamlog_out(MESSAGE) << "  High invariant mass:  "
    << n_m_invariant_too_high_ << std::endl;
  streamlog_out(MESSAGE) << "  Negative inv. mass :  "
    << n_m_invariant_negative_ << std::endl;
  streamlog_out(MESSAGE) << "  No/too many tracks :  "
    << n_wrong_track_number_ << std::endl;
  streamlog_out(MESSAGE) << "  No isolation       :  "
    << n_not_isolated_ << std::endl;
  streamlog_out(MESSAGE) << "  Tried to merge     :  "
    << n_tried_to_merge_ << std::endl;
}

// ----------------------------------------------------------------------------
bool TauConesProcessor::FindTau(std::vector<RP*> &charged_rps,
    std::vector<RP*> &neutral_rps, std::vector<std::vector<RP*> > &taus) {
  std::vector<RP*>  tau;
  if (charged_rps.size() == 0) return false;  // Need charged p to seed the tau.
  double max_opening_angle = 0;
  // Look for a seed.
  RP* tau_seed = nullptr;
  std::vector<RP*>::iterator seed_rp_ptr = charged_rps.begin();
  while (seed_rp_ptr != charged_rps.end()) {
    tau_seed = *seed_rp_ptr;
    double p_t = sqrt(pow(tau_seed->getMomentum()[0], 2)
                    + pow(tau_seed->getMomentum()[1], 2));
    if (p_t > min_p_t_seed_) break; // We found our tau_seed.
	++seed_rp_ptr;
	tau_seed = nullptr; // This particle was no good seed. Try a new one.
  }
  if (!tau_seed) return false; // No particle satisfies the seed conditions.
  tau.push_back(tau_seed);
  charged_rps.erase(seed_rp_ptr);
  Tlv seed_tlv = ref_util::getTlv(tau_seed);
  if (n_events_total_== kPrintInfoOfEventNumber) {  // Just for printing info.
    streamlog_out(DEBUG) << "Seeding: " << tau_seed->getType() << "\t"
      << seed_tlv.E() << "\t" << seed_tlv.P() << "\t"
      << seed_tlv.Theta() << "\t" << seed_tlv.Phi() << std::endl;
  }

  // Assign the charged particles to the tau candidate.
  std::vector<RP*>::iterator charged_rp_ptr=charged_rps.begin();
  while (charged_rp_ptr != charged_rps.end()) {
    RP* track = *charged_rp_ptr;
    Tlv charged_tlv = ref_util::getTlv(track);
    double angle_to_seed = acos(ROOT::Math::VectorUtil::CosTheta(
        seed_tlv, charged_tlv));
    // The only condition for being part of the tau: Be inside the search cone.
    if (angle_to_seed < search_cone_angle_) {
      tau.push_back(track);
      charged_rps.erase(charged_rp_ptr);
      if (n_events_total_== kPrintInfoOfEventNumber) {// Just for printing info.
        streamlog_out(DEBUG) << "Adding charged: " << track->getType() << "\t"
          << charged_tlv.E() << "\t" << charged_tlv.P() << "\t"
          << charged_tlv.Theta() << "\t" << charged_tlv.Phi() << std::endl;
	  }
      if (angle_to_seed > max_opening_angle) max_opening_angle = angle_to_seed;
	  // Since we erased the element at position charged_rp_ptr, the pointer
      // already points to a new element. We should not invoke ++charged_rp_ptr.
	} else {  // This particle was to far away from the seed. Try the next.
      ++charged_rp_ptr;
    }
  }

  // Assign the neutral particles to the tau candidate.
  std::vector<RP*>::iterator neutral_rp_ptr=neutral_rps.begin();
  while (neutral_rp_ptr != neutral_rps.end()) {
    RP* deposit = *neutral_rp_ptr;  // Track would be a lie for a neutral p.
    Tlv neutral_tlv = ref_util::getTlv(deposit);
    double angle_to_seed = acos(ROOT::Math::VectorUtil::CosTheta(
        seed_tlv, neutral_tlv));
    // The only condition for being part of the tau: Be inside the search cone.
    if (angle_to_seed < search_cone_angle_) {
      tau.push_back(deposit);
      neutral_rps.erase(neutral_rp_ptr);
      if (n_events_total_== kPrintInfoOfEventNumber) {// Just for printing info.
        streamlog_out(DEBUG) << "Adding neutral: " << deposit->getType() << "\t"
          << neutral_tlv.E() << "\t" << neutral_tlv.P() << "\t"
          << neutral_tlv.Theta() << "\t" << neutral_tlv.Phi() << std::endl;
	  }
      if (angle_to_seed > max_opening_angle) max_opening_angle = angle_to_seed;
	  // Since we erased the element at position neutral_rp_ptr, the pointer
      // already points to a new element. We should not invoke ++neutral_rp_ptr.
	} else {  // This particle was to far away from the seed. Try the next.
      ++neutral_rp_ptr;
    }
  }

  // We successfully reconstructed a tau.
  taus.push_back(tau);
  return true;
}