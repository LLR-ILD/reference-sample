/**
 *    @author A. Muennich, CERN (original TauFinder. Here: FindTau).
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP (adaptation).
 *
 *  Collect those methods that are used in multiple modes of the TauCones
 *  processor.
 */
// -- C++ STL headers.

// -- ROOT headers.

// -- LCIO headers.
#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"

// -- Marlin headers.
#include "marlin/Exceptions.h"  // For streamlogout.

// -- Header for this processor and other project-specific headers.
#include "ref_util.h"
#include "tau_cones_util.h"

// -- Using-declarations and global constants.
// Only in .cc files, never in .h header files!
using RP = EVENT::ReconstructedParticle;
using RPImpl = IMPL::ReconstructedParticleImpl;
using Tlv = ROOT::Math::XYZTVector;
// ----------------------------------------------------------------------------

tau_cones_util::TauCandidate tau_cones_util::TauFromParts(
    std::vector<RP*> tau_parts, RP* seed) {
  TauCandidate tau;
  RPImpl* rp_impl = new RPImpl();  // SegFault if we do not provide a new obj.
  tau.seed = seed;
  tau.parts = tau_parts;
  int charge = 0;
  for (RP* part: tau.parts) {
    Tlv part_tlv = ref_util::getTlv(part);
    tau.tlv += part_tlv;
    int partial_charge = part->getCharge();
    charge += partial_charge;
    if (partial_charge) {
      ++tau.n_charged;
    }
    rp_impl->addParticle(part);
  }
  // Set the global parameters of the tau RP.
  int tau_pdg = 15;
  if (charge > 0) tau_pdg = -15;
  // TODO: Maybe replace the 2 lines above by int tau_pdg = charge*15. For this we would first have to investigate which charge values are actually taken by our constructed taus.
  rp_impl->setType(tau_pdg);
  rp_impl->setCharge(charge);
  rp_impl->setEnergy(tau.tlv.E());
  double momentum[3];
  tau.tlv.GetCoordinates(momentum);
  rp_impl->setMomentum(momentum);
  tau.rp_impl = rp_impl;
  return tau;
}
// ----------------------------------------------------------------------------
namespace {

struct ConeRPVectors {
  std::vector<RP*> charged_in_tau_cone = {};
  //std::vector<RP*> charged_in_iso_cone = {};  // As of now, this is a dummy.
  std::vector<RP*> neutral_in_tau_cone = {};
  std::vector<RP*> neutral_in_iso_cone = {};
  ConeRPVectors(RP* seed_rp) {
    charged_in_tau_cone = {seed_rp};
    neutral_in_tau_cone = {};
    neutral_in_iso_cone = {};
  }
};  // namespace.

// TODO: A cleaner approach for these values.
double strict_min_pt_seed = 5.0;

bool IsStrictPiGamma(ConeRPVectors cone_rps) {
  // Look at charged.
  double cone_energy = 0;
  if (cone_rps.charged_in_tau_cone.size() != 1) return false;
  bool pt_is_too_low = true;
  int n_charged_pions = 0;
  for (RP* charged_rp : cone_rps.charged_in_tau_cone) {
    cone_energy += charged_rp->getEnergy();
    double p_t = sqrt(pow(charged_rp->getMomentum()[0], 2)
                    + pow(charged_rp->getMomentum()[1], 2));
    if (p_t > strict_min_pt_seed) pt_is_too_low = false;
    if (abs(charged_rp->getType()) == 211) ++n_charged_pions;
  }
  if (pt_is_too_low) streamlog_out(DEBUG) << "    Seed pt too low " << std::endl;
  if (pt_is_too_low) return false;
  if (n_charged_pions != 1) streamlog_out(DEBUG) << "    Seed pdg " << cone_rps.charged_in_tau_cone[0]->getType() << std::endl;
  if (n_charged_pions != 1) return false;
  // Look at neutral.
  double neutral_tau_cone_energy = 0;
  for (RP* neutral_rp : cone_rps.neutral_in_tau_cone) {
    cone_energy += neutral_rp->getEnergy();
    neutral_tau_cone_energy += neutral_rp->getEnergy();
  }
  streamlog_out(DEBUG) << "    Neutral tau cone energy: " << neutral_tau_cone_energy << std::endl;
  //if (neutral_tau_cone_energy < 3) return false;
  if (neutral_tau_cone_energy < 1) return false;
  double neutral_iso_cone_energy = 0;
  for (RP* neutral_rp : cone_rps.neutral_in_iso_cone) {
    neutral_iso_cone_energy += neutral_rp->getEnergy();
  }
  streamlog_out(DEBUG) << "    Neutral iso cone energy: " << neutral_iso_cone_energy << std::endl;
  streamlog_out(DEBUG) << "    Tau cone energy: " << cone_energy << std::endl;
  if (neutral_iso_cone_energy > 2) return false;
  ////if (cone_energy < 25) return false;
  ////if (cone_energy > 40) return false;
  return true;  // Candidate passed all tests.
}

bool IsThreeChargedPions(ConeRPVectors cone_rps) {
  // Look at charged.
  double cone_energy = 0;
  if (cone_rps.charged_in_tau_cone.size() != 3) return false;
  bool pt_is_too_low = true;
  int n_charged_pions = 0;
  for (RP* charged_rp : cone_rps.charged_in_tau_cone) {
    cone_energy += charged_rp->getEnergy();
    double p_t = sqrt(pow(charged_rp->getMomentum()[0], 2)
                    + pow(charged_rp->getMomentum()[1], 2));
    if (p_t > strict_min_pt_seed) pt_is_too_low = false;
    if (abs(charged_rp->getType()) == 211) ++n_charged_pions;
  }
  if (pt_is_too_low) streamlog_out(DEBUG) << "Seed pt too low " << n_charged_pions << std::endl;
  if (pt_is_too_low) return false;
  streamlog_out(DEBUG) << "    Charged pions: " << n_charged_pions << std::endl;
  if (n_charged_pions < 3) return false;
  // Look at neutral.
  for (RP* neutral_rp : cone_rps.neutral_in_tau_cone) {
    cone_energy += neutral_rp->getEnergy();
  }
  double neutral_iso_cone_energy = 0;
  for (RP* neutral_rp : cone_rps.neutral_in_iso_cone) {
    neutral_iso_cone_energy += neutral_rp->getEnergy();
  }
  streamlog_out(DEBUG) << "    Neutral iso cone energy: " << neutral_iso_cone_energy << std::endl;
  streamlog_out(DEBUG) << "    Tau cone energy: " << cone_energy << std::endl;
  if (neutral_iso_cone_energy > 2) return false;
  ////if (cone_energy < 25) return false;
  ////if (cone_energy > 60) return false;
  return true;  // Candidate passed all tests.
}
}  // namespace


bool tau_cones_util::HasStrictTau(tau_cones_util::EventVector &rpv,
    tau_cones_util::CandidateDefinition cd, EVENT::LCEvent* event) {
  if (rpv.charged.size() == 0) return false;  // Need charged p to seed the tau.
 // Look for a seed.
  for (RP* tau_seed : rpv.charged) {
    bool not_strict_tau = false;
    ////if (abs(ref_util::getPrimaryPdg(tau_seed, event)) != 15) continue;// TODO: Remove
    Tlv seed_tlv = ref_util::getTlv(tau_seed);
    streamlog_out(DEBUG) << "  Tau seed found (pT = " << seed_tlv.Pt() << ")" << std::endl;
    // If the energy is that low, the particle can not fullfill the p_t
    // requirement for any of the cases. As the particles are sorted by energy,
    // the same applies for all following particles.
    if (seed_tlv.E() < strict_min_pt_seed) return false;
    if (seed_tlv.Pt() < strict_min_pt_seed) continue;  // The seed doesn't work.
   // Fill lists for the different cones.
   ConeRPVectors cone_rps(tau_seed);
   // Assign the charged particles to the tau candidate.
    for (RP* track : rpv.charged) {
      if (track  == tau_seed) continue;
      Tlv charged_tlv = ref_util::getTlv(track);
      double angle_to_seed = acos(ROOT::Math::VectorUtil::CosTheta(
          seed_tlv, charged_tlv));
      if (angle_to_seed < cd.search_cone_angle) {
        cone_rps.charged_in_tau_cone.push_back(track);
      // Finding a charged particle in the iso cone is a show stopper (for now).
      } else if (angle_to_seed < cd.isolation_cone_angle) {
        streamlog_out(DEBUG) << "    Charged track in iso cone! " << angle_to_seed << "  " << charged_tlv.Pt() << std::endl;
        not_strict_tau = true;
        break;
      }
    }
    streamlog_out(DEBUG) << "    Charged tracks: " << cone_rps.charged_in_tau_cone.size() << std::endl;
    if (not_strict_tau) continue;
   // And now the neutral particles. Same procedure as for charged.
    for (RP* neutr_rp : rpv.neutral) {
      Tlv neutral_tlv = ref_util::getTlv(neutr_rp);
      double angle_to_seed = acos(ROOT::Math::VectorUtil::CosTheta(
          seed_tlv, neutral_tlv));
      if (angle_to_seed < cd.search_cone_angle) {
        cone_rps.neutral_in_tau_cone.push_back(neutr_rp);
      // Finding a charged particle in the iso cone is a show stopper (for now).
      } else if (angle_to_seed < cd.isolation_cone_angle) {
        cone_rps.neutral_in_iso_cone.push_back(neutr_rp);
      }
    }
   // Having collected the relevant particles, let's check if the set of
   // particles satisfies one of the strict selections.
    if (IsStrictPiGamma(cone_rps) || IsThreeChargedPions(cone_rps)) {
      ////rpv.strict_seed = tau_seed;
      return true;
    }
  }
  return false;
}


bool tau_cones_util::FindATau(tau_cones_util::EventVector &rpv,
    double min_p_t_seed, double search_cone_angle,
    bool print_info) {
  std::vector<RP*>  tau_parts;
  if (rpv.charged.size() == 0) return false;  // Need charged p to seed the tau.
 // Look for a seed.
  RP* tau_seed = nullptr;
  std::vector<RP*>::iterator seed_rp_ptr = rpv.charged.begin();
  while (seed_rp_ptr != rpv.charged.end()) {
    tau_seed = *seed_rp_ptr;
    // double p_t = sqrt(pow(tau_seed->getMomentum()[0], 2) // TODO: Maybe use pT again
    //                 + pow(tau_seed->getMomentum()[1], 2));
    // if (p_t > min_p_t_seed) break; // We found our tau_seed.
    if (tau_seed->getEnergy() > min_p_t_seed) break;  // We found our tau_seed.
    ////if (p_t > min_p_t_seed && abs(tau_seed->getType()) == 211) break; // We found our tau_seed.
    ++seed_rp_ptr;
    tau_seed = nullptr; // This particle was no good seed. Try a new one.
  }
  if (!tau_seed) return false; // No particle satisfies the seed conditions.
  tau_parts.push_back(tau_seed);
  rpv.charged.erase(seed_rp_ptr);
  Tlv seed_tlv = ref_util::getTlv(tau_seed);
  if (print_info) {  // Just for printing info.
    streamlog_out(MESSAGE) << "Seeding: " << tau_seed->getType() << "\t"
      << seed_tlv.E() << "\t" << seed_tlv.P() << "\t"
      << seed_tlv.Theta() << "\t" << seed_tlv.Phi() << std::endl;
  }

 // Assign the charged particles to the tau candidate.
  std::vector<RP*>::iterator charged_rp_ptr = rpv.charged.begin();
  while (charged_rp_ptr != rpv.charged.end()) {
    RP* track = *charged_rp_ptr;
    Tlv charged_tlv = ref_util::getTlv(track);
    double angle_to_seed = acos(ROOT::Math::VectorUtil::CosTheta(
        seed_tlv, charged_tlv));
    // The only condition for being part of the tau: Be inside the search cone.
    if (angle_to_seed < search_cone_angle) {
      tau_parts.push_back(track);
      rpv.charged.erase(charged_rp_ptr);
      if (print_info) {// Just for printing info.
        streamlog_out(MESSAGE) << "  Adding charged: " << track->getType()
          << "\t" << charged_tlv.E() << "\t" << charged_tlv.P() << "\t"
          << charged_tlv.Theta() << "\t" << charged_tlv.Phi() << std::endl;
    }
    // Since we erased the element at position charged_rp_ptr, the pointer
    // already points to a new element. We should not invoke ++charged_rp_ptr.
  } else {  // This particle was too far away from the seed. Try the next.
      ++charged_rp_ptr;
    }
  }

 // And now the neutral particles. Same procedure as for charged.
  std::vector<RP*>::iterator neutral_rp_ptr = rpv.neutral.begin();
  while (neutral_rp_ptr != rpv.neutral.end()) {
    RP* neutr_rp = *neutral_rp_ptr;
    Tlv neutral_tlv = ref_util::getTlv(neutr_rp);
    double angle_to_seed = acos(ROOT::Math::VectorUtil::CosTheta(
        seed_tlv, neutral_tlv));
    // The only condition for being part of the tau: Be inside the search cone.
    if (angle_to_seed < search_cone_angle) {
      tau_parts.push_back(neutr_rp);
      rpv.neutral.erase(neutral_rp_ptr);
      if (print_info) {// Just for printing info.
        streamlog_out(MESSAGE) << "  Adding neutral: " << neutr_rp->getType()
          << "\t" << neutral_tlv.E() << "\t" << neutral_tlv.P() << "\t"
          << neutral_tlv.Theta() << "\t" << neutral_tlv.Phi() << std::endl;
    }
  } else {  // This particle was to far away from the seed. Try the next.
      ++neutral_rp_ptr;
    }
  }



 // We successfully reconstructed a tau.
  rpv.taus.push_front(tau_cones_util::TauFromParts(tau_parts, tau_seed));
  return true;
}
// ----------------------------------------------------------------------------

bool tau_cones_util::FindAllTaus(tau_cones_util::EventVector &rpv,
    double min_p_t_seed, double search_cone_angle,
    bool print_info) {
  bool found_a_tau = false;
  bool look_for_more_taus = true;
  while(rpv.charged.size() && look_for_more_taus) {
    found_a_tau = look_for_more_taus;
    look_for_more_taus = FindATau(rpv,
        min_p_t_seed, search_cone_angle, print_info);
  }
  return found_a_tau;
}
// ----------------------------------------------------------------------------

int tau_cones_util::MergeCloseByTaus(tau_cones_util::EventVector &rpv,
    double search_cone_angle,
    bool print_info) {
  int n_merged = 0;
  std::vector<TauCandidate> tc_to_remove;
  if (! rpv.taus.empty()) {
    for (auto iter1 = rpv.taus.begin(); iter1 != rpv.taus.end(); iter1++) {
      TauCandidate tau1 = *iter1;
      for (auto iter2 = std::next(iter1); iter2 != rpv.taus.end(); iter2++) {
        TauCandidate tau2 = *iter2;
        Tlv tau1_tlv = ref_util::getTlv(tau1.rp_impl);
        Tlv tau2_tlv = ref_util::getTlv(tau2.rp_impl);
        double angle_between_taus = acos(ROOT::Math::VectorUtil::CosTheta(
        tau1_tlv, tau2_tlv));
        // The angular criterion for merging: The reconstructed tracks of both
        // taus are closer together than the search_cone_angle_.
        if (angle_between_taus > search_cone_angle) continue;
        ++n_merged;
        if (print_info) {
          tau2_tlv += tau1_tlv;
          streamlog_out(MESSAGE) << "Tau Merging: " << std::endl;
          streamlog_out(MESSAGE) << "  Angle was: " << angle_between_taus
            << ". Combined kinematics: " << "E: " << tau2_tlv.E()
            << ", theta: " << tau2_tlv.Theta() << ", phi: " << tau2_tlv.Phi()
            << "." << std::endl
            << "  Charged tracks: " << tau1.n_charged << "+"
            << tau2.n_charged << ", particles: "
            << tau1.rp_impl->getParticles().size() << "+"
            << tau2.rp_impl->getParticles().size() << std::endl;
        }
        *iter2 += tau1;  // tau2 += tau1 would update a copy instead of the
        // object inside the container.
        tc_to_remove.push_back(tau1);
        break;
      }
    }
  }
  for (auto tc : tc_to_remove) rpv.taus.remove(tc);
  return n_merged;
}
// ----------------------------------------------------------------------------

tau_cones_util::EventVector::EventVector(EVENT::LCCollection* pfo_collection,
    tau_cones_util::CandidateDefinition cd) {
  for (int e = 0; e < pfo_collection->getNumberOfElements(); ++e) {
    RP* particle = static_cast<RP*>(pfo_collection->getElementAt(e));
    Tlv tlv = ref_util::getTlv(particle);
	if ((tlv.Pt() > cd.p_t_cut) && (cos(tlv.Theta()) < cd.max_cos_theta)) {
        int charge = particle->getCharge();
        if (charge == 0) {
          neutral.push_back(particle);
        } else if (fabs(charge) == 1) {
          charged.push_back(particle);
        } else {
          streamlog_out(ERROR) << "The reconstructed particle has charge="
            << charge << ". How that? Only charges 0, +-1 should appear."
            << std::endl;
        }
    } else {
        ++n_background_suppressed;
    }
  }

  // Algorithms using the event vector might depend on the order of the charged
  // RPs (the possible seeds). To have a reliable result, order them by energy.
  std::sort(charged.begin(), charged.end(), ref_util::rpEnergySort);
}


std::vector<EVENT::ReconstructedParticle*>
    tau_cones_util::EventVector::EventVector::chargedAndNeutral() {
  std::vector<EVENT::ReconstructedParticle*> all = {};
  all.reserve(charged.size() + neutral.size());
  all.insert(all.begin(), charged.begin(), charged.end());
  all.insert(all.end(),   neutral.begin(), neutral.end());
  return all;
}