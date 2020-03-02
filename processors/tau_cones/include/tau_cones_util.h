/**
 *  Since the TauConesProcessor now uses a BDT, its functionality was split into
 *  multiple parts (preparing the training data, applying a trained model).
 *  Functionality used by multiple (sub)-processors is provided here (e.g. the
 *  algorithm for building a tau candidate.)
 *
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
 */
#ifndef _TAU_CONES_UTIL_H_
#define _TAU_CONES_UTIL_H_
// -- C++ STL headers.
#include <forward_list>

// -- ROOT headers.

// -- LCIO headers.
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/ReconstructedParticleImpl.h"

// -- Marlin headers.

// -- Header for this processor and other project-specific headers.

namespace tau_cones_util {

struct CandidateDefinition {
  float search_cone_angle = 0.f;
  float isolation_cone_angle = 0.f;
  float p_t_cut = 0.f;
  float max_cos_theta = 0.f;
  float min_p_t_seed = 0.f;
};

struct TauParameters {
  int n_charged = 0;
  int n_remnants = 0;
  float isolation_energy = 0.f;
  int n_in_isolation = 0.f;
  float m_invariant = 0.f;
};

struct CandidateCounts {
  int n_m_invariant_too_high = 0;
  int n_m_invariant_negative = 0;
  int n_wrong_track_number = 0;
  int n_not_isolated = 0;
  int n_tried_to_merge = 0;
  int n_taus_identified = 0;
  int n_background_suppressed = 0;

  CandidateCounts& operator+=(const CandidateCounts& rhs) {
    this->n_m_invariant_too_high += rhs.n_m_invariant_too_high;
    this->n_m_invariant_negative += rhs.n_m_invariant_negative;
    this->n_wrong_track_number += rhs.n_wrong_track_number;
    this->n_not_isolated += rhs.n_not_isolated;
    this->n_tried_to_merge += rhs.n_tried_to_merge;
    this->n_taus_identified += rhs.n_taus_identified;
    this->n_background_suppressed
        += rhs.n_background_suppressed;
    return *this;
  }
};

struct TauCandidate;
// From a vector of RPs, build a new RPImpl as well as additional variables
// (e.g. Lorentz vector) and collect them in a struct.
TauCandidate TauFromParts(std::vector<EVENT::ReconstructedParticle*> tau_parts,
    EVENT::ReconstructedParticle* seed);

struct TauCandidate {
  int n_charged = 0;
  IMPL::ReconstructedParticleImpl* rp_impl{};
  ROOT::Math::XYZTVector tlv{};
  std::vector<EVENT::ReconstructedParticle*> parts{};
  EVENT::ReconstructedParticle* seed{};

  // Needed for remove method of forward_list.
  bool operator==(const TauCandidate& rhs) {
    return this->rp_impl == rhs.rp_impl;
  }
  TauCandidate& operator+=(const TauCandidate& rhs) {
    this->parts.reserve(this->parts.size() + rhs.parts.size());
    this->parts.insert(this->parts.end(), rhs.parts.begin(), rhs.parts.end());
    if (rhs.seed->getEnergy() > this->seed->getEnergy()) this->seed = rhs.seed;
    TauCandidate dummy = TauFromParts(this->parts, this->seed);
    this->rp_impl = dummy.rp_impl;
    this->tlv = dummy.tlv;
    this->n_charged = dummy.n_charged;
    return *this;
  }
};

struct EventVector {
  int n_background_suppressed = 0;
  // Combines charged and neutral.
  std::vector<EVENT::ReconstructedParticle*> all{};
  std::vector<EVENT::ReconstructedParticle*> charged{};
  std::vector<EVENT::ReconstructedParticle*> neutral{};
  // Particles can be moved out of the above vectors into one of the taus.
  // Since we will want to remove candidates from any position inside the
  // container, forward_list is more appropriate than vector.
  std::forward_list<TauCandidate> taus{};

  EventVector() {all = {}; charged = {}; neutral = {}; taus = {};
    n_background_suppressed = 0;};
  EventVector(EVENT::LCCollection* pfo_collection, CandidateDefinition cd) ;
};

// In each iteration of this function, build a new tau as a vector of RPs.
// Any particle that is identified as part of the tau is removed from the
// particle vector of its type (charged_rps or neutral_rps).
// The new tau is appended to the taus vector.
// True is returned if a tau was found. If no tau can be found/seeded (any
// more), the function returns false. Note that in that case the three vectors
// which were passed by adresse remain unchanged.
bool FindATau(EventVector &rpv,
    double min_p_t_seed_, double search_cone_angle_,
    bool print_info=false);

// True if at least one tau candidate was found. Fills the taus vectors in rpv
// with RPs that are removed from the charged and neutral vectors.
bool FindAllTaus(EventVector &rpv,
    double min_p_t_seed, double search_cone_angle,
    bool print_info=false);

// Sometimes the FindTau algorithm splits a tau. Improve the performance by
// merging taus that are very close together.
int MergeCloseByTaus(EventVector &rpv,
    double search_cone_angle,
    bool print_info);
}  // namespace tau_cones_util

#endif