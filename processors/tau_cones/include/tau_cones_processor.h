/**
 *  Identify the reconstructed particles steming from a tau with the tau cone/
 *  reconstruction cone approach.
 *
 *  This processor was written starting from the TauFinder processor as found at
 *    /cvmfs/ilc.desy.de/sw/x86_64_gcc49_sl6/v02-00-02/MarlinReco/v01-25/Analysis/TauFinder.
 *
 *    @author A. Muennich, CERN (original TauFinder).
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP (adaptation).
 */
#ifndef _TAU_CONES_PROCESSOR_H_
#define _TAU_CONES_PROCESSOR_H_
// -- C++ STL headers.

// -- ROOT headers.
#include "TFile.h"
#include "TNtuple.h"

// -- LCIO headers.
#include "EVENT/ReconstructedParticle.h"

// -- Marlin headers.
#include "marlin/Processor.h"

// -- Header for this processor and other project-specific headers.

class TauConesProcessor : public marlin::Processor {
 public:
  marlin::Processor* newProcessor() { return new TauConesProcessor(); }
  TauConesProcessor();

  // These two lines avoid frequent compiler warnings when using -Weffc++.
  TauConesProcessor(const TauConesProcessor&) = delete;
  TauConesProcessor& operator=(const TauConesProcessor&) = delete;

  void init();
  void processRunHeader(EVENT::LCRunHeader* run);
  void processEvent(EVENT::LCEvent* event);
  void end();

 private:
  // -- Parameters registered in steering file.
  // Collections
  std::string pfo_collection_name_{""};
  std::string tau_collection_name_{""};
  std::string rest_collection_name_{""};
  std::string tau_relation_collection_name_{""};
  // Parameters
  float search_cone_angle_{0.f};
  float isolation_cone_angle_{0.f};
  float max_isolation_energy_{0.f};
  float p_t_cut_{0.f};
  float max_cos_theta_{0.f};
  float min_p_t_seed_{0.f};
  float max_m_invariant_{0.f};

  int n_m_invariant_too_high_{-1};
  int n_m_invariant_negative_{-1};
  int n_wrong_track_number_{-1};
  int n_not_isolated_{-1};
  int n_tried_to_merge_{-1};
  int n_taus_identified_{-1};
  int n_background_suppressed_particles_{-1};

  int n_m_invariant_too_high_per_event_{-1};
  int n_m_invariant_negative_per_event_{-1};
  int n_wrong_track_number_per_event_{-1};
  int n_not_isolated_per_event_{-1};
  int n_tried_to_merge_per_event_{-1};
  int n_taus_identified_per_event_{-1};
  int n_background_suppressed_particles_per_event_{-1};

  int n_events_total_{-1};
  int n_runs_{-1};

  // In each iteration of this function, build a new tau as a vector of RPs.
  // Any particle that is identified as part of the tau is removed from the
  // particle vector of its type (charged_rps or neutral_rps).
  // The new tau is appended to the taus vector.
  // True is returned if a tau was found. If no tau can be found/seeded (any
  // more), the function returns false. Note that in that case the three vectors
  // which were passed by adresse remain unchanged.
  bool FindTau(std::vector<EVENT::ReconstructedParticle*> &charged_rps,
               std::vector<EVENT::ReconstructedParticle*> &neutral_rps,
	           std::vector<std::vector<EVENT::ReconstructedParticle*>> &taus);

  // -- The root file
  std::string out_root_filename_{};
  TFile* root_out_{};
  TNtuple* fail_reason_tuple_{};
};
#endif