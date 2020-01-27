/**
 *  Simple processor that transforms a MCParticle collection into a
 *  ReconstructedParticle collection. This can can be helpful in identifying
 *  the source of performance loss from the MCTruth to a fully reconstructed
 *  event analysis (Particle ID, Detector acceptance,...).
 *
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
 */
#ifndef _MC_AS_RECONSTRUCTED_COLLECTION_PROCESSOR_H_
#define _MC_AS_RECONSTRUCTED_COLLECTION_PROCESSOR_H_
// -- C++ STL headers.

// -- ROOT headers.

// -- LCIO headers.

// -- Marlin headers.
#include "marlin/Processor.h"

// -- Header for this processor and other project-specific headers.

class McAsReconstructedProcessor : public marlin::Processor {
 public:
  marlin::Processor* newProcessor() { return new McAsReconstructedProcessor(); }
  McAsReconstructedProcessor();

  // These two lines avoid frequent compiler warnings when using -Weffc++.
  McAsReconstructedProcessor(const McAsReconstructedProcessor&) = delete;
  McAsReconstructedProcessor& operator=(const McAsReconstructedProcessor&)
      = delete;

  void init();
  void processRunHeader(EVENT::LCRunHeader* run);
  void processEvent(EVENT::LCEvent* event);
  void end();

 private:
  // -- Parameters registered in steering file.
  // Collections
  std::string mc_input_collection_name_{""};
  std::string now_rp_output_collection_name_{""};
  std::string mc_to_rp_relation_collection_name_{""};
  // Parameters
  bool keep_only_visible_{true};
  bool imitate_detector_acceptance_{true};
  float max_fb_cos_theta_{0.f};
  float min_p_t_for_charged_{0.f};
  float photon_ecal_cut_{0.f};
  float neutral_hcal_cut_{0.f};
  // Private helper functions
  float b_field_{0.f};
};
#endif