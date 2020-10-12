/**
 *  Run this processor on Pqqh events instead of running SplitOffZProcessor
 *  to imitate Pnnh events.
 *  Since the Pnnh events are needed to train a BDT in increasing the available
 *  statistics in this way is a great help.
 *
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
 */
#ifndef _FAKE_Z_NUNU_FROM_QQ_PROCESSOR_H_
#define _FAKE_Z_NUNU_FROM_QQ_PROCESSOR_H_
// -- C++ STL headers.

// -- ROOT headers.

// -- LCIO headers.
#include "EVENT/ReconstructedParticle.h"

// -- Marlin headers.
#include "marlin/Processor.h"

// -- Header for this processor and other project-specific headers.

class FakeZnunuFromZqqProcessor : public marlin::Processor {
 public:
  marlin::Processor* newProcessor() { return new FakeZnunuFromZqqProcessor(); }
  FakeZnunuFromZqqProcessor();

  // These two lines avoid frequent compiler warnings when using -Weffc++.
  FakeZnunuFromZqqProcessor(const FakeZnunuFromZqqProcessor&) = delete;
  FakeZnunuFromZqqProcessor& operator=(const FakeZnunuFromZqqProcessor&) = delete;

  // void init();
  // void processRunHeader(EVENT::LCRunHeader* run);
  // void end();
  void processEvent(EVENT::LCEvent* event);

 private:
  // -- Parameters registered in steering file.
  // Collections
  std::string full_pfo_collection_name_{""};
  std::string mc_collection_name_{""};
  std::string relation_collection_name_{""};
  std::string higgs_only_collection_name_{""};
  std::string z_only_collection_name_{""};
  std::string overlay_collection_name_{""};

  // Private helper functions
  void fillQQRemnants(std::vector<ReconstructedParticle*> &zqq_remnants,
                      EVENT::LCEvent* event);
};
#endif