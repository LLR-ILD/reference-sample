/**
 *  Split a ReconstructedParticle collection into a Z-boson part,
 *  a Higgs-boson part and an overlay part. The output of this processor
 *  are thus three ReconstructedParticle collections.
 *  TODO: A more detailed description.
 *
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
 */
#ifndef _EXAMPLE_PROCESSOR_H_
#define _EXAMPLE_PROCESSOR_H_
// -- C++ STL headers.

// -- ROOT headers.
////#include "TFile.h"
////#include "TH1F.h"
////#include "TH2F.h"

// -- LCIO headers.

// -- Marlin headers.
#include "marlin/Processor.h"

// -- Header for this processor and other project-specific headers.

class ExampleProcessor : public marlin::Processor {
 public:
  marlin::Processor* newProcessor() { return new ExampleProcessor(); }
  ExampleProcessor();

  // These two lines avoid frequent compiler warnings when using -Weffc++.
  ExampleProcessor(const ExampleProcessor&) = delete;
  ExampleProcessor& operator=(const ExampleProcessor&) = delete;

  void init();
  void processRunHeader(EVENT::LCRunHeader* run);
  void processEvent(EVENT::LCEvent* event);
  void end();

 private:
  // -- Parameters registered in steering file.
  // Collections
  std::string pfo_collection_name_{""};
  // Parameters
  float pfo_energy_cut_{0.f};
  // Private helper functions

  // -- The root file
  ////std::string out_root_filename_{};
  ////TFile* root_out_{};
  // And its histograms
  ////TH1F* h_example_hist_{};
};
#endif