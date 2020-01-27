/**
 *  Extract some (minimal) information from .slcio file(s):
 *    4-Momentum, PID and charge, corresponding event number and process type.
 *  Necessary for analysis outside of the Marlin framework.
 *
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
 */
#ifndef _SLCIO_TO_ROOT_PROCESSOR_H_
#define _SLCIO_TO_ROOT_PROCESSOR_H_
// -- C++ STL headers.

// -- ROOT headers.
#include "TFile.h"
#include "TTree.h"

// -- LCIO headers.

// -- Marlin headers.
#include "marlin/Processor.h"

// -- Header for this processor and other project-specific headers.

class SlcioToRootProcessor : public marlin::Processor {
 public:
  marlin::Processor* newProcessor() { return new SlcioToRootProcessor(); }
  SlcioToRootProcessor();

  // These two lines avoid frequent compiler warnings when using -Weffc++.
  SlcioToRootProcessor(const SlcioToRootProcessor&) = delete;
  SlcioToRootProcessor& operator=(const SlcioToRootProcessor&) = delete;

  void init();
  void processRunHeader(EVENT::LCRunHeader* run);
  void processEvent(EVENT::LCEvent* event);
  void end();

 private:
  // -- Parameters registered in steering file.
  // Collections
  std::string rp_collection_name_{""};
  // Parameters
  double px_{0.f};
  double py_{0.f};
  double pz_{0.f};
  double e_{0.f};
  double phi_{0.f};
  double theta_{0.f};

  int charge_{0};
  int n_event_{0};
  int n_particle_{0};
  int pid_{0};
  int mc_parent_{0};

  // -- The root file
  TFile* root_file_{};
  std::string root_file_name_ = {""};
  TTree*      tree_ = {};
  std::string tree_name_ = {""};

  // Private helper functions
  void initBranches(TTree* tree);
};
#endif