/**
 *  Build some histograms in order to get to know better the external
 *  IsolatedLeptonTagging processor.
 *
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
 */
#ifndef _CHECK_ISOLATED_LEPTON_TAGGING_PROCESSOR_H_
#define _CHECK_ISOLATED_LEPTON_TAGGING_PROCESSOR_H_
// -- C++ STL headers.

// -- ROOT headers.
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

// -- LCIO headers.

// -- Marlin headers.
#include "marlin/Processor.h"

// -- Header for this processor and other project-specific headers.

class CheckIsoLeptonProcessor : public marlin::Processor {
 public:
  marlin::Processor* newProcessor() { return new CheckIsoLeptonProcessor(); }
  CheckIsoLeptonProcessor();

  // These two lines avoid frequent compiler warnings when using -Weffc++.
  CheckIsoLeptonProcessor(const CheckIsoLeptonProcessor&) = delete;
  CheckIsoLeptonProcessor& operator=(const CheckIsoLeptonProcessor&) = delete;

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
  std::string out_root_filename_{};
  TFile* root_out_{};
  // And its histograms
  TH1F* h_n_isolated_muons_{};
  TH1F* h_n_isolated_electrons_{};
  TH1F* h_mva_isolated_muon_{};
  TH1F* h_mva_isolated_muon_possibly_z_{};
  TH1F* h_mva_isolated_electron_{};
  TH1F* h_mva_isolated_electron_possibly_z_{};
  TH2I* h_n_objects_per_type_{};
};
#endif