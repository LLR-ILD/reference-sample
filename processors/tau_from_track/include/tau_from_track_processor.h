/**
 *  Short description of the purpose of this processor.
 *
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
 */
#ifndef _TAU_FROM_TRACK_PROCESSOR_H_
#define _TAU_FROM_TRACK_PROCESSOR_H_
// -- C++ STL headers.

// -- ROOT headers.
#include "TFile.h"
#include "TTree.h"

// -- LCIO headers.

// -- Marlin headers.
#include "marlin/Processor.h"

// -- Header for this processor and other project-specific headers.

class TauFromTrackProcessor : public marlin::Processor {
 public:
  marlin::Processor* newProcessor() { return new TauFromTrackProcessor(); }
  TauFromTrackProcessor();

  // These two lines avoid frequent compiler warnings when using -Weffc++.
  TauFromTrackProcessor(const TauFromTrackProcessor&) = delete;
  TauFromTrackProcessor& operator=(const TauFromTrackProcessor&) = delete;

  void init();
  void processRunHeader(EVENT::LCRunHeader* run);
  void processEvent(EVENT::LCEvent* event);
  void end();

 private:
  // -- Parameters registered in steering file.
  // Collections
  std::string pfo_collection_name_{""};
  // -- The root file
  TFile* root_file_{};
  std::string root_file_name_ = {""};
  TTree*      track_tree_ = {};
  std::string track_tree_name_ = {""};
  double TanLambda_{0.};
  double Z0_{0.};
  double D0_{0.};
  double RadiusOfInnermostHit_{0.};
  double Phi_{0.};
  double Omega_{0.};
  double e_track_{0.};
  double Ndf_{0.};
  double dEdxError_{0.};
  double dEdx_{0.};
  double Chi2_{0.};

  double cov_d0_d0_{0.};
  double cov_phi_d0_{0.};
  double cov_phi_phi_{0.};
  double cov_omega_d0_{0.};
  double cov_omega_phi_{0.};
  double cov_omega_omega_{0.};
  double cov_z0_d0_{0.};
  double cov_z0_phi_{0.};
  double cov_z0_omega_{0.};
  double cov_z0_z0_{0.};
  double cov_tanlambda_d0_{0.};
  double cov_tanlambda_phi_{0.};
  double cov_tanlambda_omega_{0.};
  double cov_tanlambda_z0_{0.};
  double cov_tanlambda_tanlambda_{0.};

  int Type_{-1};

  int n_no_higgs_invisible_{-1};
};
#endif