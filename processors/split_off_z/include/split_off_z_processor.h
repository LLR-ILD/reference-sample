/**
 *  Split a ReconstructedParticle collection into a Z-boson part,
 *  a Higgs-boson part and an overlay part. The output of this processor
 *  are thus three ReconstructedParticle collections.
 *  TODO: A more detailed description.
 *  TODO: Add an option that allows changing the PID source (from IsoLep
 *  TODO:   tagging to MCTruthRelation, PandoraPID, ...).
 *
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
 */
#ifndef _SPLIT_OFF_Z_PROCESSOR_H_
#define _SPLIT_OFF_Z_PROCESSOR_H_
// -- C++ STL headers.

// -- ROOT headers.
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
// The ROOT::Math vector functions need the GenVector component of root.
// Make sure you have the following line in the CMakeLists.txt file:
//   `FIND_PACKAGE(ROOT REQUIRED COMPONENTS GenVector)`.
#include "Math/VectorUtil.h"
#include "Math/Vector4D.h"

// -- LCIO headers.

// -- Marlin headers.
#include "marlin/Processor.h"

// -- Header for this processor and other project-specific headers.

class SplitOffZProcessor : public marlin::Processor {
 public:
  marlin::Processor* newProcessor() { return new SplitOffZProcessor(); }
  SplitOffZProcessor();

  // These two lines avoid frequent compiler warnings when using -Weffc++.
  SplitOffZProcessor(const SplitOffZProcessor&) = delete;
  SplitOffZProcessor& operator=(const SplitOffZProcessor&) = delete;

  void init();
  void processRunHeader(EVENT::LCRunHeader* run);
  void processEvent(EVENT::LCEvent* event);
  void end();

 private:
  // -- Parameters registered in steering file.
  // Collections
  std::string full_pfo_collection_name_{""};
  std::string higgs_only_collection_name_{""};
  std::string z_only_collection_name_{""};
  std::string overlay_collection_name_{""};

  // Parameters
  std::string z_decay_channel_{""};
  float photon_recombination_cos_{0.f};

  // Private helper functions
  float identifyZ(
      const std::vector<ROOT::Math::XYZTVector> &plus_candidate_momenta,
      const std::vector<ROOT::Math::XYZTVector> &minus_candidate_momenta,
      IntVec &plus_minus_indices, FloatVec &second_best_z_mass);

  // -- The root file
  std::string out_root_filename_{};
  TFile* root_out_{};
  // And its histograms
  TH1F* h_best_z_mass_{};
  TH2F* h_two_best_z_masses_{};
};
#endif