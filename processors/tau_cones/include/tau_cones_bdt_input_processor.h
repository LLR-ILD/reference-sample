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
#ifndef _TAU_CONES_BDT_INPUT_PROCESSOR_H_
#define _TAU_CONES_BDT_INPUT_PROCESSOR_H_
// -- C++ STL headers.

// -- ROOT headers.
#include "TFile.h"
#include "TNtuple.h"

// -- LCIO headers.
#include "EVENT/ReconstructedParticle.h"

// -- Marlin headers.
#include "marlin/Processor.h"

// -- Header for this processor and other project-specific headers.
#include "tau_cones_util.h"

class TauConesBDTInputProcessor : public marlin::Processor {
 public:
  marlin::Processor* newProcessor() { return new TauConesBDTInputProcessor(); }
  TauConesBDTInputProcessor();

  // These two lines avoid frequent compiler warnings when using -Weffc++.
  TauConesBDTInputProcessor(const TauConesBDTInputProcessor&) = delete;
  TauConesBDTInputProcessor& operator=(const TauConesBDTInputProcessor&) =
      delete;

  void init();
  void processEvent(EVENT::LCEvent* event);
  void end();

 private:
  // -- Parameters registered in steering file.
  // Collections
  std::string pfo_collection_name{""};

  typedef tau_cones_util::CandidateDefinition CandidateDefinition;
  CandidateDefinition cd{};

  typedef tau_cones_util::TauParameters TauParameters;

  typedef tau_cones_util::EventVector EventVector;
  typedef tau_cones_util::TauCandidate TauCandidate;

  // Proper definition in tau_cones_utils.h.
  bool FindAllTaus(EventVector &rpv,
                double min_p_t_seed, double search_cone_angle,
                bool print_info=false) {
    return tau_cones_util::FindAllTaus(rpv,
            min_p_t_seed, search_cone_angle, print_info);
  };

  // Proper definition in tau_cones_utils.h.
  int MergeCloseByTaus(EventVector &rpv, double search_cone_angle,
                    bool print_info=false) {
      return tau_cones_util::MergeCloseByTaus(rpv,
                search_cone_angle, print_info);
    };

  bool APartFromTau(TauCandidate tau_candidate,
    UTIL::LCRelationNavigator* relation_navigator);

  // -- The root file
  std::string bdt_inputs_name{};
  TFile* bdt_inputs{};
  TNtuple* tau_parameters{};
};
#endif