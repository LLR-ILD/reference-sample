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
#include "IMPL/LCCollectionVec.h"

// -- Marlin headers.
#include "marlin/Processor.h"

// -- Header for this processor and other project-specific headers.
#include "tau_cones_util.h"

class TauConesProcessor : public marlin::Processor {
 public:
  marlin::Processor* newProcessor() { return new TauConesProcessor(); }
  TauConesProcessor();

  // These two lines avoid frequent compiler warnings when using -Weffc++.
  TauConesProcessor(const TauConesProcessor&) = delete;
  TauConesProcessor& operator=(const TauConesProcessor&) = delete;

  void init();
  void processEvent(EVENT::LCEvent* event);
  void end();

 private:
  // -- Parameters registered in steering file.
  // Collections
  std::string pfo_collection_name_{""};
  std::string tau_collection_name_{""};
  std::string rest_collection_name_{""};
  std::string tau_relation_collection_name_{""};

  typedef tau_cones_util::CandidateDefinition CandidateDefinition;
  CandidateDefinition cd{};

  typedef tau_cones_util::TauParameters TauParameters;
  TauParameters tau_cut{};

  typedef tau_cones_util::EventVector EventVector;
  typedef tau_cones_util::TauCandidate TauCandidate;
  typedef tau_cones_util::CandidateCounts CandidateCounts;
  CandidateCounts total_count{};

  int n_events_total{-1};

  // Proper definition in tau_cones_utils.h.
  bool HasStrictTau(EventVector &rpv, CandidateDefinition cad,
      EVENT::LCEvent* event) {
    return tau_cones_util::HasStrictTau(rpv, cad, event);
  };

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

  // Proper definition in tau_cones_processor.cpp.
  int MassTracksRejection(EventVector &rpv, CandidateCounts &event_count,
    IMPL::LCCollectionVec* rest_collection);

  // Proper definition in tau_cones_processor.cpp.
  int Isolation(EventVector &rpv,
    CandidateCounts &event_count, IMPL::LCCollectionVec* rest_collection);

  int getHiggsTruth(EVENT::LCEvent* event);
  // -- The root file
  std::string out_root_filename_{};
  TFile* root_out_{};
  TNtuple* fail_reason_tuple_{};
};
#endif