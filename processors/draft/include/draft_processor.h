/**
 *  Processor to test out some ideas without designing a full-fledged processor
 *  for Marlin.
 *
 *  To simplify fast changes, the inputs should be taken "hard-coded" in the
 *  .cc file (instead of getting parameters with the *register* mechanism from
 *  the Marlin steering file, define the parameter value in
 *  `DraftProcessor::init()`).
 *  Multiple (independent) drafts can be drawn up at the same time. The idea is
 *  to have all the necessary parameters in `init()` and `end()`.
 *  `processEvent()` is used to call the different drafts: The code that would
 *  usually be called in processEvent()` is outsourced into independent
 *  functions.
 *  When a draft becomes useless, its function and variables should be removed.
 *  In case a draft matures enough to be a processor in its own right, move it
 *  into a processor with the help of the `example` folder. Then remove it
 *  from here.
 *
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
 */
#ifndef _DRAFT_PROCESSOR_H_
#define _DRAFT_PROCESSOR_H_
// -- C++ STL headers.

// -- ROOT headers.
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

// -- LCIO headers.
#include "EVENT/ReconstructedParticle.h"

// -- Marlin headers.
#include "marlin/Processor.h"

// -- Header for this processor and other project-specific headers.

class DraftProcessor : public marlin::Processor {
 public:
  marlin::Processor* newProcessor() { return new DraftProcessor(); }
  DraftProcessor();

  // These two lines avoid frequent compiler warnings when using -Weffc++.
  DraftProcessor(const DraftProcessor&) = delete;
  DraftProcessor& operator=(const DraftProcessor&) = delete;

  void init();
  void processRunHeader(EVENT::LCRunHeader* run);
  void processEvent(EVENT::LCEvent* event);
  void end();

 private:
  // -- Parameters registered in steering file.
  // Collections

  // Parameters

  // Private draft functions
  void printCollectionInfo(EVENT::LCEvent* event);
  // Requires LCFIPlus to be run first. Investigates its flavor tagging output.
  // Might be useful for checking wether the LCFIPlus weights have to be
  // retrained for our usecase.
  void bbVertexPlayground(EVENT::LCEvent* event);

  // Dump a bunch of information about the reconstructed particle into the stream.
  void rpPrint(EVENT::ReconstructedParticle* rp); // bbVertexPlayground
  class MyJet {
   public :
    EVENT::ReconstructedParticle* jet;
    float b_tag;
    float c_tag;
    float bc_tag;
    bool is_assigned;
  };  // bbVertexPlayground

  // -- The root file
  std::string out_root_filename_{};
  TFile* root_out_{};
  // And its histograms
  TH2F* h_no_b_b_tags_2jets_{};  // bbVertexPlayground
  TH2F* h_b_b_tags_2jets_{};  // bbVertexPlayground
  TH2F* h_no_c__c_tags_2jets_{};  // bbVertexPlayground
  TH2F* h_cb_c_tags_2jets_{};  // bbVertexPlayground
  TH2F* h_c_no_b_c_tags_2jets_{};  // bbVertexPlayground
};
#endif