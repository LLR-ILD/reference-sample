/**
 *  Use this processor to investigate what you can do with the relations (in
 *  your .slcio files).
 *
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
 */
#ifndef _CHECK_RELATION_PROCESSOR_H_
#define _CHECK_RELATION_PROCESSOR_H_
// -- C++ STL headers.

// -- ROOT headers.

// -- LCIO headers.

// -- Marlin headers.
#include "marlin/Processor.h"

// -- Header for this processor and other project-specific headers.

class CheckRelationProcessor : public marlin::Processor {
 public:
  marlin::Processor* newProcessor() { return new CheckRelationProcessor(); }
  CheckRelationProcessor();

  // These two lines avoid frequent compiler warnings when using -Weffc++.
  CheckRelationProcessor(const CheckRelationProcessor&) = delete;
  CheckRelationProcessor& operator=(const CheckRelationProcessor&) = delete;

  void init();
  void processRunHeader(EVENT::LCRunHeader* run);
  void processEvent(EVENT::LCEvent* event);
  void end();

 private:
  // -- Parameters registered in steering file.
  // Collections
  std::string mc_collection_name_{""};
  std::string pfo_collection_name_{""};
  std::string relation_collection_name_{""};
  // Parameters
};
#endif