/**
 *  Almost trivial processor that skips the event if the specified collection
 *  has no members.
 *
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
 */
#ifndef _SKIP_EMPTY_EVENTS_PROCESSOR_H_
#define _SKIP_EMPTY_EVENTS_PROCESSOR_H_
// -- C++ STL headers.

// -- ROOT headers.

// -- LCIO headers.

// -- Marlin headers.
#include "marlin/Processor.h"

// -- Header for this processor and other project-specific headers.

class SkipEmptyEventsProcessor : public marlin::Processor {
 public:
  marlin::Processor* newProcessor() { return new SkipEmptyEventsProcessor(); }
  SkipEmptyEventsProcessor();

  // These two lines avoid frequent compiler warnings when using -Weffc++.
  SkipEmptyEventsProcessor(const SkipEmptyEventsProcessor&) = delete;
  SkipEmptyEventsProcessor& operator=(const SkipEmptyEventsProcessor&) = delete;

  void init();
  void processRunHeader(EVENT::LCRunHeader* run);
  void processEvent(EVENT::LCEvent* event);
  void end();

 private:
  // -- Parameters registered in steering file.
  // Collections
  // Parameters
  std::string collection_name_{""};
  // Private helper functions
};
#endif