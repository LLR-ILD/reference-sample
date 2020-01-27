/**
*    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
*/
// -- C++ STL headers.

// -- ROOT headers.

// -- LCIO headers.
#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
#include "UTIL/LCTOOLS.h"

// -- Marlin headers.
#include "marlin/Exceptions.h"

// -- Header for this processor and other project-specific headers.
#include "skip_empty_events_processor.h"

// -- Using-declarations and global constants.
// Only in .cc files, never in .h header files!

// This line allows to register your processor in marlin when calling
// "Marlin steering_file.xml".
SkipEmptyEventsProcessor aSkipEmptyEventsProcessor;

// ----------------------------------------------------------------------------
SkipEmptyEventsProcessor::SkipEmptyEventsProcessor() :
    marlin::Processor("SkipEmptyEventsProcessor") {
  // _description comes from marlin::Processor.
  _description = "Processor that skips events for which the specified "
    "collection is empty. ";

  registerProcessorParameter(
    "SkipCollection",
    "Name of the LCCollection for which the event should be skipped in case it "
      "is empty.",
    collection_name_,
    std::string("PandoraPFOs"));
}

// ----------------------------------------------------------------------------
void SkipEmptyEventsProcessor::init() {
  // Usually a good idea to print parameters.
  printParameters(); // method from marlin::Processor.
}

// ----------------------------------------------------------------------------
void SkipEmptyEventsProcessor::processRunHeader(EVENT::LCRunHeader* run) {
  streamlog_out(MESSAGE) << "Starting run no " << run->getRunNumber()
    << std::endl;
}

// ----------------------------------------------------------------------------
void SkipEmptyEventsProcessor::processEvent(EVENT::LCEvent* event) {
  streamlog_out(DEBUG) << "Processing event no " << event->getEventNumber()
    << " - run " << event->getEventNumber() << std::endl;
  EVENT::LCCollection* collection = nullptr;
  try {
    collection = event->getCollection(collection_name_);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "The collection " << collection_name_
      << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }
  streamlog_out(DEBUG) << "Number of elements in " <<  collection_name_ << ": "
    << collection->getNumberOfElements() <<"." << std::endl;
  if (collection->getNumberOfElements() == 0) {
    streamlog_out(MESSAGE) << "The collection " << collection_name_
      << "in event number " << event->getEventNumber()
      << "has no Pandora objects. Skip the further investigation of this event."
      << std::endl;
    throw marlin::SkipEventException(this);
  }
}

// ----------------------------------------------------------------------------
void SkipEmptyEventsProcessor::end() {
  streamlog_out(MESSAGE) << "end" << std::endl;
}