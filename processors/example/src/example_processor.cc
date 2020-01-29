/**
*    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
*/
// -- C++ STL headers.

// -- ROOT headers.

// -- LCIO headers.
#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
////#include "UTIL/LCTOOLS.h"

// -- Marlin headers.
#include "marlin/Exceptions.h"

// -- Header for this processor and other project-specific headers.
#include "example_processor.h"

// -- Using-declarations and global constants.
// Only in .cc files, never in .h header files!
using RP = EVENT::ReconstructedParticle;
using MCP = EVENT::MCParticle;

// This line allows to register your processor in marlin when calling
// "Marlin steering_file.xml".
ExampleProcessor aExampleProcessor;

// ----------------------------------------------------------------------------
ExampleProcessor::ExampleProcessor() :
    marlin::Processor("ExampleProcessor") {
  // _description comes from marlin::Processor.
  _description = "An example processor for ILD analysis. ";

  registerInputCollection(
    LCIO::RECONSTRUCTEDPARTICLE,  // The collection type.
    "PfoCollection",  // Parameter name in the steering file.
    "The Pandora PFO collection name.",  // A parameter description. Please fill this correctly.
    pfo_collection_name_,  // Your variable to store the result after steering file parsing.
    std::string("PandoraPFOs"));  // Default parameter value.

  registerProcessorParameter(
    "PfoEnergyCut",
    "A cut to apply on the pfo energy.",
    pfo_energy_cut_,
    0.f);
}

// ----------------------------------------------------------------------------
void ExampleProcessor::init() {
  // Usually a good idea to print parameters.
  printParameters(); // method from marlin::Processor.
  // This is also the place where a root file would be opened
  // and histograms would be defined.
  ////TString fnn(out_root_filename_.c_str());
  ////fnn += ".root";
  ////root_out_ = new TFile(fnn, "recreate");
  ////h_example_hist_ = new TH1F("exampleHist", "exampleHist", 100, 0, 100);
}

// ----------------------------------------------------------------------------
void ExampleProcessor::processRunHeader(EVENT::LCRunHeader* run) {
  streamlog_out(MESSAGE) << "Starting run no " << run->getRunNumber()
    << std::endl;
  // LCRunHeader objects can be printed using the LCTOOLS class.
  ////UTIL::LCTOOLS::dumpRunHeader(run);
}

// ----------------------------------------------------------------------------
void ExampleProcessor::processEvent(EVENT::LCEvent* event) {
  streamlog_out(DEBUG) << "Processing event no " << event->getEventNumber()
    << " - run " << event->getEventNumber() << std::endl;
  // Prefer accessing any collections with a try/catch approach.
  EVENT::LCCollection* collection = nullptr;
  try {
    collection = event->getCollection(pfo_collection_name_);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "RP collection" << pfo_collection_name_
      << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }
  for (int e = 0; e < collection->getNumberOfElements(); ++e) {
    // A cautious reading of a particle involves casting it with static_cast
    // (instead of dynamic_cast) and ensuring it is not a nullpointer.
    RP* particle = static_cast<RP*>(collection->getElementAt(e));
    // Start analysing your particle here! Enjoy!
    if (particle->getEnergy() < pfo_energy_cut_) {
      continue;
    }
  }
}

// ----------------------------------------------------------------------------
void ExampleProcessor::end() {
  // Write the histograms.
  ////root_out_->cd();
  ////root_out_->Write(0);
  ////root_out_->Close();
  // Cleanup your mess here!
  streamlog_out(MESSAGE) << "end" << std::endl;
}