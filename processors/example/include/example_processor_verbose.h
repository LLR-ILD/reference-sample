/**
 *  Short description of the purpose of this processor.
 *
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
 */
#ifndef _EXAMPLE_PROCESSOR_VERBOSE_H_
#define _EXAMPLE_PROCESSOR_VERBOSE_H_
// -- C++ STL headers.

// -- ROOT headers.
// -- ROOT headers.
////#include "TFile.h"
////#include "TH1F.h"
////#include "TH2F.h"

// -- LCIO headers.

// -- Marlin headers.
#include "marlin/Processor.h"

// -- Header for this processor and other project-specific headers.

class ExampleProcessor : public marlin::Processor {
 public:
  /**
   *  @brief  Factory method to create a processor.
   *
   *  This method is called internally by Marlin to create an instance of your
   *  processor.
   */
  marlin::Processor* newProcessor() { return new ExampleProcessor(); }

  /** @brief  Constructor.
   *  Register your parameters in there. See example_processor.cc.
   */
  ExampleProcessor();

  // These two lines avoid frequent compiler warnings when using -Weffc++.
  ExampleProcessor(const ExampleProcessor&) = delete;
  ExampleProcessor& operator=(const ExampleProcessor&) = delete;

  /**
   *  @brief  Called once at the begin of the job before anything is read.
   *
   *  Use to initialize the processor, e.g. book histograms. The parameters
   *  that you have registered in the constructor are initialized.
   */
  void init();

  /**
   *  @brief  Called for every run read from the file.
   *
   *  If you use simulation data output from ddsim (e.g from production data on
   *  the grid), this will be called before the first event and only once in
   *  the job.
   */
  void processRunHeader(EVENT::LCRunHeader* run);

  /**
   *  @brief  Called for every event - the working horse.
   *
   *  Analise your reconstructed particle objects from the event object here.
   *  See ExampleProcessor.cc for an example
   */
  void processEvent(EVENT::LCEvent* event);

  /**
   *  @brief  Called after data processing for clean up.
   *
   *  You have allocated memory somewhere in your code?
   *  This is the best place to clean your mess!
   */
  void end();

 private:
  /**
   *  @brief  Member initializations.
   *
   *  Initialize your members in the class definition to
   *  be more efficient and avoid compiler warnings.
   */
  // -- Parameters registered in steering file.
  // Collections
  std::string pfo_collection_name_{""};
  // Parameters
  float pfo_energy_cut_{0.f};
  // Private helper functions

  // -- The root file
  ////std::string out_root_filename_{};
  ////TFile* root_out_{};
  // And its histograms
  ////TH1F* h_example_hist_{};
};
#endif