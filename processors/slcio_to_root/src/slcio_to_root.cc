/**
*    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
*/
// -- C++ STL headers.

// -- ROOT headers.
#include "TMath.h"
#include "Math/Vector4D.h"

// -- LCIO headers.
#include "EVENT/MCParticle.h"
#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
#include "UTIL/LCRelationNavigator.h"

// -- Marlin headers.
#include "marlin/Exceptions.h"

// -- Header for this processor and other project-specific headers.
#include "ref_util.h"
#include "slcio_to_root.h"

// -- Using-declarations and global constants.
// Only in .cc files, never in .h header files!
using RP = EVENT::ReconstructedParticle;

// This line allows to register your processor in marlin when calling
// "Marlin steering_file.xml".
SlcioToRootProcessor aSlcioToRootProcessor;

// ----------------------------------------------------------------------------
SlcioToRootProcessor::SlcioToRootProcessor() :
    marlin::Processor("SlcioToRootProcessor") {
  // _description comes from marlin::Processor.
  _description = "Fill some ReconstructedParticle information into a .root "
    "file. ";

  registerInputCollection(
    LCIO::RECONSTRUCTEDPARTICLE,
    "PfoCollection",
    "The RP collection whoms (partial) information should be passed into the "
      ".root file.",
    rp_collection_name_,
    std::string("PandoraPFOs"));

  registerProcessorParameter(
    "outputRootFileName",
    "Name of the output root file.",
    root_file_name_,
    std::string("rp_simple"));

  registerProcessorParameter(
    "treeName",
    "Name of the tree. Should be the process type.",
    tree_name_,
    std::string("process?"));
}

// ----------------------------------------------------------------------------
void SlcioToRootProcessor::init() {
  // Usually a good idea to print parameters.
  printParameters(); // method from marlin::Processor.
  // This is also the place where a root file would be opened
  // and histograms would be defined.
  TString fnn(root_file_name_.c_str()); fnn+=".root";
  root_file_ = new TFile(fnn,"update");

  tree_ = new TTree(tree_name_.c_str(), "Tree with some data from "
    "reconstructed particle collections in .slcio files.");
  initBranches(tree_);

}

// ----------------------------------------------------------------------------
void SlcioToRootProcessor::processRunHeader(EVENT::LCRunHeader* run) {
  streamlog_out(MESSAGE) << "Starting run no " << run->getRunNumber()
    << std::endl;
}

// ----------------------------------------------------------------------------
void SlcioToRootProcessor::processEvent(EVENT::LCEvent* event) {
  streamlog_out(DEBUG) << "Processing event no " << event->getEventNumber()
    << " - run " << event->getEventNumber() << std::endl;
  n_event_ = event->getEventNumber();
  EVENT::LCCollection* rp_collection = nullptr;
  try {
    rp_collection = event->getCollection(rp_collection_name_);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "RP collection " << rp_collection_name_
      << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }
  UTIL::LCRelationNavigator* relation_navigator = nullptr;
  EVENT::LCCollection* relation_collection = nullptr;
  try {
    relation_collection = event->getCollection("RecoMCTruthLink");
    relation_navigator = new UTIL::LCRelationNavigator(relation_collection);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "The relation collection " << "RecoMCTruthLink"
      << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }
  for (int e = 0; e < rp_collection->getNumberOfElements(); ++e) {
    RP* particle = static_cast<RP*>(rp_collection->getElementAt(e));
    // Now collect the actual information.
    EVENT::MCParticle* earliest_mc_parent = ref_util::getMcChainFromRp(
        particle, relation_navigator).back();
    mc_parent_ = earliest_mc_parent->getPDG();
    ROOT::Math::XYZTVector rp_tlv =  ref_util::getTlv(particle);
    px_    = rp_tlv.X();
    py_    = rp_tlv.Y();
    pz_    = rp_tlv.Z();
    e_     = rp_tlv.E();
    phi_   = rp_tlv.Phi();
    theta_ = rp_tlv.Theta();
    //
    n_particle_ = e;
    charge_ = particle->getCharge();
    pid_ = particle->getType();
    tree_->Fill();
  }
}

// ----------------------------------------------------------------------------
void SlcioToRootProcessor::end() {
  // Fill the root file.
  root_file_->cd();
  tree_->Write();
  root_file_->Write(0);
  root_file_->Close();
  streamlog_out(MESSAGE) << "end()." << std::endl;
}

//-----------------------------------------------------------------------------
void SlcioToRootProcessor::initBranches(TTree* tree) {
  if (tree_ == 0) {
    streamlog_out(ERROR) << "::initBranches - Invalid tree pointer!"
      << std::endl;
    throw marlin::StopProcessingException(this);
  }
  tree->Branch(("px"),     &px_,     ("px/D"));
  tree->Branch(("py"),     &py_,     ("py/D"));
  tree->Branch(("pz"),     &pz_,     ("pz/D"));
  tree->Branch(("E"),      &e_,      ("E/D") );
  tree->Branch(("phi"),    &phi_,    ("phi/D") );
  tree->Branch(("theta"),  &theta_,  ("theta/D") );
  //
  tree->Branch(("charge"),      &charge_,      ("charge/I"));
  tree->Branch(("nEvent"),      &n_event_,      ("nEvent/I"));
  tree->Branch(("nParticle"),   &n_particle_,   ("nParticle/I"));
  tree->Branch(("pid"),         &pid_,         ("pid/I"));
  tree->Branch(("mcParent"),    &mc_parent_,    ("mcParent/I"));
}