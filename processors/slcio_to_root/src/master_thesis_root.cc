/**
*    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
*/
// -- C++ STL headers.

// -- ROOT headers.
#include "TMath.h"
#include "Math/Vector4D.h"

// -- LCIO headers.
#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"

// -- Marlin headers.
#include "marlin/Exceptions.h"

// -- Header for this processor and other project-specific headers.
#include "ref_util.h"
#include "master_thesis_root.h"

// -- Using-declarations and global constants.
// Only in .cc files, never in .h header files!
using RP = EVENT::ReconstructedParticle;
using Tlv = ROOT::Math::XYZTVector;

// This line allows to register your processor in marlin when calling
// "Marlin steering_file.xml".
MasterThesisRootProcessor aMasterThesisRootProcessor;

// ----------------------------------------------------------------------------
MasterThesisRootProcessor::MasterThesisRootProcessor() :
    marlin::Processor("MasterThesisRootProcessor") {
  _description = "Build a .root file of the format necessary for my master "
    "thesis analysis.";

  registerInputCollection(
    LCIO::RECONSTRUCTEDPARTICLE,
    "HiggsCollection",
    "Name of the new PFO collection with (only) the Higgs decay remnants.",
    higgs_only_collection_name_,
    std::string("HRemnantsPFOs"));

  registerInputCollection(
    LCIO::RECONSTRUCTEDPARTICLE,
    "ZCollection",
    "Name of the new PFO collection with (only) the Z decay remnants.",
    z_only_collection_name_,
    std::string("ZRemnantsPFOs"));

  registerProcessorParameter(
    "ZDecayChannel",
    "Can be ZDec_EL, ZDec_MU, ZDec_TAU or ZDec_NU.",
    z_decay_channel_,
    std::string("ZDec_MU"));

  registerProcessorParameter(
    "OutputRootFile",
    "Name of the output root file.",
    root_file_name_,
    std::string("master_thesis"));
}

// ----------------------------------------------------------------------------

void MasterThesisRootProcessor::initRoot() {
  const char* tree_name = z_decay_channel_.c_str();
  if (z_decay_channel_ == "ZDec_EL") {
    tree_name = "e1Tree";
  } else if (z_decay_channel_ == "ZDec_MU") {
    tree_name = "e2Tree";
  } else if (z_decay_channel_ == "ZDec_TAU") {
    tree_name = "e3Tree";
  } else if (z_decay_channel_ == "ZDec_NU") {
    tree_name = "nuTree";
  }
  tree_ = new TTree(tree_name, "Input for my master thesis analysis.");
  tv.initBranches(tree_);
}

void MasterThesisRootProcessor::endRoot() {
  TString fnn(root_file_name_.c_str()); fnn+=".root";
  root_file_ = new TFile(fnn,"update");
  root_file_->cd();
  TTree* tree_in_write_file = tree_->CloneTree();
  tree_in_write_file->Write();
  root_file_->Write();
  root_file_->Close();
}

void MasterThesisRootProcessor::init() {
  printParameters();
  initRoot();
}

void MasterThesisRootProcessor::end() {
  endRoot();
}

// ----------------------------------------------------------------------------
void MasterThesisRootProcessor::processEvent(EVENT::LCEvent* event) {
  streamlog_out(DEBUG) << "Processing event no " << event->getEventNumber()
    << " - run " << event->getEventNumber() << std::endl;

  tv.resetValues();
  EVENT::LCCollection* z_collection = nullptr;
  try {
    z_collection = event->getCollection(z_only_collection_name_);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "RP collection " << z_only_collection_name_
      << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }
  Tlv z_tlv = {0, 0, 0, 0};
  for (int e = 0; e < z_collection->getNumberOfElements(); ++e) {
    RP* z_part = static_cast<RP*>(z_collection->getElementAt(e));
    z_tlv += ref_util::getTlv(z_part);
  }

  EVENT::LCCollection* higgs_collection = nullptr;
  try {
    higgs_collection = event->getCollection(higgs_only_collection_name_);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "RP collection " << higgs_only_collection_name_
      << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }
  Tlv higgs_tlv = {0, 0, 0, 0};
  for (int e = 0; e < higgs_collection->getNumberOfElements(); ++e) {
    RP* higgs_part = static_cast<RP*>(higgs_collection->getElementAt(e));
    higgs_tlv += ref_util::getTlv(higgs_part);

    int abs_pdg = abs(higgs_part->getType());
    if (abs_pdg == 11) {
      tv.n_electrons += 1;
    } else if (abs_pdg == 22) {
      tv.n_gamma += 1;
    } else if (abs_pdg == 13) {
      tv.n_muons += 1;
    } else if (abs_pdg == 211 || abs_pdg == 321 || abs_pdg == 2212) {
      tv.n_ch_hadrons += 1;
    } else if (abs_pdg == 130 || abs_pdg == 310 || abs_pdg == 2112) {
      tv.n_n_hadrons += 1;
    } else {
      streamlog_out(WARNING) << "An unexpected PDG was found: " << abs_pdg
        << "." << std::endl;
    }
  }
  std::cout << z_tlv.mass() << std::endl;
  tv.m_z = z_tlv.mass();
  tv.m_recoil = (Tlv(0, 0, 0, 250) - z_tlv).mass();
  tv.m_vis = (z_tlv + higgs_tlv).mass();
  tv.m_miss = (Tlv(0, 0, 0, 250) - z_tlv - higgs_tlv).mass();
  tv.m_h = higgs_tlv.mass();
  tv.m_h_recoil = (Tlv(0, 0, 0, 250) - higgs_tlv).mass();

  tv.cos_theta_z = cos(z_tlv.theta());
  tv.cos_theta_miss = cos((z_tlv +  higgs_tlv).theta());

  tree_->Fill();
}