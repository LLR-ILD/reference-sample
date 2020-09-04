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
using MCP = EVENT::MCParticle;
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


  registerInputCollection(
    LCIO::MCPARTICLE,
    "MCParticleCollection",
    "MCParticle collection for Higgs decay truth information.",
    mc_collection_name,
    std::string("MCParticlesSkimmed"));

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
  if (missing_mc_collection) {
    streamlog_out(ERROR) << "At least one event did not provide "
      << "a MC Collection named " << mc_collection_name << ". " << std::endl;
  }
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
    } else if (abs_pdg == 130 || abs_pdg == 310 || abs_pdg == 2112 ||
               abs_pdg == 3122) {
      tv.n_n_hadrons += 1;
    } else {
      streamlog_out(WARNING) << "An unexpected PDG was found: " << abs_pdg
        << "." << std::endl;
    }
  }
  setIsoLeptonNumbers(event, higgs_collection);

  tv.m_z = z_tlv.mass();
  tv.m_recoil = (Tlv(0, 0, 0, 250) - z_tlv).mass();
  tv.m_vis = (z_tlv + higgs_tlv).mass();
  tv.m_miss = (Tlv(0, 0, 0, 250) - z_tlv - higgs_tlv).mass();
  tv.m_h = higgs_tlv.mass();
  tv.m_h_recoil = (Tlv(0, 0, 0, 250) - higgs_tlv).mass();

  tv.cos_theta_z = cos(z_tlv.theta());
  tv.cos_theta_miss = cos((z_tlv +  higgs_tlv).theta());

  tv.higgs_truth = getHiggsTruth(event);
  tree_->Fill();
}

void MasterThesisRootProcessor::setIsoLeptonNumbers(
    EVENT::LCEvent* event, EVENT::LCCollection* higgs_collection
  ) {
  EVENT::LCCollection* lepton_collection = nullptr;
  try {
    lepton_collection = event->getCollection("ISOLeptons");
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "RP collection " << "ISOLeptons"
      << " is not available! Remember calling the IsoLeptonTagging "
      "Processor before this one." << std::endl;
    throw marlin::StopProcessingException(this);
  }
  // principleThrustAxis.x(), .y(), .z(), majorThrustAxis, minorThrustAxis
  FloatVec principle_thrust_axis;
  higgs_collection->getParameters().getFloatVals(
      "principleThrustAxis", principle_thrust_axis);
  tv.principle_thrust_z = principle_thrust_axis[2];
  tv.principle_thrust = higgs_collection->getParameters().getFloatVal("principleThrustValue");
  tv.major_thrust = higgs_collection->getParameters().getFloatVal("majorThrustValue");
  tv.minor_thrust = higgs_collection->getParameters().getFloatVal("minorThrustValue");
  tv.oblateness = higgs_collection->getParameters().getFloatVal("Oblateness");

  tv.sphericity = higgs_collection->getParameters().getFloatVal("sphericity");
  tv.aplanarity = higgs_collection->getParameters().getFloatVal("aplanarity");

  IntVec tagged_lepton_types;
  lepton_collection->getParameters().getIntVals(
      "ISOLepType", tagged_lepton_types);

  RP* highest_e_iso_lepton = nullptr;
  double highest_e = 0;
  for (int e = 0; e < lepton_collection->getNumberOfElements(); ++e) {
    RP* pfo_iso_lepton = static_cast<RP*>(lepton_collection->getElementAt(e));
    double pfo_e = pfo_iso_lepton->getEnergy();
    for (int h = 0; h < higgs_collection->getNumberOfElements(); ++h) {
      RP* higgs_remnant = static_cast<RP*>(higgs_collection->getElementAt(h));
      if (pfo_e == higgs_remnant->getEnergy()) {
        // if (tagged_lepton_types[e] == 11) {
        //   tv.n_iso_electrons++;
        // } else if (tagged_lepton_types[e] == 13) {
        //   tv.n_iso_muons++;
        // } else {
        //   streamlog_out(ERROR) << "Unexpected ISOLepType"
        //     << tagged_lepton_types[e] << std::endl;
        // }
        if (pfo_e > highest_e) {
          highest_e = pfo_iso_lepton->getEnergy();
          highest_e_iso_lepton = pfo_iso_lepton;
        }
        tv.n_iso_leptons++;
        break;
      }
    }
  }
  if (highest_e > 0) {
    tv.e_highest_iso_lep = highest_e;
    double cos_theta = cos(ref_util::getTlv(highest_e_iso_lepton).theta());
    tv.cos_theta_iso_lep = cos_theta;
  }
}


bool isHiggsToSameParticlePair(std::vector<MCP*> remnants) {
  if (remnants.size() != 2) return false;
  return fabs(remnants[0]->getPDG()) == fabs(remnants[1]->getPDG());
}

bool isHiggsToZGamma(std::vector<MCP*> r) {
  if (r.size() != 2) return false;
  if ((fabs(r[0]->getPDG()) == 22)  & (fabs(r[1]->getPDG()) == 23)) return true;
  if ((fabs(r[0]->getPDG()) == 23)  & (fabs(r[1]->getPDG()) == 22)) return true;
  return false;
}

MasterThesisRootProcessor::HiggsTruth MasterThesisRootProcessor::getHiggsTruth(
    EVENT::LCEvent* event
  ) {
  HiggsTruth higgs_info;
  EVENT::LCCollection* mc_collection = nullptr;
  try {
    mc_collection = event->getCollection(mc_collection_name);
  } catch (DataNotAvailableException &e) {
    missing_mc_collection = true;
    return higgs_info;
  }
  for (int i = 0; i < mc_collection->getNumberOfElements(); ++i) {
    MCP* mcp = static_cast<MCP*>(mc_collection->getElementAt(i));
    if (mcp->getPDG() != 25) continue;

    std::vector<MCP*> remnants = mcp->getDaughters();
    // A Higgs can appear as 'intermediate' decay product (hadronization?).
    if (remnants[0]->getPDG() == 25) continue;

    if (isHiggsToSameParticlePair(remnants)) {
      higgs_info.decay_mode = fabs(remnants[0]->getPDG());
    } else if (isHiggsToZGamma(remnants)) {
      higgs_info.decay_mode = 20;
    } else {
      streamlog_out(ERROR) << "An unforeseen Higgs decay occurred. "
        << "The decay prodcuts are: ";
      for (auto r : remnants) streamlog_out(ERROR) << r->getPDG() << ", ";
      streamlog_out(ERROR) << "." << std::endl;
    }

    // Check wether the decay is invisible.
    std::vector<MCP*> stable_higgs_remnants;
    while (!remnants.empty()) {
      MCP* decay_product = remnants.back();
      remnants.pop_back();
      if (decay_product->getGeneratorStatus() == 1) {
        stable_higgs_remnants.push_back(decay_product);
      } else {
        for (auto daughter : decay_product->getDaughters()) {
          remnants.push_back(daughter);
        }
      }
    }
    higgs_info.decays_invisible = true;
    for (auto stable_mcp : stable_higgs_remnants) {
      int abs_pdg = fabs(stable_mcp->getPDG());
      if ((abs_pdg != 12) & (abs_pdg != 14) & (abs_pdg != 16)) {
        higgs_info.decays_invisible = false;
        break;
      }
    }
  }
  return higgs_info;
}