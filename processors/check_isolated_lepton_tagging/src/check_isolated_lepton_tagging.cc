/**
*    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
*/
// -- C++ STL headers.

// -- ROOT headers.

// -- LCIO headers.
#include "EVENT/MCParticle.h"
#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"

// -- Marlin headers.
#include "marlin/Exceptions.h"

// -- Header for this processor and other project-specific headers.
#include "check_isolated_lepton_tagging.h"
#include "ref_util.h"

// -- Using-declarations and global constants.
// Only in .cc files, never in .h header files!
using RP = EVENT::ReconstructedParticle;
using MCP = EVENT::MCParticle;

// This line allows to register your processor in marlin when calling
// "Marlin steering_file.xml".
CheckIsoLeptonProcessor aCheckIsoLeptonProcessor;

// ----------------------------------------------------------------------------
CheckIsoLeptonProcessor::CheckIsoLeptonProcessor() :
    marlin::Processor("CheckIsolatedLeptonTaggingProcessor") {
  // _description comes from marlin::Processor.
  _description = "Some investigations on the IsolatedLeptonTagging processor. ";

  registerInputCollection(
    LCIO::RECONSTRUCTEDPARTICLE,
    "PfoCollection",
    "The Pandora PFO collection name.",
    pfo_collection_name_,
    std::string("PandoraPFOs"));

  registerProcessorParameter(
    "OutputRootFile",
    "Name of the output root file.",
    out_root_filename_,
    std::string("check_iso_lepton_tagging"));
}

// ----------------------------------------------------------------------------
void CheckIsoLeptonProcessor::init() {
  // Usually a good idea to print parameters.
  printParameters(); // method from marlin::Processor.
  // This is also the place where a root file would be opened:
  TString fnn(out_root_filename_.c_str());
  fnn += ".root";
  root_out_ = new TFile(fnn, "recreate");
  // and histograms would be defined.
  h_n_isolated_muons_ = new TH1F(
      "h_n_isolated_muons_", "h_n_isolated_muons_",
      10, -0.5, 9.5);
  h_n_isolated_electrons_ = new TH1F(
      "h_n_isolated_electrons_", "h_n_isolated_electrons_",
      10, -0.5, 9.5);

  h_mva_isolated_muon_ = new TH1F(
      "h_mva_isolated_muon_", "h_mva_isolated_muon_",
      150, 0, 1.5);
  h_mva_isolated_muon_possibly_z_ = new TH1F(
      "h_mva_isolated_muon_possibly_z_", "h_mva_isolated_muon_possibly_z_",
      150, 0, 1.5);
  h_mva_isolated_electron_ = new TH1F(
      "h_mva_isolated_electron_", "h_mva_isolated_electron_",
      150, 0, 1.5);
  h_mva_isolated_electron_possibly_z_ = new TH1F(
      "h_mva_isolated_electron_possibly_z_",
      "h_mva_isolated_electron_possibly_z_",
      150, 0, 1.5);

  h_n_objects_per_type_ = new TH2I(
      "h_n_objects_per_type_", "h_n_objects_per_type_",
      ref_util::kPfoTypes, -.5, ref_util::kPfoTypes-0.5, 20, 0, 19);
  for (size_t i_type = 0; i_type < ref_util::kPfoTypes; ++i_type) {
    h_n_objects_per_type_->GetXaxis()->SetBinLabel(
        i_type+1, ref_util::getPfoLabel(i_type));
  }
}

// ----------------------------------------------------------------------------
void CheckIsoLeptonProcessor::processRunHeader(EVENT::LCRunHeader* run) {
  streamlog_out(MESSAGE) << "Starting run no " << run->getRunNumber()
    << std::endl;
}

// ----------------------------------------------------------------------------
void CheckIsoLeptonProcessor::processEvent(EVENT::LCEvent* event) {
  streamlog_out(DEBUG) << "Processing event no " << event->getEventNumber()
    << " - run " << event->getEventNumber() << std::endl;
  EVENT::LCCollection* collection = nullptr;
  try {
    collection = event->getCollection(pfo_collection_name_);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "RP collection" << pfo_collection_name_
      << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }
  // The relation_navigator is used to find out which RPs might stem from a
  // final state Z boson.
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

  // Collect the particle type occurrances in the event.
  std::vector<int> particles_per_type(ref_util::kPfoTypes);
  for (int e = 0; e < collection->getNumberOfElements(); ++e) {
    RP* particle = static_cast<RP*>(collection->getElementAt(e));
    int type_id = ref_util::particlePdgToName(abs(particle->getType()));
    particles_per_type[type_id]++;
  }
  for (size_t i_type = 0; i_type < ref_util::kPfoTypes; ++i_type) {
    h_n_objects_per_type_->Fill(i_type, particles_per_type[i_type]);
  }

  // If the "ISOLeptons" collection is available, fill a vector for each
  // electrons and muons.
  std::vector<RP*> isolated_muons;
  std::vector<RP*> isolated_electrons;
  EVENT::LCCollection* isolated_lepton_collection = nullptr;
  try {
    isolated_lepton_collection = event->getCollection("ISOLeptons");
  } catch (EVENT::DataNotAvailableException &) {
    streamlog_out(WARNING) << "Pfo collection '" << "ISOLeptons"
      << "' is not available!" << std::endl;
    streamlog_out(WARNING) << "Check that you called the"
      "IsolatedLeptonTaggingProcessor first." << std::endl;
      throw marlin::StopProcessingException(this);
  }
  streamlog_out(DEBUG) << "Number of isolated leptons: "
    << isolated_lepton_collection->getNumberOfElements() << "." << std::endl;
  if (isolated_lepton_collection->getNumberOfElements() > 0) {
    IntVec isolated_lepton_types;
    isolated_lepton_collection->getParameters().getIntVals(
        "ISOLepType", isolated_lepton_types);
    FloatVec isolated_lepton_tags;
    isolated_lepton_collection->getParameters().getFloatVals(
        "ISOLepTagging", isolated_lepton_tags);
    for (int e = 0;
         e < isolated_lepton_collection->getNumberOfElements();
         ++e) {
      RP* isolated_lepton = static_cast<RP*>(
          isolated_lepton_collection->getElementAt(e));
      // Fill some histograms only with those IsoLeptons that are possibly
      // primary Z decay remnants.
      bool might_be_from_z = true;
      for (auto mcp : ref_util::getMcChainFromRp(
          isolated_lepton, relation_navigator)) {
        int mcp_pdg = fabs(mcp->getPDG());
        if (mcp_pdg != 13 and mcp_pdg != 11 and mcp_pdg != 94 ) {
            might_be_from_z = false;
        }
      }
      if (might_be_from_z) {

      }
      if (abs(isolated_lepton_types[e]) == 13) {
        isolated_muons.push_back(isolated_lepton);
        h_mva_isolated_muon_->Fill(isolated_lepton_tags[e]);
        if (might_be_from_z) {
            h_mva_isolated_muon_possibly_z_->Fill(isolated_lepton_tags[e]);
        }
      } else if (abs(isolated_lepton_types[e]) == 11) {
        isolated_electrons.push_back(isolated_lepton);
        h_mva_isolated_electron_->Fill(isolated_lepton_tags[e]);
        if (might_be_from_z) {
            h_mva_isolated_electron_possibly_z_->Fill(isolated_lepton_tags[e]);
        }
      }
    }
  }
  h_n_isolated_muons_->Fill(isolated_muons.size());
  h_n_isolated_electrons_->Fill(isolated_electrons.size());
}

// ----------------------------------------------------------------------------
void CheckIsoLeptonProcessor::end() {
  // Write the histograms.
  root_out_->cd();
  root_out_->Write(0);
  root_out_->Close();
  streamlog_out(MESSAGE) << "end" << std::endl;
}