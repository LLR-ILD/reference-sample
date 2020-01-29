/**
*    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
*/
// -- C++ STL headers.

// -- ROOT headers.

// -- LCIO headers.
#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/LCCollectionVec.h"

// -- Marlin headers.
#include "marlin/Exceptions.h"

// -- Header for this processor and other project-specific headers.
#include "ref_util.h"
#include "split_off_z_processor.h"

// -- Using-declarations and global constants.
// Only in .cc files, never in .h header files!
using RP = EVENT::ReconstructedParticle;
using Tlv = ROOT::Math::XYZTVector;
const float kMZ = 91.19;

// This line allows to register your processor in marlin when calling
// "Marlin steering_file.xml".
SplitOffZProcessor aSplitOffZProcessor;

// ----------------------------------------------------------------------------
SplitOffZProcessor::SplitOffZProcessor() :
    marlin::Processor("SplitOffZProcessor") {
  // _description comes from marlin::Processor.
  _description = "An example processor for ILD analysis. ";

  registerInputCollection(
    LCIO::RECONSTRUCTEDPARTICLE,
    "PfoCollection",
    "The Pandora PFO collection name.",
    full_pfo_collection_name_,
    std::string("PandoraPFOs"));

  registerOutputCollection(
    LCIO::RECONSTRUCTEDPARTICLE,
	"HiggsCollection",
    "Name of the new PFO collection with (only) the Higgs decay remnants.",
    higgs_only_collection_name_,
    std::string("HRemnantsPFOs"));

  registerOutputCollection(
    LCIO::RECONSTRUCTEDPARTICLE,
	"ZCollection",
	"Name of the new PFO collection with (only) the Z decay remnants.",
    z_only_collection_name_,
    std::string("ZRemnantsPFOs"));

  registerOutputCollection(
    LCIO::RECONSTRUCTEDPARTICLE,
	"OverlayCollection",
	"Name of the new PFO collection with (only) the overlay.",
    overlay_collection_name_,
    std::string("OverlayPFOs"));

  registerProcessorParameter(
    "ZDecayChannel",
    "Can be ZDec_EL, ZDec_MU, ZDec_TAU or ZDec_NU.",
    z_decay_channel_,
    std::string("ZDec_MU"));

  registerProcessorParameter(
    "photonRecombinationCos",
    "Cosine of the angle up to which photons should be recombined with an "
      "isolated lepton for the Z reconstruction.",
    photon_recombination_cos_,
    0.99f);

  registerProcessorParameter(
    "OutputRootFile",
    "Name of the output root file.",
    out_root_filename_,
    std::string("split_of_z"));
}

// ----------------------------------------------------------------------------
void SplitOffZProcessor::init() {
  // Usually a good idea to print parameters.
  printParameters(); // method from marlin::Processor.
  // This is also the place where a root file would be opened
  // and histograms would be defined.
  TString fnn(out_root_filename_.c_str());
  fnn += ".root";
  root_out_ = new TFile(fnn, "recreate");
  h_best_z_mass_ = new TH1F(
      ("h_best_z_mass_"+z_decay_channel_).c_str(),
      ("h_best_z_mass_"+z_decay_channel_).c_str(),
      140, 0, 140);
  h_two_best_z_masses_ = new TH2F(
      ("h_two_best_z_masses_"+z_decay_channel_).c_str(),
      ("h_two_best_z_masses_"+z_decay_channel_).c_str(),
      140, 0, 140, 140, 0, 140);
}

// ----------------------------------------------------------------------------
void SplitOffZProcessor::processRunHeader(EVENT::LCRunHeader* run) {
  streamlog_out(MESSAGE) << "Starting run no " << run->getRunNumber()
    << std::endl;
}

// ----------------------------------------------------------------------------
void SplitOffZProcessor::processEvent(EVENT::LCEvent* event) {
  streamlog_out(DEBUG) << "Processing event no " << event->getEventNumber()
    << " - run " << event->getEventNumber() << std::endl;
  // Prepare all those objects that are used for all Z decay modes.
  EVENT::LCCollection* full_collection = nullptr;
  try {
    full_collection = event->getCollection(full_pfo_collection_name_);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "RP collection" << full_pfo_collection_name_
      << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }
  // Some objects that will be altered/ filled depending on the Z decay mode.
  std::vector<RP*> chosen_plus;
  std::vector<RP*> chosen_minus;
  std::vector<RP*> overlay_particles;  // TODO: Actually fill the overlay collection.
  IntVec plus_minus_indices{-1, -1};
  FloatVec second_best_z_mass{-1.};
  float best_z_mass{-1.};

  // Find the Z boson remants, dependent on the chosen Z decay mode.
  if (z_decay_channel_ == "ZDec_NU") {
    streamlog_out(DEBUG) << z_decay_channel_ << ", Z -> nu nu_bar"
      << event->getEventNumber() << std::endl;
  } else if ((z_decay_channel_ == "ZDec_EL") || z_decay_channel_ == "ZDec_MU"
           || z_decay_channel_ == "ZDec_TAU") {
    // Common steps for all charged lepton channels.
    std::vector<Tlv> plus_candidate_momenta;
    std::vector<Tlv> minus_candidate_momenta;
    // plus_candidate_rps and minus_candidate_rps:
    //   It is possible that there are multiple PFO/MC objects that this
    //   candidate is build up from (e.g. electron+photon, or 3pi for tau...).
    std::vector< std::vector<RP*> > plus_candidate_rps;
    std::vector< std::vector<RP*> > minus_candidate_rps;
    if (z_decay_channel_ == "ZDec_TAU") {
      EVENT::LCCollection* tau_collection = nullptr;
      try {
        tau_collection = event->getCollection("TauJets");
      } catch (DataNotAvailableException &e) {
        streamlog_out(ERROR) << "RP collection" << "TauJets"
          << " is not available! Remember calling the TauFinder Processor "
          "before this one." << std::endl;
        throw marlin::StopProcessingException(this);
      }
      UTIL::LCRelationNavigator* tau_navigator = nullptr;
      EVENT::LCCollection* tau_relation_collection = nullptr;
      try {
        tau_relation_collection = event->getCollection("TauRelation");
        tau_navigator = new UTIL::LCRelationNavigator(tau_relation_collection);
      } catch (DataNotAvailableException &e) {
        streamlog_out(ERROR) << "The relation collection " << "TauRelation"
          << " is not available!" << std::endl;
        throw marlin::StopProcessingException(this);
      }
      for (int e = 0; e < tau_collection->getNumberOfElements(); ++e) {
        RP* pfo_tau = static_cast<RP*>(
            tau_collection->getElementAt(e));
        if (pfo_tau->getCharge() > 0) {
          plus_candidate_momenta.push_back(ref_util::getTlv(pfo_tau));
          std::vector<RP*> candidate_vector;
          for (auto tau_object : tau_navigator->getRelatedToObjects(pfo_tau)) {
            RP* tau_part = static_cast<RP*>(tau_object);
            candidate_vector.push_back(tau_part);
          }
          plus_candidate_rps.push_back (candidate_vector);
        } else if (pfo_tau->getCharge() < 0) {
          minus_candidate_momenta.push_back(ref_util::getTlv(pfo_tau));
          std::vector<RP*> candidate_vector;
          for (auto tau_object : tau_navigator->getRelatedToObjects(pfo_tau)) {
            RP* tau_part = static_cast<RP*>(tau_object);
            candidate_vector.push_back(tau_part);
          }
          minus_candidate_rps.push_back(candidate_vector);
        }
        else {
          streamlog_out(ERROR) << "Charge zero was not expected. How does an "
          "isolated lepton object not have a charge??" << std::endl;
        }
      }
    } else { // z_decay_channel_ == "ZDec_EL" or "ZDec_MU"
      EVENT::LCCollection* lepton_collection = nullptr;
      try {
        lepton_collection = event->getCollection("ISOLeptons");
      } catch (DataNotAvailableException &e) {
        streamlog_out(ERROR) << "RP collection" << "ISOLeptons"
          << " is not available! Remember calling the IsoLeptonTagging "
          "Processor before this one." << std::endl;
        throw marlin::StopProcessingException(this);
      }
      IntVec tagged_lepton_types;
      lepton_collection->getParameters().getIntVals(
          "ISOLepType", tagged_lepton_types);
      int lepton_type{0};
      if (z_decay_channel_ == "ZDec_EL") lepton_type = 11;
      if (z_decay_channel_ == "ZDec_MU") lepton_type = 13;
      for (int e = 0; e < lepton_collection->getNumberOfElements(); ++e) {
        RP* pfo_iso_lepton = static_cast<RP*>(
            lepton_collection->getElementAt(e));
        // ISOtypes should always return positive number.
        // Must get the charge from a separate call.
        if (abs(tagged_lepton_types[e] == lepton_type)) {
          if (pfo_iso_lepton->getCharge() > 0) {
            plus_candidate_momenta.push_back(ref_util::getTlv(pfo_iso_lepton));
            // Has to be wrapped in this form for 2 reasons: We might later want
            // to recombine the lepton with a nearby photon, and the tau case
            // will need the possibility to store multiple RPs together as one
            // tau candidate (-> Aim for common interface).
            std::vector<RP*> candidate_vector;
            candidate_vector.push_back(pfo_iso_lepton);
            plus_candidate_rps.push_back (candidate_vector);
          } else if (pfo_iso_lepton->getCharge() < 0) {
            minus_candidate_momenta.push_back(ref_util::getTlv(pfo_iso_lepton));
            std::vector<RP*> candidate_vector;
            candidate_vector.push_back(pfo_iso_lepton);
            minus_candidate_rps.push_back(candidate_vector);
          }
          else {
            streamlog_out(ERROR) << "Charge zero was not expected. How does an "
            "isolated lepton object not have a charge??" << std::endl;
          }
        }
      }

      // Recombine nearby photons with the isolated leptons (FSR,
      // Bremsstrahlung).
      // Theoretically a photon could be counted twice (if it is recombined with
      // both the chosen plus and minus lepton). However, for sensible (small)
      // values of photon_recombination_cos_ this does not happen.
      if (photon_recombination_cos_ < 1.0) {
        for (int i = 0; i < full_collection->getNumberOfElements(); ++i) {
          RP* photon = static_cast<RP*>(full_collection->getElementAt(i));
          if (photon->getType() != 22) continue;
          Tlv photon_tlv = ref_util::getTlv(photon);
          for (size_t j = 0; j < plus_candidate_momenta.size(); ++j) {
            float lepton_photon_angle = ROOT::Math::VectorUtil::CosTheta(
                plus_candidate_momenta[j], photon_tlv);
            if (lepton_photon_angle > photon_recombination_cos_) {
              // Make sure that the photon is not just a lepton that was tagged
              // differently for ISOleptons and Pandora.
              if (plus_candidate_rps[j][0] == photon) continue;
              ////streamlog_out(DEBUG) << "cos:" << lepton_photon_angle
              ////  << "  Energy: " << plus_candidate_rps[j][0]->getEnergy()
              ////  <<        " , " << photon->getEnergy() << std::endl
              ////  << " Type of PFO: " << plus_candidate_rps[j][0]->getType()
              ////  <<  std::endl;
              plus_candidate_momenta[j] += photon_tlv;
              plus_candidate_rps[j].push_back(photon);
            }
          }
          for (size_t j = 0; j < minus_candidate_momenta.size(); ++j) {
            float lepton_photon_angle = ROOT::Math::VectorUtil::CosTheta(
                minus_candidate_momenta[j], photon_tlv);
            if (lepton_photon_angle > photon_recombination_cos_) {
              // Make sure that the photon is not just a lepton that was tagged
              // differently for ISOleptons and Pandora.
              if (minus_candidate_rps[j][0] == photon) continue;
              minus_candidate_momenta[j] += photon_tlv;
              minus_candidate_rps[j].push_back(photon);
            }
          }
        }
      }
      // The actual Z mass calculation step.
      best_z_mass = identifyZ(plus_candidate_momenta, minus_candidate_momenta,
          plus_minus_indices, second_best_z_mass);
      // In case that there is no lepton pair identified, the plus_minus_indices
      // variable takes the value (-1, -1). This happens if the decay mode into
      // neutrinos is specified, or if just no candidate could be found.
      if (plus_minus_indices[0] != -1) {
        chosen_plus = plus_candidate_rps[plus_minus_indices[0]];
        chosen_minus = minus_candidate_rps[plus_minus_indices[1]];
      }
    }
  } else {
    streamlog_out(ERROR) << "The specified decay mode '"
    << z_decay_channel_ << "' is not foreseen as input. Please change it to "
     "ZDec_EL, ZDec_MU, ZDec_TAU or ZDec_NU."  << std::endl;
    throw marlin::StopProcessingException(this);
  }
  // Populate the histograms.
  h_best_z_mass_->Fill(best_z_mass);
  if (second_best_z_mass[0] > 0) {
    h_two_best_z_masses_->Fill(best_z_mass, second_best_z_mass[0]);
  }

  // Fill the new (sub-)collections.
  LCCollectionVec* higgs_remnants_vec = new LCCollectionVec(
      LCIO::RECONSTRUCTEDPARTICLE);
  LCCollectionVec* z_remnants_vec = new LCCollectionVec(
      LCIO::RECONSTRUCTEDPARTICLE);
  LCCollectionVec* overlay_vec = new LCCollectionVec(
      LCIO::RECONSTRUCTEDPARTICLE);
  z_remnants_vec->setSubset(true);
  higgs_remnants_vec->setSubset(true);
  overlay_vec  ->setSubset(true);
  for (int i = 0; i < full_collection->getNumberOfElements(); ++i) {
    RP* particle = static_cast<RP*>(full_collection->getElementAt(i));
    if (std::find(chosen_plus.begin(), chosen_plus.end(), particle)
        != chosen_plus.end()) {
      z_remnants_vec->addElement(particle);
    } else if (std::find(chosen_minus.begin(), chosen_minus.end(), particle)
               != chosen_minus.end()) {
      z_remnants_vec->addElement(particle);
    } else if (std::find(overlay_particles.begin(), overlay_particles.end(),
                         particle)
               != overlay_particles.end()) {
      overlay_vec->addElement(particle);
    } else {
      higgs_remnants_vec->addElement(particle);
    }
  }
  event->addCollection(higgs_remnants_vec, higgs_only_collection_name_.c_str());
  event->addCollection(z_remnants_vec, z_only_collection_name_.c_str());
  event->addCollection(overlay_vec, overlay_collection_name_.c_str());
}

// ----------------------------------------------------------------------------
void SplitOffZProcessor::end() {
  // Write the histograms.
  root_out_->cd();
  root_out_->Write(0);
  root_out_->Close();
  streamlog_out(MESSAGE) << "end" << std::endl;
}

// ----------------------------------------------------------------------------
// Return the index pair for those candidates that form the best Z boson (mass
// closest to the Z boson mass).
// If no pair can be found, return (-1,-1).
float SplitOffZProcessor::identifyZ(
    const std::vector<Tlv> &plus_candidate_momenta,
    const std::vector<Tlv> &minus_candidate_momenta,
    IntVec &plus_minus_indices, FloatVec &second_best_z_mass) {
  // Variable preparations.
  float smallest_diff_to_z = 1000.;
  float best_z_mass = -1;
  ////second_best_z_mass.resize(1);
  ////second_best_z_mass[0] = -1;
  ////plus_minus_indices.resize(2);
  ////plus_minus_indices[0] = -1;
  ////plus_minus_indices[1] = -1;
  // Search for the best Z remant candidates.
  for (size_t i_pl = 0; i_pl < plus_candidate_momenta.size(); ++i_pl) {
    for (size_t i_min = 0; i_min < minus_candidate_momenta.size(); ++i_min) {
      float candidate_z_mass = (
          plus_candidate_momenta[i_pl] + minus_candidate_momenta[i_min]).M();
      if ((abs(candidate_z_mass - kMZ) < smallest_diff_to_z)
          || (smallest_diff_to_z < 0)) {
        smallest_diff_to_z = abs(candidate_z_mass-kMZ);
        second_best_z_mass[0] = best_z_mass;
        best_z_mass = candidate_z_mass;
        plus_minus_indices[0] = i_pl; plus_minus_indices[1] = i_min;
      }
    }
  }
  return best_z_mass;
}