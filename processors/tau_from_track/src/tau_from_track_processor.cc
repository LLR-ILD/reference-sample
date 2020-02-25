/**
*    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
*/
// -- C++ STL headers.
#include <cassert>

// -- ROOT headers.

// -- LCIO headers.
#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Track.h"
#include "EVENT/Vertex.h"

// -- Marlin headers.
#include "marlin/Exceptions.h"

// -- Header for this processor and other project-specific headers.
#include "ref_util.h"
#include "tau_from_track_processor.h"

// -- Using-declarations and global constants.
// Only in .cc files, never in .h header files!
using RP = EVENT::ReconstructedParticle;
using Tlv = ROOT::Math::XYZTVector;

// This line allows to register your processor in marlin when calling
// "Marlin steering_file.xml".
TauFromTrackProcessor aTauFromTrackProcessor;

// ----------------------------------------------------------------------------
TauFromTrackProcessor::TauFromTrackProcessor() :
    marlin::Processor("TauFromTrackProcessor") {
  // _description comes from marlin::Processor.
  _description = "An example processor for ILD analysis. ";

  registerInputCollection(
    LCIO::RECONSTRUCTEDPARTICLE,
    "PfoCollection",
    "The Pandora PFO collection name.",
    pfo_collection_name_,
    std::string("PandoraPFOs"));
}

// ----------------------------------------------------------------------------
void TauFromTrackProcessor::init() {
  // Usually a good idea to print parameters.
  printParameters(); // method from marlin::Processor.
  // This is also the place where a root file would be opened
  // and histograms would be defined.
  std::string root_file_name_ = "LCIO";
  std::string track_tree_name_ = "track";


  TString fnn(root_file_name_.c_str()); fnn+=".root";
  root_file_ = new TFile(fnn,"update");

  track_tree_ = new TTree(track_tree_name_.c_str(), track_tree_name_.c_str());
  track_tree_->Branch(("TanLambda"), &TanLambda_, ("TanLambda/D"));
  track_tree_->Branch(("Z0"), &Z0_, ("Z0/D"));
  track_tree_->Branch(("D0"), &D0_, ("D0/D"));
  track_tree_->Branch(("Type"), &Type_, ("Type/I"));
  track_tree_->Branch(("RadiusOfInnermostHit"), &RadiusOfInnermostHit_, ("RadiusOfInnermostHit/D"));
  track_tree_->Branch(("Phi"), &Phi_, ("Phi/D"));
  track_tree_->Branch(("Omega"), &Omega_, ("Omega/D"));
  track_tree_->Branch(("e_track"), &e_track_, ("e_track/D"));
  track_tree_->Branch(("Ndf"), &Ndf_, ("Ndf/D"));
  track_tree_->Branch(("dEdxError"), &dEdxError_, ("dEdxError/D"));
  track_tree_->Branch(("dEdx"), &dEdx_, ("dEdx/D"));
  track_tree_->Branch(("Chi2"), &Chi2_, ("Chi2/D"));

  track_tree_->Branch(("cov_d0_d0"), &cov_d0_d0_, ("cov_d0_d0/D"));
  track_tree_->Branch(("cov_phi_d0"), &cov_phi_d0_, ("cov_phi_d0/D"));
  track_tree_->Branch(("cov_phi_phi"), &cov_phi_phi_, ("cov_phi_phi/D"));
  track_tree_->Branch(("cov_omega_d0"), &cov_omega_d0_, ("cov_omega_d0/D"));
  track_tree_->Branch(("cov_omega_phi"), &cov_omega_phi_, ("cov_omega_phi/D"));
  track_tree_->Branch(("cov_omega_omega"), &cov_omega_omega_, ("cov_omega_omega/D"));
  track_tree_->Branch(("cov_z0_d0"), &cov_z0_d0_, ("cov_z0_d0/D"));
  track_tree_->Branch(("cov_z0_phi"), &cov_z0_phi_, ("cov_z0_phi/D"));
  track_tree_->Branch(("cov_z0_omega"), &cov_z0_omega_, ("cov_z0_omega/D"));
  track_tree_->Branch(("cov_z0_z0"), &cov_z0_z0_, ("cov_z0_z0/D"));
  track_tree_->Branch(("cov_tanlambda_d0"), &cov_tanlambda_d0_, ("cov_tanlambda_d0/D"));
  track_tree_->Branch(("cov_tanlambda_phi"), &cov_tanlambda_phi_, ("cov_tanlambda_phi/D"));
  track_tree_->Branch(("cov_tanlambda_omega"), &cov_tanlambda_omega_, ("cov_tanlambda_omega/D"));
  track_tree_->Branch(("cov_tanlambda_z0"), &cov_tanlambda_z0_, ("cov_tanlambda_z0/D"));
  track_tree_->Branch(("cov_tanlambda_tanlambda"), &cov_tanlambda_tanlambda_, ("cov_tanlambda_tanlambda/D"));

  n_no_higgs_invisible_ = 0;
}

// ----------------------------------------------------------------------------
void TauFromTrackProcessor::processRunHeader(EVENT::LCRunHeader* run) {
  streamlog_out(MESSAGE) << "Starting run no " << run->getRunNumber()
    << std::endl;
}

// ----------------------------------------------------------------------------
void TauFromTrackProcessor::processEvent(EVENT::LCEvent* event) {
  streamlog_out(DEBUG) << "Processing event no " << event->getEventNumber()
    << " - run " << event->getEventNumber() << std::endl;
  ////if (!ref_util::isHiggsInvisibleEvent(event, "MCParticlesSkimmed")) {
  ////  ++n_no_higgs_invisible_;
  ////  throw marlin::SkipEventException(this);
  ////}
  // Prefer accessing any collections with a try/catch approach.
  EVENT::LCCollection* collection = nullptr;
  try {
    collection = event->getCollection(pfo_collection_name_);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "RP collection " << pfo_collection_name_
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
  for (int e = 0; e < collection->getNumberOfElements(); ++e) {
    // A cautious reading of a particle involves casting it with static_cast
    // (instead of dynamic_cast) and ensuring it is not a nullpointer.
    RP* particle = static_cast<RP*>(collection->getElementAt(e));
    if (false) {
      // Let's only look at tau remnants for a bit.
      int primary_pdg = ref_util::getPrimaryPdg(particle, relation_navigator);
      std::cout << primary_pdg;
      if (primary_pdg != 15 and primary_pdg != -15) continue;
    }
    Tlv tlv = ref_util::getTlv(particle);
    if (tlv.Pt() < 5) {
        continue;
    }
    float charge = particle->getCharge();
    int n_tracks = particle->getTracks().size();
    std::cout << charge <<", " << n_tracks << std::endl;
    ////if ( int(n_tracks - fabs(charge)) % 2 == 0) continue;
////
    ////streamlog_out(MESSAGE)
    ////  << "New particle: " << std::endl
    ////  << "  pdg: " <<particle->getType() << std::endl
    ////  << "  charge: " << particle->getCharge() << std::endl
    ////  << "  pT: " << tlv.Pt()
    ////<< std::endl;
    ////if (particle->getStartVertex() != 0) {
    ////  EVENT::Vertex* vertex = particle->getStartVertex();
    ////  streamlog_out(MESSAGE)
    ////    << "  " << vertex->getPosition()[0] << ", "
    ////    << "  " << vertex->getPosition()[1] << ", "
    ////    << "  " << vertex->getPosition()[2] << ", "
    ////    << "  " << vertex->getPosition()[3] << std::endl
    ////    << "  vertex Prob:  " << vertex->getProbability()
    ////  << std::endl;
    ////}
    for (EVENT::Track* track : particle->getTracks()) {
      TanLambda_ = track->getTanLambda();
      Z0_ = track->getZ0();
      D0_ = track->getD0();
      Type_ = track->getType();
      RadiusOfInnermostHit_ = track->getRadiusOfInnermostHit();
      Phi_ = track->getPhi();
      Omega_ = track->getOmega();
      e_track_ = 3e-4 * 3.5 / Omega_; // 3.5 T magnetic field.
      Ndf_ = track->getNdf();
      dEdxError_ = track->getdEdxError();
      dEdx_ = track->getdEdx();
      Chi2_ = track->getChi2();
      cov_d0_d0_               = track->getCovMatrix()[0];
      cov_phi_d0_              = track->getCovMatrix()[1];
      cov_phi_phi_             = track->getCovMatrix()[2];
      cov_omega_d0_            = track->getCovMatrix()[3];
      cov_omega_phi_           = track->getCovMatrix()[4];
      cov_omega_omega_         = track->getCovMatrix()[5];
      cov_z0_d0_               = track->getCovMatrix()[6];
      cov_z0_phi_              = track->getCovMatrix()[7];
      cov_z0_omega_            = track->getCovMatrix()[8];
      cov_z0_z0_               = track->getCovMatrix()[9];
      cov_tanlambda_d0_        = track->getCovMatrix()[10];
      cov_tanlambda_phi_       = track->getCovMatrix()[11];
      cov_tanlambda_omega_     = track->getCovMatrix()[12];
      cov_tanlambda_z0_        = track->getCovMatrix()[13];
      cov_tanlambda_tanlambda_ = track->getCovMatrix()[14];
      track_tree_->Fill();
      streamlog_out(MESSAGE)
        << "  New track: " << std::endl
        << "    TanLambda_: " << TanLambda_ << std::endl
        << "    Z0_: " << Z0_ << std::endl
        << "    D0_: " << D0_ << std::endl
        << "    Type_: " << Type_ << std::endl
        << "    RadiusOfInnermostHit_: " << RadiusOfInnermostHit_ << std::endl
        << "    Phi_: " << Phi_ << std::endl
        << "    Omega_: " << Omega_ << std::endl
        << "    Ndf_: " << Ndf_ << std::endl
        << "    dEdxError_: " << dEdxError_ << std::endl
        << "    dEdx_: " << dEdx_ << std::endl
        << "    Chi2_: " << Chi2_ << std::endl
        << "    cov {d0, phi0, Omega, z0, tanLambda}: " << std::endl
          << "          " << cov_d0_d0_ << std::endl
          << "          " << cov_phi_d0_ << " " << cov_phi_phi_ << std::endl
          << "          " << cov_omega_d0_ << " " << cov_omega_phi_ << " " << cov_omega_omega_ << std::endl
          << "          " << cov_z0_d0_ << " " << cov_z0_phi_ << " " << cov_z0_omega_ << " " << cov_z0_z0_ << std::endl
          << "          " << cov_tanlambda_d0_ << " " << cov_tanlambda_phi_ << " " << cov_tanlambda_omega_ << " " << cov_tanlambda_z0_ << " " << cov_tanlambda_tanlambda_ << std::endl;
    }
  }
}

// ----------------------------------------------------------------------------
void TauFromTrackProcessor::end() {
  // Write the histograms.
  root_file_->cd();
  root_file_->Write(0);
  root_file_->Close();
  // Cleanup your mess here!
  streamlog_out(MESSAGE) << n_no_higgs_invisible_ << " events were skipped for "
    << "not having an invisibly decaying Higgs boson." << std::endl;
}