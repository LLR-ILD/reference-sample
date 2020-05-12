/**
*    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
*/
// -- C++ STL headers.

// -- ROOT headers.

// -- LCIO headers.
#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Vertex.h"
#include "UTIL/PIDHandler.h"

// -- Marlin headers.
#include "marlin/Exceptions.h"

// -- Header for this processor and other project-specific headers.
#include "draft_processor.h"
#include "ref_util.h"

// -- Using-declarations and global constants.
// Only in .cc files, never in .h header files!
using RP = EVENT::ReconstructedParticle;
using MCP = EVENT::MCParticle;

// This line allows to register your processor in marlin when calling
// "Marlin steering_file.xml".
DraftProcessor aDraftProcessor;

// ----------------------------------------------------------------------------
DraftProcessor::DraftProcessor() : marlin::Processor("DraftProcessor") {
  // _description comes from marlin::Processor.
  _description = "Processor to quickly draft up new ideas. ";
}

// ----------------------------------------------------------------------------
void DraftProcessor::init() {
  // Usually a good idea to print parameters.
  printParameters(); // method from marlin::Processor.
  // This is also the place where a root file would be opened
  // and histograms would be defined.
  out_root_filename_ = "drafter";
  TString fnn(out_root_filename_.c_str()); fnn+=".root";
  root_out_ = new TFile(fnn,"update");

  h_no_b_b_tags_2jets_  = new TH2F(
      "no_b_b_tags_2jets_", "no_b_b_tags_2jets_",
      22, -0.1, 1, 22, -0.1, 1);
  h_b_b_tags_2jets_    = new TH2F(
      "b_b_tags_2jets_",   "b_b_tags_2jets_",
      22, -0.1, 1, 22, -0.1, 1);
  h_no_c__c_tags_2jets_  = new TH2F(
      "no_c__c_tags_2jets_", "no_c__c_tags_2jets_",
      22, -0.1, 1, 22, -0.1, 1);
  h_cb_c_tags_2jets_   = new TH2F(
      "cb_c_tags_2jets_",  "cb_c_tags_2jets_",
      22, -0.1, 1, 22, -0.1, 1);
  h_c_no_b_c_tags_2jets_ = new TH2F(
      "c_no_b_c_tags_2jets_", "c_no_b_c_tags_2jets_",
      22, -0.1, 1, 22, -0.1, 1);
}

// ----------------------------------------------------------------------------
void DraftProcessor::processRunHeader(EVENT::LCRunHeader* run) {
  streamlog_out(MESSAGE) << "Starting run no " << run->getRunNumber()
    << std::endl;
}

// ----------------------------------------------------------------------------
void DraftProcessor::processEvent(EVENT::LCEvent* event) {
  streamlog_out(DEBUG) << "Processing event no " << event->getEventNumber()
    << " - run " << event->getEventNumber() << std::endl;
  // The function calls below are independent an can be toggled as needed.
  // Each function is a draft on its own.
  //printCollectionInfo(event);
  ////bbVertexPlayground(event);

  std::string mc_col_name = "MCParticlesSkimmed";
  EVENT::LCCollection* mc_collection = event->getCollection(mc_col_name);
  ref_util::printFamilyTree(mc_collection);
}

// ----------------------------------------------------------------------------
void DraftProcessor::end() {
  // Write the histograms.
  root_out_->cd();
  root_out_->Write(0);
  root_out_->Close();
  // Cleanup your mess here!
  streamlog_out(MESSAGE) << "end" << std::endl;
}

// ----------------------------------------------------------------------------
// Print some information about the contents of the .slcio files.
// `streamlog_out(MESSAGE)` is sed as print importance level.
void DraftProcessor::printCollectionInfo(EVENT::LCEvent* event) {
  // Print some information about the PFO and collect one of the vertices for
  // further investigation.
  EVENT::Vertex* vertex = nullptr;
  std::string collection_name = "PandoraPFOs"; //"PrimaryVertex_RP";
  EVENT::LCCollection* pfo_collection = nullptr;
  try {
    pfo_collection = event->getCollection(collection_name);
  } catch (DataNotAvailableException &e) {
    streamlog_out(MESSAGE) << "RP collection " << collection_name
      << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }
  for (int e = 0; e < pfo_collection->getNumberOfElements(); ++e) {
    RP* pf_particle = static_cast<RP*>(pfo_collection->getElementAt(e));
    rpPrint(pf_particle);
    for (auto daughter_of_pf_particle : pf_particle->getParticles()) {
      RP* rp = static_cast<RP*>(daughter_of_pf_particle);
      streamlog_out(MESSAGE) << std::endl;
      rpPrint(rp);
    }
    vertex = pf_particle->getEndVertex();
  }

  // Investigate:
  // Is the entry of PrimaryVertex the EndVertex in PrimaryVertex_RP Collection?
  std::string vertex_collection_name = "RefinedVertex";
  streamlog_out(MESSAGE) << std::endl << "Vertex collection name: "
    << vertex_collection_name << std::endl;
  EVENT::LCCollection* vertex_collection = nullptr;
  try {
    vertex_collection = event->getCollection(vertex_collection_name);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "Vertex collection " << vertex_collection_name
      << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }
  for (int e = 0; e < vertex_collection->getNumberOfElements(); ++e) {
    EVENT::Vertex* vertex_object =
      static_cast<EVENT::Vertex*>(vertex_collection->getElementAt(e));
    streamlog_out(MESSAGE)
      << "  vertex_object->getPosition()          "
      <<    vertex_object->getPosition()[0] << ", "
      <<    vertex_object->getPosition()[1] << ", "
      <<    vertex_object->getPosition()[2] << ", "
      <<    vertex_object->getPosition()[3] << ", "
      <<    vertex_object->getPosition()[4] << ", "
      <<    vertex_object->getPosition()[5] << ", "
      <<    vertex_object->getPosition()[6] << std::endl;
    streamlog_out(MESSAGE)
      << "  (vertex_object == vertex)             "
      <<    (vertex_object == vertex)       << std::endl;
    streamlog_out(MESSAGE)
      << "  vertex_object->isPrimary()            "
      <<    vertex_object->isPrimary()      << std::endl;
    streamlog_out(MESSAGE)
      << "  vertex_object->getAssociatedParticle()"
      <<    vertex_object->getAssociatedParticle() << std::endl;
  }
}

// ----------------------------------------------------------------------------
// LCFIPLus has to be run before this!
// Shows some output of the flavor tagging.
void DraftProcessor::bbVertexPlayground(EVENT::LCEvent* event) {
  streamlog_out(MESSAGE) << "Hello from bbVertexPlayground" << std::endl;
  // The collections used in the function.
  std::string jet_collection_name = "RefinedJets";
  std::string mc_collection_name  = "MCParticlesSkimmed";
  std::string pfo_collection_name = "PandoraPFOs";
  EVENT::LCCollection* jet_collection = nullptr;
  try {
    jet_collection = event->getCollection(jet_collection_name);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "The RP/jet collection " << jet_collection_name
      << " is not available!" << std::endl;
  }
  const int kJetsPerEvent = jet_collection->getNumberOfElements();
  if (kJetsPerEvent == 0) throw marlin::SkipEventException(this);

  // Flavor tagging: Create jet objects consisting of a RP with b-tag, c-tag and
  // an is_assigned flag.
  UTIL::PIDHandler particle_id_handler(jet_collection);
  int algo = -1;
  std::string algo_name = "lcfiplus";
  try {
    algo = particle_id_handler.getAlgorithmID(algo_name);
  } catch (...) {
    streamlog_out(ERROR) << "The algorithm " << algo_name
      << " is not available for this collection in event no "
      << event->getEventNumber() << "." << std::endl;
    }
  std::vector<MyJet*> jet_list;
  for (int i_jet = 0; i_jet < kJetsPerEvent; ++i_jet) {
    RP* rp_jet = static_cast<RP*>(jet_collection->getElementAt(i_jet));
    const ParticleID& jet_id = particle_id_handler.getParticleID(rp_jet , algo);
    FloatVec params = jet_id.getParameters();
    float b_tag = params[particle_id_handler.getParameterIndex(algo , "BTag")];
    float c_tag = params[particle_id_handler.getParameterIndex(algo , "CTag")];
    //float bc_tag = params[particle_id_handler.getParameterIndex(algo ,"BCTag")];
    MyJet* a_jet = new MyJet;
    a_jet->jet   = rp_jet;
    a_jet->b_tag  = b_tag;
    a_jet->c_tag  = c_tag;
    ////a_jet->bc_tag = bc_tag;
    a_jet->is_assigned = false;
    jet_list.push_back(a_jet);
    streamlog_out(MESSAGE) << "kJetsPerEvent: " << kJetsPerEvent
      << ", rp_jet->getEnergy() = " << rp_jet->getEnergy() << std::endl;
  }

  // Evaluation: Print the MC true number of b quarks and the number of b-tagged
  // jets.
  if (true) {
    // Jets:
    int n_b_jets = 0;
    streamlog_out(MESSAGE) << std::endl << "Jet information: " << std::endl;
    for (MyJet* jet : jet_list) {
      streamlog_out(MESSAGE) << jet->b_tag << "  ";
      if (jet->b_tag > 0.6) ++n_b_jets;
    }
    streamlog_out(MESSAGE) << std::endl
      << "  Number of b-tagged jets: " << n_b_jets << std::endl;
    // MC:
    EVENT::LCCollection* mc_collection = nullptr;
    try {
      mc_collection = event->getCollection(mc_collection_name);
    } catch (DataNotAvailableException &e) {
      streamlog_out(ERROR) << "The MC collection " << mc_collection_name
        << " is not available!" << std::endl;
    }
    bool b_quark_is_in_mc = ref_util::pdgIsInMcCol(5, mc_collection);
    streamlog_out(MESSAGE) << "B quarks in the event: "
      << b_quark_is_in_mc << std::endl;
  }

  // Get the b_tag values of the 2 jets in the collection of the event.
  if (jet_list.size() != 2) {
    streamlog_out(ERROR) << "This code expects 2 jets. The collection gives "
      << jet_list.size() << "jets." << std::endl;
  }
  float b_tag1 = -.1;
  float b_tag2 = -.1;
  for (MyJet* jet : jet_list) {
    float jet_b_tag = jet->b_tag;
    // Next if-statement just to test my understanding of LCFIPLus.
    if (jet_b_tag > b_tag1) {
      b_tag2 = b_tag1;
      b_tag1 = jet_b_tag;
    } else if (jet_b_tag > b_tag2) {
      b_tag2 = jet_b_tag;
    } else {
      streamlog_out(MESSAGE) << "Weird. Maybe more than 2 jets in event? "
        << "b_tag1, b_tag2, jet_b_tag = "
        << b_tag1 << b_tag2 << jet_b_tag << std::endl;
    }
  }
  if (b_tag2 < 0) {
    streamlog_out(MESSAGE) << "No 2 jets found in the event." << std::endl;
  }
  // Fill the b_tags into a histogram, depending on wether the event actuallay
  // has a b quark.
  EVENT::LCCollection* mc_collection = nullptr;
  try {
    mc_collection = event->getCollection(mc_collection_name);
  } catch (DataNotAvailableException &e) {
    streamlog_out(MESSAGE) << "The MC collection " << mc_collection_name
      << " is not available!" << std::endl;
  }
  bool b_quark_is_in_mc = ref_util::pdgIsInMcCol(5, mc_collection);
  if (b_quark_is_in_mc) {
    h_b_b_tags_2jets_->Fill(b_tag1, b_tag2);
  } else {
    h_no_b_b_tags_2jets_->Fill(b_tag1, b_tag2);
  }

  // Now let's tackle the c_tags in a similar manner.
    float c_tag1 = -.1;
  float c_tag2 = -.1;
  for (MyJet* jet : jet_list) {
    float jet_c_tag = jet->c_tag;
    // Next if-statement just to test my understanding of LCFIPLus.
    if (jet_c_tag < 0 || jet_c_tag > 1) { // Should never happen.
      streamlog_out(ERROR) << "jet_c_tag: " << jet_c_tag << std::endl;
    }
    if (jet_c_tag > c_tag1) {
      c_tag2 = c_tag1;
      c_tag1 = jet_c_tag;
    } else if (jet_c_tag > c_tag2) {
      c_tag2 = jet_c_tag;
    } else {
      streamlog_out(MESSAGE) << "Weird. Maybe more than 2 jets in event? "
        << "c_tag1, c_tag2, jet_c_tag = "
        << c_tag1 << c_tag2 << jet_c_tag << std::endl;
    }
  }
  if (c_tag2 < 0) {
    streamlog_out(MESSAGE) << "No 2 jets found in the event." << std::endl;
  }
  bool c_quark_is_in_mc = (ref_util::pdgIsInMcCol(4, mc_collection)
    || ref_util::pdgIsInMcCol(-4, mc_collection)); // TODO: Why do I test for both signs here, but not above with the b_quark??
  if (c_quark_is_in_mc) {
    if (b_quark_is_in_mc) {
      h_cb_c_tags_2jets_->Fill(c_tag1, c_tag2);
    } else {
      h_c_no_b_c_tags_2jets_->Fill(c_tag1, c_tag2);
    }
  } else {
    h_no_c__c_tags_2jets_->Fill(c_tag1, c_tag2);
  }

  // Print total energy of the reconstructed particles in the event.
  float pfo_total_energy = 0.;
  EVENT::LCCollection* pfo_collection = nullptr;
  try {
    pfo_collection = event->getCollection(pfo_collection_name );
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "The MC collection " << pfo_collection_name
      << " is not available!" << std::endl;
  }
  for (int i = 0; i < pfo_collection->getNumberOfElements(); ++i) {
    RP* pfo_particle = static_cast<RP*>(pfo_collection->getElementAt(i));
    pfo_total_energy += pfo_particle->getEnergy();
    streamlog_out(MESSAGE) << "type: " << pfo_particle->getType()
      << ", E= " << pfo_particle->getEnergy()
      << ", theta: " << pfo_particle->getMomentum()[2]/pfo_particle->getEnergy()
      << std::endl;
  }
  streamlog_out(MESSAGE) << "Energy sum of the Pandora objects: "
    << pfo_total_energy << std::endl;

  // Print total energy of the Monte Carlo particles in the event.
  float mc_total_energy = 0.;
  for (int i = 0; i < mc_collection->getNumberOfElements(); ++i) {
    MCP* mc_particle = static_cast<MCP*>(mc_collection->getElementAt(i));
    int abs_pdg = fabs(mc_particle->getPDG());
    if ((mc_particle->getGeneratorStatus() == 1) && (abs_pdg != 12)
        && (abs_pdg != 14) && (abs_pdg != 16)) {
      mc_total_energy += mc_particle->getEnergy();
    }
  }
  streamlog_out(MESSAGE) << "Energy sum of the MC objects: "
    << mc_total_energy << std::endl;
  streamlog_out(MESSAGE) << "b_tag1, b_tag2 = "
    << b_tag1 << ", " << b_tag2 << std::endl;
  ref_util::printFamilyTree(mc_collection);
}

// Is (only) used inside bbVertexPlayground.
void DraftProcessor::rpPrint(EVENT::ReconstructedParticle* rp) {
    streamlog_out(MESSAGE)
      << " ReconstructedParticle      " << rp << std::endl;
    streamlog_out(MESSAGE)
      << "  rp->getType()             " << rp->getType() << std::endl;
    streamlog_out(MESSAGE)
      << "  rp->getCharge()           " << rp->getCharge() << std::endl;
    streamlog_out(MESSAGE)
      << "  rp->getEnergy()           " << rp->getEnergy() << std::endl;
    streamlog_out(MESSAGE)
      << "  rp->getMass()             " << rp->getEnergy() << std::endl;
    streamlog_out(MESSAGE)
      << "  rp->getMomentum()         "
      << rp->getMomentum()[0] << ", " << rp->getMomentum()[1] << ", "
      << rp->getMomentum()[2] << ", " << rp->getMomentum()[3] << std::endl;
    streamlog_out(MESSAGE)
      << "  rp->getParticles().size() " << rp->getParticles().size()
      << std::endl;
    streamlog_out(MESSAGE)
      << "  rp->getStartVertex()      " << rp->getStartVertex() << std::endl;
    streamlog_out(MESSAGE)
      << "  rp->getEndVertex()        " << rp->getEndVertex() << std::endl;
}

// ----------------------------------------------------------------------------
/*

500 IDR                                                                      250 DBD
---------------------------------------------------------------------------  ---------------------------------------------------------------------------
COLLECTION NAME               COLLECTION TYPE          NUMBER OF ELEMENTS    COLLECTION NAME               COLLECTION TYPE          NUMBER OF ELEMENTS
===========================================================================  ===========================================================================
BCalClusters                  Cluster                          1             BCALClusters                  Cluster                          0
BCalRecoParticle              ReconstructedParticle            1             BCALParticles                 ReconstructedParticle            0
BuildUpVertex                 Vertex                           1             BuildUpVertex                 Vertex                           0
BuildUpVertex_RP              ReconstructedParticle            1             BuildUpVertex_RP              ReconstructedParticle            0
BuildUpVertex_V0              Vertex                           0             BuildUpVertex_V0              Vertex                           1
BuildUpVertex_V0_RP           ReconstructedParticle            0             BuildUpVertex_V0_RP           ReconstructedParticle            1
ClusterMCTruthLink            LCRelation                     148
DistilledPFOs                 ReconstructedParticle           53
GammaGammaCandidateEtaPrimes  ReconstructedParticle            4
GammaGammaCandidateEtas       ReconstructedParticle            9
GammaGammaCandidatePi0s       ReconstructedParticle           10
GammaGammaParticles           ReconstructedParticle            9
MCParticle                    MCParticle                     496             MCParticlesSkimmed            MCParticle                     141
MCTruthClusterLink            LCRelation                     148
MCTruthMarlinTrkTracksLink    LCRelation                      32             MCTruthMarlinTrkTracksLink    LCRelation                      29
MCTruthRecoLink               LCRelation                     156
MarlinTrkTracks               Track                           29             MarlinTrkTracks               Track                           24
MarlinTrkTracksMCTruthLink    LCRelation                      32             MarlinTrkTracksMCTruthLink    LCRelation                      29
PandoraClusters               Cluster                         56             PandoraClusters               Cluster                         41
                                                                             PandoraPFANewReclusterMonitoringLCGenericObject                0
PandoraPFOs                   ReconstructedParticle           62             PandoraPFOs                   ReconstructedParticle           47
PrimaryVertex                 Vertex                           1             PrimaryVertex                 Vertex                           1
PrimaryVertex_RP              ReconstructedParticle            1             PrimaryVertex_RP              ReconstructedParticle            1
                                                                             ProngRecoParticles            ReconstructedParticle            1
                                                                             ProngVertices                 Vertex                           1
RecoMCTruthLink               LCRelation                     156             RecoMCTruthLink               LCRelation                      47
                                                                             V0RecoParticles               ReconstructedParticle            1
                                                                             V0Vertices                    Vertex                           1

---------------------------------------------------------------------------  --------------------------------------------------------------------
*/