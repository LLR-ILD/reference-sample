/**
*    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
*/
// -- C++ STL headers.

// -- ROOT headers.

// -- LCIO headers.
#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
#include "IMPL/LCCollectionVec.h"

// -- Marlin headers.
#include "marlin/Exceptions.h"

// -- Header for this processor and other project-specific headers.
#include "ref_util.h"
#include "fake_znunu_from_zqq_processor.h"

// -- Using-declarations and global constants.
// Only in .cc files, never in .h header files!
using MCP = EVENT::MCParticle;
using RP = EVENT::ReconstructedParticle;

// This line allows to register your processor in marlin when calling
// "Marlin steering_file.xml".
FakeZnunuFromZqqProcessor aFakeZnunuFromZqqProcessor;

// ----------------------------------------------------------------------------
FakeZnunuFromZqqProcessor::FakeZnunuFromZqqProcessor() :
    marlin::Processor("FakeZnunuFromZqqProcessor") {
  // _description comes from marlin::Processor.
  _description = "A processor to transform Pqqh into Pnnh events.";

  registerInputCollection(
    LCIO::RECONSTRUCTEDPARTICLE,
    "PfoCollection",
    "The Pandora PFO collection name.",
    full_pfo_collection_name_,
    std::string("PandoraPFOs"));

  registerInputCollection(
    LCIO::MCPARTICLE,
    "MCParticleCollection",
    "MCParticle collection used to identify the remnants from the recoiling"
    "Z -> qq.",
    mc_collection_name_,
    std::string("MCParticlesSkimmed"));

  registerInputCollection(
    LCIO::LCRELATION,
    "RelationCollection",
    "Collection of relations that link the Monte Carlo particles with the "
    "reconstructed particles.",
    relation_collection_name_,
    std::string("RecoMCTruthLink"));

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
}

// ----------------------------------------------------------------------------
void FakeZnunuFromZqqProcessor::processEvent(EVENT::LCEvent* event) {
  streamlog_out(DEBUG) << "Processing event no " << event->getEventNumber()
    << " - run " << event->getEventNumber() << std::endl;
  // Prepare all those objects that are used for all Z decay modes.
  EVENT::LCCollection* full_collection = nullptr;
  try {
    full_collection = event->getCollection(full_pfo_collection_name_);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "RP collection " << full_pfo_collection_name_
      << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }

  std::vector<RP*> zqq_remnants;
  std::vector<RP*> overlay_particles;  // TODO: Actually fill the overlay collection  as soon as it is done in SplitOffZProcessor.

  fillQQRemnants(zqq_remnants, event);

  // Fill the new (sub-)collections.
  LCCollectionVec* higgs_remnants_vec = new LCCollectionVec(
      LCIO::RECONSTRUCTEDPARTICLE);
  LCCollectionVec* z_remnants_vec = new LCCollectionVec(
      LCIO::RECONSTRUCTEDPARTICLE);
  LCCollectionVec* overlay_vec = new LCCollectionVec(
      LCIO::RECONSTRUCTEDPARTICLE);
  z_remnants_vec->setSubset(true);
  higgs_remnants_vec->setSubset(true);
  overlay_vec->setSubset(true);
  for (int i = 0; i < full_collection->getNumberOfElements(); ++i) {
    RP* particle = static_cast<RP*>(full_collection->getElementAt(i));
    if (std::find(zqq_remnants.begin(), zqq_remnants.end(), particle)
        != zqq_remnants.end()) {
      // The Zqq information should be thrown away to imitate Znunu.
      //// z_remnants_vec->addElement(particle);
      continue;
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
void FakeZnunuFromZqqProcessor::fillQQRemnants(
    std::vector<RP*> &zqq_remnants, EVENT::LCEvent* event) {
  bool has_higgs_in_event = false;
  int quark_pdg = 0;

  EVENT::LCCollection* mc_collection = nullptr;
  try {
    mc_collection = event->getCollection(mc_collection_name_);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "MC collection " << mc_collection_name_
      << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }

  UTIL::LCRelationNavigator* rel_nav = nullptr;
  try {
    EVENT::LCCollection* relation_collection = event->getCollection(
        relation_collection_name_);
    rel_nav = new UTIL::LCRelationNavigator(relation_collection);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "The relation collection "
      << relation_collection_name_ << " is not available!" << std::endl;
    throw marlin::SkipEventException(this);
  }

  for (int e = 0; e < mc_collection->getNumberOfElements(); ++e) {
   // Find a q from Z decay.
    MCP* mc_particle = static_cast <MCP*>(mc_collection->getElementAt(e));
    if (mc_particle->getParents().size() != 0) continue;

    int pdg = mc_particle->getPDG();
    if (pdg == 25) {
      has_higgs_in_event = true;
      continue;
    }
    if (abs(pdg) > 6) continue;  // Not a quark.

    if (quark_pdg == 0) { // Found the first quark.
      quark_pdg = pdg;
    } else if (quark_pdg == 99) {  // Found more than 2 quarks.
      streamlog_out(DEBUG) << "There is more quarks than the required pair "
        << "without parents in the event. We expect some low-pT quarks from "
        << "gamma-gamma overlay in some of the events." << std::endl;
      continue;
    } else if (quark_pdg != -1*pdg) {  // Found 2nd quark; does not fit to 1st.
      streamlog_out(ERROR) << "Same-charge quark pair!" << std::endl;
      throw marlin::SkipEventException(this);
    } else {  // Found 2nd quark; does fit to 1st.
      quark_pdg = 99;
    }

   // Find the stable MC daughters of a q.
    std::vector<MCP*> mc_q_remnants = mc_particle->getDaughters();
    std::vector<MCP*> stable_mc_q_remnants;
    while (!mc_q_remnants.empty()) {
      MCP* decay_product = mc_q_remnants.back();
      mc_q_remnants.pop_back();
      if (decay_product->getGeneratorStatus() < 2) {
        stable_mc_q_remnants.push_back(decay_product);
      }
      for (auto daughter : decay_product->getDaughters()) {
        mc_q_remnants.push_back(daughter);
      }
    }
    for (auto mcp_from_q : stable_mc_q_remnants) {
      for (auto object_from_q : rel_nav->getRelatedFromObjects(mcp_from_q)) {
        RP* rp_from_q = static_cast<RP*>(object_from_q);
        if (std::find(zqq_remnants.begin(), zqq_remnants.end(), rp_from_q)
            == zqq_remnants.end()) {
          zqq_remnants.push_back(rp_from_q);
        }
      }
    }
  }

  if (!has_higgs_in_event) {
    streamlog_out(ERROR) << "Event without a Higgs boson. Skip." << std::endl;
    throw marlin::SkipEventException(this);
  }
}