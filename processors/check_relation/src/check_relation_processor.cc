/**
*    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
*/
// -- C++ STL headers.

// -- ROOT headers.

// -- LCIO headers.
#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include "UTIL/LCRelationNavigator.h"

// -- Marlin headers.
#include "marlin/Exceptions.h"

// -- Header for this processor and other project-specific headers.
#include "check_relation_processor.h"

// -- Using-declarations and global constants.
// Only in .cc files, never in .h header files!
using RP = EVENT::ReconstructedParticle;
using MCP = EVENT::MCParticle;

// This line allows to register your processor in marlin when calling
// "Marlin steering_file.xml".
CheckRelationProcessor aCheckRelationProcessor;

// ----------------------------------------------------------------------------
CheckRelationProcessor::CheckRelationProcessor() :
    marlin::Processor("CheckRelationProcessor") {
  // _description comes from marlin::Processor.
  _description = "An example processor for ILD analysis. ";

  registerInputCollection(
	LCIO::MCPARTICLE,
    "MCParticleCollection",
    "This MCParticle collection should link to the specified ."
      "ReconstructedParticle collection.",
    mc_collection_name_,
    std::string("MCParticlesSkimmed"));

  registerInputCollection(
	LCIO::RECONSTRUCTEDPARTICLE,
    "PfoCollection",
    "This ReconstructedParticle collection should link to the specified ."
      "MCParticle collection.",
    pfo_collection_name_,
    std::string("PandoraPFOs"));

  registerInputCollection(
	LCIO::LCRELATION,
	"RelationCollection",
    "Collection of relations that link the Monte Carlo particles with the "
    "reconstructed particles.",
	relation_collection_name_,
    std::string("RecoMCTruthLink"));

}

// ----------------------------------------------------------------------------
void CheckRelationProcessor::init() {
  // Usually a good idea to print parameters.
  printParameters(); // method from marlin::Processor.
}

// ----------------------------------------------------------------------------
void CheckRelationProcessor::processRunHeader(EVENT::LCRunHeader* run) {
  streamlog_out(MESSAGE) << "Starting run no " << run->getRunNumber()
    << std::endl;
}

// ----------------------------------------------------------------------------
void CheckRelationProcessor::processEvent(EVENT::LCEvent* event) {
  streamlog_out(DEBUG) << "Processing event no " << event->getEventNumber()
    << " - run " << event->getEventNumber() << std::endl;
  // Get the relation collection up and running.
  UTIL::LCRelationNavigator* relation_navigator = nullptr;
  try {
    EVENT::LCCollection* relation_collection = event->getCollection(
        relation_collection_name_);
    relation_navigator = new UTIL::LCRelationNavigator(relation_collection);
    streamlog_out(DEBUG) << relation_navigator->getFromType()
      << std::endl; // ReconstructedParticle.
    streamlog_out(DEBUG) << relation_navigator->getToType()
      << std::endl; // MCParticle.
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "The relation collection "
      << relation_collection_name_ << " does not exist in event no "
      << event->getEventNumber() << "!" << std::endl;
      throw marlin::SkipEventException(this);
  }

  // Check how the relations interplay with the reconstructed particles.
  EVENT::LCCollection* rp_collection = nullptr;
  try {
    rp_collection = event->getCollection(pfo_collection_name_);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "RP collection " << pfo_collection_name_
      << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }
  for (int e = 0; e < rp_collection->getNumberOfElements(); ++e) {
    // Information from relation when starting with Pfo/reconstructed particle.
    RP* pfo_particle = static_cast<RP*>(rp_collection->getElementAt(e));
    if (relation_navigator->getRelatedToObjects(pfo_particle).size() > 1) {
      streamlog_out(MESSAGE) << "Reconstructed particle # " << e
        << ", Pandora PID: " << pfo_particle->getType() << std::endl;
      for (size_t i = 0;
           i < relation_navigator->getRelatedToObjects(pfo_particle).size();
           ++i) {
        MCP* mc_particle = static_cast<MCP*>(
            relation_navigator->getRelatedToObjects(pfo_particle)[i]);
        streamlog_out(MESSAGE) << mc_particle << "  " << mc_particle->getPDG();
      }
      streamlog_out(MESSAGE) << "relation weight: " << std::endl;
      for (auto relation_weight :
           relation_navigator->getRelatedToWeights(pfo_particle)) {
        streamlog_out(MESSAGE) << "  " << relation_weight;
      }
      streamlog_out(MESSAGE) << std::endl;
    }
    // Information from relation when starting with Monte Carlo particle.
    MCP* mc_from_pfo = static_cast <MCP*>(
        relation_navigator->getRelatedToObjects(pfo_particle)[0]);
    if (relation_navigator->getRelatedFromObjects(mc_from_pfo).size() > 1) {
      streamlog_out(MESSAGE) << "MC Particle: " << mc_from_pfo
        << " ,genStatus: " << mc_from_pfo->getGeneratorStatus()
        << " ,sim?: " << mc_from_pfo->isCreatedInSimulation() << std::endl;
      streamlog_out(MESSAGE) << " PDG: " << mc_from_pfo->getPDG()
        << " Energy: " << mc_from_pfo->getEnergy() << std::endl;
      for (size_t i = 0;
           i < relation_navigator->getRelatedFromObjects(mc_from_pfo).size();
           ++i) {
        RP* rp = static_cast<RP*>(
            relation_navigator->getRelatedFromObjects(mc_from_pfo)[i]);
        streamlog_out(MESSAGE) << "PFO, type: " << rp << "  "
          << rp->getType() << ", PFO Energy: " << rp->getEnergy() << std::endl;
      }
      streamlog_out(MESSAGE) << "relation weight: " << std::endl;
      for (auto relation_weight :
           relation_navigator->getRelatedFromWeights(mc_from_pfo)) {
        streamlog_out(MESSAGE) << "  " << relation_weight;
      }
      streamlog_out(MESSAGE) << std::endl;
      }
    if (relation_navigator->getRelatedToObjects(pfo_particle).size() == 0) {
      streamlog_out(MESSAGE) << "No MC particle related to PFO Particle # "
        << e << ", Pandora PID: " << pfo_particle->getType()
        << ", Energy: " << pfo_particle->getEnergy() << std::endl;
    }
  }
}

// ----------------------------------------------------------------------------
void CheckRelationProcessor::end() {
  streamlog_out(MESSAGE) << "end" << std::endl;
}