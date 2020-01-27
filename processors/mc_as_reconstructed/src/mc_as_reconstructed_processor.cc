/**
*    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
*/
// -- C++ STL headers.

// -- ROOT headers.

// -- LCIO headers.
#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
  // IMPL headers are needed since we will build new LCIO objects.
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCRelationImpl.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/TrackImpl.h"

// -- Marlin headers.
#include "marlin/Exceptions.h"
#include "marlinutil/GeometryUtil.h"
#include "HelixClass.h"

// -- Header for this processor and other project-specific headers.
#include "mc_as_reconstructed_processor.h"

// -- Using-declarations and global constants.
// Only in .cc files, never in .h header files!
using MCP = EVENT::MCParticle;

// This line allows to register your processor in marlin when calling
// "Marlin steering_file.xml".
McAsReconstructedProcessor aMcAsReconstructedProcessor;

// ----------------------------------------------------------------------------
McAsReconstructedProcessor::McAsReconstructedProcessor() :
    marlin::Processor("McAsReconstructedProcessor") {
  // _description comes from marlin::Processor.
  _description = "Construct a (new) ReconstructedParticle collection from an "
    "input MCParticle collection. ";

  registerInputCollection(
	LCIO::MCPARTICLE,
    "MCInCollection",
    "MCParticle collection from which we want a ReconstructedParticle version.",
    mc_input_collection_name_,
    std::string("MCParticlesSkimmed"));

  registerOutputCollection(
	LCIO::RECONSTRUCTEDPARTICLE,
    "RPOutCollection",
    "Name of the ReconstructedParticles collection created by this processor.",
    now_rp_output_collection_name_,
    std::string("MCParticlesNowRP"));

  registerOutputCollection(
	LCIO::LCRELATION,
	"MCRecLinkCollectionName",
	"The new (trivial) relations collection between the newly created "
      "ReconstructedParticles and their respective underlying MCParticle.",
	mc_to_rp_relation_collection_name_,
	std::string("MCRecLink"));

  registerProcessorParameter(
	"KeepOnlyVisible",
    "Choose that only those MC particles with a type that is visible in the "
      "detector should be kept.",
    keep_only_visible_,
    bool(true));

  registerProcessorParameter(
	"ImitateDetectorAcceptance",
    "Choose that an imitation of the detector acceptance is applied (beam pipe "
      "region rejection, calorimeter thresholds for neutral, pT thresholds for "
      "charged particles).",
    imitate_detector_acceptance_,
    bool(true));

  registerProcessorParameter(
	"MaxCosTheta",
    "Cut on |cos(Theta)| to imitate the forward/backward region without "
      "reliable reconstruction.",
    max_fb_cos_theta_,
    0.995f); // DBD: [0.799, 0.985] region of ETD (p.201)
			 //      <0.99 reach a barrel ECAL layer (p.186)
			 //      < 0.969 to reach at least one tracking layer (p.192)
			 //      < 0.996 reach FTD (forward silicon tracking discs) (p.205)
			 //      [0.997, 0.995] LumiCal
			 //      [0.9992, 0.99998] BeamCal (p.242)

  registerProcessorParameter(
	"chargedMinPT",
    "Cut on the transverse momentum [Gev]. This is applied for charged "
      "particles.",
    min_p_t_for_charged_,
    0.105f); // Motivation: 3.5 T B-field. Beampipe up to ~ 5 cm.
             //     Require 10 cm minimal radius to get some TPC hits.
			 //     p_T[GeV] >= 0.3*|q|*B[T]*R[m] >= 0.3*3.5*0.1 = 0.105.

  registerProcessorParameter(
	"PhotonECALCut",
    "Cut on the energy [Gev] applied to the photons. Above this value we assume"
      " it is possible to detect the photon in the ECAL.",
    photon_ecal_cut_,
    0.2f);

  registerProcessorParameter(
	"NeutralHCALCut",
    "Cut on the energy [Gev] applied to the neutral particles (but the photon)."
      " Above this value we assume it is possible to detect it in the HCAL.",
    neutral_hcal_cut_,
    0.5f);
}

// ----------------------------------------------------------------------------
void McAsReconstructedProcessor::init() {
  // Usually a good idea to print parameters.
  printParameters(); // method from marlin::Processor.
  b_field_ = MarlinUtil::getBzAtOrigin();
  if (b_field_ == 0) {
    streamlog_out(ERROR) <<  "The B-Field was not loaded (b_field_ == 0). "
      "Maybe you forgot to execute the InitD4hep processor (first)?"
      << std::endl;
    throw marlin::StopProcessingException(this);
  }
}

// ----------------------------------------------------------------------------
void McAsReconstructedProcessor::processRunHeader(EVENT::LCRunHeader* run) {
  streamlog_out(MESSAGE) << "Starting run no " << run->getRunNumber()
    << std::endl;
}

// ----------------------------------------------------------------------------
void McAsReconstructedProcessor::processEvent(EVENT::LCEvent* event) {
  streamlog_out(DEBUG) << "Processing event no " << event->getEventNumber()
    << " - run " << event->getEventNumber() << std::endl;
  // Load the original MC collection.
  EVENT::LCCollection* mc_in_collection = nullptr;
  try {
    mc_in_collection = event->getCollection(mc_input_collection_name_);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "The ReconstructedParticle collection "
      << mc_input_collection_name_ << " is not available!" << std::endl;
    throw marlin::StopProcessingException(this);
  }
  // Create and then fill the vectors for the new collections.
  LCCollectionVec *rp_collection = new LCCollectionVec(
      LCIO::RECONSTRUCTEDPARTICLE);
  LCCollectionVec *relation_collection = new LCCollectionVec(LCIO::LCRELATION);
  relation_collection->parameters().setValue(
      std::string("FromType"),LCIO::RECONSTRUCTEDPARTICLE);
  relation_collection->parameters().setValue(
      std::string("ToType"),LCIO::MCPARTICLE);
  for (int i = 0; i < mc_in_collection->getNumberOfElements(); ++i) {
    MCP* mcp = static_cast<MCP*>(mc_in_collection->getElementAt(i));
    if (mcp == nullptr) {
      streamlog_out(ERROR) << "Wrong object type in collection '"
        << mc_input_collection_name_ << "'" << std::endl;
      continue;
    }
	// Only keep the stable MC particles.
	if (mcp->getGeneratorStatus()!=1) continue;
	// If not vetoed by the steering file, remove invisible particles
    // (SM: neutrinos)
	if (keep_only_visible_) {
      int abs_pdg = fabs(mcp->getPDG());
	  if (abs_pdg==12 || abs_pdg==14 || abs_pdg==16 || abs_pdg==1000022) {
	  	continue;
      }
	}
	// If not vetoed by the steering file, remove particles that travel outside
    // of the acceptance region of the detector (both geometrical and
    // energy-wise).
    // The acceptance cuts can be modified from the steering file.
	if (imitate_detector_acceptance_) {
	  const double* mc_mom = mcp->getMomentum();
	  double p_t = sqrt(mc_mom[0]*mc_mom[0] + mc_mom[1]*mc_mom[1]);
	  double p = sqrt(p_t*p_t + mc_mom[2]*mc_mom[2]);
	  double abs_cos_theta = fabs(p_t / p);
	  // Directional cut on all particles.
	  if (abs_cos_theta > max_fb_cos_theta_) continue;
	  // Momentum cut on charged particles to veto centrally spiralling
      // particles.
	  if (mcp->getCharge()) {
          if (p_t < min_p_t_for_charged_) continue;
	  // Energy cut for the photons to be detected in the ECAL.
      } else if (mcp->getPDG() == 22) {
        if (mcp->getEnergy() < photon_ecal_cut_) continue;
	  // Energy cut for the remaining neutral particles to be detected in the
      // HCAL.
      } else {
        if (mcp->getEnergy() < neutral_hcal_cut_) continue;
      }
    }
	// Create an empty reconstructed particle and fill it with information from
    // the Monte Carlo particle.
    ReconstructedParticleImpl *rp_filled_from_mc =
        new ReconstructedParticleImpl();
    rp_filled_from_mc->setMomentum(mcp->getMomentum());
    rp_filled_from_mc->setType    (mcp->getPDG());
    rp_filled_from_mc->setEnergy  (mcp->getEnergy());
    rp_filled_from_mc->setMass    (mcp->getMass());
    rp_filled_from_mc->setCharge  (mcp->getCharge());

	// If charged, add the track. This makes the information for D0 present.
	// The HelixClass is needed for this.
    if (mcp->getCharge()) {
      HelixClass *helix = new HelixClass();
      TrackImpl  *track = new TrackImpl();
	  float vertex[3];
	  float momentum[3];
	  for (int xyz = 0; xyz < 3; ++xyz) {
	    vertex[xyz] = mcp->getVertex()[xyz];
	    momentum[xyz] = mcp->getMomentum()[xyz];
	  }
	  helix->Initialize_VP(vertex, momentum ,mcp->getCharge(), b_field_);
	  const float* reference_point = helix->getReferencePoint();
	  track->setReferencePoint(reference_point);
	  track->setD0(fabs(  helix->getD0()));
	  track->setPhi(      helix->getPhi0());
	  track->setOmega(    helix->getOmega());
	  track->setZ0(       helix->getZ0());
	  track->setTanLambda(helix->getTanLambda());
      // Track is set, the helix did its job.
	  rp_filled_from_mc->addTrack(track);
	  delete helix;
	}
	rp_collection->addElement(rp_filled_from_mc);
    // The relation is trivial: The reconstructed particle is directly related
    // to the Monte Carlo particle that it is created based on.
	LCRelationImpl *relation_from_rp_to_mp = new LCRelationImpl(
        rp_filled_from_mc, mcp);
	relation_collection->addElement(relation_from_rp_to_mp);
  }
  event->addCollection(rp_collection, now_rp_output_collection_name_);
  event->addCollection(relation_collection, mc_to_rp_relation_collection_name_);
}

// ----------------------------------------------------------------------------
void McAsReconstructedProcessor::end() {
  streamlog_out(MESSAGE) << "end" << std::endl;
}