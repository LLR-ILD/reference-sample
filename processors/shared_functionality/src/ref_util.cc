/**
*    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
*/
// -- C++ STL headers.
#include <cassert>
#include <iomanip>  // Has std::setPrecision().
#include <iostream>  // Has std::fixed.
#include <string>  // Has std::to_string().
#include <sstream>  // Has std::ostringstream.

// -- ROOT headers.

// -- LCIO headers.

// -- Marlin headers.
#include "marlin/Exceptions.h"

// -- Header for this processor and other project-specific headers.
#include "ref_util.h"

// -- Using-declarations and global constants.
// Only in .cc files, never in .h header files!
using RP = EVENT::ReconstructedParticle;
using MCP = EVENT::MCParticle;
using Tlv = ROOT::Math::XYZTVector;

// ----------------------------------------------------------------------------
// Local helper functions -> Unnamed namespace makes them only available in this
// compilation unit. Not part of the ref_util namespace as they are not intended
// for outside use.
namespace {

// Helper for ref_util::printFamilyTree.
// Adds a rather long string with information about the MC particle to the
// stream that is handed to the function by reference.
void particleInfoForFamilyTree(
    std::ostringstream& particle_info, MCP* mcp) {
  particle_info << mcp->getGeneratorStatus() << "  " << mcp->getPDG() << "  "
    << "(E,p)=("
          << std::fixed << std::setprecision(2) << mcp->getEnergy()      << ", "
          << std::fixed << std::setprecision(2) << mcp->getMomentum()[0] << ", "
          << std::fixed << std::setprecision(2) << mcp->getMomentum()[1] << ", "
          << std::fixed << std::setprecision(2) << mcp->getMomentum()[2] << ") "
    ////<<  mcp->id()                  << "  "
    ////<<  mcp->getDaughters().size() << "  "
    ////<<  mcp->getParents().size()   << "  "
    << "(x,y,z)=("
          << std::fixed << std::setprecision(2) << mcp->getVertex()[0] << ", "
          << std::fixed << std::setprecision(2) << mcp->getVertex()[1] << ", "
          << std::fixed << std::setprecision(2) << mcp->getVertex()[2] << ") "
    << "cos(theta)=" << std::fixed << std::setprecision(2)
          << cos(ref_util::getTlv(mcp).Theta()) << " "
    << "phi=" << std::fixed << std::setprecision(2)
          << ref_util::getTlv(mcp).Phi() << " ";
    ////<< "simulated: "            << mcp->isCreatedInSimulation()
    ////<< "status: "               << mcp->getSimulatorStatus()
    ////<< "backscatter: "          << mcp->isBackscatter()
    ////<< "notEndpoint: "          << mcp->vertexIsNotEndpointOfParent()
    ////<< "trackDecay : "          << mcp->isDecayedInTracker()
    ////<< "caloDecay : "           << mcp->isDecayedInCalorimeter()
    ////<< "left detector: "        << mcp->hasLeftDetector()
    ////<< "stopped in detector: "  << mcp->isStopped();
}

// Helper for ref_util::printFamilyTree.
// Takes care of the styling of the tree.
void downOneGenerationLoop(std::ostringstream& particle_info,
    MCP* particle, std::vector<MCP*>& printed_daughters,
    std::string indentation_level, bool last_daughter) {
  particle_info << indentation_level << "|" << std::endl
                << indentation_level << "|--- ";
  particleInfoForFamilyTree(particle_info, particle);
  // Handle some particles having more than one parent particle
  if (particle->getParents().size() > 1) {
    particle_info << std::endl << indentation_level
      << "        This Particle has " << particle->getParents().size()
      << " parents and will therefore appear multiple times in this tree! "
      << "Its descendents are only shown once. ";
    if (std::find(printed_daughters.begin(), printed_daughters.end(), particle)
        != printed_daughters.end()) {
      return;
    }
    printed_daughters.push_back(particle);
  }
  if (last_daughter) {
    indentation_level = indentation_level + "     ";
  } else {
    indentation_level = indentation_level + "|    ";
  }
  particle_info << std::endl;
  const int kDaughters = particle->getDaughters().size();
  for (int i = 0; i < kDaughters; ++i) {
    auto daughter = particle->getDaughters()[i];
    downOneGenerationLoop(particle_info,
        daughter, printed_daughters, indentation_level, (i == kDaughters-1));
  }
  return;
}
}  // namespace

// ----------------------------------------------------------------------------
TString ref_util::getPfoLabel(int code) {
  switch (code) {
    case kPhoton:
      return "GAM";
    case kElectron:
      return "EL";
    case kMuon:
      return "MU";
    case kPiCharged:
      return "PIc";
    case kKShort:
      return "Ks";
    case kNeutron:
      return "N";
    case kLambda:
      return "LAM";
    default: {
      streamlog_out(ERROR) << "Unexpected PFO object with code '"
        << code << "'" << std::endl;
      assert(0);
    }
  }
}

int ref_util::particlePdgToName(int pdg_val) {
  switch (pdg_val) {
    case 11:
      return ref_util::kElectron;
    case 13:
      return ref_util::kMuon;
    case 22:
      return ref_util::kPhoton;
    case 211:
      return ref_util::kPiCharged;
    case 310:
      return ref_util::kKShort;
    case 2112:
      return ref_util::kNeutron;
    case 3122:
      return ref_util::kLambda;
    default:
      streamlog_out(ERROR) << "Unexpected PDG value "
        << std::to_string(pdg_val) << std::endl;
      assert(0);
  }
}

TString ref_util::getPfoLabelFromPdg(int pdg_val) {
  int code = abs(ref_util::particlePdgToName(pdg_val));
  return ref_util::getPfoLabel(code);
}

// ----------------------------------------------------------------------------
Tlv ref_util::getTlv(std::vector<RP*> rps) {
  Tlv total_momentum(0,0,0,0);
  for (auto rp : rps)
    total_momentum += getTlv(rp);
  return total_momentum;
}

Tlv ref_util::getTlv(RP* rp) {
  assert(rp);
  Tlv four_vector = Tlv(rp->getMomentum()[0], rp->getMomentum()[1],
      rp->getMomentum()[2], rp->getEnergy());
  return four_vector;
}

Tlv ref_util::getTlv(MCP* mcp) {
  assert(mcp);
  Tlv four_vector = Tlv(mcp->getMomentum()[0], mcp->getMomentum()[1],
      mcp->getMomentum()[2], mcp->getEnergy());
  return four_vector;
}

Tlv ref_util::getTlv(std::vector<MCP*> mc_particles) {
  Tlv total_momentum(0,0,0,0);
  for (auto mcp : mc_particles)
    total_momentum += getTlv(mcp);
  return total_momentum;
}

// ----------------------------------------------------------------------------
// Rewriting of the python code in angularPlots.py (ZH).
void ref_util::printFamilyTree(EVENT::LCCollection* collection) {
  // Fill this stream with the information that we want to print.
  std::ostringstream particle_info;
  // Some particles can have multiple parents. After printing once, fill these
  // parents into the vector to avoid printing the decay chain twice.
  std::vector<MCP*> printed_daughters;
  for (int i = 0; i < collection->getNumberOfElements(); ++i) {
    MCP* particle = static_cast<MCP*>(collection->getElementAt(i));
    // To avoid printing the same family branch multiple times:
    if (particle->getParents().size() == 0) {
      std::string indentation_level = "";
      particle_info << std::endl << "|" << std::endl << "|" << std::endl;
      // Always true, since for the starting particle (no parents), there
      // definitely is no other daughter (hell, there is not even a parent).
      downOneGenerationLoop(particle_info,
          particle, printed_daughters, indentation_level, true);
    }
  }
  // Convert the filled stream to a string and print it out.
  std::string print_string = particle_info.str();
  streamlog_out(MESSAGE) << print_string << std::endl;
  return;
}

// Similar to printFamilyTree(EVENT::LCCollection* collection), but now there is
// a specified particle to start from. No comments in this function, as they
// would just be a repetition of those from above.
void ref_util::printFamilyTree(MCP* particle) {
  std::ostringstream particle_info;
  std::vector<MCP*> printed_daughters;
  std::string indentation_level = "";
  particle_info << std::endl
    << "|" << std::endl
    << "|" << std::endl;
  downOneGenerationLoop(particle_info,
      particle, printed_daughters, indentation_level, true);
  std::string print_string = particle_info.str();
  streamlog_out(MESSAGE) << print_string << std::endl;
  return;
}

// ----------------------------------------------------------------------------
std::vector<RP*> ref_util::getRpsFromMcParticle(
    MCP* particle, UTIL::LCRelationNavigator* relation_navigator) {
  std::vector<RP*> reconstructed_relatives_of_mc;
  std::vector<MCP*> parent_candidates;
  parent_candidates.push_back(particle);
  while (!parent_candidates.empty()) {
    // Retrieve all reconstructed particles that are related to a MC particle
    // from the list.
    MCP* mc_particle= parent_candidates.back();
    parent_candidates.pop_back();
    for (size_t i_rel = 0;
         i_rel < relation_navigator->getRelatedFromObjects(mc_particle).size();
         ++i_rel) {
      RP* pfo_particle = static_cast<RP*>(
          relation_navigator->getRelatedFromObjects(mc_particle)[i_rel]);
      reconstructed_relatives_of_mc.push_back(pfo_particle);
    }
    // Add the Monte Carlo daughters of the candidate parent particle to the
    // list of candidate particles.
    for (size_t i_daughters = 0;
         i_daughters < mc_particle->getDaughters().size();
         ++i_daughters) {
      MCP* daughter = static_cast<MCP*>(
          mc_particle->getDaughters()[i_daughters]);
      parent_candidates.push_back(daughter);
    }
  }
  return reconstructed_relatives_of_mc;
}

std::vector<MCP*> ref_util::getMcChainFromRp(
    RP* rp, UTIL::LCRelationNavigator* relation_navigator) {
  std::vector<MCP*> mc_chain;
  MCP* mcp = static_cast<MCP*>(relation_navigator->getRelatedToObjects(rp)[0]);
  mc_chain.push_back(mcp);
  while (mcp->getParents().size()) {
    mcp = mcp->getParents()[0];
    mc_chain.push_back(mcp);
  }
  return mc_chain;
}

int ref_util::getPrimaryPdg(
    RP* rp, UTIL::LCRelationNavigator* relation_navigator) {
  std::vector<MCP*> mc_chain = ref_util::getMcChainFromRp(rp,
      relation_navigator);
  int primaryPDG = 0;
  if (mc_chain.size() > 0) {
    primaryPDG = mc_chain.back()->getPDG();
  }
  return primaryPDG;
}

// ----------------------------------------------------------------------------
void ref_util::thetaToRpAndVector(EVENT::LCCollection* rp_collection,
    std::map<double, RP*> &theta_to_rp, std::map<double, Tlv> &theta_to_tlv) {
  for (int i_pfo = 0; i_pfo < rp_collection->getNumberOfElements(); ++i_pfo) {
    RP* rp = static_cast<RP*>(rp_collection->getElementAt(i_pfo));
    Tlv rp_tlv = ref_util::getTlv(rp);
    double theta = rp_tlv.Theta();
    while(theta_to_rp.find(theta) != theta_to_rp.end()) {
      // Make sure that the theta does not overwrite a previous value.
      // Should not happen, but due to limited float precision it is possible.
      // Just to be save...
      theta +=  1e-7;
    }
    theta_to_rp[theta]  = rp;
    theta_to_tlv[theta] = rp_tlv;
  }
  return;
}

std::vector<double> ref_util::energyLineupWrtSeedDistance(
    RP* seed_particle, EVENT::LCCollection* rp_collection) {
  std::map<double, RP*> angle_to_seed_map = ref_util::seedDistance(
      seed_particle, rp_collection);
  // The keys of a map are sorted by default. Go through them to get the
  // reconstructed particles sorted by their angular distance to the seed, then
  // retrieve their energy.
  std::vector<double> sorted_energies;
  for (std::map<double, RP*>::iterator it = angle_to_seed_map.begin();
       it != angle_to_seed_map.end(); ++it) {
    streamlog_out(MESSAGE) << "  " << it->first << " E="
                                   << it->second->getEnergy();
    sorted_energies.push_back(it->second->getEnergy());
  }
  return sorted_energies;
}

std::map<double, RP*> ref_util::seedDistance(
    RP* seed_particle, EVENT::LCCollection* rp_collection) {
  std::map<double, RP*> angle_to_seed_map;
  Tlv seed_tlv = getTlv(seed_particle);
  double angle_to_seed = 0.;
  for (int i_pfo = 0; i_pfo < rp_collection->getNumberOfElements(); ++i_pfo) {
    RP* rp = static_cast<RP*>(rp_collection->getElementAt(i_pfo));
    angle_to_seed = ROOT::Math::VectorUtil::CosTheta(getTlv(rp), seed_tlv);
    angle_to_seed_map[angle_to_seed] = rp;
  }
  return angle_to_seed_map;
}

// ----------------------------------------------------------------------------
bool ref_util::pdgIsInMcCol(int pdg, EVENT::LCCollection* mcCol) {
  for (int e = 0; e < mcCol->getNumberOfElements(); ++e) {
    MCP* mcp = static_cast<MCP*>(mcCol->getElementAt(e));
    if (mcp->getPDG() == pdg) {
      return true;
    }
  }
  // Went through all particles and none was of type pdg.
  return false;
}