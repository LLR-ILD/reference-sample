/**
 *  With my C++ programming knowledge so far (11/19) the best way to make these
 *  universally (== in this project) useful functions available to all
 *  classes/processors is by putting them as free functions into a namespace.
 *
 *  As this is a collection of methods, there is no global explanation. Instead,
 *  each method should be explained individually, or be clear from its name.
 *  It is advisable to check from time to time wether the methods described here
 *  are still in use in the project.
 *
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
 */
#ifndef _REF_UTIL_H_
#define _REF_UTIL_H_
// -- C++ STL headers.
#include <map>

// -- ROOT headers.
#include "TMath.h"
// The ROOT::Math vector functions need the GenVector component of root.
// Make sure you have the following line in the CMakeLists.txt file:
//   `FIND_PACKAGE(ROOT REQUIRED COMPONENTS GenVector)`.
#include "Math/VectorUtil.h"
#include "Math/Vector4D.h"
#include "TString.h"

// -- LCIO headers.
#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include "UTIL/LCRelationNavigator.h"

// -- Marlin headers.

// -- Header for this processor and other project-specific headers.

namespace ref_util {

enum PfoCodes {
  kElectron = 0,
  kMuon,
  kPhoton,
  kPiCharged,
  kKShort,
  kNeutron,
  kLambda,
  kPfoTypes,
};
// The following three functions use the PfoCodes enum.
int particlePdgToName(int pdg_val);
TString getPfoLabel(int code);
TString getPfoLabelFromPdg(int pdg_val);

// Returns the (ROOT T)Lorentz 4-vector of a RP.
ROOT::Math::XYZTVector getTlv(EVENT::ReconstructedParticle* rp);
// Returns the sum of the (ROOT T)Lorentz 4-vectors of the given RPs.
ROOT::Math::XYZTVector getTlv(std::vector<EVENT::ReconstructedParticle*> rps);
// Returns the (ROOT T)Lorentz 4-vector of a MC particle.
ROOT::Math::XYZTVector getTlv(EVENT::MCParticle* mcp);
// Returns the sum of the (ROOT T)Lorentz 4-vectors of the given MC particles.
ROOT::Math::XYZTVector getTlv(std::vector<EVENT::MCParticle*> mc_particles);

// Print a family tree for the Monte Carlo particle collection into the shell.
void printFamilyTree(EVENT::LCCollection* collection);
// Print the family tree of the descendents of a MC particle into the shell.
void printFamilyTree(EVENT::MCParticle* particle);

// Given a MC particle, return those reconstructed particles that are linked to
// the MC particle or one of its daughters.
std::vector<EVENT::ReconstructedParticle*> getRpsFromMcParticle(
    EVENT::MCParticle* particle, UTIL::LCRelationNavigator* relation_navigator);
// Return the MC Particle inheritance chain  of the ReconstructedParticle.
// In case of multiple candidates (multiple MC particles are related to the RP
// or the RP has multiple parents), the first candidate is chosen.
std::vector<EVENT::MCParticle*> getMcChainFromRp(
    EVENT::ReconstructedParticle* rp,
    UTIL::LCRelationNavigator* relation_navigator);
// Return the PDG code of the particle on top of the MC inheritance chain
// produced by getMcChainFromRp.
int getPrimaryPdg(EVENT::ReconstructedParticle* rp,
    UTIL::LCRelationNavigator* relation_navigator);

// Should be passed by reference two maps,which are then populated:
// A map from the theta values of the particles in an events reconstructed
// particle collection to a pointer to the particle (theta_rp_map) as well as a
// map from the theta values to the 4-vector of these particles (theta_tlv_map).
void thetaToRpAndVector(EVENT::LCCollection* rp_collection,
    std::map<double, EVENT::ReconstructedParticle*> &theta_to_rp,
    std::map<double, ROOT::Math::XYZTVector> &theta_to_tlv);
// Have the angular distance of a reconstructed particle to the seed_particle
// mapped to the reconstructed particle.
std::map<double, EVENT::ReconstructedParticle*> seedDistance(
    EVENT::ReconstructedParticle* seed_particle,
    EVENT::LCCollection* rp_collection);
// Return a vector of energies of all particles in the collection, sorted by
// their angle with the seed particle.
std::vector<double> energyLineupWrtSeedDistance(
    EVENT::ReconstructedParticle* seed_particle,
    EVENT::LCCollection* rp_collection);

// True if any of the members of the (Monte Carlo) collection has this pdg code.
bool pdgIsInMcCol(int pdg, EVENT::LCCollection* mc_col);
bool rpEnergySort(EVENT::ReconstructedParticle* rp1,
                  EVENT::ReconstructedParticle* rp2){
  return fabs(rp1->getEnergy()) > fabs(rp2->getEnergy());
}
}  // namespace ref_util

#endif