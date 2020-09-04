/**
 *  This processor can (only) be called after a Z candidate was chosen.
 *  It produces a .root file that then can readily be input to my master thesis
 *  analysis chain.
 *
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
 */
#ifndef _MASTER_THESIS_ROOT_PROCESSOR_H_
#define _MASTER_THESIS_ROOT_PROCESSOR_H_
// -- C++ STL headers.

// -- ROOT headers.
#include "TFile.h"
#include "TTree.h"

// -- LCIO headers.

// -- Marlin headers.
#include "marlin/Processor.h"

// -- Header for this processor and other project-specific headers.

class MasterThesisRootProcessor : public marlin::Processor {
 public:
  marlin::Processor* newProcessor() { return new MasterThesisRootProcessor(); }
  MasterThesisRootProcessor();

  // These two lines avoid frequent compiler warnings when using -Weffc++.
  MasterThesisRootProcessor(const MasterThesisRootProcessor&) = delete;
  MasterThesisRootProcessor& operator=(const MasterThesisRootProcessor&) = delete;

  void init();
  void processEvent(EVENT::LCEvent* event);
  void end();

  void initRoot();
  void endRoot();
  void setIsoLeptonNumbers(EVENT::LCEvent* event,
                           EVENT::LCCollection* higgs_collection);

 private:
  // -- Parameters registered in steering file.
  // Collections
  std::string higgs_only_collection_name_{""};
  std::string z_only_collection_name_{""};
  std::string mc_collection_name{""};

  // -- The root file
  TFile* root_file_{};
  std::string root_file_name_ = {""};
  std::string z_decay_channel_ = {""};
  TTree* tree_ = {};

  bool missing_mc_collection = false;
  struct HiggsTruth {
    bool decays_invisible = false;
    int decay_mode = -1;
  };
  HiggsTruth getHiggsTruth(EVENT::LCEvent* event);

  struct TreeVars {
    TreeVars() {higgs_truth = HiggsTruth();};
    ~TreeVars() {};

    double m_z = -1;
    double m_recoil = -1;
    double m_vis = -1;
    double m_miss = -1;
    double m_h = -1;
    double m_h_recoil = -1;
    double cos_theta_z = -1;
    double cos_theta_miss = -1;
    double cos_theta_iso_lep = -1;
    double e_highest_iso_lep = -1;
    double principle_thrust_z = -1;
    double principle_thrust = -1;
    double major_thrust = -1;
    double minor_thrust = -1;
    double oblateness = -1;
    double sphericity = -1;
    double aplanarity = -1;

    int n_ch_hadrons = -1;
    int n_n_hadrons = -1;
    int n_gamma = -1;
    int n_electrons = -1;
    int n_muons = -1;
    int n_iso_leptons = -1;

    HiggsTruth higgs_truth{};

    void initBranches(TTree* tree) {
      tree->Branch(("mZ"), &m_z, ("mZ/D"));
      tree->Branch(("mRecoil"), &m_recoil, ("mRecoil/D"));
      tree->Branch(("mVis"), &m_vis, ("mVis/D"));
      tree->Branch(("mMiss"), &m_miss, ("mMiss/D"));
      tree->Branch(("mH"), &m_h, ("mH/D"));
      tree->Branch(("mHRecoil"), &m_h_recoil, ("mHRecoil/D"));
      tree->Branch(("cosTZ"), &cos_theta_z, ("cosTZ/D"));
      tree->Branch(("cosTMiss"), &cos_theta_miss, ("cosTMiss/D"));
      tree->Branch(("cosTIsoLep"), &cos_theta_iso_lep, ("cosTIsoLep/D"));
      tree->Branch(("eHighestIsoLep"), &e_highest_iso_lep, ("eHighestIsoLep/D"));
      tree->Branch(("principleThrustZ"), &principle_thrust_z, ("principleThrustZ/D"));
      tree->Branch(("principleThrust"), &principle_thrust, ("principleThrust/D"));
      tree->Branch(("majorThrust"), &major_thrust, ("majorThrust/D"));
      tree->Branch(("minorThrust"), &minor_thrust, ("minorThrust/D"));
      tree->Branch(("sphericity"), &sphericity, ("sphericity/D"));
      tree->Branch(("aplanarity"), &aplanarity, ("aplanarity/D"));

      tree->Branch(("nChargedHadrons"), &n_ch_hadrons, ("nChargedHadrons/I"));
      tree->Branch(("nNeutralHadrons"), &n_n_hadrons, ("nNeutralHadrons/I"));
      tree->Branch(("nGamma"), &n_gamma, ("nGamma/I"));
      tree->Branch(("nElectrons"), &n_electrons, ("nElectrons/I"));
      tree->Branch(("nMuons"), &n_muons, ("nMuons/I"));
      tree->Branch(("nIsoLeptons"), &n_iso_leptons, ("nIsoLeptons/I"));

      tree->Branch(("hInvisible"), &higgs_truth.decays_invisible, ("hInvisible/I"));
      tree->Branch(("hDecay"), &higgs_truth.decay_mode, ("hDecay/I"));
    }

    void resetValues() {
      m_z = 0;
      m_recoil = 0;
      m_vis = 0;
      m_miss = 0;
      m_h = 0;
      m_h_recoil = 0;
      cos_theta_z = 0;
      cos_theta_miss = 0;
      cos_theta_iso_lep = -2;
      e_highest_iso_lep = 0;
      principle_thrust_z = 0;
      principle_thrust = 0;
      major_thrust = 0;
      minor_thrust = 0;
      oblateness = 0;
      sphericity = 0;
      aplanarity = 0;

      n_ch_hadrons = 0;
      n_n_hadrons = 0;
      n_gamma = 0;
      n_electrons = 0;
      n_muons = 0;
      n_iso_leptons = 0;

      higgs_truth = HiggsTruth();
    }
  };
  TreeVars tv{};
};
#endif