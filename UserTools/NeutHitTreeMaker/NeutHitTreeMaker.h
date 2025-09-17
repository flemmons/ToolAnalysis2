#ifndef NeutHitTreeMaker_H
#define NeutHitTreeMaker_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "TApplication.h"
#include <Math/PxPyPzE4D.h>
#include <Math/LorentzVector.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"
#include "ADCPulse.h"
#include "Waveform.h"
#include "CalibratedADCWaveform.h"
#include "Hit.h"
#include "RecoDigit.h"
#include "ANNIEalgorithms.h"
#include "TimeClass.h"
#include "BeamStatus.h"

/**
 * \class NeutHitTreeMaker
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: F.A.Lemmons $
* $Date: 2024/04/14 10:44:00 $
* Contact: franklin.lemmons@mines.sdsmt.edu
*/
class NeutHitTreeMaker: public Tool {


 public:

  NeutHitTreeMaker(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.
  void ResetVariables();
  void LoadTrueNeutronInfo();
  void LoadTankClusterHitsMC(std::vector<MCHit> cluster_hits, std::vector<unsigned long> cluster_detkeys);
  void LoadAllTankHits();


 private:
	 //Digits
	 int fNHits = 0;
	 std::vector<int> fIsFiltered;
	 std::vector<double> fHitX;
	 std::vector<double> fHitY;
	 std::vector<double> fHitZ;
	 std::vector<double> fHitT;
	 std::vector<double> fHitQ;
	 std::vector<double> fHitPE;
	 std::vector<int> fHitType;
	 std::vector<int> fHitDetID;
	 std::vector<int> fHitChankey;
	 std::vector<int> fHitChankeyMC;
	 uint32_t fEventNumber;

	 //Primary Pdg
	 std::vector<int>* fTruePrimaryPdgs = nullptr;

	 //Neutrons
	 std::vector<double>* fTrueNeutCapVtxX = nullptr;
	 std::vector<double>* fTrueNeutCapVtxY = nullptr;
	 std::vector<double>* fTrueNeutCapVtxZ = nullptr;
	 std::vector<double>* fTrueNeutCapNucleus = nullptr;
	 std::vector<double>* fTrueNeutCapTime = nullptr;
	 std::vector<double>* fTrueNeutCapGammas = nullptr;
	 std::vector<double>* fTrueNeutCapE = nullptr;
	 std::vector<double>* fTrueNeutCapGammaE = nullptr;
	 int fTrueNeutrons;
	 bool onlyNeutronEvents=0;
	 double fTrueNeutrinoEnergy;
	 double fTrueMuonEnergy;

	 double fTrueVtxX;
	 double fTrueVtxY;
	 double fTrueVtxZ;
	 double fTrueVtxTime;
	 double fTrueDirX;
	 double fTrueDirY;
	 double fTrueDirZ;

	 Geometry* geom = nullptr;
	 std::map<int, double> ChannelKeyToSPEMap;
	 std::map<int, unsigned long> pmtid_to_channelkey;
	 std::map<unsigned long, int> channelkey_to_pmtid;

	 /// \brief TTree that will be used to store output
	 TTree* fPhaseIITrigTree = nullptr;

	 /// \brief ROOT TFile that will be used to store the output from this tool
	 TFile* fOutput_tfile = nullptr;

	 // \brief Event Status flag masks
	 int fEventStatusApplied;
	 int fEventStatusFlagged;

	 /// \brief Integer that determines the level of logging to perform
	 int verbosity = 0;
	 int v_error = 0;
	 int v_warning = 1;
	 int v_message = 2;
	 int v_debug = 3;
	 std::string logmessage;
	 int get_ok;
};


#endif
