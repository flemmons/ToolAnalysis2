#ifndef VertexGeometryCheck_H
#define VertexGeometryCheck_H

#include <string>
#include <iostream>

#include "FoMCalculator.h"
#include "VertexGeometry.h"
#include "Parameters.h"
#include "Tool.h"
#include "TTree.h"
#include "TH2D.h"

class VertexGeometryCheck: public Tool {


 public:

  VertexGeometryCheck();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();


 private:
 	/// \brief ROOT TFile that will be used to store the output from this tool
  TFile* fOutput_tfile = nullptr;

  /// \brief TTree that will be used to store output
  TTree* fVertexGeometry = nullptr;
  
 	/// \brief MC entry number
  uint64_t fMCEventNum;
  
  /// \brief trigger number
  uint16_t fMCTriggerNum;
  
  /// \brief ANNIE event number
  uint32_t fEventNumber;
  
  /// \brief recodigit vector
 	std::vector<RecoDigit>* fDigitList = 0;
    std::vector<RecoCluster>* fClusterList = 0;
    std::vector<RecoDigit> tempDigitList;
 		
 	/// \brief true vertex pointer
 	RecoVertex* fTrueVertex = 0;

  /// \brief reco vertex pointer
  RecoVertex* fRecoVertex = 0;
 	
 	/// \brief select a particle event to show
 	int fShowEvent = 0;
  
  /// \brief Histogram
  double recoVtxX, recoVtxY, recoVtxZ, recoVtxT, recoDirX, recoDirY, recoDirZ, recoVtxFoM;
  double trueVtxX, trueVtxY, trueVtxZ, trueVtxT, trueDirX, trueDirY, trueDirZ;
  TH1D *flappdextendedtres; ///< lappd extended time residual
  TH1D *fpmtextendedtres; ///< pmt extended time residual
  TH1D *fpointtres; ///< point time residual
  TH1D *fdelta; ///< extended time residual
  TH1D *fmeanres; ///< mean value of the time residual distribution
  TH1D *fltrack; ///< muon track path length before producing a particular Cherenkov photon
  TH1D *flphoton; ///< Cherenkov photon track path length
  TH1D *fzenith;
  TH1D *fazimuth;
  // Reco-vertex counterparts of the above
  TH1D *flappdextendedtresReco; ///< lappd extended time residual (reco vertex)
  TH1D *fpmtextendedtresReco; ///< pmt extended time residual (reco vertex)
  TH1D *fpointtresReco; ///< point time residual (reco vertex)
  TH1D *fdeltaReco; ///< extended time residual (reco vertex)
  TH1D *fmeanresReco; ///< mean value of the time residual distribution (reco vertex)
  TH1D *fltrackReco; ///< muon track path length (reco vertex)
  TH1D *flphotonReco; ///< Cherenkov photon track path length (reco vertex)
  TH1D *fzenithReco;
  TH1D *fazimuthReco;
  TH1D *fconeangle;
  TH1D *fdigitcharge;
  TH1D *fdigittime;
  TH1D *fpmtdigittime;
  TH1D *flappddigittime;
  TH1D *flappdtimesmear; 
  TH1D *fpmttimesmear;
  TH2D* fYvsDigitTheta_all;
  TH2D* StripHits1;
  double vertheta = -999, verphi = -999;
  int StripTimePlot = -1;
  bool cleanHitsOnly = 0;
  bool cleanEventsOnly = 0;
  bool fRecoCluster;


/// verbosity levels: if 'verbosity' < this level, the message type will be logged.
  int verbosity=-1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	std::string logmessage;
	int get_ok;	
};


#endif
