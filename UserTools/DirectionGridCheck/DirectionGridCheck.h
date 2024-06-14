#ifndef DirectionGridCheck_H
#define DirectionGridCheck_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "TH1.h"
#include "TH2D.h"
#include "TFile.h"

#include "VertexGeometry.h"
#include "Parameters.h"
#include "TTree.h"


/**
 * \class DirectionGridCheck
 *
 * This is a debug tool to look at the differences between the various directions associated with a muontrack event.
*
* $Author: F. A. Lemmons $
* $Date: 2024/02/22 01:38:00 $
* Contact: franklin.lemmons@mines.sdsmt.edu
*/
class DirectionGridCheck: public Tool {


 public:

  DirectionGridCheck(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.
  RecoVertex* FindSimpleDirection(RecoVertex* myVertex);
  Direction findDirectionMRD();


 private:
	 std::string OutputFile;
	 int showEvent;
	 TFile* fOutput_tfile = nullptr;
	 TTree* fVertexGeometry = nullptr;
	 /// \brief MC entry number
	 uint64_t fMCEventNum;

	 /// \brief trigger number
	 uint16_t fMCTriggerNum;

	 /// \brief ANNIE event number
	 uint32_t fEventNumber;

	 std::vector<RecoDigit>* fDigitList = 0;
	 RecoVertex* fTrueVertex = 0;
	 RecoVertex* fSDVertex = 0;
	 int verbosity = -1;
	 int v_error = 0;
	 int v_warning = 1;
	 int v_message = 2;
	 int v_debug = 3;

	 int eventCount = 0;

	 // \brief Event Status flag masks
	 int fEventStatusApplied;
	 int fEventStatusFlagged;
	 std::vector<BoostStore>* theMRDTracks;

	 std::vector<BoostStore>* theMrdTracks;   // the MRD tracks
	 TGraph* AngularPlot = 0;
	 TGraph* MRDDirPlot = 0;
	 TGraph* TrueDirPlot = 0;
	 TGraph* SimpDirPlot = 0;
	 TGraph* RecoDirPlot = 0;
	 TTree* DirectionTree = 0;

};


#endif
