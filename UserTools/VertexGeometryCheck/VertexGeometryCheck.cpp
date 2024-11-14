#include "VertexGeometryCheck.h"

    VertexGeometryCheck::VertexGeometryCheck() :Tool() {}


    bool VertexGeometryCheck::Initialise(std::string configfile, DataModel& data) {

        /////////////////// Usefull header ///////////////////////
        if (configfile != "")  m_variables.Initialise(configfile); //loading config file
        //m_variables.Print();
        std::string output_filename;
        m_variables.Get("verbosity", verbosity);
        m_variables.Get("OutputFile", output_filename);
        m_variables.Get("ShowEvent", fShowEvent);
        m_variables.Get("Theta", vertheta);
        m_variables.Get("Phi", verphi);
	m_variables.Get("Steps", steps);
	m_variables.Get("phimod", phimod);
	m_variables.Get("thetamod", thetamod);
        fOutput_tfile = new TFile(output_filename.c_str(), "recreate");

        // Histograms
        flappdextendedtres = new TH1D("lappdextendedtres", "lappd Extended time residual", 1000, -10, 30);
        fpmtextendedtres = new TH1D("pmtextendedtres", "pmt Extended time residual", 1000, -10, 30);
        fpointtres = new TH1D("pointres", "Point time residual", 1000, -10, 30);
        fdelta = new TH1D("delta", "delta", 1000, -10, 30);
        fmeanres = new TH1D("meanres", "Mean time residual", 1000, -10, 30);
        fltrack = new TH1D("ltrack", "track path length", 1000, 0, 1000);
        flphoton = new TH1D("lphoton", "photon path length", 1000, 0, 1000);
        fzenith = new TH1D("zenith", "zenith angle", 180, 0, 180);
        fazimuth = new TH1D("azimuth", "azimuth angle", 180, 0, 360);
        fconeangle = new TH1D("coneangle", "cone angle", 90, 0, 90);
        fdigitcharge = new TH1D("digitcharge", "digit charge", 150, 0, 150);
        fdigittime = new TH1D("digittime", "digit time", 200, -10, 10);
        flappdtimesmear = new TH1D("lappdtimesmear", "lappdtimesmear", 100, 0, 0.1);
        fpmttimesmear = new TH1D("pmttimesmear", "pmttimesmear", 100, 0, 1.0);
        fYvsDigitTheta_all = new TH2D("YvsDigitTheta_all", "Y vs DigitTheta", 400, -200, 200, 400, -200, 200);
        fYvsDigitTheta_all->GetXaxis()->SetTitle("DigitTheta [deg]");
        fYvsDigitTheta_all->GetYaxis()->SetTitle("Digit Y [cm]");
        m_data = &data; //assigning transient data pointer
        /////////////////////////////////////////////////////////////////

        return true;
    }


    bool VertexGeometryCheck::Execute() {
        Log("===========================================================================================", v_debug, verbosity);

        Log("VertexGeometryCheck Tool: Executing", v_debug, verbosity);

        // Get a pointer to the ANNIEEvent Store
        auto* reco_event = m_data->Stores["RecoEvent"];
        if (!reco_event) {
            Log("Error: The PhaseITreeMaker tool could not find the ANNIEEvent Store",
                0, verbosity);
            return false;
        }

        // MC entry number
        m_data->Stores.at("ANNIEEvent")->Get("MCEventNum", fMCEventNum);

        // MC trigger number
        m_data->Stores.at("ANNIEEvent")->Get("MCTriggernum", fMCTriggerNum);

        // ANNIE Event number
        m_data->Stores.at("ANNIEEvent")->Get("EventNumber", fEventNumber);

        // Only check this event
        if (fShowEvent > 0 && (int)fEventNumber != fShowEvent) return true;

        // check if event passes the cut
        bool EventCutstatus = false;
        auto get_evtstatus = m_data->Stores.at("RecoEvent")->Get("EventCutStatus", EventCutstatus);
        if (!get_evtstatus) {
            Log("Error: The VertexGeometryCheck tool could not find the Event selection status", v_error, verbosity);
            return false;
        }
        if (!EventCutstatus) {
            Log("Message: This event doesn't pass the event selection. ", v_message, verbosity);
            return true;
        }

        // Read True Vertex   
        RecoVertex* truevtx = 0;
        auto get_vtx = m_data->Stores.at("RecoEvent")->Get("TrueVertex", fTrueVertex);  ///> Get digits from "RecoEvent" 
        if (!get_vtx) {
            Log("VertexGeometryCheck  Tool: Error retrieving TrueVertex! ", v_error, verbosity);
            return true;
        }

        // Retrive digits from RecoEvent
        auto get_digit = m_data->Stores.at("RecoEvent")->Get("RecoDigit", fDigitList);  ///> Get digits from "RecoEvent" 
        if (!get_digit) {
            Log("VertexGeometryCheck  Tool: Error retrieving RecoDigits,no digit from the RecoEvent!", v_error, verbosity);
            return true;
        }
    

	Log("This Works!", v_message, verbosity);
        double recoVtxX, recoVtxY, recoVtxZ, recoVtxT, recoDirX, recoDirY, recoDirZ;
        double trueVtxX, trueVtxY, trueVtxZ, trueVtxT, trueDirX, trueDirY, trueDirZ;
        double digitX, digitY, digitZ, digitT;
        double dx, dy, dz, px, py, pz, ds, cosphi, sinphi, phi, phideg;

        Position vtxPos = fTrueVertex->GetPosition();
        Direction vtxDir = fTrueVertex->GetDirection();
        trueVtxX = vtxPos.X();
        trueVtxY = vtxPos.Y();
        trueVtxZ = vtxPos.Z();
        trueVtxT = fTrueVertex->GetTime();
        trueDirX = vtxDir.X();
        trueDirY = vtxDir.Y();
        trueDirZ = vtxDir.Z();
	double seedT = trueVtxT;
	double seedDirX = trueDirX;
	double seedDirY = trueDirY;
	double seedDirZ = trueDirZ;
	Log("This Works!", v_message, verbosity);

	double PathDirX, PathDirY, PathDirZ, MagR, MagDirR, PathX, PathY, PathZ;
	PathDirX = -trueDirX;
	PathDirY = -trueDirY;
	PathDirZ = -trueDirZ;
	MagR = 10;
	MagDirR = sqrt(((PathDirX * PathDirX) + (PathDirY * PathDirY) + (PathDirZ * PathDirZ)));
	PathX = PathDirX * MagR / MagDirR;
	PathY = PathDirY * MagR / MagDirR;
	PathZ = PathDirZ * MagR / MagDirR;
	Log("This Works!", v_message, verbosity);
	std::ofstream outputFile2;
	outputFile2.open("reco_points.txt");
if (outputFile2.is_open()) {
  outputFile2 <<" Debug1" << std::endl;

  for (int m = 0; m < 10; m++) {
    double NewX, NewY, NewZ;
    NewX = trueVtxX + PathX * m;
    NewY = trueVtxY + PathY * m;
    NewZ = trueVtxZ + PathZ * m;
    outputFile2 << NewX << " " << NewY << " " << NewZ << std::endl;
  }
 }
 else {
   // Handle error if the file cannot be opened
   std::cerr << "Error: Unable to open file for writing." << std::endl;
 }

outputFile2.close();
 std::ofstream outputFile;
 trueVtxX = trueVtxX + (PathX*steps);
 trueVtxY = trueVtxY + (PathY*steps);
 trueVtxZ = trueVtxZ + (PathZ*steps);
        if (vertheta != -999 && verphi != -999) {
            Log("overriding direction", v_debug, verbosity);
            trueDirX = cos(vertheta) * sin(verphi);
            trueDirY = sin(vertheta) * sin(verphi);
            trueDirZ = cos(verphi);
        }

        double ConeAngle = Parameters::CherenkovAngle();

        FoMCalculator* myFoMCalculator = new FoMCalculator();
        VertexGeometry* myvtxgeo = VertexGeometry::Instance();
        myvtxgeo->LoadDigits(fDigitList);
        myFoMCalculator->LoadVertexGeometry(myvtxgeo); //Load vertex geometry
        int nhits = myvtxgeo->GetNDigits();
        myvtxgeo->CalcExtendedResiduals(trueVtxX, trueVtxY, trueVtxZ, trueVtxT, trueDirX, trueDirY, trueDirZ);
        double meantime = myFoMCalculator->FindSimpleTimeProperties(ConeAngle);
        fmeanres->Fill(meantime);
        double fom = -999.999 * 100;
        myFoMCalculator->TimePropertiesLnL(meantime, fom);
        for (int n = 0; n < nhits; n++) {
            digitX = fDigitList->at(n).GetPosition().X();
            digitY = fDigitList->at(n).GetPosition().Y();
            digitZ = fDigitList->at(n).GetPosition().Z();
            digitT = fDigitList->at(n).GetCalTime();
            dx = digitX - trueVtxX;
            dy = digitY - trueVtxY;
            dz = digitZ - trueVtxZ;
            ds = sqrt(dx * dx + dy * dy + dz * dz);
            px = dx / ds;
            py = dy / ds;
            pz = dz / ds;
            cosphi = 1.0;
            sinphi = 1.0;
            phi = 0.0;
            phideg = 0.0;
	    if( trueDirX*trueDirX + trueDirY*trueDirY + trueDirZ*trueDirZ>0.0 ){
	      // zenith angle
	      cosphi = px*trueDirX+py*trueDirY+pz*trueDirZ;
	      phi = acos(cosphi); // radians
	      phideg = phi/(TMath::Pi()/180.0); // radians->degrees
	      phideg = phideg + phimod;
	    }
            // Y vs theta
            double theta = 0.0;
            double thetadeg = 0.0;
            if (digitZ != 0.0) {
                theta = atan(digitX / digitZ);
            }
            if (digitZ <= 0.0) {
                if (digitX > 0.0) theta += TMath::Pi();
                if (digitX < 0.0) theta -= TMath::Pi();
            }
            thetadeg = theta / (TMath::Pi() / 180.0); // radians->degrees
	    thetadeg = thetadeg + thetamod;
            fYvsDigitTheta_all->Fill(thetadeg, digitY);
            fdelta->Fill(myvtxgeo->GetDelta(n));
            fpointtres->Fill(myvtxgeo->GetPointResidual(n));
            if (myvtxgeo->GetDigitType(n) == RecoDigit::lappd_v0) flappdextendedtres->Fill(myvtxgeo->GetExtendedResidual(n));
            if (myvtxgeo->GetDigitType(n) == RecoDigit::PMT8inch) fpmtextendedtres->Fill(myvtxgeo->GetExtendedResidual(n));
            fltrack->Fill(myvtxgeo->GetDistTrack(n)); //cm
            flphoton->Fill(myvtxgeo->GetDistPhoton(n)); //cm
            fzenith->Fill(phideg,myvtxgeo->GetDigitQ(n)); //
            fazimuth->Fill(myvtxgeo->GetAzimuth(n)); //
            fconeangle->Fill(myvtxgeo->GetConeAngle(n)); //
            fdigitcharge->Fill(myvtxgeo->GetDigitQ(n));
            fdigittime->Fill(digitT);
            if (myvtxgeo->GetDigitType(n) == RecoDigit::lappd_v0) flappdtimesmear->Fill(Parameters::TimeResolution(RecoDigit::lappd_v0, myvtxgeo->GetDigitQ(n)));
            if (myvtxgeo->GetDigitType(n) == RecoDigit::PMT8inch) fpmttimesmear->Fill(Parameters::TimeResolution(RecoDigit::PMT8inch, myvtxgeo->GetDigitQ(n)));
        }
	return true;
    }
        bool VertexGeometryCheck::Finalise() {
            fOutput_tfile->cd();
            fOutput_tfile->Write();
            fOutput_tfile->Close();
            Log("VertexGeometryCheck exitting", v_debug, verbosity);
            return true;
        }
