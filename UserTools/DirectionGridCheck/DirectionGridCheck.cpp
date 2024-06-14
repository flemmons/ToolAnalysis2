#include "DirectionGridCheck.h"

DirectionGridCheck::DirectionGridCheck():Tool(){}


bool DirectionGridCheck::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_variables.Get("OutputFile", OutputFile);
  m_variables.Get("ShowEvent", showEvent);
  m_variables.Get("verbosity", verbosity);

  fOutput_tfile = new TFile("DGCheck.root", "recreate");
  AngularPlot = new TGraph();
  AngularPlot->SetTitle("Seed Directions in Angular Space");
  AngularPlot->GetXaxis()->SetTitle("Polar Angle");
  AngularPlot->GetYaxis()->SetTitle("ZenithAngle");
  AngularPlot->SetMarkerStyle(3);
  AngularPlot->SetMarkerColor(1);
  AngularPlot->SetLineColor(0);
  TrueDirPlot = new TGraph();
  TrueDirPlot->SetTitle("TrueDirection");
  TrueDirPlot->SetMarkerColor(3);
  TrueDirPlot->SetMarkerStyle(25);
  TrueDirPlot->SetLineColor(0);
  MRDDirPlot = new TGraph();
  MRDDirPlot->SetTitle("MRDDirection");
  MRDDirPlot->SetMarkerColor(2);
  MRDDirPlot->SetMarkerStyle(26);
  MRDDirPlot->SetLineColor(0);
  SimpDirPlot = new TGraph();
  SimpDirPlot->SetTitle("SimpleDirection");
  SimpDirPlot->SetMarkerColor(4);
  SimpDirPlot->SetMarkerStyle(27);
  SimpDirPlot->SetLineColor(0);
  RecoDirPlot = new TGraph();
  RecoDirPlot->SetMarkerColor(5);
  RecoDirPlot->SetMarkerStyle(28);
  RecoDirPlot->SetLineColor(0);

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  return true;
}


bool DirectionGridCheck::Execute(){

    auto* annie_event = m_data->Stores["RecoEvent"];
    if (!annie_event) {
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

    auto get_flags = m_data->Stores.at("RecoEvent")->Get("EventFlagged", fEventStatusFlagged);

    if ((fEventStatusFlagged) != 0) {
        //  if (!fEventCutStatus){
        Log("PhaseIITreeMaker Tool: Event was flagged with one of the active cuts.", v_debug, verbosity);
        return true;
    }

    // Only check this event
    //if (showEvent > 0 && (int)fEventNumber != showEvent) return true;
    if (eventCount > 9)return true;

    VertexGeometry* myvtxgeo = VertexGeometry::Instance();
    Direction tempdir;


    if (eventCount == 0) {
        for (int l = 0; l < 100; l++) {
            double phi = (6 * TMath::Pi() / 50) * l;
            double theta = (TMath::Pi() / 200) * l;

            while (phi > 2 * TMath::Pi()) {
                std::cout << "correcting theta" << endl;
                phi -= 2 * TMath::Pi();
            }

            tempdir.SetTheta(theta);
            tempdir.SetPhi(phi);
            AngularPlot->SetPoint(l, phi, theta);
            std::cout << "Angular Plot set:" << l << ", " << phi << ", " << theta << endl;
        }
    }
    RecoVertex* truevtx = 0;
    RecoVertex* recovtx = 0;
    auto get_muonMC = m_data->Stores.at("RecoEvent")->Get("TrueVertex", truevtx);
    auto get_muonReco = m_data->Stores.at("RecoEvent")->Get("ExtendedVertex", recovtx);
    double truePhi = truevtx->GetDirection().GetPhi();
    double trueTheta = truevtx->GetDirection().GetTheta();
    double RecoDirX = recovtx->GetDirection().X();
    double RecoDirY = recovtx->GetDirection().Y();
    double RecoDirZ = recovtx->GetDirection().Z();
    double recoPhi = atan(RecoDirY/RecoDirX);
    double recoTheta = asin(RecoDirX * RecoDirX + RecoDirY * RecoDirY);
    if (trueTheta < 0) trueTheta += 2*TMath::Pi();
    if (trueTheta > 2 * TMath::Pi()) trueTheta -= 2 * TMath::Pi();
    TrueDirPlot->SetPoint(eventCount, truePhi, trueTheta);
    std::cout << "True Plot set: " << eventCount << ", " << truePhi << ", " << trueTheta << endl;
    Direction MRDDirection = this->findDirectionMRD();
    MRDDirPlot->SetPoint(eventCount, MRDDirection.GetPhi(), MRDDirection.GetTheta());
    std::cout << "MRD Plot set: " << eventCount << ", " << MRDDirection.GetPhi() << ", " << MRDDirection.GetTheta() << endl;
    RecoDirPlot->SetPoint(eventCount, recoPhi, recoTheta);
    std::cout << "recoPlot set: " << eventCount << ", " << recoPhi << ", " << recoTheta << endl;

    eventCount++;
  return true;
}


bool DirectionGridCheck::Finalise(){
    fOutput_tfile->cd();
    AngularPlot->Write("Angular Grid");
    TrueDirPlot->Write("trueDirection");
    MRDDirPlot->Write("MRDDirection");
    RecoDirPlot->Write("recoDirection");
    fOutput_tfile->Write();
    fOutput_tfile->Close();
  return true;
}

Direction DirectionGridCheck::findDirectionMRD() {
	std::vector<BoostStore>* Tracks;
    int numtracksinev;
	m_data->Stores["MRDTracks"]->Get("MRDTracks", Tracks);
	m_data->Stores["MRDTracks"]->Get("NumMrdTracks", numtracksinev);

	if (numtracksinev > 1) Log("Multiple tracks need work; just using first for now", v_debug, verbosity);
	double gradx, grady, theta, phi;
	Direction startVertex, endVertex, result;
	BoostStore* thisTrack = &(Tracks->at(0));

	thisTrack->Get("VTrackGradient", gradx);
	thisTrack->Get("HTrackGradient", grady);
	phi = atan(grady / gradx);
	theta = asin(pow((gradx * gradx + grady * grady), 0.5));
	/*TRandom3 smear;
	Direction vtxDir = fTrueVertex->GetDirection();
	Direction result;
	result.SetTheta(smear.Gaus(vtxDir.GetTheta(), 0.4));
	result.SetPhi(smear.Gaus(vtxDir.GetPhi(), 0.4));*/
	result.SetTheta(theta);
	result.SetPhi(phi);

	return result;
}

RecoVertex* DirectionGridCheck::FindSimpleDirection(RecoVertex* myVertex) {

    /// get vertex position
    double vtxX = myVertex->GetPosition().X();
    double vtxY = myVertex->GetPosition().Y();
    double vtxZ = myVertex->GetPosition().Z();
    double vtxTime = myVertex->GetTime();

    std::cout << "Simple Direction Input Position: (" << vtxX << "," << vtxY << "," << vtxZ << ")\n";
    // current status
    // ==============
    int status = myVertex->GetStatus();

    /// loop over digits
    /// ================
    double Swx = 0.0;
    double Swy = 0.0;
    double Swz = 0.0;
    double Sw = 0.0;
    double digitq = 0.;
    double dx, dy, dz, ds, px, py, pz, q;

    RecoDigit digit;
    for (int idigit = 0; idigit < fDigitList->size(); idigit++) {
        digit = fDigitList->at(idigit);
        if (digit.GetFilterStatus()) {
            q = digit.GetCalCharge();
            dx = digit.GetPosition().X() - vtxX;
            dy = digit.GetPosition().Y() - vtxY;
            dz = digit.GetPosition().Z() - vtxZ;
            ds = sqrt(dx * dx + dy * dy + dz * dz);
            px = dx / ds;
            py = dx / ds;
            pz = dz / ds;
            Swx += q * px;
            Swy += q * py;
            Swz += q * pz;
            Sw += q;
        }
    }

    /// average direction
    /// =================
    double dirX = 0.0;
    double dirY = 0.0;
    double dirZ = 0.0;

    int itr = 0;
    bool pass = 0;
    double fom = 0.0;

    if (Sw > 0.0) {
        double qx = Swx / Sw;
        double qy = Swy / Sw;
        double qz = Swz / Sw;
        double qs = sqrt(qx * qx + qy * qy + qz * qz);

        dirX = qx / qs;
        dirY = qy / qs;
        pass = 1;
    }

    // set vertex and direction
    // ========================
    RecoVertex* newVertex = new RecoVertex(); // Note: pointer must be deleted by the invoker

    if (pass) {
        newVertex->SetVertex(vtxX, vtxY, vtxZ, vtxTime);
        newVertex->SetDirection(dirX, dirY, dirZ);
        newVertex->SetFOM(fom, itr, pass);
    }

    // set status
    // ==========
    if (!pass) status |= RecoVertex::kFailSimpleDirection;
    newVertex->SetStatus(status);
    cout << "simple direction: " << dirX << ", " << dirY << ", " << dirZ << endl;

        // return vertex
        // =============
        return newVertex;
}