#include "LikelihoodFitterCheck.h"
#include "TVector3.h"
#include <TFile.h>
#include <TH1D.h>
LikelihoodFitterCheck::LikelihoodFitterCheck():Tool(){}


bool LikelihoodFitterCheck::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();
  std::string output_filename;
  mode = "Direction";
  m_variables.Get("verbosity", verbosity);
  m_variables.Get("OutputFile", output_filename);
  m_variables.Get("ifPlot2DFOM", ifPlot2DFOM);
  m_variables.Get("ShowEvent", fShowEvent);
  m_variables.Get("EventRange", fEventRange);
  m_variables.Get("CleanEventsOnly", ifCleanEventsOnly);
  m_variables.Get("UsePDFFile", fUsePDFFile);
  m_variables.Get("PDFFile", pdffile);
  m_variables.Get("2DMode", mode);
  m_variables.Get("DrawTrueDir", DrawTrueDir);
  fOutput_tfile = new TFile(output_filename.c_str(), "recreate");
  
  // Histograms
  if (mode == "Position"){
      Likelihood2D = new TH2D("Likelihood2D","Figure of merit 2D", 200, 0, 200, 100, -100, 100);
  Likelihood2D_pdf = new TH2D("Likelihood2D_pdf", "pdf-based figure of merit 2D", 200, 0, 200, 100, -100, 100);

  }
  else if (mode == "Direction"){
  Likelihood2D = new TH2D("Likelihood2D","Figure of merit 2D", 200, 0, 2*TMath::Pi(), 100, 0, TMath::Pi());
  Likelihood2D_pdf = new TH2D("Likelihood2D_pdf", "pdf-based figure of merit 2D", 200, 0, 2*TMath::Pi(), 100, 0, TMath::Pi());
  }
  gr_parallel = new TGraph();
  gr_parallel->SetTitle("Figure of merit parallel to the track direction");
	gr_transverse = new TGraph();
  gr_transverse->SetTitle("Figure of merit transverse to the track direction");
  pdf_parallel = new TGraph();
  pdf_parallel->SetTitle("PDF-based Figure of merit parallel to track direction");
  pdf_transverse = new TGraph();
  pdf_transverse->SetTitle("PDF-based figure of merit transverse to track direction");
  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  return true;
}


bool LikelihoodFitterCheck::Execute(){
  TH1D* Chi2Values = new TH1D("Chi2Values", "Chi2Values", 100, -100, 100);
	Log("===========================================================================================",v_debug,verbosity);
	
	Log("LikelihoodFitterCheck Tool: Executing",v_debug,verbosity);	
	// Get a pointer to the ANNIEEvent Store
  auto* annie_event = m_data->Stores["RecoEvent"];
  if (!annie_event) {
    Log("Error: The PhaseITreeMaker tool could not find the ANNIEEvent Store",
      0, verbosity);
    return false;
  }

  if (ifPlot2DFOM && (mode != "Position" && mode != "Direction")) {
      Log("Error: invalid 2d-plot mode setting.  Set configvariable '2DMode' to either Position or Direction", v_error);
      return false;
  }

  // MC entry number
  m_data->Stores.at("ANNIEEvent")->Get("MCEventNum",fMCEventNum);  
  
  // MC trigger number
  m_data->Stores.at("ANNIEEvent")->Get("MCTriggernum",fMCTriggerNum); 
  
  // ANNIE Event number
  m_data->Stores.at("ANNIEEvent")->Get("EventNumber",fEventNumber);
  
  // check if event passes the cut
  bool EventCutstatus = false;
  if (ifCleanEventsOnly) {
      auto get_evtstatus = m_data->Stores.at("RecoEvent")->Get("EventCutStatus", EventCutstatus);
      if (!get_evtstatus) {
          Log("Error: The LikelihoodFitterCheck tool could not find the Event selection status", v_error, verbosity);
          return false;
      }
      if (!EventCutstatus) {
          Log("Message: This event doesn't pass the event selection. ", v_message, verbosity);
          return true;
      }
  }

  // Only check this event, or this range of events
  if ((fShowEvent > 0 && (int)fEventNumber <= fShowEvent) || fEventsShown > fEventRange) return true;
  fEventsShown++;
  std::cout << "LFCheck event " << fEventsShown << " out of " << fEventRange << endl;
  	
  logmessage = "Likelihood check for MC Entry Number " + to_string(fMCEventNum) 
               + " , MC Trigger Number" + to_string(fMCTriggerNum) 
               + " and ANNIIE Event Number " + to_string(fEventNumber);
	Log(logmessage,v_message,verbosity);
  
  // Read True Vertex   
  RecoVertex* truevtx = 0;
  auto get_vtx = m_data->Stores.at("RecoEvent")->Get("TrueVertex",fTrueVertex);  ///> Get digits from "RecoEvent" 
  if(!get_vtx){ 
  	Log("LikelihoodFitterCheck  Tool: Error retrieving TrueVertex! ",v_error,verbosity); 
  	return false;
  }
	
	// Retrive digits from RecoEvent
	auto get_digit = m_data->Stores.at("RecoEvent")->Get("RecoDigit",fDigitList);  ///> Get digits from "RecoEvent" 
  if(!get_digit){
  	Log("LikelihoodFitterCheck  Tool: Error retrieving RecoDigits,no digit from the RecoEvent!",v_error,verbosity); 
  	return false;
  }

  if (fUsePDFFile) {
      bool pdftest = this->GetPDF(pdf);
      if (!pdftest) {
          Log("LikelihoodFitterCheck  Tool: Error retrieving pdffile; running without!", v_error, verbosity);
          fUsePDFFile = 0;
          return false;
      }
  }
	
	double recoVtxX, recoVtxY, recoVtxZ, recoVtxT, recoDirX, recoDirY, recoDirZ;
  double trueVtxX, trueVtxY, trueVtxZ, trueVtxT, trueDirX, trueDirY, trueDirZ;
  double seedX, seedY, seedZ, seedT, seedDirX, seedDirY, seedDirZ;
  double ConeAngle = Parameters::CherenkovAngle();

  // Get true Vertex information
  Position vtxPos = fTrueVertex->GetPosition();
	Direction vtxDir = fTrueVertex->GetDirection();
	trueVtxX = vtxPos.X();
  trueVtxY = vtxPos.Y();
  trueVtxZ = vtxPos.Z();
  trueVtxT = fTrueVertex->GetTime();
  trueDirX = vtxDir.X();
  trueDirY = vtxDir.Y();
  trueDirZ = vtxDir.Z();
      seedDirX = trueDirX;
    seedDirY = trueDirY;
    seedDirZ = trueDirZ;
    double Seedtheta = 0;
    double Seedphi = 0;
    double Truethetadeg = 0;
    double Truephideg = 0;
    double Truethetarad = 0;
    double Truephirad = 0;

    
  std::string plotname;
  
  if(verbosity>0) cout<<"True vertex  = ("<<trueVtxX<<", "<<trueVtxY<<", "<<trueVtxZ<<", "<<trueVtxT<<", "<<trueDirX<<", "<<trueDirY<<", "<<trueDirZ<<")"<<endl;
  
  FoMCalculator * myFoMCalculator = new FoMCalculator();
  VertexGeometry* myvtxgeo = VertexGeometry::Instance();
  myvtxgeo->LoadDigits(fDigitList);
  myFoMCalculator->LoadVertexGeometry(myvtxgeo); //Load vertex geometry
  //parallel direction
  double dl = 1.0; // step size  = 1 cm along the track
  double dx = dl * trueDirX;
  double dy = dl * trueDirY;
  double dz = dl * trueDirZ;
  int nbins = 400;
  double dlpara[400], dlfom[400];
  double minphi, maxphi;
  for(int j=0;j<400;j++) {
    seedX = trueVtxX - 350*dx + j*dx;
    seedY = trueVtxY - 350*dy + j*dy;
    seedZ = trueVtxZ - 350*dz + j*dz;
    seedT = trueVtxT;
    seedDirX = trueDirX;
    seedDirY = trueDirY;
    seedDirZ = trueDirZ;
    myvtxgeo->CalcExtendedResiduals(seedX, seedY, seedZ, 0.0, seedDirX, seedDirY, seedDirZ);
    int nhits = myvtxgeo->GetNDigits();
    double meantime = myFoMCalculator->FindSimpleTimeProperties(ConeAngle);
    Double_t fom = -999.999*100;
    double timefom = -999.999*100;
    double timefomlikelihood = -999.999*100;
    double conefom = -999.999*100;
    double conefomlnl = -999.999 * 100;
    Double_t fompdf = -999.999 * 100;
    myFoMCalculator->TimePropertiesLnL(meantime,timefom);
    myFoMCalculator->ConePropertiesFoM(ConeAngle,conefom);
     fom = timefom*0.5+conefom*0.5;
    // fom = conefom;
    cout<<"timeFOM, coneFOM, fom = "<<timefom<<", "<<conefom<<", "<<fom<<endl;
    // fom = timefom;
    dlpara[j] = - 350*dl + j*dl;
    dlfom[j] = fom;
    gr_parallel->SetPoint(j, dlpara[j], dlfom[j]);

    if (fUsePDFFile) {
      //myFoMCalculator->ConePropertiesFOM(seedX, seedY, seedZ, seedDirX, seedDirY, seedDirZ, ConeAngle, conefomlnl, pdf, maxphi, minphi);
      myFoMCalculator->ConePropertiesNotMine(seedX, seedY, seedZ, seedDirX, seedDirY, seedDirZ, ConeAngle, conefomlnl);
        cout << "conefomlnl: " << conefomlnl << endl;
        myFoMCalculator->TimePropertiesLnL(meantime,timefomlikelihood);
	//fompdf = 0.5*timefomlikelihood + 0.5*conefomlnl;
	fompdf = conefomlnl;
	//fompdf = timefomlikelihood;
        pdf_parallel->SetPoint(j, dlpara[j], fompdf);
    }
  } 
  
  //transverse direction
  double dltrans[200];
  //first find the projection vector of x axis on the plane perpendicular to the track direction
  TVector3 n(trueDirX, trueDirY, trueDirZ);
  TVector3 v = n.Orthogonal();
  dl = 1.0; // step size  = 1 cm along the track
  dx = dl * v.X();
  dy = dl * v.Y();
  dz = dl * v.Z();
  for(int j=0;j<100;j++) {
    seedX = trueVtxX - 50*dx + j*dx;
    seedY = trueVtxY - 50*dy + j*dy;
    seedZ = trueVtxZ - 50*dz + j*dz;
    seedT = trueVtxT;
    seedDirX = trueDirX;
    seedDirY = trueDirY;
    seedDirZ = trueDirZ;
    double Seedtheta = std::acos(seedDirZ / std::sqrt(seedDirX*seedDirX + seedDirY*seedDirY + seedDirZ*seedDirZ)*180*(1/TMath::Pi()));
                               double Seedphi = std::atan2(seedDirY,seedDirX)*180*(1/TMath::Pi());

    myvtxgeo->CalcExtendedResiduals(seedX, seedY, seedZ, 0.0, seedDirX, seedDirY, seedDirZ);
    int nhits = myvtxgeo->GetNDigits();
    double meantime = myFoMCalculator->FindSimpleTimeProperties(ConeAngle);
    Double_t fom = -999.999*100;
    Double_t fompdf = -999.999 * 100;
    double timefom = -999.999*100;
    double timefomlikelihood = -999.999*100;
    double conefom = -999.999*100;
    double conefomlnl = -999.999 * 100;
    double coneAngle = 42.0;
    double phimax,phimin;
    myFoMCalculator->TimePropertiesLnL(meantime,timefom);
    myFoMCalculator->ConePropertiesFoM(ConeAngle,conefom);
    fom = timefom*0.5+conefom*0.5;
    //fom = timefom;
    //fom = conefom;
    cout<<"timeFOM, coneFOM, fom = "<<timefom<<", "<<conefom<<", "<<fom<<endl;
    dltrans[j] = - 50*dl + j*dl;
    dlfom[j] = fom;
    gr_transverse->SetPoint(j, dltrans[j], dlfom[j]);
    if (fUsePDFFile) {
        cout << "pdf fom coming\n";
        //myFoMCalculator->ConePropertiesFOM(seedX, seedY, seedZ, seedDirX, seedDirY, seedDirZ, ConeAngle, conefomlnl, pdf, maxphi, minphi);                                                              
      myFoMCalculator->ConePropertiesNotMine(seedX, seedY, seedZ, seedDirX, seedDirY, seedDirZ, ConeAngle, conefomlnl);
        myFoMCalculator->TimePropertiesLnL(meantime,timefomlikelihood);
	fompdf = 0.5*timefomlikelihood + 0.5*conefomlnl;
	//fompdf = conefomlnl;
	//fompdf = timefomlikelihood;
        pdf_transverse->SetPoint(j, dltrans[j], conefomlnl);
    }
  }
  
    if(ifPlot2DFOM) {
      //2D scan around the true vertex position
        cout << "2DPlot starting now" << endl;
      double dl_para = 1.0, dl_trans = 1.0;
      double dx_para = dl_para * trueDirX;
      double dy_para = dl_para * trueDirY;
      double dz_para = dl_para * trueDirZ;
      double dx_trans = dl_trans * v.X();
      double dy_trans = dl_trans * v.Y();
      double dz_trans = dl_trans * v.Z();
      double phimax, phimin;
      for (int k = 0; k < 100; k++) {
          for (int m = 0; m < 200; m++) {
              seedX = trueVtxX - 50 * dx_trans + k * dx_trans - 50 * dx_para + m * dx_para;
              seedY = trueVtxY - 50 * dy_trans + k * dy_trans - 50 * dy_para + m * dy_para;
              seedZ = trueVtxZ - 50 * dz_trans + k * dz_trans - 50 * dz_para + m * dz_para;
              seedT = trueVtxT;
              seedDirX = cos(m * TMath::Pi() / 100) * sin(k * TMath::Pi() / 100);
              seedDirY = sin(m * TMath::Pi() / 100) * sin(k * TMath::Pi() / 100);
              seedDirZ = cos(k * TMath::Pi() / 100);
cout<<"Check1"<<endl;
              if (mode == "Position") myvtxgeo->CalcExtendedResiduals(seedX, seedY, seedZ, seedT, trueDirX, trueDirY, trueDirZ);
              else if (mode == "Direction") myvtxgeo->CalcExtendedResiduals(trueVtxX, trueVtxY, trueVtxZ, seedT, seedDirX, seedDirY, seedDirZ);
              cout<<"Check2"<<endl;
int nhits = myvtxgeo->GetNDigits();
              double meantime = myFoMCalculator->FindSimpleTimeProperties(ConeAngle);
              Double_t fom = -999.999 * 100;
              Double_t fompdf = -999.999 * 100;
              double timefom = -999.999 * 100;
              double timefomlikelihood = -999.999 * 100;
              double conefom = -999.999 * 100;
              double conefomlnl = -999.999 * 100;
              double coneAngle = 42.0;
              myFoMCalculator->TimePropertiesLnL(meantime, timefom);
              myFoMCalculator->ConePropertiesFoM(coneAngle, conefom);
              fom = timefom * 0.5 + conefom * 0.5;
              //fom = timefom;
	      // fom = conefom;
              cout << "k,m, timeFOM, coneFOM, fom = " << k << ", " << m << ", " << timefom << ", " << conefom << ", " << fom << endl;
              Likelihood2D->SetBinContent(m, k, fom);
              if (fUsePDFFile) {
cout<<"Check3"<<endl;
 //myFoMCalculator->ConePropertiesFOM(seedX, seedY, seedZ, seedDirX, seedDirY, seedDirZ, ConeAngle, conefomlnl, pdf, maxphi, minphi);                                                                 
                  if (mode == "Position")myFoMCalculator->ConePropertiesNotMine(seedX, seedY, seedZ, trueDirX, trueDirY, trueDirZ, coneAngle, conefomlnl);
                  else if (mode == "Direction")myFoMCalculator->ConePropertiesNotMine(trueVtxX, trueVtxY, trueVtxZ, seedDirX, seedDirY, seedDirZ, coneAngle, conefomlnl);

                  myFoMCalculator->TimePropertiesLnL(meantime, timefomlikelihood);
		   fompdf = 0.5 * timefomlikelihood + 0.5 * conefomlnl;
		  //fompdf = conefomlnl;
		  // fompdf = timefomlikelihood;
                  cout << "coneFOMlnl: " << conefomlnl << endl;
                  if (k == 50 && m == 50) {
                      std::cout << "!!!OUTPUT!!! at true:\n";
                  }
                  std::cout << "conefomlnl, timefom, fompdf: " << conefomlnl << ", " << timefom << ", " << fompdf << endl;
                  std::cout << "phimax, phimin: " << phimax << ", " << phimin << endl;
		  
                  Likelihood2D_pdf->SetBinContent(m, k, fompdf);
		  Chi2Values->Fill(fompdf);
              }
          }
      }
    }

cout<<"CheckW1"<<endl;
    if (seedDirX != 0){
      Seedtheta = acos(seedDirZ /sqrt(seedDirX*seedDirX + seedDirY*seedDirY + seedDirZ*seedDirZ))*180*(1/TMath::Pi());
                               Seedphi =atan2(seedDirY,seedDirX)*180*(1/TMath::Pi());
      } else{
  Seedtheta = acos(seedDirZ /sqrt(seedDirX*seedDirX + seedDirY*seedDirY + seedDirZ*seedDirZ))*180*(1/TMath::Pi());
                                Seedphi = 0;

    }
    if (trueDirX != 0){
      Truephideg = acos(trueDirZ /sqrt(trueDirX*trueDirX + trueDirY*trueDirY + trueDirZ*trueDirZ))*180*(1/TMath::Pi());
      Truethetadeg = atan2(trueDirY,trueDirX)*180*(1/TMath::Pi());
      Truephirad = acos(trueDirZ /sqrt(trueDirX*trueDirX + trueDirY*trueDirY + trueDirZ*trueDirZ));
      Truethetarad = atan2(trueDirY,trueDirX);
   
    }
    fOutput_tfile->cd();
    plotname = "FoM_parallel" + std::to_string(fEventNumber);
    gr_parallel->Write(plotname.c_str());
    plotname = "FoM_transvers" + std::to_string(fEventNumber);
    gr_transverse->Write(plotname.c_str());
    if (DrawTrueDir == 1){
    TCanvas *c1  = new TCanvas("c1", "Histogram", 800,600);
    c1->Range(-2*TMath::Pi(),-2*TMath::Pi(),2*TMath::Pi(),2*TMath::Pi());
    c1->DrawFrame(-2*TMath::Pi(),-2*TMath::Pi(),2*TMath::Pi(),2*TMath::Pi());
    Likelihood2D_pdf->Draw("COLZ");
    TMarker *m1 = new TMarker(Truethetarad, Truephirad, 29);
    m1 -> Draw();
    TMarker *m2 =  new TMarker(Truethetarad+2*TMath::Pi(), Truephirad, 32);
      m2 -> Draw();
      c1->Update();
      plotname = "True Direction" + std::to_string(fEventNumber);
      c1->Write(plotname.c_str());
      plotname = "True Direction" + std::to_string(fEventNumber)+".png";
      c1->SaveAs(plotname.c_str());
      delete m1;
      delete m2;
      delete c1;
    }
   plotname = "FoM_2D" + to_string(fEventNumber);
    if(ifPlot2DFOM) Likelihood2D->Write(plotname.c_str());
    	plotname = "Chi2Values" + std::to_string(fEventNumber);
	Chi2Values->Write(plotname.c_str());
cout<<"CheckW2"<<endl;


    if (fUsePDFFile) {
        plotname = "pdfFoM_parallel" + std::to_string(fEventNumber);
        pdf_parallel->Write(plotname.c_str());
        plotname = "pdfFoM_transverse" + std::to_string(fEventNumber);
        pdf_transverse->Write(plotname.c_str());
        plotname = "pdfFoM_2D" + std::to_string(fEventNumber);
        if (ifPlot2DFOM) Likelihood2D_pdf->Write(plotname.c_str());
	plotname = "Chi2Values" + std::to_string(fEventNumber);
	Chi2Values->Write(plotname.c_str());

cout<<"CheckW3"<<endl;
    }
    Likelihood2D->Reset();

    if (fUsePDFFile) {
      Likelihood2D_pdf->Reset();
    }
    
  delete myFoMCalculator;
  return true;
}


bool LikelihoodFitterCheck::Finalise() {
    

    fOutput_tfile->Write();
    fOutput_tfile->Close();

    Log("LikelihoodFitterCheck exitting", v_debug, verbosity);
    return true;
}


  bool LikelihoodFitterCheck::GetPDF(TH1D & pdf) {
      TFile f1(pdffile.c_str(), "READ");
      if (!f1.IsOpen()) {
          Log("VtxExtendedVertexFinder: pdffile does not exist", v_error, verbosity);
          return false;
      }
      pdf = *(TH1D*)f1.Get("zenith");
      cout << "pdf entries: " << pdf.GetEntries() << endl;
      return true;
  }
