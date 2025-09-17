#include "NeutHitTreeMaker.h"

NeutHitTreeMaker::NeutHitTreeMaker():Tool(){}


bool NeutHitTreeMaker::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  m_variables.Get("verbose", verbosity);
  std::string output_filename;
  m_variables.Get("OutputFile", output_filename);
  m_variables.Get("OnlyNeutrons", onlyNeutronEvents);
  fOutput_tfile = new TFile(output_filename.c_str(), "recreate");
  
  m_data->CStore.Get("pmt_tubeid_to_channelkey_data", pmtid_to_channelkey);
  m_data->CStore.Get("channelkey_to_pmtid", channelkey_to_pmtid);
  m_data->CStore.Get("ChannelNumToTankPMTSPEChargeMap", ChannelKeyToSPEMap);

  auto get_geometry = m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry", geom);
  if (!get_geometry) {
      Log("PhaseIITreeMaker Tool: Error retrieving Geometry from ANNIEEvent!", v_error, verbosity);
      return false;
  }

  fPhaseIITrigTree = new TTree("phaseIINeutTree", "ANNIE Phase II Neutron Hit Tree");
  fPhaseIITrigTree->Branch("filter", &fIsFiltered);
  fPhaseIITrigTree->Branch("hitX", &fHitX);
  fPhaseIITrigTree->Branch("hitY", &fHitY);
  fPhaseIITrigTree->Branch("hitZ", &fHitZ);
  fPhaseIITrigTree->Branch("hitT", &fHitT);
  fPhaseIITrigTree->Branch("hitQ", &fHitQ);
  fPhaseIITrigTree->Branch("hitPE", &fHitPE);
  fPhaseIITrigTree->Branch("hitType", &fHitType);
  fPhaseIITrigTree->Branch("hitDetID", &fHitDetID);
  //fPhaseIITrigTree->Branch("hitChankey", &fHitChankey);
  //fPhaseIITrigTree->Branch("hitChankeyMC", &fHitChankeyMC);

  fTrueNeutCapVtxX = new std::vector<double>;
  fTrueNeutCapVtxY = new std::vector<double>;
  fTrueNeutCapVtxZ = new std::vector<double>;
  fTrueNeutCapNucleus = new std::vector<double>;
  fTrueNeutCapTime = new std::vector<double>;
  fTrueNeutCapGammas = new std::vector<double>;
  fTrueNeutCapE = new std::vector<double>;
  fTrueNeutCapGammaE = new std::vector<double>;
  fTruePrimaryPdgs = new std::vector<int>;

  fPhaseIITrigTree->Branch("eventNumber", &fEventNumber, "eventNumber/I");
  fPhaseIITrigTree->Branch("truePrimaryPdgs", &fTruePrimaryPdgs);
  fPhaseIITrigTree->Branch("trueNeutCapVtxX", &fTrueNeutCapVtxX);
  fPhaseIITrigTree->Branch("trueNeutCapVtxY", &fTrueNeutCapVtxY);
  fPhaseIITrigTree->Branch("trueNeutCapVtxZ", &fTrueNeutCapVtxZ);
  fPhaseIITrigTree->Branch("trueNeutCapNucleus", &fTrueNeutCapNucleus);
  fPhaseIITrigTree->Branch("trueNeutCapTime", &fTrueNeutCapTime);
  fPhaseIITrigTree->Branch("trueNeutCapGammas", &fTrueNeutCapGammas);
  fPhaseIITrigTree->Branch("trueNeutCapE", &fTrueNeutCapE);
  fPhaseIITrigTree->Branch("trueNeutCapGammaE", &fTrueNeutCapGammaE);
  fPhaseIITrigTree->Branch("trueNeutrinoEnergy", &fTrueNeutrinoEnergy, "trueNeutrinoEnergy/D");
  fPhaseIITrigTree->Branch("trueNeutrons", &fTrueNeutrons, "trueNeutrons/I");
  fPhaseIITrigTree->Branch("trueMuonEnergy", &fTrueMuonEnergy, "trueMuonEnergy/D");
  fPhaseIITrigTree->Branch("trueVtxX", &fTrueVtxX, "trueVtxX/D");
  fPhaseIITrigTree->Branch("trueVtxY", &fTrueVtxY, "trueVtxY/D");
  fPhaseIITrigTree->Branch("trueVtxZ", &fTrueVtxZ, "trueVtxZ/D");
  fPhaseIITrigTree->Branch("trueVtxTime", &fTrueVtxTime, "trueVtxTime/D");
  fPhaseIITrigTree->Branch("trueDirX", &fTrueDirX, "trueDirX/D");
  fPhaseIITrigTree->Branch("trueDirY", &fTrueDirY, "trueDirY/D");
  fPhaseIITrigTree->Branch("trueDirZ", &fTrueDirZ, "trueDirZ/D");


  return true;
}


bool NeutHitTreeMaker::Execute(){
	Log("===========================================================================================", v_debug, verbosity);
	Log("NeutHitTreeMaker Tool: Executing", v_debug, verbosity);
    
    auto get_flagsapp = m_data->Stores.at("RecoEvent")->Get("EventFlagApplied", fEventStatusApplied);
    auto get_flags = m_data->Stores.at("RecoEvent")->Get("EventFlagged", fEventStatusFlagged);
    //auto get_cutstatus = m_data->Stores.at("RecoEvent")->Get("EventCutStatus",fEventCutStatus);
    if (!get_flagsapp || !get_flags) {
        Log("NeutHitTreeMaker tool: No Event status applied or flagged bitmask!!", v_error, verbosity);
        return false;
    }
    // check if event passes the cut
    if (fEventStatusFlagged != 0) {
        //  if (!fEventCutStatus){
        Log("NeutHitTreeMaker Tool: Event was flagged with one of the active cuts.", v_debug, verbosity);
        std::cout << "flagged Status: " << fEventStatusFlagged << endl;
        return true;
    }

    // Reset variables
    this->ResetVariables();

    Log("Event is clean; writing tree", v_debug, verbosity);
    //std::map<double, std::vector<MCHit>>::iterator it_cluster_pair_mc;
    m_data->Stores.at("RecoEvent")->Get("NeutronCount", fTrueNeutrons);
    if (onlyNeutronEvents && !fTrueNeutrons) {
        Log("NeutHitTreeMaker Tool: This event contains no neutrons.  Skipping.", v_message, verbosity);
    }

    auto get_muonMCEnergy = m_data->Stores.at("RecoEvent")->Get("TrueMuonEnergy", fTrueMuonEnergy);
	
    m_data->Stores.at("ANNIEEvent")->Get("EventNumber", fEventNumber);
    this->LoadAllTankHits();
    this->LoadTrueNeutronInfo();

    fPhaseIITrigTree->Fill();

  return true;
}


bool NeutHitTreeMaker::Finalise(){

    fOutput_tfile->cd();
    fPhaseIITrigTree->Write();
    fOutput_tfile->Close();
    if (verbosity > 0) cout << "NeutHitTreeMaker exitting" << endl;

  return true;
}


void NeutHitTreeMaker::ResetVariables() {
    Log("Resetting variables", v_debug, verbosity);
    fEventNumber = -9999;
    fTruePrimaryPdgs->clear();
    fTrueNeutCapVtxX->clear();
    fTrueNeutCapVtxY->clear();
    fTrueNeutCapVtxZ->clear();
    fTrueNeutCapNucleus->clear();
    fTrueNeutCapTime->clear();
    fTrueNeutCapGammas->clear();
    fTrueNeutCapE->clear();
    fTrueNeutCapGammaE->clear();
    fTrueNeutrons = -9999;
    fTrueMuonEnergy = -9999;


    fIsFiltered.clear();
    fHitX.clear();
    fHitY.clear();
    fHitZ.clear();
    fHitT.clear();
    fHitQ.clear();
    fHitPE.clear();
    fHitType.clear();
    fHitDetID.clear();
    fHitChankey.clear();
    fHitChankeyMC.clear();
    Log("Variables reset; starting tree", v_debug, verbosity);
}

void NeutHitTreeMaker::LoadTankClusterHitsMC(std::vector<MCHit> cluster_hits, std::vector<unsigned long> cluster_detkeys) {
    Position detector_center = geom->GetTankCentre();
    double tank_center_x = detector_center.X();
    double tank_center_y = detector_center.Y();
    double tank_center_z = detector_center.Z();

    //fClusterCharge = 0;
    //fClusterPE = 0;
    //fClusterHits = 0;
    for (int i = 0; i < (int)cluster_hits.size(); i++) {
        unsigned long detkey = cluster_detkeys.at(i);
        int channel_key = (int)detkey;
        int tubeid = cluster_hits.at(i).GetTubeId();
        unsigned long utubeid = (unsigned long)tubeid;
        int wcsimid = channelkey_to_pmtid.at(utubeid);
        unsigned long detkey_data = pmtid_to_channelkey[wcsimid];
        int channel_key_data = (int)detkey_data;
        std::map<int, double>::iterator it = ChannelKeyToSPEMap.find(channel_key_data);
        if (it != ChannelKeyToSPEMap.end()) { //Charge to SPE conversion is available
            Detector* this_detector = geom->ChannelToDetector(tubeid);
            Position det_position = this_detector->GetDetectorPosition();
            unsigned long detkey = this_detector->GetDetectorID();
            double hit_PE = cluster_hits.at(i).GetCharge();
            double hit_charge = hit_PE * ChannelKeyToSPEMap.at(channel_key_data);
            fHitX.push_back((det_position.X() - tank_center_x));
            fHitY.push_back((det_position.Y() - tank_center_y));
            fHitZ.push_back((det_position.Z() - tank_center_z));
            fHitQ.push_back(hit_charge);
            fHitPE.push_back(hit_PE);
            fHitT.push_back(cluster_hits.at(i).GetTime());
            //if (cluster_hits.at(i).GetTime()>2000.) std::cout <<"found cluster hit >  2us! Time: "<<cluster_hits.at(i).GetTime()<<std::endl;
            //fHitDetID.push_back(tubeid);
            //fHitChankey.push_back(channel_key_data);
            //fHitDetID.push_back(detkey);
            //fHitChankey.push_back(channel_key_data);
            //fHitChankeyMC.push_back(channel_key);
            //fHitType.push_back(RecoDigit::PMT8inch);
            //fClusterCharge += hit_charge;
            //fClusterPE += hit_PE;
            //fClusterHits += 1;
        }
        else {
            if (verbosity > 4) {
                std::cout << "FOUND A HIT FOR CHANNELKEY " << channel_key_data << "(MC detkey: " << channel_key << ", chankey = " << utubeid << ", wcsimid = " << wcsimid << ") BUT NO CONVERSION " <<
                    "TO PE AVAILABLE.  SKIPPING PE." << std::endl;
            }
        }
    }
    return;
}

void NeutHitTreeMaker::LoadTrueNeutronInfo() {

    std::map<std::string, std::vector<double>> MCNeutCap;
    bool get_neutcap = m_data->Stores.at("ANNIEEvent")->Get("MCNeutCap", MCNeutCap);
    if (!get_neutcap) {
        Log("PhaseIITreeMaker: Did not find MCNeutCap in ANNIEEvent Store!", v_warning, verbosity);
    }
    std::map<std::string, std::vector<std::vector<double>>> MCNeutCapGammas;
    bool get_neutcap_gammas = m_data->Stores.at("ANNIEEvent")->Get("MCNeutCapGammas", MCNeutCapGammas);
    if (!get_neutcap_gammas) {
        Log("PhaseIITreeMaker: Did not find MCNeutCapGammas in ANNIEEvent Store!", v_warning, verbosity);
    }
    m_data->Stores.at("RecoEvent")->Get("NeutronCount", fTrueNeutrons);

    for (std::map<std::string, std::vector<std::vector<double>>>::iterator it = MCNeutCapGammas.begin(); it != MCNeutCapGammas.end(); it++) {
        std::vector<std::vector<double>> mcneutgammas = it->second;
        for (int i_cap = 0; i_cap < (int)mcneutgammas.size(); i_cap++) {
            std::vector<double> capgammas = mcneutgammas.at(i_cap);
            for (int i_gamma = 0; i_gamma < (int)capgammas.size(); i_gamma++) {
                std::cout << "gamma # " << i_gamma << ", energy: " << capgammas.at(i_gamma) << std::endl;
            }
        }
    }

    if (MCNeutCap.count("CaptVtxX") > 0) {
        std::vector<double> n_vtxx = MCNeutCap["CaptVtxX"];
        std::vector<double> n_vtxy = MCNeutCap["CaptVtxY"];
        std::vector<double> n_vtxz = MCNeutCap["CaptVtxZ"];
        std::vector<double> n_parent = MCNeutCap["CaptParent"];
        std::vector<double> n_ngamma = MCNeutCap["CaptNGamma"];
        std::vector<double> n_totale = MCNeutCap["CaptTotalE"];
        std::vector<double> n_time = MCNeutCap["CaptTime"];
        std::vector<double> n_nuc = MCNeutCap["CaptNucleus"];

        for (int i_cap = 0; i_cap < (int)n_vtxx.size(); i_cap++) {
            fTrueNeutCapVtxX->push_back(n_vtxx.at(i_cap));
            fTrueNeutCapVtxY->push_back(n_vtxy.at(i_cap));
            fTrueNeutCapVtxZ->push_back(n_vtxz.at(i_cap));
            fTrueNeutCapNucleus->push_back(n_nuc.at(i_cap));
            fTrueNeutCapTime->push_back(n_time.at(i_cap));
            fTrueNeutCapGammas->push_back(n_ngamma.at(i_cap));
            fTrueNeutCapE->push_back(n_totale.at(i_cap));
        }
    }

    std::cout << "MCNeutCapGammas count CaptGammas: " << MCNeutCapGammas.count("CaptGammas") << std::endl;
    if (MCNeutCapGammas.count("CaptGammas") > 0) {
        std::vector<std::vector<double>> cap_energies = MCNeutCapGammas["CaptGammas"];
        std::cout << "cap_energies size: " << cap_energies.size() << std::endl;
        for (int i_cap = 0; i_cap < (int)cap_energies.size(); i_cap++) {
            for (int i_gamma = 0; i_gamma < cap_energies.at(i_cap).size(); i_gamma++) {
                std::cout << "gamma energy: " << cap_energies.at(i_cap).at(i_gamma) << std::endl;
                fTrueNeutCapGammaE->push_back(cap_energies.at(i_cap).at(i_gamma));
            }
        }
    }
}

void NeutHitTreeMaker::LoadAllTankHits() {
    bool IsData = 0; //Hard coded for simplicity; MUST change before use on data.
    std::map<unsigned long, std::vector<Hit>>* Hits = nullptr;
    std::map<unsigned long, std::vector<MCHit>>* MCHits = nullptr;
    bool got_hits = false;
    if (IsData) got_hits = m_data->Stores["ANNIEEvent"]->Get("Hits", Hits);
    else got_hits = m_data->Stores["ANNIEEvent"]->Get("MCHits", MCHits);
    if (!got_hits) {
        std::cout << "No Hits store in ANNIEEvent. Continuing to build tree " << std::endl;
        return;
    }
    Position detector_center = geom->GetTankCentre();
    double tank_center_x = detector_center.X();
    double tank_center_y = detector_center.Y();
    double tank_center_z = detector_center.Z();
    fNHits = 0;

    std::map<unsigned long,std::vector<Hit>>::iterator it_tank_data;
    std::map<unsigned long,std::vector<MCHit>>::iterator it_tank_mc;
    if (IsData) it_tank_data = (*Hits).begin();
    else it_tank_mc = (*MCHits).begin();
    bool loop_tank = true;
    //int hits_size = fDigitList->size();
    int hits_size = (IsData)? Hits->size() : MCHits->size();
    //std::cout << "fDigitList: Viable: " << hits_size << "hits\n";
    if (hits_size == 0) {
        std::cout << "Hey! no hits..." << endl;
        loop_tank = false;
    }

    Log("Starting Loop", v_debug, verbosity);
    while (loop_tank) {
        //for (std::pair<unsigned long, std::vector<MCHit>>&& apair : *MCHits) {
            //Log("into second layer; hits list defined", v_debug, verbosity);
            unsigned long channel_key;
            if (IsData) channel_key = it_tank_data->first;
            else channel_key = it_tank_mc->first;
            Detector* this_detector = geom->ChannelToDetector(channel_key);
            Position det_position = this_detector->GetDetectorPosition();
            unsigned long detkey = this_detector->GetDetectorID();
            unsigned long channel_key_data = channel_key;
            Log("Found detector; good so far", v_debug, verbosity);
            if (!IsData) {
                int wcsimid = channelkey_to_pmtid.at(channel_key);
                channel_key_data = pmtid_to_channelkey[wcsimid];
            }
            Log("setting up channels", v_debug, verbosity);
            std::map<int, double>::iterator it = ChannelKeyToSPEMap.find(channel_key);
            std::map<int, double>::iterator it_mc = ChannelKeyToSPEMap.find(channel_key_data);
            bool SPE_available = true;
            if (IsData) SPE_available = (it != ChannelKeyToSPEMap.end());
            else SPE_available = (it_mc != ChannelKeyToSPEMap.end());
            //if (SPE_available) { //Charge to SPE conversion is available
                if (IsData) {
                    std::vector<Hit> ThisPMTHits = it_tank_data->second;
                    fNHits += ThisPMTHits.size();
                    Log("Loading hits", v_debug, verbosity);
                    for (Hit& ahit : ThisPMTHits) {
                        double hit_charge = ahit.GetCharge();
                        //double hit_PE = hit_charge / ChannelKeyToSPEMap.at(channel_key);
                        fHitX.push_back((det_position.X() - tank_center_x));
                        fHitY.push_back((det_position.Y() - tank_center_y));
                        fHitZ.push_back((det_position.Z() - tank_center_z));
                        fHitT.push_back(ahit.GetTime());
                        fHitQ.push_back(hit_charge);
                        //fHitPE.push_back(hit_PE);
                        fHitDetID.push_back(detkey);
                        fHitChankey.push_back(channel_key);
                        fHitChankeyMC.push_back(channel_key);
                        fHitType.push_back(RecoDigit::PMT8inch); // 0 For PMTs
                    }
                }
                else {
                    std::vector<MCHit> ThisPMTHits = it_tank_mc->second;
                    fNHits += ThisPMTHits.size();
                    //fNHits = fDigitList->size();
                    std::cout << "Nhits: " << fNHits << endl;
                    for (/*int i = 0; i < fDigitList->size(); i++ */ MCHit & ahit : ThisPMTHits) {
                        //double hit_PE = ahit.GetCharge();
                        double hit_charge = ahit.GetCharge();
                        //double hit_charge = fDigitList->at(i).GetCalCharge();
                        fHitX.push_back((det_position.X() - tank_center_x));
                        //fHitX.push_back(fDigitList->at(i).GetPosition().X());
                        fHitY.push_back((det_position.Y() - tank_center_y));
                        //fHitY.push_back(fDigitList->at(i).GetPosition().Y());
                        fHitZ.push_back((det_position.Z() - tank_center_z));
                        //fHitZ.push_back(fDigitList->at(i).GetPosition().Z());
                        fHitT.push_back(ahit.GetTime());
                        std::cout << "Hit Time: " << ahit.GetTime() << endl;
                        //fHitT.push_back(fDigitList->at(i).GetCalTime());
                        fHitQ.push_back(hit_charge);
                        //fHitPE.push_back(hit_PE);
                        //fHitDetID.push_back(detkey);
                        //fHitChankey.push_back(channel_key_data);
                        //fHitChankeyMC.push_back(channel_key);
                        fHitType.push_back(RecoDigit::PMT8inch); // 0 For PMTs
                    }
                }
            //}
            if (IsData) {
                it_tank_data++;
                if (it_tank_data == (*Hits).end()) loop_tank = false;
            }
            else {
                it_tank_mc++;
                if (it_tank_mc == (*MCHits).end()) loop_tank = false;
            }
        //}
        Log("tank loop processed.", v_debug, verbosity);
    }

    RecoVertex* truevtx = 0;
    auto get_muonMC = m_data->Stores.at("RecoEvent")->Get("TrueVertex", truevtx);
    fTrueVtxX = truevtx->GetPosition().X();
    fTrueVtxY = truevtx->GetPosition().Y();
    fTrueVtxZ = truevtx->GetPosition().Z();
    fTrueVtxTime = truevtx->GetTime();
    fTrueDirX = truevtx->GetDirection().X();
    fTrueDirY = truevtx->GetDirection().Y();
    fTrueDirZ = truevtx->GetDirection().Z();

    return;
}
