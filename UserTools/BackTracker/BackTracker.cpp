#include "BackTracker.h"

BackTracker::BackTracker():Tool(){}

// To sort
struct sort_by_charge {
    bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
        return left.second < right.second;
    }
};


bool BackTracker::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  // Load my config parameters
  bool gotVerbosity = m_variables.Get("verbosity",verbosity);
  bool gotStore = m_variables.Get("EventStore",eventStore);
  if (!gotVerbosity) {
    verbosity = 0;
    logmessage = "BackTracker::Initialize: \"verbosity\" not set in the config, defaulting to 0";
    Log(logmessage, v_error, verbosity);
  }
  if (!(eventStore=="ANNIEEvent"||eventStore=="RecoEvent")) {
      eventStore="ANNIEEvent";
      Log("BackTracker::Initialize: \"EventStore\" not set in the config, defaulting to ANNIEEvent",v_error,verbosity);
  }


  // Set up the pointers we're going to save. No need to 
  // delete them at Finalize, the store will handle it
  fClusterToBestParticleID  = new std::map<double, int>;
  fClusterToBestParticlePDG = new std::map<double, int>;
  fClusterEfficiency        = new std::map<double, double>;
  fClusterPurity            = new std::map<double, double>;
  fClusterTotalCharge       = new std::map<double, double>;
  
  return true;
}

//------------------------------------------------------------------------------
bool BackTracker::Execute()
{
  if (!LoadFromStores())
    return false;

  fClusterToBestParticleID ->clear();
  fClusterToBestParticlePDG->clear();
  fClusterEfficiency       ->clear();
  fClusterPurity           ->clear();
  fClusterTotalCharge      ->clear();

  fParticleToTankTotalCharge.clear();
  SumParticleTankCharge();

  
  // Loop over the clusters and do the things
  if(eventStore=="ANNIEEvent"){
  for (std::pair<double, std::vector<MCHit>>&& apair : *fClusterMapMC) {
    int prtId = -5;
    int prtPdg = -5;
    double eff = -5;
    double pur = -5;
    double totalCharge = 0;

    MatchMCParticle(apair.second, prtId, prtPdg, eff, pur, totalCharge);

    fClusterToBestParticleID ->emplace(apair.first, prtId);
    fClusterToBestParticlePDG->emplace(apair.first, prtPdg);
    fClusterEfficiency       ->emplace(apair.first, eff);
    fClusterPurity           ->emplace(apair.first, pur);
    fClusterTotalCharge      ->emplace(apair.first, totalCharge);
  }

  m_data->Stores.at("ANNIEEvent")->Set("ClusterToBestParticleID",  fClusterToBestParticleID );
  m_data->Stores.at("ANNIEEvent")->Set("ClusterToBestParticlePDG", fClusterToBestParticlePDG);
  m_data->Stores.at("ANNIEEvent")->Set("ClusterEfficiency",        fClusterEfficiency       );
  m_data->Stores.at("ANNIEEvent")->Set("ClusterPurity",            fClusterPurity           );
  m_data->Stores.at("ANNIEEvent")->Set("ClusterTotalCharge",       fClusterTotalCharge      );

  }
  else if (eventStore == "RecoEvent") {
      RecoCluster* i_cluster = nullptr;
      for (int i=0; i<fClusterList->size(); i++) {
          i_cluster=fClusterList->at(i);

          MatchMCRecoParticle(i_cluster);
      }
  }

  return true;
}

//------------------------------------------------------------------------------
bool BackTracker::Finalise()
{

  return true;
}

//------------------------------------------------------------------------------
void BackTracker::SumParticleTankCharge()
{

    if(eventStore=="ANNIEEvent"){
  for (auto mcHitsIt : *fMCHitsMap) {
    std::vector<MCHit> mcHits = mcHitsIt.second;
    for (uint mcHitIdx = 0; mcHitIdx < mcHits.size(); ++mcHitIdx) {

      // technically a MCHit could have multiple parents, but they don't appear to in practice
      // skip any cases we come across
      std::vector<int> parentIdxs = *(mcHits[mcHitIdx].GetParents());
      if (parentIdxs.size() != 1) continue;
      
      int particleId = -5;
      for (auto it : *fMCParticleIndexMap) {
	if (it.second == parentIdxs[0]) particleId = it.first;
      }
      if (particleId == -5) continue;
	
      double depositedCharge = mcHits[mcHitIdx].GetCharge();      
      if (!fParticleToTankTotalCharge.count(particleId)) 
	fParticleToTankTotalCharge.emplace(particleId, depositedCharge);
      else 
	fParticleToTankTotalCharge.at(particleId) += depositedCharge;
    }    
  }
    }
    else if (eventStore == "RecoEvent") {
        RecoDigit i_digit;
        for (int i=0;i<fDigitList->size();i++) {
            i_digit=fDigitList->at(i);
            std::vector<int> parentIdxs = i_digit.GetParents();
            if (parentIdxs.size() != 1) {
                Log("BackTracker: Digit has more than one parent!",v_debug,verbosity);
                continue;
            }
            int particleId = -5;
            particleId=parentIdxs.at(0);
            if(particleId==-5) continue;
            double depositedCharge = i_digit.GetCalCharge();
            if (!fParticleToTankTotalCharge.count(particleId)) {
                fParticleToTankTotalCharge.emplace(particleId,depositedCharge);
            }
            else {
                fParticleToTankTotalCharge.at(particleId)+=depositedCharge;
            }
        }
    }
}

//------------------------------------------------------------------------------
void BackTracker::MatchMCParticle(std::vector<MCHit> const &mchits, int &prtId, int &prtPdg, double &eff, double &pur, double &totalCharge)
{
  // Loop over the hits and get all of their parents and the energy that each one contributed
  //  be sure to bunch up all neutronic contributions
  std::map<int, double> mapParticleToTotalClusterCharge;
  totalCharge = 0;

  for (auto mchit : mchits) {    
    std::vector<int> parentIdxs = *(mchit.GetParents());
    if (parentIdxs.size() != 1) {
      logmessage = "BackTracker::MatchMCParticle: this MCHit has ";
      logmessage += std::to_string(parentIdxs.size()) + " parents!";
      Log(logmessage, v_debug, verbosity);
      continue;
    }
    
    int particleId = -5;
    for (auto it : *fMCParticleIndexMap) {
      if (it.second == parentIdxs[0]) particleId = it.first;
    }
    if (particleId == -5) continue;
    
    double depositedCharge = mchit.GetCharge();
    totalCharge += depositedCharge;
    
    if (mapParticleToTotalClusterCharge.count(particleId) == 0) 
      mapParticleToTotalClusterCharge.emplace(particleId, depositedCharge);
    else
      mapParticleToTotalClusterCharge[particleId] += depositedCharge;    
  }       

  // Loop over the particleIds to find the primary contributer to the cluster
  double maxCharge = 0;
  for (auto apair : mapParticleToTotalClusterCharge) {
    if (apair.second > maxCharge) {
      maxCharge = apair.second;
      prtId = apair.first;
    }
  }

  // Check that we have some charge, if not then something is wrong so pass back all -5
  if (totalCharge > 0) {
      Log("Check 2c: prtId, fMCParticleIndexMap, fParticleToTankTotalCharge, mapParticleToTotalClusterCharge size: " + to_string(prtId) + " " + to_string(fMCParticleIndexMap->size()) + " " + to_string(fParticleToTankTotalCharge.size()) + " " + to_string(mapParticleToTotalClusterCharge.size()), v_debug, verbosity);
    eff = maxCharge/fParticleToTankTotalCharge.at(prtId);
    pur = maxCharge/totalCharge;
    prtPdg = (fMCParticles->at(fMCParticleIndexMap->at(prtId))).GetPdgCode();
  } else {
    prtId = -5;
    eff = -5;
    pur = -5;
    totalCharge = -5;
  }

  logmessage = "BackTracker::MatchMCParticle: best particleId is : ";
  logmessage += std::to_string(prtId) + " which has PDG: " + std::to_string(prtPdg);
  Log(logmessage, v_message, verbosity);

}

//------------------------------------------------------------------------------
bool BackTracker::LoadFromStores()
{
  // Grab the stuff we need from the stores

    if(eventStore=="ANNIEEvent"){
  bool goodMCClusters = m_data->CStore.Get("ClusterMapMC", fClusterMapMC);
  if (!goodMCClusters) {
    std::cerr<<"BackTracker: no ClusterMapMC in the CStore!"<<endl;
    return false;
  }

  bool goodAnnieEvent = m_data->Stores.count("ANNIEEvent");
  if (!goodAnnieEvent) {
    std::cerr<<"BackTracker: no ANNIEEvent store!"<<endl;
    return false;
  }
    
  bool goodMCHits = m_data->Stores.at("ANNIEEvent")->Get("MCHits", fMCHitsMap);
  if (!goodMCHits) {
    std::cerr<<"BackTracker: no MCHits in the ANNIEEvent!"<<endl;
    return false;
  }
    }
  bool goodMCParticles = m_data->Stores.at("ANNIEEvent")->Get("MCParticles", fMCParticles);
  if (!goodMCParticles) {
    std::cerr<<"BackTracker: no MCParticles in the ANNIEEvent!"<<endl;
    return false;
  }

  bool goodMCParticleIndexMap = m_data->Stores.at("ANNIEEvent")->Get("TrackId_to_MCParticleIndex", fMCParticleIndexMap);
  if (!goodMCParticleIndexMap) {
    std::cerr<<"BackTracker: no TrackId_to_MCParticleIndex in the ANNIEEvent!"<<endl;
    return false;
  }

    else if (eventStore == "RecoEvent") {
        bool goodRecoDigits = m_data->Stores.at("RecoEvent")->Get("RecoDigit",fDigitList);
        if (!goodRecoDigits) {
            Log("Backtracker: no RecoDigits in the RecoEvent!",v_error,verbosity);
        }
        bool goodRecoClusters = m_data->Stores.at("RecoEvent")->Get("RecoClusters",fClusterList);
        if (!goodRecoClusters) {
            Log("Backtracker: no RecoClusters in the RecoEvent!",v_error,verbosity);
            return false;
        }
    }

  return true;
}

void BackTracker::MatchMCRecoParticle(RecoCluster* cluster) {
    // Loop over the hits and get all of their parents and the energy that each one contributed
  //  be sure to bunch up all neutronic contributions
    std::map<int, double> mapParticleToTotalClusterCharge;
    double totalCharge = 0;
    double eff=-5;
    double pur=-5;
    int prtId=-5;
    int prtPdg=-5;

    std::vector<RecoDigit*> cluster_digits=cluster->GetDigitList();
    RecoDigit* i_digit=nullptr;
    for (int i = 0; i<cluster_digits.size();i++) {

        i_digit=cluster_digits.at(i);
        std::vector<int> parentIdxs = i_digit->GetParents();
        if (parentIdxs.size() != 1) {
            logmessage = "BackTracker::MatchMCParticle: this MCHit has ";
            logmessage += std::to_string(parentIdxs.size()) + " parents!";
            Log(logmessage, v_debug, verbosity);
            continue;
        }

        int particleId = -5;
        Log("Check 1: fMCParticleIndexMap size: "+to_string(fMCParticleIndexMap->size()));
        for (auto it : *fMCParticleIndexMap) {
            if (it.second == parentIdxs[0]) particleId = it.first;
        }
        if (particleId == -5) {
            Log("BackTracker: parent identification failed!",v_debug,verbosity);
            continue;
        }

        double depositedCharge = i_digit->GetCalCharge();
        totalCharge += depositedCharge;

        if (mapParticleToTotalClusterCharge.count(particleId) == 0)
            mapParticleToTotalClusterCharge.emplace(particleId, depositedCharge);
        else
            mapParticleToTotalClusterCharge[particleId] += depositedCharge;
    }
    Log("Check 2", v_debug, verbosity);
    i_digit=nullptr;
        delete i_digit;

    // Loop over the particleIds to find the primary contributer to the cluster
    double maxCharge = 0;
    for (auto apair : mapParticleToTotalClusterCharge) {
        if (apair.second > maxCharge) {
            maxCharge = apair.second;
            prtId = apair.first;
        }
    }

    // Check that we have some charge, if not then something is wrong so pass back all -5
    if (totalCharge > 0) {
        Log("Check 2c: prtId, fMCParticleIndexMap, fParticleToTankTotalCharge, mapParticleToTotalClusterCharge size: "+to_string(prtId)+" "+to_string(fMCParticleIndexMap->size())+ " " + to_string(fParticleToTankTotalCharge.size()) + " " + to_string(mapParticleToTotalClusterCharge.size()), v_debug, verbosity);
        eff = maxCharge / fParticleToTankTotalCharge.at(prtId);
        pur = maxCharge / totalCharge;
        prtPdg = (fMCParticles->at(fMCParticleIndexMap->at(prtId))).GetPdgCode();
        Log("Check 2a: eff, pur: "+to_string(eff)+" "+to_string(pur), v_debug, verbosity);
    }
    else {
        prtId = -5;
        eff = -5;
        pur = -5;
        totalCharge = -5;
        Log("Check 2b",v_debug,verbosity);

    }

    cluster->SetParticle(prtId,prtPdg,eff,pur);

    logmessage = "BackTracker::MatchMCParticle: best particleId is : ";
    logmessage += std::to_string(prtId) + " which has PDG: " + std::to_string(prtPdg);
    Log(logmessage, v_message, verbosity);
}
