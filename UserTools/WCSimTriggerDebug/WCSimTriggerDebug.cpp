#include "WCSimTriggerDebug.h"

WCSimTriggerDebug::WCSimTriggerDebug():Tool(){}


bool WCSimTriggerDebug::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  return true;
}


bool WCSimTriggerDebug::Execute(){

  return true;
}


bool WCSimTriggerDebug::Finalise(){

  return true;
}
