#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <calotreegen/caloTreeGen.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffarawobjects.so)
R__LOAD_LIBRARY(libcaloTreeGen.so)
#endif

void Fun4All_CaloTreeGen(const int nEvents = 0, const char *listFile = "fileList_withGeo_timingCut_Template.list", const char *inName = "commissioning.root")
{
    std::cout << "working ... " << std::endl;
  Fun4AllServer *se = Fun4AllServer::instance();
  recoConsts *rc = recoConsts::instance();

  gSystem->Load("libg4dst");

  caloTreeGen *calo = new caloTreeGen(inName);
  //What subsystems do you want?
  /*calo -> doEMCal(1,"TOWERINFO_CALIB_CEMC");
  //Store EMCal clusters?
  calo -> doClusters(1,"CLUSTERINFO_CEMC");

  //Store tower information for each EMCal cluster?
  calo -> doClusterDetails(1);
  
  //Store HCal information?
  calo -> doHCals(0,"TOWERINFO_CALIB_HCALOUT","TOWERINFO_CALIB_HCALIN");

  //Store ZDC information?
  calo -> doZDC(0,"TOWERINFO_CALIB_ZDC");

  //Store GL1 Information?
  calo -> doTrig(0,"GL1Packet");*/
  
  
  se->registerSubsystem(calo);

  Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTcalo");

  //in->AddFile(listFile);    // condor
  in->AddListFile(listFile, 1);      // local

  se->registerInputManager(in);

  se->run(nEvents);
  se->End();
  se->PrintTimer();
  std::cout << "All done!" << std::endl;

  gSystem->Exit(0);
}
