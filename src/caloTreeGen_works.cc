#include "caloTreeGen.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

//Fun4All
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <ffaobjects/EventHeader.h>

//ROOT stuff
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>

//For clusters and geometry
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerGeomContainer.h>

//Tower stuff
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoDefs.h>

//GL1 Information
#include <ffarawobjects/Gl1Packet.h>

//for cluster vertex correction
#include <CLHEP/Geometry/Point3D.h>

//for the vertex
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

// hot towers 

//____________________________________________________________________________..
// THIS IS THE CODE FOR ETHANS DDPDC

caloTreeGen::caloTreeGen(const std::string &name):
SubsysReco("CaloTreeGen")
  ,T(nullptr)
  ,Outfile(name)
{
  std::cout << "caloTreeGen::caloTreeGen(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
caloTreeGen::~caloTreeGen()
{
  std::cout << "caloTreeGen::~caloTreeGen() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int caloTreeGen::Init(PHCompositeNode *topNode)
{
  
  out = new TFile(Outfile.c_str(),"RECREATE");
  h_etaTow = new TH1F("h_etaTow", "eta", 8*96, -0.5, 95.5);
  h_phiTow = new TH1F("h_phiTow", "phi", 16, -0.5, 1.5); 
  h_position = new TH2F("h_position", "Position; Eta; Phi", 8*96, -0.5, 95.5, 16, -0.5, 1.5);

  for (int i = 0; i < 16; i++) {      // just sets up the 16 2D hisotgrams (one for each phi bin)
    h_pi0[i] = new TH2F(Form("h_pi0%i", i), Form("pi0 candidates in phi bin %i", i), bins_etabin, low_etabin, high_etabin, bins_pi0, low_pi0, high_pi0);
    h_pi0[i]->GetXaxis()->SetTitle("Eta");
    h_pi0[i]->GetYaxis()->SetTitle("Pi0 Mass (GeV)");
    h_pi0_onblock[i] = new TH2F(Form("h_pi0_onblock%i", i), Form("pi0 candidates in phi bin %i, conflated eta", i), 16, 0, 16, bins_pi0, low_pi0, high_pi0);
    h_pi0_onblock[i]->GetXaxis()->SetTitle("Eta");
    h_pi0_onblock[i]->GetYaxis()->SetTitle("Pi0 Mass (GeV)");
  }

  
  T = new TTree("T","T");
  
  //Electromagnetic Calorimeter
  if(storeEMCal)
    {
      /*T -> Branch("emcTowE",&m_emcTowE);
      T -> Branch("emcTowiEta",&m_emciEta);
      T -> Branch("emcTowiPhi",&m_emciPhi);
      T -> Branch("emcTime",&m_emcTime);
      T -> Branch("emcChi2",&m_emcChi2);
      T -> Branch("emcPed",&m_emcPed);*/
      //EMCal Cluster information
      if(storeEMCal && storeClusters)
	{
	  /*T -> Branch("clusterE",&m_clusterE);*/
	  T -> Branch("clusterPhi",&m_clusterPhi);
	  T -> Branch("clusterEta", &m_clusterEta);
	  T -> Branch("clusterEta_act", &m_clusterEta_act);
	  T -> Branch("clusterPhi_act", &m_clusterPhi_act);
	  T -> Branch("clusterPt_Ecore", &m_clusterPt_Ecore);
	  T -> Branch("clusterPt_Ecore", &m_clusterPt_Ecore);
	  //T -> Branch("clusterPt_E", &m_clusterPt_E);
	  T -> Branch("clusterChi2", &m_clusterChi);
	  //T -> Branch("clusterNtow",&m_clusterNtow);
	 // T -> Branch("clusterTowMaxE",&m_clusterTowMaxE);
	  T -> Branch("clusterECore",&m_clusterECore);
      //T -> Branch("clusterNtow",&m_clusterNtow);
      T->Branch("pi0_phi",   &pi0_phi_vec);
      T->Branch("pi0_eta",   &pi0_eta_vec);
      T->Branch("pi0_pt",    &pi0_pt_vec);
      T->Branch("pi0_mass",  &pi0_mass_vec);
      T->Branch("asym",      &asym_vec);
      T->Branch("chi2_max",  &chi2_max_vec);
      T->Branch("pi0clus1_etatow", &pi0clus1_etatow_vec);
      T->Branch("pi0clus1_phitow", &pi0clus1_phitow_vec);
      T->Branch("pi0clus1_E", &pi0clus1_E_vec);
      T->Branch("pi0clus2_etatow", &pi0clus2_etatow_vec);
      T->Branch("pi0clus2_phitow", &pi0clus2_phitow_vec);
      T->Branch("pi0clus2_E", &pi0clus2_E_vec);

	  //Information for towers within clusters
	  //Enabled by setting "DoFineClusters" in the macro
	  if(storeEMCal && storeClusters && storeClusterDetails)
	    {
	      T -> Branch("clusTowPhi","vector<vector<int> >",&m_clusTowPhi);
	      T -> Branch("clusTowEta","vector<vector<int> >",&m_clusTowEta);
	      T -> Branch("clusTowE","vector<vector<float> >",&m_clusTowE);
	    }
	}
    }
  //Outer Hadronic Calorimeter
  if(storeHCals)
    {
      T -> Branch("ohcTowE",&m_ohcTowE);
      T -> Branch("ohcTowiEta",&m_ohciTowEta);
      T -> Branch("ohcTowiPhi",&m_ohciTowPhi);
      T -> Branch("ohcTime",&m_ohcTime);
      T -> Branch("ohcChi2",&m_ohcChi2);
      T -> Branch("ohcPed",&m_ohcPed);
  
      //Inner Hadronic Calorimeter
      T -> Branch("ihcTowE",&m_ihcTowE);
      T -> Branch("ihcTowiEta",&m_ihciTowEta);
      T -> Branch("ihcTowiPhi",&m_ihciTowPhi);
      T -> Branch("ihcTime",&m_ihcTime);
      T -> Branch("ihcChi2",&m_ihcChi2);
      T -> Branch("ihcPed",&m_ihcPed);
    }
  //ZDC information
  if(storeZDC)
    {
      T -> Branch("zdcTowE",&m_zdcTowE);
      T -> Branch("zdcTowside",&m_zdcSide);
  
      //SMD information
      T -> Branch("smdE",&m_smdE);
      T -> Branch("smdSide",&m_smdSide);
    }
  //Total 
  T -> Branch("totalCaloEEMCal",&totalCaloEEMCal);
  /*
  T -> Branch("totalCaloEOHCal",&totalCaloEOHCal);
  T -> Branch("totalCaloEIHCal",&totalCaloEIHCal);
  T -> Branch("totalCaloEZDC",&totalCaloEZDC);
  */
  T -> Branch("isminbias",&isminbias);
  T -> Branch("zvertex",&m_vertex);

  //GL1 Information
  if(storeTrig)T -> Branch("trigscaledvec",&gl1_scaledvec);
  
  

 //so that the histos actually get written out
  Fun4AllServer *se = Fun4AllServer::instance();
  se -> Print("NODETREE"); 
  
  std::cout << "caloTreeGen::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int caloTreeGen::InitRun(PHCompositeNode *topNode)
{
  

  std::cout << "caloTreeGen::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

static float getAvgEta(const std::vector<float> &toweretas, const std::vector<float> &towerenergies) {                  
    float etamult = 0;                                                                                                  
    float etasum = 0;                                                                                                   
                                                                                                                                    
    for (UInt_t j = 0; j < towerenergies.size(); j++) {                                                                 
        float energymult = towerenergies.at(j) * toweretas.at(j);                                                       
        etamult += energymult;                                                                                      
        etasum += towerenergies.at(j);                                                                                                                                      
    }                                                                                                                   
    return etamult / etasum;                                                                                            
}    

static float getAvgPhi(const std::vector<float> &towerphis, const std::vector<float> &towerenergies, Int_t nphibin) {   
    float phimult = 0;                                                                                                  
    float phisum  = 0;                                                                                                  
                                                                                                                                    
    for (UInt_t j = 0; j < towerenergies.size(); j++) {                                                                 
        int phibin = towerphis.at(j);                                                                                   
        if (phibin - towerphis.at(0) < -nphibin / 2.0) {                                                                
            phibin += nphibin;                                                                                      
        }                                                                                                               
        else if (phibin - towerphis.at(0) > +nphibin / 2.0) {                                                           
            phibin -= nphibin;
        }                                                                                                               
        assert(std::abs(phibin - towerphis.at(0)) <= nphibin / 2.0);                                                    
        float energymult = towerenergies.at(j) * phibin;                                                                
        phimult += energymult;                                                                                          
        phisum += towerenergies.at(j);                                                                                  
    }                                                                                                                   
    float avgphi = phimult / phisum;
        if (avgphi < 0) {                                                                                                   
            avgphi += nphibin;                                                                                              
        }    
        avgphi = fmod(avgphi, nphibin);
        if (avgphi >= 255.5) {                                                                                              
            avgphi -= nphibin;                                                                                              
        }
        avgphi = fmod(avgphi + 0.5, 2) - 0.5;
        return avgphi;
}


//____________________________________________________________________________..
int caloTreeGen::process_event(PHCompositeNode *topNode)
{
    isminbias = 1;
   _gl1_packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    // check if trigger is #11 (MBD coincidence)
    if (!_gl1_packet) isminbias = 0; 
    gl1_scaledvec = _gl1_packet->lValue(0, "ScaledVector");
    std::bitset<64> bits(gl1_scaledvec);
    if (!bits.test(10) and !bits.test(12)) isminbias = 0;
    if (bits.test(10)) isminbias =10;
    if (bits.test(12)) isminbias =12;

    //if (isminbias == 1) std::cout << "minbias, gl1_scaledvec (bits): " << bits.to_string() << std::endl;
    //if (isminbias != 1) std::cout << "not minbias, gl1_scaledvec (bits): " << bits.to_string() << std::endl;

  m_vertex = -9999;
  GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
    {
      std::cout << "GlobalVertexMap node is missing" << std::endl;
    }
  if (vertexmap && !vertexmap->empty())
    {
      GlobalVertex* vtx = vertexmap->begin()->second;
      if (vtx)
	{
	  m_vertex = vtx->get_z();
	  //zVertex -> Fill(m_vertex);
	}
    }
  

   //Information on clusters
  RawClusterContainer *clusterContainer = findNode::getClass<RawClusterContainer>(topNode,m_clusterNode.c_str());
  if(!clusterContainer && storeClusters)
    {
      std::cout << PHWHERE << "caloTreeGen::process_event: "<<  m_clusterNode << " node is missing. Output related to this node will be empty" << std::endl;
      return 0;
    }
  
  //tower information
  TowerInfoContainer *emcTowerContainer;
  emcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode,m_emcTowerNode.c_str());
  if(!emcTowerContainer && storeEMCal)
    {
      std::cout << PHWHERE << "caloTreeGen::process_event: "<< m_emcTowerNode << " node is missing. Output related to this node will be empty" << std::endl;
    }
RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");                                                                      
if (!towergeom)                                                                                                       
{                                                                                                                     
    std::cout << PHWHERE << "caloTreeGen::process_event Could not find node TOWERGEOM_CEMC"  << std::endl;              
    return Fun4AllReturnCodes::ABORTEVENT;                                                                              
}
const int nphibin = towergeom->get_phibins(); 

  
  //grab all the towers and fill their energies. 
  unsigned int tower_range = 0;
  if(storeEMCal && emcTowerContainer)
    {
      tower_range = emcTowerContainer->size();
      totalCaloEEMCal = 0;
      for(unsigned int iter = 0; iter < tower_range; iter++)
	{
	  unsigned int towerkey = emcTowerContainer->encode_key(iter);
	  unsigned int ieta = TowerInfoDefs::getCaloTowerEtaBin(towerkey);
	  unsigned int iphi = TowerInfoDefs::getCaloTowerPhiBin(towerkey);
	  double energy = emcTowerContainer -> get_tower_at_channel(iter) -> get_energy();
	  int time =  emcTowerContainer -> get_tower_at_channel(iter) -> get_time();
	  float chi2 = emcTowerContainer -> get_tower_at_channel(iter) -> get_chi2();
	  float pedestal = emcTowerContainer -> get_tower_at_channel(iter) -> get_pedestal();
	  totalCaloEEMCal += energy;
	  m_emciPhi.push_back(iphi);
	  m_emciEta.push_back(ieta);
	  m_emcTowE.push_back(energy);
	  m_emcTime.push_back(time);
	  m_emcChi2.push_back(chi2);
	  m_emcPed.push_back(pedestal);
	} 
    }
  


  if(storeClusters && storeEMCal)
    {
      //TH1F *h_invmass = new TH1F("h_invmass", "invariant diphoton mass", 100, 0, 1);
      RawClusterContainer::ConstRange clusterEnd = clusterContainer -> getClusters();
      RawClusterContainer::ConstIterator clusterIter;
      RawClusterContainer::ConstIterator clusterIter2;
      for(clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; clusterIter++)
	{
	  RawCluster *recoCluster = clusterIter -> second;

	  CLHEP::Hep3Vector vertex(0,0,0);
	  if(m_vertex != -9999)vertex.setZ(m_vertex);
	  CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*recoCluster, vertex);
	  CLHEP::Hep3Vector E_vec_cluster_Full = RawClusterUtility::GetEVec(*recoCluster, vertex);

	  float clusE = E_vec_cluster_Full.mag();
	  float clusEcore = E_vec_cluster.mag();
	  float clus_eta = E_vec_cluster.pseudoRapidity();
	  float clus_phi = E_vec_cluster.phi();
	  float clus_pt = E_vec_cluster.perp();
	  float clus_pt_E = E_vec_cluster_Full.perp();
	  float clus_chi = recoCluster -> get_chi2();
	  float nTowers = recoCluster ->getNTowers(); 
	  float maxTowerEnergy = getMaxTowerE(recoCluster,emcTowerContainer);

      if (clusE < clus1min) continue;
	  m_clusterE.push_back(clusE);
	  m_clusterECore.push_back(clusEcore);
	  m_clusterPhi.push_back(clus_phi);
	  m_clusterEta.push_back(clus_eta);
	  m_clusterPt_Ecore.push_back(clus_pt);
	  m_clusterPt_E.push_back(clus_pt_E);
	  m_clusterChi.push_back(clus_chi);
	  m_clusterNtow.push_back(nTowers);
	  m_clusterTowMaxE.push_back(maxTowerEnergy);
      TLorentzVector photon1;
      photon1.SetPtEtaPhiE(clus_pt, clus_eta, clus_phi, clusEcore);

      std::vector<float> toweretas; 
      std::vector<float> towerphis; 
      std::vector<float> towerenergies; 

        RawCluster::TowerConstRange towers = recoCluster->get_towers();                                                     
        RawCluster::TowerConstIterator toweriter;                                                                           

        for (toweriter = towers.first; toweriter != towers.second; ++toweriter)                                             
        {                                                                                                                   
            int iphi = RawTowerDefs::decode_index2(toweriter->first);  //index2 is phi in CYL                               
            int ieta = RawTowerDefs::decode_index1(toweriter->first);

            unsigned int towerkey = iphi + (ieta << 16U);                                                                   
            unsigned int towerindex = emcTowerContainer->decode_key(towerkey);                                              
            TowerInfo* towinfo = emcTowerContainer->get_tower_at_channel(towerindex);                                      
            double towerenergy = towinfo->get_energy(); 

            toweretas.push_back(ieta);                                                                                      
            towerphis.push_back(iphi);                                                                                      
            towerenergies.push_back(towerenergy);
        }                                                                                                                   
        float avgeta1 = getAvgEta(toweretas, towerenergies);                                                                    
        float avgphi1 = getAvgPhi(towerphis, towerenergies, nphibin);

        m_clusterEta_act.push_back(avgeta1);
        m_clusterPhi_act.push_back(avgphi1);


        



      for (clusterIter2 = clusterEnd.first; clusterIter2 != clusterEnd.second; clusterIter2++) {
        if (clusterIter == clusterIter2) continue;
        if (sqrt(pow(m_vertex, 2)) > vertexcut)continue;
        if (isminbias == 0)continue; 
        RawCluster* recoCluster2 = clusterIter2->second;
        CLHEP::Hep3Vector vertex2(0, 0, 0);
        CLHEP::Hep3Vector E_vec_cluster2 = RawClusterUtility::GetECoreVec(*recoCluster2, vertex2);
        float clus2E = E_vec_cluster2.mag();
        float clus2_eta = E_vec_cluster2.pseudoRapidity();
        float clus2_phi = E_vec_cluster2.phi();
        float clus2_pt = E_vec_cluster2.perp();
        float clus2_chisq = recoCluster2->get_chi2();
        //if (clus2_chisq > 4) continue;
        if (std::max(clusE, clus2E) < clus1min) continue;
        if (std::min(clusE, clus2E) < clus2min) continue;
        TLorentzVector photon2;
        photon2.SetPtEtaPhiE(clus2_pt, clus2_eta, clus2_phi, clus2E);
        float asym = sqrt(pow(clusEcore - clus2E, 2)) / (clusEcore + clus2E);
        //if (asym > maxalpha) continue;
        TLorentzVector pi0 = photon1 + photon2;
        if (pi0.M() > 1) continue;

        RawCluster::TowerConstRange towers2 = recoCluster2->get_towers();                                                   
        RawCluster::TowerConstIterator toweriter;                                                                           

        std::vector<float> toweretas2;                                                                                      
        std::vector<float> towerphis2;                                                                                      
        std::vector<float> towerenergies2;

        for (toweriter = towers2.first; toweriter != towers2.second; ++toweriter)                                           
        {                                                                                                                   
            int iphi = RawTowerDefs::decode_index2(toweriter->first);  //index2 is phi in CYL                               
            int ieta = RawTowerDefs::decode_index1(toweriter->first);  //index1 is eta in CYL

            unsigned int towerkey = iphi + (ieta << 16U);                                                                   
            unsigned int towerindex = emcTowerContainer->decode_key(towerkey);                                              
            TowerInfo* towinfo = emcTowerContainer->get_tower_at_channel(towerindex);                                      
            double towerenergy = towinfo->get_energy();

            toweretas2.push_back(ieta);                                                                                     
            towerphis2.push_back(iphi);                                                                                     
            towerenergies2.push_back(towerenergy);                                                                          
        }                                        
        float avgeta2 = getAvgEta(toweretas2, towerenergies2);                                                                    
        float avgphi2 = getAvgPhi(towerphis2, towerenergies2, nphibin);

        if (std::max(clusE, clus2E) < 0.5) continue;
        if (std::min(clusE, clus2E) < 0.5) continue;
        float chi2max = std::max(clus_chi, clus2_chisq);
        if (chi2max > 4) continue;
        if (asym > 0.7) continue;
        pi0_phi_vec.push_back(pi0.Phi());
        pi0_mass_vec.push_back(pi0.M());
        pi0_eta_vec.push_back(pi0.Eta());
        pi0_pt_vec.push_back(pi0.Pt());
        asym_vec.push_back(asym);
        chi2_max_vec.push_back(chi2max);
        pi0clus1_etatow_vec.push_back(avgeta1);
        pi0clus2_etatow_vec.push_back(avgeta2);
        pi0clus1_E_vec.push_back(clusE);
        pi0clus2_E_vec.push_back(clus2E);
        pi0clus1_phitow_vec.push_back(avgphi1);
        pi0clus2_phitow_vec.push_back(avgphi2);

        h_etaTow->Fill(avgeta1);
        h_etaTow->Fill(avgeta2);
        h_phiTow->Fill(avgphi1);
        h_phiTow->Fill(avgphi2);

        phibin1 = floor((8*(avgphi1+.5))); 
        phibin2 = floor((8*(avgphi2+.5)));
        etabin1 = floor((8*(avgeta1+.5)));
        etabin2 = floor((8*(avgeta2+.5)));

        h_position->Fill(avgeta1, avgphi1);
        h_position->Fill(avgeta2, avgphi2);

        h_pi0[phibin1]->Fill(etabin1, pi0.M());
        h_pi0[phibin2]->Fill(etabin2, pi0.M());
        h_pi0_onblock[phibin1]->Fill(etabin1%16, pi0.M());
        h_pi0_onblock[phibin2]->Fill(etabin2%16, pi0.M());




        //h_invmass->Fill(pi0.M());
        } // end pi0 stuff 


      
      if(storeClusterDetails)
	    {
	      m_clusTowPhi.push_back(returnClusterTowPhi(recoCluster,emcTowerContainer));
	      m_clusTowEta.push_back(returnClusterTowEta(recoCluster,emcTowerContainer));
	      m_clusTowE.push_back(returnClusterTowE(recoCluster,emcTowerContainer));
	    }		     
	}
    }

  
  //tower information
  TowerInfoContainer *ohcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode,m_ohcTowerNode);
  TowerInfoContainer *ihcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode,m_ihcTowerNode.c_str());
	
  if(!ohcTowerContainer || !ihcTowerContainer)
    {
      std::cout << PHWHERE << "caloTreeGen::process_event: " << m_ohcTowerNode << " or " << m_ohcTowerNode << " node is missing. Output related to this node will be empty" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  if(storeHCals && ohcTowerContainer && ihcTowerContainer)
    {
      tower_range = ohcTowerContainer->size();
      totalCaloEOHCal = 0;
      for(unsigned int iter = 0; iter < tower_range; iter++)
	{
	  unsigned int towerkey = ohcTowerContainer->encode_key(iter);
	  unsigned int ieta = TowerInfoDefs::getCaloTowerEtaBin(towerkey);
	  unsigned int iphi = TowerInfoDefs::getCaloTowerPhiBin(towerkey);
	  int time =  ohcTowerContainer -> get_tower_at_channel(iter) -> get_time();
	  float chi2 = ohcTowerContainer -> get_tower_at_channel(iter) -> get_chi2();
	  double energy = ohcTowerContainer -> get_tower_at_channel(iter) -> get_energy();
	  float pedestal = ohcTowerContainer -> get_tower_at_channel(iter) -> get_pedestal();

	  totalCaloEOHCal += energy;
	  m_ohcTowE.push_back(energy);
	  m_ohciTowEta.push_back(ieta);
	  m_ohciTowPhi.push_back(iphi);
	  m_ohcTime.push_back(time);
	  m_ohcChi2.push_back(chi2);
      	  m_ohcPed.push_back(pedestal);

	} 
    
      tower_range = ihcTowerContainer->size();
      totalCaloEIHCal = 0;
      for(unsigned int iter = 0; iter < tower_range; iter++)
	{
	  unsigned int towerkey = ihcTowerContainer->encode_key(iter);
	  unsigned int ieta = TowerInfoDefs::getCaloTowerEtaBin(towerkey);
	  unsigned int iphi = TowerInfoDefs::getCaloTowerPhiBin(towerkey);
	  int time =  ohcTowerContainer -> get_tower_at_channel(iter) -> get_time();
	  float chi2 = ohcTowerContainer -> get_tower_at_channel(iter) -> get_chi2();
	  double energy = ihcTowerContainer -> get_tower_at_channel(iter) -> get_energy();
	  float pedestal = ihcTowerContainer -> get_tower_at_channel(iter) -> get_pedestal();

	  totalCaloEIHCal += energy;
	  m_ihcTowE.push_back(energy);
	  m_ihciTowEta.push_back(ieta);
	  m_ihciTowPhi.push_back(iphi);
	  m_ihcTime.push_back(time);
	  m_ihcChi2.push_back(chi2);
	  m_ihcPed.push_back(pedestal);

	} 
    }
  
  TowerInfoContainer *zdcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode,m_zdcTowerNode.c_str());
  if(!zdcTowerContainer)
    {
      std::cout << PHWHERE << "caloTreeGen::process_event: " << m_emcTowerNode << " node is missing. Output related to this node will be empty" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  if(storeZDC && zdcTowerContainer)
    {
      tower_range = zdcTowerContainer ->size();
      totalCaloEZDC = 0;
      for(unsigned int iter = 0; iter < tower_range; iter++)
	{
	  if(iter < 16)
	    {
	      float energy = zdcTowerContainer -> get_tower_at_channel(iter) -> get_energy();
	      m_zdcTowE.push_back(energy);
	      unsigned int towerkey = zdcTowerContainer->encode_key(iter);
	      unsigned int side = TowerInfoDefs::get_zdc_side(towerkey);
	      m_zdcSide.push_back(side);
	      totalCaloEZDC += energy;
	    }
	  if(iter > 15 && iter < 48)
	    {
	      //smd north stuff
	      float energy = zdcTowerContainer -> get_tower_at_channel(iter) -> get_energy();
	      m_smdE.push_back(energy);
	      unsigned int towerkey = zdcTowerContainer->encode_key(iter);
	      unsigned int side = TowerInfoDefs::get_zdc_side(towerkey);
	      m_smdSide.push_back(side);
	    }
	}
    }

  
  /*Gl1Packet *gl1PacketInfo = findNode::getClass<Gl1Packet>(topNode,m_trigNode.c_str());
  if(!gl1PacketInfo && storeTrig)
    {
      std::cout << PHWHERE << "caloTreeGen::process_event: " << m_trigNode<< " node is missing. Output related to this node will be empty" << std::endl;
    }

  if(storeTrig && gl1PacketInfo)
    {
      uint64_t triggervec =  gl1PacketInfo->getTriggerVector();
      for(int i = 0; i < 64; i++)
	{
	  bool trig_decision = ((triggervec & 0x1U) == 0x1U);
	  m_triggerVector.push_back(trig_decision);
	  triggervec = (triggervec >> 1U) & 0xffffffffU;
	}
    }*/
    if (isminbias != 0) T -> Fill();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int caloTreeGen::ResetEvent(PHCompositeNode *topNode)
{
  m_clusterE.clear();
  m_clusterPhi.clear();
  m_clusterEta.clear();
  m_clusterPhi_act.clear();
  m_clusterEta_act.clear();
  m_clusterPt_E.clear();
  m_clusterPt_Ecore.clear();
  m_clusterChi.clear();
  m_clusterTowMaxE.clear();
  m_clusterNtow.clear();
  m_clusterECore.clear();

  m_emcTowE.clear();
  m_emciEta.clear();
  m_emciPhi.clear();
  m_emcTime.clear();
  m_emcChi2.clear();
  m_emcPed.clear();

  m_ihcTowE.clear();
  m_ihciTowEta.clear();
  m_ihciTowPhi.clear();
  m_ihcTime.clear();
  m_ihcChi2.clear();
  m_ihcPed.clear();

  m_ohcTowE.clear();
  m_ohciTowEta.clear();
  m_ohciTowPhi.clear();
  m_ohcTime.clear();
  m_ohcChi2.clear();
  m_ohcPed.clear();


  m_clusTowPhi.clear();
  m_clusTowEta.clear();
  m_clusTowE.clear();
  pi0_phi_vec.clear();
  pi0_mass_vec.clear();
  pi0_eta_vec.clear();
  pi0_pt_vec.clear();
  asym_vec.clear();
  chi2_max_vec.clear();
  pi0clus1_etatow_vec.clear();
  pi0clus2_etatow_vec.clear();
  pi0clus1_phitow_vec.clear();
  pi0clus2_phitow_vec.clear();
  pi0clus1_E_vec.clear();
  pi0clus2_E_vec.clear();
 
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int caloTreeGen::EndRun(const int runnumber)
{
  std::cout << "caloTreeGen::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int caloTreeGen::End(PHCompositeNode *topNode)
{
  std::cout << "caloTreeGen::End(PHCompositeNode *topNode) This is the End..." << std::endl;

  out -> cd();
  h_etaTow->Write();
  h_phiTow->Write();
  h_position->Write();
  for (int i = 0; i < 16; i++) h_pi0[i]->Write();
  for (int i = 0; i < 16; i++) h_pi0_onblock[i]->Write();

  //T -> Write();
  //zVertex -> Write();
  out -> Close();
  delete out;

  //hm -> dumpHistos(Outfile.c_str(), "UPDATE");
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int caloTreeGen::Reset(PHCompositeNode *topNode)
{
 std::cout << "caloTreeGen::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void caloTreeGen::Print(const std::string &what) const
{
  std::cout << "caloTreeGen::Print(const std::string &what) const Printing info for " << what << std::endl;
}
//____________________________________________________________________________..
float caloTreeGen::getMaxTowerE(RawCluster *cluster, TowerInfoContainer *towerContainer)
{
  RawCluster::TowerConstRange towers = cluster -> get_towers();
  RawCluster::TowerConstIterator toweriter;
  
  float maxEnergy = 0;
  for(toweriter = towers.first; toweriter != towers.second; toweriter++)
    {
      float towE = toweriter -> second;
   
      if( towE > maxEnergy)  maxEnergy = towE;
    }
  return maxEnergy;
}
//____________________________________________________________________________..
std::vector<int> caloTreeGen::returnClusterTowEta(RawCluster *cluster, TowerInfoContainer *towerContainer)
{
  RawCluster::TowerConstRange towers = cluster -> get_towers();
  RawCluster::TowerConstIterator toweriter;
  
  std::vector<int> towerIDsEta;
  for(toweriter = towers.first; toweriter != towers.second; toweriter++) towerIDsEta.push_back(RawTowerDefs::decode_index1(toweriter -> first));

  return towerIDsEta;
}
//____________________________________________________________________________..
std::vector<int> caloTreeGen::returnClusterTowPhi(RawCluster *cluster, TowerInfoContainer *towerContainer)
{
  RawCluster::TowerConstRange towers = cluster -> get_towers();
  RawCluster::TowerConstIterator toweriter;
  
  std::vector<int> towerIDsPhi;
  for(toweriter = towers.first; toweriter != towers.second; toweriter++) towerIDsPhi.push_back(RawTowerDefs::decode_index2(toweriter -> first));
  return towerIDsPhi;
}
//____________________________________________________________________________..
std::vector<float> caloTreeGen::returnClusterTowE(RawCluster *cluster, TowerInfoContainer *towerContainer)
{
  RawCluster::TowerConstRange towers = cluster -> get_towers();
  RawCluster::TowerConstIterator toweriter;
  
  std::vector<float> towerE;
  for(toweriter = towers.first; toweriter != towers.second; toweriter++) towerE.push_back(toweriter -> second);
  
  return towerE;
}
