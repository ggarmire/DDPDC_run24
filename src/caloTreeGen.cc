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
  h_position = new TH2F("h_position", "Position of all clusters going into diphoton spectra ; Eta; Phi", 8*96, -0.5625, 95.4375, 16, -0.5, 1.5);
  h_position_clus = new TH2F("h_position_clus", "Position of all clusters E>1GeV; Eta; Phi", 8*96, -0.5, 95.5, 16, -0.5, 1.5);
  h_position_clus_g1 = new TH2F("h_position_clus_g1", "Position of all clusters E>1GeV, Ntow>1; Eta; Phi", 8*96, -0.5625, 95.4375, 16, -0.5, 1.5);
  h_position_clus_n1 = new TH2F("h_position_clus_n1", "Position of single-tower clusters E>1GeV; Eta; Phi", 8*96, -0.5625, 95.4375, 16, -0.5, 1.5);
  h_fullpi0 = new TH1F("h_fullpi0", "Full Diphoton spectrum; mass [GeV]; counts", bins_pi0, low_pi0, high_pi0);

  for (int i = 0; i < 16; i++) {      // just sets up the 16 2D hisotgrams (one for each phi bin)
    h_pi0[i] = new TH2F(Form("h_pi0%i", i), Form("pi0 candidates in phi bin %i", i), bins_etabin, low_etabin, high_etabin, bins_pi0, low_pi0, high_pi0);
    h_pi0[i]->GetXaxis()->SetTitle("Eta");
    h_pi0[i]->GetYaxis()->SetTitle("Pi0 Mass (GeV)");
    h_pi0_onblock[i] = new TH2F(Form("h_pi0_onblock%i", i), Form("pi0 candidates in phi bin %i, conflated eta", i), 16, 0, 16, bins_pi0, low_pi0, high_pi0);
    h_pi0_onblock[i]->GetXaxis()->SetTitle("Eta");
    h_pi0_onblock[i]->GetYaxis()->SetTitle("Pi0 Mass (GeV)");
  }


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
 /*   isminbias = 1;
   _gl1_packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    // check if trigger is #11 (MBD coincidence)
    if (!_gl1_packet) isminbias = 0; 
    gl1_scaledvec = _gl1_packet->lValue(0, "ScaledVector");
    std::bitset<64> bits(gl1_scaledvec);
    if (!bits.test(10) and !bits.test(12)) isminbias = 0;
    if (bits.test(10)) isminbias =10;
    if (bits.test(12)) isminbias =12;*/

   _gl1_packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    if(_gl1_packet) b_gl1_scaledvec = _gl1_packet->lValue(0, "ScaledVector");

    isminbias = 0;
    std::vector<int> trig_bits;
    std::bitset<64> bits(b_gl1_scaledvec);
    for (unsigned int b =0; b < 64; b++) {
        if ((b_gl1_scaledvec >> b) & 0x1U) trig_bits.push_back(b);
        }

    for (const int &bit : trig_bits){
        if ((bit == 10) or (bit ==12)) isminbias = bit;
    }

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
    if (abs(m_vertex) > 20) isminbias = 0;

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
	  //unsigned int towerkey = emcTowerContainer->encode_key(iter);
	  //unsigned int ieta = TowerInfoDefs::getCaloTowerEtaBin(towerkey);
	  //unsigned int iphi = TowerInfoDefs::getCaloTowerPhiBin(towerkey);
	  double energy = emcTowerContainer -> get_tower_at_channel(iter) -> get_energy();
	  //int time =  emcTowerContainer -> get_tower_at_channel(iter) -> get_time();
	  //float chi2 = emcTowerContainer -> get_tower_at_channel(iter) -> get_chi2();
	  //float pedestal = emcTowerContainer -> get_tower_at_channel(iter) -> get_pedestal();
	  totalCaloEEMCal += energy;
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
      if (isminbias == 0) continue;

	  CLHEP::Hep3Vector vertex(0,0,0);
	  if(m_vertex != -9999)vertex.setZ(m_vertex);
	  CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*recoCluster, vertex);
	  CLHEP::Hep3Vector E_vec_cluster_Full = RawClusterUtility::GetEVec(*recoCluster, vertex);

	  float clusE = E_vec_cluster_Full.mag();
	  float clusEcore = E_vec_cluster.mag();
	  float clus_eta = E_vec_cluster.pseudoRapidity();
	  float clus_phi = E_vec_cluster.phi();
	  float clus_pt = E_vec_cluster.perp();
	  //float clus_pt_E = E_vec_cluster_Full.perp();
	  float clus_chi = recoCluster -> get_chi2();
      float nTowers = recoCluster ->getNTowers(); 
	  //float maxTowerEnergy = getMaxTowerE(recoCluster,emcTowerContainer);

      if (clusE < clus1min) continue;
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
        h_position_clus ->Fill(avgeta1, avgphi1);
        if (nTowers > 1) h_position_clus_g1->Fill(avgeta1, avgphi1);
        if (nTowers == 1) h_position_clus_n1->Fill(avgeta1, avgphi1);

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

        h_fullpi0->Fill(pi0.M());
        //std::cout<< pi0.M()<<std::endl;



        //h_invmass->Fill(pi0.M());
        } // end pi0 stuff 


      
	}
    }

  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int caloTreeGen::ResetEvent(PHCompositeNode *topNode)
{
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
  h_position_clus->Write();
  h_position_clus_g1->Write();
  h_position_clus_n1->Write();
  h_fullpi0->Write();
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
