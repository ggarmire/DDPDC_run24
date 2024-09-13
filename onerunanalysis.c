#include <complex>
#include <numeric>
#include <stdlib.h>
#include <cmath>
void onerunanalysis(string infile, string outfile, int event_start=0, int event_stop=0)
{
    TFile *f = TFile::Open(infile.c_str());     // file containing the Tree with pi0 info 
    TTree *T = (TTree*) f->Get("T");
    TFile *fout = TFile::Open(outfile.c_str(), "recreate");

    vector<float> *pt = {0};
    vector<float> *asym = {0};
    vector<float> *mass = {0};
    vector<float> *chi = {0};
    vector<float> *etaTow1 = {0};
    vector<float> *phiTow1 = {0};
    vector<float> *etaTow2 = {0};
    vector<float> *phiTow2 = {0};
    vector<float> *clusE = {0};
    vector<float> *pclus1E = {0};
    vector<float> *pclus2E = {0};
    vector<float> *clusphi = {0};
    vector<float> *cluseta = {0};
    //vector<float> *clusNTow = {0};
    float zVertex;
    float totalEEmcal;


// get info out of the tree: (check that the names are right)
    T->SetBranchAddress("pi0_pt", &pt);
    T->SetBranchAddress("pi0_mass", &mass);
    T->SetBranchAddress("asym", &asym);
    T->SetBranchAddress("clusterChi2", &chi);
    //T->SetBranchAddress("clusterECore", &clusE);
    T->SetBranchAddress("clusterEta_act", &cluseta);
    T->SetBranchAddress("clusterPhi_act", &clusphi);
    T->SetBranchAddress("pi0clus1_etatow", &etaTow1);
    T->SetBranchAddress("pi0clus1_phitow", &phiTow1);
    T->SetBranchAddress("pi0clus2_etatow", &etaTow2);
    T->SetBranchAddress("pi0clus2_phitow", &phiTow2);
    T->SetBranchAddress("pi0clus1_E", &pclus1E);
    T->SetBranchAddress("pi0clus1_E", &pclus2E);
    T->SetBranchAddress("totalCaloEEMCal", &totalEEmcal);
    T->SetBranchAddress("zvertex", &zVertex);
    //T->SetBranchAddress("clusterNtow", &clusNTow);
// set up histograms
    TH1F *h_clusE = new TH1F("h_clusE", "Energy Spectrum of clusters", 200, 0, 20);       //  genuinely dont know what this is 
    TH1F *h_clusE_bin385 = new TH1F("h_clusE_bin385", "Energy Spectrum of clusters in bin 385", 200, 0, 20);       //  genuinely dont know what this is 
    TH1F *h_clusE_mod1 = new TH1F("h_clusE_mod1", "Energy Spectrum of clusters in bins in same position as 385", 200, 0, 20);       //  genuinely dont know what th
    const int bins_etabin = 768;  
    const float high_etabin = 767.5;
    const float low_etabin = -0.5;
    
    TH2F *h_clusE_byeta = new TH2F("h_clusE_byeta", "Cluster energy spectra per eta bin", bins_etabin, low_etabin, high_etabin, 200, 0, 20);
    
// get number of events for the for loop 
    event_stop = T->GetEntries();
    cout << "running " << event_stop - event_start << " events" << endl;

    for (Long64_t entry=event_start; entry<event_stop; entry++)     // loop through events 
    //for (Long64_t entry=event_start; entry<10000; entry++)
    {
        T->GetEntry(entry);
        //cout << entry <<endl;
        int size = pclus1E->size();
        //cout << "size" << size << endl;
        for(int j = 0; j < size; j++){  // loop through clusters in the event
            if(!(abs(zVertex) <= 20)) {
                continue;       // vertex cut (since there is only one vertex per event you can move this to outside the cluster loop 
            }
            if(chi->at(j) > 4){
                continue;       // chi2 and asymmetry cut 
            }
            // cout << zVertex << ", " << chi->at(j) << ", " << asym->at(j) << endl;

           // cout << cluseta->size() << "etasize"  <<endl;
            h_clusE->Fill(pclus1E->at(j));
            h_clusE->Fill(pclus2E->at(j));
            int eta1int = floor((8*(etaTow1->at(j)+.5)));
            int eta2int = floor((8*(etaTow2->at(j)+.5)));
            if (eta1int==385) h_clusE_bin385->Fill(pclus1E->at(j));
            if (eta2int==385) h_clusE_bin385->Fill(pclus2E->at(j));
            //cout << "a"  <<endl;
            if (eta1int%16==1 and eta1int !=385) h_clusE_mod1->Fill(pclus1E->at(j));
            if (eta2int%16==1 and eta2int !=385) h_clusE_mod1->Fill(pclus2E->at(j));
        
        }// end pi0 candidate loop
    }//end event loop
    

    fout->Write();      // write histograms to outfile 
    fout->Close();

}
