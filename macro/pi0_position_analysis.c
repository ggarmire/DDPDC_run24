#include <complex>
#include <numeric>
#include <stdlib.h>
#include <cmath>
void pi0_position_analysis(string infile, string outfile, int event_start=0, int event_stop=0)
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
    //vector<float> *clusNTow = {0};
    float zVertex;
    float totalEEmcal;


// get info out of the tree: (check that the names are right)
    T->SetBranchAddress("pi0_pt", &pt);
    T->SetBranchAddress("pi0_mass", &mass);
    T->SetBranchAddress("asym", &asym);
    T->SetBranchAddress("chi2_max", &chi);
    T->SetBranchAddress("pi0clus1_etatow", &etaTow1);
    T->SetBranchAddress("pi0clus1_phitow", &phiTow1);
    T->SetBranchAddress("pi0clus2_etatow", &etaTow2);
    T->SetBranchAddress("pi0clus2_phitow", &phiTow2);
    T->SetBranchAddress("totalCaloEEMCal", &totalEEmcal);
    T->SetBranchAddress("zvertex", &zVertex);
    //T->SetBranchAddress("clusterNtow", &clusNTow);
// set up histograms
    TH2F *h_pi0[16];        // 16 2D histograms contaning (Eta, pi0mass) of each pi0
    TH1F *h_etaTow = new TH1F("h_etaTow", "eta", 8*96, -0.5, 95.5);     // distribution of all the eta tower ids 
    TH1F *h_phiTow = new TH1F("h_phiTow", "phi", 16, -0.5, 1.5);       //distribution of all the phi tower ids
    TH1F *h_etaEntries = new TH1F("h_etaEntries", "Eta Entries; Invarient Mass (GeV); Counts", 60, -0.5, 10);       //  genuinely dont know what this is 
    TH2F *h_position = new TH2F("h_position", "Position; Eta; Phi", 8*96, -0.5, 95.5, 16, -0.5, 1.5);      // 2D of eta phi position of every cluster 

    const int bins_etabin = 768;  
    const float high_etabin = 767.5;
    const float low_etabin = -0.5;
    const int bins_pi0 = 100;       // can change this to get less bins if needed (might be helpful with low statistics_
    const float high_pi0 = 1;
    const float low_pi0 = 0.0;

    for (int i = 0; i < 16; i++) {      // just sets up the 16 2D hisotgrams (one for each phi bin)
        h_pi0[i] = new TH2F(Form("h_pi0%i", i), Form("pi0 candidates in phi bin %i", i), bins_etabin, low_etabin, high_etabin, bins_pi0, low_pi0, high_pi0);
        h_pi0[i]->GetXaxis()->SetTitle("Eta");
        h_pi0[i]->GetYaxis()->SetTitle("Pi0 Mass (GeV)");
    //    h_pi0_onBlock[i] = new TH2F(Form("h_pi0_onBlock%i", i), Form("pi0 candidates in phi bin %i", i), 16, low_etabin, 15.5, bins_pi0, low_pi0, high_pi0);
    }
    // not sure this next one gets used: 
    TH2F *h_mass_mean = new TH2F("h_mass_mean", "Mean mass of pi0 candidates; Eta Bin ID; Phi Bin ID; Mass (GeV)", 768, -0.5, 767.5, 16, -0.5, 15.5);
    
// get number of events for the for loop 
    event_stop = T->GetEntries();
    cout << "running " << event_stop - event_start << " events" << endl;

    for (Long64_t entry=event_start; entry<event_stop; entry++)     // loop through events 
    //for (Long64_t entry=event_start; entry<10000; entry++)
    {
        T->GetEntry(entry);
        int size = mass->size();
        for(int j = 0; j < size; j++){  // loop through pi0 candidates in the event
            if(!(abs(zVertex) <= 20)) {
                continue;       // vertex cut (since there is only one vertex per event you can move this to outside the cluster loop 
            }
//	    if(clusNTow->at(j) > 1) {
//		continue;
//	    }
            if(chi->at(j) > 4 || asym->at(j) > 0.7){
                continue;       // chi2 and asymmetry cut 
            }
            // cout << zVertex << ", " << chi->at(j) << ", " << asym->at(j) << endl;

            int phiBin1 = floor((8*(phiTow1->at(j)+.5)));       // going from tower id to mass bin # here 
            int phiBin2 = floor((8*(phiTow2->at(j)+.5)));
            float etaBin1 = floor((8*(etaTow1->at(j)+.5)));
            float etaBin2 = floor((8*(etaTow2->at(j)+.5)));
            float pi0Mass = mass->at(j);

            h_etaTow->Fill(etaTow1->at(j));
            h_etaTow->Fill(etaTow2->at(j));
            h_phiTow->Fill(phiTow1->at(j));
            h_phiTow->Fill(phiTow2->at(j));
            
            if(phiBin1 < 0 || phiBin1 >= 16 || phiBin2 < 0 || phiBin2 >= 16){
                continue;       // I dont know why this would ever happen
            } 

            // cout << phiBin1 << ", " << phiBin2 << endl;
            h_pi0[phiBin1]->Fill(etaBin1, pi0Mass); // for each cluster in the diphoton pair, fill the 2D histogram corresponding to that phi with the eta and mass
            h_pi0[phiBin2]->Fill(etaBin2, pi0Mass);

            h_position->Fill(etaTow1->at(j), phiTow1->at(j));
            h_position->Fill(etaTow2->at(j), phiTow2->at(j));
        }// end pi0 candidate loop
    }//end event loop
    

    fout->Write();      // write histograms to outfile 
    fout->Close();

}
