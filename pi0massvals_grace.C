#include <complex>
#include <numeric>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
void pi0massvals_grace(std::string infile, std::string  outfile)
{
    //std::string infile = "/sphenix/user/ggarmire/DDPDC_run24/DDPDC_histograms/46xxx.root";
    //std::string outfile = "output46xxx.root";
    // set binning values 
    int nfit =0;
    int etabins = 768; 
    double etamin = -0.5;
    double etamax = 767.5;
    int phibins = 16;
    double phimin = -.5; 
    double phimax = 15.5;
    int pi0bins = 100; 
    double pi0min = 0;
    double pi0max = 0.7;
    double fitStart = 0.1;
    double fitEnd = 0.4;//<.5 eta meson
    double lowerSignalBoundEstimate = 0.1;
    double upperSignalBoundEstimate = 0.2;
    double sigmaEstimate = 0.025; // sigmaEstimate value

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(5000);
    
    // open infile 
    TFile *f = TFile::Open(infile.c_str());
    // TFile *f = new TFile(infile.c_str(), "read");
    // std::cout <<  << std::cout;
    // read in histograms 
    TH2F *h_pi0candidates[phibins];
    for (int iphihist = 0; iphihist < 16; iphihist++){
        // Form("h_pi0candidates%i", iphihist)
        // h_pi0candidates[iphihist] = (TH2F*)f->Get("h_pi0candidates0");
        f->GetObject(Form("h_pi0%i", iphihist), h_pi0candidates[iphihist]);
        //h_pi0_block[iphihist] = new TH2F(Form("h_pi0_byblock%i", iphihist), Form("pi0 masses on block for phi %i", iphihist), 16, -0.5, 15.5, pi0bins, pi0min, pi0max); 
    }
    //TCanvas c1;
    //h_pi0candidates[0]->Draw();

    // make outfile/ out histogram
    TFile *fout = TFile::Open(outfile.c_str(), "RECREATE");

    TH2F *h_pi0mass_bybin = new TH2F("h_pi0mass_bybin", "Pi0 Mass from fit for each bin", etabins, etamin, etamax, phibins, phimin, phimax);
    // one histogram for mass, one for peak width (sigma of gaussian)
    
    TH2F *h_diphotons_bybin = new TH2F("h_diphotons_bybin", "# of diphoton pairs for each bin", etabins, etamin, etamax, phibins, phimin, phimax);
    TH2F *h_pi0s_bybin = new TH2F("h_pi0s_bybin", "# of pi0s for each bin", etabins, etamin, etamax, phibins, phimin, phimax);
    TH2F *h_pi0mass_error = new TH2F("h_pi0mass_error", "Pi0 Mass from error for each bin", etabins, etamin, etamax, phibins, phimin, phimax);
    TH2F *h_pi0mass_negative = new TH2F("h_pi0mass_negative", "bins with negative pi0 mass", etabins, etamin, etamax, phibins, phimin, phimax);
    TH2F *h_pi0mass_bad = new TH2F("h_pi0mass_bad", "bins with pi0 mass greater than .3 or under .1", etabins, etamin, etamax, phibins, phimin, phimax);
    TH1F *h_pi0mass_distribution = new TH1F("h_pi0mass_distribution", "distribution of pi0 mass values from fits", 2500, -2.5, 2.5);
    TH2F *h_blockcount = new TH2F("h_pi0mass_blockcount", "should fill 1 time per block! ", etabins, etamin, etamax, phibins, phimin, phimax);
    TH1F *h_badfitcounts_eta = new TH1F("h_failedfitcounts_eta", "# of bad fits in each eta bin; eta bin; count", etabins, etamin, etamax);
    TH1F *h_badfitcounts_eta_conflated = new TH1F("h_failedfitcounts_eta_conflated", "# of bad fits across block in eta; eta bin; count", 16, etamin, phimax);
    TH1F *h_badfitcounts_phi = new TH1F("h_failedfitcounts_phi", "# of bad fits in each phi bin; eta phi; count", phibins, phimin, phimax);
    TH1F *h_pi0candidates_eta = new TH1F("h_pi0candidates_eta", "# of diphoton candidates in all fits in each eta bin; eta bin; pi0 candidates", etabins, etamin, etamax);
    TH1F *h_pi0candidates_phi = new TH1F("h_pi0candidates_phi", "# of diphoton candidates in all fits in each phi bin; phi bin; pi0 candidates", phibins, phimin, phimax);

    TH1F *h_meanmass_phi = new TH1F("h_meanmass_phi", "mean all mass values for each phi bin; phi bin; mass", phibins, phimin, phimax);


    TH1F *h_pi0mass_dist_byphibin[16];

    for (int i = 0; i < 16; i++) {      // just sets up the 16 2D hisotgrams (one for each phi bin)                    
        h_pi0mass_dist_byphibin[i] = new TH1F(Form("h_pi0mass_dist_byphibin%i", i), Form("pi0 mass ditribution in phi bin %i", i), 150, 0,0.3);
        h_pi0mass_dist_byphibin[i]->GetXaxis()->SetTitle("Pi0 Mass (GeV)");                                                           
        h_pi0mass_dist_byphibin[i]->GetYaxis()->SetTitle("");                                                           
    }        

    // loop thru phi histograms 
    int bin1, bin2; 

    for (int iphi = 0; iphi < 16; iphi++){
        cout << "starting Phi " <<iphi<< endl;
    // loop thru eta bins 
        //for (int ieta = 260; ieta < 100; ieta++){
            //if (ieta%10==0) cout << " on phi bin "<<iphi << " eta bin " << ieta << endl;
        for (int ieta = 0; ieta < etabins; ieta++){
    // build new histogram for the eta bin
            /*TH1F *h_pi0tofit = new TH1F("h_pi0tofit", "pi0mass", pi0bins, pi0min, pi0max);
            for (int ipi0bin = 0; ipi0bin < pi0bins; ipi0bin++){
                double val = h_pi0candidates[iphi]->GetBinContent(ieta, ipi0bin);
                h_pi0tofit->Fill(ieta, val);
            }*/
            TH1D *h_pi0tofit = h_pi0candidates[iphi]->ProjectionY("h_pi0tofit",ieta, ieta);
            /*for (int ipi0 = 0; ipi0 < pi0bins; ipi0++){
                int ietablock = ieta % 16;
                h_pi0_block[iphi]->Fill(ietablock, ipi0, h_pi0tofit->GetBinContent(ipi0));
            }*/

    //fit the bin with justins fitting - fit h_pi0tofit, put the value in h_pi0mass_bybin
            if (h_pi0tofit->Integral()==0) continue;
            nfit ++;
            TF1 *totalFit = new TF1("totalFit", "gaus(0) + pol2(3)", fitStart, fitEnd);
            
            
            
            double maxBinContent = h_pi0tofit->GetBinContent(h_pi0tofit->GetXaxis()->FindBin(0.15));
           // double maxBinCenter = h_pi0tofit->GetXaxis()->GetBinCenter(maxBin);                                                       
            totalFit->SetParameter(0, maxBinContent);                                                                   
            totalFit->SetParameter(1, 0.15); //.15                                                               
            totalFit->SetParameter(2, sigmaEstimate);
            TFitResultPtr fitResult = h_pi0tofit->Fit("totalFit", "QSR+");
            //return fitResult;
            //if (fitResult == nullptr) continue;
            double fitMean = totalFit->GetParameter(1);
            double fitsig = totalFit->GetParameter(2);
            //double fitMeanError = totalFit->GetParError(1);
            //if (fitMean > 0 && fitMean < 0.7) i
            bin1 = h_pi0tofit->FindBin(fitMean-2*fitsig);
            bin2 = h_pi0tofit->FindBin(fitMean+2*fitsig);
            if (fitMean > 0.12 && fitMean < 0.18) {
                double pi0yield = h_pi0tofit->Integral(bin1, bin2);
                h_pi0s_bybin->Fill(ieta, iphi, pi0yield);
            }

            h_pi0mass_bybin->Fill(ieta, iphi, fitMean);
            h_diphotons_bybin->Fill(ieta, iphi, h_pi0tofit->GetEntries());
            h_pi0mass_dist_byphibin[iphi]->Fill(fitMean);
            if (fitMean < 0) h_pi0mass_negative->Fill(ieta, iphi);
            h_pi0candidates_eta->Fill(ieta, h_pi0tofit->Integral());
            h_pi0candidates_phi->Fill(iphi, h_pi0tofit->Integral());
            //if (fitMeanError)h_pi0mass_error->Fill(ieta, iphi, fitMeanError);
            if (fitMean < 0.1 || fitMean > 0.3){ 
                h_pi0mass_bad->Fill(ieta, iphi);
                h_badfitcounts_eta->Fill(ieta);
                h_badfitcounts_eta_conflated->Fill(ieta%16);
                h_badfitcounts_phi->Fill(iphi);
            }
            h_pi0mass_distribution->Fill(fitMean);
         

    //fill the pi0mass hist
    

    //end loops 
        }
    }
    
    for (int iphi = 0; iphi < 16; iphi++) {
        h_meanmass_phi->Fill(iphi, h_pi0mass_dist_byphibin[iphi]->GetMean());
        h_meanmass_phi->SetBinError(h_pi0mass_dist_byphibin[iphi]->FindBin(iphi), (1/sqrt(768))*h_pi0mass_dist_byphibin[iphi]->GetStdDev());
        //cout << "std for bin " << iphi << ": " << h_pi0mass_dist_byphibin[iphi]->GetStdDev()  << endl;
        //cout << iphi << endl;
       // cout << "overflow: " << h_pi0mass_dist_byphibin[iphi]->GetBinContent(0) << "+" << h_pi0mass_dist_byphibin[iphi]->GetBinContent(2001) << endl; 
    }

    //h_meanmass_phi->SetBinError(15, h_pi0mass_dist_byphibin[15]->GetStdDev());

    cout <<nfit << endl;
    //cout <<h_pi0s_bybin->Integral() << endl;

    //close files 
    f->Close();
    fout->Write();
    fout->Close();
}

