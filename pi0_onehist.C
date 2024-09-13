#include <complex>
#include <numeric>
#include <stdlib.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TLatex.h>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TList.h"
#include "TSystemFile.h"
#include <sys/stat.h>
void pi0_onehist(string infile, int phibin, int etabin)
{

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(5000);
// files and trees in
    TFile *f = TFile::Open(infile.c_str());     // open the infile

    //TFile *fout = TFile::Open(outfile.c_str(), "RECREATE");     // make the outfile
// read in the histogram
    TH2F *h_pi0_byphi;
    f->GetObject(Form("h_pi0%i", phibin), h_pi0_byphi);

// make 1D hist to fit 
    h_pi0_byphi->SetTitle(Form("pi0 Mass in Phi bin %i, Eta bin %i", phibin, etabin));
    TH1D *h_pi0_foreta = h_pi0_byphi->ProjectionY("", etabin, etabin);

// initial fit parameters 
    
    double fitStart = 0.1;
    double fitEnd = 0.4;
    double meanest = 0.15;

    double sigmaEstimate = 0.025; // sigmaEstimate value

                        
// Define totalFit
//
//
//
//
//
//
// Alex Fitting code copied here:
    TF1 *totalFit = new TF1("totalFit", "gaus(0) + pol2(3)", fitStart, fitEnd);

    double maxBinContent = h_pi0_foreta->GetBinContent(h_pi0_foreta->GetXaxis()->FindBin(meanest));
    totalFit->SetParameter(0, maxBinContent);                                                               
    totalFit->SetParameter(1, meanest); //.15                                                               
    totalFit->SetParameter(2, sigmaEstimate);
    TFitResultPtr fitResult = h_pi0_foreta->Fit("totalFit","SR+");
    //return fitResult;
    //if (fitResult == nullptr) cout << "fit didn't work! " << endl;


    auto c1 = new TCanvas("c1");
    h_pi0_foreta->SetMarkerStyle(20);
    h_pi0_foreta->SetMarkerSize(1.0);
    h_pi0_foreta->SetMarkerColor(kBlack);
    h_pi0_foreta->Draw("PE");

    double fitMean = totalFit->GetParameter(1);
    double fitMeanError = totalFit->GetParError(1);
    double fitSigma = totalFit->GetParameter(2);
    double fitSigmaError = totalFit->GetParError(2);
    double numEntries = h_pi0_foreta->GetEntries();
    
    TF1 *gaussFit = new TF1("gaussFit", "gaus", fitStart, .3);
    TF1 *polyFit = new TF1("polyFit", "pol2", 0, .3);
    gaussFit->SetParameter(0, totalFit->GetParameter(0));
    gaussFit->SetParameter(1, totalFit->GetParameter(1));
    gaussFit->SetParameter(2, totalFit->GetParameter(2));

    for (int i = 3; i < 6; i++) {
        polyFit->SetParameter(i - 3, totalFit->GetParameter(i));
        polyFit->SetParError(i - 3, totalFit->GetParError(i));
    }
    
    TLatex *t = new TLatex(0.4, 0.67, Form("Gaussian Mean: %f GeV", fitMean));
    t->SetNDC(kTRUE);
    t->Draw("same");

    gaussFit->SetLineColor(kOrange+2);
    gaussFit->SetLineStyle(2);
    gaussFit->Draw("SAME");
    polyFit->SetLineColor(kBlue);
    polyFit->SetLineWidth(3);
    polyFit->SetLineStyle(2);
    polyFit->Draw("SAME");
    totalFit->SetLineColor(kBlack);
    totalFit->SetLineStyle(2);
    totalFit->Draw("SAME");



    c1->SaveAs(Form("pi0fit_phi%i_eta%i.pdf", phibin, etabin));


}

