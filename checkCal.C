#include <iostream>
using std::cout;
using std::endl;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "TROOT.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TError.h"

// Fitting routines
#include "calibration/fitFunctions.h"

// Gain match routine
#include "calibration/gainMatch.h"

// CSV I/O routine
#include "calibration/csvIO.h"

TH1D *HBACK;



void checkCal(){

    double low1173 = 1030;
    double high1173 = 1110;

    const char *fileLoc = "calibration/60Co/24Mg_run/run0419.root";
	const char *detect = "h0-0";
	const char *fileBackground = "calibration/background/run0422.root";

    TFile *fyield = new TFile(fileLoc);
    TFile *fbackground = new TFile(fileBackground);

    // Get runtime info from root file
    TVectorD *run_t0 = static_cast<TVectorD*>(fyield->Get("lastTimeStampSeconds-0"));
    TVectorD *back_t0 = static_cast<TVectorD*>(fbackground->Get("lastTimeStampSeconds-0"));


    // Run time and background time for each board
    double runTime0 = (*run_t0)[0];
    double backTime0 = (*back_t0)[0];


    // Scale background spectra to ratio of run time
    double scale0 = runTime0/backTime0;


    // Initial guess of centroid position
    // double centroid = (low+high)/2;


    // Get histograms from root file, prepare bg subtracted histogram
    TH1D *hyield = static_cast<TH1D*>(fyield->Get(detect));
    TH1D *ysubtracted = static_cast<TH1D*>(hyield->Clone("ysubtracted"));
    HBACK = static_cast<TH1D*>(fbackground->Get(detect));

    HBACK->Scale(scale0);		// Scale the background
    HBACK->SetDirectory(0);
    gStyle->SetOptFit(1111);


    // Create the canvas
    TCanvas *c0 = new TCanvas("c0","c0",600,800);
    c0->Divide(1,2);
    c0->Update();
    c0->cd(1);


    ////////////////////////
    ///    GAIN MATCH    ///
    ////////////////////////


    // New gain matched BG histogram
    TH1D *h2 = new TH1D("h2","h2",8192,0,8192);
    gain_match(1322,1787,1323,1789,h2,HBACK);


    TH1D *h3 = new TH1D("h3","h3",8192,0,8192);

    // Draw results
    hyield->Draw();
    HBACK->Draw("SAME");
    HBACK->SetLineColor(kGreen);
    h2->SetLineColor(kRed);
    h2->Draw("SAME");
    // gPad->SetLogy();


    hyield->GetXaxis()->SetRangeUser(800,1600);


    // Background subtracted spectrum
    c0->cd(2);


    // Recalculate errors manually
    for(int i = 0; i<=ysubtracted->GetNbinsX(); i++){
        double yval = ysubtracted->GetBinContent(i);
        double yval2 = h2->GetBinContent(i);	//fback->Eval(ysubtracted->GetBinCenter(i));
        double yerr = ysubtracted->GetBinError(i);
        double yerr2 = h2->GetBinError(i);	// Takes the error from the bg spectrum

        ysubtracted->SetBinContent(i,yval-yval2);
        ysubtracted->SetBinError(i,TMath::Sqrt(yerr*yerr+yerr2*yerr2));
    }

    // Rebinning, currently does nothing (keeps bins same)
    const int nbins = ysubtracted->GetXaxis()->GetNbins();
    double new_bins[nbins+1];
    for(int i=0; i <= nbins; i++){
        new_bins[i] = ysubtracted->GetBinLowEdge(i+1);
    }
    ysubtracted->SetBins(nbins, new_bins);

    h3 = calibrate2(1.12507,-32.3309,h3,ysubtracted);
    h3->Draw();
    h3->GetXaxis()->SetRangeUser(1000,1400);
    // gPad->SetLogy();
    h3->SetStats(kFALSE);

    TLine *line1 = new TLine(1173,0,1173,600);
    TLine *line2 = new TLine(1332,0,1332,600);
    line1->Draw("SAME");
    line2->Draw("SAME");

    // c0->Clear();
    // fyield->Close();
    // fbackground->Close();
    //
    // delete fyield;
    // delete fbackground;
    //
    // delete c0;
    //
    // delete HBACK;
    // delete h2; // This breaks code idk why, should be deleted though

    gROOT->Reset();
}
