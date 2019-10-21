#include <iostream>
using std::cout;
using std::endl;

#include <vector>
using std::vector;

#include "TROOT.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"

// #include <iostream>
// using std::cout;
// using std::endl;


//vector gloabl for a background histogram
TH1D *HBACK;



///////////////////////
// Fitting Functions //
///////////////////////

double bgfunc(double *x, double *par){
    // Quadratic background

	return par[0]*(HBACK->Interpolate(par[1]+par[2]*x[0]+par[3]*x[0]*x[0]));
}

double func(double *x, double *par){
    // Fit a Gaussian with a Quadratic bakground

	double norm = par[0];
	double mean = par[1];
	double sigma = par[2];

    // If norm=kTRUE (default is kFALSE) the result is divided \
    by sqrt(2*Pi)*sigma.
	double result = norm*TMath::Gaus(x[0],mean,sigma,kTRUE);


	result += par[3]*(HBACK->Interpolate(par[4]+par[5]*x[0]+par[6]*x[0]*x[0]));
	return result;
}



//////////
// MAIN //
//////////

int peakFitter(){
	TFile *fyield = new TFile("rootHistograms/run11284.root");      // 60Co
    // TFile *fyield = new TFile("rootHistograms/run11282.root");   // 137Cs
	TFile *fbackground = new TFile("rootHistograms/run10427.root");

	TF1 *ffit = new TF1("ffit",func,0,50000,7);
	ffit->SetParNames("norm","mean","sigma","Nbg","a0","a1","a2");
	TF1 *fback = new TF1("fback",bgfunc,0,50000,6);
	fback->SetParNames("Nbg","a0","a1","a2");

	TH1D *hyield = (TH1D*)fyield->Get("h45deg");
	HBACK = (TH1D*)fbackground->Get("h45deg");
	HBACK->SetDirectory(0);
	gStyle->SetOptFit(1111);

	TCanvas *c0 = new TCanvas("c0","c0");
	hyield->Draw();
	// hyield->GetXaxis()->SetRangeUser(1000,6000);     // 137Cs
	hyield->GetXaxis()->SetRangeUser(2600,6000);        // 60Co
    // ffit->SetParameters(5000,1500,30,0.2,0,1.0,0);   // 137Cs
    ffit->SetParameters(5000,2900,30,0.2,0,0.9,0);      // 60Co

    /* Setting parameter limits (parameter[i],min,max) */
    ffit->SetParLimits(0,0,1e6);
    ffit->SetParLimits(1,2700,3000);// 60Co
    // ffit->SetParLimits(1,1200,1600);// 137Cs
    ffit->SetParLimits(2,10,50);
	ffit->SetParLimits(4,-10,10);
    ffit->SetParLimits(5,0.85,1.15);

    // The specified value will be used in a fit operation
    ffit->FixParameter(6,0); // Makes it a linear background (0)*x^2

	hyield->Fit("ffit","R"); // R =  Use the range specified in the function \
                                range
	ffit->GetXaxis()->SetRangeUser(1000,6000);

	for(int i =0; i<4; i++){
		fback->SetParameter(i,ffit->GetParameter(i+3));
	}

    ffit->SetLineColor(kAzure);
    fback->SetLineColor(kBlack);
    fback->SetNpx(1e4); // Set the number of points used to draw the function
    ffit->SetNpx(1e5);
	ffit->Draw("SAME");
	fback->Draw("LSAME");

	TH1D *ysubtracted = new TH1D();
	ysubtracted = (TH1D*)hyield->Clone("ysubtracted");
	c0->Update();

    // Calculate errors
	for(int i = 0; i<=ysubtracted->GetNbinsX(); i++){
		double yval = ysubtracted->GetBinContent(i);
		double yval2 = fback->Eval(ysubtracted->GetBinCenter(i));
		double yerr = ysubtracted->GetBinError(i);
		double yerr2 = HBACK->GetBinError(i)*fback->GetParameter(0);

		ysubtracted->SetBinContent(i,yval-yval2);
		ysubtracted->SetBinError(i,TMath::Sqrt(yerr*yerr+yerr2*yerr2));
	}



	const int nbins = ysubtracted->GetXaxis()->GetNbins();
	double new_bins[nbins+1];

	for(int i=0; i <= nbins; i++){
	new_bins[i] = ysubtracted->GetBinLowEdge(i+1);
	}

	ysubtracted->SetBins(nbins, new_bins);

	TCanvas *c1 = new TCanvas("c1","c1");
	c1->cd();
	c1->Update();
	// HBACK->Draw();
    ysubtracted->SetLineColor(kRed);
	ysubtracted->Draw();
    ffit->SetLineColor(kAzure);
	ffit->Draw("LSAME");

	return 0;
}
