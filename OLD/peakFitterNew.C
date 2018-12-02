
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
#include "TVectorD.h" //
#include "TCanvas.h"
#include "TStyle.h"


//vector global for a background histogram
TH1D *HBACK;



/*=============================FITTING=======================================*/

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

/*===========================================================================*/




/*==============================MAIN=========================================*/

int calibration(){

	const char* fileName = "137Cs/24Mg_run/run0420.root";
	const char *detector = "h1-4";
	int low = 625;
	int high = 690;

	TFile *fyield = new TFile(fileName);      // 60Co
	TFile *fbackground = new TFile("background/run0422.root");    // Intrinsic \
	 																 background

	// Get runtime info from root file
	TVectorD *run_t0 = static_cast<TVectorD*>(fyield->Get("lastTimeStampSeconds-0"));
	TVectorD *back_t0 = static_cast<TVectorD*>(fbackground->Get("lastTimeStampSeconds-0"));

	// Run time and background time for each board
	double runTime0 = (*run_t0)[0];
	double backTime0 = (*back_t0)[0];

	// Scale background spectra to ratio of run time
	double scale0 = runTime0/backTime0;

	// Fit definitions
	TF1 *ffit = new TF1("ffit",func,low,high,6);
	TF1 *fback = new TF1("fback",bgfunc,low,high,3);
	ffit->SetParNames("norm","mean","sigma","a0","a1","a2");
	fback->SetParNames("a0","a1","a2");

	// Get histograms from root file
	TH1D *hyield = static_cast<TH1D*>(fyield->Get(detector));
	HBACK = static_cast<TH1D*>(fbackground->Get(detector));
	HBACK->Scale(scale0);
	HBACK->SetDirectory(0);
	gStyle->SetOptFit(1111);

	// Create the canvas and draw
	TCanvas *c0 = new TCanvas("c0","c0");
	HBACK->SetLineColor(kViolet);
	ffit->SetLineColor(kGreen);
	hyield->Draw();
	HBACK->Draw("SAME");

	// hyield->GetXaxis()->SetRangeUser(1000,6000);     // 137Cs
	hyield->GetXaxis()->SetRangeUser(500,2500);			// 60Co
    ffit->SetParameters(500,650,10,1,0,1.0,0);   // 137Cs
    // ffit->SetParameters(450,1075,25,1,0,0.9,0);			// 60Co
	// ffit->SetParameters(400,1300,25,1,0,1,0);			// Other 60Co peak

    // Setting parameter limits (parameter[i],min,max)
    ffit->SetParLimits(0,0,1e6);
    ffit->SetParLimits(1,500,800);	// 60Co
    ffit->SetParLimits(2,2,55);
	ffit->SetParLimits(4,-5,5);
    // ffit->SetParLimits(5,0,1.15);	//
	ffit->FixParameter(3,1);			// Nbg scales background
    ffit->FixParameter(6,0);			// Makes it a linear background (0)*x^2


	hyield->Fit("ffit","R"); // R =  Use the range specified in the function \
                                range
	ffit->GetXaxis()->SetRangeUser(500,6000); // 1000,6000 for 60Co

	for(int i =0; i<4; i++){
		fback->SetParameter(i,ffit->GetParameter(i+3));
	}

	ffit->SetNpx(1e5);
    fback->SetNpx(1e5);
	ffit->Draw("SAME");
	fback->Draw("LSAME"); //


	TH1D *ysubtracted = new TH1D();
	ysubtracted = (TH1D*)hyield->Clone("ysubtracted");
	c0->Update();
    c0->Flush();


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
	// c1->Update();

	double par[7];
	ffit->GetParameters(par);
    ysubtracted->SetLineColor(kRed);
	ysubtracted->GetXaxis()->SetRangeUser(par[1] - 3*par[2] , par[1] + 3*par[2]);

	// ysubtracted->Fit("ffit","0");	// "0" Don't draw the previous fit
	ysubtracted->Draw();


	TF1 *gauss = new TF1("gauss",func,low,high,7);
	gauss->SetParNames("norm","mean","sigma","gar","gar0","gar1","gar2");
	gauss->SetParameters(800,600,10,0,0,0,0);
    gauss->SetParLimits(0,0,1e6);
    gauss->SetParLimits(1,500,800);
    gauss->SetParLimits(2,2,45);
	gauss->FixParameter(3,1);
	gauss->FixParameter(4,0);
    gauss->FixParameter(5,0);
    gauss->FixParameter(6,0);
	gauss->SetLineColor(kGreen);
	// ysubtracted->SetLineColor(kBlue);
	ysubtracted->Fit("gauss","R");
	gauss->Draw("SAME");


	double ppar[7];
	gauss->GetParameters(ppar);
	double intg = gauss->Integral(ppar[1] - 3*ppar[2] , ppar[1] + 3*ppar[2],1.e-8);
	double sumError = gauss->IntegralError(ppar[1] - 3*ppar[2] , ppar[1] + 3*ppar[2]);
	double binw = ysubtracted->GetBinWidth(1);
	double sumOfPeak = intg/binw;
	double mean = ppar[1];
	double sigma = ppar[2];


	cout << "Peak position is: " << mean << " +/- " << sigma << " Sum is: " << sumOfPeak << " +/- " << sumError << endl;
	// cout << "Peak position is: " << mean << " +/- " << sigma << endl;

	return 0;

}
