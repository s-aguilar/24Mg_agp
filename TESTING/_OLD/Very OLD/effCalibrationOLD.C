
#include <iostream>
using std::cout;
using std::endl;

#include <vector>
using std::vector;

#include <list>

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

void calibration(const char* fileName,const char* detector,int low, int high, int loop){

	TFile *fyield = new TFile(fileName);
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
	TF1 *ffit = new TF1("ffit",func,low,high,7);
	ffit->SetParNames("norm","mean","sigma","Nbg","a0","a1","a2");
	ffit->SetLineColor(kGreen);
	ffit->SetNpx(1e5);

	TF1 *fback = new TF1("fback",bgfunc,low,high,4);
	fback->SetParNames("Nbg","a0","a1","a2");
	fback->SetLineColor(kCyan);
	fback->SetNpx(1e5);


	// Get histograms from root file
	TH1D *hyield = static_cast<TH1D*>(fyield->Get(detector));
	HBACK = static_cast<TH1D*>(fbackground->Get(detector));
	HBACK->Scale(scale0);
	HBACK->SetDirectory(0);
	gStyle->SetOptFit(1111);


	// Create the canvas
	TCanvas *c0 = new TCanvas("c0","c0",600,800);
	c0->Divide(1,2);
	c0->Update();
	c0->cd(1);


	// Draw
	hyield->Draw();			// spectrum
	HBACK->Draw("SAME");	// background
	HBACK->SetLineColor(kViolet);
	hyield->GetXaxis()->SetRangeUser(500,2500);


	// Initial guess of parameters
	double centroid = (low+high)/2;
	if(loop==0) ffit->SetParameters(450,centroid,25,1,0,0.9,0);	  // 1173 peak
	else if(loop==1) ffit->SetParameters(450,centroid,25,1,0,0.9,0); // 1332 peak
	else if(loop==2) ffit->SetParameters(500,centroid,10,1,0,0.9,0); // 661 peak


    // Set parameter limits (parameter[i],min,max) and fix parameters
    ffit->SetParLimits(0,0,1e6);		// peak amplitude
	ffit->SetParLimits(1,low,high);		// Peak centroid
    ffit->SetParLimits(2,2,55);			// Peak Width
	// ffit->FixParameter(3,1);			// Nbg (scales background)
	ffit->SetParLimits(4,-5,5);			// y-intercept (linear background)
    ffit->FixParameter(6,0);			// Makes it a linear background (0)*x^2


	hyield->Fit("ffit","QR"); // R =  Use the range specified in the function \
                                range										  \
								Q = suppress fit print statements

	ffit->GetXaxis()->SetRangeUser(500,6000);


	// Reuse background parameters from fit for background fit
	for(int i =0; i<4; i++){
		fback->SetParameter(i,ffit->GetParameter(i+3));
	}


	// Draw fits
	ffit->Draw("SAME");
	fback->Draw("LSAME");


	// Background subtracted spectrum
	c0->cd(2);
	TH1D *ysubtracted = new TH1D();
	ysubtracted = (TH1D*)hyield->Clone("ysubtracted");

    // Recalculate errors manually
	for(int i = 0; i<=ysubtracted->GetNbinsX(); i++){
		double yval = ysubtracted->GetBinContent(i);
		double yval2 = fback->Eval(ysubtracted->GetBinCenter(i));
		double yerr = ysubtracted->GetBinError(i);
		double yerr2 = HBACK->GetBinError(i)*fback->GetParameter(0);

		ysubtracted->SetBinContent(i,yval-yval2);
		ysubtracted->SetBinError(i,TMath::Sqrt(yerr*yerr+yerr2*yerr2));
	}

	// Rebinning
	const int nbins = ysubtracted->GetXaxis()->GetNbins();
	double new_bins[nbins+1];
	for(int i=0; i <= nbins; i++){
		new_bins[i] = ysubtracted->GetBinLowEdge(i+1);
	}
	ysubtracted->SetBins(nbins, new_bins);

/*
	double par[7];
	ffit->GetParameters(par);
    ysubtracted->SetLineColor(kRed);
	ysubtracted->GetXaxis()->SetRangeUser(par[1] - 3*par[2] , par[1] + 3*par[2]);

	// ysubtracted->Fit("ffit","0");	// "0" Don't draw the previous fit
	ysubtracted->Draw();


	TF1 *gauss = new TF1("gauss",func,low,high,7);
	gauss->SetParNames("norm","mean","sigma","gar","gar0","gar1","gar2");

	if(loop==0) gauss->SetParameters(450,1075,25,1,0,0.9,0);		// 1173 peak
	else if(loop==1) gauss->SetParameters(450,1300,25,1,0,0.9,0);	// 1332 peak
	else if(loop==2) gauss->SetParameters(500,650,10,1,0,0.9,0);	// 661 peak

    gauss->SetParLimits(0,0,1e6);
    gauss->SetParLimits(1,low,high);
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
*/
	//
	// // double chi2NDF = ffit->GetChisquare()/ffit->GetNDF();
	// double A_err = ffit->GetParError(3);
	// double sig_err = ffit->GetParError(5);
	// double area = par[3]*TMath::Sqrt(2*TMath::Pi()*par[5]*par[5]);
	// double area_err = TMath::Sqrt(2*TMath::Pi()*((par[5]*A_err)*(par[5]*A_err)+(par[3]*sig_err)*(par[3]*sig_err)));
	//
	// double yield = area/charge;
	// double yield_err = area_err/charge;
	// // cout << chi2NDF << endl;
	// int goodFit;
	//
	// // if (chi2NDF <= 2) goodFit = 0;
	// // else goodFit = 1;

	ysubtracted->GetXaxis()->SetRangeUser(low,high);
	// ysubtracted->Fit("ffit","0Q");	// "0" Don't draw the previous fit
	// ysubtracted->Fit("fback","0Q");	// "0" Don't draw the previous fit
	delete ffit;
	delete fback;
	ysubtracted->Draw();

	string detNum = detector;

	if (loop==0) c0->SaveAs(Form("calPlots/60Co_1173peak/det_%s_Fit.pdf",detNum.c_str()));
	else if (loop==1) c0->SaveAs(Form("calPlots/60Co_1332peak/det_%s_Fit.pdf",detNum.c_str()));
	else if (loop==2) c0->SaveAs(Form("calPlots/60Co_1332peak/det_%s_Fit.pdf",detNum.c_str()));

	// cout << "Peak position is: " << mean << " +/- " << sigma << " Sum is: " << sumOfPeak << " +/- " << sumError << endl;
	// cout << "Peak position is: " << mean << " +/- " << sigma << endl;
	c0->Clear();
	fyield->Close();
	fbackground->Close();
	delete c0;

	gROOT->Reset();
}
/*===========================END OF FITTING==================================*/



/*============================START OF MAIN==================================*/
int effCalibration(){

	const char *fileLoc[] = {"60Co/24Mg_run/run0419.root","60Co/24Mg_run/run0419.root","137Cs/24Mg_run/run0420.root"};
	const char *detect[] = {"h0-0","h0-1","h0-2","h0-3","h0-4","h0-5","h0-6","h0-7","h1-0","h1-1","h1-2","h1-3","h1-4"};

	int *low, *high;

	for (int i=0;i<=2;i++){

		const char *file = fileLoc[i];

		int low1173[] = {1045,1130,1145,1125,1090,1150,1130,1140,1120,1140,1150,1120,1130};
		int high1173[] = {1095,1190,1200,1185,1150,1220,1200,1210,1190,1200,1210,1190,1190};
		int low1332[] = {1180,1280,1295,1270,1240,1305,1280,1290,1270,1290,1300,1270,1280};
		int high1332[] = {1240,1350,1355,1335,1305,1375,1360,1360,1350,1360,1380,1350,1350};
		int low661[] = {595,630,645,630,615,645,640,640,635,640,640,620,630};
		int high661[] = {625,675,685,672,675,685,680,680,675,675,680,670,675};

		// 1173 keV gamma
		if (strcmp(file,"60Co/24Mg_run/run0419.root") == 0 && i == 0){
            low = low1173;
			high = high1173;
		}

		// 1332 keV gamma
		else if (strcmp(file,"60Co/24Mg_run/run0419.root") == 0 && i == 1){
			low = low1332;
			high = high1332;
		}

		// 661 keV gamma
		else if (strcmp(file,"137Cs/24Mg_run/run0420.root")==0){
			low = low661;
			high = high661;
		}

		for (int j=0;j<=12;j++){
			const char *det = detect[j];
			calibration(file,det,low[j],high[j],i);
		}
		break;
	}

	return 0;
}
