
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
	/*
	COEFFICIENTS:
	===============
	constant = par[0]
	linear = par[1]
	quadratic = par[2]
	// */
	return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
}

double gauss(double *x, double *par){
	/*
	norm = par[0]
	mean = par[1]
	sigma = par[2]
	// */
	return par[0]*TMath::Gaus(x[0],par[1],par[2],kTRUE);
}

double func(double *x, double *par){
    // Fit a Gaussian with a Quadratic bakground

	return bgfunc(x,par) + gauss(x,&par[3]);
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
	TF1 *ffit = new TF1("ffit",func,low,high,6);
	ffit->SetParNames("a0","a1","a2","norm","mean","sigma");
	ffit->SetLineColor(kGreen);
	ffit->SetNpx(1e5);

	TF1 *fback = new TF1("fback",bgfunc,low,high,3);
	fback->SetParNames("a0","a1","a2");
	fback->SetLineColor(kCyan);
	fback->SetNpx(1e5);


	// Get histograms from root file
	TH1D *hyield = static_cast<TH1D*>(fyield->Get(detector));
	TH1D *ysubtracted = new TH1D();
	ysubtracted = (TH1D*)hyield->Clone("ysubtracted");

	HBACK = static_cast<TH1D*>(fbackground->Get(detector));
	HBACK->SetLineColor(kViolet);
	HBACK->Scale(scale0);	// Scale the background
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
	hyield->GetXaxis()->SetRangeUser(500,2500);
	ffit->GetXaxis()->SetRangeUser(500,2500);
	fback->GetXaxis()->SetRangeUser(500,2500);


	// Initial guess of parameters
	double centroid = (low+high)/2;
	ffit->SetParameters(0,0.9,0,500,centroid,10);


    // Set parameter limits (parameter[i],min,max) and fix parameters
	ffit->FixParameter(2,0);			// Makes it a linear background (0)*x^2
    ffit->SetParLimits(3,0,1e6);		// Peak Area
	ffit->SetParLimits(4,low,high);		// Peak Centroid
    ffit->SetParLimits(5,2,55);			// Peak Width



	hyield->Fit("ffit","QR"); // R =  Use the range specified in the function \
                                range										  \
								Q = suppress fit print statements


	// Draw fits
	ffit->Draw("SAME");
	fback->Draw("LSAME");


	// Store fit parameters in an array
	double par[6];
	ffit->GetParameters(par);
	fback->SetParameters(par[0],par[1],par[2]);


	// Background subtracted spectrum
	c0->cd(2);


    // Recalculate errors manually
	for(int i = 0; i<=ysubtracted->GetNbinsX(); i++){
		double yval = ysubtracted->GetBinContent(i);
		double yval2 = fback->Eval(ysubtracted->GetBinCenter(i));
		double yerr = ysubtracted->GetBinError(i);
		double yerr2 = HBACK->GetBinError(i);

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


	// Draw the background subtracted histogram
	ysubtracted->GetXaxis()->SetRangeUser(low,high);
	ysubtracted->Draw();


	// Calculate the peak area and error
	string detNum = detector;
	double area_err;
	double area = ysubtracted->IntegralAndError(par[4]-3*par[5],par[4]+3*par[5],area_err,"width");
	double mean = par[4];
	double sigma = par[5];

	// Save fits and calibration information
	if (loop==0) {
		c0->SaveAs(Form("calPlots/60Co_1173peak/det_%s_Fit.pdf",detNum.c_str()));

		ofstream myfile;
		myfile.open ("60Co_1173cal.csv",std::ios::app);
		myfile<< Form("det_%s",detNum.c_str())<<","<<mean<<","<<sigma<<","<<area
				<<","<<area_err<<","<<runTime0<<"\n";
		myfile.close();
	}
	else if (loop==1) {
		c0->SaveAs(Form("calPlots/60Co_1332peak/det_%s_Fit.pdf",detNum.c_str()));

		ofstream myfile;
		myfile.open ("60Co_1332cal.csv",std::ios::app);
		myfile<< Form("det_%s",detNum.c_str())<<","<<mean<<","<<sigma<<","<<area
				<<","<<area_err<<","<<runTime0<<"\n";
		myfile.close();
	}
	else if (loop==2) {
		c0->SaveAs(Form("calPlots/137Cs_661peak/det_%s_Fit.pdf",detNum.c_str()));

		ofstream myfile;
		myfile.open ("137Cs_661cal.csv",std::ios::app);
		myfile<< Form("det_%s",detNum.c_str())<<","<<mean<<","<<sigma<<","<<area
				<<","<<area_err<<","<<runTime0<<"\n";
		myfile.close();
	}

	c0->Clear();
	fyield->Close();
	fbackground->Close();
	delete ffit;
	delete fback;
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

		int low1173[] = {1040,1125,1140,1116,1090,1145,1130,1135,1120,1138,1140,1113,1123};
		int high1173[] = {1100,1200,1205,1183,1155,1220,1200,1210,1190,1205,1215,1190,1190};
		int low1332[] = {1180,1280,1290,1265,1235,1300,1280,1285,1270,1290,1295,1270,1280};
		int high1332[] = {1240,1350,1360,1343,1310,1380,1360,1365,1350,1360,1380,1350,1350};
		int low661[] = {590,630,645,634,615,647,638,638,632,638,638,619,630};
		int high661[] = {630,675,685,675,662,690,687,682,675,677,685,672,675};

		// 1173 keV gamma
		if (strcmp(file,"60Co/24Mg_run/run0419.root") == 0 && i == 0){
            low = low1173;
			high = high1173;

			ofstream myfile;
			myfile.open ("60Co_1173cal.csv",std::ios::app);
			myfile<<"Detector"<<","<<"Centroid"<<","<<"Width"<<","<<"Area"<<","
					<<"Area err"<<","<<"Runtime"<<"\n";
			myfile.close();
		}

		// 1332 keV gamma
		else if (strcmp(file,"60Co/24Mg_run/run0419.root") == 0 && i == 1){
			low = low1332;
			high = high1332;

			ofstream myfile;
			myfile.open ("60Co_1332cal.csv",std::ios::app);
			myfile<<"Detector"<<","<<"Centroid"<<","<<"Width"<<","<<"Area"<<","
					<<"Area err"<<","<<"Runtime"<<"\n";
			myfile.close();
		}

		// 661 keV gamma
		else if (strcmp(file,"137Cs/24Mg_run/run0420.root")==0 && i == 2){
			low = low661;
			high = high661;

			ofstream myfile;
			myfile.open ("137Cs_661cal.csv",std::ios::app);
			myfile<<"Detector"<<","<<"Centroid"<<","<<"Width"<<","<<"Area"<<","
					<<"Area err"<<","<<"Runtime"<<"\n";
			myfile.close();
		}

		for (int j=0;j<=12;j++){
			const char *det = detect[j];
			calibration(file,det,low[j],high[j],i);
		}

	}

	return 0;
}
