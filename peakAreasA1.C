#include <iostream>
using std::cout;
using std::endl;

#include <string>
using std::string;

#include <unistd.h>


#include "TROOT.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TStyle.h"


//vector global for a background histogram
TH1D *HBACK;


/*=============================FITTING=======================================*/

double background(double *x, double *par){
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


double func(double *x, double *par) {
   return background(x,par) + gauss(x,&par[3]);
}


void peakFitter(const char* fileName,const char* detector,int low, int high){

	TFile *fyield = new TFile(fileName);
	TFile *fbackground = new TFile("calibration/background/run0422.root");
	// TFile *fbackground = new TFile("background/run0422.root");
	TH1D *qcharge = static_cast<TH1D*>(fyield->Get("h1-7"));
	double charge = qcharge->GetEntries();
	gStyle->SetOptFit(1111);

	// Get runtime info from root file
	TVectorD *run_t0 = static_cast<TVectorD*>(fyield->Get("lastTimeStampSeconds-0"));
	TVectorD *back_t0 = static_cast<TVectorD*>(fbackground->Get("lastTimeStampSeconds-0"));

	// Run time and background time for board
	double runTime0 = (*run_t0)[0];
	double backTime0 = (*back_t0)[0];

	// Scale background spectra to ratio of run time
	double scale0 = runTime0/backTime0;

	TH1D *hyield = static_cast<TH1D*>(fyield->Get(detector));
	HBACK = static_cast<TH1D*>(fbackground->Get(detector));
	HBACK->SetLineColor(kBlack);
	HBACK->Scale(scale0);
	HBACK->SetDirectory(0);

	// Create the canvas
	TCanvas *c0 = new TCanvas("c0","c0",600,800);
	c0->Divide(1,2);
	c0->Update();
	c0->cd(1);

	// Draw experimental and intrinsic background spectrum
	hyield->Draw();
	hyield->GetXaxis()->SetRangeUser(500,1800);
	HBACK->Draw("SAME");

	c0->cd(2);

	// Subtract the intrinsic background spectrum
	TH1D *ysubtracted = new TH1D();
	ysubtracted = (TH1D*)hyield->Clone("ysubtracted");


	// Recalculate errors manually
	for(int i = 0; i<=ysubtracted->GetNbinsX(); i++){
		double yval = ysubtracted->GetBinContent(i);
		double yval2 = HBACK->GetBinContent(i);	//fback->Eval(ysubtracted->GetBinCenter(i));
		double yerr = ysubtracted->GetBinError(i);
		double yerr2 = HBACK->GetBinError(i);	// Takes the error from the bg spectrum

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


	TF1 *ffit = new TF1("ffit",func,low,high,6);
	ffit->SetParNames("a0","a1","a2","norm","mean","sigma");
	ffit->SetNpx(500);
	ffit->SetParameters(1,1,0,500,(low+high)/2,10);   //// give it good range

	// ffit->SetParLimits(0,-1e3,1e3)
	// ffit->SetParLimits(1,-15,15);
	ffit->FixParameter(2,0);
	ffit->SetParLimits(3,0,1e6);
    ffit->SetParLimits(4,low,high);	///// give it good range
    ffit->SetParLimits(5,2.5,55);

	ffit->SetLineColor(kRed);

	TF1 *fback = new TF1("fback",background,low,high,3);
	fback->SetParNames("a0","a1","a2");
	fback->SetLineColor(kCyan);
	fback->SetNpx(1e5);


	// Perform the fit
	TFitResultPtr r = ysubtracted->Fit("ffit","SQR");  // TFitResultPtr contains the TFitResult
	TMatrixDSym cov = r->GetCovarianceMatrix();   //  to access the covariance matrix

	// Draw the background subtracted histogram
	ysubtracted->GetXaxis()->SetRangeUser(low-20,high+20);
	ysubtracted->Draw();


	// Draw fits
	ffit->Draw("SAME");

	// Store fit parameters in an array
	double par[9];
	ffit->GetParameters(par);
	fback->SetParameters(par[0],par[1],par[2]);
	fback->Draw("SAME");


	double A = par[3];
	double A_err = ffit->GetParError(3);

	double yield = A/charge;
	// cout << yield << " " <<typeid(area).name() << " " <<typeid(charge).name()<<endl;
	double yield_err = A_err/charge;

	double chi2NDF = ffit->GetChisquare()/ffit->GetNDF();
	double goodFit;
	if (chi2NDF <= 2 && chi2NDF >=.5 ) goodFit = 0;
	else goodFit = 1;


	string runNum = fileName;
	cout << "________________________________________________________" << endl;
	runNum = runNum.substr(4,3);
	// runNum = runNum.substr(73,3);
	cout << "Run: "<<runNum << endl;
	cout << "Chi2: "<<chi2NDF << endl;
	cout << "________________________________________________________" << endl;
	string detNum = detector;

	c0->SaveAs(Form("peakAreasA1/run0%s/det_%s_Fit.png",runNum.c_str(),detNum.c_str()));

	ofstream myfile;
	myfile.open ("peakAreasA1.csv",std::ios::app);
	myfile<<Form("run0%s",runNum.c_str())<<","<< Form("det_%s",detNum.c_str())<<
				","<<yield<<","<<yield_err<<","<<goodFit<<"\n";
	myfile.close();



	c0->Clear();
	fyield->Close();
	fbackground->Close();
	delete c0;
	delete ffit;
	gROOT->Reset();

}

/*===========================================================================*/




/*==============================MAIN=========================================*/

int peakAreasA1(){

	// const char *path = "/afs/crc.nd.edu/user/s/saguilar/Group/24Mg_ap/";
	// chdir(path);
	gSystem->Exec(Form("mkdir peakAreasA1"));

	ofstream myfile;
	myfile.open ("peakAreasA1.csv",std::ios::out);
	myfile<<"Run"<<","<<"Detector"<<","<<"Yield"<<","<<"Yield err"<<","<<"Fit Status"<<"\n";
	myfile.close();

	// 159 to 410
	for(int i=159;i<162;i++){///////////////////////

		if(i==163) continue;
		else if(i==164) continue;
		else if(i==166) continue;
		else if(i==182) continue;
		else if(i==244) continue;
		else if(i==164) continue;
		else if(i==166) continue;
		else if(i==182) continue;
		else if(i>=244 && i<=255) continue;
		else if(i==276) continue;
		else if(i==277) continue;
		else if(i==285) continue;
		else if(i==289) continue;
		else if(i==290) continue;
		else if(i==294) continue;
		else if(i==255) continue;

		double a1;
		gSystem->Exec(Form("mkdir peakAreasA1/run0%d",i));



		for(int j=0;j<8;j++){
			const char *files = Form("run0%d.root",i);
			// const char *files = Form("/afs/crc.nd.edu/group/nsl/activetarget/data/24Mg_alpha_gamma/spectra/run0%d.root",i);
			const char *detect = Form("h0-%d",j);
			if(j==0){
				a1 = 1242;
			}
			else if(j==1){
				a1 = 1349;
			}
			else if(j==2){
				a1 = 1360;
			}
			else if(j==3){
				a1 = 1342;
			}
			else if(j==4){
				a1 = 1308;
			}
			else if(j==5){
				a1 = 1374;
			}
			else if(j==6){
				a1 = 1355;
			}
			else if(j==7){
				a1 = 1364;
			}
			peakFitter(files,detect,a1-60,a1+60);

		}
		for(int k=0;k<5;k++){
			const char *files = Form("run0%d.root",i);
			// const char *files = Form("/afs/crc.nd.edu/group/nsl/activetarget/data/24Mg_alpha_gamma/spectra/run0%d.root",i);
			const char*detect = Form("h1-%d",k);
			if(k==0){
				a1 = 1345;
			}
			else if(k==1){
				a1 = 1361;
			}
			else if(k==2){
				a1 = 1375;
			}
			else if(k==3){
				a1 = 1344;
			}
			else if(k==4){
				a1 = 1347;
			}
			peakFitter(files,detect,a1-60,a1+60);
		}
	}

	return 0;
}
