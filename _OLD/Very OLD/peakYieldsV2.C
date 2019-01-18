#include <iostream>
using std::cout;
using std::endl;

#include <string>
using std::string;


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

double bgfunc(double *x, double *par){

	return par[0]*(HBACK->Interpolate(par[1]+par[2]*x[0]+par[3]*x[0]*x[0]));
}


double func(double *x, double *par){

	double norm = par[0];
	double mean = par[1];
	double sigma = par[2];
	double result = norm*TMath::Gaus(x[0],mean,sigma,kTRUE);
	result += par[3]*(HBACK->Interpolate(par[4]+par[5]*x[0]+par[6]*x[0]*x[0]));

	return result;
}


void peakFitter(const char* fileName,const char* detector,int low, int high){

	TFile *fyield = new TFile(fileName);
	TFile *fbackground = new TFile("background/run0422.root");

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
	HBACK->Scale(scale0);
	HBACK->SetDirectory(0);

	gStyle->SetOptFit(0111);

	TCanvas *c0 = new TCanvas("c0","c0",600,800);
	c0->Divide(1,2);
	c0->Update();
	c0->cd(1);

	hyield->Draw();
	HBACK->Draw("SAME");
	HBACK->SetLineColor(kViolet);
	HBACK->GetXaxis()->SetRangeUser(500,1400);
	hyield->GetXaxis()->SetRangeUser(500,1400);

	TF1 *ffit = new TF1("ffit",func,low,high,7);
	ffit->SetParNames("norm","mean","sigma","Nbg","a0","a1","a2");
	ffit->SetParameters(500,(low+high)/2,10,1,1,1,1);
	// ffit->SetParLimits(0,0,1e6);			// Amplitude
	// ffit->SetParLimits(1,low,high);			// Mean
	// ffit->SetParLimits(2,2.5,25);			// Sigma
	ffit->FixParameter(3,1);				// Nbg scales background
    ffit->FixParameter(6,0);				// Makes it a linear background
	ffit->SetLineColor(kGreen);

	TF1 *fback = new TF1("fback",bgfunc,low,high,6);
	fback->SetParNames("Nbg","a0","a1","a2");
	fback->FixParameter(0,1);
	fback->SetLineColor(kRed);

	ffit->SetNpx(1e5);
    fback->SetNpx(1e5);

	for(int i =0; i<4; i++){
		fback->SetParameter(i,ffit->GetParameter(i+3));
	}

	// hyield->Fit("ffit","R");
	// hyield->Fit("fback");
	// ffit->Draw("SAME");
	// fback->Draw("SAME"); // L



	c0->cd(2);
	TH1D *ysubtracted = new TH1D();
	ysubtracted = (TH1D*)hyield->Clone("ysubtracted");

	// Calculate errors
	// double peakArea;
	// double peakErr;
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

	double par[7];
	ffit->GetParameters(par);
    ysubtracted->SetLineColor(kRed);
	ysubtracted->GetXaxis()->SetRangeUser(par[1] - 3*par[2] , par[1] + 3*par[2]);

	ysubtracted->Fit("ffit","0");	// "0" Don't draw the previous fit
	ysubtracted->Draw();


	string runNum = fileName;
	// cout << "________________________________________________________" << endl;
	// cout << runNum << endl;
	runNum = runNum.substr(4,3);
	// runNum = runNum.substr(73,3);
	// cout << runNum << endl;
	string detNum = detector;

	c0->SaveAs(Form("peakAreas/run0%s/det_%s_Fit.pdf",runNum.c_str(),detNum.c_str()));


	// ofstream myfile;
	// myfile.open ("peakArea.csv",std::ios::app);
	// myfile << Form("run0%s",runNum.c_str()) <<","<< Form("det_%s",detNum.c_str())<<","<< yield <<","<< yield_err << "\n";
	// myfile.close();


	c0->Clear();
	fyield->Close();
	delete c0;
	delete ffit;
	gROOT->Reset();

}

/*===========================================================================*/




/*==============================MAIN=========================================*/

int peakYieldsV2(){

	// 159 to 410
	for(int i=159;i<160;i++){

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

		const char *files = Form("run0%d.root",i); //"/afs/crc.nd.edu/group/nsl/activetarget/data/24Mg_alpha_gamma/spectra/run0%d.root"

		double p1;
		double p2;
		// int low;
		// int high;

		// const char *files = Form("run0%d.root",i);

		string runNum = files;
		runNum = runNum.substr(4,3);
		gSystem->Exec(Form("mkdir peakAreas/run0%d",i));

		for(int j=0;j<8;j++){
			// const char *files = Form("run0%d.root",i);
			const char *files = Form("run0%d.root",i); //"/afs/crc.nd.edu/group/nsl/activetarget/data/24Mg_alpha_gamma/spectra/run0%d.root"
			const char *detect = Form("h0-%d",j); //   run0158.root  -  run0409.root
			// cout <<files << "  "<< detect <<endl;
			if(j==0){
				p2 = 927;
			}
			else if(j==1){
				p2 = 1002;
			}
			else if(j==2){
				p2 = 1014;
			}
			else if(j==3){
				p2 = 999;
			}
			else if(j==4){
				p2 = 973;
			}
			else if(j==5){
				p2 = 1022;
			}
			else if(j==6){
				p2 = 1008;
			}
			else if(j==7){
				p2 = 1013;
			}
			peakFitter(files,detect,p2-35,p2+35);

		}
		for(int k=0;k<5;k++){
			// const char *files = Form("run0%d.root",i);
			const char *files = Form("run0%d.root",i); //"/afs/crc.nd.edu/group/nsl/activetarget/data/24Mg_alpha_gamma/spectra/run0%d.root"
			const char*detect = Form("h1-%d",k);
			// cout <<files << "  "<< detect <<endl;
			if(k==0){
				p2 = 999;
			}
			else if(k==1){
				p2 = 1010;
			}
			else if(k==2){
				p2 = 1019;
			}
			else if(k==3){
				p2 = 995;
			}
			else if(k==4){
				p2 = 1000;
			}
			peakFitter(files,detect,p2-35,p2+35);

		}
	}

	return 0;
}
