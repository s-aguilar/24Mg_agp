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

double double_gauss(double*x, double *par){
	/*
	norm1 = par[0]
	mean1 = par[1]
	sigma1 = par[2]
	norm2 = par[3]
	mean2 = par[4]
	sigma2 = par[5]
	// */
	double g1 = par[0]*TMath::Gaus(x[0],par[1],par[2],kTRUE);
	double g2 = par[3]*TMath::Gaus(x[0],par[4],par[5],kTRUE);
	return g1+g2;

}

double func(double *x, double *par) {

   return background(x,par) + double_gauss(x,&par[3]);
}


void peakFitter(const char* fileName,const char* detector,int low, int high){

	TFile *fyield = new TFile(fileName);
	TFile *fbackground = new TFile("background/run0422.root");
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
	HBACK->Scale(scale0);
	HBACK->SetDirectory(0);

	// Create the canvas
	TCanvas *c0 = new TCanvas("c0","c0",600,800);
	c0->Divide(1,2);
	c0->Update();
	c0->cd(1);

	// Draw experimental and intrinsic background spectrum
	hyield->Draw();
	hyield->GetXaxis()->SetRangeUser(500,1400);
	HBACK->Draw("SAME");
	HBACK->SetLineColor(kViolet);


	// Subtract the intrinsic background spectrum
	TH1D *ysubtracted = new TH1D();
	ysubtracted = (TH1D*)hyield->Clone("ysubtracted");
	c0->Update();
    c0->Flush();





	TF1 *ffit = new TF1("ffit",func,low,high,9);
	ffit->SetParNames("a0","a1","a2","norm1","mean1","sigma1","norm2","mean2","sigma2");
	ffit->SetNpx(500);
	ffit->SetParameters(1,1,1,500,low+30,10,500,low+55,10);   //// give it good range

	// ffit->SetParLimits(0,-1e3,1e3)
	// ffit->SetParLimits(1,-15,15);
    ffit->FixParameter(2,0);			// Makes it a linear background (0)*x^2
    ffit->SetParLimits(3,0,1e6);
    ffit->SetParLimits(4,low+20,low+40);	///// give it good range // first peak
    ffit->SetParLimits(5,2.5,12);
	ffit->SetParLimits(6,0,1e6);
    ffit->SetParLimits(7,low+25,high);	///// give it good range
    ffit->SetParLimits(8,2.5,15);

	ffit->SetLineColor(kRed);
	hyield->Fit("ffit","QR");


	// int fitStatus1 = hyield->Fit("ffit","QR");
	TF1 *back = new TF1("back",background,low,high,3);
	back->SetLineColor(kMagenta);
	back->SetNpx(500);

	TF1 *peak = new TF1("peak",double_gauss,low,high,6);
	peak->SetLineColor(kBlue);
	peak->SetNpx(500);

	double par[9];
	ffit->GetParameters(par);
	back->SetParameters(par[0],par[1],par[2]);
	peak->SetParameters(par[3],par[4],par[5],par[6],par[7],par[8]);

	back->Draw("SAME");
	peak->Draw("SAME");


	c0->cd(2);
	// c0->Flush();


	// Calculate errors for whole histogram
	for(int i = 0; i<=ysubtracted->GetNbinsX(); i++){
		double yval = ysubtracted->GetBinContent(i);
		double yval2 = back->Eval(ysubtracted->GetBinCenter(i));
		double yerr = ysubtracted->GetBinError(i);
		double yerr2 = hyield->GetBinError(i);

		ysubtracted->SetBinContent(i,yval-yval2);
		ysubtracted->SetBinError(i,TMath::Sqrt(yerr*yerr+yerr2*yerr2));
	}

	const int nbins = ysubtracted->GetXaxis()->GetNbins();
	double new_bins[nbins+1];

	for(int i=0; i <= nbins; i++){
	new_bins[i] = ysubtracted->GetBinLowEdge(i+1);
	}

	ysubtracted->SetBins(nbins, new_bins);
	ysubtracted->GetXaxis()->SetRangeUser(low,high);


/*
	int lower = par[4]-3*par[5];
	int upper = par[4]+3*par[5];
	int area = 0;
	double err = 0;
	double area_err = 0;

	// Calculate area and its error
	for(int i = lower; i<=upper; i++){
		area += ysubtracted->GetBinContent(i);
		err += (ysubtracted->GetBinError(i))*(ysubtracted->GetBinError(i));
	}
	area_err = TMath::Sqrt(err);
*/


	double chi2NDF = ffit->GetChisquare()/ffit->GetNDF();
	double A_err = ffit->GetParError(3);
	double sig_err = ffit->GetParError(5);
	double area = par[3]*TMath::Sqrt(2*TMath::Pi()*par[5]*par[5]);
	double area_err = TMath::Sqrt(2*TMath::Pi()*((par[5]*A_err)*(par[5]*A_err)+(par[3]*sig_err)*(par[3]*sig_err)));

	double yield = area/charge;
	double yield_err = area_err/charge;
	// cout << chi2NDF << endl;
	int goodFit;

	if (chi2NDF <= 2) goodFit = 0;
	else goodFit = 1;


	// ysubtracted->Fit("ffit","0Q");	// "0" Don't draw the previous fit
	ysubtracted->Draw();

	string runNum = fileName;
	// cout << "________________________________________________________" << endl;
	// cout << runNum << endl;
	// runNum = runNum.substr(4,3);
	runNum = runNum.substr(73,3);
	// cout << runNum << endl;
	string detNum = detector;

	c0->SaveAs(Form("peakAreasP1/run0%s/det_%s_Fit.pdf",runNum.c_str(),detNum.c_str()));

// /*
	ofstream myfile;
	myfile.open ("peakAreasP1.csv",std::ios::app);
	myfile<<Form("run0%s",runNum.c_str())<<","<< Form("det_%s",detNum.c_str())<<
				","<<yield<<","<<yield_err<<","<<goodFit<<"\n";
	myfile.close();

// */
	c0->Clear();
	fyield->Close();
	fbackground->Close();
	delete c0;
	delete ffit;
	gROOT->Reset();

}

/*===========================================================================*/




/*==============================MAIN=========================================*/

int peakAreasP1(){

	const char *path = "/afs/crc.nd.edu/user/s/saguilar/Group/24Mg_ap/";
	chdir(path);
	gSystem->Exec(Form("mkdir peakAreasP1"));

	ofstream myfile;
	myfile.open ("peakAreasP1.csv",std::ios::app);
	myfile<<"Run"<<","<<"Detector"<<","<<"Yield"<<","<<"Yield err"<<","<<"Fit Status"<<"\n";
	myfile.close();

	// 159 to 410
	for(int i=159;i<410;i++){

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

		double p1;
		gSystem->Exec(Form("mkdir peakAreasP1/run0%d",i));



		for(int j=0;j<8;j++){
			// const char *files = Form("run0%d.root",i);
			const char *files = Form("/afs/crc.nd.edu/group/nsl/activetarget/data/24Mg_alpha_gamma/spectra/run0%d.root",i);
			const char *detect = Form("h0-%d",j);
			// cout <<files << "  "<< detect <<endl;
			if(j==0){
				p1 = 773;
			}
			else if(j==1){
				p1 = 833;
			}
			else if(j==2){
				p1 = 845;
			}
			else if(j==3){
				p1 = 832;
			}
			else if(j==4){
				p1 = 810;
			}
			else if(j==5){
				p1 = 851;
			}
			else if(j==6){
				p1 = 840;
			}
			else if(j==7){
				p1 = 843;
			}
			peakFitter(files,detect,p1-30,p1+75);

		}
		for(int k=0;k<5;k++){
			// const char *files = Form("run0%d.root",i);
			const char *files = Form("/afs/crc.nd.edu/group/nsl/activetarget/data/24Mg_alpha_gamma/spectra/run0%d.root",i);
			const char*detect = Form("h1-%d",k);
			// cout <<files << "  "<< detect <<endl;
			if(k==0){
				p1 = 831;
			}
			else if(k==1){
				p1 = 840;
			}
			else if(k==2){
				p1 = 846;
			}
			else if(k==3){
				p1 = 827;
			}
			else if(k==4){
				p1 = 832;
			}
			peakFitter(files,detect,p1-30,p1+75);
		}
	}

	return 0;
}
