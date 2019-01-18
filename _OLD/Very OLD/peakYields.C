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


/*=============================FITTING=======================================*/

double background(double *x, double *par){
	return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
}


double gauss(double *x, double *par){
	// double arg = 0;
	// if (par[2]!=0) arg = (x[0] - par[1])/par[2];
	// double fitval = par[0]*TMath::Exp(-0.5*arg*arg);
	// return fitval;
	double norm = par[0];
	double mean = par[1];
	double sigma = par[2];
	return norm*TMath::Gaus(x[0],mean,sigma,kTRUE);
}


double func(double *x, double *par) {
   return background(x,par) + gauss(x,&par[3]);
}


void peakFitter(const char* fileName,const char* detector,int low, int high){

	TFile *fyield = new TFile(fileName);
	TH1D *hyield = static_cast<TH1D*>(fyield->Get(detector));
	gStyle->SetOptFit(1111);

	TCanvas *c0 = new TCanvas("c0","c0",600,800);
	c0->Divide(1,2);
	c0->Update();
	c0->cd(1);

	hyield->Draw();
	hyield->GetXaxis()->SetRangeUser(230,390);


	TF1 *ffit = new TF1("ffit",func,low,high,6);
	ffit->SetParNames("a0","a1","a2","norm","mean","sigma");
	ffit->SetNpx(500);
	ffit->SetParameters(3000,0,0,5000,290,5);   //// give it good range

	// ffit->SetParLimits(0,-1e3,1e3)
	// ffit->SetParLimits(1,-15,15);
    // ffit->FixParameter(2,0);			// Makes it a linear background (0)*x^2
    ffit->SetParLimits(3,1e2,1e6);
    ffit->SetParLimits(4,low,high);	///// give it good range
    ffit->SetParLimits(5,2.5,15);
	// ffit->SetParLimits(3,-5,5);

	ffit->SetLineColor(kRed);
	hyield->Fit("ffit","QR");
	TF1 *back = new TF1("back",background,low,high,3);
	back->SetLineColor(kMagenta);
	back->SetNpx(500);

	TF1 *peak = new TF1("peak",gauss,low,high,3);
	peak->SetLineColor(kBlue);
	peak->SetNpx(500);

	double par[6];
	ffit->GetParameters(par);
	back->SetParameters(par[0],par[1],par[2]);
	peak->SetParameters(par[3],par[4],par[5]);

	back->Draw("SAME");
	peak->Draw("SAME");


	c0->cd(2);
	// c0->Flush();


	TH1D *ysubtracted = new TH1D();
	ysubtracted = (TH1D*)hyield->Clone("ysubtracted");
	c0->Update();
    c0->Flush();

	// Calculate errors
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

	ysubtracted->GetXaxis()->SetRangeUser(par[4] - 3*par[5] , par[4] + 3*par[5]);
	ysubtracted->Draw();


	// ysubtracted->SetLineColor(kBlue);
	ysubtracted->Fit("peak","QR");
	peak->Draw("SAME");


	string runNum = fileName;
	// cout << "________________________________________________________" << endl;
	// cout << runNum << endl;
	// runNum = runNum.substr(4,3);
	runNum = runNum.substr(73,3);
	// cout << runNum << endl;
	string detNum = detector;


	double intg = peak->Integral(par[4] - 3*par[5] , par[4] + 3*par[5],1.e-8);
	double binw = hyield->GetBinWidth(1);
	double yield = intg/binw;

	double intg_err = peak->IntegralError(par[4] - 3*par[5] , par[4] + 3*par[5]);
	double yield_err = intg_err/binw;

	c0->SaveAs(Form("YIELDS/run0%s/det_%s_Fit.pdf",runNum.c_str(),detNum.c_str()));


	ofstream myfile;
	myfile.open ("yields.csv",std::ios::app);
	myfile << Form("run0%s",runNum.c_str()) <<","<< Form("det_%s",detNum.c_str())<<","<< yield <<","<< yield_err << "\n";
	myfile.close();


	c0->Clear();
	fyield->Close();
	delete c0;
	delete ffit;
	gROOT->Reset();

}

/*===========================================================================*/




/*==============================MAIN=========================================*/

int peakYields(){

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

		const char *files = Form("/afs/crc.nd.edu/group/nsl/activetarget/data/24Mg_alpha_gamma/spectra/run0%d.root",i);

		int low = 230;
		int high = 330;

		// const char *files = Form("run0%d.root",i);

		string runNum = files;
		runNum = runNum.substr(4,3);
		gSystem->Exec(Form("mkdir YIELDS/run0%d",i));

		for(int j=0;j<8;j++){
			// const char *files = Form("run0%d.root",i);
			const char *files = Form("/afs/crc.nd.edu/group/nsl/activetarget/data/24Mg_alpha_gamma/spectra/run0%d.root",i);
			const char *detect = Form("h0-%d",j); //   run0158.root  -  run0409.root
			cout <<files << "  "<< detect <<endl;
			peakFitter(files,detect,low,high);

		}
		for(int k=0;k<5;k++){
			// const char *files = Form("run0%d.root",i);
			const char *files = Form("/afs/crc.nd.edu/group/nsl/activetarget/data/24Mg_alpha_gamma/spectra/run0%d.root",i);
			const char*detect = Form("h1-%d",k);
			cout <<files << "  "<< detect <<endl;
			peakFitter(files,detect,low,high);

		}
	}

	return 0;
}
