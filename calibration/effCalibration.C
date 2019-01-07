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
#include "TVectorD.h"
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


	// Initial guess of centroid position
	double centroid = (low+high)/2;


	// Get histograms from root file
	TH1D *hyield = static_cast<TH1D*>(fyield->Get(detector));
	TH1D *ysubtracted = static_cast<TH1D*>(hyield->Clone("ysubtracted"));
	HBACK = static_cast<TH1D*>(fbackground->Get(detector));
	HBACK->SetLineColor(kBlack);
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


	// Background subtracted spectrum
	c0->cd(2);


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


	// Fit definitions
	TF1 *ffit = new TF1("ffit",func,low,high,6);
	ffit->SetParNames("a0","a1","a2","norm","mean","sigma");
	ffit->SetLineColor(kGreen);
	ffit->SetNpx(1e5);

	TF1 *fback = new TF1("fback",bgfunc,low,high,3);
	fback->SetParNames("a0","a1","a2");
	fback->SetLineColor(kCyan);
	fback->SetNpx(1e5);


	// Initial guess of parameters
	ffit->SetParameters(0,0.9,0,5000,centroid,10);

	ffit->FixParameter(2,0);			// Makes it a linear background (0)*x^2
    ffit->SetParLimits(3,0,1e6);		// Peak Area
	ffit->SetParLimits(4,low,high);		// Peak Centroid
    ffit->SetParLimits(5,2,55);			// Peak Width


	// Perform the fit
	TFitResultPtr r = ysubtracted->Fit("ffit","SQR");  // TFitResultPtr contains the TFitResult
	TMatrixDSym cov = r->GetCovarianceMatrix();   //  to access the covariance matrix

	// cout << TMath::Sqrt(cov[3][3]) << " " << r->ParError(3) << endl;


	// Draw the background subtracted histogram
	ysubtracted->GetXaxis()->SetRangeUser(low-20,high+20);
	ysubtracted->Draw();


	// Draw fits
	ffit->Draw("SAME");

	// Store fit parameters in an array
	double par[6];
	ffit->GetParameters(par);
	fback->SetParameters(par[0],par[1],par[2]);
	fback->Draw("SAME");

/*
	// Plot markers showing interval of integration
	TBox *box1 = new TBox(par[4]-3*par[5]-.25,-5,par[4]-3*par[5]+.25,100);
	TBox *box2 = new TBox(par[4]+3*par[5]-.25,-5,par[4]+3*par[5]+.25,100);
	box1->SetFillColor(2);
	box2->SetFillColor(2);
	box1->Draw("SAME");
	box2->Draw("SAME");
// */

	// Calculate the peak area and error
	string detNum = detector;
	double area_err;
	double area = ysubtracted->IntegralAndError(par[4]-3*par[5],par[4]+3*par[5],area_err,"width");
	double mean = par[4];
	double sigma = par[5];
	double A = par[3];
	double A_err = ffit->GetParError(3);


	// cout<<area<<"\t"<<A<<endl;
	// cout <<area_err<<"\t"<<area_err/area<<"\t"<<A_err<<"\t"<<A_err/A<<endl;
	area = A;
	area_err = A_err;


	// Save fits and calibration information
	if (loop==0) {
		c0->SaveAs(Form("calPlots/60Co_1173peak/det_%s_Fit.png",detNum.c_str()));

		ofstream myfile;
		myfile.open ("60Co_1173cal.csv",std::ios::app);
		myfile<< Form("det_%s",detNum.c_str())<<","<<mean<<","<<sigma<<","<<area
				<<","<<area_err<<","<<runTime0<<"\n";
		myfile.close();
	}
	else if (loop==1) {
		c0->SaveAs(Form("calPlots/60Co_1332peak/det_%s_Fit.png",detNum.c_str()));

		ofstream myfile;
		myfile.open ("60Co_1332cal.csv",std::ios::app);
		myfile<< Form("det_%s",detNum.c_str())<<","<<mean<<","<<sigma<<","<<area
				<<","<<area_err<<","<<runTime0<<"\n";
		myfile.close();
	}
	else if (loop==2) {
		c0->SaveAs(Form("calPlots/137Cs_661peak/det_%s_Fit.png",detNum.c_str()));

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

		// Wider fits to capture more background
		int low1173[] = {1030,1115,1136,1100,1070,1140,1120,1120,1100,1120,1140,1100,1110};
		int high1173[] = {1110,1200,1210,1200,1170,1230,1220,1220,1200,1220,1240,1200,1200};
		int low1332[] = {1170,1260,1275,1250,1220,1280,1270,1270,1260,1280,1270,1260,1260};
		int high1332[] = {1250,1360,1375,1360,1320,1400,1370,1380,1360,1370,1410,1360,1360};
		int low661[] = {570,623,637,625,600,636,638,630,617,620,630,610,620};
		int high661[] = {650,690,700,690,675,705,690,700,690,690,700,695,685};

		// // Tight fits
		// int low1173[] = {1040,1123,1137,1120,1090,1145,1130,1135,1120,1138,1136,1113,1123};
		// int high1173[] = {1100,1194,1205,1190,1155,1220,1200,1210,1190,1205,1222,1190,1190};
		// int low1332[] = {1180,1280,1288,1269,1237,1298,1283,1290,1270,1290,1297,1268,1275};
		// int high1332[] = {1243,1350,1362,1346,1312,1380,1357,1367,1350,1360,1380,1350,1350};
		// int low661[] = {590,628,646,631,615,645,638,638,630,638,638,619,630};
		// int high661[] = {630,681,690,680,662,695,682,688,680,677,690,677,675};

		// 1173 keV gamma
		if (strcmp(file,"60Co/24Mg_run/run0419.root") == 0 && i == 0){
            low = low1173;
			high = high1173;

			ofstream myfile;
			myfile.open ("60Co_1173cal.csv",std::ios::out);
			myfile<<"Detector"<<","<<"Centroid"<<","<<"Width"<<","<<"Area"<<","
					<<"Area err"<<","<<"Runtime"<<"\n";
			myfile.close();
		}

		// 1332 keV gamma
		else if (strcmp(file,"60Co/24Mg_run/run0419.root") == 0 && i == 1){
			low = low1332;
			high = high1332;

			ofstream myfile;
			myfile.open ("60Co_1332cal.csv",std::ios::out);
			myfile<<"Detector"<<","<<"Centroid"<<","<<"Width"<<","<<"Area"<<","
					<<"Area err"<<","<<"Runtime"<<"\n";
			myfile.close();
		}

		// 661 keV gamma
		else if (strcmp(file,"137Cs/24Mg_run/run0420.root")==0 && i == 2){
			low = low661;
			high = high661;

			ofstream myfile;
			myfile.open ("137Cs_661cal.csv",std::ios::out);
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
