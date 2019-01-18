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
#include "TError.h"

// Fitting routines
#include "calibration/fitFunctions.h"

// Gain match routine
#include "calibration/gainMatch.h"

//vector global for a background histogram
TH1D *HBACK;


void peakFitter(const char *fileName,const char *fileBack,const char *detector,int low,int high,int detectorLoop){

	TFile *fyield = new TFile(fileName);
	TFile *fbackground = new TFile(fileBack);	// bg spectra


	// Get runtime info from root file
	TVectorD *run_t0 = static_cast<TVectorD*>(fyield->Get("lastTimeStampSeconds-0"));
	TVectorD *back_t0 = static_cast<TVectorD*>(fbackground->Get("lastTimeStampSeconds-0"));


	// Run time and background time for each board
	double runTime0 = (*run_t0)[0];
	double backTime0 = (*back_t0)[0];


	// Scale background spectra to ratio of run time
	double scale0 = runTime0/backTime0;


	// Get histograms from root file
	TH1D *hyield = static_cast<TH1D*>(fyield->Get(detector));
	HBACK = static_cast<TH1D*>(fbackground->Get(detector));


	TH1D *qcharge = static_cast<TH1D*>(fyield->Get("h1-7"));
	double charge = qcharge->GetEntries();
	gStyle->SetOptFit(1111);


	// Create the canvas
	TCanvas *c0 = new TCanvas("c0","c0",600,800);
	c0->Divide(1,2);
	c0->Update();
	c0->cd(1);


	// Only gain match longer runs where background peaks appear
	if(runTime0 > 500){

		/////////////////////////
		///    GAIN MATCH     ///
		/////////////////////////

		HBACK->Scale(scale0);	// Scale the background
		HBACK->SetDirectory(0);

		// // Background peaks from run0422.root to be gain matched
		// int peak1Back[] = {1323,1439,1451,1431,1397,1466,1447,1456,1434,1452,1463,1431,1437};
		// int peak2Back[] = {1789,1969,1980,1968,1930,2025,1981,2000,1954,1991,2026,1987,1978};


		// Background peak ranges
		int peak1BackLow[] = {1270,1380,1385,1370,1330,1390,1380,1380,1370,1390,1390,1350,1380};
		int peak1BackHigh[] = {1380,1500,1510,1485,1450,1540,1510,1520,1490,1510,1530,1500,1490};

		int peak2BackLow[] = {1750,1927,1935,1930,1887,1980,1942,1946,1912,1948,1980,1942,1942};
		int peak2BackHigh[] = {1830,2015,2025,2010,1972,2064,2017,2050,1993,2029,2073,2027,2020};

		// Find BG peak positions using same procedure for finding in spectrum
		double peak1Back = iterative_double_gauss_peak(peak1BackLow[detectorLoop],peak1BackHigh[detectorLoop],HBACK);
		double peak2Back = iterative_single_gauss_peak(peak2BackLow[detectorLoop],peak2BackHigh[detectorLoop],HBACK);

		// Estimated peak position range MAYBE READ FILE FOR THESE VALUES
		int peak1SpecLow[] = {1285,1400,1400,1390,1360,1430,1410,1420,1410,1430,1420,1410,1410};
		int peak1SpecHigh[] = {1360,1480,1485,1480,1440,1510,1480,1500,1480,1510,1510,1490,1490};

		int peak2SpecLow[] = {1740,1885,1923,1890,1841,1934,1895,1913,1890,1924,1936,1895,1898};
		int peak2SpecHigh[] = {1840,2034,2048,2036,2007,2090,2037,2065,2018,2059,2085,2097,2035};


		// Find the peak positions in runs given an estimated range
		double peak1Position = iterative_double_gauss_peak(peak1SpecLow[detectorLoop],peak1SpecHigh[detectorLoop],hyield);
		double peak2Position = iterative_single_gauss_peak(peak2SpecLow[detectorLoop],peak2SpecHigh[detectorLoop],hyield);


		// New gain matched BG histogram
		TH1D *h2 = new TH1D("h2","h2",8192,0,8192);

		h2 = gain_match(peak1Position,peak2Position,peak1Back,peak2Back,h2,HBACK);

		// Draw results
		// hyield->Draw();
		// hyield->GetXaxis()->SetRangeUser(1200,1600);
		// hyield->GetXaxis()->SetRangeUser(1500,2400);
		// HBACK->Draw("SAME");			// Not gain matched BG
		HBACK->Draw();
		HBACK->GetXaxis()->SetRangeUser(1500,2400);
		HBACK->SetLineColor(kGreen);
		h2->SetLineColor(kBlue);
		h2->Draw("SAME");				// Gain matched BG


		// Prepare bg subtracted histogram
		TH1D *ysubtracted = static_cast<TH1D*>(hyield->Clone("ysubtracted"));


		// Recalculate errors manually
		for(int i = 0; i<=ysubtracted->GetNbinsX(); i++){
			double yval = ysubtracted->GetBinContent(i);
			double yval2 = HBACK->GetBinContent(i);	//fback->Eval(ysubtracted->GetBinCenter(i));
			double yerr = ysubtracted->GetBinError(i);
			double yerr2 = HBACK->GetBinError(i);	// Takes the error from the bg spectrum

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


	}
	else{
		cout << Form("Runtime too short: %f seconds. Not performing gain match",runTime0);

		// Draw results
		hyield->Draw();
		hyield->GetXaxis()->SetRangeUser(700,2000);
	}

	c0->cd(2);

	// Calculate the peak area and error
	// double area_err;
	// double area = ysubtracted->IntegralAndError(par[4]-3*par[5],par[4]+3*par[5],area_err,"width");
	// double mean = par[4];
	// double sigma = par[5];
	double A = par[3];
	double A_err = ffit->GetParError(3);

	double yield = A/charge;
	double yield_err = A_err/charge;

	double chi2NDF = ffit->GetChisquare()/ffit->GetNDF();
	double goodFit;
	if (chi2NDF <= 2 && chi2NDF >=.5 ) goodFit = 0;
	else goodFit = 1;

	string runNum = fileName;
	runNum = runNum.substr(4,3);
	// runNum = runNum.substr(73,3);

	string detNum = detector;

	c0->SaveAs(Form("peakAreasP2/run0%s/det_%s_Fit.png",runNum.c_str(),detNum.c_str()));


	// ofstream myfile;
	// myfile.open ("peakAreasP2.csv",std::ios::app);
	// myfile<<Form("run0%s",runNum.c_str())<<","<< Form("det_%s",detNum.c_str())<<
	// 			","<<yield<<","<<yield_err<<","<<goodFit<<"\n";
	// myfile.close();


	c0->Clear();
	fyield->Close();
	fbackground->Close();

	// delete->fyield;
	// delete->fbackground;
	delete c0;
	// delete ffit;
	// delete fback;
	gROOT->Reset();

}



/*==============================MAIN=========================================*/
void peakAreasP2(){

	// Suppress certain root print statements
	gErrorIgnoreLevel = kWarning;

	// // When running on CRC
	// const char *path = "/afs/crc.nd.edu/user/s/saguilar/Group/24Mg_ap/";
	// chdir(path);

	// Create file directory for output incase its not made
	try {
		cout << "\nATTEMPTING TO CREATE OUTPUT FILE DIRECTORIES" << endl;
		gSystem->Exec(Form("mkdir peakAreasP2"));
	}catch(...){}

	// Background Spectrum file
	const char *fileBackground = "calibration/background/run0422.root";

	const char *detect;
	const char *files;

	// Prepare structure of data output in CSV file
	ofstream myfile;
	myfile.open ("peakAreasP2.csv",std::ios::out);
	myfile<<"Run"<<","<<"Detector"<<","<<"Area"<<","<<"Area err"<<","<<"Fit Status"<<"\n";
	myfile.close();

	// Loop through runs: 159-410
	cout << "\nBEGINNING PEAK FITTING:" << endl;
	int fileNum = 1;
	for(int i=159;i<162;i++){ // 410

		// Skip bad runs
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

		double p2;

		// Create file directory for output
		try {
			gSystem->Exec(Form("mkdir peakAreasP2/run0%d",i));
		}catch(...){}


		// Loop through detectors on board 1 (0-7) and board 2 (8-12)
		for(int j=0;j<13;j++){

			files = Form("run0%d.root",i);

			// const char *files = Form("/afs/crc.nd.edu/group/nsl/activetarget/data/24Mg_alpha_gamma/spectra/run0%d.root",i);
			if (j<8){
				detect = Form("h0-%d",j);
			}
			else{
				detect = Form("h1-%d",j-8);
			}

			// Estimated peak positions from calibration
			int peakPos[] = {927,1002,1014,999,973,1022,1009,1014,1000,1011,1020,996,1000};
			p2 = peakPos[j];

			peakFitter(files,fileBackground,detect,p2-60,p2+60,j);
		}
		cout << Form("Fitting %d complete",fileNum) << endl;
		fileNum+=1;
	}
}
