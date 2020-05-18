#include <iostream>
using std::cout;
using std::endl;
using std::ofstream;

#include <ostream>
#include <fstream>

#include <unistd.h>

#include <chrono>  // for high_resolution_clock

#include <string>
using std::string;

#include <vector>
using std::vector;

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
#include "TFitResult.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TGaxis.h"


// Fitting routines
// #include "/Users/sebastian/Desktop/24Mg_agp/calibration/fitFunctions.h"
#include "calibration/fitFunctions.h"

// Gain match routine
// #include "/Users/sebastian/Desktop/24Mg_agp/calibration/gainMatch.h"
#include "calibration/gainMatch.h"



// if 1, use local, if 0 use CRC
int loc = 1;

// if 1, save plots, if 0 don't save
int plot = 0;


void eCalibrate(TFile *TFitOut,
	const vector < double > &peak1BackLow, const vector < double > &peak1BackHigh,
	const vector < double > &peak2BackLow, const vector < double > &peak2BackHigh,
	const char *detector, int detLoop){

	// Reset global variables
	gROOT->Reset();

	// Get root file for raw BG spectra
	TFile *fback = new TFile("calibration/background/run0422.root","READ");

	// Change output directory to TFitOut
	TFitOut->cd();


	// Get histograms from root file
	TH1D *hback = static_cast<TH1D*>(fback->Get(detector));

	gStyle->SetOptFit(1111);

	// Create the canvas
	TCanvas *c0 = new TCanvas("c0","c0",1920,1080);
	// c0->Divide(1,2);
	c0->Update();


	// Find the BG peak positions in runs
	vector < double > peak1Position;
	peak1Position = BGdouble_gauss_peak(peak1BackLow[detLoop],peak1BackHigh[detLoop],hback);

	vector < double > peak2Position;
	peak2Position = BGsingle_gauss_peak(peak2BackLow[detLoop],peak2BackHigh[detLoop],hback);

	vector < double > calibrators;
	// hback->Draw();
	TH1D *h3 = new TH1D(Form("h3 - det-%i",detLoop),Form("h3 - ECAL - det-%i",detLoop),8192,0,8192);
	calibrators = calibrate(1468,2615,peak1Position[0],peak2Position[0],h3,hback);


	h3->GetXaxis()->SetRangeUser(0,3500);
	gPad->SetLogy();
	h3->Draw();
	h3->SetStats(kFALSE);
	h3->GetXaxis()->SetTitle("Energy (keV)");
	h3->GetYaxis()->SetTitle("Counts / Bin");
	h3->GetXaxis()->CenterTitle();
	h3->GetYaxis()->CenterTitle();
	h3->SetTitle("");

	h3->Write(Form("det-%i",detLoop));

    int n = h3->GetNbinsX();

	string outFileName = Form("E_cal_spectras/Spectra_det-%i",detLoop);
	cout << outFileName << endl;
	ofstream myfile (outFileName.c_str(), ios::out);

	if (myfile.is_open()) {
		for (int i=1; i<=n; i++) {
			myfile << h3->GetBinLowEdge(i)+h3->GetBinWidth(i)/2<<"\t"<<h3->GetBinContent(i)<<"\n";
			// myfile << Form("%g\t%g\n",h3->GetBinLowEdge(i)+h3->GetBinWidth(i)/2,h3->GetBinContent(i));
	    }
		myfile.close();
	}
	else cout << "something bad\n";




	c0->Clear();
	delete h3;
	fback->Close();
	delete c0;
}



/*==============================MAIN=========================================*/
void eCalBGSpectra(){

	// Suppress certain root print statements
	gErrorIgnoreLevel = kWarning;

	// When running on CRC
	// const char *path = "/afs/crc.nd.edu/user/s/saguilar/Group/24Mg_ap/";
	// chdir(path);

	// Background Spectrum file
	const char *fileBackground = "calibration/background/run0422.root";

	// BG spectra peak ranges
	// 1468 keV: 138La + Ba EC
	const vector < double >  peak1BackLow ({1270,1380,1385,1370,1330,1390,1380,1380,1370,1390,1400,1350,1380});
	const vector < double >  peak1BackHigh ({1380,1500,1510,1485,1450,1540,1510,1520,1490,1510,1530,1520,1500});

	// 2615 keV: 208Ti
	const vector < double >  peak2BackLow ({2250,2450,2450,2450,2400,2500,2450,2500,2450,2480,2500,2500,2480});
	const vector < double >  peak2BackHigh ({2350,2600,2600,2600,2530,2650,2600,2620,2550,2600,2650,2620,2580});

	const char *detect;

	// Record start time
	auto start = std::chrono::high_resolution_clock::now();

	// Save TFitResult results.
	TFile *ff = new TFile("E_cal_spectras/eCalBGSpectra.root","RECREATE");

	// Loop through detectors on board
	for(int j=0;j<13;j++){ // 13

		if (j<8){
			detect = Form("h0-%d",j);
		}
		else{
			detect = Form("h1-%d",j-8);
		}

		// Perform calibration
		eCalibrate(ff,peak1BackLow,peak1BackHigh,peak2BackLow,peak2BackHigh,
					detect,j);
	}

	// // Write all TObjects in memory (TFitResult) to TFile
	ff->Write();
	ff->Close();
	delete ff;

	// Record end time
	auto finish = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> elapsed = finish - start;

	cout << "Elapsed time: " << elapsed.count() << " s\n";
}
