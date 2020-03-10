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


// Flourine line from first run, these get updated run by run with the new ranges
double fLow[] = {113,123,129,125,124,129,126,127,126,128,127,126,123};
double fHigh[] = {140,143,144,144,141,147,146,146,143,143,147,143,147};


// if 1, use local, if 0 use CRC
int loc = 0;

// if 1, save plots, if 0 don't save
int plot = 0;


void peakFitter(TFile *TFitOut, const char *fileName, const char *detector,
	int detLoop){

	// Reset global variables
	gROOT->Reset();

	// Get root file
	TFile *fyield = new TFile(fileName,"READ");

	// Change output directory to TFitOut
	TFitOut->cd();

	// Get integrated charge from root file
	TH1D *qcharge = static_cast<TH1D*>(fyield->Get("h1-7"));
	double charge = qcharge->GetEntries();

	// Get histograms from root file
	TH1D *hyield = static_cast<TH1D*>(fyield->Get(detector));

	gStyle->SetOptFit(1111);

	// Create the canvas
	TCanvas *c0 = new TCanvas("c0","c0",1920,1080);
	// c0->Divide(1,2);
	c0->Update();

	// Energy calibrate the spectra
	vector < double >  peak2SpecLow ({1285,1400,1400,1390,1360,1430,1410,1420,1410,1410,1420,1410,1400});
	vector < double >  peak2SpecHigh ({1360,1480,1510,1485,1440,1510,1480,1500,1480,1510,1530,1490,1490});

	// Find the BG peak positions in runs
	vector < double > peak2Position;
	peak2Position = iterative_double_gauss_peak(peak2SpecLow[detLoop],peak2SpecHigh[detLoop],hyield);

	vector < double > fPosition;
	fPosition = single_gauss_peak(fLow[detLoop],fHigh[detLoop],hyield);

	vector < double > calibrators;

	TH1D *h3 = new TH1D(Form("h3 - det-%i",detLoop),Form("h3 - ECAL - det-%i",detLoop),8192,0,8192);
	calibrators = calibrate(109.9,1468,fPosition[0],peak2Position[0],h3,hyield);


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

	string runNum = fileName;
		if (loc==1) runNum = runNum.substr(9,3);
		else runNum = runNum.substr(73,3);

    int n = h3->GetNbinsX();

	string outFileName = Form("E_cal_spectras/Spectra_run_%s_det-%i",runNum.c_str(),detLoop);
	cout << outFileName << endl;
	ofstream myfile (outFileName.c_str(), ios::out);
	// ofstream myfile (Form("run_%s_det-%i_Spectra.dat",runNum.c_str(),detLoop),ios::out);
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
	fyield->Close();
	delete c0;
}



/*==============================MAIN=========================================*/
void paperPlotter(){

	// Suppress certain root print statements
	gErrorIgnoreLevel = kWarning;

	// When running on CRC
	const char *path = "/afs/crc.nd.edu/user/s/saguilar/Group/24Mg_ap/";
	chdir(path);

	const char *detect;
	const char *files;

	// Loop through runs:
	cout << "\nBEGINNING PEAK FITTING:" << endl;

	// Record start time
	auto start = std::chrono::high_resolution_clock::now();

	int fileNum = 1;
	int runStart = 159 ;
	int upToRun;

	if (loc==1) upToRun = 160;
	else upToRun = 410;

	for(int i=runStart;i<upToRun;i++){

		// Skip bad runs
		if(i>=163 && i<=166) continue;
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

		TString runNum_TString;

		if(i < 100) runNum_TString = "00";
		else if((i >= 100) && (i < 1000)) runNum_TString = "0";
		else runNum_TString = "";

		runNum_TString += i;	// Should be format 0001 -> 9999
		const char *runNum_String = (const char*)runNum_TString;


		// // Save TFitResult results.
		TFile *ff = new TFile(Form("E_cal_spectras/run_%s.root",runNum_String),"RECREATE");

		// Loop through detectors on board
		for(int j=0;j<13;j++){ // 13

			if (loc==1) files = Form("runs/run%s.root",runNum_String);
			else files = Form("/afs/crc.nd.edu/group/nsl/activetarget/data/24Mg_alpha_gamma/spectra/run0%d.root",i);

			if (j<8){
				detect = Form("h0-%d",j);
			}
			else{
				detect = Form("h1-%d",j-8);
			}

			// Perform peak fitting
			peakFitter(ff,files,detect,j);
		}

		// // Write all TObjects in memory (TFitResult) to TFile
		ff->Write();
		ff->Close();

		cout << Form("Fitting run_%s complete",runNum_String) << endl;
		fileNum+=1;
		delete ff;
	}

	// Record end time
	auto finish = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> elapsed = finish - start;

	cout << "Elapsed time: " << elapsed.count() << " s\n";
}
