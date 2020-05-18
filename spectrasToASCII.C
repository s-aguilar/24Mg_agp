#include <iostream>
using std::cout;
using std::endl;

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

// Fitting routines
#include "calibration/fitFunctions.h"

// Gain match routine
#include "calibration/gainMatch.h"



// if 1, use local, if 0 use CRC
int loc = 0;


// Flourine line from first run, these get updated run by run with the new ranges
double fLow[] = {113,123,129,125,124,129,126,127,126,128,127,126,123};
double fHigh[] = {140,143,144,144,141,147,146,146,143,143,147,143,147};

// global array for gain match constants
double _a_gain[13];
double _b_gain[13];

// global array for calibration constants
double _a_calibrator[13];
double _b_calibrator[13];


void peakFitter(const char *fileName, const char *fileBack, const char *detector,
	const vector < double > &peak1BackLow, const vector < double > &peak1BackHigh,
	const vector < double > &peak2BackLow, const vector < double > &peak2BackHigh,
	vector < double > &peak1SpecLow, vector < double > &peak1SpecHigh,
	vector < double > &peak2SpecLow, vector < double > &peak2SpecHigh, int detLoop){


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
	TH1D *HBACK = static_cast<TH1D*>(fbackground->Get(detector));

	// Scale the background
	HBACK->Scale(scale0);
	HBACK->SetDirectory(0);

	// Get integrated charge information
	TH1D *qcharge = static_cast<TH1D*>(fyield->Get("h1-7"));
	double charge = qcharge->GetEntries();
	gStyle->SetOptFit(1111);

	// Create the canvas
	TCanvas *c0 = new TCanvas("c0","c0",1920,1080);
	c0->Divide(1,2);
	c0->Update();

	c0->cd(1);

	// Prepare bg subtracted histogram and other histograms
	TH1D *ysubtracted = static_cast<TH1D*>(hyield->Clone("ysubtracted"));
	TH1D *h2 = new TH1D(Form("h2 - det-%i",detLoop),Form("h2 - det-%i",detLoop),8192,0,8192);
	TH1D *h3 = new TH1D(Form("h3 - det-%i",detLoop),Form("h3 - ECAL - det-%i",detLoop),8192,0,8192);
	TH1D *h4 = new TH1D(Form("h4 - det-%i",detLoop),Form("h4 - ECAL - det-%i",detLoop),8192,0,8192);

	double area;
	double area_err;
	double chi2NDF;
	double sig1;
	double sig2;

	double a = 0;
	double b = 0;
	double a_gm = 0;
	double b_gm = 0;
	double a_cal = 0;
	double b_cal = 0;

	double linear = 0;
	double offset = 0;


	// Only gain match longer runs where background peaks appear
	if(runTime0 > 850){

		/////////////////////////
		///    GAIN MATCH     ///
		/////////////////////////

		// Find BG peak positions
		vector < double > peak1Back;
		vector < double > peak2Back;

		peak1Back = iterative_single_gauss_peak(peak1BackLow[detLoop],peak1BackHigh[detLoop],HBACK);
		peak2Back = iterative_double_gauss_peak(peak2BackLow[detLoop],peak2BackHigh[detLoop],HBACK);

		// Find the BG peak positions in runs
		vector < double > peak1Position;
		vector < double > peak2Position;

		peak1Position = iterative_single_gauss_peak(peak1SpecLow[detLoop],peak1SpecHigh[detLoop],hyield);
		peak2Position = iterative_double_gauss_peak(peak2SpecLow[detLoop],peak2SpecHigh[detLoop],hyield);

		// Add new ranges for each detector
		peak1SpecLow[detLoop] = peak1Position[1];
		peak1SpecHigh[detLoop] = peak1Position[2];
		peak2SpecLow[detLoop] = peak2Position[1];
		peak2SpecHigh[detLoop] = peak2Position[2];

		// New gain matched BG histogram -> h2
		vector < double > gain;
		gain = gain_match(peak1Position[0],peak2Position[0],peak1Back[0],peak2Back[0],h2,HBACK);
		a = gain[0];
		b = gain[1];
		a_gm = gain[0];
		b_gm = gain[1];

		_a_gain[detLoop] = a_gm;
		_b_gain[detLoop] = b_gm;

		// Draw results
		hyield->Draw();
		hyield->GetXaxis()->SetRangeUser(600,1200);
		hyield->SetStats(kFALSE);


		// Recalculate errors manually
		for(int i = 0; i<=ysubtracted->GetNbinsX(); i++){
			double yval = ysubtracted->GetBinContent(i);
			double yval2 = h2->GetBinContent(i);	//fback->Eval(ysubtracted->GetBinCenter(i));
			double yerr = ysubtracted->GetBinError(i);
			double yerr2 = h2->GetBinError(i);	// Takes the error from the gain matched bg spectrum

			ysubtracted->SetBinContent(i,yval-yval2);
			ysubtracted->SetBinError(i,TMath::Sqrt(yerr*yerr+yerr2*yerr2));
		}

		ysubtracted->SetLineColor(kRed);
		ysubtracted->Draw("SAME");

		c0->cd(2);

		// Now I've gainmatched and subtracted the BG -> ysubtracted, find flourine line position
		vector < double > fPosition;
		fPosition = single_gauss_peak(fLow[detLoop],fHigh[detLoop],ysubtracted);

		// Found the lines update their ranges run by run
		fLow[detLoop] = fPosition[1];
		fLow[detLoop] = fPosition[2];

		vector < double > calibrators;


		calibrators = calibrate(109.9,1468,fPosition[0],peak2Position[0],h4,hyield);


		calibrators = calibrate(109.9,1468,fPosition[0],peak2Position[0],h3,ysubtracted);
		a_cal = calibrators[0];
		b_cal = calibrators[1];

		// Use these calibrators for when the runtime is too short and no gainmatch fitting will be performed
		_a_calibrator[detLoop] = a_cal;
		_b_calibrator[detLoop] = b_cal;
	}
	else{
		// cout << Form("Runtime too short: %f seconds. Not performing gain match",runTime0);
		a_cal = _a_calibrator[detLoop];
		b_cal = _b_calibrator[detLoop];
		a_gm = _a_gain[detLoop];
		b_gm = _b_gain[detLoop];

		h2 = gain_match2(a_gm,b_gm,h2,HBACK);

		// Draw results
		hyield->Draw();
		hyield->GetXaxis()->SetRangeUser(600,1200);
		hyield->SetStats(kFALSE);

		// Recalculate errors manually
		for(int i = 0; i<=ysubtracted->GetNbinsX(); i++){
			double yval = ysubtracted->GetBinContent(i);
			double yval2 = h2->GetBinContent(i);	//fback->Eval(ysubtracted->GetBinCenter(i));
			double yerr = ysubtracted->GetBinError(i);
			double yerr2 = h2->GetBinError(i);	// Takes the error from the gain matched bg spectrum

			ysubtracted->SetBinContent(i,yval-yval2);
			ysubtracted->SetBinError(i,TMath::Sqrt(yerr*yerr+yerr2*yerr2));
		}

		ysubtracted->SetLineColor(kRed);
		ysubtracted->Draw("SAME");

		c0->cd(2);

		h3 = calibrate2(a_cal,b_cal,h3,ysubtracted);
		h4 = calibrate2(a_cal,b_cal,h4,hyield);
	}


	// h3->GetXaxis()->SetRangeUser(600,1200);
	// h3->Draw();
	// h3->SetStats(kFALSE);

	string runNum = fileName;
	if (loc==1) runNum = runNum.substr(9,3);
	else runNum = runNum.substr(73,3);

	string detNum = detector;


	int n = h3->GetNbinsX();
	string outFileName = Form("E_cal_spectras/BGSubSpectra_run_%s_det-%i",runNum.c_str(),detLoop);
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

	n = h4->GetNbinsX();
	string outFileName1 = Form("E_cal_spectras/Spectra_run_%s_det-%i",runNum.c_str(),detLoop);
	cout << outFileName1 << endl;
	ofstream myfile1 (outFileName1.c_str(), ios::out);
	if (myfile1.is_open()) {
		for (int i=1; i<=n; i++) {
			myfile1 << h4->GetBinLowEdge(i)+h4->GetBinWidth(i)/2<<"\t"<<h4->GetBinContent(i)<<"\n";
			// myfile << Form("%g\t%g\n",h3->GetBinLowEdge(i)+h3->GetBinWidth(i)/2,h3->GetBinContent(i));
	    }
		myfile1.close();
	}
	else cout << "something bad\n";


	c0->Clear();
	delete ysubtracted;
	delete h2;
	delete h3;
	delete h4;
	fyield->Close();
	fbackground->Close();
	delete c0;

	gROOT->Reset();
}



/*==============================MAIN=========================================*/
void spectrasToASCII(){

	// Suppress certain root print statements
	gErrorIgnoreLevel = kWarning;

	// When running on CRC
	const char *path = "/afs/crc.nd.edu/user/s/saguilar/Group/24Mg_ap/";
	chdir(path);

	// Background Spectrum file
	const char *fileBackground = "calibration/background/run0422.root";

	const char *detect;
	const char *files;


	// BG spectra peak ranges
	const vector < double >  peak1BackLow ({310,330,330,330,320,340,330,330,330,330,320,310,330});
	const vector < double >  peak1BackHigh ({350,370,380,370,360,380,380,380,380,380,380,375,380});
	const vector < double >  peak2BackLow ({1270,1380,1385,1370,1330,1390,1380,1380,1370,1390,1400,1350,1380});
	const vector < double >  peak2BackHigh ({1380,1500,1510,1485,1450,1540,1510,1520,1490,1510,1530,1520,1500});

	// BG peak ranges starting from run0159 \
		These get updated as runs progress
	vector < double >  peak1SpecLow ({300,320,330,330,320,330,330,330,330,330,330,330,330});
	vector < double >  peak1SpecHigh ({350,380,390,380,360,380,380,380,380,380,380,380,380});
	vector < double >  peak2SpecLow ({1285,1400,1400,1390,1360,1430,1410,1420,1410,1410,1420,1410,1400});
	vector < double >  peak2SpecHigh ({1360,1480,1510,1485,1440,1510,1480,1500,1480,1510,1530,1490,1490});


	// Loop through runs: 159-410
	cout << "\nBEGINNING PEAK FITTING:" << endl;
	int fileNum = 1;

	int upToRun;
	if (loc==1) upToRun = 175;
	else upToRun = 410;
	for(int i=160;i<upToRun;i++){

		// Skip bad runs
		if(i>=163 && i<=166) continue;
		else if(i>=168 && i<=171) continue;
		else if(i==182) continue;
		else if(i>=244 && i<=255) continue;
		else if(i==276) continue;
		else if(i==277) continue;
		else if(i==285) continue;
		else if(i==289) continue;
		else if(i==290) continue;
		else if(i==294) continue;
		else if(i==406) continue;


		// Loop through detectors on board 1 (0-7) and board 2 (8-12)
		for(int j=0;j<13;j++){ // 13

			if (loc==1) files = Form("runs/run0%d.root",i);
			else files = Form("/afs/crc.nd.edu/group/nsl/activetarget/data/24Mg_alpha_gamma/spectra/run0%d.root",i);

			if (j<8){
				detect = Form("h0-%d",j);
			}
			else{
				detect = Form("h1-%d",j-8);
			}


			// Perform peak fitting
			peakFitter(files,fileBackground,detect,peak1BackLow,peak1BackHigh,
				peak2BackLow,peak2BackHigh,peak1SpecLow,peak1SpecHigh,
				peak2SpecLow,peak2SpecHigh,j);
		}

		cout << Form("Fitting %d complete",fileNum) << endl;
		fileNum+=1;
	}
}



// #include <iostream>
// using std::cout;
// using std::endl;
// using std::ofstream;
//
// #include <ostream>
// #include <fstream>
//
// #include <unistd.h>
//
// #include <chrono>  // for high_resolution_clock
//
// #include <string>
// using std::string;
//
// #include <vector>
// using std::vector;
//
// #include "TROOT.h"
// #include "TF1.h"
// #include "TFile.h"
// #include "TH1D.h"
// #include "TMath.h"
// #include "TString.h"
// #include "TVectorD.h"
// #include "TCanvas.h"
// #include "TStyle.h"
// #include "TError.h"
// #include "TFitResult.h"
// #include "TArrow.h"
// #include "TLatex.h"
// #include "TAxis.h"
// #include "TGaxis.h"
//
//
// // Fitting routines
// // #include "/Users/sebastian/Desktop/24Mg_agp/calibration/fitFunctions.h"
// #include "calibration/fitFunctions.h"
//
// // Gain match routine
// // #include "/Users/sebastian/Desktop/24Mg_agp/calibration/gainMatch.h"
// #include "calibration/gainMatch.h"
//
//
// // Flourine line from first run, these get updated run by run with the new ranges
// double fLow[] = {113,123,129,125,124,129,126,127,126,128,127,126,123};
// double fHigh[] = {140,143,144,144,141,147,146,146,143,143,147,143,147};
//
//
// // if 1, use local, if 0 use CRC
// int loc = 0;
//
// // if 1, save plots, if 0 don't save
// int plot = 0;
//
//
// void peakFitter(TFile *TFitOut, const char *fileName, const char *fileBack,
// 	const char *detector, int detLoop){
//
// 	// Reset global variables
// 	gROOT->Reset();
//
// 	// Get root file
// 	TFile *fyield = new TFile(fileName,"READ");
// 	TFile *fbackground = new TFile(fileBack);	// bg spectra
//
// 	// Change output directory to TFitOut
// 	TFitOut->cd();
//
// 	// Get integrated charge from root file
// 	TH1D *qcharge = static_cast<TH1D*>(fyield->Get("h1-7"));
// 	double charge = qcharge->GetEntries();
//
// 	// Get histograms from root file
// 	TH1D *hyield = static_cast<TH1D*>(fyield->Get(detector));
//
// 	gStyle->SetOptFit(1111);
//
// 	// Create the canvas
// 	TCanvas *c0 = new TCanvas("c0","c0",1920,1080);
// 	// c0->Divide(1,2);
// 	c0->Update();
//
// 	// Energy calibrate the spectra
// 	vector < double >  peak2SpecLow ({1285,1400,1400,1390,1360,1430,1410,1420,1410,1410,1420,1410,1400});
// 	vector < double >  peak2SpecHigh ({1360,1480,1510,1485,1440,1510,1480,1500,1480,1510,1530,1490,1490});
//
// 	// Find the BG peak positions in runs
// 	vector < double > peak2Position;
// 	peak2Position = iterative_double_gauss_peak(peak2SpecLow[detLoop],peak2SpecHigh[detLoop],hyield);
//
// 	vector < double > fPosition;
// 	fPosition = single_gauss_peak(fLow[detLoop],fHigh[detLoop],hyield);
//
// 	vector < double > calibrators;
//
// 	TH1D *h3 = new TH1D(Form("h3 - det-%i",detLoop),Form("h3 - ECAL - det-%i",detLoop),8192,0,8192);
// 	calibrators = calibrate(109.9,1468,fPosition[0],peak2Position[0],h3,hyield);
//
//
// 	h3->GetXaxis()->SetRangeUser(0,3500);
// 	gPad->SetLogy();
// 	h3->Draw();
// 	h3->SetStats(kFALSE);
// 	h3->GetXaxis()->SetTitle("Energy (keV)");
// 	h3->GetYaxis()->SetTitle("Counts / Bin");
// 	h3->GetXaxis()->CenterTitle();
// 	h3->GetYaxis()->CenterTitle();
// 	h3->SetTitle("");
//
// 	h3->Write(Form("det-%i",detLoop));
//
// 	string runNum = fileName;
// 		if (loc==1) runNum = runNum.substr(9,3);
// 		else runNum = runNum.substr(73,3);
//
//     int n = h3->GetNbinsX();
//
// 	string outFileName = Form("E_cal_spectras/Spectra_run_%s_det-%i",runNum.c_str(),detLoop);
// 	cout << outFileName << endl;
// 	ofstream myfile (outFileName.c_str(), ios::out);
// 	// ofstream myfile (Form("run_%s_det-%i_Spectra.dat",runNum.c_str(),detLoop),ios::out);
// 	if (myfile.is_open()) {
// 		for (int i=1; i<=n; i++) {
// 			myfile << h3->GetBinLowEdge(i)+h3->GetBinWidth(i)/2<<"\t"<<h3->GetBinContent(i)<<"\n";
// 			// myfile << Form("%g\t%g\n",h3->GetBinLowEdge(i)+h3->GetBinWidth(i)/2,h3->GetBinContent(i));
// 	    }
// 		myfile.close();
// 	}
// 	else cout << "something bad\n";
//
//
//
//
// 	c0->Clear();
// 	delete h3;
// 	fyield->Close();
// 	delete c0;
// }
//
//
//
// /*==============================MAIN=========================================*/
// void spectrasToASCII(){
//
// 	// Suppress certain root print statements
// 	gErrorIgnoreLevel = kWarning;
//
// 	// When running on CRC
// 	const char *path = "/afs/crc.nd.edu/user/s/saguilar/Group/24Mg_ap/";
// 	chdir(path);
//
// 	const char *detect;
// 	const char *files;
//
// 	// Background Spectrum file
// 	const char *fileBackground = "calibration/background/run0422.root";
//
// 	// Loop through runs:
// 	cout << "\nBEGINNING PEAK FITTING:" << endl;
//
// 	// Record start time
// 	auto start = std::chrono::high_resolution_clock::now();
//
// 	int fileNum = 1;
// 	int runStart = 159 ;
// 	int upToRun;
//
// 	if (loc==1) upToRun = 160;
// 	else upToRun = 410;
//
// 	for(int i=runStart;i<upToRun;i++){
//
// 		// Skip bad runs
// 		if(i>=163 && i<=166) continue;
// 		else if(i==182) continue;
// 		else if(i==244) continue;
// 		else if(i==164) continue;
// 		else if(i==166) continue;
// 		else if(i==182) continue;
// 		else if(i>=244 && i<=255) continue;
// 		else if(i==276) continue;
// 		else if(i==277) continue;
// 		else if(i==285) continue;
// 		else if(i==289) continue;
// 		else if(i==290) continue;
// 		else if(i==294) continue;
// 		else if(i==255) continue;
//
// 		TString runNum_TString;
//
// 		if(i < 100) runNum_TString = "00";
// 		else if((i >= 100) && (i < 1000)) runNum_TString = "0";
// 		else runNum_TString = "";
//
// 		runNum_TString += i;	// Should be format 0001 -> 9999
// 		const char *runNum_String = (const char*)runNum_TString;
//
//
// 		// // Save TFitResult results.
// 		TFile *ff = new TFile(Form("E_cal_spectras/run_%s.root",runNum_String),"RECREATE");
//
// 		// Loop through detectors on board
// 		for(int j=0;j<13;j++){ // 13
//
// 			if (loc==1) files = Form("runs/run%s.root",runNum_String);
// 			else files = Form("/afs/crc.nd.edu/group/nsl/activetarget/data/24Mg_alpha_gamma/spectra/run0%d.root",i);
//
// 			if (j<8){
// 				detect = Form("h0-%d",j);
// 			}
// 			else{
// 				detect = Form("h1-%d",j-8);
// 			}
//
// 			// Perform peak fitting
// 			peakFitter(ff,files,fileBackground,detect,j);
// 		}
//
// 		// // Write all TObjects in memory (TFitResult) to TFile
// 		ff->Write();
// 		ff->Close();
//
// 		cout << Form("Fitting run_%s complete",runNum_String) << endl;
// 		fileNum+=1;
// 		delete ff;
// 	}
//
// 	// Record end time
// 	auto finish = std::chrono::high_resolution_clock::now();
//
// 	std::chrono::duration<double> elapsed = finish - start;
//
// 	cout << "Elapsed time: " << elapsed.count() << " s\n";
// }
