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

// CSV I/O routine
#include "calibration/csvIO.h"

// global histogram for background
TH1D *HBACK;

// global array for new peak ranges
double _peak1SpecLow[13];
double _peak1SpecHigh[13];
double _peak2SpecLow[13];
double _peak2SpecHigh[13];

// global array for gain match constants
double _a_gain[13];
double _b_gain[13];

// global array for calibration constants
double _a_calibrator[13];
double _b_calibrator[13];


int peakFitter(const char *fileName,const char *fileBack,const char *detector,double low,double high,int detectorLoop){

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


	// Get integrated charge information
	TH1D *qcharge = static_cast<TH1D*>(fyield->Get("h1-7"));
	double charge = qcharge->GetEntries();
	gStyle->SetOptFit(1111);


	// Create the canvas
	TCanvas *c0 = new TCanvas("c0","c0",600,800);
	c0->Divide(1,3); // 2
	c0->Update();

	c0->cd(1);

	// Scale the background
	HBACK->Scale(scale0);
	HBACK->SetDirectory(0);

	// Prepare bg subtracted histogram
	TH1D *ysubtracted = static_cast<TH1D*>(hyield->Clone("ysubtracted"));
	TH1D *h2 = new TH1D("h2","h2",8192,0,8192);
	TH1D *h3 = new TH1D("h3","h3",8192,0,8192);

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

	vector < double > a1Peak;


	// Only gain match longer runs where background peaks appear
	if(runTime0 > 500){

		/////////////////////////
		///    GAIN MATCH     ///
		/////////////////////////

		ifstream backfileREAD;
		backfileREAD.open("Yields/A1/_backgroundPeakRanges.csv",std::ios::in);
		double bgPeaks[4][13];
		// Read in background peak ranges	\
			-1st and 2nd row are low and high of peak1 \
			-3rd and 4th row are low and high of peak2
		readIn(backfileREAD,bgPeaks);
		backfileREAD.close();

		// Find BG peak positions using same procedure for finding in spectrum
		vector < double > peak1Back;
		vector < double > peak2Back;
		HBACK->Draw();
		peak1Back = single_gauss_peak(bgPeaks[0][detectorLoop],bgPeaks[1][detectorLoop],HBACK);
		peak2Back = iterative_double_gauss_peak(bgPeaks[2][detectorLoop],bgPeaks[3][detectorLoop],HBACK);
		HBACK->GetXaxis()->SetRangeUser(1200,1500);
		HBACK->SetStats(kFALSE);
		c0->cd(2);

		// Read in estimated peak ranges from previous run
		ifstream rangefileREAD;
		rangefileREAD.open("Yields/A1/_ranges/peakRanges.csv",std::ios::in);
		double runPeaks[4][13];
		// Read in background peak ranges for runs\
			-1st and 2nd row are low and high of peak1 \
			-3rd and 4th row are low and high of peak2
		readIn(rangefileREAD,runPeaks);
		rangefileREAD.close();

		// Find the peak positions in runs given an estimated range
		vector < double > peak1Position;
		vector < double > peak2Position;


		// hyield->Draw(); ///////////////////////////////
		peak1Position = iterative_double_gauss_peak(runPeaks[0][detectorLoop],runPeaks[1][detectorLoop],hyield);
		peak2Position = iterative_single_gauss_peak(runPeaks[2][detectorLoop],runPeaks[3][detectorLoop],hyield);


		// Add new ranges for each detector
		_peak1SpecLow[detectorLoop] = peak1Position[1];
		_peak1SpecHigh[detectorLoop] = peak1Position[2];
		_peak2SpecLow[detectorLoop] = peak2Position[1];
		_peak2SpecHigh[detectorLoop] = peak2Position[2];


		// New gain matched BG histogram
		vector < double > gain;
		gain = gain_match(peak1Position[0],peak2Position[0],peak1Back[0],peak2Back[0],h2,HBACK);
		a = gain[0];
		b = gain[1];
		a_gm = gain[0];
		b_gm = gain[1];

		_a_gain[detectorLoop] = a_gm;
		_b_gain[detectorLoop] = b_gm;

		// Draw results
		// hyield->GetXaxis()->SetRangeUser(1000,1500);
		// hyield->GetXaxis()->SetRangeUser(50,200);////////////////////////
		// HBACK->Draw("SAME");			// Not gain matched BG
		// HBACK->SetLineColor(kOrange);
		// h2->Draw("SAME");				// Gain matched BG
		// h2->SetLineColor(kRed);
		gPad->SetLogy();
		hyield->SetStats(kFALSE);

		c0->cd(3); //2

		// Recalculate errors manually
		for(int i = 0; i<=ysubtracted->GetNbinsX(); i++){
			double yval = ysubtracted->GetBinContent(i);
			double yval2 = h2->GetBinContent(i);	//fback->Eval(ysubtracted->GetBinCenter(i));
			double yerr = ysubtracted->GetBinError(i);
			double yerr2 = h2->GetBinError(i);	// Takes the error from the gain matched bg spectrum

			ysubtracted->SetBinContent(i,yval-yval2);
			ysubtracted->SetBinError(i,TMath::Sqrt(yerr*yerr+yerr2*yerr2));
		}

		vector < double > calibrators;
		calibrators = calibrate(1468,109.9,peak1Position[0],peak2Position[0],h3,ysubtracted);
		a_cal = calibrators[0];
		b_cal = calibrators[1];

		// cout << a_cal << "\t" << b_cal << endl;

		_a_calibrator[detectorLoop] = a_cal;
		_b_calibrator[detectorLoop] = b_cal;

		h3->Draw();
		// h3->GetXaxis()->SetRangeUser(1000,1500);
		h3->GetXaxis()->SetRangeUser(0,1500);
		// gPad->SetLogy();
		h3->SetStats(kFALSE);
	}
	else{
		// cout << Form("Runtime too short: %f seconds. Not performing gain match",runTime0);
		a_cal = _a_calibrator[detectorLoop];
		b_cal = _b_calibrator[detectorLoop];
		a_gm = _a_gain[detectorLoop];
		b_gm = _b_gain[detectorLoop];

		h2 = gain_match2(a_gm,b_gm,h2,HBACK);

		// Draw results
		// hyield->Draw();
		hyield->GetXaxis()->SetRangeUser(1000,1500);
		// HBACK->GetXaxis()->SetRangeUser(1000,1500);
		// HBACK->Draw("SAME");			// Not gain matched BG
		// HBACK->SetLineColor(kOrange);
		// h2->Draw("SAME");				// Gain matched BG
		// h2->SetLineColor(kRed);
		gPad->SetLogy();
		hyield->SetStats(kFALSE);

		c0->cd(2);

		// Recalculate errors manually
		for(int i = 0; i<=ysubtracted->GetNbinsX(); i++){
			double yval = ysubtracted->GetBinContent(i);
			double yval2 = h2->GetBinContent(i);	//fback->Eval(ysubtracted->GetBinCenter(i));
			double yerr = ysubtracted->GetBinError(i);
			double yerr2 = h2->GetBinError(i);	// Takes the error from the gain matched bg spectrum

			ysubtracted->SetBinContent(i,yval-yval2);
			ysubtracted->SetBinError(i,TMath::Sqrt(yerr*yerr+yerr2*yerr2));
		}

		h3 = calibrate2(a_cal,b_cal,h3,ysubtracted);
		h3->Draw();
		h3->GetXaxis()->SetRangeUser(1000,1500);
		// gPad->SetLogy();
		h3->SetStats(kFALSE);
	}

	c0->cd(3);

	h3->GetXaxis()->SetRangeUser(1000,1500);
	h3->Draw();
	gPad->SetLogy();
	h3->SetStats(kFALSE);

	a1Peak = single_gauss_area(1368.6,h3); // 1368.63

	area = a1Peak[0];
	area_err = a1Peak[1];
	chi2NDF = a1Peak[2];
	sig1 = a1Peak[3];
	linear = a1Peak[4];
	offset = a1Peak[5];

	// The charge is integrated charge of alpha which is 2+
	double yield = area/(charge);
	double yield_err = area_err/(charge);


	double goodFit;
	if (chi2NDF <= 1.4 && chi2NDF >=.6 ) goodFit = 0;
	else goodFit = 1;

	string runNum = fileName;
	runNum = runNum.substr(9,3);
	// runNum = runNum.substr(73,3);

	string detNum = detector;

	c0->SaveAs(Form("Yields/A1/run0%s/det_%s_Fit.png",runNum.c_str(),detNum.c_str()));


	ofstream myfile;
	myfile.open ("Yields/A1/_A1.csv",std::ios::app);
	myfile<<Form("run0%s",runNum.c_str())<<","<< Form("det_%s",detNum.c_str())<<","<<
				yield<<","<<yield_err<<","<<area<<","<<area_err<<","<<runTime0<<","<<
				goodFit<<","<<a<<","<<b<<","<<sig1<<","<<chi2NDF<<","<<linear<<","<<
				offset<<","<<charge<<"\n";
	myfile.close();


	c0->Clear();
	fyield->Close();
	fbackground->Close();
	delete c0;

	gROOT->Reset();

	// Test for peak area containing minimum number of counts
	int test = 0;
	if (area<1000){
		test = 1;	// FAILS
	}
	return test;
}



/*==============================MAIN=========================================*/
void a1Yields(){

	// Suppress certain root print statements
	gErrorIgnoreLevel = kWarning;

	// When running on CRC
	const char *path = "/afs/crc.nd.edu/user/s/saguilar/Group/24Mg_ap/";
	chdir(path);

	// Create file directory for output incase its not made
	// try {
	// 	cout << "\nATTEMPTING TO CREATE OUTPUT FILE DIRECTORIES" << endl;
	// 	gSystem->Exec(Form("mkdir Yields"));
	// 	gSystem->Exec(Form("mkdir Yields/A1"));
	// 	gSystem->Exec(Form("mkdir Yields/A1/_ranges"));
	// }catch(...){}

	// Background Spectrum file
	const char *fileBackground = "calibration/background/run0422.root";

	const char *detect;
	const char *files;

	// Prepare structure of data output in CSV file
	ofstream myfile;
	myfile.open ("Yields/A1/_A1.csv",std::ios::out);
	myfile<<"Run"<<","<<"Detector"<<","<<"Yield"<<","<<"Yield err"<<","<<"Area"<<","
			<<"Area err"<<","<<"Time"<<","<<"Fit Status"<<","<<"a"<<","<<"b"<<","
			<<"sig1"<<","<<"X2NDF"<<","<<"Linear"<<","<<"Offset"<<","<<"Q_int"<<"\n";
	myfile.close();

	// BG spectra peak ranges
	ofstream backfileOUT;
	backfileOUT.open ("Yields/A1/_backgroundPeakRanges.csv",std::ios::out);
	double peak1BackLow[] = {210,220,225,210,209,226,220,220,218,215,220,210,200};
	double peak1BackHigh[] = {235,255,264,260,250,267,260,260,270,280,260,270,270};
	double peak2BackLow[] = {1270,1380,1385,1370,1330,1390,1380,1380,1370,1390,1390,1350,1380};
	double peak2BackHigh[] = {1380,1500,1510,1485,1450,1540,1510,1520,1490,1510,1530,1500,1490};
	writeOut(backfileOUT,peak1BackLow,peak1BackHigh,peak2BackLow,peak2BackHigh);
	backfileOUT.close();

	// Spectra peak ranges from run0159 \
		These get updated as runs progress
	ofstream rangefileOUT;
	rangefileOUT.open ("Yields/A1/_ranges/peakRanges.csv",std::ios::out);
	double peak1SpecLow[] = {1285,1400,1400,1390,1360,1430,1410,1420,1410,1430,1420,1410,1410};
	double peak1SpecHigh[] = {1360,1480,1485,1480,1440,1510,1480,1500,1480,1510,1510,1490,1490};
	// double peak2SpecLow[] = {1740,1885,1923,1890,1841,1934,1895,1913,1890,1924,1936,1895,1898};
	// double peak2SpecHigh[] = {1840,2034,2048,2036,2007,2090,2037,2065,2018,2059,2085,2097,2035};
	double peak2SpecLow[] = {113,123,129,125,124,129,126,127,126,128,127,126,123};
	double peak2SpecHigh[] = {140,143,144,144,141,147,146,146,143,143,147,143,147};
	writeOut(rangefileOUT,peak1SpecLow,peak1SpecHigh,peak2SpecLow,peak2SpecHigh);
	rangefileOUT.close();

	// // BG spectra peak ranges
	// ofstream backfileOUT;
	// backfileOUT.open ("Yields/A1/_backgroundPeakRanges.csv",std::ios::out);
	// double peak1BackLow[] = {1270,1380,1385,1370,1330,1390,1380,1380,1370,1390,1390,1350,1380};
	// double peak1BackHigh[] = {1380,1500,1510,1485,1450,1540,1510,1520,1490,1510,1530,1500,1490};
	// double peak2BackLow[] = {1750,1927,1935,1930,1887,1980,1942,1946,1912,1948,1980,1942,1942};
	// double peak2BackHigh[] = {1830,2015,2025,2010,1972,2064,2017,2050,1993,2029,2073,2027,2020};
	// writeOut(backfileOUT,peak1BackLow,peak1BackHigh,peak2BackLow,peak2BackHigh);
	// backfileOUT.close();
	//
	// // Spectra peak ranges from run0159 \
	// 	These get updated as runs progress
	// ofstream rangefileOUT;
	// rangefileOUT.open ("Yields/A1/_ranges/peakRanges.csv",std::ios::out);
	// double peak1SpecLow[] = {1285,1400,1400,1390,1360,1430,1410,1420,1410,1430,1420,1410,1410};
	// double peak1SpecHigh[] = {1360,1480,1485,1480,1440,1510,1480,1500,1480,1510,1510,1490,1490};
	// // double peak2SpecLow[] = {1740,1885,1923,1890,1841,1934,1895,1913,1890,1924,1936,1895,1898};
	// // double peak2SpecHigh[] = {1840,2034,2048,2036,2007,2090,2037,2065,2018,2059,2085,2097,2035};
	// double peak2SpecLow[] = {113,123,129,125,124,129,126,127,126,128,127,122,123};
	// double peak2SpecHigh[] = {140,143,144,144,141,147,146,146,143,143,147,143,147};
	// writeOut(rangefileOUT,peak1SpecLow,peak1SpecHigh,peak2SpecLow,peak2SpecHigh);
	// rangefileOUT.close();

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

		double a1;
		int test;

		// Create file directory for output
		// try {
		// 	gSystem->Exec(Form("mkdir Yields/A1/run0%d",i));
		// }catch(...){}


		// Loop through detectors on board 1 (0-7) and board 2 (8-12)
		for(int j=0;j<13;j++){

			files = Form("runs/run0%d.root",i);
			// files = Form("/afs/crc.nd.edu/group/nsl/activetarget/data/24Mg_alpha_gamma/spectra/run0%d.root",i);
			if (j<8){
				detect = Form("h0-%d",j);
			}
			else{
				detect = Form("h1-%d",j-8);
			}

			// Estimated peak positions from calibration
			// double peakPos[] = {1247,1352,1366,1347,1315,1379,1358,1371,1351,1372,1379,1355,1356};
			double peakPos[] = {1247,1352,1366,1347,1315,1379,1358,1371,1351,1372,1379,1355,1356};
			a1 = peakPos[j];

			// Perform peak fitting
			test = peakFitter(files,fileBackground,detect,a1-50,a1+50,j);
		}

		if (test==0){
			// Write out new spectra peak ranges to gain match to
			ofstream rangefile;
			rangefile.open ("Yields/A1/_ranges/peakRanges.csv",std::ios::out);
			writeOut(rangefile,_peak1SpecLow,_peak1SpecHigh,_peak2SpecLow,_peak2SpecHigh);
			rangefile.close();
		}

		cout << Form("Fitting %d complete",fileNum) << endl;
		fileNum+=1;
	}
}
