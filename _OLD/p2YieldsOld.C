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

//vector global for a background histogram
TH1D *HBACK;

// global array for new peak ranges
double _peak1SpecLow[13];
double _peak1SpecHigh[13];
double _peak2SpecLow[13];
double _peak2SpecHigh[13];


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


	// Get integrated charge information
	TH1D *qcharge = static_cast<TH1D*>(fyield->Get("h1-7"));
	double charge = qcharge->GetEntries();
	gStyle->SetOptFit(1111);


	// Create the canvas
	TCanvas *c0 = new TCanvas("c0","c0",600,800);
	c0->Divide(1,3); // 2
	c0->Update();
	c0->cd(1);

	// Prepare bg subtracted histogram
	TH1D *ysubtracted = static_cast<TH1D*>(hyield->Clone("ysubtracted"));


	// Only gain match longer runs where background peaks appear
	if(runTime0 > 500){

		/////////////////////////
		///    GAIN MATCH     ///
		/////////////////////////

		HBACK->Scale(scale0);	// Scale the background
		HBACK->SetDirectory(0);

		ifstream backfileREAD;
		backfileREAD.open("Yields/P2/_backgroundPeakRanges.csv",std::ios::in);
		int bgPeaks[4][13];
		// Read in background peak ranges	\
			-1st and 2nd row are low and high of peak1 \
			-3rd and 4th row are low and high of peak2
		string num = "";
		int row = 0;
		int column = 0;
		while (!backfileREAD.eof()) {
			char curr = backfileREAD.get();

			if (curr == ',') {
				bgPeaks[row][column] = std::stoi(num);
				// Reset
				num = "";
				column++;
				if (column > 12) {
					row++;
					column = 0;
				}
			}
			else {
				num += curr;
			}
		}
		backfileREAD.close();

		// Find BG peak positions using same procedure for finding in spectrum
		vector < double > peak1Back;
		vector < double > peak2Back;
		peak1Back = iterative_double_gauss_peak(bgPeaks[0][detectorLoop],bgPeaks[1][detectorLoop],HBACK);
		peak2Back = iterative_single_gauss_peak(bgPeaks[2][detectorLoop],bgPeaks[3][detectorLoop],HBACK);


		// Read in estimated peak ranges from previous run
		ifstream rangefileREAD;
		rangefileREAD.open("Yields/P2/_ranges/peakRanges.csv",std::ios::in);
		int runPeaks[4][13];
		// Read in background peak ranges for runs\
			-1st and 2nd row are low and high of peak1 \
			-3rd and 4th row are low and high of peak2
		num = "";
		row = 0;
		column = 0;
		while (!rangefileREAD.eof()) {
			char curr = rangefileREAD.get();

			if (curr == ',') {
				runPeaks[row][column] = std::stoi(num);
				// Reset
				num = "";
				column++;
				if (column > 12) {
					row++;
					column = 0;
				}
			}
			else {
				num += curr;
			}
		}
		rangefileREAD.close();

		// Find the peak positions in runs given an estimated range
		vector < double > peak1Position;
		vector < double > peak2Position;
		peak1Position = iterative_double_gauss_peak(runPeaks[0][detectorLoop],runPeaks[1][detectorLoop],hyield);
		peak2Position = iterative_single_gauss_peak(runPeaks[2][detectorLoop],runPeaks[3][detectorLoop],hyield);


		// Add new ranges for each detector
		_peak1SpecLow[detectorLoop] = peak1Position[1];
		_peak1SpecHigh[detectorLoop] = peak1Position[2];
		_peak2SpecLow[detectorLoop] = peak2Position[1];
		_peak2SpecHigh[detectorLoop] = peak2Position[2];


		// New gain matched BG histogram
		TH1D *h2 = new TH1D("h2","h2",8192,0,8192);

		h2 = gain_match(peak1Position[0],peak2Position[0],peak1Back[0],peak2Back[0],h2,HBACK);

		// Draw results
		hyield->Draw();
		hyield->GetXaxis()->SetRangeUser(1200,1600);
		// hyield->GetXaxis()->SetRangeUser(0,2200);
		HBACK->Draw("SAME");			// Not gain matched BG
		HBACK->SetLineColor(kOrange);
		h2->Draw("SAME");				// Gain matched BG
		h2->SetLineColor(kRed);
		// gPad->SetLogy();



		c0->cd(2);

		// Recalculate errors manually
		for(int i = 0; i<=ysubtracted->GetNbinsX(); i++){
			double yval = ysubtracted->GetBinContent(i);
			double yval2 = h2->GetBinContent(i);	//fback->Eval(ysubtracted->GetBinCenter(i));
			double yerr = ysubtracted->GetBinError(i);
			double yerr2 = h2->GetBinError(i);	// Takes the error from the gain matchedbg spectrum

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

		ysubtracted->Draw();
		ysubtracted->GetXaxis()->SetRangeUser(1200,1600);
		gPad->SetLogy();
	}
	else{
		cout << Form("Runtime too short: %f seconds. Not performing gain match",runTime0);

		// Draw results
		ysubtracted->Draw();
		ysubtracted->GetXaxis()->SetRangeUser(700,2000);
		gPad->SetLogy();
	}


	c0->cd(3);

	// At this point, yusbtracted histogram has or has not been BG subtracted \
		now we fit the peak of interest:									  \
			P2 : single gaussian

	ysubtracted->GetXaxis()->SetRangeUser(0,2000);
	ysubtracted->Draw();
	// gPad->SetLogy();
	vector < double > p2Peak;
	p2Peak = iterative_single_gauss_area(low,high,ysubtracted);

	ysubtracted->GetXaxis()->SetRangeUser(600,1300);

	double area = p2Peak[0];
	double area_err = p2Peak[1];
	double chi2NDF = p2Peak[2];

	double yield = area/charge;
	double yield_err = area_err/charge;


	double goodFit;
	if (chi2NDF <= 1.25 && chi2NDF >=.75 ) goodFit = 0;
	else goodFit = 1;

	string runNum = fileName;
	runNum = runNum.substr(4,3);
	// runNum = runNum.substr(73,3);

	string detNum = detector;

	c0->SaveAs(Form("Yields/P2/run0%s/det_%s_Fit.png",runNum.c_str(),detNum.c_str()));


	ofstream myfile;
	myfile.open ("Yields/P2/_P2.csv",std::ios::app);
	myfile<<Form("run0%s",runNum.c_str())<<","<< Form("det_%s",detNum.c_str())<<
				","<<yield<<","<<yield_err<<","<<area<<","<<area_err<<","<<goodFit<<"\n";
	myfile.close();


	c0->Clear();
	fyield->Close();
	fbackground->Close();

	delete c0;
	gROOT->Reset();

}



/*==============================MAIN=========================================*/
void p2Yields(){

	// Suppress certain root print statements
	gErrorIgnoreLevel = kWarning;

	// // When running on CRC
	// const char *path = "/afs/crc.nd.edu/user/s/saguilar/Group/24Mg_ap/";
	// chdir(path);

	// Create file directory for output incase its not made
	try {
		cout << "\nATTEMPTING TO CREATE OUTPUT FILE DIRECTORIES" << endl;
		gSystem->Exec(Form("mkdir Yields"));
		gSystem->Exec(Form("mkdir Yields/P2"));
		gSystem->Exec(Form("mkdir Yields/P2/_ranges"));
	}catch(...){}

	// Background Spectrum file
	const char *fileBackground = "calibration/background/run0422.root";

	const char *detect;
	const char *files;

	// Prepare structure of data output in CSV file
	ofstream myfile;
	myfile.open ("Yields/P2/_P2.csv",std::ios::out);
	myfile<<"Run"<<","<<"Detector"<<","<<"Yield"<<","<<"Yield err"<<","<<"Area"<<","<<"Area err"<<","<<"Fit Status"<<"\n";
	myfile.close();


	// BG spectra peak ranges
	ofstream backfileIN;
	backfileIN.open ("Yields/P2/_backgroundPeakRanges.csv",std::ios::out);
	int peak1BackLow[] = {1270,1380,1385,1370,1330,1390,1380,1380,1370,1390,1390,1350,1380};
	int peak1BackHigh[] = {1380,1500,1510,1485,1450,1540,1510,1520,1490,1510,1530,1500,1490};
	int peak2BackLow[] = {1750,1927,1935,1930,1887,1980,1942,1946,1912,1948,1980,1942,1942};
	int peak2BackHigh[] = {1830,2015,2025,2010,1972,2064,2017,2050,1993,2029,2073,2027,2020};
	for(int ii = 0;ii<13;ii++){
		backfileIN<<peak1BackLow[ii]<<",";
	}backfileIN<<endl;
	for(int ii = 0;ii<13;ii++){
		backfileIN<<peak1BackHigh[ii]<<",";
	}backfileIN<<endl;
	for(int ii = 0;ii<13;ii++){
		backfileIN<<peak2BackLow[ii]<<",";
	}backfileIN<<endl;
	for(int ii = 0;ii<13;ii++){
		backfileIN<<peak2BackHigh[ii]<<",";
	}backfileIN<<endl;
	backfileIN.close();

	// Run spectra peak ranges from run0159 \
		These get updated as runs progress
	ofstream rangefileIN;
	rangefileIN.open ("Yields/P2/_ranges/peakRanges.csv",std::ios::out);
	int peak1SpecLow[] = {1285,1400,1400,1390,1360,1430,1410,1420,1410,1430,1420,1410,1410};
	int peak1SpecHigh[] = {1360,1480,1485,1480,1440,1510,1480,1500,1480,1510,1510,1490,1490};
	int peak2SpecLow[] = {1740,1885,1923,1890,1841,1934,1895,1913,1890,1924,1936,1895,1898};
	int peak2SpecHigh[] = {1840,2034,2048,2036,2007,2090,2037,2065,2018,2059,2085,2097,2035};
	for(int ii = 0;ii<13;ii++){
		rangefileIN<<peak1SpecLow[ii]<<",";
	}rangefileIN<<endl;
	for(int ii = 0;ii<13;ii++){
		rangefileIN<<peak1SpecHigh[ii]<<",";
	}rangefileIN<<endl;
	for(int ii = 0;ii<13;ii++){
		rangefileIN<<peak2SpecLow[ii]<<",";
	}rangefileIN<<endl;
	for(int ii = 0;ii<13;ii++){
		rangefileIN<<peak2SpecHigh[ii]<<",";
	}rangefileIN<<endl;
	rangefileIN.close();


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
			gSystem->Exec(Form("mkdir Yields/P2/run0%d",i));
		}catch(...){}


		// Loop through detectors on board 1 (0-7) and board 2 (8-12)
		for(int j=0;j<13;j++){

			files = Form("run0%d.root",i);
			// files = Form("/afs/crc.nd.edu/group/nsl/activetarget/data/24Mg_alpha_gamma/spectra/run0%d.root",i);
			if (j<8){
				detect = Form("h0-%d",j);
			}
			else{
				detect = Form("h1-%d",j-8);
			}

			// Estimated peak positions from calibration
			int peakPos[] = {927,1002,1014,999,973,1022,1009,1014,1000,1011,1020,996,1000};
			p2 = peakPos[j];


			// Perform peak fitting
			peakFitter(files,fileBackground,detect,p2-50,p2+50,j);
		}

		// Write out new spectra peak ranges to gain match to
		ofstream rangefileOUT;
		rangefileOUT.open ("Yields/P2/_ranges/peakRanges.csv",std::ios::out);
		for(int ii = 0;ii<13;ii++){
			rangefileOUT<<_peak1SpecLow[ii]<<",";
		}rangefileOUT<<endl;
		for(int ii = 0;ii<13;ii++){
			rangefileOUT<<_peak1SpecHigh[ii]<<",";
		}rangefileOUT<<endl;
		for(int ii = 0;ii<13;ii++){
			rangefileOUT<<_peak2SpecLow[ii]<<",";
		}rangefileOUT<<endl;
		for(int ii = 0;ii<13;ii++){
			rangefileOUT<<_peak2SpecHigh[ii]<<",";
		}rangefileOUT<<endl;
		rangefileOUT.close();


		cout << Form("Fitting %d complete",fileNum) << endl;
		fileNum+=1;
	}
}
