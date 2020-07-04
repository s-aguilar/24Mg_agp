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
	vector < double > &peak2SpecLow, vector < double > &peak2SpecHigh,
	double low, double high, int detLoop){


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

	vector < double > O17_1Peak;


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

		// Found the lines, update their ranges run by run
		fLow[detLoop] = fPosition[1];
		fHigh[detLoop] = fPosition[2];

		vector < double > calibrators;

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
	}


	h3->GetXaxis()->SetRangeUser(600,1200);
	h3->Draw();
	h3->SetStats(kFALSE);

	O17_1Peak = double_gauss_area_p1(843.8,870.7,h3,detLoop); // 843.76, 867.7

	area = O17_1Peak[6];
	area_err = O17_1Peak[7];
	chi2NDF = O17_1Peak[2];
	sig1 = O17_1Peak[3];
	linear = O17_1Peak[4];
	offset = O17_1Peak[5];

	// The charge is integrated charge of alpha which is 2+
	double yield = area/(charge);
	double yield_err = area_err/(charge);

	double goodFit;
	if (chi2NDF <= 1.4 && chi2NDF >=.6 ) goodFit = 1;
	else goodFit = 0;

	string runNum = fileName;
	if (loc==1) runNum = runNum.substr(9,3);
	else runNum = runNum.substr(73,3);

	string detNum = detector;

  // nt n = h3->GetNbinsX();
  // FILE *fptr = fopen(Form("histSubCal_0%s_%s.dat",runNum.c_str(),detNum.c_str()), "w");
  // or (int i=1; i<=n; i++) {
  // fprintf(fptr,"%g\t%g\n",
  // 		h3->GetBinLowEdge(i)+h3->GetBinWidth(i)/2,
  // 		h3->GetBinContent(i));
  //  }
  //  fclose(fptr);
  //
  //  h3 = calibrate2(a_cal,b_cal,h3,hyield);
  //
  //  n = h3->GetNbinsX();
  //  fptr = fopen(Form("histRawCal_0%s_%s.dat",runNum.c_str(),detNum.c_str()), "w");
  //  for (int i=1; i<=n; i++) {
  //   fprintf(fptr,"%g\t%g\n",
  // 	   h3->GetBinLowEdge(i)+h3->GetBinWidth(i)/2,
  // 	   h3->GetBinContent(i));
  // }
  // fclose(fptr);

	c0->SaveAs(Form("Yields/O17_1/run0%s/det_%s_Fit.png",runNum.c_str(),detNum.c_str()));
	c0->SaveAs(Form("Yields/O17_1/det-%i/run0%s_Fit.png",detLoop,runNum.c_str()));


	ofstream myfile;
	myfile.open ("Yields/O17_1/_O17_1.csv",std::ios::app);
	myfile<<Form("run0%s",runNum.c_str())<<","<< Form("det_%s",detNum.c_str())<<","<<
				yield<<","<<yield_err<<","<<area<<","<<area_err<<","<<runTime0<<","<<
				goodFit<<","<<a<<","<<b<<","<<sig1<<","<<chi2NDF<<","<<linear<<","<<
				offset<<","<<charge<<"\n";
	myfile.close();


	c0->Clear();
	delete ysubtracted;
	delete h2;
	delete h3;
	fyield->Close();
	fbackground->Close();
	delete c0;


	gROOT->Reset();
}



/*==============================MAIN=========================================*/
void O17_1Yields(){

	// Suppress certain root print statements
	gErrorIgnoreLevel = kWarning;

	// When running on CRC
	const char *path = "/afs/crc.nd.edu/user/s/saguilar/Group/24Mg_ap/";
	chdir(path);

	// Background Spectrum file
	const char *fileBackground = "calibration/background/run0424.root"; //422 short BG run

	const char *detect;
	const char *files;

	// // Prepare structure of data output in CSV file
	ofstream myfile;
	myfile.open ("Yields/O17_1/_O17_1.csv",std::ios::out);
	myfile<<"Run"<<","<<"Detector"<<","<<"Yield"<<","<<"Yield err"<<","<<"Area"<<","
			<<"Area err"<<","<<"Time"<<","<<"Fit Status"<<","<<"a"<<","<<"b"<<","
			<<"sig1"<<","<<"X2NDF"<<","<<"Linear"<<","<<"Offset"<<","<<"Q_int"<<"\n";
	myfile.close();

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


	// Make directory to visually inspect the fits
	for(int ii = 0; ii<13; ii++){
		try {
			gSystem->Exec(Form("mkdir Yields/O17_1/det-%i",ii));
		}catch(...){}
	}


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

		try {
			gSystem->Exec(Form("mkdir Yields/O17_1/run0%d",i));
		}catch(...){}

		double O17_1;

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

			// Estimated peak positions from calibration
			double peakPos[] = {778.4,836.3,850.7,838.1,815.1,859.8,844.8,851.2,839.2,849.2,850.3,837.7,839.2};
			O17_1 = peakPos[j];

			// Perform peak fitting
			peakFitter(files,fileBackground,detect,peak1BackLow,peak1BackHigh,
				peak2BackLow,peak2BackHigh,peak1SpecLow,peak1SpecHigh,
				peak2SpecLow,peak2SpecHigh,O17_1-40,O17_1+90,j);
		}

		cout << Form("Fitting %d complete",fileNum) << endl;
		fileNum+=1;
	}
}
