#include <iostream>
using std::cout;
using std::endl;

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
#include "TLine.h"

// Fitting routines
#include "fitFunctions.h"

// Gain match routine
#include "gainMatch.h"

//vector global for a background histogram
TH1D *HBACK;


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void calibration(const char *fileName,const char *fileBack,const char *detector,int low, int high, int fileNameLoop, int detectorLoop){

	TFile *fyield = new TFile(fileName);
	TFile *fbackground = new TFile(fileBack);    // bg spectra

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


	// Get histograms from root file, prepare bg subtracted histogram
	TH1D *hyield = static_cast<TH1D*>(fyield->Get(detector));
	TH1D *ysubtracted = static_cast<TH1D*>(hyield->Clone("ysubtracted"));
	HBACK = static_cast<TH1D*>(fbackground->Get(detector));

	HBACK->Scale(scale0);		// Scale the background
	HBACK->SetDirectory(0);
	gStyle->SetOptFit(1111);


	// Create the canvas
	TCanvas *c0 = new TCanvas("c0","c0",600,800);
	c0->Divide(1,2);
	c0->Update();
	c0->cd(1);


	////////////////////////
	///    GAIN MATCH    ///
	////////////////////////

	// Background peaks from run0422.root to be gain matched
	int peak1Back[] = {1323,1439,1451,1431,1397,1466,1447,1456,1434,1452,1463,1431,1437};
	int peak2Back[] = {1789,1969,1980,1968,1930,2025,1981,2000,1954,1991,2026,1987,1978};
	// int peak1Back[] = {1332,1450,1460,1440,1406,1466,1456,1464,1444,1462,1470,1431,1438};
	// int peak2Back[] = {1791,1970,1982,1969,1932,2026,1984,2005,1956,1994,2032,1985,1979};

	int peak1Spec[13];
	int peak2Spec[13];

	// Cs137
	if (fileNameLoop==2){
		int p1[] = {1322,1431,1443,1421,1394,1453,1435,1444,1426,1444,1455,1424,1426};
		int p2[] = {1787,1959,1970,1955,1927,2008,1965,1984,1944,1980,2015,1977,1962};
		// int p1[] = {1319,1433,1443,1422,1392,1455,1436,1444,1425,1444,1455,1424,1429};
		// int p2[] = {1784,1959,1970,1955,1927,2008,1965,1984,1944,1980,2015,1977,1962};
		for (int i = 0; i<=12; i++){
			peak1Spec[i] = p1[i];
			peak2Spec[i] = p2[i];
		}
	}
	else{ // Its Co60
		int p1[] = {1321,1436,1446,1426,1390,1460,1440,1450,1430,1451,1460,1429,1433};
		int p2[] = {1785,1967,1975,1960,1922,2018,1973,1992,1952,1987,2017,1429,1969};
		for (int i = 0; i<=12; i++){
			peak1Spec[i] = p1[i];
			peak2Spec[i] = p2[i];
		}
	}

	// New gain matched BG histogram
	TH1D *h2 = new TH1D("h2","h2",8192,0,8192);
	gain_match(peak1Spec[0],peak2Spec[0],peak1Back[0],peak2Back[0],h2,HBACK);

	// Draw results
	hyield->Draw();
	HBACK->Draw("SAME");
	HBACK->SetLineColor(kGreen);
	h2->SetLineColor(kRed);
	h2->Draw("SAME");
	// gPad->SetLogy();


	if(fileNameLoop==0) hyield->GetXaxis()->SetRangeUser(800,1600);
	// if(fileNameLoop==0) hyield->GetXaxis()->SetRangeUser(0,1600);
	else if(fileNameLoop==1) hyield->GetXaxis()->SetRangeUser(1000,1600);
	else if(fileNameLoop==2) hyield->GetXaxis()->SetRangeUser(500,800);


	// Background subtracted spectrum
	c0->cd(2);


    // Recalculate errors manually
	for(int i = 0; i<=ysubtracted->GetNbinsX(); i++){
		double yval = ysubtracted->GetBinContent(i);
		double yval2 = h2->GetBinContent(i);	//fback->Eval(ysubtracted->GetBinCenter(i));
		double yerr = ysubtracted->GetBinError(i);
		double yerr2 = h2->GetBinError(i);	// Takes the error from the bg spectrum

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


	// Fit definitions
	TF1 *ffit = new TF1("ffit",fit_single_gauss_func,low,high,6);
	ffit->SetParNames("a0","a1","a2","norm","mean","sigma");
	ffit->SetLineColor(kGreen);
	ffit->SetNpx(1e5);

	TF1 *fback = new TF1("fback",background_func,low,high,3);
	fback->SetParNames("a0","a1","a2");
	fback->SetLineColor(kCyan);
	fback->SetNpx(1e5);


	// Initial guess of parameters
	ffit->SetParameters(0,0,0,5000,centroid,10);
	ffit->FixParameter(2,0);			// Makes it a linear background (0)*x^2
    ffit->SetParLimits(3,0,1e8);		// Peak Area
	ffit->SetParLimits(4,low,high);		// Peak Centroid
    ffit->SetParLimits(5,2,20);			// Peak Width


	// Perform the fit
	TFitResultPtr r = ysubtracted->Fit("ffit","SQR");  // TFitResultPtr contains the TFitResult
	// TMatrixDSym cov = r->GetCovarianceMatrix();   //  to access the covariance matrix
	// cout << TMath::Sqrt(cov[3][3]) << " " << r->ParError(3) << endl;


	// Draw the background subtracted histogram
	ysubtracted->GetXaxis()->SetRangeUser(low-20,high+20);
	ysubtracted->Draw();


	// Store fit parameters in an array and Draw
	double par[6];
	double par_err[6];
	ffit->GetParameters(par);
	fback->SetParameters(par[0],par[1],par[2]);	// Set parameters on BG from S+B fit
	fback->SetParErrors(ffit->GetParErrors());	// Set parameter errors on BG terms from previous fit
	ysubtracted->GetXaxis()->SetRangeUser(par[4]-3*par[5]-200,par[4]+3*par[5]+200);
	ffit->Draw("SAME");
	fback->Draw("SAME");


	// Plot markers showing interval of integration
	TLine *line1 = new TLine(par[4]-3*par[5],ysubtracted->GetMinimum(),par[4]-3*par[5],ysubtracted->GetMaximum());
	TLine *line2 = new TLine(par[4]+3*par[5],ysubtracted->GetMinimum(),par[4]+3*par[5],ysubtracted->GetMaximum());
	line1->SetFillColor(2);
	line2->SetFillColor(2);
	line1->Draw("SAME");
	line2->Draw("SAME");
	// gPad->SetLogy();


	// Calculate the peak area and error
	string detNum = detector;

	double mean = par[4];
	double sigma = par[5];

	// Peak area using fit parameter (OLD)
	double A = par[3];
	double A_err = ffit->GetParError(3);


	// Histogram integration
	double area1_err;
	double area1 = ysubtracted->IntegralAndError(par[4]-3*par[5],par[4]+3*par[5],area1_err,"width");

	// BG fit function integration
	double area2 = fback->Integral(par[4]-3*par[5],par[4]+3*par[5]);
	double area2_err = TMath::Sqrt(abs(area2));	// THIS IS WRONG, TEMPORARY

	// Peak area and its error (NEW)
	double area = area1 - area2;
	double area_err = TMath::Sqrt(area1_err * area1_err + area2_err * area2_err);


	// Check effects of new area method vs old
	double percentChange;
	percentChange = area/A*100-100;


	// Save fits and calibration information
	if (fileNameLoop==0) {
		c0->SaveAs(Form("calPlots/60Co_1173peak/det_%s_Fit.png",detNum.c_str()));

		ofstream myfile;
		myfile.open ("csv/60Co_1173cal_.csv",std::ios::app);
		myfile<< Form("det_%s",detNum.c_str())<<","<<mean<<","<<sigma<<","<<area
				<<","<<area_err<<","<<runTime0<<","<<A<<","<<A_err<<","<<percentChange<<"\n";
		myfile.close();
	}
	else if (fileNameLoop==1) {
		c0->SaveAs(Form("calPlots/60Co_1332peak/det_%s_Fit.png",detNum.c_str()));

		ofstream myfile;
		myfile.open ("csv/60Co_1332cal_.csv",std::ios::app);
		myfile<< Form("det_%s",detNum.c_str())<<","<<mean<<","<<sigma<<","<<area
				<<","<<area_err<<","<<runTime0<<","<<A<<","<<A_err<<","<<percentChange<<"\n";
		myfile.close();
	}
	else if (fileNameLoop==2) {
		c0->SaveAs(Form("calPlots/137Cs_661peak/det_%s_Fit.png",detNum.c_str()));

		ofstream myfile;
		myfile.open ("csv/137Cs_661cal_.csv",std::ios::app);
		myfile<< Form("det_%s",detNum.c_str())<<","<<mean<<","<<sigma<<","<<area
				<<","<<area_err<<","<<runTime0<<","<<A<<","<<A_err<<","<<percentChange<<"\n";
		myfile.close();
	}

	c0->Clear();
	delete h2;

	fyield->Close();
	fbackground->Close();

	delete fyield;
	delete fbackground;
	delete ffit;
	delete fback;
	delete c0;

	delete HBACK;
	// delete h2; // This breaks code idk why, should be deleted though


	gROOT->Reset();
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


/*============================START OF MAIN==================================*/
void effCalibration(){

	// Suppress certain root print statements
	gErrorIgnoreLevel = kWarning;


	// Create file directory for output incase its not made
	try {
		cout << "\nATTEMPTING TO CREATE OUTPUT FILE DIRECTORIES" << endl;
		gSystem->Exec(Form("mkdir csv"));
		gSystem->Exec(Form("mkdir calPlots"));
		gSystem->Exec(Form("mkdir calPlots/60Co_1173peak"));
		gSystem->Exec(Form("mkdir calPlots/60Co_1332peak"));
		gSystem->Exec(Form("mkdir calPlots/137Cs_661peak"));
	}catch(...){}


	// List of relative root file directories and detector names\
	   These will then be fed into actual calibration function

	// const char *fileLoc[] = {"60Co/15N_run/run0417.root","60Co/15N_run/run0417.root","137Cs/15N_run/run0416.root"};
	const char *fileLoc[] = {"60Co/24Mg_run/run0419.root","60Co/24Mg_run/run0419.root","137Cs/24Mg_run/run0420.root"};
	const char *detect[] = {"h0-0","h0-1","h0-2","h0-3","h0-4","h0-5","h0-6","h0-7","h1-0","h1-1","h1-2","h1-3","h1-4"};
	const char *fileBackground = "background/run0422.root";
	int *low, *high;

	// Outer loop: Runs over the calibration spectra
	cout << "\nBEGINNING CALIBRATION FITTING:" << endl;
	for (int i=0;i<=2;i++){

		// Estimate of the peak ranges [low,high]
		int low1173[] = {1030,1115,1136,1100,1070,1140,1120,1120,1100,1120,1140,1100,1110};
		int high1173[] = {1110,1200,1210,1200,1170,1230,1220,1220,1200,1220,1220,1200,1200};

		int low1332[] = {1170,1260,1275,1250,1220,1280,1270,1270,1260,1280,1270,1260,1260};
		int high1332[] = {1250,1360,1375,1360,1320,1400,1370,1380,1360,1370,1410,1360,1360};

		int low661[] = {570,623,637,625,600,636,638,630,617,620,630,610,620};
		int high661[] = {650,690,700,700,675,705,710,700,690,690,700,695,685};

		// 1173 keV gamma
		if (i == 0){
            low = low1173;
			high = high1173;

			ofstream myfile;
			myfile.open ("csv/60Co_1173cal.csv",std::ios::out);
			myfile<<"Detector"<<","<<"Centroid"<<","<<"Width"<<","<<"Area"<<","
					<<"Area err"<<","<<"Runtime"<<","<<"A"<<','<<"A err"<<","<<"percentChange"<<"\n";
			myfile.close();
		}

		// 1332 keV gamma
		else if (i == 1){
			low = low1332;
			high = high1332;

			ofstream myfile;
			myfile.open ("csv/60Co_1332cal.csv",std::ios::out);
			myfile<<"Detector"<<","<<"Centroid"<<","<<"Width"<<","<<"Area"<<","
					<<"Area err"<<","<<"Runtime"<<","<<"A"<<','<<"A err"<<","<<"percentChange"<<"\n";
			myfile.close();
		}

		// 661 keV gamma
		else if (i == 2){
			low = low661;
			high = high661;

			ofstream myfile;
			myfile.open ("csv/137Cs_661cal.csv",std::ios::out);
			myfile<<"Detector"<<","<<"Centroid"<<","<<"Width"<<","<<"Area"<<","
					<<"Area err"<<","<<"Runtime"<<","<<"A"<<','<<"A err"<<","<<"percentChange"<<"\n";
			myfile.close();
		}

		// Inner loop: Runs over all 13 detectors
		for (int j=0;j<=12;j++){

			calibration(fileLoc[i],fileBackground,detect[j],low[j],high[j],i,j);
		}
		cout << Form("Files %i/3 complete",i+1) << endl;

	}

	cout << "\nDONE!" << endl;
}
