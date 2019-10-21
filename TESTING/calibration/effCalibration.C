#include <iostream>
using std::cout;
using std::endl;

#include "TROOT.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TLine.h"

// Fitting routines
#include "fitFunctions.h"

// Gain match routine
#include "gainMatch.h"

//vector global for a background histogram
TH1D *HBACK;


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void calibration(const char *fileName,const char *fileBack,const char *detector,
	const vector < double > &peak1BackLow, const vector < double > &peak1BackHigh,
	const vector < double > &peak2BackLow, const vector < double > &peak2BackHigh,
	vector < double > &peak1SpecLow, vector < double > &peak1SpecHigh,
	vector < double > &peak2SpecLow, vector < double > &peak2SpecHigh,
	int low, int high, int fileNameLoop, int detLoop){

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
	TCanvas *c0 = new TCanvas("c0","c0",1920,1080);
	c0->Divide(1,2);
	c0->Update();
	c0->cd(1);


	////////////////////////
	///    GAIN MATCH    ///
	////////////////////////

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


	// New gain matched BG histogram
	TH1D *h2 = new TH1D("h2","h2",8192,0,8192);


	// New gain matched BG histogram -> h2
	vector < double > gain;
	gain = gain_match(peak1Position[0],peak2Position[0],peak1Back[0],peak2Back[0],h2,HBACK);

	// Draw results
	hyield->Draw();
	HBACK->Draw("SAME");
	HBACK->SetLineColor(kGreen); // not gain matched
	h2->SetLineColor(kRed);
	h2->Draw("SAME");			// gain matched
	hyield->SetStats(kFALSE);
	// hyield->GetXaxis()->SetRangeUser(200,1550);

	// gPad->SetLogy();


	if(fileNameLoop==0) hyield->GetXaxis()->SetRangeUser(800,1600);
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
	TFitResultPtr r = ysubtracted->Fit("ffit","SQRN0");  // TFitResultPtr contains the TFitResult      dont draw fit
	// TMatrixDSym cov = r->GetCovarianceMatrix();   //  to access the covariance matrix
	// cout << TMath::Sqrt(cov[3][3]) << " " << r->ParError(3) << endl;


	// Draw the background subtracted histogram
	ysubtracted->GetXaxis()->SetRangeUser(low-20,high+20);
	ysubtracted->Draw();


	// Store fit parameters in an array and Draw
	double par[6];
	double par_err[6];
	ffit->GetParameters(par);
	// cout << par[5] << " ";
	fback->SetParameters(par[0],par[1],par[2]);	// Set parameters on BG from S+B fit
	fback->SetParErrors(ffit->GetParErrors());	// Set parameter errors on BG terms from previous fit
	ysubtracted->GetXaxis()->SetRangeUser(par[4]-3*par[5]-200,par[4]+3*par[5]+200);
	// ffit->Draw("SAME");
	// fback->Draw("SAME");


	// Plot markers showing interval of integration
	TLine *line1 = new TLine(par[4]-3*par[5],ysubtracted->GetMinimum(),par[4]-3*par[5],ysubtracted->GetMaximum());
	TLine *line2 = new TLine(par[4]+3*par[5],ysubtracted->GetMinimum(),par[4]+3*par[5],ysubtracted->GetMaximum());
	TLine *line3 = new TLine(par[4],ysubtracted->GetMinimum(),par[4],1.1*ysubtracted->GetMaximum());
	line1->SetLineColor(2);
	line2->SetLineColor(2);
	line3->SetLineColor(kBlack);
	// line1->Draw("SAME");
	// line2->Draw("SAME");
	// line3->Draw("SAME");
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
	double area2_err = TMath::Sqrt(abs(area2));	// TEMPORARY

	// Peak area and its error (NEW)
	double area = area1 - area2;
	double area_err = TMath::Sqrt(area1_err * area1_err + area2_err * area2_err);

	if(fileNameLoop==0){
		cout << "Old area parameter: "<<area << "\t Old histo int - bg int: " << A << endl;
	}

	if(fileNameLoop==0) {





	// Fit definitions
	TF1 *fitGaus = new TF1("fitGaus",single_gauss_func,par[4]-3*par[5],par[4]+3*par[5],3);
	fitGaus->SetParNames("norm","mean","sigma");
	fitGaus->SetParameters(4000,(low+high)/2,10);
	fitGaus->SetParLimits(2,9,13);
	fitGaus->SetLineColor(kViolet);
	fitGaus->SetNpx(1e5);

	ysubtracted->Fit("fitGaus","SQR");

	// Store fit parameters in an array and Draw
	double par1[3];
	double par1_err[3];
	fitGaus->GetParameters(par1);
	A = par1[0];
	A_err = fitGaus->GetParError(0);

	cout<< "Area of single gauss param: "<<A << endl;


	TF1 *fbackOnly = new TF1("fbackOnly",background_func,par1[1]+3*par1[2],par1[1]+9*par1[2],3);
	fbackOnly->SetParNames("a0","a1","a2");
	fbackOnly->FixParameter(2,0);		// Makes it a linear background
	fbackOnly->SetLineColor(kRed);
	fbackOnly->SetNpx(1e5);

	ysubtracted->Fit("fbackOnly","SQR");
	fitGaus->Draw("SAME");


	// Fit definitions
	TF1 *fit = new TF1("fit",fit_single_gauss_func,par[4]-3*par[5],par[4]+3*par[5],6);
	fit->SetParNames("a0","a1","a2","norm","mean","sigma");
	fit->SetParameters(fbackOnly->GetParameter(0),fbackOnly->GetParameter(1),fbackOnly->GetParameter(2),par1[0],par1[1],par1[2]);
	fit->FixParameter(0,fbackOnly->GetParameter(0));
	fit->FixParameter(1,fbackOnly->GetParameter(1));
	fit->FixParameter(2,fbackOnly->GetParameter(2));
	fit->SetParLimits(5,9,13);
	// fit->SetLineColor(kViolet);
	fit->SetNpx(1e5);
	fit->SetLineColor(kBlue);
	fit->Draw("SAME");


	// fitGaus->SetLineColor(kViolet);

	// // Store fit parameters in an array and Draw
	// double par1[3];
	// double par1_err[3];
	// fitGaus->GetParameters(par1);
	// area = par1[0];

	TLine *line4 = new TLine(par1[1]-3*par1[2],ysubtracted->GetMinimum(),par1[1]-3*par1[2],ysubtracted->GetMaximum());
	TLine *line5 = new TLine(par1[1]+3*par1[2],ysubtracted->GetMinimum(),par1[1]+3*par1[2],ysubtracted->GetMaximum());
	TLine *line6 = new TLine(par1[1],ysubtracted->GetMinimum(),par1[1],1.1*ysubtracted->GetMaximum());
	line4->SetLineColor(2);
	line5->SetLineColor(2);
	line6->SetLineColor(kBlack);

	line4->Draw("SAME");
	line5->Draw("SAME");
	line6->Draw("SAME");

	// Histogram integration
	area1 = ysubtracted->IntegralAndError(par1[1]-3*par1[2],par1[1]+3*par1[2],area1_err,"width");

	// BG fit function integration
	area2 = fbackOnly->Integral(par1[1]+3*par1[2],par1[1]+9*par1[2]);
	area2_err = TMath::Sqrt(abs(area2));	// TEMPORARY

	// Peak area and its error (NEW)
	double area = area1 - area2;
	double area_err = TMath::Sqrt(area1_err * area1_err + area2_err * area2_err);

	cout <<"histo int - BG int of new: "<<area <<"\n"<<endl;

	// cout << A << "\t" << area << endl;
	// cout << A_err << "\t" << area_err << "\n" <<endl;

	////////////////////
	// area = A;
	// area_err = A_err;
	////////////////////

	}

	ysubtracted->SetStats(kFALSE);

	// Check effects of new area method vs old
	double percentChange;
	percentChange = area/A*100-100;


	// Save fits and calibration information
	if (fileNameLoop==0) {
		c0->SaveAs(Form("calPlots/60Co_1173peak/det_%s_Fit.png",detNum.c_str()));

		ofstream myfile;
		myfile.open ("csv/60Co_1173cal.csv",std::ios::app);
		myfile<< Form("det_%s",detNum.c_str())<<","<<mean<<","<<sigma<<","<<area
				<<","<<area_err<<","<<runTime0<<","<<A<<","<<A_err<<","<<percentChange<<"\n";
		myfile.close();
	}
	else if (fileNameLoop==1) {
		c0->SaveAs(Form("calPlots/60Co_1332peak/det_%s_Fit.png",detNum.c_str()));

		ofstream myfile;
		myfile.open ("csv/60Co_1332cal.csv",std::ios::app);
		myfile<< Form("det_%s",detNum.c_str())<<","<<mean<<","<<sigma<<","<<area
				<<","<<area_err<<","<<runTime0<<","<<A<<","<<A_err<<","<<percentChange<<"\n";
		myfile.close();
	}
	else if (fileNameLoop==2) {
		c0->SaveAs(Form("calPlots/137Cs_661peak/det_%s_Fit.png",detNum.c_str()));

		ofstream myfile;
		myfile.open ("csv/137Cs_661cal.csv",std::ios::app);
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

	const char *fileLoc[] = {"60Co/15N_run/run0417.root","60Co/15N_run/run0417.root","137Cs/15N_run/run0416.root"};
	// const char *fileLoc[] = {"60Co/24Mg_run/run0419.root","60Co/24Mg_run/run0419.root","137Cs/24Mg_run/run0420.root"};
	const char *detect[] = {"h0-0","h0-1","h0-2","h0-3","h0-4","h0-5","h0-6","h0-7","h1-0","h1-1","h1-2","h1-3","h1-4"};
	const char *fileBackground = "background/run0422.root";
	int *low, *high;


	// BG spectra peak ranges
	const vector < double >  peak1BackLow ({310,330,330,330,320,340,330,330,330,330,320,310,330});
	const vector < double >  peak1BackHigh ({350,370,380,370,360,380,380,380,380,380,380,375,380});
	const vector < double >  peak2BackLow ({1270,1380,1385,1370,1330,1390,1380,1380,1370,1390,1400,1350,1380});
	const vector < double >  peak2BackHigh ({1380,1500,1510,1485,1450,1540,1510,1520,1490,1510,1530,1520,1500});

	// BG peak ranges in calibration runs
	vector < double >  peak1SpecLow;
	vector < double >  peak1SpecHigh;
	vector < double >  peak2SpecLow;
	vector < double >  peak2SpecHigh;

	peak1SpecLow={300,300,320,310,310,325,320,310,320,320,320,312,310};
	peak1SpecHigh={350,400,400,380,380,385,380,395,385,390,375,375,380};
	peak2SpecLow={1250,1360,1370,1355,1320,1385,1375,1360,1385,1385,1355,1360,1380};
	peak2SpecHigh={1400,1530,1530,1500,1480,1540,1520,1510,1540,1540,1520,1520,1500};

	// Outer loop: Runs over the calibration spectra
	cout << "\nBEGINNING CALIBRATION FITTING:" << endl;
	for (int i=0;i<=2;i++){

		// Estimate of the peak ranges [low,high]
		int low1173[] = {1040,1125,1136,1115,1090,1140,1130,1135,1120,1135,1140,1110,1120};
		int high1173[] = {1110,1200,1210,1200,1170,1230,1220,1220,1200,1220,1220,1200,1200};

		int low1332[] = {1170,1260,1275,1250,1220,1280,1270,1270,1260,1280,1270,1250,1260};
		int high1332[] = {1250,1360,1375,1360,1320,1400,1370,1380,1360,1370,1410,1360,1360};

		int low661[] = {580,623,637,625,600,630,630,630,610,620,630,610,620};
		int high661[] = {660,700,720,700,685,720,710,715,700,700,700,695,700};

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

			calibration(fileLoc[i],fileBackground,detect[j],peak1BackLow,peak1BackHigh,
				peak2BackLow,peak2BackHigh,peak1SpecLow,peak1SpecHigh,
				peak2SpecLow,peak2SpecHigh,low[j],high[j],i,j);
		}
		cout << Form("Files %i/3 complete",i+1) << endl;

	}

	cout << "\nDONE!" << endl;
}
