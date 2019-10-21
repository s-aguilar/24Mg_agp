#include <iostream>
using std::cout;
using std::endl;

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
// #include "/Users/sebastian/Desktop/27Al_pga/calibration/fitFunctions.h"
#include "/Users/sebastian/Desktop/24Mg_agp/calibration/fitFunctions.h"

// Gain match routine
#include "/Users/sebastian/Desktop/24Mg_agp/calibration/gainMatch.h"


// Flourine line from first run, these get updated run by run with the new ranges
double fLow[] = {113,123,129,125,124,129,126,127,126,128,127,126,123};
double fHigh[] = {140,143,144,144,141,147,146,146,143,143,147,143,147};


// if 1, use local, if 0 use CRC
int loc = 1;

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
	c0->SetRightMargin(0.09);
	c0->SetLeftMargin(0.13);
	c0->SetBottomMargin(0.15);

	// c0->cd(1);

	// // Prepare histograms and name them
	// TH1D *h1 = new TH1D(Form("h1 - det%i",detLoop),Form("h1 - det%i",detLoop),4096,0,4096);
	// TH1D *h2 = new TH1D(Form("h2 - det%i",detLoop),Form("h2 - ECAL - det%i",detLoop),4096,0,4096);

	double area;
	double area_err;
	double chi2NDF;
	double sig1;
	double sig2;
	bool isValid;
	int status;

	// Draw results
	// hyield->Draw();
	// hyield->GetXaxis()->SetRangeUser(2,2500);
	// hyield->SetStats(kFALSE);
	// gPad->SetLogy();
	// c0->cd(2);
	// hyield->SetBinContent(0,0);
	// hyield->GetXaxis()->SetRangeUser(0,2500);
	//
	// gPad->SetLogy();
	// hyield->GetYaxis()->SetTitle("Counts/Channel");
	// hyield->Draw();
	// hyield->SetStats(kFALSE);
	// hyield->GetXaxis()->SetTitle("Channel");
	// hyield->GetYaxis()->SetTitle("Counts/Channel");
	// hyield->GetXaxis()->CenterTitle();
	// hyield->GetYaxis()->CenterTitle();
	// hyield->SetTitle("");

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

	// h3->SetBinContent(0,0);
	h3->GetXaxis()->SetRangeUser(0,2500);

	gPad->SetLogy();
	// h3->GetYaxis()->SetTitle("Counts/bin");
	h3->Draw();
	h3->SetStats(kFALSE);
	h3->GetXaxis()->SetTitle("Energy (keV)");
	h3->GetXaxis()->SetTitleSize(0.06);
	h3->GetXaxis()->SetLabelSize(0.06);
	h3->GetYaxis()->SetTitle("Counts / Bin");
	h3->GetYaxis()->SetTitleSize(0.06);
	h3->GetYaxis()->SetLabelSize(0.06);
	// cout << h3->GetBinWidth(1000)<<endl;  // 1 keV bins
	h3->GetXaxis()->CenterTitle();
	h3->GetYaxis()->CenterTitle();
	h3->SetTitle("");


	// Tlatex (x ,y)
	// TArrow (x1,y1,x2,y2)

	TLatex tt(1830,600000,"E_{#bf{#alpha}} = 5.441 MeV");
	tt.SetTextSize(0.05);
	tt.Draw();

	TLatex t(1450,3200,"#bf{#alpha}_{1} \t1.369 MeV");
	t.SetTextSize(0.05);
	TArrow ar(1600,2200,1379,600,0.01,"|>");
	ar.Draw();
	t.Draw("SAME");

	TLatex t1(550,2300,"p_{1} \t0.844 MeV");
	t1.SetTextSize(0.05);
	TArrow ar1(750,2100,840,400,0.01,"|>");
	ar1.Draw();
	t1.Draw("SAME");

	TLatex t2(900,1100,"p_{2} \t1.015 MeV");
	t2.SetTextSize(0.05);
	TArrow ar2(1015,800,1015,300,0.01,"|>");
	ar2.Draw();
	t2.Draw("SAME");

	h3->SetLineColor(kBlack);
	TGaxis *axis1 = new TGaxis(2500,1,2500,323000,1,323000,10,"+GLB"); //GLB
  	axis1->SetName("axis3");
	axis1->SetLabelSize(0);
  	axis1->Draw();

	TGaxis *axis2 = new TGaxis(0,323000,2500,323000,0,2500,25,"-LB"); //GLB
  	axis2->SetName("axis2");
	axis2->SetLabelSize(0);
  	axis2->Draw();

	string runNum = fileName;
		if (loc==1) runNum = runNum.substr(9,3);
	c0->SaveAs(Form("run_%s_det-%i_Spectra.pdf",runNum.c_str(),detLoop));
	// TFile *aa = new TFile("hist.dat","RECREATE");


	int n = h3->GetNbinsX();
	FILE *fptr = fopen("hist.dat", "w");
   	for (int i=1; i<=n; i++) {
      	fprintf(fptr,"%g\t%g\n",
				h3->GetBinLowEdge(i)+h3->GetBinWidth(i)/2,
             	h3->GetBinContent(i));
   }
   fclose(fptr);

	// h3->Print("all"); > "hist.dat"
	// aa->Write();
	// aa->Close();

	// // ATTEMPT TO FIT
	// try {
	//
	// 	string runNum = fileName;
	// 	if (loc==1) runNum = runNum.substr(9,3);
	// 	else runNum = runNum.substr(73,3);
	//
	//
	// 	// .Get() returns the contained pointer to TFitResult. Dereference it with "*""
	// 	TFitResult fitResults = static_cast<TFitResult>(*single_gauss_area_p1(hyield).Get());
	//
	// 	// Write out fitResults to current TFile (TFitOut)
	// 	fitResults.Write(Form("det-%i",detLoop));
	//
	// 	area = fitResults.Parameter(3);
	// 	area_err = fitResults.Error(3);
	// 	chi2NDF = fitResults.Chi2()/fitResults.Ndf();
	// 	sig1 = fitResults.Parameter(5);
	// 	isValid = fitResults.IsValid();
	// 	status = fitResults.Status();
	//
	//
	// 	// The charge is integrated charge of proton which is 1
	// 	double yield = area/(charge);
	// 	double yield_err = area_err/(charge);
	//
	// 	string detNum = detector;
	//
	// 	// Save time: svg << pdf<< eps << jpg << png
	// 	// c0->SaveAs(Form("Yields/P1/run_%s/%s_Fit.pdf",runNum.c_str(),detNum.c_str()));
	// 	c0->SaveAs(Form("run_%s_det-%i_Spectra.pdf",runNum.c_str(),detLoop));
	//
	// 	// ofstream myfile;
	// 	// myfile.open ("Yields/P1/_P1.csv",std::ios::app);
	// 	// myfile<<Form("run_%s",runNum.c_str())<<","<< Form("%s",detNum.c_str())
	// 	// 	<<","<<yield<<","<<yield_err<<","<<area<<","<<area_err<<","<<sig1
	// 	// 	<<","<<chi2NDF<<","<<isValid<<","<<status<<","<<charge<<"\n";
	// 	// myfile.close();
	// }catch(...){}

	c0->Clear();
	fyield->Close();
	delete c0;
}



/*==============================MAIN=========================================*/
void paperPlotter(){

	// Suppress certain root print statements
	gErrorIgnoreLevel = kWarning;

	// When running on CRC
	// const char *path = "/afs/crc.nd.edu/user/s/saguilar/Group/24Mg_ap/";
	// chdir(path);

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
	else upToRun = 160;

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
		TFile *ff = new TFile(Form("run_%s.root",runNum_String),"RECREATE");

		// Loop through detectors on board
		for(int j=0;j<1;j++){ // 13

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
	}

	// Record end time
	auto finish = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> elapsed = finish - start;

	cout << "Elapsed time: " << elapsed.count() << " s\n";
}
