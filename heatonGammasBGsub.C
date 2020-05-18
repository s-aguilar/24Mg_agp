// From bg subtracted summed spectras file, calculate yield for Heaton gamma rays
// for the top 100keV energy runs and compare with the all data summed spectra

#include "calibration/fitFunctions.h"

void heatonGammasBGsub(){


    TFile *sumBGsubSpec = new TFile("E_cal_spectras/summedSpectrasBGsub.root","READ");
    TFile *sumBGsubSpecALL = new TFile("E_cal_spectras/summedSpectrasBGsubALL.root","READ");


    TH1D *sumBGsubSpecDet = static_cast<TH1D*>(sumBGsubSpec->Get("h0"));
    TH1D *sumBGsubSpecDetALL = static_cast<TH1D*>(sumBGsubSpecALL->Get("h0"));


    // Create the canvas
	TCanvas *c0 = new TCanvas("c0","c0",1920,1080);
	c0->Divide(1,3);
	c0->Update();


	c0->cd(1);
    gPad->SetLogy();
    sumBGsubSpecDet->Draw();
    sumBGsubSpecDet->GetXaxis()->SetRangeUser(1600,3500);

    vector < double > p3PeakNoBG;
    p3PeakNoBG = single_gauss_peak_heaton(2120,2265,sumBGsubSpecDet,false);

    double cent1 = p3PeakNoBG[0];
    double area1 = p3PeakNoBG[4];
	double area1_err = p3PeakNoBG[5];
    double q1 = 8925578.0;  // summed charge for relevant runs
    double yield1 = area1/q1;


    c0->cd(2);
    gPad->SetLogy();
    sumBGsubSpecDetALL->Draw();
    sumBGsubSpecDetALL->GetXaxis()->SetRangeUser(1600,3500);

    vector < double > p3PeakNoBGALL;
    p3PeakNoBGALL = single_gauss_peak_heaton(2130,2284,sumBGsubSpecDetALL);

    double cent2 = p3PeakNoBGALL[0];
    double area2 = p3PeakNoBGALL[4];
	double area2_err = p3PeakNoBGALL[5];
    double q2 = 55601009.0;  // summed charge for all runs
    double yield2 = area2/q2;


    // c0->cd(3);
    // gPad->SetLogy();
    // sumBGsubSpecDetALL->Draw();
    // sumBGsubSpecDetALL->GetXaxis()->SetRangeUser(3000,3500);
    //
    // vector < double > otherPeak;
    // otherPeak = single_gauss_peak_heaton(3250,3400,sumBGsubSpecDetALL);
    //
    // double cent3 = otherPeak[0];
    // double area3 = otherPeak[4];
	// double area3_err = otherPeak[5];


    std::cout << "\nBG SUB:\nCentroid = " << cent1 << "\tArea = " << area1 << " +/- " << area1_err << "\n";
    std::cout << "\nBG SUB ALL:\nCentroid = " << cent2 << "\tArea = " << area2 << " +/- " << area2_err << "\n";
    // std::cout << "\nOther peak:\nCentroid = " << cent3 << "\tArea = " << area3 << " +/- " << area3_err << "\n";
    std::cout << "\np3 count ratio:\nArea = " << area1/area2 <<"\n";
    std::cout << "\np3 yield ratio:\nArea = " << yield1/yield2 <<"\n";

}
