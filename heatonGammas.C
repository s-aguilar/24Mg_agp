// From summed spectras file, calculate yield for Heaton gamma rays

#include "calibration/fitFunctions.h"

void heatonGammasBGsub(){


    TFile *sumSpec = new TFile("E_cal_spectras/summedSpectrasALL.root","READ");
    TFile *sumBGsubSpec = new TFile("E_cal_spectras/summedSpectrasBGsubALL.root","READ");


    TH1D *sumSpecDet = static_cast<TH1D*>(sumSpec->Get("h0"));
    TH1D *sumBGsubSpecDet = static_cast<TH1D*>(sumBGsubSpec->Get("h0"));


    // Create the canvas
	TCanvas *c0 = new TCanvas("c0","c0",1920,1080);
	c0->Divide(1,3);
	c0->Update();


	c0->cd(1);
    gPad->SetLogy();
    sumSpecDet->Draw();
    sumSpecDet->GetXaxis()->SetRangeUser(1600,3500);

    vector < double > p3PeakWithBG;
    p3PeakWithBG = single_gauss_peak_heaton(2120,2265,sumSpecDet,false);

    double cent1 = p3PeakWithBG[0];
    double area1 = p3PeakWithBG[4];
	double area1_err = p3PeakWithBG[5];


    c0->cd(2);
    gPad->SetLogy();
    sumBGsubSpecDet->Draw();
    sumBGsubSpecDet->GetXaxis()->SetRangeUser(1600,3500);

    vector < double > p3PeakNoBG;
    p3PeakNoBG = single_gauss_peak_heaton(2130,2284,sumBGsubSpecDet);

    double cent2 = p3PeakNoBG[0];
    double area2 = p3PeakNoBG[4];
	double area2_err = p3PeakNoBG[5];


    c0->cd(3);
    gPad->SetLogy();
    sumBGsubSpecDet->Draw();
    sumBGsubSpecDet->GetXaxis()->SetRangeUser(3000,3500);

    vector < double > otherPeak;
    otherPeak = single_gauss_peak_heaton(3250,3400,sumSpecDet);

    double cent3 = otherPeak[0];
    double area3 = otherPeak[4];
	double area3_err = otherPeak[5];


    std::cout << "\nNO BG SUB:\nCentroid = " << cent1 << "\tArea = " << area1 << " +/- " << area1_err << "\n";
    std::cout << "\nBG SUB:\nCentroid = " << cent2 << "\tArea = " << area2 << " +/- " << area2_err << "\n";
    std::cout << "\nOther peak:\nCentroid = " << cent3 << "\tArea = " << area3 << " +/- " << area3_err << "\n";

}
