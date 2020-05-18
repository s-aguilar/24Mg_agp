#ifndef FITFUNCTIONS_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define FITFUNCTIONS_H

# include <iostream>

#include <vector>
using std::vector;

#include "TF1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TLine.h"


double background_func(double *x, double *par){
	/*
	COEFFICIENTS:
	===============
	constant = par[0]
	linear = par[1]
	quadratic = par[2]
	*/

	return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
}


double single_gauss_func(double *x, double *par){
	/*
	norm = par[0]
	mean = par[1]
	sigma = par[2]
	*/

	return par[0]*TMath::Gaus(x[0],par[1],par[2],kTRUE);
}


double double_gauss_func(double*x, double *par){
	/*
	norm1 = par[0]
	mean1 = par[1]
	sigma1 = par[2]
	norm2 = par[3]
	mean2 = par[4]
	sigma2 = par[5]
	*/

	double g1 = par[0]*TMath::Gaus(x[0],par[1],par[2],kTRUE);
	double g2 = par[3]*TMath::Gaus(x[0],par[4],par[5],kTRUE);
	return g1+g2;
}


double double_gauss_same_width_func(double*x, double *par){
	/*
	norm1 = par[0]
	mean1 = par[1]
	sigma = par[2]
	norm2 = par[3]
	mean2 = par[4]
	*/

	double g1 = par[0]*TMath::Gaus(x[0],par[1],par[2],kTRUE);
	double g2 = par[3]*TMath::Gaus(x[0],par[4],par[2],kTRUE); // par[5]->par[2]
	return g1+g2;
}


double double_gauss_same_width_func_p1(double*x, double *par){
	/*
	norm1 = par[0]
	mean1 = par[1]
	sigma = par[2]
	norm2 = par[3]
	mean2 = par[4]
	*/

	double g1 = par[0]*TMath::Gaus(x[0],par[1],par[2],kTRUE);
	double g2 = par[3]*TMath::Gaus(x[0],par[1]+26.9,par[2],kTRUE); // par[4]->par[1]+24 , par[5]->par[2]
	return g1+g2;
}


double fit_single_gauss_func(double *x, double *par){
    // Fit a Gaussian with a Quadratic bakground

	return background_func(x,par) + single_gauss_func(x,&par[3]);
}


double fit_double_gauss_func(double *x, double *par) {
    // Fit a double Gaussian with a Quadratic bakground, return array of parameters

   return background_func(x,par) + double_gauss_func(x,&par[3]);
}


double fit_double_gauss_same_width_func(double *x, double *par) {
    // Fit a double Gaussian with a Quadratic bakground, return array of parameters

   return background_func(x,par) + double_gauss_same_width_func(x,&par[3]);
}


double fit_double_gauss_same_width_func_p1(double *x, double *par) {
    // Fit a double Gaussian with a Quadratic bakground, return array of parameters

   return background_func(x,par) + double_gauss_same_width_func_p1(x,&par[3]);
}


/*===========================================================================*/
						///////////////////
						// FINDING PEAKS //
						///////////////////

vector < double > single_gauss_peak_heaton(double low, double high, TH1D *H0, bool isQuadratic = false){

	TF1 *ffit1 = new TF1("ffit1",fit_single_gauss_func,low,high,6);
	ffit1->SetParNames("a0","a1","a2","norm","mean","sigma");
	ffit1->SetNpx(500);
	ffit1->SetParameters(1,1,1,4000,(low+high)/2,12);

	if (isQuadratic) ffit1->SetParLimits(2,-100,100);
	else ffit1->FixParameter(2,0);			// Makes it a linear background

	ffit1->SetParLimits(3,0,1e8);
	ffit1->SetParLimits(4,low,high);		// Peak centroid range
	ffit1->SetParLimits(5,8,18);			// Std dev range

	ffit1->SetLineColor(kBlack);
	H0->Fit("ffit1","SQR");

	double area = ffit1->GetParameter(3);
	double area_err = ffit1->GetParError(3);
	double peakPos = ffit1->GetParameter(4);
	double sigma = ffit1->GetParameter(5);

	// Plot line marking peak position
	TLine *line1 = new TLine(peakPos,H0->GetMinimum(),peakPos,1.05*H0->GetMaximum());
	line1->Draw("SAME");

	TF1 *bgfit = new TF1("bgfit",background_func,low,high,3);
	bgfit->SetNpx(500);
	bgfit->SetParameter(0,ffit1->GetParameter(0));
	bgfit->SetParameter(1,ffit1->GetParameter(1));
	bgfit->SetParameter(2,ffit1->GetParameter(2));
	// bgfit->SetLineColor(3);
	bgfit->Draw("SAME");

	// LOW AND HIGH ARE NOT BEING CALIBRATED RUN BY RUN
	// results format: [centroid,low,high,sigma]
	vector < double > results;
	results.push_back(peakPos);
	results.push_back(low);
	results.push_back(high);
	results.push_back(sigma);
	results.push_back(area);
	results.push_back(area_err);

	return results;
	}

vector < double > single_gauss_peak(double low, double high, TH1D *H0){
	/* FOR FINDING PEAK POSITIONS:
		Fit using fit_single_gauss_func, hone in on parameters and then
		fit it again. Return peak position and its range.
	*/

	TF1 *ffit1 = new TF1("ffit1",fit_single_gauss_func,low,high,6);
	ffit1->SetParNames("a0","a1","a2","norm","mean","sigma");
	ffit1->SetNpx(500);
	ffit1->SetParameters(1,1,1,4000,(low+high)/2,12);

	ffit1->FixParameter(2,0);				// Makes it a linear background

    ffit1->SetParLimits(3,0,1e8);
    ffit1->SetParLimits(4,low,high);		// Peak centroid range
    ffit1->SetParLimits(5,8,15);			// Std dev range

	ffit1->SetLineColor(kBlack);
	H0->Fit("ffit1","SQR");


	double peakPos = ffit1->GetParameter(4);
	double sigma = ffit1->GetParameter(5);

	// Plot line marking peak position
	TLine *line1 = new TLine(peakPos,H0->GetMinimum(),peakPos,1.05*H0->GetMaximum());
	line1->Draw("SAME");

	// LOW AND HIGH ARE NOT BEING CALIBRATED RUN BY RUN
	// results format: [centroid,low,high,sigma]
	vector < double > results;
	results.push_back(peakPos);
	results.push_back(low);
	results.push_back(high);
	results.push_back(sigma);

	return results;
}


vector < double > double_gauss_peak(double low, double high, TH1D *H0){
	/* FOR FINDING PEAK POSITIONS:
		Fit using fit_double_gauss_same_width_func
	*/

	TF1 *ffit1 = new TF1("ffit1",fit_double_gauss_same_width_func,low,high,8);
	ffit1->SetParNames("a0","a1","a2","norm1","mean1","sigma1","norm2","mean2");
	ffit1->SetNpx(500);
	ffit1->SetParameters(1,1,0,4000,(low+high)/2-15,12,4000,(low+high)/2+15);

	// Turn off linear background terms
	ffit1->FixParameter(0,0);
	ffit1->FixParameter(1,0);
	ffit1->FixParameter(2,0);		// Makes it a linear background

	ffit1->SetParLimits(3,0,1e8);	// Normalization
	ffit1->SetParLimits(4,low,high);	// centroid 1
	ffit1->SetParLimits(6,0,1e8);	// Normalization
	ffit1->SetParLimits(7,low,high);	// centroid 2

	ffit1->SetLineColor(kBlack);
	H0->Fit("ffit1","SQR");

	// Store fit parameters and sort
	double par1[8];
	double temp[8];
	ffit1->GetParameters(par1);
	ffit1->GetParameters(temp);

	// I want Peak 2 to be the higher energy peak
	// Case: Peak 1 > Peak2  -> Sort
	if(par1[4]>par1[7]){
		// Lower energy peak is first peak
		par1[3] = temp[6];
		par1[4] = temp[7];

		// Higher energy peak is second peak
		par1[6] = temp[3];
		par1[7] = temp[4];
	}

	double peakPos = par1[7];
	double sigma = par1[5];
	double otherPeak = par1[4];

	// Plot line marking peak position
	TLine *line1 = new TLine(peakPos,H0->GetMinimum(),peakPos,1.05*H0->GetMaximum());
	line1->Draw("SAME");

	//Draw the individual peak
	TF1 *peak1fit = new TF1("peak1fit",fit_single_gauss_func,low,high,6);
	peak1fit->SetNpx(500);
	peak1fit->SetParameter(0,par1[0]);
	peak1fit->SetParameter(1,par1[1]);
	peak1fit->SetParameter(2,par1[2]);
	peak1fit->SetParameter(3,par1[6]);
	peak1fit->SetParameter(4,par1[7]);
	peak1fit->SetParameter(5,par1[5]);
	peak1fit->SetLineColor(3);
	peak1fit->Draw("SAME");

	TF1 *peak2fit = new TF1("peak2fit",fit_single_gauss_func,low,high,6);
	peak2fit->SetNpx(500);
	peak2fit->SetParameter(0,par1[0]);
	peak2fit->SetParameter(1,par1[1]);
	peak2fit->SetParameter(2,par1[2]);
	peak2fit->SetParameter(3,par1[3]);
	peak2fit->SetParameter(4,par1[4]);
	peak2fit->SetParameter(5,par1[5]);
	peak2fit->SetLineColor(3);
	peak2fit->Draw("SAME");

	// Calibrate new ranges
	// low = peakPos-3*sigma;
	// high = peakPos+3*sigma;
	low = peakPos-75;
	high = peakPos+85;

	// LOW AND HIGH ARE NOT BEING CALIBRATED RUN BY RUN
	// results format: [centroid,low,high,sigma,other centroid]
	vector < double > results;
	results.push_back(peakPos);
	results.push_back(low);
	results.push_back(high);
	results.push_back(sigma);
	results.push_back(otherPeak);

	return results;
}


vector < double > BGsingle_gauss_peak(double low, double high, TH1D *H0){
	/* FOR FINDING PEAK POSITIONS:
		Fit using fit_single_gauss_func, hone in on parameters and then
		fit it again. Return peak position and its range.
	*/

	TF1 *ffit1 = new TF1("ffit1",fit_single_gauss_func,low,high,6);
	ffit1->SetParNames("a0","a1","a2","norm","mean","sigma");
	ffit1->SetNpx(500);
	ffit1->SetParameters(1,1,1,4000,(low+high)/2,12);

	ffit1->FixParameter(2,0);				// Makes it a linear background

    ffit1->SetParLimits(3,0,1e8);
    ffit1->SetParLimits(4,low,high);		// Peak centroid range
    ffit1->SetParLimits(5,8,15);			// Std dev range

	ffit1->SetLineColor(kBlack);
	H0->Fit("ffit1","SQR");


	double peakPos = ffit1->GetParameter(4);
	double sigma = ffit1->GetParameter(5);

	// Plot line marking peak position
	TLine *line1 = new TLine(peakPos,H0->GetMinimum(),peakPos,1.05*H0->GetMaximum());
	line1->Draw("SAME");

	// LOW AND HIGH ARE NOT BEING CALIBRATED RUN BY RUN
	// results format: [centroid,low,high,sigma]
	vector < double > results;
	results.push_back(peakPos);
	results.push_back(low);
	results.push_back(high);
	results.push_back(sigma);

	return results;
}


vector < double > BGdouble_gauss_peak(double low, double high, TH1D *H0){
	/* FOR FINDING PEAK POSITIONS:
		Fit using fit_double_gauss_same_width_func
	*/

	TF1 *ffit1 = new TF1("ffit1",fit_double_gauss_same_width_func,low,high,8);
	ffit1->SetParNames("a0","a1","a2","norm1","mean1","sigma1","norm2","mean2");
	ffit1->SetNpx(500);
	ffit1->SetParameters(1,1,0,4000,(low+high)/2-15,12,4000,(low+high)/2+15);

	// Turn off linear background terms
	ffit1->FixParameter(0,0);
	ffit1->FixParameter(1,0);
	ffit1->FixParameter(2,0);		// Makes it a linear background

	ffit1->SetParLimits(3,0,1e8);	// Normalization
	ffit1->SetParLimits(4,low,high);	// centroid 1
	ffit1->SetParLimits(6,0,1e8);	// Normalization
	ffit1->SetParLimits(7,low,high);	// centroid 2

	ffit1->SetLineColor(kBlack);
	H0->Fit("ffit1","SQR");

	// Store fit parameters and sort
	double par1[8];
	double temp[8];
	ffit1->GetParameters(par1);
	ffit1->GetParameters(temp);

	// I want Peak 2 to be the higher energy peak
	// Case: Peak 1 > Peak2  -> Sort
	if(par1[4]>par1[7]){
		// Lower energy peak is first peak
		par1[3] = temp[6];
		par1[4] = temp[7];

		// Higher energy peak is second peak
		par1[6] = temp[3];
		par1[7] = temp[4];
	}

	double peakPos = par1[7];
	double sigma = par1[5];
	double otherPeak = par1[4];

	// Plot line marking peak position
	TLine *line1 = new TLine(peakPos,H0->GetMinimum(),peakPos,1.05*H0->GetMaximum());
	line1->Draw("SAME");

	//Draw the individual peak
	TF1 *peak1fit = new TF1("peak1fit",fit_single_gauss_func,low,high,6);
	peak1fit->SetNpx(500);
	peak1fit->SetParameter(0,par1[0]);
	peak1fit->SetParameter(1,par1[1]);
	peak1fit->SetParameter(2,par1[2]);
	peak1fit->SetParameter(3,par1[6]);
	peak1fit->SetParameter(4,par1[7]);
	peak1fit->SetParameter(5,par1[5]);
	peak1fit->SetLineColor(3);
	peak1fit->Draw("SAME");

	TF1 *peak2fit = new TF1("peak2fit",fit_single_gauss_func,low,high,6);
	peak2fit->SetNpx(500);
	peak2fit->SetParameter(0,par1[0]);
	peak2fit->SetParameter(1,par1[1]);
	peak2fit->SetParameter(2,par1[2]);
	peak2fit->SetParameter(3,par1[3]);
	peak2fit->SetParameter(4,par1[4]);
	peak2fit->SetParameter(5,par1[5]);
	peak2fit->SetLineColor(3);
	peak2fit->Draw("SAME");

	// Calibrate new ranges
	// low = peakPos-3*sigma;
	// high = peakPos+3*sigma;
	low = peakPos-75;
	high = peakPos+85;

	// LOW AND HIGH ARE NOT BEING CALIBRATED RUN BY RUN
	// results format: [centroid,low,high,sigma,other centroid]
	vector < double > results;
	results.push_back(peakPos);
	results.push_back(low);
	results.push_back(high);
	results.push_back(sigma);
	results.push_back(otherPeak);

	return results;
}


vector < double > iterative_single_gauss_peak(double low, double high, TH1D *H0){
	/* FOR FINDING PEAK POSITIONS:
		Fit using fit_single_gauss_func, hone in on parameters and then
		fit it again. Return peak position and its range.
	*/

	// FIRST ITERATION FIT
	TF1 *ffit1 = new TF1("ffit1",fit_single_gauss_func,low,high,6);
	ffit1->SetParNames("a0","a1","a2","norm","mean","sigma");
	ffit1->SetNpx(500);
	ffit1->SetParameters(0,0,0,4000,(low+high)/2,10);

	// Turn off linear background terms
	ffit1->FixParameter(0,0);
	ffit1->FixParameter(1,0);
	ffit1->FixParameter(2,0);		// Makes it a linear background

    ffit1->SetParLimits(3,0,1e6);
    ffit1->SetParLimits(4,low,high);	// Peak centroid range
    ffit1->SetParLimits(5,5,10);	// Std dev range

	H0->Fit("ffit1","SQRN0");

	// Store fit parameters
	double par1[6];
	ffit1->GetParameters(par1);

	// Recalibrate range of integration
	low = par1[4]-2*par1[5]; // 3
	high = par1[4]+2.5*par1[5]; // 3


	// SECOND ITERATION FIT
	TF1 *ffit2 = new TF1("ffit2",fit_single_gauss_func,low,high,6);
	ffit2->SetParNames("a0","a1","a2","norm","mean","sigma");
	ffit2->SetNpx(500);
	ffit2->SetParameters(par1);

	ffit2->FixParameter(0,0);
	ffit2->FixParameter(1,0);
	ffit2->FixParameter(2,0);		// Makes it a linear background

	ffit2->SetParLimits(3,0,1e6);
	ffit2->SetParLimits(4,low,high);	// Peak centroid range
	ffit2->SetParLimits(5,5,15);	// Std dev range

	ffit2->SetLineColor(kBlack);
	H0->Fit("ffit2","SQRN0");

	double peakPos = ffit2->GetParameter(4);
	double sigma = ffit2->GetParameter(5);

	// Plot line marking peak position
	TLine *line1 = new TLine(peakPos,H0->GetMinimum(),peakPos,1.0*H0->GetMaximum());
	line1->Draw("SAME");

	// results format: [centroid,low,high,sigma]
	vector < double > results;
	results.push_back(peakPos);
	results.push_back(low);
	results.push_back(high);
	results.push_back(sigma);

	return results;
}


vector < double > iterative_double_gauss_peak(double low, double high, TH1D *H0){
	/* FOR FINDING PEAK POSITIONS:
		Fit using fit_double_gauss_same_width_func, hone in on parameters and then
		fit it again.
	*/

	// FIRST ITERATION FIT
	TF1 *ffit1 = new TF1("ffit1",fit_double_gauss_func,low,high,9);
	ffit1->SetParNames("a0","a1","a2","norm1","mean1","sigma1","norm2","mean2","sigma2");
	ffit1->SetNpx(500);
	ffit1->SetParameters(0,0,0,4000,(low+high)/2-15,12,4000,(low+high)/2+15,12);

	// Turn off linear background terms
	// ffit1->FixParameter(0,0);
	// ffit1->FixParameter(1,0);
	ffit1->FixParameter(2,0);		// Makes it a linear background

	ffit1->SetParLimits(3,0,1e8);	// Normalization
	ffit1->SetParLimits(4,low,high);	// centroid 1
	ffit1->SetParLimits(5,8,13);
	ffit1->SetParLimits(6,0,1e8);	// Normalization
	ffit1->SetParLimits(7,low,high);	// centroid 2
	ffit1->SetParLimits(8,8,18);

	H0->Fit("ffit1","SQRN0");

	// Store fit parameters
	double par1[9];
	ffit1->GetParameters(par1);

	// Sort peaks, I want the higher energy peak
	double temp[9];
	for(int ii = 0;ii<9;ii++){
		temp[ii] = par1[ii];
	}

	// Case: Peak 1 > Peak2  -> Sort
	if(par1[4]>par1[7]){
		// Lower energy peak is first peak
		par1[3] = temp[6];
		par1[4] = temp[7];
		par1[5] = temp[8];

		// Higher energy peak is second peak
		par1[6] = temp[3];
		par1[7] = temp[4];
		par1[8] = temp[5];
	}


	// Recalibrate range of integration
	low = par1[4]-35;
	high = par1[7]+4*par1[5];


	// SECOND ITERATION FIT
	TF1 *ffit2 = new TF1("ffit2",fit_double_gauss_func,low,high,9);
	ffit2->SetParNames("a0","a1","a2","norm1","mean1","sigma1","norm2","mean2","sigma2");
	ffit2->SetNpx(500);
	ffit2->SetParameters(par1);

	ffit2->FixParameter(2,0);		// Makes it a linear background
	ffit2->SetParLimits(3,0,1e8);
	ffit2->SetParLimits(4,low,high);	// Centroid 1
	ffit2->SetParLimits(5,8,13);	// Std dev range
	ffit2->SetParLimits(6,0,1e8);
	ffit2->SetParLimits(7,low,high);	// centroid 2
	ffit2->SetParLimits(8,8,18);	// Std dev range

	ffit2->SetLineColor(kBlack);
	H0->Fit("ffit2","SQRN0");

	// Store fit parameters
	double par2[9];
	ffit1->GetParameters(par2);

	double peakPos = ffit2->GetParameter(4);
	double sigma = ffit2->GetParameter(5);
	double otherPeak = ffit2->GetParameter(7);

	// Sort peak
	if (ffit2->GetParameter(7) > peakPos){
		otherPeak = peakPos;
		peakPos = ffit2->GetParameter(7);
		sigma = ffit2->GetParameter(8);
	}

	// Plot line marking peak position
	TLine *line1 = new TLine(peakPos,H0->GetMinimum(),peakPos,1.05*H0->GetMaximum());
	line1->Draw("SAME");

	//Draw the individual peak
	TF1 *peak1fit = new TF1("peak1fit",fit_single_gauss_func,low,high,6);
	peak1fit->SetNpx(500);
	peak1fit->SetParameter(0,par2[0]);
	peak1fit->SetParameter(1,par2[1]);
	peak1fit->SetParameter(2,par2[2]);
	peak1fit->SetParameter(3,par2[6]);
	peak1fit->SetParameter(4,par2[7]);
	peak1fit->SetParameter(5,par2[8]);
	peak1fit->SetLineColor(3);
	peak1fit->Draw("SAME");

	TF1 *peak2fit = new TF1("peak2fit",fit_single_gauss_func,low,high,6);
	peak2fit->SetNpx(500);
	peak2fit->SetParameter(0,par2[0]);
	peak2fit->SetParameter(1,par2[1]);
	peak2fit->SetParameter(2,par2[2]);
	peak2fit->SetParameter(3,par2[3]);
	peak2fit->SetParameter(4,par2[4]);
	peak2fit->SetParameter(5,par2[5]);
	peak2fit->SetLineColor(3);
	peak2fit->Draw("SAME");

	// results format: [centroid,low,high,sigma,other centroid]
	vector < double > results;
	results.push_back(peakPos);
	results.push_back(low);
	results.push_back(high);
	results.push_back(sigma);
	results.push_back(otherPeak);

	return results;
}


vector < double > iterative_double_gauss_peak_same_width(double low, double high, TH1D *H0){
	/* FOR FINDING PEAK POSITIONS:
		Fit using fit_double_gauss_same_width_func, hone in on parameters and then
		fit it again.
	*/

	// FIRST ITERATION FIT
	TF1 *ffit1 = new TF1("ffit1",fit_double_gauss_same_width_func,low,high,8); //fit_double_gauss_func
	ffit1->SetParNames("a0","a1","a2","norm1","mean1","sigma1","norm2","mean2");
	ffit1->SetNpx(500);
	ffit1->SetParameters(0,0,0,4000,(low+high)/2-15,12,4000,(low+high)/2+15);

	// Turn off linear background terms
	ffit1->FixParameter(0,0);
	ffit1->FixParameter(1,0);
	ffit1->FixParameter(2,0);		// Makes it a linear background

	ffit1->SetParLimits(3,0,1e8);	// Normalization
	ffit1->SetParLimits(4,low,high);	// centroid 1
	ffit1->SetParLimits(5,10,20);
	ffit1->SetParLimits(6,0,1e8);	// Normalization
	ffit1->SetParLimits(7,low,high);	// centroid 2

	H0->Fit("ffit1","SQRN0");

	// Store fit parameters
	double par1[8];
	ffit1->GetParameters(par1);

	// Sort peaks, I want the higher energy peak
	double temp[8];
	for(int ii = 0;ii<8;ii++){
		temp[ii] = par1[ii];
	}

	// Case: Peak 1 > Peak2  -> Sort
	if(par1[4]>par1[7]){
		// Lower energy peak is first peak
		par1[3] = temp[6];
		par1[4] = temp[7];

		// Higher energy peak is second peak
		par1[6] = temp[3];
		par1[7] = temp[4];
	}


	// Recalibrate range of integration
	low = par1[4]-35;
	high = par1[7]+3*par1[5];


	// SECOND ITERATION FIT
	TF1 *ffit2 = new TF1("ffit2",fit_double_gauss_func,low,high,8);
	ffit2->SetParNames("a0","a1","a2","norm1","mean1","sigma1","norm2","mean2");
	ffit2->SetNpx(500);
	ffit2->SetParameters(par1);

	ffit2->FixParameter(2,0);		// Makes it a linear background
	ffit2->SetParLimits(3,0,1e8);
	ffit2->SetParLimits(4,low,high);	// Centroid 1
	ffit2->SetParLimits(5,10,20);	// Std dev range
	ffit2->SetParLimits(6,0,1e8);
	ffit2->SetParLimits(7,low,high);	// centroid 2

	ffit2->SetLineColor(kBlack);
	H0->Fit("ffit2","SQR");

	// Store fit parameters
	double par2[8];
	ffit1->GetParameters(par2);

	double peakPos = ffit2->GetParameter(4);
	double sigma = ffit2->GetParameter(5);
	double otherPeak = ffit2->GetParameter(7);

	// Sort peak
	if (ffit2->GetParameter(7) > peakPos){
		otherPeak = peakPos;
		peakPos = ffit2->GetParameter(7);
		sigma = ffit2->GetParameter(5);
	}

	// Plot line marking peak position
	TLine *line1 = new TLine(peakPos,H0->GetMinimum(),peakPos,1.05*H0->GetMaximum());
	line1->Draw("SAME");

	//Draw the individual peak
	TF1 *peak1fit = new TF1("peak1fit",fit_single_gauss_func,low,high,6);
	peak1fit->SetNpx(500);
	peak1fit->SetParameter(0,par2[0]);
	peak1fit->SetParameter(1,par2[1]);
	peak1fit->SetParameter(2,par2[2]);
	peak1fit->SetParameter(3,par2[6]);
	peak1fit->SetParameter(4,par2[7]);
	peak1fit->SetParameter(5,par2[5]);
	peak1fit->SetLineColor(3);
	peak1fit->Draw("SAME");

	TF1 *peak2fit = new TF1("peak2fit",fit_single_gauss_func,low,high,6);
	peak2fit->SetNpx(500);
	peak2fit->SetParameter(0,par2[0]);
	peak2fit->SetParameter(1,par2[1]);
	peak2fit->SetParameter(2,par2[2]);
	peak2fit->SetParameter(3,par2[3]);
	peak2fit->SetParameter(4,par2[4]);
	peak2fit->SetParameter(5,par2[5]);
	peak2fit->SetLineColor(3);
	peak2fit->Draw("SAME");

	// results format: [centroid,low,high,sigma,other centroid]
	vector < double > results;
	results.push_back(peakPos);
	results.push_back(low);
	results.push_back(high);
	results.push_back(sigma);
	results.push_back(otherPeak);

	return results;
}


/*===========================================================================*/
					//////////////////////////////////////
					// FINDING PEAKS HEALTHY BACKGROUND //
					//////////////////////////////////////

vector < double > iterative_single_gauss_peakWide(double low, double high, TH1D *H0){
	/* FOR FINDING PEAK POSITIONS:
		Fit using fit_single_gauss_func, hone in on parameters and then
		fit it again. Return peak position and its range.
	*/

	// FIRST ITERATION FIT
	TF1 *ffit1 = new TF1("ffit1",fit_single_gauss_func,low,high,6);
	ffit1->SetParNames("a0","a1","a2","norm","mean","sigma");
	ffit1->SetNpx(500);
	ffit1->SetParameters(1,1,1,4000,(low+high)/2,12);

	// Turn off linear background terms
	ffit1->FixParameter(0,0);
	ffit1->FixParameter(1,0);
	ffit1->FixParameter(2,0);		// Makes it a linear background

    ffit1->SetParLimits(3,0,1e6);
    ffit1->SetParLimits(4,low,high);	// Peak centroid range
    ffit1->SetParLimits(5,2.5,55);	// Std dev range

	H0->Fit("ffit1","SQRN0");

	// Store fit parameters
	double par1[6];
	ffit1->GetParameters(par1);

	// Recalibrate range of integration
	low = par1[4]-3*par1[5];
	high = par1[4]+3*par1[5];


	// SECOND ITERATION FIT
	TF1 *ffit2 = new TF1("ffit2",fit_single_gauss_func,low,high,6);
	ffit2->SetParNames("a0","a1","a2","norm","mean","sigma");
	ffit2->SetNpx(500);
	ffit2->SetParameters(par1);
	ffit2->FixParameter(2,0);		// Makes it a linear background
	ffit2->SetParLimits(3,0,1e6);
	ffit2->SetParLimits(4,low,high);	// Peak centroid range
	ffit2->SetParLimits(5,2.5,55);	// Std dev range
	H0->Fit("ffit2","SQR");


	double peakPos = ffit2->GetParameter(4);
	double sigma = ffit2->GetParameter(5);

	// results format: [centroid,low,high,sigma]
	vector < double > results;
	results.push_back(peakPos);
	results.push_back(low);
	results.push_back(high);
	results.push_back(sigma);

	return results;
}


vector < double > iterative_double_gauss_peakWide(double low, double high, TH1D *H0){
	/* FOR FINDING PEAK POSITIONS:
		Fit using fit_double_gauss_func, hone in on parameters and then
		fit it again.
	*/

	// FIRST ITERATION FIT
	TF1 *ffit1 = new TF1("ffit1",fit_double_gauss_func,low,high,9);
	ffit1->SetParNames("a0","a1","a2","norm1","mean1","sigma1","norm2","mean2","sigma2");
	ffit1->SetNpx(500);
	ffit1->SetParameters(1,1,1,4000,(low+high)/2-20,12,4000,(low+high)/2+20,12);

	// Turn off linear background terms
	ffit1->FixParameter(0,0);
	ffit1->FixParameter(1,0);
	ffit1->FixParameter(2,0);		// Makes it a linear background

	ffit1->SetParLimits(3,0,1e8);	// Normalization
	ffit1->SetParLimits(5,2.5,15);	// Std dev range
	ffit1->SetParLimits(6,0,1e8);
	ffit1->SetParLimits(8,2.5,15);

	H0->Fit("ffit1","SQRN0");

	// Store fit parameters
	double par1[9];
	ffit1->GetParameters(par1);

	// Sort peaks, I want the lower energy peak
	double temp[9];
	for(int ii = 0;ii<6;ii++){
		temp[ii] = par1[ii];
	}

	// Case: Peak 1 > Peak2  -> Sort
	if(par1[4]>par1[7]){
		// Lower energy peak is first peak
		par1[3] = temp[6];
		par1[4] = temp[7];
		par1[5] = temp[8];

		// Higher energy peak is second peak
		par1[6] = temp[3];
		par1[7] = temp[4];
		par1[8] = temp[5];
	}


	// Recalibrate range of integration
	if(par1[4]<par1[7]){
		// low = par1[4]-2.2*par1[5]; // Constrain left bound
		low = par1[4]-22; // 30
		high = par1[7]+3*par1[8]+30;
	}
	else{
		low = par1[7]-22; // 30;
		high = par1[4]+3*par1[4]+30;
	}


	// SECOND ITERATION FIT
	TF1 *ffit2 = new TF1("ffit2",fit_double_gauss_func,low,high,9);
	ffit2->SetParNames("a0","a1","a2","norm1","mean1","sigma1","norm2","mean2","sigma2");
	ffit2->SetNpx(500);
	ffit2->SetParameters(par1);

	ffit2->FixParameter(2,0);		// Makes it a linear background
	ffit2->SetParLimits(3,0,1e8);
	ffit2->SetParLimits(5,2.5,15);	// Std dev range
	ffit2->SetParLimits(6,0,1e8);
	ffit2->SetParLimits(8,2.5,15);

	H0->Fit("ffit2","SQR");

	double peakPos = ffit2->GetParameter(4);
	double sigma = ffit2->GetParameter(5);
	double otherPeak = ffit2->GetParameter(7);

	// Sort peak
	if (ffit2->GetParameter(7) > peakPos){
		otherPeak = peakPos;
		peakPos = ffit2->GetParameter(7);
		sigma = ffit2->GetParameter(8);
	}

	// results format: [centroid,low,high,sigma,other centroid]
	vector < double > results;
	results.push_back(peakPos);
	results.push_back(low);
	results.push_back(high);
	results.push_back(sigma);
	results.push_back(otherPeak);

	return results;
}


/*===========================================================================*/
						///////////////////
						// FITTING PEAKS //
						///////////////////

vector < double > iterative_single_gauss_area(double low, double high, TH1D *H0){

	/* FOR FINDING PEAK POSITIONS:
		Fit using fit_single_gauss_func, hone in on parameters and then
		fit it again. Return peak area.
		*/

	// FIRST ITERATION FIT
	TF1 *ffit1 = new TF1("ffit1",fit_single_gauss_func,low,high,6);
	ffit1->SetParNames("a0","a1","a2","norm","mean","sigma");
	ffit1->SetNpx(500);
	ffit1->SetParameters(0,0,0,4000,(low+high)/2,12);
	ffit1->FixParameter(2,0);		// Makes it a linear background
	ffit1->SetParLimits(3,0,1e8);
	ffit1->SetParLimits(4,low,high);	// Peak centroid range
	ffit1->SetParLimits(5,2.5,18);		// Std dev range

	H0->Fit("ffit1","SQRN0");

	// Store fit parameters
	double par1[6];
	ffit1->GetParameters(par1);

	// Recalibrate range of integration
	low = par1[4]-4*par1[5];
	high = par1[4]+4*par1[5];


	// SECOND ITERATION FIT
	TF1 *ffit2 = new TF1("ffit2",fit_single_gauss_func,low,high,6);
	ffit2->SetParNames("a0","a1","a2","norm","mean","sigma");
	ffit2->SetNpx(500);
	ffit2->SetParameters(par1);
	ffit2->FixParameter(2,0);		// Makes it a linear background
	ffit2->SetParLimits(3,0,1e6);
	ffit2->SetParLimits(4,low,high);	// Peak centroid range
	ffit2->SetParLimits(5,2.5,18);	// Std dev range
	H0->Fit("ffit2","SQRN0");


	double peakPos = ffit2->GetParameter(4);
	double sigma = ffit2->GetParameter(5);

	// Store fit parameters
	double par2[6];
	ffit2->GetParameters(par2);

	TF1 *fback = new TF1("fback",background_func,low,high,3);
	fback->SetParNames("a0","a1","a2");
	fback->SetLineColor(kCyan);
	fback->SetNpx(1e5);

	fback->SetParameters(par2[0],par2[1],par2[2]);
	fback->SetParErrors(ffit2->GetParErrors());	// Set parameter errors on BG terms from previous fit
	fback->Draw("SAME");


	// Calculate peak area
	// Histogram integration
	double area1_err;
	double area1 = H0->IntegralAndError(par2[4]-3*par2[5],par2[4]+3*par2[5],area1_err,"width");

	// BG fit function integration
	double area2 = fback->Integral(par2[4]-3*par2[5],par2[4]+3*par2[5]);
	double area2_err = TMath::Sqrt(abs(area2));	// THIS IS WRONG, TEMPORARY

	// Peak area and its error (NEW)
	double area = area1 - area2;
	double area_err = TMath::Sqrt(area1_err * area1_err + area2_err * area2_err);

	double chi2NDF = ffit2->GetChisquare()/ffit2->GetNDF();


	// results format: [area, area err, chi2NDF]
	vector < double > results;
	results.push_back(area);
	results.push_back(area_err);
	results.push_back(chi2NDF);

	return results;
}


vector < double > iterative_double_gauss_area(double low, double high, TH1D *H0){
	/* FOR FINDING PEAK POSITIONS:
		Fit using fit_double_gauss_func, hone in on parameters and then
		fit it again. Return peak 1 area
	*/

	// FIRST ITERATION FIT
	TF1 *ffit1 = new TF1("ffit1",fit_double_gauss_func,low,high,9);
	ffit1->SetParNames("a0","a1","a2","norm1","mean1","sigma1","norm2","mean2","sigma2");
	ffit1->SetNpx(500);
	ffit1->SetParameters(0,0,0,4000,low+40,12,4000,high-50,12);

	// Turn off linear background terms
	// ffit1->FixParameter(0,0);
	// ffit1->FixParameter(1,0);
	ffit1->FixParameter(2,0);		// Makes it a linear background

	ffit1->SetParLimits(3,0,1e6);	// Normalization
	ffit1->SetParLimits(5,2.5,15);	// Std dev range
	ffit1->SetParLimits(6,0,1e6);
	ffit1->SetParLimits(8,2.5,15);

	H0->Fit("ffit1","SQRN0");

	// Store fit parameters
	double par1[9];
	ffit1->GetParameters(par1);

	// Sort peaks, I want the lower energy peak
	double temp[9];
	for(int ii = 0;ii<9;ii++){
		temp[ii] = par1[ii];
	}

	// Case: Peak 1 > Peak2  -> Sort
	if(par1[4]>par1[7]){
		// Lower energy peak is first peak
		par1[3] = temp[6];
		par1[4] = temp[7];
		par1[5] = temp[8];

		// Higher energy peak is second peak
		par1[6] = temp[3];
		par1[7] = temp[4];
		par1[8] = temp[5];
	}


	// Recalibrate range of integration
	if(par1[4]<par1[7]){
		// low = par1[4]-2.2*par1[5]; // Constrain left bound
		low = par1[4]-35;
		// high = par1[7]+3*par1[8];
		high = low +100;
	}
	else{
		low = par1[7]-35;
		// high = par1[4]+3*par1[4];
		high = low +100;
	}


	// SECOND ITERATION FIT
	TF1 *ffit2 = new TF1("ffit2",fit_double_gauss_func,low,high,9);
	ffit2->SetParNames("a0","a1","a2","norm1","mean1","sigma1","norm2","mean2","sigma2");
	ffit2->SetNpx(500);
	ffit2->SetParameters(par1);

	ffit2->FixParameter(2,0);		// Makes it a linear background
	ffit2->SetParLimits(3,0,1e6);
	ffit2->SetParLimits(4,low,high);
	ffit2->SetParLimits(5,2,15);	// Std dev range
	ffit2->SetParLimits(6,0,1e6);
	ffit2->SetParLimits(7,low,high);
	ffit2->SetParLimits(8,2,15);

	H0->Fit("ffit2","SQR");
	ffit2->Draw("SAME");

	// Store fit parameters
	double par2[9];
	ffit2->GetParameters(par2);

	double area = par2[3];
	double area_err = ffit2->GetParError(3);

	// Sort peak
	if (par2[7] < par2[4]){
		area = par2[6];
		area_err = ffit2->GetParError(6);
	}

	// Plot line marking peak position
	TLine *line1 = new TLine(par2[4],0,par2[4],600);
	TLine *line2 = new TLine(par2[7],0,par2[7],600);
	line1->Draw("SAME");
	line2->Draw("SAME");

	//Draw the individual peaks
	TF1 *peak1fit = new TF1("peak1fit",fit_single_gauss_func,low,high,6);
	peak1fit->SetNpx(500);
	// peak1fit->SetParameters(par2);
	peak1fit->SetParameter(0,par2[0]);
	peak1fit->SetParameter(1,par2[1]);
	peak1fit->SetParameter(2,par2[2]);
	peak1fit->SetParameter(3,par2[3]);
	peak1fit->SetParameter(4,par2[4]);
	peak1fit->SetParameter(5,par2[5]);
	peak1fit->SetLineColor(3);
	peak1fit->Draw("SAME");
	TF1 *peak2fit = new TF1("peak2fit",fit_single_gauss_func,low,high,6);
	peak2fit->SetNpx(500);
	peak2fit->SetParameter(0,par2[0]);
	peak2fit->SetParameter(1,par2[1]);
	peak2fit->SetParameter(2,par2[2]);
	peak2fit->SetParameter(3,par2[6]);
	peak2fit->SetParameter(4,par2[7]);
	peak2fit->SetParameter(5,par2[8]);
	peak2fit->SetLineColor(3);
	peak2fit->Draw("SAME");


	double chi2NDF = ffit2->GetChisquare()/ffit2->GetNDF();

	// results format: [area, area err, chi2NDF]
	vector < double > results;
	results.push_back(area);
	results.push_back(area_err);
	results.push_back(chi2NDF);

	return results;
}


vector < double > single_gauss_area_noBG(double peak1, TH1D *H0){
	/* FOR FINDING PEAK POSITIONS:
		Fit using fit_double_gauss_func, hone in on parameters and then
		fit it again. Return peak 1 area
	*/

	double low = peak1-70;
	double high = peak1+60;

	TF1 *ffit1 = new TF1("ffit1",fit_single_gauss_func,low,high,6);
	ffit1->SetParNames("a0","a1","a2","norm","mean","sigma");
	ffit1->SetNpx(500);
	ffit1->SetParameters(0,0,0,4000,peak1,10);

	// Turn off linear background terms, fix peak positions
	ffit1->FixParameter(0,0);
	ffit1->FixParameter(1,0);
	ffit1->FixParameter(2,0);


	// ffit1->FixParameter(4,peak1);
	ffit1->SetParLimits(4,peak1-10,peak1+10);

	ffit1->SetParLimits(3,0,1e6);	// Normalization

	ffit1->SetParLimits(5,5,14);	// Std dev range

	// H0->Fit("ffit1","SQRE");
	TFitResultPtr r = H0->Fit("ffit1","SQRE");

	// Store fit parameters
	double par1[6];
	ffit1->GetParameters(par1);

	// Plot line marking peak position
	TLine *line1 = new TLine(peak1,H0->GetMinimum(),peak1,1.05*H0->GetMaximum());
	line1->Draw("SAME");

	TF1 *fback = new TF1("fback",background_func,low,high,3);
	fback->SetParNames("a0","a1","a2");
	fback->SetNpx(500);
	fback->SetParameters(par1[0],par1[1],par1[2]);
	// fback->SetParErrors(ffit1->GetParErrors());	// Set parameter errors on BG terms from previous fit
	fback->SetLineColor(kCyan);
	fback->Draw("SAME");

	double chi2NDF = ffit1->GetChisquare()/ffit1->GetNDF();

	// results format: [area, area err, chi2NDF, sigma1, linear, offset]
	vector < double > results;
	results.push_back(ffit1->GetParameter(3));
	if (r->IsValid()) results.push_back(ffit1->GetParError(3));
	else results.push_back(TMath::Sqrt(ffit1->GetParameter(3)));
	// results.push_back(ffit1->GetParError(3));
	results.push_back(chi2NDF);
	results.push_back(ffit1->GetParameter(5));
	results.push_back(ffit1->GetParameter(1));
	results.push_back(ffit1->GetParameter(0));

	return results;

}
vector < double > single_gauss_area(double peak1, TH1D *H0){
	/* FOR FINDING PEAK POSITIONS:
		Fit using fit_double_gauss_func, hone in on parameters and then
		fit it again. Return peak 1 area
	*/

	double low = peak1-70;
	double high = peak1+60;

	TF1 *ffit1 = new TF1("ffit1",fit_single_gauss_func,low,high,6);
	ffit1->SetParNames("a0","a1","a2","norm","mean","sigma");
	ffit1->SetNpx(500);
	ffit1->SetParameters(0,0,0,4000,peak1,10);

	// Turn off linear background terms, fix peak positions
	ffit1->FixParameter(2,0);		// Makes it a linear background


	// ffit1->FixParameter(4,peak1);
	ffit1->SetParLimits(4,peak1-10,peak1+10);

	ffit1->SetParLimits(3,0,1e6);	// Normalization

	ffit1->SetParLimits(5,7,14);	// Std dev range

	// H0->Fit("ffit1","SQRE");
	TFitResultPtr r = H0->Fit("ffit1","SQRE");

	// Store fit parameters
	double par1[6];
	ffit1->GetParameters(par1);

	// Plot line marking peak position
	TLine *line1 = new TLine(peak1,H0->GetMinimum(),peak1,1.05*H0->GetMaximum());
	line1->Draw("SAME");

	TF1 *fback = new TF1("fback",background_func,low,high,3);
	fback->SetParNames("a0","a1","a2");
	fback->SetNpx(500);
	fback->SetParameters(par1[0],par1[1],par1[2]);
	// fback->SetParErrors(ffit1->GetParErrors());	// Set parameter errors on BG terms from previous fit
	fback->SetLineColor(kCyan);
	fback->Draw("SAME");

	double chi2NDF = ffit1->GetChisquare()/ffit1->GetNDF();

	// results format: [area, area err, chi2NDF, sigma1, linear, offset]
	vector < double > results;
	results.push_back(ffit1->GetParameter(3));
	if (r->IsValid()) results.push_back(ffit1->GetParError(3));
	else results.push_back(TMath::Sqrt(ffit1->GetParameter(3)));
	// results.push_back(ffit1->GetParError(3));
	results.push_back(chi2NDF);
	results.push_back(ffit1->GetParameter(5));
	results.push_back(ffit1->GetParameter(1));
	results.push_back(ffit1->GetParameter(0));

	return results;

}


vector < double > double_gauss_area(double peak1, double peak2, TH1D *H0){
	/* FOR FINDING PEAK POSITIONS:
		Fit using fit_double_gauss_func, hone in on parameters and then
		fit it again. Return peak 1 area
	*/

	double low = peak1-60;
	double high = peak2+60;

	TF1 *ffit1 = new TF1("ffit1",fit_double_gauss_same_width_func,low,high,8); //9
	ffit1->SetParNames("a0","a1","a2","norm1","mean1","sigma1","norm2","mean2","sigma2");
	ffit1->SetNpx(500);
	ffit1->SetParameters(0,0,0,4000,peak1,11,4000,peak2,11);

	// Turn off linear background terms, fix peak positions
	ffit1->FixParameter(2,0);		// Makes it a linear background


	ffit1->FixParameter(4,peak1);  //////// WORK ON THIS
	ffit1->FixParameter(7,peak2); /////////

	ffit1->SetParLimits(3,0,1e8);	// Normalization
	ffit1->SetParLimits(5,9,13);	// width range
	ffit1->SetParLimits(6,0,1e8);
	// ffit1->SetParLimits(8,2.5,20);

	H0->Fit("ffit1","SQRE");

	// Store fit parameters
	double par1[8]; //9
	ffit1->GetParameters(par1);

	// Plot line marking peak position
	TLine *line1 = new TLine(peak1,H0->GetMinimum(),peak1,1.05*H0->GetMaximum());
	TLine *line2 = new TLine(peak2,H0->GetMinimum(),peak2,1.05*H0->GetMaximum());
	line1->Draw("SAME");
	line2->Draw("SAME");

	//Draw the individual peaks
	TF1 *peak1fit = new TF1("peak1fit",fit_single_gauss_func,low,high,6);
	peak1fit->SetNpx(500);
	peak1fit->SetParameter(0,par1[0]);
	peak1fit->SetParameter(1,par1[1]);
	peak1fit->SetParameter(2,par1[2]);
	peak1fit->SetParameter(3,par1[3]);
	peak1fit->SetParameter(4,par1[4]);
	peak1fit->SetParameter(5,par1[5]);
	peak1fit->SetLineColor(3);
	peak1fit->Draw("SAME");
	TF1 *peak2fit = new TF1("peak2fit",fit_single_gauss_func,low,high,6);
	peak2fit->SetNpx(500);
	peak2fit->SetParameter(0,par1[0]);
	peak2fit->SetParameter(1,par1[1]);
	peak2fit->SetParameter(2,par1[2]);
	peak2fit->SetParameter(3,par1[6]);
	peak2fit->SetParameter(4,par1[7]);
	peak2fit->SetParameter(5,par1[5]);//8
	peak2fit->SetLineColor(3);
	peak2fit->Draw("SAME");

	TF1 *fback = new TF1("fback",background_func,low,high,3);
	fback->SetParNames("a0","a1","a2");
	fback->SetNpx(500);
	fback->SetParameters(par1[0],par1[1],par1[2]);
	// fback->SetParErrors(ffit1->GetParErrors());	// Set parameter errors on BG terms from previous fit
	fback->SetLineColor(kCyan);
	fback->Draw("SAME");

	double chi2NDF = ffit1->GetChisquare()/ffit1->GetNDF();

	// results format: [area, area err, chi2NDF, sigma1, sigma2, linear, offset]
	vector < double > results;
	results.push_back(ffit1->GetParameter(3));
	if (ffit1->GetParameter(3) < ffit1->GetParError(3)) results.push_back(TMath::Sqrt(ffit1->GetParameter(3)));
	else results.push_back(ffit1->GetParError(3));
	results.push_back(chi2NDF);
	results.push_back(par1[5]);
	results.push_back(par1[5]); //8
	results.push_back(ffit1->GetParameter(1));
	results.push_back(ffit1->GetParameter(0));

	return results;

}


vector < double > double_gauss_area_p1(double peak1, double peak2, TH1D *H0, int detLoop){
	/* FOR FINDING PEAK POSITIONS:
		Fit using fit_double_gauss_same_width_func_p1.
		Return fit parameters including peak 1 and peak 2 area
	*/

	double low = peak1-58;
	double high = peak2+60;

	double E_delta_peaks = 26.9;	// The energy difference between the two peaks

	TF1 *ffit1 = new TF1("ffit1",fit_double_gauss_same_width_func_p1,low,high,7); //9
	ffit1->SetParNames("a0","a1","a2","norm1","mean1","sigma1","norm2");
	ffit1->SetNpx(500);
	ffit1->SetParameters(0,0,0,4000,peak1,11,200);

	// Turn off linear background terms, fix peak positions
	ffit1->FixParameter(2,0);		// Makes it a linear background


	// ffit1->FixParameter(4,peak1);  //////// WORK ON THIS
	// ffit1->FixParameter(7,peak2); /////////

	ffit1->SetParLimits(3,0,1e8);	// peak 1 area
	ffit1->SetParLimits(4,peak1-14,peak1+5); // peak 1 centroid
	if(detLoop>6) ffit1->SetParLimits(4,peak1-17,peak1+5);

	ffit1->SetParLimits(5,9,13);	// peak 1 width
	ffit1->SetParLimits(6,0,1e8);	// peak 2 area

	H0->Fit("ffit1","SQR"); // E

	// Store fit parameters
	double par1[7];
	ffit1->GetParameters(par1);

	// Plot line marking peak position
	TLine *line1 = new TLine(peak1,H0->GetMinimum(),peak1,1.05*H0->GetMaximum());
	TLine *line2 = new TLine(peak2,H0->GetMinimum(),peak2,1.05*H0->GetMaximum());
	line1->Draw("SAME");
	line2->Draw("SAME");

	// Plot line marking expected peak position
	TLine *line3 = new TLine(par1[4],H0->GetMinimum(),par1[4],1.05*H0->GetMaximum());
	TLine *line4 = new TLine(par1[4]+E_delta_peaks,H0->GetMinimum(),par1[4]+E_delta_peaks,1.05*H0->GetMaximum());
	line3->SetLineColor(kViolet);
	line4->SetLineColor(kViolet);
	line3->Draw("SAME");
	line4->Draw("SAME");

	//Draw the individual peaks
	TF1 *peak1fit = new TF1("peak1fit",fit_single_gauss_func,low,high,6);
	peak1fit->SetNpx(500);
	peak1fit->SetParameter(0,par1[0]);
	peak1fit->SetParameter(1,par1[1]);
	peak1fit->SetParameter(2,par1[2]);
	peak1fit->SetParameter(3,par1[3]);
	peak1fit->SetParameter(4,par1[4]);
	peak1fit->SetParameter(5,par1[5]);
	peak1fit->SetLineColor(3);
	peak1fit->Draw("SAME");
	TF1 *peak2fit = new TF1("peak2fit",fit_single_gauss_func,low,high,6);
	peak2fit->SetNpx(500);
	peak2fit->SetParameter(0,par1[0]);
	peak2fit->SetParameter(1,par1[1]);
	peak2fit->SetParameter(2,par1[2]);
	peak2fit->SetParameter(3,par1[6]);
	peak2fit->SetParameter(4,par1[4]+E_delta_peaks); // par1[7]->par1[4]+E_delta_peaks
	peak2fit->SetParameter(5,par1[5]);    // par1[8]
	peak2fit->SetLineColor(3);
	peak2fit->Draw("SAME");

	TF1 *fback = new TF1("fback",background_func,low,high,3);
	fback->SetParNames("a0","a1","a2");
	fback->SetNpx(500);
	fback->SetParameters(par1[0],par1[1],par1[2]);
	// fback->SetParErrors(ffit1->GetParErrors());	// Set parameter errors on BG terms from previous fit
	fback->SetLineColor(kCyan);
	fback->Draw("SAME");

	double chi2NDF = ffit1->GetChisquare()/ffit1->GetNDF();

	// results format: [area, area err, chi2NDF, sigma1, linear, offset,
	//					area2, area2 err]
	vector < double > results;
	results.push_back(ffit1->GetParameter(3));
	results.push_back(ffit1->GetParError(3));
	// if(ffit1->GetParError(3) < 0.1) results.push_back(ffit1->GetParError(3));
	// else results.push_back(TMath::Sqrt(ffit1->GetParameter(3)));
	results.push_back(chi2NDF);
	results.push_back(par1[5]);
	results.push_back(ffit1->GetParameter(1));
	results.push_back(ffit1->GetParameter(0));
	results.push_back(ffit1->GetParameter(6));
	results.push_back(ffit1->GetParError(6));

	return results;
}

/*===========================================================================*/
#endif
