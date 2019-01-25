#ifndef FITFUNCTIONS_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define FITFUNCTIONS_H

# include <iostream>

#include <vector>
using std::vector;

#include "TF1.h"
#include "TH1D.h"
#include "TMath.h"


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


double fit_single_gauss_func(double *x, double *par){
    // Fit a Gaussian with a Quadratic bakground

	return background_func(x,par) + single_gauss_func(x,&par[3]);
}


double fit_double_gauss_func(double *x, double *par) {
    // Fit a double Gaussian with a Quadratic bakground, return array of parameters

   return background_func(x,par) + double_gauss_func(x,&par[3]);
}


/*===========================================================================*/
						///////////////////
						// FINDING PEAKS //
						///////////////////

vector < double > iterative_single_gauss_peak(int low, int high, TH1D *H0){
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
	H0->Fit("ffit2","SQRN0");


	double peakPos = ffit2->GetParameter(4);
	double sigma = ffit2->GetParameter(5);

	// results format: [centroid,low,high]
	vector < double > results;
	results.push_back(peakPos);
	results.push_back(low);
	results.push_back(high);
	results.push_back(sigma);

	return results;
}


vector < double > iterative_double_gauss_peak(int low,int high,TH1D *H0){
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
		high = par1[7]+3*par1[8];
	}
	else{
		low = par1[7]-22; // 30;
		high = par1[4]+3*par1[4];
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

	H0->Fit("ffit2","SQRN0");

	double peakPos = ffit2->GetParameter(4);
	double sigma = ffit2->GetParameter(5);

	// Sort peak
	if (ffit2->GetParameter(7) < peakPos){
		peakPos = ffit2->GetParameter(7);
		sigma = ffit2->GetParameter(8);
	}

	// results format: [centroid,low,high,sigma]
	vector < double > results;
	results.push_back(peakPos);
	results.push_back(low);
	results.push_back(high);
	results.push_back(sigma);

	return results;
}


/*===========================================================================*/




/*===========================================================================*/
						///////////////////
						// FITTING PEAKS //
						///////////////////

vector < double > iterative_single_gauss_area(int low, int high, TH1D *H0){

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

	// // Plot markers showing interval of integration
	// TBox *box1 = new TBox(low-5,-5,low,140);
	// TBox *box2 = new TBox(high,-5,high+5,140);
	// box1->SetFillColor(2);
	// box2->SetFillColor(2);
	// box1->Draw("SAME");
	// box2->Draw("SAME");

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
	H0->Fit("ffit2","SQR");


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


	// Plot markers showing interval of integration
	TBox *box3 = new TBox(low-5,-5,low,140);
	TBox *box4 = new TBox(high,-5,high+5,140);
	box3->SetFillColor(3);
	box4->SetFillColor(3);
	box3->Draw("SAME");
	box4->Draw("SAME");


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



vector < double > iterative_double_gauss_area(int low,int high,TH1D *H0){
	/* FOR FINDING PEAK POSITIONS:
		Fit using fit_double_gauss_func, hone in on parameters and then
		fit it again. Return peak 1 area
	*/

	// FIRST ITERATION FIT
	TF1 *ffit1 = new TF1("ffit1",fit_double_gauss_func,low,high,9);
	ffit1->SetParNames("a0","a1","a2","norm1","mean1","sigma1","norm2","mean2","sigma2");
	ffit1->SetNpx(500);
	ffit1->SetParameters(0,0,0,4000,(low+high)/2,12,4000,(low+high)/2+30,12); // -20 , +20

	// Turn off linear background terms
	// ffit1->FixParameter(0,0);
	// ffit1->FixParameter(1,0);
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
		low = par1[4]-30;
		high = par1[7]+3*par1[8];
	}
	else{
		low = par1[7]-30;
		high = par1[4]+3*par1[4];
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

	double area = ffit2->GetParameter(3);
	double area_err = ffit2->GetParError(3);

	// Sort peak
	if (ffit2->GetParameter(7) < ffit2->GetParameter(4)){
		area = ffit2->GetParameter(6);
		area_err = ffit2->GetParError(6);
	}

	double chi2NDF = ffit2->GetChisquare()/ffit2->GetNDF();

	// results format: [area, area err, chi2NDF]
	vector < double > results;
	results.push_back(area);
	results.push_back(area_err);
	results.push_back(chi2NDF);

	return results;
}



/*===========================================================================*/
#endif
