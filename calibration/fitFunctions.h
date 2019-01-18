#ifndef FITFUNCTIONS_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define FITFUNCTIONS_H

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

double iterative_single_gauss_peak(int low, int high, TH1D *H0){
	/* FOR FINDING PEAK POSITIONS:
		Fit using fit_single_gauss_func, hone in on parameters and then
		fit it again.
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

	return peakPos;
}


double iterative_double_gauss_peak(int low,int high,TH1D *H0){
	/* FOR FINDING PEAK POSITIONS:
		Fit using fit_single_gauss_func, hone in on parameters and then
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

	// Sort to higher energy peak I as the second peak
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
	low = par1[4]-2.5*par1[5]; //
	high = par1[7]+3*par1[8];


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

	// Sort peak
	if (ffit2->GetParameter(7) > peakPos){
		peakPos = ffit2->GetParameter(7);
	}

	return peakPos;
}
/*===========================================================================*/




// /*===========================================================================*/
// 						///////////////////
// 						// FITTING PEAKS //
// 						///////////////////
//
// double *iterative_single_gauss_fit(int low, int high, TH1D *H0){
// 	/* FOR FINDING PEAK POSITIONS:
// 		Fit using fit_single_gauss_func, hone in on parameters and then
// 		fit it again.
// 		*/
//
// 	// FIRST ITERATION FIT
// 	TF1 *ffit1 = new TF1("ffit1",fit_single_gauss_func,low,high,6);
// 	ffit1->SetParNames("a0","a1","a2","norm","mean","sigma");
// 	ffit1->SetNpx(500);
// 	ffit1->SetParameters(1,1,1,4000,(low+high)/2,12);
// 	ffit1->FixParameter(2,0);		// Makes it a linear background
//     ffit1->SetParLimits(3,0,1e6);
//     ffit1->SetParLimits(4,low,high);	// Peak centroid range
//     ffit1->SetParLimits(5,2.5,55);	// Std dev range
//
// 	H0->Fit("ffit1","SQRN0");
//
// 	// Store fit parameters
// 	double par1[6];
// 	ffit1->GetParameters(par1);
//
// 	// Recalibrate range of integration
// 	low = par1[4]-3*par1[5];
// 	high = par1[4]+3*par1[5];
//
//
// 	// SECOND ITERATION FIT
// 	TF1 *ffit2 = new TF1("ffit2",fit_single_gauss_func,low,high,6);
// 	ffit2->SetParNames("a0","a1","a2","norm","mean","sigma");
// 	ffit2->SetNpx(500);
// 	ffit2->SetParameters(par1);
// 	ffit2->FixParameter(2,0);		// Makes it a linear background
// 	ffit2->SetParLimits(3,0,1e6);
// 	ffit2->SetParLimits(4,low,high);	// Peak centroid range
// 	ffit2->SetParLimits(5,2.5,55);	// Std dev range
// 	H0->Fit("ffit2","SQRN0");
//
// 	double par2[6];
// 	double par2_err[6];
//
// 	ffit2->GetParameters(par2);
// 	ffit2->GetParErrors(par2_err);
//
// 	double nested_results[2];
// 	nested_results[0] = {par2};
// 	nested_results[1] = {par2_err};
// 	return nested_results;
// }
//

// double *iterative_double_gauss_fit(int low, int high, TH1D *H0){
// 	// Fit using fit_double_gauss_func, hone in on parameters and then\
// 		fit it again.
//
// 	double blah[1];
// 	blah[0]=1;
// 	return blah;
// }




/*===========================================================================*/
#endif
