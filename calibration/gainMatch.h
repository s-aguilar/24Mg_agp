#ifndef GAINMATCH_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define GAINMATCH_H

#include <iostream>

#include <vector>
using std::vector;

#include "TRandom3.h"
#include "TH1D.h"


vector < double > gain_match(Double_t E1, Double_t E2, Double_t E1p, Double_t E2p, TH1D* H2, TH1D* H2p){ // TH1D *gain_match
	// H2p is the original spectrum, H2 will be the new gain matched background spectrum

	TRandom3 *rand = new TRandom3(0); // Random number [0,1)

	// "linear" map
	// E(Ep) = a*Ep + b
	// E2 = 0;
	// E2p = 0;
	double a = (E2 - E1)/(E2p - E1p);
	double b = E2 - a*E2p;

	// cout << "a: " << a << "\t" << "b: " << b << endl;

	// get histograms

	Int_t ch_max = H2p->GetNbinsX();
	// cout << "Number of bins: " << ch_max << endl;
 	// check for correct number of bins

	for(int ch=0; ch<ch_max; ch++){
		int bin_value;
		bin_value = H2p->GetBinContent(ch);

		for(int N=0; N<bin_value; N++){
			H2->Fill(a*(ch+(rand->Rndm())) + b); // fill new histogram with gain matched energies
		}
	}

	// cout << "E1: " << E1 << "\t" << "E2: " << E2 << endl;
	vector < double > results;
	results.push_back(a);
	results.push_back(b);

	return results;
}


TH1D *gain_match2(double a, double b, TH1D* H2, TH1D* H2p){
	// H2p is the original spectrum, H2 will be the new calibrated spectrum

	TRandom3 *rand = new TRandom3(0); // Random number [0,1)

	Int_t ch_max = H2p->GetNbinsX();
	// cout << "Number of bins: " << ch_max << endl;
 	// check for correct number of bins

	for(int ch=0; ch<ch_max; ch++){
		int bin_value;
		bin_value = H2p->GetBinContent(ch);

		for(int N=0; N<bin_value; N++){
			H2->Fill(a*(ch+(rand->Rndm())) + b); // fill new histogram with gain matched energies
		}
	}

	return H2;
}


vector < double > calibrate(Double_t E1, Double_t E2, Double_t E1p, Double_t E2p, TH1D* H2, TH1D* H2p){
	// H2p is the original spectrum, H2 will be the new calibrated spectrum

	TRandom3 *rand = new TRandom3(0); // Random number [0,1)

	// "linear" map
	// E(Ep) = a*Ep + b
	double a = (E2 - E1)/(E2p - E1p);
	double b = E2 - a*E2p;

	// cout << "a: " << a << "\t" << "b: " << b << endl;

	// get histograms

	Int_t ch_max = H2p->GetNbinsX();
	// cout << "Number of bins: " << ch_max << endl;
 	// check for correct number of bins

	for(int ch=0; ch<ch_max; ch++){
		int bin_value;
		bin_value = H2p->GetBinContent(ch);

		for(int N=0; N<bin_value; N++){
			H2->Fill(a*(ch+(rand->Rndm())) + b); // fill new histogram with gain matched energies
		}
	}

	// cout << "E1: " << E1 << "\t" << "E2: " << E2 << endl;
	vector < double > results;
	results.push_back(a);
	results.push_back(b);

	return results;
}


TH1D *calibrate2(double a, double b, TH1D* H2, TH1D* H2p){
	// H2p is the original spectrum, H2 will be the new calibrated spectrum

	TRandom3 *rand = new TRandom3(0); // Random number [0,1)

	// "linear" map
	// E(Ep) = a*Ep + b

	Int_t ch_max = H2p->GetNbinsX();
	// cout << "Number of bins: " << ch_max << endl;
 	// check for correct number of bins

	for(int ch=0; ch<ch_max; ch++){
		int bin_value;
		bin_value = H2p->GetBinContent(ch);

		for(int N=0; N<bin_value; N++){
			H2->Fill(a*(ch+(rand->Rndm())) + b); // fill new histogram with gain matched energies
		}
	}

	return H2;
}

#endif
