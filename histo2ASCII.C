#include <iostream>
using std::cout;
using std::endl;
using std::ofstream;

#include <ostream>
#include <fstream>

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



/*==============================MAIN=========================================*/
void histo2ASCII(){

	// Suppress certain root print statements
	gErrorIgnoreLevel = kWarning;

	// When running on CRC
	// const char *path = "/afs/crc.nd.edu/user/s/saguilar/Group/24Mg_ap/";
	// chdir(path);

	// Record start time
	auto start = std::chrono::high_resolution_clock::now();

	TFile *f = new TFile("for James/run0093.root","READ");

	const char* detect;

	// Loop through detectors on board
	for(int j=0;j<13;j++){ // 13

		if (j<8){
			detect = Form("h0-%d",j);
		}
		else{
			detect = Form("h1-%d",j-8);
		}

		// Get histograms from root file
        TH1D *h = static_cast<TH1D*>(f->Get(detect));

		// First file
        int n = h->GetNbinsX();
        string outFileName = Form("for James/run0093_det-%d",j);
    	ofstream myfile (outFileName.c_str(), ios::out);
    	if (myfile.is_open()) {
    		for (int i=1; i<=n; i++) {
    			myfile << h->GetBinLowEdge(i)+h->GetBinWidth(i)/2<<"\t"<<h->GetBinContent(i)<<"\n";
    	    }
    		myfile.close();
    	}
    	else cout << Form("run0093_det-%d Messed up\n",j);


	}

	// // Write all TObjects in memory (TFitResult) to TFile
	f->Close();
	delete f;

	// Record end time
	auto finish = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> elapsed = finish - start;

	cout << "Elapsed time: " << elapsed.count() << " s\n";
}
