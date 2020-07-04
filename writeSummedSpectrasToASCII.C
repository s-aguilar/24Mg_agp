// Read in summed spectras file and write them to ASCII

#include <iostream>
using std::cout;
using std::endl;
using std::ofstream;

#include <ostream>
#include <string>
using std::string;


void writeSummedSpectrasToASCII(){

    TFile *summedAll = new TFile("E_cal_spectras/summedSpectrasALL.root","READ");
    TFile *summedBGsubAll = new TFile("E_cal_spectras/summedSpectrasBGsubALL.root","READ");

    // Loop through detectors in files
    for (int detLoop = 0; detLoop <= 12; detLoop += 1){

        // Get histograms from root file
        TH1D *h = static_cast<TH1D*>(summedAll->Get(Form("h%d",detLoop)));

        // First file
        int n = h->GetNbinsX();
        string outFileName = Form("E_cal_spectras/summedSpectrasALL_det-%d",detLoop);
    	ofstream myfile (outFileName.c_str(), ios::out);
    	if (myfile.is_open()) {
    		for (int i=1; i<=n; i++) {
    			myfile << h->GetBinLowEdge(i)+h->GetBinWidth(i)/2<<"\t"<<h->GetBinContent(i)<<"\n";
    	    }
    		myfile.close();
    	}
    	else cout << Form("summedSpectrasALL_det-%d: Messed up\n",detLoop);

        // Other file
        TH1D *hh = static_cast<TH1D*>(summedBGsubAll->Get(Form("h%d",detLoop)));
        int nn = hh->GetNbinsX();
        outFileName = Form("E_cal_spectras/summedSpectrasBGsubALL_det-%d",detLoop);
    	ofstream myfile1 (outFileName.c_str(), ios::out);
    	if (myfile1.is_open()) {
    		for (int i=1; i<=nn; i++) {
    			myfile1 << hh->GetBinLowEdge(i)+hh->GetBinWidth(i)/2<<"\t"<<hh->GetBinContent(i)<<"\n";
    	    }
    		myfile1.close();
    	}
    	else cout << Form("summedSpectrasBGsubALL_det-%d: Messed up\n",detLoop);
    }

    summedAll->Close();
    summedBGsubAll->Close();
}
