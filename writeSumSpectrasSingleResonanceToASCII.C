// Read in summed spectras file and write them to ASCII

#include <iostream>
using std::cout;
using std::endl;
using std::ofstream;

#include <ostream>
#include <string>
using std::string;


void writeSumSpectrasSingleResonanceToASCII(){

    TFile *summed = new TFile("E_cal_spectras/summedSpectrasSingleResonance.root","READ");

    // Loop through detectors in files
    for (int detLoop = 0; detLoop <= 12; detLoop += 1){

        // Get histograms from root file
        TH1D *h = static_cast<TH1D*>(summed->Get(Form("h%d",detLoop)));

        // First file
        int n = h->GetNbinsX();
        string outFileName = Form("E_cal_spectras/summedSpectrasSingleResonance_det-%d",detLoop);
    	ofstream myfile (outFileName.c_str(), ios::out);
    	if (myfile.is_open()) {
    		for (int i=1; i<=n; i++) {
    			myfile << h->GetBinLowEdge(i)+h->GetBinWidth(i)/2<<"\t"<<h->GetBinContent(i)<<"\n";
    	    }
    		myfile.close();
    	}
    	else cout << Form("summedSpectrasSingleResonance_det-%d: Messed up\n",detLoop);
    }

    summed->Close();
}
