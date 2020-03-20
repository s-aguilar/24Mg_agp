// Sum all the energy calibrated spectra per detector to address referee
// comment to check for other higher energy gamma-ray transitions

#include <iostream>
#include <list>

void sumSpectrasALL(){

    const char *detect;
    const char *files;

    TFile *ff = new TFile("E_cal_spectras/summedSpectrasALL.root","RECREATE");

    TH1D *h0 = new TH1D("h0","120 deg",8192,0,8192);
    TH1D *h1 = new TH1D("h1","105 deg",8192,0,8192);
    TH1D *h2 = new TH1D("h2","90 deg",8192,0,8192);
    TH1D *h3 = new TH1D("h3","45 deg",8192,0,8192);
    TH1D *h4 = new TH1D("h4","30 deg",8192,0,8192);
    TH1D *h5 = new TH1D("h5","15 deg",8192,0,8192);
    TH1D *h6 = new TH1D("h6","0 deg",8192,0,8192);
    TH1D *h7 = new TH1D("h7","-15 deg",8192,0,8192);
    TH1D *h8 = new TH1D("h8","-30 deg",8192,0,8192);
    TH1D *h9 = new TH1D("h9","-45 deg",8192,0,8192);
    TH1D *h10 = new TH1D("h10","-90 deg",8192,0,8192);
    TH1D *h11 = new TH1D("h11","-105 deg",8192,0,8192);
    TH1D *h12 = new TH1D("h12","-120 deg",8192,0,8192);

    int startRun = 160;
    int upToRun = 410;
	for(int run=startRun;run<upToRun;run++){

        // Skip bad runs
        if(run>=163 && run<=166) continue;
        else if(run>=168 && run<=171) continue;
        else if(run==182) continue;
        else if(run>=244 && run<=255) continue;
        else if(run==276) continue;
        else if(run==277) continue;
        else if(run==285) continue;
        else if(run==289) continue;
        else if(run==290) continue;
        else if(run==294) continue;
        else if(run==406) continue;

        // Create the canvas
        TCanvas *c0 = new TCanvas("c0","c0",1920,1080);
        c0->Update();

        TString runNum_TString;

        // Pad the values
		if(run < 100) runNum_TString = "00";
		else if((run >= 100) && (run < 1000)) runNum_TString = "0";
		else runNum_TString = "";

        runNum_TString += run;	// Should be format 0001 -> 9999

        // Cast to const char*
		const char *runNum_String = (const char*)runNum_TString;

        // Open file
        TFile *fyield = new TFile(Form("E_cal_spectras/run_%s.root",runNum_String),"READ");

        // Loop through detectors and sum their spectras with the previous runs
        for (int det = 0; det <= 12; det++){
            detect = Form("det-%d",det);
            if(det == 0) h0->Add(static_cast<TH1D*>(fyield->Get(detect)),1);
            else if(det == 1) h1->Add(static_cast<TH1D*>(fyield->Get(detect)),1);
            else if(det == 2) h2->Add(static_cast<TH1D*>(fyield->Get(detect)),1);
            else if(det == 3) h3->Add(static_cast<TH1D*>(fyield->Get(detect)),1);
            else if(det == 4) h4->Add(static_cast<TH1D*>(fyield->Get(detect)),1);
            else if(det == 5) h5->Add(static_cast<TH1D*>(fyield->Get(detect)),1);
            else if(det == 6) h6->Add(static_cast<TH1D*>(fyield->Get(detect)),1);
            else if(det == 7) h7->Add(static_cast<TH1D*>(fyield->Get(detect)),1);
            else if(det == 8) h8->Add(static_cast<TH1D*>(fyield->Get(detect)),1);
            else if(det == 9) h9->Add(static_cast<TH1D*>(fyield->Get(detect)),1);
            else if(det == 10) h10->Add(static_cast<TH1D*>(fyield->Get(detect)),1);
            else if(det == 11) h11->Add(static_cast<TH1D*>(fyield->Get(detect)),1);
            else if(det == 12) h12->Add(static_cast<TH1D*>(fyield->Get(detect)),1);
        }


        c0->Clear();
        fyield->Close();
        delete fyield;
        delete c0;
    }

    // // // Write all TObjects in memory (TFitResult) to TFile
    ff->Write();
    delete h0;
    delete h1;
    delete h2;
    delete h3;
    delete h4;
    delete h5;
    delete h6;
    delete h7;
    delete h8;
    delete h9;
    delete h10;
    delete h11;
    delete h12;
    ff->Close();


}
