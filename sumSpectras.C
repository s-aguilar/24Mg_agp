// Sum the highest energy spectras to address referee comment to check for other
// higher energy gamma-ray transitions

#include <iostream>
#include <list>

void sumSpectras(){

    std::list<int> listOfRuns({159,160,161,162,167,172,
                            260,261,262,263,264,265,266,
                            267,268,269,270,271,272,273,
                            274,275,278});

    TFile *ff = new TFile("E_cal_spectras/summedSpectras.root","RECREATE");

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

    const char *detect;
    const char *files;

    for (int val : listOfRuns){

        // Create the canvas
        TCanvas *c0 = new TCanvas("c0","c0",1920,1080);
        c0->Update();

        TString runNum_TString;

        // Pad the values
		if(val < 100) runNum_TString = "00";
		else if((val >= 100) && (val < 1000)) runNum_TString = "0";
		else runNum_TString = "";

        runNum_TString += val;	// Should be format 0001 -> 9999

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

    // Create the canvas
    // TCanvas *c0 = new TCanvas("c0","c0",1920,1080);
    // c0->Update();
    // gPad->SetLogy();
    // h0->Draw();
    //
    // TCanvas *c1 = new TCanvas("c1","c1",1920,1080);
    // c1->Update();
    // gPad->SetLogy();
    // h1->Draw();
    //
    // TCanvas *c2 = new TCanvas("c2","c2",1920,1080);
    // c2->Update();
    // gPad->SetLogy();
    // h2->Draw();

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
