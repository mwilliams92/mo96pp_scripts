{

    TFile *mo96data = TFile::Open("../root_outputs/Kinematics_Histograms_Run0060.root");
    TFile *CD2data = TFile::Open("../root_outputs/Kinematics_Histograms_Run0080.root");

    TCanvas *c1 = new TCanvas();
    
    TH1D *hExc_mo96 = (TH1D*)mo96data->Get("hExcitation_US");
    hExc_mo96->Rebin(2);
    hExc_mo96->Draw();

    TH1D *hExc_CD2 = (TH1D*)CD2data->Get("hExcitation_US");
    hExc_CD2->Rebin(2);
    //hExc_CD2->Scale(0.0653);
    hExc_CD2->Draw("HIST+same");

    TCanvas *c2 = new TCanvas();

    TH1F *hSub = new TH1F("hSub","Subtraction",1250,-5,20);

    hSub->Add(hExc_mo96);
    hSub->Add(hExc_CD2,-1);
    hSub->Draw();

    TCanvas *c3 = new TCanvas();
    
    TH1D *hExc_mo96_DS = (TH1D*)mo96data->Get("hExcitation_DS");
    hExc_mo96_DS->Rebin(2);
    hExc_mo96_DS->Draw();

    TH1D *hExc_CD2_DS = (TH1D*)CD2data->Get("hExcitation_DS");
    hExc_CD2_DS->Rebin(2);
    //hExc_CD2_DS->Scale(0.07);
    hExc_CD2_DS->Draw("HIST+same");


}
