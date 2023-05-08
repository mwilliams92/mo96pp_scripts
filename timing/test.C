{

    //TFile *raw = TFile::Open("/mnt/sandisk/files/sorted/Run0060_combined.root");
    TFile *raw = TFile::Open("/mnt/d/sorted/Run0060_combined.root");
    TTree *Traw = (TTree*)raw->Get("data");

    TCanvas *c1 = new TCanvas();

    Traw->Draw("timeStamp-xtals_timestamp>>h1(2000,-1000,1000)");

    TFile *sorted = TFile::Open("Timing_Histograms.root");
    TTree *Tsorted = (TTree*)sorted->Get("t1");

    TH1F *h2 = (TH1F*)sorted->Get("hTS_Diff");
    h2->SetLineColor(2);
    h2->Draw("same");

    TH1F *h3 = new TH1F("h3","ORRUBA-GRETINA TS Difference",2000,-1000,1000);

    //Tsorted->SetAlias("delta_TS","timeStamp-timeStamp_GT");

    Tsorted->Draw("TimeStampDiff>>h3","","same");
    h3->SetLineColor(6);
    h3->Draw("same");
    h3->SetLineColor(6);



}
