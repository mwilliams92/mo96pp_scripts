{

    TFile *infile = TFile::Open("../root_outputs/Run0060_Coinc.root");

    TH1D *hExc = (TH1D*)infile->Get("hExc");
    hExc->Draw();

    TF1 *f1 = new TF1("f1","gaus",0,5);
    f1->FixParameter(0,253);
    f1->FixParameter(1,2.121);
    f1->FixParameter(2,0.053);
    f1->Draw("same");

    cout << f1->Integral(0,2500)/0.01 << endl;


}
