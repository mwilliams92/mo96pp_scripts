{

	TChain* Chain_Coinc = new TChain ("t1");
	TChain* Chain_Singles = new TChain ("t1");

	int firstRun = 60;
	int lastRun = 78;
	char fname[64];

	for (int i=firstRun; i<=lastRun; i++) {

		sprintf(fname,"../root_outputs/Run00%d_Coinc.root", i);
		Chain_Coinc->Add(fname);

		sprintf(fname,"../root_outputs/Run00%d_Singles.root", i);
		Chain_Singles->Add(fname);
	}


	TFile *OutFile = new TFile("Exc_v_Gamma.root", "recreate");

	TCut Proton_Cut = "(Range_QQQ5>2700 && Range_QQQ5<4200 && Energy>Range_QQQ5) || (Range_BB10>2300 && Range_BB10<2900 && Energy>Range_BB10)";
	TCut Proton_Cut_Barrel = "Range_BB10>2300 && Range_BB10<2900";
	TCut Proton_Cut_QQQ5 = "Range_QQQ5>2700 && Range_QQQ5<4200 && Energy>Range_QQQ5 && Energy<12000";

	TCanvas *c1 = new TCanvas();
	TH2D *hExc_Egam = new TH2D("hExc_Egam","Excitation vs Gamma Energy",10000,0,10000,2000,-5,15);
	Chain_Coinc->Draw("Excitation:Egamma>>hExc_Egam",Proton_Cut_Barrel,"colz");

	TCanvas *c2 = new TCanvas();
	TH1D *hExc = new TH1D("hExc","Excitation Energy Barrel",2000,-5,15);
	Chain_Singles->Draw("Excitation>>hExc",Proton_Cut_Barrel);
	hExc->SetLineColor(1);

	hExc->Write();
	hExc_Egam->Write();
	OutFile->Write();


}
