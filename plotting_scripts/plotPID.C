{

	TChain* Chain_Coinc = new TChain ("t1");
	TChain* Chain_Singles = new TChain ("t1");

	int firstRun = 60;
	int lastRun = 60;
	char fname[64];

	for (int i=firstRun; i<=lastRun; i++) {

		sprintf(fname,"../root_outputs/Run00%d_Coinc.root", i);
		Chain_Coinc->Add(fname);
	
		sprintf(fname,"../root_outputs/Run00%d_Singles.root", i);
		Chain_Singles->Add(fname);
	}


	TFile *OutFile = new TFile("ParticleID.root", "recreate");

	TCanvas *c1 = new TCanvas("c1");
	TH2D *hPID_Barrel = new TH2D("hPID_Barrel","Particle ID Barrel",2000,0,20000,2000,0,20000);
	Chain_Coinc->Draw("Range_BB10:Energy>>hPID_Barrel","Range_BB10>0","colz");

	TCanvas *c2 = new TCanvas("c2");
	TH2D *hPID_QQQ5 = new TH2D("hPID_QQQ5","Particle ID QQQ5",2000,0,20000,2000,0,20000);
	Chain_Coinc->Draw("Range_QQQ5:Energy>>hPID_QQQ5","Range_QQQ5>0","colz");

	c1->Write();
	c2->Write();
	hPID_Barrel->Write();
	hPID_QQQ5->Write();

}
