{

	//gROOT->SetBatch(kTRUE);

	TChain* Chain = new TChain ("data");

	Chain->Add("../OutputFolder/Run0080.root"); // CD2

	//TFile *root_outfile = new TFile("QQQ5down.root", "RECREATE");

	TCanvas *c1 = new TCanvas();
	TH2D *hE2 = new TH2D("hE2","E2 Layer Detector",20000,0,20000,10,0,10);
	Chain->Draw("SX3Mul:SX3SectorEnergy>>hE2","QQQ5Upstream==0 && QQQ5Det==5", "colz");

	TCanvas *c2 = new TCanvas();
	TH2D *hE1 = new TH2D("hE1","E1 Layer Detector",20000,0,20000,10,0,10);
	Chain->Draw("SX3Mul:SX3SectorEnergy>>hE1","QQQ5Upstream==0 && QQQ5Det==3", "colz");

	TCanvas *c3 = new TCanvas();
	TH2D *hdE = new TH2D("hdE","dE Layer Detector",20000,0,20000,10,0,10);
	Chain->Draw("SX3Mul:SX3SectorEnergy>>hdE","QQQ5Upstream==0 && QQQ5Det==1", "colz");



}
