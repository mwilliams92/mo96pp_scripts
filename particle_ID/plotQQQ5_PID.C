{

	TChain* Chain = new TChain ("data");
    Chain->Add("../OutputFolder/Run0031.root");

	TCanvas *c1 = new TCanvas();

	TH2D *hPID;

	//char cut[4096];
	//char draw[4096];
	//char hname[4096];

	//for(int i=0; i<12; i++) {

		//c1->cd(i+1);

		sprintf(hname,"hPID[%d]", i);
		sprintf(draw,"BB10Energy:SX3SectorEnergy>>hPID[%d]", i);
		sprintf(cut,"BB10Mul>0 && BB10Det==%d", i);

		hdeltaE = new TH1D(hname, hname,5000,0,50000);

		deltaE

		Chain->Draw(draw,cut,"colz");
		Chain->Draw("QQQ5SectorEnergy>>hdeltaE","QQQ5Detector==0");


	//}


}
