{

	gROOT->SetBatch(kTRUE);

	TChain* Chain = new TChain ("data");
/*  Chain->Add("../OutputFolder/Run0030.root");
	Chain->Add("../OutputFolder/Run0031.root");
	Chain->Add("../OutputFolder/Run0032.root");
	Chain->Add("../OutputFolder/Run0033.root");
	Chain->Add("../OutputFolder/Run0034.root");
	Chain->Add("../OutputFolder/Run0035.root");
	Chain->Add("../OutputFolder/Run0037.root");
	Chain->Add("../OutputFolder/Run0038.root");
	Chain->Add("../OutputFolder/Run0039.root");
	Chain->Add("../OutputFolder/Run0040.root");
	Chain->Add("../OutputFolder/Run0041.root");
	Chain->Add("../OutputFolder/Run0042.root");
	Chain->Add("../OutputFolder/Run0043.root");
	Chain->Add("../OutputFolder/Run0044.root");
	Chain->Add("../OutputFolder/Run0045.root");
	Chain->Add("../OutputFolder/Run0046.root");
	Chain->Add("../OutputFolder/Run0047.root");
	Chain->Add("../OutputFolder/Run0048.root");
	Chain->Add("../OutputFolder/Run0049.root");
*/
	Chain->Add("../../OutputFolder/Run0064.root");

	Chain->SetAlias("Etotal","BB10Energy+SX3SectorEnergy");
	Chain->SetAlias("Range","pow((Etotal-SX3SectorEnergy),1.73)");

	TFile *root_outfile = new TFile("particleID.root", "RECREATE");

	TH2D *hPID_Indv[12][4];
	TH2D *hPID[4];

	char cut[4096];
	char draw[4096];
	char hname[4096];

	// Plot each individual sector for each detector

	for(int i=0; i<12; i++) {
		for (int j=0; j<4; j++) {

			sprintf(hname,"hPID_Indv[%d][%d]", i, j);
			sprintf(draw,"Range:Etotal>>hPID_Indv[%d][%d]", i, j);
			sprintf(cut,"BB10Mul>0 && SX3Upstream==0 && BB10Energy>250 && BB10Det==%d && SX3Det==%d && SX3Sector==%d", i, i, j);

			hPID_Indv[i][j] = new TH2D(hname,hname,3000,0,30000,1000,0,10000);
			Chain->Draw(draw,cut,"goff");
		}
	}

	// Plot each individual sector including all detectors
	for (int j=0; j<4; j++) {

			sprintf(hname,"hPID[%d]", j);
			sprintf(draw,"Range:Etotal>>hPID[%d]",j);
			sprintf(cut,"BB10Mul>0 && SX3Upstream==0 && BB10Energy>250 && SX3Sector==%d", j);

			hPID[j] = new TH2D(hname,hname,3000,0,30000,1000,0,10000);
			Chain->Draw(draw,cut,"goff");
			
	}
	

/*
	TCanvas *c1 = new TCanvas();
	c1->Divide(4,3);

	for(int i=0; i<12; i++) {

		c1->cd(i+1);

		sprintf(hname,"hPID[%d]", i);
		sprintf(draw,"BB10Energy:SX3SectorEnergy>>hPID[%d]", i);
		sprintf(cut,"BB10Mul>0 && BB10Det==%d && SX3Upstream==0 && SX3Det==%d && SX3Sector==1", i);

		hPID[i] = new TH2D(hname, hname,2000,0,20000,1000,0,10000);

		Chain->Draw(draw,cut,"colz");


	}

	TCanvas *c2 = new TCanvas();

	TH2D *hPID_2 = new TH2D("hPID_2","BB10 vs SX3Sector Detector 1",3000,0,30000,1000,0,10000);
	Chain->Draw("BB10Energy:Etot>>hPID_2","BB10Mul>0 && BB10Det<5 && SX3Upstream==0 && SX3Sector==1","colz");


	TCanvas *c3 = new TCanvas();
	c3->Divide(4,3);

	for(int i=0; i<12; i++) {

		c3->cd(i+1);
		sprintf(hname,"hPID_Sector0[%d]", i);
		sprintf(draw,"BB10Energy:Etot>>hPID_Sector0[%d]", i);
		sprintf(cut,"BB10Mul>0 && BB10Det==%d && SX3Upstream==0 && SX3Det==%d && SX3Sector==0", i, i);

		hPID_Sector0[i] = new TH2D(hname,hname,3000,0,30000,1000,0,10000);
		Chain->Draw(draw,cut,"colz");
		root_outfile->Append(hPID_Sector0[i]);
	}

	TCanvas *c4 = new TCanvas();
	c4->Divide(4,3);

	for(int i=0; i<12; i++) {

		c4->cd(i+1);
		sprintf(hname,"hPID_Sector1[%d]", i);
		sprintf(draw,"BB10Energy:Etot>>hPID_Sector1[%d]", i);
		sprintf(cut,"BB10Mul>0 && BB10Det==%d && SX3Upstream==0 && SX3Det==%d && SX3Sector==1", i, i);

		hPID_Sector1[i] = new TH2D(hname,hname,3000,0,30000,1000,0,10000);
		Chain->Draw(draw,cut,"colz");
		root_outfile->Append(hPID_Sector1[i]);
	}

	TCanvas *c5 = new TCanvas();
	c5->Divide(4,3);

	for(int i=0; i<12; i++) {

		c5->cd(i+1);
		sprintf(hname,"hPID_Sector2[%d]", i);
		sprintf(draw,"BB10Energy:Etot>>hPID_Sector2[%d]", i);
		sprintf(cut,"BB10Mul>0 && BB10Det==%d && SX3Upstream==0 && SX3Det==%d && SX3Sector==2", i, i);

		hPID_Sector2[i] = new TH2D(hname,hname,3000,0,30000,1000,0,10000);
		Chain->Draw(draw,cut,"colz");
		root_outfile->Append(hPID_Sector2[i]);
	}

	TCanvas *c6 = new TCanvas();
	c6->Divide(4,3);

	for(int i=0; i<12; i++) {

		c6->cd(i+1);
		sprintf(hname,"hPID_Sector3[%d]", i);
		sprintf(draw,"BB10Energy:Etot>>hPID_Sector3[%d]", i);
		sprintf(cut,"BB10Mul>0 && BB10Det==%d && SX3Upstream==0 && SX3Det==%d && SX3Sector==3", i, i);

		hPID_Sector3[i] = new TH2D(hname,hname,3000,0,30000,1000,0,10000);
		Chain->Draw(draw,cut,"colz");
		root_outfile->Append(hPID_Sector3[i]);
	}
*/
	root_outfile->Write();
	root_outfile->Close();	


}
