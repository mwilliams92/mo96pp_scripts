{
	gROOT->SetBatch(kTRUE);

	TFile *infile1 = TFile::Open("../OutputFolder/cal228th_315deg_US_BR.root");
	TTree *data1 = (TTree*)infile1->Get("data");

	TFile *infile2 = TFile::Open("../OutputFolder/cal228th_45deg_US_BL.root");
	TTree *data2 = (TTree*)infile2->Get("data");

	TFile *infile3 = TFile::Open("../OutputFolder/cal228th_270deg_BR_noBB10.root");
	TTree *data3 = (TTree*)infile3->Get("data");

	TFile *infile4 = TFile::Open("../OutputFolder/cal228th_90deg_BL_noBB10.root");
	TTree *data4 = (TTree*)infile4->Get("data");

	TFile *infile5 = TFile::Open("../OutputFolder/cal228th_180deg_DS_noBB10.root");
	TTree *data5 = (TTree*)infile5->Get("data");

	char cut[4096];
	char draw[4096];
	char hname[4096];

	TFile *outfile = TFile::Open("FrontEnergyVsStrip.root","RECREATE");
	outfile->cd();

	TH2D *hFrontEnergyUpstream[12];
	TH2D *hFrontEnergyDownstream[12];

	TDirectory *dirUp = outfile->mkdir("Upstream");
	TDirectory *dirDown = outfile->mkdir("Downstream");

	dirUp->cd();

	for( int i=0; i<6; i++ )	{

		cout << "Detector " << i << endl;

		for( int i=0; i<4; j++ )	{

			sprintf(hname,"hFrontEnergyUpstream[%d]", i);
			sprintf(draw,"SX3StripEnergy:SX3Strip>>hFrontEnergyUpstream[%d]", i);
			sprintf(cut,"SX3Upstream==1 && SX3Det==%d", i);
			
			TCut Cut = cut;

			hFrontEnergyUpstream[i] = new TH2D(hname, hname, 4, 0, 4, 4096, 0, 4096);

			data1->Draw(draw,cut,"same");

		}

	}

	for( int i=6; i<12; i++ )	{

		cout << "Detector " << i << endl;

			sprintf(hname,"hFrontEnergyUpstream[%d]", i);
			sprintf(draw,"SX3StripEnergy:SX3Strip>>hFrontEnergyUpstream[%d]", i);
			sprintf(cut,"SX3Upstream==1 && SX3Det==%d", i);
			
			TCut Cut = cut;

			hFrontEnergyUpstream[i] = new TH2D(hname, hname, 4, 0, 4, 4096, 0, 4096);

			data2->Draw(draw,cut,"same");

	}

	dirDown->cd();

	for( int i=0; i<6; i++ )	{

		cout << "Detector " << i << endl;

			sprintf(hname,"hFrontEnergyDownstream[%d]", i);
			sprintf(draw,"SX3StripEnergy:SX3Strip>>hFrontEnergyDownstream[%d]", i);
			sprintf(cut,"SX3Upstream==0 && SX3Det==%d", i);
			
			TCut Cut = cut;

			hFrontEnergyDownstream[i] = new TH2D(hname, hname, 4, 0, 4, 4096, 0, 4096);

			data3->Draw(draw,cut,"same");

	}

	for( int i=6; i<12; i++ )	{

		cout << "Detector " << i << endl;

			sprintf(hname,"hFrontEnergyDownstream[%d]", i);
			sprintf(draw,"SX3StripEnergy:SX3Strip>>hFrontEnergyDownstream[%d]", i);
			sprintf(cut,"SX3Upstream==0 && SX3Det==%d", i);
			
			TCut Cut = cut;

			hFrontEnergyDownstream[i] = new TH2D(hname, hname, 4, 0, 4, 4096, 0, 4096);

			data4->Draw(draw,cut,"same");

	}

	outfile->cd();
	outfile->Write();
	outfile->Close();	

}
