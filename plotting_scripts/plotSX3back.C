{

	TFile *infile1 = TFile::Open("../OutputFolder/cal228th_315deg_US_BR.root");
	TTree *data1 = (TTree*)infile1->Get("data");

	TFile *infile2 = TFile::Open("../OutputFolder/cal228th_45deg_US_BL.root");
	TTree *data2 = (TTree*)infile2->Get("data");

	TFile *infile3 = TFile::Open("../OutputFolder/cal228th_270deg_BR_noBB10.root");
	TTree *data3 = (TTree*)infile3->Get("data");

	TFile *infile4 = TFile::Open("../OutputFolder/cal228th_90deg_BL_noBB10.root");
	TTree *data4 = (TTree*)infile4->Get("data");

	char cut[64];
	char draw[64];
	char hname[64];

	TFile *outfile = TFile::Open("SX3_Back_Energy.root","RECREATE");
	outfile->cd();
	TH1D *hSX3backUpstream[12][4];
	TH1D *hSX3backDownstream[12][4];

	TDirectory *dirUp = outfile->mkdir("Upstream");
	TDirectory *dirDown = outfile->mkdir("Downstream");

	dirUp->cd();

	for( int i=0; i<6; i++ )	{

		cout << "Detector " << i << endl;

		for ( int j=0; j<4; j++ )	{
		
			sprintf(hname,"hSX3backUpstream[%d][%d]", i, j);
			sprintf(draw,"SX3SectorEnergy>>hSX3backUpstream[%d][%d]", i, j);
			sprintf(cut,"SX3Upstream==1 && SX3Det==%d && SX3Sector==%d", i, j);
			TCut Cut = cut;

			hSX3backUpstream[i][j] = new TH1D(hname, hname, 1000, 0, 10000);
	
			data1->Draw(draw,cut,"goff");


		}

	}

	for( int i=6; i<12; i++ )	{

		cout << "Detector " << i << endl;

		for ( int j=0; j<4; j++ )	{
		
			sprintf(hname,"hSX3backUpstream[%d][%d]", i, j);
			sprintf(draw,"SX3SectorEnergy>>hSX3backUpstream[%d][%d]", i, j);
			sprintf(cut,"SX3Upstream==1 && SX3Det==%d && SX3Sector==%d", i, j);
			TCut Cut = cut;

			hSX3backUpstream[i][j] = new TH1D(hname, hname, 1000, 0, 10000);
	
			data2->Draw(draw,cut,"goff");


		}

	}

	dirDown->cd();

	for( int i=0; i<5; i++ )	{

		cout << "Detector " << i << endl;

		for ( int j=0; j<4; j++ )	{
		
			sprintf(hname,"hSX3backDownstream[%d][%d]", i, j);
			sprintf(draw,"SX3SectorEnergy>>hSX3backDownstream[%d][%d]", i, j);
			sprintf(cut,"SX3Upstream==0 && SX3Det==%d && SX3Sector==%d", i, j);
			TCut Cut = cut;

			hSX3backDownstream[i][j] = new TH1D(hname, hname, 1000, 0, 10000);
	
			data3->Draw(draw,cut,"goff");


		}

	}

	for( int i=5; i<12; i++ )	{

		cout << "Detector " << i << endl;

		for ( int j=0; j<4; j++ )	{
		
			sprintf(hname,"hSX3backDownstream[%d][%d]", i, j);
			sprintf(draw,"SX3SectorEnergy>>hSX3backDownstream[%d][%d]", i, j);
			sprintf(cut,"SX3Upstream==0 && SX3Det==%d && SX3Sector==%d", i, j);
			TCut Cut = cut;

			hSX3backDownstream[i][j] = new TH1D(hname, hname, 1000, 0, 10000);
	
			data4->Draw(draw,cut,"goff");

		}

	}

	outfile->cd();
	outfile->Write();
	outfile->Close();	

}
