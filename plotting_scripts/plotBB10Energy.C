{
	gROOT->SetBatch(kTRUE);

	TFile *infile1 = TFile::Open("../OutputFolder/cal228th_225deg_DS_BR.root");
	TTree *data1 = (TTree*)infile1->Get("data");

	TFile *infile2 = TFile::Open("../OutputFolder/cal228th_135deg_DS_BL.root");
	TTree *data2 = (TTree*)infile2->Get("data");

	char cut[4096];
	char draw[4096];
	char hname[4096];

	TFile *outfile = TFile::Open("BB10Energy.root","RECREATE");
	outfile->cd();

	TH1D *hBB10Energy[12][8];

	for( int i=0; i<5; i++ )	{

		cout << "Detector " << i << endl;

		for ( int j=0; j<8; j++ )	{

			sprintf(hname,"hBB10Energy[%d][%d]", i, j);
			sprintf(draw,"BB10Energy>>hBB10Energy[%d][%d]", i, j);
			sprintf(cut,"BB10Det==%d && BB10Strip==%d", i, j);
			
			TCut Cut = cut;

			hBB10Energy[i][j] = new TH1D(hname, hname, 100, 0, 10000);

			data1->Draw(draw,cut,"same");

		}

	}

	for( int i=7; i<12; i++ )	{

		cout << "Detector " << i << endl;

			for ( int j=0; j<8; j++ )	{

			sprintf(hname,"hBB10Energy[%d][%d]", i, j);
			sprintf(draw,"BB10Energy>>hBB10Energy[%d][%d]", i, j);
			sprintf(cut,"BB10Det==%d && BB10Strip==%d", i, j);
			
			TCut Cut = cut;

			hBB10Energy[i][j] = new TH1D(hname, hname, 100, 0, 10000);

			data2->Draw(draw,cut,"same");

		}

	}

	outfile->cd();
	outfile->Write();
	outfile->Close();	

}
