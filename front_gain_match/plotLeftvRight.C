{
	gROOT->SetBatch(kTRUE);

	TFile *infile1 = TFile::Open("../../OutputFolder/cal228th_315deg_US_BR.root");
	TTree *data1 = (TTree*)infile1->Get("data");

	TFile *infile2 = TFile::Open("../../OutputFolder/cal228th_45deg_US_BL.root");
	TTree *data2 = (TTree*)infile2->Get("data");

	TFile *infile3 = TFile::Open("../../OutputFolder/cal228th_180deg_DS_noBB10.root");
	TTree *data3 = (TTree*)infile3->Get("data");

	TFile *infile4 = TFile::Open("../../OutputFolder/cal228th_90deg_BL_noBB10.root");
	TTree *data4 = (TTree*)infile4->Get("data");

	char cut[4096];
	char draw[4096];
	char hname[4096];

	TFile *outfile = TFile::Open("LeftvRight.root","RECREATE");
	outfile->cd();

	ofstream txt_outfile;
	txt_outfile.open ("SX3_Strip_Gains.txt");

	TH2D *hLeftvRightUpstream[12][4][3];
	TH2D *hLeftvRightDownstream[12][4][3];

	TDirectory *dirUp = outfile->mkdir("Upstream");
	TDirectory *dirDown = outfile->mkdir("Downstream");

	int alpha_cut_low[3] = {5550, 6600, 8600};
	int alpha_cut_high[3] = {5800, 7000, 9000};

	int alpha_range_low[3][4];
	int alpha_range_high[3][4];

	for (int i=0; i<4; i++) {

		alpha_range_low[0][i] = 400 + i*20;
		alpha_range_low[1][i] = 450 + i*20;
		alpha_range_low[2][i] = 600 + i*20;

		alpha_range_high[0][i] = 1300; //+ i*20;
		alpha_range_high[1][i] = 1550; //+ i*20;
		alpha_range_high[2][i] = 2100; //+ i*50;

	}

	TF1 *alpha_fit[12][4][3];

	Color_t colour[3] = {kRed, kBlue, kGreen};

	dirUp->cd();

	double par[2];
	double sum_par;

	for( int i=2; i<3; i++ )	{

		cout << "Detector " << i << endl;

		for ( int j=0; j<4; j++ )	{

			sum_par=0.;

			for (int k=0; k<3; k++)	{

				sprintf(hname,"hLeftvRightUpstream[%d][%d][%d]", i, j, k);
				sprintf(draw,"SX3StripLeftEnergy:SX3StripRightEnergy>>hLeftvRightUpstream[%d][%d][%d]", i, j, k);
				sprintf(cut,"SX3Upstream==1 && SX3Det==%d && SX3Strip==%d && SX3SectorEnergy>%d && SX3SectorEnergy<%d", i, j, alpha_cut_low[k], alpha_cut_high[k]);
			
				TCut Cut = cut;

				hLeftvRightUpstream[i][j][k] = new TH2D(hname, hname, 4096, 0, 4096, 4096, 0, 4096);
				data1->Draw(draw,cut,"same");

				alpha_fit[i][j][k] = new TF1("f1","[0]+[1]*x",alpha_range_low[k][j],alpha_range_high[k][j]);
   				alpha_fit[i][j][k]->SetParameters(2500.,-1.);
   				alpha_fit[i][j][k]->SetLineColor(colour[k]);
   				hLeftvRightUpstream[i][j][k]->Fit(alpha_fit[i][j][k],"QR+");
   				alpha_fit[i][j][k]->Draw("same");

				alpha_fit[i][j][k]->GetParameters(par);

				sum_par += par[1];

			}

			txt_outfile << i << "\t" << j << "\t" << abs(sum_par/3.) << endl;

		}

	}
/*
	for( int i=6; i<12; i++ )	{

		cout << "Detector " << i << endl;

			for ( int j=0; j<4; j++ )	{

				sum_par=0.;

				for (int k=0; k<3; k++)	{

					sprintf(hname,"hLeftvRightUpstream[%d][%d][%d]", i, j, k);
					sprintf(draw,"SX3StripLeftEnergy:SX3StripRightEnergy>>hLeftvRightUpstream[%d][%d][%d]", i, j, k);
					sprintf(cut,"SX3Upstream==1 && SX3Det==%d && SX3Strip==%d && SX3SectorEnergy>%d && SX3SectorEnergy<%d", i, j, alpha_cut_low[k], alpha_cut_high[k]);
				
					TCut Cut = cut;

					hLeftvRightUpstream[i][j][k] = new TH2D(hname, hname, 4096, 0, 4096, 4096, 0, 4096);
					data2->Draw(draw,cut,"same");

					alpha_fit[i][j][k] = new TF1("f1","[0]+[1]*x",alpha_range_low[k],alpha_range_high[k]);
   					alpha_fit[i][j][k]->SetParameters(2500.,-1.);
   					alpha_fit[i][j][k]->SetLineColor(colour[k]);
   					hLeftvRightUpstream[i][j][k]->Fit(alpha_fit[i][j][k],"QR+");
   					alpha_fit[i][j][k]->Draw("same");

					alpha_fit[i][j][k]->GetParameters(par);

					sum_par += par[1];

			}

			txt_outfile << i << "\t" << j << "\t" << abs(sum_par/3.) << endl;

		}

	}
*/

/*

	dirDown->cd();

	for( int i=0; i<6; i++ )	{ // loop over detector

		cout << "Detector " << i << endl;

		for ( int j=0; j<4; j++ )	{ // loop over front strips
		
				sprintf(hname,"hLeftvRightDownstream[%d][%d]", i, j);
				sprintf(draw,"SX3StripLeftEnergy:SX3StripRightEnergy>>hLeftvRightDownstream[%d][%d]", i, j);
				sprintf(cut,"SX3Upstream==0 && SX3Det==%d && SX3Strip==%d", i, j);
				TCut Cut = cut;

				hLeftvRightDownstream[i][j] = new TH2D(hname, hname, 4096, 0, 4096, 4096, 0, 4096);
	
				data3->Draw(draw,cut,"goff");

		}

	}

	for( int i=6; i<12; i++ )	{

		cout << "Detector " << i << endl;

		for ( int j=0; j<4; j++ )	{
		
			sprintf(hname,"hLeftvRightDownstream[%d][%d]", i, j);
			sprintf(draw,"SX3StripLeftEnergy:SX3StripRightEnergy>>hLeftvRightDownstream[%d][%d]", i, j);
			sprintf(cut,"SX3Upstream==0 && SX3Det==%d && SX3Strip==%d", i, j);
			TCut Cut = cut;

			hLeftvRightDownstream[i][j] = new TH2D(hname, hname, 4096, 0, 4096, 4096, 0, 4096);
	
			data4->Draw(draw,cut,"goff");


		}

	}
*/
	outfile->cd();
	outfile->Write();
	outfile->Close();	

}
