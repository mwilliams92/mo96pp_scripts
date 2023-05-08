{
	gROOT->SetBatch(kTRUE);

	TFile *infile = TFile::Open("../../OutputFolder/Run0080.root");
	TTree *data = (TTree*)infile->Get("data");

	TFile *infile1 = TFile::Open("../../OutputFolder/cal228th_315deg_US_BR.root");
	TTree *data1 = (TTree*)infile1->Get("data");

	TFile *infile2 = TFile::Open("../../OutputFolder/cal228th_45deg_US_BL.root");
	TTree *data2 = (TTree*)infile2->Get("data");

	TFile *infile3 = TFile::Open("../../OutputFolder/cal228th_270deg_BR_noBB10.root");
	TTree *data3 = (TTree*)infile3->Get("data");

	TFile *infile4 = TFile::Open("../../OutputFolder/cal228th_90deg_BL_noBB10.root");
	TTree *data4 = (TTree*)infile4->Get("data");

	TFile *infile5 = TFile::Open("../../OutputFolder/cal228th_180deg_DS_noBB10.root"); // good for detector 4
	TTree *data5 = (TTree*)infile5->Get("data");

	TFile *infile6 = TFile::Open("../../OutputFolder/cal228th_135deg_DS_BL.root");
	TTree *data6 = (TTree*)infile6->Get("data");

	char cut[4096];
	char draw[4096];
	char hname[4096];

	TFile *outfile = TFile::Open("FrontPosition.root","RECREATE");
	outfile->cd();

	TH1D *hFrontPositionUpstream[12][4];
	TH1D *hFrontPositionDownstream[12][4];

	TH1 *hFrontPositionUpstreamCumulative[12][4];
	TH1 *hFrontPositionDownstreamCumulative[12][4];

	TH1 *hFrontPositionUpstreamInvCumulative[12][4];
	TH1 *hFrontPositionDownstreamInvCumulative[12][4];

	TDirectory *dirUp = outfile->mkdir("Upstream");
	TDirectory *dirDown = outfile->mkdir("Downstream");

	dirUp->cd();

	bool foundRisingEdge;
	double RisingEdge, FallingEdge, Length;

	ofstream outfile_txt;
	outfile_txt.open ("SX3_PositionCal.txt");

	TLine *l1[12][4];
	TLine *l2[12][4];

	TLine *l3[12][4];
	TLine *l4[12][4];

	for( int i=0; i<12; i++ )	{

		if(i==11) continue;

		cout << "Detector " << i << endl;

		for ( int j=0; j<4; j++ )	{

			sprintf(hname,"hFrontPositionUpstream[%d][%d]", i, j);
			sprintf(draw,"SX3StripPosition>>hFrontPositionUpstream[%d][%d]", i, j);
			sprintf(cut,"SX3Upstream==0 && SX3Det==%d && SX3Strip==%d && SX3SectorEnergy>5000", i, j);
			
			TCut Cut = cut;

			hFrontPositionUpstream[i][j] = new TH1D(hname, hname, 300, -1, 2);

			data->Draw(draw,cut);

			hFrontPositionUpstreamCumulative[i][j] = hFrontPositionUpstream[i][j]->GetCumulative();

           	foundRisingEdge = false; 
			
			for (int k=1; k <= hFrontPositionUpstreamCumulative[i][j]->GetNbinsX(); k++)	{

				if ( foundRisingEdge <1 && hFrontPositionUpstreamCumulative[i][j]->GetBinContent(k) > 0.001*hFrontPositionUpstream[i][j]->GetEntries() ) {

					RisingEdge = k*0.01 - 1.;

					l1[i][j] = new TLine(RisingEdge,0.,RisingEdge,80.);
					l1[i][j]->SetLineColor(2);
					hFrontPositionUpstream[i][j]->GetListOfFunctions()->Add(l1[i][j]);
			
					foundRisingEdge = 1;

					outfile_txt << i << "\t" << j << "\t" << RisingEdge << "\t";

				}

				if ( hFrontPositionUpstreamCumulative[i][j]->GetBinContent(k) > 0.995*hFrontPositionUpstream[i][j]->GetEntries() ) {

					FallingEdge = k*0.01 - 1.;

					l2[i][j] = new TLine(FallingEdge,0.,FallingEdge,80.);
					l2[i][j]->SetLineColor(2);
					hFrontPositionUpstream[i][j]->GetListOfFunctions()->Add(l2[i][j]);

					Length = FallingEdge - RisingEdge;			

					outfile_txt << FallingEdge << "\t" << Length << endl;

					break;

				}

			}

		}

	}
/*
		for( int i=6; i<12; i++ )	{

		cout << "Detector " << i << endl;

		for ( int j=0; j<4; j++ )	{

			sprintf(hname,"hFrontPositionUpstream[%d][%d]", i, j);
			sprintf(draw,"SX3StripPosition>>hFrontPositionUpstream[%d][%d]", i, j);
			sprintf(cut,"SX3Upstream==1 && SX3Det==%d && SX3Strip==%d && SX3SectorEnergy>6600 && SX3SectorEnergy<7000", i, j);
			
			TCut Cut = cut;

			hFrontPositionUpstream[i][j] = new TH1D(hname, hname, 200, -1, 1);

			data2->Draw(draw,cut,"same");

			hFrontPositionUpstreamCumulative[i][j] = hFrontPositionUpstream[i][j]->GetCumulative();

           	foundRisingEdge = false; 

			for (int k=1; k <= hFrontPositionUpstreamCumulative[i][j]->GetNbinsX(); k++)	{

				if ( foundRisingEdge <1 && hFrontPositionUpstreamCumulative[i][j]->GetBinContent(k) > 0.005*hFrontPositionUpstream[i][j]->GetEntries() ) {

					RisingEdge = k*0.01 - 1.;

					l1[i][j] = new TLine(RisingEdge,0.,RisingEdge,80.);
					l1[i][j]->SetLineColor(2);
					hFrontPositionUpstream[i][j]->GetListOfFunctions()->Add(l1[i][j]);
			
					foundRisingEdge = 1;

					outfile_txt << i << "\t" << j << "\t" << RisingEdge << "\t";

				}

				if ( hFrontPositionUpstreamCumulative[i][j]->GetBinContent(k) > 0.995*hFrontPositionUpstream[i][j]->GetEntries() ) {

					FallingEdge = k*0.01 - 1.;

					l2[i][j] = new TLine(FallingEdge,0.,FallingEdge,80.);
					l2[i][j]->SetLineColor(2);
					hFrontPositionUpstream[i][j]->GetListOfFunctions()->Add(l2[i][j]);

					Length = FallingEdge - RisingEdge; 			

					outfile_txt << FallingEdge << "\t" << Length << endl;

					break;

				}

			}

		}

	}
*/
/*
	dirDown->cd();

	for( int i=0; i<6; i++ )	{

	cout << "Detector " << i << endl;

		for ( int j=0; j<4; j++ )	{

			sprintf(hname,"hFrontPositionDownstream[%d][%d]", i, j);
			sprintf(draw,"SX3StripPosition>>hFrontPositionDownstream[%d][%d]", i, j);
			sprintf(cut,"SX3Upstream==0 && SX3Det==%d && SX3Strip==%d && SX3SectorEnergy>6600 && SX3SectorEnergy<7000", i, j);
			
			TCut Cut = cut;

			hFrontPositionDownstream[i][j] = new TH1D(hname, hname, 200, -1, 1);

			data3->Draw(draw,cut,"same");

			hFrontPositionDownstreamCumulative[i][j] = hFrontPositionDownstream[i][j]->GetCumulative();

           	foundRisingEdge = false; 

			for (int k=1; k <= hFrontPositionDownstreamCumulative[i][j]->GetNbinsX(); k++)	{

				if ( foundRisingEdge <1 && hFrontPositionDownstreamCumulative[i][j]->GetBinContent(k) > 0.005*hFrontPositionDownstream[i][j]->GetEntries() ) {

					RisingEdge = k*0.01 - 1.;

					l3[i][j] = new TLine(RisingEdge,0.,RisingEdge,80.);
					l3[i][j]->SetLineColor(2);
					hFrontPositionDownstream[i][j]->GetListOfFunctions()->Add(l3[i][j]);
			
					foundRisingEdge = 1;

					outfile_txt << i << "\t" << j << "\t" << RisingEdge << "\t";

				}

				if ( hFrontPositionDownstreamCumulative[i][j]->GetBinContent(k) > 0.995*hFrontPositionDownstream[i][j]->GetEntries() ) {

					FallingEdge = k*0.01 - 1.;

					l4[i][j] = new TLine(FallingEdge,0.,FallingEdge,80.);
					l4[i][j]->SetLineColor(2);
					hFrontPositionDownstream[i][j]->GetListOfFunctions()->Add(l4[i][j]);

					Length = FallingEdge - RisingEdge;		

					outfile_txt << FallingEdge << "\t" << Length << endl;

					break;

				}

			}

		}

	}

	for( int i=6; i<12; i++ )	{

	cout << "Detector " << i << endl;

		for ( int j=0; j<4; j++ )	{

			sprintf(hname,"hFrontPositionDownstream[%d][%d]", i, j);
			sprintf(draw,"SX3StripPosition>>hFrontPositionDownstream[%d][%d]", i, j);
			sprintf(cut,"SX3Upstream==0 && SX3Det==%d && SX3Strip==%d && SX3SectorEnergy>6600 && SX3SectorEnergy<7000", i, j);
			
			TCut Cut = cut;

			hFrontPositionDownstream[i][j] = new TH1D(hname, hname, 200, -1, 1);

			data4->Draw(draw,cut,"same");

			hFrontPositionDownstreamCumulative[i][j] = hFrontPositionDownstream[i][j]->GetCumulative();

           	foundRisingEdge = false; 

			for (int k=1; k <= hFrontPositionDownstreamCumulative[i][j]->GetNbinsX(); k++)	{

				if ( foundRisingEdge <1 && hFrontPositionDownstreamCumulative[i][j]->GetBinContent(k) > 0.005*hFrontPositionDownstream[i][j]->GetEntries() ) {

					RisingEdge = k*0.01 - 1.;

					l3[i][j] = new TLine(RisingEdge,0.,RisingEdge,80.);
					l3[i][j]->SetLineColor(2);
					hFrontPositionDownstream[i][j]->GetListOfFunctions()->Add(l3[i][j]);
			
					foundRisingEdge = 1;

					outfile_txt << i << "\t" << j << "\t" << RisingEdge << "\t";

				}

				if ( hFrontPositionDownstreamCumulative[i][j]->GetBinContent(k) > 0.995*hFrontPositionDownstream[i][j]->GetEntries() ) {

					FallingEdge = k*0.01 - 1.;

					l4[i][j] = new TLine(FallingEdge,0.,FallingEdge,80.);
					l4[i][j]->SetLineColor(2);
					hFrontPositionDownstream[i][j]->GetListOfFunctions()->Add(l4[i][j]);

					Length = FallingEdge - RisingEdge;			

					outfile_txt << FallingEdge << "\t" << Length << endl;

					break;

				}

			}

		}

	}
*/
	outfile->cd();
	outfile->Write();
	outfile->Close();	

}
