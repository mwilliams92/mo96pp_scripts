{

			TFile *file=TFile::Open("../HistFiles/Sum_Nov2020_Target2.root");
			TDirectory *dir = (TDirectory*)file->Get("S3_EMMA");

			double high_theta;
			double low_theta;

			const int n=14;
			double theta[n] = {159,157,155,153,151,149,147,145,143,141,139,137,135,133};
			double dummy[n] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};

			double counts[n] = {0};
			
			
			double error[n] = {0};

			for(int i=0; i<1; i++) {

				TCanvas *c1 = new TCanvas("c1");
				c1->SetTicky();
				c1->SetTickx();
				TH2F *h_exc_theta = (TH2F*) dir->Get("h_theta_v_excite_gated");

				high_theta = 160. - 2*i;
				low_theta = 158. - 2*i;

				TH1F *hExc = (TH1F*)h_exc_theta->ProjectionY("h_theta_v_excite_gated",low_theta,high_theta);
				hExc->SetLineColor(1);
				hExc->Rebin(10);
				hExc->Draw();

				hExc->SetTitle("Excite 158 to 160 deg");
				hExc->GetYaxis()->SetTitle("Counts [80 keV/bin]");
				hExc->GetXaxis()->SetTitle("Excitation Energy [MeV]");
				hExc->GetXaxis()->SetRangeUser(1.0,8.5);

				TSpectrum *s = new TSpectrum();
   				Int_t nfound = s->Search(hExc,0.1,"",0.05);
				printf("Found %d candidate peaks to fit\n",nfound);

				//if (nfound != 8) { continue; }

				int npeaks = 0;
				double par[24] = {0};

				Double_t *xpeaks;
   				xpeaks = s->GetPositionX();

   				for (int p=0;p<nfound;p++) {
      					Double_t xp = xpeaks[p];
      					Int_t bin = hExc->GetXaxis()->FindBin(xp);
      					Double_t yp = hExc->GetBinContent(bin);
      					par[3*npeaks] = yp; // "height"
      					par[3*npeaks+1] = xp; // "mean"
      					par[3*npeaks+2] = 0.15; // "sigma"
      					npeaks++;
   				}

				//TF1 *total = new TF1("total","gaus(0)+gaus(3)+pol0(6)",6.7,8.4);
				//TF1 *total = new TF1("total","gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)+gaus(15)+gaus(18)+gaus(21)+gaus(24)+pol2(27)",1.0,10.0);
				TF1 *total = new TF1("total","gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)+gaus(15)+gaus(18)+gaus(21)",1.0,8.5);
   				// Use the parameters on the sum.
   				total->SetParameters(par);
				total->SetNpx(1000);
   				auto fitResult = hExc->Fit(total,"RS");

				total->GetParameters(par);

/*	
				TF1 *bg = new TF1("bg","pol2",1.0,11.0);
				bg->FixParameter(0,par[27]);
				bg->FixParameter(1,par[28]);
				bg->FixParameter(2,par[29]);
				bg->SetLineColor(4);
				bg->Draw("same");			
*/
				TF1 *g1 = new TF1("g1","gaus",1.0,11.0);
				g1->FixParameter(0,par[0]);
				g1->FixParameter(1,par[1]);
				g1->FixParameter(2,par[2]);
				g1->SetLineColor(1);
				g1->SetLineStyle(9);
				g1->Draw("same");
				
				counts[i] = 10.* g1->Integral(0,100);
				error[i] = sqrt(counts[i]);

				cout << "Integral = " << counts[i] << endl;

				//delete c1;
				//delete h_exc_theta;
				//delete hExc;
				//delete 2178.87s;
				//delete total;
				//delete g1;
			}

			counts[0] = 1698.14;
			error[0] = sqrt(1698.14); 

			counts[11] = 2178.87;

			error[11] = sqrt(2178.87);

			counts[13] = 1419.1;
			error[13] = sqrt(1419.1); 

			TCanvas *c2 = new TCanvas();
			TGraphErrors *gr1 = new TGraphErrors(n,theta,counts,dummy,error);
			gr1->SetMarkerStyle(20);
			gr1->Draw("AP");

}
