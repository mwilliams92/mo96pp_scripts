#include <iostream>
#include <algorithm>

using namespace std;

void calqqq5(){

	TFile *file=TFile::Open("../protons/cal228th_QQQ5_E1_a.root");
	
	Int_t start_chan = 609;
	Int_t end_chan = 624;

	char histname[10];
	
	TCanvas *c1 = new TCanvas("c1");
	c1->Divide(4,4);
	
	TH1F *hist[16];
	
	Int_t j = 0;
	Int_t nfound;
	
	// Guess required threshhold and sigma for TSpectrum Search
	Double_t thresh = 0.08;
	Double_t sigma = 2.0;
	
	// Actual Th228 alpha energies:
    	Double_t alpha_energies[7] = {5340.36, 5423.15, 5685.37, 6050.78, 6288.08, 6778.30, 8784.86}; // the 6050.78 one should not be used
    	
    	ofstream outfile;
  	outfile.open ("qqq5cal.txt");
	
	for (int i=start_chan; i<=end_chan; i++)	{	
		sprintf(histname,"d%d",i);
		hist[j] = (TH1F*) file->Get(histname);
		c1->cd(j+1);
		hist[j]->Draw();
		hist[j]->GetXaxis()->SetRangeUser(400,1000);
		
		Double_t centroid[7] = {0};
		TSpectrum *s = new TSpectrum();
   		Int_t nfound = s->Search(hist[j],sigma,"",thresh);
   		j++;
		if (nfound != 7) { printf("Channel %d found %d peaks, skipping\n",i,nfound);	continue; }

		Double_t *xpeaks;
   		xpeaks = s->GetPositionX();
   		
   		// allocate elements of pointer into an array:
   		for (int p=0;p<nfound;p++) {
      			centroid[p] = xpeaks[p];
   		}
   	
   	
   		// Sort array in ascending order:
   		size_t len = sizeof(centroid) / sizeof(centroid[0]);
    		sort(centroid, centroid + len);
    		
    		TGraph *gr = new TGraph(7,centroid,alpha_energies);
    		gr->Draw("goff");
    		
    		TF1 *fit = new TF1("fit","pol1",0,1000);
    		gr->Fit(fit);
    		
    		Double_t offset = fit->GetParameter(0);
    		Double_t slope = fit->GetParameter(1);
    		
    		outfile << i << "\t" << offset << "\t" << slope << endl; 
    		
    		delete gr;
    		delete fit;
    	
    		
	}
   	


}
