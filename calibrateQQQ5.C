#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"

// Feed this a histogram and it will spit out the peak locations as vector in ascending order
vector < Double_t > findpeaks(TH1* h){

 	vector< Double_t > peaks;
 	int npeaks = 6;

 	h->GetXaxis()->SetRange(300,1500);
 
 	TSpectrum * s = new TSpectrum(npeaks);
 	Int_t nfound = s->Search(h, 6., "", 0.3);

    Double_t * xpeaks = s->GetPositionX();  
    for (int i=0; i<=nfound; i++){
       	if (xpeaks[i]>400 and xpeaks[i]<10000){
          	peaks.push_back(xpeaks[i]);  
      	}
    }

	sort(peaks.begin(), peaks.end()); 
 	return peaks;
}

// Feed this a vector of peak locations and it will spit out slopes and offset at a pointer to an array
double * findcalib(vector < Double_t > peaks)	{
	// Th228 alpha lines:
	vector < Double_t > energies = {5423., 5685., 6288., 6778, 8784.};
	
	TGraph *gr = new TGraph(5, &peaks[0], &energies[0]);
	gr->SetMarkerStyle(20);
	gr->Draw("goff");

	TF1 *fit = new TF1("fit","pol1");
	gr->Fit(fit,"Q");
	static double par[2];
	fit->GetParameters(par);

	return par;
}

// Feed this a layer to tell it which root file to use and the channels corresponding to the first ring and sector
void calibrateQQQ5(TString layer, Int_t FirstRing){

			TString fname;
			Int_t ringStart, sectorStart;

			if (layer == "E1") {	
				fname = "../protons/cal228th_QQQ5_E1_a.root";	// file name
			}
			else if (layer == "E2")	{	
				fname =	"../protons/cal228th_QQQ5_E1_a.root";	
			}
			else if (layer == "dE")	{	
				fname =	"../protons/cal228th_QQQ5_E1_a.root";	
			}
			else { printf("Invalid detector, choose either dE, E1, or E2"); return 0;	}

			TFile *file=TFile::Open(fname);

			TCanvas *cRings = new TCanvas("Rings"); // canvas for rings
			cRings->Divide(6,6);
			//TCanvas *cSectors = new TCanvas("Sectors"); // canvas for sectors
			//cSectors->Divide(2,2);
			
			const int nRings = 32;
			const int nSectors = 4;

			int nchan;
			char hname[8];
			TH1F *hRing[32];

			ofstream outfile;
			outfile.open ("ringcalib.dat");

			for(int i=0; i<nRings; i++) {

				cRings->cd(i+1);
				nchan = FirstRing + i;
				sprintf(hname,"d%i", nchan);

				hRing[i] = (TH1F*) file->Get(hname);	
				hRing[i]->Draw();

				vector < Double_t > peaks = findpeaks(hRing[i]);

				if (peaks.size() != 5) {
					outfile << i << "\t" << -1 << "\t" << 0 << endl;
				} else {
					double *par = findcalib(peaks);
					outfile << i << "\t" << std::fixed << setprecision(5) << *(par + 0) << "\t" << *(par + 1) << endl;
				}

				hRing[i]->GetXaxis()->SetRangeUser(300,1000);	

			}
			outfile.close();
}
