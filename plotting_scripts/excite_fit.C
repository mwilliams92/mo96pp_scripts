#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"

void findpeaks(TH1* hist, double thresh, double sigma, const int npeaks){

	TSpectrum *s = new TSpectrum();
   	Int_t nfound = s->Search(hist,thresh,"",sigma);
	printf("Found %d peaks\n",nfound);

	if (nfound == npeaks) {

		Double_t *xpeaks;
		Double_t xp[npeaks];
   		xpeaks = s->GetPositionX();

		for( int p=0; p<npeaks; p++)	{

			xp[p] = xpeaks[p];

		}

		size_t len = sizeof(xp) / sizeof(xp[0]);
		sort(xp, xp + len);

		return xp;

	} else { 

		printf("Found %d peaks != %d, skipping\n",nfound,npeaks);
		return 0; 
	
	}

}


void calibrateQQQ5(TString layer, Int_t detector){

			TString fname;
			Int_t ringStart, sectorStart;

			if (layer == "E1") {	
				fname = "../protons/cal228th_QQQ5_E1_a.root";	// file name
				ringStart = 578; // first ring channel in ldf file
				sectorStart = 578; // NEED TO CHANGE // first sector channel in ldf file
			}
			else if (layer == "E2")	{	
				fname =	"../protons/cal228th_QQQ5_E1_a.root";	
				ringStart = 578; // NEED TO CHANGE
				sectorStart = 578; // NEED TO CHANGE
			}
			else if (layer == "dE")	{	
				fname =	"../protons/cal228th_QQQ5_E1_a.root";	
				ringStart = 578; // NEED TO CHANGE
				sectorStart = 578; // NEED TO CHANGE
			}
			else { printf("Invalid detector, choose either dE, E1, or E2"); return 0;	}

			TFile *file=TFile::Open(fname);

			TCanvas *cRings = new TCanvas("Rings"); // canvas for rings
			//cRings->Divide(6,6);
			//TCanvas *cSectors = new TCanvas("Sectors"); // canvas for sectors
			//cSectors->Divide(2,2);
			
			const int nRings = 1;
			const int nSectors = 4;

			int chan;

			for(int i=0; i<nRings; i++) {

				//cRings->cd(i+1);

				nchan = ringStart + i;
				sprintf(hname,"d%i", nchan);
				TH1F *hChan = (TH1F*) file->GetObject(hname);
				hChan->Draw();
				findpeaks(hChan, 0.1, 0.05, 5);

			}

}
