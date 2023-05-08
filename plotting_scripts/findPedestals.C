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

// Feed this a histogram and it will spit out the peak locations as vector in ascending order. Also set sigma and thresh to desired value
vector < Double_t > findped(TH1* h, Double_t sigma, Double_t thresh){

 	vector< Double_t > peaks;
 	int npeaks = 10;

 	h->GetXaxis()->SetRange(0,500);
 
 	TSpectrum * s = new TSpectrum(npeaks);
 	Int_t nfound = s->Search(h, sigma, "", thresh);

    Double_t * xpeaks = s->GetPositionX();  
    for (int i=0; i<=nfound; i++){
       	//if (xpeaks[i]>400 and xpeaks[i]<10000){
          	peaks.push_back(xpeaks[i]);  
      	//}
    }
	// sorts vector into ascending order
	sort(peaks.begin(), peaks.end()); 
 	return peaks;
}

void pedestalSX3(TFile *root_outfile, Int_t FirstChan, TString outputFile)	{

	TFile *file=TFile::Open("~/experiments/mo96pp/protons/cal228th/pedestal03/pedestal03.root");
	
	const int nDets = 12;		
	const int nStrips = 4;

	int chanID = FirstChan;
	char dname[32], hname[32];

	TDirectory *dir[nDets];
	
	TH1F *hStrip[nDets][nStrips];

	ofstream outfile;
	outfile.open (outputFile);
	
	for(int j=0; j<nDets; j++)	{

		sprintf(dname,"dirSX3_%d",j);
		dir[j] = root_outfile->mkdir(dname);
		dir[j]->cd();	

		for(int i=0; i<nStrips; i++) {
		
			sprintf(hname,"d%d", chanID);	// Name of hist inside input root file
			
			hStrip[j][i] = (TH1F*) file->Get(hname);	// Grabs hist from input root file
			hStrip[j][i]->Draw("goff");

			vector < Double_t > pedestal = findped(hStrip[j][i],5.,0.2);	// Uses findpeaks function to return vector of peak loactions in ascending order.

			outfile << j << "\t" << i << "\t" << std::fixed << setprecision(5) << pedestal.at(1) << endl;

			hStrip[j][i]->GetXaxis()->SetRangeUser(0,200);

			chanID++;

			dir[j]->Append(hStrip[j][i]);

		}

		file->cd();

	}

}


void findPedestals()	{

	gROOT->SetBatch(kTRUE); // Set Batch mode so you don't get a shit load of canvases popping up

	// Create ROOT file
	TFile *root_outfile = new TFile("pedestalORRUBA.root", "RECREATE");

	// SX3 Front Downstream Pedestals
	//pedestalSX3(root_outfile,289,"SX3d_pedestals.dat");
	// SX3 Front Upstream Pedestals
	//pedestalSX3(root_outfile,193,"SX3u_pedestals.dat");
	// SX3 Back Downstream Pedestals
	//pedestalSX3(root_outfile,385,"SX3dBack_pedestals.dat");
	// SX3 Back Upstream Pedestals
	pedestalSX3(root_outfile,145,"SX3uBack_pedestals.dat");

	// Write root file
	root_outfile->cd();
	root_outfile->Write();
	root_outfile->Close();	

}
