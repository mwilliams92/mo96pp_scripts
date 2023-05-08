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
vector < Double_t > findpeaks(TH1* h, Double_t sigma, Double_t thresh){

 	vector< Double_t > peaks;
 	int npeaks = 100;

 	h->GetXaxis()->SetRange(100,2000);
 
 	TSpectrum * s = new TSpectrum(npeaks);
 	Int_t nfound = s->Search(h, sigma, "", thresh);

    Double_t * xpeaks = s->GetPositionX();  
    for (int i=0; i<=nfound; i++){
       	if (xpeaks[i]>100 and xpeaks[i]<2000){
          	peaks.push_back(xpeaks[i]);  
      	}
    }
	// sorts vector into ascending order
	sort(peaks.begin(), peaks.end()); 
 	return peaks;
}

// Feed this a vector of peak locations and it will spit out slopes and offset at a pointer to an array
double * findcalib(vector < Double_t > peaks)	{
	// Eu152 lines:
	vector < Double_t > energies = {121.782, 344.279, 778.905, 964.057, 1408.018};
	
	TGraph *gr = new TGraph(5, &peaks[0], &energies[0]);
	gr->SetMarkerStyle(20);
	gr->Draw("goff");

	TF1 *fit = new TF1("fit","pol1");
	//fit->FixParameter(0,0);
	gr->Fit(fit,"Q");
	static double par[2];
	fit->GetParameters(par);

	return par;
}

// Feed this a layer to tell it which root file to use and the channels corresponding to the first ring and sector
void calibrate(TString inputFile, TString outputCalFile, TString outputRootFile){

	TFile *file=TFile::Open(inputFile);

    ofstream txtoutfile;
	txtoutfile.open (outputCalFile);
			
    TH2D *hEu152 = (TH2D*)file->Get("hEgam1_v_Crys[0]");

	const int nCrys = 50;
    int chan;
    
    TFile *outfile = new TFile(outputRootFile,"RECREATE");
    TH1D *hEu152_Crys[nCrys];

    for (int i=0; i<nCrys; i++) {

        hEu152_Crys[i] = (TH1D*) hEu152->ProjectionX(hname_eu152, i, i+1);

        hEu152_Crys[i]->Draw();
		vector < Double_t > peaks = findpeaks(hEu152_Crys[i],2.,0.3);	// Uses findpeaks function to return vector of peak loactions in ascending order.

        chan = i*40 + 9;
/*
        if (peaks.size() != 5) {
			printf("Channel %d (Ring %d) not calibrated! Found %lu peaks instead of 5.\n", i, rID, peaks.size());
			txtoutfile << chan << "\t" << -1 << "\t" << 0 << endl;	// Turns off channel if the number of peaks found is not 5.
		} else {
			double *par = findcalib(/peaks);	// Uses findcalib function to return a pointer to an array containing the calibration offset and slope for each channel
			txtoutfile << chan << "\t" << std::fixed << setprecision(5) << *(par + 0) << "\t" << *(par + 1) << endl;
		}
*/
        hEu152_Crys[i]->Write();

    }

    outfile->Write();
    outfile->Close();

    txtoutfile.close();
}

void calibrateGRETINA()	{

	gROOT->SetBatch(kTRUE); // Set Batch mode so you don't get a shit load of canvases popping up

	// Calibrate
	calibrate("../CheckGretina.root","Eu152_Calib.dat","Eu152.root");

}
