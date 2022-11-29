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
 	int npeaks = 6;

 	h->GetXaxis()->SetRange(300,4000);
 
 	TSpectrum * s = new TSpectrum(npeaks);
 	Int_t nfound = s->Search(h, sigma, "", thresh);

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
void calibrateQQQ5(TString inputFile, Int_t FirstRing, Int_t FirstSector, TString outputFile){

	TFile *file=TFile::Open(inputFile);

	//if (vis)	{
		TCanvas *cQQQ5 = new TCanvas("QQQ5"); // canvas for rings
		cQQQ5->Divide(6,6);
	//}
			
	const int nRings = 32;
	const int nSectors = 4;

	int cID, rID, sID;
	char hname[8];
	TH1F *hRing[nRings], *hSector[nSectors];

	ofstream outfile;
	outfile.open (outputFile);

	for(int i=FirstRing; i<nRings+FirstRing; i++) {

		rID = i - FirstRing;		// Ring number within QQQ5 0->31
		
		sprintf(hname,"d%i", i);	// Name of hist inside input root file

		hRing[rID] = (TH1F*) file->Get(hname);	// Grabs hist from input root file

		cID = rID + 1;	// Subcanvas ID
		cQQQ5->cd(cID);

		hRing[rID]->Draw();

		vector < Double_t > peaks = findpeaks(hRing[rID],5.,0.3);	// Uses findpeaks function to return vector of peak loactions in ascending order.

		if (peaks.size() != 5) {
			printf("Channel %d (Ring %d) not calibrated!\n", i, rID);
			outfile << i << "\t" << -1 << "\t" << 0 << endl;	// Turns off channel if the number of peaks found is not 5.
		} else {
			double *par = findcalib(peaks);	// Uses findcalib function to return a pointer to an array containing the calibration offset and slope for each channel
			outfile << i << "\t" << std::fixed << setprecision(5) << *(par + 0) << "\t" << *(par + 1) << endl;
		}

		hRing[rID]->GetXaxis()->SetRangeUser(300,4000);

	}

	for(int i=FirstSector; i<nSectors+FirstSector; i++) {

		sID = i - FirstSector;		// Sector number (within QQQ5) 0->3
		sprintf(hname,"d%i", i);		// Name of hist inside input root file

		hSector[sID] = (TH1F*) file->Get(hname);	// Grabs hist from input root file
		
		cID = sID + nRings + 1;	// Subcanvas number (starts after rings)
		cQQQ5->cd(cID);

		hSector[sID]->Draw();

		vector < Double_t > peaks = findpeaks(hSector[sID], 5., 0.3); // Uses findpeaks function to return vector of peak loactions in ascending order.

		if (peaks.size() != 5) {
			printf("Channel %d (Sector %d) not calibrated!\n", i, sID);
			outfile << i << "\t" << -1 << "\t" << 0 << endl;	// Turns off channel if the number of peaks found is not 5.
		} else {
			double *par = findcalib(peaks);	// Uses findcalib function to return a pointer to an array containing the calibration offset and slope for each channel
			outfile << i << "\t" << std::fixed << setprecision(5) << *(par + 0) << "\t" << *(par + 1) << endl;
		}

		hSector[sID]->GetXaxis()->SetRangeUser(300,4000);	

		}

		outfile.close();
}

// Modified version of calibrateQQQ5 to account for rings not being instrumented in the downstream QQQ5 A detector in the E1 layer.
void calibrateQQQ5mod(TString inputFile, Int_t FirstSector, TString outputFile){

	TFile *file=TFile::Open(inputFile);

	TCanvas *cQQQ5 = new TCanvas("QQQ5"); // canvas for rings
	cQQQ5->Divide(2,2);
			
	const int nSectors = 4;

	int cID, sID;
	char hname[8];
	TH1F *hSector[nSectors];

	ofstream outfile;
	outfile.open (outputFile);

	for(int i=FirstSector; i<nSectors+FirstSector; i++) {

		sID = i - FirstSector;		// Sector number (within QQQ5) 0->3
		sprintf(hname,"d%i", i);		// Name of hist inside input root file

		hSector[sID] = (TH1F*) file->Get(hname);	// Grabs hist from input root file
		
		cID = sID + 1;	// Subcanvas number (starts after rings)
		cQQQ5->cd(cID);

		hSector[sID]->Draw();

		vector < Double_t > peaks = findpeaks(hSector[sID], 5., 0.3); // Uses findpeaks function to return vector of peak loactions in ascending order.

		if (peaks.size() != 5) {
			printf("Channel %d (Sector %d) not calibrated!\n", i, sID);
			outfile << i << "\t" << -1 << "\t" << 0 << endl;	// Turns off channel if the number of peaks found is not 5.
		} else {
			double *par = findcalib(peaks);	// Uses findcalib function to return a pointer to an array containing the calibration offset and slope for each channel
			outfile << i << "\t" << std::fixed << setprecision(5) << *(par + 0) << "\t" << *(par + 1) << endl;
		}

		hSector[sID]->GetXaxis()->SetRangeUser(300,4000);	

		}

		outfile.close();
}

void calibrateIndv(TString inputFile, Int_t chan, Double_t sigma, Double_t thresh)	{

	TFile *file=TFile::Open(inputFile);

	TCanvas *cChan = new TCanvas("Chan"); // canvas for rings

	char hname[8];
	sprintf(hname,"d%i", chan);

	TH1F *hChan = (TH1F*) file->Get(hname);
	hChan->Draw();

	vector < Double_t > peaks = findpeaks(hChan, sigma, thresh);
	cout << "Found " << peaks.size() << " peaks" << endl;

	double *par = findcalib(peaks);	// Uses findcalib function to return a pointer to an array containing the calibration offset and slope for each channel
	cout << "Offset = " << *(par + 0) << " Slope = " << *(par + 1) << endl;

}

void calibrateORRUBA()	{

	// At the moment visualisation does not work correctly, since the next instance of calibrateQQQ5 overwrites the previous canvas.
	// Ultimately it would be good to create an output root file in this function, which is then passed as an input for the calibrateQQQ5
	// then write the canvases to that output rootfile. One could even write a canvas for each channel and group them into directories
	// so each calibrateQQQ5 will create a new directory in the output root file to dump the canvases. 

	//calibrateQQQ5("../protons/cal228th_180deg_DS.root",497,561,"QQQ5dAdEcalib.dat");	// Downstream dE A
	calibrateQQQ5("../protons/cal228th_180deg_DS.root",529,569,"QQQ5dBdEcalib.dat");	// Downstream dE B
	//calibrateQQQ5("../protons/cal228th_QQQ5_E1_a.root",577,673,"QQQ5dAE1calib.dat");	// Downstream E1 A
	//calibrateQQQ5("../protons/cal228th_QQQ5_E1_a.root",609,681,"QQQ5dBE1calib.dat");	// Downstream E1 B
	//calibrateQQQ5mod("../protons/cal228th_QQQ5_E2.root",677,"QQQ5dAE2calib.dat");		// Downstream E2 A // Modified to only calibrate sectors (rings not instrumented)
	//calibrateQQQ5("../protons/cal228th_QQQ5_E2.root",641,685,"QQQ5dBE2calib.dat");	// Downstream E2 B

}
