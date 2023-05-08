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

 	h->GetXaxis()->SetRange(200,4000);
 
 	TSpectrum * s = new TSpectrum(npeaks);
 	Int_t nfound = s->Search(h, sigma, "", thresh);

    Double_t * xpeaks = s->GetPositionX();  
    for (int i=0; i<=nfound; i++){
       	if (xpeaks[i]>100 and xpeaks[i]<10000){
          	peaks.push_back(xpeaks[i]);  
      	}
    }
	// sorts vector into ascending order
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
	//fit->FixParameter(0,0);
	gr->Fit(fit,"Q");
	static double par[2];
	fit->GetParameters(par);

	return par;
}

// Feed this a layer to tell it which root file to use and the channels corresponding to the first ring and sector
void calibrateQQQ5(TString inputFile, Int_t FirstRing, Int_t FirstSector, TString outputFile, TDirectory *dir){

	TFile *file=TFile::Open(inputFile);
	dir->cd();

	//TCanvas *cQQQ5 = new TCanvas("QQQ5"); // canvas for rings
	//cQQQ5->Divide(6,6);
			
	const int nRings = 32;
	const int nSectors = 4;

	int cID, rID, sID;
	char hname[8];
	TH1F *hRing[nRings], *hSector[nSectors];

	ofstream outfile;
	outfile.open (outputFile);

	TList *listRings = new TList();

	TCanvas *cRing[nRings];  

	for(int i=FirstRing; i<nRings+FirstRing; i++) {

		rID = i - FirstRing;		// Ring number within QQQ5 0->31
		
		sprintf(hname,"d%i", i);	// Name of hist inside input root file

		hRing[rID] = (TH1F*) file->Get(hname);	// Grabs hist from input root file

		cID = rID + 1;	// Subcanvas ID
		//cQQQ5->cd(cID);
		cRing[cID-1] = new TCanvas(Form("c%1.1i",i),Form("Ring %1.1i",rID),30,113,800,600);
		hRing[rID]->Rebin(2);
		hRing[rID]->Draw();

		vector < Double_t > peaks = findpeaks(hRing[rID],3.5,0.32);	// Uses findpeaks function to return vector of peak loactions in ascending order.

		listRings->Add(cRing[cID-1]);  

		if (peaks.size() != 5) {
			printf("Channel %d (Ring %d) not calibrated! Found %lu peaks instead of 5.\n", i, rID, peaks.size());
			outfile << i << "\t" << -1 << "\t" << 0 << endl;	// Turns off channel if the number of peaks found is not 5.
		} else {
			double *par = findcalib(peaks);	// Uses findcalib function to return a pointer to an array containing the calibration offset and slope for each channel
			outfile << i << "\t" << std::fixed << setprecision(5) << *(par + 0) << "\t" << *(par + 1) << endl;
		}

		hRing[rID]->GetXaxis()->SetRangeUser(300,4000);

	}

	dir->Append(listRings);

	TList *listSectors = new TList();

	TCanvas *cSector[nSectors];

	for(int i=FirstSector; i<nSectors+FirstSector; i++) {

		sID = i - FirstSector;		// Sector number (within QQQ5) 0->3
		sprintf(hname,"d%i", i);		// Name of hist inside input root file

		hSector[sID] = (TH1F*) file->Get(hname);	// Grabs hist from input root file
		
		cID = sID + nRings + 1;	// Subcanvas number (starts after rings)
		//cQQQ5->cd(cID);

		cSector[cID-1] = new TCanvas(Form("c%1.1i",i),Form("Sector%1.1i",sID),30,113,800,600);

		hSector[sID]->Draw();

		vector < Double_t > peaks = findpeaks(hSector[sID], 5., 0.4); // Uses findpeaks function to return vector of peak loactions in ascending order.

		listSectors->Add(cSector[cID-1]);


		if (peaks.size() != 5) {
			printf("Channel %d (Sector %d) not calibrated! Found %lu peaks instead of 5.\n", i, sID, peaks.size());
			outfile << i << "\t" << -1 << "\t" << 0 << endl;	// Turns off channel if the number of peaks found is not 5.
		} else {
			double *par = findcalib(peaks);	// Uses findcalib function to return a pointer to an array containing the calibration offset and slope for each channel
			outfile << i << "\t" << std::fixed << setprecision(5) << *(par + 0) << "\t" << *(par + 1) << endl;
		}

		hSector[sID]->GetXaxis()->SetRangeUser(300,4000);	

	}

	dir->Append(listSectors);
	outfile.close();
}

// Modified version of calibrateQQQ5 to account for rings not being instrumented in the downstream QQQ5 A detector in the E1 layer.
void calibrateQQQ5mod(TString inputFile, Int_t FirstSector, TString outputFile, TDirectory *dir){

	TFile *file=TFile::Open(inputFile);
	dir->cd();

	//TCanvas *cQQQ5 = new TCanvas("QQQ5"); // canvas for rings
	//cQQQ5->Divide(2,2);
			
	const int nSectors = 4;

	int cID, sID;
	char hname[8];
	TH1F *hSector[nSectors];

	ofstream outfile;
	outfile.open (outputFile);

	TList *listSectors = new TList();

	TCanvas *cSector[nSectors];

	for(int i=FirstSector; i<nSectors+FirstSector; i++) {

		sID = i - FirstSector;		// Sector number (within QQQ5) 0->3
		sprintf(hname,"d%i", i);		// Name of hist inside input root file

		hSector[sID] = (TH1F*) file->Get(hname);	// Grabs hist from input root file
		
		cID = sID + 1;	// Subcanvas number (starts after rings)
		//cQQQ5->cd(cID);

		cSector[cID-1] = new TCanvas(Form("c%1.1i",i),Form("Sector%1.1i",sID),30,113,800,600);

		hSector[sID]->Draw();

		vector < Double_t > peaks = findpeaks(hSector[sID], 5., 0.34); // Uses findpeaks function to return vector of peak loactions in ascending order.

		listSectors->Add(cSector[cID-1]);

		if (peaks.size() != 5) {
			printf("Channel %d (Sector %d) not calibrated!\n", i, sID);
			outfile << i << "\t" << -1 << "\t" << 0 << endl;	// Turns off channel if the number of peaks found is not 5.
		} else {
			double *par = findcalib(peaks);	// Uses findcalib function to return a pointer to an array containing the calibration offset and slope for each channel
			outfile << i << "\t" << std::fixed << setprecision(5) << *(par + 0) << "\t" << *(par + 1) << endl;
		}

		hSector[sID]->GetXaxis()->SetRangeUser(300,4000);	

	}
	dir->Append(listSectors);
	outfile.close();
}

void calibrateIndv(TString inputFile, Int_t chan, Double_t sigma, Double_t thresh)	{

	TFile *file=TFile::Open(inputFile);

	TCanvas *cChan = new TCanvas("Chan"); // canvas for rings

	char hname[8];
	sprintf(hname,"d%i", chan);

	TH1F *hChan = (TH1F*) file->Get(hname);
	//hChan->Rebin(2);
	hChan->Draw();

	vector < Double_t > peaks = findpeaks(hChan, sigma, thresh);
	cout << "Found " << peaks.size() << " peaks" << endl;

	double *par = findcalib(peaks);	// Uses findcalib function to return a pointer to an array containing the calibration offset and slope for each channel
	cout << "Offset = " << *(par + 0) << " Slope = " << *(par + 1) << endl;

}

void calibrateBB10(TString inputFile, Int_t FirstStrip, TString outputFile, TDirectory *dir)	{

	TFile *file=TFile::Open(inputFile);
	dir->cd();

	//TCanvas *cBB10 = new TCanvas("BB10"); // canvas for rings
	//cBB10->Divide(8,1);
			
	const int nStrips = 8;

	int cID, sID;
	char hname[8];
	TH1F *hStrip[nStrips];

	ofstream outfile;
	outfile.open (outputFile);

	TList *listStrips = new TList();

	TCanvas *cStrip[nStrips];  

	for(int i=FirstStrip; i<nStrips+FirstStrip; i++) {

		sID = i - FirstStrip;		// Ring number within QQQ5 0->31
		
		sprintf(hname,"d%i", i);	// Name of hist inside input root file

		hStrip[sID] = (TH1F*) file->Get(hname);	// Grabs hist from input root file

		cID = sID + 1;	// Subcanvas ID
		//cBB10->cd(cID);
		cStrip[cID-1] = new TCanvas(Form("c%1.1i",i),Form("Strip %1.1i",sID),30,113,800,600);
		//hRing[rID]->Rebin(2);
		hStrip[sID]->Draw();

		vector < Double_t > peaks = findpeaks(hStrip[sID],5.,0.4);	// Uses findpeaks function to return vector of peak loactions in ascending order.

		listStrips->Add(cStrip[cID-1]);  

		if (peaks.size() != 5) {
			printf("Channel %d (Strip %d) not calibrated! Found %lu peaks instead of 5.\n", i, sID, peaks.size());
			outfile << i << "\t" << -1 << "\t" << 0 << endl;	// Turns off channel if the number of peaks found is not 5.
		} else {
			double *par = findcalib(peaks);	// Uses findcalib function to return a pointer to an array containing the calibration offset and slope for each channel
			outfile << i << "\t" << std::fixed << setprecision(5) << *(par + 0) << "\t" << *(par + 1) << endl;
		}

		hStrip[sID]->GetXaxis()->SetRangeUser(300,4000);

	}

	dir->Append(listStrips);

}

void calibrateSX3(TString inputFile, TString outputFile)	{

	TFile *infile=TFile::Open(inputFile);
	TTree *data = (TTree*)infile->Get("data");

	TFile *outfile=TFile::Open("test.root","RECREATE");
	outfile->cd();

	char cut[4096];
	char draw[4096];
	char hname[4096];
	char fname[4096];
	char gname[4096];

	const int nDets = 8;
	const int nStrips = 4;
	int detID;
			
	TH1D *hBackStrip[nDets][nStrips];
	TGraphErrors *gr[nDets][nStrips];
	TF1 *fit[nDets][nStrips];
	TCanvas *cStrip[nDets][nStrips];

	ofstream txt_outfile;
	txt_outfile.open (outputFile);

	for( int i=7; i<nDets; i++ )	{

		detID = i;

		cout << "Detector " << detID << endl;

		for ( int j=0; j<nStrips; j++ )	{

			sprintf(hname,"hBackStrip[%d][%d]", i, j);
			sprintf(draw,"SX3SectorEnergy>>hBackStrip[%d][%d]", i, j);
			sprintf(cut,"SX3Upstream==1 && SX3Det==%d && SX3Sector==%d && SX3Mul==1", detID, j);
			TCut Cut = cut;

			hBackStrip[i][j] = new TH1D(hname, hname, 4096, 0, 4096);
	
			data->Draw(draw,cut,"goff");
            //hBackStrip[i][j]->Rebin(2);

			vector < Double_t > peaks = findpeaks(hBackStrip[i][j],3.,0.25);

			if (peaks.size() != 5) {
				printf("Strip %d on Detector %d not calibrated! Found %lu peaks instead of 5.\n", j, detID, peaks.size());
				txt_outfile << detID << "\t" << j << "\t" << -1 << "\t" << 0 << endl;	// Turns off channel if the number of peaks found is not 5.
			} else {
				//double *par = findcalib(peaks);	// Uses findcalib function to return a pointer to an array containing the calibration offset and slope for each channel

				vector < Double_t > energies = {5423., 5685., 6288., 6778, 8784.};
				vector < Double_t > energiesError = {1., 1., 1., 1., 1.};
				vector < Double_t > peaksError = {2., 2., 2., 2., 2.};
	
				//cStrip[i][j] = new TCanvas(Form("c%1.1i_%1.1i",i,j),Form("Strip[%1.1i][%1.1i]",i,j),30,113,800,600);

				//cStrip[i][j]->cd();

				gr[i][j] = new TGraphErrors(5, &peaks[0], &energies[0], &peaksError[0], &energiesError[0]);
				gr[i][j]->SetMarkerStyle(20);
				gr[i][j]->Draw();
				TAxis *Xaxis = gr[i][j]->GetXaxis();
   				Xaxis->SetLimits(0,4096);
				TAxis *Yaxis = gr[i][j]->GetYaxis();
   				Yaxis->SetLimits(0,10000);

				sprintf(gname,"gr[%d][%d]", i, j);
				gr[i][j]->Write(gname);

				sprintf(fname,"fit[%d][%d]", i, j);
				fit[i][j] = new TF1(fname,"pol1",0,4096);
				fit[i][j]->FixParameter(0,0);
				gr[i][j]->Fit(fit[i][j],"QR");
				fit[i][j]->Write(fname);
				static double par[2];
				fit[i][j]->GetParameters(par);

				double chisq = fit[i][j]->GetChisquare();
				int ndf = fit[i][j]->GetNDF();

				txt_outfile << detID << "\t" << j << "\t" << std::fixed << setprecision(5) << par[0] << "\t" << par[1] << "\t" << chisq/ndf << endl;
			}

		}

	}

	outfile->Write();
	outfile->Close();

}


void calibrateORRUBA()	{

	gROOT->SetBatch(kTRUE); // Set Batch mode so you don't get a shit load of canvases popping up

	// Create ROOT file
//	TFile *root_outfile = new TFile("calibrateORRUBA.root", "RECREATE");
/*
	// Make a TDirectory for each detector
	// QQQ5s Downstream
	TDirectory *dirQQQ5dAdE = root_outfile->mkdir("QQQ5dAdE");
	TDirectory *dirQQQ5dBdE = root_outfile->mkdir("QQQ5dBdE");
	TDirectory *dirQQQ5dAE1 = root_outfile->mkdir("QQQ5dAE1");
	TDirectory *dirQQQ5dBE1 = root_outfile->mkdir("QQQ5dBE1");
	TDirectory *dirQQQ5dAE2 = root_outfile->mkdir("QQQ5dAE2");
	TDirectory *dirQQQ5dBE2 = root_outfile->mkdir("QQQ5dBE2");
	// QQQ5s Upstream
	TDirectory *dirQQQ5uA = root_outfile->mkdir("QQQ5uA");
	TDirectory *dirQQQ5uB = root_outfile->mkdir("QQQ5uB");
	TDirectory *dirQQQ5uC = root_outfile->mkdir("QQQ5uC");
	TDirectory *dirQQQ5uD = root_outfile->mkdir("QQQ5uD");
	// BB10s Beam-Right
	TDirectory *dirBB10_1 = root_outfile->mkdir("BB10_1");
	TDirectory *dirBB10_2 = root_outfile->mkdir("BB10_2");
	TDirectory *dirBB10_3 = root_outfile->mkdir("BB10_3");
	TDirectory *dirBB10_4 = root_outfile->mkdir("BB10_4");
	// BB10s Beam-Left
	TDirectory *dirBB10_5 = root_outfile->mkdir("BB10_5");
	TDirectory *dirBB10_6 = root_outfile->mkdir("BB10_6");
	TDirectory *dirBB10_7 = root_outfile->mkdir("BB10_7");
	TDirectory *dirBB10_8 = root_outfile->mkdir("BB10_8");
	// SX3s Upstream
	TDirectory *dirSX3u_1 = root_outfile->mkdir("SX3u_1");
	TDirectory *dirSX3u_2 = root_outfile->mkdir("SX3u_2");
	TDirectory *dirSX3u_3 = root_outfile->mkdir("SX3u_3");
	TDirectory *dirSX3u_4 = root_outfile->mkdir("SX3u_4");
	TDirectory *dirSX3u_5 = root_outfile->mkdir("SX3u_5");
	TDirectory *dirSX3u_6 = root_outfile->mkdir("SX3u_6");
	TDirectory *dirSX3u_7 = root_outfile->mkdir("SX3u_7");
	TDirectory *dirSX3u_8 = root_outfile->mkdir("SX3u_8");
	TDirectory *dirSX3u_9 = root_outfile->mkdir("SX3u_9");
	TDirectory *dirSX3u_10 = root_outfile->mkdir("SX3u_10");
	TDirectory *dirSX3u_11 = root_outfile->mkdir("SX3u_11");
	TDirectory *dirSX3u_12 = root_outfile->mkdir("SX3u_12");
	// BB10s Downstream
	TDirectory *dirSX3d_1 = root_outfile->mkdir("SX3d_1");
	TDirectory *dirSX3d_2 = root_outfile->mkdir("SX3d_2");
	TDirectory *dirSX3d_3 = root_outfile->mkdir("SX3d_3");
	TDirectory *dirSX3d_4 = root_outfile->mkdir("SX3d_4");
	TDirectory *dirSX3d_5 = root_outfile->mkdir("SX3d_5");
	TDirectory *dirSX3d_6 = root_outfile->mkdir("SX3d_6");
	TDirectory *dirSX3d_7 = root_outfile->mkdir("SX3d_7");
	TDirectory *dirSX3d_8 = root_outfile->mkdir("SX3d_8");
	TDirectory *dirSX3d_9 = root_outfile->mkdir("SX3d_9");
	TDirectory *dirSX3d_10 = root_outfile->mkdir("SX3d_10");
	TDirectory *dirSX3d_11 = root_outfile->mkdir("SX3d_11");
	TDirectory *dirSX3d_12 = root_outfile->mkdir("SX3d_12");

	// Calibrate Downstream
	calibrateQQQ5("../protons/cal228th/cal228th_180deg_DS_noBB10/cal228th_180deg_DS_noBB10.root",497,561,"QQQ5dAdEcalib.dat", dirQQQ5dAdE);	// Downstream dE A
	calibrateQQQ5("../protons/cal228th/cal228th_180deg_DS_noBB10/cal228th_180deg_DS_noBB10.root",529,569,"QQQ5dBdEcalib.dat", dirQQQ5dBdE);	// Downstream dE B
	calibrateQQQ5("../protons/cal228th/cal228th_QQQ5_E1_a/cal228th_QQQ5_E1_a.root",577,673,"QQQ5dAE1calib.dat", dirQQQ5dAE1); 	// Downstream E1 A
	calibrateQQQ5("../protons/cal228th/cal228th_QQQ5_E1_a/cal228th_QQQ5_E1_a.root",609,681,"QQQ5dBE1calib.dat", dirQQQ5dBE1);	// Downstream E1 B
	calibrateQQQ5mod("../protons/cal228th/cal228th_QQQ5_E2/cal228th_QQQ5_E2.root",677,"QQQ5dAE2calib.dat", dirQQQ5dAE2);		// Downstream E2 A (only sectors instrumented)
	calibrateQQQ5("../protons/cal228th/cal228th_QQQ5_E2/cal228th_QQQ5_E2.root",641,685,"QQQ5dBE2calib.dat", dirQQQ5dBE2);		// Downstream E2 B

	// Calibrate Upstream

	calibrateQQQ5("../../protons/cal228th/cal228th_0deg_US/cal228th_0deg_US.root",1,129,"QQQ5uAcalib.dat", dirQQQ5uA);	// Downstream dE A
	calibrateQQQ5("../../protons/cal228th/cal228th_0deg_US/cal228th_0deg_US.root",33,133,"QQQ5uBcalib.dat", dirQQQ5uB);	// Downstream dE B
	calibrateQQQ5("../../protons/cal228th/cal228th_0deg_US/cal228th_0deg_US.root",65,137,"QQQ5uCcalib.dat", dirQQQ5uC); 	// Downstream E1 A
	calibrateQQQ5("../../protons/cal228th/cal228th_0deg_US/cal228th_0deg_US.root",97,141,"QQQ5uDcalib.dat", dirQQQ5uD);	// Downstream E1 B

	// Calibrate Beam-Right BB10s
	calibrateBB10("../protons/cal228th/cal228th_225deg_DS_BR/cal228th_225deg_DS_BR.root",433,"BB10_1_calib.dat", dirBB10_1);	// BB10 1
	calibrateBB10("../protons/cal228th/cal228th_225deg_DS_BR/cal228th_225deg_DS_BR.root",441,"BB10_2_calib.dat", dirBB10_2);	// BB10 2
	calibrateBB10("../protons/cal228th/cal228th_225deg_DS_BR/cal228th_225deg_DS_BR.root",449,"BB10_3_calib.dat", dirBB10_3); 	// BB10 3
	calibrateBB10("../protons/cal228th/cal228th_225deg_DS_BR/cal228th_225deg_DS_BR.root",457,"BB10_4_calib.dat", dirBB10_4);	// BB10 4

	// Calibrate Beam-Left BB10s
	calibrateBB10("../protons/cal228th/cal228th_135deg_DS_BL/cal228th_135deg_DS_BL.root",465,"BB10_5_calib.dat", dirBB10_5);	// BB10 1
	calibrateBB10("../protons/cal228th/cal228th_135deg_DS_BL/cal228th_135deg_DS_BL.root",473,"BB10_6_calib.dat", dirBB10_6);	// BB10 2
	calibrateBB10("../protons/cal228th/cal228th_135deg_DS_BL/cal228th_135deg_DS_BL.root",481,"BB10_7_calib.dat", dirBB10_7); 	// BB10 3
	calibrateBB10("../protons/cal228th/cal228th_135deg_DS_BL/cal228th_135deg_DS_BL.root",489,"BB10_8_calib.dat", dirBB10_8);	// BB10 4

	// Calibrate Downstream Beam-Right SX3s
	calibrateSX3("../protons/cal228th/cal228th_270deg_BR_noBB10/cal228th_270deg_BR_noBB10.root",385,"SX3d_1_calib.dat", dirSX3d_1);	// SX3d 1
	calibrateSX3("../protons/cal228th/cal228th_270deg_BR_noBB10/cal228th_270deg_BR_noBB10.root",389,"SX3d_2_calib.dat", dirSX3d_2);	// SX3d 2
	calibrateSX3("../protons/cal228th/cal228th_270deg_BR_noBB10/cal228th_270deg_BR_noBB10.root",393,"SX3d_3_calib.dat", dirSX3d_3); // SX3d 3
	calibrateSX3("../protons/cal228th/cal228th_270deg_BR_noBB10/cal228th_270deg_BR_noBB10.root",397,"SX3d_4_calib.dat", dirSX3d_4);	// SX3d 4
*/
	// Calibrate Upstream Beam-Right SX3s
//	calibrateSX3("../../OutputFolder/cal228th_315deg_US_BR.root","SX3_US_BR_calib.dat");	// SX3 Upstream Beam Right
	calibrateSX3("../../OutputFolder/cal228th_45deg_US_BL.root","SX3_US_BL_calib.dat");	// SX3 Upstream Beam Left
	//calibrateSX3("../OutputFolder/cal228th_270deg_BR_noBB10.root","SX3_DS_BR_calib.dat");	// SX3 Downstream Beam Right
	//calibrateSX3("../OutputFolder/cal228th_90deg_BL_noBB10.root","SX3_DS_BL_calib.dat");	// SX3 Downstream Beam Left

	// Write root file
//	root_outfile->cd();
//	root_outfile->Write();
//	root_outfile->Close();	

}
