#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>
#include <string>
#include <vector>
#include <map>
#include <random>
#include <cassert>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <list>
#include <functional>
#include <algorithm>
#include <iomanip>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TChain.h"
#include "analysis_functions.cxx"

using namespace std;

void absolute_efficiency() {

	// === GRETINA ===
    const Int_t NMAX = 44;
	bool  foundGRETINA = 0;
	int   xtalsMul = 0;
    float xtals_xlab[NMAX] = {0};
	float xtals_ylab[NMAX] = {0};
	float xtals_zlab[NMAX] = {0};
	float xtals_cc[NMAX] = {0};
    float xtals_cc1[NMAX] = {0};
    float xtals_cc2[NMAX] = {0};
    float xtals_cc3[NMAX] = {0};
    float xtals_cc4[NMAX] = {0};
	float xtals_edop[NMAX] = {0};
	float xtals_edopMaxInt[NMAX] = {0};
	float xtals_edopSeg[NMAX] = {0};
	float xtals_edopXtal[NMAX] = {0};
	int   xtals_crystalNum[NMAX] = {0};
	int   xtals_quadNum[NMAX] = {0};
	float xtals_t0[NMAX] = {0};
	long long  xtals_timestamp[NMAX] = {0};
    
	//TFile *inputFile = TFile::Open("/mnt/sandisk/files/sorted/Run0098_combined.root");
    TFile *inputFile = TFile::Open("/mnt/d/sorted/Run0098_combined.root");
	TTree *data = (TTree*)inputFile->Get("data_combined");

    TFile *outFile = new TFile("Efficiency.root","RECREATE");
    TH2F *hCrys_v_TS = new TH2F("hCrys_v_TS","Crystal vs TS",5000,0,5000,50,0,50);
    TH1F *hEgam = new TH1F("hEgam","All Singles (except Crystal 2)",3000,0,3000);
 /*   
    TH1F *hEgam1 = new TH1F("hEgam1","Crystal 1 Spectrum",3000,0,3000);
    TH2F *hCrys_v_Egam1 = new TH2F("hCrys_v_Egam1","Crystal vs Egam1",5000,0,5000,50,0,50);
    TH1F *hEgam2 = new TH1F("hEgam2","Coincidences in Every Other Crystal",3000,0,3000);
    TH2F *hCrys_v_Egam2 = new TH2F("hCrys_v_Egam2","Crystal vs Egam2",5000,0,5000,50,0,50);
    TH1F *TS_Diff = new TH1F("hTS_Diff","Time Stamp Diff Between 60Co gammas",2000,-1000,1000);
*/

    TH1F *hEgam1[49];
    TH2F *hCrys_v_Egam1[49];
    TH1F *hEgam2[49];
    TH2F *hCrys_v_Egam2[49];
    TH1F *TS_Diff[49];
    TH1F *TS_Diff2[49];
    TH1F *TS_Diff3[49];

    char hname[32];

    for(int n=0; n<=48; n++) {

        sprintf(hname,"hEgam1[%d]",n);
        hEgam1[n] = new TH1F(hname,"Crystal 1 Spectrum",3000,0,3000);
        sprintf(hname,"hCrys_v_Egam1[%d]",n);
        hCrys_v_Egam1[n] = new TH2F(hname,"Crystal vs Egam1",5000,0,5000,50,0,50);
        sprintf(hname,"hEgam2[%d]",n);
        hEgam2[n] = new TH1F(hname,"Coincidences in Every Other Crystal",3000,0,3000);
        sprintf(hname,"hCrys_v_Egam2[%d]",n);
        hCrys_v_Egam2[n] = new TH2F(hname,"Crystal vs Egam2",5000,0,5000,50,0,50);
        sprintf(hname,"hTS_Diff[%d]",n);
        TS_Diff[n] = new TH1F(hname,"Time Stamp Diff Between gammas",2000,-1000,1000);
        sprintf(hname,"hTS_Diff2[%d]",n);
        TS_Diff2[n] = new TH1F(hname,"Time Stamp Diff Between gammas Hitting Different Crsytal",2000,-1000,1000);
        sprintf(hname,"hTS_Diff3[%d]",n);
        TS_Diff3[n] = new TH1F(hname,"Time Stamp Diff Between 60Co gammas",2000,-1000,1000);


    }

	// ================= GRETINA Branch Address ================
    data->SetBranchAddress("xtalsMul",&xtalsMul);
    data->SetBranchAddress("xtals_xlab",xtals_xlab);
    data->SetBranchAddress("xtals_ylab",xtals_ylab);
    data->SetBranchAddress("xtals_zlab",xtals_zlab);
    data->SetBranchAddress("xtals_cc",xtals_cc);
    data->SetBranchAddress("xtals_cc1",xtals_cc1);
    data->SetBranchAddress("xtals_cc2",xtals_cc2);
    data->SetBranchAddress("xtals_cc3",xtals_cc3);
    data->SetBranchAddress("xtals_cc4",xtals_cc4);
    data->SetBranchAddress("xtals_edop",xtals_edop);
    data->SetBranchAddress("xtals_edopMaxInt",xtals_edopMaxInt);
    data->SetBranchAddress("xtals_edopSeg",xtals_edopSeg);
    data->SetBranchAddress("xtals_edopXtal",xtals_edopXtal);
    data->SetBranchAddress("xtals_crystalNum",xtals_crystalNum);
    data->SetBranchAddress("xtals_quadNum",xtals_quadNum);
    data->SetBranchAddress("xtals_t0",xtals_t0);
    data->SetBranchAddress("xtals_timestamp",xtals_timestamp);
	
    //Getting the number of entries to loop through
	unsigned long long int nEntries = data->GetEntries();

    cout << "\n" << "Sorting..." << endl;

    float counter1, counter2;
    long long deltaTS;

     counter1 = 0.;
     counter2 = 0.;

    // 1165 -> 1180
    // 1320 -> 1340

	//Looping through each event:
	for ( unsigned long long int i=0; i<nEntries; i++ )	{
    	
        data->GetEntry(i);

		for (int l=0; l<xtalsMul; l++) { // Loop over all gamma-ray events

            hCrys_v_TS->Fill(xtals_timestamp[l]/1e9, xtals_crystalNum[l]);
            if (xtals_crystalNum[l] != 2 ) hEgam->Fill(xtals_cc1[l]);

            for (int n=1; n<=48; n++) {

                if (xtals_crystalNum[l] == n ) { // find an event in crystal 1

                    hEgam1[n]->Fill(xtals_cc1[l]); // fill spectrum for crystal 1
                    hCrys_v_Egam1[n]->Fill(xtals_cc1[l],xtals_crystalNum[l]);

                    if (xtals_cc1[l]>1165 && xtals_cc1[l]<1180.) { // gate on 1333-keV peak in crystal 1

                        counter1++;

                        for (int k=0; k<xtalsMul; k++) { // If the above condition is satisfied then loop over all gamma entries again 

                            deltaTS = abs(xtals_timestamp[k] - xtals_timestamp[l]);
                            TS_Diff[n]->Fill(deltaTS);

                            if (xtals_crystalNum[k] != n && xtals_crystalNum[k] != 2) { 

                                hEgam2[n]->Fill(xtals_cc1[k]); // plot the sum of all events in every other crystal
                                hCrys_v_Egam2[n]->Fill(xtals_cc1[k],xtals_crystalNum[k]);
                                TS_Diff2[n]->Fill(deltaTS);


                                if (xtals_cc1[k]>1320. && xtals_cc1[k]<1340.) { 
                                    TS_Diff3[n]->Fill(deltaTS);
                                    counter2++; 
                                } 
                            }
                        }
                    }
                }
            }
        }
		    
        if (i % 10000 == 0)
            cout << setiosflags(ios::fixed) << "Entry " << i << " of " << nEntries << ", " << 100 * i / nEntries << "% complete" << "\r" << flush; // Event counter
    }
	
	cout << "\n" << "Finished" << endl;
    
    double efficiency;

    efficiency = (43./42.) * counter2 / counter1; // rough efficiency. correction for 1 crystal used to gate on 1333 keV.

    cout << "Efficiency Estimate = " << efficiency << endl;

    hCrys_v_TS->Write();
    hEgam->Write();

    for(int n=0; n<=48; n++) {

        hEgam1[n]->Write();
        hCrys_v_Egam1[n]->Write();
        hEgam2[n]->Write();
        hCrys_v_Egam2[n]->Write();
        TS_Diff[n]->Write();
        TS_Diff2[n]->Write();
        TS_Diff3[n]->Write();


    }
    
    outFile->Close();
    inputFile->Close();

	return;
}
