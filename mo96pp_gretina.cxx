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

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TChain.h"
#include "analysis_functions.cxx"

void mo96pp_gretina() {

	int firstRun = 54;
	int lastRun = 59;
    int nruns = 1 + lastRun - firstRun;

	char inputFileName[512];

    //  === BB10 ===
	int   BB10Mul = 0;
    int   BB10Det[512] = {0};
	int   BB10Strip[512] = {0};
	int   BB10Channel[512] = {0};
	int   BB10ADC[512] = {0};
	float BB10Energy[512] = {0};
	float En_BB10 = 0;
	//
	// === QQQ5 ===
	int   QQQ5Mul = 0;
	int   QQQ5Upstream[512] = {0};
	int   QQQ5Det[512] = {0};
	int   QQQ5Ring[512] = {0};
	int   QQQ5RingChannel[512] = {0};
	int   QQQ5Sector[512] = {0};
	int   QQQ5SectorChannel[512] = {0};
	int   QQQ5RingADC[512] = {0};
	float QQQ5RingEnergy[512] = {0};
	int   QQQ5SectorADC[512] = {0};
	float QQQ5SectorEnergy[512] = {0};
	float QQQ5Angle[512] = {0};
	//
	// === SX3 ===
	int   SX3Mul = 0;
	int   SX3Upstream[36] = {0};
	int   SX3Det[36] = {0};
	int   SX3Sector[36] = {0};
	int   SX3SectorChannel[36] = {0};
	int   SX3SectorADC[36] = {0};
	float SX3SectorEnergy[36] = {0};
	int   SX3Strip[36] = {0};
	int   SX3StripLeftChannel[36] = {0};
	int   SX3StripRightChannel[36] = {0};
	int   SX3StripLeftADC[36] = {0};
	int   SX3StripRightADC[36] = {0};
	float SX3StripLeftEnergy[36] = {0};
	float SX3StripRightEnergy[36] = {0};
	float SX3StripEnergy[36] = {0};
	float SX3StripPosition[36] = {0};
	float SX3StripPositionCal[36] = {0};
	
	// === TDC ===
	int   tdcSilicon = 0;
	int   tdcGRETINA = 0;
	int   tdcRF = 0;
    long long timeStamp = 0;
	unsigned long long GRETINATimeStamp = 0;

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

    int n;

    long long TS_Diff;

    TFile *outFile = new TFile("CheckGretina.root","RECREATE");

    TH2F *hCrys_v_TS[nruns];
    TH2F *hDiffTS_v_TS[nruns];
    TH2F *hDiffTS_v_TS_ORRUBA[nruns];
    TH2F *hEgam1_v_Crys[nruns];
    
    //TH1D *hMultiRaw = new TH1D("hMultiRaw","Multiplicity (xtalsMul)",100,0,100);
    //TH1D *hMulti = new TH1D("hMulti","Multiplicity (sorted)",100,0,100);    
    //TH2F *hEgam1_v_Crys = new TH2F("hEgam_v_Crys","Gamma vs Crystal",10000,0,10000,500,0,500);
    //TH2D *hQuad_v_Crys = new TH2D("hQuad_v_Crys","Quad vs Crystal",50,0,50,500,0,500);


    char hname[32], hname2[32], hname3[32], hname4[32], htitle[32];

    for (int i=0; i<nruns; i++) {

        sprintf(hname, "hCrys_v_TS[%d]", i);
        sprintf(hname2, "hDiffTS_v_TS[%d]", i);
        sprintf(hname3, "hDiffTS_v_TS_ORRUBA[%d]", i);
        sprintf(htitle, "Run %d", i+firstRun);
        hCrys_v_TS[i] = new TH2F(hname,htitle,10000,0,10000,500,0,500);
        hDiffTS_v_TS[i] = new TH2F(hname2,htitle,10000,0,100000,1000,0,1000);
        hDiffTS_v_TS_ORRUBA[i] = new TH2F(hname3,htitle,10000,0,100000,1000,0,1000);

        sprintf(hname4, "hEgam1_v_Crys[%d]", i);
        hEgam1_v_Crys[i] = new TH2F(hname4,htitle,10000,0,10000,50,0,50);

    }

    for(int r=firstRun; r<=lastRun; r++) { // loop over runs

        n = r - firstRun; // n runs from first run

        if ( r == 79 || r == 80 ) { continue; } // skip CD2 runs

        //sprintf(inputFileName,"/mnt/sandisk/files/sorted/Run00%d_combined.root",r); // open sorted root file.
        sprintf(inputFileName,"/mnt/d/sorted/Run00%d_combined.root",r); // open sorted root file.


	    TFile *inputFile = TFile::Open(inputFileName);
	    TTree *data = (TTree*)inputFile->Get("data_combined");

	    cout << "\n" << "Sorting Run " << r << endl;

        //============================================================
        //   Allocating the branch addresses of the "raw" variables
        //============================================================
/*
        // ================ BB10 Branch Addresses ===================
        data->SetBranchAddress("BB10Mul",&BB10Mul);
        data->SetBranchAddress("BB10Det",&BB10Det);
        data->SetBranchAddress("BB10Strip",&BB10Strip);
        data->SetBranchAddress("BB10Channel",&BB10Channel);
        data->SetBranchAddress("BB10ADC",&BB10ADC);
        data->SetBranchAddress("BB10Energy",&BB10Energy);

        // ================ QQQ5 Branch Addresses ===================
        data->SetBranchAddress("QQQ5Mul",&QQQ5Mul);
        data->SetBranchAddress("QQQ5Upstream",&QQQ5Upstream);
        data->SetBranchAddress("QQQ5Det",&QQQ5Det);
        data->SetBranchAddress("QQQ5Ring",&QQQ5Ring);
        data->SetBranchAddress("QQQ5RingChannel",&QQQ5RingChannel);
        data->SetBranchAddress("QQQ5Sector",&QQQ5Sector);
        data->SetBranchAddress("QQQ5SectorChannel",&QQQ5SectorChannel);
        data->SetBranchAddress("QQQ5RingADC",&QQQ5RingADC);
        data->SetBranchAddress("QQQ5RingEnergy",&QQQ5RingEnergy);
        data->SetBranchAddress("QQQ5SectorADC",&QQQ5SectorADC);
        data->SetBranchAddress("QQQ5SectorEnergy",&QQQ5SectorEnergy);
        data->SetBranchAddress("QQQ5Angle",&QQQ5Angle);
    
        // =================== SX3 Branch Address ==================
        data->SetBranchAddress("SX3Mul",&SX3Mul);
        data->SetBranchAddress("SX3Upstream",SX3Upstream);
        data->SetBranchAddress("SX3Det",SX3Det);
        data->SetBranchAddress("SX3Sector",SX3Sector);
        data->SetBranchAddress("SX3SectorChannel",SX3SectorChannel);
        data->SetBranchAddress("SX3SectorADC",SX3SectorADC);
        data->SetBranchAddress("SX3SectorEnergy",SX3SectorEnergy);
        data->SetBranchAddress("SX3Strip",SX3Strip);
        data->SetBranchAddress("SX3StripLeftChannel",SX3StripLeftChannel);
        data->SetBranchAddress("SX3StripRightChannel",SX3StripRightChannel);
        data->SetBranchAddress("SX3StripLeftADC",SX3StripLeftADC);
        data->SetBranchAddress("SX3StripRightADC",SX3StripRightADC);
	    data->SetBranchAddress("SX3StripLeftEnergy",SX3StripLeftEnergy);
        data->SetBranchAddress("SX3StripRightEnergy",SX3StripRightEnergy);
        data->SetBranchAddress("SX3StripEnergy",SX3StripEnergy);
	    data->SetBranchAddress("SX3StripPosition",SX3StripPosition);
	    data->SetBranchAddress("SX3StripPositionCal",SX3StripPositionCal);
*/
	    // =================== TDC Branch Address ==================
        data->SetBranchAddress("tdcSilicon",&tdcSilicon);
        data->SetBranchAddress("tdcGRETINA",&tdcGRETINA);
        data->SetBranchAddress("tdcRF",&tdcRF);
        data->SetBranchAddress("timeStamp",&timeStamp);
//        data->SetBranchAddress("GRETINATimeStamp",&GRETINATimeStamp);
        
	    // ================= GRETINA Branch Address ================
        data->SetBranchAddress("xtalsMul",&xtalsMul);
        data->SetBranchAddress("xtals_xlab",&xtals_xlab);
        data->SetBranchAddress("xtals_ylab",&xtals_ylab);
        data->SetBranchAddress("xtals_zlab",&xtals_zlab);
        data->SetBranchAddress("xtals_cc",&xtals_cc);
        data->SetBranchAddress("xtals_cc1",&xtals_cc1);
        data->SetBranchAddress("xtals_cc2",&xtals_cc2);
        data->SetBranchAddress("xtals_cc3",&xtals_cc3);
        data->SetBranchAddress("xtals_cc4",&xtals_cc4);
        data->SetBranchAddress("xtals_edop",&xtals_edop);
        data->SetBranchAddress("xtals_edopMaxInt",&xtals_edopMaxInt);
        data->SetBranchAddress("xtals_edopSeg",&xtals_edopSeg);
        data->SetBranchAddress("xtals_edopXtal",&xtals_edopXtal);
        data->SetBranchAddress("xtals_crystalNum",&xtals_crystalNum);
        data->SetBranchAddress("xtals_quadNum",&xtals_quadNum);
        data->SetBranchAddress("xtals_t0",&xtals_t0);
        data->SetBranchAddress("xtals_timestamp",&xtals_timestamp);
	
	    //Getting the number of entries to loop through
	    unsigned long long int nEntries = data->GetEntries();

	    //Looping through each event:
	    for ( unsigned long long int i=0; i<nEntries; i++ )	{
    	
		    data->GetEntry(i);
               
            //hMulti->Fill(xtalsMul);

		    for (int l=0; l<xtalsMul; l++) {


                hEgam1_v_Crys[n]->Fill(xtals_cc1[l],xtals_crystalNum[l]);
                //hQuad_v_Crys->Fill(xtals_quadNum[l],xtals_crystalNum[l]);

                TS_Diff = timeStamp - xtals_timestamp[l];

               // cout << xtals_t0[l] << endl;

                hCrys_v_TS[n]->Fill(xtals_timestamp[l]/1e8,xtals_crystalNum[l]);
                hDiffTS_v_TS[n]->Fill(xtals_timestamp[l]/1e8,TS_Diff);
                hDiffTS_v_TS_ORRUBA[n]->Fill(timeStamp/1e8,TS_Diff);
                                        
                
            }
                  
		    
            if (i % 10000 == 0)
      		    cout << setiosflags(ios::fixed) << "Entry " << i << " of " << nEntries << ", " << 100 * i / nEntries << "% complete" << "\r" << flush; // Event counter
    }
	
	cout << "\n" << "Finished Run " << r << endl;

       
    //hMulti->Write();        
    //hEgam_v_Crys->Write();
    //hQuad_v_Crys->Write();

    inputFile->Close();
              

	}

    outFile->cd();

    for (int i=0; i<nruns; i++) {

        hEgam1_v_Crys[i]->Write();
        hCrys_v_TS[i]->Write();
        hDiffTS_v_TS[i]->Write();
        hDiffTS_v_TS_ORRUBA[i]->Write();

    }
    
    outFile->Close();
	return;
}
