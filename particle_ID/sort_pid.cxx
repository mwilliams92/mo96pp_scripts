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

#include "../analysis_functions.cxx"

void sort_pid() {

	TChain* Chain = new TChain ("data");
  	//Chain->Add("../../OutputFolder/Run0080.root");
	
	Chain->Add("../../OutputFolder/Run0060_combined.root");
	//Chain->Add("../../OutputFolder/Run0061_combined.root");
	//Chain->Add("../../OutputFolder/Run0062_combined.root");
	//Chain->Add("../../OutputFolder/Run0063_combined.root");
	//Chain->Add("../../OutputFolder/Run0064_combined.root");
	
	// ===================== Data variables =====================
	// 
	//  === BB10 ===
	int   BB10Mul = 0;
    int   BB10Det[512] = {0};
	int   BB10Strip[512] = {0};
	int   BB10Channel[512] = {0};
	int   BB10ADC[512] = {0};
	float BB10Energy[512] = {0};
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

	// === GRETINA ===
    const Int_t NMAX = 44;
	bool  foundGRETINA = 0;
	int   xtalsMul = 0;
    float xtals_xlab[NMAX] = {0};
	float xtals_ylab[NMAX] = {0};
	float xtals_zlab[NMAX] = {0};
	float xtals_cc[NMAX] = {0};
	float xtals_edop[NMAX] = {0};
	float xtals_edopMaxInt[NMAX] = {0};
	float xtals_edopSeg[NMAX] = {0};
	float xtals_edopXtal[NMAX] = {0};
	int   xtals_crystalNum[NMAX] = {0};
	int   xtals_quadNum[NMAX] = {0};
	float xtals_t0[NMAX] = {0};
	long long  xtals_timestamp[NMAX] = {0};

	// === Analysis parameters ===
	double r2d = 180./ TMath::Pi();
	//double angle;
    double LabTheta;
    vector<double> hit_pos;
	vector<double> spherical_polar_coord;
    double initial_energy;
	double angle_IC_corrected = 0;
	double excitation = 0.0;
	float Barrel_Etotal = 0.0;
	float BB10_Range = 0.0;
	float BB10_Range_Shifted = 0.0;
	double QQQ5_Exponent = 1.7;
	double Barrel_Exponent = 1.6;

	double Etotal, dE, Ethick, Range, E1, E2;

    //============================================================
    //   Allocating the branch addresses of the "raw" variables
    //============================================================

    // ================ BB10 Branch Addresses ===================
    Chain->SetBranchAddress("BB10Mul",&BB10Mul);
    Chain->SetBranchAddress("BB10Det",&BB10Det);
    Chain->SetBranchAddress("BB10Strip",&BB10Strip);
    Chain->SetBranchAddress("BB10Channel",&BB10Channel);
    Chain->SetBranchAddress("BB10ADC",&BB10ADC);
    Chain->SetBranchAddress("BB10Energy",&BB10Energy);

    // ================ QQQ5 Branch Addresses ===================
    Chain->SetBranchAddress("QQQ5Mul",&QQQ5Mul);
    Chain->SetBranchAddress("QQQ5Upstream",&QQQ5Upstream);
    Chain->SetBranchAddress("QQQ5Det",&QQQ5Det);
    Chain->SetBranchAddress("QQQ5Ring",&QQQ5Ring);
    Chain->SetBranchAddress("QQQ5RingChannel",&QQQ5RingChannel);
    Chain->SetBranchAddress("QQQ5Sector",&QQQ5Sector);
    Chain->SetBranchAddress("QQQ5SectorChannel",&QQQ5SectorChannel);
    Chain->SetBranchAddress("QQQ5RingADC",&QQQ5RingADC);
    Chain->SetBranchAddress("QQQ5RingEnergy",&QQQ5RingEnergy);
    Chain->SetBranchAddress("QQQ5SectorADC",&QQQ5SectorADC);
    Chain->SetBranchAddress("QQQ5SectorEnergy",&QQQ5SectorEnergy);
    Chain->SetBranchAddress("QQQ5Angle",&QQQ5Angle);

    // =================== SX3 Branch Address ==================
    Chain->SetBranchAddress("SX3Mul",&SX3Mul);
    Chain->SetBranchAddress("SX3Upstream",&SX3Upstream);
    Chain->SetBranchAddress("SX3Det",&SX3Det);
    Chain->SetBranchAddress("SX3Sector",&SX3Sector);
    Chain->SetBranchAddress("SX3SectorChannel",&SX3SectorChannel);
    Chain->SetBranchAddress("SX3SectorADC",&SX3SectorADC);
    Chain->SetBranchAddress("SX3SectorEnergy",&SX3SectorEnergy);
    Chain->SetBranchAddress("SX3Strip",&SX3Strip);
    Chain->SetBranchAddress("SX3StripLeftChannel",&SX3StripLeftChannel);
    Chain->SetBranchAddress("SX3StripRightChannel",&SX3StripRightChannel);
    Chain->SetBranchAddress("SX3StripLeftADC",&SX3StripLeftADC);
    Chain->SetBranchAddress("SX3StripRightADC",&SX3StripRightADC);
	Chain->SetBranchAddress("SX3StripLeftEnergy",&SX3StripLeftEnergy);
    Chain->SetBranchAddress("SX3StripRightEnergy",&SX3StripRightEnergy);
    Chain->SetBranchAddress("SX3StripEnergy",&SX3StripEnergy);
	Chain->SetBranchAddress("SX3StripPosition",&SX3StripPosition);
	Chain->SetBranchAddress("SX3StripPositionCal",&SX3StripPositionCal);

	// ================= GRETINA Branch Address ================
    Chain->SetBranchAddress("foundGRETINA",&foundGRETINA);
    Chain->SetBranchAddress("xtalsMul",&xtalsMul);
    Chain->SetBranchAddress("xtals_xlab",xtals_xlab);
    Chain->SetBranchAddress("xtals_ylab",xtals_ylab);
    Chain->SetBranchAddress("xtals_zlab",xtals_zlab);
    Chain->SetBranchAddress("xtals_cc",xtals_cc);
    Chain->SetBranchAddress("xtals_edop",xtals_edop);
    Chain->SetBranchAddress("xtals_edopMaxInt",xtals_edopMaxInt);
    Chain->SetBranchAddress("xtals_edopSeg",xtals_edopSeg);
    Chain->SetBranchAddress("xtals_edopXtal",xtals_edopXtal);
    Chain->SetBranchAddress("xtals_crystalNum",xtals_crystalNum);
    Chain->SetBranchAddress("xtals_quadNum",xtals_quadNum);
    Chain->SetBranchAddress("xtals_t0",xtals_t0);
    Chain->SetBranchAddress("xtals_timestamp",xtals_timestamp);

    //Output root file for histograms
    TFile write("PID_Histograms.root", "recreate");

	TH2F* hBarrel_PID[12]; // each sector and detector combination
	//TH2F* hBarrel_PID_2[4];		// All detectors for a given sector
	//TH2F* hBarrel_PID_3[12];	// All sectors for each detector (likely to be shit)


	TH2F* hBarrel_PID_Corrected = new TH2F("hBarrel_PID_Corrected", "Corrected PID", 1000, 0, 20000, 750, 0, 15000);		// Everything
	TH2F* hBarrel_PID_NonCorrected = new TH2F("hBarrel_PID_NonCorrected", "Non-Corrected PID", 1000, 0, 20000, 750, 0, 15000);

	TH2F* hQQQ5_PID = new TH2F("hQQQ5_PID", "QQQ5 PID", 3000, 0, 30000, 3000, 0, 30000);

	TH1D* hBarrel_Range[12];
	TH1D* hBarrel_Range_Corrected = new TH1D("hBarrel_Range_Corrected", "Corrected Range", 750, 0, 15000);
	TH1D* hBarrel_Range_NonCorrected = new TH1D("hBarrel_Range_NonCorrected", "Non-Corrected Range", 750, 0, 15000);

	TH2D* hQQQ5a_PID = new TH2D("hQQQ5a_PID","QQQ5 A-stack PID",4000,0,40000,2000,0,20000);
	TH2D* hQQQ5b_PID = new TH2D("hQQQ5b_PID","QQQ5 B-stack PID",4000,0,40000,2000,0,20000);

	TH2D* hQQQ5a_dE_E = new TH2D("hQQQ5a_dE_E","QQQ5 A-stack #DeltaE-E",4000,0,40000,2000,0,20000);
	TH2D* hQQQ5b_dE_E = new TH2D("hQQQ5b_dE_E","QQQ5 B-stack #DeltaE-E",4000,0,40000,2000,0,20000);

	char hname[4096], hname2[4096];

	for (int i=0; i<12; i++) { 

		//for (int j=0; j<4; j++) {

			sprintf(hname,"hBarrel_PID[%d]", i);
			sprintf(hname2,"hBarrel_Range[%d]", i);
			hBarrel_PID[i] = new TH2F(hname, hname, 1000, 0, 20000, 750, 0, 15000);
			hBarrel_Range[i] = new TH1D(hname2, hname2, 750, 0, 15000);

		//}

	}

	double shift[12] = { 0.0, -140., -81, 0.0, 118., 0.0, 0.0, -117., -99., -159., 161., 0.0 };
	double gain[12] = { 0.0, 0.949146, 0.969933, 1.0, 1.047295, 0.0, 0.0, 0.957143, 0.963496, 0.942641, 1.065661, 0.0 };
	double gain_CD2[4] = {0.939202, 1.0, 1.033238, 1.029436};

	//Getting the number of entries to loop through
	unsigned long long int nEntries = Chain->GetEntries();

	//Looping through each event:
	for ( unsigned long long int i=0; i<nEntries; i++ )	{

	Chain->GetEntry(i);

		if(SX3Mul >= 1 && SX3Mul <= 1 && BB10Mul >= 1 && BB10Mul <= 1) {

			for(int j=0; j<SX3Mul; j++) {  

				if(SX3Upstream[j] == 0) {

					for( int k=0; k<BB10Mul; k++ ) {

						if ( SX3Det[j] == BB10Det[k] ) {

							Barrel_Etotal = BB10Energy[k] + SX3SectorEnergy[j];
							BB10_Range = pow((pow(Barrel_Etotal,Barrel_Exponent) - pow(SX3SectorEnergy[j],Barrel_Exponent)),1/Barrel_Exponent);

							BB10_Range_Shifted = gain[SX3Det[j]]*BB10_Range;

							hBarrel_PID[SX3Det[j]]->Fill(Barrel_Etotal ,BB10_Range);

							hBarrel_PID_Corrected->Fill(Barrel_Etotal ,BB10_Range_Shifted);
			
							hBarrel_PID_NonCorrected->Fill(Barrel_Etotal ,BB10_Range);

							hBarrel_Range[SX3Det[j]]->Fill(BB10_Range);
							hBarrel_Range_Corrected->Fill(BB10_Range_Shifted);
							hBarrel_Range_NonCorrected->Fill(BB10_Range);

						}
					}	
				}
			}
		}

        if(QQQ5Mul >= 1 && QQQ5Mul<= 3) {

			Etotal = 0.0;
			dE = 0.0;
			Ethick = 0.0;
			Range = 0.0;
			E1 = 0.0;
			E2 = 0.0;
								

            for(int j=0; j<QQQ5Mul; j++) {
				
				if ( QQQ5Upstream[j] == 0) {

					if ( QQQ5Det[j] == 0 || QQQ5Det[j] == 1 ) {
						dE = QQQ5SectorEnergy[j];
					}
					if ( QQQ5Det[j] == 2 || QQQ5Det[j] == 3 ) {
						E1 = QQQ5SectorEnergy[j];
					}
					if ( QQQ5Det[j] == 4 || QQQ5Det[j] == 5 ) {
						E2 = QQQ5SectorEnergy[j];
					}
				}
			}

			Ethick = E1 + E2;
			Etotal = dE + E1 + E2;
			Range = pow((pow(Etotal,QQQ5_Exponent) - pow(Ethick,QQQ5_Exponent)),1/QQQ5_Exponent);
			
			hQQQ5a_PID->Fill(Etotal, Range);

			hQQQ5a_dE_E->Fill(Etotal, dE);

		} 
		
		if (i % 10000 == 0)
      		cout << setiosflags(ios::fixed) << "Entry " << i << " of " << nEntries << ", " << 100 * i / nEntries << "% complete" << "\r" << flush; // Event counter

	}

	for (int i=0; i<12; i++) {
		//for (int j=0; j<4; j++) {
			hBarrel_PID[i]->Write();
			hBarrel_Range[i]->Write();
		//}
	}

	hBarrel_PID_Corrected->Write();
	hBarrel_PID_NonCorrected->Write();
	hBarrel_Range_Corrected->Write();
	hBarrel_Range_NonCorrected->Write();

	hQQQ5a_PID->Write();
	hQQQ5b_PID->Write();

	hQQQ5a_dE_E->Write();
	hQQQ5b_dE_E->Write();

	cout << "\n" << "Done!" << endl;

	return;
}
