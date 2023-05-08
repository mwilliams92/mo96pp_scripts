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

void mo96pp_DS_QQQ5s() {

	TChain* Chain = new TChain ("data");
  
	//Chain->Add("../OutputFolder/cal228th_180deg_DS_noBB10.root");
	//Chain->Add("../OutputFolder/cal228th_QQQ5_E1_a.root");
	//Chain->Add("../OutputFolder/cal228th_QQQ5_E2.root");

	//Chain->Add("../OutputFolder/Run0080.root");

	Chain->Add("../OutputFolder/Run0060.root");
	Chain->Add("../OutputFolder/Run0061.root");
	Chain->Add("../OutputFolder/Run0062.root");
	Chain->Add("../OutputFolder/Run0063.root");
	Chain->Add("../OutputFolder/Run0064.root");
	
	// ===================== Data variables ===================
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

	// === Analysis parameters ===
	double r2d = 180./ TMath::Pi();
	//double angle;
    double LabTheta, LabTheta_dE, LabTheta_E1, LabTheta_E2;
    vector<double> hit_pos;
	vector<double> spherical_polar_coord;
    double initial_energy;
	double excitation = 0.0;

	//double Etotal[512] = {0};
	//double dE[512] = {0};
	//double E1[512] = {0};
	//double E2[512] = {0};

	double Etotal = 0;
	double Ethick = 0;
	double Range = 0;
	double dE = 0;
	double E1 = 0;
	double E2 = 0;
	double exponent = 1.7;
	
    //============================================================
    //   Allocating the branch addresses of the "raw" variables
    //============================================================

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

    //Output root file for histograms
    TFile write("QQQ5_Hists.root", "recreate");

	TH2D *hdE = new TH2D("hdE","LabEnergy vs Angle",90,0,90,2000,0,20000);
	TH2D *hE1 = new TH2D("hE1","LabEnergy vs Angle",90,0,90,2000,0,20000);
	TH2D *hE2 = new TH2D("hE2","LabEnergy vs Angle",90,0,90,2000,0,20000);
	TH2D *hEtotal = new TH2D("hEtotal","LabEnergy vs Angle",90,0,90,2000,0,20000);

	TH2D *hPID1 = new TH2D("hPID1","#DeltaE vs E total",2000,0,20000,2000,0,20000);
	TH2D *hPID2 = new TH2D("hPID2","Range vs E total",2000,0,20000,2000,0,20000);

	TH2D *hPID1_gated = new TH2D("hPID1_gated","#DeltaE vs E total",2000,0,20000,2000,0,20000);
	TH2D *hPID2_gated = new TH2D("hPID2_gated","Range vs E total",2000,0,20000,2000,0,20000);

	TH2D *hExcite = new TH2D("hExcite","Excitation vs LabAngle",90,0,90,2500,-5,20);

	//Getting the number of entries to loop through
	unsigned long long int nEntries = Chain->GetEntries();

	//Looping through each event:
	for ( unsigned long long int i=0; i<nEntries; i++ )	{
    	
		Chain->GetEntry(i);

		Etotal=0.0;
		dE=0.0;
		Ethick=0.0;
	
        if(QQQ5Mul >= 1 && QQQ5Mul<= 6) {

		//if(QQQ5Mul == 6) {
	
            for(int j=0; j<QQQ5Mul; j++) {
	
				if ( QQQ5Upstream[j] == 0 && (QQQ5Det[j] % 2 != 0) ) {

					Etotal += QQQ5SectorEnergy[j];

				}

				if (QQQ5Upstream[j] == 0 && QQQ5Det[j] == 1) {

					spherical_polar_coord = hit_position_r_theta_phi("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);
					LabTheta_dE = spherical_polar_coord.at(1)*r2d;
					dE += QQQ5SectorEnergy[j];

				}

				if (QQQ5Upstream[j] == 0 && (QQQ5Det[j] == 3 || QQQ5Det[j] == 5)) {

					spherical_polar_coord = hit_position_r_theta_phi("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);
					LabTheta_dE = spherical_polar_coord.at(1)*r2d;
					Ethick += QQQ5SectorEnergy[j];

				}

			}
		
			//Ethick = Etotal - dE;
			Range = pow((pow(Etotal,exponent) - pow(Ethick,exponent)),1/exponent);

			hEtotal->Fill(LabTheta_dE, Etotal);
			
			
			if (Ethick > 0) {
				hPID2->Fill(Etotal,Range);
				hPID1->Fill(Etotal,dE);
			}

		}
/*
		if(SX3Mul >= 1 && SX3Mul <= 1) {

			for(int j=0; j<SX3Mul; j++) {

				hit_pos = hit_position_3D("SX3", SX3Upstream[j], SX3Det[j], SX3Strip[j], SX3StripPositionCal[j]);     
				initial_energy = initial_proton_energy((SX3SectorEnergy[j]/1000.0), proton_distance_through_target(hit_pos)); 

				spherical_polar_coord = hit_position_r_theta_phi("SX3", SX3Upstream[j], SX3Det[j], SX3Strip[j], SX3StripPositionCal[j]);
				LabTheta = spherical_polar_coord.at(1)*r2d;	           

				if (SX3Upstream[j]==1 && SX3StripRightEnergy[j]>0.0)	{
            		
            		//Filling histograms
					hkinematics->Fill(LabTheta, SX3SectorEnergy[j]);
					hkinematics_US->Fill(LabTheta, SX3SectorEnergy[j]);
					hkinematics_US_Det[SX3Det[j]][SX3Sector[j]]->Fill(LabTheta, SX3SectorEnergy[j]);

                }

				if(SX3Upstream[j] == 0) { // Detector 4 and 0 are weird!

					
					if (SX3Det[j] == 0 || SX3Det[j] == 5 || SX3Det[j] == 6)	{ // non telescope detectors

						hkinematics->Fill(LabTheta, SX3SectorEnergy[j]);
						hkinematics_DS->Fill(LabTheta, SX3SectorEnergy[j]);
						hkinematics_DS_Det[SX3Det[j]][SX3Sector[j]]->Fill(LabTheta, SX3SectorEnergy[j]);

					}
					
					if(BB10Mul >= 1 && BB10Mul <= 1 && SX3Det[j]!=0 && SX3Det[j]!=5 && SX3Det[j]!=6) { // telescope detectors

						for ( int k=0; k<=BB10Mul; k++ )	{ // do not want to include BB10Mul = 0 events

							if ( SX3Det[j] == BB10Det[k] ) {
							
								Etotal = SX3SectorEnergy[j] + BB10Energy[k];
								BB10_Range = gain_CD2[SX3Sector[j]]*pow((pow(Etotal,exponent) - pow(SX3SectorEnergy[j],exponent)),1/exponent);

								if (protons_Barrel->IsInside(Etotal,BB10_Range))	{

									hkinematics->Fill(LabTheta, Etotal);
									hkinematics_DS->Fill(LabTheta, Etotal);
									hkinematics_DS_Det[SX3Det[j]][SX3Sector[j]]->Fill(LabTheta, Etotal);
									hEtotal[SX3Det[j]][SX3Sector[j]]->Fill(Etotal);
								}

							}
						} 
					}
				}
			}
		}*/
		
		if (i % 10000 == 0)
      		cout << setiosflags(ios::fixed) << "Entry " << i << " of " << nEntries << ", " << 100 * i / nEntries << "% complete" << "\r" << flush; // Event counter

	}

	hdE->Write();
	hE1->Write();	
	hE2->Write();
	hEtotal->Write();
	hExcite->Write();

	hPID1->Write();
	hPID2->Write();

	hPID1_gated->Write();
	hPID2_gated->Write();

	cout << "\n" << "Done!" << endl;

	return;
}
