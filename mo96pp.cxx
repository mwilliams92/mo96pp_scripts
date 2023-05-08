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
#include "TH2.h"
#include "TCanvas.h"
#include "TChain.h"
#include "analysis_functions.cxx"
#include "TCutG.h"
#include "TMath.h"

using namespace std;

void mo96pp() {

	TCutG *CD2_Barrel = new TCutG("CD2_Barrel",6);
  	CD2_Barrel->SetPoint(0,14252.14,2041.359);
   	CD2_Barrel->SetPoint(1,14188.03,3006.408);
   	CD2_Barrel->SetPoint(2,12008.55,3006.408);
   	CD2_Barrel->SetPoint(3,10512.82,2693.786);
   	CD2_Barrel->SetPoint(4,10448.72,2000.583);
   	CD2_Barrel->SetPoint(5,14252.14,2041.359);

	TCutG *protons_Barrel = new TCutG("protons_Barrel",6);
   	protons_Barrel->SetPoint(0,7931.658,2234.933);
   	protons_Barrel->SetPoint(1,2131.658,2284.698);
   	protons_Barrel->SetPoint(2,2816.583,2898.472);
   	protons_Barrel->SetPoint(3,10423.62,2837.647);
   	protons_Barrel->SetPoint(4,10394.47,2655.174);
   	protons_Barrel->SetPoint(5,7931.658,2234.933);

	TChain* Chain = new TChain("data");
  
	//Chain->Add("/home/williams328/experiments/mo96pp/OutputFolder/cal228th_315deg_US_BR.root");
    Chain->Add("/mnt/sandisk/files/sorted/Run0036.root");
    //Chain->Add("/mnt/d/sorted/Run0080.root");
	
	// ===================== Data variables =====================
	// 
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

	// === Analysis parameters ===
	double r2d = 180./ TMath::Pi();
	//double angle;
    double LabTheta;
	double LabTheta_dE_A;
	double LabTheta_dE_B;
    vector<double> hit_pos;
	vector<double> hit_pos_A;
	vector<double> hit_pos_B;	
	vector<double> spherical_polar_coord;
    double initial_energy;
	double initial_energy_A;
	double initial_energy_B;
	double angle_IC_corrected = 0;
	double Excitation = 0.0;
    double Excitation_A = 0.0;
    double Excitation_B = 0.0;
	double Etotal = 0.0;
	double Etotal_A = 0.0;
	double dE_A = 0.0;
	double E1_A = 0.0;
	double E2_A = 0.0;
	double Ethick_A = 0.0;
	double Etotal_B = 0.0;
	double dE_B = 0.0;
	double E1_B = 0.0;
	double E2_B = 0.0;
	double Ethick_B = 0.0;

	double BB10_Range = 0.0;
	double exponent = 1.6;
	double Range_A = 0.0;
	double Range_B = 0.0;

	double gain_CD2[4] = {0.939202, 1.0, 1.033238, 1.029436};

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

    //Output root file for histograms
    TFile write("./root_outputs/Kinematics_Histograms.root", "recreate");

	TH2D* hkinematics = new TH2D("kinematics", "kinematics", 360, 0, 180, 1000, 0, 40000);
	TH2D* hkinematics_gated = new TH2D("kinematics_gated", "kinematics gated", 360, 0, 180, 1000, 0, 40000);
	TH2D* hkinematics_US = new TH2D("hkinematics_US", "Kinematics Upstream", 360, 0, 180, 500, 0, 20000);
	TH2D* hkinematics_DS = new TH2D("hkinematics_DS", "Kinematics Downstream", 360, 0, 180, 500, 0, 20000);
	TH2D *hkinematics_DS_SX3 = new TH2D("hkinematics_DS_SX3", "Kinematics Downstream SX3", 360, 0, 180, 500, 0, 20000);

    TH1D* hExcitation = new TH1D("hExcitation", "Excitation Energy", 2500,-5,20);
	TH1D* hExcitation_DS = new TH1D("hExcitation_DS", "Excitation Energy (Downstream)", 2500,-5,20);
	TH1D* hExcitation_US = new TH1D("hExcitation_US", "Excitation Energy (Upstream)", 2500,-5,20);
	TH1D *hExcitation_US_SX3 = new TH1D("hExcitation_US_SX3", "Excitation Energy (Upstream SX3)", 2500,-5,20);
	TH1D *hExcitation_DS_SX3 = new TH1D("hExcitation_DS_SX3", "Excitation Energy (Downstream SX3)", 2500,-5,20);
	TH1D *hExcitation_US_Q5 = new TH1D("hExcitation_US_Q5", "Excitation Energy (Upstream QQQ5s)", 2500,-5,20);

    TH2D* hExcitation_vs_LabAngle_DS = new TH2D("hExcitation_vs_angle_DS", "Excitation Energy vs Angle (Downstream)", 360,0,180,2500,-5,20);
	TH2D* hExcitation_vs_LabAngle_US = new TH2D("hExcitation_vs_angle_US", "Excitation Energy vs Angle (Upstream)", 360,0,180,2500,-5,20);

	TH2D* hPID_Sector0 = new TH2D("hPID_Sector0", "Sector 0 #DeltaE-E", 2000, 0, 20000, 500, 0, 5000);
	TH2D* hPID_Sector1 = new TH2D("hPID_Sector1", "Sector 1 #DeltaE-E", 2000, 0, 20000, 500, 0, 5000);
	TH2D* hPID_Sector2 = new TH2D("hPID_Sector2", "Sector 2 #DeltaE-E", 2000, 0, 20000, 500, 0, 5000);
	TH2D* hPID_Sector3 = new TH2D("hPID_Sector3", "Sector 3 #DeltaE-E", 2000, 0, 20000, 500, 0, 5000);

	TH2D* gammaEx_matrix = new TH2D("hgammaEx_matrix", "Ex vs Gamma", 20000, 0, 20000, 2000, 0, 20000);

	TH2D* hBarrel_PID = new TH2D("hBarrel_PID", "Barrel PID", 2000, 0, 20000, 2000, 0, 20000);
    TH2D* hQQQ5a_PID = new TH2D("hQQQ5a_PID", "QQQ5a PID", 2000, 0, 20000, 2000, 0, 20000);
    TH2D* hQQQ5b_PID = new TH2D("hQQQ5b_PID", "QQQ5b PID", 2000, 0, 20000, 2000, 0, 20000);

	TH2D* hkinematics_DS_Det[12][4];
	TH2D* hkinematics_US_Det[12][4];
	TH2D* hkinematics_US_Q5[4];
	TH2D* hkinematics_DS_Q5[2];
	
	

    TH2D* hExcitationAng_US_Det[12][4];

	TH1D* hEtotal[12][4];
    TH1D* hExcitation_US_Det[12][4];
	TH1D* hExcitation_DS_Det[12][4];

	TH1D *hExcitation_DS_Q5[2];

	char hname[4096];
	//char hname2[4096];
	//char hname3[4096];

	for (int i=0; i<2; i++) {

		sprintf(hname,"hkinematics_DS_Q5[%d]", i);
		hkinematics_DS_Q5[i] = new TH2D(hname, hname, 360, 0, 180, 500, 0, 20000);
		sprintf(hname,"hExcitation_DS_Q5[%d]", i);
		hExcitation_DS_Q5[i] = new TH1D(hname, hname, 2500,-5,20);

	}

	for (int i=0; i<4; i++) {

		sprintf(hname,"hkinematics_US_Q5[%d]", i);
		hkinematics_US_Q5[i] = new TH2D(hname, hname, 360, 0, 180, 500, 0, 20000);

	}


	for (int i=0; i<12; i++) { 

		for (int j=0; j<4; j++) {

			sprintf(hname,"hkinematics_DS_Det[%d][%d]", i, j);
			hkinematics_DS_Det[i][j] = new TH2D(hname, hname, 360, 0, 180, 3000, 0, 30000);

			sprintf(hname,"hkinematics_US_Det[%d][%d]", i, j);
			hkinematics_US_Det[i][j] = new TH2D(hname, hname, 360, 0, 180, 3000, 0, 30000);

			sprintf(hname,"hEtotal[%d][%d]", i, j);
			hEtotal[i][j] = new TH1D(hname, hname, 500, 0, 20000);

            sprintf(hname,"hExcitationAng_US_Det[%d][%d]", i, j);
			hExcitationAng_US_Det[i][j] = new TH2D(hname, hname, 360, 0, 180, 2500, -5, 20);

            sprintf(hname,"hExcitation_US_Det[%d][%d]", i, j);
			hExcitation_US_Det[i][j] = new TH1D(hname, hname, 2500, -5, 20);

			sprintf(hname,"hExcitation_DS_Det[%d][%d]", i, j);
			hExcitation_DS_Det[i][j] = new TH1D(hname, hname, 2500, -5, 20);


		} 

	}

	//Getting the number of entries to loop through
	unsigned long long int nEntries = Chain->GetEntries();

	//Looping through each event:
	for ( unsigned long long int i=0; i<nEntries; i++ )	{
    	
		Chain->GetEntry(i);

		if(SX3Mul == 1) {

			for(int j=0; j<SX3Mul; j++) {

				hit_pos = hit_position_3D("SX3", SX3Upstream[j], SX3Det[j], SX3Strip[j], SX3StripPositionCal[j]);     
				initial_energy = initial_proton_energy((SX3SectorEnergy[j]/1000.0), proton_distance_through_target(hit_pos)); 

				spherical_polar_coord = hit_position_r_theta_phi("SX3", SX3Upstream[j], SX3Det[j], SX3Strip[j], SX3StripPositionCal[j]);
				LabTheta = spherical_polar_coord.at(1)*r2d;

                //Excitation = -rel_q_value_C12(LabTheta, SX3SectorEnergy[j]/1000);   
				Excitation = -rel_q_value_C12(LabTheta, initial_energy);      

				if (SX3Upstream[j]==1 && SX3StripRightEnergy[j]>0.0 && SX3StripLeftEnergy[j]>0.0)	{
            		
            		//Filling histograms
					hkinematics->Fill(LabTheta, SX3SectorEnergy[j]);
					hkinematics_US->Fill(LabTheta, SX3SectorEnergy[j]);
					hkinematics_US_Det[SX3Det[j]][SX3Sector[j]]->Fill(LabTheta, SX3SectorEnergy[j]);

                    hExcitation->Fill(Excitation);
                    hExcitation_US->Fill(Excitation);
                    hExcitation_US_Det[SX3Det[j]][SX3Sector[j]]->Fill(Excitation);
                    hExcitationAng_US_Det[SX3Det[j]][SX3Sector[j]]->Fill(LabTheta, Excitation);
					hExcitation_US_SX3->Fill(Excitation);

					hExcitation_vs_LabAngle_US->Fill(LabTheta,Excitation);

                }

				if(SX3Upstream[j] == 0) { // Detector 4 and 0 are weird!

					
					if (SX3Det[j] == 0 || SX3Det[j] == 5 || SX3Det[j] == 6)	{ // non telescope detectors

						hkinematics->Fill(LabTheta, SX3SectorEnergy[j]);
						hkinematics_DS->Fill(LabTheta, SX3SectorEnergy[j]);
						hkinematics_DS_Det[SX3Det[j]][SX3Sector[j]]->Fill(LabTheta, SX3SectorEnergy[j]);
						hkinematics_DS_SX3->Fill(LabTheta, SX3SectorEnergy[j]);
						hExcitation_DS_Det[SX3Det[j]][SX3Sector[j]]->Fill(Excitation);
                        hExcitation_DS_SX3->Fill(Excitation);
						hExcitation_DS->Fill(Excitation);
					}
					
					if(BB10Mul >= 1 && BB10Mul <= 1 && SX3Det[j]!=0 && SX3Det[j]!=5 && SX3Det[j]!=6) { // telescope detectors

						for ( int k=0; k<BB10Mul; k++ )	{ // do not want to include BB10Mul = 0 events

							if ( SX3Det[j] == BB10Det[k] ) {
							
								Etotal = SX3SectorEnergy[j] + BB10Energy[k];
								BB10_Range = pow((pow(Etotal,exponent) - pow(SX3SectorEnergy[j],exponent)),1/exponent);
    
								initial_energy = initial_proton_energy((Etotal/1000.0), proton_distance_through_target(hit_pos));
                               // Excitation = -rel_q_value_C12(LabTheta, Etotal/1000);
								Excitation = -rel_q_value_C12(LabTheta, initial_energy);

								hkinematics->Fill(LabTheta, Etotal);
								hBarrel_PID->Fill(Etotal,BB10_Range);

								//if (protons_Barrel->IsInside(Etotal,BB10_Range))	{
								if (BB10_Range > 2200. && BB10_Range<3000. && (Etotal > BB10_Range) )	{

									hkinematics_gated->Fill(LabTheta, Etotal);
									hkinematics_DS->Fill(LabTheta, Etotal);
									hkinematics_DS_SX3->Fill(LabTheta, Etotal);
									hkinematics_DS_Det[SX3Det[j]][SX3Sector[j]]->Fill(LabTheta, Etotal);
									hEtotal[SX3Det[j]][SX3Sector[j]]->Fill(Etotal);

                                    hExcitation->Fill(Excitation);
									hExcitation_DS->Fill(Excitation);
                                    hExcitation_DS_SX3->Fill(Excitation);
									hExcitation_DS_Det[SX3Det[j]][SX3Sector[j]]->Fill(Excitation);
                                    hExcitation_vs_LabAngle_DS->Fill(LabTheta,Excitation);

									
								}
							}
						} 
					}
				}
			}
		}
	
        if(QQQ5Mul >= 1 && QQQ5Mul<= 3) {

			Etotal_A = 0.0;
			Ethick_A = 0.0;
			dE_A = 0.0;
			E1_A = 0.0;
			E2_A = 0.0;
			LabTheta_dE_A = 0.0;
			Etotal_B = 0.0;
			Ethick_B = 0.0;
			dE_B = 0.0;
			E1_B = 0.0;
			E2_B = 0.0;	
			LabTheta_dE_B = 0.0;

            for(int j=0; j<QQQ5Mul; j++) {

				if (QQQ5Upstream[j] == 1) {

					spherical_polar_coord = hit_position_r_theta_phi("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);
					LabTheta = spherical_polar_coord.at(1)*r2d;

					hit_pos = hit_position_3D("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);     
					initial_energy = initial_proton_energy((QQQ5SectorEnergy[j]/1000.0), proton_distance_through_target(hit_pos)); 
                    //Excitation = -rel_q_value_C12(LabTheta, QQQ5SectorEnergy[j]/1000);
					Excitation = -rel_q_value_C12(LabTheta, initial_energy);

                	//Filling histograms
					hkinematics->Fill(LabTheta, QQQ5SectorEnergy[j]);
					hkinematics_US->Fill(LabTheta, QQQ5SectorEnergy[j]);
					hkinematics_US_Q5[QQQ5Det[j]]->Fill(LabTheta, QQQ5SectorEnergy[j]);

                    hExcitation->Fill(Excitation);
                    hExcitation_US->Fill(Excitation);
					hExcitation_US_Q5->Fill(Excitation);

					hExcitation_vs_LabAngle_US->Fill(LabTheta,Excitation);

				} 

	
				if ( QQQ5Upstream[j] == 0 ) {

					if ( QQQ5Det[j] == 0 && QQQ5Sector[j]!=0 && QQQ5Sector[j]!=3 ) {

						spherical_polar_coord = hit_position_r_theta_phi("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);
						LabTheta_dE_A = spherical_polar_coord.at(1)*r2d;
						dE_A = QQQ5SectorEnergy[j];
						hit_pos_A = hit_position_3D("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);     
						

					}

					if ( QQQ5Det[j] == 1 ) {

						spherical_polar_coord = hit_position_r_theta_phi("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);
						LabTheta_dE_B = spherical_polar_coord.at(1)*r2d;
						dE_B = QQQ5SectorEnergy[j];
						hit_pos_B = hit_position_3D("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);

					}

					if ( QQQ5Det[j] == 2 ) E1_A = QQQ5SectorEnergy[j];
					if ( QQQ5Det[j] == 3 ) E1_B = QQQ5SectorEnergy[j];
					
					if ( QQQ5Det[j] == 4 ) E2_A = QQQ5SectorEnergy[j];
					if ( QQQ5Det[j] == 5 ) E2_B = QQQ5SectorEnergy[j];

				}

			}

			Ethick_A = E1_A + E2_A;
			Etotal_A = 1.005*(dE_A + Ethick_A);
			Range_A = pow((pow(Etotal_A,exponent) - pow(Ethick_A,exponent)),1/exponent);
			
			Ethick_B = E1_B + E2_B;
			Etotal_B = dE_B + Ethick_B;
			Range_B = pow((pow(Etotal_B,exponent) - pow(Ethick_B,exponent)),1/exponent);
		

			if ( LabTheta_dE_A > 25 && LabTheta_dE_A < 45 && dE_A>0 && Ethick_A > 0 && Range_A > 2400 && Range_A < 3900) {

				initial_energy_A = initial_proton_energy((Etotal_A/1000.0), proton_distance_through_target(hit_pos_A));
            	//Excitation_A = -rel_q_value_C12(LabTheta_dE_A, Etotal_A/1000);
				Excitation_A = -rel_q_value_C12(LabTheta_dE_A, initial_energy_A);

				hkinematics->Fill(LabTheta_dE_A, Etotal_A);
				hkinematics_DS->Fill(LabTheta_dE_A, Etotal_A);
				hkinematics_DS_Q5[0]->Fill(LabTheta_dE_A, Etotal_A);
				hExcitation_DS->Fill(Excitation_A);
                hExcitation->Fill(Excitation_A);
                hExcitation_DS_Q5[0]->Fill(Excitation_A);
                hExcitation_vs_LabAngle_DS->Fill(LabTheta_dE_A,Excitation_A);
				

                hQQQ5a_PID->Fill(Etotal,Range_A);

			}

			if ( LabTheta_dE_B > 0 && dE_B>0 && Ethick_B > 0 && Range_B > 2400 && Range_B < 3900) {

				initial_energy_B = initial_proton_energy((Etotal_B/1000.0), proton_distance_through_target(hit_pos_B));
            	//Excitation_B = -rel_q_value_C12(LabTheta_dE_B, Etotal_B/1000);
				Excitation_B = -rel_q_value_C12(LabTheta_dE_B, initial_energy_B);

				hkinematics->Fill(LabTheta_dE_B, Etotal_B);
				hkinematics_DS->Fill(LabTheta_dE_B, Etotal_B);
				hkinematics_DS_Q5[1]->Fill(LabTheta_dE_B, Etotal_B);
				hExcitation_DS->Fill(Excitation_B);
                hExcitation->Fill(Excitation_B);
                hExcitation_DS_Q5[1]->Fill(Excitation_B);
                hExcitation_vs_LabAngle_DS->Fill(LabTheta_dE_B,Excitation_B);

                hQQQ5b_PID->Fill(Etotal,Range_B);
			}

		} 
		
		if (i % 10000 == 0)
      		cout << setiosflags(ios::fixed) << "Entry " << i << " of " << nEntries << ", " << 100 * i / nEntries << "% complete" << "\r" << flush; // Event counter

	}


	hkinematics->Write();
	hBarrel_PID->Write();
    hQQQ5a_PID->Write();
    hQQQ5b_PID->Write();
	hkinematics_gated->Write();

	hkinematics_US->Write();
	hExcitation_US_Q5->Write();
	hkinematics_DS->Write();
	hkinematics_DS_SX3->Write();
    hExcitation->Write();
	hExcitation_US->Write();
	hExcitation_DS->Write();
	hExcitation_US_SX3->Write();
	hExcitation_DS_SX3->Write();
	hExcitation_vs_LabAngle_US->Write();
	hExcitation_vs_LabAngle_DS->Write();


	for(int i=0; i<2; i++)	{

		hkinematics_DS_Q5[i]->Write();
		hExcitation_DS_Q5[i]->Write();

	}

	for(int i=0; i<4; i++)	{

		hkinematics_US_Q5[i]->Write();

	}

	for (int i=0; i<12; i++) {
		for (int j=0; j<4; j++) {
			hkinematics_DS_Det[i][j]->Write();
			hkinematics_US_Det[i][j]->Write();
			hEtotal[i][j]->Write();

            hExcitationAng_US_Det[i][j]->Write();
            hExcitation_US_Det[i][j]->Write();
			hExcitation_DS_Det[i][j]->Write();
		}	
	}

	cout << "\n" << "Done!" << endl;

	return;
}
