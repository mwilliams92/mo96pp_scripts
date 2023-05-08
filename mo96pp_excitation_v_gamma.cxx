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

void mo96pp_excitation_v_gamma() {

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


	TChain* Chain = new TChain ("data");
   
	//Chain->Add("../OutputFolder/Run0080.root");
	Chain->Add("../OutputFolder/Run0060_combined.root");
	Chain->Add("../OutputFolder/Run0061_combined.root");
	Chain->Add("../OutputFolder/Run0062_combined.root");
	Chain->Add("../OutputFolder/Run0063_combined.root");
	Chain->Add("../OutputFolder/Run0064_combined.root");
	Chain->Add("../OutputFolder/Run0065_combined.root");
	//Chain->Add("../OutputFolder/Run0066_combined.root");
	//Chain->Add("../OutputFolder/Run0067_combined.root");
	//Chain->Add("../OutputFolder/Run0068_combined.root");
	//Chain->Add("../OutputFolder/Run0069_combined.root");
	//Chain->Add("../OutputFolder/Run0070_combined.root");
	//Chain->Add("../OutputFolder/Run0071_combined.root");
	//Chain->Add("../OutputFolder/Run0072_combined.root");
	//Chain->Add("../OutputFolder/Run0073_combined.root");
	//Chain->Add("../OutputFolder/Run0074_combined.root");
	//Chain->Add("../OutputFolder/Run0075_combined.root");
	//Chain->Add("../OutputFolder/Run0076_combined.root");
	//Chain->Add("../OutputFolder/Run0077_combined.root");
	//Chain->Add("../OutputFolder/Run0078_combined.root");

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

	// === TDC ===
	int   tdcSilicon = 0;
	int   tdcGRETINA = 0;
	int   tdcRF = 0;
	unsigned long long timeStamp = 0;
	unsigned long long GRETINATimeStamp = 0;

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
	double LabTheta_QQQ5;
    vector<double> hit_pos;
	vector<double> spherical_polar_coord;
    double initial_energy;
	double angle_IC_corrected = 0;
	double excitation = 0.0;
	double Etotal = 0.0;
	double dE = 0.0;
	double E1 = 0.0;
	double E2 = 0.0;
	double Ethick = 0.0;
	double BB10_Range = 0.0;
	double Range = 0.0;

	double Barrel_Exponent = 1.6;
	double QQQ5_Exponent = 1.7;

	double gain_CD2[4] = {0.939202, 1.0, 1.033238, 1.029436};

    //============================================================
    //   Allocating the branch addresses of the "raw" variables
    //============================================================

    // ================ BB10 Branch Addresses ===================
    Chain->SetBranchAddress("BB10Mul",&BB10Mul);
    Chain->SetBranchAddress("BB10Det",BB10Det);
    Chain->SetBranchAddress("BB10Strip",BB10Strip);
    Chain->SetBranchAddress("BB10Channel",BB10Channel);
    Chain->SetBranchAddress("BB10ADC",BB10ADC);
    Chain->SetBranchAddress("BB10Energy",BB10Energy);

    // ================ QQQ5 Branch Addresses ===================
    Chain->SetBranchAddress("QQQ5Mul",&QQQ5Mul);
    Chain->SetBranchAddress("QQQ5Upstream",QQQ5Upstream);
    Chain->SetBranchAddress("QQQ5Det",QQQ5Det);
    Chain->SetBranchAddress("QQQ5Ring",QQQ5Ring);
    Chain->SetBranchAddress("QQQ5RingChannel",QQQ5RingChannel);
    Chain->SetBranchAddress("QQQ5Sector",QQQ5Sector);
    Chain->SetBranchAddress("QQQ5SectorChannel",QQQ5SectorChannel);
    Chain->SetBranchAddress("QQQ5RingADC",QQQ5RingADC);
    Chain->SetBranchAddress("QQQ5RingEnergy",QQQ5RingEnergy);
    Chain->SetBranchAddress("QQQ5SectorADC",QQQ5SectorADC);
    Chain->SetBranchAddress("QQQ5SectorEnergy",QQQ5SectorEnergy);
    Chain->SetBranchAddress("QQQ5Angle",QQQ5Angle);

    // =================== SX3 Branch Address ==================
    Chain->SetBranchAddress("SX3Mul",&SX3Mul);
    Chain->SetBranchAddress("SX3Upstream",SX3Upstream);
    Chain->SetBranchAddress("SX3Det",SX3Det);
    Chain->SetBranchAddress("SX3Sector",SX3Sector);
    Chain->SetBranchAddress("SX3SectorChannel",SX3SectorChannel);
    Chain->SetBranchAddress("SX3SectorADC",SX3SectorADC);
    Chain->SetBranchAddress("SX3SectorEnergy",SX3SectorEnergy);
    Chain->SetBranchAddress("SX3Strip",SX3Strip);
    Chain->SetBranchAddress("SX3StripLeftChannel",SX3StripLeftChannel);
    Chain->SetBranchAddress("SX3StripRightChannel",SX3StripRightChannel);
    Chain->SetBranchAddress("SX3StripLeftADC",SX3StripLeftADC);
    Chain->SetBranchAddress("SX3StripRightADC",SX3StripRightADC);
	Chain->SetBranchAddress("SX3StripLeftEnergy",SX3StripLeftEnergy);
    Chain->SetBranchAddress("SX3StripRightEnergy",SX3StripRightEnergy);
    Chain->SetBranchAddress("SX3StripEnergy",SX3StripEnergy);
	Chain->SetBranchAddress("SX3StripPosition",SX3StripPosition);
	Chain->SetBranchAddress("SX3StripPositionCal",SX3StripPositionCal);

	// =================== TDC Branch Address ==================
    Chain->SetBranchAddress("tdcSilicon",&tdcSilicon);
    Chain->SetBranchAddress("tdcGRETINA",&tdcGRETINA);
    Chain->SetBranchAddress("tdcRF",&tdcRF);
    Chain->SetBranchAddress("timeStamp",&timeStamp);
    Chain->SetBranchAddress("GRETINATimeStamp",&GRETINATimeStamp);

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
    TFile write("Excitation_v_Gamma_Histograms.root", "recreate");

	CD2_Barrel->Write();

	TH2D* hExcite = new TH2D("hExcite", "Excitation Energy", 360, 0, 180, 2000,-5,15);
	TH2D* hExcite_US = new TH2D("hExcite_US", "Excitation Energy Upstream", 360, 0, 180, 2000,-5,15);
	TH2D* hExcite_DS = new TH2D("hExcite_DS", "Excitation Energy Downstream", 360, 0, 180, 2000,-5,15);

	TH2D* hExcite_v_Gamma = new TH2D("hExcite_v_Gamma", "Excitation vs Egamma", 8000, 0, 16000, 2000,-5,15);

	char hname[4096];

	//Getting the number of entries to loop through
	unsigned long long int nEntries = Chain->GetEntries();

	//Looping through each event:
	for ( unsigned long long int i=0; i<nEntries; i++ )	{
    	
		Chain->GetEntry(i);

		if ( tdcSilicon>750 && tdcSilicon<850 ) {

				//if (xtals_cc[k]>773 && xtals_cc[k]<783) {

					if(SX3Mul >= 1 && SX3Mul <= 1) {

						for(int j=0; j<SX3Mul; j++) {

							hit_pos = hit_position_3D("SX3", SX3Upstream[j], SX3Det[j], SX3Strip[j], SX3StripPositionCal[j]);     
							initial_energy = initial_proton_energy((SX3SectorEnergy[j]/1000.0), proton_distance_through_target(hit_pos));

							spherical_polar_coord = hit_position_r_theta_phi("SX3", SX3Upstream[j], SX3Det[j], SX3Strip[j], SX3StripPositionCal[j]);
							LabTheta = spherical_polar_coord.at(1)*r2d;

							/*
							if (SX3Upstream[j]==1 && SX3StripRightEnergy[j]>0.0 && SX3SectorEnergy[j]>0.0)	{
            		
								excitation = -rel_q_value(LabTheta, initial_energy);

            					//Filling histograms
								hExcite->Fill(LabTheta, excitation);
								hExcite_US->Fill(LabTheta, excitation);

               				}
							*/
				
							if(SX3Upstream[j] == 0) { // Detector 4 and 0 are weird!

							/*
								if (SX3Det[j] == 5 || SX3Det[j] == 6)	{ // non telescope detectors

								//excitation = -rel_q_value(LabTheta, 0.001*SX3SectorEnergy[j]);
								excitation = -rel_q_value(LabTheta, initial_energy);

								//hExcite->Fill(LabTheta, excitation);
								//hExcite_DS->Fill(LabTheta, excitation);

							}
							*/
					
							if(BB10Mul >= 1 && BB10Mul <= 1 && SX3Det[j]!=0 && SX3Det[j]!=5 && SX3Det[j]!=6) { // telescope detectors

								for ( int l=0; l<=BB10Mul; l++ )	{ // do not want to include BB10Mul = 0 events

										if ( SX3Det[j] == BB10Det[l] ) {

											hit_pos = hit_position_3D("SX3", SX3Upstream[j], SX3Det[j], SX3Strip[j], SX3StripPositionCal[j]);     

											spherical_polar_coord = hit_position_r_theta_phi("SX3", SX3Upstream[j], SX3Det[j], SX3Strip[j], SX3StripPositionCal[j]);
											LabTheta = spherical_polar_coord.at(1)*r2d;
							
											Etotal = SX3SectorEnergy[j] + BB10Energy[l];

											initial_energy = initial_proton_energy((Etotal/1000.0), proton_distance_through_target(hit_pos));
											excitation = -rel_q_value(LabTheta, initial_energy);

											BB10_Range = pow((pow(Etotal,Barrel_Exponent) - pow(SX3SectorEnergy[j],Barrel_Exponent)),1/Barrel_Exponent);

											if (BB10_Range > 2200. && BB10_Range<3000. && (Etotal > BB10_Range) && Etotal> 0 && Etotal<8000.)	{

												hExcite->Fill(LabTheta, excitation);
												hExcite_DS->Fill(LabTheta, excitation);

												if(xtalsMul >= 1 && xtalsMul <= 10)	{
                        	
													for(int k=0; k<xtalsMul; k++)	{

														if ( tdcGRETINA > 1400 ) {
															hExcite_v_Gamma->Fill(xtals_cc[k], excitation);
														}
													}
												}
											}
										}
									} 
								}
							}
						}
					}
		
        			if(QQQ5Mul >= 1 && QQQ5Mul<= 3) {

						Etotal = 0.0;
						dE = 0.0;
						E1 = 0.0;
						E2 = 0.0;
						Ethick = 0.0;

            			for(int j=0; j<QQQ5Mul; j++) {
						/*
							if (QQQ5Upstream[j] == 1) {

								spherical_polar_coord = hit_position_r_theta_phi("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);
								LabTheta = spherical_polar_coord.at(1)*r2d;

								excitation = -rel_q_value(LabTheta, 0.001*QQQ5SectorEnergy[j]);

                				//Filling histograms
								hExcite->Fill(LabTheta, excitation);
								hExcite_US->Fill(LabTheta, excitation);
							} 
						*/
				
							if ( QQQ5Upstream[j] == 0) {

								if ( QQQ5Det[j] == 1) { // Use front QQQ5s to set the angle

									hit_pos = hit_position_3D("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);
									spherical_polar_coord = hit_position_r_theta_phi("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);
									LabTheta_QQQ5 = spherical_polar_coord.at(1)*r2d;

									dE = QQQ5SectorEnergy[j];

								}

								if ( QQQ5Det[j] == 3) {
									E1 = QQQ5SectorEnergy[j];
								}

								if ( QQQ5Det[j] == 5) {
									E2 = QQQ5SectorEnergy[j];
								}
							}
						}

						Ethick = E1 + E2;
						Etotal = dE + Ethick;
						Range = pow((pow(Etotal,QQQ5_Exponent) - pow(Ethick,QQQ5_Exponent)),1/QQQ5_Exponent);

						initial_energy = initial_proton_energy((Etotal/1000.0), proton_distance_through_target(hit_pos));
						excitation = -rel_q_value(LabTheta, initial_energy);

						if ( Range > 2700. && Range < 4200. && ( Etotal > Range ) && Etotal> 0 && Etotal<8000.) {

							hExcite->Fill(LabTheta_QQQ5, excitation);
							hExcite_DS->Fill(LabTheta_QQQ5, excitation);

							if(xtalsMul >= 1 && xtalsMul <= 10)	{
                        	
								for(int k=0; k<xtalsMul; k++)	{

									if ( tdcGRETINA > 1400 ) {
										hExcite_v_Gamma->Fill(xtals_cc[k], excitation);

									}

								}
							}
						}
					}
				} 
		
		if (i % 10000 == 0)
      		cout << setiosflags(ios::fixed) << "Entry " << i << " of " << nEntries << ", " << 100 * i / nEntries << "% complete" << "\r" << flush; // Event counter
	}

	hExcite->Write();
	//hExcite_US->Write();
	hExcite_DS->Write();
	hExcite_v_Gamma->Write();
	
	cout << "\n" << "Done!" << endl;

	return;
}
