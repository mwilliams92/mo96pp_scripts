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

void plotKinematics() {

	TChain* Chain = new TChain ("data");
    Chain->Add("../OutputFolder/Run0030.root");
	Chain->Add("../OutputFolder/Run0031.root");
	Chain->Add("../OutputFolder/Run0032.root");
	Chain->Add("../OutputFolder/Run0033.root");
	Chain->Add("../OutputFolder/Run0034.root");
	Chain->Add("../OutputFolder/Run0035.root");
	Chain->Add("../OutputFolder/Run0037.root");
	Chain->Add("../OutputFolder/Run0038.root");
	Chain->Add("../OutputFolder/Run0039.root");
	Chain->Add("../OutputFolder/Run0040.root");
	Chain->Add("../OutputFolder/Run0041.root");
	Chain->Add("../OutputFolder/Run0042.root");
	Chain->Add("../OutputFolder/Run0043.root");
	Chain->Add("../OutputFolder/Run0044.root");
	Chain->Add("../OutputFolder/Run0045.root");
	Chain->Add("../OutputFolder/Run0046.root");
	Chain->Add("../OutputFolder/Run0047.root");
	Chain->Add("../OutputFolder/Run0048.root");
	Chain->Add("../OutputFolder/Run0049.root");

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
	int   QQQ5Upstream = 0;
	int   QQQ5Det[512] = {0};
	int   QQQ5Ring = 0;
	int   QQQ5RingChannel = 0;
	int   QQQ5Sector = 0;
	int   QQQ5SectorChannel = 0;
	int   QQQ5RingADC = 0;
	float QQQ5RingEnergy = 0;
	int   QQQ5SectorADC = 0;
	float QQQ5SectorEnergy[512] = {0};
	float QQQ5Angle = 0;
	float deltaE = 0;
	float E1 = 0;
	float E2 = 0;
	float Etot = 0;
	//
	// === SX3 ===
	int   SX3Mul = 0;
	int   SX3Upstream = 0;
	int   SX3Det = 0;
	int   SX3Sector = 0;
	int   SX3SectorChannel = 0;
	int   SX3SectorADC = 0;
	float SX3SectorEnergy[512] = {0};
	int   SX3Strip = 0;
	int   SX3StripLeftChannel = 0;
	int   SX3StripRightChannel = 0;
	int   SX3StripLeftADC = 0;
	int   SX3StripRightADC = 0;
	float SX3StripLeftEnergy = 0;
	float SX3StripRightEnergy = 0;
	float SX3StripEnergy = 0;
	float SX3StripPosition = 0;
	float SX3StripPositionCal = 0;
	float En_SX3 = 0;
	
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
    TFile write("Kinematics_Histograms.root", "recreate");

	//Getting the number of entries to loop through
	unsigned long long int nEntries = Chain->GetEntries();

	//Looping through each event:
	for ( unsigned long long int i=0; i<nEntries; i++ )	{
    	
		Chain->GetEntry(i);

		if(SX3Mul >= 1 && SX3Mul <= 5) {

			for(int j=0; j<SX3Mul; j++) {

        		//Applying the proton energy loss through target correcttion.
            	hit_pos = hit_position_3D("SX3", SX3Upstream[j], SX3Det[j], SX3Strip[j], SX3Position[j]);
            	initial_energy = initial_proton_energy((BSX3_En/1000.0), proton_distance_through_target(hit_pos));
                    
            	//Calculating the ejectile angle
            	angle = getSX3Angle(SX3Upstream[j], SX3Position[j]);
                    
            	//Filling histograms
            	kinematics->Fill(angle, SX3SectorEnergy[j]);
                
				if(SX3Upstream[j] == 1) {
                      
					//Testing function to calculate r, theta, phi
                	spherical_polar_coord = hit_position_r_theta_phi("SX3", SX3Upstream[j], SX3Det[j], SX3Strip[j], SX3Position[j]);
                	angle_vs_strip->Fill(spherical_polar_coord.at(2)*(180./3.14159), spherical_polar_coord.at(1)*(180./3.14159));
                    	
					// Coincidences with the dE layer (downstream barrel ONLY!)
                	if(BB10Mul >= 1) {
                        	
						for(int l=0; l<BB10Mul; l++) {

                    		if(DS_protons->IsInside(BSX3_En, BB10_En)) {
                            
								kinematics_proton_PID->Fill(angle, BSX3_En+BB10_En);

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
