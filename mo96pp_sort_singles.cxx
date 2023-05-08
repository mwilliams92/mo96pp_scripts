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
#include "TH2D.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TCutG.h"
#include "analysis_functions.cxx"

void mo96pp_sort_singles() {

	int firstRun = 59;	// Min = 59
	int lastRun = 68;	// Max = 93
    int nruns = lastRun - firstRun;

	char inputFileName[512];

    //Double_t TS_Gate[35] = {5000., 50000., 50000., 50000., 27700., 1570., 50000., 50000., 50000., 24360., 600., 50000., 50000., 50000., 50000., 50000., 1850., 210., 50000., 50000., 50000., 50000., 50000., 50000., 50000., 9200., 50000., 50000., 50000., 25900., 0., 550., 50000., 50000., 50000.};
	//Double_t TS_Gate[23] = {50000., 50000., 50000., 50000., 1850., 210., 50000., 50000., 50000., 50000., 50000., 50000., 50000., 9200., 50000., 50000., 50000., 25900., 0., 550., 50000., 50000., 50000.};
	Double_t TS_Gate[10] = {5000., 50000., 50000., 50000., 27700., 1570., 50000., 50000., 50000., 24360.};

    int r;

	for(int n=0; n<=nruns; n++) {

    r = n + firstRun;

    if ( r==79 || r==80 || r == 89 ) {continue;}

	//sprintf(inputFileName,"/mnt/d/sorted/Run00%d_combined.root",r);
    sprintf(inputFileName,"/mnt/sandisk/files/sorted/Run00%d_combined.root",r);

	TFile *inputFile = TFile::Open(inputFileName);
	TTree *data = (TTree*)inputFile->Get("data_combined");

	cout << "\n" << "Sorting Run " << r << endl;

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
    int tdcSilicon_Div = 0;
    int tdcSilicon_GRETINA = 0;
    int tdcSilicon_Delay = 0;
    int tdcSilicon_Upstream = 0;

	Long64_t timeStamp = 0;
	unsigned long long GRETINATimeStamp = 0;

	// === GRETINA ===
    const Int_t NMAX = 128;
	bool  foundGRETINA = 0;
	int   xtalsMul = 0;
    float xtals_xlab[NMAX] = {0};
	float xtals_ylab[NMAX] = {0};
	float xtals_zlab[NMAX] = {0};
	float xtals_cc[NMAX] = {0};
    float xtals_cc1[NMAX] = {0};
	float xtals_edop[NMAX] = {0};
	float xtals_edopMaxInt[NMAX] = {0};
	float xtals_edopSeg[NMAX] = {0};
	float xtals_edopXtal[NMAX] = {0};
	int   xtals_crystalNum[NMAX] = {0};
	int   xtals_quadNum[NMAX] = {0};
	float xtals_t0[NMAX] = {0};
	Long64_t  xtals_timestamp[NMAX] = {0};

	// === Analysis parameters ===
	double r2d = 180./ TMath::Pi();

    vector<double> hit_pos;
	vector<double> spherical_polar_coord;

    double initial_energy;
	double LabTheta, LabTheta_A, LabTheta_B;
	double Etotal;
	double dEa, E1a, E2a, Etotal_A, Ethick_A, Range_A;
	double dEb, E1b, E2b, Etotal_B, Ethick_B, Range_B;
	double BB10_Range;
	double distance;
    double Excitation, Excitation_Cor, Excitation_O16, Excitation_C12;

	double Barrel_Exponent = 1.6;
	double QQQ5_Exponent = 1.7;

	double Ring_A, Ring_B;


    //============================================================
    //   Allocating the branch addresses of the "raw" variables
    //============================================================

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
    data->SetBranchAddress("SX3Upstream",&SX3Upstream);
    data->SetBranchAddress("SX3Det",&SX3Det);
    data->SetBranchAddress("SX3Sector",&SX3Sector);
    data->SetBranchAddress("SX3SectorChannel",&SX3SectorChannel);
    data->SetBranchAddress("SX3SectorADC",&SX3SectorADC);
    data->SetBranchAddress("SX3SectorEnergy",&SX3SectorEnergy);
    data->SetBranchAddress("SX3Strip",&SX3Strip);
    data->SetBranchAddress("SX3StripLeftChannel",&SX3StripLeftChannel);
    data->SetBranchAddress("SX3StripRightChannel",&SX3StripRightChannel);
    data->SetBranchAddress("SX3StripLeftADC",&SX3StripLeftADC);
    data->SetBranchAddress("SX3StripRightADC",&SX3StripRightADC);
	data->SetBranchAddress("SX3StripLeftEnergy",&SX3StripLeftEnergy);
    data->SetBranchAddress("SX3StripRightEnergy",&SX3StripRightEnergy);
    data->SetBranchAddress("SX3StripEnergy",&SX3StripEnergy);
	data->SetBranchAddress("SX3StripPosition",&SX3StripPosition);
	data->SetBranchAddress("SX3StripPositionCal",&SX3StripPositionCal);

	// =================== TDC Branch Address ==================
    data->SetBranchAddress("tdcSilicon",&tdcSilicon);
    data->SetBranchAddress("tdcGRETINA",&tdcGRETINA);
    data->SetBranchAddress("tdcRF",&tdcRF);
    data->SetBranchAddress("timeStamp",&timeStamp);
    //data->SetBranchAddress("GRETINATimeStamp",&GRETINATimeStamp);

    data->SetBranchAddress("tdcSilicon_Div", &tdcSilicon_Div);
    data->SetBranchAddress("tdcSilicon_GRETINA", &tdcSilicon_GRETINA);
    data->SetBranchAddress("tdcSilicon_Delay", &tdcSilicon_Delay);
    data->SetBranchAddress("tdcSilicon_Upstream", &tdcSilicon_Upstream);

	// ================= GRETINA Branch Address ================
    data->SetBranchAddress("foundGRETINA",&foundGRETINA);
    data->SetBranchAddress("xtalsMul",&xtalsMul);
    data->SetBranchAddress("xtals_xlab",&xtals_xlab);
    data->SetBranchAddress("xtals_ylab",&xtals_ylab);
    data->SetBranchAddress("xtals_zlab",&xtals_zlab);
    data->SetBranchAddress("xtals_cc",&xtals_cc);
    data->SetBranchAddress("xtals_cc1",&xtals_cc1);
    data->SetBranchAddress("xtals_edop",&xtals_edop);
    data->SetBranchAddress("xtals_edopMaxInt",&xtals_edopMaxInt);
    data->SetBranchAddress("xtals_edopSeg",&xtals_edopSeg);
    data->SetBranchAddress("xtals_edopXtal",&xtals_edopXtal);
    data->SetBranchAddress("xtals_crystalNum",&xtals_crystalNum);
    data->SetBranchAddress("xtals_quadNum",&xtals_quadNum);
    data->SetBranchAddress("xtals_t0",&xtals_t0);
    data->SetBranchAddress("xtals_timestamp",&xtals_timestamp);

    TCutG *PID_Barrel = new TCutG("PID_Barrel",6);
    PID_Barrel->SetPoint(0,10854.95,2955.739);
    PID_Barrel->SetPoint(1,3077.65,2955.739);
    PID_Barrel->SetPoint(2,2351.979,2263.173);
    PID_Barrel->SetPoint(3,7762.963,2181.694);
    PID_Barrel->SetPoint(4,10791.85,2724.884);
    PID_Barrel->SetPoint(5,10854.95,2955.739);

    //Output root file for histograms
	char outFileName[512];
	sprintf(outFileName,"./root_outputs/Run00%d_Singles.root",r);

	const int MaxMulti = 128;
	double LabAngle[MaxMulti], Energy[MaxMulti], Excite[MaxMulti], ExciteC12[MaxMulti], Range_Barrel[MaxMulti], Range_EndCap[MaxMulti];
	int Mul;

    TFile *outFile = new TFile(outFileName, "recreate");
	TTree *t1 = new TTree("t1","sorted");
	t1->Branch("Mul", &Mul, "Mul/I");
	t1->Branch("LabAngle", &LabAngle, "LabAngle[Mul]/D");
	t1->Branch("Energy", &Energy, "Energy[Mul]/D");
	t1->Branch("Excite", &Excite, "Excite[Mul]/D");
	t1->Branch("ExciteC12", &ExciteC12, "ExciteC12[Mul]/D");
	t1->Branch("Range_Barrel", &Range_Barrel, "Range_Barrel[Mul]/D");
	t1->Branch("Range_EndCap", &Range_EndCap, "Range_EndCap[Mul]/D");

	// Kinematics
	TH2D *hKin = new TH2D("hKin","Energy v Angle",180,0,180,2000,0,20000);
	// Exc v Ang
	TH2D *hExc_Ang = new TH2D("hExc_Ang","Excitation v Angle",180,0,180,2000,-5,15);
	TH2D *hExc_Ang_Cor = new TH2D("hExc_Ang_cor","Excitation v Angle (w/ Target Correction)",180,0,180,2000,-5,15);
	TH2D *hExc_Ring = new TH2D("hExc_Ring","Excitation v Ring",32,0,32,2000,-5,15);
	TH2D *hExcCor_Ring = new TH2D("hExcCor_Ring","Excitation v Ring (w/ Target Correction)",32,0,32,2000,-5,15);
	TH2D *hExc_Ang_C12 = new TH2D("hExc_Ang_C12","Excitation v Angle",180,0,180,2000,-5,15);
	// PID Hists
	TH2D *hPID = new TH2D("hPID","Barrel Particle ID",2000,0,20000,2000,0,20000);
	TH2D *hPID_Q5_A = new TH2D("hPID_Q5_A","QQQ5 A Particle ID",2000,0,20000,2000,0,20000); 
	TH2D *hPID_Q5_B = new TH2D("hPID_Q5_B","QQQ5 B Particle ID",2000,0,20000,2000,0,20000);
	TH2D *hPID2_Q5_A = new TH2D("hPID2_Q5_A","QQQ5 A Particle ID",2000,0,20000,2000,0,20000); 
	TH2D *hPID2_Q5_B = new TH2D("hPID2_Q5_B","QQQ5 B Particle ID",2000,0,20000,2000,0,20000);

	TH2D *hTS = new TH2D("hTS","TimeStamp",2000,0,20000,50,0,50);
	TH2D *hTS_cut = new TH2D("hTS_cut","TimeStamp Cut",2000,0,20000,50,0,50);
	
	//Getting the number of entries to loop through
	unsigned long long int nEntries = data->GetEntries();

    int counter1 = 0;
    int counter2 = 0;
    int counter3 = 0;

	//Looping through each event:
	for ( unsigned long long int i=0; i<nEntries; i++ )	{
    	
        if (i % 10000 == 0)
      		cout << setiosflags(ios::fixed) << "Entry " << i << " of " << nEntries << ", " << 100 * i / nEntries << "% complete" << "\r" << flush; // Event counter

		data->GetEntry(i);

		Excitation = 0; Excitation_Cor = 0; Etotal = 0;	initial_energy = 0; BB10_Range = 0; LabTheta = 0;

		for(int l=0; l<xtalsMul; l++) {
			hTS->Fill(timeStamp/1e8, xtals_crystalNum[l]);
		}

        if ( timeStamp/1e8 > TS_Gate[n] ) { continue; }

		else {

			for(int l=0; l<xtalsMul; l++) {
				hTS_cut->Fill(timeStamp/1e8, xtals_crystalNum[l]);
			}

		    if( SX3Mul == 1 && BB10Mul==1 && QQQ5Mul==0) {

		        for(int j=0; j<SX3Mul; j++) {

					hit_pos = hit_position_3D("SX3", SX3Upstream[j], SX3Det[j], SX3Strip[j], SX3StripPositionCal[j]);
					initial_energy = initial_proton_energy((SX3SectorEnergy[j]/1000.0), proton_distance_through_target(hit_pos)); 

					spherical_polar_coord = hit_position_r_theta_phi("SX3", SX3Upstream[j], SX3Det[j], SX3Strip[j], SX3StripPositionCal[j]);
					LabTheta = spherical_polar_coord.at(1)*r2d;

					if (SX3Upstream[j]==1 && SX3StripRightEnergy[j]>0.0) {

						hKin->Fill(LabTheta, SX3SectorEnergy[j]);

					}

				    if (SX3Det[j]==BB10Det[j] && SX3Upstream[j]==0 && SX3Det[j]!=0 && SX3Det[j]!=5 && SX3Det[j]!=6 && SX3StripRightEnergy[j]>0.0) {
							         
				        Etotal = SX3SectorEnergy[j] + BB10Energy[j];
						initial_energy = initial_proton_energy((Etotal/1000.0), proton_distance_through_target(hit_pos));

                        Excitation = -rel_q_value(LabTheta, Etotal/1000);
                        Excitation_Cor = -rel_q_value(LabTheta, initial_energy);

						Excitation_C12 = -rel_q_value_C12(LabTheta, Etotal/1000);

					    BB10_Range = pow((pow(Etotal,Barrel_Exponent) - pow(SX3SectorEnergy[j],Barrel_Exponent)),1/Barrel_Exponent);         

                        hPID->Fill(Etotal,BB10_Range);

						LabAngle[j] = LabTheta; Energy[j] = Etotal; Excite[j] = Excitation; ExciteC12[j] = Excitation_C12;
						Range_Barrel[j] = BB10_Range; Range_EndCap[j] = 0;

						//if (PID_Barrel->IsInside(Etotal,BB10_Range)) {
						if ( BB10_Range > 2200. && BB10_Range<3000. ) {

							hKin->Fill(LabTheta, Etotal);
							hExc_Ang->Fill(LabTheta,Excitation);
					    	hExc_Ang_Cor->Fill(LabTheta,Excitation_Cor);

							hExc_Ang_C12->Fill(LabTheta,Excitation_C12);
						}
			    	}
		   	 	}
			}
			t1->Fill();
	    
	        // Now just look at QQQ5s
			Excitation = 0; Excitation_Cor = 0; Etotal = 0;	initial_energy = 0; LabTheta_A = 0; LabTheta_B = 0;
		    dEa = 0.; E1a = 0.; E2a = 0; Etotal_A=0; Ethick_A=0; Range_A=0; // reset dE, E1 and E2
            dEb = 0.; E1b = 0.; E2b = 0; Etotal_A=0; Ethick_A=0; Range_B=0; // reset dE, E1 and E2

		    if( SX3Mul == 0 && BB10Mul==0 && QQQ5Mul>0 && QQQ5Mul<4) {

			    for(int j=0; j<QQQ5Mul; j++) {

					if (QQQ5Upstream[j]==1) {

						hit_pos = hit_position_3D("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);     
			    		distance = proton_distance_through_target(hit_pos);
			    		spherical_polar_coord = hit_position_r_theta_phi("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);
			    		LabTheta = spherical_polar_coord.at(1)*r2d;

						hKin->Fill(LabTheta, QQQ5SectorEnergy[j]);

					}

			        if ( QQQ5Upstream[j] == 0 ) {
			    					
			    	    if (QQQ5Det[j]==0 && QQQ5Sector[j]!=0 && QQQ5Sector[j]!=3) {

			    		    hit_pos = hit_position_3D("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);     
			    			distance = proton_distance_through_target(hit_pos);
			    			spherical_polar_coord = hit_position_r_theta_phi("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);
			    			LabTheta = spherical_polar_coord.at(1)*r2d;

			    			dEa = QQQ5SectorEnergy[j];
							Ring_A = QQQ5Ring[j];
    
			    		}
                        
                        if (QQQ5Det[j]==1) {

			    		    hit_pos = hit_position_3D("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);     
			    			distance = proton_distance_through_target(hit_pos);
			    			spherical_polar_coord = hit_position_r_theta_phi("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);
			    			LabTheta = spherical_polar_coord.at(1)*r2d;

			    			dEb = QQQ5SectorEnergy[j];
							Ring_B = QQQ5Ring[j];

			    		}

                        if (QQQ5Det[j]==2) { E1a = QQQ5SectorEnergy[j]; }
			    		if (QQQ5Det[j]==3) { E1b = QQQ5SectorEnergy[j]; }
                        if (QQQ5Det[j]==4) { E2a = QQQ5SectorEnergy[j]; }
        	    		if (QQQ5Det[j]==5) { E2b = QQQ5SectorEnergy[j]; }
                    
			    	}
			    }

                Etotal_A = 1.005*(dEa + E1a + E2a);	// Stack A needs a small gain correction to match Stack B. Perhaps related to the bad sectors on the dE layer.
			    Ethick_A = E1a + E2a;
				Etotal_B = dEb + E1b + E2b;
			    Ethick_B = E1b + E2b;

			    Range_A = pow((pow(Etotal_A,QQQ5_Exponent) - pow(Ethick_A,QQQ5_Exponent)),1/QQQ5_Exponent);
				Range_B = pow((pow(Etotal_B,QQQ5_Exponent) - pow(Ethick_B,QQQ5_Exponent)),1/QQQ5_Exponent);

                hPID_Q5_A->Fill(Etotal_A,Range_A);
				hPID_Q5_B->Fill(Etotal_B,Range_B);

				hPID2_Q5_A->Fill(Etotal_A,dEa);
				hPID2_Q5_B->Fill(Etotal_B,dEb);

                if (LabTheta_A > 25 && Range_A > 2800 && Range_A<4200 && Etotal_A>Range_A) {

					initial_energy = initial_proton_energy((Etotal_A/1000.0), distance); 
                	Excitation = -rel_q_value(LabTheta_A, Etotal_A/1000);
                	Excitation_Cor = -rel_q_value(LabTheta_A, initial_energy);

					Excitation_C12 = -rel_q_value_C12(LabTheta_A, Etotal_A/1000);

					hKin->Fill(LabTheta_A, Etotal_A);
                    hExc_Ang->Fill(LabTheta_A,Excitation);
			        hExc_Ang_Cor->Fill(LabTheta_A,Excitation_Cor);

					hExc_Ring->Fill(Ring_A,Excitation);
					hExcCor_Ring->Fill(Ring_A,Excitation_Cor);

					hExc_Ang_C12->Fill(LabTheta_A,Excitation_C12);

                }

				if (LabTheta_B > 25 && Range_B > 2800 && Range_B<4200 && Etotal_B>Range_B) {

					initial_energy = initial_proton_energy((Etotal_B/1000.0), distance); 
                	Excitation = -rel_q_value(LabTheta_B, Etotal_B/1000);
                	Excitation_Cor = -rel_q_value(LabTheta_B, initial_energy);
					Excitation_C12 = -rel_q_value_C12(LabTheta_A, Etotal_A/1000);

					hKin->Fill(LabTheta_B, Etotal_B);
                    hExc_Ang->Fill(LabTheta_B,Excitation);
			        hExc_Ang_Cor->Fill(LabTheta_B,Excitation_Cor);

					hExc_Ring->Fill(Ring_B,Excitation);
					hExcCor_Ring->Fill(Ring_B,Excitation_Cor);

					hExc_Ang_C12->Fill(LabTheta_B,Excitation_C12);

                }
		    }
		}
	}

	t1->Write();
	
	cout << "\n" << "Finished Run " << r << endl;

	hExc_Ang->Write();
	hExc_Ang_Cor->Write();
	hExc_Ang_C12->Write();

	hPID->Write();
	hPID_Q5_A->Write();
	hPID_Q5_B->Write();
	hPID2_Q5_A->Write();
	hPID2_Q5_B->Write();
	hTS->Write();
	hTS_cut->Write();
	hKin->Write();
	hExc_Ring->Write();
	hExcCor_Ring->Write();

    outFile->Close();
	}
	return;
}