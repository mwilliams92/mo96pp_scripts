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

void mo96pp_sort_coinc() {

	int firstRun = 60;	// Low Energy Data Min = 59 // High Energy Data Min = 13
	int lastRun = 60;	// Low Energy Data Max = 93 // High Energy Data Max = 52
    int nruns = lastRun - firstRun;

	char inputFileName[512];

    double TS_Gate[35] = {5000., 50000., 50000., 50000., 27700., 1570., 50000., 50000., 50000., 24360., 600., 50000., 
	50000., 50000., 50000., 50000., 1850., 210., 50000., 50000., 50000., 50000., 50000., 50000., 50000., 9200., 50000., 50000., 50000., 25900., 0., 550., 50000., 50000., 50000.};

    int r;

	for(int n=0; n<=nruns; n++) {

    r = n + firstRun;

    if ( (r>13 && r<16) || (r>22 && r<26) || r==30 || (r>31 && r<34) || (r>35 && r<39) || r==40  || r==43 || r==36 || r==79 || r==80 || r == 89 ) {continue;}

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
	Long64_t GRETINATimeStamp = 0;

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
	double LabTheta;
	double Etotal;
	double dE, E1, E2, DeltaE, TotalEnergy, Ethick, Range; 
	double distance;
    double Excitation, Excitation_Cor, Excitation_O16, Excitation_C12;

	double Barrel_Exponent = 1.6;
	double QQQ5_Exponent = 1.6;
	int left, upstream, dE_det, right_ring, left_ring;

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
    data->SetBranchAddress("GRETINA_TimeStamp",&GRETINATimeStamp);

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

	TString outFileName = Form("./root_outputs/Run00%d_sorted.root", r);
	
	TFile *outFile = new TFile(outFileName, "recreate");

	TH2D *hTS = new TH2D("hTS","TimeStamp",5000,0,500000,50,0,50);
	TH2D *hTS_cut = new TH2D("hTS_cut","TimeStamp Cut",5000,0,500000,50,0,50);

	TH2D *hEnergy_v_Ring_Left = new TH2D("hEnergy_v_Ring_Left","Energy vs Ring Left",32,0,32,2000,0,20000);
	TH2D *hEnergy_v_Ring_Right = new TH2D("hEnergy_v_Ring_Right","Energy vs Ring Right",32,0,32,2000,0,20000);

	TH1F *hDeltaTS = new TH1F("hDeltaTS","ORRUBA - GRETINA Timestamp",2000,-1000,1000);
	TH1F *hTDC_Silicon_GRETINA = new TH1F("hTDC_Silicon_GRETINA","TDC Silicon-GRETINA",4097,-1,4096);
	TH2F *hDeltaTS_TDC_COINC = new TH2F("hDeltaTS_TDC_COINC","ORRUBA - GRETINA Timestamp vs TDC SILICON GRETINA",4097,-1,4096,2000,-1000,1000);		
	TH2F *hDeltaTS_TDC_GRETINA = new TH2F("hDeltaTS_TDC_GRETINA","ORRUBA - GRETINA Timestamp vs TDC GRETINA",4097,-1,4096,2000,-1000,1000);
	TH2F *hDeltaTS_TDC_Silicon = new TH2F("hDeltaTS_TDC_Silicon","ORRUBA - GRETINA Timestamp vs TDC Silicon",4097,-1,4096,2000,-1000,1000);	
						

	Int_t fEndCap, fLeft, fUpstream, fE1_Det, fdE_Det;
	Int_t fEndCapPG, fLeftPG, fUpstreamPG, fCrys, fE1_DetPG, fdE_DetPG;
    Float_t fEtot, fThetaLab, fDeltaE, fRange, fExc, fExc_C12, fExc_O16, fE1_Energy; 
	Float_t fEtotPG, fThetaLabPG, fDeltaEPG, fRangePG, fExcPG, fExcPG_C12, fExcPG_O16, fE1_EnergyPG, fEgam;
	Long64_t fTimeStamp;
	Long64_t fdeltaTS;

	TTree* SinglesTree = new TTree("t1","Singles Tree");
    TTree* CoincidenceTree = new TTree("t2","Coicidence Tree");

	SinglesTree->Branch("Etot", &fEtot, "Etot/F");
    SinglesTree->Branch("ThetaLab", &fThetaLab, "ThetaLab/F");
	SinglesTree->Branch("DeltaE", &fDeltaE, "DeltaE/F");
	SinglesTree->Branch("E1_Energy", &fE1_Energy, "E1_Energy/F");
	SinglesTree->Branch("Range", &fRange, "Range/F");
	SinglesTree->Branch("Upstream", &fUpstream, "Upstream/I");
    SinglesTree->Branch("EndCap", &fEndCap, "EndCap/I");
	SinglesTree->Branch("Left", &fLeft, "Left/I");
	SinglesTree->Branch("dE_Det", &fdE_Det, "dE_Det/I");
	SinglesTree->Branch("E1_Det", &fE1_Det, "E1_Det/I");
    SinglesTree->Branch("Exc", &fExc, "Exc/F");
	SinglesTree->Branch("Exc_C12", &fExc_C12, "Exc_C12/F");
	SinglesTree->Branch("Exc_O16", &fExc_O16, "Exc_O16/F");
	SinglesTree->Branch("tdcSilicon_Div", &tdcSilicon_Div,"tdcSilicon_Div/I");
    SinglesTree->Branch("tdcSilicon_GRETINA", &tdcSilicon_GRETINA,"tdcSilicon_GRETINA/I");
	SinglesTree->Branch("deltaTS", &fdeltaTS, "deltaTS/L");
	SinglesTree->Branch("tdcGRETINA",&tdcGRETINA,"tdcGRETINA/I");


	CoincidenceTree->Branch("Etot", &fEtotPG, "EtotPG/F");
    CoincidenceTree->Branch("ThetaLab", &fThetaLabPG, "ThetaLabPG/F");
	CoincidenceTree->Branch("DeltaE", &fDeltaEPG, "DeltaEPG/F");
	CoincidenceTree->Branch("E1_Energy", &fE1_EnergyPG, "E1_EnergyPG/F");
	CoincidenceTree->Branch("Range", &fRangePG, "RangePG/F");
	CoincidenceTree->Branch("Upstream", &fUpstreamPG, "UpstreamPG/I");
    CoincidenceTree->Branch("EndCap", &fEndCapPG, "EndCapPG/I");
	CoincidenceTree->Branch("Left", &fLeftPG, "LeftPG/I");
	CoincidenceTree->Branch("dE_Det", &fdE_DetPG, "dE_DetPG/I");
	CoincidenceTree->Branch("E1_Det", &fE1_DetPG, "E1_DetPG/I");
    CoincidenceTree->Branch("Exc", &fExcPG, "ExcPG/F");
	CoincidenceTree->Branch("Exc_C12", &fExcPG_C12, "ExcPG_C12/F");
	CoincidenceTree->Branch("Exc_O16", &fExcPG_O16, "ExcPG_O16/F");
	CoincidenceTree->Branch("Egam", &fEgam, "Egam/F");
	CoincidenceTree->Branch("TimeStamp", &fTimeStamp, "TimeStamp/L");
	CoincidenceTree->Branch("Crys", &fCrys, "Crys/I");
	CoincidenceTree->Branch("deltaTS", &fdeltaTS, "deltaTS/L");
	CoincidenceTree->Branch("tdcSilicon_Div", &tdcSilicon_Div,"tdcSilicon_Div/I");
    CoincidenceTree->Branch("tdcSilicon_GRETINA", &tdcSilicon_GRETINA,"tdcSilicon_GRETINA/I");
	CoincidenceTree->Branch("tdcGRETINA",&tdcGRETINA,"tdcGRETINA/I");

	double correction_factor = 1.0167;
	double correction_factor2 = 1.008;
	//double correction_factor_US_SX3 = 1.0262;

	//Getting the number of entries to loop through
	unsigned long long int nEntries = data->GetEntries();

	//Looping through each event:
	for ( unsigned long long int i=0; i<nEntries; i++ )	{
    	
        if (i % 10000 == 0)
      		cout << setiosflags(ios::fixed) << "Entry " << i << " of " << nEntries << ", " << 100 * i / nEntries << "% complete" << "\r" << flush; // Event counter

		data->GetEntry(i);
		
		for(int l=0; l<xtalsMul; l++) {
			hTS->Fill(timeStamp/1e8, xtals_crystalNum[l]);
		}
		
        //if ( timeStamp/1e8 > TS_Gate[n] ) { continue; }
			
		for(int l=0; l<xtalsMul; l++) {
			hTS_cut->Fill(timeStamp/1e8, xtals_crystalNum[l]);

		}

		if (SX3Mul == 1 && QQQ5Mul == 0) {

			for(int j=0; j<SX3Mul; j++) {

				Excitation = 0; Excitation_Cor = 0; Excitation_O16 = 0; Excitation_C12 = 0;
				Etotal = 0;	initial_energy = 0; Range = 0; LabTheta = 0; DeltaE = 0; dE_det = -1;

				hit_pos = hit_position_3D("SX3", SX3Upstream[j], SX3Det[j], SX3Strip[j], SX3StripPositionCal[j]);
				spherical_polar_coord = hit_position_r_theta_phi("SX3", SX3Upstream[j], SX3Det[j], SX3Strip[j], SX3StripPositionCal[j]);
				LabTheta = spherical_polar_coord.at(1)*r2d;

				if (SX3Upstream[j] == 1 && SX3StripRightEnergy[j]>0) {

					Etotal = SX3SectorEnergy[j] * correction_factor;
					DeltaE = 0.;
					dE_det = -1;
					Range = 0.;	     
					initial_energy = initial_proton_energy((Etotal/1000.0), proton_distance_through_target(hit_pos)); 

					Excitation = -rel_q_value(LabTheta, Etotal/1000);
               	 	Excitation_Cor = -rel_q_value(LabTheta, initial_energy);

					Excitation_O16 = -rel_q_value_O16(LabTheta, Etotal/1000);
					Excitation_C12 = -rel_q_value_C12(LabTheta, Etotal/1000);

				} 

				if (SX3Upstream[j] == 0 && BB10Mul == 1 && SX3Det[j]==BB10Det[j] && SX3Det[j]!=0 && SX3Det[j]!=5 && SX3Det[j]!=6 && SX3StripRightEnergy[j]>0) {
							         
				    Etotal = SX3SectorEnergy[j] + BB10Energy[j];
					DeltaE = BB10Energy[j];
					dE_det = BB10Det[j];
					initial_energy = initial_proton_energy((Etotal/1000.0), proton_distance_through_target(hit_pos)); 

					Excitation = -rel_q_value(LabTheta, Etotal/1000);
                    Excitation_Cor = -rel_q_value(LabTheta, initial_energy);

					Excitation_O16 = -rel_q_value_O16(LabTheta, Etotal/1000);
					Excitation_C12 = -rel_q_value_C12(LabTheta, Etotal/1000);				

					Range = pow((pow(Etotal,Barrel_Exponent) - pow(SX3SectorEnergy[j],Barrel_Exponent)),1/Barrel_Exponent);
				}
                    	
				fEtot = Etotal;
            	fThetaLab = LabTheta;
				fDeltaE = DeltaE;
				fE1_Energy = SX3SectorEnergy[j];
				fE1_Det = SX3Det[j];
				fdE_Det = dE_det;
				fRange = Range;
				fUpstream = SX3Upstream[j];
				fEndCap = 0;
				fLeft = -1;
            	fExc = Excitation;
				fExc_C12 = Excitation_C12;
				fExc_O16 = Excitation_O16;
				fdeltaTS = timeStamp - GRETINATimeStamp;
				SinglesTree->Fill();

            	// particle-gamma tree
            	for (int l=0; l<xtalsMul; l++){
					//if (xtals_crystalNum[l] != 2) {
						fEtotPG = Etotal;
						fThetaLabPG = LabTheta;
						fDeltaEPG = DeltaE;
						fE1_EnergyPG = SX3SectorEnergy[j];
						fE1_DetPG = SX3Det[j];
						fdE_DetPG = dE_det;
						fRangePG = Range;
						fUpstreamPG = SX3Upstream[j];
						fEndCapPG = 0;
						fLeftPG = -1;
                    	fExcPG = Excitation;
						fExcPG_C12 = Excitation_C12;
						fExcPG_O16 = Excitation_O16;
                    	fEgam = xtals_cc1[l];
						fCrys = xtals_crystalNum[l];
						fTimeStamp = xtals_timestamp[l];
						fdeltaTS = timeStamp - GRETINATimeStamp;
                    	CoincidenceTree->Fill();

						hDeltaTS->Fill(fdeltaTS);
						hTDC_Silicon_GRETINA->Fill(tdcSilicon_GRETINA);
						hDeltaTS_TDC_COINC->Fill(tdcSilicon_GRETINA,fdeltaTS);
						hDeltaTS_TDC_GRETINA->Fill(tdcGRETINA,fdeltaTS);
						hDeltaTS_TDC_Silicon->Fill(tdcSilicon,fdeltaTS);
						
					//}
            	}
			}
		}	

		if( QQQ5Mul>1 && QQQ5Mul<4 && SX3Mul == 0) {

			dE = 0.; E1 = 0.; E2 = 0.; Ethick = 0.; TotalEnergy = 0.; left = -1;

			for(int j=0; j<QQQ5Mul; j++) {
				
				if ( QQQ5Upstream[j] == 1 ) {

					hit_pos = hit_position_3D("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);     
			    	distance = proton_distance_through_target(hit_pos);
			    	spherical_polar_coord = hit_position_r_theta_phi("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);
			    	LabTheta = spherical_polar_coord.at(1)*r2d;

					dE = QQQ5SectorEnergy[j];
					E1 = 0.; E2 = 0.; left = -1;
					upstream = 1;

				} 
				else {

					upstream = 0;

					if (QQQ5Det[j] < 2 ) {

			    		hit_pos = hit_position_3D("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);     
			    		distance = proton_distance_through_target(hit_pos);
			    		spherical_polar_coord = hit_position_r_theta_phi("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);
			    		LabTheta = spherical_polar_coord.at(1)*r2d;

						if (QQQ5Det[j] == 0) { dE = QQQ5RingEnergy[j]; left = 0; right_ring = QQQ5Ring[j]; } // 2 bad sectors
						if (QQQ5Det[j] == 1) { dE = QQQ5RingEnergy[j]; left = 1; left_ring = QQQ5Ring[j]; } 
    
			    	}

					if ( QQQ5Det[j] == 2 ) E1 = QQQ5RingEnergy[j]; 
					if ( QQQ5Det[j] == 3 ) E1 = QQQ5RingEnergy[j];
					if ( QQQ5Det[j] == 4 ) E2 = QQQ5RingEnergy[j]; // only instrumented sectors
					if ( QQQ5Det[j] == 5 ) E2 = QQQ5RingEnergy[j];
				}
			}

			TotalEnergy = (dE + E1 + E2) * correction_factor;
			if (left == 0) TotalEnergy = TotalEnergy;
			Excitation = -rel_q_value(LabTheta, TotalEnergy/1000);
			Excitation_O16 = -rel_q_value_O16(LabTheta, TotalEnergy/1000);
			Excitation_C12 = -rel_q_value_C12(LabTheta, TotalEnergy/1000);

			if (!upstream) {
				Ethick = E1 + E2;
				Range = pow((pow(TotalEnergy,QQQ5_Exponent) - pow(Ethick,QQQ5_Exponent)),1/QQQ5_Exponent);

				if (left == 0) { hEnergy_v_Ring_Right->Fill(right_ring,TotalEnergy); }
				if (left == 1) { hEnergy_v_Ring_Left->Fill(left_ring,TotalEnergy); }
			}
			else { 
				Ethick = 0.; 
				Range = 0.; 
			}

			if ( (upstream == 0 && left==1 && left_ring>6 && left_ring <28 && left_ring != 10 && left_ring != 11) || (upstream == 0 && right_ring>6 && right_ring <26 && right_ring != 8) ) {

				fEtot = TotalEnergy;
				fThetaLab = LabTheta;
				fDeltaE = dE;
				fRange = Range;
				fUpstream = upstream;
				fEndCap = 1;
				fLeft = left;
				fExc = Excitation;
				fExc_C12 = Excitation_C12;
				fExc_O16 = Excitation_O16;
				SinglesTree->Fill();

				for (int l=0; l<xtalsMul; l++){
					if (xtals_crystalNum[l] != 2) {
						fEtotPG = TotalEnergy;
						fThetaLabPG = LabTheta;
						fDeltaEPG = dE;
						fRangePG = Range;
						fUpstreamPG = upstream;
						fEndCapPG = 1;
						fLeftPG = left;
                    	fExcPG = Excitation;
						fExcPG_C12 = Excitation_C12;
						fExcPG_O16 = Excitation_O16;
                    	fEgam = xtals_cc1[l];
						fCrys = xtals_crystalNum[l];
						fTimeStamp = xtals_timestamp[l];
                    	CoincidenceTree->Fill();
					}
				}
			}
		}
	}
	
	cout << "\n" << "Finished Run " << r << endl;

	hTS->Write();
	hTS_cut->Write();

	hEnergy_v_Ring_Right->Write();
	hEnergy_v_Ring_Left->Write();

	hDeltaTS->Write();
	hTDC_Silicon_GRETINA->Write();
	hDeltaTS_TDC_COINC->Write();
	hDeltaTS_TDC_GRETINA->Write();
	hDeltaTS_TDC_Silicon->Write();

	SinglesTree->Write();
    CoincidenceTree->Write();
	
	outFile->Close();

	}
	return;
}