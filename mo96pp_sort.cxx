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

	int firstRun = 59;	// Min = 59
	int lastRun = 93;	// Max = 93
    int nruns = lastRun - firstRun;

	char inputFileName[512];

    Double_t TS_Gate[35] = {5000., 50000., 50000., 50000., 27700., 1570., 50000., 50000., 50000., 24360., 600., 50000., 50000., 50000., 50000., 50000., 1850., 210., 50000., 50000., 50000., 50000., 50000., 50000., 50000., 9200., 50000., 50000., 50000., 25900., 0., 550., 50000., 50000., 50000.};

    //Double_t TS_Gate[6] = {25900., 0., 550., 50000., 50000., 50000.};

    int r;

	for(int n=0; n<=nruns; n++) {

    r = n + firstRun;

    if ( r==79 || r==80 || r == 89 ) {continue;}

	sprintf(inputFileName,"/mnt/d/sorted/Run00%d_combined.root",r);
    //sprintf(inputFileName,"/mnt/sandisk/files/sorted/Run00%d_combined.root",r);

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
	double distance = 0.0;

    float Egam_Cal;

    double Excitation, Excitation_Cor, Excitation_O16, Excitation_Cor_O16, Excitation_C12, Excitation_Cor_C12;
    bool upsteam;

	double Barrel_Exponent = 1.6;
	double QQQ5_Exponent = 1.7;

	double gain_CD2[4] = {0.939202, 1.0, 1.033238, 1.029436};

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
    data->SetBranchAddress("GRETINATimeStamp",&GRETINATimeStamp);

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
	sprintf(outFileName,"./root_outputs/Run00%d_Coinc.root",r);

	int fMulti = 0;
	int fSX3Det[64] = {0};
	float fEnergy[64] = {0};
	float fAngle[64] = {0};
	float fDeltaE_BB10[64] = {0};
	float fDeltaE_QQQ5[64] = {0};
	float fRange_BB10[64] = {0};
	float fRange_QQQ5[64] = {0};
	float fExcitation[64] = {0};
	float fExcitation_Cor[64] = {0};
	float fExcitation_O16[64] = {0};
	float fExcitation_Cor_O16[64] = {0};
	float fExcitation_C12[64] = {0};
	float fExcitation_Cor_C12[64] = {0};
	float fEgamma[64] = {0};

    int fTDC_Silicon = 0;
    int fTDC_GRETINA = 0;
    int fTDC_RF = 0;
    int fTDC_Silicon_Div = 0;
    int fTDC_Silicon_GRETINA = 0;
    int fTDC_Silicon_Delay = 0;
    int fTDC_Silicon_Upstream = 0;
    unsigned long long fTimeStamp[64] = { 0 };
    unsigned long long fTimeStampGRETINA = 0;
    unsigned long long fGRET_timestamp[64] = { 0 };

    TFile *outFile = new TFile(outFileName, "recreate");
	TTree *t1 = new TTree("t1","sorted data");
	t1->Branch("Multi", &fMulti, "Multi/I");
	t1->Branch("Energy", &fEnergy, "Energy[Multi]/F");
	t1->Branch("Angle", &fAngle, "Angle[Multi]/F");
	t1->Branch("DeltaE_BB10", &fDeltaE_BB10, "DeltaE_BB10[Multi]/F");
	t1->Branch("Range_BB10", &fRange_BB10, "Range_BB10[Multi]/F");
	t1->Branch("DeltaE_QQQ5", &fDeltaE_QQQ5, "DeltaE_QQQ5[Multi]/F");
	t1->Branch("Range_QQQ5", &fRange_QQQ5, "Range_QQQ5[Multi]/F");
	t1->Branch("SX3Det", &fSX3Det, "SX3Det[Multi]/I");
	t1->Branch("Excitation", &fExcitation, "Excitation[Multi]/F");
	t1->Branch("Excitation_Cor", &fExcitation_Cor, "Excitation_Cor[Multi]/F");
	t1->Branch("Excitation_O16", &fExcitation_O16, "Excitation_O16[Multi]/F");
	t1->Branch("Excitation_C11", &fExcitation_C12, "Excitation_C12[Multi]/F");
	t1->Branch("Excitation_Cor_O16", &fExcitation_Cor_O16, "Excitation_Cor_O16[Multi]/F");
	t1->Branch("Excitation_Cor_C12", &fExcitation_Cor_C12, "Excitation_Cor_C12[Multi]/F");
	t1->Branch("Egamma", &fEgamma, "Egamma[Multi]/F");

    t1->Branch("tdcSilicon",&fTDC_Silicon,"tdcSilicon/I");
    t1->Branch("tdcGRETINA",&fTDC_GRETINA,"tdcGRETINA/I");
    t1->Branch("tdcRF",&fTDC_RF,"tdcRF/I");
    t1->Branch("tdcSilicon_Div", &fTDC_Silicon_Div, "tdcSilicon_Div/I");
    t1->Branch("tdcSilicon_GRETINA", &fTDC_Silicon_GRETINA, "tdcSilicon_GRETINA/I");
    t1->Branch("tdcSilicon_Delay", &fTDC_Silicon_Delay,"tdcSilicon_Delay/I");
    t1->Branch("tdcSilicon_Upstream", &fTDC_Silicon_Upstream,"tdcSilicon_Upstream/I");
    t1->Branch("timeStamp",&fTimeStamp, "timeStamp[Multi]/L");
    //t1->Branch("timeStamp",&fTimeStampGRETINA);
    
    t1->Branch("GRET_timestamp",&fGRET_timestamp,"GRET_timestamp[Multi]/L");

    TH1F *hTS_Diff = new TH1F("hTS_Diff","ORRUBA-GRETINA TS Difference",2000,-1000,1000);

	TH2D *hExc_Egam = new TH2D("hExc_Egam","Excitation vs Gamma Energy",10000,0,10000,2000,-5,15);
	TH2D *hExc_Egam_O16 = new TH2D("hExc_Egam_O16","Excitation vs Gamma Energy 16O",10000,0,10000,2000,-5,15);
	TH2D *hExc_Egam_C12 = new TH2D("hExc_Egam_C12","Excitation vs Gamma Energy 12C",10000,0,10000,2000,-5,15);
	TH2D *hExc_Egam_Q5 = new TH2D("hExc_Egam_Q5","Excitation vs Gamma Energy",10000,0,10000,2000,-5,15);
	TH2D *hExc_Egam_O16_Q5 = new TH2D("hExc_Egam_O16_Q5","Excitation vs Gamma Energy 16O",10000,0,10000,2000,-5,15);
	TH2D *hExc_Egam_C12_Q5 = new TH2D("hExc_Egam_C12_Q5","Excitation vs Gamma Energy 12C",10000,0,10000,2000,-5,15);

	TH1D *hExc = new TH1D("hExc","Excitation Coincidences",2000,-5,15);
	TH1D *hExc_O16 = new TH1D("hExc_O16","Excitation Coincidences 16O",2000,-5,15);
	TH1D *hExc_C12 = new TH1D("hExc_C12","Excitation Coincidences 12C",2000,-5,15);
	TH1D *hExc_Q5 = new TH1D("hExc_Q5","Excitation Coincidences",2000,-5,15);
	TH1D *hExc_O16_Q5 = new TH1D("hExc_O16_Q5","Excitation Coincidences 16O",2000,-5,15);
	TH1D *hExc_C12_Q5 = new TH1D("hExc_C12_Q5","Excitation Coincidences 12C",2000,-5,15);

	TH1D *hExc_Cor = new TH1D("hExc_cor","Excitation Coincidences (Target Correction)",2000,-5,15);
	TH1D *hExc_Cor_O16 = new TH1D("hExc_Cor_O16","Excitation Coincidences 16O (Target Correction)",2000,-5,15);
	TH1D *hExc_Cor_C12 = new TH1D("hExc_Cor_C12","Excitation Coincidences 12C (Target Correction)",2000,-5,15);
	TH1D *hExc_Cor_Q5 = new TH1D("hExc_Cor_Q5","Excitation Coincidences (Target Correction)",2000,-5,15);
	TH1D *hExc_Cor_O16_Q5 = new TH1D("hExc_Cor_O16_Q5","Excitation Coincidences 16O (Target Correction)",2000,-5,15);
	TH1D *hExc_Cor_C12_Q5 = new TH1D("hExc_Cor_C12_Q5","Excitation Coincidences 12C (Target Correction)",2000,-5,15);

    TH1D *hExc778 = new TH1D("hExc778","Excitation Coincidences",2000,-5,15);
    TH1D *hExc778_Q5 = new TH1D("hExc778_Q5","Excitation Coincidences",2000,-5,15);
    TH1D *hExc778_O16 = new TH1D("hExc778_O16","Excitation Coincidences",2000,-5,15);
    TH1D *hExc778_O16_Q5 = new TH1D("hExc778_O16_Q5","Excitation Coincidences",2000,-5,15);

	TH2D *hKin = new TH2D("hKin","Kinematics",180,0,180,2000,0,20000);
	TH2D *hPID = new TH2D("hPID","Barrel Particle ID",2000,0,20000,2000,0,20000);
    TH2D *hPID_778 = new TH2D("hPID_778","Barrel Particle ID Gated on 778-keV",2000,0,20000,2000,0,20000);
	TH2D *hExc_Ang = new TH2D("hExc_Ang","Excitation v Angle",180,0,180,2000,-5,15);
    TH2D *hExc_Ang_778 = new TH2D("hExc_Ang_778","Excitation v Angle w/ 778-keV gamma",180,0,180,2000,-5,15);
	TH2D *hExc_Ang_Cor = new TH2D("hExc_Ang_cor","Excitation v Angle (w/ Target Correction)",180,0,180,2000,-5,15);
	TH2D *hKin_Q5 = new TH2D("hKin_Q5","Kinematics",180,0,180,2000,0,20000);
	TH2D *hPID_Q5 = new TH2D("hPID_Q5","Barrel Particle ID ",2000,0,20000,2000,0,20000);
    TH2D *hPID_Q5_778 = new TH2D("hPID_Q5_778","QQQ5 Particle ID ",2000,0,20000,2000,0,20000);
	TH2D *hPID2_Q5 = new TH2D("hPID2_Q5","#DeltaE vs Etotal (QQQ5)",2000,0,20000,2000,0,20000);
	TH2D *hExc_Ang_Q5 = new TH2D("hExc_Ang_Q5","Excitation v Angle",180,0,180,2000,-5,15);
	TH2D *hExc_Ang_Cor_Q5 = new TH2D("hExc_Ang_Cor_Q5","Excitation v Angle (w/ Target Correction)",180,0,180,2000,-5,15);

    TH2D *hCrys_v_TS = new TH2D("hCrys_v_TS","Crystal vs TS",10000,0,10000,50,0,50);
    TH2D *hCrys_v_TS_Cut = new TH2D("hCrys_v_TS_Cut","Crystal vs TS",10000,0,10000,50,0,50);

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

		fMulti=0;

        fTDC_Silicon = tdcSilicon;
		fTDC_GRETINA = tdcGRETINA;
        fTDC_RF = tdcRF;
        fTDC_Silicon_Div = tdcSilicon_Div;
        fTDC_Silicon_GRETINA = tdcSilicon_GRETINA;
        fTDC_Silicon_Delay = tdcSilicon_Delay;
        fTDC_Silicon_Upstream = tdcSilicon_Upstream;

        //fTimeStamp = timeStamp;
        fTimeStampGRETINA = GRETINATimeStamp;

        for (int l=0; l<xtalsMul; l++) {

            hTS_Diff->Fill(timeStamp-xtals_timestamp[l]);

            fGRET_timestamp[l] = xtals_timestamp[l];
            fTimeStamp[l] = timeStamp;
            hCrys_v_TS->Fill(xtals_timestamp[l]/1e8,xtals_crystalNum[l]);
            fMulti++;
        }

        fMulti=0;

        if ( timeStamp/1e8 > TS_Gate[n] ) { continue; }

        else {
    
            for (int l=0; l<xtalsMul; l++) {

                hCrys_v_TS_Cut->Fill(xtals_timestamp[l]/1e8,xtals_crystalNum[l]);
            }

		    if( SX3Mul == 1 && BB10Mul==1 && QQQ5Mul==0) {

		        for(int j=0; j<SX3Mul; j++) {

				    if (SX3Det[j]==BB10Det[j] && SX3Upstream[j]==0 && SX3Det[j]!=4 && SX3Det[j]!=5 && SX3Det[j]!=6 && SX3StripRightEnergy[j]>0.0) {
							         
				        Etotal = SX3SectorEnergy[j] + BB10Energy[j];					

					    hit_pos = hit_position_3D("SX3", SX3Upstream[j], SX3Det[j], SX3Strip[j], SX3StripPositionCal[j]);     
					    initial_energy = initial_proton_energy((Etotal/1000.0), proton_distance_through_target(hit_pos)); 

					    spherical_polar_coord = hit_position_r_theta_phi("SX3", SX3Upstream[j], SX3Det[j], SX3Strip[j], SX3StripPositionCal[j]);
					    LabTheta = spherical_polar_coord.at(1)*r2d;

                        Excitation = -rel_q_value(LabTheta, Etotal/1000);
                        Excitation_O16 = -rel_q_value_O16(LabTheta, Etotal/1000);
                        Excitation_C12 = -rel_q_value_C12(LabTheta, Etotal/1000);

                        Excitation_Cor = -rel_q_value(LabTheta, initial_energy);
                        Excitation_Cor_O16 = -rel_q_value_O16(LabTheta, initial_energy);
                        Excitation_Cor_C12 = -rel_q_value_C12(LabTheta, initial_energy);

					    BB10_Range = pow((pow(Etotal,Barrel_Exponent) - pow(SX3SectorEnergy[j],Barrel_Exponent)),1/Barrel_Exponent);         

					    fEnergy[fMulti] = Etotal;
					    fAngle[fMulti] = LabTheta;
					    fDeltaE_BB10[fMulti] = BB10Energy[j];
					    fRange_BB10[fMulti] = BB10_Range;
					    fSX3Det[fMulti] = SX3Det[j] + 1;
					    fExcitation[fMulti] = Excitation;
					    fExcitation_Cor[fMulti] = Excitation_Cor;
					    fExcitation_O16[fMulti] = Excitation_O16;
					    fExcitation_Cor_O16[fMulti] = Excitation_Cor_O16;
					    fExcitation_C12[fMulti] = Excitation_C12;
					    fExcitation_Cor_C12[fMulti] = Excitation_Cor_C12;
					    fRange_QQQ5[fMulti] = 0.0;
                        fDeltaE_QQQ5[fMulti] = 0.0;

                        hPID->Fill(Etotal,BB10_Range);

                        for (int l=0; l<xtalsMul; l++) {

                            if (xtals_cc1[l]>770. && xtals_cc1[l]<786.) {

                                hPID_778->Fill(Etotal,BB10_Range);

                            }
                        }

						if (PID_Barrel->IsInside(Etotal,BB10_Range)) {

							hExc_Ang->Fill(LabTheta,Excitation);
					    	hExc_Ang_Cor->Fill(LabTheta,Excitation_Cor);
						}

                        if (PID_Barrel->IsInside(Etotal,BB10_Range) && LabTheta>65 && LabTheta<87) {

					        hExc->Fill(Excitation);
					        hExc_Cor->Fill(Excitation_Cor);

					        hExc_O16->Fill(Excitation_O16);
					        hExc_Cor_O16->Fill(Excitation_Cor_O16);

					        hExc_C12->Fill(Excitation_C12);
					        hExc_Cor_C12->Fill(Excitation_Cor_C12);

					        hKin->Fill(LabTheta,Etotal);
					    
                            counter1++;

                            for (int l=0; l<xtalsMul; l++) {

                                counter3++;

                                fEgamma[fMulti] = xtals_cc1[l];
                                hExc_Egam->Fill(xtals_cc1[l],Excitation);
				                hExc_Egam_O16->Fill(xtals_cc1[l],Excitation_O16);
					            hExc_Egam_C12->Fill(xtals_cc1[l],Excitation_C12);

                                if (xtals_cc1[l]>770. && xtals_cc1[l]<786.) {

                                    hExc778->Fill(Excitation);
                                    hExc778_O16->Fill(Excitation_O16);
                                    hExc_Ang_778->Fill(LabTheta,Excitation);

                                    counter2++;
                                }
						        fMulti++;
                            }
                        }
		            }
			    }
		    }
	    
	        // Now just look at QQQ5s
		    //dEa = 0.; E1a = 0.; E2a = 0; // reset dE, E1 and E2
            dE = 0.; E1 = 0.; E2 = 0; // reset dE, E1 and E2

		    if( SX3Mul == 0 && BB10Mul==0 && QQQ5Mul>0 && QQQ5Mul<4) {

			    for(int j=0; j<QQQ5Mul; j++) {

			        if ( QQQ5Upstream[j] == 0 ) {
			    					
			    	    if (QQQ5Det[j]==0) {

			    		    hit_pos = hit_position_3D("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);     
			    			distance = proton_distance_through_target(hit_pos);
			    			spherical_polar_coord = hit_position_r_theta_phi("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);
			    			LabTheta = spherical_polar_coord.at(1)*r2d;

			    			dE = QQQ5RingEnergy[j];
    
			    		}
                        
                        if (QQQ5Det[j]==1) {

			    		    hit_pos = hit_position_3D("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);     
			    			distance = proton_distance_through_target(hit_pos);
			    			spherical_polar_coord = hit_position_r_theta_phi("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);
			    			LabTheta = spherical_polar_coord.at(1)*r2d;

			    			dE = QQQ5SectorEnergy[j];

			    		}

                        if (QQQ5Det[j]==2) { E1 = QQQ5SectorEnergy[j]; }
			    		if (QQQ5Det[j]==3) { E1 = QQQ5SectorEnergy[j]; }
                        if (QQQ5Det[j]==4) { E2 = QQQ5SectorEnergy[j]; }
        	    		if (QQQ5Det[j]==5) { E2 = QQQ5SectorEnergy[j]; }
                    
			    	}
                    /*
                    if ( QQQ5Upstream[j] == 1 && QQQ5Mul==1) {

                        hit_pos = hit_position_3D("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);     
			    	    distance = proton_distance_through_target(hit_pos);
			    	    spherical_polar_coord = hit_position_r_theta_phi("QQQ5", QQQ5Upstream[j], QQQ5Det[j], QQQ5Ring[j], QQQ5Sector[j]);
			    		LabTheta = spherical_polar_coord.at(1)*r2d;

                        dE = 0;
                        E1 = QQQ5SectorEnergy[j];
                        E2 = 0;

                    }
                */
			    }
			          
                Etotal = dE + E1 + E2;
			    Ethick = E1 + E2;
			    Range = pow((pow(Etotal,QQQ5_Exponent) - pow(Ethick,QQQ5_Exponent)),1/QQQ5_Exponent);
			    initial_energy = initial_proton_energy((Etotal/1000.0), distance); 
                Excitation = -rel_q_value(LabTheta, Etotal/1000);
                Excitation_Cor = -rel_q_value(LabTheta, initial_energy);
                Excitation_O16 = -rel_q_value_O16(LabTheta, Etotal/1000);
                Excitation_Cor_O16 = -rel_q_value_O16(LabTheta, initial_energy);
                Excitation_C12 = -rel_q_value_C12(LabTheta, Etotal/1000);
                Excitation_Cor_C12 = -rel_q_value_C12(LabTheta, initial_energy);

			    fEnergy[fMulti] = Etotal;						
			    fAngle[fMulti] = LabTheta;
			    fDeltaE_QQQ5[fMulti] = dE;
			    fRange_QQQ5[fMulti] = Range;
						
			    fSX3Det[fMulti] = 0;
			    fExcitation[fMulti] = Excitation;
			    fExcitation_Cor[fMulti] = Excitation_Cor;
			    fExcitation_O16[fMulti] = Excitation_O16;
			    fExcitation_Cor_O16[fMulti] = Excitation_Cor_O16;
			    fExcitation_C12[fMulti] = Excitation_C12;
			    fExcitation_Cor_C12[fMulti] = Excitation_Cor_C12;
			    fRange_BB10[fMulti] = 0.0;

                hPID_Q5->Fill(Etotal,Range);
			    hPID2_Q5->Fill(Etotal,dE);

                for (int l=0; l<xtalsMul; l++) {

                    if (xtals_cc1[l]>770. && xtals_cc1[l]<786.) {

                        hPID_Q5_778->Fill(Etotal,Range);

                    }
                }

                if (Range > 2800 && Range<4200 && E1>0.) {

                    //cout << "Test 2" << endl;

                    //hExc->Fill(Excitation);
                    hExc_Q5->Fill(Excitation);
			        hExc_Cor_Q5->Fill(Excitation_Cor);

                    //hExc_O16->Fill(Excitation_O16);
			        hExc_O16_Q5->Fill(Excitation_O16);
			        hExc_Cor_O16_Q5->Fill(Excitation_Cor_O16);

			        hExc_C12_Q5->Fill(Excitation_C12);
			        hExc_Cor_C12_Q5->Fill(Excitation_Cor_C12);

			        hKin_Q5->Fill(LabTheta,Etotal);
			        hPID_Q5->Fill(Etotal,Range);
			        hPID2_Q5->Fill(Etotal,dE);
                    hExc_Ang->Fill(LabTheta,Excitation);
			        hExc_Ang_Q5->Fill(LabTheta,Excitation);
			        hExc_Ang_Cor_Q5->Fill(LabTheta,Excitation_Cor);

                    for (int l=0; l<xtalsMul; l++) {

			            fEgamma[fMulti] = xtals_cc1[l];
						
                        //hExc_Egam->Fill(xtals_cc1[l],Excitation);
			            hExc_Egam_Q5->Fill(xtals_cc1[l],Excitation);
		                hExc_Egam_O16_Q5->Fill(xtals_cc1[l],Excitation_O16);
			        	hExc_Egam_C12_Q5->Fill(xtals_cc1[l],Excitation_C12);

                        if (xtals_cc1[l]>770. && xtals_cc1[l]<786.) {

                            //hExc778->Fill(Excitation);
                            hExc778_Q5->Fill(Excitation);
                            hExc_Ang_778->Fill(LabTheta,Excitation);
                            //hExc778_O16->Fill(Excitation_O16);
                            hExc778_O16_Q5->Fill(Excitation_O16);
                        
                        }
		
			        	fMulti++;

			        }
                }
                    
		    }
		}

		t1->Fill();

	}
	
	cout << "\n" << "Finished Run " << r << endl;
    cout << "\n" << "Singles Counter " << counter1 << endl;
    cout << "\n" << "Coinc Counter " << counter3 << endl;
    cout << "\n" << "Coinc 778 Counter " << counter2 << endl;

    PID_Barrel->Write();

    hCrys_v_TS->Write();
    hCrys_v_TS_Cut->Write();

    hTS_Diff->Write();

	hExc_Egam->Write();
	hExc_Egam_O16->Write();
	hExc_Egam_C12->Write();

	hExc->Write();

    hExc778->Write();
    hExc778_Q5->Write();
    hExc778_O16->Write();
    hExc778_O16_Q5->Write();

	hExc_Cor->Write();
	hExc_O16->Write();
	hExc_Cor_O16->Write();
	hExc_C12->Write();
	hExc_Cor_C12->Write();

	hKin->Write();
	hPID->Write();
    hPID_778->Write();
	hExc_Ang->Write();
	hExc_Ang_Cor->Write();

    hExc_Ang_778->Write();

	hExc_Egam_Q5->Write();
	hExc_Egam_O16_Q5->Write();
	hExc_Egam_C12_Q5->Write();

	hExc_Q5->Write();
	hExc_Cor_Q5->Write();
	hExc_O16_Q5->Write();
	hExc_Cor_O16_Q5->Write();
	hExc_C12_Q5->Write();
	hExc_Cor_C12_Q5->Write();

	hKin_Q5->Write();
	hPID_Q5->Write();
	hPID2_Q5->Write();
	hExc_Ang_Q5->Write();
	hExc_Ang_Cor_Q5->Write();

	t1->Write();
    outFile->Close();

	}
	return;
}
