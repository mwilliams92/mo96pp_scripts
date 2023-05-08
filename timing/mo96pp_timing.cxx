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

void mo96pp_timing() {

	//TFile *inputFile = TFile::Open("/mnt/sandisk/files/sorted/Run0060_combined.root");
    TFile *inputFile = TFile::Open("/mnt/d/sorted/Run0060_combined.root");
	TTree *Chain = (TTree*)inputFile->Get("data_combined");

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

	// === TDC ===
	int   tdcSilicon = 0;
	int   tdcGRETINA = 0;
	int	  tdcDiff = 0;
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
	float Excitation = 0.0;
	double Etotal = 0.0;
	double dE = 0.0;
	double E1 = 0.0;
	double E2 = 0.0;
	double Ethick = 0.0;
	double BB10_Range = 0.0;
	double Range = 0.0;

	double Barrel_Exponent = 1.6;
	double QQQ5_Exponent = 1.7;

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


    int fMulti;
    long long fTimeStamp;
    long long fTimeStamp_GT;
    long long fTimeStampDiff[128] = {0};
    float fExcite[128] = {0} ;
    float fExcite_Coinc[128] = {0};
    float fEgamma[128] = {0};

    //Output root file for histograms
    TFile write("Timing_Histograms.root", "recreate");
    TTree *t1 = new TTree("t1","sorted data");
    t1->Branch("Multi",&fMulti, "Multi/I");
    t1->Branch("timeStamp",&fTimeStamp, "timeStamp/L");
    t1->Branch("timeStamp_GT",&fTimeStamp_GT, "timeStamp_GT/L");
    t1->Branch("Excite",&fExcite, "Excite[Multi]/F");
    t1->Branch("Egamma",&fEgamma, "Egamma[Multi]/F");
    t1->Branch("Excite_Coinc",&fExcite_Coinc, "Excite_Coinc[Multi]/F");
    t1->Branch("TimeStampDiff",&fTimeStampDiff, "TimeStampDiff[Multi]/L");

	// Histograms
	TH1D* silicon_TDC_hist = new TH1D("silicon_TDC_hist", "silicon_TDC_hist", 4048, 0, 4048);
    TH1D* gretina_TDC_hist = new TH1D("gretina_TDC_hist", "grertina_TDC_hist", 4048, 0, 4048);
	TH1D* tdcDiff_hist = new TH1D("tdcDiff_hist", "tdcDiff_hist", 4048, 0, 4048);
	TH1D* RF_TDC_hist = new TH1D("RF_TDC_hist", "RF_TDC_hist", 4048, 0, 4048);

	TH2D* delta_timestamp_vs_Run_hist = new TH2D("delta_timestamp_vs_Run_hist", "delta_timestamp_vs_Run_hist", 120, 0, 120, 500, -100, 400);
    TH2D* tdcGRETINA_vs_Run_hist = new TH2D("tdcGRETINA_vs_Run_hist", "tdcGRETINA_vs_Run_hist", 120, 0, 120, 4096, 0, 4096);

    TH1F *hTS_Diff = new TH1F("hTS_Diff","ORRUBA-GRETINA TS Difference",2000,-1000,1000);
    TH1F *hExc = new TH1F("hExc","Excitation Singles",2000,-5,15);
    TH1F *hExc_Coinc = new TH1F("hExc_Coinc","Excitation Coincidences",2000,-5,15);
    TH1F *hExc_Coinc2 = new TH1F("hExc_Coinc2","Excitation Coincidences",2000,-5,15);
    TH1F *hEgam = new TH1F("hEgam","Excitation Coincidences",10000,0,10000);
    TH2F *hExcEgam = new TH2F("hExcEgam","Excitation Coincidences",10000,0,10000,2000,-5,15);
    TH2F *hExcAng = new TH2F("hExcAng","Excitation vs Angle",180,0,180,2000,-5,15);
    TH2F *hExcAng_Coinc = new TH2F("hExcAng_Coinc","Excitation vs Angle",180,0,180,2000,-5,15);
    

	string runNumber_str;
    int runNumber;

    float counter1, counter2;

	unsigned long long int nEntries = Chain->GetEntries();

    long long TS_Diff;

	for ( unsigned long long int i=0; i<nEntries; i++ )	{
    		
		Chain->GetEntry(i);

        fMulti = 0;
        //counter1 = 0;
        //counter2 = 0;

        for (int l=0; l<xtalsMul; l++) {

            TS_Diff = timeStamp-xtals_timestamp[l];

            hTS_Diff->Fill(TS_Diff);

            //fTimeStamp = timeStamp;
            //fTimeStamp_GT = xtals_timestamp[l];

            //fTimeStampDiff[l] = timeStamp-xtals_timestamp[l];

            //cout 

            
            
            //fTimeStamp[l] = timeStamp;
            //fMulti++;
        }
        
        //if (foundGRETINA) {

             if( SX3Mul == 1 && BB10Mul==1) {

                for(int j=0; j<SX3Mul; j++) {
            		
                    if (SX3Det[j]==BB10Det[j] && SX3Upstream[j]==0 && SX3Det[j]!=4 && SX3Det[j]!=5 && SX3Det[j]!=6 && SX3StripRightEnergy[j]>0.0) {

                        Etotal = SX3SectorEnergy[j] + BB10Energy[j];					
                        BB10_Range = pow((pow(Etotal,1.6) - pow(SX3SectorEnergy[j],1.6)),1/1.6);

                        if (BB10_Range > 2200 && BB10_Range<3000) {
    
                            spherical_polar_coord = hit_position_r_theta_phi("SX3", SX3Upstream[j], SX3Det[j], SX3Strip[j], SX3StripPositionCal[j]);
					        LabTheta = spherical_polar_coord.at(1)*r2d;
                            Excitation = -rel_q_value(LabTheta, Etotal/1000);
                            fExcite[fMulti] = Excitation;

                            //cout << "fMulti" << fMulti << endl;

                            hExc->Fill(Excitation);
                            hExcAng->Fill(LabTheta, Excitation);

                            //if (foundGRETINA) {

                                hExc_Coinc->Fill(Excitation);

                                for (int k=0; k<xtalsMul; k++) {

                                    hEgam->Fill(xtals_cc[k]);
                                    hExcEgam->Fill(xtals_cc[k],Excitation);
                                    fExcite_Coinc[fMulti] = Excitation;
                                    fEgamma[fMulti] = xtals_cc[k];
                                    fMulti++;
            
                                    if (xtals_cc[k]>770 && xtals_cc[k]<785) {

                                        hExc_Coinc2->Fill(Excitation);
                                        hExcAng_Coinc->Fill(LabTheta, Excitation);
                                    }

                                }      

                           // }

                        }

                    }

                }
            }

        //}
      
        t1->Fill();

		if (i % 10000 == 0)
      		cout << setiosflags(ios::fixed) << "Entry " << i << " of " << nEntries << ", " << 100 * i / nEntries << "% complete" << "\r" << flush; // Event counter
	}

    cout << "\n" << counter1 << " " << counter2 << endl;

    hEgam->Write();
    hExcAng->Write();
    hExcAng_Coinc->Write();
    hTS_Diff->Write();
    hExc->Write();
    hExc_Coinc->Write();
    hExc_Coinc2->Write();
    hExcEgam->Write();
    t1->Write();

	silicon_TDC_hist->Write();
	gretina_TDC_hist->Write();
	tdcDiff_hist->Write();
	RF_TDC_hist->Write();

	delta_timestamp_vs_Run_hist->Write();
	tdcGRETINA_vs_Run_hist->Write();

	cout << "\n" << "Done!" << endl;

	return;
}
