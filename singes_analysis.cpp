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
#include "analysis.h"

using namespace std;

void processInputFile(const string& inputFileName) {

    // Open the input file
    inputFile = TFile::Open(inputFileName.c_str());
    if (!inputFile) {
        cout << "Could not open input file " << inputFileName << "!" << endl;
        return;
    }

    // Get the tree from the input file
    inputTree = (TTree*)inputFile->Get("data_combined");
    if (!inputTree) {
        cout << "Could not get input tree from file " << inputFileName << "!" << endl;
        inputFile->Close();
        return;
    }

    // Create the output root file
    string outputFileName = "output_" + inputFileName;
    outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    if (!outputFile) {
        cout << "Could not create output file for input file " << inputFileName << "!" << endl;
        inputFile->Close();
        return;
    }

    // Create the output tree
    outputTree = new TTree("output_data", "Output Data Tree");

	// === Analysis parameters ===
	double r2d = 180./ TMath::Pi();
    vector<double> hit_pos;
	vector<double> spherical_polar_coord;
    double initial_energy;
	
    float Theta = 0.0;
	float Exc = 0.0;
	float Etot = 0.0;
	float dE = 0.0;
	float E1 = 0.0;
	float E2 = 0.0;
	float Ethick = 0.0;
    float Range = 0.0;
	double Barrel_Exponent = 1.6;
	double QQQ5_Exponent = 1.7;

    //============================================================
    //   Allocating the branch addresses of the "raw" variables
    //============================================================
    // ================ BB10 Branch Addresses ===========================
    inputTree->SetBranchAddress("BB10Mul",&BB10Mul);
    inputTree->SetBranchAddress("BB10Det",&BB10Det);
    inputTree->SetBranchAddress("BB10Strip",&BB10Strip);
    inputTree->SetBranchAddress("BB10Energy",&BB10Energy);
    // ================ QQQ5 Branch Addresses ===========================
    inputTree->SetBranchAddress("QQQ5Mul",&QQQ5Mul);
    inputTree->SetBranchAddress("QQQ5Upstream",&QQQ5Upstream);
    inputTree->SetBranchAddress("QQQ5Det",&QQQ5Det);
    inputTree->SetBranchAddress("QQQ5Ring",&QQQ5Ring);
    inputTree->SetBranchAddress("QQQ5Sector",&QQQ5Sector);
    inputTree->SetBranchAddress("QQQ5RingEnergy",&QQQ5RingEnergy);
    inputTree->SetBranchAddress("QQQ5SectorEnergy",&QQQ5SectorEnergy);
    // =================== SX3 Branch Address ===========================
    inputTree->SetBranchAddress("SX3Mul",&SX3Mul);
    inputTree->SetBranchAddress("SX3Upstream",&SX3Upstream);
    inputTree->SetBranchAddress("SX3Det",&SX3Det);
    inputTree->SetBranchAddress("SX3Sector",&SX3Sector);
    inputTree->SetBranchAddress("SX3SectorEnergy",&SX3SectorEnergy);
    inputTree->SetBranchAddress("SX3Strip",&SX3Strip);
	inputTree->SetBranchAddress("SX3StripLeftEnergy",&SX3StripLeftEnergy);
    inputTree->SetBranchAddress("SX3StripRightEnergy",&SX3StripRightEnergy);
    inputTree->SetBranchAddress("SX3StripEnergy",&SX3StripEnergy);
	inputTree->SetBranchAddress("SX3StripPositionCal",&SX3StripPositionCal);
	// =================== TDC Branch Address ===========================
    inputTree->SetBranchAddress("tdcSilicon",&tdcSilicon);
    inputTree->SetBranchAddress("tdcGRETINA",&tdcGRETINA);
    inputTree->SetBranchAddress("tdcRF",&tdcRF);
    inputTree->SetBranchAddress("tdcSilicon_Div", &tdcSilicon_Div);
    inputTree->SetBranchAddress("tdcSilicon_GRETINA", &tdcSilicon_GRETINA);
    inputTree->SetBranchAddress("tdcSilicon_Delay", &tdcSilicon_Delay);
    inputTree->SetBranchAddress("tdcSilicon_Upstream", &tdcSilicon_Upstream);
    // =================== Timestamp Branch Address =====================
    inputTree->SetBranchAddress("timeStamp",&timeStamp);
    inputTree->SetBranchAddress("GRETINATimeStamp",&GRETINATimeStamp);
	// ================= GRETINA Branch Address ===================
    inputTree->SetBranchAddress("foundGRETINA",&foundGRETINA);
    inputTree->SetBranchAddress("xtalsMul",&xtalsMul);
    inputTree->SetBranchAddress("xtals_cc",&xtals_cc);
    inputTree->SetBranchAddress("xtals_cc1",&xtals_cc1);
    inputTree->SetBranchAddress("xtals_crystalNum",&xtals_crystalNum);
    inputTree->SetBranchAddress("xtals_timestamp",&xtals_timestamp);

    //============================================================
    //   Assigning New Branches to output ROOT file
    //============================================================
    outputTree->Branch("Multi", &Multi, "Multi/I");
    outputTree->Branch("Excitation", &Excitation, "Excitation[NMAX]/F");
    outputTree->Branch("Excitation_Cor", &Excitation_Cor, "Excitation_Cor[NMAX]/F");
    outputTree->Branch("Etotal", &Etotal, "Etotal[NMAX]/F");
    outputTree->Branch("Labtheta", &LabTheta, "LabTheta[NMAX]/F");
    outputTree->Branch("DeltaE_Barrel", &DeltaE_Barrel, "DeltaE_Barrel[NMAX]/F");
    outputTree->Branch("DeltaE_Endcap", &DeltaE_Endcap, "DeltaE_Endcap[NMAX]/F");
    outputTree->Branch("Ethick_Barrel", &Ethick_Barrel, "Ethick_Barrel[NMAX]/F");
    outputTree->Branch("Ethick_Endcap", &Ethick_Endcap, "Ethick_Endcap[NMAX]/F");
    outputTree->Branch("Range_Barrel", &Range_Barrel, "Range_Barrel[NMAX]/F");
    outputTree->Branch("Range_Endcap", &Range_Endcap, "Range_Endcap[NMAX]/F");
    // Timing Branches to Copy
    outputTree->Branch("tdcSilicon",&tdcSilicon, "tdcSilicon/I");
    outputTree->Branch("tdcGRETINA",&tdcGRETINA, "tdcGRETINA/I");
    outputTree->Branch("tdcRF",&tdcRF, "tdcRF/I");
    outputTree->Branch("tdcSilicon_Div", &tdcSilicon_Div, "tdcSilicon_Div/I");
    outputTree->Branch("tdcSilicon_GRETINA", &tdcSilicon_GRETINA, "tdcSilicon_GRETINA/I");
    outputTree->Branch("tdcSilicon_Delay", &tdcSilicon_Delay, "tdcSilicon_Delay/I");
    outputTree->Branch("tdcSilicon_Upstream", &tdcSilicon_Upstream, "tdcSilicon_Upstream/I");
    // Timestamp Branches to Copy:
    outputTree->Branch("timeStamp",&timeStamp, "timeStamp[NMAX]/L");
    outputTree->Branch("GRETINATimeStamp",&GRETINATimeStamp, "GRETINATimeStamp[NMAX]/L");
    outputTree->Branch("xtals_timestamp",&xtals_timestamp, "xtals_timestamp[NMAX]/L");
    // GRETINA Branches to copy:
    outputTree->Branch("foundGRETINA",&foundGRETINA, "foundGRETINA[NMAX]/B");
    outputTree->Branch("xtals_cc", &xtals_cc, "xtals_cc[NMAX]/F");
    outputTree->Branch("xtals_cc1", &xtals_cc1, "xtals_cc1[NMAX]/F");
    outputTree->Branch("Egam", &xtals_cc1, "Egam[NMAX]/F"); // testing ability to rename a copied value
    outputTree->Branch("xtals_crystalNum",&xtals_crystalNum, "xtals_crystalNum[NMAX]/I");

    unsigned long long int nEntries = inputTree->GetEntries();
    // Loop over events in the input tree
    for (unsigned long long int i = 0; i < nEntries; i++) {

        // Get the event from the input tree
        inputTree->GetEntry(i);

        // Select events with a hit in an SX3 and BB10
        if( SX3Mul == 1 && BB10Mul==1 ) {

            // Loop over SX3 Multiplicity (should just be 1)
            for(int j=0; j<SX3Mul; j++) {
                
                // Select events in telescope detectors and where the hit in the SX3 and BB10 are in the same detector
			    if (SX3Det[j]==BB10Det[j] && SX3Upstream[j]==0 && SX3Det[j]!=0 && SX3Det[j]!=5 && SX3Det[j]!=6) {

                    DeltaE_Barrel[j] = BB10Energy[j];
                    Ethick_Barrel[j] = SX3SectorEnergy[j];
                    Etotal[j] = SX3SectorEnergy[j] + BB10Energy[j];					

				    hit_pos = hit_position_3D("SX3", SX3Upstream[j], SX3Det[j], SX3Strip[j], SX3StripPositionCal[j]);     
					initial_energy = initial_proton_energy((Etot/1000.0), proton_distance_through_target(hit_pos)); 

				    spherical_polar_coord = hit_position_r_theta_phi("SX3", SX3Upstream[j], SX3Det[j], SX3Strip[j], SX3StripPositionCal[j]);
				    LabTheta[j] = spherical_polar_coord.at(1)*r2d;

                    Range_Barrel[j] = pow((pow(Etot,Barrel_Exponent) - pow(SX3SectorEnergy[j],Barrel_Exponent)),1/Barrel_Exponent);

                    Excitation[j] = -rel_q_value(LabTheta[j], Etotal[j]/1000);

                }
            }        
        }
        outputTree->Fill();
    }    

    // Write the output tree to the output file and close the files
    outputFile->cd();
    outputTree->Write();
    outputFile->Close();
    inputFile->Close();
}

int main(int argc, char* argv[]) {

    if (argc == 2) {
        // Process single input file
        processInputFile(argv[1]);
    } else if (argc == 3 && std::string(argv[1]) == "-f") {
        // Process list of input files
        std::ifstream inputListFile(argv[2]);
        std::vector<std::string> inputFiles;
        std::string line;
        while (std::getline(inputListFile, line)) {
            inputFiles.push_back(line);
        }
        inputListFile.close();

        for (const auto& inputFile : inputFiles) {
            processInputFile(inputFile);
        }
    } else {
        // Print usage message
        std::cerr << "Usage: " << argv[0] << " [-f input_file_list] [input_file]" << std::endl;
        std::cerr << "       -f: process input files from list" << std::endl;
        std::cerr << "       input_file: process single input file" << std::endl;
        return 1;
    }

    return 0

}