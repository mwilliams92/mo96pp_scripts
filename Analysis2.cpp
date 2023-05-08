#define Analysis_cxx
#include "Analysis.h"
#include "Utilities.h"
#include "json/json.h"
#include <string>
#include <algorithm>
#include <iomanip>
#include <fstream>

void Analysis::Loop() {
    std::cout << blue << "GODDESS Mo-96 Analysis..." << reset << std::endl;
    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();

    // read and parse config.json
    /*Json::Value config;
    std::ifstream config_stream("config.json");
    std::cout << red << "Opened" << reset << std::endl;
    ASSERT_WITH_MESSAGE(config_stream.is_open(),
                        "Could not find 'config.json'\n");
    config_stream >> config;
    config_stream.close();
    TString OutputPath = config["OutputPath"].asString();
    TString OutputFilePrefix = config["OutputFilePrefix"].asString();*/
    TString OutputPath = "./root_outputs/";
    TString OutputFilePrefix = "goddess96Mo";
    std::cout << green << "\t Finished reading config.json" << reset << std::endl;

    // move these to config.json
    //Reaction Parameters
	Float_t Ma = 1.007276466; //Mass of beam (proton)
	Float_t Mx = 95.9046795; //Mass of target (96Mo)
	Float_t Mb = 1.007276466; //Mass of Ejectile (proton) 
	Float_t My = 95.9046795; //Mass of Recoil (96Mo)
	Float_t Qgs = 0; //Q-value of Ground State
    Float_t Ta = 21.3;
    //Float_t Ta = 14.94;

    // read barrel PID cut file
    //TFile *barrelPIDCutFile = TFile::Open("/Users/ghimire1/Desktop/96Mo/Analysis/96MoAnalysisCode/graphicalCuts/barrelPIDCut.root");
    //TCutG *barrelCut = static_cast<TCutG*>(barrelPIDCutFile->Get("barrelPIDCut"));
    //barrelPIDCutFile->Close();
    //std::cout << green << "\t Finished reading barrelPIDCut.root" << reset << std::endl;

    // create output file   
    TFile* outputHistFile = new TFile(OutputPath + OutputFilePrefix + "_hist.root", "recreate");
    TFile* outputTreeFile = new TFile(OutputPath + OutputFilePrefix + "_tree.root", "recreate");

    // tree
    Int_t fBB10Det, fSX3Det;
    Float_t fBB10Energy, fSX3Energy, fTotalEnergy, fSX3Angle, fExcitationEnergy, fExcitationEnergyPG,fGammaEnergy;
    TTree* barrelTree = new TTree("barrelTree","Barrel Tree");
    TTree* barrelPGTree = new TTree("barrelPGTree","Barrel PG Tree");
    barrelTree->Branch("BB10Energy", &fBB10Energy, "BB10Energy/F");
    barrelTree->Branch("BB10Det", &fBB10Det, "BB10Det/I");
    barrelTree->Branch("SX3Energy", &fSX3Energy, "SX3Energy/F");
    barrelTree->Branch("SX3Det", &fSX3Det, "SX3Det/I");
    barrelTree->Branch("TotalEnergy", &fTotalEnergy, "TotalEnergy/F");
    barrelTree->Branch("SX3Angle", &fSX3Angle, "SX3Angle/F");
    barrelTree->Branch("ExcitationEnergy", &fExcitationEnergy, "ExcitationEnergy/F");

    barrelPGTree->Branch("ExcitationEnergy_PG", &fExcitationEnergyPG, "ExcitationEnergy_PG/F");
    barrelPGTree->Branch("GammaEnergy", &fGammaEnergy, "GammaEnergy/F");

    // barrel histograms 
    TH1F *hBB10Energy = new TH1F("hBB10Energy","BB10 Energy",1000,0,10000);
    TH2F *hBB10DetSX3Det = new TH2F("hBB10DetSX3Det","BB10 Detector Versus SX3 Detector",15,0,15,15,0,15);
    TH2F *hBB10EnergySX3Energy = new TH2F("hBB10EnergySX3Energy","BB10 Energy Versus SX3 Energy",1000,0,20000,1000,0,10000);
    TH2F *hBB10EnergyTotEnergy = new TH2F("hBB10EnergyTotEnergy","BB10 Energy Versus Total Energy",1000,0,30000,1000,0,10000);
    TH2F *hTotEnergyAngle = new TH2F("hTotEnergyAngle","Total Barrel Energy Versus Angle",90,45,95,1800,0,18000);
    TH1F *hExEnergy = new TH1F("hExEnergy","Excitation Energy",2000,0,30000);    
    TH2F *hExEnergyAngle = new TH2F("hExEnergyAngle","Excitation Energy Versus Angle",90,45,95,2000,0,20000);
    TH2F *hExEnergyGammaEnergy = new TH2F("hExEnergyGammaEnergy","Excitation Energy Versus Gamma Energy",2000,0,20000,2000,0,20000);

    // endcap histograms
    TH2F *hdEresE = new TH2F("hdEresE","QQQ5 PID delEnergy vs residualEnergy",1000,0,20000,1000,0,10000);

    // Analysis variables
    Float_t totalEnergy, range, position, angle, excitationEnergy, delE, resE, e1, e2;

    // Fixed Analysis Variable:
    Float_t Exponent = 1.6;

    // data processing loop
    int prevRunNumber = -1;
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {   
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        // get run number for run by run cuts (if any)
        std::string fileName = fChain->GetCurrentFile()->GetName();
        fileName.erase(fileName.begin(), fileName.end() - 18);
        fileName.erase(fileName.end() - 14, fileName.end());
        TString runNumberStr = fileName;
        int runNumber = std::stoi(fileName);

        if(runNumber != prevRunNumber) {
            std::cout << cyan << " \t\t Processing Run: " << magenta << runNumberStr << reset << std::endl;
            std::cout << yellow << " \t\t\t Data Processing Loop" << reset << std::endl;
            prevRunNumber = runNumber;
        }

        // loop over BB10 multiplicity
        for (Int_t j=0; j<BB10Mul; j++){
            for (Int_t k=0; k<SX3Mul; k++){

                //removing the upstrem stuffs
                if (SX3Upstream[k])continue;

                // making telescopes
                if (SX3Det[k]!=BB10Det[j]) continue;

                // making PID cut in BB10Energy vs SX3Energy
                //if(!barrelCut->IsInside(SX3SectorEnergy[k],BB10Energy[j])) continue;

                // totalE = dE + E
                totalEnergy = SX3SectorEnergy[k] + BB10Energy[j];

                // Range
                range = pow((pow(totalEnergy,Exponent) - pow(SX3SectorEnergy[k],Exponent)),1/Exponent);

                // Exclude events outside PID range
                if(range < 2200. || range > 3000.) continue;

                // sx3 Angle
                position = 75*SX3StripPositionCal[k];
                angle = (atan2(100,position)*(180./M_PI)); 

                // excitation energy calculation
                excitationEnergy = (Qgs) -  (1./My) * ( (My+Mb)*(totalEnergy/1000.) - (My-Ma)*Ta - (2.*sqrt(Ma*Mb*Ta*(totalEnergy/1000.))* cos( angle * (M_PI/180.) ) ) ) ;
                excitationEnergy = 1000 * excitationEnergy;

                // filling histograms
                if(!SX3Upstream[k]) {
                    hBB10Energy->Fill(BB10Energy[j]);
                    hBB10DetSX3Det->Fill(SX3Det[j],BB10Det[k]);
                    hBB10EnergySX3Energy->Fill(SX3SectorEnergy[k],BB10Energy[j]);
                    hBB10EnergyTotEnergy->Fill(totalEnergy,BB10Energy[j]);  
                    hTotEnergyAngle->Fill(angle,totalEnergy);
                    hExEnergy->Fill(excitationEnergy);
                    hExEnergyAngle->Fill(angle,excitationEnergy);

                    // particle-gamma matrix
                    for (Int_t i=0; i<xtalsMul; i++){
                        hExEnergyGammaEnergy->Fill(excitationEnergy,xtals_cc1[i]);
                    }
                }
                // fill barrel tree
                if(!SX3Upstream[k]){
                    fBB10Det = BB10Det[j];
                    fSX3Det = SX3Det[k];
                    fBB10Energy = BB10Energy[j];
                    fSX3Energy = SX3SectorEnergy[k];
                    fTotalEnergy = totalEnergy;
                    fSX3Angle = angle;
                    fExcitationEnergy = excitationEnergy;
                    // particle-gamma matrix
                    for (Int_t i=0; i<xtalsMul; i++){
                        fExcitationEnergyPG = excitationEnergy;
                        fGammaEnergy = xtals_cc1[i];
                        barrelPGTree->Fill();
                    }
                    barrelTree->Fill();
                }          
            }

        }
        
        // loop over QQQ5 multiplicity
        for (Int_t j=0; j<QQQ5Mul; j++){
            if (QQQ5Upstream[j]){
                continue;
            }else{
                if (QQQ5Det[j]==0 || QQQ5Det[j]==1) delE = QQQ5RingEnergy[j];
                if (QQQ5Det[j]==2 || QQQ5Det[j]==3) e1 = QQQ5RingEnergy[j];
                if (QQQ5Det[j]==4 || QQQ5Det[j]==4) e2 = QQQ5RingEnergy[j];
                resE = e1+e2;
            }
        }
        hdEresE->Fill(delE,resE);

    }    
        
    // writing histograms
    outputHistFile->cd(); 
    hBB10Energy->Write();
    hBB10DetSX3Det->Write();
    hBB10EnergySX3Energy->Write();
    hBB10EnergyTotEnergy->Write();
    hTotEnergyAngle->Write();
    hExEnergy->Write();
    hExEnergyAngle->Write();
    hExEnergyGammaEnergy->Write();
    outputHistFile->Close(); 
    std::cout << green << "\t Finished Writing Histograms " << reset << std::endl;

    // writing tree
    outputTreeFile->cd();
    barrelTree->Write();
    barrelPGTree->Write();
    outputTreeFile->Close();
    std::cout << green << "\t Finished Writing 'barrelTree' " << reset << std::endl;
    std::cout << blue << "Finished GODDESS Mo-96 Analysis" << reset << std::endl;

}