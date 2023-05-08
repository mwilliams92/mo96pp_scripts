#define Analysis_cxx
#include "Analysis.h"
#include <string>
#include <algorithm>
#include <iomanip>
#include <fstream>

void Analysis::Loop() {
    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();

    // histograms
    TH1F *htDiff = new TH1F("htDiff","Time Difference",1000,-500,500);
	TH1F *hGamma = new TH1F("hGamma","Gamma Histogram",3000,0,3000);
	TH2F *hGammaGamma = new TH2F("hGammaGamma","Gamma-Gamma Histogram",3000,0,3000,3000,0,3000);
    TH1F *hGammaProjection = new TH1F("hGammaProjection","Gamma Histogram Gated on First Peak of G-G Histogram",3000,0,3000);

	// create output file   
	TFile* outputFile = new TFile("/Users/ghimire1/Desktop/96Mo/Analysis/60Co_efficiency/output/60Co_Efficiency.root", "recreate");

    // create .dat file to save counts
    std::ofstream pfile;
    pfile.open("/Users/ghimire1/Desktop/96Mo/Analysis/60Co_efficiency/output/60Co_Efficiency.dat", std::ofstream::out | std::ofstream::trunc);
    pfile.close();

    // making an array of crystals
    Int_t nCrystals=48;
    std::array<int,48> allCrystals({0});
    for (Int_t x=0; x<nCrystals; x++ ){
        allCrystals[x]=x+1;
    }

    // loop over allCrystals
    for (Int_t x=0; x<nCrystals; x++){ 
        std:: cout << "Crystal Number: " << allCrystals[x] << std::endl;
        htDiff->Reset();
        hGamma->Reset();
        hGammaGamma->Reset();
        hGammaProjection->Reset();

        // loop over all entries
        Long64_t nbytes = 0, nb = 0;
        for (Long64_t jentry=0; jentry<nentries;jentry++) {   
            Long64_t ientry = LoadTree(jentry);
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);   nbytes += nb;

            // loop over Gamma Mult
            for (Int_t k=0; k<xtalsMul; k++) {
                Int_t crystal = allCrystals[x];
                if(xtals_crystalNum[k]==crystal){
                    hGamma->Fill(xtals_cc1[k]);
                    for (Int_t l=0; l<xtalsMul; l++) {
                        if(xtals_crystalNum[l]!=crystal) {
                            Float_t tDiff = fabs(xtals_timestamp[k] - xtals_timestamp[l]);
                            htDiff->Fill(tDiff);

                            // without time cut
                            hGammaGamma->Fill(xtals_cc1[k],xtals_cc1[l]);
                            if (xtals_cc1[k]>=1165. && xtals_cc1[k] <= 1180.) hGammaProjection->Fill(xtals_cc1[l]);
                        }
                    }
                }
            }
    
        }    
        // fitting
        std::cout << "Fitting..." << std::endl;
        TF1* fit_gaus = new TF1("fit_gaus","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2])) + [3]*x+[4]",1000,4000);
        TF1* bkgrd = new TF1("bkgrd","[3]*x+[4]",1000,4000);
        fit_gaus->SetParLimits (0,400,50000);
        fit_gaus->SetParLimits(2,0.8,50);

        // fit 1173 keV gamma line
        Float_t guessX = 1173.; 
        fit_gaus->SetParLimits (1,guessX-20,guessX+20);
        hGamma->Fit("fit_gaus","Q","",guessX-50,guessX+50); 
        bkgrd->FixParameter(3,fit_gaus->GetParameter(3));
        bkgrd->FixParameter(4,fit_gaus->GetParameter(4));
        hGamma->Fit("bkgrd","R+","",guessX-50,guessX+50);
        Float_t counts_1173 = (fit_gaus->Integral(guessX-20,guessX+20)/hGamma->GetBinWidth(0)) - (bkgrd->Integral(guessX-20,guessX+20)/hGamma->GetBinWidth(0));
        Float_t binLower = (hGamma->GetXaxis())->FindBin(guessX-20);
        Float_t binUpper = (hGamma->GetXaxis())->FindBin(guessX+20);
        Float_t counts_1173_new = hGamma->Integral(binLower,binUpper) - (bkgrd->Integral(guessX-20,guessX+20)/hGamma->GetBinWidth(0));
        std::cout << counts_1173 << '\t' << counts_1173_new << std::endl;

        // fit 1333 keV gamma line
        guessX = 1333.; 
        fit_gaus->SetParLimits (1,guessX-20,guessX+20);
        fit_gaus->FixParameter(3,0);
        fit_gaus->FixParameter(4,0);
        hGammaProjection->Fit("fit_gaus","Q","",guessX-50,guessX+50);
        Float_t counts_1333 = fit_gaus->Integral(guessX-20,guessX+20)/hGammaProjection->GetBinWidth(0);
        binLower = (hGammaProjection->GetXaxis())->FindBin(guessX-20);
        binUpper = (hGammaProjection->GetXaxis())->FindBin(guessX+20);
        Float_t counts_1333_new = hGammaProjection->Integral(binLower,binUpper);
        std::cout << counts_1333 << '\t' << counts_1333_new << std::endl;

        // efficiency
        Float_t eff = (counts_1333_new/counts_1173_new * 100)*(nCrystals/(nCrystals-1));
        std::cout << "crystal= " << allCrystals[x] << std::endl;
        std::cout << "Efficiency(counts_1333 / counts_1173 * 100 %) is: " << std::setprecision(4) << eff << " %" << std::endl;

        // writing to .dat file
        pfile.open("/Users/ghimire1/Desktop/96Mo/Analysis/60Co_efficiency/output/60Co_Efficiency.dat", std::ofstream::app);
		pfile << allCrystals[x] << '\t' << counts_1173_new << '\t' << counts_1333_new << std::endl;
		pfile.close();
    } // end of loop over x

    // writing histograms
    outputFile->cd(); 
    htDiff->Write();
    hGamma->Write();
    hGammaGamma->Write();
    hGammaProjection->Write();
    outputFile->Close();  
}