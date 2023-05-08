// g++ -std=c++11 AlphaCalib.cpp `root-config --cflags --libs` -lSpectrum -o calibrate
#include <iostream>
#include <fstream>
#include <vector>
#include "TH1D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TSpectrum.h"
#include "TTree.h"

using namespace std;

vector<double> findPeaksAndFitGaussians(TH1F* hist, double threshold, double sigma) {
    
    hist->GetXaxis()->SetRangeUser(300,4096);
    int nbins = hist->GetNbinsX();
    int npeaks = 0;
    double* peaks;
    TSpectrum *spectrum = new TSpectrum(1000);
    npeaks = spectrum->Search(hist, sigma, "", threshold);
    int counter = 0;
    while (npeaks != 5) {
        if (npeaks > 5) {
            threshold += 0.01;
        } else  {
            threshold -= 0.01;
        }
        if (threshold >= 1.0 || threshold <= 0.0 || counter==50) {
            cerr << "Error: Unable to find appropriate threshold." << endl;
            break;
        }
        npeaks = spectrum->Search(hist, sigma, "", threshold);
        counter++; 
    }
    peaks = spectrum->GetPositionX();
    printf("Found %d peaks \n",npeaks); 
    vector<double> peakpos(npeaks);
    for (int i=0; i<npeaks; i++) {
        peakpos[i] = peaks[i];
    }
    sort(peakpos.begin(), peakpos.end());

    vector<double> peakamp(npeaks);
    for (int i=0; i<npeaks; i++) {
        peakamp[i] = hist->GetBinContent(hist->GetXaxis()->FindBin(peakpos[i]));
    }
    TF1* fitFunc = new TF1("fitFunc","gaus");
    vector<double> fitResults(npeaks*3);
    for (int i=0; i<npeaks; i++) {
        double lowlim = peakpos[i] - 2*sigma*hist->GetBinWidth(hist->GetXaxis()->FindBin(peakpos[i]));
        double highlim = peakpos[i] + 2*sigma*hist->GetBinWidth(hist->GetXaxis()->FindBin(peakpos[i]));
        fitFunc->SetRange(lowlim, highlim);
        hist->Fit(fitFunc,"QR+");
        fitResults[i*3] = fitFunc->GetParameter(0);
        fitResults[i*3+1] = fitFunc->GetParameter(1);
        fitResults[i*3+2] = fitFunc->GetParError(1);
    }
    
    delete fitFunc;
    delete spectrum;
    return fitResults;
}

TH1F* getHistogramFromTree(TFile* file, TString treeName, TString branchName, TString cut="", TString histName="hist") {
    
    TTree* tree = (TTree*) file->Get(treeName);
    if (!tree) {
        std::cout << "Error: Could not find tree " << treeName << " in file " << file->GetName() << std::endl;
        return nullptr;
    }
    TH1F* hist = new TH1F(histName, histName, 2048, 0, 4096);
    hist->GetXaxis()->SetTitle(branchName);
    hist->GetYaxis()->SetTitle("counts");
    tree->Draw(branchName + ">>" + histName, cut);
    return hist;
}

TGraphErrors* getMeansGraph(vector<double> fitresults) {

    const int n = 5;
    double xValues[n];
    double xErrors[n];
    double yValues[5] = {5423., 5685., 6288., 6778, 8784.};;
    double yErrors[5] = {1., 1., 1., 1., 1.};

    for (int i = 0; i < 5; i++) {
        xValues[i] = fitresults[i*3+1];
        //xErrors[i] = fitresults[i*3+2];
        xErrors[i] = 1.0;
    }

    TGraphErrors* graph = new TGraphErrors(n, xValues, yValues, xErrors, yErrors);

    graph->SetTitle("Means of Gaussian Fits");

    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1);
    graph->SetMarkerColor(kBlack);

    graph->GetXaxis()->SetTitle("ADC Value");
    graph->GetYaxis()->SetTitle("Source Energy");

    return graph;
}

TF1* fitMeansGraph(TGraphErrors* graph) {
    
    TF1* fitFunc = new TF1("fitFunc", "[0] + [1]*x");
    fitFunc->SetParameter(0, 0);
    fitFunc->SetParameter(1, 1);
    graph->Fit(fitFunc, "Q");

    return fitFunc;
}

int main(int argc, char* argv[]) {
  // Define the default values for the user options
  double threshold = 0.5;
  double sigma = 2.0;

  // Parse the user options
  int opt;
  while ((opt = getopt(argc, argv, "t:s:")) != -1) {
    switch (opt) {
      case 't':
        threshold = std::atof(optarg);
        break;
      case 's':
        sigma = std::atof(optarg);
        break;
      default: /* '?' */
        std::cerr << "Usage: " << argv[0] << " [-t threshold] [-s sigma] input.root\n";
        return 1;
    }
  }

  // Check that the correct number of non-option arguments were provided
  if (optind != argc - 1) {
    std::cerr << "Usage: " << argv[0] << " [-t threshold] [-s sigma] input.root\n";
    return 1;
  }

  // Open the input file
  TFile* inputFile = TFile::Open(argv[optind], "READ");
  if (!inputFile) {
    std::cerr << "Error: failed to open input file '" << argv[optind] << "'\n";
    return 1;
  }

  TFile *outFile = TFile::Open("calib_output.root","RECREATE");
  outFile->Close();

  TString hName;
  TString gName;
  TString cutString;

  int n=0;

    ofstream txtoutFile("Calibration.dat");
    if (!txtoutFile.is_open()) {
        cout << "Error opening file!" << endl;
        return 1;
    }

  bool SX3 = false;

  if (SX3) {
  // Loop over SX3Det
    for (int i = 5; i < 11; i++) {
        // Loop over SX3Sector
        for (int j = 0; j < 4; j++) {
            // Construct cut string
            cutString = "SX3Upstream == 0 && SX3Det == ";
            cutString += i;
            cutString += " && SX3Sector == ";
            cutString += j;
            // Get histogram with cut
            hName = Form("hist_%d_%d", i, j);
            TH1F* hist = getHistogramFromTree(inputFile, "data", "SX3SectorEnergy", cutString, hName);
            
            // Do something with histogram
            std::cout << "Got histogram for SX3Det = " << i << " and SX3Sector = " << j << std::endl;
            // Find and fit the peaks in the histogram
            std::vector<double> fitresult = findPeaksAndFitGaussians(hist, threshold, sigma);
            
            TFile* outFile = TFile::Open("calib_output.root", "UPDATE");
            hist->Write();

            if (fitresult.size() == 15) { 
                TGraphErrors *Graph = getMeansGraph(fitresult);
                gName = Form("Graph_%d_%d", i, j);
                Graph->SetName(gName);
                TF1* fitFunc = fitMeansGraph(Graph);
                txtoutFile << i << "\t" << j+8 << "\t" << fitFunc->GetParameter(1) << "\t" << fitFunc->GetChisquare() / fitFunc->GetNDF() << "\n";

                Graph->Write();
                delete Graph;
                delete fitFunc;
             }
            else {
                txtoutFile << i << "\t" << j+8 << "\t" << 0 << "\t" << "n/a" << "\n";
            }
            
            
            
            //fitFunc->Write();
            outFile->Write();
            outFile->Close();

            // Clean up memory
            delete outFile;
            delete hist;
        }
    }
  }

  if (!SX3) {
  // Loop over SX3Det
    for (int i = 0; i < 2; i++) {
        // Loop over SX3Sector
        for (int j = 0; j < 4; j++) {
            // Construct cut string
            cutString = "QQQ5Upstream == 0 && QQQ5Det == ";
            cutString += i;
            cutString += " && QQQ5Sector == ";
            cutString += j;
            // Get histogram with cut
            hName = Form("hist_%d_%d", i, j);
            TH1F* hist = getHistogramFromTree(inputFile, "data", "QQQ5SectorEnergy", cutString, hName);
            
            // Do something with histogram
            std::cout << "Got histogram for QQQ5Det = " << i << " and QQQ5Sector = " << j << std::endl;
            // Find and fit the peaks in the histogram
            std::vector<double> fitresult = findPeaksAndFitGaussians(hist, threshold, sigma);
            
            TFile* outFile = TFile::Open("calib_output.root", "UPDATE");
            hist->Write();

            if (fitresult.size() == 15) { 
                TGraphErrors *Graph = getMeansGraph(fitresult);
                gName = Form("Graph_%d_%d", i, j);
                Graph->SetName(gName);
                TF1* fitFunc = fitMeansGraph(Graph);
                txtoutFile << i << "\t" << j << "\t" << fitFunc->GetParameter(0) << "\t" << fitFunc->GetParameter(1) << "\t" << fitFunc->GetChisquare() / fitFunc->GetNDF() << "\n";

                Graph->Write();
                delete Graph;
                delete fitFunc;
             }
            else {
                txtoutFile << i << "\t" << j << "\t" << 0 << "\t" << "n/a" << "\n";
            }
            
            //fitFunc->Write();
            outFile->Write();
            outFile->Close();

            // Clean up memory
            delete outFile;
            delete hist;
        }
    }
  }

  txtoutFile.close();  
 
  inputFile->Close();
  outFile->Close();

  delete inputFile;
  delete outFile;

  return 0;
}