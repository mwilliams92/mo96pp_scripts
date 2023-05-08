// g++ -std=c++11 fitpedestals.cpp `root-config --cflags --libs` -o fitpedestals
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"

using namespace std;

void findGaussianPeak(TH1F* h, double &mean, double &sigma) {
  
  h->GetXaxis()->SetRangeUser(25,200);

  // Find the bin with the maximum content
  int maxBin = h->GetMaximumBin();

  // Fit the range around the maximum bin to a Gaussian
  TF1* fitFunc = new TF1("fitFunc", "gaus",1,200);
  fitFunc->SetParameter(0, h->GetBinContent(maxBin));
  fitFunc->SetParameter(1, h->GetBinCenter(maxBin));
  fitFunc->SetParameter(2, h->GetRMS());
  h->Fit(fitFunc,"QR0");

  // Get the fit parameters
  mean = fitFunc->GetParameter(1);
  sigma = fitFunc->GetParameter(2);

  // Delete the fit function
  delete fitFunc;
}

double analyseHistogram(const char* inputFileName, int Channel) {

  // Open the input file
  TFile* inFile = TFile::Open(inputFileName);

  // Open the output file or create a new one if it doesn't exist
  TFile* outFile = TFile::Open("output.root", "UPDATE");

  // Get the original histogram from the input file
  TH1F* hist = (TH1F*)inFile->Get(Form("d%d", Channel));

  // Fit the histogram to a Gaussian
  Double_t mean, sigma;
  findGaussianPeak(hist, mean, sigma);

  // Create a fit function from the Gaussian fit results
  TF1* fitFunc = new TF1(Form("fitFunc%d", Channel), "gaus", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
  fitFunc->SetParameter(0, hist->GetMaximum());
  fitFunc->SetParameter(1, mean);
  fitFunc->SetParameter(2, sigma);

  // Write the original histogram and the fit function to the output file
  outFile->cd();
  hist->Write(Form("d%d", Channel));
  fitFunc->Write(Form("fitFunc_%d", Channel));

  // Clean up the original histogram and fit function
  delete hist;
  delete fitFunc;

  // Close the input and output files
  inFile->Close();
  outFile->Close();
  
  // Delete the TTree and the input file
  delete inFile;
  delete outFile;
  
  return mean;
}

// Define a struct to hold a row of data
struct DataRow {
    int col1;
    int col2;
    double col3;
};

// Define a comparison function to sort the rows by the first and second columns
bool compareRows(const DataRow& row1, const DataRow& row2) {
    if (row1.col1 < row2.col1) {
        return true;
    } else if (row1.col1 > row2.col1) {
        return false;
    } else {
        return row1.col2 < row2.col2;
    }
}

void orderfile(string input_txt) {
    // Open the input file
    ifstream infile(input_txt);
    if (!infile.is_open()) {
        cerr << "Error opening input file." << endl;
    }
    
    // Read the data from the input file into a vector of DataRow objects
    vector<DataRow> data;
    string line;
    while (getline(infile, line)) {
        istringstream iss(line);
        DataRow row;
        if (iss >> row.col1 >> row.col2 >> row.col3) {
            data.push_back(row);
        } else {
            cerr << "Error parsing input line: " << line << endl;
        }
    }
    
    // Sort the rows by the first and second columns
    sort(data.begin(), data.end(), compareRows);
    
    infile.close();

    // Open the output file
    ofstream outfile(input_txt);
    if (!outfile.is_open()) {
        cerr << "Error opening output file." << endl;
    }
    
    // Write the sorted data to the output file
    for (const auto& row : data) {
        outfile << row.col1 << "\t" << row.col2 << "\t" << row.col3 << endl;
    }

    outfile.close();
}


int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Error: expected input ROOT file name as argument.\n";
        return 1;
    }

    const char* inputFileName = argv[1];
    TFile* inputFile = TFile::Open(inputFileName);

    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error opening file " << inputFileName << "\n";
        return 1;
    }

    const char* outputFileName = "output.root";

    if (access(outputFileName, F_OK) != -1) {
        remove(outputFileName);
    }

    TFile* outputFile = TFile::Open(outputFileName, "RECREATE");

    if (!outputFile) {
        std::cerr << "Error creating output file " << outputFileName << "\n";
        return 1;
    }

    TTree *tree = new TTree("fitResults", "Fit Results");
    tree->Write();

    outputFile->Close();
  
    ofstream uSX3_outfile("SX3uPedestals.dat");
    
    // Check if the file was successfully opened
    if (!uSX3_outfile.is_open()) {
        cout << "Error opening file!" << endl;
        return 1;
    }

    int Chan, nDets, nStrips;
    double pedestal;

    // Upstream SX3 Front Strips
    // Detectors 0 -> 11
    Chan = 193;
    nDets = 12;
    nStrips = 8;

    for (int i=0; i<nDets; i++) {

        for (int j=0; j<nStrips; j++) {

            pedestal = analyseHistogram(argv[1], Chan);
            uSX3_outfile << i << "\t" << j << "\t" << pedestal << endl;

            Chan++;
        }
    }
  
    // Upstream SX3 Back Strips
    // Detectors 0 -> 3  
    Chan = 145;
    nDets = 4;
    nStrips = 4;

    for (int i=0; i<nDets; i++) {

        for (int j=8; j<nStrips+8; j++) {

            pedestal = analyseHistogram(argv[1], Chan);

            uSX3_outfile << i << "\t" << j << "\t" << pedestal << endl;

            Chan++;

        }
    }

    // Upstream SX3 Back Strips
    // Detectors 4 -> 5 (reversed)
    Chan = 177;
    nDets = 2;
    nStrips = 4;

    for (int i=4; i<6; i++) {

        for (int j=11; j>=8; --j) {

            pedestal = analyseHistogram(argv[1], Chan);
            
            uSX3_outfile << i << "\t" << j << "\t" << pedestal << endl;

            Chan++;

        }
    }

    // Upstream SX3 Back Strips
    // Detectors 6 -> 9
    Chan = 161;
    nDets = 4;
    nStrips = 4;

    for (int i=6; i<10; i++) {

        for (int j=8; j<nStrips+8; j++) {

            pedestal = analyseHistogram(argv[1], Chan);

            uSX3_outfile << i << "\t" << j << "\t" << pedestal << endl;

            Chan++;

        }
    }

    // Upstream SX3 Back Strips
    // Detectors 10 -> 11 (reversed)
    Chan = 185;
    nDets = 2;
    nStrips = 4;

    for (int i=10; i<12; i++) {

        for (int j=11; j>=8; --j) {

            pedestal = analyseHistogram(argv[1], Chan);
            
            uSX3_outfile << i << "\t" << j << "\t" << pedestal << endl;

            Chan++;

        }
    }
    
    // close Upstream SX3 Outfile
    uSX3_outfile.close();
    orderfile("SX3uPedestals.dat");
    cout << "Finished Upstream SX3 Pedestals" << endl;

    // Create Downstream SX3 Output txt file
    ofstream dSX3_outfile("SX3dPedestals.dat");
    
    // Check if the file was successfully opened
    if (!dSX3_outfile.is_open()) {
        cout << "Error opening downstream SX3 file!" << endl;
        return 1;
    }

    // Downstream SX3 Front Strips
    // Detectors 0 -> 11
    Chan = 289;

    for (int i=0; i<12; i++) {

        for (int j=0; j<8; j++) {

            pedestal = analyseHistogram(argv[1], Chan);

            dSX3_outfile << i << "\t" << j << "\t" << pedestal << endl;

            Chan++;

        }
    }

    // Dowmstream SX3 Back Strips
    // Detectors 0 -> 3  
    Chan = 385;
    
    for (int i=0; i<4; i++) {

        for (int j=8; j<12; j++) {

            pedestal = analyseHistogram(argv[1], Chan);

            dSX3_outfile << i << "\t" << j << "\t" << pedestal << endl;

            Chan++;

        }
    }

  // Downstream SX3 Back Strips
  // Detectors 4 -> 5 (reversed)
  Chan = 417;
  
  for (int i=4; i<6; i++) {

        for (int j=11; j>=8; --j) {

            pedestal = analyseHistogram(argv[1], Chan);
            
            dSX3_outfile << i << "\t" << j << "\t" << pedestal << endl;

            Chan++;

        }
    }

  // Upstream SX3 Back Strips
  // Detectors 6 -> 9
  Chan = 401;
  
  for (int i=6; i<10; i++) {

    for (int j=8; j<12; j++) {

      pedestal = analyseHistogram(argv[1], Chan);

      dSX3_outfile << i << "\t" << j << "\t" << pedestal << endl;

      Chan++;

    }
  }

  // Upstream SX3 Back Strips
  // Detectors 10 -> 11 (reversed)
  Chan = 425;

  for (int i=10; i<12; i++) {

        for (int j=11; j>=8; --j) {

            pedestal = analyseHistogram(argv[1], Chan);
            
            dSX3_outfile << i << "\t" << j << "\t" << pedestal << endl;

            Chan++;

        }
    }
    
    dSX3_outfile.close();
    orderfile("SX3dPedestals.dat");
    cout << "Finished Downstream SX3 Pedestals" << endl;

    inputFile->Close();

    delete inputFile;
    delete outputFile;

    return 0;
  
}
