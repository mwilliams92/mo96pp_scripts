// g++ -std=c++11 -o merge mergeRootFiles.cpp `root-config --cflags --libs`
#include "TChain.h"
#include "TFile.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TString.h"
#include "TSystemFile.h"
#include <iostream>

using namespace std;

void mergeRootFiles(const char* outputFileName, int firstRun, int lastRun) 
{
  int nFiles = lastRun - firstRun; // Replace with the actual number of files
  //int firstRun = 59;
  int run;
 
 std::vector<std::string> fileNames;
  for (int i = 0; i <= nFiles; i++) {
    run = firstRun + i;
    fileNames.push_back(Form("Run%04d_sorted.root", run));
    cout << "Adding Run " << run << " to merge" << endl;
  }

  TChain chain1("t1");
  TChain chain2("t2");

 for (const auto& fileName : fileNames) {

   TFile* inFile = TFile::Open(fileName.c_str(), "READ");
    if (!inFile || inFile->IsZombie()) {
      std::cerr << "Could not open file: " << fileName << std::endl;
      continue;
    }
    TTree* inTree1 = (TTree*)inFile->Get("t1");
    if (!inTree1) {
      std::cerr << "Could not find tree in file: " << fileName << std::endl;
      inFile->Close();
      continue;
    }
    TTree* inTree2 = (TTree*)inFile->Get("t2");
    if (!inTree2) {
      std::cerr << "Could not find tree in file: " << fileName << std::endl;
      inFile->Close();
      continue;
    }
    
    Long64_t entries1 = inTree1->GetEntries();
    Long64_t entries2 = inTree2->GetEntries();
   
    chain1.Add(inFile->GetName(), entries1);
    chain2.Add(inFile->GetName(), entries2);
    
    inFile->Close();
  }

    TFile* outFile = TFile::Open(outputFileName, "RECREATE");
    if (!outFile) {
        printf("Error creating output file %s\n", outputFileName);
        return;
    }

    TTree* outTree1 = chain1.CloneTree(-1, "fast");
    outTree1->Write();
    TTree* outTree2 = chain2.CloneTree(-1, "fast");
    outTree2->Write();

    outFile->Close();
    printf("Merged TTree saved in %s\n", outputFileName);

}

int main(int argc, char* argv[])
{
    if (argc != 4) {
        printf("Usage: mergeRootFiles outputFileName.root firstRun lastRun\n");
        return 1;
    }

    const char* outputFileName = argv[1];
    int firstRun = std::atof(argv[2]);
    int lastRun = std::atof(argv[3]);

    mergeRootFiles(outputFileName, firstRun, lastRun);

    return 0;
}