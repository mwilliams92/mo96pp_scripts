#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <string>
#include <stdio.h>
#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"

// Feed this a histogram and it will spit out the peak locations as vector in ascending order. Also set sigma and thresh to desired value
vector < Double_t > findpeaks(TH1* h, Double_t sigma, Double_t thresh){

 	vector< Double_t > peaks;
 	int npeaks = 100;

 	h->GetXaxis()->SetRange(50,2000);
 
 	TSpectrum * s = new TSpectrum(npeaks);
 	Int_t nfound = s->Search(h, sigma, "", thresh);

    Double_t * xpeaks = s->GetPositionX();  
    for (int i=0; i<=nfound; i++){
       	if (xpeaks[i]>100 and xpeaks[i]<2000){
          	peaks.push_back(xpeaks[i]);  
      	}
    }
	// sorts vector into ascending order
	sort(peaks.begin(), peaks.end()); 
 	return peaks;
}

// Feed this a vector of peak locations and it will spit out slopes and offset at a pointer to an array
double * findcalib(vector < Double_t > peaks)	{
	// Eu152 lines:
	vector < Double_t > energies = {121.782, 244.697, 344.279, 778.905, 964.057, 1408.018};
	
	TGraph *gr = new TGraph(5, &peaks[0], &energies[0]);
	gr->SetMarkerStyle(20);
	gr->Draw("goff");

	TF1 *fit = new TF1("fit","pol1");
	//fit->FixParameter(0,0);
	gr->Fit(fit,"Q");
	static double par[2];
	fit->GetParameters(par);

	return par;
}

double * findcalibration ( const int n, vector < Double_t > Raw, vector < Double_t > Energies  ) {

    TGraph *gr = new TGraph(n, &Raw[0], &Energies[0]);
    gr->SetMarkerStyle(20);
    gr->Draw("goff");

    TF1 *fit = new TF1("fit","pol1");
	gr->Fit(fit,"Q");
	static double par[2];
	fit->GetParameters(par);

	return par;
}

// Feed this a layer to tell it which root file to use and the channels corresponding to the first ring and sector
void calibrate_Eu152(TString inputFile, TString outputCalFile, TString outputRootFile){

	TFile *file=TFile::Open(inputFile);

    ofstream txtoutfile;
	txtoutfile.open (outputCalFile);
			
    TH2D *hEu152 = (TH2D*)file->Get("hEgam1_v_Crys[0]");

	const int nCrys = 50;
    int chan;
    char hname[64];
    
    TFile *outfile = new TFile(outputRootFile,"RECREATE");
    TH1F *hEu152_Crys[nCrys];

    for (int i=0; i<nCrys; i++) {

        sprintf(hname, "Eu152 Crystal %d", i);
        hEu152_Crys[i] = (TH1F*) hEu152->ProjectionX(hname, i, i);

        //hEu152_Crys[i]->Draw();
        hEu152_Crys[i]->GetXaxis()->SetRangeUser(0,2000);
		vector < Double_t > peaks = findpeaks(hEu152_Crys[i],2.,0.09);	// Uses findpeaks function to return vector of peak loactions in ascending order. //0.09 is good for finding 8

        //chan = i*40 + 9;

        cout <<  peaks.size() << " Peaks found in crystal " << i << endl; 

/*
        if (peaks.size() != 5) {
			printf("Channel %d (Ring %d) not calibrated! Found %lu peaks instead of 5.\n", i, rID, peaks.size());
			txtoutfile << chan << "\t" << -1 << "\t" << 0 << endl;	// Turns off channel if the number of peaks found is not 5.
		} else {
			double *par = findcalib(/peaks);	// Uses findcalib function to return a pointer to an array containing the calibration offset and slope for each channel
			txtoutfile << chan << "\t" << std::fixed << setprecision(5) << *(par + 0) << "\t" << *(par + 1) << endl;
		}
*/
        hEu152_Crys[i]->Write();

    }

    outfile->Write();
    outfile->Close();

    txtoutfile.close();
}

void FindPeaks_Eu152(TString inputFile, TString PeakFile, TString outputRootFile){

    ofstream txtoutfile;                         // text file to write peak locations
	txtoutfile.open (PeakFile, ios::app);   // open in append mode
    if (!txtoutfile) {

        cout << "No such file found";
        return 0;

    } else {

        const int Eu152_peaks = 6;                           // number of Co60 peaks
        double Eu152_energies[Eu152_peaks] = {121.7817, 244.6974, 344.2785, 778.9045, 964.057, 1408.018};  // Energies of Co60 lines
        double Eu152_error[Eu152_peaks] = { [0 ... 5] = 0.1 };       // Error in known Co60 line energies

	    TFile *file = TFile::Open(inputFile);       // open sorted gretina file

        if (!file) { 

            cout << "Eu152 data file not found!" << endl;
            return 0;

        } else { 

	        const int nCrys = 48;   // number of crystals
            int chan;               // channel for sort code (cc1 chan = crystal*40 + 9)
            char hname[64];         // name of indv crystal hists
    
            TFile *outfile = new TFile(outputRootFile,"RECREATE");  // Create output root file for Co60 spectra
           
            TH2F *hEu152 = (TH2F*)file->Get("hEu152");    // grabs 2D hist of crystal vs energy
            hEu152->Write();                             // saves to current root file
            TH1F *hEu152_Crys[nCrys];                    // hists for indv crystals

            for (int i=1; i<=nCrys; i++) {   // loop over crystals indexed from 1

                sprintf(hname, "Eu152 Crystal %d", i);                           // separate name for each cystal hist
                hEu152_Crys[i] = (TH1F*) hEu152->ProjectionX(hname, i+1, i+1);    // project 2D onto energy axis for each crystal (crystal 1 = bin 2)

                hEu152_Crys[i]->GetXaxis()->SetRangeUser(0,2000);                // sets X-axis range optimised to find Co60 peaks.               
		        
                vector < Double_t > peaks;
                
                double thresh;
                bool found = false;
                int j=0;

                //for(int j=0; j<20; j++) {
                while (!found) {
                    thresh = 0.5 - 0.005*j;
                    peaks = findpeaks(hEu152_Crys[i],2.,thresh);	// Uses findpeaks function to return vector of peak loactions in ascending order.
                    //cout << "Peaks = " << peaks.size() << endl;
                    if (peaks.size() == 6 || peaks.size() == 0 || j == 99) { found = true; }
                    j++;
                }
                cout << j << ", " << peaks.size() << endl;
                
                chan = (i-1)*40 + 9;    // Channel ID is indexed with crsytal 1 = 0

                if (peaks.size() != 6) {
			        printf("Skipping Crystal %d (channel %d). Found %lu peaks instead of 6.\n", i, chan, peaks.size());   // Print statement if 2 peaks are not found
			        continue;	                                                                                          // Skips crystal
		        }
                
                for (int j=0; j<Eu152_peaks; j++) {
                    txtoutfile << i << "," << Eu152_energies[j] << "," << peaks.at(j) << "," << Eu152_error[j] << "," << 0.5 << endl;
                }
                
                hEu152_Crys[i]->Write(); // writes indv crystal hist to output root file

            } // end loop over crystals

            outfile->Write();   // write Co60 root file
            outfile->Close();   // close Co60 root file

        } // end if found input root file

        txtoutfile.close();

    }   // end if text file open
}   // end Eu152 calib


void FindPeaks_Co60(TString inputFile, TString PeakFile, TString outputRootFile){

    ofstream txtoutfile;                         // text file to write peak locations
	txtoutfile.open (PeakFile, ios::app);   // open in append mode
    if (!txtoutfile) {

        cout << "No such file found";
        return 0;

    } else {

        const int Co60_peaks = 2;                           // number of Co60 peaks
        double Co60_energies[Co60_peaks] = {1113., 1173.};  // Energies of Co60 lines
        double Co60_error[Co60_peaks] = {0.1, 0.1};         // Error in known Co60 line energies

	    TFile *file = TFile::Open(inputFile);       // open sorted gretina file

        if (!file) { 

            cout << "Co60 data file not found!" << endl;
            return 0;

        } else { 

	        const int nCrys = 48;   // number of crystals
            int chan;               // channel for sort code (cc1 chan = crystal*40 + 9)
            char hname[64];         // name of indv crystal hists
    
            TFile *outfile = new TFile(outputRootFile,"RECREATE");  // Create output root file for Co60 spectra
            
            TH2F *hCo60 = (TH2F*)file->Get("hCo60");    // grabs 2D hist of crystal vs energy
            hCo60->Write();                             // saves to current root file
            TH1F *hCo60_Crys[nCrys];                    // hists for indv crystals

            for (int i=1; i<nCrys; i++) {   // loop over crystals indexed from 1

                sprintf(hname, "Co60 Crystal %d", i);                           // separate name for each cystal hist
                hCo60_Crys[i] = (TH1F*) hCo60->ProjectionX(hname, i+1, i+1);    // project 2D onto energy axis for each crystal (crystal 1 = bin 2)

                hCo60_Crys[i]->GetXaxis()->SetRangeUser(0,2000);                // sets X-axis range optimised to find Co60 peaks.               
		        vector < Double_t > peaks = findpeaks(hCo60_Crys[i],2.,0.5);	// Uses findpeaks function to return vector of peak loactions in ascending order.

                chan = (i-1)*40 + 9;    // Channel ID is indexed with crsytal 1 = 0

                if (peaks.size() != 2) {
			        printf("Skipping Crystal %d (channel %d). Found %lu peaks instead of 2.\n", i, chan, peaks.size());   // Print statement if 2 peaks are not found
			        continue;	                                                                                          // Skips crystal
		        }
        
                for (int j=0; j<Co60_peaks; j++) {
                    txtoutfile << i << "," << Co60_energies[j] << "," << peaks.at(j) << "," << Co60_error[j] << "," << 0.5 << endl;
                }

                hCo60_Crys[i]->Write(); // writes indv crystal hist to output root file

            } // end loop over crystals

            outfile->Write();   // write Co60 root file
            outfile->Close();   // close Co60 root file

        } // end if found input root file

        txtoutfile.close();

    }   // end if text file open
}   // end Co60 calib

void generate2Dmatrix(TString inputFile, TString source) {

    TFile *file = TFile::Open(inputFile); 
    TTree *data = (TTree*) file->Get("teb");

    ostringstream fname;
    fname << source << ".root";
    ostringstream hname;
    hname << "h" << source;
    ostringstream htitle;
    htitle << source << " Cystal vs Energy";
    ostringstream hdraw;
    hdraw << "xtals.crystalNum:xtals.cc1>>h" << source;

    TFile *outputfile = new TFile(fname.str().c_str(),"RECREATE");   
    TH2F *Crys_v_CC1 = new TH2F(hname.str().c_str(),htitle.str().c_str(),10000,0,10000,50,0,50);    // create 2D hist for crystal vs cc1
    cout << "Generating 2D Histogram. Be Patient, this may take a while..." << endl;                // a big root file will take a while to draw
    data->Draw(hdraw.str().c_str(), "","goff");                                                      // drawing 2D hist from sorted data
    cout << "Writing 2D histogram...";
    outputfile->cd();
    Crys_v_CC1->Write();                                                                            // writing 2D hist to output root file
    cout << " Done!" << endl;
    outputfile->Close();

}

void calibrateGRETINA()	{

	//gROOT->SetBatch(kTRUE); // set Batch mode to stop unwanted canvases appearing

    //ofstream txtoutfile;                         // text file to write peak locations
    //txtoutfile.open("Peaks.csv", ios::trunc);

    //txtoutfile << "crys,source,location,source_error,location_error" << endl;

    //generate2Dmatrix("/mnt/d/sorted/Run0098_gretina.root","Co60"); // Comment-out if already have the 2D hist needed to find peaks
    //generate2Dmatrix("/mnt/d/sorted/Run0097_gretina.root","Eu152");

	// Find peak locations for each source in each crystal, appends peak locations to Peaks.dat and saves histograms a .root file 
    //FindPeaks_Co60("Co60.root","Peaks.csv","Co60_IndvCrystals.root");           // comment out if already have peak locations in Peaks.txt
	//FindPeaks_Eu152("Eu152.root","Peaks.csv","Eu152_IndvCrystals.root");

    //txtoutfile.close();

    
    //gROOT->SetBatch(kTRUE);
    
    
    auto rdf = ROOT::RDF::MakeCsvDataFrame("Calibration_Data.csv",true,',');

    auto colNames = rdf.GetColumnNames();
    // Print columns' names
    for (auto &&colName : colNames) std::cout << colName << std::endl;

    ofstream txtoutfile;
	txtoutfile.open ("calpars.txt");
    

    TFile *file = TFile::Open("cal.root","RECREATE");

    /*
    for (int i=1; i<=48; i++) {

        sprintf(cut,"crys==%d",i);
        auto filter = rdf.Filter(cut);   
        // auto gr = filter.GraphAsymmErrors("peak", "source","peak_er","peak_er","source_er","source_er");
        using f = float;
        auto myGAE2 = myDf.GraphAsymmErrors<f, f, f, f, f, f>("peak", "source","peak_er","peak_er","source_er","source_er");    


        gr->DrawClone("AP");
        TF1 *f1 = new TF1("f1","pol1");
        gr->Fit("f1","Q");

        cout << f1->GetChisquare()/f1->GetNDF() << endl;
        txtoutfile << (i-1)*40 + 9 << " " << f1->GetParameter(0) << " " << f1->GetParameter(1) << endl;

        gr->SetMarkerStyle(20);
        gr->Write();
        f1->Write();
        
    }
    */
    FILE *fp;
    fp = fopen("Gretina_Calibration_Manual.txt","rb");

    Int_t Crystal;
    Float_t Peak, Source, Peak_Error, Source_Error;

    TFile *hfile = TFile::Open("cal.root","RECREATE");
    TTree *tree = new TTree("tree","Eu152 Calibration");
    tree->Branch("Crystal",&Crystal,"Crystal/I");
    tree->Branch("Peak",&Peak,"Peak/F");
    tree->Branch("Source",&Source,"Source/F");
    tree->Branch("Peak_Error",&Peak_Error,"Peak_Error/F");
    tree->Branch("Source_Error",&Source_Error,"Source_Error/F");
      
    char line[80];
    while (fgets(line,80,fp)) {
      sscanf(&line[0],"%d %f %f %f %f", &Crystal, &Peak, &Source, &Peak_Error, &Source_Error);
      tree->Fill();
    }

    tree->Print();
    tree->Write();
    fclose(fp);

    char cut[16];

    for (int i=1; i<=48; i++) {

        sprintf(cut,"Crystal==%d",i);
        int n = tree->Draw("Peak:Source:Peak_Error:Source_Error",cut,"goff"); 
        TGraphErrors *g = new TGraphErrors(n,tree->GetV1(),tree->GetV2(),tree->GetV3(),tree->GetV4());
        //g->SetMarkerStyle(20);
        //g->Draw("ape");
        TF1 *f1 = new TF1("f1","pol1");
        f1->SetNpx(10000);
        g->Fit("f1","Q0");

        txtoutfile << (i-1)*40 + 9 << " " << f1->GetParameter(0) << " " << f1->GetParameter(1) << endl;

        cout << i << " Chi-Sq = " << f1->GetChisquare()/f1->GetNDF() << endl;

    }

   
    
}
