#include "GRootFunctions.h"
#include "GRootCommands.h"
#include "GCanvas.h"
#include "GPeak.h"

void Calibrate_Gretina(){

    gROOT->SetBatch(true);

    // outfile
    ofstream outfile;
    outfile.open("peak_locations.dat");

    // Grabs root files for each source  
    TFile *Co60 = TFile::Open("../root_outputs/Run0098_Gretina.root");
    TFile *Eu152 = TFile::Open("../root_outputs/Run0097_Gretina.root");
    TFile *Co56 = TFile::Open("../root_outputs/Run0096_Gretina.root");

    TFile *Gretina_Calb = TFile::Open("./Gretina_Calib.root","RECREATE");

    // Grabs 2D histogram of crystal vs energy 
    TH2D *hCo60 = (TH2D*)Co60->Get("hEgam_v_Crys");
    TH2D *hEu152 = (TH2D*)Eu152->Get("hEgam_v_Crys");
    TH2D *hCo56 = (TH2D*)Co56->Get("hEgam_v_Crys");

    // Array of Histograms for each source and crystal
    TH1D *hCo60_Crys[24];
    TH1D *hEu152_Crys[24];
    TH1D *hCo56_Crys[24];
   
    // Number of peaks in each source
    Int_t num_peaks[] = {2, 11, 13};

    // Source Energies
    Double_t** source_energy = new Double_t*[5];
    source_energy[0] = new Double_t[num_peaks[0]]{1173.240, 1332.508}; // Co60
    source_energy[1] = new Double_t[num_peaks[1]]{121.7817, 244.6974, 344.2785, 367.7891, 443.9606, 656.489, 778.9045, 964.057, 1112.076, 1212.948, 1408.018}; // Eu152
    source_energy[2] = new Double_t[num_peaks[2]]{846.770, 1037.843, 1175.101, 1238.2883, 1360.212, 1771.3567, 2015.215, 2034.791, 2598.500, 3009.645, 3202.029, 3451.232, 3548.05}; // Co56

    // Source Energy Errors    
    Double_t** source_energy_er = new Double_t*[5];
    source_energy_er[0] = new Double_t[num_peaks[0]]{0.003, 0.004}; // Co60
    source_energy_er[1] = new Double_t[num_peaks[1]]{0.0003, 0.0008, 0.0012, 0.0016, 0.0020, 0.005, 0.0024, 0.005, 0.003, 0.011, 0.003}; // Eu152
    source_energy_er[2] = new Double_t[num_peaks[2]]{0.002, 0.004, 0.004, 0.003, 0.004, 0.004, 0.005, 0.005, 0.004, 0.004, 0.004, 0.004}; // Co56

    // Fit range minimums
    Int_t** fit_min_range = new Int_t*[5];
    fit_min_range[0] = new Int_t[num_peaks[0]]{30, 30}; //co60
    fit_min_range[1] = new Int_t[num_peaks[1]]{5, 6, 8, 5, 5, 5, 6, 5, 8, 10, 10}; //eu152
    fit_min_range[2] = new Int_t[num_peaks[2]]{10, 6, 8, 10, 10, 10, 6, 10, 8, 10, 10, 10, 10}; //eu152
    // Fit range maximums
    Int_t** fit_max_range = new Int_t*[5];
    fit_max_range[0] = new Int_t[num_peaks[0]]{15, 15}; //co60
    fit_max_range[1] = new Int_t[num_peaks[1]]{5, 6, 8, 5, 5, 5, 6, 10, 8, 10, 10}; //eu152
    fit_max_range[2] = new Int_t[num_peaks[2]]{5, 6, 8, 5, 5, 5, 6, 10, 8, 10, 10, 10, 10}; //eu152

    // Containers
    Double_t** centroids = new Double_t*[5];
    centroids[0] = new Double_t[num_peaks[0]]{}; // co60
    Double_t** centroids_er = new Double_t*[5];
    centroids_er[0] = new Double_t[num_peaks[0]]{}; // co60

    int Crystal;
    int FirstCrystal = 21;
    double centroid, centroid_error;
    char hname_co60[64]; 
    char hname_eu152[32]; 
    char hname_co56[32]; 

    TH1D *hist[3][24];      

    for (int i=0; i<24; i++) { // loop over crystals

        Crystal = FirstCrystal + i;

        sprintf(hname_co60,"hCo60_Crys[%d]",i);
        sprintf(hname_eu152,"hEu152_Crys[%d]",i);
        sprintf(hname_co56,"hCo56_Crys[%d]",i);

        hist[0][i] = (TH1D*) hCo60->ProjectionX(hname_co60, Crystal, Crystal + 1);
        hist[1][i] = (TH1D*) hEu152->ProjectionX(hname_eu152, Crystal, Crystal + 1);
        hist[2][i] = (TH1D*) hCo56->ProjectionX(hname_co56, Crystal, Crystal + 1);
        
        for (int k=0; k<3; k++) { // loop over sources

            for (int j=0; j<num_peaks[k]; j++ ) { // loop over peaks

                //if (j==10) {continue;}

                GPeak *g = PhotoPeakFit(hist[k][i],source_energy[k][j] - fit_min_range[k][j],source_energy[k][j] + fit_max_range[k][j],"QR");

                centroid = g->GetCentroid();

                centroid_error = g->GetCentroidErr();      

                outfile << Crystal << " " <<  centroid << " " << centroid_error  << endl;
                //g->Delete();
            }
        } 
    }        

   // Gretina_Calb->Close();
 
}
