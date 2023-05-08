#include <iostream>
#include <iomanip>
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TFitResult.h"
#include <cmath>
#include "Math/PdfFuncMathCore.h"
#include "Math/MinimizerOptions.h"
//#include "/opt/root_v6.14.06/include/Math/MinimizerOptions.h"

using namespace std;

Double_t RadwareEfficiency(Double_t *dim, Double_t *par) {

	if (dim[0] < 1) {
		return 0;
	}

	Double_t x = dim[0];

	Double_t A = par[0];
	Double_t B = par[1];
	Double_t C = par[2];
	Double_t D = par[3];
	Double_t E = par[4];
	Double_t F = par[5];
	Double_t G = par[6];
	Double_t H = par[7];

	return H * exp(pow((pow(A + B * log(x / 100.) + C * log(x / 100.) * log(x / 100.), -G) + pow(D + E * log(x / 1000.) + F * log(x / 1000.) * log(x / 1000.), -G)), -1/G));
}

Double_t CarlottaEfficiency(Double_t *dim, Double_t *par) {

        if (dim[0] < 1) {
                return 0;
        }

        Double_t x = dim[0];

        Double_t A = par[0];
        Double_t B = par[1];
        Double_t C = par[2];
        Double_t D = par[3];
        Double_t E = par[4];
        Double_t F = par[5];

        return exp(A + B*x*0.05 + C*log(x*0.05)/(x*0.05)+ D/(x*0.05) + E/pow((x*0.05),2) + F/pow((x*0.05),3));
}

Double_t DirkEfficiency( Double_t *dim, Double_t *par) {

    if (dim[0] < 1) {

        return 0;
    }

    Double_t x = dim[0];
    
    return par[0]*pow(x+100.,par[1]);


}


void Singles(){

	const int n=2;
	Double_t x[n] = {0,8000};
	Double_t y[n] = {0,1};

	gStyle->SetLineWidth(2);
	TCanvas *c1 = new TCanvas();
	c1->SetTicky();
	c1->SetTickx();
	c1->SetLeftMargin(0.17); c1->SetBottomMargin(0.17);
	TGraph *g1 = new TGraph(n,x,y);
    g1->SetTitle("GRETINA Singles Efficiency");
	g1->GetXaxis()->SetTitle("Energy [keV]");
	g1->GetYaxis()->SetTitle("Efficiency");
	g1->GetXaxis()->SetTitleSize(0.07); g1->GetXaxis()->SetLabelSize(0.05); g1->GetXaxis()->CenterTitle();
	g1->GetYaxis()->SetTitleSize(0.07); g1->GetYaxis()->SetLabelSize(0.05); g1->GetYaxis()->CenterTitle();
	g1->GetYaxis()->SetRangeUser(0,0.24);
	g1->GetXaxis()->SetRangeUser(0,4000);
	g1->Draw("AP");

	TGraphErrors *g_Eu152 = new TGraphErrors("Eu152andCo56_Singles.txt");
	g_Eu152->SetMarkerStyle(20);
	g_Eu152->SetMarkerColor(1);
	g_Eu152->Draw("Psame");

	TF1 *func = new TF1("fit",RadwareEfficiency,100,4600,8);
	
    double par[8] = {7.04,0.70,0,5.273,-0.863,0.01,11,1.25e-4};
    double par_lower[8] = {5.0, 0.0, 0.0, 5.0, -5, -1 , -1 , 1.0e-6};
    double par_upper[8] = {25.0, 20.0, 0.0, 20.0, 5, 5 , 10, 1.0e-3};

	//double par[2] = {3.6,-0.62};
	//double par_lower[3] = {0.,-1.0};
	//double par_upper[3] = {10.,0.0};

	for (int i=0; i<8; i++) {
		func->FixParameter(i, par[i]);
        func->SetParLimits(i, par_lower[i], par_upper[i]);
	}	

    //func->Draw();
	g_Eu152->Fit(func,"R");
	func->SetLineColor(2);
    cout << func->GetChisquare()/func->GetNDF() << endl;
    cout << func->GetProb() << endl;

	cout << "Co60: " << func->Eval(1332.508) << endl;

	
/*
	TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
	legend->AddEntry(g_Singles,"Add-Back Data","lep");
	legend->AddEntry(func2,"Fit to Data","l");
	legend->AddEntry(g_Eu152,"GEANT4","l"); 
	legend->AddEntry(g_G4Scaled,"GEANT4 (Scaled 81.92%)","l");
	legend->Draw();
*/

    cout << "608.7 keV =  " << func->Eval(608.7) << endl;
	cout << "778.2 keV = " << func->Eval(778.2) << endl;
	cout << "719.6 keV = " << func->Eval(719.6) << endl;
	cout << "849.0 keV = " << func->Eval(849.0) << endl;
	cout << "1091.3 keV = " << func->Eval(1091.3) << endl;
	cout << "812.6 keV = " << func->Eval(812.6) << endl;

	cout << "2741.5 keV = " << func->Eval(2741.5) << endl;


	cout << "204.0 keV = " << func->Eval(204.0) << endl;
	cout << "765.0 keV = " << func->Eval(765.0) << endl;
	cout << "948.0 keV = " << func->Eval(948.0) << endl;


	cout << "4439 keV = " << func->Eval(4439) << endl;

	double Energy, Eff;
/*
	for(int i=0; i<36 ; i++) {

		Energy = i*100.;
		Eff = func->Eval(Energy);

		cout << Energy << "," << Eff << endl;
	}
	*/
}
