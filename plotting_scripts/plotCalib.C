#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"

// Plot QQQ5 rings
void plotRings( TFile *file, TCut cut )	{

	TTree *data = (TTree*)file->Get("data");
	TCanvas *cRings = new TCanvas("cRings");
	TH2D *hRings = new TH2D("hRings","Ring Number vs Energy",32,0,32,1000,0,10000);
	data->Draw("QQQ5RingEnergy:QQQ5Ring>>hRings",cut,"colz");
}
// Plot QQQ5 sectors
void plotSectors(TFile *file, TCut cut)	{

	TTree *data = (TTree*)file->Get("data");
	TCanvas *cSectors = new TCanvas("cSectors");
	TH2D *hSectors = new TH2D("hSectors","Sector Number vs Energy",4,0,4,500,0,10000);
	data->Draw("QQQ5SectorEnergy:QQQ5Sector>>hSectors",cut,"colz");

}
// Plot BB10 strips
void plotStrips(TFile *file, TCut cut)	{

	TTree *data = (TTree*)file->Get("data");
	TCanvas *cStrips = new TCanvas("cStrips");
	TH2D *hStrips = new TH2D("hStrips","Strip Number vs Energy",8,0,8,1000,0,10000);
	data->Draw("BB10Energy:BB10Strip>>hStrips",cut,"colz");

}

void plotSX3backStrips(TFile *file, TCut cut)	{

	TTree *data = (TTree*)file->Get("data");
	TCanvas *cStrips = new TCanvas("cStrips");
	TH2D *hStrips = new TH2D("hStrips","Strip Number vs Energy",4,0,4,500,0,10000);
	data->Draw("SX3SectorEnergy:SX3Sector>>hStrips",cut,"colz");

}

void LeftvRight(TFile *file, TCut cut)	{

	TTree *data = (TTree*)file->Get("data");
	TCanvas *cLeftvRight = new TCanvas("cLeftvRight");
	TH2D *hLeftvRight = new TH2D("hLeftvRight","Left v Right",2000,0,2000,2000,0,2000);
	data->Draw("SX3StripLeftADC:SX3StripRightADC>>hLeftvRight",cut,"colz");

}

// Plot QQQ5 - give it a detector number
void plotQQQ5(Int_t detector){

	Int_t upstream=0;

	// if detecotr < 4 then ask whether upstream or downstream
	if (detector < 4)	{
		char response;
  		cout << "Upstream? (U) or Downstream (D)? ";
  		cin >> response;
  		if (toupper( response ) == 'U') upstream = 1;
	}

	char c[32];
	sprintf(c,"QQQ5Det==%d && QQQ5Upstream==%d", detector, upstream);
	TCut cut = c;

	if (upstream==1)	{

		TFile *file = TFile::Open("../OutputFolder/cal228th_0deg_US.root");
		plotRings(file, cut);
		plotSectors(file, cut);

	} else if (upstream!=1 && (detector==0 || detector==1))	{

		TFile *file = TFile::Open("../OutputFolder/cal228th_180deg_DS_noBB10.root");
		plotRings(file, cut);
		plotSectors(file, cut);
		
	} else if (upstream!=1 && (detector==2 || detector==3))	{

		TFile *file = TFile::Open("../OutputFolder/cal228th_QQQ5_E1_a.root");
		plotRings(file, cut);
		plotSectors(file, cut);

	} else if (upstream!=1 && (detector==4 || detector==5) )	{

		TFile *file = TFile::Open("../OutputFolder/cal228th_QQQ5_E2.root");
		plotRings(file, cut);
		plotSectors(file, cut);

	} else if (detector<0 || detector>5) { printf("Invalid detector number. Needs to be 0, 1, 2, 3, 4 or 5\n"); }

}

void plotBB10(Int_t detector){

	char c[32];
	sprintf(c,"BB10Det==%d", detector);
	TCut cut = c;

	if (detector>=0 && detector<4)	{

		TFile *file = TFile::Open("../OutputFolder/cal228th_225deg_DS_BR.root");
		plotStrips(file, cut);
		
	} else if (detector>=4 && detector<8)	{

		TFile *file = TFile::Open("../OutputFolder/cal228th_135deg_DS_BL.root");
		plotStrips(file, cut);

	} else if (detector<0 || detector>7) { printf("Invalid detector number. Needs to be 0, 1, 2, 3, 4, 5, 6 or 7\n"); }

}

void plotSX3uBack(Int_t detector){

	char c[32];
	sprintf(c,"SX3Det==%d && SX3Upstream==1 && SX3Strip==0", detector);
	TCut cut = c;

	if (detector>=0 && detector<5)	{

		TFile *file = TFile::Open("../../OutputFolder/cal228th_315deg_US_BR.root");
		plotSX3backStrips(file, cut);
		
	} else if (detector>5 && detector<12)	{

		TFile *file = TFile::Open("../../OutputFolder/cal228th_45deg_US_BL.root");
		plotSX3backStrips(file, cut);

	} else if (detector<0 || detector>11 || detector==5) { printf("Invalid detector number. Needs to be 0->11 except 5\n"); }

}

void plotSX3dBack(Int_t detector){

	char c[32];
	sprintf(c,"SX3Upstream==0 && SX3Det==%d", detector);
	TCut cut = c;

	if (detector>=0 && detector<5)	{

		TFile *file = TFile::Open("../../OutputFolder/cal228th_270deg_BR_noBB10.root");
		plotSX3backStrips(file, cut);
		
	} else if (detector>4 && detector<11)	{

		TFile *file = TFile::Open("../../OutputFolder/cal228th_90deg_BL_noBB10.root");
		plotSX3backStrips(file, cut);

	} else if (detector<0 || detector>10) { printf("Invalid detector number. Needs to be 0->10\n"); }

}



void plotSX3uLeftvRight(Int_t detector, Int_t back){

	char c[32];
	sprintf(c,"SX3Upstream==1 && SX3Det==%d && SX3Sector==%d", detector, back);
	TCut cut = c;

	if (detector>=0 && detector<5)	{

		TFile *file = TFile::Open("../OutputFolder/cal228th_315deg_US_BR.root");
		LeftvRight(file, cut);
		
	} else if (detector>5 && detector<12)	{

		TFile *file = TFile::Open("../OutputFolder/cal228th_45deg_US_BL.root");
		LeftvRight(file, cut);

	} else if (detector<0 || detector>11 || detector==5) { printf("Invalid detector number. Needs to be 0->11 except 5\n"); }

}
