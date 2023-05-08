#ifndef MYHEADER_H
#define MYHEADER_H

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
#include <iomanip>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TMath.h"
#include "CubicSpline.h"
//#include "analysis_functions.cxx"

using namespace std;

const int NMAX = 128;

TFile* inputFile;
TTree* inputTree;
TFile* outputFile;
TTree* outputTree;

// ===================== Data variables =====================
	//  === BB10 ===
	int   BB10Mul = 0;
    int   BB10Det[NMAX] = {0};
	int   BB10Strip[NMAX] = {0};
	int   BB10Channel[NMAX] = {0};
	int   BB10ADC[NMAX] = {0};
	float BB10Energy[NMAX] = {0};
	//
	// === QQQ5 ===
	int   QQQ5Mul = 0;
	int   QQQ5Upstream[NMAX] = {0};
	int   QQQ5Det[NMAX] = {0};
	int   QQQ5Ring[NMAX] = {0};
	int   QQQ5Sector[NMAX] = {0};
	float QQQ5RingEnergy[NMAX] = {0};
	float QQQ5SectorEnergy[NMAX] = {0};
	//
	// === SX3 ===
	int   SX3Mul = 0;
	int   SX3Upstream[NMAX] = {0};
	int   SX3Det[NMAX] = {0};
	int   SX3Sector[NMAX] = {0};
	float SX3SectorEnergy[NMAX] = {0};
	int   SX3Strip[NMAX] = {0};
	float SX3StripLeftEnergy[NMAX] = {0};
	float SX3StripRightEnergy[NMAX] = {0};
	float SX3StripEnergy[NMAX] = {0};
	float SX3StripPositionCal[NMAX] = {0};
	//
	// === TDC ===
	int tdcSilicon = 0;
	int tdcGRETINA = 0;
	int tdcRF = 0;
    int tdcSilicon_Div = 0;
    int tdcSilicon_GRETINA = 0;
    int tdcSilicon_Delay = 0;
    int tdcSilicon_Upstream = 0;
    //
    // Timestamps
	Long64_t timeStamp = 0;
	unsigned long long GRETINATimeStamp = 0;
    //
	// === GRETINA ===
	bool  foundGRETINA = 0;
	int   xtalsMul = 0;
	float xtals_cc[NMAX] = {0};
    float xtals_cc1[NMAX] = {0};
	int   xtals_crystalNum[NMAX] = {0};
	Long64_t  xtals_timestamp[NMAX] = {0};


    // Declare the output variables
    int Multi;
    int GMulti;
    float Excitation[NMAX];
    float Excitation_Cor[NMAX];
    float Etotal[NMAX];
    float LabTheta[NMAX];
    float DeltaE_Barrel[NMAX];
    float DeltaE_Endcap[NMAX];
    float Ethick_Barrel[NMAX];
    float Ethick_Endcap[NMAX];
    float Range_Barrel[NMAX];
    float Range_Endcap[NMAX];

    float Egamma[NMAX]; // test

    double rel_q_value ( double psi, double T3 ) {

        double q_val;
	    double M1 = 938.272, M2 = 89334.643, M3 = 938.272; // protons on Mo96 //938.783 //  938.27208816
        double P1, P3;
        double E1, E3;
        double T1 = 14.77;     //after energy loss of beam through half of the target, keV (right now just the beam energy) //14.7 is correct MeV


        E1 = T1 + M1; // total energy beam
        E3 = T3 + M3; // total energy ejectile
        P1 = sqrt ( ( E1*E1 ) - ( M1*M1 ) ); // beam momentum
        P3 = sqrt ( ( E3*E3 ) - ( M3*M3 ) ); // 3.1415

        q_val = M1 + M2 - M3 - sqrt ( ( M1*M1 ) + ( M2*M2 ) + ( M3*M3 ) + 2*M2*E1 - 2*E3* ( E1+M2 ) + 2*P1*P3*cos ( psi* ( TMath::Pi()/180.0 ) ) );

        return q_val;
    }

    double initial_proton_energy(double detected_energy, double material_length)    {

    //double E_p[38] = {0.005,0.01,0.02,0.03,0.04,0.05,0.1,0.15,0.2,0.25,0.5,
    //0.75,1,1.25,1.5,1.75,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
    //19,20};

    //double E_loss[38] = {5,5.1,6.23,6.96,7.5,7.89,8.44,7.93,7.17,6.41,4.15,3.19,
    //2.61,2.21,1.95,1.76,1.61,1.37,1.18,1.05,0.953,0.871,0.803,0.696,0.615,0.552,
    //0.501,0.46,0.423,0.393,0.367,0.345,0.326,0.309,0.293,0.28,0.267,0.256};

    double E_p[36] = {0.005,0.01,0.02,0.03,0.04,0.05,0.1,0.15,0.2,0.25,0.5,
	    0.7,1,1.2,1.5,1.7,2,2.5,3,3.5,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20};

	//double E_loss[36] = { 0.7529, 1.088, 1.601, 1.922, 2.13, 2.263, 2.412, 2.287, 
	//2.122, 1.966, 1.437, 1.2, 0.9787, 0.9037, 0.871, 0.7552, 0.6903, 0.6066, 0.5436, 0.4942, 
	//0.4543, 0.3933, 0.3486, 0.3141, 0.2866, 0.2642, 0.2454, 0.2294, 0.2156, 0.2036, 0.193, 0.1836, 0.1752, 0.1677, 0.1608, 0.1488};
/*
    0.0270, 0.0385, 0.0562, 0.0673, 0.0745, 0.0791, 0.0842, 0.0798, 0.0740, 0.0686
0.0501
0.0418
0.0340
0.0313
0.0278
0.0260
0.0237
0.0208
0.0187
0.0170
0.0156
0.0135
0.0120
0.0113
0.0108
0.0098
0.0091
0.0084
0.0079
0.0074
0.0070
0.0066
0.0063
0.0060
0.0058
0.0055
0.0051
*/
   
    // protons in CD2 target
       double E_loss[36] = { 0.1506, 0.2511, 0.3740, 0.4262, 0.46971093, 0.5001, 0.5157, 0.4608, 0.4064, 0.3619, 0.2376, 0.1900, 0.1489, 
        0.1348, 0.1146, 0.1044, 0.4064, 0.0779, 0.0676, 0.0599, 0.0539, 0.0452, 0.0390, 0.0345, 0.0309, 0.0281, 0.0258, 0.0239, 0.0222, 
        0.0208, 0.0196, 0.0185, 0.0176,0.0167, 0.0160, 0.0146};

    CubicSpline E_loss_spline(E_p, E_loss);
    double updated_energy = detected_energy;

    //Material length is in um!!!

    //incrementing through the material in 0.1um steps
    for(double dx = 0; dx < material_length; dx=dx+0.01)
    {
        updated_energy = updated_energy + (E_loss_spline(updated_energy)/1000.0);
    }

    return updated_energy;

}


std::vector<double>  hit_position_3D(string SX3_or_QQQ5, int upstream, int det, int strip, double pos_or_segment)
{
    vector<double> hit_pos;
    double x, y, z;


    if(SX3_or_QQQ5 == "SX3")
    {
        //correcting the weird detector numbering
        //Now it SHOULD be like a clock
        if (det != 0) det = det + 1;
        if (det == 0) det = 12;

        double rDetCenterD = 98.0;

        //How far from the center of the detector was the hit? (in the x-y plane)
        double strip_adj_mag = ( strip * 10.0 ) - 15.0; //mm

        //Angle from the 12 oclock detector to the hit detector
        double det_angle = 0.524 * det; //rads

        //angle the FACE of the detector makes in the x-y plane
        double perp_det_angle = det_angle + 3.14159/2.0; //rads

        //x-y position of the center of the hit detector
        double xDetCenter = rDetCenterD * sin ( det_angle );
        double yDetCenter = rDetCenterD * cos ( det_angle );

        //Adjusting the x-y position to match the strip
        double xAdj = strip_adj_mag * sin ( perp_det_angle );
        double yAdj = strip_adj_mag * cos ( perp_det_angle );

        //constructing the final x-y of the hit
        x = xDetCenter + xAdj;
        y = yDetCenter + yAdj;


        //now need to add z from the other function
        if(upstream == 1)
        {
            //z = -(8.2 + pos_or_segment*75.);
			z = pos_or_segment*75. - 5.9;  // was -8.2
        }
        if(upstream == 0)
        {
            z = z = pos_or_segment*75.0 + 9.6;
            //(pos_or_segment*75.0 - 3.0); // was -3.0
        }

        // cout << "SX3:  " << x << "  " << y << "  " << z << endl;

        //fill hit_pos vector mm
        hit_pos.push_back ( x );
        hit_pos.push_back ( y );
        hit_pos.push_back ( z );


        

    }

    if(SX3_or_QQQ5 == "QQQ5")
    {
        double radius[32] = {2.6475,2.9,3.1475,3.39,3.6275,3.86,4.0875,4.31,
        4.5275,4.74,4.9475,5.15,5.3475,5.54,5.7275,5.919,6.0875,6.26,6.4275,
        6.59,6.7475,6.9,7.0475,7.19,7.3275,7.46,7.5875,7.71,7.8275,7.94,8.0475,8.15};

        double theta_detCenter = 45. + (90. * det);
        double theta_segmentCenter = (theta_detCenter - 30.) + (pos_or_segment*20.);

        x = sin(theta_segmentCenter*(3.14159/180.))*radius[strip]*10.;
        y = cos(theta_segmentCenter*(3.14159/180.))*radius[strip]*10.;

        if(upstream==1)
        {
            z = -86.;
        }
        if(upstream==0)
        {
            z = 75.;
        }

        // cout << "QQQ5:  " << x << "  " << y << "  " << z << endl;

        //fill hit_pos vector
        hit_pos.push_back ( x );
        hit_pos.push_back ( y );
        hit_pos.push_back ( z );


    }
    
    return hit_pos;

}


vector<double>  hit_position_r_theta_phi(string SX3_or_QQQ5, int upstream, int det, int strip, double pos_or_segment)
{
    vector<double> hit_pos;
    double x, y, z;


    if(SX3_or_QQQ5 == "SX3")
    {
        //correcting the weird detector numbering
        //Now it SHOULD be like a clock
        if (det != 0) det = det + 1;
        if (det == 0) det = 12;

        double rDetCenterD = 98.8; //98

        //How far from the center of the detector was the hit? (in the x-y plane)
        double strip_adj_mag = -( ( strip * 10.0 ) - 15.0 ); //mm

        //Angle from the 12 oclock detector to the hit detector
        double det_angle = 0.524 * det; //rads

        //angle the FACE of the detector makes in the x-y plane
        double perp_det_angle = det_angle + 3.14159/2.0; //rads

        //x-y position of the center of the hit detector
        double xDetCenter = rDetCenterD * sin ( det_angle );
        double yDetCenter = rDetCenterD * cos ( det_angle );

        //Adjusting the x-y position to match the strip
        double xAdj = strip_adj_mag * sin ( perp_det_angle );
        double yAdj = strip_adj_mag * cos ( perp_det_angle );

        //constructing the final x-y of the hit
        x = xDetCenter + xAdj;
        y = yDetCenter + yAdj;


        //now need to add z from the other function
        if(upstream == 1)
        {
            //z = -(8.2 + pos_or_segment*75.);
			z = pos_or_segment*75. - 5.6; // my code already gives negative positions for upstream (8.2) was -3.3
        }
        if(upstream == 0)
        {
            z = pos_or_segment*75.0 + 9.6; // was +3.3 // 5.6
        }

        // cout << "SX3:  " << x << "  " << y << "  " << z << endl;

        // Converting to spherical polar coordinates
        double r, theta, phi;
        r = sqrt( pow(x,2) + pow(y,2) + pow(z,2) );
        phi = atan(-y/x);
        theta = acos(  z / sqrt( pow(x,2) + pow(y,2) + pow(z,2) ));

        // cout << r << "  " << theta << "  " << phi << endl;



        //fill hit_pos vector mm
        hit_pos.push_back ( r );
        hit_pos.push_back ( theta );
        hit_pos.push_back ( phi );

    }

    if(SX3_or_QQQ5 == "QQQ5")
    {
        double radius[32] = {2.6475,2.9,3.1475,3.39,3.6275,3.86,4.0875,4.31,
        4.5275,4.74,4.9475,5.15,5.3475,5.54,5.7275,5.919,6.0875,6.26,6.4275,
        6.59,6.7475,6.9,7.0475,7.19,7.3275,7.46,7.5875,7.71,7.8275,7.94,8.0475,8.15};

        double theta_detCenter = 45. + (90. * det);
        double theta_segmentCenter = (theta_detCenter - 30.) + (pos_or_segment*20.);

        x = sin(theta_segmentCenter*(3.14159/180.))*radius[strip]*10.;
        y = cos(theta_segmentCenter*(3.14159/180.))*radius[strip]*10.;

        if(upstream==1)
        {
            z = -86.;
        }
        if(upstream==0)
        {
            z = 75.;
        }

        // cout << "QQQ5:  " << x << "  " << y << "  " << z << endl;

		// Converting to spherical polar coordinates
        double r, theta, phi;
        r = sqrt( pow(x,2) + pow(y,2) + pow(z,2) );
        phi = atan(-y/x);
        theta = acos(  z / sqrt( pow(x,2) + pow(y,2) + pow(z,2) ));

        //fill hit_pos vector
        hit_pos.push_back ( r );
        hit_pos.push_back ( theta );
        hit_pos.push_back ( phi );


    }
    
    return hit_pos;
}

double proton_distance_through_target(vector<double> hit_pos) {

    double target_norm[3] = {0, 0.454, -0.891};
   // double half_targ = 0.343; //um
   double half_targ = 6.45; //um

    double proton_vector[3] = {hit_pos.at(0), hit_pos.at(1), hit_pos.at(2)};

    double numerator = (target_norm[0]*proton_vector[0]) + (target_norm[1]*proton_vector[1]) + (target_norm[2]*proton_vector[2]);
    double denominator = sqrt( pow(proton_vector[0],2) + pow(proton_vector[1],2) + pow(proton_vector[2],2));
    double angle = acos( (numerator*1.)/(denominator*1.) );

    double material_length = half_targ / cos(angle);

    return abs(material_length);
}

#endif