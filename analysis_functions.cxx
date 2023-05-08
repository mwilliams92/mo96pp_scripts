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

#include "CubicSpline.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TChain.h"

double rel_q_value ( double psi, double T3 )
{
    double q_val;
	double M1 = 938.272, M2 = 89334.643, M3 = 938.272; // protons on Mo96 //938.783 //  938.27208816
    double P1, P3;
    double E1, E3;
    double T1 = 14.937;     //after energy loss of beam through half of the target, keV (right now just the beam energy) //14.7 is correct MeV
    //double T1 = 21.1;

	// BeamE = 14940, 14030, 21300

    E1 = T1 + M1; // total energy beam
    E3 = T3 + M3; // total energy ejectile
    P1 = sqrt ( ( E1*E1 ) - ( M1*M1 ) ); // beam momentum
    P3 = sqrt ( ( E3*E3 ) - ( M3*M3 ) ); // 3.1415

    q_val = M1 + M2 - M3 - sqrt ( ( M1*M1 ) + ( M2*M2 ) + ( M3*M3 ) + 2*M2*E1 - 2*E3* ( E1+M2 ) + 2*P1*P3*cos ( psi* ( TMath::Pi()/180.0 ) ) );

    return q_val;
}

double rel_q_value_O16 ( double psi, double T3 )
{
    double q_val;
	double M1 = 938.272, M2 = 14899.3, M3 = 938.272; // protons on O16 
    double P1, P3;
    double E1, E3;
    double T1 = 14.937;     //after energy loss of beam through half of the target, keV (right now just the beam energy) //14.7 is correct MeV
    //double T1 = 21.1;


    E1 = T1 + M1; // total energy beam
    E3 = T3 + M3; // total energy ejectile
    P1 = sqrt ( ( E1*E1 ) - ( M1*M1 ) ); // beam momentum
    P3 = sqrt ( ( E3*E3 ) - ( M3*M3 ) ); // 3.1415

    q_val = M1 + M2 - M3 - sqrt ( ( M1*M1 ) + ( M2*M2 ) + ( M3*M3 ) + 2*M2*E1 - 2*E3* ( E1+M2 ) + 2*P1*P3*cos ( psi* ( TMath::Pi()/180.0 ) ) );

    return q_val;
}

double rel_q_value_C12 ( double psi, double T3 )
{
    double q_val;
	double M1 = 938.272, M2 = 11178, M3 = 938.272; // protons on carbon

    double P1, P3;
    double E1, E3;
    double T1 = 14.937;     //after energy loss of beam through half of the target, keV (right now just the beam energy) //14.7 is correct MeV
    //double T1 = 21.1;        // 14.77;
            //
	// BeamE = 14940, 14030, 21300

    E1 = T1 + M1; // total energy beam
    E3 = T3 + M3; // total energy ejectile
    P1 = sqrt ( ( E1*E1 ) - ( M1*M1 ) ); // beam momentum
    P3 = sqrt ( ( E3*E3 ) - ( M3*M3 ) ); // 3.1415

    q_val = M1 + M2 - M3 - sqrt ( ( M1*M1 ) + ( M2*M2 ) + ( M3*M3 ) + 2*M2*E1 - 2*E3* ( E1+M2 ) + 2*P1*P3*cos ( psi* ( TMath::Pi()/180.0 ) ) );

    return q_val;
}

double rel_q_value_D2 ( double psi, double T3 )
{
    double q_val;
	double M1 = 938.272, M2 = 1876.124, M3 = 938.272; // protons on carbon

    double P1, P3;
    double E1, E3;
    double T1 = 21.1;     //after energy loss of beam through half of the target, keV (right now just the beam energy) //14.7 is correct MeV
            // 14.77;
            //
	// BeamE = 14940, 14030, 21300

    E1 = T1 + M1; // total energy beam
    E3 = T3 + M3; // total energy ejectile
    P1 = sqrt ( ( E1*E1 ) - ( M1*M1 ) ); // beam momentum
    P3 = sqrt ( ( E3*E3 ) - ( M3*M3 ) ); // 3.1415

    q_val = M1 + M2 - M3 - sqrt ( ( M1*M1 ) + ( M2*M2 ) + ( M3*M3 ) + 2*M2*E1 - 2*E3* ( E1+M2 ) + 2*P1*P3*cos ( psi* ( TMath::Pi()/180.0 ) ) );

    return q_val;
}


double integrate_func (vector<double> xValue, vector<double> yValue, double xStart, double xStop)
{
    int found_start_flag = 0;
    int nStart, nStop;

    // Finding index of start of function
    for(int i=0; i<xValue.size(); i++)
    {
        if(xValue.at(i)>=xStart && found_start_flag==0){
            nStart = i;
            found_start_flag = 1;
        }
        if(xValue.at(i)<=xStop){
            nStop = i;
        }
    }

    double nPoints = nStop-nStart;
    double coefficient, sum;
    int loop_number = 0;

    for(int n = nStart; n <= nStop; n++){

        if( loop_number % 2 == 0) coefficient = 2;
        if( loop_number % 2 != 0) coefficient = 4;
        if(loop_number == 0 || loop_number == nPoints) coefficient = 1.;

        sum = sum + ( coefficient*yValue.at(n) );

        loop_number++;
    }

    double integral = ( ( (1.0*(xStop - xStart) ) / (1.0*nPoints) ) / 3. ) * sum;

    return integral;
}

double getSX3Angle ( int upstream, double pos )
{
    double angle;

    //NEW
    //if(upstream==1) angle = 180. - (180./3.14159)* atan(98.8/(8.2 + pos*75.)); //nominal offset is ~1.7mm
	if(upstream==1) angle = 180. - (180./3.14159)* atan(98.8/(8.2 - pos*75.)); //nominal offset is ~1.7mm
    if(upstream==0) angle = ((180./3.14159)* atan(98./(pos*75.0 -3.0)));

    //OLD
    // if(upstream==1) angle = 180. - (180./3.14159)* atan(98.8/(6.7 + pos*75.)); //nominal offset is ~1.7mm
    // if(upstream==0) angle = ((180./3.14159)* atan(98./(pos*75.0 -3.5)));



    return angle;
}



double proton_energy_loss(double initial_energy, double material_length)
{
    //double E_p[38] = {0.005,0.01,0.02,0.03,0.04,0.05,0.1,0.15,0.2,0.25,0.5,
   // 0.75,1,1.25,1.5,1.75,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
    //19,20};

    //double E_loss[38] = {5,5.1,6.23,6.96,7.5,7.89,8.44,7.93,7.17,6.41,4.15,3.19,
    //2.61,2.21,1.95,1.76,1.61,1.37,1.18,1.05,0.953,0.871,0.803,0.696,0.615,0.552,
    //0.501,0.46,0.423,0.393,0.367,0.345,0.326,0.309,0.293,0.28,0.267,0.256};

	double E_p[36] = {0.005,0.01,0.02,0.03,0.04,0.05,0.1,0.15,0.2,0.25,0.5,
	0.7,1,1.2,1.5,1.7,2,2.5,3,3.5,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20};

	double E_loss[36] = { 0.7529, 1.088, 1.601, 1.922, 2.13, 2.263, 2.412, 2.287, 
	2.122, 1.966, 1.437, 1.2, 0.9787, 0.9037, 0.871, 0.7552, 0.6903, 0.6066, 0.5436, 0.4942, 
	0.4543, 0.3933, 0.3486, 0.3141, 0.2866, 0.2642, 0.2454, 0.2294, 0.2156, 0.2036, 0.193, 0.1836, 0.1752, 0.1677, 0.1608, 0.1488};

    CubicSpline E_loss_spline(E_p, E_loss);
    double updated_energy = initial_energy;

    //Material length is in um!!!

    //incrementing through the material in 0.01um steps
    for(double dx = 0; dx < material_length; dx=dx+0.01)
    {
        updated_energy = updated_energy - (E_loss_spline(updated_energy)/1000.0);
    }

    return updated_energy;

}

double initial_proton_energy(double detected_energy, double material_length)
{
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

    double E_loss[36] = { 0.1506, 0.2511, 0.3740, 0.4262, 0.46971093, 0.5001, 0.5157, 0.4608, 0.4064, 0.3619, 0.2376, 0.1900, 0.1489, 
    0.1348, 0.1146, 0.1044, 0.4064, 0.0779, 0.0676, 0.0599, 0.0539, 0.0452, 0.0390, 0.0345, 0.0309, 0.0281, 0.0258, 0.0239, 0.0222, 
    0.0208, 0.0196, 0.0185, 0.0176,0.0167, 0.0160, 0.0146};

    //cout << "test" <<endl;

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


vector<double>  hit_position_3D(string SX3_or_QQQ5, int upstream, int det, int strip, double pos_or_segment)
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
			z = pos_or_segment*75. - 3.3;  // was -8.2, -5.6
        }
        if(upstream == 0)
        {
            z = z = pos_or_segment*75.0 + 3.3;
            //(pos_or_segment*75.0 - 3.0); // was -3.0 // + 9.6
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

        double rDetCenterD; //98.8 is for the 1-layer barrel

        if (upstream) { rDetCenterD = 98.8; }
        else { rDetCenterD = 100.2; }

        //How far from the center of the detector was the hit? (in the x-y plane)
        double strip_adj_mag = -( ( strip * 10.0 ) - 15.0 ); //mm

        //Angle from the 12 oclock detector to the hit detector
        double det_angle = 0.524 * det; //rads

        //angle the FACE of the detector makes in the x-y plane
        double perp_det_angle = det_angle + TMath::Pi()/2.0; //rads

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
			z = pos_or_segment*75. - 3.3; // my code already gives negative positions for upstream (8.2) was -5.6 for CD2
        }
        if(upstream == 0)
        {
            z = pos_or_segment*75.0 + 3.3; // was +3.3 // 9.6 for the CD2
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
        //double radius[32] = {2.6475,2.9,3.1475,3.39,3.6275,3.86,4.0875,4.31,
        //4.5275,4.74,4.9475,5.15,5.3475,5.54,5.7275,5.919,6.0875,6.26,6.4275,
        //6.59,6.7475,6.9,7.0475,7.19,7.3275,7.46,7.5875,7.71,7.8275,7.94,8.0475,8.15};

        // From Micron
        //double radius[32] = {26.136, 27.908, 29.680, 31.452, 33.223, 34.995, 36.767, 
        //38.539, 40.311, 42.083, 43.855, 45.627, 47.398, 49.170, 50.942, 52.714,
        //54.486, 56.258, 58.030, 59.802, 61.573, 63.345, 65.117, 66.889, 68.661,
        //70.433, 72.205, 73.977, 75.748, 77.520, 79.292, 81.064};

        // From NIM
        double radius[32] = {26.5, 29.0, 31.5, 33.9, 36.3, 38.6, 40.9, 43.1, 45.3, 47.4, 49.5, 51.5, 53.5, 55.4, 57.3, 59.1, 
        60.9, 62.6, 64.3, 65.9, 67.5, 69.0, 70.5, 71.9, 73.3, 74.6, 75.9, 77.1, 78.3, 79.4, 80.5, 81.5};


        double theta_detCenter = 45. + (90. * det);
        double theta_segmentCenter = (theta_detCenter - 30.) + (pos_or_segment*20.);

        x = sin(theta_segmentCenter*(TMath::Pi()/180.))*radius[strip];
        y = cos(theta_segmentCenter*(TMath::Pi()/180.))*radius[strip];

        if(upstream==1)
        {
            z = -80.; // -86
        }
        if(upstream==0)
        {
            z = 80.; // 63
        }

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

double proton_distance_through_target(vector<double> hit_pos)
{

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


double IC_corrected_angle(double protonX, double protonY, double protonZ, float IC_X, float IC_Y )
{


    // Fixing interaction point and changing the angle
    // ===============================================
    double beam_vector[3];
    double proton_vector[3];

    beam_vector[0] = 372.11;
    beam_vector[1] = IC_Y;
    beam_vector[2] = IC_X;

    proton_vector[0] = protonZ;
    proton_vector[1] = protonX;
    proton_vector[2] = protonY;

    double numerator = (beam_vector[0]*proton_vector[0]) + (beam_vector[1]*proton_vector[1]) + (beam_vector[2]*proton_vector[2]);
    double denominator = sqrt( pow(proton_vector[0],2) + pow(proton_vector[1],2) + pow(proton_vector[2],2)) * sqrt(pow(beam_vector[0],2) + pow(beam_vector[1],2) + pow(beam_vector[2],2));
    
    double angle = (180./3.14159) * acos( (numerator*1.)/(denominator*1.) );

    return angle;


    // Fixing the anlge and changing the interaction point
    // ===================================================


    // double beam_vector[3];
    // double proton_vector[3];

    // cout << IC_X << "   " << IC_Y << endl;

    // beam_vector[0] = 372.11;
    // beam_vector[1] = 0;
    // beam_vector[2] = 0;

    // IC_Y = IC_Y*0.0;
    // IC_X = IC_X*0.0;

    // proton_vector[0] = protonZ - (IC_Y*0.5095);
    // proton_vector[1] = protonX - IC_Y;
    // proton_vector[2] = protonY - IC_X;

    // // cout << endl;
    // // cout << protonZ << "  " << protonZ - (IC_Y*0.5095) << endl;
    // // cout << protonX << "  " << protonX - IC_X << endl;
    // // cout << protonY << "  " << protonY - IC_Y << endl;

    // double numerator = (beam_vector[0]*proton_vector[0]) + (beam_vector[1]*proton_vector[1]) + (beam_vector[2]*proton_vector[2]);
    // double denominator = sqrt( pow(proton_vector[0],2) + pow(proton_vector[1],2) + pow(proton_vector[2],2)) * sqrt(pow(beam_vector[0],2) + pow(beam_vector[1],2) + pow(beam_vector[2],2));
    
    // double angle = (180./3.14159) * acos( (numerator*1.)/(denominator*1.) );

    // return angle;
}

double get_solid_angle(string det_group, double bin_center, double bin_width, double QQQ5_z_offset)
{
    double liveStripsInner, liveStripsOuter;
    double Ephi;
    double solid_angle;



    double bin_begin = bin_center - (bin_width/2.);
    double bin_end = bin_center + (bin_width/2.);

    //For a radius of 98.8mm, inner strips have a phi coverage of 5.779 degrees
    //and outer strips have a phi coverage of 5.665 degrees.

    if (det_group == "Upstream_SX3")
    {
        liveStripsInner = 22;
        liveStripsOuter = 22;

        Ephi = ((liveStripsInner*5.779) + (liveStripsOuter*5.665)) / 360.;

        solid_angle = Ephi * 2 * 3.14159 * (cos(bin_begin*(3.14159/180.)) - cos(bin_end*(3.14159/180.)) );
    }
    if (det_group == "Downstream_SX3")
    {
        liveStripsInner = 16;
        liveStripsOuter = 16;

        Ephi = ((liveStripsInner*5.779) + (liveStripsOuter*5.665)) / 360.;

        solid_angle = Ephi * 2 * 3.14159 * (cos(bin_begin*(3.14159/180.)) - cos(bin_end*(3.14159/180.)) );
    }

    
    if (det_group == "Upstream_QQQ5")
    {
            double ring_radius[32] = {2.6475,2.9,3.1475,3.39,3.6275,3.86,
                4.0875,4.31,4.5275,4.74,4.9475,5.15,5.3475,5.54,5.7275,5.919,
                6.0875,6.26,6.4275,6.59,6.7475,6.9,7.0475,7.19,7.3275,7.46,
                7.5875,7.71,7.8275,7.94,8.0475,8.15}; //cm
            double ring_inner_radius[32] = {2.53, 2.78, 3.03, 3.28, 3.52, 3.75, 3.98,
                4.21, 4.43, 4.64, 4.85, 5.06, 5.26, 5.45, 5.64, 5.83, 6.01, 6.18, 6.35, 
                6.52, 6.68, 6.83, 6.98, 7.13, 7.27, 7.40, 7.53, 7.66, 7.78, 7.89, 8.00, 8.11}; //cm
            double ring_outer_radius[32] = {2.77, 3.02, 3.27, 3.51, 3.74, 3.97, 4.20, 4.42, 4.63,
                4.84, 5.05, 5.25, 5.44, 5.63, 5.82, 6.00, 6.17, 6.34, 6.51, 6.67, 6.82, 6.97, 7.12, 7.26, 7.39, 7.52, 
                7.65, 7.77, 7.88, 7.99, 8.10, 8.20}; //cm
            double QQQ5_Ephi[32] = {0.749,0.771,0.789,0.805,0.818,0.829,0.838,0.847,0.854,
                0.861,0.866,0.872,0.877,0.881,0.885,0.888,0.892,0.895,0.897,0.9,0.902,0.904,0.906,0.908,
                0.91,0.912,0.913,0.914,0.916,0.917,0.918,0.919};
            double angleStart, angleStop;
            double strip_SA;


        //NOTE!!!! For QQQ5's, the bin_center variable acutally represents the starting strip
        // and bin_width represents the number of strips in the bin.
        int start_strip = round(bin_center);
        int number_of_strips = round(bin_width);


        for(int strip = start_strip; strip< start_strip+number_of_strips; strip++)
        {


            angleStart = 180.0 - (atan(ring_inner_radius[strip]/QQQ5_z_offset))*(180./3.14159); //degrees
            angleStop = 180.0 - (atan(ring_outer_radius[strip]/QQQ5_z_offset))*(180./3.14159); //degrees



            strip_SA = QQQ5_Ephi[strip] * 2 * 3.14159 * (cos(angleStop*(3.14159/180.)) - cos(angleStart*(3.14159/180.)) );

            solid_angle = solid_angle + strip_SA;
        }

    }
    return solid_angle;
}


double get_solid_angle_per_det(string det_group, double bin_center, double bin_width, double QQQ5_z_offset)
{
    double liveStripsInner, liveStripsOuter;
    double Ephi;
    double solid_angle;



    double bin_begin = bin_center - (bin_width/2.);
    double bin_end = bin_center + (bin_width/2.);

    //For a radius of 98.8mm, inner strips have a phi coverage of 5.779 degrees
    //and outer strips have a phi coverage of 5.665 degrees.

    if (det_group == "Upstream_SX3")
    {
        liveStripsInner = 2;
        liveStripsOuter = 2;

        Ephi = ((liveStripsInner*5.779) + (liveStripsOuter*5.665)) / 360.;

        solid_angle = Ephi * 2 * 3.14159 * (cos(bin_begin*(3.14159/180.)) - cos(bin_end*(3.14159/180.)) );
    }
    if (det_group == "Downstream_SX3")
    {
        liveStripsInner = 2;
        liveStripsOuter = 2;

        Ephi = ((liveStripsInner*5.779) + (liveStripsOuter*5.665)) / 360.;

        solid_angle = Ephi * 2 * 3.14159 * (cos(bin_begin*(3.14159/180.)) - cos(bin_end*(3.14159/180.)) );
    }

    
    if (det_group == "Upstream_QQQ5")
    {
            double ring_radius[32] = {2.6475,2.9,3.1475,3.39,3.6275,3.86,
                4.0875,4.31,4.5275,4.74,4.9475,5.15,5.3475,5.54,5.7275,5.919,
                6.0875,6.26,6.4275,6.59,6.7475,6.9,7.0475,7.19,7.3275,7.46,
                7.5875,7.71,7.8275,7.94,8.0475,8.15}; //cm
            double ring_inner_radius[32] = {2.53, 2.78, 3.03, 3.28, 3.52, 3.75, 3.98,
                4.21, 4.43, 4.64, 4.85, 5.06, 5.26, 5.45, 5.64, 5.83, 6.01, 6.18, 6.35, 
                6.52, 6.68, 6.83, 6.98, 7.13, 7.27, 7.40, 7.53, 7.66, 7.78, 7.89, 8.00, 8.11}; //cm
            double ring_outer_radius[32] = {2.77, 3.02, 3.27, 3.51, 3.74, 3.97, 4.20, 4.42, 4.63,
                4.84, 5.05, 5.25, 5.44, 5.63, 5.82, 6.00, 6.17, 6.34, 6.51, 6.67, 6.82, 6.97, 7.12, 7.26, 7.39, 7.52, 
                7.65, 7.77, 7.88, 7.99, 8.10, 8.20}; //cm
            double QQQ5_Ephi[32] = {0.749,0.771,0.789,0.805,0.818,0.829,0.838,0.847,0.854,
                0.861,0.866,0.872,0.877,0.881,0.885,0.888,0.892,0.895,0.897,0.9,0.902,0.904,0.906,0.908,
                0.91,0.912,0.913,0.914,0.916,0.917,0.918,0.919};
            double angleStart, angleStop;
            double strip_SA;


        //NOTE!!!! For QQQ5's, the bin_center variable actually represents the starting strip
        // and bin_width represents the number of strips in the bin.
        int start_strip = round(bin_center);
        int number_of_strips = round(bin_width);


        for(int strip = start_strip; strip< start_strip+number_of_strips; strip++)
        {


            angleStart = 180.0 - (atan(ring_inner_radius[strip]/QQQ5_z_offset))*(180./3.14159); //degrees
            angleStop = 180.0 - (atan(ring_outer_radius[strip]/QQQ5_z_offset))*(180./3.14159); //degrees



            strip_SA = QQQ5_Ephi[strip] * 2 * 3.14159 * (cos(angleStop*(3.14159/180.)) - cos(angleStart*(3.14159/180.)) );

            solid_angle = solid_angle + strip_SA;
        }
        solid_angle = solid_angle/4.;

    }

    return solid_angle;

}

double get_transverse_momentum(double angle, double energy)
{
    // subtract 90 deg from the angle to get it relative to he transverse plane
    double transverseAngle = abs(angle-90.);

    double transverseMom =  cos( (3.14159/180.)*transverseAngle ) * sqrt(2. * 1. * energy );

    return transverseMom;

}

double calculate_recoil_angle(double angle, double energy, string particle)
{
    // Here need to take the energy and angle of the proton/deuteron and 
    // get the transverse momentum. Then need to use this to calculate the 
    // recoil angle.
    double recoilMom_Z = 10324.; //units are sqrt(AMU.keV) 

    // subtract 90 deg from the angle to get it relative to he transverse plane
    double transverseAngle = abs(angle-90.);
    double transverseMom;

    if(particle == "proton") transverseMom =  cos( (3.14159/180.)*transverseAngle ) * sqrt(2. * 1. * energy );
    if(particle == "deuteron") transverseMom =  cos( (3.14159/180.)*transverseAngle ) * sqrt(2. * 2. * energy );

    double recoilAngle = (180.0/3.14159) * atan(transverseMom / recoilMom_Z);

    return recoilAngle;

}
