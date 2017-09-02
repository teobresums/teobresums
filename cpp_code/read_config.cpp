/**
 * Copyright (C) 2017 Alessandro Nagar, Gregorio Carullo, Ka Wa Tsang, Philipp Fleig, Sebastiano Bernuzzi, Walter Del Pozzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#include <iostream>
#include <list>
#include <fstream>
#include <stdio.h>
#include <string>
#include <math.h>

#include "TEOBResum.h"

using namespace::std;

/** Sets the dynamics controlling flags to their default value */
void SetDefaultFlagsValues(TEOBResumParams *p)
{
    p->flags.NQC        = 1;
    p->flags.spin       = 1;
    p->flags.tidal      = 0;
    p->flags.RWZ        = 0;
    p->flags.speedy     = 1;
    p->flags.dynamics   = 0;
    p->flags.waveform   = 0;
}

double logQ(double x)
{
    /** logQ-vs-log(lambda) fit of Table I of Yunes-Yagi
        here x = log(lambda) and the output is the log of the coefficient
        that describes the quadrupole deformation due to spin. */
    double ai = 0.194;
    double bi = 0.0936;
    double ci = 0.0474;
    double di = -4.21e-3;
    double ei = 1.23e-4;
    double x2 = x*x;
    double x3 = x*x2;
    double x4 = x*x3;
    
    return ai + bi*x + ci*x2 + di*x3 + ei*x4;
}

double Yagi13_fit_barlamdel(double barlam2, int ell)
{
    /*
     Yagi 2013 fits for multipolar
     $\bar{\lambda}_\ell$ = 2 k_\ell/(C^{2\ell+1} (2\ell-1)!!)$
     Eq.(10),(61); Tab.I; Fig.8 http://arxiv.org/abs/1311.0872
     */
    double lnx = log(barlam2);
    double coeffs[5] = {0.0};
    if (ell == 3)
    {
        coeffs[0] = 2.52e-5;
        coeffs[1] = -1.31e-3;
        coeffs[2] = 2.51e-2;
        coeffs[3] = 1.18;
        coeffs[4] = -1.15;
    }
    else if (ell == 4)
    {
        coeffs[0] = 2.8e-5;
        coeffs[1] =-1.81e-3;
        coeffs[2] =3.95e-2;
        coeffs[3] =1.43;
        coeffs[4] =-2.45;
    }
    else return 0.0;
    
    double lny = coeffs[0]*lnx*lnx*lnx*lnx+coeffs[1]*lnx*lnx*lnx+coeffs[2]*lnx*lnx+coeffs[3]*lnx+coeffs[4];
    return exp(lny);
}

double radius0(double M, double f_start)
{
    double MSUN_SEC = 4.925490949141889e-06;
    double x = (M*f_start*MSUN_SEC*2.*M_PI)/2.;
    return cbrt(1/(x*x));
}

double time_units_conversion(double M, double Srate)
{
    double MSUN_SEC = 4.925490949141889e-06;
    return (1./Srate)/(M*MSUN_SEC);
}

TEOBResumParams read_config(char *fname)
{
    TEOBResumParams params;
 
    SetDefaultFlagsValues(&params);
    
    string param_name;
    double param_value;

    ifstream fin (fname);
    if (!fin) //checks to see if file opens properly
    {
        cerr << "error: Could not find the parameters file.";
    }

    int i = 0;
    
    while ( fin >> param_name >> param_value )
    {
        //param_values.push_back(param_value);
        cout << param_name <<"\t"<< param_value << endl;
        switch (i) {
        case 0:
            params.q    = param_value;
            break;
        case 1:
            params.chi1 = param_value;
            break;
        case 2:
            params.chi2 = param_value;
            break;
        case 3:
            params.r0   = param_value;
            break;
        case 4:
            if (param_value==0) {
                params.flags.NQC = 0;
            } else {
                params.flags.NQC = 1;
            }
            break;
        case 5:
            if (param_value==0) {
                params.flags.tidal = 0;
            } else {
                params.flags.tidal = 1;
            }
            break;
        case 6:
            if (param_value==0) {
                params.flags.RWZ = 0;
            } else {
                params.flags.RWZ = 1;
            }
            break;
        case 7:
            if (param_value==0) {
                params.flags.speedy = 0;
            } else {
                params.flags.speedy = 1;
            }
            break;
        case 8:
            if (param_value==0) {
                params.flags.dynamics = 0;
            } else {
                params.flags.dynamics = 1;
            }
            break;
        case 9:
            if (param_value==0) {
                params.flags.waveform = 0;
            } else {
                params.flags.waveform = 1;
            }
            break;
        case 10:
            params.lm = param_value;
            break;
        case 11:
            params.dt = param_value;
            break;
        case 12:
            params.solver_scheme = param_value;
            break;
        case 13:
            params.LambdaAl2 = param_value;
            break;
        case 14:
            params.LambdaBl2 = param_value;
            break;
        default:
            break;
        }
        
        i++;
    }
    fin.close();

    // calculate parameters

    double chi1 = params.chi1;
    double chi2 = params.chi2;
    double q = params.q;
    
    params.LambdaAl3 = Yagi13_fit_barlamdel(params.LambdaAl2, 3);
    params.LambdaBl3 = Yagi13_fit_barlamdel(params.LambdaBl2, 3);
    params.LambdaAl4 = Yagi13_fit_barlamdel(params.LambdaAl2, 4);
    params.LambdaBl4 = Yagi13_fit_barlamdel(params.LambdaBl2, 4);
    
    /** Override spin settings if spins are given in input */
    if (chi1 != .0 || chi2 != .0) params.flags.spin = 1;

    double nu = q/((q+1.)*(q+1.));
    params.nu = nu;
    
    double X1 = 0.5*(1.+sqrt(1.-4.*nu));
    double X2 = 1. - X1;
    double XA = X1; // a different notation used in tidal part, keep here for simplicity
    double XB = X2;
    params.X1 = X1;
    params.X2 = X2;
    
    double S1 = pow(params.X1, 2.) * params.chi1;
    double S2 = pow(params.X2, 2.) * params.chi2;
    params.S1 = S1;
    params.S2 = S2;
    
    double a1  = X1*chi1;
    double a2  = X2*chi2;
    double aK  = a1 + a2;
    double aK2 = aK*aK;
    
    params.S     = S1 + S2;        // => in the EMRL this becomes the spin of the BH
    params.Sstar = X2*a1 + X1*a2;  // => in the EMRL this becomes the spin of the particle
    
    params.a1  = X1*chi1;
    params.a2  = X2*chi2;
    params.aK  = a1 + a2;
    params.aK2 = aK2;
    params.rLR = 0.;
    params.cN3LO = c3_fit_global(nu,chi1,chi2,X1,X2,a1,a2,params.flags.tidal);
    
    /** if tidal corrections are requested, unset the NQC corrections */
    if (params.flags.tidal==1) params.flags.NQC = 0;

    // tidal params
    double LambdaAl2 = params.LambdaAl2;
    double LambdaAl3 = params.LambdaAl3;
    double LambdaAl4 = params.LambdaAl4;
    
    double LambdaBl2 = params.LambdaBl2;
    double LambdaBl3 = params.LambdaBl3;
    double LambdaBl4 = params.LambdaBl4;
    
    /** Computing the tidal coupling constants */
    
    double kapA2 = 3.   * LambdaAl2 * XA*XA*XA*XA*XA / q; //Note: kap stands for kappa; see eqn(1) of REF
    double kapA3 = 15.  * LambdaAl3 * XA*XA*XA*XA*XA*XA*XA / q;
    double kapA4 = 105. * LambdaAl4 * XA*XA*XA*XA*XA*XA*XA*XA*XA / q;
    
    double kapB2 = 3.   * LambdaBl2 * XB*XB*XB*XB*XB * q;
    double kapB3 = 15.  * LambdaBl3 * XB*XB*XB*XB*XB*XB*XB * q;
    double kapB4 = 105. * LambdaBl4 * XB*XB*XB*XB*XB*XB*XB*XB*XB * q;
    
    double kapT2 = kapA2 + kapB2;
    double kapT3 = kapA3 + kapB3;
    double kapT4 = kapA4 + kapB4;

    params.kappaAl2 = kapA2;
    params.kappaAl3 = kapA3;
    params.kappaAl4 = kapA4;

    params.kappaBl2 = kapB2;
    params.kappaBl3 = kapB3;
    params.kappaBl4 = kapB4;

    params.kappaTl2 = kapT2;
    params.kappaTl3 = kapT3;
    params.kappaTl4 = kapT4;

    //-----------------------------------------------------------------------------------
    // Definition of the conservative tidal coefficients \bar{\alpha}_n^{(\ell)}, Eq.(37)
    // of Damour&Nagar, PRD 81, 084016 (2010)
    //-----------------------------------------------------------------------------------
    params.bar_alph2_1 = (5./2.*XA*kapA2 + 5./2.*XB*kapB2)/kapT2;
    params.bar_alph2_2 = ((3.+XA/8.+ 337./28.*XA*XA)*kapA2 + (3.+XB/8.+ 337./28.*XB*XB)*kapB2)/kapT2; 
    params.bar_alph3_1 = ((-2.+15./2.*XA)*kapA3 + (-2.+15./2.*XB)*kapB3)/kapT3;			     			   
    params.bar_alph3_2 = ((8./3.-311./24.*XA+110./3.*XA*XA)*kapA3 + (8./3.-311./24.*XB+110./3.*XB*XB)*kapB3)/kapT3;
    

    //double lambda1 = params.LambdaAl2;
    //double lambda2 = params.LambdaBl2;
    double logC_Q1 = logQ(log(LambdaAl2));
    double logC_Q2 = logQ(log(LambdaBl2));
    double C_Q1    = exp(logC_Q1);
    double C_Q2    = exp(logC_Q2);    
    params.C_Q1 = C_Q1;
    params.C_Q2 = C_Q2;
    
    return params;
}

TEOBResumParams process_input_parameters(
                                 double m1,
                                 double m2,
                                 double chi1,
                                 double chi2,
                                 double f_min,
                                 double sampling_rate,
                                 double LambdaAl2,
                                 double LambdaBl2,
                                 int    NQC,
                                 int    tidal,
                                 int    speedy,
                                 int    RWZ,
                                 int    dynamics,
                                 int    waveform,
                                 int    lm,
                                 int    solver_scheme
                                 )
{
    TEOBResumParams params;
    
    SetDefaultFlagsValues(&params);
    
    double mtot = m1+m2;
    double q = m1/m2;
    params.mtot = mtot;
    params.q = q;
    params.chi1 = chi1;
    params.chi2 = chi2;
    params.r0 = radius0(mtot, f_min);
    params.flags.NQC = NQC;
    params.flags.tidal = tidal;
    params.flags.RWZ = RWZ;
    params.flags.speedy = speedy;
    params.flags.dynamics = dynamics;
    params.flags.waveform = waveform;
    params.lm = lm;
    params.dt = time_units_conversion(mtot, sampling_rate);
    params.solver_scheme = solver_scheme;
    
    if (params.flags.tidal == 1)
    {
        params.LambdaAl2 = LambdaAl2;
        params.LambdaBl2 = LambdaBl2;
        params.LambdaAl3 = Yagi13_fit_barlamdel(LambdaAl2, 3);
        params.LambdaBl3 = Yagi13_fit_barlamdel(LambdaBl2, 3);
        params.LambdaAl4 = Yagi13_fit_barlamdel(LambdaAl2, 4);
        params.LambdaBl4 = Yagi13_fit_barlamdel(LambdaBl2, 4);
    }

    /** Override spin settings if spins are given in input */
    if (chi1 != .0 || chi2 != .0) params.flags.spin = 1;

    double nu = q/((q+1.)*(q+1.));
    params.nu = nu;
    
    double X1 = 0.5*(1.+sqrt(1.-4.*nu));
    double X2 = 1. - X1;
    double XA = X1; // a different notation used in tidal part, keep here for simplicity
    double XB = X2;
    params.X1   = X1;
    params.X2   = X2;
    
    double S1 = params.X1*params.X1 * params.chi1;
    double S2 = params.X2*params.X2 * params.chi2;
    params.S1 = S1;
    params.S2 = S2;
    
    double a1  = X1*chi1;
    double a2  = X2*chi2;
    double aK  = a1 + a2;
    double aK2 = aK*aK;

    params.S     = S1 + S2;        // => in the EMRL this becomes the spin of the BH
    params.Sstar = X2*a1 + X1*a2;  // => in the EMRL this becomes the spin of the particle
    
    params.a1  = X1*chi1;
    params.a2  = X2*chi2;
    params.aK  = a1 + a2;
    params.aK2 = aK2;
    
    params.rLR = 0.;

    params.cN3LO = c3_fit_global(nu,chi1,chi2,X1,X2,a1,a2,params.flags.tidal);

    /** if tidal corrections are requested, unset the NQC corrections */
    if (params.flags.tidal==1) params.flags.NQC = 0;

    // tidal params
    double LambdaAl3 = params.LambdaAl3;
    double LambdaAl4 = params.LambdaAl4;
    
    double LambdaBl3 = params.LambdaBl3;
    double LambdaBl4 = params.LambdaBl4;
    
    /** Computing the tidal coupling constants */

    double kapA2 = 3.   * LambdaAl2 * XA*XA*XA*XA*XA / q; //Note: kap stands for kappa; see eqn(1) of REF
    double kapA3 = 15.  * LambdaAl3 * XA*XA*XA*XA*XA*XA*XA / q;
    double kapA4 = 105. * LambdaAl4 * XA*XA*XA*XA*XA*XA*XA*XA*XA / q;
    
    double kapB2 = 3.   * LambdaBl2 * XB*XB*XB*XB*XB * q;
    double kapB3 = 15.  * LambdaBl3 * XB*XB*XB*XB*XB*XB*XB * q;
    double kapB4 = 105. * LambdaBl4 * XB*XB*XB*XB*XB*XB*XB*XB*XB * q;
    
    double kapT2 = kapA2 + kapB2;
    double kapT3 = kapA3 + kapB3;
    double kapT4 = kapA4 + kapB4;

    params.kappaAl2 = kapA2;
    params.kappaAl3 = kapA3;
    params.kappaAl4 = kapA4;

    params.kappaBl2 = kapB2;
    params.kappaBl3 = kapB3;
    params.kappaBl4 = kapB4;

    params.kappaTl2 = kapT2;
    params.kappaTl3 = kapT3;
    params.kappaTl4 = kapT4;

    //-----------------------------------------------------------------------------------
    // Definition of the conservative tidal coefficients \bar{\alpha}_n^{(\ell)}, Eq.(37)
    // of Damour&Nagar, PRD 81, 084016 (2010)
    //-----------------------------------------------------------------------------------
    params.bar_alph2_1 = (5./2.*XA*kapA2 + 5./2.*XB*kapB2)/kapT2;
    params.bar_alph2_2 = ((3.+XA/8.+ 337./28.*XA*XA)*kapA2 + (3.+XB/8.+ 337./28.*XB*XB)*kapB2)/kapT2; 
    params.bar_alph3_1 = ((-2.+15./2.*XA)*kapA3 + (-2.+15./2.*XB)*kapB3)/kapT3;			     			   
    params.bar_alph3_2 = ((8./3.-311./24.*XA+110./3.*XA*XA)*kapA3 + (8./3.-311./24.*XB+110./3.*XB*XB)*kapB3)/kapT3;
    
    double logC_Q1 = logQ(log(LambdaAl2));
    double logC_Q2 = logQ(log(LambdaBl2));
    double C_Q1    = exp(logC_Q1);
    double C_Q2    = exp(logC_Q2);    
    params.C_Q1 = C_Q1;
    params.C_Q2 = C_Q2;
    
    return params;
}
