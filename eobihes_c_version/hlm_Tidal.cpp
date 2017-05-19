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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#include <gsl/gsl_math.h>
#include <vector>
#include "hlm_Tidal.h"

using namespace::std;

vector<double> hlm_Tidal(double x,double nu)
{

/** Calculate tidal correction to multipolar waveform.
 *
 *   Reference(s)
 *   Damour, Nagar & Villain, Phys.Rev. D85 (2012) 123007
 */

    int kmax  = 35;
    int x5    = x*x*x*x*x
    double CA = 0.1738106852;
    double CB = 0.1738106852;
    double XA = 0.5;
    double XB = 0.5;
    double khatA_2;
    double khatB_2;
    
    vector<double> kAl(3);
    vector<double> kBl(3);
    vector<double> hTidallm(kmax);
    vector<double> hA(kmax);
    vector<double> hB(kmax);
    vector<double> betaA1(kmax);
    vector<double> betaB1(kmax);

    kAl[0] = 0.093330885635;
    kAl[1] = 0.025545679620;
    kAl[2] = 0.009495642804;
    kBl[0] = 0.093330885635;
    kBl[1] = 0.025545679620;
    kBl[2] = 0.009495642804;
    
    khatA_2 = kAl[0] * XB/XA * gsl_pow_int(XA/CA,5);
    khatB_2 = kBl[0] * XA/XB * gsl_pow_int(XB/CB,5);
    
    for (int i=kmax; i--;)
    {
        hA[i]=0.;
    }
    
    hB=hA;
    betaA1=hA;
    betaB1=hA;
    
    /** l=2 ------------------------------------------------------------------ */
    
    hA[1]     = 2 * khatA_2 *(XA/XB+3);
    hB[1]     = 2 * khatB_2 *(XB/XA+3);
    betaA1[1] = (-202. + 560*XA - 340*XA*XA + 45*XA*XA*XA)/(42*(3-2*XA));
    betaB1[1] = (-202. + 560*XB - 340*XB*XB + 45*XB*XB*XB)/(42*(3-2*XB));
    hA[0]     = 3 * khatA_2 * XB * (3-4*XA)/XA;
    hB[0]     = 3 * khatB_2 * XA * (3-4*XB)/XB;
    
    
    /** l=3 ------------------------------------------------------------------ */
    
    hA[4] = 12 * khatA_2 * XB*XB/XA;
    hB[4] = 12 * khatB_2 * XA*XA/XB;
    hA[2] = hA[4];
    hB[2] = hB[4];
    

    /** l=2 ------------------------------------------------------------------
     
     *  (2,2) */
    hTidallm[1] = ( hA[1]*(1. + betaA1[1]*x) + hB[1]*(1. + betaB1[1]*x) )*x5;
    /** (2,1) */
    hTidallm[0] = ( hA[0] - hB[0] )*x5;


    /** l=3 ------------------------------------------------------------------
     
     *  (3,3) */
    hTidallm[4] = ( hA[4] - hB[4] )*x5;
    /** (3,1) */
    hTidallm[2] = ( hA[2] - hB[2] )*x5;

    return hTidallm;

}
