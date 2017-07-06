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
#include "input_struc.h"

using namespace::std;

vector<double> hlm_Tidal(double x,void *params)
{

    //EOBhlmTidal Calculate tidal correction to multipolar waveform.
    //   hTidallm = EOBhlmTidal( x, nu, Topt, EOBopt )
    //
    //   Reference(s)
    //    Damour, Nagar & Villain, Phys.Rev. D85 (2012) 123007
    //

    int kmax   = 35;
    double x5    = gsl_pow_int(x,5);

    vector<double> kAl(3);
    vector<double> kBl(3);
    vector<double> hA(kmax);

    kAl[0] = (*(input *)params).kAl1;
    kAl[1] = (*(input *)params).kAl2;
    kAl[2] = (*(input *)params).kAl3;
    kBl[0] = (*(input *)params).kBl1;
    kBl[1] = (*(input *)params).kBl2;
    kBl[2] = (*(input *)params).kBl3;
    
    double CA = (*(input *)params).CA;
    double CB = (*(input *)params).CB;
    double XA = (*(input *)params).X1;
    double XB = (*(input *)params).X2;
    
    double khatA_2 = kAl[0] * XB/XA *gsl_pow_int(XA/CA,5);
    double khatB_2 = kBl[0] * XA/XB *gsl_pow_int(XB/CB,5);
    
    for (int i=kmax; i--; )
    {
        hA[i]=0.;
    }
    vector<double> hB=hA;
    vector<double> betaA1=hA;
    //vector<double> betaA2(kmax)=hA;
    vector<double> betaB1=hA;
    //vector<double> betaB2(kmax)=hA;
    vector<double> hTidallm(kmax);

    
    // l=2 ------------------------------------------------------------------
    
    hA[1] = 2 * khatA_2 *(XA/XB+3);
    hB[1] = 2 * khatB_2 *(XB/XA+3);
    
    betaA1[1] = (-202. + 560*XA - 340*XA*XA + 45*XA*XA*XA)/(42*(3-2*XA));
    betaB1[1] = (-202. + 560*XB - 340*XB*XB + 45*XB*XB*XB)/(42*(3-2*XB));
    
    hA[0] = 3 * khatA_2 * XB * (3-4*XA)/XA;
    hB[0] = 3 * khatB_2 * XA * (3-4*XB)/XB;
    
    
    // l=3 ------------------------------------------------------------------
  
    hA[2] = hA[4];
    hB[2] = hB[4];
    
    hA[4] = 12 * khatA_2 * XB*XB/XA;
    hB[4] = 12 * khatB_2 * XA*XA/XB;
    
    /** l=2 ------------------------------------------------------------------
     * (2,1) */
    hTidallm[0] = ( -hA[0] + hB[0] )*x5;
    /* (2,2) */
    hTidallm[1] = ( hA[1]*(1. + betaA1[1]*x) + hB[1]*(1. + betaB1[1]*x) )*x5;

    /** l=3 ------------------------------------------------------------------
     * (3,1) */
    hTidallm[2] = ( -hA[2] + hB[2] )*x5;
    /* (3,3) */
    hTidallm[4] = ( -hA[4] + hB[4] )*x5;

    return hTidallm;

}

//double nu = (*(input *)params).nu;
/*
 kAl[0]      =     0.093330885635; //k.A.2
 kAl[1]      =     0.025545679620; //k.A.3
 kAl[2]      =     0.009495642804; //k.A.4
 
 kBl[0]      =     0.093330885635; //k.B.2
 kBl[1]      =     0.025545679620; //k.B.3
 kBl[2]      =     0.009495642804; //k.B.4
 
 double CA = 0.1738106852;
 double CB = 0.1738106852;
 double XA =0.5;//= MA/M;
 double XB =0.5;// MB/M;
 */
