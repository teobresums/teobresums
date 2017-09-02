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

#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <vector>

#include "TEOBResum.h"

using namespace::std;

vector<gsl_complex> hlmNQC(double nu, double r, double prstar, double  Omega, double ddotr)
{

    /** This file computes the NQC corrections to the RWZ multipolar waveform. */
    /*
     EOBNQCabFit Coefficients for NQC as fitted function of nu.
     
     [a,b] = EOBInitNQCFit(nu, EOBopt )
     
     Reference(s)
     Nagar, Damour, Reisswig, Pollney
     http://arxiv.org/abs/1506.08457
     */
    
    vector<double> n(6);
    const int kmax = 35;
    double a1;
    double a2;
    double a3;
    double b1;
    double b2;
    double b3;
    
    /** FITS: possibly to be improved further. Current fits: 9/02/2016 */
    const double xnu  = 1-4*nu;
    const double xnu2 = (1-4*nu)*(1-4*nu);

    /** NQC multipolar correction factor */
    vector<gsl_complex> psilmnqc(kmax);

    for (int i=kmax;i--;)
    {
        switch (i)
        {
            /** l=2 -------------------------------------------------------------------
              *(2,1) */
            case 0:
                
                a1 = 0.0162387198*(7.32653082*xnu2 + 1.19616248*xnu + 0.73496656);
                a2 = -1.80492460*xnu2 + 1.78172686*xnu + 0.30865284;
                a3 = 0.0;
                
                b1 =  -0.0647955017*(3.59934444*xnu2 - 4.08628784*xnu + 1.37890907);
                b2 =  1.3410693180*(0.38491989*xnu2 + 0.10969453*xnu + 0.97513971);
                b3 =  0.0;
                
                n  = NQC(r,prstar, Omega,ddotr,0);
                
                break;
                
            /* (2,2) */
            case 1:
                
                a1 = -0.0805236959*( 1 - 2.00332326*xnu2)/( 1 + 3.08595088*xnu2);
                a2 =  1.5299534255*( 1 + 1.16438929*xnu2)/( 1 + 1.92033923*xnu2);
                a3 =  0.0;
                
                b1 = 0.146768094955*( 0.07417121*xnu + 1.01691256);
                b2 = 0.896911234248*(-0.61072011*xnu + 0.94295129);
                b3 = 0.0;
                
                n  = NQC(r,prstar, Omega,ddotr,1);
                
                break;
                
            /** l=3 -------------------------------------------------------------------
              * (3,3) */
            case 4:
                
                a1 = -0.0377680000*(1 - 14.61548907*xnu2)/( 1 + 2.44559263*xnu2);
                a2 =  1.9898000000*(1 + 2.09750346 *xnu2)/( 1 + 2.57489466*xnu2);
                a3 =  0.0;
                
                b1 = 0.1418400000*(1.07430512 - 1.23906804*xnu + 4.44910652*xnu2);
                b2 = 0.6191300000*(0.80672432 + 4.07432829*xnu - 7.47270977*xnu2);
                b3 = 0.0;
                
                n  = NQC(r,prstar, Omega,ddotr,4);
                
                break;
                
            default:
                
                a1 = 0.;
                a2 = 0.;
                a3 = 0.;
                
                b1 = 0.;
                b2 = 0.;
                b3 = 0.;
                break;
        }
        
        psilmnqc[i].dat[0] = 1. + a1*n[0] + a2*n[1] + a3*n[2];
        psilmnqc[i].dat[1] = 0. + b1*n[3] + b2*n[4] + b3*n[5];
    }
        
    return psilmnqc;
}
