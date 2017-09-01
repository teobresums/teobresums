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
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <cmath>

#include "hlmNewt.h"
#include "hhatlmtail.h"
#include "hlmNQC.h"
#include "deltalm.h"
#include "SpeedyTail.h"
#include "s_flm.h"
#include "ringdown_match.h"
#include "hlm_Tidal.h"
#include "multipole_index.h"

#include "hlm.h"
#include "input_struc.h"
#include "f_lm.h"

using namespace::std;

vector<gsl_complex> hlm(double t, const double phi, const double r, const double pph, const double prstar, double Omega, const double ddotr, const double H, const double Heff,const double jhat, const double rw,void *params)
{
    int kmax = 35; 

    vector<gsl_complex> hlm(kmax);
    double nu            = (*(input *)params).nu;
    bool tidal_flag      = (*(input *)params).tidal;
    bool spin_flag       = (*(input *)params).spin;
    bool NQC_flag        = (*(input *)params).NQC;
    bool speedytail_flag = (*(input *)params).speedy;
    
    double source[] = {
        jhat,Heff,
        Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,
        Heff,jhat,Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,jhat,Heff,
        Heff,jhat,Heff,jhat,Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,jhat,Heff,jhat,Heff};
    
    /** Newtonian waveform */
    vector<gsl_complex> hNewt = hlmNewt( rw,Omega,phi, nu,tidal_flag);
    
    /** Compute corrections */
    double x = gsl_pow_int(rw*Omega,2);
    vector<double> flm(35);
    if (spin_flag==true)
    {
        flm = s_flm(x,params);
    }
    else
    {
        flm = f_lm(x,nu);
    }
    
    /** Computing the tail */
    const double r0    = 1.213061319425267e+00;   // 2/sqrt(e);
    const double Hreal = H * nu;
    vector<gsl_complex> tlm(kmax);
    if (speedytail_flag==false)
    {
        tlm = hhatlmtail(Omega,Hreal, r0, L, M);
    }
    else if (speedytail_flag==true)
    {
        tlm = speedyTail(Omega,Hreal, r0, L, M);
    }
    
    /** Residual phase corrections delta_{lm} */
    const vector<double> EOBdeltalm = deltalm(Hreal, Omega, nu);

    vector<gsl_complex> h_NQC(kmax);
    if (NQC_flag==true)
    {
        h_NQC = hlmNQC(nu,r,prstar,Omega,ddotr);
    }
    
    for (int k=35; k--;)
    {
        tlm[k].dat[1] += EOBdeltalm[k];
        
        /** Compute \hat{h}_lm */
        hlm[k].dat[0] =   hNewt[k].dat[0] * flm[k] * source[k] * tlm[k].dat[0] ;
        hlm[k].dat[1] =   hNewt[k].dat[1] + tlm[k].dat[1] ;
        hlm[k].dat[1] = - hlm[k].dat[1]; /** Minus sign by convention */
    
        /** NQC correction */
        if (NQC_flag==true)
        {
            hlm[k].dat[0] *= h_NQC[k].dat[0];
            hlm[k].dat[1] -= h_NQC[k].dat[1];
        }
    }
    
    if (tidal_flag==true)
    {
            
        /** Compute tidal contribution */
        vector<double> hlmtidal = hlm_Tidal(x, params);

        /** Update waveform */
        double p2 = sqrt(1-4*nu);
        for (int k=35; k--;)
        {
            switch (k)
            {
                case 0:
                    hlm[0].dat[0] *= p2;
                    break;
                case 2:
                    hlm[2].dat[0] *= p2;
                    break;
                case 4:
                    hlm[4].dat[0] *= p2;
                    break;
                default:
                    break;
            }

            double rad_tid = hNewt[k].dat[0] * tlm[k].dat[0] * hlmtidal[k] ;
            hlm[k].dat[0] += rad_tid;
        }
    }
    
    return hlm;
}
