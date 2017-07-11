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
#include <stdbool.h>

#include "flux.h"
#include "f_lm.h"
#include "NQC.h"
#include "FlmNewt.h"
#include "Tlm.h"
#include "HorizonFlux.h"
#include "hlmNQC.h"
#include "hlm_Tidal.h"
#include "multipole_index.h"
#include "input_struc.h"

double flux(const double x,const double Omega,const double r_omega,const double E, const double Heff,const double jhat,const double r,const double prstar, const double ddotr, double source[],void *params)
{

    bool tidal_flag = (*(input *)params).tidal;
    bool NQC_flag   = (*(input *)params).NQC;
    double nu       = (*(input *)params).nu;
    double Flm;
    double Modhhatlm;
    double sqrt_one_4nu = sqrt(1.-4.*nu);
    const vector<double> flm     = f_lm(x,nu);
    const vector<double> FNewtlm = FlmNewt(x,params);
    const double FNewt22         = FNewtlm[1];
    vector<double> MTlm          = Tlm(E*Omega);
    vector<double> hlmTidal      = hlm_Tidal(x,params);

    /** Compute NQC correction to the modulus of the (l,m) waveform */
    vector<gsl_complex> hlm_NQC = hlmNQC(nu,r,prstar,Omega,ddotr);
    double SFlm=0.;
    for (int k=35; k--;)
    {
        /** Compute modulus of hhat_lm */
        Modhhatlm = source[k] * MTlm[k] * flm[k];
    
        /** NQCs */
        if (NQC_flag==true)
        {
            switch (k)
            {
                case 0:
                    Modhhatlm *= hlm_NQC[k].dat[0];
                    break;
                case 1:
                    Modhhatlm *= hlm_NQC[k].dat[0];
                    break;
                case 4:
                    Modhhatlm *= hlm_NQC[k].dat[0];
                    break;
                default:
                    Modhhatlm *= 1.;
                    break;
            }
        }
    
        //Make tidal corrections
        if (tidal_flag==true)
        {
            switch (k)
            {
                case 0: // (2,1)
                    Modhhatlm *=sqrt_one_4nu;
                    break;
                case 2: // (3,1)
                    Modhhatlm *=sqrt_one_4nu;
                    break;
                case 4: // (3,3)
                    Modhhatlm *=sqrt_one_4nu;
                    break;
                default: Modhhatlm *= 1.;
                    break;
            }
            
            Modhhatlm += MTlm[k]*hlmTidal[k];
            
        }
        //if (k==1)
        //{
        //    printf(" %s %.8e \n","hlmTidal",hlmTidal[k]);
        //}
        /** Total flux multipoles */
        Flm = (Modhhatlm * Modhhatlm) * FNewtlm[k];
        SFlm += Flm;
    }
    
    /** Sum over multipoles and normalize to the 22 Newtonian multipole */
    double hatf = SFlm/(FNewt22);

    if (tidal_flag==false)
    {
        double hatFH = HorizonFlux(x,Heff,jhat,nu);
        hatf += hatFH;
    }
    
    double Fphi = -32.0/5.0* nu * gsl_pow_int(r_omega,4) * gsl_pow_int(Omega,5) * hatf;

    return Fphi;
}
