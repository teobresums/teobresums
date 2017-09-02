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

#include <gsl/gsl_math.h>
#include <cmath>

#include "TEOBResum.h"

double s_Flux(double x, double Omega, double r_omega, double E, double Heff, double jhat, double r, double pr_star, double ddotr, void *params){
/*
% DINFLUX This function computes the Newton.Normalized energy flux according to
%         the DIN resummation procedure. It is also designed so to add non-QC
%         and non-K corrections to  (2,2) partial flux.
%
%         USAGE:
%
%         [Flm F hatF hatF_resum]=DINFlux(x,Omega,E,Heff,jhat,nu,lmax,r,pr_star,ddotr)
%
%         where:
%
%         x       :: PN argument
%         Omega   :: Orbital frequency
%         E       :: Energy
%         Heff    :: Effective energy
%         jhat    :: Newton-Normalized angular momentum
%         nu      :: symmetric mass ratio
%         lmax    :: maximum l
%         r       :: EOB radius
%         pr_star :: radial momentum
%         ddotr   :: \ddot{r}
*/
    
    double nu       = (*(input *)params).nu;
    double chi1     = (*(input *)params).chi1;
    double chi2     = (*(input *)params).chi2;
    double X1       = (*(input *)params).X1;
    double X2       = (*(input *)params).X2;
    bool tidal_flag = (*(input *)params).tidal;
    bool NQC_flag   = (*(input *)params).NQC;
    
    double prefact[] = {
        jhat,Heff,
        Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,
        Heff,jhat,Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,jhat,Heff,
        Heff,jhat,Heff,jhat,Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,jhat,Heff,jhat,Heff};

    
    double hnqclm = 1.;
    double Flm;
    double Modhhatlm;
    double sqrt_one_4nu = sqrt(1.-4.*nu);
    
    vector<double> flm     = s_flm(x,params);
    vector<double> FNewtlm = FlmNewt(x,params);
    double FNewt22         = FNewtlm[1];

    vector<double> MTlm     = Tlm(E*Omega);
    vector<double> hlmTidal = hlm_Tidal(x,params);

    double SFlm = 0.;
    for (int k=35; k--;) {
        // Compute modulus of hhat_lm
        Modhhatlm = prefact[k] * MTlm[k] * flm[k];
        
        //Include NQC with flag
        if (NQC_flag==true) {
            Modhhatlm *= hnqclm;
        }

    //Make tidal corrections
    if (tidal_flag==true) {
        switch (k) {
            case 0: // (2,1)
                Modhhatlm *= sqrt_one_4nu;
                break;
            case 2: // (3,1)
                Modhhatlm *= sqrt_one_4nu;
                break;
            case 4: // (3,3)
                Modhhatlm *= sqrt_one_4nu;
                break;
                
            default: Modhhatlm *= 1.;
                break;
        }
        
        Modhhatlm += MTlm[k]*hlmTidal[k];
        
    }

	
        // Total flux multipoles
        Flm = (Modhhatlm * Modhhatlm) * FNewtlm[k];
        
        SFlm += Flm;
    }
    
    // Sum over multipoles and normalize to the 22 Newtonian multipole
    double hatf = SFlm/(FNewt22);
    
    if (tidal_flag==false) {
        double hatFH = s_HorizonFlux(x, Heff, jhat, nu, X1, X2, chi1, chi2);
        hatf += hatFH;
    }
    
    double Fphi = -32.0/5.0* nu * gsl_pow_int(r_omega,4) * gsl_pow_int(Omega,5) * hatf;

    return Fphi;
}

