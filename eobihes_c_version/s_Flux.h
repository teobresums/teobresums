//
//  EOB_Flux.h
//
//  Created by Philipp Fleig on 25/04/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef EOB_Flux_h
#define EOB_Flux_h

#include "s_flm.h"
#include "FlmNewt.h"
#include "s_HorizonFlux.h"
#include "hlm_Tidal.h"
#include "Tlm.h"
#include "input_struc.h"
#include "multipole_index.h"

double s_Flux(double x,double Omega,double r_omega,double E,double Heff,double jhat,double r,double pr_star,double ddotr,void *params){
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
    
    vector<double> flm = s_flm(x,params);
    vector<double> FNewtlm =FlmNewt(x,params);
    double FNewt22 = FNewtlm[1];

    vector<double> MTlm = Tlm(E*Omega);
    vector<double> hlmTidal = hlm_Tidal(x,nu);

    double SFlm=0.;
    for (int k=35; k--;) {
        // Compute modulus of hhat_lm
        Modhhatlm = prefact[k] * MTlm[k] * flm[k];
        
        //Include NQC with flag
        if (NQC_flag==true) {
            Modhhatlm*=hnqclm;
        }

    //Make tidal corrections
    if (tidal_flag==true) {
        switch (k) {
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

	
        // Total flux multipoles
        Flm = (Modhhatlm * Modhhatlm) * FNewtlm[k];
        
        SFlm += Flm;
    }
    
    // Sum over multipoles and normalize to the 22 Newtonian multipole
    double hatf = SFlm/(FNewt22);
    
    if (tidal_flag==false) {
        double hatFH =s_HorizonFlux(x,Heff,jhat,nu,X1,X2,chi1,chi2);
        hatf += hatFH;
    }
    
    double Fphi = -32.0/5.0* nu * gsl_pow_int(r_omega,4) * gsl_pow_int(Omega,5) * hatf;

    return Fphi;
}


#endif /* EOB_Flux_h */
