//
//  hlm_hybrid.h
//
//  Created by Philipp Fleig on 02/07/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef hlm_hybrid_h
#define hlm_hybrid_h

//#include "A5pnP15.h"
#include "hlmNewt.h"
#include "hhatlmtail.h"
#include "hlmNQC.h"
#include "deltalm.h"
#include "SpeedyTail.h"
#include "s_flm.h"
#include "ringdown_match.h"
#include "hlm_Tidal.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "cmath"
#include "multipole_index.h"

vector<gsl_complex> hlm(double t, const double phi, const double r, const double pph, const double prstar, double Omega, const double ddotr, const double H, const double Heff,const double jhat, const double rw,void *params){
    
    double nu = (*(input *)params).nu;
    bool tidal_flag = (*(input *)params).tidal;
    bool spin_flag = (*(input *)params).spin;
    bool NQC_flag = (*(input *)params).NQC;
    bool speedytail_flag = (*(input *)params).speedy;
    
    double source[] = {
        jhat,Heff,
        Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,
        Heff,jhat,Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,jhat,Heff,
        Heff,jhat,Heff,jhat,Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,jhat,Heff,jhat,Heff};
    
    int kmax = 35;
    
    // Newtonian waveform
    vector<gsl_complex> hNewt = hlmNewt( rw,Omega,phi, nu,tidal_flag);
    
    // Compute corrections
    double x = gsl_pow_int(rw*Omega,2);
    vector<double> flm(35);
    if (spin_flag==true) {
        flm = s_flm(x,params);
    }else if (spin_flag==false)
    {
        flm = f_lm(x,nu);
    }
    
    // Computing the tail
    const double r0  = 1.213061319425267e+00;   // 2/sqrt(e);
    const double Hreal = H * nu;
    vector<gsl_complex> tlm(kmax);
    if (speedytail_flag==false) {
        tlm = hhatlmTail(Omega,Hreal, r0, L, M);
    }else if (speedytail_flag==true)
    {
        tlm = speedyTail(Omega,Hreal, r0, L, M);
    }
    
    // Residual phase corrections delta_{lm}
    const vector<double> deltalm = EOBdeltalm(Hreal,Omega, nu);

    vector<gsl_complex> h_NQC(kmax);
    if (NQC_flag==true) {
        h_NQC = hlmNQC(nu,r,prstar,Omega,ddotr);
    }
    
    vector<gsl_complex> hlm(kmax);
    for (int k=35; k--;) {
    
            tlm[k].dat[1] += deltalm[k];
            
            // Compute \hat{h}_lm
            hlm[k].dat[0] = hNewt[k].dat[0] * flm[k] * source[k] * tlm[k].dat[0] ;
            hlm[k].dat[1] = hNewt[k].dat[1] + tlm[k].dat[1] ;
            
            //NQC correction
            if (NQC_flag==true) {
                hlm[k].dat[0] *= h_NQC[k].dat[0];
                hlm[k].dat[1] += h_NQC[k].dat[1];
            }

    }
    
    if (tidal_flag==true) {
            
        // Compute tidal contribution
        vector<double> hlmtidal = hlm_Tidal(x, nu);

        // Update waveform
        double p2 = sqrt(1-4*nu); //see in EOBRun file
        for (int k=35; k--; ) {
            
            switch (k) {
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

#endif /* hlm_hybrid_h */
