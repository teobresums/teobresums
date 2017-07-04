//
//  initial_Spin.h
//
//  Created by Philipp Fleig on 26/04/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef initial_Spin_h
#define initial_Spin_h

#include "s_Flux.h"
#include <gsl/gsl_math.h>
#include "s_Metric.h"
#include "s_get_rc.h"
#include "s_bisec.h"
#include "s_D1.h"
#include "input_struc.h"
#include "multipole_index.h"

using namespace std;

vector<double> s_initial(input *params){

/*
% EOB_ModinSpin
% This function computes post-post-circular (ppc) initial
% data at relative separation r0.
%
% USAGE:
%
% [pph0, pr0, prstar0, j0] = EOB_ModinSpin(r0)
%
% where
%
% r0       => relative separation
% pph0     => post-post-circular angular momentum
% pr       => post-circular radial momentum
% pprstar0 => post-circular r*-conjugate radial momentum
% j0       => circular angular momentum
%
% It consists of three step:
%
% 0. Compute j0                         =>           circular ID, j!=0, pr =0
% 1. From j0, compute pr*               =>      post circular ID, j!=0, pr !=0
% 2. From pr* and j0, re-compute pph0   => post-post-circular ID, pph0!=j!=0, pr !=0
%
*/
    
    double nu = (*params).nu;
    double r0 = (*params).r0;
    double chi1 = (*params).chi1;
    double chi2 = (*params).chi2;
    double S1 = (*params).S1;
    double S2 = (*params).S2;
    double c3 = (*(input *)params).cN3LO;
    bool tidal_flag = (*params).tidal;

//-----------------------------------------------------------------
// Build  a small  grid (2*N points) around the initial position r0
// This is used later to compute a spatial derivative numerically
//-----------------------------------------------------------------
    
    vector<double> y_init(7);
    
    int N  = 10;
    const double dr = 1.e-8;
    vector<double> r(2*N);
    
    // Angular momentum for circular orbit: circular ID
    vector<double> j2(2*N);
    vector<double> j(2*N);
    vector<double> djdr(2*N);
    
    // For circular orbit at r0=r(N)
    vector<double> E0(2*N);         // real Hamiltonian      H_0
    vector<double> Omega_j(2*N);                  // Orbital frequency (from Hamilton's equation)
    
    
    vector<double> Fphi(2*N);
    
    vector<double> prstar(2*N);
    vector<double> pr(2*N);
    vector<double> pph(2*N);
    
    vector<double> dprstardt(2*N);  // NOTE: Fr* here
    
    double X1  = (*(input *)params).X1;
    double X2  = (*(input *)params).X2;
    // Kerr parameter
    double a1  = (*(input *)params).a1;
    double a2  = (*(input *)params).a2;
    double aK2 = (*(input *)params).aK2;

    // spin variable (with dimensions)
    double S     = S1 + S2;        // => in the EMRL this becomes the spin of the BH
    double Ss = X2*a1 + X1*a2;  // => in the EMRL this becomes the spin of the particle

    vector<double> rc(2*N);
    vector<double> drc(2*N);
    vector<double> d2rc(2*N);
    
    double rorb;
    double pphorb;
    
    vector<double> metric(5);
    vector<double> A(2*N);
    vector<double> dA(2*N);
    vector<double> B(2*N);
    vector<double> d2A(2*N);
    for (int i=2*N; i--;) {
        
        r[i]  = r0+(i-N+1)*dr;

	/*if (tidal_flag==true) {
            metric = Metric(r[i], params,false);
            A[i] = metric[0];
            B[i] = metric[3];
            dA[i] = metric[1];
            d2A[i] = A5pnP15_dd(r[i],params)[0];
        } else {
            metric = s_Metric(r[i], params,false);
            A[i] = metric[0];
            B[i] = metric[1];
            dA[i] = metric[2];
            d2A[i] = metric[3];
	    }*/

	// metric is here. The tidal parameters are within this routine 
	metric = s_Metric(r[i], params,false);
        A[i]   = metric[0];
        B[i]   = metric[1];
        dA[i]  = metric[2];
        d2A[i] = metric[3];
		
        vector<double> rc_rad = s_get_rc(r[i],params);
        rc[i] = rc_rad[0];
        drc[i] = rc_rad[1];
        
        //Compute minimum of Heff0 using bisection method
        rorb   = r[i];
        pphorb = rorb/sqrt(rorb-3.);
        pph[i] = s_bisec(pphorb,rorb,A[i],dA[i],rc[i],drc[i],aK2,S,Ss,params);
    }

    vector<double> dpph_dr = s_D1(pph,r,12-1); // derivative is computed on a grid with 12 points
    
    for (int i=2*N; i--;) {
        
        double sqrtAbyB     = sqrt(A[i]/B[i]);
        
        double uc  = 1./rc[i];
        double uc2 = uc*uc;
        
        // circular angular momentum
        double pph2 = pph[i]*pph[i];
        
        // Orbital effective Hamiltonian
        double Horbeff0   = sqrt(A[i]*(1. + pph2*uc2));
        
        // Compute gyro-gravitomagnetic coupling functions
        vector<double> ggm0 = s_GS(r[i],rc[i],drc[i],aK2,0.,pph[i],nu,chi1,chi2,X1,X2,c3);
        double GS_0 = ggm0[2];
        double GSs_0 = ggm0[3];
        double dGS_dr_0 = ggm0[6];
        double dGSs_dr_0 = ggm0[7];
        double dGSs_dpph_0 = ggm0[9];
        double dGS_dprstarbyprstar_0 = ggm0[10];
        double dGSs_dprstarbyprstar_0 = ggm0[11];
        
        // Final effective Hamiltonian
        double Heff0    = (GS_0*S + GSs_0*Ss)*pph[i] + Horbeff0;
        
        // Real Hamiltonian: beware that this is NOT divided by nu
        double H0     = sqrt( 1. + 2.*nu*(Heff0 - 1.));
        double one_H0 = 1./H0;
        
        // get gyro-gravitomagnetic (derivative) functions
        double dHeff_dprstarbyprstar = pph[i]*(dGS_dprstarbyprstar_0*S + dGSs_dprstarbyprstar_0*Ss) + 1./Horbeff0;
        
        double C0 = sqrtAbyB*one_H0*dHeff_dprstarbyprstar;
        
        // orbital frequency for circular orbit
        double dHeff_dpph = GS_0*S + (GSs_0 + pph[i]*dGSs_dpph_0)*Ss + pph[i]*A[i]*uc2/Horbeff0;
        double Omg        = one_H0*dHeff_dpph;
        
        //================
        // computing F_phi
        //================
        
        // get the psi correctly with the spin contribution
        double Gtilde     =  GS_0*S     + GSs_0*Ss;
        double dGtilde_dr =  dGS_dr_0*S + dGSs_dr_0*Ss;
        double duc_dr     = -uc2*drc[i];
        double psic       = (duc_dr + dGtilde_dr*rc[i]*sqrt(A[i]/pph2 + A[i]*uc2)/A[i])/(-0.5*dA[i]);
        double r_omg      =  pow((pow(rc[i]*rc[i]*rc[i]*psic,-1./2)+Gtilde)*one_H0,-2./3.);
        double v_phi      =  r_omg*Omg;
        double x          =  v_phi*v_phi;
        double jhat       =  pph[i]/(r_omg*v_phi); // Newton-normalized angular momentum

        Fphi[i] = s_Flux(x,Omg,r_omg,H0,Heff0,jhat,r[i],0.,0.,params);
        
        //=================================================================
        // 1. The radial momentum conjugate to r*: post-circular correction
        //=================================================================
    
        prstar[i] = Fphi[i]/(dpph_dr[i]*C0);
        pr[i]     = prstar[i]* sqrt(B[i]/A[i]);
        
        j[i]       = pph[i];
        E0[i]      = H0;
        Omega_j[i] = Omg;
        
    }

    y_init[0] = r[N-1];
    y_init[1] = pph[N-1];
    y_init[2] = prstar[N-1];
    y_init[3] = pr[N-1];
    y_init[4] = j[N-1];
    y_init[5] = E0[N-1];
    y_init[6] = Omega_j[N-1];
    
    return y_init;
}

#endif /* initial_Spin_h */
