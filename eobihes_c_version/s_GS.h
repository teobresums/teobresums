//
//  EOB_GS.h
//
//  Created by Philipp Fleig on 21/04/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef EOB_GS_h
#define EOB_GS_h

//double c3_fit_global(void *params){
double c3_fit_global(double nu,double chi1,double chi2,double X1,double X2,double a1,double a2,bool tidal_flag){

    double c3 =0.;
    if (tidal_flag==true) {
        c3 =0.;
    } else {
    double nu2      = nu*nu;
    double nu3      = nu2*nu;
    
    // equal-mass, equal-spin coefficients
    //----------------------------------------------------
    // NEW values used in the paper Nagar et al. The value
    // in Eq. (12) was kept, by mistake, to c0 = 44.786477
    // that is an old value obtained with Fitting_c3.m
    //----------------------------------------------------
    double c0 =  44.822889;
    double n1 =  -1.879350;
    double n2 =   0.894242;
    double d1 =  -0.797702;
    
    
    // the two pieces
    double c3_eq   = c0*(1 + n1*(a1+a2) + n2*(a1+a2)*(a1+a2))/(1.+d1*(a1+a2));
    
    //-----------------------------------
    // New fit: different functional form
    //-----------------------------------
    double cnu    = 1222.36;
    double cnu2   = -12764.4;
    double cnu3   =  36689.6;
    double ca1_a2 = -358.086;
    
    double c3_uneq = cnu*(a1+a2)*nu*sqrt(1-4*nu) + cnu2*(a1+a2)*nu2*sqrt(1-4*nu) + cnu3*(a1+a2)*nu3*sqrt(1-4*nu) + ca1_a2*(a1-a2)*nu2;
    
    c3 = c3_eq + c3_uneq;
    }
    
    return c3;
}

vector<double> s_GS(double r,double rc,double drc_dr,double aK2,double prstar,double pph,double nu,double chi1,double chi2,double X1,double X2,double cN3LO){
    
/*
% EOB_GetGSs(r,aK2,prstar,nu). This function computes the
% gyro-gravitomagnetic functions GS and GS*, that are called GS and GSs.
%
% Usage: ggm=EOB_GetGSs(r,aK2,prstar,nu)
%
% where ggm is the output structure. Then we have:
%
% r      => BL radius
% aK2    => squared Kerr parameter
% prstar => r* conjugate momentum
% nu     => symmetric mass ratio
%
% the CN3LO parameter is hard-coded in this routine and can be modified
% here
%
% (c) nagar@ihes.fr, January 2013
*/

double cN4LO = 0.;

double nu2   = nu*nu;

// Boyer-Lindquist radius
double u   = 1./r;
double u2  = u*u;

double uc      = 1./rc;
double uc2     = uc*uc;
double uc3     = uc2*uc;
double uc4     = uc3*uc;
double prstar2 = prstar*prstar;
double prstar4 = prstar2*prstar2;

double GS0       = 2.*u*uc2;
double dGS0_duc  = 2*u2/drc_dr + 4*u*uc;

double GSs0          =  3./2.*uc3;
double dGSs0_duc     =  9./2.*uc2;
double dGSs0_dprstar =  0.0;
double dGSs0_dpph    =  0.0;

// coefficients of hat{GS}
double c10 =  5./16.*nu;
double c20 =  51./8.*nu + 41./256.*nu2;
double c30 =  nu*cN3LO;
double c40 =  nu*cN4LO;
double c02 =  27./16.*nu;
double c12 =  12.*nu - 49./128.*nu2;
double c04 = -5./16.*nu + 169./256.*nu2;

// coefficients of hat{GS*}
double cs10 = 3./4.   + nu/2.;
double cs20 = 27./16. + 29./4.*nu + 3./8.*nu2;
double cs02 = 5./4.   + 3./2.*nu;
double cs12 = 4.   + 11.*nu     - 7./8.*nu2;
double cs04 = 5./48.  + 25./12.*nu + 3./8.*nu2;
double cs30 = nu*cN3LO + 135./32.;
double cs40 = nu*cN4LO + 2835./256.;

double hGS  =  1./(1.  + c10*uc + c20*uc2 + c30*uc3 + c02*prstar2      + c12*uc*prstar2 + c04*prstar4); //write

double hGSs = 1./(1.  + cs10*uc + cs20*uc2  + cs30*uc3 + cs40*uc4 + cs02*prstar2 + cs12*uc*prstar2 + cs04*prstar4); //write

// complete gyro-gravitomagnetic functions
double GS  =  GS0*hGS; //write
double GSs = GSs0*hGSs; //write

//--------------------------------------------------
// Get derivatives of gyro-gravitomagnetic functions
//--------------------------------------------------
double dhGS_dprstar  = -2.*prstar*hGS*hGS *( c02 +  c12*uc +  2.*c04*prstar2);
double dhGSs_dprstar = -2.*prstar*hGSs*hGSs*(cs02 + cs12*uc + 2.*cs04*prstar2);


double dGS_dprstar  = GS0 *dhGS_dprstar; //write
double dGSs_dprstar = GSs0*dhGSs_dprstar + dGSs0_dprstar*hGSs; //write

// derivatives of hat{G} with respect to uc
double dhGS_duc  = -hGS*hGS*(  c10 + 2.*c20*uc  + 3.*c30*uc2  + 4.*c40*uc3);
double dhGSs_duc = -hGSs*hGSs*(cs10 + 2.*cs20*uc + 3.*cs30*uc2 + 4.*cs40*uc3);

// derivatives of G with respect to uc
double dGS_duc  =  dGS0_duc*hGS  +  GS0*dhGS_duc;
double dGSs_duc = dGSs0_duc*hGSs + GSs0*dhGSs_duc;

// derivatives of (G,G*) with respect to r
double dGS_dr  = -drc_dr*uc2*dGS_duc; //write
double dGSs_dr = -drc_dr*uc2*dGSs_duc; //write

// derivatives of (G,G*) with respect to pph
double dGS_dpph  = 0.; //write
double dGSs_dpph = dGSs0_dpph*hGSs; //write


// For initial data: compute the two ratios of ggm.dG_dprstar/prstar for GS
// and GSs
double dGS_dprstarbyprstar  = -2*GS0*hGS*hGS *( c02  +  c12*uc +  2.*c04*prstar2);
double dGSs_dprstarbyprstar = -2*GSs0*hGSs*hGSs*(cs02 + cs12*uc + 2.*cs04*prstar2);

// --- for NQC --
// Second derivatives neglecting all pr_star^2 terms
double d2GS_dprstar20  =  GS0*(-2*hGS*hGS *( c02 +  c12*uc +  2.*c04*prstar2));
double d2GSs_dprstar20 =  GSs0*(-2*hGSs*hGSs*(cs02 + cs12*uc + 2.*cs04*prstar2));

    return {hGS,hGSs,GS,GSs,dGS_dprstar,dGSs_dprstar,dGS_dr,dGSs_dr,dGS_dpph,dGSs_dpph,dGS_dprstarbyprstar,dGSs_dprstarbyprstar,d2GS_dprstar20,d2GSs_dprstar20};
}

#endif /* EOB_GS_h */
