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

#include "TEOBResumS.h"

/** Point-mass PN coefs for A function */
void acoeffs(double r, double nu, double *a)
{    
  const double u    = 1./r;
  const double logu = log(u);
  const double a5l  =  64./5.;
  const double a5c0 = -4237./60.+2275./512.*pi*pi+256./5.*Log2+128./5.*EulerGamma;
  const double a5c1 = -221./6.+41./32.*pi*pi;
  const double a5c  =  a5c0+nu*a5c1;
  const double a6c = 3097.3*nu*nu-1330.6*nu+81.38;
  const double a6l = (-7004./105.-144./5.*nu);
    
  a[0] = 2.;
  a[1] = (94./3. - 41./32.*pi*pi);  
  a[2] = a5c;
  a[3] = a5l;
  a[4] = a5c+a5l*logu;
  a[5] = a6c;
  a[6] = a6l;
  a[7] = a6c+a6l*logu;
}

/** Pade P^1_5 resummed expression for A function at 5P,
    return also derivatives w.r.t to u */
void A_NumDenom(double r, double *acoefs, double nu, 
		double *N, double *D, double *dN, dounle *dD, double *d2N, double *d2D)
{
    
  /** Shorthands */
  const double u  = 1./r;
  const double u2 = u*u;
  const double u3 = u2*u;
  const double u4 = u3*u;
  const double u5 = u4*u;
  
  /**  Point-mass PN coefs */
  const double a3  = acoefs[0];
  const double a4  = acoefs[1];
  const double a5  = acoefs[4];
  const double a6  = acoefs[7];
  const double a5l = acoefs[3];
  const double a6l = acoefs[6];
  
  /**  Coefficients of the denominator and numerator of the Pade' A function */
  const double C1  = (8.*a3 + 3.*a4 + a5);
  const double C2  = (4.*a4 + 2.*a5 + a6 + a3*(8. - a3*nu));
  const double C3  = (12.*a3 + 4.*a4 + a5);
  const double C4  = (-32. + C3*nu);
  const double C6  = (4.*a3 + a4);
  const double C5  = (16.*a6 + (2.*a4 + a5)*(8.*a3 + 2.*a4 + a5)*nu - C6*a6*nu + a3*a3*C6*nu*nu);
  const double C7  = (64. - (4.*C1 + a6)*nu + a3*a3*nu*nu);
  const double C8  = a6*(8 - a3*nu) + nu*(8.*a3*a3 + 16.*a3*a4 + 4.*a4*a4 + a3*a3*a3*nu) + a5*(16. + a4*nu);
  const double C9  = 4.*(4.*a4 + 2.*a5 + a6) + a3*(8.*a3 + 4.*a4 + a5)*nu;
  const double C10 = (-8 + a3*nu);
  const double C11 = (-32. + C3*nu);
  const double C12 = (-4.*a4 - 2.*a5 - a6 + a3*C10);
  
  /**  WIT */
  *N = 1. + (C7*u)/C11;
  *D = (-32. + nu*(12.*a3 + 4.*a4 + a5 - C2*u - 2.*C2*u2 - (C9)*u3 - (C8)*u4 - C5*u5))/C11;

  *dN = (-((4.*a5l + a6l)*nu*C4) + C4*C7 + a5l*nu*(-64 + nu*(12.*a4 + 4.*a5 + a6 + a3*(32 - a3*nu))))/(C4*C4);

  *dD = (nu*(-((2.*a5l + a6l)*C4) + a5l*nu*C2 + C4*C12 - 2.*(2.*a5l + a6l)*C4*u + 2.*a5l*nu*C2*u + 4.*C4*C12*u + a5l*nu*(C9)*u2 - 3.*(C9)*C4*u2 - C4*(4.*a6l + a5l*(8 + a3*nu))*u2 + a5l*nu*(C8)*u3 - 4.*C4*(C8)*u3 + C4*(a6l*C10 - a5l*(16 + a4*nu))*u3 + a5l*nu*C5*u4 - 5.*C4*C5*u4 + C4*(-2.*(4.*a3 + 2.*a4 + a5)*a5l*nu + a6l*(-16 + 4.*a3*nu + a4*nu))*u4))/(C11*C11);

  *d2N = -((nu*(-32. + (12.*a3 + 4.*a4 + a5 - 2.*a5l)*nu)*(a6l*(-32. + (12.*a3 + 4.*a4 + a5)*nu) + a5l*(-64. + nu*(4.*a4 - a6 + a3*(16. + a3*nu)))))/(gsl_pow_int(-32. + (12.*a3 + 4.*a4 + a5)*nu,3)*u));
    
   *d2D = (nu*(-2048.*(8.*a4 + 4.*a5 + 6.*a5l + 2.*a6 + 3.*a6l) - (a6l*gsl_pow_int(-32 + (12.*a3 + 4.*a4 + a5)*nu,2) + 2.*(16384.*a3 - 64.*(224.*gsl_pow_int(a3,2) + (4.*a4 + a5)*(8.*a4 + 4.*a5 + 2.*a6 + 3.*a6l) + 4.*a3*(40.*a4 + 16.*a5 + 6.*a6 + 9.*a6l))*nu + (12.*a3 + 4.*a4 + a5)*(320.*gsl_pow_int(a3,2) + (4.*a4 + a5)*(8.*a4 + 4.*a5 + 2.*a6 + 3.*a6l) + 4.*a3*(40.*a4 + 16.*a5 + 6.*a6 + 9.*a6l))*gsl_pow_int(nu,2) - 2.*gsl_pow_int(a3,2)*gsl_pow_int(12.*a3 + 4.*a4 + a5,2)*gsl_pow_int(nu,3))*u + 2.*(12.*(4.*a4 + 2.*a5 + a6) + 10.*a6l + 3.*a3*(8.*a3 + 4.*a4 + a5)*nu)*gsl_pow_int(-32 + (12.*a3 + 4.*a4 + a5)*nu,2)*u2 + gsl_pow_int(-32 + (12.*a3 + 4.*a4 + a5)*nu,2)*(a3*nu*(-12.*a6 - 7.*a6l + 12.*gsl_pow_int(a3,2)*nu) + 12.*a5*(16 + a4*nu) + 8.*(12.*a6 + 7.*a6l + 6.*(2.*gsl_pow_int(a3,2) + 4.*a3*a4 + gsl_pow_int(a4,2))*nu))*u3 + gsl_pow_int(-32 + (12.*a3 + 4.*a4 + a5)*nu,2)*(144.*a6l + 20.*(2.*a4 + a5)*(8.*a3 + 2.*a4 + a5)*nu + (4.*a3 + a4)*nu*(-9.*a6l + 20.*gsl_pow_int(a3,2)*nu) - 20.*a6*(-16 + 4.*a3*nu + a4*nu))*u4 + a5l*((-32 + (12.*a3 + 4.*a4 + a5)*nu)*(-64 + nu*(4.*a4 - a6 - 2.*a6l + a3*(16 + a3*nu))) + 2.*nu*(-32.*(120.*a3 + 36.*a4 + 6.*a5 - 3.*a6 - 2.*a6l) + (480.*gsl_pow_int(a3,2) + (4.*a4 + a5)*(12.*a4 - 3.*a6 - 2.*a6l) + 12.*a3*(28.*a4 + 4.*a5 - 3.*a6 - 2.*a6l))*nu + 3.*gsl_pow_int(a3,2)*(12.*a3 + 4.*a4 + a5)*gsl_pow_int(nu,2))*u + 4.*(-32 + (12.*a3 + 4.*a4 + a5)*nu)*(-320 + nu*(20.*a4 - 5.*a6 - 2.*a6l + 5.*a3*(16 + a3*nu)))*u2 - (-8 + a3*nu)*(-32 + (12.*a3 + 4.*a4 + a5)*nu)*(-448 + nu*(28.*a4 - 7.*a6 - 2.*a6l + 7.*a3*(16 + a3*nu)))*u3 - nu*(-32 + (12.*a3 + 4.*a4 + a5)*nu)*(144.*(16.*a3 + 8.*a4 + 4.*a5 + a6) + 32.*a6l - (9.*(96.*gsl_pow_int(a3,2) + 12.*gsl_pow_int(a4,2) + 8.*a4*a5 + gsl_pow_int(a5,2) + a4*a6 + 4.*a3*(16.*a4 + 6.*a5 + a6)) + 2.*(4.*a3 + a4)*a6l)*nu + 9.*gsl_pow_int(a3,2)*(4.*a3 + a4)*gsl_pow_int(nu,2))*u4) + 2.*gsl_pow_int(a5l,2)*nu*(-64 + nu*(4.*a4 - a6 + a3*(16 + a3*nu)))*(-1 + u*(-2 + u*(-4 + u*(-8 + a3*nu + (-16 + 4.*a3*nu + a4*nu)*u)))))/u))/gsl_pow_int(-32 + (12.*a3 + 4.*a4 + a5)*nu,3);
    
}

/** Tidal term of A potential, Resummed expression [Bini&Damour 2014, Bernuzzi+ 2015] */
void A_TResum(double r, double nu, )
{

    //Missing: b3NR (not needed), rlR (is calculated), kTl (yes, this has to be passed).
    
    /** Shorthands */
    const double u  = 1./r;
    const double u2 = u*u;
    const double u3 = u2*u;
    const double u4 = u3*u;
    const double u5 = u3*u2;
    const double u6 = u5*u;
    const double u7 = u6*u;
    const double u8 = u5*u3;
    
    
    /** Compactness of the star */
    double XA = (*(TEOBResumParams *)params).X1;
    double XB = (*(TEOBResumParams *)params).X2;
    
    //-----------------------------------------------------------------------------------
    // Definition of the conservative tidal coefficients \bar{\alpha}_n^{(\ell)}, Eq.(37)
    // of Damour&Nagar, PRD 81, 084016 (2010)
    //-----------------------------------------------------------------------------------
    double kapA2 = (*(TEOBResumParams *)params).kappaAl2; // 3.   * lambdaAl2 * pow(XA, 2.*2 +1.) / q; //Note: kap stands for kappa; see eqn(1) of REF
    
    double kapB2 = (*(TEOBResumParams *)params).kappaBl2; // 3.   * lambdaBl2 * pow(XB, 2.*2 +1.) * q;
    
    double kapT3 = (*(TEOBResumParams *)params).kappaTl3;//kapA3 + kapB3;
    double kapT4 = (*(TEOBResumParams *)params).kappaTl4;//kapA4 + kapB4;
    
    double bar_alph3_1 = (*(TEOBResumParams *)params).bar_alph3_1;//((-2.+15./2.*XA)*kapA3 + (-2.+15./2.*XB)*kapB3)/kapT3;
    double bar_alph3_2 = (*(TEOBResumParams *)params).bar_alph3_2;//((8./3.-311./24.*XA+110./3.*XA*XA)*kapA3 + (8./3.-311./24.*XB+110./3.*XB*XB)*kapB3)/kapT3;
    
    double p      = 4.;// % 4<p<6


    double c1     = 8.53353;
    double c2     = 3.04309;
    double Acub   = 5./2.*u*(1. -  (c1+c2)*u +   c1*c2*u2);
    double dAcub  = 5./2.*   (1. -2*(c1+c2)*u + 3*c1*c2*u2);
    double d2Acub = 5   *   (  -  (c1+c2)   + 3*c1*c2*u);
    double n1     = 0.840058;
    double d2     = 17.73239;
    double Den    = 1./(1. + d2*u2);
    double f23    = (1. + n1*u)*Den;
    double df23   = (n1 - 2*d2*u - n1*d2*u2)*(Den*Den);
    double A1SF   = Acub*f23;
    double dA1SF  = dAcub*f23 + Acub*df23;
    double A2SF   = 337./28.*u2;
    double dA2SF  = 674./28.*u;
    double oom3u  = 1./(1.-rLR*u);
    double f0     = 1 + 3*u2*oom3u;
    double f1     = A1SF *pow(oom3u,7./2.);
    double f2     = A2SF *pow(oom3u,p);
    double df0    = 3*u*(2.-rLR*u)*(oom3u*oom3u);
    double df1    = 0.5*(7*rLR*A1SF + 2*(1.-rLR*u)*dA1SF)*pow(oom3u,9./2.);
    double df2    = (rLR*p*A2SF + (1.-rLR*u)*dA2SF)*pow(oom3u,p+1);
    
    double d2f23  = 2*d2*(-1 + 3*d2*u2 + n1*(-3+d2*u2))*(Den*Den*Den);
    double d2A1SF = d2Acub*f23 + 2*dAcub*df23 + Acub*d2f23;
    double d2A2SF = 674./28.;
    
    double d2f0   = 6*(oom3u*oom3u*oom3u);
    double d2f1   = 0.25*(63*(rLR*rLR)*A1SF + 4*(-1+rLR*u)*(-7*rLR*dA1SF + (-1+rLR*u)*d2A1SF))*pow(oom3u,11./2.);
    double d2f2   = (rLR*p*((1+p)*rLR*A2SF - 2*(-1+rLR*u)*dA2SF +(-1.+rLR*u)*(-1.+rLR*u) *d2A2SF))*pow(oom3u,p+2);
    
    double d2AT2  = - kapA2*30*u4*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*30*u4*( f0 + XB*f1 + XB*XB*f2 ) - 2*kapA2*6*u5*( df0 + XA*df1 + XA*XA*df2 ) - 2*kapB2*6*u5*( df0 + XB*df1 + XB*XB*df2 ) - kapA2*u6*( d2f0 + XA*d2f1 + XA*XA*d2f2 ) - kapB2*u6*( d2f0 + XB*d2f1 + XB*XB*d2f2 );
    
    double d2AT3  = - kapT3*(56*u6 + 72*bar_alph3_1*u7 + 90*bar_alph3_2*u8);
    double d2AT4  = - kapT4*90*u8;
    double A_ddu  =   d2AT2 + d2AT3 + d2AT4;
    
    A_ddu = A_ddu + d2A0_u;
        
}




/** Metric function A at 5PN + Pade P^1_5 resummation, 
    return also derivatives w.r.t to u and r */
void A5pnP15_dd(double r, void *params, double *A, )
{
    
  // fixme:
  double nu         = (*(TEOBResumParams *)params).nu;
  bool   tidal_flag = (*(TEOBResumParams *)params).flags.tidal;
  double rLR        = (*(TEOBResumParams *)params).rLR;
    
  /** Point-mass PN coefs */
  double acoefs[8];
  acoeffs(r, nu, acoefs);
 
  /** Num & Den of A metric function */
  double N,D,dN,dD,d2N,d2D;
  A_NumDenom(r,a,nu, &N,&D, &dN,&dD, &d2N, &d2D);    

  double ooD = 1./D;
    
    
  double A_ddu = (2.*N*dD*dD-2.*D*dD*dN-D*N*d2D+D*D*d2N)*ooD*ooD*ooD; // A_ddu
  //A_ddr=-2.*A5pnP15(r,nu,tidal_flag)[1]/r+A_dd[1]/(r*r*r*r); //A_ddr
    
  double d2A0_u = A_ddu;
  
  if (tidal_flag==true) {
    
    //Missing: b3NR (not needed), rlR (is calculated), kTl (yes, this has to be passed).
    
    /** Shorthands */
    const double u  = 1./r;
    const double u2 = u*u;
    const double u3 = u2*u;
    const double u4 = u3*u;
    const double u5 = u3*u2;
    const double u6 = u5*u;
    const double u7 = u6*u;
    const double u8 = u5*u3;
    
    
    /** Compactness of the star */
    double XA = (*(TEOBResumParams *)params).X1;
    double XB = (*(TEOBResumParams *)params).X2;
    
    //-----------------------------------------------------------------------------------
    // Definition of the conservative tidal coefficients \bar{\alpha}_n^{(\ell)}, Eq.(37)
    // of Damour&Nagar, PRD 81, 084016 (2010)
    //-----------------------------------------------------------------------------------
    double kapA2 = (*(TEOBResumParams *)params).kappaAl2; // 3.   * lambdaAl2 * pow(XA, 2.*2 +1.) / q; //Note: kap stands for kappa; see eqn(1) of REF
    
    double kapB2 = (*(TEOBResumParams *)params).kappaBl2; // 3.   * lambdaBl2 * pow(XB, 2.*2 +1.) * q;
    
    double kapT3 = (*(TEOBResumParams *)params).kappaTl3;//kapA3 + kapB3;
    double kapT4 = (*(TEOBResumParams *)params).kappaTl4;//kapA4 + kapB4;
    
    double bar_alph3_1 = (*(TEOBResumParams *)params).bar_alph3_1;//((-2.+15./2.*XA)*kapA3 + (-2.+15./2.*XB)*kapB3)/kapT3;
    double bar_alph3_2 = (*(TEOBResumParams *)params).bar_alph3_2;//((8./3.-311./24.*XA+110./3.*XA*XA)*kapA3 + (8./3.-311./24.*XB+110./3.*XB*XB)*kapB3)/kapT3;
    
    double p      = 4.;// % 4<p<6
    double c1     = 8.53353;
    double c2     = 3.04309;
    double Acub   = 5./2.*u*(1. -  (c1+c2)*u +   c1*c2*u2);
    double dAcub  = 5./2.*   (1. -2*(c1+c2)*u + 3*c1*c2*u2);
    double d2Acub = 5   *   (  -  (c1+c2)   + 3*c1*c2*u);
    double n1     = 0.840058;
    double d2     = 17.73239;
    double Den    = 1./(1. + d2*u2);
    double f23    = (1. + n1*u)*Den;
    double df23   = (n1 - 2*d2*u - n1*d2*u2)*(Den*Den);
    double A1SF   = Acub*f23;
    double dA1SF  = dAcub*f23 + Acub*df23;
    double A2SF   = 337./28.*u2;
    double dA2SF  = 674./28.*u;
    double oom3u  = 1./(1.-rLR*u);
    double f0     = 1 + 3*u2*oom3u;
    double f1     = A1SF *pow(oom3u,7./2.);
    double f2     = A2SF *pow(oom3u,p);
    double df0    = 3*u*(2.-rLR*u)*(oom3u*oom3u);
    double df1    = 0.5*(7*rLR*A1SF + 2*(1.-rLR*u)*dA1SF)*pow(oom3u,9./2.);
    double df2    = (rLR*p*A2SF + (1.-rLR*u)*dA2SF)*pow(oom3u,p+1);
    
    double d2f23  = 2*d2*(-1 + 3*d2*u2 + n1*(-3+d2*u2))*(Den*Den*Den);
    double d2A1SF = d2Acub*f23 + 2*dAcub*df23 + Acub*d2f23;
    double d2A2SF = 674./28.;
    
    double d2f0   = 6*(oom3u*oom3u*oom3u);
    double d2f1   = 0.25*(63*(rLR*rLR)*A1SF + 4*(-1+rLR*u)*(-7*rLR*dA1SF + (-1+rLR*u)*d2A1SF))*pow(oom3u,11./2.);
    double d2f2   = (rLR*p*((1+p)*rLR*A2SF - 2*(-1+rLR*u)*dA2SF +(-1.+rLR*u)*(-1.+rLR*u) *d2A2SF))*pow(oom3u,p+2);
    
    double d2AT2  = - kapA2*30*u4*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*30*u4*( f0 + XB*f1 + XB*XB*f2 ) - 2*kapA2*6*u5*( df0 + XA*df1 + XA*XA*df2 ) - 2*kapB2*6*u5*( df0 + XB*df1 + XB*XB*df2 ) - kapA2*u6*( d2f0 + XA*d2f1 + XA*XA*d2f2 ) - kapB2*u6*( d2f0 + XB*d2f1 + XB*XB*d2f2 );
    
    double d2AT3  = - kapT3*(56*u6 + 72*bar_alph3_1*u7 + 90*bar_alph3_2*u8);
    double d2AT4  = - kapT4*90*u8;
    double A_ddu  =   d2AT2 + d2AT3 + d2AT4;
    
    A_ddu = A_ddu + d2A0_u;
        
  }
  
  double A_ddr = -2.*Metric(r,params,false)[1]/r+A_ddu/(r*r*r*r);
  
  A_dd[1] = A_ddu;
  A_dd[0] = A_ddr;
  
  return A_dd;
  
}


void s_A5PNlog(double r, void *params, bool nnlo_flag){
    
    /*
     %EOB_A5PNlog function EOB_A5PNlog(r,nu,a5,a6)
     %   This function computes the Pade' (1,5) resummed A function (with its
     %   derivatives) starting from the 5PN-expanded version of the A function
     %   including 4PN and 5PN log terms.
     %
     %   This represents the current, stable, most accurate implementation of
     %   the EOB effective potential
     %
     %   Usage:
     %
     %   [A dA d2A] = EOB_A5PNlog(r,a5,a6)
     %
     %   where a5 and a6 are the nonlog contributions to the 4PN and 5PN terms.
     %   In practice, a5 is fixed to its GSF value computed in Akcay et al,
     %
     %   a5 \equiv a5_GSF = +23.50190(5) \approx +23.5
     %
     %   and a6 \equiv a6(nu) = (-110.5 + 129*(1-4*nu)).*(1-1.5e-5/((0.26-nu)^2)
     %   as obtained from comparison with the Caltech-Cornell-CITA numerical data.
     %   [ the effective case]
     %   These values are used as default.
     %
     %   (c) anagar, IHES, November 12, 2012
     %               revised  March 13, 2013
     */
    
    
    // parameters
    bool tidal_flag = (*(TEOBResumParams *)params).flags.tidal;
    double rLR      = (*(TEOBResumParams *)params).rLR;
    double nu       = (*(TEOBResumParams *)params).nu;
    vector<double> data(5);
    
    // constants
    double nu2 = nu*nu;
    double pi2 = pi*pi;
    double pi4 = pi2*pi2;
    
    double a5c0 = -4237./60. + 2275./512.*pi2 + 256./5.*log(2) + 128./5.*EulerGamma;
    double a5c1 = -221./6.   + 41./32.*pi2;
    double a5   =  a5c0 + nu*a5c1;
    double a6   =  3097.3*nu2 - 1330.6*nu + 81.38;
    
    // shortcuts
    double u    = 1./r;
    double u2   = u*u;
    double u3   = u*u2;
    double u4   = u2*u2;
    double u5   = u4*u;
    double u6   = u5*u;
    double u7   = u6*u;
    double u10  = u5*u5;
    double u8   = u5*u3;
    double u9   = u8*u;
    double logu = log(u);
    
    // 4PN and 5PN coefficients including all known log terms
    double a5tot  = a5  + 64./5.*logu;
    double a6tot  = a6  + (-7004./105. - 144./5.*nu)*logu;
    double a5tot2 = a5tot*a5tot;
    
    // Coefficients of the Padeed function
    double N1 = (-3*(-512 - 32*nu2 + nu*(3520 + 32*a5tot + 8*a6tot - 123*pi2)))/(-768 + nu*(3584 + 24*a5tot - 123*pi2));
    double D1 = (nu*(-3392 - 48*a5tot - 24*a6tot + 96*nu + 123*pi2))/(-768 + nu*(3584 + 24*a5tot - 123*pi2));
    double D2 = (2*nu*(-3392 - 48*a5tot - 24*a6tot + 96*nu + 123*pi2))/(-768 + nu*(3584 + 24*a5tot - 123*pi2));
    double D3 = (-2*nu*(6016 + 48*a6tot + 3392*nu + 24*a5tot*(4 + nu) - 246*pi2 - 123*nu*pi2))/(-768 + nu*(3584 + 24*a5tot - 123*pi2));
    double D4 = -(nu*(-4608*a6tot*(-4 + nu) + a5tot*(36864 + nu*(72192 - 2952*pi2)) + nu*(2048*(5582 + 9*nu) - 834432*pi2 + 15129*pi4)))/(96.*(-768 + nu*(3584 + 24*a5tot - 123*pi2)));
    double D5 = (nu*(-24*a6tot*(1536 + nu*(-3776 + 123*pi2)) + nu*(-2304*a5tot2 + 96*a5tot*(-3392 + 123*pi2) - (-3776 + 123*pi2)*(-3008 - 96*nu + 123*pi2))))/(96.*(-768 + nu*(3584 + 24*a5tot - 123*pi2)));
    
    // First derivatives
    
    double dN1 = (160*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2)*u);
    double dD1 = (160*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2)*u);
    double dD2 = (320*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2)*u);
    double dD3 = (640*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2)*u);
    double dD4 = (-320*(-4 + nu)*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2)*u);
    double dD5 = (nu*(-8400*nu*(-24*(a6 - (4*logu*(1751 + 756*nu))/105.)*(1536 + nu*(-3776 + 123*pi2)) + nu*(-2304*gsl_pow_int(a5 + (64*logu)/5.,2) + 96*(a5 + (64*logu)/5.)*(-3392 + 123*pi2) - (-3776 + 123*pi2)*(-32*(94 + 3*nu) + 123*pi2))) - (1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)))*(4128768*logu*nu + 5*(-2689536 + nu*(11170624 + 64512*a5 - 380685*pi2) - 756*nu*(1536 + nu*(-3776 + 123*pi2))))))/(2625.*gsl_pow_int(-768 + nu*(3584 + 24*(a5 + (64*logu)/5.) - 123*pi2),2)*u);
    
    // Numerator and denominato of the Pade
    double Num = 1 + N1*u;
    double Den = 1 + D1*u + D2*u2 + D3*u3 + D4*u4 + D5*u5;
    double A   = Num/Den;
    
    // First derivative
    double dNum  = dN1*u + N1;
    double dDen  = D1 + u*(dD1 + 2*D2) + u2*(dD2 + 3*D3) + u3*(dD3 + 4*D4) + u4*(dD4 + 5*D5) + dD5*u5;
    
    
    // derivative of A function with respect to u
    double prefactor = A/(Num*Den);
    double dA_u      = prefactor*(dNum*Den - dDen*Num);
    
    //
    if (tidal_flag==true)
    {
        double A0    = A;
        double A0_du = dA_u;
        // CHECK IF THE SAME (BE AWARE ANOTHE DEF BELOW IN OTHER ROUTINE):
        double XA = (*(TEOBResumParams *)params).X1;
        double XB = (*(TEOBResumParams *)params).X2;
        
        double kapA2 = (*(TEOBResumParams *)params).kappaAl2; // 3.   * lambdaAl2 * pow(XA, 2.*2 +1.) / q; //Note: kap stands for kappa; see eqn(1) of REF
        double kapB2 = (*(TEOBResumParams *)params).kappaBl2; // 3.   * lambdaBl2 * pow(XB, 2.*2 +1.) * q;
        
        double kapT2 = (*(TEOBResumParams *)params).kappaTl2; // kapA2 + kapB2;
        double kapT3 = (*(TEOBResumParams *)params).kappaTl3;//kapA3 + kapB3;
        double kapT4 = (*(TEOBResumParams *)params).kappaTl4;//kapA4 + kapB4;
        
        double bar_alph2_1 = (*(TEOBResumParams *)params).bar_alph2_1;//(5./2.*XA*kapA2 + 5./2.*XB*kapB2)/kapT2;
        double bar_alph2_2 = (*(TEOBResumParams *)params).bar_alph2_2;//((3.+XA/8.+ 337./28.*XA*XA)*kapA2 + (3.+XB/8.+ 337./28.*XB*XB)*kapB2)/kapT2;
        double bar_alph3_1 = (*(TEOBResumParams *)params).bar_alph3_1;//((-2.+15./2.*XA)*kapA3 + (-2.+15./2.*XB)*kapB3)/kapT3;
        double bar_alph3_2 = (*(TEOBResumParams *)params).bar_alph3_2;//((8./3.-311./24.*XA+110./3.*XA*XA)*kapA3 + (8./3.-311./24.*XB+110./3.*XB*XB)*kapB3)/kapT3;
        
        //case 'nnlo'
        if (nnlo_flag==true) { //Used for calculating the rLR
            //case 'nnlo'
            A     = -(kapT4*u10) - kapT2*u6*(1. + bar_alph2_1*u + bar_alph2_2*u2) - kapT3*u8*(1. + bar_alph3_1*u + bar_alph3_2*u2);
            dA_u = -10.*kapT4*u9 - kapT2*u6*(bar_alph2_1 + 2.*bar_alph2_2*u) - kapT3*u8*(bar_alph3_1 + 2.*bar_alph3_2*u)
            - 6.*kapT2*u5*(1. + bar_alph2_1*u + bar_alph2_2*u2) - 8.*kapT3*u7*(1. + bar_alph3_1*u + bar_alph3_2*u2);
        } else { //Used for calculting the dynamcis
            //case 'nnlo_gsfLR'; Bini & Damour, 1409.6933 + free light-ring
            // Tidal PN coefs
            double p      =  4.;// % 4<p<6
            double c1     =  8.53353;
            double c2     =  3.04309;
            double Acub   =  5./2.* u * (1. -  (c1+c2)*u +   c1*c2*u2);
            double n1     =  0.840058;
            double d2     =  17.73239;
            double DenI   =  1./(1. + d2*u2);
            double f23    =  (1. + n1*u)*DenI;
            double A1SF   =  Acub*f23;
            double A2SF   =  337./28.*u2;
            double oom3u  =  1./(1.-rLR*u);
            double f0     =  1. + 3.*u2*oom3u;
            double f1     =  A1SF *pow(oom3u,7./2.);
            double f2     =  A2SF *pow(oom3u,p);
            double AT2    = - kapA2*u6*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*u6*( f0 + XB*f1 + XB*XB*f2 );
            double AT3    = - kapT3*u8*(1. + bar_alph3_1*u + bar_alph3_2*u2);
            double AT4    = - kapT4*u10;
            
            A = AT2 + AT3 + AT4;
            
            //Derivative of potential w.r.t. u
            double dAcub = 5./2.*   (1. -2.*(c1+c2)*u + 3.*c1*c2*u2);
            double df23  = (n1 - 2.*d2*u - n1*d2*u2)*(DenI*DenI);
            double dA1SF = dAcub*f23 + Acub*df23;
            double dA2SF = 674./28.*u;
            double df0   = 3.*u*(2.-rLR*u)*(oom3u*oom3u);
            double df1   = 0.5*(7.*rLR*A1SF + 2.*(1.-rLR*u)*dA1SF)*pow(oom3u,9./2.);
            double df2   = (rLR*p*A2SF + (1.-rLR*u)*dA2SF)*pow(oom3u,p+1.);
            double dAT2  = - kapA2*6.*u5*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*6.*u5*( f0 + XB*f1 + XB*XB*f2 ) - kapA2*u6*( df0 + XA*df1 + XA*XA*df2 ) - kapB2*u6*( df0 + XB*df1 + XB*XB*df2 );
            double dAT3  = - kapT3*(8.*u7 + 9*bar_alph3_1*u8 + 10*bar_alph3_2*u9);
            double dAT4  = - kapT4*10.*u9;
            
            dA_u =  dAT2 + dAT3 + dAT4;
            
        }
        A    = A+A0;
        dA_u = dA_u+A0_du;
    }
    
    
    // Derivative of A with respect to r
    double dA    = -u2*dA_u;
    
    // Second derivatives of Pade coefficients
    
    double d2N1 = (160*nu*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*pi2))*(828672 + nu*(-42024*a5 - 8064*a6 + 3584*(-1397 + 9*nu) + 174045*pi2) + 756*nu*(768 + nu*(-3584 - 24*a5 + 123*pi2))))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),3)*u2);
    double d2D1 = (160*nu*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*pi2))*(828672 + nu*(-42024*a5 - 8064*a6 + 3584*(-1397 + 9*nu) + 174045*pi2) + 756*nu*(768 + nu*(-3584 - 24*a5 + 123*pi2))))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),3)*u2);
    double d2D2 = (320*nu*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*pi2))*(828672 + nu*(-42024*a5 - 8064*a6 + 3584*(-1397 + 9*nu) + 174045*pi2) + 756*nu*(768 + nu*(-3584 - 24*a5 + 123*pi2))))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),3)*u2);
    double d2D3 = (640*nu*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*pi2))*(828672 + nu*(-42024*a5 - 8064*a6 + 3584*(-1397 + 9*nu) + 174045*pi2) + 756*nu*(768 + nu*(-3584 - 24*a5 + 123*pi2))))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),3)*u2);
    double d2D4 = (320*(-4 + nu)*nu*(-828672 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 32256*nu - 174045*pi2))*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*pi2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),3)*u2);
    double d2D5 = (nu*(gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2)*(4128768*logu*nu - 7680*(1751 + 756*nu) + nu*(64*(808193 + 5040*a5 + 223020*nu) - 615*(3095 + 756*nu)*pi2)) + 3072*nu*(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)))*(4128768*logu*nu - 7680*(1751 + 756*nu) + 5*nu*(64*(174541 + 1008*a5 + 44604*nu) - 123*(3095 + 756*nu)*pi2)) + 25804800*nu2*(-24*(a6 - (4*logu*(1751 + 756*nu))/105.)*(1536 + nu*(-3776 + 123*pi2)) + nu*(-2304*gsl_pow_int(a5 + (64*logu)/5.,2) + 96*(a5 + (64*logu)/5.)*(-3392 + 123*pi2) - (-3776 + 123*pi2)*(-32*(94 + 3*nu) + 123*pi2))) + 42000*nu*(-768 + nu*(3584 + 24*(a5 + (64*logu)/5.) - 123*pi2))*(-24*(a6 - (4*logu*(1751 + 756*nu))/105.)*(1536 + nu*(-3776 + 123*pi2)) + nu*(-2304*gsl_pow_int(a5 + (64*logu)/5.,2) + 96*(a5 + (64*logu)/5.)*(-3392 + 123*pi2) - (-3776 + 123*pi2)*(-32*(94 + 3*nu) + 123*pi2)))))/(13125.*gsl_pow_int(-768 + nu*(3584 + 24*(a5 + (64*logu)/5.) - 123*pi2),3)*u2);
    
    // Second derivative of numerator and denominator
    double d2Num = 2.*dN1 + d2N1*u;
    double d2Den = 2.*(D2 + dD1) + u*(6.*D3 + 4.*dD2 + d2D1) + u2*(12.*D4 + 6.*dD3 + d2D2) + u3*(20.*D5 + 8.*dD4 + d2D3) + u4*(10.*dD5 + d2D4) + u5*d2D5;
    
    // second derivative with respect of u
    double d2A_u = prefactor*(2.*dDen*dDen*A - 2.*dNum*dDen + Den*d2Num - d2Den*Num);
    
    // second derivative with respect of r
    double d2A = u4*d2A_u + 2.*u3*dA_u;
    
    //-----------------------
    // The D and dD functions
    //-----------------------
    double Dp  = 1.0 + 6.*nu*u2 - 2.*(3.0*nu-26.0)*nu*u3; // Pade' resummation of D
    double D   = 1./Dp;
    double dD  = 6.*u2*(2.*nu*u-(3.*nu-26.)*nu*u2)*D*D;
    
    return {A,dA,d2A,D,dD};
}










/** */
vector<double> Metric(const double r, void *params, bool nnlo_flag)
{
    
    double nu         = (*(TEOBResumParams *)params).nu;
    bool   tidal_flag = (*(TEOBResumParams *)params).flags.tidal;
    double rLR        = (*(TEOBResumParams *)params).rLR;
    vector<double> data(5);
    
    const double u   = 1./r;
    const double u2  = u*u;
    const double u3  = u2*u;
    const double u5  = u3*u2;
    const double u6  = u5*u;
    const double u7  = u6*u;
    const double u8  = u5*u3;
    const double u9  = u8*u;
    const double u10 = u5*u5;
    const vector<double> a = acoeffs(r,nu);
    
    /** Unpack numerator and denominator and compute A and its derivatives */
    const vector<double> frac = A_NumDenom(r, a, nu);
    const double Num  = frac[0];
    const double Den  = frac[1];
    const double dNum = frac[2];
    const double dDen = frac[3];
    
    const double ooD = 1./Den;
    double A    = Num*ooD;
    double A_du = (-Num*dDen+Den*dNum)*ooD*ooD; //A_du
    
    double A0    = A;
    double A0_du = A_du;
    
    if (tidal_flag==true)
    {
        
        //Missing: b3NR (not needed), rlR (is calculated), kTl (yes, this has to be passed).
        
        double XA = (*(TEOBResumParams *)params).X1;
        double XB = (*(TEOBResumParams *)params).X2;
        
        /** Computing the tidal coupling constants */
        double kapA2 = (*(TEOBResumParams *)params).kappaAl2; // 3.   * lambdaAl2 * pow(XA, 2.*2 +1.) / q; //Note: kap stands for kappa; see eqn(1) of REF
        
        double kapB2 = (*(TEOBResumParams *)params).kappaBl2; // 3.   * lambdaBl2 * pow(XB, 2.*2 +1.) * q;
        
        double kapT2 = (*(TEOBResumParams *)params).kappaTl2; // kapA2 + kapB2;
        double kapT3 = (*(TEOBResumParams *)params).kappaTl3;//kapA3 + kapB3;
        double kapT4 = (*(TEOBResumParams *)params).kappaTl4;//kapA4 + kapB4;
        
        double bar_alph2_1 = (*(TEOBResumParams *)params).bar_alph2_1;//(5./2.*XA*kapA2 + 5./2.*XB*kapB2)/kapT2;
        double bar_alph2_2 = (*(TEOBResumParams *)params).bar_alph2_2;//((3.+XA/8.+ 337./28.*XA*XA)*kapA2 + (3.+XB/8.+ 337./28.*XB*XB)*kapB2)/kapT2;
        double bar_alph3_1 = (*(TEOBResumParams *)params).bar_alph3_1;//((-2.+15./2.*XA)*kapA3 + (-2.+15./2.*XB)*kapB3)/kapT3;
        double bar_alph3_2 = (*(TEOBResumParams *)params).bar_alph3_2;//((8./3.-311./24.*XA+110./3.*XA*XA)*kapA3 + (8./3.-311./24.*XB+110./3.*XB*XB)*kapB3)/kapT3;
        
	       
        /** Case 'nnlo' */
        if (nnlo_flag==true)
        { //Used for calculating the rLR
            //case 'nnlo'
            A    = -(kapT4*u10) - kapT2*u6*(1. + bar_alph2_1*u + bar_alph2_2*u2) - kapT3*u8*(1. + bar_alph3_1*u + bar_alph3_2*u2);
            A_du = -10.*kapT4*u9 - kapT2*u6*(bar_alph2_1 + 2.*bar_alph2_2*u) - kapT3*u8*(bar_alph3_1 + 2.*bar_alph3_2*u)
            - 6.*kapT2*u5*(1. + bar_alph2_1*u + bar_alph2_2*u2) - 8.*kapT3*u7*(1. + bar_alph3_1*u + bar_alph3_2*u2);
        }
        else
        {   //Used for calculting the dynamcis
            //case 'nnlo_gsfLR'; Bini & Damour, 1409.6933 + free light-ring
            // Tidal PN coefs
            double p      =  4.;// % 4<p<6
            double c1     =  8.53353;
            double c2     =  3.04309;
            double Acub   =  5./2.* u * (1. -  (c1+c2)*u +   c1*c2*u2);
            double n1     =  0.840058;
            double d2     =  17.73239;
            double DenI   =  1./(1. + d2*u2);
            double f23    =  (1. + n1*u)*DenI;
            double A1SF   =  Acub*f23;
            double A2SF   =  337./28.*u2;
            double oom3u  =  1./(1.-rLR*u);
            double f0     =  1. + 3.*u2*oom3u;
            double f1     =  A1SF *pow(oom3u,7./2.);
            double f2     =  A2SF *pow(oom3u,p);
            double AT2    = - kapA2*u6*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*u6*( f0 + XB*f1 + XB*XB*f2 );
            double AT3    = - kapT3*u8*(1. + bar_alph3_1*u + bar_alph3_2*u2);
            double AT4    = - kapT4*u10;
            
            
            A = AT2 + AT3 + AT4;
            
            /** Derivative of potential w.r.t. u */
            double dAcub = 5./2.*   (1. -2.*(c1+c2)*u + 3.*c1*c2*u2);
            double df23  = (n1 - 2.*d2*u - n1*d2*u2)*(DenI*DenI);
            double dA1SF = dAcub*f23 + Acub*df23;
            double dA2SF = 674./28.*u;
            double df0   = 3.*u*(2.-rLR*u)*(oom3u*oom3u);
            double df1   = 0.5*(7.*rLR*A1SF + 2.*(1.-rLR*u)*dA1SF)*pow(oom3u,9./2.);
            double df2   = (rLR*p*A2SF + (1.-rLR*u)*dA2SF)*pow(oom3u,p+1.);
            double dAT2  = - kapA2*6.*u5*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*6.*u5*( f0 + XB*f1 + XB*XB*f2 ) - kapA2*u6*( df0 + XA*df1 + XA*XA*df2 ) - kapB2*u6*( df0 + XB*df1 + XB*XB*df2 );
            double dAT3  = - kapT3*(8.*u7 + 9.*bar_alph3_1*u8 + 10.*bar_alph3_2*u9);
            double dAT4  = - kapT4*10.*u9;
            
            A_du =  dAT2 + dAT3 + dAT4;
            
        }
        A    = A+A0;
        A_du = A_du+A0_du;
    }
    
    double A_dr = - A_du/(r*r); //A_dr
    
    /** Compute D and B functions */
    const double Dp   = 1. + 6.*nu*u2 - 2.*(3.*nu-26.)*nu*u3;
    const double D    = 1./Dp;
    const double D_dr = 6.*u2*(2.*nu*u-(3.*nu-26.)*nu*u2)*D*D;
    
    const double B    = D/A; // B
    const double B_dr = (D_dr*A - D*A_dr)/(A*A); // dB
    
    data[0] = A; data[1] = A_dr; data[2] = A_du; data[3] = B; data[4] = B_dr;
    
    return data;
}


vector<double> s_Metric(double r, void *params, bool nnlo_flag){
    /*
     %                       This function computes the EOB metric potentials
     %                       A(r), B(r), and their derivatives, as functions
     %                       of the EOB (Boyer-Lindquist) radius r.
     %                       USAGE:
     %
     %                       [A dA d2A B dB] = EOB_Metric(r,d2A_flag)
     %
     %                       where r is the EOB radius; dA,dB and d2A are the
     %                       first and second derivatives of A and B. The
     %                       default value for d2A_flag is 'no'. One must set
     %                       it up to 'yes' if the computation of the second
     %                       radial derivative of A is needed.
     */
    
    vector<double> rc_vec = s_get_rc(r,params); //[rc, drc, d2rc]
    double rc   = rc_vec[0];
    double drc  = rc_vec[1];
    double d2rc = rc_vec[2];
    
    double r2  = r*r;
    double u   = 1./r;
    double u2  = u*u;
    double u3  = u2*u;
    double u4  = u2*u2;
    double uc  = 1./rc;
    double uc2 = uc*uc;
    double uc3 = uc2*uc;
    
    
    vector<double> metric = s_A5PNlog(rc,params,nnlo_flag);
    double Aorb   = metric[0];
    double dAorb  = metric[1];
    double d2Aorb = metric[2];
    double Dorb   = metric[3];
    
    double AKerr_Multipole = (1.+2.*uc)/(1.+2.*u);
    
    double fss = 1.;
    
    // Correction factor
    double A   = Aorb*AKerr_Multipole*fss;
    
    double dA  = dAorb*drc*(1.+2.*uc)/(1.+2.*u) - 2.*Aorb*drc*uc2/(1.+2.*u) + 2.*Aorb*(1.+2.*uc)*u2/((1.+2.*u)*(1.+2.*u));
    
    double d2A = d2Aorb*(1.+2.*uc)/(1.+2.*u) + 4.*dAorb*( u2*(1.+2.*uc)/((1.+2.*u)*(1.+2.*u)) - uc2/(1.+2.*u)*drc) + Aorb*(-4.*u3*(1.+2.*uc)/((1.+2.*u)*(1.+2.*u)) + 8.*u4*(1.+2.*uc)/((1.+2.*u)*(1.+2.*u)*(1.+2.*u))+4.*uc3*(1.+2.*u)*drc*drc - 2.*uc2/(1.+2.*u)*d2rc);
    
    
    // The B function
    double B  = r2*uc2*Dorb/A;
    
    return {A,B,dA,d2A};
}


