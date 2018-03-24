/**
 * Copyright (C) 2017 Sebastiano Bernuzzi, Gregorio Carullo, Walter Del Pozzo, Alessandro Nagar, Ka Wa Tsang
 * This file is part of TEOBResumS
 *
 * TEOBResumS is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * TEOBResumS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 */

#include "TEOBResumS.h"

/** EOB Metric A function 5PN log 
    This function computes the Pade' (1,5) resummed A function (with its
    derivatives) starting from the 5PN-expanded version of the A function
    including 4PN and 5PN log terms.
    This represents the current, stable, most accurate implementation of
    the EOB effective potential
    
    Coefficients a5 and a6 are the nonlog contributions to the 4PN and 5PN terms.
    In practice, a5 is fixed to its GSF value computed in Akcay et al,
    
    a5 \equiv a5_GSF = +23.50190(5) \approx +23.5
    
    and a6 \equiv a6(nu) = (-110.5 + 129*(1-4*nu)).*(1-1.5e-5/((0.26-nu)^2)
    as obtained from comparison with the Caltech-Cornell-CITA numerical data.
   These values are used as default. */
void eob_metric_A5PNlog(double r, double nu, double *A, double *dA, double *d2A)
{

  /* shortcuts */
  double nu2 = nu*nu;
  double pi2 = Pi*Pi;
  double pi4 = pi2*pi2;
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

  double a5c0 = -4237./60. + 2275./512.*pi2 + 256./5.*Log2 + 128./5.*EulerGamma;
  double a5c1 = -221./6.   + 41./32.*pi2;
  double a5   =  a5c0 + nu*a5c1;
  double a6   =  3097.3*nu2 - 1330.6*nu + 81.38;
  
  /* 4PN and 5PN coefficients including all known log terms */
  double a5tot  = a5  + 64./5.*logu;
  double a6tot  = a6  + (-7004./105. - 144./5.*nu)*logu;
  double a5tot2 = a5tot*a5tot;
  
  /* Coefficients of the Padeed function */
  double N1 = (-3*(-512 - 32*nu2 + nu*(3520 + 32*a5tot + 8*a6tot - 123*pi2)))/(-768 + nu*(3584 + 24*a5tot - 123*pi2));
  double D1 = (nu*(-3392 - 48*a5tot - 24*a6tot + 96*nu + 123*pi2))/(-768 + nu*(3584 + 24*a5tot - 123*pi2));
  double D2 = (2*nu*(-3392 - 48*a5tot - 24*a6tot + 96*nu + 123*pi2))/(-768 + nu*(3584 + 24*a5tot - 123*pi2));
  double D3 = (-2*nu*(6016 + 48*a6tot + 3392*nu + 24*a5tot*(4 + nu) - 246*pi2 - 123*nu*pi2))/(-768 + nu*(3584 + 24*a5tot - 123*pi2));
  double D4 = -(nu*(-4608*a6tot*(-4 + nu) + a5tot*(36864 + nu*(72192 - 2952*pi2)) + nu*(2048*(5582 + 9*nu) - 834432*pi2 + 15129*pi4)))/(96.*(-768 + nu*(3584 + 24*a5tot - 123*pi2)));
  double D5 = (nu*(-24*a6tot*(1536 + nu*(-3776 + 123*pi2)) + nu*(-2304*a5tot2 + 96*a5tot*(-3392 + 123*pi2) - (-3776 + 123*pi2)*(-3008 - 96*nu + 123*pi2))))/(96.*(-768 + nu*(3584 + 24*a5tot - 123*pi2)));
  
  /* First derivatives */
  double dN1 = (160*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2)*u);
  double dD1 = (160*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2)*u);
  double dD2 = (320*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2)*u);
  double dD3 = (640*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2)*u);
  double dD4 = (-320*(-4 + nu)*nu*(-828672 - 32256*nu2 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 174045*pi2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2)*u);
  double dD5 = (nu*(-8400*nu*(-24*(a6 - (4*logu*(1751 + 756*nu))/105.)*(1536 + nu*(-3776 + 123*pi2)) + nu*(-2304*gsl_pow_int(a5 + (64*logu)/5.,2) + 96*(a5 + (64*logu)/5.)*(-3392 + 123*pi2) - (-3776 + 123*pi2)*(-32*(94 + 3*nu) + 123*pi2))) - (1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)))*(4128768*logu*nu + 5*(-2689536 + nu*(11170624 + 64512*a5 - 380685*pi2) - 756*nu*(1536 + nu*(-3776 + 123*pi2))))))/(2625.*gsl_pow_int(-768 + nu*(3584 + 24*(a5 + (64*logu)/5.) - 123*pi2),2)*u);
  
  /* Numerator and denominato of the Pade */
  double Num = 1 + N1*u;
  double Den = 1 + D1*u + D2*u2 + D3*u3 + D4*u4 + D5*u5;
  *A = Num/Den;
    
  /* First derivative */
  double dNum  = dN1*u + N1;
  double dDen  = D1 + u*(dD1 + 2*D2) + u2*(dD2 + 3*D3) + u3*(dD3 + 4*D4) + u4*(dD4 + 5*D5) + dD5*u5;
  
  /* Derivative of A function with respect to u */
  double prefactor = (*A)/(Num*Den);
  double dA_u      = prefactor*(dNum*Den - dDen*Num);

  /* Derivative of A with respect to r */
  /* *dA = -u2*dA_u; */

  *dA = dA_u;

  if (d2A != NULL) {
    
    /* Second derivatives of Pade coefficients */
    double d2N1 = (160*nu*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*pi2))*(828672 + nu*(-42024*a5 - 8064*a6 + 3584*(-1397 + 9*nu) + 174045*pi2) + 756*nu*(768 + nu*(-3584 - 24*a5 + 123*pi2))))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),3)*u2);
    double d2D1 = (160*nu*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*pi2))*(828672 + nu*(-42024*a5 - 8064*a6 + 3584*(-1397 + 9*nu) + 174045*pi2) + 756*nu*(768 + nu*(-3584 - 24*a5 + 123*pi2))))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),3)*u2);
    double d2D2 = (320*nu*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*pi2))*(828672 + nu*(-42024*a5 - 8064*a6 + 3584*(-1397 + 9*nu) + 174045*pi2) + 756*nu*(768 + nu*(-3584 - 24*a5 + 123*pi2))))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),3)*u2);
    double d2D3 = (640*nu*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*pi2))*(828672 + nu*(-42024*a5 - 8064*a6 + 3584*(-1397 + 9*nu) + 174045*pi2) + 756*nu*(768 + nu*(-3584 - 24*a5 + 123*pi2))))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),3)*u2);
    double d2D4 = (320*(-4 + nu)*nu*(-828672 + 756*nu*(-768 + nu*(3584 + 24*a5 - 123*pi2)) + nu*(5006848 + 42024*a5 + 8064*a6 - 32256*nu - 174045*pi2))*(-3840 + 1536*logu*nu + nu*(20992 + 120*a5 - 615*pi2)))/(7.*gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),3)*u2);
    double d2D5 = (nu*(gsl_pow_int(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)),2)*(4128768*logu*nu - 7680*(1751 + 756*nu) + nu*(64*(808193 + 5040*a5 + 223020*nu) - 615*(3095 + 756*nu)*pi2)) + 3072*nu*(1536*logu*nu + 5*(-768 + nu*(3584 + 24*a5 - 123*pi2)))*(4128768*logu*nu - 7680*(1751 + 756*nu) + 5*nu*(64*(174541 + 1008*a5 + 44604*nu) - 123*(3095 + 756*nu)*pi2)) + 25804800*nu2*(-24*(a6 - (4*logu*(1751 + 756*nu))/105.)*(1536 + nu*(-3776 + 123*pi2)) + nu*(-2304*gsl_pow_int(a5 + (64*logu)/5.,2) + 96*(a5 + (64*logu)/5.)*(-3392 + 123*pi2) - (-3776 + 123*pi2)*(-32*(94 + 3*nu) + 123*pi2))) + 42000*nu*(-768 + nu*(3584 + 24*(a5 + (64*logu)/5.) - 123*pi2))*(-24*(a6 - (4*logu*(1751 + 756*nu))/105.)*(1536 + nu*(-3776 + 123*pi2)) + nu*(-2304*gsl_pow_int(a5 + (64*logu)/5.,2) + 96*(a5 + (64*logu)/5.)*(-3392 + 123*pi2) - (-3776 + 123*pi2)*(-32*(94 + 3*nu) + 123*pi2)))))/(13125.*gsl_pow_int(-768 + nu*(3584 + 24*(a5 + (64*logu)/5.) - 123*pi2),3)*u2);
    
    /* Second derivative of numerator and denominator */
    double d2Num = 2.*dN1 + d2N1*u;
    double d2Den = 2.*(D2 + dD1) + u*(6.*D3 + 4.*dD2 + d2D1) + u2*(12.*D4 + 6.*dD3 + d2D2) + u3*(20.*D5 + 8.*dD4 + d2D3) + u4*(10.*dD5 + d2D4) + u5*d2D5;
    
    /* Second derivative with respect of u */
    double d2A_u = prefactor*(2.*dDen*dDen*(*A) - 2.*dNum*dDen + Den*d2Num - d2Den*Num);

    *d2A = d2A_u;
    
    /* Second derivative with respect of r */
    /* *d2A = u4*d2A_u + 2.*u3*dA_u; */
    
  }
}

/** Tidal potential, two version implemented: 
    1. TEOB NNLO, Bernuzzi+ 1205.3403
    2. TEOBResum, Bini&Damour, 1409.6933, Bernuzzi+ 1412.4553 */
void eob_metric_Atidal(double r, Dynamics *dyn, double *AT, double *dAT, double *d2AT)
{
  
  double A, dA_u, d2A_u, dA, d2A;

  const double nu    = dyn->nu;
  const double rLR   = dyn->rLR_tidal;
  const double XA    = dyn->X1;
  const double XB    = dyn->X2;
  const double kapA2 = dyn->kapA2;
  const double kapB2 = dyn->kapB2;
  const double kapT2 = dyn->kapT2;
  const double kapT3 = dyn->kapT3;
  const double kapT4 = dyn->kapT4;

  /* Definition of the conservative tidal coefficients \bar{\alpha}_n^{(\ell)}, 
     Eq.(37) of Damour&Nagar, PRD 81, 084016 (2010) */
  const double bar_alph2_1 = dyn->bar_alph2_1;
  const double bar_alph2_2 = dyn->bar_alph2_2;
  const double bar_alph3_1 = dyn->bar_alph3_1;
  const double bar_alph3_2 = dyn->bar_alph3_2;

  const double p = dyn->pGSF_tidal; 
  
  /* shortcuts */
  double nu2  = nu*nu;
  double pi2  = Pi*Pi;
  double pi4  = pi2*pi2;
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

  if (dyn->use_tidal==TIDES_NNLO) {

    A    = -(kapT4*u10) - kapT2*u6*(1. + bar_alph2_1*u + bar_alph2_2*u2) - kapT3*u8*(1. + bar_alph3_1*u + bar_alph3_2*u2);
    dA_u = -10.*kapT4*u9 - kapT2*u6*(bar_alph2_1 + 2.*bar_alph2_2*u) - kapT3*u8*(bar_alph3_1 + 2.*bar_alph3_2*u)
      - 6.*kapT2*u5*(1. + bar_alph2_1*u + bar_alph2_2*u2) - 8.*kapT3*u7*(1. + bar_alph3_1*u + bar_alph3_2*u2);
    d2A_u = -90.*kapT4*u8
      - kapT2*(2*bar_alph2_2*u6 + 12.*u5*(bar_alph2_1 + 2*bar_alph2_2*u)
	     + 30.*u4*(1 + bar_alph2_1*u + bar_alph2_2*u2))
      - kapT3*(2.*bar_alph3_2*u8 + 16.*u7*(bar_alph3_1 + 2*bar_alph3_2*u) + 56.*u6*(1 + bar_alph3_1*u + bar_alph3_2*u2));
  
  } else { 

    const double c1  =  8.53353;
    const double c2  =  3.04309;
    const double n1  =  0.840058;
    const double d2  =  17.73239;

    double Acub   = 5./2.* u * (1. -  (c1+c2)*u +   c1*c2*u2);
    double dAcub  = 5./2.*     (1. -2*(c1+c2)*u + 3*c1*c2*u2);
    double d2Acub = 5    *     (   -  (c1+c2)   + 3*c1*c2*u);
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
    
    double d2f23  = 2*d2*(-1 + 3*d2*u2 + n1*u*(-3+d2*u2))*(Den*Den*Den);
    double d2A1SF = d2Acub*f23 + 2*dAcub*df23 + Acub*d2f23;
    double d2A2SF = 674./28.;
    double d2f0   = 6*(oom3u*oom3u*oom3u);
    double d2f1   = 0.25*(63*(rLR*rLR)*A1SF + 4*(-1+rLR*u)*(-7*rLR*dA1SF + (-1+rLR*u)*d2A1SF))*pow(oom3u,11./2.);
    double d2f2   = (rLR*p*((1+p)*rLR*A2SF - 2*(-1+rLR*u)*dA2SF +(-1.+rLR*u)*(-1.+rLR*u) *d2A2SF))*pow(oom3u,p+2);

    double AT2    = - kapA2*u6*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*u6*( f0 + XB*f1 + XB*XB*f2 );
    double AT3    = - kapT3*u8*(1. + bar_alph3_1*u + bar_alph3_2*u2);
    double AT4    = - kapT4*u10;

    double dAT2  = - kapA2*6.*u5*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*6.*u5*( f0 + XB*f1 + XB*XB*f2 ) - kapA2*u6*( df0 + XA*df1 + XA*XA*df2 ) - kapB2*u6*( df0 + XB*df1 + XB*XB*df2 );
    double dAT3  = - kapT3*(8.*u7 + 9*bar_alph3_1*u8 + 10*bar_alph3_2*u9);
    double dAT4  = - kapT4*10.*u9;

    double d2AT2  = - kapA2*30*u4*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*30*u4*( f0 + XB*f1 + XB*XB*f2 ) - 2*kapA2*6*u5*( df0 + XA*df1 + XA*XA*df2 ) - 2*kapB2*6*u5*( df0 + XB*df1 + XB*XB*df2 ) - kapA2*u6*( d2f0 + XA*d2f1 + XA*XA*d2f2 ) - kapB2*u6*( d2f0 + XB*d2f1 + XB*XB*d2f2 );
    
    double d2AT3  = - kapT3*(56*u6 + 72*bar_alph3_1*u7 + 90*bar_alph3_2*u8);
    double d2AT4  = - kapT4*90*u8;

    A     = AT2   + AT3   + AT4;
    dA_u  = dAT2  + dAT3  + dAT4;
    d2A_u = d2AT2 + d2AT3 + d2AT4;
    
  }

  *AT   = A;
  *dAT  = dA_u;
  *d2AT = d2A_u;
}

/** EOB Metric potentials A(r), B(r), and their derivatives, no spin version */
void eob_metric(double r, Dynamics *dyn, double *A, double *B, double *dA, double *d2A, double *dB)
{
  const double nu    = dyn->nu;
  const double u     = 1./r;
  const double u2    = u*u;
  const double u3    = u2*u;
  const double u4    = u2*u2;

  /* A potential and derivative with respect to u */  
  double Atmp, dAtmp_u, d2Atmp_u;
  eob_metric_A5PNlog(r, nu, &Atmp, &dAtmp_u, &d2Atmp_u);

  /* Add here tides if needed */
  if (dyn->use_tidal) {
    double AT, dAT_u, d2AT_u;
    eob_metric_Atidal(r, dyn, &AT, &dAT_u, &d2AT_u);
    Atmp     += AT;
    dAtmp_u  += dAT_u;
    d2Atmp_u += d2AT_u;
  }

  /* A potential and derivative with respect to r */  
  *A   = Atmp;
  *dA  = -dAtmp_u*u2;
  *d2A = 2.*dAtmp_u*u3 + d2Atmp_u*u4;

  /* D potential and derivative with respect to r */
  const double Dp  = 1.0 + 6.*nu*u2 - 2.*(3.0*nu-26.0)*nu*u3; // Pade' resummation of D
  const double D   = 1./Dp;
  const double dD  = 6.*u2*(2.*nu*u-(3.*nu-26.)*nu*u2)*D*D;

  /* B potential and derivative with respect to r */
  *B   = D/(Atmp);
  *dB  = (dD*(Atmp) - D*(*dA))/((Atmp)*(Atmp));

}
 
/** EOB Metric potentials A(r), B(r), and their derivatives, spin version */
void eob_metric_s(double r, Dynamics *dyn, double *A, double *B, double *dA, double *d2A, double *dB)
{

  const double nu    = dyn->nu;
  const double a1    = dyn->a1;
  const double a2    = dyn->a2;
  const double aK2   = dyn->aK2;
  const double C_Q1  = dyn->C_Q1;
  const double C_Q2  = dyn->C_Q2;
  const int usetidal = dyn->use_tidal;

  const double u   = 1./r;
  const double u2  = u*u;
  const double u3  = u2*u;
  const double u4  = u2*u2;
  
  double rc, drc, d2rc;
  eob_dyn_s_get_rc(r, nu, a1, a2, aK2, C_Q1, C_Q2, usetidal, &rc, &drc, &d2rc);

  /* A potential and derivative with respect to u */  
  double Aorb, dAorb_u, d2Aorb_u;
  eob_metric_A5PNlog(rc, nu, &Aorb, &dAorb_u, &d2Aorb_u);

  /* Add here tides if needed */
  if (usetidal) {
    double AT, dAT_u, d2AT_u;
    eob_metric_Atidal(rc, dyn, &AT, &dAT_u, &d2AT_u);
    Aorb     += AT;
    dAorb_u  += dAT_u;
    d2Aorb_u += d2AT_u;
  }

  /* A potential and derivative with respect to r */  
  double uc  = 1./rc;
  double uc2 = uc*uc;
  double uc3 = uc2*uc;
  double uc4 = uc2*uc2;

  double dAorb  = -dAorb_u*uc2;
  double d2Aorb = 2.*dAorb_u*uc3 + d2Aorb_u*uc4;

  /* Correct A for spin */
  double AKerr_Multipole = (1.+2.*uc)/(1.+2.*u);
  double fss = 1.;

  *A   = Aorb*AKerr_Multipole*fss;    
  *dA  = dAorb*drc*(1.+2.*uc)/(1.+2.*u) - 2.*Aorb*drc*uc2/(1.+2.*u) + 2.*Aorb*(1.+2.*uc)*u2/((1.+2.*u)*(1.+2.*u));
  *d2A = d2Aorb*(1.+2.*uc)/(1.+2.*u) + 4.*dAorb*( u2*(1.+2.*uc)/((1.+2.*u)*(1.+2.*u)) - uc2/(1.+2.*u)*drc) + Aorb*(-4.*u3*(1.+2.*uc)/((1.+2.*u)*(1.+2.*u)) + 8.*u4*(1.+2.*uc)/((1.+2.*u)*(1.+2.*u)*(1.+2.*u))+4.*uc3*(1.+2.*u)*drc*drc - 2.*uc2/(1.+2.*u)*d2rc);

  /* D potential and derivative with respect to r */
  double Dp = 1.0 + 6.*nu*uc2 - 2.*(3.0*nu-26.0)*nu*uc3; // Pade' resummation of D
  double D  = 1./Dp;
  double dD = 6.*uc2*(2.*nu*uc-(3.*nu-26.)*nu*uc2)*D*D;

  /* B potential and derivative with respect to r */
  *B   = r*r*uc2*D/(*A);
  *dB  = (dD*(*A) - D*(*dA))/((*A)*(*A));

}


