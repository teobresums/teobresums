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
#include <stdbool.h>

#include "constants.h"
#include "input_struc.h"
#include "s_A5PNlog.h"

void s_A5PNlog(
    double result[],             /** OUTPUT Dimension: 5*/
    double r,
    void *params,
    bool nnlo_flag){

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
    bool tidal_flag = (*(input *)params).tidal;
    double rLR      = (*(input *)params).rLR;
    double nu       = (*(input *)params).nu;
    double data[5];

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

    if (tidal_flag==true)
    {

        double A0    = A;
        double A0_du = dA_u;
   
        //Missing: b3NR (not needed), rlR (is calculated), kTl (yes, this has to be passed).
        
        // Tidal PN coefs
        double q   =  (1.+sqrt(1-4*nu)-2.*nu)/(2.*nu);
        double XA  =  0.5*(1.+sqrt(1.-4.*nu));
        double XB  =  1.-XA;
			   
        //dimensionless Love numbers (apsidal constants)
        double kAl[3];
        double kBl[3];

        kAl[0] = (*(input *)params).kAl1;
        kAl[1] = (*(input *)params).kAl2;
        kAl[2] = (*(input *)params).kAl3;
        
        kBl[0] = (*(input *)params).kBl1;
        kBl[1] = (*(input *)params).kBl2;
        kBl[2] = (*(input *)params).kBl3;
        
        //Compactness of the star
        double CA = (*(input *)params).CA;
        double CB = (*(input *)params).CB;
        
        //Computing the tidal coupling constants
        double kapA2 = 2. * kAl[0] * pow(XA/CA, 2.*2 +1.) * q; //Note: kap stands for kappa; see eqn(1) of REF
        double kapA3 = 2. * kAl[1] * pow(XA/CA, 2.*3 +1.) * q;
        double kapA4 = 2. * kAl[2] * pow(XA/CA, 2.*4 +1.) * q;
        
        double kapB2 = 2. * kBl[0] * pow(XB/CB, 2.*2 +1.) * q;
        double kapB3 = 2. * kBl[1] * pow(XB/CB, 2.*3 +1.) * q;
        double kapB4 = 2. * kBl[2] * pow(XB/CB, 2.*4 +1.) * q;
        
        double kapT2 = kapA2 + kapB2;
        double kapT3 = kapA3 + kapB3;
        double kapT4 = kapA4 + kapB4;

        double bar_alph2_1 = (5/2.*XA*kapA2 + 5/2.*XB*kapB2)/kapT2;
        double bar_alph2_2 = ((3.+XA/8.+ 337./28.*XA*XA)*kapA2 + (3.+XB/8.+ 337./28.*XB*XB)*kapB2)/kapT2; 
        double bar_alph3_1 = ((-2.+15./2.*XA)*kapA3 + (-2.+15./2.*XB)*kapB3)/kapT3;
        double bar_alph3_2 = ((8./3.-311./24.*XA+110./3.*XA*XA)*kapA3 + (8./3.-311./24.*XB+110./3.*XB*XB)*kapB3)/kapT3;

			       
        //case 'nnlo'
        if (nnlo_flag==true)
        { //Used for calculating the rLR
            //case 'nnlo'
            A     = -(kapT4*u10) - kapT2*u6*(1. + bar_alph2_1*u + bar_alph2_2*u2) - kapT3*u8*(1. + bar_alph3_1*u + bar_alph3_2*u2);
            dA_u = -10.*kapT4*u9 - kapT2*u6*(bar_alph2_1 + 2.*bar_alph2_2*u) - kapT3*u8*(bar_alph3_1 + 2.*bar_alph3_2*u)
	         - 6.*kapT2*u5*(1. + bar_alph2_1*u + bar_alph2_2*u2) - 8.*kapT3*u7*(1. + bar_alph3_1*u + bar_alph3_2*u2);
        }
        else
        { //Used for calculting the dynamcis
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
            double df23  = (n1 - 2.*d2*u - n1*d2*u2)*pow(DenI,2.);
            double dA1SF = dAcub*f23 + Acub*df23;
            double dA2SF = 674./28.*u;
            double df0   = 3.*u*(2.-rLR*u)*pow(oom3u,2.);
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

    result[0] = A;
    result[1] = dA;
    result[2] = d2A;
    result[3] = D;
    result[4] = dD;

}

