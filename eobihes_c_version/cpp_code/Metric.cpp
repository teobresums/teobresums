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

#include <math.h>
#include <vector>

#include "A_NumDenom.h"
#include "input_struc.h"

#include <gsl/gsl_math.h>
#include "Metric.h"

#include "constants.h"

using namespace::std;

vector<double> acoeffs(const double r, const double nu)
{
    
    vector<double> a(8);
    const double u    = 1./r;
    const double logu = log(u);

    
    /** Point-mass PN coefs */
    a[0] = 2.;
    
    a[1] = (94./3. - 41./32.*pi*pi);
    
    const double a5l  =  64./5.;
    const double a5c0 = -4237./60.+2275./512.*pi*pi+256./5.*Log2+128./5.*EulerGamma;
    const double a5c1 = -221./6.+41./32.*pi*pi;
    const double a5c  =  a5c0+nu*a5c1;
    a[2] = a5c;
    a[3] = a5l;
    a[4] = a5c+a5l*logu;
    
    const double a6c = 3097.3*nu*nu-1330.6*nu+81.38;
    const double a6l = (-7004./105.-144./5.*nu);
    a[5] = a6c;
    a[6] = a6l;
    a[7] = a6c+a6l*logu;
    
    return a;
}

vector<double> Metric(const double r, void *params, bool nnlo_flag)
{

    double nu         = (*(input *)params).nu;
    bool   tidal_flag = (*(input *)params).tidal;
    double rLR        = (*(input *)params).rLR;
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
    
    if (tidal_flag==true) {
        
        //Missing: b3NR (not needed), rlR (is calculated), kTl (yes, this has to be passed).
        
        // Tidal PN coefs
        double q   =  (1.+sqrt(1.-4.*nu)-2.*nu)/(2.*nu);
        //double XA  =  0.5*(1+sqrt(1.-4.*nu));
        //double XB  =  1.-XA;
	// CHECK IF THE SAME (BE AWARE ANOTHE DEF BELOW IN OTHER ROUTINE):
	double XA = (*(input *)params).X1;
	double XB = (*(input *)params).X2;
			   

        /** \bar{lambda} from Yagi et al. REF: */
        double lambdaAl2, lambdaAl3, lambdaAl4;
        double lambdaBl2, lambdaBl3, lambdaBl4;

        lambdaAl2 = (*(input *)params).LambdaAl2;
        lambdaAl3 = (*(input *)params).LambdaAl3;
        lambdaAl4 = (*(input *)params).LambdaAl4;
        
        lambdaBl2 = (*(input *)params).LambdaBl2;
        lambdaBl3 = (*(input *)params).LambdaBl3;
        lambdaBl4 = (*(input *)params).LambdaBl4; 
       
        /** Computing the tidal coupling constants */
        double kapA2 = 3.   * lambdaAl2 * pow(XA, 2.*2 +1.) * q; //Note: kap stands for kappa; see eqn(1) of REF
        double kapA3 = 15.  * lambdaAl3 * pow(XA, 2.*3 +1.) * q;
        double kapA4 = 105. * lambdaAl4 * pow(XA, 2.*4 +1.) * q;
 
        double kapB2 = 3.   * lambdaBl2 * pow(XB, 2.*2 +1.) / q;
        double kapB3 = 15.  * lambdaBl3 * pow(XB, 2.*3 +1.) / q;
        double kapB4 = 105. * lambdaBl4 * pow(XB, 2.*4 +1.) / q;
        
        double kapT2 = kapA2 + kapB2;
        double kapT3 = kapA3 + kapB3;
        double kapT4 = kapA4 + kapB4;

        double bar_alph2_1 = (5./2.*XA*kapA2 + 5./2.*XB*kapB2)/kapT2;
      	double bar_alph2_2 = ((3.+XA/8.+ 337./28.*XA*XA)*kapA2 + (3.+XB/8.+ 337./28.*XB*XB)*kapB2)/kapT2; 
      	double bar_alph3_1 = ((-2.+15./2.*XA)*kapA3 + (-2.+15./2.*XB)*kapB3)/kapT3;			     			   
      	double bar_alph3_2 = ((8./3.-311./24.*XA+110./3.*XA*XA)*kapA3 + (8./3.-311./24.*XB+110./3.*XB*XB)*kapB3)/kapT3;

			       
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
            double df23  = (n1 - 2.*d2*u - n1*d2*u2)*pow(DenI,2.);
            double dA1SF = dAcub*f23 + Acub*df23;
            double dA2SF = 674./28.*u;
            double df0   = 3.*u*(2.-rLR*u)*pow(oom3u,2.);
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

vector<double> A5pnP15_dd(const double r, void *params)
{
    
    double nu         = (*(input *)params).nu;
    bool   tidal_flag = (*(input *)params).tidal;
    double rLR        = (*(input *)params).rLR;

    
        
    vector<double> A_dd(2);
    
    /** Shorthands */
    const double u  = 1./r;
    const double u2 = u*u;
    const double u3 = u2*u;
    const double u4 = u3*u;
    const double u5 = u3*u2;
    const double u6 = u5*u;
    const double u7 = u6*u;
    const double u8 = u5*u3;        
    const double sm = nu;
    
    /** Point-mass PN coefs */
    const vector<double> a = acoeffs(r, nu);
    const double a3  = a[0];
    const double a4  = a[1];
    const double a5l = a[3];
    const double a5  = a[4];
    const double a6l = a[6];
    const double a6  = a[7];
    
    const vector<double> frac = A_NumDenom(r,a,nu);
    const double N  = frac[0];
    const double D  = frac[1];
    const double dN = frac[2];
    const double dD = frac[3];
    
    const double ooD = 1./D;
    
    const double d2N = -((sm*(-32. + (12.*a3 + 4.*a4 + a5 - 2.*a5l)*sm)*(a6l*(-32. + (12.*a3 + 4.*a4 + a5)*sm) + a5l*(-64. + sm*(4.*a4 - a6 + a3*(16. + a3*sm)))))/(gsl_pow_int(-32. + (12.*a3 + 4.*a4 + a5)*sm,3)*u));
    
    const double d2D = (sm*(-2048.*(8.*a4 + 4.*a5 + 6.*a5l + 2.*a6 + 3.*a6l) - (a6l*gsl_pow_int(-32 + (12.*a3 + 4.*a4 + a5)*sm,2) + 2.*(16384.*a3 - 64.*(224.*gsl_pow_int(a3,2) + (4.*a4 + a5)*(8.*a4 + 4.*a5 + 2.*a6 + 3.*a6l) + 4.*a3*(40.*a4 + 16.*a5 + 6.*a6 + 9.*a6l))*sm + (12.*a3 + 4.*a4 + a5)*(320.*gsl_pow_int(a3,2) + (4.*a4 + a5)*(8.*a4 + 4.*a5 + 2.*a6 + 3.*a6l) + 4.*a3*(40.*a4 + 16.*a5 + 6.*a6 + 9.*a6l))*gsl_pow_int(sm,2) - 2.*gsl_pow_int(a3,2)*gsl_pow_int(12.*a3 + 4.*a4 + a5,2)*gsl_pow_int(sm,3))*u + 2.*(12.*(4.*a4 + 2.*a5 + a6) + 10.*a6l + 3.*a3*(8.*a3 + 4.*a4 + a5)*sm)*gsl_pow_int(-32 + (12.*a3 + 4.*a4 + a5)*sm,2)*u2 + gsl_pow_int(-32 + (12.*a3 + 4.*a4 + a5)*sm,2)*(a3*sm*(-12.*a6 - 7.*a6l + 12.*gsl_pow_int(a3,2)*sm) + 12.*a5*(16 + a4*sm) + 8.*(12.*a6 + 7.*a6l + 6.*(2.*gsl_pow_int(a3,2) + 4.*a3*a4 + gsl_pow_int(a4,2))*sm))*u3 + gsl_pow_int(-32 + (12.*a3 + 4.*a4 + a5)*sm,2)*(144.*a6l + 20.*(2.*a4 + a5)*(8.*a3 + 2.*a4 + a5)*sm + (4.*a3 + a4)*sm*(-9.*a6l + 20.*gsl_pow_int(a3,2)*sm) - 20.*a6*(-16 + 4.*a3*sm + a4*sm))*u4 + a5l*((-32 + (12.*a3 + 4.*a4 + a5)*sm)*(-64 + sm*(4.*a4 - a6 - 2.*a6l + a3*(16 + a3*sm))) + 2.*sm*(-32.*(120.*a3 + 36.*a4 + 6.*a5 - 3.*a6 - 2.*a6l) + (480.*gsl_pow_int(a3,2) + (4.*a4 + a5)*(12.*a4 - 3.*a6 - 2.*a6l) + 12.*a3*(28.*a4 + 4.*a5 - 3.*a6 - 2.*a6l))*sm + 3.*gsl_pow_int(a3,2)*(12.*a3 + 4.*a4 + a5)*gsl_pow_int(sm,2))*u + 4.*(-32 + (12.*a3 + 4.*a4 + a5)*sm)*(-320 + sm*(20.*a4 - 5.*a6 - 2.*a6l + 5.*a3*(16 + a3*sm)))*u2 - (-8 + a3*sm)*(-32 + (12.*a3 + 4.*a4 + a5)*sm)*(-448 + sm*(28.*a4 - 7.*a6 - 2.*a6l + 7.*a3*(16 + a3*sm)))*u3 - sm*(-32 + (12.*a3 + 4.*a4 + a5)*sm)*(144.*(16.*a3 + 8.*a4 + 4.*a5 + a6) + 32.*a6l - (9.*(96.*gsl_pow_int(a3,2) + 12.*gsl_pow_int(a4,2) + 8.*a4*a5 + gsl_pow_int(a5,2) + a4*a6 + 4.*a3*(16.*a4 + 6.*a5 + a6)) + 2.*(4.*a3 + a4)*a6l)*sm + 9.*gsl_pow_int(a3,2)*(4.*a3 + a4)*gsl_pow_int(sm,2))*u4) + 2.*gsl_pow_int(a5l,2)*sm*(-64 + sm*(4.*a4 - a6 + a3*(16 + a3*sm)))*(-1 + u*(-2 + u*(-4 + u*(-8 + a3*sm + (-16 + 4.*a3*sm + a4*sm)*u)))))/u))/gsl_pow_int(-32 + (12.*a3 + 4.*a4 + a5)*sm,3);
    
    double A_ddu = (2.*N*dD*dD-2.*D*dD*dN-D*N*d2D+D*D*d2N)*ooD*ooD*ooD; // A_ddu
    //A_ddr=-2.*A5pnP15(r,nu,tidal_flag)[1]/r+A_dd[1]/(r*r*r*r); //A_ddr
    
    double d2A0_u = A_ddu;
    
    if (tidal_flag==true) {

        //Missing: b3NR (not needed), rlR (is calculated), kTl (yes, this has to be passed).

       /** Compactness of the star */
      
        /** Tidal PN coefs */
        double q    = (1.+sqrt(1.-4.*nu)-2.*nu)/(2.*nu);
        //double XA   =  0.5*(1.+sqrt(1.-4.*nu));
        //double XB   =  1.-XA;
	// CHECK IF THE SAME (BE AWARE ANOTHE DEF ABOVE IN OTHER ROUTINE):
	double XA = (*(input *)params).X1;
	double XB = (*(input *)params).X2;

        /** \bar{lambda} from Yagi et al. REF: */
        double lambdaAl2;
        double lambdaAl3;
        double lambdaAl4;

        double lambdaBl2;
        double lambdaBl3;
        double lambdaBl4;

        lambdaAl2 = (*(input *)params).LambdaAl2;
        lambdaAl3 = (*(input *)params).LambdaAl3;
        lambdaAl4 = (*(input *)params).LambdaAl4;

        lambdaBl2 = (*(input *)params).LambdaBl2;
        lambdaBl3 = (*(input *)params).LambdaBl3;
        lambdaBl4 = (*(input *)params).LambdaBl4;

        /** Computing the tidal coupling constants: Eq. (31) of D&N, PRD 81, 084016 (2010) */
       
        double kapA2 = 3.   * lambdaAl2 * pow(XA, 2.*2 +1.) * q; //Note: kap stands for kappa; see eqn(1) of REF
        double kapA3 = 15.  * lambdaAl3 * pow(XA, 2.*3 +1.) * q;
        double kapA4 = 105. * lambdaAl4 * pow(XA, 2.*4 +1.) * q;

        double kapB2 = 3.   * lambdaBl2 * pow(XB, 2.*2 +1.) / q;
        double kapB3 = 15.  * lambdaBl3 * pow(XB, 2.*3 +1.) / q;
        double kapB4 = 105. * lambdaBl4 * pow(XB, 2.*4 +1.) / q;
 
        double kapT2 = kapA2 + kapB2;
        double kapT3 = kapA3 + kapB3;
        double kapT4 = kapA4 + kapB4;

	//-----------------------------------------------------------------------------------
	// Definition of the conservative tidal coefficients \bar{\alpha}_n^{(\ell)}, Eq.(37)
	// of Damour&Nagar, PRD 81, 084016 (2010)
	//-----------------------------------------------------------------------------------
	double bar_alph2_1 = (5/2.*XA*kapA2 + 5/2.*XB*kapB2)/kapT2;
	double bar_alph2_2 = ((3.+XA/8.+ 337./28.*XA*XA)*kapA2 + (3.+XB/8.+ 337./28.*XB*XB)*kapB2)/kapT2; 
	double bar_alph3_1 = ((-2.+15./2.*XA)*kapA3 + (-2.+15./2.*XB)*kapB3)/kapT3;			     			   
	double bar_alph3_2 = ((8./3.-311./24.*XA+110./3.*XA*XA)*kapA3 + (8./3.-311./24.*XB+110./3.*XB*XB)*kapB3)/kapT3;
	                
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
        double df23   = (n1 - 2*d2*u - n1*d2*u2)*pow(Den,2);
        double A1SF   = Acub*f23;
        double dA1SF  = dAcub*f23 + Acub*df23;
        double A2SF   = 337./28.*u2;
        double dA2SF  = 674./28.*u;
        double oom3u  = 1./(1.-rLR*u);
        double f0     = 1 + 3*u2*oom3u;
        double f1     = A1SF *pow(oom3u,7./2.);
        double f2     = A2SF *pow(oom3u,p);
        double df0    = 3*u*(2.-rLR*u)*pow(oom3u,2);
        double df1    = 0.5*(7*rLR*A1SF + 2*(1.-rLR*u)*dA1SF)*pow(oom3u,9./2.);
        double df2    = (rLR*p*A2SF + (1.-rLR*u)*dA2SF)*pow(oom3u,p+1);
        
        double d2f23  = 2*d2*(-1 + 3*d2*u2 + n1*(-3+d2*u2))*pow(Den,3);
        double d2A1SF = d2Acub*f23 + 2*dAcub*df23 + Acub*d2f23;
        double d2A2SF = 674./28.;
        
        double d2f0   = 6*pow(oom3u,3);
        double d2f1   = 0.25*(63*pow(rLR,2)*A1SF + 4*(-1+rLR*u)*(-7*rLR*dA1SF + (-1+rLR*u)*d2A1SF))*pow(oom3u,11./2.);
        double d2f2   = (rLR*p*((1+p)*rLR*A2SF - 2*(-1+rLR*u)*dA2SF +pow(-1+rLR*u,2) *d2A2SF))*pow(oom3u,p+2);
        
        double d2AT2  = - kapA2*30*u4*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*30*u4*( f0 + XB*f1 + XB*XB*f2 ) - 2*kapA2*6*u5*( df0 + XA*df1 + XA*XA*df2 ) - 2*kapB2*6*u5*( df0 + XB*df1 + XB*XB*df2 ) - kapA2*u6*( d2f0 + XA*d2f1 + XA*XA*d2f2 ) - kapB2*u6*( d2f0 + XB*d2f1 + XB*XB*d2f2 );
        
        double d2AT3  = - kapT3*(56*u6 + 72*bar_alph3_1*u7 + 90*bar_alph3_2*u8);
        double d2AT4  = - kapT4*90*u8;
        double A_ddu  =   d2AT2 + d2AT3 + d2AT4;
        
        A_ddu=A_ddu+d2A0_u;
        
    }
    double A_ddr = -2.*Metric(r,params,false)[1]/r+A_ddu/(r*r*r*r);

    A_dd[1] = A_ddu;
    A_dd[0] = A_ddr;
    
return A_dd;

}


