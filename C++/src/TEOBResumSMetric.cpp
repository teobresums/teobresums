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
#include <ios>
#include <cmath>
#include <limits>
#include <vector>
#include <fstream>
#include <stdio.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include "TEOBResumS.h"

using namespace::std;

vector<double> acoeffs(const double r, const double nu)
{
  /** This function precomputes the  coefficients of the PN-expanded EOB interaction potential,
       A(u;nu). See Eq. (1) of Nagar, Damour, Reisswig, Pollney, PRD 93 (2016), 044046 [NDRP]
       We give here the analyticall known 2PN, 3PN and 4PN coefficients plus the effective
       5PN coefficient informed by nonspinning NR simulations in NDRP */
  
    
    vector<double> a(8);
    const double u    = 1./r;
    const double logu = log(u);
    
    
    /** Point-mass PN coefs */
    /** 2PN coefficient */
    a[0] = 2.;
    /** 3PN coefficient*/
    a[1] = (94./3. - 41./32.*pi*pi);

    /** 4PN coefficient */
    const double a5l  =  64./5.;
    const double a5c0 = -4237./60.+2275./512.*pi*pi+256./5.*Log2+128./5.*EulerGamma;
    const double a5c1 = -221./6.+41./32.*pi*pi;
    const double a5c  =  a5c0+nu*a5c1;
    a[2] = a5c;
    a[3] = a5l;
    a[4] = a5c+a5l*logu;

    /** effective 5PN coefficient. The fit comes from Eq.(5) of NDRP */
    const double a6c = 3097.3*nu*nu-1330.6*nu+81.38;
    const double a6l = (-7004./105.-144./5.*nu);
    a[5] = a6c;
    a[6] = a6l;
    a[7] = a6c+a6l*logu;
    
    return a;
}

vector<double> Metric(const double r, void *params, bool nnlo_flag)
{

  /** Computation of the effective metric, aka as EOB potentials.
     This functions gives the A function as a (1,5) Pade' approximant
     as well as the D function as a (0,3) Pade' approximant.
     Derivatives of A are also analytically computed here. 
     These derivatives are used in the EoM */

  double nu         = (*(TEOBResumParams *)params).nu;
  bool   tidal_flag = (*(TEOBResumParams *)params).flags.tidal;
  double rLR        = (*(TEOBResumParams *)params).rLR;
  vector<double> data(5);\

    
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
      /*The tidal part of the potential stars here */
    {

      /* kapA2 stands for the \kappa^A_\ell of Eq.(26) in Damour & Nagar, PRD 81, 084016*/
        double XA = (*(TEOBResumParams *)params).X1;
        double XB = (*(TEOBResumParams *)params).X2;
        
        /** Computing the tidal coupling constants */
        double kapA2 = (*(TEOBResumParams *)params).kappaAl2; // 3.*lambdaAl2 * pow(XA, 2.*2 +1.) / q; 
        double kapB2 = (*(TEOBResumParams *)params).kappaBl2; // 3.*lambdaBl2 * pow(XB, 2.*2 +1.) * q;
        
        double kapT2 = (*(TEOBResumParams *)params).kappaTl2; //kapA2 + kapB2;
        double kapT3 = (*(TEOBResumParams *)params).kappaTl3; //kapA3 + kapB3;
        double kapT4 = (*(TEOBResumParams *)params).kappaTl4; //kapA4 + kapB4;
        
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
        {   /*Used for calculating the dynamics
              The l=2 part of the potential  Bini & Damour, Eq.(7.33) of  PRD 90, 124037 (2014),
	      with the free coefficent p fixed to p=4
	      The l=3 and l=4 corrections are given in the simple (nonresummed) nnlo form.*/
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
            
            /** the tidal part of the potential */
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

vector<double> A5pnP15_dd(const double r, void *params)
{
  /** The A potential. This routine computes explcitly the (1,5) 
      Pade' approximant of the A function */

    double nu         = (*(TEOBResumParams *)params).nu;
    bool   tidal_flag = (*(TEOBResumParams *)params).flags.tidal;
    double rLR        = (*(TEOBResumParams *)params).rLR;
        
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
        
                
        double XA = (*(TEOBResumParams *)params).X1;
        double XB = (*(TEOBResumParams *)params).X2;
        
        //-----------------------------------------------------------------------------------
        // Definition of the conservative tidal coefficients \bar{\alpha}_n^{(\ell)}, Eq.(37)
        // of Damour&Nagar, PRD 81, 084016 (2010)
        //-----------------------------------------------------------------------------------
        double kapA2 = (*(TEOBResumParams *)params).kappaAl2; // 3.   * lambdaAl2 * pow(XA, 2.*2 +1.) / q; 
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
        
        double d2f23  = 2*d2*(-1 + 3*d2*u2 + n1*u*(-3+d2*u2))*(Den*Den*Den);
        double d2A1SF = d2Acub*f23 + 2*dAcub*df23 + Acub*d2f23;
        double d2A2SF = 674./28.;
        
        double d2f0   = 6*(oom3u*oom3u*oom3u);
        double d2f1   = 0.25*(63*(rLR*rLR)*A1SF + 4*(-1+rLR*u)*(-7*rLR*dA1SF + (-1+rLR*u)*d2A1SF))*pow(oom3u,11./2.);
        double d2f2   = (rLR*p*((1+p)*rLR*A2SF - 2*(-1+rLR*u)*dA2SF +(-1.+rLR*u)*(-1.+rLR*u) *d2A2SF))*pow(oom3u,p+2);
        
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

int rhs(double t, const double y[], double f[], void *params){
    (void)(t); /* avoid unused parameter warning */

    /** This routine gives the RHS of the EOB equation of motion for the nonspinning case. 
	These RHS are written explicitly in Eq.(6a)-(6d) of Damour, Nagar and Bernuzzi, 
	PRD 87, 084035 (2013) [DNB]. Note however that the radial flux F_r*=0 here.
	Here we work with rescaled variables, see e.g. Sec IIIA in https://arxiv.org/pdf/1406.6913.pdf”*/
    
    double nu = (*(TEOBResumParams *)params).nu;
    
    const double z3 = 2.0*nu*(4.0-3.0*nu);
    const double prstar  = y[2];
    const double prstar2 = prstar*prstar;
    const double prstar3 = prstar2*prstar;
    const double prstar4 = prstar3*prstar;
    
    const double pphi  = y[3];
    const double pphi2 = pphi*pphi;
    
    const double r  = y[0];
    const double r2 = r*r;
    const double u  = 1./r;
    const double u2 = u*u;
    const double u3 = u2*u;
    
    const vector<double> metric = Metric(r, params, false);
    const double A     = metric[0];
    const double dA    = metric[1];
    const double B     = metric[3];
    const double dB    = metric[4];
    const double one_A = 1./A;
    const double one_B = 1./B;
    
    //const double j = pphi;
    
    double Heff = sqrt(prstar2+A*(1. + pphi2*u2 +  z3*prstar4*u2));
    const double H = sqrt( 1. + 2.*nu*(Heff - 1.) )/nu;
    const double E = H*nu;
    
    const double sqrAB = sqrt(A/B);
    
    //r evol eqn rhs
    f[0]  = (prstar+2.0*z3*A*prstar3*u2)/Heff;
    f[0] *= (sqrAB / E);
    
    //phi evol eqn rhs
    double Omega   = A*pphi*u2/Heff;
    Omega *= 1./E;
    f[1] = Omega;
    
    //prstar evol eqn rhs
    f[2]  = (dA + ( pphi2 + z3*prstar4 )*(dA*u2-2.0*A*u3))/Heff;
    f[2] *= - 0.5 * sqrAB / E;
    
    //pphi evol eqn rhs
    const double sqrW = sqrt(A*(1. + pphi2*u2));

    /* computation of psi in the nonspinning limit, See Eq.(69) of
       PRD 90, 044018 (2014) and references therein, notably
       Damour & Gopakumar, Phys. Rev. D73 (2006) 124006 */
    double psi = 0.;
    psi = 2.*(1.0 + 2.0*nu*(sqrW - 1.0))/(r2*dA);
    
    /* The code below is not used. The standard choice of psi
       is the one above*/
    /*bool psi_flag = false; // flag for an alternative way of computing psi
    if (psi_flag==false) {
        psi = 2.*(1.0 + 2.0*nu*(sqrW - 1.0))/(r2*dA);
    } if (psi_flag==true) {
        psi = 2.*(1.0 + 2.0*nu*(Heff - 1.0))/(r2*dA);
    }*/

    const double r_omega = r*cbrt(psi);
    const double v_phi   = r_omega*Omega;
    const double x       = v_phi * v_phi;
    const double jhat    = pphi/(r_omega*v_phi);
    
    const double E2       = E*E;
    const double tmpE     = 1./Heff+nu/E2;
    const double denE     = E*Heff;
    const double one_denE = 1./denE;
    
    const double dHeff_dr      = 0.5*(dA + (pphi2 + z3*prstar4)*(dA*u2 - 2.*A*u3))/Heff;
    const double dHeff_dprstar = (prstar + z3*2.0*A*u2*prstar3)/Heff;
    const double dprstar_dt    = f[2];
    const double dr_dt         = f[0];
    const double ddotr_dr      = sqrAB*( (prstar + z3*2.*A*u2*prstar3)*(0.5*(dA*one_A-dB*one_B)-dHeff_dr*tmpE)+ 2.0*z3*(dA*u2 - 2.*A*u3)*prstar3)*one_denE;
    const double ddotr_dprstar = sqrAB*( 1.+z3*6.*A*u2*prstar2-(prstar + z3*2.*A*u2*prstar3)*dHeff_dprstar*tmpE)*one_denE;
    
    
    /* Approximate ddot(r) without the flux. See Appendix of Damour, Nagar and Bernuzzi (DNB), 
       arXiv:1212.4357 where this approximation is explicitly introduced and discussed */
    const double ddotr = dprstar_dt*ddotr_dprstar + dr_dt*ddotr_dr;

    /* Vector of souces for the multipolar flux, Eq.(16) of DNB */
    double source[] = {
        jhat,Heff,
        Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,
        Heff,jhat,Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,jhat,Heff,
        Heff,jhat,Heff,jhat,Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,jhat,Heff,jhat,Heff};
    
    f[3] = flux(x,Omega,r_omega,E,Heff,jhat,r, prstar,ddotr,source,params);
    
    
    return GSL_SUCCESS;
}

int s_RHS(double t, const double y[], double f[], void *params)
{
    /*
      This function provides the rhs of the full EOB equations for
      spinning binaries. The reference Hamiltonian is given in 
      Eq. (24) and following in Damour&Nagar, PRD 90, 044018.
      Due to the simple stracture of the effective Hamiltonian, that
      is written as

      H_eff = H_SO + H_orb

      writing this r.h.s. essentially amounts to adding to the orbital
      pieces written in rhs above the derivatives of the spin-orbit term.
      Spin-spin terms do not appear explictly here because they are
      incorporated in special resummed form within H_SO and H_orb.

     */
    
    /* y content
     %---------------------------
     %      y(0) = r
     %      y(1) = phi
     %      y(2) = pr*
     %      y(3) = pphi
     %---------------------------
     */
    
    (void)(t); /* avoid unused parameter warning */
    
    double nu       = (*(TEOBResumParams *)params).nu;
    double S        = (*(TEOBResumParams *)params).S;
    double Sstar    = (*(TEOBResumParams *)params).Sstar;
    double chi1     = (*(TEOBResumParams *)params).chi1;
    double chi2     = (*(TEOBResumParams *)params).chi2;
    double X1       = (*(TEOBResumParams *)params).X1;
    double X2       = (*(TEOBResumParams *)params).X2;
    double c3       = (*(TEOBResumParams *)params).cN3LO;
    bool tidal_flag = (*(TEOBResumParams *)params).flags.tidal;
    double aK2      = (*(TEOBResumParams *)params).aK2;
    
    
    double r      = y[0];
    double prstar = y[2];
    double pph    = y[3];
    double pphi2  = pph*pph;
    
    // shorthands
    double prstar2 = prstar*prstar;
    double prstar3 = prstar2*prstar;
    double prstar4 = prstar3*prstar;
    
    double z3 = 2.*nu*(4.-3.*nu);
    
    double A, B, dA;
    vector<double> metric;

    /** Call to metric functions. The centrifugal radius, redefined below, is also called inside.
        The tidal parameters are also defined within this routine */
    metric = s_Metric(r, params, false); //{A,B,dA,d2A} data[0]=A; data[1]=A_dr; data[2]=A_du; data[3]=B; data[4]=B_dr;
    A      = metric[0];
    B      = metric[1];
    dA     = metric[2];
    
    /** Introduce here the "centrifugal radius", Eq.(58) of Damour&Nagar, PRD 90, 044018 (2014)
        Spin-spin terms are all incorporated within this particular variable */
    vector<double> rc_vec;
    rc_vec = s_get_rc(r,params); //[rc, drc, d2rc]
    double rc     = rc_vec[0];
    double drc_dr = rc_vec[1];
    double uc     = 1./rc;
    double uc2    = uc*uc;
    double uc3    = uc2*uc;
    
    double Heff_orb = sqrt( prstar2+A*(1. + pphi2*uc2 +  z3*prstar4*uc2) );
    
    vector<double> ggm = s_GS(r, rc, drc_dr, aK2, prstar, pph, nu, chi1, chi2, X1, X2, c3);
    
    double GS              = ggm[2];
    double GSs             = ggm[3];
    double dGS_dprstar     = ggm[4];
    double dGSs_dprstar    = ggm[5];
    double dGS_dr          = ggm[6];
    double dGSs_dr         = ggm[7];
    double dGSs_dpph       = ggm[9];
    double d2GS_dprstar20  = ggm[12];
    double d2GSs_dprstar20 = ggm[13];
    
    
    double Heff  = Heff_orb + (GS*S + GSs*Sstar)*pph;
    double H     = sqrt( 1. + 2.*nu*(Heff - 1.) );
    double one_H = 1./H;
    
    double sqrtAbyB = sqrt(A/B);
    
    double dHeff_dr = pph*(dGS_dr*S + dGSs_dr*Sstar) + 1./(2.*Heff_orb)*( dA*(1. + pphi2*uc2 + z3*prstar4*uc2) - 2.*A*uc3*drc_dr*(pphi2 + z3*prstar4) );
    
    double dp_rstar_dt_0 = - sqrtAbyB*one_H*dHeff_dr;
    
    double dHeff_dprstar = pph*(dGS_dprstar*S + dGSs_dprstar*Sstar) + (prstar/Heff_orb)*(1. + 2.*A*uc2*z3*prstar2);
    
    // second derivative of Heff wrt to pr_star neglecting all pr_star^2 terms
    double d2Heff_dprstar20 = pph*(d2GS_dprstar20*S + d2GSs_dprstar20*Sstar) +  (1./Heff_orb)*(1. + 2.*A*uc2*z3*prstar2);
    
    double ddotr_dp_rstar = sqrtAbyB*one_H*d2Heff_dprstar20;
    
    //-------------------------------------------
    // 0.th -- approximate ddot(r)_0 without Fphi
    //-------------------------------------------
    double ddotr = dp_rstar_dt_0*ddotr_dp_rstar; //order pr_star^2 neglected
    
    
    //------------------ dr/dt ------------------
    f[0] = sqrtAbyB*one_H*dHeff_dprstar;
    
    //----------------- d\phi/dt ----------------
    double dHeff_dpph = GS*S + (GSs + pph*dGSs_dpph)*Sstar + pph*A*uc2/Heff_orb;
    f[1] = one_H*dHeff_dpph;
    
    //----------------- dp_{r*}/dt --------------
    f[2] = -sqrtAbyB*one_H*dHeff_dr;
    
    //------------------ dp_{\phi}/dt -----------
    double Omg = f[1];
    
    
    /*----------------------------------------------------
      Compute here the new r_omg radius as in Eq. (70)
      of Damour&Nagar, arXiv:1406.6913. One starts by
      computing some quantities with prstar=0, consitently
      with the nonspinning case. See Eq.(71) of the above
      reference.
      ----------------------------------------------------*/
    vector<double> ggm0 = s_GS(r, rc, drc_dr, aK2, 0., pph, nu, chi1, chi2, X1, X2, c3);
    
    double GS_0      = ggm0[2];
    double GSs_0     = ggm0[3];
    double dGS_dr_0  = ggm0[6];
    double dGSs_dr_0 = ggm0[7];
    double Heff_orb_0 = sqrt(A*(1.0 + pphi2*uc2));   // circularized effective Hamiltonian H_0^eff
    double Heff_0     = Heff_orb_0 + (GS_0*S + GSs_0*Sstar)*pph;
    double H0         = sqrt(1.0 + 2.0*nu*(Heff_0 - 1.0) );
    double one_H0     = 1./H0;
    double Gtilde     = GS_0*S     + GSs_0*Sstar;
    double dGtilde_dr = dGS_dr_0*S + dGSs_dr_0*Sstar;
    double duc_dr     = -uc2*drc_dr;
    double psic       = (duc_dr + dGtilde_dr*rc*sqrt(A/pphi2 + A*uc2)/A)/(-0.5*dA);
    double r_omg      = pow( ((1./sqrt(rc*rc*rc*psic))+Gtilde)*one_H0, -2./3. );
    double v_phi = r_omg*Omg;
    double x     = v_phi*v_phi;
    double jhat  = pph/(r_omg*v_phi);
    double Fphi  = s_Flux(x,Omg,r_omg,H,Heff,jhat,r,prstar,ddotr,params);
    f[3] = Fphi;
    
    return GSL_SUCCESS;
}

vector<double> s_A5PNlog(double r, void *params, bool nnlo_flag){
    
    /*  This function computes the Pade' (1,5) resummed A function 
        with the 5PN NR-informed coefficient.
        The function is defined (for example) in Damour,Nagar PRD 90, 044018 (2014). The most up to date value of a6c can be found here: https://arxiv.org/abs/1506.08457, while the c3 value can be found here: https://arxiv.org/abs/1703.06814. */
    
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
        } else { /*Used for calculting the dynamics
                   case 'nnlo_gsfLR'; Bini & Damour, 1409.6933 + free light-ring
                   Tidal PN coefs*/
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
     %                       first and second derivatives of A and B.
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

double c3_fit_global(double nu, double chi1, double chi2, double X1, double X2, double a1, double a2, bool tidal_flag)

/** This function computed the next-to-next-to-next-to-leading order effective spin-orbit parameter informed by
   SXS NR simulations. Reference is Nagar, Riemenschneider and Pratten, PRD 96, 084045 (2017) mentioned below.
   Here, the spin variables are a1 = X1 chi1 and a2 = X2 chi2 that are addressed as tilde{a}_1 and tilde{a}_2
   in the reference paper and are called at1 and at2 in other parts of the code. This will need cleaning and
   uniformity when porting the code to LAL.*/
{
    
    double c3 = 0.;
    if (tidal_flag==true)
      /* NO effective spin-orbit calibration in the BNS case. This is doen on purpose because the BBH calibration
         of c3 effectively incorporates effects, like the "missing" spin-spin information, that is different,
         and EOS dependent, in the BNS case. So the NR-calibrated parameter is set to zero to avoid
         systematics */
    {
        c3 = 0.;
    }
    else
    {
      /* Effective NNNLO spin-orbit parameter informed by NR simulations.
         The functional form is given in Eq.(26) of Nagar, Riemenschneider
         and Pratten, PRD 96, 084045 (2017), while the numerical values are
         found in Eqs.(27)-(34) of the same reference */

      
        double nu2 = nu*nu;
        double nu3 = nu2*nu;

	/*------------------------------------
          equal-mass, equal-spin coefficients
	  ------------------------------------*/
	// OLD fits with the incorrect treatment of the l=5, m=odd modes
	// Before 10/05/2018
	// -- TO BE DELETED EVENTUALLY
	/*
        double c0 =  44.822889;
        double n1 =  -1.879350;
        double n2 =   0.894242;
        double d1 =  -0.797702;
	*/
	/**********************************************************/
	/* New coefficient with the correct l=5 modes: 10/05/2018 */
	/**********************************************************/
	double c0 =  43.371638;
        double n1 =  -1.174839;
        double n2 =   0.354064;
        double d1 =  -0.151961;
	
        double c3_eq = c0*(1. + n1*(a1+a2) + n2*(a1+a2)*(a1+a2))/(1.+d1*(a1+a2));

	/*----------------------------------------
          unequal-mass, unequal-spin coefficients
	  ----------------------------------------*/
	// OLD fits with the incorrect treatment of the l=5, m=odd modes
	// Before 10/05/2018
	// -- TO BE DELETED EVENTUALLY
	/*
        double cnu    = 1222.36;
        double cnu2   = -12764.4;
        double cnu3   =  36689.6;
        double ca1_a2 = -358.086;
	*/
	/**********************************************************/
	/* New coefficient with the correct l=5 modes: 10/05/2018 */
	/**********************************************************/
	double cnu    =  929.579;
        double cnu2   = -9178.87;
        double cnu3   =  23632.3;
        double ca1_a2 = -104.891;
	        
        double c3_uneq = cnu*(a1+a2)*nu*sqrt(1.-4.*nu) + cnu2*(a1+a2)*nu2*sqrt(1.-4.*nu) + cnu3*(a1+a2)*nu3*sqrt(1.-4.*nu) + ca1_a2*(a1-a2)*nu2;
        
        c3 = c3_eq + c3_uneq;
    }
    
    return c3;
}

vector<double> s_GS(double r, double rc, double drc_dr, double aK2, double prstar, double pph, double nu, double chi1, double chi2, double X1, double X2, double cN3LO)
{
    
    /*
      EOB_GetGSs(r,aK2,prstar,nu). This function computes the
      gyro-gravitomagnetic functions GS and GS*, that are called GS and GSs.
      The functions (and coefficients) are written explicitly in Sec.IIIB of
      Damour & Nagar, PRD 90, 044018 (2014)
          
     */
    
    
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
    double dGS0_duc  = 2.*u2/drc_dr + 4.*u*uc;
    
    double GSs0          =  3./2.*uc3;
    double dGSs0_duc     =  9./2.*uc2;
    double dGSs0_dprstar =  0.0;
    double dGSs0_dpph    =  0.0;
    
    // coefficients of hat{GS}
    double c10 =  5./16.*nu;
    double c20 =  51./8.*nu + 41./256.*nu2;
    double c30 =  nu*cN3LO;
    double c02 =  27./16.*nu;
    double c12 =  12.*nu - 49./128.*nu2;
    double c04 = -5./16.*nu + 169./256.*nu2;
    
    // coefficients of hat{GS*}
    // Could be precomputed and put in a matrix
    double cs10 = 3./4.   + nu/2.;
    double cs20 = 27./16. + 29./4.*nu + 3./8.*nu2;
    double cs02 = 5./4.   + 3./2.*nu;
    double cs12 = 4.   + 11.*nu     - 7./8.*nu2;
    double cs04 = 5./48.  + 25./12.*nu + 3./8.*nu2;
    double cs30 = nu*cN3LO + 135./32.;
    double cs40 = 2835./256.;
    
    double hGS  =  1./(1.  + c10*uc + c20*uc2 + c30*uc3 + c02*prstar2 + c12*uc*prstar2 + c04*prstar4); //write
    
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
    double dhGS_duc  = -hGS*hGS*(  c10 + 2.*c20*uc  + 3.*c30*uc2);
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
    double dGS_dprstarbyprstar  = -2.*GS0*hGS*hGS *( c02  +  c12*uc +  2.*c04*prstar2);
    double dGSs_dprstarbyprstar = -2.*GSs0*hGSs*hGSs*(cs02 + cs12*uc + 2.*cs04*prstar2);
    
    // --- for NQC --
    // Second derivatives neglecting all pr_star^2 terms
    double d2GS_dprstar20  =  GS0*(-2.*hGS*hGS *( c02 +  c12*uc +  2.*c04*prstar2));
    double d2GSs_dprstar20 =  GSs0*(-2.*hGSs*hGSs*(cs02 + cs12*uc + 2.*cs04*prstar2));
    
    return {hGS,hGSs,GS,GSs,dGS_dprstar,dGSs_dprstar,dGS_dr,dGSs_dr,dGS_dpph,dGSs_dpph,dGS_dprstarbyprstar,dGSs_dprstarbyprstar,d2GS_dprstar20,d2GSs_dprstar20};
}

vector <double> s_get_rc(double r, void *params)
{
    /** The centrifugal radius, introduced in eq. 8-22 of Damour,Nagar PRD 90, 044018 (2014), is needed to factorize out the metric A potential.
        This function computes the centrifugal radius in the spinning case. */
    
    
    double nu         = (*(TEOBResumParams *)params).nu;
    double at1        = (*(TEOBResumParams *)params).a1;
    double at2        = (*(TEOBResumParams *)params).a2;
    bool   tidal_flag = (*(TEOBResumParams *)params).flags.tidal;
    double aK2        = (*(TEOBResumParams *)params).aK2;
    
    double C_Q1       = (*(TEOBResumParams *)params).C_Q1;
    double C_Q2       = (*(TEOBResumParams *)params).C_Q2;
    double C_Oct1     = (*(TEOBResumParams *)params).C_Oct1;
    double C_Oct2     = (*(TEOBResumParams *)params).C_Oct2;
    double C_Hex1     = (*(TEOBResumParams *)params).C_Hex1;
    double C_Hex2     = (*(TEOBResumParams *)params).C_Hex2;
    
    double rc, drc_dr, d2rc_dr2;
    
    double u   = 1./r;
    double u2  = u*u;
    double u3  = u*u2;
    double u4  = u*u3;
    double r2  = r*r;
    double X12 = sqrt(1.-4.*nu);   //(X1-X2) will be defined at the beginning and not redefined several times
    
    if (tidal_flag==true)
    {
        
      /* Inclusion of LO spin-square coupling. The S1*S1 term coincides with the BBH one, no effect of structure.
         The self-spin couplings, S1*S1 and S2*S2 get a EOS-dependent coefficient, CQ, that describe the quadrupole
         deformation due to spin. Notation of Levi-Steinhoff, JCAP 1412 (2014), no.12, 003. Notation analogous to
         the parameter a of Poisson, PRD 57, (1998) 5287-5290 or C_ES^2 in Porto & Rothstein, PRD 78 (2008), 044013

	 NS quadrupoles due to rotation
         These are parameters that should be specified in the parameter file
         
         Inclusion of LO spin-square coupling. The S1*S1 term coincides with the BBH one, no effect of structure.
         The self-spin couplings, S1*S1 and S2*S2 get a EOS-dependent coefficient, CQ, that describe the quadrupole
         deformation due to spin. Notation of Levi-Steinhoff, JCAP 1412 (2014), no.12, 003. Notation analogous to
         the parameter a of Poisson, PRD 57, (1998) 5287-5290 or C_ES^2 in Porto & Rothstein, PRD 78 (2008), 044013
                           
         The implementation uses the I-Love-Q fits of Table I of Yunes-Yagi paper, PRD 88, 023009, the 
         bar{Q}(bar{\lambda)^{tid}) relation, line 3 of the table. The dimensionless bar{\lambda} Love number 
         is related to our apsidal constant as lambda = 2/3 k2/(C^5) so that both quantities have to appear here.*/
        
        //BNS effective spin parameter
        double a02      = C_Q1*at1*at1 + 2.*at1*at2 + C_Q2*at2*at2;
        
	double delta_a2 = X12*(at1*at1*(C_Q1+0.25) - at2*at2*(C_Q2+0.25))
	                + at1*at1*(-17./4.+3.*C_Q1-0.5*nu)
	                + at2*at2*(-17./4.+3.*C_Q2-0.5*nu)
	                + at1*at2*(nu-2.0);

        double delta_a2_nnlo  =
	    (  387./28.  - 207./28.*nu)*(C_Q1*at1*at1 + C_Q2*at2*at2)
          + (-2171./212. - 269./28.*nu + 0.375*nu*nu)*(at1*at1 + at2*at2)
	  + (- 25./2    - 145./8.*nu - 0.75 *nu*nu)*at1*at2
	  + 163./28.                               *X12*(C_Q1*at1*at1 - C_Q2*at2*at2)
	  + (  -29./112. - 2.625   *nu ) *X12*(at1*at1 - at2*at2);

	double delta_a4_lo = 0.75*(C_Hex1 - C_Q1*C_Q1)*at1*at1*at1*at1
	                     + 3.*(C_Oct1 - C_Q1)     *at1*at1*at1*at2
	                     + 3.*(C_Q1*C_Q2 - 1)     *at1*at1*at2*at2
	                     + 3.*(C_Oct2 - C_Q2)     *at1*at2*at2*at2
	                   + 0.75*(C_Hex2 - C_Q2*C_Q2)*at2*at2*at2*at2;

	delta_a4_lo = 0.;
	
	double alphanu2 = 1. + 0.5/a02*delta_a2;
        
	double rc2   = r2 + a02*(1. + 2.*alphanu2*u) + (delta_a2_nnlo+delta_a4_lo)*u2;
	rc          = sqrt(rc2);
	drc_dr      = r/rc*(1. -alphanu2*a02*u3 - (delta_a2_nnlo+delta_a4_lo)*u4);
	d2rc_dr2    = 1./rc*(1. -drc_dr*drc_dr +2.*alphanu2*a02*u3 +3.*(delta_a2_nnlo+delta_a4_lo)*u4);

	//NO spin-spin-tidal couplings
        /*double rc2 = r2;
         rc = r;
         drc_dr = 1;
         d2rc_dr2 = 0;*/
        
    }
    else
    {
      /* NLO spin-spin coupling introduced via the centrifugal radius, see Damour & Nagar, PRD 90, 044018 (2014).
         The quantity here addressed as u*c_ss_nlo corresponds to the \delta a^2 defined in Eq.(59) of this paper
         and then detailed in Eqs.(60) to (65). One actually finds that these complicates equations simplify when
         they are expressed in terms of the total-mass scaled Kerr parameters of the two black holes, 
         at1 = X1 chi1 and at2 = X2 chi2. This calculation is also mentioned in the Conclusions section
         of arXiv:1801.02366.
          
        In the first version of the code, the NLO spin-orbit correction was written with the following 
        piece of code:
     
        double alphanu2 = 1. + 0.5/aK2*(- at2*at2*(5./4. + 5./4.*X12 + nu/2.) - at1*at1*(5./4. - 5./4.*X12 +nu/2.) + at1*at2*(-2.+nu));	
	double rc2      = r2 + aK2*(1. + 2.*alphanu2/r);
	rc         = sqrt(rc2);
	drc_dr     = r/rc*(1.+aK2*(-alphanu2*u3 ));
	d2rc_dr2   = 1./rc*(1.-drc_dr*r/rc*(1.-alphanu2*aK2*u3)+ 2.*alphanu2*aK2*u3);

        This was actually formally singular and creating nans because of the 0.5/aK2, that actually simplifies
        with the aK2 at numerator giving a regular expression. This was corrected on 25/01/2018 explicitly 
        writing the expression above */
	  
      double c_ss_nlo = (- at2*at2*(1.25 + 1.25*X12 + 0.5*nu) - at1*at1*(1.25 - 1.25*X12 + 0.5*nu) + at1*at2*(-2.+nu));
      double rc2      = r2 + aK2*(1. + 2.*u) + u*c_ss_nlo;
	  	
      rc         = sqrt(rc2);
      drc_dr     = r/rc*(1-(aK2 + 0.5*c_ss_nlo)*u3);	
      d2rc_dr2   = 1./rc*(1.-drc_dr*r/rc*(1.-(aK2+0.5*c_ss_nlo)*u3)+ (2.*aK2 + c_ss_nlo)*u3);
    }
    
    return {rc, drc_dr, d2rc_dr2};
}

double DHeff0(double x, void *DHeff_params){

  /* This function computed the radial derivative of the circularized effective 
     Hamiltonian (pr*=0). This is used in the computation of the derivatives of
     the gyro-gravitomagnetic ratios */
  
    struct energy_params *p
    = (struct energy_params *) DHeff_params;
    
    double rorb   = p->rorb;
    double A      = p->A;
    double dA     = p->dA;
    double rc     = p->rc;
    double drc_dr = p->drc_dr;
    double ak2    = p->ak2;
    double S      = p->S;
    double Ss     = p->Ss;
    double nu     = p->nu;
    double chi1   = p->chi1;
    double chi2   = p->chi2;
    double X1     = p->X1;
    double X2     = p->X2;
    double c3     = p->c3;
    
    vector<double> ggm0 = s_GS(rorb, rc, drc_dr, ak2, 0., x, nu, chi1, chi2, X1, X2, c3);
    double dGS_dr  = ggm0[6];
    double dGSs_dr = ggm0[7];
    
    double x2 = x*x;
    
    double uc  = 1./rc;
    double uc2 = uc*uc;
    double uc3 = uc2*uc;
    
    // Orbital circular effective Hamiltonian
    double Horbeff0 = sqrt(A*(1. + x2*uc2));
    
    double dHeff_dr = x*(dGS_dr*S + dGSs_dr*Ss) + 1./(2.*Horbeff0)*( dA*(1. + x2*uc2) - 2.*A*uc3*drc_dr*x2);
    
    return dHeff_dr;
}

double s_bisec(double pph, double rorb, double A, double dA, double rc, double drc_dr, double ak2, double S, double Ss, void *params){
      /* This function solves the equation dHeff0=/dr=0 for the angular momentum.
         That is it computed the angular momentum of an EOB circular orbit.
         Spin is included. */

    double nu   = (*(TEOBResumParams *)params).nu;
    double chi1 = (*(TEOBResumParams *)params).chi1;
    double chi2 = (*(TEOBResumParams *)params).chi2;
    double X1   = (*(TEOBResumParams *)params).X1;
    double X2   = (*(TEOBResumParams *)params).X2;
    double c3   = (*(TEOBResumParams *)params).cN3LO;
    
    int status;
    int iter = 0, max_iter = 200;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    
    double r;
    double x_lo = 0.5*pph, x_hi = 1.5*pph;
    gsl_function F;
    struct energy_params DHeff_params = {rorb,A,dA,rc,drc_dr,ak2,S,Ss,nu,chi1,chi2,X1,X2,c3};
    F.function = &DHeff0;
    
    F.params = &DHeff_params;
    T = gsl_root_fsolver_bisection;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        r      = gsl_root_fsolver_root (s);
        x_lo   = gsl_root_fsolver_x_lower (s);
        x_hi   = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi,
                                         0, 0.000000000000001);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (s);
    
    return r;
}

vector<double> get_Omg_orb(vector<double> r,vector<double> pph,vector<double> pr_star,vector<double> A,vector<double> B,void *params)

/* This function computed the orbital frequency of the relative dynamics, 
   that is the derivative of the effective Hamiltonian with respect to
   the angular momentum. This is done here with the full spin dependence */ 

{
    
    double nu   = (*(TEOBResumParams *)params).nu;
    double aK2  = (*(TEOBResumParams *)params).aK2;
    double S1   = (*(TEOBResumParams *)params).S1;
    double S2   = (*(TEOBResumParams *)params).S2;
    double a1   = (*(TEOBResumParams *)params).a1;
    double a2   = (*(TEOBResumParams *)params).a2;
    double X1   = (*(TEOBResumParams *)params).X1;
    double X2   = (*(TEOBResumParams *)params).X2;
    double chi1 = (*(TEOBResumParams *)params).chi1;
    double chi2 = (*(TEOBResumParams *)params).chi2;
    double c3   = (*(TEOBResumParams *)params).cN3LO;
    
    double S     = S1 + S2;
    double Sstar = X2*a1 + X1*a2;
    double z3    = 2.*nu*(4.-3.*nu);
    
    long int r_length = r.size();
    vector<double> Omg_orb(r_length);
    
    for (int i=r_length; i--;)
    {
        
        vector<double> rc_vec = s_get_rc(r[i],params);
        double rc             = rc_vec[0];
        double drc_dr         = rc_vec[1];
        double uc             = 1./rc;
        double uc2            = uc*uc;
        
        vector<double> ggm    = s_GS(r[i],rc,drc_dr,aK2,pr_star[i],pph[i],nu,chi1,chi2,X1,X2,c3);
        double GS             = ggm[2];
        double GSs            = ggm[3];
        
        double prstar2        = pr_star[i]*pr_star[i];
        double prstar4        = prstar2*prstar2;
        double pphi2          = pph[i]*pph[i];
        double Horbeff        = sqrt(prstar2+A[i]*(1. + pphi2*uc2 +  z3*prstar4*uc2));
        
        double Heff           = Horbeff + (GS*S + GSs*Sstar)*pph[i];
        double H              = sqrt( 1. + 2.*nu*(Heff - 1.) );
        double one_H          = 1./H;
        
        Omg_orb[i]            = one_H*pph[i]*A[i]*uc2/Horbeff;
    }
    
    return Omg_orb;
    
}

vector<double> A_NumDenom(const double r, const vector<double> a, const double nu)
/* This function computed the A function as Pade(1,5). The numerator and denominator
   are computed separately. This is a duplication of another function and the derivatives
   are not included. It will have to be eventually removed, as this comes directly from
   the original version of the public Matlab code */

{
    
    /** Shorthands */
    const double u  = 1./r;
    const double u2 = u*u;
    const double u3 = u2*u;
    const double u4 = u3*u;
    const double u5 = u4*u;
    
    /**  Point-mass PN coefs */
    const double a3  = a[0];
    const double a4  = a[1];
    const double a5  = a[4];
    const double a6  = a[7];
    const double a5l = a[3];
    const double a6l = a[6];
    
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
    vector<double> frac(4);
    frac[0] = 1. + (C7*u)/C11;
    frac[1] = (-32. + nu*(12.*a3 + 4.*a4 + a5 - C2*u - 2.*C2*u2 - (C9)*u3 - (C8)*u4 - C5*u5))/C11;
    frac[2] = (-((4.*a5l + a6l)*nu*C4) + C4*C7 + a5l*nu*(-64 + nu*(12.*a4 + 4.*a5 + a6 + a3*(32 - a3*nu))))/(C4*C4);
    frac[3] = (nu*(-((2.*a5l + a6l)*C4) + a5l*nu*C2 + C4*C12 - 2.*(2.*a5l + a6l)*C4*u + 2.*a5l*nu*C2*u + 4.*C4*C12*u + a5l*nu*(C9)*u2 - 3.*(C9)*C4*u2 - C4*(4.*a6l + a5l*(8 + a3*nu))*u2 + a5l*nu*(C8)*u3 - 4.*C4*(C8)*u3 + C4*(a6l*C10 - a5l*(16 + a4*nu))*u3 + a5l*nu*C5*u4 - 5.*C4*C5*u4 + C4*(-2.*(4.*a3 + 2.*a4 + a5)*a5l*nu + a6l*(-16 + 4.*a3*nu + a4*nu))*u4))/(C11*C11);
    
    return frac;
}
