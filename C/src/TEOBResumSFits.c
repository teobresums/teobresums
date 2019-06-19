/**
 * This file is part of TEOBResumS
 *
 * Copyright (C) 2017-2018 See AUTHORS file
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
 * along with this program. If not, see http://www.gnu.org/licenses/.       
 *
 */

#include "TEOBResumS.h"

/** Fit of a6c, TEOBResumS paper Nagar et al. (2018) */
double eob_a6c_fit(double nu)
{
  return 3097.3*nu*nu - 1330.6*nu + 81.38;
}

/** Fit of a6c, higher multipoles: REF TO BE UPDATED */
double eob_a6c_fit_HM(double nu)
{
  double nu2 = nu*nu;
  double nu3 = nu2*nu;
  
  const double n0 =   5.9951;
  const double n1 = -34.4844;
  const double n2 = -79.2997;
  const double n3 = 713.4451;
  const double d1 =  -3.167;
  
  return n0*(1 + n1*nu + n2*nu2 + n3*nu3)/(1 + d1*nu);
}

/** Fit of c3, TEOBResumS paper Nagar et al. (2018) 
    Note: c3 = 0 with tides*/
double eob_c3_fit_global(double nu, double a1, double a2)
{  
  const double nu2 = nu*nu;
  const double nu3 = nu2*nu;
  const double X12 = sqrt(1.-4.*nu);
  const double a12 = a1+a2;
    
  /* Equal-mass, equal-spin coefficients */
  const double c0 =  43.371638;
  const double n1 =  -1.174839;
  const double n2 =   0.354064;
  const double d1 =  -0.151961;
	
  const double c3_eq = c0*(1. + n1*a12 + n2*a12*a12)/(1.+d1*a12);
  
  /* Coefficients 10/05/2018 */
  const double cnu    =  929.579;
  const double cnu2   = -9178.87;
  const double cnu3   =  23632.3;
  const double ca1_a2 = -104.891;
  
  const double c3_uneq = cnu*a12*nu*X12 + cnu2*a12*nu2*X12 + cnu3*a12*nu3*X12 + ca1_a2*(a1-a2)*nu2;
  
  return c3_eq + c3_uneq;
}

/** Fit of c3, HM paper: REF TO BE UPDATED 
    Note: c3 = 0 with tides*/
double eob_c3_fit_HM(double nu, double a1, double a2)
{  
  const double nu2 = nu*nu;
  const double X12 = sqrt(1.-4.*nu);
  const double a0  = a1+a2;
  const double a02 = a0*a0;
  const double a03 = a02*a0;
  const double a04 = a03*a0;
  
  /* Equal-mass, equal-spin coefficients */
  const double p0 =  46.323097;
  const double n1 =  -1.418623;
  const double n2 =   0.457118;
  const double n3 =   0.000435;
  const double n4 =   0.097791;
  const double d1 =  -0.519003;
	  
  /* Other coefficients */
  const double p1 =   46.2914;
  const double p2 = -103.0997;
  
  const double c3 = p0*(1 + n1*a0 + n2*a02 + n3*a03 + n4*a04)/(1 + d1*a0)
    + p1*nu*X12*a0 + p2*nu2*(a1 - a2);
  
  return c3;
}

/** Function providing a fit of Deltat_NQC vs chi, via a simple rational function. */
double eob_nqc_dtfit(const double chi, const double chi0)
{
  const double n1 = -16.06288206;
  const double d1 = -4.04266459;

  double x     = chi-chi0;

  double dtnqc = (1.+n1*x)/(1.+d1*x);

  return dtnqc;
}

/** Fits for NR point used to determine NQC corrections */
void eob_nqc_point(Dynamics *dyn, double *A_tmp, double *dA_tmp, double *omg_tmp, double *domg_tmp)
{

  const double nu   = dyn->nu;
  const double X1   = dyn->X1;
  const double X2   = dyn->X2;
  const double chi1 = dyn->chi1;
  const double chi2 = dyn->chi2;
  const double aK   = dyn->a1 + dyn->a2;

  const double nu2  = SQ(nu);
  const double nu3  = nu2*nu;
  const double X12  = X1 - X2;  
  const double aK2  = SQ(aK);
  const double aK3  = aK2*aK;
  const double aK4  = aK2*aK2;
  const double a12  = X1*chi1 - X2*chi2;
  const double aeff     = aK + 1./3.*a12*X12;
  const double aeff_omg = aK + a12*X12;

  double pA[5], pdA[5], pomg[5], pdomg[5];
  double c_p1,     c_p2,     c_p3,   c_p4;
  double c_pdA1,   c_pdA2,   c_pdA3, c_pdA4;
  double c_pdomg1, c_pdomg2;
  double n0, d1;

    double a0_omg_tmp, a1_omg_tmp, a2_omg_tmp, b0_omg_tmp, b1_omg_tmp, b2_omg_tmp, a0_domg_tmp, a1_domg_tmp, a2_domg_tmp, b0_domg_tmp, b1_domg_tmp, b2_domg_tmp, a0_A_tmp, a1_A_tmp , a2_A_tmp, b0_A_tmp, b1_A_tmp, b2_A_tmp, a0_dA_tmp, a1_dA_tmp, a2_dA_tmp, b0_dA_tmp, b1_dA_tmp, b2_dA_tmp, omg_tmp_nu, omg_tmp_equal, domg_tmp_nu, domg_tmp_equal,  A_tmp_scale_nu, A_tmp_scale_equal, dA_tmp_scale_nu, dA_tmp_scale_equal ;

    double p1[2], p2[2], p3[2], p4[2]; 
    double pn0[2], pd1[2], ppdomg1[2], ppdomg2[2], pdA1[2],pdA2[2],pdA3[2],pdA4[2];

      
  if (DEQUAL(nu,0.25,1e-9)) {

    pA[0]    =  0.00178195;
    pA[1]    =  0.00435589;
    pA[2]    =  0.00344489;
    pA[3]    = -0.00076165;
    pA[4]    =  0.31973334;
    *A_tmp    =  pA[0]*aK4    + pA[1]*aK3   + pA[2]*aK2    + pA[3]*aK     + pA[4];
    
    pdA[0]   =  0.00000927;
    pdA[1]   = -0.00024550;
    pdA[2]   =  0.00012469;
    pdA[3]   =  0.00123845;
    pdA[4]   = -0.00195014;
    *dA_tmp   =  pdA[0]*aK4   + pdA[1]*aK3   + pdA[2]*aK2   + pdA[3]*aK   + pdA[4];
    
    pomg[0]  =  0.00603482;
    pomg[1]  =  0.01604555;
    pomg[2]  =  0.02290799;
    pomg[3]  =  0.07084587;
    pomg[4]  =  0.38321834;
    *omg_tmp  =  pomg[0]*aK4  + pomg[1]*aK3  + pomg[2]*aK2  + pomg[3]*aK  + pomg[4];
    
    pdomg[0] =  0.00024066;
    pdomg[1] =  0.00038123;
    pdomg[2] = -0.00049714;
    pdomg[3] =  0.00041219;
    pdomg[4] =  0.01190548;
    *domg_tmp =  pdomg[0]*aK4 + pdomg[1]*aK3 + pdomg[2]*aK2 + pdomg[3]*aK + pdomg[4];
  
  }  else if( nu > 0.16) {
    
    p1[0]      =  0.04680896;
    p1[1]      = -0.00632114;
    p2[0]      =  0.06586192;
    p2[1]      = -0.01180039;
    p3[0]      = -0.11617413;
    p3[1]      =  0.02704959;
    p4[0]      =  0.15597465;
    p4[1]      =  0.28034978;
    c_p1       =  p1[0]*nu + p1[1];
    c_p2       =  p2[0]*nu + p2[1];
    c_p3       =  p3[0]*nu + p3[1];
    c_p4       =  p4[0]*nu + p4[1];
    *A_tmp      =  c_p1*aK3 + c_p2*aK2 + c_p3*aK + c_p4;
    
    pdA1[0]    = -0.00130824;
    pdA1[1]    =  0.00006202;
    pdA2[0]    =  0.00199855;
    pdA2[1]    = -0.00027474;
    pdA3[0]    =  0.00218838;
    pdA3[1]    =  0.00071540;
    pdA4[0]    = -0.00362779;
    pdA4[1]    = -0.00105397;
    c_pdA1     =  pdA1[0]*nu + pdA1[1];
    c_pdA2     =  pdA2[0]*nu + pdA2[1];
    c_pdA3     =  pdA3[0]*nu + pdA3[1];
    c_pdA4     =  pdA4[0]*nu + pdA4[1];
    *dA_tmp     =  c_pdA1*aK3   + c_pdA2*aK2 + c_pdA3*aK+ c_pdA4;
    
    pn0[0]     =  0.46908067;
    pn0[1]     =  0.27022141;
    pd1[0]     =  0.64131115;
    pd1[1]     = -0.37878384;
    n0         =  pn0[0]*nu + pn0[1];
    d1         =  pd1[0]*nu + pd1[1];
    *omg_tmp    =  n0/(1 + d1*aK);
    
    ppdomg1[0] =  0.00061175;
    ppdomg1[1] =  0.00074001;
    ppdomg2[0] =  0.02504442;
    ppdomg2[1] =  0.00548217;
    c_pdomg1   =  ppdomg1[0]*nu + ppdomg1[1];
    c_pdomg2   =  ppdomg2[0]*nu + ppdomg2[1];
    *domg_tmp   =  c_pdomg1*aK   + c_pdomg2;
  
  }  else {

    /* Fit by G.Riemanschneider incorporating the test-particle NQC point
       obtained from the most-recent Teukolsky waveforms done by
       M. Colleoni using the 6PN-accurare iResum-radiation reaction.
       These points assure a smooth connection between merger and
       ringdown also outside the "calibration" domain, notably for
       large-mass ratios (though q<=20) and large (negative) spins
       Updated, 28/09/2017 */

    a0_A_tmp 	= -0.2750516062;
    b0_A_tmp 	= -0.4693776065;
    a1_A_tmp 	=  0.143066;
    a2_A_tmp 	= -0.0425947;
    b1_A_tmp 	=  0.176955;
    b2_A_tmp 	= -0.111902;
    
    A_tmp_scale_nu    = -0.9862040409*nu3 +0.8167558040*nu2 -0.0427442282*nu+0.2948879452;
    A_tmp_scale_equal = ((a2_A_tmp*X12*X12 + a1_A_tmp*X12 +a0_A_tmp)*aeff+1)/((b2_A_tmp*X12*X12 + b1_A_tmp*X12 +b0_A_tmp)*aeff+1);
    *A_tmp             = A_tmp_scale_nu*A_tmp_scale_equal*(1-0.5*(*omg_tmp)*aeff);
    
    a0_dA_tmp 	= +0.0037461628;
    b0_dA_tmp 	= +0.0636082543;
    a1_dA_tmp 	=  0.00129393;
    a2_dA_tmp 	= -0.00239069;
    b1_dA_tmp 	= -0.0534209;
    b2_dA_tmp 	= -0.186101;
    
    dA_tmp_scale_nu    = ( -0.0847947167*nu -0.0042142765)/( +16.1559461812*nu+1);
    dA_tmp_scale_equal = ((a2_dA_tmp*X12*X12 + a1_dA_tmp*X12+ a0_dA_tmp)*aeff)/((b2_dA_tmp*X12*X12 + b1_dA_tmp*X12 + b0_dA_tmp)*aeff+1);
    *dA_tmp             = (dA_tmp_scale_nu +dA_tmp_scale_equal)*(*omg_tmp);
    
    a0_omg_tmp    = -0.1460961247;
    a1_omg_tmp    =  0.0998056;
    a2_omg_tmp    = -0.118098;
    b0_omg_tmp    = -0.3430184009;
    b1_omg_tmp    =  0.0921551;
    b2_omg_tmp    = -0.0740285;
    omg_tmp_nu    = +0.5427169903*nu2 +0.2512395608*nu +0.2863992248;
    omg_tmp_equal =((a2_omg_tmp*X12*X12 + a1_omg_tmp*X12 + a0_omg_tmp)*aeff_omg+1)/((b2_omg_tmp*X12*X12 +b1_omg_tmp*X12 + b0_omg_tmp)*aeff_omg+1);
    *omg_tmp       = omg_tmp_nu*omg_tmp_equal;

    a0_domg_tmp    = +0.0604556289;
    b0_domg_tmp    = -0.0299583285;
    a1_domg_tmp    = 0.0711715;
    a2_domg_tmp    = -0.0500886;
    b1_domg_tmp    = 0.0461239;
    b2_domg_tmp    = -0.0153068;
    
    domg_tmp_nu    = ( +0.0045213831*nu +0.0064934920)/( -1.4466409969*nu+1);
    domg_tmp_equal = (a2_domg_tmp*X12*X12 +a1_domg_tmp*X12 +b0_domg_tmp)*aeff_omg*aeff_omg +(b2_domg_tmp*X12*X12 +b1_domg_tmp*X12+a0_domg_tmp)*aeff_omg+1;
    *domg_tmp       = domg_tmp_nu*domg_tmp_equal;

  }

}

/** Fits for NR point used to determine NQC corrections 
    New fits for higer modes paper
    Ref TO BE UPDATED                                   */
void eob_nqc_point_HM(Dynamics *dyn, double *A_tmp, double *dA_tmp, double *omg_tmp, double *domg_tmp)
{

  const double nu   = dyn->nu;
  const double X1   = dyn->X1;
  const double X2   = dyn->X2;
  const double chi1 = dyn->chi1;
  const double chi2 = dyn->chi2;
  const double aK   = dyn->a1 + dyn->a2;

  const double nu2  = SQ(nu);
  const double nu3  = nu2*nu;
  const double X12  = X1 - X2;  
  const double aK2  = SQ(aK);
  const double aK3  = aK2*aK;
  const double aK4  = aK2*aK2;
  const double a12  = X1*chi1 - X2*chi2;
  const double Shat = 0.5*(aK + a12*X12);

  double pA[5], pdA[5], pomg[5], pdomg[5];
  double c_p1,     c_p2,     c_p3,   c_p4;
  double c_pdA1,   c_pdA2,   c_pdA3, c_pdA4;
  double c_pdomg1, c_pdomg2;
  double n0, d1;

  double a0_omg_tmp, a1_omg_tmp, a2_omg_tmp, b0_omg_tmp, b1_omg_tmp, b2_omg_tmp, c11_omg_tmp, c12_omg_tmp, c21_omg_tmp, c22_omg_tmp, c31_omg_tmp, c32_omg_tmp;
  double a0_domg_tmp, a1_domg_tmp, a2_domg_tmp, b1_domg_tmp, b2_domg_tmp, c1_domg_tmp, c2_domg_tmp;
  double scale, a0_A_tmp, a1_A_tmp , a2_A_tmp, b0_A_tmp, b1_A_tmp, b2_A_tmp, c11_A_tmp, c12_A_tmp, c21_A_tmp, c22_A_tmp, c31_A_tmp, c32_A_tmp;
  double a0_dA_tmp, a1_dA_tmp, a2_dA_tmp, b1_dA_tmp, b2_dA_tmp, c1_dA_tmp, c2_dA_tmp;
  double omg_tmp_nu, omg_tmp_equal, omg_tmp_den, domg_tmp_nu, domg_tmp_equal, A_tmp_scale_nu, A_tmp_scale_equal, A_tmp_scale_den, dA_tmp_scale_nu, dA_tmp_scale_equal;

  double p1[2], p2[2], p3[2], p4[2]; 
  double pn0[2], pd1[2], ppdomg1[2], ppdomg2[2], pdA1[2],pdA2[2],pdA3[2],pdA4[2];

      
  if (DEQUAL(nu,0.25,1e-9)) {

    pA[0]    =  0.00178195;
    pA[1]    =  0.00435589;
    pA[2]    =  0.00344489;
    pA[3]    = -0.00076165;
    pA[4]    =  0.31973334;
    *A_tmp   =  pA[0]*aK4    + pA[1]*aK3   + pA[2]*aK2    + pA[3]*aK     + pA[4];
      
    pdA[0]   =  0.00000927;
    pdA[1]   = -0.00024550;
    pdA[2]   =  0.00012469;
    pdA[3]   =  0.00123845;
    pdA[4]   = -0.00195014;
    *dA_tmp  =  pdA[0]*aK4   + pdA[1]*aK3   + pdA[2]*aK2   + pdA[3]*aK   + pdA[4];
    
    pomg[0]  =  0.00603482;
    pomg[1]  =  0.01604555;
    pomg[2]  =  0.02290799;
    pomg[3]  =  0.07084587;
    pomg[4]  =  0.38321834;
    *omg_tmp =  pomg[0]*aK4  + pomg[1]*aK3  + pomg[2]*aK2  + pomg[3]*aK  + pomg[4];
    
    pdomg[0]  =  0.00024066;
    pdomg[1]  =  0.00038123;
    pdomg[2]  = -0.00049714;
    pdomg[3]  =  0.00041219;
    pdomg[4]  =  0.01190548;
    *domg_tmp =  pdomg[0]*aK4 + pdomg[1]*aK3 + pdomg[2]*aK2 + pdomg[3]*aK + pdomg[4];
    
  }  else if( nu > 0.16) {
      
    p1[0]      =  0.04680896;
    p1[1]      = -0.00632114;
    p2[0]      =  0.06586192;
    p2[1]      = -0.01180039;
    p3[0]      = -0.11617413;
    p3[1]      =  0.02704959;
    p4[0]      =  0.15597465;
    p4[1]      =  0.28034978;
    c_p1       =  p1[0]*nu + p1[1];
    c_p2       =  p2[0]*nu + p2[1];
    c_p3       =  p3[0]*nu + p3[1];
    c_p4       =  p4[0]*nu + p4[1];
    *A_tmp     =  c_p1*aK3 + c_p2*aK2 + c_p3*aK + c_p4;
      
    pdA1[0]    = -0.00130824;
    pdA1[1]    =  0.00006202;
    pdA2[0]    =  0.00199855;
    pdA2[1]    = -0.00027474;
    pdA3[0]    =  0.00218838;
    pdA3[1]    =  0.00071540;
    pdA4[0]    = -0.00362779;
    pdA4[1]    = -0.00105397;
    c_pdA1     =  pdA1[0]*nu + pdA1[1];
    c_pdA2     =  pdA2[0]*nu + pdA2[1];
    c_pdA3     =  pdA3[0]*nu + pdA3[1];
    c_pdA4     =  pdA4[0]*nu + pdA4[1];
    *dA_tmp    =  c_pdA1*aK3   + c_pdA2*aK2 + c_pdA3*aK+ c_pdA4;
    
    pn0[0]     =  0.46908067;
    pn0[1]     =  0.27022141;
    pd1[0]     =  0.64131115;
    pd1[1]     = -0.37878384;
    n0         =  pn0[0]*nu + pn0[1];
    d1         =  pd1[0]*nu + pd1[1];
    *omg_tmp   =  n0/(1 + d1*aK);
    
    ppdomg1[0] =  0.00061175;
    ppdomg1[1] =  0.00074001;
    ppdomg2[0] =  0.02504442;
    ppdomg2[1] =  0.00548217;
    c_pdomg1   =  ppdomg1[0]*nu + ppdomg1[1];
    c_pdomg2   =  ppdomg2[0]*nu + ppdomg2[1];
    *domg_tmp  =  c_pdomg1*aK   + c_pdomg2;
  
    }  else {

    a0_omg_tmp    = 0.285588;
    a1_omg_tmp    = 0.91704;
    a2_omg_tmp    = 1.7912;
    b0_omg_tmp    = -0.46550;
    b1_omg_tmp    = -0.078787;
    b2_omg_tmp    = -0.852284;
    c11_omg_tmp   = -0.338008;
    c12_omg_tmp   = 1.077812;
    c21_omg_tmp   = 0.0555533;
    c22_omg_tmp   = -0.312861;
    c31_omg_tmp   = 0.289185;
    c32_omg_tmp   = -0.195838;

    omg_tmp_nu    = a0_omg_tmp*(1. + a1_omg_tmp*nu + a2_omg_tmp*nu2);
    omg_tmp_equal = 1. + (b0_omg_tmp + c11_omg_tmp*X12)/(1. + c12_omg_tmp*X12)*Shat + (b1_omg_tmp + c21_omg_tmp*X12)/(1. + c22_omg_tmp*X12)*Shat*Shat;
    omg_tmp_den   = 1. + (b2_omg_tmp + c31_omg_tmp*X12)/(1. + c32_omg_tmp*X12)*Shat;
    *omg_tmp      = omg_tmp_nu*omg_tmp_equal/omg_tmp_den;
    
    a0_domg_tmp    = 0.00628027;
    a1_domg_tmp    = 2.4351;
    a2_domg_tmp    = 4.4928;
    b1_domg_tmp    = 0.001425242;
    b2_domg_tmp    = -0.00096073;
    c1_domg_tmp    = -0.000063766;
    c2_domg_tmp    = 0.000513197;

    domg_tmp_nu    = a0_domg_tmp*(1. + a1_domg_tmp*nu + a2_domg_tmp*nu2);
    domg_tmp_equal = (b1_domg_tmp + c1_domg_tmp*X12)*Shat
      + (b2_domg_tmp + c2_domg_tmp*X12)*Shat*Shat;
    *domg_tmp      = domg_tmp_nu + domg_tmp_equal;

    scale       = 1 - (*omg_tmp)*Shat;
    a0_A_tmp 	= 0.294773;
    a1_A_tmp 	= -0.052697;
    a2_A_tmp 	= 1.6088;
    b0_A_tmp 	= -0.705226;
    b1_A_tmp 	= -0.0953944;
    b2_A_tmp 	= -1.087280; 
    c11_A_tmp   = 0.009335;
    c12_A_tmp   = 0.582869;
    c31_A_tmp   = -0.140747;
    c32_A_tmp   = 0.505807;
    
    A_tmp_scale_nu    = a0_A_tmp*(1. + a1_A_tmp*nu + a2_A_tmp*nu2)*scale;
    A_tmp_scale_equal = 1. + (b0_A_tmp + c11_A_tmp*X12)/(1. + c12_A_tmp*X12)*Shat + b1_A_tmp*Shat*Shat;
    A_tmp_scale_den   = 1. + (b2_A_tmp + c31_A_tmp*X12)/(1. + c32_A_tmp*X12)*Shat;
    *A_tmp            = A_tmp_scale_nu*A_tmp_scale_equal/A_tmp_scale_den;
    
    scale       = 0.5*(*omg_tmp)/sqrt(6);
    a0_dA_tmp 	= -0.0011936600;
    a1_dA_tmp 	=  2.86637;
    a2_dA_tmp 	= -1.3667;
    b1_dA_tmp 	=  0.02679530;
    b2_dA_tmp 	= -0.0064409;
    c1_dA_tmp   = -0.015395218;
    c2_dA_tmp   = 0.008732589;

    dA_tmp_scale_nu    = a0_dA_tmp*(1. + a1_dA_tmp*nu + a2_dA_tmp*nu2)*scale;
    dA_tmp_scale_equal = (b1_dA_tmp + c1_dA_tmp*X12)*Shat + (b2_dA_tmp + c2_dA_tmp*X12)*Shat*Shat;
    *dA_tmp            = dA_tmp_scale_nu*dA_tmp_scale_equal;

    }

}

/** Time-shift for NQC */
double eob_nqc_timeshift(double nu, double chi1)
{

  double DeltaT_nqc = 1.0;  
  
  /* Additional time-shift only needed ONLY for large, negative, spins.
     This change from 1.0 to 4.0 eliminates unphysical features in the
     frequency related to the imperfect behavior of the NQC functions */

  /* Old Delta_T NQC
  if ((chi1 <-0.85) && (nu <= 14./225.)) {
    DeltaT_nqc = 4.0;
  } else {
    DeltaT_nqc = 1.0; // standard choice inspired by test-particle results
  }  
  */

  /* New Delta_T NQC determined in TEOBResumS paper (arXiv:1806.01772) */
  if (((chi1 < -0.9) && (nu < 8./81.)) || ((chi1 < -0.8) && (nu < 11./144.))) {
    DeltaT_nqc = 4.0;
  } else {
    DeltaT_nqc = 1.0; // standard choice inspired by test-particle results
  }
    
  return DeltaT_nqc;  
}

/** Set NQC coefficients */
void eob_nqc_setcoefs(NQCdata *nqc)
{

  nqc->flx->add = 1;
  nqc->hlm->add = 1;

  if (STREQUAL(par_get_s("nqc_coefs_flx"),"none"))
    nqc->flx->add = 0;
  if (STREQUAL(par_get_s("nqc_coefs_hlm"),"none"))
    nqc->hlm->add = 0;
      
  /* Init NQC coefs to zero */
  for (int k = 0; k < KMAX; k++) {
    for (int j = 0; j < 6; j++) {
      nqc->flx->n[k][j] = 0.;
      nqc->hlm->n[k][j] = 0.;
    }
    nqc->flx->a1[k] = 0.;
    nqc->flx->a2[k] = 0.;
    nqc->flx->a3[k] = 0.;
    nqc->flx->b1[k] = 0.;
    nqc->flx->b2[k] = 0.;
    nqc->flx->b3[k] = 0.;
    nqc->flx->activemode[k] = 0; 
    nqc->hlm->a1[k] = 0.;
    nqc->hlm->a2[k] = 0.;
    nqc->hlm->a3[k] = 0.;
    nqc->hlm->b1[k] = 0.;
    nqc->hlm->b2[k] = 0.;
    nqc->hlm->b3[k] = 0.;
    nqc->hlm->activemode[k] = 0; 
  }
  nqc->flx->maxk = -1;
  nqc->hlm->maxk = -1; 
  
  if (nqc->flx->add + nqc->hlm->add == 0) 
    return;

  if (STREQUAL(par_get_s("nqc_coefs_flx"),"nrfit_nospin201602")) 
    eob_nqc_setcoefs_nospin201602(nqc->flx);
  if (STREQUAL(par_get_s("nqc_coefs_flx"),"fromfile")) 
    eob_nqc_setcoefs_fromfile(nqc->flx, par_get_s("nqc_coefs_flx_file"));
  // TODO: ADD HERE YOUR LATEST FITS
  //else if (STREQUAL(par_get_s("nqc_coefs_flx"),"nrfit_spin_202001")) 
  
  if (STREQUAL(par_get_s("nqc_coefs_hlm"),"nrfit_nospin201602")) 
    eob_nqc_setcoefs_nospin201602(nqc->hlm);
   if (STREQUAL(par_get_s("nqc_coefs_hlm"),"fromfile")) 
    eob_nqc_setcoefs_fromfile(nqc->hlm, par_get_s("nqc_coefs_hlm_file"));
   // TODO: ADD HERE YOUR LATEST FITS
   //else if (STREQUAL(par_get_s("nqc_coefs_hlm"),"nrfit_spin_202001")) 
  
}

/** Set NQC coefficients 
    NR fits for nonspinning case 2016/02/09 
    Hardcoded in eob_wav_hlmNQC_nospin201602() */
void eob_nqc_setcoefs_nospin201602(NQCcoefs *nqc)
{

  const double nu = par_get_d("nu");  
  const double xnu  = 1-4*nu;
  const double xnu2 = SQ(xnu);

  const int k21 = 0;
  const int k22 = 1;
  const int k33 = 4;

  nqc->activemode[k21]=1;
  nqc->activemode[k22]=1;
  nqc->activemode[k33]=1;
  
  /* (2,1) */
  nqc->a1[k21] =  0.0162387198*(7.32653082*xnu2 + 1.19616248*xnu + 0.73496656);
  nqc->a2[k21] = -1.80492460*xnu2 + 1.78172686*xnu + 0.30865284;
  nqc->a3[k21] =  0.0;
    
  nqc->b1[k21] = -0.0647955017*(3.59934444*xnu2 - 4.08628784*xnu + 1.37890907);
  nqc->b2[k21] =  1.3410693180*(0.38491989*xnu2 + 0.10969453*xnu + 0.97513971);
  nqc->b3[k21] =  0.0;
  
  /* (2,2) */
  nqc->a1[k22] = -0.0805236959*( 1 - 2.00332326*xnu2)/( 1 + 3.08595088*xnu2);
  nqc->a2[k22] =  1.5299534255*( 1 + 1.16438929*xnu2)/( 1 + 1.92033923*xnu2);
  nqc->a3[k22] =  0.0;
    
  nqc->b1[k22] = 0.146768094955*( 0.07417121*xnu + 1.01691256);
  nqc->b2[k22] = 0.896911234248*(-0.61072011*xnu + 0.94295129);
  nqc->b3[k22] = 0.0;
    
  /* (3,3) */
  nqc->a1[k33] = -0.0377680000*(1 - 14.61548907*xnu2)/( 1 + 2.44559263*xnu2);
  nqc->a2[k33] =  1.9898000000*(1 + 2.09750346 *xnu2)/( 1 + 2.57489466*xnu2);
  nqc->a3[k33] =  0.0;
    
  nqc->b1[k33] = 0.1418400000*(1.07430512 - 1.23906804*xnu + 4.44910652*xnu2);
  nqc->b2[k33] = 0.6191300000*(0.80672432 + 4.07432829*xnu - 7.47270977*xnu2);
  nqc->b3[k33] = 0.0;

  nqc->add = 1;
  nqc->maxk = k33;

}

/** Set NQC coefficients from file */
//TODO: how portable is 'getline'
void eob_nqc_setcoefs_fromfile(NQCcoefs *nqc, const char *fname)
{
  FILE *fp;
  char *line = NULL;
  size_t len = 0;
  size_t read;
  int nl = 0;
  
  int k, maxk=-1;
  double a1k,a2k,a3k, b1k,b2k,b3k;
  
  if ((fp = fopen(fname, "r")) == NULL)
    errorexits("error opening file",fname);
  
  while ((read = getline(&line, &len, fp)) != -1) {
    /* printf("Retrieved line of length %zu :\n", read);  */
    /* printf("%s", line); */
    if (line[0]=='#') continue; /* skip comment */
    nl++;
    if (nl>KMAX) break;
    sscanf(line, "%d %*d %*d %lf %lf %lf %lf %lf %lf", 
	   &k, &a1k,&a2k,&a3k, &b1k,&b2k,&b3k);
    /* printf("%d: %d %e %e %e %e %e %e\n",nl, k,a1k,a2k,a3k, b1k,b2k,b3k); */
    if (k>=0 && k<KMAX) {
      nqc->a1[k] = a1k;
      nqc->a2[k] = a2k;
      nqc->a3[k] = a3k;
      nqc->b1[k] = b1k;
      nqc->b2[k] = b2k;
      nqc->b3[k] = b3k;
      if (k>maxk) maxk = k;
      nqc->activemode[k]=1;
    }
  }
  
  if (line) free(line);
  fclose(fp);

  nqc->add = 1;
  nqc->maxk = maxk;
}

/** logQ-vs-log(lambda) fit of Table I of Yunes-Yagi
    here x = log(lambda) and the output is the log of the coefficient
    that describes the quadrupole deformation due to spin. */
double logQ(double x)
{
  const double ai = 0.194;
  const double bi = 0.0936;
  const double ci = 0.0474;
  const double di = -4.21e-3;
  const double ei = 1.23e-4;
  const double x2 = x*x;
  const double x3 = x*x2;
  const double x4 = x*x3;
  return ai + bi*x + ci*x2 + di*x3 + ei*x4;
}

/** Yagi 2013 fits for NS multipolar
    $\bar{\lambda}_\ell$ = 2 k_\ell/(C^{2\ell+1} (2\ell-1)!!)$
    Eq.(9,10),(61); Tab.I; Fig.8 http://arxiv.org/abs/1311.0872 */
double Yagi13_fit_barlamdel(double barlam2, int ell)
{
  if (barlam2<=0.) return 0.;
  double lnx = log(barlam2);
  double coeffs[5];
  if (ell == 3) {
    coeffs[0] = 2.52e-5;
    coeffs[1] = -1.31e-3;
    coeffs[2] = 2.51e-2;
    coeffs[3] = 1.18;
    coeffs[4] = -1.15;
  } else if (ell == 4) {
    coeffs[0] = 2.8e-5;
    coeffs[1] =-1.81e-3;
    coeffs[2] =3.95e-2;
    coeffs[3] =1.43;
    coeffs[4] =-2.45;
  }
  else 
    errorexit("Yagi fits are for ell=3,4.");
  double lny = coeffs[0]*lnx*lnx*lnx*lnx+coeffs[1]*lnx*lnx*lnx+coeffs[2]*lnx*lnx+coeffs[3]*lnx+coeffs[4];
  return exp(lny);
}

/** Yagi 2013 fits for NS multipolar
    $\bar{\sigma_2}( \bar{\lambda}_2 )$
    Eq.(9,11),(61); Tab.I; Fig.9 http://arxiv.org/abs/1311.0872 
    See also later erratum */
double Yagi13_fit_barsigmalambda(double barlam2)
{
  if (barlam2<=0.) return 0.;
  double lnx = log(barlam2);
  double coeffs[5];
  /*
  coeffs[4] = 0.126;
  coeffs[3] = 0.617;
  coeffs[2] = 2.81e-2;
  coeffs[1] = 3.59e-4;
  coeffs[0] = -3.61e-5;
  */
  coeffs[4] = -2.01;
  coeffs[3] =  0.462;
  coeffs[2] =  1.68e-2;
  coeffs[1] = -1.58e-4;
  coeffs[0] = -6.03e-6;
  double lny = coeffs[0]*lnx*lnx*lnx*lnx+coeffs[1]*lnx*lnx*lnx+coeffs[2]*lnx*lnx+coeffs[3]*lnx+coeffs[4];

  return -1.0*exp(lny);
}

/* Yagi et al. fits for C_Oct
   Eq. (90) and Table I of https://arxiv.org/abs/1403.6243 */
double Yagi14_fit_Coct(double C_Q)
{
  double A0  = -0.925;
  double B1  =  1.98;
  double nu1 =  0.273;

  double cubrootCoct = A0 + B1*pow(C_Q,nu1);

  return cubrootCoct*cubrootCoct*cubrootCoct;
}

/* Yagi et al. fits for C_Hex
   Eq. (90) and Table I of https://arxiv.org/abs/1403.6243 */
double Yagi14_fit_Chex(double C_Q)
{
  double A0  = -0.413;
  double B1  =  1.5;
  double nu1 =  0.466;

  double fourthrootChex = A0 + B1*pow(C_Q,nu1);

  return SQ(SQ(fourthrootChex));
}

double JFAPG_fit_Sigma_Irrotational(double barlam2)
{
  if (barlam2<=0.) return 0.;
  double lnx = log(barlam2);
  double coeffs[6];
 
  coeffs[5] = -2.03;
  coeffs[4] =  0.487;
  coeffs[3] =  9.69e-3;
  coeffs[2] =  1.03e-3;
  coeffs[1] = -9.37e-5;
  coeffs[0] =  2.24e-6;
  double lny = coeffs[0]*lnx*lnx*lnx*lnx*lnx+coeffs[1]*lnx*lnx*lnx*lnx+coeffs[2]*lnx*lnx*lnx+coeffs[3]*lnx*lnx+coeffs[4]*lnx+coeffs[5];

  return -1.0*exp(lny);
}

double JFAPG_fit_Sigma_Static(double barlam2)
{
  if (barlam2<=0.) return 0.;
  double lnx = log(barlam2);
  double coeffs[6];
 
  coeffs[5] = -2.66;
  coeffs[4] =  0.786;
  coeffs[3] =  -0.01;
  coeffs[2] =  1.28e-3;
  coeffs[1] = -6.37e-5;
  coeffs[0] =  1.18e-6;
  double lny = coeffs[0]*lnx*lnx*lnx*lnx*lnx+coeffs[1]*lnx*lnx*lnx*lnx+coeffs[2]*lnx*lnx*lnx+coeffs[3]*lnx*lnx+coeffs[4]*lnx+coeffs[5];

  return exp(lny);
}
 
/** Mass and angular momentum of the final black hole
  Healey, Lousto and Zochlower (HLZ),
  arXiv: 1406.7295, published as PRD 90, 104004 (2014)
  WARNING: the formula uses the convention that M2 > M1, so that
  chi2 should refer to the black hole with the largest
  mass. In the EOB code, this is given by chi1, since
  in EOB code we use the convention that M1 > M2
  
  Here it is q=M2/M1, with M2>M1
  
  Improved with (Eisco, Jisco) + iterative procedure 23/02/2016
  parameters (TABLE VI)
*/
void HealyBBHFitRemnant(double chi1,double chi2, double q, double *mass, double *spin)
{
       
  /** Final mass:                    Angular momentum: */
  
  double M0  =  0.951507;            double L0  =  0.686710;
  double K1  = -0.051379;            double L1  =  0.613247;
  double K2a = -0.004804;            double L2a = -0.145427;
  double K2b = -0.054522;            double L2b = -0.115689;
  double K2c = -0.000022;            double L2c = -0.005254;
  double K2d =  1.995246;            double L2d =  0.801838;
  double K3a =  0.007064;            double L3a = -0.073839;
  double K3b = -0.017599;            double L3b =  0.004759;
  double K3c = -0.119175;            double L3c = -0.078377;
  double K3d =  0.025000;            double L3d =  1.585809;
  double K4a = -0.068981;            double L4a = -0.003050;
  double K4b = -0.011383;            double L4b = -0.002968;
  double K4c = -0.002284;            double L4c =  0.004364;
  double K4d = -0.165658;            double L4d = -0.047204;
  double K4e =  0.019403;            double L4e = -0.053099;
  double K4f =  2.980990;            double L4f =  0.953458;
  double K4g =  0.020250;            double L4g = -0.067998;
  double K4h = -0.004091;            double L4h =  0.001629;
  double K4i =  0.078441;            double L4i = -0.066693;
      
  /** Parameters */
  double nu      = q/((1.+q)*(1.+q));
  
  /** Masses: convention here is that m2>m1 */
  double X2      = 0.5*(1.+sqrt(1.-4*nu));
  double X1      = 1.-X2;
  
  /** Spin variables */
  double s1      = X1*X1*chi1;
  double s2      = X2*X2*chi2;
  double S       = s1 + s2;
  double S2      = S*S;
  double S3      = S*S2;
  double S4      = S2*S2;
  double Delta   = X1/X2*s2 - X2/X1*s1 + s2 - s1;
  double Delta2  = Delta*Delta;
  double Delta3  = Delta*Delta2;
  double Delta4  = Delta2*Delta2;
  
  /** Mass ratio variables */
  double deltam  = -sqrt(1-4*nu); // X1 - X2
  double deltam2 =  deltam*deltam;
  double deltam3 =  deltam*deltam2;
  double deltam4 =  deltam*deltam3;
  double deltam6 =  deltam2*deltam4;
  
  /** Initialize the angular momentum */
  double a0 = s1 + s2;
  int a0_sign = 0.;
  
  if (a0==0) {
    a0_sign=0;
  } else if (a0>0) {
    a0_sign=1;
  } else { // if (a0<0) {
    a0_sign=-1;
  }
  
  /** Set-up an interative procedure to compute properly the "isco" quantities */
  double a2;
  double Z1;
  double Z2;
  double risco;
  double uisco;
  double Eisco;
  double Jisco;
  double abh;
  double Mbh=0.;
  
  int i;
  for(i=0; i<20; i++) {
    a2     = a0*a0;
    Z1     = 1 + cbrt(1-a2)*(cbrt(1+a0) + cbrt(1-a0));
    Z2     = sqrt(3*a2 + Z1*Z1);
    risco  = 3 + Z2 - a0_sign*sqrt((3-Z1)*(3+Z1+2.*Z2));
    uisco  = 1./risco;
    Eisco  = (1 - 2.*uisco + a0*sqrt(uisco*uisco*uisco))/sqrt(1-3*uisco + 2*a0*sqrt(uisco*uisco*uisco));
    Jisco  = 2./(sqrt(3.*risco))*(3.*sqrt(risco)-2.*a0);
    
    /** Dimensionless spin: J/Mbh^2 */
    abh = (4*nu)*(4*nu)*(L0 + L1*S + L2a*Delta*deltam + L2b*S2 + L2c*Delta2 + L2d*deltam2 + L3a*Delta*S*deltam + L3b*S*Delta2 + L3c*S3 + L3d*S*deltam2 + L4a*Delta*S2*deltam + L4b*Delta3*deltam + L4c*Delta4 + L4d*S4 + L4e*Delta2*S2 + L4f*deltam4 + L4g*Delta*deltam3 + L4h*Delta2*deltam2 + L4i*S2*deltam2) + S*(1+8*nu)*deltam4 + nu*Jisco*deltam6;
    
    Mbh = (4*nu)*(4*nu)*(M0 + K1*S + K2a*Delta*deltam + K2b*S2 + K2c*Delta2 + K2d*deltam2 + K3a*Delta*S*deltam + K3b*S*Delta2 + K3c*S3 + K3d*S*deltam2 + K4a*Delta*S2*deltam + K4b*Delta3*deltam + K4c*Delta4 + K4d*S4 + K4e*Delta2*S2 + K4f*deltam4 + K4g*Delta*deltam3 + K4h*Delta2*deltam2 + K4i*S2*deltam2) + (1 + nu*(Eisco + 11))*deltam6;
    
    a0 = abh;
    
  }

  *mass = Mbh;
  *spin = abh;
}

/** Final spin fit of */
double JimenezFortezaRemnantSpin(double nu, double X1, double X2, double chi1, double chi2)
{

  const double xnu     = sqrt(1.0-4.0*nu);
  const double Dchi    = chi1-chi2;
  const double S       = (X1*X1*chi1+X2*X2*chi2)/(X1*X1+X2*X2);
  const double a2      = 3.833;
  const double a3      = -9.49;
  const double a5      = 2.513;
  
  /** The functional form is taken from eq. (7), page 5. */
  double Lorb_spin_zero  = (1.3*a3*nu*nu*nu + 5.24*a2*nu*nu + 2.*sqrt(3)*nu)/(2.88*a5*nu + 1);
  
  /** Coeffcients taken from Table II, page 6: */
  double b1      = 1.00096;
  double b2      = 0.788;
  double b3      = 0.654;
  double b5      = 0.840;
  
  /** These values are taken from Table III, page 7: */
  double f21     = 8.774;
  double f31     = 22.83;
  double f50     = 1.8805;
  double f11     = 0.345225*f21 + 0.0321306*f31 - 3.66556*f50 + 7.5397;
  
  /** These values are taken from Table IV, page 10 */
  double f12     = 0.512;
  double f22     = -32.1;
  double f32     = -154;
  double f51     = -4.77;
  
  /** The following quantities were taken from the relation given in eq. (11), */
  /** page 7: fi3 = 64 - 64.*fi0 - 16.*fi1 - 4.*fi2; */
  double f13     = 64 - 16.*f11 - 4.*f12;
  double f23     = 64 - 16.*f21 - 4.*f22;
  double f33     = 64 - 16.*f31 - 4.*f32;
  double f53     = 64 - 64.*f50 - 16.*f51;
  
  /** this transformation is given in eq. (9), page (7) */
  double b1t     = b1*(f11*nu + f12*nu*nu + f13*nu*nu*nu);
  double b2t     = b2*(f21*nu + f22*nu*nu + f23*nu*nu*nu);
  double b3t     = b3*(f31*nu + f32*nu*nu + f33*nu*nu*nu);
  double b5t     = b5*(f50 + f51*nu + f53*nu*nu*nu);
  
  /** The functional form is taken from eq. (8), page 6. */
  double Lorb_eq_spin  = (0.00954*b3t*S*S*S + 0.0851*b2t*S*S - 0.194*b1t*S)/(1 - 0.579*b5t*S);
  
  /** These values are taken from Table IV, page 10: */
  double d10     = 0.322;
  double d11     = 9.33;
  double d20     = -0.0598;
  double d30     = 2.32;
  double d31     = -3.26;
  
  /** The functional form is taken from eq. (19a-c), page 10.*/
  double A1      = d10*xnu*nu*nu*(d11*nu+1);
  double A2      = d20*nu*nu*nu;
  double A3      = d30*xnu*nu*nu*nu*(d31*nu+1);
  
  /** The functional form is taken from eq. (15), page 9. */
  double Lorb_uneq_mass  = A1*Dchi + A2*Dchi*Dchi + A3*S*Dchi;
  
  return X1*X1*chi1+X2*X2*chi2 + Lorb_spin_zero + Lorb_eq_spin + Lorb_uneq_mass;
}

/** QNM fits for the 22 mode for spinning systems */
void QNMHybridFitCab(double nu, double X1, double X2, double chi1, double chi2, double aK, 
		     double Mbh, double abh,  
		     double *a1, double *a2, double *a3, double *a4, double *b1, double *b2, double *b3, double *b4, 
		     double *sigmar, double *sigmai)
{

  const double a12        = X1*chi1 - X2*chi2;
  const double X12        = X1 - X2;
  const double aeff       = aK + 1./3.*a12*X12;
  const double aeff_omg   = aK + a12*X12;
  const double af         = abh; 
  const double nu2        = SQ(nu);
  const double nu3        = nu2*nu;
  const double aeff2      = SQ(aeff);
  const double aeff3      = aeff2*aeff;
  const double af2        = SQ(af);
  const double af3        = af2*af;
  const double aeff_omg2  = SQ(aeff_omg); 
  const double aeff_omg3  = aeff_omg2*aeff_omg;
  const double aeff_omg4  = SQ(aeff_omg2);
  const double X12_2      = SQ(X12);

  double alpha21[KMAX], alpha1[KMAX], omega1[KMAX], c3A[KMAX], c3phi[KMAX], c4phi[KMAX], Domg[KMAX], Amrg[KMAX], c2A[KMAX];
      
  const int usespins = par_get_i("use_spins");  

  int modeon[KMAX];
  const int k21 = 0;
  const int k22 = 1;
  const int k33 = 4;
  const int k44 = 8;
 
  for (int k=0; k<KMAX; k++) {
    modeon[k] = 0; /* off */
    sigmar[k] = sigmai[k] = 0.;
    a1[k] = a2[k] = a3[k] = a4[k] = 0.;
    b1[k] = b2[k] = b3[k] = b4[k] = 0.;
  }

  if (!(usespins)) {

    modeon[k21]=1;
    modeon[k22]=1;
    modeon[k33]=1;
    
    /** Last updates: 05/09/2017 from CoM extrapolated SXS data */
    
    // l=2 -------------------------------------------------------------------
    
    /* (l=2, m=2)*/
    alpha21[k22] = -0.3025985041156393 *nu2 +  0.0032794155172817 *nu +  0.1828276903682022;
    alpha1[k22]  = -0.1615300454109702 *nu2 +  0.0147030662812516 *nu +  0.0878204175700328;
    c3A[k22]     =  0.8118901739129283 *nu  -  0.5584875090785957;
    c3phi[k22]   =  0.7156419884962878 *nu  +  3.8436474282409803;
    c4phi[k22]   =  2.2336960710670901 *nu  +  1.4736119175780844;
    Domg[k22]    =  0.8846304360111242 *nu2 +  0.0872792137250448 *nu +  0.1058414813686749;
    Amrg[k22]     = 1.4935750287318139 *nu2 +  0.2157497669089671 *nu +  1.4292027468283439;
    
    /* (l=2, m=1)*/
    alpha21[k21] = -0.2741607253846813 *nu2 +  0.0079342900879431 *nu +  0.1835522430667348;
    alpha1[k21]  = -0.1277546304610336 *nu2 +  0.0093615534859368 *nu +  0.0882855170502398;
    c3A[k21]     = -0.9431151070942140 *nu  +  0.2569989171628133;
    c3phi[k21]   = -3.4479482376671666 *nu  +  2.4755856452648359;
    c4phi[k21]   = -3.4024504071619841 *nu  +  1.0650118588151427;
    Domg[k21]    =  0.2660644668923829 *nu2 +  0.2276854484140649 *nu +  0.0884880283627388;
    Amrg[k21]    = -5.7236432632743952 *nu2 +  0.0390010969627653 *nu +  0.4291847351869338;
    
    // l=3 ------------------------------------------------------------------
    /* (l=3,m=3)*/
    alpha21[k33] = -0.3620553934265325 *nu2 +  0.0171973908686402 *nu +  0.1865364041200878;
    alpha1[k33]  = -0.1821867653548689 *nu2 +  0.0134440240947561 *nu +  0.0916720214797975;
    c3A[k33]     =  2.7565431398030675 *nu  -  0.5506682334306747;
    c3phi[k33]   = -0.2497526471104979 *nu  +  2.3737675006958683;
    c4phi[k33]   = -2.9538823110315420 *nu  +  1.4483501341373066;
    Domg[k33]    =  1.3341439550896721 *nu2 -  0.1717105341058959 *nu +  0.1694617455660599;
    Amrg[k33]    = -9.3034388918614841 *nu2 +  1.0189351143222705 *nu +  0.4533252110436300;
    
    // l=4 ------------------------------------------------------------------
    /* (l=4,m=4)*/
    /* alpha21[k44] = -0.3991680748908423 *nu2 +   0.0287698202159666 *nu +  0.1880112530796091; */
    /* alpha1[k44]  = -0.2003781755488581 *nu2 +   0.0171888841352427 *nu +  0.0930836242032652; */
    /* c3A[k44]     =  3.1899853343683140 *nu  +  -0.4131730594856833; */
    /* c3phi[k44]   = 31.5753575286023747 *nu  +  -1.0375600524681363; */
    /* c4phi[k44]   = 25.4170586178559716 *nu  +  -0.4151371540505313; */
    /* Domg[k44]    = -1.5342842283421341 *nu2 +   1.5224173843877831 *nu +  0.0897013049238634; */
    /* Amrg[k44]    =  0.9438333992719329 *nu2 +  -1.0464153920266663 *nu +  0.2897769169572948; */
    
    //sigma[k][0] = -0.208936*nu3 - 0.028103*nu2 - 0.005383*nu + 0.08896;
    //sigma[k][1] =  0.733477*nu3 + 0.188359*nu2 + 0.220659*nu + 0.37367;
    sigmar[k21] = -0.208936*nu3 - 0.028103*nu2 - 0.005383*nu + 0.08896;
    sigmai[k21] =  0.733477*nu3 + 0.188359*nu2 + 0.220659*nu + 0.37367;

    sigmar[k22] = -0.364177*nu3 + 0.010951*nu2 - 0.010591*nu + 0.08896;
    sigmai[k22] =  2.392808*nu3 + 0.051309*nu2 + 0.449425*nu + 0.37365;

    sigmar[k33] = -0.319703*nu3 - 0.030076*nu2-0.009034*nu + 0.09270;
    sigmai[k33] =  2.957425*nu3 + 0.178146*nu2 + 0.709560*nu + 0.59944;
    
  } else {

    modeon[k22]=1;
    
    /** Setting up coefficients from the phenomenological description of the ringdown.
	For notation: Damour&Nagar, PRD 90 (2015), 024054 and Del Pozzo & Nagar, PRD 95 (2017), 124034
	Current global fits are new. See Nagar+ 2017 (in preparation) for a global performance
	and Riemenschneider& Nagar (2017) in preparation for the description of the fits */
    
    /* omg1 - imaginary part of the fundamental mode */
    double omega1_c    = -0.0598837831 * af3 + 0.8082136788 * af2 - 1.7408467418 * af + 1;
    double omega1_d    = -0.2358960279 * af3 + 1.3152369374 * af2 - 2.0764065380 * af + 1;
    omega1[k22]        =  0.3736716844 * (omega1_c/omega1_d);
    
        /* alpha1 - real part (damping time) of the fundamental mode */
    double alpha1_c    =  0.1211263886 * af3 + 0.7015835813 * af2 - 1.8226060896 * af + 1;
    double alpha1_d    =  0.0811633377 * af3 + 0.7201166020 * af2 - 1.8002031358 * af + 1;
    alpha1[k22]        =  0.0889623157 * (alpha1_c/alpha1_d);
    
    /* alpha2 - alpha1 */
    double alpha21_c   =  0.4764196512 * af3 - 0.0593165805 * af2 - 1.4168096833 * af + 1;
    double alpha21_d   =  0.4385578151 * af3 - 0.0763529088 * af2 - 1.3595491146 * af + 1;
    alpha21[k22]       =  0.1849525596 * (alpha21_c/alpha21_d);
    
    /* c3A */
    double a_c3A 	=  0.0169543;
    double b_c3A 	= -0.0799343;
    double c_c3A 	= -0.115928;
    double c3A_nu       =  0.8298678603 * nu - 0.5615838975;
    double c3A_eq       =  (c_c3A * X12 + 0.0907476903) * aeff3 + (b_c3A * X12 + 0.0227344099) * aeff2 + (a_c3A * X12 - 0.1994944332)*aeff;
    c3A[k22]            =  c3A_nu + c3A_eq;
    
    /* c3_phi */
    double a_c3phi      = -0.462321;
    double b_c3phi      = -0.904512;
    double c_c3phi      =  0.437747;
    double d_c3phi      =  1.8275;
    double c3phi_nu     =  0.4558467286 * nu + 3.8883812141;
    double c3phi_equal  =  (d_c3phi*X12-2.0575868122) * aeff_omg4 +(c_c3phi*X12-0.5051534498)*aeff_omg3 +(b_c3phi*X12+2.5742292762)*aeff_omg2 +(a_c3phi*X12+2.5599640181)*aeff_omg;
    c3phi[k22]          = c3phi_nu + c3phi_equal;
    
    /* c4_phi */
    double a_c4phi      = -0.449976;
    double b_c4phi      = -0.980913;
    double c4phi_nu     =  2.0822327682 * nu + 1.4996868401;
    double c4phi_equal  =  (b_c4phi*X12+3.5695199109) * aeff_omg2 + (a_c4phi * X12 + 4.1312404030) * aeff_omg;
    c4phi[k22]          =  c4phi_nu + c4phi_equal;
    
    /* omg_mrg: the "merger frequency", i.e. the frequency at the peak of |h22| */
    /* Special scaling and independent variables used for the fit. AN&GR 2017 */
    double a2_omgmx     = -0.122735;
    double a1_omgmx     =  0.0857478;
    double b2_omgmx     = -0.0760023;
    double b1_omgmx     =  0.0826514;
    double omgmx_eq_c   =  (a2_omgmx*X12_2 +a1_omgmx*X12 -0.1416002395) * aeff_omg + 1;
    double omgmx_eq_d   =  (b2_omgmx*X12_2 +b1_omgmx*X12 -0.3484804901) * aeff_omg + 1;
    double omgmx_eq     =  omgmx_eq_c/omgmx_eq_d;
    double omgmx        =  (0.481958619443355 * nu2 + 0.223976694441952 * nu + 0.273813064427363) * omgmx_eq;
    
    /* the peak of the h22 metric (strain) waveform.*/
    /* Special scaling and independent variables used for the fit. AN& GR 2017*/	
    double a2_A_scaled = -0.0820894;
    double a1_A_scaled = 0.176126;
    double b2_A_scaled = -0.150239;
    double b1_A_scaled = 0.20491;
    double A_scaled_eq = ((a2_A_scaled*X12*X12 + a1_A_scaled*X12 -0.2935238329)*aeff + 1)/((b2_A_scaled*X12*X12 + b1_A_scaled*X12 -0.4728707630)*aeff + 1);
    double A_scaled    = (+1.826573640739664*nu2 +0.100709438291872*nu +1.438424467327531)*A_scaled_eq;
    
    Amrg[k22]      = A_scaled*(1-0.5*omgmx*aeff);
    Domg[k22]      = omega1[k22] - Mbh*omgmx;
    
    /* renaming real & imaginary part of the QNM complex frequency sigma */
    //sigma[k22][0] = alpha1[k22];
    //sigma[k22][1] = omega1[k22];
    sigmar[k22] = alpha1[k22];
    sigmai[k22] = omega1[k22];
  
  }

  for (int k=0; k<KMAX; k++) {
    if (modeon[k]) {
      c2A[k] = 0.5*alpha21[k];
      double cosh_c3A = cosh(c3A[k]);  
      a1[k] = Amrg[k] * alpha1[k] * cosh_c3A * cosh_c3A / c2A[k];
      a2[k] = c2A[k];
      a3[k] = c3A[k];
      a4[k] = Amrg[k] - a1[k] * tanh(c3A[k]);
      b2[k] = alpha21[k];
      b3[k] = c3phi[k];
      b4[k] = c4phi[k];
      b1[k] = Domg[k] * (1+c3phi[k]+c4phi[k]) / (b2[k]*(c3phi[k] + 2.*c4phi[k]));
    }
  }
  
}

/** Light-ring approximation, vacuum nonspinning case */
double eob_approxLR(const double nu) 
{
  const double x = 1.-4.*nu;
  const double r0 = 1.818461553848201e+00; // x=0 (q=1)
  const double r1 = 1.821043720041472e+00; // x=1 (3M)
  return r0*x + r1;
}

/** Compute optimized timestep after merger */
double get_mrg_timestep(double q, double chi1, double chi2)
{
  double dt = 0.1;
  // ...
  return dt;
}

