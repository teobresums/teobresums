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
  double a6   =  EOBPars->a6c;
  
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
  
  /* Numerator and denominator of the Pade */
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

/** EOB Metric A function GSF-informed */

void eob_metric_AGSF(double r, double nu, double *A, double *dA, double *d2A)
{

  /* shortcuts */
  double nu2    = nu*nu;
  double pi2    = Pi*Pi;
  double pi4    = pi2*pi2;
  double u      = 1./r;
  double u2     = u*u;
  double u3     = u*u2;
  double u4     = u2*u2;
  double u5     = u4*u;
  double u6     = u5*u;
  double u7     = u6*u;
  double u10    = u5*u5;
  double u8     = u5*u3;
  double u9     = u8*u;
  double logu   = log(u);
  double logu2  = logu*logu;
  double logu3  = logu2*logu;
  double logu4  = logu3*logu;
  double logu5  = logu4*logu;
  double logu6  = logu5*logu;
  double logu7  = logu6*logu;
  double logu8  = logu7*logu;
  double logu9  = logu8*logu;
  double logu10 = logu9*logu;
  double logu11 = logu10*logu;
  double logu12 = logu11*logu;
  double logu13 = logu12*logu;
  double logu14 = logu13*logu;
  double logu15 = logu14*logu;
  double logu16 = logu15*logu;
  double logu17 = logu16*logu;
  double logu18 = logu17*logu;
  double logu19 = logu18*logu;
  double logu20 = logu19*logu;
  double logu21 = logu20*logu;
  double logu22 = logu21*logu;
  double logu23 = logu22*logu;
  double u3by2  = sqrt(u3);
  double u5by2  = sqrt(u5);
  double u7by2  = sqrt(u7);
  double log2   = log(2.);
  double log3   = log(3.);

  /* The A function is factorized into two pieces:
  integer and tail part (semi-integer powers of u) */

  double den   = 2051.88197956336 + 857.7557817914457*logu + 148.1965571937127*logu2 - 27.90426097620334*logu3 - 1.*logu4;

  /* Integer part */
  double cNUM1 = (21180.72520413395 + 8187.222705098468*logu + 1103.683057935346*logu2 - 358.7178617813733*logu3 - 9.367468814953554*logu4)/den;   
  double cNUM2 = (57871.66044469341 + 26711.92378791032*logu + 4122.408176113611*logu2 - 771.266342096333*logu3 - 184.8981051714863*logu4 - 6.4*logu5)/den;     
  double cNUM3 = (52990.49627106789 + 128389.8685809894*logu + 52506.56659502195*logu2 + 10166.34892393265*logu3 - 422.3415307067122*logu4 - 29.12567868268803*logu5)/(-den); 

  double cDEN1 = (982.0988270772413 - 256.4834786521463*logu - 355.1566409048896*logu2 - 84.02967435153155*logu3 + 0.4764825322652425*logu4)/den; 
  double cDEN2 = (14760.97950381165 + 2124.411322949478*logu - 286.4909581835134*logu2 - 437.7379326226858*logu3 + 5.297364505707555*logu4)/den; 
  double cDEN3 = (82438.45336374756 + 34414.29393768581*logu + 6868.976797992799*logu2 - 460.4929791137166*logu3 + 5.427034905560372*logu4 + 4.076190476190476*logu5)/(-den); 

  double NumInt = 1. + cNUM1*u + cNUM2*u2 + cNUM3*u3;
  double DenInt = 1. + cDEN1*u + cDEN2*u2 + cDEN3*u3;

  double ResumInt = NumInt/DenInt;

  /* Tail part */

  double NumTail = 1. + 14.68963100342707*u + 2.442129652761482e-14*(3.679688615621096e15 + 2.62066348228608e14*logu)*u2 + 1.187543790303097e-22*(2.329829442816534e24 + 7.589067344988391e22*logu)*u3 + 40.97833617482458*u7by2;
  double DenTail = 1. + 14.68963100342707*u + 2.442129652761482e-14*(3.679688615621096e15 + 2.62066348228608e14*logu)*u2 + 1.187543790303097e-22*(2.329829442816534e24 + 7.589067344988391e22*logu)*u3;

  double ResumTail = NumTail/DenTail;

  /* a1SF function purely analytical */
  double a1SF_tmp = ResumInt*ResumTail;

  /* Schwarzschild Hamiltonian */
  double Hschw = (1. - 2.*u)/sqrt(1. - 3.*u);

  /* First derivatives */ 

  /* of the Schwarzschild Hamiltonian */
  double dHschw_u = -2./sqrt(1. - 3.*u) + (3.*(1. - 2.*u))/(2.*pow(1. - 3.*u, 1.5));

  /* of the integer part */
  double NumDResumInt = -5.118204693470436e9 - 2.517552739621222e10*u + 2.169722891287838e10*u2 - 4.413372601505725e11*u3 - 2.464425852212004e11*u4 - 1.081879419281513e12*u5 + 
                        logu*(-6.418745881013987e9 - 3.490011790031027e10*u + 1.028325989896825e11*u2 - 5.186443254518964e11*u3 - 1.920423399415217e11*u4 - 1.431457710853209e12*u5) + 
                        logu2*(-3.792234354872924e9 - 2.209171382512806e10*u + 1.162712175515147e11*u2 - 3.070776803447059e11*u3 - 6.838930005618407e10*u4 - 9.212539531809577e11*u5) + 
                        logu3*(-1.0922677718373e9 - 7.111695102247094e9*u + 6.766204471277715e10*u2 - 1.032992477934339e11*u3 - 8.010683626671883e9*u4 - 3.210177223952128e11*u5) + 
                        logu4*(-9.182825467188085e7 - 7.929907976461177e8*u + 2.109470329557746e10*u2 - 2.301568792721061e10*u3 - 4.755998812616166e9*u4 - 5.681544205408146e10*u5) + 
                        logu5*(3.942701724425052e7 + 2.299118683070849e8*u + 2.806216066588674e9*u2 - 3.75534035039741e9*u3 - 4.11942353797047e9*u4 + 5.469918530074278e7*u5) + 
                        logu6*(1.022972472584576e7 + 8.871918762561751e7*u - 3.379061074805489e8*u2 - 6.765317369960402e8*u3 - 1.661481447616641e9*u4 + 1.815288220018964e9*u5) +
                        logu7*(150491.7723380209 + 6.709806062411051e6*u - 1.576201448773277e8*u2 - 6.239066492363331e7*u3 - 2.099718235240898e8*u4 + 2.750763213323691e8*u5) +
                        logu8*(-254793.0308863475 - 1.572310528895074e6*u - 1.241586966402864e7*u2 + 9.574241801152974e6*u3 + 1.795680321608006e7*u4 - 6.291657800031504e6*u5) + 
                        logu9*(-3351.894734917506 - 267138.1016475637*u + 2.22013888531215e6*u2 + 4.55445574698243e6*u3 + 8.928747190607484e6*u4 - 2.183876868770208e6*u5) + 
                        logu10*(1120.555225658835 + 5786.04278981747*u + 157624.6619596899*u2 + 19120.10231997654*u3 - 59090.42397756783*u4 - 174468.4943790704*u5) + 
                        logu11*(49.59664489934719 + 1721.587023312448*u - 3607.272313296524*u2 - 17697.9147801814*u3 - 27515.5517263019*u4 - 747.9778266824749*u5)  + 
                        logu12*(0.5924620250845447 + 66.27253867316881*u - 389.0492017114477*u2 - 537.6433047794397*u3 - 1571.172723544291*u4 + 94.0985561765124*u5) + 
                        logu13*u*(0.7703729583369727 - 6.178336699201579*u - 2.925697147745601*u2 - 81.97727104890009*u3) - 1.570093457943925*logu14*u4;

   
   double fDInt = -503.3822613414784 - 240.9354599138092*u - 3621.268336215476*u2 + 20224.38692334929*u3 +
                  logu*(-210.4307408600509 + 62.92234873475551*u - 521.1756750226522*u2 + 8442.759026768716*u3) + 
                  logu2*(-36.35663202182205 + 87.12954975470422*u + 70.2839967506283*u2 + 1685.146177077673*u3) + 
                  logu3*(6.84567150117138 + 20.61475655820283*u + 107.3889787976215*u2 - 112.9714084274305*u3) + 
                  logu4*(0.2453271028037383 - 0.1168940791772207*u - 1.299587086680592*u2 + 1.331398750195886*u3) + 
                  logu5*u3;

   double DenDResumInt = -den*fDInt*fDInt;

  double dResumInt = NumDResumInt/DenDResumInt;

  double fDTail = 8.420742107916455e21 + 1.236975943403134e23*u + 7.567111844223912e23*u2 + 5.621323169503642e17*logu*(95872. + 135005.*u)*u2 + 2.329829442816534e24*u3;

  double dResumTail = (1.725340004698507e23*(5.894519475541519e22 + 6.184879717015671e23*u + 2.162348054285843e24*u2 + 5.621323169503642e17*logu*(287616. + 135005.*u)*u2 + 2.178048095916766e24*u3)*u5by2)/(fDTail*fDTail);

  /* derivative of a1SF with respect to u */
  double da1SF_u_tmp = dResumInt*ResumTail + ResumInt*dResumTail; 

  /* GSF-informed correction and improved functions */
  double f1SF    = 1 + u7*(963.433 + 811.383*logu - 205.256*logu2);
  double df1SF_u = u6*(7555.41 + 5269.17*logu - 1436.79*logu2); 

  double a1SF    = a1SF_tmp*f1SF;
  double da1SF_u = da1SF_u_tmp*f1SF + a1SF_tmp*df1SF_u;

  *A          = 1. - 2.*u + 2.*nu*u3*Hschw*a1SF;
  double dA_u = -2. + 6.*nu*u2*a1SF*Hschw + 2.*nu*u3*dHschw_u*a1SF + 2.*nu*u3*Hschw*da1SF_u;
  *dA         = dA_u;

  /* Second derivatives */

  /* of the Schwarzschild Hamiltonian */
  double d2Hschw_u = -6./pow(1. - 3.*u, 1.5) + (27.*(1. - 2.*u))/(4.*pow(1. - 3.*u, 2.5));


  /* of the integer part */
  
   double d2ResumInt_Num = -2.781672130920104e16 + 2.480967515433256e17*u - 2.119819245405735e18*u2 - 6.151188461573498e18*u3 - 3.812978833531141e18*u4 - 5.887074391136803e19*u5 - 3.829883408018222e19*u6 - 4.174514962629836e19*u7 + 
   logu22*u3*(-5.392610708358809 + 25.08041605478338*u + 8.777091443236803*u2 + 255.7599782479884*u3) + 3.14018691588785*logu23*u6 + 
   logu2*(-7.302337491971747e16 + 7.749856020668733e17*u - 4.536430813976029e18*u2 - 1.521334974495676e19*u3 + 1.542177195408716e19*u4 - 1.145429786298632e20*u5 - 5.52331248729711e19*u6 - 9.650131997859638e19*u7) + 
   logu*(-6.50156152291766e16 + 6.329783360595484e17*u - 4.402043640276021e18*u2 - 1.403027319090653e19*u3 + 2.812197677454418e18*u4 - 1.170096731707718e20*u5 - 6.624386740349145e19*u6 - 9.020299949920489e19*u7) + 
   logu3*(-4.972114285445602e16 + 5.818629455571703e17*u - 2.898929144906484e18*u2 - 9.90882999704736e18*u3 + 1.998325016128144e19*u4 - 6.89799884883713e19*u5 - 2.635650951596912e19*u6 - 6.376112794968176e19*u7) + 
   logu4*(-2.170550138562532e16 + 2.875480348965157e17*u - 1.241536781424846e18*u2 - 4.059941051327444e18*u3 + 1.419885553092043e19*u4 - 2.768767847985642e19*u5 - 7.249235340601194e18*u6 - 2.790211691369108e19*u7) + 
   logu5*(-5.679823459453885e15 + 9.178877487631606e16*u - 3.557426366324769e17*u2 - 9.380974411479764e17*u3 + 6.322785438227277e18*u4 - 7.508274425432074e18*u5 - 9.114550213721848e17*u6 - 7.946731176429188e18*u7) + 
   logu6*(-5.047690792268042e14 + 1.554521608676324e16*u - 6.243043602261672e16*u2 - 2.666356729071885e16*u3 + 1.724242911971985e18*u4 - 1.371971242179132e18*u5 - 5.634033933868715e16*u6 - 1.256784884847882e18*u7) + 
   logu11*(-5.306847305704753e11 + 6.253279252719161e12*u + 5.772758721537055e12*u2 - 1.083228391778925e14*u3 + 1.518639264300058e14*u4 + 3.062621664460186e14*u5 + 4.106007319348497e14*u6 - 1.325028820802392e14*u7) + 
   logu12*(-2.282501982123623e10 + 8.309007855476562e11*u - 9.867237828462509e11*u2 - 6.272287227604082e12*u3 + 6.096444114076848e13*u4 + 3.778639630599216e13*u5 + 5.495825133154462e13*u6 - 1.304032852625348e13*u7) + 
   logu15*(-6.734754074389435e7 - 1.403634723835746e8*u + 5.234804348968131e9*u2 - 7.382838334102685e9*u3 - 1.208134438862878e10*u4 + 3.97484756995695e10*u5 + 6.376355152343071e10*u6 - 6.87050945676072e9*u7) + 
   logu16*(-9.556742362152208e6 + 1.09893048651619e8*u + 4.691809801858975e8*u2 - 1.256658246746765e9*u3 + 5.8429498350217e9*u4 + 5.151313952315446e9*u5 + 4.385815183742413e9*u6 - 6.996676469152076e8*u7) + 
   logu20*(27.13339429214482 - 364.2153704188829*u - 825.5506382242851*u2 - 24545.77289755097*u3 + 50071.16293057652*u4 + 98923.90322166853*u5 + 106161.4458583436*u6 - 2191.252654792529*u7) + 
   logu21*(0.1889933659471545 - 2.941374847534438*u - 2.70453031957207*u2 - 627.663562075684*u3 + 2246.124935000245*u4 + 1876.540616673677*u5 + 7266.195831026774*u6 - 94.0985561765124*u7) + 
   logu19*(1438.256772439385 - 15437.04173106886*u - 66150.72026509648*u2 - 268876.3899561888*u3 - 733873.231412491*u4 + 1.140013125257407e6*u5 + 1.008334501442573e6*u6 + 217033.9097392946*u7) + 
   logu18*(30010.55813123302 - 178113.168214963*u - 1.968865318923377e6*u2 + 3.7869306215258e6*u3 - 3.364552565504599e7*u4 - 3.067273835106073e7*u5 - 2.065160143492527e7*u6 + 8.077212576591474e6*u7) + 
   logu17*(-32593.92532632617 + 4.382743427464224e6*u - 1.061485556435447e7*u2 + 4.167949610238702e7*u3 - 2.785523652360581e7*u4 - 4.824174221204235e8*u5 - 5.786057537162013e8*u6 + 4.826363507307306e7*u7) + 
   logu14*(1.008661833309523e9 - 1.674250207602155e10*u - 3.258437262950869e10*u2 + 1.890757007272333e11*u3 - 7.82767074686612e11*u4 - 8.986832409444984e11*u5 - 1.120826633549568e12*u6 + 4.804082438956363e10*u7) + 
   logu13*(9.180436686067664e9 - 4.808439555916074e10*u - 3.6063716274591e11*u2 + 1.400194307958538e12*u3 + 1.964532005334606e12*u4 - 2.642410768352027e12*u5 - 3.367971355170659e12*u6 + 9.117076236751579e11*u7) + 
   logu10*(-1.56199755143428e12 - 2.145999411703856e12*u + 7.033290830585901e13*u2 - 3.208755969038362e14*u3 - 2.038267230103109e15*u4 + 4.978646696144575e14*u5 - 1.118487175397344e14*u6 + 2.168235817748831e14*u7) + 
   logu7*(2.118995835105034e14 - 6.597666702003277e14*u - 3.504236206040361e15*u2 + 5.860005166500923e16*u3 + 2.130066103019003e17*u4 - 1.766475571189837e17*u5 - 8.210724855563997e16*u6 + 1.939839886244514e15*u7) + 
   logu9*(9.253859190079912e12 - 2.115059207990373e14*u + 4.35078393347501e14*u2 + 1.875870855765536e15*u3 - 1.605032939210946e16*u4 - 2.590498743288691e15*u5 - 1.089805281078113e16*u6 + 8.910128460870063e15*u7) + 
   logu8*(8.492202298842821e13 - 1.029092880540039e15*u + 1.298086876633344e15*u2 + 1.87885792794751e16*u3 - 2.915049315910143e16*u4 - 2.276509023314673e16*u5 - 4.972860774439563e16*u6 + 4.704998193388095e16*u7);

   double d2ResumInt_Den = den*den*fDInt*fDInt*fDInt; 

   double d2ResumInt = d2ResumInt_Num/d2ResumInt_Den;

  /* of the tail part */

  /*  some factors */
  double fD2Tail1 = 5.894519475541519e22 + 6.184879717015671e23*u + 2.162348054285843e24*u2 + 5.621323169503642e17*logu*(287616. + 135005.*u)*u2 + 2.178048095916766e24*u3;
  double fD2Tail2 = 8.420742107916455e21 + 1.236975943403134e23*u + 7.567111844223912e23*u2 + 5.621323169503642e17*logu*(95872. + 135005.*u)*u2 + 2.329829442816534e24*u3;

  double d2ResumTail = (8.626700023492536e22*(-12.*u*(5.897984982238551e22 + 5.890863629721913e23*u + 1.873774389834547e17*logu*u*(191744. + 405015.*u) - 2.188678318287006e22*u2 + 1.585655577881217e8*(-1.11924169555968e14 - 4.20318202284288e14*u + 1.499072781193489e16*u2))*(fD2Tail1) + 
                                                6.*u*(2.948992491119276e23 + 1.695402089595687e24*u + 5.621323169503642e17*logu*u*(191744. + 135005.*u) - 1.736681300826379e23*u2 + 1.585655577881217e8*(-5.5962084777984e14 - 1.260954606852864e15*u + 1.499072781193489e16*u2))*(fD2Tail2) + 
                                                5.*(fD2Tail1)*(fD2Tail2))*u3by2)/(fDTail*fDTail*fDTail);
  
  double d2a1SF_u_tmp = d2ResumInt*ResumTail + 2*dResumInt*dResumTail + ResumInt*d2ResumTail;

  /* adding the GSF correction */
  double d2f1SF_u = u5*(50601.6 + 28741.4*logu - 8620.74*logu2);

  double d2a1SF_u     = d2a1SF_u_tmp*f1SF + 2.*df1SF_u*da1SF_u_tmp + d2f1SF_u*a1SF_tmp;
  double d2A_u = da1SF_u*(12.*u2*nu*Hschw + 2.*nu*u3*dHschw_u) + 2.*nu*u3*Hschw*d2a1SF_u + 2.*nu*u3*dHschw_u*da1SF_u + a1SF*(12.*nu*u*Hschw + 12.*nu*u2*dHschw_u + 2.*nu*u3*d2Hschw_u);
  *d2A         = d2A_u;

}

/** EOB Metric D function at 3PN, resummed */

void eob_metric_D3PN(double r, double nu, double *D, double *dD, double *d2D)
{

  /* shortcuts */
  double u  = 1./r;
  double u2 = u*u;
  double u3 = u2*u;

  double Dp       = 1.0 + 6.*nu*u2 - 2.*(3.*nu-26.)*nu*u3; // Pade' resummation of D
  double dDp_du   = 6.*nu*u*(2. - (3.*nu-26.)*u);
  double d2Dp_du2 = 12.*nu*(1. - (3.*nu-26.)*u);
  double D_tmp    = 1./Dp;
  
  /* derivatives wrt to u */
  *D   = D_tmp;
  *dD  = -SQ(D_tmp)*dDp_du;
  *d2D = 2.*SQ(D_tmp)*D_tmp*SQ(dDp_du) - SQ(D_tmp)*d2Dp_du2;

}

/** EOB Metric D function GSF-informed */

void eob_metric_DGSF(double r, double nu, double *D, double *dD, double *d2D)

{

  /* shortcuts */
  double nu2    = nu*nu;
  double pi2    = Pi*Pi;
  double pi4    = pi2*pi2;
  double u      = 1./r;
  double u2     = u*u;
  double u3     = u*u2;
  double u4     = u2*u2;
  double u5     = u4*u;
  double u6     = u5*u;
  double u7     = u6*u;
  double u10    = u5*u5;
  double u8     = u5*u3;
  double u9     = u8*u;
  double u3by2  = sqrt(u3);
  double u5by2  = sqrt(u5);
  double u7by2  = sqrt(u7);
  double u11by2 = sqrt(u10*u);
  double u13by2 = sqrt(u10*u3);
  double u15by2 = sqrt(u10*u5);
  double logu   = log(u);
  double log2   = log(2.);
  double log3   = log(3.);
  double logu2  = logu*logu;
  double logu3  = logu2*logu;
  double logu4  = logu3*logu;
  double logu5  = logu4*logu;
  double logu6  = logu5*logu;
  double logu7  = logu6*logu;
  double logu8  = logu7*logu;
  double logu9  = logu8*logu;
  double logu10 = logu9*logu;
  double logu11 = logu10*logu;
  double logu12 = logu11*logu;
  double logu13 = logu12*logu;
  double logu14 = logu13*logu;
  double logu15 = logu14*logu;
  double logu16 = logu15*logu;
  double logu17 = logu16*logu;
  double logu18 = logu17*logu;
  double logu19 = logu18*logu;
  double logu20 = logu19*logu;
  double logu21 = logu20*logu;
  double logu22 = logu21*logu;
  double logu23 = logu22*logu; 

  /* The d function is factorized into two pieces:
  integer and tail part (semi-integer powers of u) */

  /* Integer part */

  double Num_d1SFInt = -4.899034652014102e35 - 5.371367422700039e36*u - 3.819509358031804e37*u2 - 1.441829657233304e37*u3 + 
                        logu*(-1.891500947547505e35 - 1.906893527277763e36*u - 1.594035352919957e37*u2 + 1.061680081899536e37*u3) +
                        logu2*(-2.897217981124678e34 - 2.588448404744833e35*u - 2.805757109454611e36*u2 + 6.278884354865547e36*u3) +
                        logu3*(2.319421673577702e33 + 3.475562446777183e34*u + 1.166947884436259e35*u2 + 1.52323780670844e36*u3) + 
                        logu4*(1.161291662015839e32 + 1.290586003749631e33*u + 2.286657612744122e34*u2 + 3.228978016726249e34*u3) + 
                        logu5*(7.638718487926409e32 - 1.0324503849308e31*u)*u2;

  double Den_d1SFInt = -4.899034652014102e35 - 1.125537390954484e36*u - 1.034895514594423e37*u2 + 7.769486956780932e37*u3 + 
                        logu*(-1.891500947547505e35 - 2.675927060699253e35*u - 3.41368009510312e36*u2 + 2.581126241057948e37*u3) + 
                        logu2*(-2.897217981124678e34 - 7.752615443677862e33*u - 4.244765573412268e35*u2 + 3.294242589442275e36*u3) + 
                        logu3*(2.319421673577702e33 + 1.465396996343174e34*u + 9.461312952527907e34*u2 - 5.811173613250686e35*u3) + 
                        logu4*(1.161291662015839e32 + 2.841332300025704e32*u + 8.589526938685824e32*u2 + 5.659350379307624e33*u3) + 
                        logu5*2.046982059584161e33*u3;                       

  double d1SFInt = Num_d1SFInt/Den_d1SFInt;

  /* Tail part */

  double Num_d1SFTail = 3.467367151833754e27 + 5.492753720117619e28*u + 5.290545158749652e29*u2 + 4.154845306616869e19*logu*(5.48939104e8 + 3.918241855e9*u)*u2 + 1.829586553431019e30*u3 + 3.053882426133018e29*u7by2;
  double Den_d1SFTail = 3.467367151833754e27 + 5.492753720117619e28*u + 5.290545158749652e29*u2 + 4.154845306616869e19*logu*(5.48939104e8 + 3.918241855e9*u)*u2 + 1.829586553431019e30*u3;
  double d1SFTail     = Num_d1SFTail/Den_d1SFTail;

  /* d1SF */

  double d1SF_tmp = d1SFInt*d1SFTail;

  /* GSF-informed correction and improved function */

  double f1SF = 1.*(1. + (-9.898639062741917e6 + 2.6012818800187605e6*logu)*u15by2 + 
                         (1.0018072583960408e7 + 2.6340749138477915e6*logu + 198577.79085898417*logu2)*u7);
  
  double d1SF = d1SF_tmp*f1SF;

  /* First derivatives */

  /* of the integer part */

  double Dd1SF_Int = (-1.019022157248837e108 - 1.41398591658835e109*u - 5.972962305277044e109*u2 - 4.293368527973855e110*u3 - 1.325220602340017e111*u4 - 5.864268870835916e110*u5 + 
                    logu*(-1.180323173596444e108 - 1.723025531544031e109*u - 5.39019335672999e109*u2 - 4.42007147492885e110*u3 - 1.411692164939391e111*u4 - 7.50941887625607e110*u5) + 
                    logu2*(-6.365093162041004e107 - 9.842775233739604e108*u - 1.909856460903128e109*u2 - 2.088225359054651e110*u3 - 6.907763538653848e110*u4 - 4.610828293268436e110*u5) + 
                    logu3*(-1.837821898027187e107 - 3.102912220003404e108*u - 4.798257691855606e107*u2 - 4.745231850755357e109*u3 - 1.650915498648733e110*u4 - 1.643284330819324e110*u5) + 
                    logu4*(-2.574122122956699e106 - 5.155748105722514e107*u + 1.752901463013083e108*u2 - 2.950281688263012e108*u3 - 1.137699573781095e109*u4 - 3.523801850520969e109*u5) + 
                    logu5*(3.010170247358785e104 - 1.605345012586149e106*u + 5.930679631349931e107*u2 + 1.012695147735136e108*u3 + 3.908120786967131e108*u4 - 3.883153364054411e108*u5) + 
                    logu6*(5.754024584321218e104 + 9.193768830359903e105*u + 6.901177454847088e106*u2 + 1.571272052824783e107*u3 + 8.342105678692772e107*u4 - 1.215942066268117e106*u5) + 
                    logu7*(5.039325597214167e103 + 1.365062148855902e105*u - 2.254151277414218e105*u2 - 1.40506432118166e106*u3 + 1.300920119261987e105*u4 + 5.157021248913295e106*u5) +
                    logu8*(-4.339080797143522e102 - 1.240071562740252e103*u - 1.269773512540493e105*u2 - 4.657412698053918e105*u3 - 1.522130274577523e106*u4 + 4.540214939989015e105*u5) +  
                    logu9*(-3.639754585836174e101 - 1.132048599724896e103*u - 3.213189206629774e103*u2 + 8.045408013366297e103*u3 - 4.055332885668775e104*u4 + 4.588559283895894e103*u5) + 
                    logu10*(6.084616370125772e99 - 3.676502495475943e101*u + 4.323847172137712e102*u2 + 2.501052663247333e103*u3 + 7.214491425303449e103*u4 - 2.018936626349478e103*u5) + 
                    logu11*(8.132717273867246e98 + 2.286266731488472e100*u + 1.945468008307626e101*u2 + 2.888188284581554e101*u3 + 4.105889405095191e102*u4 - 8.762387987027005e101*u5) + 
                    logu12*(1.357300524129673e97 + 1.426893628578643e99*u + 1.791343186981558e98*u2 - 2.862393807419481e100*u3 - 6.635782346409161e100*u4 - 7.68252855465966e99*u5)  + 
                    logu13*u*(2.060312590473761e97 - 5.802957152660475e97*u - 6.142628859334867e98*u2 - 9.565492862182726e99*u3) - 
                    logu14*1.815832769570697e98*u4)/((-4.899034652014102e35 - 1.891500947547505e35*logu - 2.897217981124678e34*logu2 + 2.319421673577702e33*logu3 + 1.161291662015839e32*logu4)*Den_d1SFInt*Den_d1SFInt);

  /* of the tail part */
  
  double Dd1SF_Tail = (1.526941213066509e29*(2.427157006283628e28 + 2.746376860058809e29*u + 1.541548406427438e30*u2 + 4.154845306616869e19*logu*(1.646817312e9 + 3.918241855e9*u)*u2 + 1.503992777802289e30*u3)*u5by2)/(Den_d1SFTail*Den_d1SFTail);

  /* D(d1SF) */

  double Dd1SF_tmp = d1SFInt*Dd1SF_Tail + d1SFTail*Dd1SF_Int;

  /* adding the GSF correction */

  double Df1SF = 1.*(2.6012818800187605e6*u13by2 + 
                 7.5*(-9.898639062741917e6 + 2.6012818800187605e6*logu)*u13by2 + 
                 7.*(1.0018072583960408e7 + 2.6340749138477915e6*logu + 198577.79085898417*logu2)*u6 + 
                 (2.6340749138477915e6/u + (397155.58171796834*logu)/u)*u7);

  double Dd1SF = d1SF_tmp*Df1SF + f1SF*Dd1SF_tmp;

  /* Second derivatives */ 

  /* of the integer part */

  /* some factors */

  double fD2Int = 4.899034652014102e35 + 1.891500947547505e35*logu + 2.897217981124678e34*logu2 - 2.319421673577702e33*logu3 - 1.161291662015839e32*logu4;

  double D2d1SF_Int = (-2.641101451551516e180 + 9.515208764625697e178*u - 3.230643783622936e182*u2 - 3.470652067064428e183*u3 - 9.549640873937975e183*u4 - 5.060505333409733e184*u5 - 1.038118798424198e185*u6 - 1.718696019161502e184*u7 +  
                      logu*(-5.469848301307015e180 + 6.049977967745608e180*u - 5.703138110142176e182*u2 - 6.712294550366805e183*u3 - 1.584946569302385e184*u4 - 8.867370385165948e184*u5 - 1.865080762101521e185*u6 - 3.402679321376194e184*u7) + 
                      logu2*(-5.434745415956582e180 + 1.178284018448924e181*u - 4.782118833780379e182*u2 - 6.23738337988015e183*u3 - 1.217182289418095e184*u4 - 7.408920636948101e184*u5 - 1.600865114263083e185*u6 - 3.259169384961876e184*u7) + 
                      logu3*(-3.326872248969142e180 + 1.114523011273131e181*u - 2.380740435588645e182*u2 - 3.543760691773076e183*u3 - 5.179621170723038e183*u4 - 3.682947277388125e184*u5 - 8.224430998891984e184*u6 - 1.934523175589101e184*u7) + 
                      logu4*(-1.359059879933684e180 + 6.504755519519597e180*u - 7.340588405886695e181*u2 - 1.321349047118826e183*u3 - 1.076015721628102e183*u4 - 1.135655615419909e184*u5 - 2.652655707375446e184*u6 - 7.723437875763073e183*u7) + 
                      logu5*(-3.723056700913299e179 + 2.530025690483013e180*u - 1.247664113742667e181*u2 - 3.175874031680214e182*u3 + 6.723039694836914e181*u4 - 1.923745755727437e183*u5 - 4.871090456390997e183*u6 - 2.09622322599808e183*u7) + 
                      logu6*(-6.253102804277507e178 + 6.572991520874566e179*u - 3.017115838764146e179*u2 - 4.101325568370858e181*u3 + 1.107363537142588e182*u4 - 3.512442720863993e181*u5 - 1.907450956796042e182*u6 - 3.626516356932299e182*u7) + 
                      logu7*(-3.553959604720957e177 + 1.038451193952306e179*u + 3.203333214433537e179*u2 + 7.452754618772544e179*u3 + 3.107046432677259e181*u4 + 5.737348291773185e181*u5 + 1.335781704072417e182*u6 - 2.833829594641033e181*u7) + 
                      logu8*(1.001275185030998e177 + 5.308328459180423e177*u + 4.269607161360426e178*u2 + 1.271388520296124e180*u3 + 3.313948039810947e180*u4 + 9.501457914152601e180*u5 + 2.896784540030839e181*u6 + 2.822879836019659e180*u7) + 
                      logu9*(2.557624650506344e176 - 1.52447658776954e177*u - 6.93138636912931e177*u2 + 1.795591997415555e179*u3 - 2.643839464911647e179*u4 - 4.370974734175302e179*u5 + 5.421997288258592e179*u6 + 9.442105044367457e179*u7) + 
                      logu10*(1.657192797268684e175 - 3.483451703321782e176*u - 2.276849022974393e177*u2 - 2.517759223662924e177*u3 - 1.105995973972376e179*u4 - 2.822438605417451e179*u5 - 5.937108488899216e179*u6 + 6.120853659802775e178*u7) + 
                      logu11*(-2.199549299343665e174 - 1.855780694184339e175*u - 8.000389864913327e175*u2 - 3.104755703637823e177*u3 - 6.473230229858546e177*u4 - 1.39379462439791e178*u5 - 6.008013855414782e178*u6 - 7.93363313883039e177*u7) + 
                      logu12*(-4.120452435673366e173 + 2.813787820633901e174*u + 3.623854640207053e175*u2 - 1.835200955485576e176*u3 + 1.199753189782112e177*u4 + 3.926096750861723e177*u5 + 5.694402562682274e177*u6 - 1.343709744486429e177*u7) + 
                      logu13*(-5.302660932744038e171 + 4.081959154917619e173*u + 3.686878938792337e174*u2 + 2.6381520517586e175*u3 + 1.500303895332389e176*u4 + 3.658040985394538e176*u5 + 1.079913152859212e177*u6 + 7.582261410937595e174*u7) + 
                      logu14*(3.012031949980434e171 - 1.19191517520802e171*u - 2.496037034411261e173*u2 + 2.49466254053343e174*u3 - 7.779709515003898e174*u4 - 4.060047777625348e175*u5 - 4.003487759935744e175*u6 + 1.165313354094316e175*u7) + 
                      logu15*(1.466674388126535e170 - 2.315296459918937e171*u - 3.221337956476751e172*u2 - 1.273042027717591e173*u3 - 1.241775925973065e174*u4 - 2.96125686887704e174*u5 - 8.806523427884157e174*u6 + 6.0983247000273e173*u7) + 
                      logu16*(-7.853013852910665e168 - 6.713055347530957e169*u + 4.063014752969727e170*u2 - 1.661206026837197e172*u3 + 1.974987673765557e172*u4 + 2.171434214355e173*u5 + 9.346577358562122e172*u6 - 4.829902671134612e172*u7) + 
                      logu17*(-6.777957330395917e167 + 4.436929149406857e168*u + 1.13649688486896e170*u2 + 2.045131957139145e170*u3 + 5.145572568217222e171*u4 + 1.37577650016579e172*u5 + 4.275032505756639e172*u6 - 4.464178761368798e171*u7) + 
                      logu18*(-2.395762731314108e165 + 2.327796603512565e167*u + 1.777920771928429e168*u2 + 5.90350187032057e169*u3 + 4.857273276226472e169*u4 - 4.793042804570739e170*u5 + 5.581319425544714e170*u6 + 1.599927204313149e169*u7) + 
                      logu19*(9.591923547625109e164 - 1.316441609367546e163*u - 1.303498056799728e167*u2 + 9.003188094060316e167*u3 - 8.796189767806308e168*u4 - 3.562173918229893e169*u5 - 1.12966419355446e170*u6 + 9.323532095546503e168*u7) + 
                      logu20*(2.97242101502897e163 - 1.589298102464979e164*u - 4.994772576534919e165*u2 - 7.932447668923482e166*u3 - 2.655979043542925e167*u4 + 5.113683727713408e166*u5 - 3.872050210420418e168*u6 + 2.534713029370771e167*u7) + 
                      logu21*(2.778534106980166e161 - 2.24499559607378e162*u - 5.037625527001666e163*u2 - 3.178534212586987e165*u3 - 3.412456879494412e164*u4 + 3.048308166392644e166*u5 + 9.48461315740492e166*u6 + 1.826247049784983e165*u7) + 
                      logu22*u3*(-3.428360642270303e163 + 4.319474056046067e163*u + 4.380572523409785e164*u2 + 6.510635448543724e165*u3) + 
                      logu23*8.632989034111797e163*u6)/(fD2Int*fD2Int*Den_d1SFInt*Den_d1SFInt*Den_d1SFInt);

  double D2d1SF_Tail = (-7.634706065332545e28*u3by2*(-4.207922237965502e56 - 7.999065405344364e57*u - 5.309878379975745e58*u2 + 7.630824296487991e57*u3 + 1.347723391402615e60*u4 + 6.792525608293873e60*u5 + 3.135364343039252e60*u6 +
                        1.726273952191622e39*logu2*(-9.040024197009684e17 + 1.290525703883399e19*u + 1.535261923427384e19*u2)*u4 + 
                        8.309690613233738e19*logu*u2*(-5.710120252699959e36 + 1.073743136376743e38*u + 1.115798111024323e39*u2 + 9.589347763949992e39*u3 + 7.168762610998614e39*u4)))/(Den_d1SFTail*Den_d1SFTail*Den_d1SFTail);

  // D2(d1SF)   

  double D2d1SF = D2d1SF_Int*d1SFTail + 2.*Dd1SF_Int*Dd1SF_Tail + d1SFInt*D2d1SF_Tail;

  // adding the GSF correction 
  
  double D2f1SF = 1.*(3.641794632026264e7*u11by2 + 48.75*(-9.898639062741917e6 + 2.6012818800187605e6*logu)*u11by2 + 
                  42.*(1.0018072583960408e7 + 2.6340749138477915e6*logu + 198577.79085898417*logu2)*u5 + 
                  14.*(2.6340749138477915e6/u + (397155.58171796834*logu)/u)*u6 + 
                  (-2.236919332129823e6/u2 - (397155.58171796834*logu)/u2)*u7);

  D2d1SF = D2f1SF*d1SF + 2.*Dd1SF*Df1SF + f1SF*D2d1SF; 

  /* D potential and derivatives wrt u */

  double barD   = 1. + 6.*nu*u2*d1SF;
  double barD2  = barD*barD;
  double barD3  = barD2*barD;
  *D            = 1./barD;

  double DbarD  = nu*(12.*u*d1SF + 6.*u2*Dd1SF);
  *dD           = - DbarD/barD2;

  double D2barD = nu*(12.*d1SF + 24.*u*Dd1SF + 6.*u2*D2d1SF);
  *d2D          = 2.*(DbarD*DbarD)/barD3 - D2barD/barD2;
  
}

/** EOB Metric Q function at 3PN */

void eob_metric_Q3PN(double r, double prstar, double nu, double *Q, double *dQ_du, double *dQ_dprstar, double *ddQ_drdprstar, double *d2Q_dprstar2)
{
  const double z3 = 2.*nu*(4. - 3.*nu);
  double u  = 1./r;
  double u2 = u*u;
  double u3 = u2*u;
  double prstar2 = prstar*prstar;
  double prstar3 = prstar2*prstar;
  double prstar4 = prstar2*prstar2;

  *Q             = z3*u2*prstar4;
  *dQ_du         = 2.*z3*u*prstar4;
  *dQ_dprstar    = 4.*z3*u2*prstar3;
  *ddQ_drdprstar = -8.*z3*u3*prstar3;
  *d2Q_dprstar2  = 12.*z3*u2*prstar2;

}

/** EOB Metric Q function GSF-informed */

void eob_metric_QGSF(double r, double prstar, double nu, double *Q, double *dQ_du, double *dQ_dprstar, double *ddQ_drdprstar, double *d2Q_dprstar2)
{

  double u = 1./r;
  double u2     = u*u;
  double u3     = u*u2;
  double u4     = u3*u;
  double u5     = u4*u;
  double u6     = u3*u3;
  double u7     = u6*u;
  double sqrtu  = sqrt(u);
  double u3by2  = sqrt(u3);
  double u5by2  = sqrt(u5);
  double u7by2  = sqrt(u7);
  double u13by2 = sqrt(u7*u6);
  double u15by2 = sqrt(u7*u7*u);
  double logu   = log(u);    
  double logu2  = logu*logu;
  double logu3  = logu2*logu;
  double logu4  = logu3*logu;
  double logu5  = logu4*logu;
  double logu6  = logu5*logu;
  double logu7  = logu6*logu;
  double logu8  = logu7*logu;
  double logu9  = logu8*logu;
  double logu10 = logu9*logu;
  double logu11 = logu10*logu;
  double logu12 = logu11*logu;
  double logu13 = logu12*logu;
  double logu14 = logu13*logu;
  double logu15 = logu14*logu;
  double prstar2 = prstar*prstar;
  double prstar3 = prstar2*prstar;
  double prstar4 = prstar2*prstar2;

  /* Integer part */

  double den = 1.5436057887949119e47 + 4.0197014785265873e46*logu + 7.336891510046999e45*logu2 + 1.93993845760767e44*logu3 - 6.650397583926134e43*logu4;

  double cNUM1 = (8.671883671883672e-8*(2.9922854126728384e55 + 9.373903063776709e54*logu + 1.4902772115762003e54*logu2 - 1.1870358464650615e53*logu3 - 2.216517843338154e52*logu4 + 7.2559118997243416e50*logu5))/den;
  double cNUM2 = (5.865722180657246e-13*(4.1904196783613106e61 + 9.328728689754337e60*logu + 4.215054515446438e59*logu2 - 2.540526550671673e59*logu3 - 6.01282301840797e57*logu4 - 4.732256858711584e56*logu5))/den;
  double cNUM3 = (7.273961037521385e-18*(-1.8884724609902733e67 - 1.333202857385806e67*logu - 1.7612067328905822e66*logu2 - 6.936063514271904e63*logu3 + 1.1620685429316606e63*logu4 - 3.627093949696908e63*logu5 + 1.639400804564506e62*logu6))/den;

  double cDEN1 = (2.6015651015651016e-7*(7.844869763556561e54 + 2.5701117634388443e54*logu + 3.955460520146066e53*logu2 - 4.2244022634971293e52*logu3 - 6.470964973997371e51*logu4 + 2.418637299908114e50*logu5))/den;
  double cDEN2 = (2.6395749812957606e-11*(7.247289853081318e59 + 9.709926313193797e58*logu - 1.1089732280061834e58*logu2 - 5.8582461503798106e57*logu3 + 1.71864535381638e55*logu4 - 2.7905911010816658e54*logu5))/den;
  double cDEN3 = (2.182188311256416e-17*(-8.209334537517876e66 - 4.205835322920629e66*logu - 4.49329548985677e65*logu2 + 2.78600018346364e64*logu3 + 3.678533385659995e63*logu4 - 1.036914394270767e63*logu5 + 3.6014080126105997e61*logu6))/den;

  double q1SF_Int = (1. + cNUM1*u + cNUM2*u2 + cNUM3*u3)/(1. + cDEN1*u + cDEN2*u2 + cDEN3*u3);

  /* Tail part */

  double den2 = 1.310215755290071e32 + 2.244932235664865e31*logu;
  double fDenTail = 1.290137297437158e48 + 4.936261309852659e49*u + 9.019433892978062e50*u2 + 8.775970408136903e51*u3 - 
                    1.049332668650155e26*u*(-1.160818436398093e21 - 4.764125797365334e22*u - 6.444964435529441e23*u2) + 
                    7.87366594294329e25*logu2*u*(1.160818436398093e21 + 4.764125797365334e22*u + 6.704749518001747e23*u2) + 
                    2.70343066432167e8*(1.176703680646096e39 - 5.883218038414272e39*u - 2.392394015070348e41*u2 - 7.319235654333975e42*u3) + 
                    2.3417856e8*logu*(1.176703680646096e39 - 4.982078072992192e39*u - 2.022556319381654e41*u2 - 6.808831425479551e42*u3) +
                    5.906114128924534e47*logu3*u3;

  double NumTail = (8.146836362023952e-17*(1.290137297437158e48 + 4.936261309852659e49*u + 9.019433892978062e50*u2 + 8.775970408136903e51*u3 + 
       5.906114128924534e47*logu3*u3 + 2.377783562583731e49*u5by2 + 
       1.049332668650155e26*u*(1.160818436398093e21 + 4.764125797365334e22*u + 6.444964435529441e23*u2 + 2.139442835033594e22*u5by2) + 
       7.87366594294329e25*logu2*u*(1.160818436398093e21 + 4.764125797365334e22*u + 6.704749518001747e23*u2 + 2.139442835033594e22*u5by2) + 
       2.3417856e8*logu*(1.176703680646096e39 - 4.982078072992192e39*u - 2.022556319381654e41*u2 - 6.808831425479551e42*u3 + 2.168720085397228e40*u5by2 - 5.59378451258966e41*u7by2) + 3.971475085647266e50*u7by2 - 
       2.70343066432167e8*(-1.176703680646096e39 + 5.883218038414272e39*u + 2.392394015070348e41*u2 + 7.319235654333975e42*u3 - 2.168720085397228e40*u5by2 + 5.759868835678712e41*u7by2)))/den2;

  double DenTail = (8.146836362023952e-17*fDenTail)/den2;

  double q1SF_Tail = NumTail/DenTail;
  
  /* q1SF */

  double q1SF_tmp = 8.*u2*q1SF_Int*q1SF_Tail;

  /* First derivatives */

  /* of the integer part */

  double Dq1SF_IntNum = 1.3199825654293425e142 + 2.8348649025312635e143*u + 3.245039844757882e144*u2 + 6.130359198669476e145*u3 + 2.3896554647022278e146*u4 + 7.35572253603403e146*u5 + 
                        logu*(1.0312100703801561e142 + 2.9245942428764234e143*u + 2.2197128114101742e144*u2 + 4.566907139165532e145*u3 + 1.9372083987866826e146*u4 + 4.842434943648264e146*u5) + 
                        logu2*(4.567571229846957e141 + 1.4845624306397553e143*u + 8.314281796430755e143*u2 + 1.600989223669338e145*u3 + 8.766049967734823e145*u4 + 1.8875747846445416e146*u5) +
                        logu3*(1.2631524163419016e141 + 4.5604107553221665e142*u + 1.4050244771641895e143*u2 + 2.506339001280541e144*u3 + 2.313991157788905e145*u4 + 6.054211547951584e145*u5) + 
                        logu4*(2.2595940214460124e140 + 8.860563816919735e141*u - 1.1244818548823825e142*u2 - 1.0242776666076194e143*u3 + 3.7117653875046295e144*u4 + 1.6152863072672808e145*u5) + 
                        logu5*(2.2517090729794913e139 + 9.344747621535008e140*u - 1.239671831222866e142*u2 - 1.3059857504603616e143*u3 + 3.0640095535887064e143*u4 + 3.0212659043188435e144*u5) + 
                        logu6*(-6.685044684924607e136 - 7.526976063423245e137*u - 3.174622735004624e141*u2 - 2.4433871862143373e142*u3 + 1.55111906626354e142*u4 + 3.214026545306396e143*u5) + 
                        logu7*(-3.3650454658844735e137 - 1.5252004857854732e139*u - 3.952418091399067e140*u2 - 1.0175819766186294e141*u3 + 5.9533817738080455e141*u4 - 4.223473550249853e141*u5) + 
                        logu8*(-3.9387443570702907e136 - 1.868363853801939e138*u - 1.094058671243926e139*u2 + 3.413930741474212e140*u3 + 1.7600214332644416e141*u4 - 5.491098215895841e141*u5) + 
                        logu9*(-9.792861762019303e133 - 1.2880965107519904e136*u + 3.800371941149868e138*u2 + 6.176344638411114e139*u3 + 1.7983541858727996e140*u4 - 4.501379719748879e140*u5) + 
                        logu10*(3.224249860768159e134 + 1.35562117016215e136*u + 4.769704565480593e137*u2 + 2.2026102840779024e138*u3 - 7.158389714574527e138*u4 + 1.7626270002695224e139*u5) + 
                        logu11*(9.237690916560419e132 + 6.543506505561098e134*u + 4.272557355107218e135*u2 - 4.28039763506431e137*u3 - 2.3764806969959435e138*u4 + 3.844701193337184e138*u5) + 
                        logu12*(-1.0556059541541925e132 - 2.6587869852767965e133*u - 2.596139588141419e135*u2 - 3.490263636027969e136*u3 - 4.72134301181584e136*u4 - 3.5710330372423043e136*u5) + 
                        logu13*u*(-1.8037975677301894e132 - 8.520542891846365e133*u + 1.3075015126643744e135*u2 + 9.479548786624334e135*u3 - 1.3060715480187995e136*u4) + 
                        logu14*u2*(6.24821576991909e132 + 1.057532549359379e134*u + 8.476415750945398e133*u2 + 4.155518512476097e134*u3) + 
                        logu15*(-3.402898905009507e132 - 8.666184902472648e132*u)*u3;

  double fDIntDen = 1.5436057887949119e47 + 2.040893940319202e48*u + 1.9129764978392074e49*u2 - 1.7914313870965105e50*u3 + 
                    logu*(4.0197014785265873e46 + 6.686313070884439e47*u + 2.563007856653173e48*u2 - 9.177924680746751e49*u3) + 
                    logu2*(7.336891510046999e45 + 1.0290388049830552e47*u - 2.9272179875719206e47*u2 - 9.805216896986615e48*u3) + 
                    logu3*(1.93993845760767e44 - 1.0990057503686754e46*u - 1.546327997281475e47*u2 + 6.079577035512586e47*u3) + 
                    logu4*(-6.650397583926134e43 - 1.6834636649801683e45*u + 4.536493277653917e44*u2 + 8.027252556753731e46*u3) + 
                    logu5*u*(6.292242392784595e43 - 7.365974453441754e43*u - 2.2627424709511944e46*u2) + 
                    logu6*7.85895046918405e44*u3;

  double Dq1SF_IntDen = den*fDIntDen*fDIntDen;

  double Dq1SF_Int = Dq1SF_IntNum/Dq1SF_IntDen;

  /* of the tail part */

  double Dq1SF_TailNum = 2.189564415845942e6*u3by2*(3.170594381937833e27*den2*den2 - 
     2.277847717464087e39*(-5.292034219762144e19 + 8.291560259986376e18*logu)*(1.310215755290071e32 + 2.244932235664865e31*logu)*u + 
     4.276666692501831e35*(-1.71712036698414e24 + 4.41626925136341e22*logu)*(1.310215755290071e32 + 2.244932235664865e31*logu)*u + 
     5.511175469608051e47*(-5.564541370244304e43 + 2.99133561138563e43*logu - 3.846714457261093e41*logu2)*u + 
     3.80222731529389e52*(2.921698136444077e39 - 7.117254389988845e37*logu + 5.575657791466199e36*logu2)*u - 
     2.540540370142495e21*(-1.310215755290071e32 - 2.244932235664865e31*logu)*(2.921698136444077e39 - 7.117254389988845e37*logu + 5.575657791466199e36*logu2)*u + 
     8.592292380672e16*(1.310215755290071e32 + 2.244932235664865e31*logu)*(5.564541370244304e43 - 2.99133561138563e43*logu + 3.846714457261093e41*logu2)*u - 
     1.20784547883346e24*(4.051595191069018e36 - 6.417538657180621e35*logu)*(-1.310215755290071e32 - 2.244932235664865e31*logu)*u2 + 
     5.711360433362489e29*(-1.71712036698414e24 + 4.41626925136341e22*logu)*(2.921698136444077e39 - 7.117254389988845e37*logu + 5.575657791466199e36*logu2)*u2 - 
     4.409252928048242e28*(-5.292034219762144e19 + 8.291560259986376e18*logu)*(5.564541370244304e43 - 2.99133561138563e43*logu + 3.846714457261093e41*logu2)*u2 + 
     8.1962496e10*(2.921698136444077e39 - 7.117254389988845e37*logu + 5.575657791466199e36*logu2)*(5.564541370244304e43 - 2.99133561138563e43*logu + 3.846714457261093e41*logu2)*u2 + 
     1.54311173510304e48*(1.266032089826438e45 - 7.119398244223421e43*logu + 5.638398123889636e42*logu2)*u2 - 
     3.4368782063616e16*(-1.310215755290071e32 - 2.244932235664865e31*logu)*(1.266032089826438e45 - 7.119398244223421e43*logu + 5.638398123889636e42*logu2)*u2 - 
     2.33803874304e13*(4.051595191069018e36 - 6.417538657180621e35*logu)*(-5.564541370244304e43 + 2.99133561138563e43*logu - 3.846714457261093e41*logu2)*u3 + 
     2.419562257278566e19*(-1.310215755290071e32 - 2.244932235664865e31*logu)*(-6.808831425479551e42 + 4.508607272660297e41*logu + 7.566167623019631e39*logu2)*u3 + 
     2.317922253799711e25*(-1.71712036698414e24 + 4.41626925136341e22*logu)*(1.266032089826438e45 - 7.119398244223421e43*logu + 5.638398123889636e42*logu2)*u3 + 
     1.99584e6*(5.564541370244304e43 - 2.99133561138563e43*logu + 3.846714457261093e41*logu2)*(1.266032089826438e45 - 7.119398244223421e43*logu + 5.638398123889636e42*logu2)*u3 + 
     2.319492146318903e42*(6.864894914709503e51 - 1.594482338501548e51*logu + 5.279095793585579e49*logu2 + 5.906114128924534e47*logu3)*u3 + 
     5.16606272e10*(-1.310215755290071e32 - 2.244932235664865e31*logu)*(6.864894914709503e51 - 1.594482338501548e51*logu + 5.279095793585579e49*logu2 + 5.906114128924534e47*logu3)*u3 - 
     4.6835712e8*(-6.808831425479551e42 + 4.508607272660297e41*logu + 7.566167623019631e39*logu2)*(5.564541370244304e43 - 2.99133561138563e43*logu + 3.846714457261093e41*logu2)*u4 + 
     3.484130371873063e19*(-1.71712036698414e24 + 4.41626925136341e22*logu)*(6.864894914709503e51 - 1.594482338501548e51*logu + 5.279095793585579e49*logu2 + 5.906114128924534e47*logu3)*u4 - 
     1.*(-5.564541370244304e43 + 2.99133561138563e43*logu - 3.846714457261093e41*logu2)*(6.864894914709503e51 - 1.594482338501548e51*logu + 5.279095793585579e49*logu2 + 5.906114128924534e47*logu3)*u4);

  double Dq1SF_TailDen = fDenTail*fDenTail;

  double Dq1SF_Tail = Dq1SF_TailNum/Dq1SF_TailDen;

  /* First derivative of q1SF wrt u */

  double Dq1SF_u_tmp = 16.*u*q1SF_Int*q1SF_Tail + 8.*u2*(Dq1SF_Int*q1SF_Tail + q1SF_Int*Dq1SF_Tail);

  /* GSF-informed correction and improved functions */

  double f1SF = 1.*(1. + (6.473498944710975e9 - 1.6413472171574914e9*logu)*u15by2 + (-6.555554279481237e9 - 1.7756296887290301e9*logu - 1.4618867042954516e8*logu2)*u7);
  double Df1SF = (-4.766450964509769e10 + logu*(-1.2721785161962301e10 - 1.0233206930068161e9*logu - 1.2310104128681187e10*sqrtu) + 4.690989486817482e10*sqrtu)*u6;

  double q1SF    = q1SF_tmp*f1SF;
  double Dq1SF_u = Dq1SF_u_tmp*f1SF + q1SF_tmp*Df1SF;

  /* complete Q function and derivatives */

  *Q             = nu*q1SF*prstar4;
  *dQ_du         = nu*Dq1SF_u*prstar4;
  *dQ_dprstar    = 4.*nu*q1SF*prstar3;
  *ddQ_drdprstar = -4.*u2*nu*Dq1SF_u*prstar3;
  *d2Q_dprstar2  = 12.*nu*q1SF*prstar2;

}

/** Tidal potential, three version implemented: 
    1. TEOB NNLO, Bernuzzi+ 1205.3403
    2. TEOBResum, Bini&Damour, 1409.6933, Bernuzzi+ 1412.4553 
    3. TEOBResum3, Akcay+ 1812.02744
    f-mode resonance model from Hinderer&Steinhoff can be added to all three.
*/

/* Macro for the bar_alpha coefs */
#define bar_alph(cA,kapA,cB,kapB,kap)( ((cA)*(kapA)+(cB)*(kapB))/(kap) )

void eob_metric_Atidal(double r, Dynamics *dyn, double *AT, double *dAT, double *d2AT)
{
  
  double A, dA_u, d2A_u, dA, d2A;

  const double elsix = 1.833333333333333333333;  // 11/6
  const double eightthird = 2.6666666666666666667; // 8/3

  const double nu    = EOBPars->nu;
  const double rLR   = EOBPars->rLR_tidal;
  const double XA    = EOBPars->X1;
  const double XB    = EOBPars->X2;

  double kapA2 = EOBPars->kapA2; 
  double kapA3 = EOBPars->kapA3;
  double kapA4 = EOBPars->kapA4;
  
  double kapB2 = EOBPars->kapB2;  
  double kapB3 = EOBPars->kapB3;
  double kapB4 = EOBPars->kapB4;
  
  double kapT2 = EOBPars->kapT2;
  double kapT3 = EOBPars->kapT3;
  double kapT4 = EOBPars->kapT4;
  const double kapT5 = EOBPars->kapT5;
  const double kapT6 = EOBPars->kapT6;
  const double kapT7 = EOBPars->kapT7;
  const double kapT8 = EOBPars->kapT8;
  
  const double kapA2j = EOBPars->kapA2j;
  const double kapB2j = EOBPars->kapB2j;
  const double kapT2j = EOBPars->kapT2j;  
  
  const double p = EOBPars->pGSF_tidal;
  
  /* Definition of the conservative tidal coefficients \bar{\alpha}_n^{(\ell)}, 
     Eq.(37) of Damour&Nagar, PRD 81, 084016 (2010) */
  double bar_alph2_1 = EOBPars->bar_alph2_1;
  double bar_alph2_2 = EOBPars->bar_alph2_2;
  double bar_alph3_1 = EOBPars->bar_alph3_1;
  double bar_alph3_2 = EOBPars->bar_alph3_2;
  double bar_alph2j_1 = EOBPars->bar_alph2j_1;

  double kapA2_u = kapA2;
  double kapA3_u = kapA3;
  double kapA4_u = kapA4;
  double kapB2_u = kapB2;
  double kapB3_u = kapB3;
  double kapB4_u = kapB4;
  double kapT2_u=0, kapT3_u=0, kapT4_u=0;
  double bar_alph2_1_u=0, bar_alph2_2_u=0, bar_alph3_1_u=0, bar_alph3_2_u=0;
  
  if (EOBPars->use_tidal_fmode_model) {
    
    /* Dress tidal coupling constants for ell=2,3,4, 
       and recompute coefficients as needed */

     kapA2 *= dyn->dress_tides_fmode_A[2]; 
     kapA3 *= dyn->dress_tides_fmode_A[3];
     kapA4 *= dyn->dress_tides_fmode_A[4];
     
     kapB2 *= dyn->dress_tides_fmode_B[2];
     kapB3 *= dyn->dress_tides_fmode_B[3];
     kapB4 *= dyn->dress_tides_fmode_B[4];

     kapT2 = kapA2 + kapB2;
     kapT3 = kapA3 + kapB3;
     kapT4 = kapA4 + kapB4;

     kapA2_u *= dyn->dress_tides_fmode_A_u[2];
     kapA3_u *= dyn->dress_tides_fmode_A_u[3];
     kapA4_u *= dyn->dress_tides_fmode_A_u[4];

     kapB2_u *= dyn->dress_tides_fmode_B_u[2];
     kapB3_u *= dyn->dress_tides_fmode_B_u[3];
     kapB4_u *= dyn->dress_tides_fmode_B_u[4];

     kapT2_u = kapA2_u + kapB2_u;
     kapT3_u = kapA3_u + kapB3_u;
     kapT4_u = kapA4_u + kapB4_u;      
    
    /* Tidal coefficients cons dynamics
       \bar{\alpha}_n^{(\ell)}, Eq.(37) of Damour&Nagar, PRD 81, 084016 (2010) 
       The structure is always
       bar_alpha = (cA * kapA + cB * kapB ) / kap
       hence
       bar_alpha' =  - bar_alpha kap'/kap + (cA * kapA' + cB * kapB' ) / kap
    */
    
    //const double bar_alph2_1 = (5./2.*XA*kapA2 + 5./2.*XB*kapB2)/kapT2;
    //const double bar_alph2_2 = ((3.+XA/8.+ 337./28.*XA*XA)*kapA2 + (3.+XB/8.+ 337./28.*XB*XB)*kapB2)/kapT2;
    //const double bar_alph3_1 = ((-2.+15./2.*XA)*kapA3 + (-2.+15./2.*XB)*kapB3)/kapT3;
    //const double bar_alph3_2 = ((8./3.-311./24.*XA+110./3.*XA*XA)*kapA3 + (8./3.-311./24.*XB+110./3.*XB*XB)*kapB3)/kapT3;

    const double cA21 = 2.5*XA;
    const double cB21 = 2.5*XB;
    const double cA22 = (3.+XA/8.+ 337./28.*XA*XA);
    const double cB22 = (3.+XB/8.+ 337./28.*XB*XB);
    const double cA31 = (-2.+15./2.*XA);
    const double cB31 = (-2.+15./2.*XB);
    const double cA32 = (8./3.-311./24.*XA+110./3.*XA*XA);
    const double cB32 = (8./3.-311./24.*XB+110./3.*XB*XB);
    
    bar_alph2_1 = bar_alph( cA21, kapA2, cB21, kapB2, kapT2 );
    bar_alph2_2 = bar_alph( cA22, kapA2, cB22, kapB2, kapT2 );
    bar_alph3_1 = bar_alph( cA31, kapA3, cB31, kapB3, kapT3 );
    bar_alph3_2 = bar_alph( cA32, kapA3, cB32, kapB3, kapT3 );
    
    bar_alph2_1_u =
      - bar_alph2_1 * kapT2_u/kapT2
      + bar_alph( cA21, kapA2_u, cB21, kapB2_u, kapT2 );      
    bar_alph2_2_u =
      - bar_alph2_2 * kapT2_u/kapT2
      + bar_alph( cA22, kapA2_u, cB22, kapB2_u, kapT2 );
    bar_alph3_1_u =
      - bar_alph3_1 * kapT3_u/kapT3
      + bar_alph( cA31, kapA3_u, cB31, kapB3_u, kapT3 );
    bar_alph3_2_u =
      - bar_alph3_2 * kapT3_u/kapT3
      + bar_alph( cA32, kapA3_u, cB32, kapB3_u, kapT3 );    
  } 


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
  double u12  = u6*u6;
  double u14  = u7*u7;
  double u16  = u10*u6;
  double u18  = u12*u6;
  double u8   = u5*u3;
  double u9   = u8*u;
  double u11  = u10*u;
  double u13  = u12*u;
  double u15  = u14*u;
  double u17  = u16*u;
  double logu = log(u);
  double oom3u  = 1./(1.-rLR*u);

  if (EOBPars->use_tidal==TIDES_NNLO) {

    A  = - kapT8*u18;
    A -= kapT7*u16;
    A -= kapT6*u14;
    A -= kapT5*u12;
    A -= kapT4*u10;
    A -= kapT3*u8*(1. + bar_alph3_1*u + bar_alph3_2*u2);
    A -= kapT2*u6*(1. + bar_alph2_1*u + bar_alph2_2*u2);
    
    dA_u  = -18.*kapT8*u17;
    dA_u -= 16.*kapT7*u15;
    dA_u -= 14.*kapT6*u13;
    dA_u -= 12.*kapT5*u11;
    dA_u -= 10.*kapT4*u9;
    dA_u -= kapT3*u8*(bar_alph3_1 + 2.*bar_alph3_2*u);
    dA_u -= 8.*kapT3*u7*(1. + bar_alph3_1*u + bar_alph3_2*u2);
    dA_u -= kapT2*u6*(bar_alph2_1 + 2.*bar_alph2_2*u);
    dA_u -= 6.*kapT2*u5*(1. + bar_alph2_1*u + bar_alph2_2*u2);
    
    if (d2AT != NULL) {
      d2A_u  = -306.*kapT8*u16;
      d2A_u -= 240.*kapT7*u14;
      d2A_u -= 182.*kapT6*u12;
      d2A_u -= 132.*kapT5*u10;
      d2A_u -= 90.*kapT4*u8;
      d2A_u -= kapT3*(2.*bar_alph3_2*u8 + 16.*u7*(bar_alph3_1 + 2*bar_alph3_2*u) + 56.*u6*(1 + bar_alph3_1*u + bar_alph3_2*u2));
      d2A_u -= kapT2*(2*bar_alph2_2*u6 + 12.*u5*(bar_alph2_1 + 2*bar_alph2_2*u) + 30.*u4*(1 + bar_alph2_1*u + bar_alph2_2*u2));
    }

    if (EOBPars->use_tidal_fmode_model) {
      /* Adding missing derivative terms from ell=4,3,2 */
      dA_u -= kapT4_u*u10;
      
      dA_u -= kapT3_u*u8*(1. + bar_alph3_1*u + bar_alph3_2*u2);
      dA_u -= kapT3*u8*(bar_alph3_1_u*u + bar_alph3_2_u*u2);      

      dA_u -= kapT2_u*u6*(1.+ bar_alph2_1*u + bar_alph2_2*u2);
      dA_u -= kapT2*u6*(bar_alph2_1_u*u + bar_alph2_2_u*u2);
    }    

  } else if (EOBPars->use_tidal==TIDES_TEOBRESUM) { 

    const double c1  =  8.533515908;  	// OLD value 8.53353;
    const double c2  =  3.043093411;	// OLD value 3.04309;
    const double n1  =  0.8400636422; 	// OLD value 0.840058;
    const double d2  = 17.7324036;	// OLD value 17.73239

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
    
    double f0     = 1 + 3*u2*oom3u;
    double f1     = A1SF *pow(oom3u,7./2.);
    double f2     = A2SF *pow(oom3u,p); 
    
    double df0    = 3*u*(2.-rLR*u)*(oom3u*oom3u);
    double df1    = 0.5*(7*rLR*A1SF + 2*(1.-rLR*u)*dA1SF)*pow(oom3u,9./2.);
    double df2    = (rLR*p*A2SF + (1.-rLR*u)*dA2SF)*pow(oom3u,p+1);

    /** Gravito-electric tides for el = 2, 3, 4 */
    double AT2    = - kapA2*u6*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*u6*( f0 + XB*f1 + XB*XB*f2 );
    double AT3    = - kapT3*u8*(1. + bar_alph3_1*u + bar_alph3_2*u2);
    double AT4    = - kapT4*u10;
    double AT5    = - kapT5*u12;
    double AT6    = - kapT6*u14;
    double AT7    = - kapT7*u16;
    double AT8    = - kapT8*u18;

    double dAT2  = - kapA2*6.*u5*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*6.*u5*( f0 + XB*f1 + XB*XB*f2 ) - kapA2*u6*( df0 + XA*df1 + XA*XA*df2 ) - kapB2*u6*( df0 + XB*df1 + XB*XB*df2 );
    double dAT3  = - kapT3*(8.*u7 + 9*bar_alph3_1*u8 + 10*bar_alph3_2*u9);
    double dAT4  = - kapT4*10.*u9;
    double dAT5  = - kapT5*12.*u11;
    double dAT6  = - kapT6*14.*u13;
    double dAT7  = - kapT7*16.*u15;
    double dAT8  = - kapT8*18.*u17;

    A     = AT2 + AT3 + AT4 + AT5 + AT6 + AT7 + AT8;
    dA_u  = dAT2 + dAT3  + dAT4 + dAT5 + dAT6 + dAT7 + dAT8;

    if (d2AT != NULL) {
      double d2f23  = 2*d2*(-1 + 3*d2*u2 + n1*u*(-3+d2*u2))*(Den*Den*Den);
      double d2A1SF = d2Acub*f23 + 2*dAcub*df23 + Acub*d2f23;
      double d2A2SF = 674./28.;
      double d2f0   = 6*(oom3u*oom3u*oom3u);
      double d2f1   = 0.25*(63*(rLR*rLR)*A1SF + 4*(-1+rLR*u)*(-7*rLR*dA1SF + (-1+rLR*u)*d2A1SF))*pow(oom3u,11./2.);
      double d2f2   = (  rLR*p*(1+p)*rLR*A2SF +(-1+rLR*u)*( -2.*p*rLR*dA2SF +(-1.+rLR*u)*d2A2SF )  )*pow(oom3u,p+2);
      
      double d2AT2  = - kapA2*30.*u4*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*30.*u4*( f0 + XB*f1 + XB*XB*f2 ) - 2.*kapA2*6.*u5*( df0 + XA*df1 + XA*XA*df2 ) - 2.*kapB2*6.*u5*( df0 + XB*df1 + XB*XB*df2 ) - kapA2*u6*( d2f0 + XA*d2f1 + XA*XA*d2f2 ) - kapB2*u6*( d2f0 + XB*d2f1 + XB*XB*d2f2 );
      double d2AT3  = - kapT3*(56.*u6 + 72.*bar_alph3_1*u7 + 90.*bar_alph3_2*u8);
      double d2AT4  = - kapT4*90.*u8;
      double d2AT5  = - kapT5*132.*u10; 
      double d2AT6  = - kapT6*182.*u12;
      double d2AT7  = - kapT7*240.*u14;
      double d2AT8  = - kapT8*306.*u16;

      d2A_u = d2AT2 + d2AT3 + d2AT4 + d2AT5 + d2AT6 + d2AT7 + d2AT8;
    }

    if (EOBPars->use_tidal_fmode_model) {
      /* Adding missing derivative terms from ell=4,3,2 */
      dA_u -= kapT4_u*u10;

      dA_u -= kapT3_u*u8*(1. + bar_alph3_1*u + bar_alph3_2*u2);
      dA_u -= kapT3*u8*(bar_alph3_1_u*u + bar_alph3_2_u*u2);
      
      dA_u -= kapA2_u*u6*( f0 + XA*f1 + XA*XA*f2 ) + kapB2_u*u6*( f0 + XB*f1 + XB*XB*f2 );      
    }
    
  } else if (EOBPars->use_tidal==TIDES_TEOBRESUM3) { 

    const double c1  =  8.533515908;  
    const double c2  = 3.043093411;
    const double n1  =  0.8400636422; 
    const double d2  =  17.7324036;	

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
    
    double f0     = 1 + 3*u2*oom3u;
    double f1     = A1SF *pow(oom3u,7./2.);
    double f2     = A2SF *pow(oom3u,p); 
    
    double df0    = 3*u*(2.-rLR*u)*(oom3u*oom3u);
    double df1    = 0.5*(7*rLR*A1SF + 2*(1.-rLR*u)*dA1SF)*pow(oom3u,9./2.);
    double df2    = (rLR*p*A2SF + (1.-rLR*u)*dA2SF)*pow(oom3u,p+1);

    /** Gravito-electric tides for el = 2, 4; el = 3 added below as a GSF series */
    double AT2    = - kapA2*u6*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*u6*( f0 + XB*f1 + XB*XB*f2 );
    double AT4    = - kapT4*u10;
    double AT5    = - kapT5*u12;
    double AT6    = - kapT6*u14;
    double AT7    = - kapT7*u16;
    double AT8    = - kapT8*u18;

    double dAT2  = - kapA2*6.*u5*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*6.*u5*( f0 + XB*f1 + XB*XB*f2 ) - kapA2*u6*( df0 + XA*df1 + XA*XA*df2 ) - kapB2*u6*( df0 + XB*df1 + XB*XB*df2 );
    double dAT4  = - kapT4*10.*u9;
    double dAT5  = - kapT5*12.*u11;
    double dAT6  = - kapT6*14.*u13;
    double dAT7  = - kapT7*16.*u15;
    double dAT8  = - kapT8*18.*u17;

    /** el = 3+, i.e.,  even parity tidal potential **/

    /* 1GSF fitting parameters */
    const double C1 = -3.6820949997216643;
    const double C2 = 5.171003322924513;
    const double C3 = -7.639164165720986;
    const double C4 = -8.63278143009751;
    const double C5 = 12.319646912775516;
    const double C6 = 16.36009385150114;

    /* 0SF -- el = 3+, i.e.,  even parity terms */
    double A3hat_Sch    = (1.0 - 2.0*u)*( 1.0 + eightthird*u2*oom3u );
    double dA3hat_Sch   = (1.0 - 2.0*u)*( eightthird*rLR*u2*oom3u*oom3u + 2.0*eightthird*u*oom3u ) - 2.0*( 1.0 + eightthird*u2*oom3u );
    double d2A3hat_Sch  = (1.0 - 2.0*u)*( 2.0*eightthird*rLR*rLR*u2*oom3u*oom3u*oom3u + 4.0*eightthird*rLR*u*oom3u*oom3u + 2.0*eightthird*oom3u ) - 4.0*( eightthird*rLR*u2*oom3u*oom3u + 2.0*eightthird*u*oom3u );

    /* 1SF -- el = 3+, i.e.,  even parity terms */
    double Denom3    = 1./(1. + C5*u2);
    double A3tilde   = 7.5*u*( 1 + C1*u + C2*u2 + C3*u3 )*( 1 + C4*u + C6*u2 )*Denom3;
    double dA3tilde  = 7.5*( 1 + 3*C2*u2 + 3*C6*u2 + 4*C3*u3 + 5*C2*C6*u4 + 6*C3*C6*u5 + C1*u*(2 + 3*C4*u + 4*C6*u2) + C4*u*(2 + 4*C2*u2 + 5*C3*u3) )*Denom3 + ( -15.*C5*u2*(1. + C4*u + C6*u2)*(1. + C1*u + C2*u2 + C3*u3) )*Denom3*Denom3;
    double d2A3tilde = 15.*( C1*(1 + 3*C4*u - 3*C5*pow(u,2) + 6*C6*pow(u,2) - C4*C5*pow(u,3) + 3*C5*C6*pow(u,4) + pow(C5,2)*C6*pow(u,6)) + C4*(1 - 3*C5*pow(u,2) + 10*C3*pow(u,3) +  9*C3*C5*pow(u,5) + 3*C3*pow(C5,2)*pow(u,7) +  C2*pow(u,2)*(6 + 3*C5*pow(u,2) + pow(C5,2)*pow(u,4))) + u*(3*(C6 + 2*C3*u + 5*C3*C6*pow(u,3)) + C5*(-3 - C6*pow(u,2) + 3*C3*pow(u,3) +17*C3*C6*pow(u,5)) + pow(C5,2)*(pow(u,2) + C3*pow(u,5) + 6*C3*C6*pow(u,7)) + C2*(3 + 10*C6*pow(u,2) + 3*pow(C5,2)*C6*pow(u,6) + C5*pow(u,2)*(-1 + 9*C6*pow(u,2)))) )*Denom3*Denom3*Denom3;
    double A3hat1GSFfit = A3tilde*pow(oom3u, 3.5);
    double dA3hat1GSFfit = 3.5*rLR*A3tilde*pow(oom3u, 4.5) + dA3tilde*pow(oom3u, 3.5);
    double d2A3hat1GSFfit = 15.75*rLR*rLR*A3tilde*pow(oom3u, 5.5) + 7.0*rLR*dA3tilde*pow(oom3u, 4.5) + d2A3tilde*pow(oom3u, 3.5);
    
    /* 2SF -- el = 3+, i.e.,  even parity terms */
    double A3hat2GSF     =  36.666666666666666667*u2*pow(oom3u,p);
    double dA3hat2GSF    =  36.666666666666666667*u*( 2. + (p - 2.)*rLR*u ) * pow(oom3u, p+1);
    double d2A3hat2GSF   =  36.666666666666666667*( 2. + 4.*(p - 1.)*rLR*u + (2. - 3.*p + 1.*p*p)*rLR*rLR*u2 ) * pow(oom3u, p+2);

    /* Hatted el = 3+ potential as a GSF series */
    double A3hatA   = A3hat_Sch + XA*A3hat1GSFfit + XA*XA*A3hat2GSF;
    double dA3hatA  = dA3hat_Sch + XA*dA3hat1GSFfit + XA*XA*dA3hat2GSF;
    double A3hatB   = A3hat_Sch + XB*A3hat1GSFfit + XB*XB*A3hat2GSF;
    double dA3hatB  = dA3hat_Sch + XB*dA3hat1GSFfit + XB*XB*dA3hat2GSF;
    
    /* Total el = 3+ tidal potential */
    double AT3      = - kapA3*u8*A3hatA - kapB3*u8*A3hatB;
    double dAT3     = - kapA3*u7*( 8.*A3hatA + u*dA3hatA ) - kapB3*u7*( 8.*A3hatB + u*dA3hatB );

    A     = AT2   + AT3   + AT4 + AT5 + AT6 + AT7 + AT8; 
    dA_u  = dAT2  + dAT3  + dAT4 + dAT5 + dAT6 + dAT7 + dAT8;;

    if (d2AT != NULL) {
      double d2f23  = 2*d2*(-1 + 3*d2*u2 + n1*u*(-3+d2*u2))*(Den*Den*Den);
      double d2A1SF = d2Acub*f23 + 2*dAcub*df23 + Acub*d2f23;
      double d2A2SF = 674./28.;
      double d2f0   = 6*(oom3u*oom3u*oom3u);
      double d2f1   = 0.25*(63*(rLR*rLR)*A1SF + 4*(-1+rLR*u)*(-7*rLR*dA1SF + (-1+rLR*u)*d2A1SF))*pow(oom3u,11./2.);
      double d2f2   = (  rLR*p*(1+p)*rLR*A2SF +(-1+rLR*u)*( -2.*p*rLR*dA2SF +(-1.+rLR*u)*d2A2SF )  )*pow(oom3u,p+2);
      
      double d2AT2  = - kapA2*30*u4*( f0 + XA*f1 + XA*XA*f2 ) - kapB2*30*u4*( f0 + XB*f1 + XB*XB*f2 ) - 2*kapA2*6*u5*( df0 + XA*df1 + XA*XA*df2 ) - 2*kapB2*6*u5*( df0 + XB*df1 + XB*XB*df2 ) - kapA2*u6*( d2f0 + XA*d2f1 + XA*XA*d2f2 ) - kapB2*u6*( d2f0 + XB*d2f1 + XB*XB*d2f2 );
      double d2AT4  = - kapT4*90*u8;
      double d2AT5  = - kapT5*132.*u10; 
      double d2AT6  = - kapT6*182.*u12;
      double d2AT7  = - kapT7*240.*u14;
      double d2AT8  = - kapT8*306.*u16;
      
      double d2A3hatA = d2A3hat_Sch + XA*d2A3hat1GSFfit + XA*XA*d2A3hat2GSF;
      double d2A3hatB = d2A3hat_Sch + XB*d2A3hat1GSFfit + XB*XB*d2A3hat2GSF;
      double d2AT3 = -1.*kapA3 * ( 56.*u6*A3hatA + 16.*u7*dA3hatA + 1.*u8*d2A3hatA ) - 1.*kapB3 * ( 56.*u6*A3hatB + 16.*u7*dA3hatB + 1.*u8*d2A3hatB );
      
      d2A_u += d2AT2  + d2AT3  + d2AT4 + d2AT5 + d2AT6 + d2AT7 + d2AT8;
    }

    if (EOBPars->use_tidal_fmode_model) {
      /* Adding missing derivative terms from ell=4,3,2 */
      dA_u -= kapT4_u*u10;
      
      dA_u -= kapA3_u*u8*A3hatA + kapB3_u*u8*A3hatB;
      
      dA_u -= kapA2_u*u6*( f0 + XA*f1 + XA*XA*f2 ) + kapB2_u*u6*( f0 + XB*f1 + XB*XB*f2 );      
    }
    
  } // EOBPars->use_tidal


#if(USEGRAVITOMAGNETICTERMS)

  if (EOBPars->use_tidal_gravitomagnetic==TIDES_GM_PN) {

    /* PN series for the (2-) tidal potential */
    A    +=-kapT2j*u7*(1. +  bar_alph2j_1*u);
    dA_u += -kapT2j*u7*bar_alph2j_1 - 7.*kapT2j*u6*(1. +  bar_alph2j_1*u);
    
    if (d2AT != NULL) {
      d2A_u += - 14.*kapT2j*u5*(3. + 4.*bar_alph2j_1*u);
    }
    
  } else if (EOBPars->use_tidal_gravitomagnetic==TIDES_GM_GSF) {

    /** GSF series for the (2-) tidal potential */
    const double a1j =  0.728591192;
    const double a2j =  3.100367557;	
    const double n1j = -15.04421708;
    const double d2j =  12.55229698;
    // Schwarzschild gravito-magnetic term
    double Ahat_Schj     =  (1.-2.*u)*oom3u;	
    double dAhat_Schj    =  (rLR-2.)*oom3u*oom3u;
    double d2Ahat_Schj   =  2.*rLR*(rLR-2.)*pow(oom3u, 3.);

    /* 1SF -- el = 2 gravitomagnetic terms */
    double Denomj = 1./(1. + d2j*u2);
    double Ahat1GSFfitj = elsix*u*(1. - a1j*u)*(1. - a2j*u)*(1. + n1j*u)*Denomj*pow(oom3u, 3.5);
    double dAhat1GSFfitj = 0.5*elsix * Denomj * Denomj * (2 + 4*n1j*u + 5*rLR*u - 2*d2j*pow(u,2) + 3*n1j*rLR*pow(u,2) + 9*d2j*rLR*pow(u,3) + 7*d2j*n1j*rLR*pow(u,4) -    a2j*u*(4 + rLR*u*(3 + 7*d2j*pow(u,2)) + n1j*u*(6 + 2*d2j*pow(u,2) + rLR*(u + 5*d2j*pow(u,3)))) + a1j*u*(-4 - 3*rLR*u - 7*d2j*rLR*pow(u,3) - n1j*u*(6 + rLR*u + 2*d2j*pow(u,2) + 5*d2j*rLR*pow(u,3)) + a2j*u*(6 + rLR*u + 2*d2j*pow(u,2) + 5*d2j*rLR*pow(u,3) + n1j*u*(8 - rLR*u + 4*d2j*pow(u,2) + 3*d2j*rLR*pow(u,3)))) ) * pow(oom3u, 4.5);
    double d2Ahat1GSFfitj = 0.25*elsix * Denomj * Denomj * Denomj * ( 8*(1 + n1j*u)*pow(-1 + rLR*u,2)*pow(1 + d2j*pow(u,2),2)*(-a2j + a1j*(-1 + 3*a2j*u)) +    4*(1 - rLR*u)*(1 + d2j*pow(u,2))*(1 - 2*a2j*u + a1j*u*(-2 + 3*a2j*u))*(-4*d2j*u + rLR*(7 + 11*d2j*pow(u,2)) + n1j*(2 - 2*d2j*pow(u,2) + rLR*u*(5 + 9*d2j*pow(u,2)))) + u*(-1 + a1j*u)*(-1 + a2j*u)*(7*rLR*(9*rLR + n1j*(4 + 5*rLR*u)) + 2*d2j*(-4 - 20*rLR*u + 87*pow(rLR,2)*pow(u,2) +3*n1j*u*(-4 + 8*rLR*u + 17*pow(rLR,2)*pow(u,2))) + pow(d2j,2)*pow(u,2)*(24 - 104*rLR*u + 143*pow(rLR,2)*pow(u,2) + n1j*u*(8 - 44*rLR*u + 99*pow(rLR,2)*pow(u,2)))) ) * pow(oom3u, 5.5); 
    
    /* 2SF -- el = 2 gravitomagnetic terms */
    double Ahat2GSFj    =  u*pow(oom3u,p);
    double dAhat2GSFj    =  ( 1.+ (p-1.)*rLR*u ) * pow(oom3u, p+1);
    double d2Ahat2GSFj   =  p*rLR * ( 2.+ (p-1.)*rLR*u ) * pow(oom3u, p+2);

    /* Total el = 2 gravitomagnetic potential as a GSF series */
    double AhatjA   = Ahat_Schj + XA*Ahat1GSFfitj + XA*XA*Ahat2GSFj;
    double dAhatjA  = dAhat_Schj + XA*dAhat1GSFfitj + XA*XA*dAhat2GSFj;
    double AhatjB   = Ahat_Schj + XB*Ahat1GSFfitj + XB*XB*Ahat2GSFj;
    double dAhatjB  = dAhat_Schj + XB*dAhat1GSFfitj + XB*XB*dAhat2GSFj;
    
    /* el = 2 gravitomagnetic total contribution */
    double ATj_2      = - kapA2j*u7*( AhatjA ) - kapB2j*u7*( AhatjB );
    double dATj_2     = - kapA2j * ( 7.*u6*AhatjA + u7*dAhatjA ) - kapB2j * ( 7.*u6*AhatjB + u7*dAhatjB );
    
    A    += ATj_2;
    dA_u += dATj_2;

    if (d2AT != NULL) {
      double d2AhatjA = d2Ahat_Schj + XA*d2Ahat1GSFfitj + XA*XA*d2Ahat2GSFj;
      double d2AhatjB = d2Ahat_Schj + XB*d2Ahat1GSFfitj + XB*XB*d2Ahat2GSFj;
      double d2ATj_2    = - kapA2j * ( 42.*u5*AhatjA + 14.*u6*dAhatjA + u7*d2AhatjA ) - kapB2j * ( 42.*u5*AhatjB + 14.*u6*dAhatjB + u7*d2AhatjB );
      d2A_u += d2ATj_2;
    }
    
  } // EOBPars->use_tidal_gravitomagnetic

#endif
    
  *AT   = A;
  *dAT  = dA_u;
  if (d2AT != NULL) *d2AT = d2A_u;

}


/** Tidal B potential
    Vines, Flanagan 1PN term */
void eob_metric_Btidal(double r, Dynamics *dyn, double *BT, double *dBT, double *d2BT)
{
  const double nu = EOBPars->nu;

  const double u  = 1./r;
  const double u2 = u*u;
  const double u4 = u2*u2;
  const double u5 = u4*u;
  const double u6 = u5*u;
  
  double kapA2 = EOBPars->kapA2; 
  double kapB2 = EOBPars->kapB2;
  double kapT2 = EOBPars->kapT2;
  double kapA2_u = kapA2;
  double kapB2_u = kapB2;
  double kapT2_u = 0;
  
  if (EOBPars->use_tidal_fmode_model) {
    kapA2 *= dyn->dress_tides_fmode_A[2]; 
    kapB2 *= dyn->dress_tides_fmode_B[2];
    kapT2 = kapA2 + kapB2;
    
    kapA2_u *= dyn->dress_tides_fmode_A_u[2]; 
    kapB2_u *= dyn->dress_tides_fmode_B_u[2]; 
    kapT2_u = kapA2_u + kapB2_u;    
  }

  const double cnu = (8. - 15.*nu);
  const double tmp = kapT2*cnu; 
  const double tmp_u = kapT2_u*cnu;

  double B, dB_u, d2B_u;    
  B     = tmp*u6; 
  dB_u  = 6*tmp*u5 + tmp_u*u6;
  d2B_u = 30*tmp*u4; 
  
  *BT   = B;
  *dBT  = dB_u;
  if (d2BT != NULL) *d2BT = d2B_u;
}


/** EOB Metric potentials A(r), B(r), and their derivatives, no spin version */
void eob_metric(double r, double prstar, Dynamics *dyn, double *A, double *B, double *dA, double *d2A, double *dB, double *d2B,
                double *Q, double *dQ, double *dQ_dprstar, double *ddQ_drdprstar, double *d2Q_dprstar2)
{
  const double nu    = EOBPars->nu;
  const double u     = 1./r;
  const double u2    = u*u;
  const double u3    = u2*u;
  const double u4    = u2*u2;
  const double u6    = u2*u4;

  double Atmp=0., dAtmp_u=0., d2Atmp_u=0.;
  double D=0., dD_u=0., d2D_u=0.;
  double Btmp=0., dBtmp_r=0., d2Btmp_r=0.;
  double Qtmp=0., dQtmp_du=0., dQtmp_dprstar=0., ddQtmp_drdprstar=0., d2Qtmp_dprstar2=0.;

  /* A potential and derivative with respect to u */  
  eob_metric_Apotential(r, nu, &Atmp, &dAtmp_u, &d2Atmp_u);

  /* Add here tides if needed */
  if (EOBPars->use_tidal) {
    double AT, dAT_u, d2AT_u;
    double BT, dBT_u, d2BT_u;
    eob_metric_Atidal(r, dyn, &AT, &dAT_u, &d2AT_u);
    Atmp     += AT;
    dAtmp_u  += dAT_u;
    d2Atmp_u += d2AT_u;
#if (USEBTIDALPOTENTIAL)
    eob_metric_Btidal(r, dyn, &BT, &dBT_u, &d2BT_u);
    Btmp     += BT;
    dBtmp_r  += -dBT_u*u2; 
#endif
  }

  /* A potential and derivative with respect to r */  
  *A   = Atmp;
  *dA  = -dAtmp_u*u2;
  *d2A = 2.*dAtmp_u*u3 + d2Atmp_u*u4;

  /* D potential and derivative with respect to r */
  eob_metric_Dpotential(r, nu, &D, &dD_u, &d2D_u); // this gives dD wtr to u
  double dD  = - u2*dD_u; // derivative wrt to r
  double d2D = 2.*dD_u*u3 + d2D_u*u4;

  /* B potential and derivative with respect to r */
  Btmp     = D/(Atmp);
  dBtmp_r  = (Btmp)*(dD/D - (*dA)/(*A));
  d2Btmp_r = SQ(dBtmp_r)/(Btmp) + (Btmp)*(d2D/D - SQ(dD/D) - (*d2A)/(*A) + SQ((*dA)/(*A)));

  *B   = Btmp;
  *dB  = dBtmp_r;
  *d2B = d2Btmp_r;

  /* Q potential and derivatives */
  eob_metric_Qpotential(r, prstar, nu, &Qtmp, &dQtmp_du, &dQtmp_dprstar, &ddQtmp_drdprstar, &d2Qtmp_dprstar2);
  *Q             = Qtmp;
  *dQ            = -u2*dQtmp_du; // derivative wrt to r
  *dQ_dprstar    = dQtmp_dprstar;
  *ddQ_drdprstar = ddQtmp_drdprstar;
  *d2Q_dprstar2  = d2Qtmp_dprstar2;

}
 
/** EOB Metric potentials A(r), B(r), and their derivatives, spin version */
void eob_metric_s(double r, double prstar, Dynamics *dyn, double *A, double *B, double *dA, double *d2A, double *dB, double *d2B,
                  double *Q, double *dQ, double *dQ_dprstar, double *ddQ_drdprstar, double *d2Q_dprstar2)
{

  const double nu    = EOBPars->nu;
  const double a1    = EOBPars->a1;
  const double a2    = EOBPars->a2;
  const double aK2   = EOBPars->aK2;
  const double C_Q1  = EOBPars->C_Q1;
  const double C_Q2  = EOBPars->C_Q2;
  const double C_Oct1 = EOBPars->C_Oct1;
  const double C_Oct2 = EOBPars->C_Oct2;
  const double C_Hex1 = EOBPars->C_Hex1;
  const double C_Hex2 = EOBPars->C_Hex2;
  const int usetidal = EOBPars->use_tidal;

  const double u   = 1./r;
  const double u2  = u*u;
  const double u3  = u2*u;
  const double u4  = u2*u2;
  
  double rc, drc, d2rc;
  eob_dyn_s_get_rc(r, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc, &drc, &d2rc);

  /* A potential and derivative with respect to u */  
  double Aorb, dAorb_u, d2Aorb_u;
  eob_metric_Apotential(rc, nu, &Aorb, &dAorb_u, &d2Aorb_u);

  /* Add here tides if needed */
  if (usetidal) {
    double AT, dAT_u, d2AT_u;
    double BT, dBT_u, d2BT_u;
    eob_metric_Atidal(rc, dyn, &AT, &dAT_u, &d2AT_u);
    Aorb     += AT;
    dAorb_u  += dAT_u;
    d2Aorb_u += d2AT_u;
#if (USEBTIDALPOTENTIAL)
    /* eob_metric_Btidal(rc, dyn, &BT, &dBT_u, &d2BT_u); */
#endif    
  }

  /* A potential and derivative with respect to r */  
  double uc  = 1./rc;
  double uc2 = uc*uc;
  double uc3 = uc2*uc;
  double uc4 = uc2*uc2;

  double dAorb  = -dAorb_u*uc2*drc;
  double d2Aorb = 2.*dAorb_u*uc3*SQ(drc) + d2Aorb_u*uc4*SQ(drc) - uc2*dAorb_u*d2rc;

  /* Correct A for spin */
  double AKerr_Multipole   = (1.+2.*uc)/(1.+2.*u);
  double dAKerr_Multipole  = -2.*uc2/(1.+2.*u)*drc + 2.*u2*(1.+2.*uc)/SQ(1.+2.*u);
  double d2AKerr_Multipole = -2.*uc2/(1.+2.*u)*d2rc - 8.*u2*uc2/SQ(1.+2.*u)*drc + 4.*uc3/(1.+2.*u)*SQ(drc)
    - 4.*u3*(1+2.*uc)/SQ(1.+2.*u) + 8.*u4*(1+2.*uc)/SQ(1.+2.*u)/(1.+2.*u);
  double fss = 1.;
  
  *A   = Aorb*AKerr_Multipole*fss;
  *dA  = dAorb*AKerr_Multipole + Aorb*dAKerr_Multipole;
  *d2A = d2Aorb*AKerr_Multipole + 2.*dAorb*dAKerr_Multipole + Aorb*d2AKerr_Multipole;
  //*d2A = d2Aorb*(1.+2.*uc)/(1.+2.*u) + 4.*dAorb*( u2*(1.+2.*uc)/((1.+2.*u)*(1.+2.*u)) - uc2/(1.+2.*u)*drc) + Aorb*(-4.*u3*(1.+2.*uc)/((1.+2.*u)*(1.+2.*u)) + 8.*u4*(1.+2.*uc)/((1.+2.*u)*(1.+2.*u)*(1.+2.*u))+4.*uc3/(1.+2.*u)*drc*drc - 2.*uc2/(1.+2.*u)*d2rc - 8.*u2*uc2/SQ(1.+2.*u)*drc); /* expanded *correct* form */

  /* D potential and derivative with respect to r */
  double D=0., dD_uc=0., d2D_uc=0.;
  eob_metric_Dpotential(rc, nu, &D, &dD_uc, &d2D_uc); // this gives dD wtr to uc

  // derivatives wrt to r
  double dD  = -uc2*drc*dD_uc; 
  double d2D = (2.*uc3*SQ(drc) - uc2*d2rc)*dD_uc + SQ(uc2*drc)*d2D_uc;

  /* B potential and derivative with respect to r */
  double fact   = r*r*uc2;
  double dfact  = 2.*r*uc2 - 2.*r*r*uc3*drc;
  double d2fact = 2.*uc2 - 8.*r*uc3*drc + 6.*r*r*uc4*SQ(drc) - 2.*r*r*uc3*d2rc;
  
  *B   = fact*D/(*A);
  *dB  = (*B)*(dfact/fact + dD/D - (*dA)/(*A));
  *d2B = SQ(*dB)/(*B) + (*B)*(d2fact/fact - SQ(dfact/fact) + d2D/D - SQ(dD/D)
			 - (*d2A)/(*A) + SQ((*dA)/(*A)));

  /* Q potential and derivatives */
  double Qtmp=0., dQtmp_duc=0., dQtmp_dprstar=0., ddQtmp_drcdprstar=0., d2Qtmp_dprstar2=0.;
  eob_metric_Qpotential(rc, prstar, nu, &Qtmp, &dQtmp_duc, &dQtmp_dprstar, &ddQtmp_drcdprstar, &d2Qtmp_dprstar2);
  *Q             = Qtmp;
  *dQ            = - uc2*drc*dQtmp_duc; // derivative wrt to r
  *dQ_dprstar    = dQtmp_dprstar;
  *ddQ_drdprstar = drc*ddQtmp_drcdprstar;
  *d2Q_dprstar2  = d2Qtmp_dprstar2;
  
}
