/** \file TEOBResumSMetric.c
 *  \brief TEOBResumS metric functions
 * 
 *  This file contains the functions that compute the EOB metric potentials
 *  A,D,Q and their derivatives.
 */


#include "TEOBResumS.h"

/**
 *  Function : eob_metric_A5PNlog
 *  ---------------
 *   EOB Metric A function 5PN log  
 *   This function computes the Pade' (1,5) resummed A function (with its
 *   derivatives) starting from the 5PN-expanded version of the A function
 *   including 4PN and 5PN log terms.
 *   This represents the current, stable, most accurate implementation of
 *   the EOB effective potential
 *   
 *   Coefficients a5 and a6 are the nonlog contributions to the 4PN and 5PN terms.
 *   In practice, a5 is fixed to its GSF value computed in Akcay et al,
 *   
 *   \f$a^5 \equiv a^5_{GSF} = +23.50190(5) \approx +23.5 \f$
 *   
 *   and \f$a^6 \equiv a^6(\nu) = (-110.5 + 129*(1-4*\nu)).*(1-1.5\times 10^{-5}/((0.26-\nu)^2)\f$
 *   as obtained from comparison with the Caltech-Cornell-CITA numerical data.
 *  These values are used as default.
 *
 *    @param[in]  r    : radial separation  
 *    @param[in]  nu   : symmetric mass ratio 
 *    @param[out] A    : A potential evaluated at r  
 *    @param[out] dA   : dA/du,   with u=1/r  
 *    @param[out] d2A  : d2A/du2, with u=1/r   
 */
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

/**
 *  Function : eob_metric_A5PNlog
 *  ---------------
 *   EOB Metric A function GSF-informed  
 *
 *    @param[in]  r    : radial separation  
 *    @param[in]  nu   : symmetric mass ratio 
 *    @param[out] A    : A potential evaluated at r  
 *    @param[out] dA   : dA/du,   with u=1/r  
 *    @param[out] d2A  : d2A/du2, with u=1/r   
 */
void eob_metric_AGSF(double r, double nu, double *A, double *dA, double *d2A)
{

  /* shortcuts */
  double u       = 1./r;
  double u2      = u*u;
  double u3      = u*u2;
  double u4      = u2*u2;
  double u5      = u4*u;
  double u6      = u5*u;
  double u7      = u6*u;
  double logu    = log(u);
  double logu2   = logu*logu;
  double onebyu  = r;
  double onebyu2 = r*r;

  /* The a function is factorized into two pieces:
  integer and tail part (semi-integer powers of u) */

  /* Integer part */
  double ResumInt, dResumInt, d2ResumInt;
  double coeffs_int[6], Dcoeffs_int[5], D2coeffs_int[5];

  // Coefficients of the Taylor-expanded function (+ derivatives), to be fed to the Padé

  double k1  =  9.843951347218796;
  double k2a =  16.29867029491112;
  double k2b =  6.4;
  double k3a = -64.26554415558926;
  double k3b = -30.15238095238095;
  double k4a =  309.0094399882301;
  double k4b = -20.58977072310406;
  double k5a =  969.247628430881;
  double k5b =  457.2474034136427;
  double k5c = -26.08761904761905;
  double k6a = -5268.885414280484;
  double k6b = -356.236938531493;
  double k6c =  110.4765532879819;

  coeffs_int[0] =  k1;
  coeffs_int[1] =  k2a + k2b*logu;
  coeffs_int[2] =  k3a + k3b*logu;
  coeffs_int[3] =  k4a + k4b*logu;
  coeffs_int[4] =  k5a + k5b*logu + k5c*logu2;
  coeffs_int[5] =  k6a + k6b*logu + k6c*logu2;

  Dcoeffs_int[0] =  k2b*onebyu;
  Dcoeffs_int[1] =  k3b*onebyu;
  Dcoeffs_int[2] =  k4b*onebyu;
  Dcoeffs_int[3] =  k5b*onebyu + 2.*k5c*logu*onebyu;
  Dcoeffs_int[4] =  k6b*onebyu + 2.*k6c*logu*onebyu;

  D2coeffs_int[0] = -k2b*onebyu2;
  D2coeffs_int[1] = -k3b*onebyu2;
  D2coeffs_int[2] = -k4b*onebyu2;
  D2coeffs_int[3] = -k5b*onebyu2 + k5c*(2.*onebyu2 - 2.*logu*onebyu2);
  D2coeffs_int[4] = -k6b*onebyu2 + k6c*(2.*onebyu2 - 2.*logu*onebyu2);

  Pade33_forGSF(coeffs_int, Dcoeffs_int, D2coeffs_int, u, &ResumInt, &dResumInt, &d2ResumInt);

  /* Tail part */
  double ResumTail, dResumTail, d2ResumTail;
  double coeffs_tail[4], Dcoeffs_tail[2], D2coeffs_tail[2];

  double j1  =  40.97833617482458;
  double j2  = -601.9566375425603;
  double j3a =  5160.094217567668;
  double j3b = -262.2613515188773;
  double j4a = -33044.17255491994;
  double j4b =  7335.733860726202;
 
  // Coefficients of the Taylor-expanded function, to be fed to the Padé
  coeffs_tail[0] = j1;
  coeffs_tail[1] = j2;
  coeffs_tail[2] = j3a + j3b*logu;
  coeffs_tail[3] = j4a + j4b*logu;

  Dcoeffs_tail[0] = j3b*onebyu;
  Dcoeffs_tail[1] = j4b*onebyu;

  D2coeffs_tail[0] = -j3b*onebyu2;
  D2coeffs_tail[1] = -j4b*onebyu2;

  Pade76v1_forGSF(coeffs_tail, Dcoeffs_tail, D2coeffs_tail, u, &ResumTail, &dResumTail, &d2ResumTail);

  /* a1SF function purely analytical + derivatives wrt to u */
  double a1SF_tmp     = ResumInt*ResumTail;
  double da1SF_u_tmp  = dResumInt*ResumTail + ResumInt*dResumTail; 
  double d2a1SF_u_tmp = d2ResumInt*ResumTail + 2*dResumInt*dResumTail + ResumInt*d2ResumTail;

  /* GSF-informed correction and improved functions */
  double f1SF     = 1 + (963.4328658734496 + (811.3827589325526 - 205.25579381038455*logu)*logu)*u7;
  double df1SF_u  = (7555.4128200467 + (5269.167724907099 - 1436.7905566726918*logu)*logu)*u6; 
  double d2f1SF_u = (50601.6446451873 + (28741.42523609721 - 8620.743340036151*logu)*logu)*u5;

  double a1SF     = a1SF_tmp*f1SF;
  double da1SF_u  = da1SF_u_tmp*f1SF + a1SF_tmp*df1SF_u;
  double d2a1SF_u = d2a1SF_u_tmp*f1SF + 2.*df1SF_u*da1SF_u_tmp + d2f1SF_u*a1SF_tmp;

  /* Schwarzschild Hamiltonian & derivatives */
  double Hschw     = (1. - 2.*u)/sqrt(1. - 3.*u);
  double dHschw_u  = -2./sqrt(1. - 3.*u) + (3.*(1. - 2.*u))/(2.*pow(1. - 3.*u, 1.5));
  double d2Hschw_u = -6./pow(1. - 3.*u, 1.5) + (27.*(1. - 2.*u))/(4.*pow(1. - 3.*u, 2.5));

  /* A potential and derivatives wrt to u */
  *A   = 1. - 2.*u + 2.*nu*u3*Hschw*a1SF;
  *dA  = -2. + 6.*nu*u2*a1SF*Hschw + 2.*nu*u3*dHschw_u*a1SF + 2.*nu*u3*Hschw*da1SF_u;
  *d2A = da1SF_u*(12.*u2*nu*Hschw + 2.*nu*u3*dHschw_u) + 2.*nu*u3*Hschw*d2a1SF_u + 2.*nu*u3*dHschw_u*da1SF_u + 
         a1SF*(12.*nu*u*Hschw + 12.*nu*u2*dHschw_u + 2.*nu*u3*d2Hschw_u);

}

/**
 *  Function : eob_metric_A5PNlogP33
 *  --------------------------------
 *    EOB Metric A function 5PN log resummed with (3, 3) Pade 
 *    See 2108.02043 and refs. therein
 *
 *    @param[in]  r    : radial separation  
 *    @param[in]  nu   : symmetric mass ratio 
 *    @param[out] A    : A potential evaluated at r  
 *    @param[out] dA   : dA/du,   with u=1/r  
 *    @param[out] d2A  : d2A/du2, with u=1/r   
 */
void eob_metric_A5PNlogP33(double r, double nu, double *A, double *dA, double *d2A)
{

  /* shortcuts */
  double nu2 = nu*nu;
  double pi2 = Pi*Pi;
  double pi4 = pi2*pi2;
  double pi6 = pi4*pi2;
  double ln2  = log(2);
  double ln3  = log(3);

  double a5c0       = -4237./60 + 2275./512*pi2 + 256./5*ln2 + 128./5*EulerGamma;
  double a5c1       = -221./6   + 41./32*pi2;
  double a5         =  a5c0 + nu*a5c1;
  double a6         =  EOBPars->a6c;

  double u       = 1./r;
  double u2      = u*u;
  double u3      = u*u2;
  double u4      = u2*u2;
  double logu    = log(u);

  double a5tot = a5  + 64./5*logu;
  double a6tot = a6  + (-7004./105 - 144./5*nu)*logu;

  

  /*----------------------------------------------------------
   preliminary: coefficients of the Pade and its derivatives
  ----------------------------------------------------------*/
  double N0     = -1737228288.  + 3538944.*a5tot + 7077888.*nu + 142073856.*pi2 - 2904768.*pi4;
  double N1     = -27216576512. - 3538944.*a6tot + 3338735616.*pi2 - 136524096.*pi4 + 1860867.*pi6 + a5tot*(110886912. - 4534272.*pi2);
  double D0     =  N0;
  double D1     = -18432.*(192.*a6tot + nu*(6016.- 246.*pi2) + a5tot*(-3008. + 123.*pi2));

  double D2     = -192.*(9216.*a5tot*a5tot + 18432.*a5tot*nu + (-3008. + 123.*pi2)*(96.*a6tot + nu*(3008. - 123.*pi2)));
  double D3     =  nu*(-3538944.*a6tot - 36864.*a5tot*(-3008. + 123.*pi2) + (-3008. + 123.*pi2)*SQ(-3008. + 123.*pi2) );

  /*1st derivatives of the above coefficients entering danu/du*/
  double dN0 = 226492416./(5*u);
  double dN1 = (-7077888*(-8186 - 504*nu + 287*pi2))/(35.*u);
  double dD0 = dN0;
  double dD1 = (-1179648*(-4*(7015 + 756*nu) + 861*pi2))/(35.*u);
  double dD2 = (-24576*(4128768*logu + 5*(64512*a5 + 64*(82297 + 36540*nu) \
                -123*(1751 + 756*nu)*pi2)))/(175.*u);
  double dD3 = (-7077888*nu*(-8186 - 504*nu + 287*pi2))/(35.*u);

  /* 2nd derivatives of the coefficients */
  double d2N0 = -226492416./(5.*u2);
  double d2N1 = -(-7077888*(-8186 - 504*nu + 287*pi2))/(35.*u2);
  double d2D0 = d2N0;
  double d2D1 = -(-1179648*(-4*(7015 + 756*nu) + 861*pi2))/(35.*u2);
  double d2D2 = -(-24576*(22206272 + 322560*a5 + 4128768*logu - 1076865*pi2 \
                + nu*(11692800 - 464940*pi2) ))/(175.*u2);
  double d2D3 = -(-7077888*nu*(-8186 - 504*nu + 287*pi2))/(35.*u2);

  /* the small a(u,nu) function*/
  double Num  = N0 + N1*u;
  double Den  = D0 + D1*u + D2*u2 + D3*u3;
  double anu  = Num/Den;

  /* derivative of anu with respect to u*/
  double dNum    = dN0 + N1 + u*dN1;
  double dDen    = dD0 + D1 + u*dD1 + 2.*u*D2 + u2*dD2 + 3*u2*D3 + u3*dD3;
  double danu_du = anu*(dNum/Num - dDen/Den);

  /* second derivative of anu with respect to u*/
  double d2Num = d2N0 + 2*dN1 + u*d2N1;
  double d2Den = u3*d2D3 + u2*(d2D2+6*dD3) + u*(d2D1 + 4*dD2 + 6*D3) + d2D0 + 2*dD1 + 2*D2;
  double d2anu_du2 = anu/(Num*Den)*(2.*dDen*dDen*anu - 2.*dNum*dDen + Den*d2Num - d2Den*Num);

  /* derivatives of A with respect to u */
  double dA_du  = -2 + 6*nu*u2*anu + 2*nu*u3*danu_du;
  double d2A_du = 12*nu*u*anu + 12*nu*u2*danu_du + 2*nu*u3*d2anu_du2;

  /* the A function and its derivatives with respect to r */
  /*
  *dA  =  -u2*dA_du;
  *d2A =   u4*d2A_du + 2*u3*dA_du;
  */

  /* output derivatives with respect to u */
  *A   =  1-2*u + 2*nu*u3*anu;
  *dA  =  dA_du;
  *d2A =  d2A_du;
}


/** EOB Metric D function at 3PN, resummed */


/**
 *  Function : eob_metric_D3PN
 *  --------------------------
 *    EOB Metric D function @3PN resummed with Pade 03
 *
 *    @param[in]  r    : radial separation  
 *    @param[in]  nu   : symmetric mass ratio 
 *    @param[out] D    : D potential evaluated at r  
 *    @param[out] dD   : dD/du,   with u=1/r  
 *    @param[out] d2D  : d2D/du2, with u=1/r   
 */
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

/**
 *  Function : eob_metric_DGSF
 *  --------------------------
 *    EOB Metric D function GSF-informed
 * 
 *    @param[in]  r    : radial separation  
 *    @param[in]  nu   : symmetric mass ratio 
 *    @param[out] D    : D potential evaluated at r  
 *    @param[out] dD   : dD/du,   with u=1/r  
 *    @param[out] d2D  : d2D/du2, with u=1/r   
 */
void eob_metric_DGSF(double r, double nu, double *D, double *dD, double *d2D)
{

  /* shortcuts */
  double u      = 1./r;
  double u2     = u*u;
  double u3     = u*u2;
  double u4     = u2*u2;
  double u5     = u4*u;
  double u6     = u5*u;
  double u7     = u6*u;
  double u11    = u5*u6;
  double u11by2 = sqrt(u11);
  double u13by2 = sqrt(u11*u2);
  double u15by2 = sqrt(u11*u4);
  double logu   = log(u);
  double logu2  = logu*logu;
  double onebyu  = r;
  double onebyu2 = r*r;

  /* The d function is factorized into two pieces:
  integer and tail part (semi-integer powers of u) */

  /* Integer part */
  double ResumInt, dResumInt, d2ResumInt;
  double coeffs_int[6], Dcoeffs_int[5], D2coeffs_int[5];

  // Coefficients of the Taylor-expanded function (+ derivatives), to be fed to the Padé

  double k1  =  8.666666666666667;
  double k2a =  36.92866519869136;
  double k2b =  6.577777777777778;
  double k3a = -79.89813242312009;
  double k3b = -33.80952380952381;
  double k4a =  777.9303456579357;
  double k4b = -118.6962962962963;
  double k5a =  5757.134871176116;
  double k5b =  1884.740681815503;
  double k5c = -115.9449735449735;
  double k6a = -42331.4249981053;
  double k6b = -3615.490356674048;
  double k6c =  879.6352532123961;

  coeffs_int[0] =  k1;
  coeffs_int[1] =  k2a + k2b*logu;
  coeffs_int[2] =  k3a + k3b*logu;
  coeffs_int[3] =  k4a + k4b*logu;
  coeffs_int[4] =  k5a + k5b*logu + k5c*logu2;
  coeffs_int[5] =  k6a + k6b*logu + k6c*logu2;

  Dcoeffs_int[0] =  k2b*onebyu;
  Dcoeffs_int[1] =  k3b*onebyu;
  Dcoeffs_int[2] =  k4b*onebyu;
  Dcoeffs_int[3] =  k5b*onebyu + 2.*k5c*logu*onebyu;
  Dcoeffs_int[4] =  k6b*onebyu + 2.*k6c*logu*onebyu;

  D2coeffs_int[0] = -k2b*onebyu2;
  D2coeffs_int[1] = -k3b*onebyu2;
  D2coeffs_int[2] = -k4b*onebyu2;
  D2coeffs_int[3] = -k5b*onebyu2 + k5c*(2.*onebyu2 - 2.*logu*onebyu2);
  D2coeffs_int[4] = -k6b*onebyu2 + k6c*(2.*onebyu2 - 2.*logu*onebyu2);

  Pade33_forGSF(coeffs_int, Dcoeffs_int, D2coeffs_int, u, &ResumInt, &dResumInt, &d2ResumInt);

  /* Tail part */
  double ResumTail, dResumTail, d2ResumTail;
  double coeffs_tail[4], Dcoeffs_tail[2], D2coeffs_tail[2];
 
  // Coefficients of the Taylor-expanded function (+ derivatives), to be fed to the Padé

  double j1  =  88.07496559797366;
  double j2  = -1395.220274500357;
  double j3a =  8663.505746839271;
  double j3b = -579.3375514888934;
  double j4a =  29169.61933412867;
  double j4b =  14219.676593263639;

  coeffs_tail[0] =  j1;
  coeffs_tail[1] =  j2;
  coeffs_tail[2] =  j3a + j3b*logu;
  coeffs_tail[3] =  j4a + j4b*logu;

  Dcoeffs_tail[0] = j3b*onebyu;
  Dcoeffs_tail[1] = j4b*onebyu;

  D2coeffs_tail[0] = -j3b*onebyu2;
  D2coeffs_tail[1] = -j4b*onebyu2;

  Pade76v1_forGSF(coeffs_tail, Dcoeffs_tail, D2coeffs_tail, u, &ResumTail, &dResumTail, &d2ResumTail);

  /* a1SF function purely analytical + derivatives wrt to u */
  double d1SF_tmp     = ResumInt*ResumTail;
  double Dd1SF_u_tmp  = dResumInt*ResumTail + ResumInt*dResumTail; 
  double D2d1SF_u_tmp = d2ResumInt*ResumTail + 2*dResumInt*dResumTail + ResumInt*d2ResumTail;

  /* GSF-informed correction and improved function */

  double f1SF     = 1.*(1. + (-9.898639062741917e6 + 2.6012818800187605e6*logu)*u15by2 + 
                    (1.0018072583960408e7 + 2.6340749138477915e6*logu + 198577.79085898417*logu2)*u7);
  double df1SF_u  = 1.*(2.6012818800187605e6*u13by2 + 
                    7.5*(-9.898639062741917e6 + 2.6012818800187605e6*logu)*u13by2 + 
                    7.*(1.0018072583960408e7 + 2.6340749138477915e6*logu + 198577.79085898417*logu2)*u6 + 
                    (2.6340749138477915e6/u + (397155.58171796834*logu)/u)*u7);
  double d2f1SF_u = 1.*(3.641794632026264e7*u11by2 + 48.75*(-9.898639062741917e6 + 2.6012818800187605e6*logu)*u11by2 + 
                    42.*(1.0018072583960408e7 + 2.6340749138477915e6*logu + 198577.79085898417*logu2)*u5 + 
                    14.*(2.6340749138477915e6/u + (397155.58171796834*logu)/u)*u6 + 
                    (-2.236919332129823e6/u2 - (397155.58171796834*logu)/u2)*u7);

  double d1SF     = d1SF_tmp*f1SF;
  double Dd1SF_u  = Dd1SF_u_tmp*f1SF + d1SF_tmp*df1SF_u;
  double D2d1SF_u = D2d1SF_u_tmp*f1SF + 2.*df1SF_u*Dd1SF_u_tmp + d2f1SF_u*d1SF_tmp;

  /* D potential and derivatives wrt u */
  double barD   = 1. + 6.*nu*u2*d1SF;
  double barD2  = barD*barD;
  double barD3  = barD2*barD;
  *D            = 1./barD;

  double DbarD  = nu*(12.*u*d1SF + 6.*u2*Dd1SF_u);
  *dD           = - DbarD/barD2;

  double D2barD = nu*(12.*d1SF + 24.*u*Dd1SF_u + 6.*u2*D2d1SF_u);
  *d2D          = 2.*(DbarD*DbarD)/barD3 - D2barD/barD2;
  
}


/*
  EOB metric D potential rewritten with (3, 2) Pade 
  (Pade coefficients already computed in matlab)
  (2108.02043 and refs. therein)
*/

/**
 *  Function : eob_metric_DGSF
 *  --------------------------
 *    EOB metric D potential rewritten with (3, 2) Pade 
 *    2108.02043 and refs. therein
 * 
 *    @param[in]  r    : radial separation  
 *    @param[in]  nu   : symmetric mass ratio 
 *    @param[out] D    : D potential evaluated at r  
 *    @param[out] dD   : dD/du,   with u=1/r  
 *    @param[out] d2D  : d2D/du2, with u=1/r   
 */
void eob_metric_D5PNP32(double r, double nu, double *D, double *dD, double *d2D)
{

  /* shortcuts */
  double nu2 = nu*nu;
  double pi2 = Pi*Pi;
  double pi4 = pi2*pi2;
  double pi6 = pi4*pi2;
  double ln2  = log(2);
  double ln3  = log(3);

  double u       = 1./r;
  double u2      = u*u;
  double u3      = u*u2;
  double u4      = u2*u2;
  double logu    = log(u);


  // only analytically uncalculated 5PN coefficient set to zero 
  double d5nu2 = 0.;        
  double d2    = -6.*nu;
  double d3    = -52.*nu + 6.*nu2;

  double c1    = 533./45 - 1184./15*EulerGamma + 23761./1536*pi2  + 6496./15*ln2 - 2916./5*ln3;
  double c2    = 296. - 123./16*pi2;
  double d4c   = c1 + nu*c2;
  double d4log = -592./15;
  double d4    = nu*(d4c + d4log*logu);
  // its derivative
  double Dd4   = nu*d4log/u;
  double D2d4  = -nu*d4log/u2;

  double d5c   = (-294464./175) + (2840./7)*EulerGamma + (-120648./35)*ln2 + (19683./7)* \
                ln3 + ((-2216./105) - d5nu2 + (6784./15)*EulerGamma + (326656./21)* \
                ln2 + (-58320./7)*ln3)*nu + (63707./512)*pi2 + nu2*((-1285./3) + (205./16) \
                *pi2);
  double d5log = (1420./7) + (3392./15)*nu;
  double d5    = nu*(d5c  +  d5log*logu);
  double Dd5   = nu*d5log/u;
  double D2d5  = -nu*d5log/u2;
  

  // define powers of coefficients
  double d2_2 = d2*d2;
  double d2_3 = d2_2*d2;
  double d3_2 = d3*d3;
  double d3_3 = d3_2*d3;
  double d3_4 = d3_2*d3_2;
  double d4_2 = d4*d4;
  double Dd4_2 = Dd4*Dd4;

  // we define here this variable to then comfortably 
  // write powers of it
  double factor = 1./(d3_2 - d2*d4);
  double factor2 = factor*factor;
  double factor3 = factor2*factor;

  double N1   = factor*((-1)*d3*d4 + d2*d5);
  double N2   = -factor*((-1)*d2*d3_2 + d2_2*d4 + (-1)*d4_2 + d3*d5);
  double N3   = factor*(d3_3 + (-2)*d2*d3*d4 + d2_2*d5);

  double dN1  = factor2*(((-1)*d3_3 + d2_2*d5)*Dd4 + d2*( \
                d3_2 + (-1)*d2*d4)*Dd5);
  double dN2  = factor2*((-1)*d2*d4_2*Dd4 + d3*d4*(2.*d3* \
                Dd4 + d2*Dd5) + (-1)*d3*(d2*d5*Dd4 + d3_2*Dd5));
  double dN3  = d2*factor2*(((-1)*d3_3 + d2_2*d5)*Dd4 + d2* \
                (d3_2 + (-1)*d2*d4)*Dd5);

  double d2N1 = factor3*((-1)*(d3_2 + (-1)*d2*d4)*(d2* \
                D2d5*((-1)*d3_2 + d2*d4) + D2d4*(d3_3 + (-1)*d2_2*d5)) + 2.*((-1) \
                *d2*d3_3 + d2_3*d5)*Dd4_2 + 2.*d2_2*(d3_2 + (-1)*d2*d4)* \
                Dd4*Dd5);
  double d2N2 = factor3*((-1)*(d3_2 + (-1)*d2*d4)*(D2d5* \
                d3_3 + (-1)*d3*(d2*D2d5 + 2.*D2d4*d3)*d4 + d2*D2d4*d4_2 + d2* \
                D2d4*d3*d5) + 2.*(d3_4 + (-1)*d2_2*d3*d5)*Dd4_2 + 2.*d2*d3*(( \
                -1)*d3_2 + d2*d4)*Dd4*Dd5);
  double d2N3 = d2*factor3*((-1)*(d3_2 + (-1)*d2*d4)*(d2* \
                D2d5*((-1)*d3_2 + d2*d4) + D2d4*(d3_3 + (-1)*d2_2*d5)) + 2.*((-1) \
                *d2*d3_3 + d2_3*d5)*Dd4_2 + 2.*d2_2*(d3_2 + (-1)*d2*d4)* \
                Dd4*Dd5);

  double D1   = factor*((-1)*d3*d4 + d2*d5);
  double D2   = -factor*((-1)*d4_2 + d3*d5);

  double dD1  = factor2*(((-1)*d3_3 + d2_2*d5)*Dd4 + d2*( \
                d3_2 + (-1)*d2*d4)*Dd5);
  double dD2  = factor2*(d2*(d4_2 + (-1)*d3*d5)*Dd4 + ((-1)* \
                d3_2 + d2*d4)*((-2)*d4*Dd4 + d3*Dd5));

  double d2D1 = factor3*((-1)*(d3_2 + (-1)*d2*d4)*(d2* \
                D2d5*((-1)*d3_2 + d2*d4) + D2d4*(d3_3 + (-1)*d2_2*d5)) + 2.*((-1) \
                *d2*d3_3 + d2_3*d5)*Dd4_2 + 2.*d2_2*(d3_2 + (-1)*d2*d4)* \
                Dd4*Dd5);
  double d2D2 = factor3*((-1)*(d3_2 + (-1)*d2*d4)*(D2d5* \
                d3_3 + (-1)*d3*(d2*D2d5 + 2.*D2d4*d3)*d4 + d2*D2d4*d4_2 + d2* \
                D2d4*d3*d5) + 2.*(d3_4 + (-1)*d2_2*d3*d5)*Dd4_2 + 2.*d2*d3*(( \
                -1)*d3_2 + d2*d4)*Dd4*Dd5);
  
  double Num = 1. + N1*u + N2*u2 + N3*u3;
  double Den = 1. + D1*u + D2*u2;
  
  double dNum =  N1 +u*dN1 + 2*N2*u + dN2*u2 + 3*N3*u2 + dN3*u3;
  double dDen =  D1 +u*dD1 + 2*D2*u + dD2*u2;
  double d2Num = 2*dN1 +u*d2N1 + 2*N2  + 4*dN2*u + d2N2*u2 + 6*N3*u + 6*dN3*u2 + d2N3*u3;
  double d2Den = 2*dD1 +u*d2D1 + 2*D2  + 4*dD2*u + d2D2*u2;
  
  // Finally, the D function and its derivative in r
  // *dD = -u2*(dNum*Den - Num*dDen)/(Den*Den);
  // *d2D = -2*u*dD + u4*( (d2Num*Den-Num*d2Den)/(Den*Den)\
  //     -2.*dDen*(dNum*Den-Num*dDen)/(Den*Den*Den));

  // Output the derivatives in u
  *D   =  Num/Den;
  *dD  = (dNum*Den - Num*dDen)/(Den*Den);
  *d2D = (d2Num*Den - Num*d2Den)/(Den*Den)\
        -2.*dDen*(dNum*Den - Num*dDen)/(Den*Den*Den);

}

/**
  *  Function : eob_metric_Q3PN
  *  --------------------------
  *    EOB Metric Q function @3PN
  *
  *    @param[in]  r    : radial separation  
  *    @param[in]  nu   : symmetric mass ratio 
  *    @param[out] Q    : Q potential evaluated at r  
  *    @param[out] dQ   : dQ/du,   with u=1/r  
  *    @param[out] dQ_dprstar  : dQ/dprstar
  *    @param[out] d2Q_du2   : d2Q/du2, with u=1/r
  *    @param[out] ddQ_drdprstar : d2Q/drdprstar
  *    @param[out] d2Q_dprstar2  : d2Q/dprstar^2
  *    @param[out] d3Q_dr2dprstar : d3Q/dr^2dprstar
  *    @param[out] d3Q_drdprstar2 : d3Q/drdprstar^2
  *    @param[out] d3Q_dprstar3   : d3Q/dprstar^3
*/
void eob_metric_Q3PN(double r, double prstar, double nu, double *Q, double *dQ_du, double *dQ_dprstar, 
                     double *d2Q_du2, double *ddQ_drdprstar, double *d2Q_dprstar2,
                     double *d3Q_dr2dprstar, double *d3Q_drdprstar2, double *d3Q_dprstar3)
{
  const double z3 = 2.*nu*(4. - 3.*nu);
  double u  = 1./r;
  double u2 = u*u;
  double u3 = u2*u;
  double u4 = u3*u;
  double prstar2 = prstar*prstar;
  double prstar3 = prstar2*prstar;
  double prstar4 = prstar2*prstar2;

  *Q              =  z3*u2*prstar4;
  *dQ_du          =  2.*z3*u*prstar4;
  *dQ_dprstar     =  4.*z3*u2*prstar3;
  *d2Q_du2        =  2.*z3*prstar4;
  *ddQ_drdprstar  = -8.*z3*u3*prstar3;
  *d2Q_dprstar2   =  12.*z3*u2*prstar2;
  *d3Q_dr2dprstar =  24.*z3*u4*prstar3;
  *d3Q_drdprstar2 = -24.*z3*u3*prstar2;
  *d3Q_dprstar3   =  24.*z3*u2*prstar;

}

/**
  *  Function : eob_metric_Q3PN
  *  --------------------------
  *     EOB Metric Q function GSF-informed
  *
  *    @param[in]  r    : radial separation  
  *    @param[in]  nu   : symmetric mass ratio 
  *    @param[out] Q    : Q potential evaluated at r  
  *    @param[out] dQ   : dQ/du,   with u=1/r  
  *    @param[out] dQ_dprstar  : dQ/dprstar
  *    @param[out] d2Q_du2   : d2Q/du2, with u=1/r
  *    @param[out] ddQ_drdprstar : d2Q/drdprstar
  *    @param[out] d2Q_dprstar2  : d2Q/dprstar^2
  *    @param[out] d3Q_dr2dprstar : d3Q/dr^2dprstar
  *    @param[out] d3Q_drdprstar2 : d3Q/drdprstar^2
  *    @param[out] d3Q_dprstar3   : d3Q/dprstar^3
*/
void eob_metric_QGSF(double r, double prstar, double nu, double *Q, double *dQ_du, double *dQ_dprstar, 
                     double *d2Q_du2, double *ddQ_drdprstar, double *d2Q_dprstar2,
                     double *d3Q_dr2dprstar, double *d3Q_drdprstar2, double *d3Q_dprstar3)
{

  double u       = 1./r;
  double u2      = u*u;
  double u3      = u*u2;
  double u4      = u2*u2;
  double u5      = u4*u;
  double u6      = u5*u;
  double u7      = u6*u;
  double sqrtu   = sqrt(u);
  double u15by2  = sqrt(u5*u5*u5);
  double logu    = log(u);
  double logu2   = logu*logu;
  double onebyu  = r;
  double onebyu2 = r*r;
  double prstar2 = prstar*prstar;
  double prstar3 = prstar2*prstar;
  double prstar4 = prstar2*prstar2;

  /* The q function is factorized into two pieces:
  integer and tail part (semi-integer powers of u) */

  /* Integer part */
  double ResumInt, dResumInt, d2ResumInt;
  double coeffs_int[6], Dcoeffs_int[5], D2coeffs_int[5];

  // Coefficients of the Taylor-expanded function (+ derivatives), to be fed to the Padé

  double k1  =  3.588880535619494;
  double k2a = -12.14327341036888;
  double k2b =  6.461904761904762;
  double k3a = -13.57177599359104;
  double k3b = -122.5626984126984;
  double k4a =  5849.419891636204;
  double k4b =  1690.691389474068;
  double k4c = -80.70857142857143;
  double k5a = -89749.62943631351;
  double k5b = -8625.740422940584;
  double k5c =  1436.31537414966;
  double k6a =  445969.7971893991;
  double k6b = -70286.9426913128;
  double k6c = -5402.297846034037;

  coeffs_int[0] =  k1;
  coeffs_int[1] =  k2a + k2b*logu;
  coeffs_int[2] =  k3a + k3b*logu;
  coeffs_int[3] =  k4a + k4b*logu + k4c*logu2;
  coeffs_int[4] =  k5a + k5b*logu + k5c*logu2;
  coeffs_int[5] =  k6a + k6b*logu + k6c*logu2;

  Dcoeffs_int[0] =  k2b*onebyu;
  Dcoeffs_int[1] =  k3b*onebyu;
  Dcoeffs_int[2] =  k4b*onebyu + 2.*k4c*logu*onebyu;
  Dcoeffs_int[3] =  k5b*onebyu + 2.*k5c*logu*onebyu;
  Dcoeffs_int[4] =  k6b*onebyu + 2.*k6c*logu*onebyu;

  D2coeffs_int[0] = -k2b*onebyu2;
  D2coeffs_int[1] = -k3b*onebyu2;
  D2coeffs_int[2] = -k4b*onebyu2 + k4c*(2.*onebyu2 - 2.*logu*onebyu2);
  D2coeffs_int[3] = -k5b*onebyu2 + k5c*(2.*onebyu2 - 2.*logu*onebyu2);
  D2coeffs_int[4] = -k6b*onebyu2 + k6c*(2.*onebyu2 - 2.*logu*onebyu2);

  Pade33_forGSF(coeffs_int, Dcoeffs_int, D2coeffs_int, u, &ResumInt, &dResumInt, &d2ResumInt);

  /* Tail part */
  double ResumTail, dResumTail, d2ResumTail;
  double coeffs_tail[5], Dcoeffs_tail[3], D2coeffs_tail[3];

  // Coefficients of the Taylor-expanded function, to be fed to the Padé

  double j1  =  18.43046912376822;
  double j2  = -397.3442276527582;
  double j3a =  2180.644722347771;
  double j3b = -119.09593619501652;
  double j4a =  64484.07468368325;
  double j4b =  650.5464012653043;
  double j5a = -1.366299312664882e6;
  double j5b = -5839.988023296008;
  double j5c =  2257.0834309598104;

  coeffs_tail[0] =  j1;
  coeffs_tail[1] =  j2;
  coeffs_tail[2] =  j3a + j3b*logu;
  coeffs_tail[3] =  j4a + j4b*logu;
  coeffs_tail[4] =  j5a + j5b*logu + j5c*logu2;

  Dcoeffs_tail[0] = j3b*onebyu;
  Dcoeffs_tail[1] = j4b*onebyu;
  Dcoeffs_tail[2] = j5b*onebyu + 2.*j5c*logu*onebyu;

  D2coeffs_tail[0] = -j3b*onebyu2;
  D2coeffs_tail[1] = -j4b*onebyu2;
  D2coeffs_tail[2] = -j5b*onebyu2 + j5c*(2.*onebyu2 - 2.*logu*onebyu2);
  
  Pade76v2_forGSF(coeffs_tail, Dcoeffs_tail, D2coeffs_tail, u, &ResumTail, &dResumTail, &d2ResumTail);

  double dIntTail = dResumInt*ResumTail + ResumInt*dResumTail;

  /* q1SF function purely analytical + derivatives wrt to u */
  double q1SF_tmp     = 8.*u2*ResumInt*ResumTail;
  double dq1SF_u_tmp  = 16.*u*ResumInt*ResumTail + 8.*u2*dIntTail;
  double d2q1SF_u_tmp = 16.*ResumInt*ResumTail + 32.*u*dIntTail + 8.*u2*(d2ResumInt*ResumTail + 2*dResumInt*dResumTail + ResumInt*d2ResumTail);

  /* GSF-informed correction and improved functions */

  double f1SF     = 1.*(1. + (6.473498944710975e9 - 1.6413472171574914e9*logu)*u15by2 + (-6.555554279481237e9 - 1.7756296887290301e9*logu - 1.4618867042954516e8*logu2)*u7);
  double df1SF_u  = (-4.766450964509769e10 + logu*(-1.2721785161962301e10 - 1.0233206930068161e9*logu - 1.2310104128681187e10*sqrtu) + 4.690989486817482e10*sqrtu)*u6;
  double d2f1SF_u = (-2.987088430325484e11 + logu*(-7.837735235778743e10 - 6.139924158040897e9*logu - 8.001567683642772e10*sqrtu) + 2.926042125144552e11*sqrtu)*u5;

  double q1SF     = q1SF_tmp*f1SF;
  double dq1SF_u  = dq1SF_u_tmp*f1SF + q1SF_tmp*df1SF_u;
  double d2q1SF_u = d2q1SF_u_tmp*f1SF + 2.*df1SF_u*dq1SF_u_tmp + d2f1SF_u*q1SF_tmp;

  /* complete Q function and derivatives */

  *Q              =  nu*q1SF*prstar4;
  *dQ_du          =  nu*dq1SF_u*prstar4;
  *dQ_dprstar     =  4.*nu*q1SF*prstar3;
  *d2Q_du2        =  nu*d2q1SF_u*prstar4;
  *ddQ_drdprstar  = -4.*u2*nu*dq1SF_u*prstar3;
  *d2Q_dprstar2   =  12.*nu*q1SF*prstar2;
  *d3Q_dr2dprstar =  4.*u3*nu*(2.*dq1SF_u + u*d2q1SF_u)*prstar3;
  *d3Q_drdprstar2 = -12.*nu*u2*dq1SF_u*prstar2;
  *d3Q_dprstar3   =  24.*nu*q1SF*prstar;

}

/**
  *  Function : eob_metric_Q3PN
  *  --------------------------
  *    EOB Metric Q function at 5PN, only local part 
  *    2108.02043 and refs. therein
  *    @param[in]  r    : radial separation  
  *    @param[in]  nu   : symmetric mass ratio 
  *    @param[out] Q    : Q potential evaluated at r  
  *    @param[out] dQ   : dQ/du,   with u=1/r  
  *    @param[out] dQ_dprstar  : dQ/dprstar
  *    @param[out] d2Q_du2   : d2Q/du2, with u=1/r
  *    @param[out] ddQ_drdprstar : d2Q/drdprstar
  *    @param[out] d2Q_dprstar2  : d2Q/dprstar^2
  *    @param[out] d3Q_dr2dprstar : d3Q/dr^2dprstar
  *    @param[out] d3Q_drdprstar2 : d3Q/drdprstar^2
  *    @param[out] d3Q_dprstar3   : d3Q/dprstar^3
*/
void eob_metric_Q5PNloc(double r, double prstar, double nu, double *Q, double *dQ_du, double *dQ_dprstar, 
                     double *d2Q_du2, double *d2Q_drdprstar, double *d2Q_dprstar2,
                     double *d3Q_dr2dprstar, double *d3Q_drdprstar2, double *d3Q_dprstar3)
{

  /* shortcuts */
  double pi2     = Pi*Pi;
  double nu2     = nu*nu;
  double nu3     = nu2*nu;
  double nu4     = nu3*nu;

  double prstar2 = prstar*prstar;
  double prstar3 = prstar2*prstar;
  double prstar4 = prstar2*prstar2;
  double prstar6 = prstar4*prstar2;
  double prstar8 = prstar4*prstar4;

  double u     = 1./r;
  double u2    = u*u;
  double u3    = u2*u;
  double u4    = u3*u;
  double du_dr = -u2;

  double z3  = 2.0*nu*(4.0-3.0*nu);
  double q42 = z3;


  double q62    = -2.783007636952232489*nu - 5.4*nu2 + 6.*nu3;
  double q43    = 92.711044284955949757*nu - 131.*nu2 + 10.*nu3;

  double q44loc = 453.70625272357658631*nu - 1081.8908931907688518*nu2 + \
                  602.31854041656388904*nu3;
  double q63loc = -21.096091643426789868*nu - 78.6*nu2 + 188.*nu3 - 14.*nu4;
  double q82loc = 6./7*nu + 18./7*nu2 + 24./7*nu3 - 6.*nu4;

  double q82    = -1.5175121231629802709e6*nu + 3.3384202338272640486*nu2 + \
                  3.4285714285714285714*nu3 - 6.*nu4;
  double q63    = -217866.45869782848474*nu - 89.529832736260960467*nu2 + 188.*nu3 - 14.*nu4;
  double q44c   = 602.31854041656388904*nu3 + -1796.1366049802412531*nu2 + \
                  452.54216699669657073*nu;
  double q44log = 51.695238095238095238*nu - 118.4*nu2;
  double q44    = q44c + q44log*log(u);




  /* Q potential and all its derivatives */

  *Q = u2*prstar4*( q42 + q43*u + q44loc*u2 \
       + prstar2*( q62 + q63loc*u + q82loc*prstar2 ));
  
  *dQ_du = u*prstar4*( 2.*q42 + 3.*q43*u + 4.*q44loc*u2 \
          + prstar2*( 2.*q62 + 3.*q63loc*u + 2.*q82loc*prstar2 ));
  
  double dQ_dprstar2 = u2*prstar2*( 2.*q42 + 2.*q43*u + 2.*q44loc*u2 \
                      + prstar2*( 3.*q62 + 3.*q63loc*u + 4.*q82loc*prstar2 ));
  
  double d2Q_dudprstar2 = u*prstar2*( 4.*q42 + 6.*q43*u + 8.*q44loc*u2 \
                      + prstar2*(6.*q62 + 9.*q63loc*u + 8.*q82loc*prstar2 ));   

  *d2Q_du2 = prstar4*( 2.*q42 + 6.*q43*u + 12.*q44loc*u2  \
            + prstar2*( 2.* q62 + 2.*q82loc*prstar2 + 6.*q63loc*u ));               

  
  double d2Q_dprstar22 = 2.*q42*u2 + 2.*q43*u3 + 2.*q44loc*u4 \
                      + prstar2*( 6.*q62*u2 + 6.*q63loc*u3 + 12.*q82loc*u2*prstar2 ); 

  
  double d3Q_dprstar23 = 6.*q62*u2+ 6.*q63loc*u3 + 24.*q82loc*u2*prstar2; 
  

  double d3Q_du2dprstar2 = prstar2*( 4.*q42 + 12.*q43*u + 24.*q44loc*u2 \
                        + prstar2*( 6.*q62+ 8.*q82loc*prstar2 + 18.*q63loc*u )); 

  double d3Q_dudprstar22 = 4.*q42*u + 6.*q43*u2 + 8.*q44loc*u3 \
                          + prstar2*( 12.*q62*u + 18.*q63loc*u2 + 24.*q82loc*u*prstar2);  



  /* We also translate derivatives from prstar2 to prstar
     and some from u to r as needed for output*/


    *d2Q_drdprstar   = -2*prstar*u2*d2Q_dudprstar2;

    *d3Q_dr2dprstar  = prstar*( 4*u3*d2Q_dudprstar2 + 2*u4*d3Q_du2dprstar2 );

    *d3Q_drdprstar2  = -2*u2*d2Q_dudprstar2 - 4*prstar2*u2*d3Q_dudprstar22;

    *dQ_dprstar   = 2.*prstar*dQ_dprstar2;

    *d2Q_dprstar2 = 2.*dQ_dprstar2 + 4.*prstar2*d2Q_dprstar22;

    *d3Q_dprstar3 = prstar*( 12.*d2Q_dprstar22 + 8.*prstar2*d3Q_dprstar23 );

}

/* Macro for the bar_alpha coefs */
#define bar_alph(cA,kapA,cB,kapB,kap)( ((cA)*(kapA)+(cB)*(kapB))/(kap) )

/**
 *  Function : eob_metric_Atidal
 *  ---------------
 *   Tidal potential, three version implemented: 
 *   1. TEOB NNLO, Bernuzzi+ 1205.3403
 *   2. TEOBResum, Bini&Damour, 1409.6933, Bernuzzi+ 1412.4553 
 *   3. TEOBResum3, Akcay+ 1812.02744
 *   
 *   f-mode resonance model from Hinderer&Steinhoff can be added to all three.
 *
 *   @param[in]  r     : radial separation
 *   @param[in]  dyn   : EOB dynamics 
 *   @param[out] AT    : tidal A potential evaluated at r  
 *   @param[out] dAT   : tidal dA/du,   with u=1/r  
 *   @param[out] d2AT  : tidal d2A/du2, with u=1/r 
 */ 
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
  
  const double kapA2j = EOBPars->japA2;
  const double kapB2j = EOBPars->japB2;
  const double kapT2j = EOBPars->japT2;  
  
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
    
  *AT   = A;
  *dAT  = dA_u;
  if (d2AT != NULL) *d2AT = d2A_u;

}

/**
 *  Function : eob_metric_Btidal
 *  ---------------
 *   Tidal B potential
 *   Vines, Flanagan 1PN term
 *
 *   @param[in]  r     : radial separation
 *   @param[in]  dyn   : EOB dynamics 
 *   @param[out] BT    : tidal B potential evaluated at r  
 *   @param[out] dBT   : tidal dB/du,   with u=1/r  
 *   @param[out] d2BT  : tidal d2B/du2, with u=1/r 
 */ 
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


/**
 *  Function : eob_metric
 *  ---------------------
 *   EOB Metric potentials A(r), B(r), and their derivatives, no spin version
 *
 *   @param[in]  r     : radial separation
 *   @param[in]  dyn   : EOB dynamics 
 *   @param[out] A     : A potential evaluated at r  
 *   @param[out] B     : B potential evaluated at r    
 *   @param[out] dA    : dA/dr
 *   @param[out] d2A   : d2A/dr2
 *   @param[out] dB    : dB/dr 
 *   @param[out] d2B   : d2B/dr2
 *   @param[out] Q     : Q potential evaluated at r
 *   @param[out] dQ    : dQ/dr
 *   @param[out] dQ_dprstar : dQ/dprstar
 *   @param[out] d2Q   : d2Q/dr2
 *   @param[out] ddQ_drdprstar : ddQ/drdprstar
 *   @param[out] d2Q_dprstar2 : d2Q/dprstar2
 *   @param[out] d3Q_dr2dprstar : d3Q/dr2dprstar
 *   @param[out] d3Q_drdprstar2 : d3Q/drdprstar2
 *   @param[out] d3Q_dprstar3 : d3Q/dprstar3
 */ 
void eob_metric(double r, double prstar, Dynamics *dyn, double *A, double *B, double *dA, double *d2A, double *dB, double *d2B,
                double *Q, double *dQ, double *dQ_dprstar, double *d2Q, double *ddQ_drdprstar, double *d2Q_dprstar2,
                double *d3Q_dr2dprstar, double *d3Q_drdprstar2, double *d3Q_dprstar3)
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
  double Qtmp=0., dQtmp_du=0., dQtmp_dprstar=0., d2Qtmp_du2=0., ddQtmp_drdprstar=0., d2Qtmp_dprstar2=0.,
         d3Qtmp_dr2dprstar= 0., d3Qtmp_drdprstar2=0., d3Qtmp_dprstar3=0.;

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
  eob_metric_Qpotential(r, prstar, nu, &Qtmp, &dQtmp_du, &dQtmp_dprstar, &d2Qtmp_du2, &ddQtmp_drdprstar, &d2Qtmp_dprstar2, &d3Qtmp_dr2dprstar, &d3Qtmp_drdprstar2, &d3Qtmp_dprstar3);
  *Q              = Qtmp;
  *dQ             = -u2*dQtmp_du; // derivative wrt to r
  *dQ_dprstar     = dQtmp_dprstar;
  *d2Q            = u3*(2.*dQtmp_du + u*d2Qtmp_du2); // derivative wrt to r
  *ddQ_drdprstar  = ddQtmp_drdprstar;
  *d2Q_dprstar2   = d2Qtmp_dprstar2;
  *d3Q_dr2dprstar = d3Qtmp_dr2dprstar;
  *d3Q_drdprstar2 = d3Qtmp_drdprstar2;
  *d3Q_dprstar3   = d3Qtmp_dprstar3;

}
 
/**
 *  Function : eob_metric_s
 *  ---------------------
 *   EOB Metric potentials A(r), B(r), and their derivatives, spin version
 *
 *   @param[in]  r     : radial separation
 *   @param[in]  dyn   : EOB dynamics 
 *   @param[out] A     : A potential evaluated at r  
 *   @param[out] B     : B potential evaluated at r    
 *   @param[out] dA    : dA/dr
 *   @param[out] d2A   : d2A/dr2
 *   @param[out] dB    : dB/dr 
 *   @param[out] d2B   : d2B/dr2
 *   @param[out] Q     : Q potential evaluated at r
 *   @param[out] dQ    : dQ/dr
 *   @param[out] dQ_dprstar : dQ/dprstar
 *   @param[out] d2Q   : d2Q/dr2
 *   @param[out] ddQ_drdprstar : ddQ/drdprstar
 *   @param[out] d2Q_dprstar2 : d2Q/dprstar2
 *   @param[out] d3Q_dr2dprstar : d3Q/dr2dprstar
 *   @param[out] d3Q_drdprstar2 : d3Q/drdprstar2
 *   @param[out] d3Q_dprstar3 : d3Q/dprstar3
 */
 void eob_metric_s(double r, double prstar, Dynamics *dyn, double *A, double *B, double *dA, double *d2A, double *dB, double *d2B,
                  double *Q, double *dQ, double *dQ_dprstar, double *d2Q, double *ddQ_drdprstar, double *d2Q_dprstar2,
                  double *d3Q_dr2dprstar, double *d3Q_drdprstar2, double *d3Q_dprstar3)
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
  double d2AKerr_Multipole = -2.*uc2/(1.+2.*u)*d2rc - 8.*u2*uc2/SQ(1.+2.*u)*drc + 4.*uc3/(1.+2.*u)*SQ(drc) - 4.*u3*(1+2.*uc)/SQ(1.+2.*u) + 8.*u4*(1.+2.*uc)/SQ(1.+2.*u)/(1.+2.*u);
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
  double d2fact = 2.*uc2   - 8.*r*uc3*drc + 6.*r*r*uc4*SQ(drc) - 2.*r*r*uc3*d2rc;
  
  *B   = fact*D/(*A);
  *dB  = (*B)*(dfact/fact + dD/D - (*dA)/(*A));
  *d2B = SQ(*dB)/(*B) + (*B)*(d2fact/fact - SQ(dfact/fact) + d2D/D - SQ(dD/D)
			 - (*d2A)/(*A) + SQ((*dA)/(*A)));

  /* Q potential and derivatives */
  double Qtmp=0., dQtmp_duc=0., dQtmp_dprstar=0.,  d2Qtmp_duc2=0., ddQtmp_drcdprstar=0., d2Qtmp_dprstar2=0.,
         d3Qtmp_drc2dprstar=0., d3Qtmp_drcdprstar2=0., d3Qtmp_dprstar3=0.;
  eob_metric_Qpotential(rc, prstar, nu, &Qtmp, &dQtmp_duc, &dQtmp_dprstar, &d2Qtmp_duc2, &ddQtmp_drcdprstar, &d2Qtmp_dprstar2,
                        &d3Qtmp_drc2dprstar, &d3Qtmp_drcdprstar2, &d3Qtmp_dprstar3);
  *Q              = Qtmp;
  *dQ             = - uc2*drc*dQtmp_duc; // derivative wrt to r
  *dQ_dprstar     = dQtmp_dprstar;
  *d2Q            = uc2*(-d2rc*dQtmp_duc + SQ(drc)*(2.*uc*dQtmp_duc + uc2*d2Qtmp_duc2)); // derivative wrt to r
  *ddQ_drdprstar  = drc*ddQtmp_drcdprstar;
  *d2Q_dprstar2   = d2Qtmp_dprstar2;
  *d3Q_dr2dprstar = d2rc*ddQtmp_drcdprstar + SQ(drc)*d3Qtmp_drc2dprstar;
  *d3Q_drdprstar2 = drc*d3Qtmp_drcdprstar2;
  *d3Q_dprstar3   = d3Qtmp_dprstar3;
}
