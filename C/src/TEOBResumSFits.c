/**  \file TEOBResumSFits.c
 *   \brief Various fits employed throughout the code
 * 
 *   This file contains all of the fits employed throughout the code,
 *   grouped depending on their purpose.
 * 
 *   Note that may of the functions are not employed in the current
 *   version of the code, but are kept for historical reasons and backwards
 *   compatibility.
 */

#include "TEOBResumS.h"

/** \defgroup a6cfits a6c fits
 *  Various fits of the pseudo-5PN NR-informed parameter a6c
 *  @{
 */

/**
 *  Function: eob_a6c_fit
 *  ---------------
 *   Fit of a6c from TEOBResumS paper Nagar et al. (2018) 
 *   
 *   @param[in] nu: symmetric mass ratio
 *
 *   @return   a6c
 */
double eob_a6c_fit(double nu)
{
  return 3097.3*nu*nu - 1330.6*nu + 81.38;
}

/**
 *  Function: eob_a6c_fit_HM
 *  ---------------
 *   Fit of a6c from higher multipoles paper: arXiv:2001.09082
 *   
 *   @param[in] nu: symmetric mass ratio
 *
 *   @return  a6c
 */
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

/**
 *  Function: eob_a6c_fit_ecc
 *  ------------------------
 *   Fit of a6c from arXiv:2101.08624
 *   
 *   @param[in] nu: symmetric mass ratio
 *
 *   @return  a6c
 */
double eob_a6c_fit_ecc(double nu)
{
  double nu2 = nu*nu;
  
  const double n0 = -0.50395;
  const double n1 =-4.8547;
  const double n2 = 52.96;
  const double n3 = 20.7013;
    
  return (n0 + n1*nu + n2*nu2)*exp(n3*nu);
}

/**
 * Function: eob_a6c_fit_next
 * --------------------------
 *   Fit of a6c from 2108.02043
 *   @note Currently not interfaced
 * 
 *   @param[in] nu: symmetric mass ratio
 *   
 *   @return  a6c
*/
double eob_a6c_fit_next(double nu)
{
  
  const double n0 = -4.6416;
  const double n1 = -503.57;
  const double n2 = 599.96;
    
  return n0 + nu*(n1 + n2*nu);
}

/**
 *  Function: eob_a6c_fit_ecc_P33_4PNh22
 *  ------------------------------------
 *   New fit for a6c obtained using the 4PN term
 *   in the rho22 with P(2,2) resummation.
 *   Eq. (34) of arXiv:2404.05288
 *   
 *   @param[in] nu: symmetric mass ratio
 *
 *   @return  a6c
 */
double eob_a6c_fit_ecc_P33_4PNh22(double nu)
{
  return -32.5953 - 269.4331*nu;
}

/**
 *  Function: eob_a6c_fit_HM_2023
 *  -----------------------------
 *   Fit of a6c from Tab. II of arXiv:2304.09662 
 *   
 *   @param[in] nu: symmetric mass ratio
 *
 *   @return  a6c
 */
double eob_a6c_fit_HM_2023(double nu)
{
  const double nu2 = nu*nu;
  const double nu3 = nu2*nu;

  const double n0 = 46.5524;
  const double n1 =-24.2516;
  const double n2 = 120.9594;
  const double n3 = -167.2242;
  const double d1 =-3.3998;
  
  return n0*(1 + n1*nu + n2*nu2 + n3*nu3)/(1 + d1*nu);
}

/**
 *  Function: eob_a6c_fit_ecc_P33_newlogs
 *  ------------------------------------
 *   New fit for a6c obtained using the new resummation
 *   scheme for the log terms
 *   Eq. (34) of arXiv:2407.04762
 *   
 *   @param[in] nu: symmetric mass ratio
 *
 *   @return  a6c
 */
double eob_a6c_fit_ecc_P33_newlogs(double nu)
{
  return 34.85 + nu*(-318.26 + nu*208.19);
}

/** @} */ // end of a6cfits

/** \defgroup c3fits cN3LO fits
 *  Various fits of the spin-orbit N3LO NR-informed parameter cN3LO
 *  @{
 */

/**
 *  Function: eob_c3_fit_global
 *  ---------------------------
 *   Fit of c3, TEOBResumS paper Nagar et al. (2018)
 *   Note: c3 = 0 with tides
 *
 *    @param[in] nu: symmetric mass ratio
 *    @param[in] a1: spin of body 1
 *    @param[in] a2: spin of body 2
 *
 *   @return  cN3LO
 */
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

/**
 *  Function: eob_c3_fit_HM
 *  -----------------------
 *   Fit of c3, HM paper: arXiv:2001.09082
 *   Note: c3 = 0 with tides
 *    @param[in] nu: symmetric mass ratio
 *    @param[in] a1: spin of body 1
 *    @param[in] a2: spin of body 2
 *
 *   @return  cN3LO
 */
double eob_c3_fit_HM(double nu, double a1, double a2)
{  
  const double nu2 = nu*nu;
  const double X12 = sqrt(1.-4.*nu);
  const double a0  = a1+a2;
  const double a02 = a0*a0;
  const double a03 = a02*a0;
  const double a04 = a03*a0;
  
  /* Equal-mass, equal-spin coefficients */
  const double p0 =  45.235903;//46.323097;
  const double n1 =  -1.688708;//-1.418623;
  const double n2 =   0.787959;//0.457118;
  const double n3 =  -0.018080;//0.000435;
  const double n4 =   -0.001906;//0.097791;
  const double d1 =  -0.751479;//-0.519003;
		  
  /* Other coefficients */
  const double p1 =  47.3756;// 46.2914;
  const double p2 = -36.1964;//-103.0997;
	
  
  const double c3 = p0*(1 + n1*a0 + n2*a02 + n3*a03 + n4*a04)/(1 + d1*a0)
    + p1*nu*X12*a0 + p2*nu2*(a1 - a2);
  
  return c3;
}

/**
 *  Function: eob_c3_fit_ecc
 *  ------------------------
 *   Fit of c3, all in one paper: arXiv:2101.08624
 *   Note: c3 = 0 with tides
 *   @param[in] nu: symmetric mass ratio
 *   @param[in] a1: spin of body 1
 *   @param[in] a2: spin of body 2
 *
 *   @return  cN3LO
 */
double eob_c3_fit_ecc(double nu, double a1, double a2)
{
  const double nu2 = nu*nu;
  const double X12 = sqrt(1.-4.*nu);
  const double a0  = a1+a2;
  const double a02 = a0*a0;
  const double a03 = a02*a0;
  const double a04 = a03*a0;
  
  /* Equal-mass, equal-spin coefficients */
  const double p0 = 35.482253; //42.720039;
  const double n1 = -1.730483; //-0.710937;
  const double n2 = 1.144438; //-0.142503;
  const double n3 = 0.098420; //0.183764;
  const double n4 = -0.329288; //0.025030;
  const double d1 = -0.345207; //0.331651;

  /* Other coefficients */
  const double p1 = 244.505; //100.743;
  const double p2 = 148.184; //-5.01934; 
  const double p3 = -1085.35; //-61.6306; (for different functional form: nu*a0^2)
  
  const double c3 = p0*(1 + n1*a0 + n2*a02 + n3*a03 + n4*a04)/(1 + d1*a0)
    + p1*nu*X12*a0 + p2*nu2*(a1 - a2) + p3*nu2*X12*a0; //p3*nu*X12*a02;
  
  return c3;
}

/**
 *  Function: eob_c3_fit_ecc_P33_4PNh22
 *  -----------------------------------
 *   New fit for c3 obtained using the 4PN term
 *   in the rho22 with P(2,2) resummation 
 *   Note: c3 = 0 with tides
 *   From Tab. III of arXiv:2404.0528
 *
 *   @param[in] nu: symmetric mass ratio
 *   @param[in] a1: spin of body 1
 *   @param[in] a2: spin of body 2
 *
 *   @return  cN3LO
 */
double eob_c3_fit_ecc_P33_4PNh22(double nu, double a1, double a2)
{
  const double nu2 = nu*nu;
  const double X12 = sqrt(1.-4.*nu);
  const double a0  = a1+a2;
  const double a02 = a0*a0;
  const double a03 = a02*a0;
  const double a04 = a03*a0;

  /* Equal mass part */
  const double p0 =  38.625208;
  const double n1 = -0.105187;
  const double n2 = -0.758427;
  const double n3 =  0.183613;
  const double n4 =  0.057817;
  const double d1 =  0.905420;
  const double c3_eq = p0*(1. + n1*a0 + n2*a02 + n3*a03 + n4*a04)/(1.+ d1*a0);

  /* Unequal mass, unequal-spin part */
  const double c1 =  23.0579594801987;
  const double c2 =  12.544546974770155;
  const double c3 = -0.015710117458387414;
  const double c4 = -119.25957682313306;
  const double c5 =  64.47094882671716;
  const double c6 =  54.65682440916807;

  const double c3_neq = c1*a0*X12 + c2*a02*X12 + c3*a03*X12 + c4*a0*nu*X12 + c5*(a1 - a2)*nu2 + c6*SQ((a1 - a2))*nu2;
  return c3_eq + c3_neq;
}

/**
 *  Function: eob_c3_fit_HM_2023
 *  ----------------------------
 *   Fit of c3, Table IV of arXiv:2304.09662
 *   Note: c3 = 0 with tides
 *    @param[in] nu: symmetric mass ratio
 *    @param[in] a1: spin of body 1
 *    @param[in] a2: spin of body 2
 *
 *   @return  cN3LO
 */
double eob_c3_fit_HM_2023(double nu, double a1, double a2)
{
  const double nu2 = nu*nu;
  const double X12 = sqrt(1.-4.*nu);
  const double a0  = a1+a2;
  const double a02 = a0*a0;
  const double a03 = a02*a0;
  const double a04 = a03*a0;
  double *coeff_eq, *coeff_neq;

  double coeff_eq_420[6]  = {43.872788,-1.849495, 1.011208, -0.086453, -0.038378, -0.888154};
  double coeff_eq_431[6]  = {42.195044,-2.010717, 1.258034, -0.129593, -0.106295, -0.966525};
  
  double coeff_neq_420[6] = {26.553,    -8.65836,  0,        -84.7473,   24.0418,   0};
  double coeff_neq_430[6] = {16.695689, 2.025017, -6.600956, -53.146114, 34.097885, -101.003721};
  double coeff_neq_431[6] = {20.99561492580911, 1.5806110388013492, -10.428048396877976, -61.19803748038313, 37.11341346539421, -37.668146004322274};
  double coeff_neq_432[6] = {18.80033140113715, 0.6175269833645626, -10.397710332484472, -47.16963298818921, 33.4449328928248,  -32.51572099183539};

  switch(EOBPars->use_cN3LO_fits)
  {
    case(cN3LO_fits_HM_2023_420):
      coeff_eq  = coeff_eq_420;
      coeff_neq = coeff_neq_420;
      break;
    case(cN3LO_fits_HM_2023_430):
      coeff_eq  = coeff_eq_420;
      coeff_neq = coeff_neq_430;
      break;
    case(cN3LO_fits_HM_2023_431):
      coeff_eq  = coeff_eq_431;
      coeff_neq = coeff_neq_431;
      break;
    case(cN3LO_fits_HM_2023_432):
      coeff_eq  = coeff_eq_431;
      coeff_neq = coeff_neq_431;
  }

  double c3_eq  = coeff_eq[0]*(1 + coeff_eq[1]*a0 + coeff_eq[2]*a02 
                                 + coeff_eq[3]*a03 + coeff_eq[4]*a04)/(1 + coeff_eq[5]*a0);
  double c3_neq = coeff_neq[0]*a0*X12    + coeff_neq[1]*a02*X12       + coeff_neq[2]*a03*X12 
                + coeff_neq[3]*a0*nu*X12 + coeff_neq[4]*(a1 - a2)*nu2 + coeff_neq[5]*SQ(a1 - a2)*nu2;

  return c3_eq + c3_neq;
}

/**
 *  Function: eob_c3_fit_ecc_P33_newlogs
 *  -----------------------------------
 *   New fit for c3 obtained using the new
 *   resummation scheme for the log terms.
 *   Note: c3 = 0 with tides
 *   From Tab. II of arXiv:2407.04762
 *
 *   @param[in] nu: symmetric mass ratio
 *   @param[in] a1: spin of body 1
 *   @param[in] a2: spin of body 2
 *
 *   @return  cN3LO
 */
double eob_c3_fit_ecc_P33_newlogs(double nu, double a1, double a2)
{
  const double nu2 = nu*nu;
  const double X12 = sqrt(1.-4.*nu);
  const double a0  = a1+a2;
  const double a02 = a0*a0;
  const double a03 = a02*a0;
  const double a04 = a03*a0;

  /* Equal mass part */
  const double p0 =  43.588561;
  const double n1 = -1.623117;
  const double n2 =  0.912729;
  const double n3 = -0.103066;
  const double n4 = -0.076885;
  const double d1 = -0.653058;
  const double c3_eq = p0*(1. + n1*a0 + n2*a02 + n3*a03 + n4*a04)/(1.+ d1*a0);

  /* Unequal mass, unequal-spin part */
  const double c1 =  12.020173836641153;
  const double c2 = -2.7102976036661226;
  const double c3 = -0.35955651905612807;
  const double c4 = -36.64350436933381;
  const double c5 =  34.80670732430782;
  const double c6 = -85.97333184971158;

  const double c3_neq = c1*a0*X12 + c2*a02*X12 + c3*a03*X12 + c4*a0*nu*X12 + c5*(a1 - a2)*nu2 + c6*SQ((a1 - a2))*nu2;
  return c3_eq + c3_neq;
}
/** @} */ // end of c3fits


/**
 *  Function: eob_nqc_dtfit
 *  -----------------------
 *   Function providing a fit of Deltat_NQC vs chi, via a simple rational function. 
 *   
 *    @param[in] chi : 
 *    @param[in] chi0: 
 *
 *   @return  dtnqc
 */
double eob_nqc_dtfit(const double chi, const double chi0)
{
  const double n1 = -16.06288206;
  const double d1 = -4.04266459;

  double x     = chi-chi0;

  double dtnqc = (1.+n1*x)/(1.+d1*x);

  return dtnqc;
}

/**
 *  Function: eob_mrg_momg
 *  -----------------------
 *   Fit of GW frequency of NR merger. 
 *   
 *    @param[in] nu   : symmetric mass ratio
 *    @param[in] X1   : mass fraction of body 1
 *    @param[in] X2   : mass fraction of body 2
 *    @param[in] chi1 : dimensionless spin of body 1
 *    @param[in] chi2 : dimensionless spin of body 2
 *
 *   @return  omg_mrg
 */
double eob_mrg_momg(double nu, double X1, double X2, double chi1, double chi2)
{
  const double nu2 = SQ(nu);
  const double X12 = sqrt(1-4*nu);
  const double a1  = X1*chi1;
  const double a2  = X2*chi2;
  const double a0  = a1+a2;
  const double a12 = a1-a2;
  const double Shat = 0.5*(a0 + X12*a12);
  const double Shat2 = SQ(Shat);
  const double b[4] = {0.066045, -0.23876, 0.76819, -0.9201};
  return( 0.273356*(1+0.84074*nu+1.6976*nu2)*
	  (1+((-0.42311+b[0]*X12)/(1+b[1]*X12))*Shat
	   +((-0.066699))*Shat2)/(1+((-0.83053+b[2]*X12)/(1+b[3]*X12))*Shat) );
}

/** \defgroup nqcpoint NQC point
 *  Various functions to compute (A, dA, omg, domg) at the NR point 
 *  used to determine NQC corrections.
 *  Depending on the function, this point is assumed to be 2M after the 
 *  peak of the corresponding multipole,\f$ t^{NQC}_{lm} = t^{peak}_{lm} +2\f$,
 *  or at the peak of the 22 mode
 *
 *  @{
 */

/** 
 *  Function: eob_nqc_point
 *  -----------------------
 *   Fits for the NR point used to determine (2,2) NQC corrections, from 
 *   https://arxiv.org/abs/1806.01772 (Table X)
 *   https://arxiv.org/abs/1703.06814 (Table II)
 *   
 *   @param[out]  A_tmp   : RWZ-normalized Amplitude at the NQC point 
 *   @param[out]  dA_tmp  : first time derivative of the RWZ-normalized amplitude at the NQC point 
 *   @param[out]  omg_tmp : omega at the NQC point 
 *   @param[out]  domg_tmp: first time derivative of the frequency at the NQC point 
 *
 */
void eob_nqc_point(Dynamics *dyn, double *A_tmp, double *dA_tmp, double *omg_tmp, double *domg_tmp)
{

  const double nu   = EOBPars->nu;
  const double X1   = EOBPars->X1;
  const double X2   = EOBPars->X2;
  const double chi1 = EOBPars->chi1;
  const double chi2 = EOBPars->chi2;
  const double aK   = EOBPars->a1 + EOBPars->a2;

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

      
  if ( nu > 0.2485) {

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

  }

}

/** 
 *  Function: eob_nqc_point_HM
 *  -----------------------
 *   Fits for the NR point used to determine NQC corrections, from 
 *   https://arxiv.org/abs/2001.09082 (App C)
 *   
 *   @param[in, out]  A_tmp   : RWZ-normalized Amplitude at the NQC point 
 *   @param[in, out]  dA_tmp  : first time derivative of the RWZ-normalized amplitude at the NQC point 
 *   @param[in, out]  omg_tmp : omega at the NQC point 
 *   @param[in, out]  domg_tmp: first time derivative of the frequency at the NQC point 
 *
 */
void eob_nqc_point_HM(Dynamics *dyn, double *A_tmp, double *dA_tmp, double *omg_tmp, double *domg_tmp)
{

  const double nu   = EOBPars->nu;
  const double X1   = EOBPars->X1;
  const double X2   = EOBPars->X2;
  const double chi1 = EOBPars->chi1;
  const double chi2 = EOBPars->chi2;
  const double aK   = EOBPars->a1 + EOBPars->a2;

  const double nu2  = SQ(nu);
  const double nu3  = nu2*nu;
  const double X12  = X1 - X2;  
  const double aK2  = SQ(aK);
  const double aK3  = aK2*aK;
  const double aK4  = aK2*aK2;
  const double a12  = X1*chi1 - X2*chi2;
  const double Shat = 0.5*(aK + a12*X12);
  const double Shat2 = Shat*Shat;

  double pA[5], pdA[5], pomg[5], pdomg[5];
  double c_p1,     c_p2,     c_p3,   c_p4;
  double c_pdA1,   c_pdA2,   c_pdA3, c_pdA4;
  double c_pdomg1, c_pdomg2;
  double n0, d1;

  double a0_omg_tmp, a1_omg_tmp, a2_omg_tmp, b0_omg_tmp, b1_omg_tmp, b2_omg_tmp, b3_omg_tmp, c11_omg_tmp, c12_omg_tmp, c21_omg_tmp, c22_omg_tmp, c31_omg_tmp, c32_omg_tmp;
  double a0_domg_tmp, a1_domg_tmp, a2_domg_tmp, b1_domg_tmp, b2_domg_tmp, b3_domg_tmp, b4_domg_tmp, c1_domg_tmp, c2_domg_tmp;
  double scale, a0_A_tmp, a1_A_tmp , a2_A_tmp, b0_A_tmp, b1_A_tmp, b2_A_tmp, b3_A_tmp, b4_A_tmp, b5_A_tmp, b6_A_tmp, c11_A_tmp, c12_A_tmp, c21_A_tmp, c22_A_tmp, c31_A_tmp, c32_A_tmp;
  double a0_dA_tmp, a1_dA_tmp, a2_dA_tmp, b1_dA_tmp, b2_dA_tmp, b3_dA_tmp, c1_dA_tmp, c2_dA_tmp;
  double omg_tmp_nu, omg_tmp_equal, omg_tmp_den, domg_tmp_nu, domg_tmp_equal, A_tmp_scale_nu, A_tmp_scale_equal, A_tmp_scale_den, dA_tmp_scale_nu, dA_tmp_scale_equal;

  double p1[2], p2[2], p3[2], p4[2]; 
  double pn0[2], pd1[2], ppdomg1[2], ppdomg2[2], pdA1[2],pdA2[2],pdA3[2],pdA4[2];

  /* l=2, m=2 */
  if (nu>0.2485) {
	    
    pA[0]    =  0.00178195;
    pA[1]    =  0.00435589;
    pA[2]    =  0.00344489;
    pA[3]    = -0.00076165;
    pA[4]    =  0.31973334;
    A_tmp[1] =  pA[0]*aK4    + pA[1]*aK3   + pA[2]*aK2    + pA[3]*aK     + pA[4];
	      
    pdA[0]   =  0.00000927;
    pdA[1]   = -0.00024550;
    pdA[2]   =  0.00012469;
    pdA[3]   =  0.00123845;
    pdA[4]   = -0.00195014;
    dA_tmp[1] =  pdA[0]*aK4   + pdA[1]*aK3   + pdA[2]*aK2   + pdA[3]*aK   + pdA[4];
	    
    pomg[0]  =  0.00603482;
    pomg[1]  =  0.01604555;
    pomg[2]  =  0.02290799;
    pomg[3]  =  0.07084587;
    pomg[4]  =  0.38321834;
    omg_tmp[1] =  pomg[0]*aK4  + pomg[1]*aK3  + pomg[2]*aK2  + pomg[3]*aK  + pomg[4];
	    
    pdomg[0]  =  0.00024066;
    pdomg[1]  =  0.00038123;
    pdomg[2]  = -0.00049714;
    pdomg[3]  =  0.00041219;
    pdomg[4]  =  0.01190548;
    domg_tmp[1] =  pdomg[0]*aK4 + pdomg[1]*aK3 + pdomg[2]*aK2 + pdomg[3]*aK + pdomg[4];
	    
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
    A_tmp[1]   =  c_p1*aK3 + c_p2*aK2 + c_p3*aK + c_p4;
	      
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
    dA_tmp[1]  =  c_pdA1*aK3   + c_pdA2*aK2 + c_pdA3*aK+ c_pdA4;
	    
    pn0[0]     =  0.46908067;
    pn0[1]     =  0.27022141;
    pd1[0]     =  0.64131115;
    pd1[1]     = -0.37878384;
    n0         =  pn0[0]*nu + pn0[1];
    d1         =  pd1[0]*nu + pd1[1];
    omg_tmp[1] =  n0/(1 + d1*aK);
	    
    ppdomg1[0] =  0.00061175;
    ppdomg1[1] =  0.00074001;
    ppdomg2[0] =  0.02504442;
    ppdomg2[1] =  0.00548217;
    c_pdomg1   =  ppdomg1[0]*nu + ppdomg1[1];
    c_pdomg2   =  ppdomg2[0]*nu + ppdomg2[1];
    domg_tmp[1] =  c_pdomg1*aK   + c_pdomg2;
	  
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
    omg_tmp[1]    = omg_tmp_nu*omg_tmp_equal/omg_tmp_den;
	    
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
    domg_tmp[1]    = domg_tmp_nu + domg_tmp_equal;
	
    scale       = 1 - omg_tmp[1]*Shat;
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
    A_tmp[1]          = A_tmp_scale_nu*A_tmp_scale_equal/A_tmp_scale_den;
  
    scale       = 0.5*omg_tmp[1]/sqrt(6);
    a0_dA_tmp 	= -0.0011936600;
    a1_dA_tmp 	=  2.86637;
    a2_dA_tmp 	= -1.3667;
    b1_dA_tmp 	=  0.02679530;
    b2_dA_tmp 	= -0.0064409;
    c1_dA_tmp   = -0.015395218;
    c2_dA_tmp   = 0.008732589;
	
    dA_tmp_scale_nu    = a0_dA_tmp*(1. + a1_dA_tmp*nu + a2_dA_tmp*nu2)*scale;
    dA_tmp_scale_equal = (b1_dA_tmp + c1_dA_tmp*X12)*Shat + (b2_dA_tmp + c2_dA_tmp*X12)*Shat*Shat;
    dA_tmp[1]          = dA_tmp_scale_nu*dA_tmp_scale_equal;
	    
  }
  
  /* l=3, m=1 */
  A_tmp[2]    = 0.00520201*X12*(1 - 4.9441*nu + 8.9339*nu2);
  dA_tmp[2]   = -0.00043382*(1 - 9.0479*nu + 23.054*nu2)/(1 + 88.626*nu2);
  omg_tmp[2]  = 0.485186*(1 - 0.4421*nu - 6.8184*nu2);
  domg_tmp[2] = 0.0673118*(1 + 13.318*nu)/(1 + 70.552*nu);
	  
  /* l=3, m=3 */
  b1_A_tmp = +0.0067063;
  b2_A_tmp = +4.814781;
  b3_A_tmp = +0.0111876;
  b4_A_tmp = -1.079532;
  b5_A_tmp = +2.967227;
  b6_A_tmp = -2.571783;
	  
  A_tmp[4] = 0.0512928*X12*(1 + 0.09537*nu + 3.7217*nu2);
  A_tmp[4] += ( (0.0092449444 + b1_A_tmp*nu)/(1 + b2_A_tmp*nu)*a12 + (-0.000087052+b3_A_tmp*nu)/(1 + b4_A_tmp*nu)*SQ(a12))/(1 + (-0.4796687 + b5_A_tmp*nu)/(1 + b6_A_tmp*nu)*a12);
	
  b1_dA_tmp = -0.0004671176;
  b2_dA_tmp = -4.0270198;
	  
  dA_tmp[4] = -0.00039568*(1 + 1.0985*nu - 13.458*nu2);
  dA_tmp[4] += (0.0001415043*Shat + (+0.0001046803 + b1_dA_tmp*nu)/(1 + b2_dA_tmp*nu)*Shat2)/(1-0.422066*Shat);
	
  omg_tmp[4] = 0.476647*(1 + 1.0886*nu + 3.0658*nu2)*(1 + (-0.2362713 - 0.085544*nu)*Shat)/(1 + (-0.5828922 - 0.523365*nu)*Shat);        
	
  b1_domg_tmp = -0.00045666;
  b2_domg_tmp = +0.0290846;
  b3_domg_tmp = -0.00388909;
  b4_domg_tmp = +0.0087659;
	
  domg_tmp[4] = 0.0110394*(1 + 2.7962*nu) + (b1_domg_tmp+b2_domg_tmp*nu)*Shat + (b3_domg_tmp+b4_domg_tmp*nu)*Shat2;
	  
  /* l=4, m=1 */
  A_tmp[5]    = 0.00043987*X12*(1-2.*nu)*(1 - 8.4975*nu + 27.31*nu2)/(1 - 1.2002*nu);
  dA_tmp[5]   = -0.00001827*(1 - 2.8242*nu - 3.1871*nu2);
  omg_tmp[5]  = 0.735051*(1 - 8.3628*nu + 20.529*nu2)/(1 - 7.4883*nu + 18.695*nu2);
  domg_tmp[5] = 0.0739078*(1 + 0.99186*nu - 19.435*nu2);
	
  /* l=4, m=4 (only dA) */
  b1_dA_tmp = -0.000152614*(1. -7.63783*nu + 15.8089*nu2)/(1. - 5.88951*nu + 11.1555*nu2);  
  b2_dA_tmp = (0.0000376236 - 0.0000645958*nu)/(1. - 2.35613*nu);
  b3_dA_tmp = (-0.819379    - 298.678*nu)     /(1. - 1063.08*nu);
  dA_tmp[8] = b1_dA_tmp + b2_dA_tmp*Shat/(1. + b3_dA_tmp*Shat);
		
  /* l=5, m=5 */
  b1_A_tmp = -0.0829447;
  b2_A_tmp = +8.434479;
  b3_A_tmp = +1.197758;
  b4_A_tmp = -0.260704;
	
  A_tmp[13] = 0.00516272*X12*(1-2*nu)*(1 - 0.68569*nu + 7.7267*nu2) + 0.034503277967117711167*(0.04402930 + b1_A_tmp*nu)/(1 + b2_A_tmp*nu)*a12/(1 + (-0.5749409+b3_A_tmp*nu)/(1 + b4_A_tmp*nu)*a12);
	
  b1_dA_tmp = +0.00171106;
  b2_dA_tmp = +0.00475246;
  b3_dA_tmp = +6.508129;
	  
  dA_tmp[13] = -0.0000658001*(1 - 3.3212*nu) + 0.034503277967117711167*((0.0002250231 + b1_dA_tmp*nu)*a12 + (0.0002322527 + b2_dA_tmp*nu)*SQ(a12))/(1 + (-0.8223640 + b3_dA_tmp*nu)*a12);
	
  b1_omg_tmp = -1.711526;
  b2_omg_tmp = -1.546737;
  b3_omg_tmp = -2.023890;
	
  omg_tmp[13] =  0.855016*(1 + 1.1639*nu)*(1 + (-0.3577869 + b1_omg_tmp*nu)*Shat + (-0.0627735 + b2_omg_tmp*nu)*Shat2)/(1 + (-0.7053108 + b3_omg_tmp*nu)*Shat);
	
  b1_domg_tmp = -0.01040098;
  b2_domg_tmp = +0.1384369;
  b3_domg_tmp = -0.03244646;
  b4_domg_tmp = +0.2250689;
	  
  domg_tmp[13] =  0.0178326*(1 + 3.1304*nu) + (b1_domg_tmp + b2_domg_tmp*nu)*Shat + (b3_domg_tmp + b4_domg_tmp*nu)*Shat2;
	  
}

/** 
 *  Function: eob_nqc_point_HM_peak22
 *  ---------------------------------
 *   Fits for the NR point used to determine NQC corrections, available for 21,33,44
 *   Here the NQC point is assumed to be the peak of the 22 mode.
 *   From Pompili et. al. arXiv:2303.18039
 *   
 *   @param[in, out]  A_tmp   : RWZ-normalized Amplitude at the NQC point 
 *   @param[in, out]  dA_tmp  : first time derivative of the RWZ-normalized amplitude at the NQC point 
 *   @param[in, out]  omg_tmp : omega at the NQC point 
 *   @param[in, out]  domg_tmp: first time derivative of the frequency at the NQC point 
 *
 */
void eob_nqc_point_HM_peak22(Dynamics *dyn, double *A_tmp, double *dA_tmp, double *omg_tmp, double *domg_tmp)
{

  int knqcpeak22_size = EOBPars->knqcpeak22_size;  
  int *knqcpeak22     = EOBPars->knqcpeak22;

  const double nu   = EOBPars->nu;
  const double chi1 = EOBPars->chi1;
  const double chi2 = EOBPars->chi2;

  const double nu2  = SQ(nu);
  const double nu3  = nu2*nu;
  const double delta = sqrt(1-4.*nu);

  const double chi_S = (chi1+chi2)/2;
  const double chi_A = (chi1-chi2)/2;

  const double chi   = chi_S + chi_A*delta/(1.-2.*nu);
  const double chisq = chi*chi;
  const double chi3  = chisq*chi;
  const double chi4  = chi3*chi;
  const double chi_21A  = delta*chi_S/(1-1.3*nu) + chi_A;
  const double chi_21A2 = SQ(chi_21A);
  const double chi_21A3 = chi_21A2*chi_21A;
  const double chi_21D  = delta*chi_S/(1-2.*nu)  + chi_A;
  const double chi_33   = chi_S*delta+chi_A;
  const double chi_44A  = (1 - 5.*nu)*chi_S + chi_A*delta;
  const double chi_44D  = (1 - 7.*nu)*chi_S + chi_A*delta;


  for(int j=0; j<knqcpeak22_size; j++){
    int k = knqcpeak22[j];
    switch(k)
    {
      case(-1):
        if (VERBOSE) printf("eob_nqc_point_HM_peak22: no modes selected in knqcpeak22\n");
        break;
      case(0):
        /* l=2, m=1 */
        A_tmp[0]    = fabs(-0.033175*chi_21A3*delta + 0.086356*chi_21A2*delta*nu - 0.049897*chi_21A2*delta + 0.012706*chi_21A*delta + 0.168668*chi_21A*nu - 0.285597*chi_21A + 1.067921*delta*nu2 - 0.189346*delta*nu + 0.431426*delta);
        dA_tmp[0]   = chi_21D*delta*(0.023534*nu - 0.008064) + delta*(0.006743 - 0.0297*nu) + 0.008256*fabs( (chi_21D - delta*(5.471011*nu2 + 1.235589*nu + 0.815482)) );
        omg_tmp[0]  = 0.01009*chi3 - 0.077343*chisq*nu + 0.02411*chisq + 0.168854*chi*nu2 - 0.159382*chi*nu + 0.047635*chi + 1.965157*nu3 - 0.53085*nu2 + 0.237904*nu + 0.176526;
        domg_tmp[0] = 0.00149*chi3 + 0.008965*chisq*nu - 0.002739*chisq - 0.033831*chi*nu2 + 0.005752*chi*nu - 0.002003*chi + 0.204368*nu3 - 0.120705*nu2 + 0.035144*nu + 0.006579;
        break;
      case(4):
        /* l=3, m=3 */
        A_tmp[4]    = fabs(- 0.088371*SQ(chi_33)*delta*nu + 0.036258*SQ(chi_33)*delta + 1.057731*chi_33*nu2 - 0.466709*chi_33*nu + 0.099543*chi_33 + 1.96267*delta*nu2 + 0.027833*delta*nu + 0.558808*delta);
        dA_tmp[4]   = SQ(chi_33)*delta*(0.004941*nu - 0.002094) + 0.001781*sqrt(fabs( (SQ(chi_33) + chi_33*delta*(39.247538*nu - 2.986889) + SQ(delta)*(85.173306*nu + 4.637906)) ));
        omg_tmp[4]  = 0.045141*chi3 - 0.346675*chisq*nu + 0.119419*chisq + 0.745924*chi*nu2 - 0.478915*chi*nu + 0.17467*chi - 8.887163*nu3 + 4.226831*nu2 + 0.427167;
        domg_tmp[4] = 0.001697*chi3 + 0.016231*chisq*nu - 0.003985*chisq - 0.154378*chi*nu2 + 0.050618*chi*nu - 0.002721*chi - 0.255402*nu3 + 0.08663*nu2 + 0.027405*nu + 0.009736;
        break;
      case(8):
        /* l=4, m=4 */
        A_tmp[8]    = fabs( 0.031483*SQ(chi_44A) - 0.180165*chi_44A*nu + 0.063931*chi_44A + 6.239418*nu3 - 1.947473*nu2 - 0.615307*nu + 0.262533 );
        dA_tmp[8]   = -0.001251*SQ(chi_44D)*chi_44D + 0.006387*SQ(chi_44D)*nu - 0.001223*SQ(chi_44D) - 0.034308*chi_44D*nu2 + 0.014373*chi_44D*nu - 0.000681*chi_44D + 1.134679*nu3 - 0.417056*nu2 + 0.024004*nu + 0.003498;
        omg_tmp[8]  = 0.042529*chi3 - 0.415864*chisq*nu + 0.155222*chisq + 0.768712*chi*nu2 - 0.592568*chi*nu + 0.244508*chi - 13.651335*nu3 + 5.490329*nu2 + 0.574041;
        domg_tmp[8] = 0.001812*chi3 + 0.024687*chisq*nu - 0.00568*chisq  - 0.162693*chi*nu2 + 0.061205*chi*nu - 0.003623*chi - 0.536664*nu3  + 0.094797*nu2 + 0.045406*nu + 0.013038;
        break;
      default:
        errorexit("A,dA,omg,domg fits only implemented for (2,1), (3,3) and (4,4) at the moment.");
    }
    /* Normalize A and dA */
    int l       = LINDEX[k];
    double fact = sqrt((l+2)*(l+1)*l*(l-1));
    A_tmp[k]    = A_tmp[k]/fact;
    dA_tmp[k]   = dA_tmp[k]/fact;
  }
}

/** 
 *  Function: eob_nqc_point_postpeak
 *  --------------------------------
 *   This function computes the NQC functioning points by evaluating 
 *   the post-peak template and its derivatives.
 *   
 *   @param[in]  tau     : Evaluation time (in units of Mbh)
 *   @param[in]  Mbh     : Mass of the final BH
 *   @param[in]  c1A, c2A,c3A,c4A,c1phi.c2phi,c3phi,c4phi : postpeak template parameters
 *   @param[in]  alpha1,omega1 : QNM frequency and inverse damping time 
 *   @param[in, out] A_tmp   : RWZ-normalized Amplitude at the NQC point 
 *   @param[in, out] dA_tmp  : first time derivative of the RWZ-normalized amplitude at the NQC point 
 *   @param[in, out] omg_tmp : omega at the NQC point 
 *   @param[in, out] domg_tmp: first time derivative of the frequency at the NQC point 
 *
 */
void eob_nqc_point_postpeak(double tau, double Mbh, double c1A, double c2A, double c3A, double c4A, 
			    double c1phi, double c2phi, double c3phi, double c4phi,
			    double alpha1, double omega1,
			    double *A_tmp, double *dA_tmp, double *omg_tmp, double *domg_tmp)
{
  
  double Mbh2, x, x2, dA_tmp1, dA_tmp2, omg_tmp1, omg_tmp2, domg_tmp_n1,domg_tmp_n2,domg_tmp_d1,domg_tmp_d2;
  
  /* The time variable in the post-peak template is given in units of Mbh */
  Mbh2 = SQ(Mbh);
  
  *A_tmp = exp(-alpha1*tau)*(c1A*tanh(c2A*tau + c3A) + c4A);
  
  double fact = cosh(c3A+ c2A*tau);
  dA_tmp1     = c1A*c2A*exp(-alpha1*tau)/(fact*fact);
  dA_tmp2     = -alpha1*exp(-alpha1*tau)*(c4A + c1A*tanh(c3A+ c2A*tau));
  *dA_tmp     = (dA_tmp1 + dA_tmp2)/Mbh;
  
  x        = exp(-c2phi*tau);
  x2       = x*x;
  omg_tmp1 = c1phi*(-2.*c2phi*c4phi*exp(-2.*c2phi*tau) - c2phi*c3phi*x);
  omg_tmp2 = 1 + c4phi*exp(-2.*c2phi*tau) + c3phi*x;
  *omg_tmp = (omg_tmp1/omg_tmp2)/Mbh + omega1/Mbh;
	  
  domg_tmp_n1 = c1phi*SQ(-2.*c2phi*c4phi*x2 - c2phi*c3phi*x);
  domg_tmp_d1 = SQ(1 + c4phi*x2 + c3phi*x);
  domg_tmp_n2 = c1phi*(4.*SQ(c2phi)*c4phi*x2 + SQ(c2phi)*c3phi*x);
  domg_tmp_d2 = 1 + c4phi*x2 + c3phi*x;
  *domg_tmp   = -(domg_tmp_n1/domg_tmp_d1 - domg_tmp_n2/domg_tmp_d2)/Mbh2;
}
/** @} */ // end of nqcpoint

/** 
 *  Function: eob_nqc_timeshift
 *  ---------------------------
 *   This function computes the time-shift for (2,2) tNQC with respect
 *   to the peak of the EOB pure orbital frequency, from
 *   https://arxiv.org/abs/1806.01772
 *   
 *   @param[in] nu   : symmmetric mass ratio
 *   @param[in] chi1 : dimensionless spin of body 1
 *
 *   @return DeltaT_nqc
 */
double eob_nqc_timeshift_bbh(double nu, double chi1)
{

  double DeltaT_nqc = 1.0;  
  /* New Delta_T NQC determined in TEOBResumS paper (arXiv:1806.01772) */
  if (((chi1 < -0.9) && (nu < 8./81.)) || ((chi1 < -0.8) && (nu < 11./144.))) {
    DeltaT_nqc = 4.0;
  } else {
    DeltaT_nqc = 1.0; // standard choice inspired by test-particle results
  }

  if (EOBPars->d_delta_t_nqc) DeltaT_nqc += EOBPars->d_delta_t_nqc;
    
  return DeltaT_nqc;  
}

/** 
 *  Function: eob_nqc_deltat_lm
 *  ---------------------------
 *   This function computes the time-shift of the peak of the modes w.r.t
 *   the (2,2) mode, as defined in Eq. 3.50 of https://arxiv.org/1904.09550
 *   
 *   @param[in, out] Dt_lm : KMAX dimensional array of time-shifts 
 *
 */
void eob_nqc_deltat_lm_bbh(double *Dt_lm)
{
  
  const double nu    = EOBPars->nu;
  const double nu2   = SQ(nu);
  const double nu3   = nu2*nu;
  const double nu4   = SQ(nu2);
  const double X1    = EOBPars->X1;
  const double X2    = EOBPars->X2;
  const double X12   = X1 - X2;
  const double X12_2 = SQ(X12);
  const double X12_3 = X12_2*X12;

  const double chi1  = EOBPars->chi1;
  const double chi2  = EOBPars->chi2;
  const double aK    = X1*chi1 + X2*chi2;
  const double aK2   = SQ(aK);
  const double a12   = X1*chi1 - X2*chi2;
  const double a12_2 = SQ(a12);
  const double Shat  = 0.5*(aK + a12*X12);
  const double Shat2 = SQ(Shat);
  
  const int usespins = EOBPars->use_spins;  

  double Dt_lm_TP[KMAX];
  for (int k=0; k<KMAX; k++) {
    Dt_lm[k] = Dt_lm_TP[k] = 0;
  }

  Dt_lm_TP[0]  = 11.75925;
  Dt_lm_TP[2]  = 12.9338;
  Dt_lm_TP[3]  = 9.16665;
  Dt_lm_TP[4]  = 3.42593;
  Dt_lm_TP[5]  = 13.1116;
  Dt_lm_TP[6]  = 11.66665;
  Dt_lm_TP[7]  = 9.53705;
  Dt_lm_TP[8]  = 5.27778;
  Dt_lm_TP[13] = 6.561811;

  if (!(usespins)) {
    
    Dt_lm[0]  = Dt_lm_TP[0]*(1 - 3.764*nu + 6.9051*nu2);
    Dt_lm[2]  = Dt_lm_TP[2]*(1 - 25.615*nu2)/(1 + 0.88803*nu + 16.292*nu2);
    Dt_lm[3]  = Dt_lm_TP[3]*(1 - 11.398*nu + 33.244*nu2)/(1 - 8.1976*nu + 19.537*nu2);
    Dt_lm[4]  = Dt_lm_TP[4]*(1 - 0.11298*nu + 5.0056*nu2);
    Dt_lm[5]  = Dt_lm_TP[5]*(1 - 9.6225*nu + 38.451*nu2)/(1 - 7.7998*nu + 32.405*nu2);
    Dt_lm[6]  = Dt_lm_TP[6]*(1 - 3.8284*nu - 12.399*nu2);
    Dt_lm[7]  = Dt_lm_TP[7]*(1 - 11.345*nu + 38.813*nu2)/(1 - 7.5049*nu + 22.399*nu2);
    Dt_lm[8]  = Dt_lm_TP[8]*(1 - 8.4686*nu + 18.006*nu2)/(1 - 6.7964*nu + 11.368*nu2);
    //Dt_lm[13] = Dt_lm_TP[13]*(1 - 12.198*nu + 40.327*nu2)/(1 - 11.501*nu + 39.431*nu2);
    //Dt_lm[13] = Dt_lm_TP[13]*(1 - 91.4401*nu +2548.5975*nu2 -11086.4884*nu3 + 27137.0063*nu4)/(1 - 67.156*nu + 1773.5942*nu2);
    Dt_lm[13] = Dt_lm_TP[13]*(1 - 91.2039*nu + 2556.5123*nu2 - 11325.217*nu3 + 27767.2164*nu4)/(1 - 66.0362*nu + 1762.4169*nu2);
    
  } else {
    
    double Dt_21_orb, Dt_32_orb, Dt_33_orb, Dt_42_orb, Dt_43_orb, Dt_44_orb;
    double Dt_32_S, Dt_33_S1, Dt_33_S2, Dt_42_S, Dt_43_S, Dt_44_S;
    double c11Dtlm, c12Dtlm, c21Dtlm, c22Dtlm, b1Dtlm, b2Dtlm;
    double a1Dt, a2Dt, a3Dt, a4Dt;
    double b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15;
    double n1, n2, d1, d2;

    /* (l = 2, m = 1) */
    Dt_21_orb	= (Dt_lm_TP[0]*X12_2 + 4*6.6264*nu)*(1-2.0728*nu*X12);

    c11Dtlm   = -1976.13; 
    c12Dtlm   = 3719.88;
    c21Dtlm   = -2545.41; 
    c22Dtlm   = 5277.62;
    b1Dtlm    = 0.0472289;
    b2Dtlm    = 0.115583;
    b1        = (b1Dtlm + c11Dtlm*X12)/(1 + c12Dtlm*X12);
    b2        = (b2Dtlm + c21Dtlm*X12)/(1 + c22Dtlm*X12);
    Dt_lm[0]  = (1 + b1*aK + b2*aK2)*Dt_21_orb;

    /* (l = 3, m = 2) */
    a1Dt      = -11.3497;
    a2Dt      = 32.9144; 
    a3Dt      = -8.36579;
    a4Dt      = 20.1017;
    Dt_32_orb = Dt_lm_TP[3]*(1 + a1Dt*nu+ a2Dt*nu2)/(1 + a3Dt*nu+ a4Dt*nu2);

    if (nu < 0.2) {

      b1 = -0.037634;
      b2 = +12.456704;
      b3 = +2.670868;
      b4 = -12.255859;
      b5 = +37.843505;
      b6 = -25.058475;
      b7 = +449.470722;
      b8 = -1413.508735;
      b9 = -11.852596;
      b10 = +41.348059;
      b11 = -5.650710;
      b12 = -9.567484;
      b13 = +173.182999;
      b14 = -10.938605;
      b15 = +35.670656;

      Dt_32_S  = (1 + (-0.34161 + b1*nu + b2*nu2 + b3*nu3)/(1 + b4*nu + b5*nu2)*Shat + (-0.46107 + b6*nu + b7*nu2 + b8*nu3)/(1 + b9*nu + b10*nu2)*Shat2)/(1 + (+0.34744 + b11*nu + b12*nu2 + b13*nu3)/(1 + b14*nu + b15*nu2)*Shat);
      
    } else {
      
      b1 = +2.497188;
      b2 = -7.532596;
      b3 = +4.645986;
      b4 = -3.652524;
      b5 = +3.398687;
      b6 = +7.054185;
      b7 = -12.260185;
      b8 = +5.724802;
      b9 = -3.242611;
      b10 = +2.714232;
      b11 = +2.614565;
      b12 = -9.507583;
      b13 = +7.321586;
      b14 = -3.937568;
      b15 = +4.584970;
      Dt_32_S  = (1 + (+0.15477 + b1*X12 + b2*X12_2 + b3*X12_3)/(1 + b4*X12 + b5*X12_2)*Shat + (-0.755639 + b6*X12 + b7*X12_2 + b8*X12_3)/(1 + b9*X12 + b10*X12_2)*Shat2)/(1 + (+0.21816 + b11*X12 + b12*X12_2 + b13*X12_3)/(1 + b14*X12 + b15*X12_2)*Shat);

    }

    Dt_lm[3] = Dt_32_orb*Dt_32_S;

    /* (l = 3, m = 3) */
    Dt_33_orb = Dt_lm_TP[4]*(1 + 0.183349*nu+ 4.22361*nu2);
    Dt_33_S1  = 1 + (-0.49791 - 1.9478*nu)/(1 + 13.9828*nu)*Shat + (-0.18754 + 1.25084*nu)/(1 - 3.41811*nu)*Shat2;
    Dt_33_S2  = 1 + (-1.07291 - 1043.15*nu)/(1 + 1033.85*nu)*Shat;
    Dt_lm[4]  = Dt_33_orb*Dt_33_S1/Dt_33_S2; 

    /* (l = 4, m = 2) */
    Dt_42_orb  = Dt_lm_TP[6]*(1 - 9.8446172795*nu + 23.3229430582*nu2)/(1 - 5.7604819848*nu + 7.1217930024*nu2);
    
    if (nu < 6./25.) {
      b1       = 0.;
      b2       = 0.;
      b3       = 0.;
      b4       = 0.;
    } else {
      b1       = +24.604717;
      b2       = -0.808279;
      b3       = +62.471781;
      b4       = +48.340961;
    }
    Dt_42_S  = (1 + (-1.3002045 + b1*nu)/(1 + b2*nu)*Shat)/(1 + (-0.9494348+b3*nu)/(1+b4*nu)*Shat);
    Dt_lm[6] = Dt_42_orb*Dt_42_S;

    /* (l = 4, m = 3) */
    a1Dt      = -11.2377;
    a2Dt      = 38.3177; 
    a3Dt      = -7.29734;
    a4Dt      = 21.4267;
    Dt_43_orb = Dt_lm_TP[7]*(1 + a1Dt*nu + a2Dt*nu2)/(1 + a3Dt*nu + a4Dt*nu2);

    b1      = +3.215984;
    b2      = +42.133767;
    b3      = -9.440398;
    b4      = +35.160776;
    b5      = +1.133942;
    b6      = -10.356311;
    b7      = -6.701429;
    b8      = +10.726960;
    b9      = -6.036207;
    b10     = +67.730599;
    b11     = -3.082275;
    b12     = +11.547917;
    Dt_43_S = (1 + (-1.371832 + b1*nu + b2*nu2)/(1 + b3*nu + b4*nu2)*Shat +((+0.362375+b5*nu+b6*nu2)/(1+b7*nu+b8*nu2))*Shat2)/(1+((-1.0808402+b9*nu+b10*nu2)/(1+b11*nu+b12*nu2))*Shat);
    Dt_lm[7] 		= Dt_43_orb*Dt_43_S;

    /* (l = 4, m = 4) */
    n1        = -8.35574;
    n2        = 17.5288;
    d1        = -6.50259;
    d2        = 10.1575;
    Dt_44_orb = Dt_lm_TP[8]*(1 + n1*nu + n2*nu2)/(1 + d1*nu + d2*nu2);
    b1        = 0.00159701 - 2.28656*X12 + 1.66532*X12_2;
    b2        = -1.14134 - 0.589331*X12 + 0.708784*X12_2;
    Dt_44_S   = (1 + b1*Shat)/(1 + b2*Shat);
    Dt_lm[8]  = Dt_44_orb*Dt_44_S;

    // These fits are taken from the testparticle limit and approixmate the spinning case
    /* (l = 3, m = 1) */
    Dt_lm[2]  = Dt_lm_TP[2]*(1 - 25.615*nu2)/(1 + 0.88803*nu + 16.292*nu2);

    /* (l = 4, m = 1) */	
    Dt_lm[5]  = Dt_lm_TP[5]*(1 - 9.6225*nu + 38.451*nu2)/(1 - 7.7998*nu + 32.405*nu2);	

    /* (l = 5, m = 5) */
    //Dt_lm[13] = Dt_lm_TP[13]*(1 - 12.198*nu + 40.327*nu2)/(1 - 11.501*nu + 39.431*nu2);
    //Dt_lm[13] = Dt_lm_TP[13]*(1 - 91.4401*nu +2548.5975*nu2 -11086.4884*nu3 + 27137.0063*nu4)/(1 - 67.156*nu + 1773.5942*nu2);
    Dt_lm[13] = Dt_lm_TP[13]*(1 - 91.2039*nu + 2556.5123*nu2 - 11325.217*nu3 + 27767.2164*nu4)/(1 - 66.0362*nu + 1762.4169*nu2);

  }

}

/** 
 *  Function: eob_nqc_setcoefs
 *  --------------------------
 *   This function sets the waveform and flux NQC coefficients according
 *   to the flags NQC_FLX_* and NQC_HLM_*
 *
 *   @param[in, out] nqc   : NQC coefficients
 *
 */
void eob_nqc_setcoefs(NQCdata *nqc)
{
  nqc->flx->add = 1;
  nqc->hlm->add = 1;

  if (EOBPars->nqc_coefs_flx == NQC_FLX_NONE) nqc->flx->add = 0;
  if (EOBPars->nqc_coefs_hlm == NQC_HLM_NONE) nqc->hlm->add = 0;
      
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

  if (EOBPars->nqc_coefs_flx == NQC_FLX_NRFIT_NOSPIN_201602)
    eob_nqc_setcoefs_nospin201602(nqc->flx);
  else if (EOBPars->nqc_coefs_flx == NQC_FLX_NRFIT_SPIN_202002)
    eob_nqc_setcoefs_spin202002(nqc->flx);
  else if (EOBPars->nqc_coefs_flx == NQC_FLX_FROMFILE) 
    eob_nqc_setcoefs_fromfile(nqc->flx, EOBPars->nqc_coefs_flx_file);

  if (EOBPars->nqc_coefs_hlm == NQC_HLM_NRFIT_NOSPIN_201602) 
    eob_nqc_setcoefs_nospin201602(nqc->hlm);
  else if (EOBPars->nqc_coefs_hlm == NQC_HLM_NRFIT_SPIN_202002)
    eob_nqc_setcoefs_spin202002(nqc->hlm);
  else if (EOBPars->nqc_coefs_hlm == NQC_HLM_FROMFILE) 
    eob_nqc_setcoefs_fromfile(nqc->hlm,EOBPars->nqc_coefs_hlm_file);
  
}

/** 
 *  Function: eob_nqc_setcoefs_nospin201602
 *  ---------------------------------------
 *  Set NQC coefficients. NR fits for nonspinning case 2016/02/09 
 *  Hardcoded in eob_wav_hlmNQC_nospin201602()
 *  https://arxiv.org/pdf/1703.06814.pdf
 *
 *   @param[in, out] nqc   : NQC coefficients
 *
 */
void eob_nqc_setcoefs_nospin201602(NQCcoefs *nqc)
{

  const double nu = EOBPars->nu;  
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

/** 
 *  Function: eob_nqc_setcoefs_fromfile
 *  -----------------------------------
 *  Set NQC coefficients from file 
 *
 *   @param[in, out] nqc  : NQC coefficients
 *   @param[in] fname : file name from which to read NQCs
 *
 */
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


/** 
 *  Function: eob_nqc_setcoefs_spin202002
 *  -------------------------------------
 *  Set NQC coefficients. NR fits for spinning case 2020/02
 *  https://arxiv.org/pdf/2104.07533.pdf
 *
 *   @param[in, out] nqc   : NQC coefficients
 *
 */
void eob_nqc_setcoefs_spin202002(NQCcoefs *nqc)
{
  const int k22 = 1;
  const double nu = EOBPars->nu;
  const double chi1 = EOBPars->chi1;  
  const double chi2 = EOBPars->chi2;  
  nqc->activemode[k22] = 1;
	  
  double a1fit = get_a1_fit_22(nu, chi1, chi2);
  double a2fit = get_a2_fit_22(nu, chi1, chi2);

  /* Limiting fits to calibrated and non-singular region */
  double X1   = 0.5*(1. + sqrt(1.-4*nu));
  double X2   = 1.-X1;
  double Shat_norm = (X1*X1*chi1 + X2*X2*chi2)/(X1*X1 + X2*X2);

  if (nu < 0.025) {
    a1fit = 0.;
    a2fit = 0.;
  } else if ((nu < 0.09) && (Shat_norm > 0.83)) {
    a1fit = 0.;
    a2fit = 0.;
  } else if ((nu < 0.13) && (Shat_norm > 0.87)) {
    a1fit = 0.;
    a2fit = 0.;
  }
  
  /* (2,2) */
  nqc->a1[k22] = a1fit;
  nqc->a2[k22] = a2fit;
  nqc->a3[k22] = 0.0;
	    
  nqc->b1[k22] = 0.0;
  nqc->b2[k22] = 0.0;
  nqc->b3[k22] = 0.0;
	
  nqc->add = 1;
  nqc->maxk = k22;
  
}

/** 
 *  Function: get_a1_fit_22
 *  -----------------------
 *  Compute the a1 coefficient for the (2,2) mode.
 *  Fit from https://arxiv.org/pdf/2104.07533.pdf
 *
 *   @param[in] nu   : symmetric mass ratio
 *   @param[in] chi1 : dimensionless spin of body 1
 *   @param[in] chi2 : dimensionless spin of body 2
 *
 *   @return a1
 */	
double get_a1_fit_22(double nu, double chi1, double chi2)
{
  const int usespins = EOBPars->use_spins;
  double res = 0.;

  double X12  = sqrt(1.-4*nu);
  double X1   = 0.5*(1.+X12);
  double X2   = 1.-X1;
  double Shat = (X1*X1*chi1 + X2*X2*chi2);
  double Shat_norm = (X1*X1*chi1 + X2*X2*chi2)/(X1*X1 + X2*X2);
    
  if (!(usespins)) 
    {
      double a1_q1 = 0.070974;
      double b1    = 0.786350;
      double b2    = -9.085105;
      double X     = 1.-4.*nu;
      res          = a1_q1/(1. + b1*X + b2*X*X);
    }
  else if (nu>0.2485)
    {
      double b1    = 0.121187;
      double b2    = -5.950663;
      double b3    = 9.420324;
      double b4    = -10.601339;
      double b5    = 17.641549;
      double b6    = -5.684777;
      double b7    = 10.910451;
      double b8    = -6.867377;
      res      = b1*(1. + b2*Shat + b3*Shat*Shat+ b4*Shat*Shat*Shat+ b5*Shat*Shat*Shat*Shat)/(1. + b6*Shat+ b7*Shat*Shat+ b8*Shat*Shat*Shat);
    }
  else if (nu>0.16)
    {
      double b1    = 0.26132647;
      double b2    = -4.90302367;
      double b3    = 20.67036124;
      double b4    		= -3.17109808;
      double a1_nospin   = b1*(1. + b2*nu + b3*nu*nu*nu)/(1. + b4*nu);
      
      double c1    = -3.082861;
      double c2    = 2.169948;
      double c3    = -0.636353;
      double c4    = 0.741419;
      double c5    = -2.843896;
      double c6    = 2.709697;
      double c7    = -0.832894;
      double a1_spin = (1. + c1*Shat_norm + c2*Shat_norm*Shat_norm+ c3*Shat_norm*Shat_norm*Shat_norm+ c4*Shat_norm*Shat_norm*Shat_norm*Shat_norm)/(1. + c5*Shat_norm+ c6*Shat_norm*Shat_norm+ c7*Shat_norm*Shat_norm*Shat_norm);
      res      = a1_nospin*a1_spin;
    }
  else
    {
      double b1    = 0.341803;
      double b2    = -1.350488;
      double b3    = -6.353357;
      double b4    = 2.216156;
      double a1_nospin   = b1*(1. + b2*nu + b3*nu*nu*nu)/(1 + b4*nu);
	
      double c1    = 7.650946;
      double c2    = 7.106992;
      double c3    = -60.630748;
      double c4    = -69.630357;
      double c5    = 47.114247;
      double c6    = 5.733002;
      double c7    = -12.905707;
      double c8    = 5.045688;
      double c9    = 3.515869;
      double c10    = 1.564146;
      double c11    = 0.642864;
      double c12    = 2.947890;
      double c13    = 31.023038;
      double c14    = 1.829543;
      double xnu    = nu-0.16;
      double a1_spin  =  (1.-2.287721*(1.+c1*xnu)/(1.+c2*xnu)*Shat_norm-0.598451*(1.+c3*xnu)/(1.+c4*xnu)*Shat_norm*Shat_norm+0.766069*(1.+c5*xnu)/(1.+c6*xnu)*Shat_norm*Shat_norm*Shat_norm+1.857169*(1.+c7*xnu)/(1.+c8*xnu)*Shat_norm*Shat_norm*Shat_norm*Shat_norm)/(1.-2.035234*(1.+c9*xnu)/(1.+c10*xnu)*Shat_norm+0.836427*(1.+c11*xnu)/(1.+c12*xnu)*Shat_norm*Shat_norm+0.297476*(1.+c13*xnu)/(1.+c14*xnu)*Shat_norm*Shat_norm*Shat_norm);
      res  = a1_nospin*a1_spin;
    }

  return res;
}

/** 
 *  Function: get_a2_fit_22
 *  -----------------------
 *  Compute the a2 coefficient for the (2,2) mode.
 *  Fit from https://arxiv.org/pdf/2104.07533.pdf
 *
 *   @param[in] nu   : symmetric mass ratio
 *   @param[in] chi1 : dimensionless spin of body 1
 *   @param[in] chi2 : dimensionless spin of body 2
 *
 *   @return a2
 */
double get_a2_fit_22(double nu, double chi1, double chi2)
{
  const int usespins = EOBPars->use_spins;
  double res = 0.;

  double X12  = sqrt(1.-4*nu);
  double X1   = 0.5*(1.+X12);
  double X2   = 1.-X1;
  double Shat = (X1*X1*chi1 + X2*X2*chi2);
  double Shat_norm = (X1*X1*chi1 + X2*X2*chi2)/(X1*X1 + X2*X2);
    
  if (!(usespins)) 
    {
      double x     	= sqrt(1.-4.*nu);
      double a2q1 	= 1.315133;
      double b1    	= -0.324849;
      double b2    	= -0.304506;
      double b3    	= -0.371614;
      res      		= a2q1*(1. + b1*x + b2*x*x)/(1. + b3*x);
    }
  else if (nu>0.2485)
    {
      double b1   = 1.331703;
      double b2   = -4.237724;
      double b3   = 1.786023;
      double b4   = 10.546205;
      double b5   = -9.698233;
      double b6   = -6.225823;
      double b7   = 13.209381;
      double b8   = -9.402513;
      res      	= b1*(1. + b2*Shat + b3*Shat*Shat+ b4*Shat*Shat*Shat+ b5*Shat*Shat*Shat*Shat)/(1. + b6*Shat+ b7*Shat*Shat+ b8*Shat*Shat*Shat);
    }
  else if (nu>0.16)
    {
      double b1    		= 1.03364144;
      double b2    		= -3.46191440;
      double b3    		= -7.86652243;
      double b4    		= -3.96268815;
      double a2_nospin 	= b1*(1. + b2*nu + b3*nu*nu*nu)/(1. + b4*nu);
	
      double c1    		= 0.036452;
      double c2    		= -64.360789;
      double c3    		= 0.275707;
      double c4    		= -34.573145;
      double c5    		= -0.113951;
      double c6    		= -2.531304;
      double c7    		= -7.691661;
      double c8    		= -1.025824;
      double c9    		= 4.237539;
      double c10    		= 0.593579;
      double c11    		= 1.661809;
      double c12    		= -0.939736;
      double c13    		= -6.333442;
      double a2_spin 		= (1. + c1*(1.+c2*nu)*Shat_norm + c3*(1.+c4*nu)*Shat_norm*Shat_norm + c5*Shat_norm*Shat_norm*Shat_norm+ c6*(1.+c7*nu)*Shat_norm*Shat_norm*Shat_norm*Shat_norm)/(1. + c8*(1.+c9*nu)*Shat_norm+ c10*(1.+c11*nu)*Shat_norm*Shat_norm + c12*(1.+c13*nu)*Shat_norm*Shat_norm*Shat_norm);
      res 	= a2_nospin*a2_spin;
    }
  else
    {
      double b1    = 0.929192;
      double b2    = 1.334263;
      double b3    = -26.389790;
      double b4    = -1.289984;
      double a2_nospin  = b1*(1. + b2*nu + b3*nu*nu*nu)/(1. + b4*nu);
      double c1    = 15.871482;
      double c2    = 5.066190;
      double c3    = 7.168498;
      double c4    = 6.709490;
      double c5    = 18.583382;
      double c6    = 5.764512;
      double c7    = -14.038564;
      double c8    = -17.126231;
      double c9    = 6.387917;
      double c10    = 3.438456;
      double c11    = 8.867098;
      double c12    = 2.910938;
      double xnu    = nu-0.16;
      double a2_spin = (1.-0.886561*(1+c1*xnu)/(1+c2*xnu)*Shat_norm-1.953955*(1+c3*xnu)/(1+c4*xnu)*Shat_norm*Shat_norm+1.366537*(1+c5*xnu)/(1+c6*xnu)*Shat_norm*Shat_norm*Shat_norm+0.950212*(1+c7*xnu)/(1+c8*xnu)*Shat_norm*Shat_norm*Shat_norm*Shat_norm)/(1.-2.531000*(1+c9*xnu)/(1+c10*xnu)*Shat_norm+1.723991*(1+c11*xnu)/(1+c12*xnu)*Shat_norm*Shat_norm);
	
      res      = a2_nospin*a2_spin;
    }

  return res;
}

/** 
 *  Function: YagiYunes13_fit_logQ_coefs
 *  ------------------------------------
 *  This function sets the coefficients for the 
 *  Yunes-Yagi logQ-vs-log(lambda) fit 
 *  Table I of https://arxiv.org/abs/1303.1528 
 *  Here x = log(lambda) and the output is the log of the coefficient
 *  that describes the quadrupole deformation due to spin.
 *
 *   @param[in, out] c   : coefficients
 *
 */
void YagiYunes13_fit_logQ_coefs(double *c)
{
  c[0] = 0.194;
  c[1] = 0.0936;
  c[2] = 0.0474;
  c[3] = -4.21e-3;
  c[4] = 1.23e-4;
}

/** 
 *  Function: YagiYunes13_fit_logQ
 *  ------------------------------------
 *  Yunes-Yagi logQ-vs-log(lambda) fit 
 *  Table I of https://arxiv.org/abs/1303.1528 
 *  Here x = log(lambda) and the output is the log of the coefficient
 *  that describes the quadrupole deformation due to spin.
 *
 *   @param[in] x   : log(Lambda2)
 *
 *   @return logQ
 */
double YagiYunes13_fit_logQ(double x)
{
  double c[5];
  YagiYunes13_fit_logQ_coefs(c);
  const double x2 = x*x;
  const double x3 = x*x2;
  const double x4 = x*x3;
  return( c[4]*x4 + c[3]*x3 + c[2]*x2 + c[1]*x + c[0] );
}

/** 
 *  Function: YagiYunes13_fit_logQ_drvts
 *  ------------------------------------
 *  Compute Yunes-Yagi logQ-vs-log(lambda) fit, as well as 
 *  first and second derivatives of Q w.r.t. u
 *  Table I of https://arxiv.org/abs/1303.1528 
 *  Here x = log(lambda) and the output is the log of the coefficient
 *  that describes the quadrupole deformation due to spin and its first
 *  and second derivatives.
 *
 *   @param[in] Lam    : Lambda_2
 *   @param[in] Lam_u  : dLam/du   (dyn tides)
 *   @param[in] Lam_uu : d2Lam/d2u (dyn tides)
 *   @param[out] Q      : Quadrupole deformation
 *   @param[out] Q_u    : first derivative of Q
 *   @param[out] Q_uu   : second derivative of Q
 *
 */
void YagiYunes13_fit_logQ_drvts(double Lam, double Lam_u, double Lam_uu,
				double *Q, double *Q_u, double *Q_uu)
{
  double c[5];
  YagiYunes13_fit_logQ_coefs(c);
  const double div_Lam = 1./Lam;
  const double x = log(Lam);
  const double x_u = Lam_u * div_Lam;
  const double x_uu = ( -SQ(Lam_u) + Lam * Lam_uu ) * SQ(div_Lam);
  const double x2 = x*x;
  const double x3 = x*x2;
  const double x4 = x*x3;
  const double logQ = c[4]*x4 + c[3]*x3 + c[2]*x2 + c[1]*x + c[0];
  const double logQ_x = 4*c[4]*x3 + 3*c[3]*x2 + 2*c[2]*x + c[1];  
  const double logQ_xx = 12*c[4]*x2 + 6*c[3]*x + 2*c[2];  
  const double logQ_u = logQ_x * x_u;
  const double logQ_uu = logQ_xx * SQ(x_u) + logQ_x * x_uu;    
  *Q = exp(logQ);
  *Q_u = *Q * logQ_u;
  *Q_uu = *Q * ( SQ(logQ_u) + logQ_uu );
}

/** 
 *  Function: Yagi13_fit_barlamdel
 *  ------------------------------
 *  Yagi 2013 fits for NS multipolar
 *  $\bar{\lambda}_\ell$ = 2 k_\ell/(C^{2\ell+1} (2\ell-1)!!)$
 *  Eq.(9,10),(61); Tab.I; Fig.8 http://arxiv.org/abs/1311.087
 *
 *   @param[in] barlam2 : Lambda_2
 *   @param[in] ell     : multipolar index
 *
 *   @return Lambda_ell
 */
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

/** 
 *  Function: Yagi13_fit_barsigmalambda
 *  -----------------------------------
 *   Yagi 2013 fits for NS multipolar
 *   $\bar{\sigma_2}( \bar{\lambda}_2 )$
 *   Eq.(9,11),(61); Tab.I; Fig.9 http://arxiv.org/abs/1311.0872 
 *   See also later erratum
 *
 *   @param[in] barlam2 : Lambda_2
 *
 *   @return Sigma_2
 */
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

/** 
 *  Function: Yagi14_fit_Coct
 *  -------------------------
 *  Yagi et al. fits for octupolar spin-induced multipole C_Oct
 *  Eq. (90) and Table I of https://arxiv.org/abs/1403.6243
 *
 *   @param[in] C_Q : Quadrupolar spin-induced parameter
 *
 *   @return C_Oct
 */
double Yagi14_fit_Coct(double C_Q)
{
  const double A0  = -0.925;
  const double B1  =  1.98;
  const double nu1 =  0.273;
  const double cubrootCoct = A0 + B1*pow(C_Q,nu1);
  return cubrootCoct*cubrootCoct*cubrootCoct;
}

/** 
 *  Function: Yagi14_fit_Chex
 *  -------------------------
 *  Yagi et al. fits for hexapolar spin-induced multipole C_Hex
 *  Eq. (90) and Table I of https://arxiv.org/abs/1403.6243
 *
 *   @param[in] C_Q : Quadrupolar spin-induced parameter
 *
 *   @return C_Hex
 */
double Yagi14_fit_Chex(double C_Q)
{
  const double A0  = -0.413;
  const double B1  =  1.5;
  const double nu1 =  0.466;
  const double fourthrootChex = A0 + B1*pow(C_Q,nu1);
  return SQ(SQ(fourthrootChex));
}

/** 
 *  Function: JFAPG_fit_Sigma_Irrotational
 *  --------------------------------------
 *  
 *
 *   @param barlam2 : Lambda_2
 *
 *   @return Sigma_2
 */
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
  return - exp(lny);
}

/** 
 *  Function: JFAPG_fit_Sigma_Static
 *  --------------------------------
 *  
 *
 *   @param barlam2 : Lambda_2
 *
 *   @return Sigma_2
 */
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

/** 
 *  Function: Godzieba20_fit_barlamdel
 *  ----------------------------------
 *  Godzieba 2020 fits for NS multipolar
 *  $\bar{\lambda}_\ell$ = 2 k_\ell/(C^{2\ell+1} (2\ell-1)!!)$
 *  Eq.(4); Tab.I; https://journals.aps.org/prd/pdf/10.1103/PhysRevD.103.063036 
 *
 *   @param[in] barlam2 : Lambda_2
 *   @param[in] ell     : multipolar index (9>ell>2)
 *
 *   @return Lambda_ell
 */
double  Godzieba20_fit_barlamdel(double barlam2, int ell)
{  
  if (barlam2<=0.) return 0.;
  double lnx = log(barlam2);
  double *coef;
  double c23[7] = {-1.052, 1.165, 6.369e-3, 5.058e-3, -7.268e-4,
		     3.749e-5, -6.803e-8};
  double c24[7] = {-2.262, 1.383, 1.662e-3, 1.225e-2, -1.752e-3,
		     9.667e-5, -1.886e-6};
  double c25[7] = {-4.511, 2.382, -0.3182, 0.08999, -0.01142,
		     6.896e-4, -1.606e-05};
  double c26[7] = {-7.662, 4.949, -1.644,  0.4715, -0.06856,
		     0.004986, -1.44e-04};
  double c27[7] = {-13.07, 9.799, -4.197, 1.242, -0.1966,
         0.016,  -0.0005265};
  double c28[7] = {-26.49, 30.03, -18.38, 6.54 , -1.291,
         0.1346, -0.005795};
  switch(ell)
  {
    case(3):
      coef = c23;
      break;
    case(4):
      coef = c24;
      break;
    case(5):
      coef = c25;
      break;
    case(6):
      coef = c26;
      break;
    case(7):
      coef = c27;
      break;
    case(8):
      coef = c28;
      break;
    default:
      errorexit("Godzieba fits are for ell=3,4,5,6,7,8.");
      break;
  }
  double lny = coef[0];
  double lnxp = lnx;
  for (int i=1; i<=6; i++) {
    lny += lnxp * coef[i];
    lnxp *= lnx;
  }
  return exp(lny);
}

/** 
 *  Function: Chang14_fit_omegaf
 *  ----------------------------
 *  Chang+ 2014 fits for f-mode frequency vs Lambda
 *  Eq.(3.5) and Tab.I of https://arxiv.org/abs/1408.3789
 *
 *   @param[in] lam  : Lambda_2
 *   @param[in] ell  : multipolar index (6>ell>1)
 *
 *   @return omegaf_ell
 */
double Chang14_fit_omegaf(double lam, int ell)
{  
  if (lam<=0.) return 0.;
  double lnx = log(lam);
  double *coef;
  double a2[5] = {+1.820*1e-1, -6.836*1e-3, -4.196*1e-3, +5.215*1e-4, -1.857*1e-5};
  double a3[5] = {+2.245*1e-1, -1.500*1e-2, -1.412*1e-3, +1.832*1e-4, -5.561*1e-6};
  double a4[5] = {+2.501*1e-1, -1.646*1e-2, -5.897*1e-4, +8.695*1e-5, -2.368*1e-6};
  double a5[5] = {+2.681*1e-1, -1.638*1e-2, -2.497*1e-4, +4.712*1e-5, -1.166*1e-6};
  switch(ell) {
    case(2):
      coef = a2;
      break;
    case(3):
      coef = a3;
      break;
    case(4):
      coef = a4;
      break;
    case(5):
      coef = a5;
      break;
    default:
      errorexit("Chang fits are for ell=2,3,4,5");
      break;
  }
  double y = coef[0];
  double lnxp = lnx;
  for (int i=1; i<=4; i++) {
    y += lnxp * coef[i];
    lnxp *= lnx;
  }
  return y;
}

/** 
 *  Function: HealyBBHFitRemnant
 *  ----------------------------
 *  Mass and angular momentum of the final black hole
 *  Healey, Lousto and Zochlower (HLZ),
 *  arXiv: 1406.7295, published as PRD 90, 104004 (2014)
 *  @note the formula uses the convention that M2 > M1, so that
 *  chi2 should refer to the black hole with the largest
 *  mass. In the EOB code, this is given by chi1, since
 *  in EOB code we use the convention that M1 > M2
 *  
 *  Here it is q=M2/M1, with M2>M1
 *  
 *  Improved with (Eisco, Jisco) + iterative procedure 23/02/2016
 *  parameters (TABLE VI)
 *
 *   
 *   @param[in] chi1  : spin of body 1
 *   @param[in] chi2  : spin of body 2
 *   @param[in] q     : mass ratio, m2/m1 > 1
 *   @param[out] mass  : remnant mass 
 *   @param[out] spin  : remanant spin 
 *
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
      
  /* Parameters */
  double nu      = q/((1.+q)*(1.+q));
  
  /* Masses: convention here is that m2>m1 */
  double X2      = 0.5*(1.+sqrt(1.-4*nu));
  double X1      = 1.-X2;
  
  /* Spin variables */
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
  
  /* Mass ratio variables */
  double deltam  = -sqrt(1-4*nu); // X1 - X2
  double deltam2 =  deltam*deltam;
  double deltam3 =  deltam*deltam2;
  double deltam4 =  deltam*deltam3;
  double deltam6 =  deltam2*deltam4;
  
  /* Initialize the angular momentum */
  double a0 = s1 + s2;
  int a0_sign = 0.;
  
  if (a0==0) {
    a0_sign=0;
  } else if (a0>0) {
    a0_sign=1;
  } else { // if (a0<0) {
    a0_sign=-1;
  }
  
  /* Set-up an interative procedure to compute properly the "isco" quantities */
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
    
    /* Dimensionless spin: J/Mbh^2 */
    abh = (4*nu)*(4*nu)*(L0 + L1*S + L2a*Delta*deltam + L2b*S2 + L2c*Delta2 + L2d*deltam2 + L3a*Delta*S*deltam + L3b*S*Delta2 + L3c*S3 + L3d*S*deltam2 + L4a*Delta*S2*deltam + L4b*Delta3*deltam + L4c*Delta4 + L4d*S4 + L4e*Delta2*S2 + L4f*deltam4 + L4g*Delta*deltam3 + L4h*Delta2*deltam2 + L4i*S2*deltam2) + S*(1+8*nu)*deltam4 + nu*Jisco*deltam6;
    
    Mbh = (4*nu)*(4*nu)*(M0 + K1*S + K2a*Delta*deltam + K2b*S2 + K2c*Delta2 + K2d*deltam2 + K3a*Delta*S*deltam + K3b*S*Delta2 + K3c*S3 + K3d*S*deltam2 + K4a*Delta*S2*deltam + K4b*Delta3*deltam + K4c*Delta4 + K4d*S4 + K4e*Delta2*S2 + K4f*deltam4 + K4g*Delta*deltam3 + K4h*Delta2*deltam2 + K4i*S2*deltam2) + (1 + nu*(Eisco + 11))*deltam6;
    
    a0 = abh;
    
  }

  if (mass) *mass = Mbh;
  if (spin) *spin = abh;
}

/** 
 *  Function: JimenezFortezaRemnantMass
 *  ----------------------------
 *  Final mass fit of Jimenez-Forteza et al. (arxiv 1611.00332)
 *
 *   @param[in] nu    : symmetric mass ratio
 *   @param[in] X1    : mass fraction of body 1
 *   @param[in] X2    : mass fraction of body 2
 *   @param[in] chi1  : spin of body 1
 *   @param[in] chi2  : spin of body 2
 *
 *   @return Mbh, mass of the remnant BH
 */
double JimenezFortezaRemnantMass(double nu, double X1, double X2, double chi1, double chi2)
{
  const double xnu     = sqrt(1.0-4.0*nu);
  const double Dchi    = chi1-chi2;
  const double S       = (X1*X1*chi1+X2*X2*chi2)/(X1*X1+X2*X2);
	
  double S2 = S*S;
  double S3 = S2*S;
	
  double Dchi2 = Dchi*Dchi;
	
  double nu2 = nu*nu;
  double nu3 = nu2*nu;
  double nu4 = nu3*nu;
	  
  /* Coeffcients taken from Table VII, page 13: */ 
  const double a2      = 0.5610;
  const double a3      = -0.847;
  const double a4      = 3.145;
	
  /* The functional form is taken from eq. (21), page 14. */
  double E_rad_spin_zero = a4*nu4 + a3*nu3 + a2*nu2 + (1-2.*sqrt(2)/3.)*nu;
	
  /* Coeffcients taken from Table VIII, page 13: */
  const double b1      = -0.209;
  const double b2      = -0.197;
  const double b3      = -0.159;
  const double b5      = 2.985;
	 
  /* Coeffcients taken from Table IX, page 14: */
  const double f20     = 4.27;
  const double f30     = 31.09;
  const double f50     = 1.56735;
  const double f10     = -0.574752*f20 - 0.280958*f30 + 64.6408*f50 - 88.3165; //1.81
	
  /* Coeffcients taken from Table X, page 17: */
  const double f11     = 15.7;
  const double f21     = 0.;
  const double f31     = -243.6;
  const double f51     = -0.58;
	  
  /* Coeffcients evaluated from eq.(24), page 15: */
  const double f12     = 16. - 16.*f10 - 4.*f11;
  const double f22     = 16. - 16.*f20 - 4.*f21;
  const double f32     = 16. - 16.*f30 - 4.*f31;
  const double f52     = 16. - 16.*f50 - 4.*f51;
	
  /* this transformation is given in eq. (9), page 7: */
  const double b1tot      = b1*(f10 + f11*nu + f12*nu2);
  const double b2tot      = b2*(f20 + f21*nu + f22*nu2);
  const double b3tot      = b3*(f30 + f31*nu + f32*nu2);
  const double b5tot      = b5*(f50 + f51*nu + f52*nu2);
	  
  /* The functional form is taken from eq. (22), page 14: */
  double E_rad_eq_spin = 0.0484161*(0.128*b3tot*S3 + 0.211*b2tot*S2 + 0.346*b1tot*S + 1)/(1 - 0.212*b5tot*S);
	
  double E_rad_eq = 0.0484161;
	  
  /* These values are taken from Table IV, page 10: */
  const double d10     = -0.098;
  const double d11     = -3.23;
  const double d20     = 0.0112;
  const double d30     = -0.0198;
  const double d31     = -4.92;
	  
  /* The functional form is taken from eq. (27a-c), page 10: */
  const double A1      = d10*sqrt(1 - 4.*nu)*nu2*(d11*nu + 1);
  const double A2      = d20*nu3;
  const double A3      = d30*sqrt(1-4.*nu)*nu*(d31*nu+1);
	  
  /* The functional form is taken from eq. (15), page 9: */
  double E_rad_uneq_mass = A1*Dchi + A2*Dchi2 + A3*S*Dchi; 
	
  return 1 - (E_rad_spin_zero*E_rad_eq_spin/E_rad_eq + E_rad_uneq_mass);
}

/** 
 *  Function: JimenezFortezaRemnantSpin
 *  -----------------------------------
 *  Final spin fit of Jimenez-Forteza et al. (arxiv 1611.00332)
 *
 *   @param[in] nu    : symmetric mass ratio
 *   @param[in] X1    : mass fraction of body 1
 *   @param[in] X2    : mass fraction of body 2
 *   @param[in] chi1  : spin of body 1
 *   @param[in] chi2  : spin of body 2
 *
 *   @return abh, spin of the remnant BH
 */
double JimenezFortezaRemnantSpin(double nu, double X1, double X2, double chi1, double chi2)
{

  const double xnu     = sqrt(1.0-4.0*nu);
  const double Dchi    = chi1-chi2;
  const double S       = (X1*X1*chi1+X2*X2*chi2)/(X1*X1+X2*X2);
  const double a2      = 3.833;
  const double a3      = -9.49;
  const double a5      = 2.513;
  
  /* The functional form is taken from eq. (7), page 5. */
  double Lorb_spin_zero  = (1.3*a3*nu*nu*nu + 5.24*a2*nu*nu + 2.*sqrt(3)*nu)/(2.88*a5*nu + 1);
  
  /* Coeffcients taken from Table II, page 6: */
  double b1      = 1.00096;
  double b2      = 0.788;
  double b3      = 0.654;
  double b5      = 0.840;
  
  /* These values are taken from Table III, page 7: */
  double f21     = 8.774;
  double f31     = 22.83;
  double f50     = 1.8805;
  double f11     = 0.345225*f21 + 0.0321306*f31 - 3.66556*f50 + 7.5397;
  
  /* These values are taken from Table IV, page 10 */
  double f12     = 0.512;
  double f22     = -32.1;
  double f32     = -154;
  double f51     = -4.77;
  
  /* The following quantities were taken from the relation given in eq. (11), */
  /* page 7: fi3 = 64 - 64.*fi0 - 16.*fi1 - 4.*fi2; */
  double f13     = 64 - 16.*f11 - 4.*f12;
  double f23     = 64 - 16.*f21 - 4.*f22;
  double f33     = 64 - 16.*f31 - 4.*f32;
  double f53     = 64 - 64.*f50 - 16.*f51;
  
  /* this transformation is given in eq. (9), page (7) */
  double b1t     = b1*(f11*nu + f12*nu*nu + f13*nu*nu*nu);
  double b2t     = b2*(f21*nu + f22*nu*nu + f23*nu*nu*nu);
  double b3t     = b3*(f31*nu + f32*nu*nu + f33*nu*nu*nu);
  double b5t     = b5*(f50 + f51*nu + f53*nu*nu*nu);
  
  /* The functional form is taken from eq. (8), page 6. */
  double Lorb_eq_spin  = (0.00954*b3t*S*S*S + 0.0851*b2t*S*S - 0.194*b1t*S)/(1 - 0.579*b5t*S);
  
  /* These values are taken from Table IV, page 10: */
  double d10     = 0.322;
  double d11     = 9.33;
  double d20     = -0.0598;
  double d30     = 2.32;
  double d31     = -3.26;
  
  /* The functional form is taken from eq. (19a-c), page 10.*/
  double A1      = d10*xnu*nu*nu*(d11*nu+1);
  double A2      = d20*nu*nu*nu;
  double A3      = d30*xnu*nu*nu*nu*(d31*nu+1);
  
  /* The functional form is taken from eq. (15), page 9. */
  double Lorb_uneq_mass  = A1*Dchi + A2*Dchi*Dchi + A3*S*Dchi;
  
  return X1*X1*chi1+X2*X2*chi2 + Lorb_spin_zero + Lorb_eq_spin + Lorb_uneq_mass;
}

/** 
 *  Function: PrecessingRemnantSpin
 *  -------------------------------
 *  Final spin fit of precessing BH, eq. 20 of 
 *  https://arxiv.org/pdf/1611.00332.pdf
 *
 *   @param[in] dyn   : EOB dynamics
 *
 *   @return abh, spin of the remnant BH
 */
double PrecessingRemnantSpin(Dynamics *dyn)
{
  double SAmrg[3];
  double SBmrg[3];
  double Lhmrg[3];
  double Sperp[3];
  double SApar, SBpar, sqSperp;
  double Mbhf;

  double omgmrg = eob_mrg_momg(EOBPars->nu, EOBPars->X1, EOBPars->X2, EOBPars->chi1, EOBPars->chi2);
  
  /*find merger*/
  int imrg = find_point_bisection(omgmrg, dyn->spins->size, dyn->spins->data[EOB_EVOLVE_SPIN_Momg], 1);

  SAmrg[0] = dyn->spins->data[EOB_EVOLVE_SPIN_SxA][imrg];
  SAmrg[1] = dyn->spins->data[EOB_EVOLVE_SPIN_SyA][imrg];
  SAmrg[2] = dyn->spins->data[EOB_EVOLVE_SPIN_SzA][imrg];

  SBmrg[0] = dyn->spins->data[EOB_EVOLVE_SPIN_SxB][imrg];
  SBmrg[1] = dyn->spins->data[EOB_EVOLVE_SPIN_SyB][imrg];
  SBmrg[2] = dyn->spins->data[EOB_EVOLVE_SPIN_SzB][imrg];

  Lhmrg[0] = dyn->spins->data[EOB_EVOLVE_SPIN_Lx][imrg];
  Lhmrg[1] = dyn->spins->data[EOB_EVOLVE_SPIN_Ly][imrg];
  Lhmrg[2] = dyn->spins->data[EOB_EVOLVE_SPIN_Lz][imrg];

  vect_dot3(SAmrg, Lhmrg, &SApar);  
  vect_dot3(SBmrg, Lhmrg, &SBpar); 
  EOBPars->abhf = JimenezFortezaRemnantSpin(EOBPars->nu, EOBPars->X1, EOBPars->X2, SApar/SQ(EOBPars->X1), SBpar/SQ(EOBPars->X2));
  Mbhf          = JimenezFortezaRemnantMass(EOBPars->nu, EOBPars->X1, EOBPars->X2, SApar/SQ(EOBPars->X1), SBpar/SQ(EOBPars->X2));
  for(int i=0; i < IN3; i++)
    Sperp[i] = (SAmrg[i] - SApar*Lhmrg[i]) + (SBmrg[i] - SBpar*Lhmrg[i]);

  vect_dot3(Sperp, Sperp, &sqSperp);
  double sqMbhf = SQ(Mbhf);
  
  return sqrt(SQ(EOBPars->abhf) + sqSperp/SQ(sqMbhf));

}

/** 
 *  Function: QNMHybridFitCab
 *  -------------------------
 *  Compute QNM fits for the 22 mode for spinning systems
 *  and the postpeak template parameters
 *  From TODO: add reference
 *
 *   
 *   @param[in] nu                              : symmetric mass ratio
 *   @param[in] X1                              : mass fraction of body 1
 *   @param[in] X2                              : mass fraction of body 2
 *   @param[in] chi1                            : dimensionless spin of body 1
 *   @param[in] chi2                            : dimensionless spin of body 2
 *   @param[in] aK                              : a1 + a2
 *   @param[in] Mbh                             : final mass of the remnant BH
 *   @param[in] abh                             : final spin of the remnant BH
 *   @param[out] ca1,ca2,ca3,ca4,cb1,cb2,cb3,cb4 : phenomenological parameters of the postpeak template
 *   @param[out] sigmar,sigmai                   : real and imaginary QNM frequency
 *
 */
void QNMHybridFitCab(double nu, double X1, double X2, double chi1, double chi2, double aK, 
		     double Mbh, double abh,  
		     double *ca1, double *ca2, double *ca3, double *ca4, double *cb1, double *cb2, double *cb3, double *cb4, 
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
      
  const int usespins = EOBPars->use_spins;  

  int modeon[KMAX];
  const int k21 = 0;
  const int k22 = 1;
  const int k33 = 4;
  const int k44 = 8;
 
  for (int k=0; k<KMAX; k++) {
    modeon[k] = 0; /* off */
    sigmar[k] = sigmai[k] = 0.;
    ca1[k] = ca2[k] = ca3[k] = ca4[k] = 0.;
    cb1[k] = cb2[k] = cb3[k] = cb4[k] = 0.;
  }

  if (!(usespins)) {

    modeon[k21]=1;
    modeon[k22]=1;
    modeon[k33]=1;
    
    /* Last updates: 05/09/2017 from CoM extrapolated SXS data */
    
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

    /* Apply user-input deviations from merger amplitude */
    if (EOBPars->delta_Alm_mrg_size > 0)
      apply_mode_deviations(Amrg, EOBPars->delta_Alm_mrg, EOBPars->delta_Alm_mrg_k, EOBPars->delta_Alm_mrg_size);
    /* Throw error if deviations from merger or QNM frequencies are in input */
    if (EOBPars->delta_Omglm_mrg_size + EOBPars->delta_alphalm0_size + EOBPars->delta_omglm0_size > 0)
      errorexit("Merger frequency and QNM frequency, damping time deviations incompatible with QNMHybridFitCab if usespins == 0.\n");
    
  } else {

    modeon[k22]=1;
    
    /* Setting up coefficients from the phenomenological description of the ringdown.
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

    /* Apply QNM deviations */
    if (EOBPars->delta_alphalm0_size + EOBPars->delta_omglm0_size > 0)
      QNM_deviations(alpha1, omega1, alpha21);
    
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

    /* Throw error if merger frequency deviation given in input */
    if (EOBPars->delta_Omglm_mrg_size > 0)
      errorexit("Merger frequency deviation incompatible with QNMHybridFitCab.\n");
    
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

    /* Apply merger amplitude deviation */
    if (EOBPars->delta_Alm_mrg_size > 0)
      apply_mode_deviations(Amrg, EOBPars->delta_Alm_mrg, EOBPars->delta_Alm_mrg_k, EOBPars->delta_Alm_mrg_size);
    
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
      ca1[k] = Amrg[k] * alpha1[k] * cosh_c3A * cosh_c3A / c2A[k];
      ca2[k] = c2A[k];
      ca3[k] = c3A[k];
      ca4[k] = Amrg[k] - ca1[k] * tanh(c3A[k]);
      cb2[k] = alpha21[k];
      cb3[k] = c3phi[k];
      cb4[k] = c4phi[k];
      cb1[k] = Domg[k] * (1+c3phi[k]+c4phi[k]) / (cb2[k]*(c3phi[k] + 2.*c4phi[k]));
    }
  }
  
}

/** 
 *  Function: QNMHybridFitCab_HM
 *  ----------------------------
 *  Compute QNM fits for HMs for spinning systems
 *  and the postpeak template parameters
 *  From arXiv:2001.09082 
 *
 *   @param[in] nu                               : symmetric mass ratio
 *   @param[in] X1                               : mass fraction of body 1
 *   @param[in] X2                               : mass fraction of body 2
 *   @param[in] chi1                             : dimensionless spin of body 1
 *   @param[in] chi2                             : dimensionless spin of body 2
 *   @param[in] aK                               : a1 + a2
 *   @param[in] Mbh                              : final mass of the remnant BH
 *   @param[in] abh                              : final spin of the remnant BH
 *   @param[out] ca1,ca2,ca3,ca4,cb1,cb2,cb3,cb4 : phenomenological parameters of the postpeak template
 *   @param[out] sigmar,sigmai                   : real and imaginary part of QNM frequency
 *
 */
void QNMHybridFitCab_HM(double nu, double X1, double X2, double chi1, double chi2, double aK, 
			double Mbh, double abh,  
			double *ca1, double *ca2, double *ca3, double *ca4, double *cb1, double *cb2, double *cb3, double *cb4, 
			double *sigmar, double *sigmai)
{
  const double a12   = X1*chi1 - X2*chi2;
  const double X12   = X1 - X2;
  const double Shat  = 0.5*(aK + a12*X12);
  const double Sbar  = 0.5*(a12 + aK*X12);
  const double af    = abh; 
  const double nu2   = SQ(nu);
  const double nu3   = nu2*nu;
  const double nu4   = SQ(nu2);
  const double af2   = SQ(af);
  const double af3   = af2*af;
  const double Shat2 = SQ(Shat); 
  const double Shat3 = Shat2*Shat;
  const double Shat4 = SQ(Shat2);
  const double Sbar2 = SQ(Sbar);
  const double X12_2 = SQ(X12);
  const double aK2   = SQ(aK);
  const double aeff  = aK + 1./3.*a12*X12;
  const double aeff2 = SQ(aeff);
  const double aeff3 = aeff2*aeff;
  
  double c3A[KMAX], c3phi[KMAX], c4phi[KMAX], Domg[KMAX], Amrg[KMAX], c2A[KMAX], omgmrg[KMAX];
      
  const int usespins = EOBPars->use_spins;  

  int modeon[KMAX];
  const int k21 = 0;
  const int k22 = 1;
  const int k31 = 2;
  const int k32 = 3;
  const int k33 = 4;
  const int k41 = 5;
  const int k42 = 6;
  const int k43 = 7;
  const int k44 = 8;
  const int k55 = 13;
  
  for (int k=0; k<KMAX; k++) {
    modeon[k] = 0; /* off */
    sigmar[k] = sigmai[k] = 0.;
    ca1[k] = ca2[k] = ca3[k] = ca4[k] = 0.;
    cb1[k] = cb2[k] = cb3[k] = cb4[k] = 0.;
  }
  
  /* Defining the test-particle data used in the fits of the peak. */
  double ATP[KMAX], omgTP[KMAX];
  
  ATP[0]    = 0.5238781992;
  ATP[1]    = 1.44959;
  ATP[2]    = 0.0623783;
  ATP[3]    = 0.1990192432;
  ATP[4]    = 0.5660165890;
  ATP[5]    = 0.00925061;
  ATP[6]    = 0.0314363901;
  ATP[7]    = 0.0941569508;
  ATP[8]    = 0.2766182761;
  ATP[13]   = 0.151492;
  omgTP[0]  = 0.2906425497;
  omgTP[1]  = 0.273356;
  omgTP[2]  = 0.411755;
  omgTP[3]  = 0.4516072248;
  omgTP[4]  = 0.4541278937;
  omgTP[5]  = 0.552201;
  omgTP[6]  = 0.6175331548;
  omgTP[7]  = 0.6361300619;
  omgTP[8]  = 0.6356586393;
  omgTP[13] = 0.818117;
  
  /* If RWZ-normalization is needed. */
  /*
  for (int k=0; k<KMAX; k++) {
    double l = LINDEX[k];
    ATP[k] = ATP[k]/sqrt((l+2)*(l+1)*l*(l-1));
  }
  */
      
  if (!(usespins)) {

    modeon[0]  = 1;
    modeon[1]  = 1;
    modeon[2]  = 1;
    modeon[3]  = 1;
    modeon[4]  = 1;
    modeon[5]  = 1;
    modeon[6]  = 1;
    modeon[7]  = 1;
    modeon[8]  = 1;
    modeon[13] = 1;
    
    // Non spinning peak fits
    /* (l=2, m=1)*/
    Amrg[0]    = ATP[0]*X12*(1+9.0912*nu+3.9331*nu2)/(1+11.108*nu);
    omgmrg[0]  = omgTP[0]*(1-0.060432*nu+1.9995*nu2)/(1+0.23248*nu);
    
    /* (l=2, m=2)*/
    Amrg[1]    = ATP[1]*(1 - 0.041285*nu + 1.5971*nu2);
    omgmrg[1]  = omgTP[1]*(1 + 0.84074*nu + 1.6976*nu2);

    /* (l=3, m=1)*/
    Amrg[2]    = ATP[2]*X12*(1 - 5.49*nu + 10.915*nu2);
    omgmrg[2]  = omgTP[2]*(1 + 7.5362*nu2)/(1 - 2.7555*nu + 38.572*nu2);

    /* (l=3, m=2)*/
    Amrg[3]    = ATP[3]*(1-3*nu)*(1 - 6.142*nu + 11.372*nu2)/(1 - 3.6448*nu);
    omgmrg[3]  = omgTP[3]*(1 - 9.0214*nu + 21.078*nu2)/(1 - 8.6636*nu + 19.493*nu2);

    /* (l=3, m=3)*/
    Amrg[4]    = ATP[4]*X12*(1 + 0.098379*nu + 3.8179*nu2);
    omgmrg[4]  = omgTP[4]*(1 + 1.1054*nu + 2.2957*nu2);

    /* (l=4, m=1)*/
    Amrg[5]    = ATP[5]*X12*(1-2*nu)*(1 - 8.4449*nu + 26.825*nu2)/(1 - 1.2565*nu);
    omgmrg[5]  = omgTP[5]*(1 - 10.876*nu + 37.904*nu2)/(1 - 11.194*nu + 42.77*nu2);

    /* (l=4, m=2)*/
    Amrg[6]    = ATP[6]*(1-3*nu)*(1 - 4.7096*nu + 7.3253*nu2);
    omgmrg[6]  = omgTP[6]*(1 - 7.0558*nu + 12.738*nu2)/(1 - 6.0595*nu + 9.3389*nu2);

    /* (l=4, m=3)*/
    Amrg[7]    = ATP[7]*X12*(1-2*nu)*(1 - 5.7791*nu + 12.589*nu2)/(1 - 3.3039*nu);
    omgmrg[7]  = omgTP[7]*(1 - 9.0124*nu + 22.011*nu2)/(1 - 8.732*nu + 20.518*nu2);

    /* (l=4, m=4)*/
    Amrg[8]    = ATP[8]*(1-3*nu)*(1 - 3.6757*nu + 0.32156*nu2)/(1 - 3.6784*nu);
    omgmrg[8]  = omgTP[8]*(1 + 3.2876*nu - 29.122*nu2)/(1 + 1.696*nu - 22.761*nu2);

    /* (l=5, m=5)*/
    Amrg[13]   = ATP[13]*X12*(1 - 2*nu)*(1 - 0.29628*nu + 6.4207*nu2);
    omgmrg[13] = omgTP[13]*(1 - 2.8918*nu - 3.2012*nu2)/(1 - 3.773*nu);
    
    
    // Non spinning postmerger parameter fits
    /* (l=2, m=1)*/
    c3A[0] 	= (0.23882 - 2.2982*nu + 5.7022*nu2)/(1 - 7.7463*nu + 27.266*nu2);
    c3phi[0]	= (2.6269 - 37.677*nu + 181.61*nu2)/(1 - 16.082*nu + 89.836*nu2);
    c4phi[0]	= (4.355 - 53.763*nu + 188.06*nu2)/(1 - 18.427*nu + 147.16*nu2);
    
    /* (l=2, m=2)*/
    c3A[1]    = -0.56187 + 0.75497*nu;
    c3phi[1]  = (+4.4414 - 63.107*nu + 296.64*nu2)/(1 - 13.299*nu + 69.129*nu2);
    c4phi[1]  = (7.1508 - 109.47*nu)/(1 + 556.34*nu + 287.42*nu2);

    /* (l=3, m=1)*/
    c3A[2]    = (3.5042 - 55.171*nu + 217*nu2)/(1 - 15.749*nu + 605.17*nu3);
    c4phi[2]  = 3.6485 + 5.4536*nu;

    if (nu>0.08271){
      c3phi[2] = (-6.1719 + 29.617*nu + 254.24*nu2)/(1 - 1.5435*nu);
    } else {
      double nu_d    = 0.08271;
      c3phi[2] = (-6.1719 + 29.617*nu_d + 254.24*nu_d*nu_d)/(1 - 1.5435*nu_d);
    }

    /* (l=3, m=2)*/
    c3A[3]    = (0.1877 - 3.0017*nu + 19.501*nu2)/(1 - 1.8199*nu) - exp(-703.67*(nu - 2./9.)*(nu - 2./9.));
    c3phi[3]  = (0.90944 - 1.8924*nu + 3.6848*nu2)/(1 - 8.9739*nu + 21.024*nu2);
    c4phi[3]  = (2.3038 - 50.79*nu + 334.41*nu2)/(1 - 18.326*nu + 99.54*nu2);

    /* (l=3, m=3)*/
    c3A[4]    = -0.39337 + 0.93118*nu;
    c3phi[4]  = 3.1017 - 6.5849*nu;
    c4phi[4]  = (3.4521 - 24.153*nu + 53.029*nu2)/(1 + 3.1413*nu);
    
    /* (l=4, m=1)*/
    c3A[5]    = 11.47 + 10.936*nu;
    c4phi[5]  = 1.6629 + 11.497*nu;

    if (nu >= 10./121.) {
      c3phi[5] = -6.0286+46.632*nu;
    } else {
      c3phi[5] = -2.1747;
    }

    /* (l=4, m=2)*/
    c3A[6]    = (0.27143 - 2.2629*nu + 4.6249*nu2)/(1 - 7.6762*nu + 15.117*nu2);
    c3phi[6]  = (2.2065 - 17.629*nu + 65.372*nu2)/(1 - 4.7744*nu + 3.1876*nu2);

    if (nu >= 2.5/12.25) {
      c4phi[6] = (132.56 - 1155.5*nu + 2516.8*nu2)/(1 - 3.8231*nu);
    } else {
      c4phi[6] = (-0.58736 + 16.401*nu)/(1 - 4.5202*nu);
    }

    /* (l=4, m=3)*/
    c3A[7]    = (-0.02833 + 2.8738*nu - 31.503*nu2 + 93.513*nu3)/(1 - 10.051*nu + 156.14*nu3);
    c3phi[7]  = (2.284 - 23.817*nu + 70.952*nu2)/(1 - 10.909*nu + 30.723*nu2);
    c4phi[7]  = (2.4966 - 6.2043*nu)/(1 - 252.47*nu4);
    
    /* (l=4, m=4)*/
    c3A[8]    = (-0.25808 + 0.84605*nu) + 1.2376*exp(-6054.7*(nu - 2.5/12.25)*(nu - 2.5/12.25));
    c3phi[8]  = (2.3328 - 9.4841*nu + 19.719*nu2)/(1 - 2.904*nu);
    c4phi[8]  = 0.94564 + 3.2761*nu;

    /* (l=5, m=5)*/
    c3A[13]   = (-0.19751 + 3.607*nu - 14.898*nu2)/(1 - 20.046*nu + 108.42*nu2);
    c3phi[13] = 0.83326 + 10.945*nu;
    c4phi[13] = (0.45082 - 9.5961*nu + 52.88*nu2)/(1 - 19.808*nu + 99.078*nu2);
    
  } else {
    
    modeon[0]  = 1;
    modeon[1]  = 1;
    modeon[2]  = 0;
    modeon[3]  = 1;
    modeon[4]  = 1;
    modeon[5]  = 0;
    modeon[6]  = 1;
    modeon[7]  = 1;
    modeon[8]  = 1;
    modeon[13] = 1;
    
    /* a_vars */
    double a1, a2, a3;
    double a1Omg, a2Omg, a3Omg, a4Omg;
    double a1Amp, a2Amp, a3Amp, a1AmpS, a2AmpS, a3AmpS;
    
    /* b_vars */
    double b1, b2, b3, b4, b5, b6, b7, b8, b9;
    double b1Omg, b2Omg, b3Omg, b4Omg;
    double b11Omg, b12Omg, b13Omg, b14Omg, b21Omg, b22Omg, b23Omg;
    double b24Omg, b31Omg, b32Omg, b33Omg, b34Omg, b41Omg, b42Omg; 
    double b43Omg, b44Omg, b3Amp, b4Amp;
    double b1Amax, b2Amax, b3Amax, b1Amp, b2Amp, b11Amp, b12Amp, b13Amp, b14Amp, b21Amp, b22Amp, b23Amp, b24Amp, b31Amp, b32Amp, b33Amp, b34Amp;
    
    /* c_vars */
    double c1, c2, c3, c4;
    double c1Omg, c2Omg, c5Omg, c6Omg;
    double c11Omg, c12Omg, c13Omg, c21Omg, c22Omg, c23Omg, c31Omg, c32Omg, c33Omg, c41Omg, c42Omg, c43Omg;
    double c1Amax, c2Amax, c3Amax, c4Amax, c1Amp, c2Amp, c3Amp, c4Amp, c5Amp, c6Amp;
    double c11Amp, c12Amp, c13Amp, c14Amp, c21Amp, c22Amp, c23Amp, c24Amp;
    double c31Amp, c32Amp, c41Amp, c42Amp;
    
    /* others vars */
    double n1Omg, n2Omg, d1Omg, d2Omg;
    double omg1, omg2, orb, num, denom, OmgOrb, Omgspin, omgOrb, omgS, omgorb, omgspin;
    double Aorb, Aspin, scale, Amax1, Amax2, num_A, denom_A;
    
    // Peak Frequency
    /* (l=2, m=1)*/
    a1Omg     = -0.563075;
    a2Omg     = 3.28677;
    b1Omg     = 0.179639;
    b2Omg     = -0.302122;
    c11Omg    = -1.20684;
    c21Omg    = 0.425645;
	
    omgOrb    = omgTP[0]*(1 + a1Omg*nu + a2Omg*nu2);
    b1        = b1Omg + c11Omg*nu;
    b2        = b2Omg + c21Omg*nu;
    omgS      = 1 + b1*Shat + b2*Shat2;
    omgmrg[0] = omgOrb*omgS;

    /* (l=2, m=2)*/
    /*
      omg1      = 0.84074;
      omg2      = 1.6976;
      orb       = omgTP[1]*(1+omg1*nu + omg2*nu2);
		
      // Equal Mass fit calibrated to the q=1 SXS data
      b1 	      = -0.42311;
      b2        = -0.066699;
      b3        = -0.83053;
		
      // Unequal Mass corrections to the q=1 fit based on SXS, BAM and TP data
      c1        = 0.15873;
      c2	      = -0.43361;
      c3        = 0.60589;
      c4        = -0.71383;
      num       = 1 + (b1+c1*X12)/(1+c2*X12)*Shat + b2*Shat2;
      denom     = 1 + (b3+c3*X12)/(1+c4*X12)*Shat;
      omgmrg[1] = (orb*num/denom);
    */
    /*
    c2 = -0.122735;
    c1 = 0.0857478;
    c4 = -0.0760023;
    c3 = 0.0826514;
    omg1 = (c2*SQ(X12) +c1*X12 -0.1416002395)*2.*Shat + 1;
    omg2 = (c4*SQ(X12) + c3*X12 -0.3484804901)*2.*Shat + 1;
    omgmrg[1] = (+0.481958619443355*nu2 +0.223976694441952*nu +0.273813064427363)*omg1/omg2;
    */
    omg1      = 0.84074;
    omg2      = 1.6976;
    orb       = omgTP[1]*(1+omg1*nu + omg2*nu2);
		
    b1 = 0.066045;
    b2 = -0.23876;
    b3 = 0.76819;
    b4 = -0.9201;

    num   = 1. + (-0.42311 + b1*X12)/(1. + b2*X12)*Shat -0.066699*Shat2;
    denom = 1. + (-0.83053 + b3*X12)/(1. + b4*X12)*Shat;
    omgmrg[1] = (orb*num/denom);
    
    /* (l=3, m=2)*/
    a1Omg     = -9.13525;
    a2Omg     = 21.488;
    a3Omg     = -8.81384;
    a4Omg     = 20.0595; 
    b1Omg     = -0.458126;
    b2Omg     = 0.0474616;
    b3Omg     = -0.486049;
    c11Omg    = 3.25319; c12Omg = 0.535555; c13Omg = -8.07905; 
    c21Omg    = 1.00066; c22Omg = -1.1333; c23Omg = 0.601572;
    
    b1        = (b1Omg + c11Omg*X12 + c12Omg*X12_2)/(1+c13Omg*X12);
    b2        = (b2Omg + c21Omg*X12 + c22Omg*X12_2)/(1+c23Omg*X12);
    b3        = b3Omg;
    omgOrb    = omgTP[3]*(1 + a1Omg*nu + a2Omg*nu2)/(1 + a3Omg*nu + a4Omg*nu2);
    omgS      = (1 + b1*aK + b2*aK2)/(1 + b3*aK);
    omgmrg[3] = omgOrb*omgS;
	
    /* (l=3, m=3)*/
    a1Omg     = 1.08224;
    a2Omg     = 2.59333;
    b1Omg     = -0.406161;
    b2Omg     = -0.0647944;
    b3Omg     = -0.748126;
    c1Omg     = 0.85777;
    c2Omg     = -0.70066;
    c5Omg     = 2.97025;
    c6Omg     = -3.96242;
    b1Omg     = (b1Omg + c1Omg*nu)/(1 + c2Omg*nu);
    b3Omg     = (b3Omg + c5Omg*nu)/(1 + c6Omg*nu);

    omgorb    = 1 + a1Omg*nu + a2Omg*nu2;
    omgspin   = (1 + b1Omg*Shat + b2Omg*Shat2)/(1 + b3Omg*Shat);
    omgmrg[4] = omgTP[4]*omgorb*omgspin;
    
    /* (l=4, m=2)*/
    a1Omg     = -7.44121;
    a2Omg     = 14.233;
    a3Omg     = -6.61754;
    a4Omg     = 11.4329;
    b1Omg     = -2.37589;
    b2Omg     = 1.97249;
    b3Omg     = -2.36107;
    b4Omg     = 2.16383;
    c11Omg    = 10.1045; 
    c12Omg    = -6.94127; 
    c13Omg    = 12.1857; 
    c21Omg    = -1.62866; 
    c22Omg    = -2.6756; 
    c23Omg    = -4.7536; 
    c31Omg    = 10.071; 
    c32Omg    = -6.7299; 
    c33Omg    = 12.0377;
    c41Omg    = -8.56139; 
    c42Omg    = -5.27136; 
    c43Omg    = 5.10653;
    
    b1        = (b1Omg + c11Omg*nu)/(1 + c12Omg*nu + c13Omg*nu2);
    b2        = (b2Omg + c21Omg*nu)/(1 + c22Omg*nu + c23Omg*nu2);
    b3        = (b3Omg + c31Omg*nu)/(1 + c32Omg*nu + c33Omg*nu2);
    b4        = (b4Omg + c41Omg*nu)/(1 + c42Omg*nu + c43Omg*nu2);
    omgOrb    = omgTP[6]*(1 + a1Omg*nu + a2Omg*nu2)/(1 + a3Omg*nu + a4Omg*nu2);
    omgS      = (1 + b1*Shat + b2*Shat2)/(1 + b3*Shat + b4*Shat2);
    omgmrg[6] = omgOrb*omgS;
    
    /* (l=4, m=3)*/
    a1Omg     = -9.02463;
    a2Omg     = 21.9802;
    a3Omg     = -8.75892;
    a4Omg     = 20.5624;
    b1Omg     = -0.973324;
    b2Omg     = -0.109921;
    b3Omg     = -1.08036;
    c11Omg    = 11.5224; 
    c12Omg    = -26.8421; 
    c13Omg    = -2.84285; 
    c21Omg    = 3.51943; 
    c22Omg    = -12.1688; 
    c23Omg    = -3.96385; 
    c31Omg    = 5.53433; 
    c32Omg    = 3.73988; 
    c33Omg    = 4.219;
    
    b1        = (b1Omg + c11Omg*nu + c12Omg*nu2)/(1 + c13Omg*nu);
    b2        = (b2Omg + c21Omg*nu + c22Omg*nu2)/(1 + c23Omg*nu);
    b3        = (b3Omg + c31Omg*nu + c32Omg*nu2)/(1 + c33Omg*nu);
    omgOrb    = omgTP[7]*(1 + a1Omg*nu + a2Omg*nu2)/(1 + a3Omg*nu + a4Omg*nu2);
    omgS      = (1 + b1*Shat + b2*Shat2)/(1+b3*Shat);
    omgmrg[7] = omgOrb*omgS;
    
    /* (l=4, m=4)*/ 
    n1Omg     = -0.964614;
    n2Omg     = -11.1828;
    d1Omg     = -2.08471;
    d2Omg     = -6.89287;
    b1Omg     = -0.445192;
    b2Omg     = -0.0985658;
    b3Omg     = -0.0307812;
    b4Omg     = -0.801552;
    b11Omg    = -0.92902;
    b12Omg    = 10.86310;
    b13Omg    = -4.44930;
    b14Omg    = 3.01808;
    b21Omg    = 0;
    b22Omg    = 1.62523;
    b23Omg    = -7.70486;
    b24Omg    = 15.06517;
    b31Omg    = 0;
    b32Omg    = 0;
    b33Omg    = 0;
    b34Omg    = 0;
    b41Omg    = 0.93790;
    b42Omg    = 8.36038;
    b43Omg    = -4.85774;
    b44Omg    = 4.80446;
    a1Omg     = (b1Omg + b11Omg*nu + b12Omg*nu2)/(1 + b13Omg*nu + b14Omg*nu2);
    a2Omg     = (b2Omg + b21Omg*nu + b22Omg*nu2)/(1 + b23Omg*nu + b24Omg*nu2);
    a3Omg     = (b3Omg + b31Omg*nu + b32Omg*nu2)/(1 + b33Omg*nu + b34Omg*nu2);
    a4Omg     = (b4Omg + b41Omg*nu + b42Omg*nu2)/(1 + b43Omg*nu + b44Omg*nu2);

    OmgOrb    = (1 + n1Omg*nu + n2Omg*nu2)/(1 + d1Omg*nu + d2Omg*nu2);
    Omgspin   = (1 + a1Omg*Shat + a2Omg*Shat2 + a3Omg*Shat3)/(1 + a4Omg*Shat);
    omgmrg[8] = omgTP[8]*OmgOrb*Omgspin;
    
    /* (l=5, m=5)*/
    b1 	       = +1.487294;
    b2 	       = -2.058537;
    b3         = +1.454248;
    b4         = -1.301284;
    OmgOrb     = omgTP[13]*(1 - 2.8918*nu - 3.2012*nu2)/(1 - 3.773*nu);
    Omgspin    = (1 + (-0.332703 + b1*nu)/(1 + b2*nu)*Shat)/(1 + (-0.675738 + b3*nu)/(1 + b4*nu)*Shat);
    omgmrg[13] = OmgOrb*Omgspin;
    
    
    // Peak Amplitude
    /* (l=2, m=1)*/
    b1    = +0.891139;
    b2 	  = -5.191702;
    b3 	  = +3.480139;
    b4 	  = +10.237782;
    b5 	  = -13.867475; 
    b6 	  = +10.525510;

    if (DEQUAL(nu,0.25,1e-9)) {
      double S_bar21 = - fabs(Sbar);
      Aspin = ((-0.4281863 + b1*nu + b2*nu2)*S_bar21 + (-0.335659 + b3*nu + b4*nu2)*S_bar21*S_bar21)/(1 + (+0.828923 + b5*nu + b6*nu2)*S_bar21);
    } else {
      Aspin = ((-0.4281863 + b1*nu + b2*nu2)*Sbar + (-0.335659 + b3*nu + b4*nu2)*Sbar2)/(1 + (+0.828923 + b5*nu + b6*nu2)*Sbar);
    }
    
    Aorb    = ATP[0]*X12*(1 + 3.3362232268*nu + 3.4708521429*nu2)/(1 + 4.7623643259*nu);
    Amrg[0] = Aorb + Aspin;
    
    /* (l=2, m=2)*/    
    /*
      scale    = 1 - Shat*omgmrg[1];
      Amax1    = -0.041285;
      Amax2    = 1.5971;
      b1Amax   = -0.74124;
      b2Amax   = -0.088705;
      b3Amax   = -1.0939;
      c1Amax   = 0.15372;
      c2Amax   = 0.32082;
      c3Amax   = -0.072578;
      c4Amax   = 0.44091;
      num_A    = 1 + (b1Amax + c1Amax*X12)/(1 + c2Amax*X12)*Shat + b2Amax*Shat2;
      denom_A  = 1 + (b3Amax + c3Amax*X12)/(1 + c4Amax*X12)*Shat;
	
      Aorb     = ATP[1]*(1 + Amax1*nu + Amax2*nu2);
      Amrg[1]  = nu*Aorb*scale*(num_A/denom_A);
    */
    /*
    scale = 1 - 0.5*aeff*omgmrg[1];
    c4Amax = -0.0820894;
    c3Amax = 0.176126;
    c2Amax = -0.150239;
    c1Amax = 0.20491;    
    num_A    = (c4Amax*SQ(X12) + c3Amax*X12 - 0.2935238329)*aeff + 1.;
    denom_A  = 1. + (c2Amax*SQ(X12) + c1Amax*X12 -0.4728707630)*aeff;    
    Aorb = +1.826573640739664*nu2 + 0.100709438291872*nu + 1.438424467327531;
    Amrg[1] = Aorb*scale*(num_A/denom_A);
    */
    scale = 1. - omgmrg[1]*Shat;

    b1 = 0.4446696;
    b2 = -0.3254310;
    b3 = 0.4582812;
    b4 = -0.2124477;

    num_A   = 1. + (-0.741 + b1*X12)/(1. + b2*X12)*Shat - 0.0887*Shat2;
    denom_A = 1. + (-1.094 + b3*X12)/(1. + b4*X12)*Shat;
      
    Aorb     = ATP[1]*(1. - 0.041285*nu + 1.5971*nu2);
    Amrg[1]  = Aorb*scale*(num_A/denom_A);
    
    /* (l=3, m=2)*/
    a1      = -6.06831;
    a2      = 10.7505;    
    a3      = -3.68883; 
    b1      = -0.258378;
    b2      = 0.679163;
    c11Amp  = 4.36263; 
    c12Amp  = -12.5897; 
    c13Amp  = -7.73233; 
    c14Amp  = 16.2082; 
    c21Amp  = 3.04724; 
    c22Amp  = 46.5711;
    c23Amp  = 2.10475; 
    c24Amp  = 56.9136;    
    b1      = (b1 + c11Amp*nu + c12Amp*nu2)/(1 + c13Amp*nu + c14Amp*nu2);
    b2      = (b2 + c21Amp*nu + c22Amp*nu2)/(1 + c23Amp*nu + c24Amp*nu2);
 
    scale   = 1 + aK*pow(fabs(omgmrg[3])/2.,1./3.);
    // FIXME: Different from Matlab. Complex for negative frequencies
    Aorb    = ATP[3]*(1-3*nu)*(1 + a1*nu + a2*nu2)/(1 + a3*nu);
    Aspin   = (1 + b1*aK)/(1 + b2*aK);
    Amrg[3] = Aorb*scale*Aspin; 

    /* (l=3, m=3)*/
    a1Amp   = -0.22523;
    a2Amp   = 3.0569;
    a3Amp   = -0.396851; 
    b1Amp   = 0.100069;
    b2Amp   = -0.455859;
    c1Amp   = -0.401156;
    c2Amp   = -0.141551;
    c3Amp   = -15.4949; 
    c4Amp   = 1.84962;
    c5Amp   = -2.03512;
    c6Amp   = -4.92334;
    b1Amp   = (b1Amp + c1Amp*nu)/(1 + c2Amp*nu + c3Amp*nu2);
    b2Amp   = (b2Amp + c4Amp*nu)/(1 + c5Amp*nu + c6Amp*nu2);

    if (DEQUAL(nu,0.25,1e-9)) {
      double a12_33 = fabs(a12);
      Aspin   = (b1Amp*a12_33)/(1 + b2Amp*a12_33);
    } else {
      Aspin   = (b1Amp*a12)/(1 + b2Amp*a12);
    }
    
    Aorb    = (1 + a1Amp*nu + a2Amp*nu2)/(1 + a3Amp*nu);
    Amrg[4] = ATP[4]*X12*Aorb + Aspin;
    
    /* (l=4, m=2)*/
    a1Amp   = -4.56243;
    a2Amp   = 6.4522;
    
    b1Amp = -1.63682;
    b2Amp = 0.854459;
    b3Amp = 0.120537;
    b4Amp = -0.399718;
    
    c11Amp  = 6.53943; 
    c12Amp  = -4.00073; 
    c21Amp  = -0.638688; 
    c22Amp  = -3.94066; 
    c31Amp  = -0.482148; 
    c32Amp  = -3.9999999923319502; 
    c41Amp  = 1.25617; 
    c42Amp  = -4.04848;
    b1      = (b1Amp + c11Amp*nu)/(1 + c12Amp*nu);
    b2      = (b2Amp + c21Amp*nu)/(1 + c22Amp*nu);
    b3      = (b3Amp + c31Amp*nu)/(1 + c32Amp*nu);
    b4      = (b4Amp + c41Amp*nu)/(1 + c42Amp*nu);

    Aorb    = ATP[6]*(1-3*nu)*(1 + a1Amp*nu + a2Amp*nu2);
    Aspin   = (1 + b1*Shat + b2*Shat2)/(1 + b3*Shat + b4*Shat2);
    Amrg[6] = Aorb*Aspin*(1 + aK*pow(omgmrg[k42]/2.,1./3.));

    /* (l=4, m=3)*/
    if (DEQUAL(nu,0.25,1e-9)) {
      b1 = 0.00452129; 
      b2 = -0.00471163; 
      b3 = 0.0291409; 
      b4 = -0.351031;
      
      Amrg[7] = (b1 + b2*a12 + b3*a12*a12)/(1 + b4*a12);
    } else {
      a1Amp = -5.74386;
      a2Amp = 12.6016;
      a3Amp = -3.27435;
      b1    = +0.249099;
      b2    = -7.345984;
      b3    = +108.923746;
      b4    = -0.104206;
      b5    = +7.073534;
      b6    = -44.374738;
      b7    = +3.545134;
      b8    = +1.341375;
      b9    = -19.552083;

      Aorb    = ATP[7]*X12*(1-2*nu)*(1 + a1Amp*nu + a2Amp*nu2)/(1 + a3Amp*nu);
      Aspin   = (+((-0.02132252+b1*nu)/(1+b2*nu+b3*nu2))*aK+((+0.02592749+b4*nu)/(1+b5*nu+b6*nu2))*aK2)/(1+((-0.826977+b7*nu)/(1+b8*nu+b9*nu2))*aK);
      Amrg[7] = Aorb + Aspin;
    }

    /* (l=4, m=4)*/
    a1Amp  = -3.7082;
    a2Amp  = 0.280906;
    a3Amp  = -3.71276; 
    b1Amp  = -0.316647;
    b2Amp  = -0.062423;
    b3Amp  = -0.852876;
    b11Amp = 1.2436;
    b12Amp = -1.60555;
    b13Amp = -4.05685; 
    b14Amp = 1.59143;
    b21Amp = 0.837418;
    b22Amp = -2.93528; 
    b23Amp = -11.5591;
    b24Amp = 34.1863;
    b31Amp = 0.950035; 
    b32Amp = 7.95168;
    b33Amp = -1.26899;
    b34Amp = -9.72147;

    a1AmpS = (b1Amp + b11Amp*nu + b12Amp*nu2)/(1 + b13Amp*nu + b14Amp*nu2);
    a2AmpS = (b2Amp + b21Amp*nu + b22Amp*nu2)/(1 + b23Amp*nu + b24Amp*nu2);
    a3AmpS = (b3Amp + b31Amp*nu + b32Amp*nu2)/(1 + b33Amp*nu + b34Amp*nu2);

    scale   = 1 - 0.5*Shat*omgmrg[8];
    Aorb    = (1-3*nu)*(1 + a1Amp*nu + a2Amp*nu2)/(1 + a3Amp*nu);
    Aspin   = (1 + a1AmpS*Shat + a2AmpS*Shat2)/(1 + a3AmpS*Shat);
    Amrg[8] = ATP[8]*Aorb*Aspin*scale;

    /* (l=5, m=5)*/
    b1 = +5.720690;
    b2 = +44.868515;
    b3 = +12.777090;
    b4 = -42.548247;

    //Aorb     = ATP[13]*sqrt(1-4*nu)*(1-2*nu)*(1 - 0.29628*nu + 6.4207*nu2);
    //Aorb     = ATP[13]*sqrt(1-4*nu)*(1-2*nu)*(11.2008774621215608 *nu2 -0.9750925916632546 *nu + 1);
    Aorb     = ATP[13]*sqrt(1-4*nu)*(1-2*nu)*(9.9813227734275785*nu2 -0.6458814376485329*nu + 0.9823996027545479);
    Aspin    = (0.04360530/(1 + b1*nu + b2*nu2)*a12)/(1 - 0.5769451/(1+b3*nu+b4*nu2)*a12);
    Amrg[13] = Aorb + Aspin;
      
    // c3A
    /* (l=2, m=2)*/
    b1 = 0.1659421;
    b2 = -0.2560047;
    b3 = -0.9418946;
    c3A[k22] = -0.5585 + 0.81196*nu + (-0.398576+b1*X12)*Shat + (0.099805+b2*X12)*Shat2 + (0.72125+b3*X12)*Shat3;
    /*
      b1 =  0.0169543;
      b2 = -0.0799343;
      b3 = -0.115928;
      double c3A_nu =  0.8298678603 * nu - 0.5615838975;
      double c3A_eq =  (b3*X12 + 0.0907476903)*aeff3 + (b2*X12 + 0.0227344099)*aeff2 + (b1*X12 - 0.1994944332)*aeff;
      c3A[k22]            =  c3A_nu + c3A_eq;
    */
      
    /* (l=3, m=3)*/
    b1 = -0.3502608;
    b2 = 1.587606;
    b3 = -1.555325;
    c3A[k33] = -0.41455 + 1.3225*nu + (b1 + b2*X12 + b3*X12_2)*a12;

    /* (l=4, m=4)*/
    b1 = -9.614738;
    b2 = 122.461125;
    c3A[k44] = -0.41591 + 3.2099*nu + b1*nu*Shat + b2*nu*Shat2;

    /* (l=5, m=5)*/
    b1 = -7.063079;
    b2 = +65.464944;
    b3 = -2.055335;
    b4 = -0.585373;
    b5 = -12.631409;
    b6 = +19.271346;
    //c3A[k55] = b1 + b2*nu + (b3 + b4*X12)*a12 + (b5 + b6*X12)*SQ(a12);
    //c3A[k55] = 9.1187519178640084 *nu +  -0.5970347579708830 + (b3 + b4*X12)*a12 + (b5 + b6*X12)*SQ(a12);
    c3A[k55] = 9.1187519178640084*nu  - 0.5970347579708830 + (b3 + b4*X12)*a12 + (b5 + b6*X12)*SQ(a12);
    
    // c3phi
    /* (l=2, m=2)*/
    b1 = -1.323643;
    b2 = -3.555007;
    b3 = 7.011267;
    b4 = 32.737824;
    c3phi[k22] = 3.8436 + 0.71565*nu + (5.12794 + b1*X12)*Shat + (9.9136 + b2*X12)*Shat2 + (-4.1075 + b3*X12)*Shat3 +(-31.5562 + b4*X12)*Shat4;
    /*
      b1 = -0.462321;
      b2 = -0.904512;
      b3 =  0.437747;
      b4 =  1.8275;
      double c3phi_nu     =  0.4558467286*nu + 3.8883812141;
      double c3phi_equal  =  (b4*X12 - 2.0575868122)*16*Shat4 +(b3*X12 - 0.5051534498)*8*Shat3 +(b2*X12 + 2.5742292762)*4*Shat2 +(b1*X12 + 2.5599640181)*2*Shat;
      c3phi[k22]          = c3phi_nu + c3phi_equal;
    */
      
    /* (l=3, m=3)*/
    b1 = -0.634377;
    b2 = 5.983525; 
    b3 = -5.881900;
    c3phi[k33] = 3.0611 - 6.1597*nu;// + (b1 + b2*X12 + b3*X12_2)*Shat;

    /* (l=4, m=4)*/
    b1 = 7.911653;
    b2 = 21.181688;
    c3phi[k44] = (3.6662 - 30.072*nu +76.371*nu2)/(1 - 3.5522*nu);// + (-4.9184 + b1*X12)*Shat + (-15.6772 + b2*X12)*Shat2;

    /* (l=5, m=5)*/
    b1 = -1.510167;
    b2 = +30.569461;
    b3 = -2.687133;
    b4 = +4.873750;
    b5 = -14.629684;
    b6 = +19.696954;
    //c3phi[k55] = b1 + b2*nu + (b3 + b4*X12)*Shat + (b5 + b6*X12)*Shat2;
    //c3phi[k55] = 373.312597*nu2 -59.69284 *nu+ 4.226238 + (b3 + b4*X12)*Shat + (b5 + b6*X12)*Shat2;
    c3phi[k55] = 373.3125969012880319*nu2 - 59.6928356312470854*nu + 4.2262379605181648 + (b3 + b4*X12)*Shat + (b5 + b6*X12)*Shat2;
    
    // c4phi
    /* (l=2, m=2)*/
    b1 = 0.779683;
    b2 = -0.069638;
    c4phi[k22] = 1.4736 + 2.2337*nu + (8.26539 + b1*X12)*Shat + (14.2053 + b2*X12)*Shat2;
    /*
      b1 = -0.449976;
      b2 = -0.980913;
      double c4phi_nu     =  2.0822327682 * nu + 1.4996868401;
      double c4phi_equal  =  (b2*X12 + 3.5695199109)*4*Shat2 + (b1*X12 + 4.1312404030)*2*Shat;
      c4phi[k22]          =  c4phi_nu + c4phi_equal;
    */
      
    /* (l=3, m=3)*/
    b1 = -3.877528;
    b2 = 12.043300; 
    b3 = -6.524665;
    c4phi[k33] = 1.789 - 5.6684*nu;// + (b1 + b2*X12 + b3*X12_2)*Shat;

    /* (l=4, m=4)*/
    b1 = 11.746452;
    b2 = 34.922883;
    c4phi[k44] = 0.21595 + 23.216*nu;// + (-3.4207 + b1*X12)*Shat + (-15.5383 + b2*X12)*Shat2;

    /* (l=5, m=5)*/
    b1 = -1.383721;
    b2 = +56.871881;
    b3 = +7.198729;
    b4 = -3.870998;
    b5 = -25.992190;
    b6 = +36.882645;
    //c4phi[k55] = b1 + b2*nu + (b3 + b4*X12)*Shat + (b5 + b6*X12)*Shat2;
    //c4phi[k55] = 14.9111373110275380 *nu  +  1.3639723340485870+ (b3 + b4*X12)*Shat + (b5 + b6*X12)*Shat2;
    c4phi[k55] = 14.9111373110275380*nu + 1.3639723340485870 + (b3 + b4*X12)*Shat + (b5 + b6*X12)*Shat2;
    
    // These fits are taken from the testparticle limit and approximate the spinning case
    /* (l=3, m=1)*/
    if (DEQUAL(nu,0.25,1e-9)){
      Amrg[2]   = 0.;
      omgmrg[2] = 0.;
    }
    else{
      Amrg[2]   = ATP[2]*X12*(1 - 5.49*nu + 10.915*nu2);
      omgmrg[2] = omgTP[2]*(1 + 7.5362*nu2)/(1 - 2.7555*nu + 38.572*nu2);
    }
	    
    /* (l=4, m=1)*/
    Amrg[5]   = 0;//ATP[5]*X12*(1-2*nu)*(1 - 8.4449*nu + 26.825*nu2)/(1 - 1.2565*nu);
    omgmrg[5] = 0;//omgTP[5]*(1 - 10.876*nu + 37.904*nu2)/(1 - 11.194*nu + 42.77*nu2);

    // Postpeak parameter fits
    /* (l=2, m=1)*/
    c3A[0]   = (0.23882 - 2.2982*nu + 5.7022*nu2)/(1 - 7.7463*nu + 27.266*nu2);
    c3phi[0] = (2.6269 - 37.677*nu + 181.61*nu2)/(1 - 16.082*nu + 89.836*nu2);
    c4phi[0] = (4.355 - 53.763*nu + 188.06*nu2)/(1 - 18.427*nu + 147.16*nu2);

    /* (l=3, m=1)*/
    c3A[2]    = (3.5042 - 55.171*nu + 217*nu2)/(1 - 15.749*nu + 605.17*nu3);
    c4phi[2]  = 3.6485 + 5.4536*nu;

    if (nu>0.08271){
      c3phi[2]  = (-6.1719 + 29.617*nu + 254.24*nu2)/(1 - 1.5435*nu);
    } else {
      double nu_d    = 0.08271;
      c3phi[2]  = (-6.1719 + 29.617*nu_d + 254.24*nu_d*nu_d)/(1 - 1.5435*nu_d);
    }
    
    /* (l=3, m=2)*/
    c3A[3]   = (0.1877 - 3.0017*nu + 19.501*nu2)/(1 - 1.8199*nu) - exp(-703.67*(nu - 2./9.)*(nu - 2./9.));
    c3phi[3] = (0.90944 - 1.8924*nu + 3.6848*nu2)/(1 - 8.9739*nu + 21.024*nu2);
    c4phi[3] = (2.3038 - 50.79*nu + 334.41*nu2)/(1 - 18.326*nu + 99.54*nu2);

    /* (l=4, m=1)*/
    c3A[5]     = 11.47 + 10.936*nu;
    c4phi[k41] = 1.6629 + 11.497*nu;

    if (nu >= 10./121.) {
      c3phi[5] = -6.0286 + 46.632*nu;
    } else {
      c3phi[5] = -2.1747;
    }

    /* (l=4, m=2)*/
    c3A[6]   = (0.27143 - 2.2629*nu + 4.6249*nu2)/(1 - 7.6762*nu + 15.117*nu2);
    c3phi[6] = (2.2065 - 17.629*nu + 65.372*nu2)/(1 - 4.7744*nu + 3.1876*nu2);

    if (nu >= 2.5/12.25) {
      c4phi[6] = (132.56 - 1155.5*nu + 2516.8*nu2)/(1 - 3.8231*nu);
    } else {
      c4phi[6] = (-0.58736 + 16.401*nu)/(1 - 4.5202*nu);
    }

    /* (l=4, m=3)*/
    c3A[7]   = (-0.02833 + 2.8738*nu - 31.503*nu2 + 93.513*nu3)/(1 - 10.051*nu + 156.14*nu3);
    c3phi[7] = (2.284 - 23.817*nu + 70.952*nu2)/(1 - 10.909*nu + 30.723*nu2);
    c4phi[7] = (2.4966 - 6.2043*nu)/(1 - 252.47*nu4);

    /* (l=5, m=5)*/
    /*
    c3A[13]   = (-0.19751 + 3.607*nu - 14.898*nu2)/(1 - 20.046*nu + 108.42*nu2);
    c3phi[13] = 0.83326 + 10.945*nu;
    c4phi[13] = (0.45082 - 9.5961*nu + 52.88*nu2)/(1 - 19.808*nu + 99.078*nu2);
    */
  }

  /* Apply deviations from merger amplitude, frequency */
  if (VERBOSE & (EOBPars->delta_Alm_mrg_size + EOBPars->delta_Omglm_mrg_size > 0))
    printf("Applying deviations from NR amplitude, frequency at merger.\n");
  if (EOBPars->delta_Alm_mrg_size > 0)
    apply_mode_deviations(Amrg, EOBPars->delta_Alm_mrg, EOBPars->delta_Alm_mrg_k, EOBPars->delta_Alm_mrg_size);
  if (EOBPars->delta_Omglm_mrg_size > 0)
    apply_mode_deviations(omgmrg, EOBPars->delta_Omglm_mrg, EOBPars->delta_Omglm_mrg_k, EOBPars->delta_Omglm_mrg_size);
  
  if (DEQUAL(nu,0.25,1e-9) && DEQUAL(chi1,chi2,1e-9)){
    modeon[0] = modeon[2] = modeon[4] = modeon[5] = modeon[7] = modeon[13] = 0;
  }

  double alpha21[KMAX], alpha1[KMAX], omega1[KMAX];
  QNM_coefs(af, alpha21, alpha1, omega1);
  
  for (int k=0; k<KMAX; k++) {
    if (modeon[k]) {
      sigmar[k] = alpha1[k];
      sigmai[k] = omega1[k];
      Domg[k] 	= omega1[k] - Mbh*omgmrg[k];
    }
  }
  
  for (int k=0; k<KMAX; k++) {
    if (modeon[k]) {
      c2A[k] = 0.5*alpha21[k];
      double cosh_c3A= cosh(c3A[k]);  
      ca1[k] = Amrg[k]*alpha1[k]*cosh_c3A*cosh_c3A/c2A[k];
      ca2[k] = c2A[k];
      ca3[k] = c3A[k];
      ca4[k] = Amrg[k] - ca1[k]*tanh(c3A[k]);
      cb2[k] = alpha21[k];
      cb3[k] = c3phi[k];
      cb4[k] = c4phi[k];
      cb1[k] = Domg[k]*(1+c3phi[k]+c4phi[k])/(cb2[k]*(c3phi[k] + 2*c4phi[k]));
    }
  }
  
}

/**
 *  Function: apply_mode_deviations
 *  ---------------------------
 *  Apply mode-by-mode fractional deviations to merger or NQC point quantities
 *
 *   @param[in] x_mrg                           : Bare array
 *   @param[in] delta_x                         : Fractional deviation array
 *   @param[in] klm                             : List of modes to deform
 *   @param[in] Nk                              : Number of modes to deform
 *
 */
void apply_mode_deviations(double *x_mrg, double *delta_x, int *klm, int Nk)
{
  int km;

  for (int k=0; k<Nk; k++) {
      km = klm[k];
      x_mrg[km]   = x_mrg[km]*(1. + delta_x[km]);
  }
}

/** 
 *  Function: QNMHybridFitCab_HM_Pompili23
 *  ----------------------------
 *  Global fits of the ringdown coefficients from the 
 *  Pompili et. al. paper, arXiv:2303.18039.
 *  Note: only available for (2,1), (3,3), (4,4)
 *
 *   @param[in] nu                              : symmetric mass ratio
 *   @param[in] X1                              : mass fraction of body 1
 *   @param[in] X2                              : mass fraction of body 2
 *   @param[in] chi1                            : dimensionless spin of body 1
 *   @param[in] chi2                            : dimensionless spin of body 2
 *   @param[in] Mbh                             : final mass of the remnant BH
 *   @param[out] ca1,ca2,ca3,ca4,cb1,cb2,cb3,cb4 : phenomenological parameters of the postpeak template
 *   @param[in] sigmar,sigmai                   : real and imaginary QNM frequency
 *
 */
void QNMHybridFitCab_HM_Pompili23(double nu, double X1, double X2, double chi1, double chi2, double aK, double Mbh, double abh, 
			double *ca1, double *ca2, double *ca3, double *ca4, double *cb1, double *cb2, double *cb3, double *cb4,
      double *sigmar, double *sigmai)
{

  int knqcpeak22_size = EOBPars->knqcpeak22_size;  
  int *knqcpeak22     = EOBPars->knqcpeak22;

  const double delta = sqrt(1-4.*nu);
  const double nu2   = SQ(nu);
  const double nu3   = nu2*nu;
  const double nu4   = nu3*nu;
  const double chi_S = (chi1+chi2)/2.;
  const double chi_A = (chi1-chi2)/2.;

  const double chi       = chi_S + chi_A*delta/(1.-2.*nu);
  const double chisq     = chi*chi;
  const double chi3      = chisq*chi;
  const double chi4      = chi3*chi;
  const double chi_21A   = delta*chi_S/(1-1.3*nu) + chi_A;
  const double chi_21A2  = SQ(chi_21A);
  const double chi_21A3  = chi_21A2*chi_21A;
  const double chi_21D   = delta*chi_S/(1-2.*nu)  + chi_A;

  double A[KMAX], dA[KMAX], omg[KMAX], domg[KMAX];
  double c1f, c2f, d1f, d2f, c1c, c2c, d1c;
  
  /* Compute A, dA, omg, domg */
  if (EOBPars->binary==BINARY_BBH){eob_nqc_point_HM_peak22(NULL,A, dA, omg, domg);
  }else{eob_nqc_point_BHNS_HM(NULL,A, dA, omg, domg);}

  /* Apply deviations from fitted merger quantities here by hand */
  if (EOBPars->delta_Alm_mrg_size > 0)
    apply_mode_deviations(A, EOBPars->delta_Alm_mrg, knqcpeak22, knqcpeak22_size);
  if (EOBPars->delta_Omglm_mrg_size > 0)
    apply_mode_deviations(omg, EOBPars->delta_Omglm_mrg, knqcpeak22, knqcpeak22_size);

  for(int j=0; j<knqcpeak22_size; j++){
    int k = knqcpeak22[j];

    /* The QNM frequencies need rescaling */
    sigmar[k] /= Mbh;
    sigmai[k] /= Mbh;

    /* Fits */
    switch(k)
    {
      case(-1):
        if (VERBOSE) printf("QNMHybridFitCab_HM_Pompili23: no modes selected in knqcpeak22\n");
        break;
      case(0):
        c1f = 0.173462*chisq*nu - 0.028873*chisq + 0.197467*chi*nu2 - 0.026139*chi - 2.934735*nu3 + 1.009106*nu2 - 0.112721*nu + 0.099889;
        c2f = 0.183489*chi3 + 0.10573*chisq - 20.792825*chi*nu2 + 6.867746*chi*nu - 0.484948*chi - 54.917585*nu3 + 16.466312*nu2 + 0.426316*nu - 0.92208;
        d1f = 0.018467*chi4 + 0.398621*chi3*nu - 0.050499*chi3 - 0.877201*chisq*nu2 + 0.414553*chisq*nu - 0.068277*chisq - 10.648526*chi*nu3 + 4.104918*chi*nu2 - 0.723576*chi*nu + 0.039227*chi + 42.715534*nu4 - 18.280603*nu3 + 2.236592*nu2 - 0.048094*nu + 0.16335;
        d2f = exp(0.814085*chi3 - 1.197363*chisq*nu + 0.560622*chisq + 6.44667*chi*nu2 - 5.630563*chi*nu + 0.949586*chi + 91.269183*nu3 - 27.329751*nu2 + 1.101262*nu + 1.040761);
        break;
      case(4):
        c1f = -0.00956*chi3 + 0.029459*chisq*nu - 0.020264*chisq - 0.494524*chi*nu2 + 0.169463*chi*nu - 0.026285*chi - 5.847417*nu3 + 1.957462*nu2 - 0.171682*nu + 0.093539;
        c2f = -0.057346*chi3 + 0.237107*chisq*nu - 0.094285*chisq - 4.250609*chi*nu2 + 1.763105*chi*nu - 0.315826*chi + 14.801916*nu3 - 7.060581*nu2 + 1.158627*nu - 0.646888;
        d1f = -0.016524*chi3 + 0.221466*chisq*nu - 0.066323*chisq + 0.678442*chi*nu2 - 0.261264*chi*nu + 0.006664*chi + 2.316434*nu3 - 2.192227*nu2 + 0.424582*nu + 0.161577;
        d2f = exp(0.275999*chi3 - 1.830695*chisq*nu + 0.512734*chisq + 29.072515*chi*nu2 - 10.581319*chi*nu + 1.310643*chi + 324.310223*nu3 - 124.681881*nu2 + 13.200426*nu + 0.410855);
        break;
      case(8):
        c1f = 4.519504*chi*nu2 - 1.489036*chi*nu + 0.068403*chi - 1656.065439*nu4 + 817.835726*nu3 - 127.055379*nu2 + 6.921968*nu + 0.009386;
        c2f = 0.964861*chi3*nu - 0.185226*chi3 - 12.647814*chisq*nu2 + 5.264969*chisq*nu - 0.539721*chisq - 254.719552*chi*nu3 + 105.698791*chi*nu2 - 12.107281*chi*nu + 0.2244*chi - 393.727702*nu4 + 145.32788*nu3 - 15.556222*nu2 + 1.592449*nu - 0.677664;
        d1f = -0.020644*chi3 + 0.494221*chisq*nu - 0.127074*chisq + 4.297985*chi*nu2 - 1.284386*chi*nu + 0.062684*chi - 44.280815*nu3 + 11.021482*nu2 - 0.162943*nu + 0.166018;
        d2f = exp(37.735116*chi*nu2 - 12.516669*chi*nu + 1.309868*chi - 528.368915*nu3 + 155.115196*nu2 - 6.612448*nu + 0.787726);
        break;
      default:
        errorexit("Ringdown fits only implemented for (2,1), (3,3) and (4,4) at the moment");
    }

    /* Correct normalization of A, dA */
    int l = LINDEX[k];
    A[k] *= sqrt((l+2)*(l+1)*l*(l-1));
    dA[k]*= sqrt((l+2)*(l+1)*l*(l-1));

    /* Constrained coefficients*/
    c1c = (dA[k]  + sigmar[k]*A[k])*cosh(c2f)*cosh(c2f)/c1f;
    c2c =  A[k]   - (dA[k] + sigmar[k]*A[k])*cosh(c2f)*sinh(c2f)/c1f;
    d1c = (omg[k] - sigmai[k])*(1+d2f)/(d1f*d2f);

    /* Finalize*/
    ca1[k] = c1c;
    ca2[k] = c1f;
    ca3[k] = c2f;
    ca4[k] = c2c;
    cb1[k] = -d1c;
    cb2[k] = d1f;
    cb3[k] = d2f;
    cb4[k] = 0.;
  }
}

/** 
 *  Function: QNM_coefs
 *  -------------------
 *  Fits of the ringdown QNM coefs, from
 *  TODO: add reference
 *
 *   @param[in]  af      : spin of the final BH
 *   @param[out] alpha21 : difference between the inverse damping times of first overtone and fundamental mode
 *   @param[out] alpha1  : real part of the fundamental QNM frequency (damping time)
 *   @param[out] omega1  : imaginary part of the fundamental QNM frequency
 *
 */void QNM_coefs(double af, double *alpha21, double *alpha1, double *omega1)
{

  const int binary = EOBPars->binary;

  for (int k=0; k<KMAX; k++) {
    alpha21[k] = alpha1[k] = omega1[k] = 0.;
  }

  double af2 = pow(af,2);
  double af3 = pow(af,3);

  /* l = 2, m = 1*/
  omega1[0]  = +0.373672*(1 - 0.79546*af - 0.1908*af2 + 0.11460*af3)/(1 - 0.96337*af - 0.1495*af2 + 0.19522*af3);
  alpha1[0]  = +0.0889623*(1 - 1.31253*af - 0.21033*af2 + 0.52502*af3)/(1 - 1.30041*af - 0.1566*af2 + 0.46204*af3);
  alpha21[0] = +0.184952*(1 - 1.1329*af - 0.3520*af2 + 0.4924*af3)/(1 - 1.10334*af - 0.3037*af2 + 0.4262*af3);
  
  /* l = 2, m = 2*/
  omega1[1]  = +0.373672*(1 - 1.5367*af + 0.5503*af2)/(1 - 1.8700*af + 0.9848*af2 - 0.10943*af3);
  alpha1[1]  = +0.08896*(1 - 1.90036*af + 0.86200*af2 + 0.0384893*af3)/(1 - 1.87933*af + 0.88062*af2);
  alpha21[1] = +0.184953*(1 - 1.89397*af + 0.88126*af2 + 0.0130256*af3)/(1 - 1.83901*af + 0.84162*af2);

  /* if(binary==BINARY_BBH){
   * omega1[1]  = +0.373672*(1 - 1.5367*af + 0.5503*af2)/(1 - 1.8700*af + 0.9848*af2 - 0.10943*af3);
   * alpha1[1]  = +0.08896*(1 - 1.90036*af + 0.86200*af2 + 0.0384893*af3)/(1 - 1.87933*af + 0.88062*af2);
   * alpha21[1] = +0.184953*(1 - 1.89397*af + 0.88126*af2 + 0.0130256*af3)/(1 - 1.83901*af + 0.84162*af2);
   * }else{ 
   * // for the BHNS model, we use the old fits
   * double omega1_c = -0.0598837831*af3 + 0.8082136788*af2 - 1.7408467418*af + 1;
   * double omega1_d = -0.2358960279*af3 + 1.3152369374*af2 - 2.0764065380*af + 1;
   * omega1[1] =  0.3736716844*(omega1_c/omega1_d);
   *
   * double alpha1_c = 0.1211263886*af3 + 0.7015835813*af2 - 1.8226060896*af + 1;
   * double alpha1_d = 0.0811633377*af3 + 0.7201166020*af2 - 1.8002031358*af + 1;
   * alpha1[1] = 0.0889623157 * (alpha1_c/alpha1_d);
   *	    
   * double alpha21_c = 0.4764196512*af3 - 0.0593165805*af2 - 1.4168096833*af + 1;
   * double alpha21_d = 0.4385578151*af3 - 0.0763529088*af2 - 1.3595491146*af + 1;
   * alpha21[1] = 0.1849525596*(alpha21_c/alpha21_d); 
   * } 
   */

  /* l = 3, m = 1*/
  omega1[2]  = +0.599443*(1 - 0.70941*af - 0.16975*af2 + 0.08559*af3)/(1 - 0.82174*af - 0.16792*af2 + 0.14524*af3);
  alpha1[2]  = +0.0927030*(1 - 1.2345*af - 0.30447*af2 + 0.5446*af3)/(1 - 1.2263*af - 0.24223*af2 + 0.47738*af3);
  alpha21[2] = +0.188595*(1 - 1.035*af - 0.3816*af2 + 0.4486*af3)/(1 - 1.023*af - 0.3170*af2 + 0.3898*af3);
  
  /* l = 3, m = 2*/
  omega1[3]  = +0.599443*(1 - 0.251*af - 0.891*af2 + 0.2706*af3)/(1 - 0.475*af - 0.911*af2 + 0.4609*af3);
  alpha1[3]  = +0.0927030*(1 - 1.58277*af + 0.2783*af2 + 0.30503*af3)/(1 - 1.56797*af + 0.3290*af2 + 0.24155*af3);
  alpha21[3] = +0.188595*(1 - 1.5212*af + 0.1563*af2 + 0.3652*af3)/(1 - 1.4968*af + 0.1968*af2 + 0.3021*af3);
  
  /* l = 3, m = 3*/
  omega1[4]  = +0.599443*(1 - 1.84922*af + 0.9294*af2 - 0.07613*af3)/(1 - 2.18719*af + 1.4903*af2 - 0.3014*af3);
  alpha1[4]  = +0.0927030*(1 - 1.8310*af + 0.7568*af2 + 0.0745*af3)/(1 - 1.8098*af + 0.7926*af2 + 0.0196*af3);
  alpha21[4] = +0.188595*(1 - 1.8011*af + 0.7046*af2 + 0.0968*af3)/(1 - 1.7653*af + 0.7176*af2 + 0.0504*af3);

  /* l = 4, m = 1*/
  omega1[5]  = +0.809178*(1 - 0.68647*af - 0.1852590*af2 + 0.0934997*af3)/(1 - 0.77272*af - 0.1986852*af2 + 0.1485093*af3);
  alpha1[5]  = +0.0941640*(1 + 1.1018882*af - 0.88643*af2 - 0.78266*af3)/(1 + 1.1065495*af - 0.80961*af2 - 0.68905*af3);
  alpha21[5] = +0.190170*(1 + 1.0590157*af - 0.8650630*af2 - 0.75222*af3)/(1 + 1.0654880*af - 0.7830051*af2 - 0.65814*af3);
  
  /* l = 4, m = 2*/
  omega1[6]  = +0.809178*(1 - 0.6644*af - 0.3357*af2 + 0.1425*af3)/(1 - 0.8366*af - 0.2921*af2 + 0.2254*af3);
  alpha1[6]  = +0.0941640*(1 - 1.44152*af + 0.0542*af2 + 0.39020*af3)/(1 - 1.43312*af + 0.1167*af2 + 0.32253*af3);
  alpha21[6] = +0.190170*(1 - 1.38840*af + 0.39333*af3)/(1 - 1.37584*af + 0.0600017*af2 + 0.32632*af3);
  
  /* l = 4, m = 3*/
  omega1[7]  = +0.809178*(1 - 1.8397*af + 0.9616*af2 - 0.11339*af3)/(1 - 2.0979*af + 1.3701*af2 - 0.2675*af3);
  alpha1[7]  = +0.0941640*(1 - 1.7177*af + 0.5320*af2 + 0.1860*af3)/(1 - 1.7065*af + 0.5876*af2 + 0.120939*af3);
  alpha21[7] = +0.190170*(1 - 1.6860*af + 0.4724*af2 + 0.2139*af3)/(1 - 1.6684*af + 0.5198*af2 + 0.1508*af3);
  
  /* l = 4, m = 4*/
  omega1[8]  = +0.809178*(1 - 1.83156*af + 0.9016*af2 - 0.06579*af3)/(1 - 2.17745*af + 1.4753*af2 - 0.2961*af3);
  alpha1[8]  = +0.0941640*(1 - 1.8662*af + 0.8248*af2 + 0.0417*af3)/(1 - 1.8514*af + 0.8736*af2 - 0.0198*af3);
  alpha21[8] = +0.190170*(1 - 1.8546*af + 0.8041*af2 + 0.0507*af3)/(1 - 1.8315*af + 0.8391*af2 - 0.0051*af3);

  /* l = 5, m = 5*/
  omega1[13]  = +1.012295*(1 - 1.5659*af + 0.5783*af2)/(1 - 1.9149*af + 1.0668*af2 - 0.14663*af3);
  alpha1[13]  = +0.0948705*(1 - 1.8845*af + 0.8585*af2 + 0.0263*af3)/(1 - 1.8740*af + 0.9147*af2 - 0.0384*af3);
  alpha21[13] = +0.190947*(1 - 1.8780*af + 0.8467*af2 + 0.0315*af3)/(1 - 1.8619*af + 0.8936*af2 - 0.0293*af3);

  /* Apply user-input deviations from fitted values of alpha, omega */
  if (EOBPars->delta_alphalm0_size > 0 || EOBPars->delta_omglm0_size > 0) {
    if (VERBOSE) printf("Applying user-input deviations from QNM frequencies.\n");
    QNM_deviations(alpha1, omega1, alpha21);
  }
  if (EOBPars->use_spins == MODE_SPINS_GENERIC && EOBPars->use_effective_QNMs == 1) {
    // For precessing spins, use effective QNM frequencies estimated by https://arxiv.org/pdf/2301.06558, Eq. (34)
    if (DEBUG) printf("Using effective QNM frequencies for precessing spins; cos(beta_final) = %.9f\n", EOBPars->cbeta_final);
    double w21_by_22 = omega1[0]/omega1[1];
    double cbeta     = fabs(EOBPars->cbeta_final);
    int m_of_k[KMAX] = {
      1, 2,
      1, 2, 3,
      1, 2, 3, 4,
      1, 2, 3, 4, 5,
      1, 2, 3, 4, 5, 6,
      1, 2, 3, 4, 5, 6, 7,
      1, 2, 3, 4, 5, 6, 7, 8
    };
    for (int k = 0; k < KMAX; k++) omega1[k] *= 1. - m_of_k[k]*(1. - w21_by_22)*(1. - cbeta);
  }
  
}

/**
 *  Function: QNM_deviations
 *  ------------------------
 *  Apply user-input deviations from fitted values of alpha, omega
 *
 *   @param[in] alpha1,omega1 : Real and imaginary part of QNM frequency
 *   @param[in] alpha21       : Difference between real parts (inverse damping times) of first overtone and fundamental QNM
 *
 */
void QNM_deviations(double *alpha1, double *omega1, double *alpha21)
{
  int km;

  if (EOBPars->delta_alphalm0_size > 0) {
    for (int k=0; k<EOBPars->delta_alphalm0_size; k++) {
      km = EOBPars->delta_alphalm0_k[k];
      alpha21[km] -= alpha1[km]*EOBPars->delta_alphalm0[km];
      alpha1[km]  += alpha1[km]*EOBPars->delta_alphalm0[km];
    }
  }
  if (EOBPars->delta_omglm0_size > 0) {
    for (int k=0; k<EOBPars->delta_omglm0_size; k++) {
      km = EOBPars->delta_omglm0_k[k];
      omega1[km]  += omega1[km]*EOBPars->delta_omglm0[km];
    }
  }
}

/** 
 *  Function: eob_approxLR
 *  ----------------------
 *  Light-ring approximation, vacuum nonspinning case
 *
 *   @param[in] nu      : symmetric mass ratio
 *
 *   @return approximate estimate for r_LR
 */
double eob_approxLR(const double nu) 
{
  const double x = 1.-4.*nu;
  const double r0 = 1.818461553848201e+00; // x=0 (q=1)
  const double r1 = 1.821043720041472e+00; // x=1 (3M)
  return r0*x + r1;
}

/** 
 *  Function: get_mrg_timestep
 *  --------------------------
 *  Compute optimized timestep after Omega_peak
 *
 *   @param[in] q      : mass ratio
 *   @param[in] chi1   : dimensionless spin of body 1
 *   @param[in] chi2   : dimensionless spin of body 2
 *
 *   @return dt to be used for merger interpolation
 */
double get_mrg_timestep(double q, double chi1, double chi2)
{
  double dt = 0.1;
  // ...
  return dt;
}

/** 
 *  Function: get_mrg_timestop
 *  --------------------------
 *  Compute optimized shift to stop integration after Omega_peak
 *
 *   @param[in] q      : mass ratio
 *   @param[in] chi1   : dimensionless spin of body 1
 *   @param[in] chi2   : dimensionless spin of body 2
 *
 *   @return time interval beyind tOmg_pk after which integration is stopped
 */
double get_mrg_timestop(double q, double chi1, double chi2)
{ 
  double tstop = 2.0;
  if (EOBPars->use_flm == USEFLM_HM) {
    if((EOBPars->binary==BINARY_BHNS_TD) || (EOBPars->binary==BINARY_BHNS) ){ // for BHNS
      tstop = 3.;
    }else{
      tstop = 10.;
    }
    
  }
  return tstop;
}

/** 
 *  Function: fmode_resonance_dress_Love
 *  ------------------------------------
 *   Effective Love number for the f-mode resonance model  
 *   Eq.(6.51) of https://arxiv.org/abs/1608.01907  
 *   Eq.(11)   of https://arxiv.org/abs/1702.02053  
 *   This routine computes the dressing factor for the star ell-th Love number
 *
 *   @param[in] nu        : mass ratio
 *   @param[in] r         : orbital separation
 *   @param[in] bomgf     : mA omega_f
 *   @param[in] ell       : multipolar index
 *   @param[out] dtides   : array of dressing factors
 *   @param[out] dtides_u : array of first derivative of dressing factors w.r.t. 1/r
 *   
 */
double fmode_resonance_dress_Love(double nu,
				  double r, // orbital radius
				  double bomgf, // = mA omega_f
				  int ell, 
				  double *dtides, double *dtides_u)
{
  if (ell<2 || ell>4) errorexit("f-mode model only for l=2,3,4");
  double a[5] = {0, 0, 1./4., 3./8., 29./64.}; // a_ell, ell = (0), (1), 2, 3, 4
  double b[5] = {0, 0, 3./4., 5./8., 35./64.}; // b_ell
  const double five_o_three = 1.66666666666666666; // = 5./3.;
  const double one_o_six = 0.166666666666666666; // = 1./6.;
  const int emm = ell;
  const double Omega = pow(r, -1.5);
  const double Omg1 = 0.375; // = 3./8.;
  double x = bomgf/(emm*Omega);
  double x_u = - 1.5*r*x;
  double x2 = SQ(x);
  double x5_3 = pow(x, five_o_three);
  const double eps = 256./5. * nu * pow(bomgf/emm,five_o_three);
  if (eps<=0.) errorexit("epsilon<=0"); 
  const double div_sqrt_eps = 1.0/sqrt(eps);
  /* Arguments of various functions */
  const double that = 1.6 * div_sqrt_eps * (1. - x5_3);
  const double that_u = - 2.6666666666666665 * div_sqrt_eps * x5_3/x * x_u;
  const double arg = Omg1 * SQ(that);// cos/sin arg
  const double arg_u = 0.75 * that * that_u;
  if (arg<=0.) errorexit("arg<=0"); 
  const double argfresnel = 0.61237243569579447*that;//sqrt(arg); 
  const double nrmfresnel = sqrt(Pi/(2.*Omg1));    
  /* Resonant term */
  const double x1 = x - 1.;
  double rterm = 0;
  double rterm_u = 0;
  if(fabs(x1)<1e-2) {
    /* Taylor expansion near x1 */
    double c0[5]  = {-0.0833333333333333,
		     -0.1157407407407407,
		     -0.006944444444444445,
		     0.01122256515775034,
		     -0.007120198902606311};
    double c1[5] = {-0.1157407407407407,
		    -0.01388888888888889,
		    0.03366769547325102,
		    -0.02848079561042524,
		    0.01919764200071127};
    const double x12 = x1*x1;
    const double x13 = x1*x12;
    const double x14 = x1*x13;
    rterm = c0[4]*x14 + c0[3]*x13 + c0[2]*x12 + c0[1]*x1 + c0[0];
    rterm_u = c1[4]*x14 + c1[3]*x13 + c1[2]*x12 + c1[1]*x1 + c1[0];
    rterm_u *= x_u;
  }
  else {
    /* Full formula */
    const double rt1 = 1./(x2 - 1.);
    const double rt2 = 1./(1. - x5_3);
    rterm   = x2*( rt1 + 5*one_o_six*rt2 );
    rterm_u = 2.*x*x_u*( -SQ(rt1) + 5*one_o_six*(1. - one_o_six*x5_3)*SQ(rt2) );
  }
  /* Fresnel term */
  const double FS = Fresnel_Sine_Integral(argfresnel);  
  const double FC = Fresnel_Cosine_Integral(argfresnel);
  const double FS1 = nrmfresnel*(0.5 + FS); 
  const double FC1 = nrmfresnel*(0.5 + FC);
  const double cosa = cos(arg);
  const double sina = sin(arg);
  const double fterm = div_sqrt_eps * x2 * (cosa*FS1 - sina*FC1);
  const double fterm_u = div_sqrt_eps * ( 2*x*x_u*(cosa*FS1 - sina*FC1) - x2*arg_u*(sina*FS1 + cosa*FC1) );
  /* Dressing factors */ 
  *dtides   = a[ell] + b[ell]*(rterm + fterm);
  *dtides_u = b[ell]*(rterm_u + fterm_u);

#if (0)//(DEBUG)
  printf("DEBUG(f-mode) Omega = %.6e x = %.6e x1 = %.6e\n",Omega,x,x1);
  printf("DEBUG(f-mode) alpha = %.6e alpha_u = %.6e\n",*dtides,*dtides_u);  
  printf("DEBUG(f-mode) that = %.6e\n",that);
  printf("DEBUG(f-mode) sqrt(eps) = %.6e\n",sqrt(eps));
  printf("DEBUG(f-mode) rt1 = %.6e rt2 = %.6e\n",rt1, rt2);
  printf("DEBUG(f-mode) r = %.6e r_u = %.6e\n",rterm, rterm_u);
  printf("DEBUG(f-mode) arg = %.6e\n",arg);
  printf("DEBUG(f-mode) cos = %.6e sin = %.6e\n",cosa,sina);  
  printf("DEBUG(f-mode) argfresnel = %.6e\n",argfresnel);  
  printf("DEBUG(f-mode) FS1 = %.6e FC1 = %.6e\n",FS1,FC1);
  printf("DEBUG(f-mode) f = %.6e f_u = %.6e\n",fterm, fterm_u);
#endif
  
}

/** 
 *  Function: fmode_resonance_dressing_factors
 *  ------------------------------------------
 *   Routine to update dressing factors for f-mode resonance model
 *   Only apply to gravitoelectric ell=2,3,4
 *   Note: if enforce_alpha_min is 1, the minimum value of the dressing
 *   factors is floored to 1.
 *  
 *   @param[in] r        : orbital separation
 *   @param[in] dyn      : EOB dynamics
 *   
 */
#define enforce_alpha_min (1) /* enforce min(alpha)=1 */
void fmode_resonance_dressing_factors(double r, Dynamics *dyn)
{
  if (!(EOBPars->use_tidal_fmode_model)) return;

  const double nu = EOBPars->nu;
  const int lmax = 4;
  
  for (int l=2; l<=lmax; l++) {
    if (EOBPars->bomgfA[l]<0. || EOBPars->bomgfB[l]<0.)
      errorexit("f-mode frequency cannot be negative");
    dyn->dress_tides_fmode_A[l] = dyn->dress_tides_fmode_B[l] = 1;
    dyn->dress_tides_fmode_A_u[l] = dyn->dress_tides_fmode_B_u[l] = 0;
  }
  
  /* Compute the dressing factors */
  if (EOBPars->LambdaAl2>0.) {
    for (int l=2; l<=lmax; l++) {
      fmode_resonance_dress_Love(nu, r, EOBPars->bomgfA[l], l,
				 &(dyn->dress_tides_fmode_A[l]),
				 &(dyn->dress_tides_fmode_A_u[l]));
      /* if (dyn->dress_tides_fmode_A[l]<=0.) errorexit("Dressing A factor<0"); */
#if(enforce_alpha_min)
      if (dyn->dress_tides_fmode_A[l]<1) {
	dyn->dress_tides_fmode_A[l] = 1;
      	dyn->dress_tides_fmode_A_u[l] = 0;
      }
#endif
    }
  }
  
  if (EOBPars->LambdaBl2>0.) {
    for (int l=2; l<=lmax; l++) {
      fmode_resonance_dress_Love(nu, r, EOBPars->bomgfB[l], l,
				 &(dyn->dress_tides_fmode_B[l]),
				 &(dyn->dress_tides_fmode_B_u[l]));
      /* if (dyn->dress_tides_fmode_B[l]<=0.) errorexit("Dressing B factor<0"); */
#if(enforce_alpha_min)
      if (dyn->dress_tides_fmode_B[l]<1) {
      	dyn->dress_tides_fmode_B[l] = 1;
	dyn->dress_tides_fmode_B_u[l] = 0;
      }
#endif
    }
  }
  
#if(DEBUG)
  if (EOBPars->output_dynamics) {
    const double Omega = pow(r, -1.5);
    FILE* fp;
    char fname[STRLEN];
    strcpy(fname,EOBPars->output_dir);
    strcat(fname,"/fmode_dressingfact.txt");
    if( access( fname, F_OK ) == 0 ) {
      if ((fp = fopen(fname, "a+")) == NULL)
	errorexits("error opening file",fname);
    } else {
      if ((fp = fopen(fname, "a+")) == NULL)
	errorexits("error opening file",fname);
      fprintf(fp, "# bomgfA = (%.6e %.6e %.6e) bomgfB = (%.6e %.6e %.6e)\n",
	      EOBPars->bomgfA[2],EOBPars->bomgfA[3],EOBPars->bomgfA[4],
	      EOBPars->bomgfB[2],EOBPars->bomgfB[3],EOBPars->bomgfB[4]);
      fprintf(fp, "# Omg alphaA2 alphaA3 alphaA4 alphaB2 alphaB3 alphaB4\n");
    }  
    fprintf(fp, "%.6e %.6e %.6e %.6e %.6e %.6e %.6e\n", Omega,
	    dyn->dress_tides_fmode_A[2],dyn->dress_tides_fmode_A[3],dyn->dress_tides_fmode_A[4],
	    dyn->dress_tides_fmode_B[2],dyn->dress_tides_fmode_B[3],dyn->dress_tides_fmode_B[4]);
    fclose(fp);
  }
#endif
  
}
 
/** 
 *  Function: fmode_resonance_dress_QOH
 *  ------------------------------------
 *   \details Routine to update the quadrupole, octupole, hexapole 
 *   using the dressed Lambdas
 *  
 *   @param[in] dyn      : EOB dynamics
 *   
 */
void fmode_resonance_dress_QOH(Dynamics *dyn)
{
  /* Dress the Lambda's */
  const double LamA2 = EOBPars->LambdaAl2 * dyn->dress_tides_fmode_A[2];
  const double LamB2 = EOBPars->LambdaBl2 * dyn->dress_tides_fmode_B[2];
  
  const double LamA2_u = EOBPars->LambdaAl2 * dyn->dress_tides_fmode_A_u[2];
  const double LamB2_u = EOBPars->LambdaBl2 * dyn->dress_tides_fmode_B_u[2];
  
  const double LamA2_uu = 0;
  const double LamB2_uu = 0;

  /* Compute the QOH */
  double C_QA = 1;
  double C_QA_u = 0;
  double C_QA_uu = 0;
  if (LamA2>0) 
    YagiYunes13_fit_logQ_drvts(LamA2,LamA2_u,LamA2_uu, &C_QA, &C_QA_u, &C_QA_uu);

  dyn->dressed_C_Q1 = C_QA;
  dyn->dressed_C_Q1_u = C_QA_u;
  dyn->dressed_C_Q1_uu = C_QA_uu;

  double C_QB = 1;
  double C_QB_u = 0;
  double C_QB_uu = 0;  
  if (LamB2>0)
    YagiYunes13_fit_logQ_drvts(LamB2,LamB2_u,LamB2_uu, &C_QB, &C_QB_u, &C_QB_uu);

  dyn->dressed_C_Q2 = C_QB;
  dyn->dressed_C_Q2_u = C_QB_u;
  dyn->dressed_C_Q2_uu = C_QB_uu;

  //TODO use other fits for OH ...
  dyn->dressed_C_Oct1 = EOBPars->C_Oct1;
  dyn->dressed_C_Oct2 = EOBPars->C_Oct2;
  dyn->dressed_C_Hex1 = EOBPars->C_Hex1; 
  dyn->dressed_C_Hex2 = EOBPars->C_Hex2;

  dyn->dressed_C_Oct1_u = dyn->dressed_C_Oct1_uu = 0;
  dyn->dressed_C_Oct2_u = dyn->dressed_C_Oct2_uu = 0;
  dyn->dressed_C_Hex1_u = dyn->dressed_C_Hex1_uu = 0;
  dyn->dressed_C_Hex2_u = dyn->dressed_C_Hex2_uu = 0;
  
}


/** 
 *  Function: eob_bhns_fit
 *  ----------------------
 *   \details Fits of BH remnant from BHNS CBC.
 *   From 
 *   https://arxiv.org/abs/1903.11622 
 *   https://arxiv.org/abs/2212.03909
 *   and updated for Gonzalez+ in prep.
 *  
 *   @param[in]  a      : spin of the BH
 *   @param[in]  nu     : symmetric mass ratio
 *   @param[out] mass   : mass of the remnant
 *   @param[out] spin   : spin of the remnant
 *   @param[in]  lambda : tidal parameter of the NS
 *   @param[in]  m_bh   : mass of the BBH remnant with same masses
 *   @param[in]  a_bh   : spin of the BBH remnant with same masses
 *    
 */
void eob_bhns_fit(double a, double nu, double *mass, double *spin, double lambda, double m_bh, double a_bh)
{
  /* m_bh & a_bh are the remnant mass & spin of Jimenez et al. fits (2017) */

     /* Mass Parameters */
  const double p102 = -0.0008738793880464456;
  const double p101 = 0.0005314645082267707; 
  const double p100 = 0.00048079143591861327;
  const double p112 = 0.0034086268737653345; 
  const double p111 = -0.0006054545388183284;
  const double p110 = -0.006365039897453894; 
  const double p122 = -0.006753889323869175; 
  const double p121 = -0.013585866660743297; 
  const double p120 = 0.02530645334339657;   
  const double p212 = -1.7262162415183805e-07;
  const double p211 = -5.40510452294782e-06;  
  const double p210 = 4.363002416591192e-05;  
  const double p222 = 7.737497619529003e-06;  
  const double p221 = 2.396600803932003e-05;  
  const double p220 = -2.1980019903488383e-05;
  const double p0   = 4.0788174743855836e-05;

  double p10 = p102*a*a + p101*a + p100;      
  double p11 = p112*a*a + p111*a + p110;      
  double p12 = p122*a*a + p121*a + p120;      
  double p21 = p212*a*a + p211*a + p210;      
  double p22 = p222*a*a + p221*a + p220;      

  double p1 = p10*nu + p11*nu*nu + p12*nu*nu*nu;
  double p2 = p21*nu*nu + p22*nu*nu*nu;         
                                                                                                              
  double m_fit = ( (1 + p1*lambda + p2*lambda*lambda) / ((1 + p0*nu*nu*lambda*lambda))  );                                                                                                                

     /* Spin Parameters */  
  const double  d112 = 0.0013574322384365632;
  const double  d111 = 0.004609266720050001;
  const double  d110 = 0.0009102105324348831;
  const double  d122 = -0.04975162622644385;
  const double  d121 = -0.02314136431572672;
  const double  d120 = -0.024890806466961492;
  const double  d132 = 0.17502248483357633;
  const double  d131 = 0.015427669282338003;
  const double  d130 = 0.12743167707858225;
  const double  d212 = 2.0730081552909773e-05;
  const double  d211 = -3.1406807993499036e-05;
  const double  d210 = 8.094258450191269e-06;
  const double  d222 = -0.0001652984994487448;
  const double  d221 = 0.00027540548306762236;
  const double  d220 = -4.4402028029482014e-05;
  const double  d232 = 0.0003856401528807317;
  const double  d231 = -0.0006293995220570712;
  const double  d230 = 2.1108847630144224e-05;
  const double  d312 =-1.561477963430114e-08;
  const double  d311 = 1.7300042605972545e-08;
  const double  d310 = -3.1476803536094817e-09;
  const double  d322 = 1.5278874946669358e-07;
  const double  d321 = -1.6399967561815403e-07;
  const double  d320 = 2.4337955107824662e-08;
  const double  d332 = -3.9321543091156496e-07;
  const double  d331 = 3.977949302502904e-07;
  const double  d330 = -2.4789449646084726e-08;
  const double  d41 = 7.867448219972702e-07;

  double  d11 = d112*a*a + d111*a + d110;
  double  d12 = d122*a*a + d121*a + d120;
  double  d13 = d132*a*a + d131*a + d130;
  double  d21 = d212*a*a + d211*a + d210;
  double  d22 = d222*a*a + d221*a + d220;
  double  d23 = d232*a*a + d231*a + d230;
  double  d31 = d312*a*a + d311*a + d310;
  double  d32 = d322*a*a + d321*a + d320;
  double  d33 = d332*a*a + d331*a + d330;

  double  d1 = d11*nu + d12*nu*nu + d13*nu*nu*nu;
  double  d2 = d21*nu + d22*nu*nu + d23*nu*nu*nu;
  double  d3 = d31*nu + d32*nu*nu + d33*nu*nu*nu;
  double  d4 = d41*nu;

  double a_fit = ( 1 + lambda*d1 + lambda*lambda*d2  + lambda*lambda*lambda*d3) / ( (1 + lambda*lambda*d4)*(1 + lambda*lambda*d4) );

  double mbh = m_bh * m_fit;
  double abh = a_bh * a_fit;

  *mass = mbh; // Mf/M                                                                                                                                                                                                                        
  *spin = abh;
}

/** 
 *  Function: bhns_cases
 *  ----------------------
 *   \details  Filter tidal disruption cases.
 *   From https://arxiv.org/abs/2212.03909 and updated for Gonzalez+ in prep.
 *  
 *   @param[in]  nu        : symmetric mass ratio
 *   @param[in]  af        : 
 *   @param[in]  chi1      : dimensionless spin of body 1
 *   @param[in]  lambda    : tidal parameter of the NS
 *   @param[out] binary    : type of system
 *   @param[out] usetidal  : tidal model
 *    
 */
void bhns_cases(double nu, double af, double chi1, double lambda, int *binary, int *usetidal)
{
  /* Amplitude peak 2024 */
  const double c102 = 0.025960639049437614; 
  const double c101 = 0.011583003041680442;
  const double c100 = -0.0016533947810772696;
  const double c112 = -0.31403656018939863;
  const double c111 = -0.12825460441414602;
  const double c110 = 0.13396228694902;
  const double c122 = 0.8917717620733512;
  const double c121 = 0.3150355510141175;
  const double c120 = -0.513304355534311;
  const double c202 = -1.984228201511826e-05;
  const double c201 = -2.804785772573806e-06;
  const double c200 = 8.018661470565014e-06;
  const double c212 = 0.0002223708961923682;
  const double c211 = 2.555424239372415e-05;
  const double c210 = -0.00010025951516653183;
  const double c222 = -0.0006022803171412016;
  const double c221 = -5.057203798367972e-05;
  const double c220 = 0.0002870181683963824;
  const double c31  = 0.009579950459173342;
  const double c32  = -0.01899088523426322;
  
  double p10 = c102*chi1*chi1 + c101*chi1 + c100;
  double p11 = c112*chi1*chi1 + c111*chi1 + c110;
  double p12 = c122*chi1*chi1 + c121*chi1 + c120;
  double p20 = c202*chi1*chi1 + c201*chi1 + c200;
  double p21 = c212*chi1*chi1 + c211*chi1 + c210;
  double p22 = c222*chi1*chi1 + c221*chi1 + c220;

  double p1 = p10*nu + p11*nu*nu + p12*nu*nu*nu;
  double p2 = p20*nu + p21*nu*nu + p22*nu*nu*nu;
  double p3 = c31*nu + c32*nu*nu;

  double Ap =  ( 1 + (p1*lambda + p2*lambda*lambda) ) / ( 1 + lambda*p3  );

  if(Ap<=0.85){
    *binary = BINARY_BHNS_TD; //Type I
    }else if(Ap<0.99){
      *binary = BINARY_BHNS; //Type III
    }else{
      *binary = BINARY_BBH; //Type II
      *usetidal = 0;
    }

    
}

/** 
 *  Function: QNMHybridFitCab_BHNS_HM
 *  ---------------------------------
 *  \details Compute QNM fits for HMs for BHNS systems
 *  and the postpeak template parameters
 *  From Gonzalez+2022 and updated for Gonzalez+ in prep.
 *
 *   @param[in] nu                               : symmetric mass ratio
 *   @param[in] X1                               : mass fraction of body 1
 *   @param[in] X2                               : mass fraction of body 2
 *   @param[in] chi1                             : dimensionless spin of body 1
 *   @param[in] chi2                             : dimensionless spin of body 2
 *   @param[in] aK                               : a1 + a2
 *   @param[in] Mbh                              : final mass of the remnant BH
 *   @param[in] abh                              : final spin of the remnant BH
 *   @param[out] ca1,ca2,ca3,ca4,cb1,cb2,cb3,cb4 : phenomenological parameters of the postpeak template
 *   @param[out] sigmar,sigmai                   : real and imaginary part of QNM frequency
 *
 */
void QNMHybridFitCab_BHNS_HM(double nu, double X1, double X2, double chi1, double chi2, double aK, 
			double Mbh, double abh,  
			double *ca1, double *ca2, double *ca3, double *ca4, double *cb1, double *cb2, double *cb3, double *cb4, 
			double *sigmar, double *sigmai)
{
  const double a12   = X1*chi1 - X2*chi2;
  const double X12   = X1 - X2;
  const double Shat  = 0.5*(aK + a12*X12);
  const double Sbar  = 0.5*(a12 + aK*X12);
  const double af    = abh; 
  const double nu2   = SQ(nu);
  const double nu3   = nu2*nu;
  const double nu4   = SQ(nu2);
  const double af2   = SQ(af);
  const double af3   = af2*af;
  const double Shat2 = SQ(Shat); 
  const double Shat3 = Shat2*Shat;
  const double Shat4 = SQ(Shat2);
  const double Sbar2 = SQ(Sbar);
  const double X12_2 = SQ(X12);
  const double aK2   = SQ(aK);
  const double aeff  = aK + 1./3.*a12*X12;
  const double aeff2 = SQ(aeff);
  const double aeff3 = aeff2*aeff;
  const double lambda = EOBPars->LambdaBl2;
  
  double c3A[KMAX], c3phi[KMAX], c4phi[KMAX], Domg[KMAX], Amrg[KMAX], c2A[KMAX], omgmrg[KMAX];
      
  const int usespins = EOBPars->use_spins;  

  int modeon[KMAX];
  const int k21 = 0;
  const int k22 = 1;
  const int k32 = 3;
  const int k33 = 4;
  const int k44 = 8;
  const int k55 = 13;
  
  for (int k=0; k<KMAX; k++) {
    modeon[k] = 0; /* off */
    sigmar[k] = sigmai[k] = 0.;
    ca1[k] = ca2[k] = ca3[k] = ca4[k] = 0.;
    cb1[k] = cb2[k] = cb3[k] = cb4[k] = 0.;
  }
    
    modeon[0]  = 1;
    modeon[1]  = 1;
    modeon[2]  = 0;
    modeon[3]  = 1;
    modeon[4]  = 1;
    modeon[5]  = 0;
    modeon[6]  = 0;
    modeon[7]  = 0;
    modeon[8]  = 1;
    modeon[13] = 1;
   
    double a1, a2, a3;
    double b1, b2, b3, b4, b5, b6, b7, b8, b9;
    double c1, c2, c3, c4;
    
  // Peak Frequency & Amplitude
    /* Modes: 21, 22, 32, 33, 44, 55 */ 
    peak_bhns(nu, lambda, chi1, X1, X2, abh, Amrg, omgmrg);

  // c3A
    /* l=2, m=1 */
    c3A[k21]   = (0.23882 - 2.2982*nu + 5.7022*nu2)/(1 - 7.7463*nu + 27.266*nu2);
    
    /* l=2, m=2 */
    b1 = 0.1659421;
    b2 = -0.2560047;
    b3 = -0.9418946;
    c3A[k22] = -0.5585 + 0.81196*nu + (-0.398576+b1*X12)*Shat + (0.099805+b2*X12)*Shat2 + (0.72125+b3*X12)*Shat3;
    /*
      b1 =  0.0169543;
      b2 = -0.0799343;
      b3 = -0.115928;
      double c3A_nu =  0.8298678603 * nu - 0.5615838975;
      double c3A_eq =  (b3*X12 + 0.0907476903)*aeff3 + (b2*X12 + 0.0227344099)*aeff2 + (b1*X12 - 0.1994944332)*aeff;
      c3A[k22]            =  c3A_nu + c3A_eq;
    */
      
    /* (l=3, m=3)*/
    b1 = -0.3502608;
    b2 = 1.587606;
    b3 = -1.555325;
    c3A[k33] = -0.41455 + 1.3225*nu + (b1 + b2*X12 + b3*X12_2)*a12;

    /* (l=4, m=4)*/
    b1 = -9.614738;
    b2 = 122.461125;
    c3A[k44] = -0.41591 + 3.2099*nu + b1*nu*Shat + b2*nu*Shat2;

    /* (l=5, m=5)*/
    b1 = -7.063079;
    b2 = +65.464944;
    b3 = -2.055335;
    b4 = -0.585373;
    b5 = -12.631409;
    b6 = +19.271346;
    //c3A[k55] = b1 + b2*nu + (b3 + b4*X12)*a12 + (b5 + b6*X12)*SQ(a12);
    //c3A[k55] = 9.1187519178640084 *nu +  -0.5970347579708830 + (b3 + b4*X12)*a12 + (b5 + b6*X12)*SQ(a12);
    c3A[k55] = 9.1187519178640084*nu  - 0.5970347579708830 + (b3 + b4*X12)*a12 + (b5 + b6*X12)*SQ(a12);
     
  // c3phi
    /* l=2, m=1 */
    c3phi[k21] = (2.6269 - 37.677*nu + 181.61*nu2)/(1 - 16.082*nu + 89.836*nu2);

    /* l=2, m=2 */
    b1 = -1.323643;
    b2 = -3.555007;
    b3 = 7.011267;
    b4 = 32.737824;
    c3phi[k22] = 3.8436 + 0.71565*nu + (5.12794 + b1*X12)*Shat + (9.9136 + b2*X12)*Shat2 + (-4.1075 + b3*X12)*Shat3 +(-31.5562 + b4*X12)*Shat4;
    /*
      b1 = -0.462321;
      b2 = -0.904512;
      b3 =  0.437747;
      b4 =  1.8275;
      double c3phi_nu     =  0.4558467286*nu + 3.8883812141;
      double c3phi_equal  =  (b4*X12 - 2.0575868122)*16*Shat4 +(b3*X12 - 0.5051534498)*8*Shat3 +(b2*X12 + 2.5742292762)*4*Shat2 +(b1*X12 + 2.5599640181)*2*Shat;
      c3phi[k22]          = c3phi_nu + c3phi_equal;
    */
      
    /* (l=3, m=3)*/
    b1 = -0.634377;
    b2 = 5.983525; 
    b3 = -5.881900;
    c3phi[k33] = 3.0611 - 6.1597*nu;// + (b1 + b2*X12 + b3*X12_2)*Shat;

    /* (l=4, m=4)*/
    b1 = 7.911653;
    b2 = 21.181688;
    c3phi[k44] = (3.6662 - 30.072*nu +76.371*nu2)/(1 - 3.5522*nu);// + (-4.9184 + b1*X12)*Shat + (-15.6772 + b2*X12)*Shat2;

    /* (l=5, m=5)*/
    b1 = -1.510167;
    b2 = +30.569461;
    b3 = -2.687133;
    b4 = +4.873750;
    b5 = -14.629684;
    b6 = +19.696954;
    //c3phi[k55] = b1 + b2*nu + (b3 + b4*X12)*Shat + (b5 + b6*X12)*Shat2;
    //c3phi[k55] = 373.312597*nu2 -59.69284 *nu+ 4.226238 + (b3 + b4*X12)*Shat + (b5 + b6*X12)*Shat2;
    c3phi[k55] = 373.3125969012880319*nu2 - 59.6928356312470854*nu + 4.2262379605181648 + (b3 + b4*X12)*Shat + (b5 + b6*X12)*Shat2;

  // c4phi
    /* l=2, m=1 */
    c4phi[k21] = (4.355 - 53.763*nu + 188.06*nu2)/(1 - 18.427*nu + 147.16*nu2);

    /* l=2, m=2 */
    b1 = 0.779683;
    b2 = -0.069638;
    c4phi[k22] = 1.4736 + 2.2337*nu + (8.26539 + b1*X12)*Shat + (14.2053 + b2*X12)*Shat2;
    /*
      b1 = -0.449976;
      b2 = -0.980913;
      double c4phi_nu     =  2.0822327682 * nu + 1.4996868401;
      double c4phi_equal  =  (b2*X12 + 3.5695199109)*4*Shat2 + (b1*X12 + 4.1312404030)*2*Shat;
      c4phi[k22]          =  c4phi_nu + c4phi_equal;
    */
      
    /* (l=3, m=3)*/
    b1 = -3.877528;
    b2 = 12.043300; 
    b3 = -6.524665;
    c4phi[k33] = 1.789 - 5.6684*nu;// + (b1 + b2*X12 + b3*X12_2)*Shat;

    /* (l=4, m=4)*/
    b1 = 11.746452;
    b2 = 34.922883;
    c4phi[k44] = 0.21595 + 23.216*nu;// + (-3.4207 + b1*X12)*Shat + (-15.5383 + b2*X12)*Shat2;

    /* (l=5, m=5)*/
    b1 = -1.383721;
    b2 = +56.871881;
    b3 = +7.198729;
    b4 = -3.870998;
    b5 = -25.992190;
    b6 = +36.882645;
    //c4phi[k55] = b1 + b2*nu + (b3 + b4*X12)*Shat + (b5 + b6*X12)*Shat2;
    //c4phi[k55] = 14.9111373110275380 *nu  +  1.3639723340485870+ (b3 + b4*X12)*Shat + (b5 + b6*X12)*Shat2;
    c4phi[k55] = 14.9111373110275380*nu + 1.3639723340485870 + (b3 + b4*X12)*Shat + (b5 + b6*X12)*Shat2;

  /* (l=3, m=2)*/
    c3A[3]   = (0.1877 - 3.0017*nu + 19.501*nu2)/(1 - 1.8199*nu) - exp(-703.67*(nu - 2./9.)*(nu - 2./9.));
    c3phi[3] = (0.90944 - 1.8924*nu + 3.6848*nu2)/(1 - 8.9739*nu + 21.024*nu2);
    c4phi[3] = (2.3038 - 50.79*nu + 334.41*nu2)/(1 - 18.326*nu + 99.54*nu2);
	    
  
  if (DEQUAL(nu,0.25,1e-9) && DEQUAL(chi1,chi2,1e-9)){
    modeon[0] = modeon[2] = modeon[4] = modeon[5] = modeon[7] = modeon[13] = 0;
  }

  double alpha21[KMAX], alpha1[KMAX], omega1[KMAX];
  QNM_coefs(af, alpha21, alpha1, omega1);  
  
  for (int k=0; k<KMAX; k++) {
    if (modeon[k]) {
      sigmar[k] = alpha1[k];
      sigmai[k] = omega1[k];
      Domg[k] 	= omega1[k] - Mbh*omgmrg[k];
    }
  }

  for (int k=0; k<KMAX; k++) {
    if (modeon[k]) {
      c2A[k] = 0.5*alpha21[k];
      double cosh_c3A= cosh(c3A[k]);  
      ca1[k] = Amrg[k]*alpha1[k]*cosh_c3A*cosh_c3A/c2A[k];
      ca2[k] = c2A[k];
      ca3[k] = c3A[k];
      ca4[k] = Amrg[k] - ca1[k]*tanh(c3A[k]);
      cb2[k] = alpha21[k];
      cb3[k] = c3phi[k];
      cb4[k] = c4phi[k];
      cb1[k] = Domg[k]*(1+c3phi[k]+c4phi[k])/(cb2[k]*(c3phi[k] + 2*c4phi[k]));
    }
  }
  
}

/** 
 *  Function: postpeak_coef
 *  -----------------------
 *   \details Equivalent to QNMHybridFitCab_BHNS_HM but for tidal
 *    disruption cases.
 *   From Gonzalez+2022 and updated for Gonzalez+ in prep.
 *
 *   @param[out] ca1,ca2,ca3,ca4,cb1,cb2,cb3,cb4 : phenomenological parameters of the postpeak template
 *   @param[in] nu                               : symmetric mass ratio
 *   @param[in] X1                               : mass fraction of body 1
 *   @param[in] X2                               : mass fraction of body 2
 *   @param[in] chi1                             : dimensionless spin of body 1
 *   @param[in] chi2                             : dimensionless spin of body 2
 *   @param[in] X1                               : dimensionless spin of body 2
 *   @param[in] X2                               : dimensionless spin of body 2
 *   @param[in] aK                               : a1 + a2
 *   @param[in] Mbh                              : final mass of the remnant BH
 *   @param[in] abh                              : final spin of the remnant BH
 *   @param[in] Apeak                            : dimensionless spin of body 2
 *   @param[in] alpha2                           : dimensionless spin of body 2
 *   @param[out] sigmar,sigmai                   : real and imaginary part of QNM frequency
 *
 */
void postpeak_coef(double *ca1, double *ca2, double *ca3, double *ca4, double *cb1, double *cb2, double *cb3, double *cb4, 
                    double *sigmar, double *sigmai, double nu, double chi1, double chi2, double X1, double X2, double aK, double Mbh, double abh,
                    double *Apeak, double *alpha2)
{
  const double a12   = X1*chi1 - X2*chi2;
  const double X12   = X1 - X2;
  const double Shat  = 0.5*(aK + a12*X12);
  const double Sbar  = 0.5*(a12 + aK*X12);
  const double af    = abh; 
  const double nu2   = SQ(nu);
  const double nu3   = nu2*nu;
  const double nu4   = SQ(nu2);
  const double af2   = SQ(af);
  const double af3   = af2*af;
  const double Shat2 = SQ(Shat); 
  const double Shat3 = Shat2*Shat;
  const double Shat4 = SQ(Shat2);
  const double Sbar2 = SQ(Sbar);
  const double X12_2 = SQ(X12);
  const double aK2   = SQ(aK);
  const double aeff  = aK + 1./3.*a12*X12;
  const double aeff2 = SQ(aeff);
  const double aeff3 = aeff2*aeff;
  const double lambda = EOBPars->LambdaBl2;
  
  double c3A[KMAX], c3phi[KMAX], c4phi[KMAX], c2A[KMAX];

  int modeon[KMAX];
  const int k21 = 0;
  const int k22 = 1; // (2,2) mode
  const int k32 = 3;
  const int k33 = 4;
  const int k44 = 8;
  const int k55 = 13;

  for (int k=0; k<KMAX; k++) {
    modeon[k] = 0; /* off */
    sigmar[k] = sigmai[k] = 0.;
    ca1[k] = ca2[k] = ca3[k] = ca4[k] = 0.;
    cb1[k] = cb2[k] = cb3[k] = cb4[k] = 0.;
  }

  modeon[1]  = 1;

  double Opeak[KMAX], Domega[KMAX];
  peak_bhns(nu, lambda, chi1, X1, X2, abh, Apeak, Opeak);

  double alpha1[KMAX], omega1[KMAX], omega2[KMAX], alpha21[KMAX];
  QNM_bhns_td(abh, alpha1, alpha2, omega1, omega2, alpha21, lambda, nu, chi1);

  /** l=2, m=2 **/
  //a3[1] = -0.56187 + 0.75497*nu;
  //b3[1] = ( 4.4414 - 63.107*nu + 296.64*nu*nu ) / ( 1. - 13.299*nu + 69.129*nu*nu );
  //b4[1] = ( 7.1508 - 109.47*nu ) / ( 1. + 556.34*nu + 287.42*nu*nu );

   double b1, b2, b3, b4, b5, b6, b7, b8, b9;

  // c3A
    /* l=2, m=1 */
    c3A[k21]   = (0.23882 - 2.2982*nu + 5.7022*nu2)/(1 - 7.7463*nu + 27.266*nu2);
    
    /* l=2, m=2 */
    b1 = 0.1659421;
    b2 = -0.2560047;
    b3 = -0.9418946;
    c3A[k22] = -0.5585 + 0.81196*nu + (-0.398576+b1*X12)*Shat + (0.099805+b2*X12)*Shat2 + (0.72125+b3*X12)*Shat3;

    /* l=3, m=3 */
    b1 = -0.3502608;
    b2 = 1.587606;
    b3 = -1.555325;
    c3A[k33] = -0.41455 + 1.3225*nu + (b1 + b2*X12 + b3*X12_2)*a12;
    c3A[k32] = c3A[k33];

    /* l=4, m=4 */
    b1 = -9.614738;
    b2 = 122.461125;
    c3A[k44] = -0.41591 + 3.2099*nu + b1*nu*Shat + b2*nu*Shat2;

    /* l=5, m=5 */
    b1 = -7.063079;
    b2 = +65.464944;
    b3 = -2.055335;
    b4 = -0.585373;
    b5 = -12.631409;
    b6 = +19.271346;
    c3A[k55] = b1 + b2*nu + (b3 + b4*X12)*a12 + (b5 + b6*X12)*SQ(a12);
     
  // c3phi
    /* l=2, m=1 */
    c3phi[k21] = (2.6269 - 37.677*nu + 181.61*nu2)/(1 - 16.082*nu + 89.836*nu2);

    /* l=2, m=2 */
    b1 = -1.323643;
    b2 = -3.555007;
    b3 = 7.011267;
    b4 = 32.737824;
    c3phi[k22] = 3.8436 + 0.71565*nu + (5.12794 + b1*X12)*Shat + (9.9136 + b2*X12)*Shat2 + (-4.1075 + b3*X12)*Shat3 +(-31.5562 + b4*X12)*Shat4;

    /* l=3, m=3 */
    b1 = -0.634377;
    b2 = 5.983525; 
    b3 = -5.881900;
    c3phi[k33] = 3.0611 - 6.1597*nu;
    c3phi[k32] = c3phi[k33];

    /* l=4, m=4 */
    b1 = 7.911653;
    b2 = 21.181688;
    c3phi[k44] = (3.6662 - 30.072*nu +76.371*nu2)/(1 - 3.5522*nu);

    /* l=5, m=5 */
    b1 = -1.510167;
    b2 = +30.569461;
    b3 = -2.687133;
    b4 = +4.873750;
    b5 = -14.629684;
    b6 = +19.696954;
    c3phi[k55] = b1 + b2*nu + (b3 + b4*X12)*Shat + (b5 + b6*X12)*Shat2;

  // c4phi
    /* l=2, m=1 */
    c4phi[k21] = (4.355 - 53.763*nu + 188.06*nu2)/(1 - 18.427*nu + 147.16*nu2);

    /* l=2, m=2 */
    b1 = 0.779683;
    b2 = -0.069638;
    c4phi[k22] = 1.4736 + 2.2337*nu + (8.26539 + b1*X12)*Shat + (14.2053 + b2*X12)*Shat2;

    /* l=3, m=3 */
    b1 = -3.877528;
    b2 = 12.043300; 
    b3 = -6.524665;
    c4phi[k33] = 1.789 - 5.6684*nu;
    c4phi[k32] = c4phi[k33];


    /* l=4, m=4 */
    b1 = 11.746452;
    b2 = 34.922883;
    c4phi[k44] = 0.21595 + 23.216*nu;

    /* l=5, m=5 */
    b1 = -1.383721;
    b2 = +56.871881;
    b3 = +7.198729;
    b4 = -3.870998;
    b5 = -25.992190;
    b6 = +36.882645;
    c4phi[k55] = b1 + b2*nu + (b3 + b4*X12)*Shat + (b5 + b6*X12)*Shat2;

	    
  /* from Nagar et al arxiv:2001.09082 */
  for (int k=0; k<KMAX; k++) {
    sigmar[k] = alpha1[k];
    sigmai[k] = omega1[k];
    Domega[k] = omega1[k]-Mbh*Opeak[k];     
    /** These are the same for all modes given the constraints
     * for physically correct behaviour of the fits */
      c2A[k] = 0.5*alpha21[k];
      double cosh_c3A= cosh(c3A[k]);  
      ca1[k] = Apeak[k]*alpha1[k]*cosh_c3A*cosh_c3A/c2A[k];
      ca2[k] = c2A[k];
      ca3[k] = c3A[k];
      ca4[k] = Apeak[k] - ca1[k]*tanh(c3A[k]);
      cb2[k] = alpha21[k];
      cb3[k] = c3phi[k];
      cb4[k] = c4phi[k];
      cb1[k] = Domega[k]*(1+c3phi[k]+c4phi[k])/(cb2[k]*(c3phi[k] + 2*c4phi[k]));
  }
}

/** 
 *  Function: peak_bhns
 *  -----------------------
 *   \details ^{peak}_{22} amplitude and omega^{peak}_{22} frequency at tpeak (\dot{A}(tpeak)=0) 
 *   ap_bbh and op_bbh is for BBH from the function QNMHybridFitCab_HM in ..Fits.c
 *   From Gonzalez+2022 and updated for Gonzalez+ in prep.
 *
 *   @param[in] nu                               : symmetric mass ratio
 *   @param[in] lambda                           : tidal polarizability
 *   @param[in] chi1                             : dimensionless spin of body 1
 *   @param[in] X1                               : mass fraction of body 1
 *   @param[in] X2                               : mass fraction of body 2
 *   @param[in] abh                              : final spin of the remnant BH
 *   @param[out] Apeak                           : Peak amplitude
 *   @param[out] Opeak                           : Omega at the peak of amplitude
 */
void peak_bhns(double nu, double lambda, double chi1, double X1, double X2, double abh, double *Apeak, double *Opeak)
{
   double op_bbh[KMAX], ap_bbh[KMAX], Ap[KMAX], Op[KMAX];
  for (int k=0; k<KMAX; k++) {
    Apeak[k] = Opeak[k] = op_bbh[k] = ap_bbh[k] = 0.;
     Ap[k] = Op[k] = 1.;
  }

  // BBH fits 
    const double chi2 = 0.;
    const double aK = X1*chi1 + X2*chi2;
    const double a12c   = X1*chi1 - X2*chi2;
    const double X12   = X1 - X2;
    const double Shat  = 0.5*(aK + a12c*X12);
    const double Sbar  = 0.5*(a12c + aK*X12);
    const double af    = abh; 
    const double nu2   = SQ(nu);
    const double nu3   = nu2*nu;
    const double nu4   = SQ(nu2);
    const double af2   = SQ(af);
    const double af3   = af2*af;
    const double Shat2 = SQ(Shat); 
    const double Shat3 = Shat2*Shat;
    const double Shat4 = SQ(Shat2);
    const double Sbar2 = SQ(Sbar);
    const double X12_2 = SQ(X12);
    const double aK2   = SQ(aK);
    const double aeff  = aK + 1./3.*a12c*X12;
    const double aeff2 = SQ(aeff);
    const double aeff3 = aeff2*aeff;

    double ATP[KMAX], omgTP[KMAX];
    for (int k=0; k<KMAX; k++) {
    ATP[k] = omgTP[k] = 0.;
    }
  
    ATP[0]    = 0.5238781992;
    ATP[1]    = 1.44959;
    ATP[3]    = 0.1990192432;
    ATP[4]    = 0.5660165890;
    ATP[8]    = 0.2766182761;
    ATP[13]   = 0.151492;
    omgTP[0]  = 0.2906425497;
    omgTP[1]  = 0.273356;
    omgTP[3]  = 0.4516072248;
    omgTP[4]  = 0.4541278937;
    omgTP[8]  = 0.6356586393;
    omgTP[13] = 0.818117;

     /* a_vars */
    double a1, a2, a3;
    double a1Omg, a2Omg, a3Omg, a4Omg;
    double a1Amp, a2Amp, a3Amp, a1AmpS, a2AmpS, a3AmpS;
    
    /* b_vars */
    double b1, b2, b3, b4, b5, b6, b7, b8, b9;
    double b1Omg, b2Omg, b3Omg, b4Omg;
    double b11Omg, b12Omg, b13Omg, b14Omg, b21Omg, b22Omg, b23Omg;
    double b24Omg, b31Omg, b32Omg, b33Omg, b34Omg, b41Omg, b42Omg; 
    double b43Omg, b44Omg, b3Amp, b4Amp;
    double b1Amax, b2Amax, b3Amax, b1Amp, b2Amp, b11Amp, b12Amp, b13Amp, b14Amp, b21Amp, b22Amp, b23Amp, b24Amp, b31Amp, b32Amp, b33Amp, b34Amp;
    
    /* c_vars */
    double c1, c2, c3, c4;
    double c1Omg, c2Omg, c5Omg, c6Omg;
    double c11Omg, c12Omg, c13Omg, c21Omg, c22Omg, c23Omg, c31Omg, c32Omg, c33Omg, c41Omg, c42Omg, c43Omg;
    double c1Amax, c2Amax, c3Amax, c4Amax, c1Amp, c2Amp, c3Amp, c4Amp, c5Amp, c6Amp;
    double c11Amp, c12Amp, c13Amp, c14Amp, c21Amp, c22Amp, c23Amp, c24Amp;
    double c31Amp, c32Amp, c41Amp, c42Amp;
    
    /* others vars */
    double n1Omg, n2Omg, d1Omg, d2Omg;
    double omg1, omg2, orb, num, denom, OmgOrb, Omgspin, omgOrb, omgS, omgorb, omgspin;
    double Aorb, Aspin, scale, Amax1, Amax2, num_A, denom_A;

    /* l=2, m=1 */
    if((chi1==0)&&(chi2==0)){
      ap_bbh[0] = ATP[0]*X12*(1+9.0912*nu+3.9331*nu2)/(1+11.108*nu);
      op_bbh[0]  = omgTP[0]*(1-0.060432*nu+1.9995*nu2)/(1+0.23248*nu);
      
    }else{
      a1Omg     = -0.563075;
      a2Omg     = 3.28677;
      b1Omg     = 0.179639;
      b2Omg     = -0.302122;
      c11Omg    = -1.20684;
      c21Omg    = 0.425645;
	
      omgOrb    = omgTP[0]*(1 + a1Omg*nu + a2Omg*nu2);
      b1        = b1Omg + c11Omg*nu;
      b2        = b2Omg + c21Omg*nu;
      omgS      = 1 + b1*Shat + b2*Shat2;
      op_bbh[0] = omgOrb*omgS;

      b1    = +0.891139;
      b2 	  = -5.191702;
      b3 	  = +3.480139;
      b4 	  = +10.237782;
      b5 	  = -13.867475; 
      b6 	  = +10.525510;

      if (DEQUAL(nu,0.25,1e-9)) {
        double S_bar21 = - fabs(Sbar);
        Aspin = ((-0.4281863 + b1*nu + b2*nu2)*S_bar21 + (-0.335659 + b3*nu + b4*nu2)*S_bar21*S_bar21)/(1 + (+0.828923 + b5*nu + b6*nu2)*S_bar21);
      } else {
        Aspin = ((-0.4281863 + b1*nu + b2*nu2)*Sbar + (-0.335659 + b3*nu + b4*nu2)*Sbar2)/(1 + (+0.828923 + b5*nu + b6*nu2)*Sbar);
      }
      Aorb    = ATP[0]*X12*(1 + 3.3362232268*nu + 3.4708521429*nu2)/(1 + 4.7623643259*nu);
      ap_bbh[0] = Aorb + Aspin;
    }
    
    /* l=2, m=2 */
    c2 = -0.122735;
    c1 = 0.0857478;
    c4 = -0.0760023;
    c3 = 0.0826514;
    omg1 = (c2*X12*X12 + c1*X12 -0.1416002395)*2.*Shat + 1;
    omg2 = (c4*X12*X12 + c3*X12 -0.3484804901)*2.*Shat + 1;
    op_bbh[1] = (0.481958619443355*nu*nu +0.223976694441952*nu +0.273813064427363)*omg1/omg2;

    scale = 1 - 0.5*aeff*op_bbh[1];
    c4Amax = -0.0820894;
    c3Amax = 0.176126;
    c2Amax = -0.150239;
    c1Amax = 0.20491;   
    num_A    = (c4Amax*X12*X12 + c3Amax*X12 - 0.2935238329)*aeff + 1.;
    denom_A  = 1. + (c2Amax*X12*X12 + c1Amax*X12 -0.4728707630)*aeff;   
    Aorb = 1.826573640739664*nu*nu + 0.100709438291872*nu + 1.438424467327531;
    ap_bbh[1] = Aorb*scale*(num_A/denom_A);

    /* l=3, m=2 */
    if((chi1==0)&&(chi2==0)&&(X12==0)){
      ap_bbh[3]    = ATP[3]*(1-3*nu)*(1 - 6.142*nu + 11.372*nu2)/(1 - 3.6448*nu);
      op_bbh[3]  = omgTP[3]*(1 - 9.0214*nu + 21.078*nu2)/(1 - 8.6636*nu + 19.493*nu2);
    }else{
      a1Omg            = -9.13525;
    a2Omg            = 21.488;
    a3Omg     = -8.81384;
    a4Omg     = 20.0595; 
    b1Omg            = -0.458126;
    b2Omg            = 0.0474616;
    b3Omg     = -0.486049;
    c11Omg           = 3.25319; 
    c12Omg    = 0.535555; 
    c13Omg    = -8.07905; 
    c21Omg           = 1.00066; 
    c22Omg    = -1.1333; 
    c23Omg    = 0.601572;
    
    b1        = (b1Omg + c11Omg*X12 + c12Omg*X12_2)/(1+c13Omg*X12);
    b2        = (b2Omg + c21Omg*X12 + c22Omg*X12_2)/(1+c23Omg*X12);
    b3        = b3Omg;
    omgOrb    = omgTP[3]*(1 + a1Omg*nu + a2Omg*nu2)/(1 + a3Omg*nu + a4Omg*nu2);
    omgS      = (1 + b1*aK + b2*aK2)/(1 + b3*aK);
    op_bbh[3] = omgOrb*omgS;

    a1      = -6.06831;
    a2      = 10.7505;    
    a3      = -3.68883; 
    b1      = -0.258378;
    b2      = 0.679163;
    c11Amp  = 4.36263; 
    c12Amp  = -12.5897; 
    c13Amp  = -7.73233; 
    c14Amp  = 16.2082; 
    c21Amp  = 3.04724; 
    c22Amp  = 46.5711;
    c23Amp  = 2.10475; 
    c24Amp  = 56.9136;    
    b1      = (b1 + c11Amp*nu + c12Amp*nu2)/(1 + c13Amp*nu + c14Amp*nu2);
    b2      = (b2 + c21Amp*nu + c22Amp*nu2)/(1 + c23Amp*nu + c24Amp*nu2);
 
    scale   = 1 + aK*pow(fabs(op_bbh[3])/2.,1./3.);
    Aorb    = ATP[3]*(1-3*nu)*(1 + a1*nu + a2*nu2)/(1 + a3*nu);
    Aspin   = (1 + b1*aK)/(1 + b2*aK);
    ap_bbh[3] = Aorb*scale*Aspin; 
    }
    

    /* l=3, m=3 */
    if((chi1==0)&&(chi2==0)&&(X12==0)){
      ap_bbh[4]    = ATP[4]*X12*(1 + 0.098379*nu + 3.8179*nu2);
      op_bbh[4]  = omgTP[4]*(1 + 1.1054*nu + 2.2957*nu2);
    }else{
      a1Omg     = 1.08224;
      a2Omg     = 2.59333;
      b1Omg     = -0.406161;
      b2Omg     = -0.0647944;
      b3Omg     = -0.748126;
      c1Omg     = 0.85777;
      c2Omg     = -0.70066;
      c5Omg     = 2.97025;
      c6Omg     = -3.96242;
      b1Omg     = (b1Omg + c1Omg*nu)/(1 + c2Omg*nu);
      b3Omg     = (b3Omg + c5Omg*nu)/(1 + c6Omg*nu);

      omgorb    = 1 + a1Omg*nu + a2Omg*nu2;
      omgspin   = (1 + b1Omg*Shat + b2Omg*Shat2)/(1 + b3Omg*Shat);
      op_bbh[4] = omgTP[4]*omgorb*omgspin/omgTP[4];

      a1Amp   = -0.22523;
      a2Amp   = 3.0569;
      a3Amp   = -0.396851; 
      b1Amp   = 0.100069;
      b2Amp   = -0.455859;
      c1Amp   = -0.401156;
      c2Amp   = -0.141551;
      c3Amp   = -15.4949; 
      c4Amp   = 1.84962;
      c5Amp   = -2.03512;
      c6Amp   = -4.92334;
      b1Amp   = (b1Amp + c1Amp*nu)/(1 + c2Amp*nu + c3Amp*nu2);
      b2Amp   = (b2Amp + c4Amp*nu)/(1 + c5Amp*nu + c6Amp*nu2);

      if (DEQUAL(nu,0.25,1e-9)) {
        double a12_33 = fabs(a12c);
        Aspin   = (b1Amp*a12_33)/(1 + b2Amp*a12_33);
      } else {
        Aspin   = (b1Amp*a12c)/(1 + b2Amp*a12c);
      }
      Aorb    = (1 + a1Amp*nu + a2Amp*nu2)/(1 + a3Amp*nu);
      ap_bbh[4] = ATP[4]*X12*Aorb + Aspin;
    }
    

    /* l=4, m=4 */
    n1Omg     = -0.964614;
    n2Omg     = -11.1828;
    d1Omg     = -2.08471;
    d2Omg     = -6.89287;
    b1Omg     = -0.445192;
    b2Omg     = -0.0985658;
    b3Omg     = -0.0307812;
    b4Omg     = -0.801552;
    b11Omg    = -0.92902;
    b12Omg    = 10.86310;
    b13Omg    = -4.44930;
    b14Omg    = 3.01808;
    b21Omg    = 0;
    b22Omg    = 1.62523;
    b23Omg    = -7.70486;
    b24Omg    = 15.06517;
    b31Omg    = 0;
    b32Omg    = 0;
    b33Omg    = 0;
    b34Omg    = 0;
    b41Omg    = 0.93790;
    b42Omg    = 8.36038;
    b43Omg    = -4.85774;
    b44Omg    = 4.80446;
    a1Omg     = (b1Omg + b11Omg*nu + b12Omg*nu2)/(1 + b13Omg*nu + b14Omg*nu2);
    a2Omg     = (b2Omg + b21Omg*nu + b22Omg*nu2)/(1 + b23Omg*nu + b24Omg*nu2);
    a3Omg     = (b3Omg + b31Omg*nu + b32Omg*nu2)/(1 + b33Omg*nu + b34Omg*nu2);
    a4Omg     = (b4Omg + b41Omg*nu + b42Omg*nu2)/(1 + b43Omg*nu + b44Omg*nu2);

    OmgOrb    = (1 + n1Omg*nu + n2Omg*nu2)/(1 + d1Omg*nu + d2Omg*nu2);
    Omgspin   = (1 + a1Omg*Shat + a2Omg*Shat2 + a3Omg*Shat3)/(1 + a4Omg*Shat);
    op_bbh[8] = omgTP[8]*OmgOrb*Omgspin;

    a1Amp  = -3.7082;
    a2Amp  = 0.280906;
    a3Amp  = -3.71276; 
    b1Amp  = -0.316647;
    b2Amp  = -0.062423;
    b3Amp  = -0.852876;
    b11Amp = 1.2436;
    b12Amp = -1.60555;
    b13Amp = -4.05685; 
    b14Amp = 1.59143;
    b21Amp = 0.837418;
    b22Amp = -2.93528; 
    b23Amp = -11.5591;
    b24Amp = 34.1863;
    b31Amp = 0.950035; 
    b32Amp = 7.95168;
    b33Amp = -1.26899;
    b34Amp = -9.72147;

    a1AmpS = (b1Amp + b11Amp*nu + b12Amp*nu2)/(1 + b13Amp*nu + b14Amp*nu2);
    a2AmpS = (b2Amp + b21Amp*nu + b22Amp*nu2)/(1 + b23Amp*nu + b24Amp*nu2);
    a3AmpS = (b3Amp + b31Amp*nu + b32Amp*nu2)/(1 + b33Amp*nu + b34Amp*nu2);

    scale   = 1 - 0.5*Shat*op_bbh[8];
    Aorb    = (1-3*nu)*(1 + a1Amp*nu + a2Amp*nu2)/(1 + a3Amp*nu);
    Aspin   = (1 + a1AmpS*Shat + a2AmpS*Shat2)/(1 + a3AmpS*Shat);
    ap_bbh[8] = ATP[8]*Aorb*Aspin*scale;

    /* l=5, m=5 */
    if((chi1==0)&&(chi2==0)&&(X12==0)){
      ap_bbh[13]   = ATP[13]*X12*(1 - 2*nu)*(1 - 0.29628*nu + 6.4207*nu2);
      op_bbh[13] = omgTP[13]*(1 - 2.8918*nu - 3.2012*nu2)/(1 - 3.773*nu);
    }else{
      b1 	       = +1.487294;
    b2 	       = -2.058537;
    b3         = +1.454248;
    b4         = -1.301284;
    OmgOrb     = omgTP[13]*(1 - 2.8918*nu - 3.2012*nu2)/(1 - 3.773*nu);
    Omgspin    = (1 + (-0.332703 + b1*nu)/(1 + b2*nu)*Shat)/(1 + (-0.675738 + b3*nu)/(1 + b4*nu)*Shat);
    op_bbh[13] = OmgOrb*Omgspin;

    b1 = +5.720690;
    b2 = +44.868515;
    b3 = +12.777090;
    b4 = -42.548247;

    Aorb     = ATP[13]*sqrt(1-4*nu)*(1-2*nu)*(1 - 0.29628*nu + 6.4207*nu2);
    Aspin    = (0.04360530/(1 + b1*nu + b2*nu2)*a12c)/(1 - 0.5769451/(1+b3*nu+b4*nu2)*a12c);
    ap_bbh[13] = Aorb + Aspin;
    }

  /* BHNS Section */

  /* Amplitude peak */
  double c102,c101,c100,c112,c111,c110,c122,c121,c120,c202,c201,c200,c212,c211,c210,c222,c221,c220,c31,c32;
  double c301,c300,c311,c310,c132,c131,c130,c232,c231,c230,c321,c320;
  double c113,c123,c213,c223,c313,c312,c323,c322,c411,c410,c421,c420;
  double c114,c124,c134,c133,c214,c224,c234,c233,c332,c331,c330;
  double p10,p11,p12,p20,p21,p22,p30,p31,p13,p23,p32,p41,p42,p33;
  double p1,p2,p3,p4;

  /* l=2, m=1 */
  c112 = 79.69190135775388;
  c111 = -63.833579388191794;
  c110 = 475.62792687058675;
  c122 = 13717.304556430048;
  c121 = -16696.974173835133;
  c120 = 187.6302117256296;
  c132 = -47483.566270772106;
  c131 = 55908.56848609845;
  c130 = -3729.579436870687;
  c212 = -3.675245190907651;
  c211 = 4.055153660963919;
  c210 = -1.810890290356132;
  c222 = 45.410498599979725;
  c221 = -41.879957729349584;
  c220 = 16.088785963718603;
  c232 = -139.34673532789955;
  c231 = 113.28977369273983;
  c230 = -36.645761242568085;
  c311 = -72.42429070293126;
  c310 = 44.48560453945776;

  p11 = c112*chi1*chi1 + c111*chi1 + c110;
  p12 = c122*chi1*chi1 + c121*chi1 + c120;
  p13 = c132*chi1*chi1 + c131*chi1 + c130;
  p21 = c212*chi1*chi1 + c211*chi1 + c210;
  p22 = c222*chi1*chi1 + c221*chi1 + c220;
  p23 = c232*chi1*chi1 + c231*chi1 + c230;
  p31 = c311*chi1 + c310;

  p1 = p11*nu + p12*nu*nu + p13*nu*nu*nu;
  p2 = p21*nu + p22*nu*nu + p23*nu*nu*nu;
  p3 = p31*nu;

  Ap[0] =  ( 1 + (p1*lambda + p2*lambda*lambda ) ) / ( (1 + lambda*p3*p3 ) );
  if(Ap[0]<0.0){Ap[0]=1.0;};

  /* l=2, m=2 */
  c102 = 0.025960639049437614; 
  c101 = 0.011583003041680442;
  c100 = -0.0016533947810772696;
  c112 = -0.31403656018939863;
  c111 = -0.12825460441414602;
  c110 = 0.13396228694902;
  c122 = 0.8917717620733512;
  c121 = 0.3150355510141175;
  c120 = -0.513304355534311;
  c202 = -1.984228201511826e-05;
  c201 = -2.804785772573806e-06;
  c200 = 8.018661470565014e-06;
  c212 = 0.0002223708961923682;
  c211 = 2.555424239372415e-05;
  c210 = -0.00010025951516653183;
  c222 = -0.0006022803171412016;
  c221 = -5.057203798367972e-05;
  c220 = 0.0002870181683963824;
  c31  = 0.009579950459173342;
  c32  = -0.01899088523426322;
  
  p10 = c102*chi1*chi1 + c101*chi1 + c100;
  p11 = c112*chi1*chi1 + c111*chi1 + c110;
  p12 = c122*chi1*chi1 + c121*chi1 + c120;
  p20 = c202*chi1*chi1 + c201*chi1 + c200;
  p21 = c212*chi1*chi1 + c211*chi1 + c210;
  p22 = c222*chi1*chi1 + c221*chi1 + c220;

  p1 = p10*nu + p11*nu*nu + p12*nu*nu*nu;
  p2 = p20*nu + p21*nu*nu + p22*nu*nu*nu;
  p3 = c31*nu + c32*nu*nu;

  Ap[1] =  ( 1 + (p1*lambda + p2*lambda*lambda) ) / ( 1 + lambda*p3  );
  if(Ap[1]<0.0){Ap[1]=1.0;};

  /* l=3, m=2 */
  c112 = 20.16870439462333;
  c111 = 11.038387003596617;
  c110 = -4.020843926856486;
  c122 = -207.72480245938715;
  c121 = -117.37710746118753;
  c120 = 40.66818784611969;
  c132 = 538.5544958151388;
  c131 = 318.28911903762025;
  c130 = -99.26273990634219;
  c212 = -0.016627367594463694;
  c211 = -0.014305531331241536;
  c210 = 0.004174534399221245;
  c222 = 0.15432035038868333;
  c221 = 0.14118798342477987;
  c220 = -0.040988466644647094;
  c232 = -0.3581036182848496;
  c231 = -0.35040084466576926;
  c230 = 0.0994625770255457;
  c312 = 2.05494620204673e-06;
  c311 = 4.055461359375113e-06;
  c310 = -6.611853315606744e-07;
  c322 = -1.0313837453840835e-05;
  c321 = -3.6352321950711845e-05;
  c320 = 5.58113544950673e-06;
  c332 = 2.421279243476082e-06;
  c331 = 8.114351204003656e-05;
  c330 = -1.1677202522832645e-05;
  c411 = -3.317069910114291;
  c410 = -0.015174553789524594;
  c421 = 22.36565160345709;
  c420 = 3.7427983470741455;

  p11 = c112*chi1*chi1 + c111*chi1 + c110;
  p12 = c122*chi1*chi1 + c121*chi1 + c120;
  p13 = c132*chi1*chi1 + c131*chi1 + c130;
  p21 = c212*chi1*chi1 + c211*chi1 + c210;
  p22 = c222*chi1*chi1 + c221*chi1 + c220;
  p23 = c232*chi1*chi1 + c231*chi1 + c230;
  p31 = c312*chi1*chi1 + c311*chi1 + c310;
  p32 = c322*chi1*chi1 + c321*chi1 + c320;
  p33 = c332*chi1*chi1 + c331*chi1 + c330;

  p41 = c411*chi1 + c410;
  p42 = c421*chi1 + c420;

  p1 = p11*nu + p12*nu*nu + p13*nu*nu*nu;
  p2 = p21*nu + p22*nu*nu + p23*nu*nu*nu;
  p3 = p31*nu + p32*nu*nu + p33*nu*nu*nu;
  p4 = p41*nu + p42*nu*nu;

  Ap[3] = ( 1 + (p1*lambda + p2*lambda*lambda + p3*lambda*lambda*lambda ) ) / ( (1 + lambda*p4*p4 ) );
  if(Ap[3]<0.0){
    Ap[3]=fabs(Ap[3]);
    if(Ap[3]>1.0){Ap[3]=1.0;}
    }

  /* l=3, m=3 */
  c114 = 98.58316168624746;
  c113 = 23.34611781664497;
  c112 = -28.887092039772345;
  c111 = 1.4659939569547824;
  c110 = 0.15438882520561592;
  c124 = -1023.2487147730186;
  c123 = -248.719272283286;
  c122 = 299.4472305323695;
  c121 = -14.353307233157919;
  c120 = -1.422004596271538;
  c134 = 2630.326171106412;
  c133 = 655.8187503615836;
  c132 = -769.2344295536586;
  c131 = 34.542048041903065;
  c130 = 3.243842025932237;
  c214 = -0.00023145853854142156;
  c213 = -8.571596581472541e-05;
  c212 = 7.124381821132649e-05;
  c211 = 1.3058237396225664e-05;
  c210 = 6.794241998472768e-07;
  c224 = 0.0024030529286666346;
  c223 = 0.0008885012949962715;
  c222 = -0.0007497610069806854;
  c221 = -0.000148744547249859;
  c220 = -6.704977980044648e-06;
  c234 = -0.006184049610372618;
  c233 = -0.002273091188747331;
  c232 = 0.001965296930731387;
  c231 = 0.0004187367371565314;
  c230 = 1.9367965673575747e-05;
  c311 = -0.054311467466514475;
  c310 = 0.009586523081306567;
  c321 = 0.2937082050407349;
  c320 = 0.009847652717538884;
    
  p11 = c114*chi1*chi1*chi1*chi1 + c113*chi1*chi1*chi1 + c112*chi1*chi1 + c111*chi1 + c110;
  p12 = c124*chi1*chi1*chi1*chi1 + c123*chi1*chi1*chi1 + c122*chi1*chi1 + c121*chi1 + c120;
  p13 = c134*chi1*chi1*chi1*chi1 + c133*chi1*chi1*chi1 + c132*chi1*chi1 + c131*chi1 + c130;
  p21 = c214*chi1*chi1*chi1*chi1 + c213*chi1*chi1*chi1 + c212*chi1*chi1 + c211*chi1 + c210;
  p22 = c224*chi1*chi1*chi1*chi1 + c223*chi1*chi1*chi1 + c222*chi1*chi1 + c221*chi1 + c220;
  p23 = c234*chi1*chi1*chi1*chi1 + c233*chi1*chi1*chi1 + c232*chi1*chi1 + c231*chi1 + c230;
  p31 = c311*chi1 + c310;
  p32 = c321*chi1 + c320;

  p1 = p11*nu + p12*nu*nu + p13*nu*nu*nu;
  p2 = p21*nu + p22*nu*nu + p23*nu*nu*nu;
  p3 = p31*nu + p32*nu*nu;       

  Ap[4] = ( 1 + (p1*lambda + p2*lambda*lambda*lambda ) ) / ( (1 + lambda*lambda*p3*p3 )*(1 + lambda*lambda*p3*p3 ) );
  if(Ap[4]>1.0){Ap[4]=1.0;};
  if(Ap[4]<0.0){
    Ap[4]=fabs(Ap[4]);
    if(Ap[4]>1.0){Ap[4]=1.0;}
    }

  /* l=4, m=4 */
  c112 = 3543.1799502856747;
  c111 = 3309.22676314805;
  c110 = 807.7348298148902;
  c122 = -40577.17930128383;
  c121 = -37096.373919250385;
  c120 = -8622.741990143393;
  c132 = 118155.67244830148;
  c131 = 106312.69486816866;
  c130 = 23776.282112663503;
  c212 = 0.26110920437332924;
  c211 = -5.627928587745965;
  c210 = -2.6550111083017383;
  c222 = -5.521939460220596;
  c221 = 56.79132247835583;
  c220 = 27.525128848923384;
  c232 = 19.288265323714658;
  c231 = -143.50684157009792;
  c230 = -70.87348176416717;
  c312 = -0.0016458349827221754;
  c311 = 0.002319476890473154;
  c310 = 0.0017785039216730385;
  c322 = 0.019459621840015342;
  c321 = -0.022402763254381477;
  c320 = -0.018327149646315592;
  c332 = -0.055329421236874586;
  c331 = 0.053749032597984206;
  c330 = 0.04672079813389738;
  c411 = -117.4017144826126;
  c410 = -10.046736035912739;
  c421 = 739.3626850513566;
  c420 = 127.27925047272853;

  p11 = c112*chi1*chi1 + c111*chi1 + c110;
  p12 = c122*chi1*chi1 + c121*chi1 + c120;
  p13 = c132*chi1*chi1 + c131*chi1 + c130;
  p21 = c212*chi1*chi1 + c211*chi1 + c210;
  p22 = c222*chi1*chi1 + c221*chi1 + c220;
  p23 = c232*chi1*chi1 + c231*chi1 + c230;
  p31 = c312*chi1*chi1 + c311*chi1 + c310;
  p32 = c322*chi1*chi1 + c321*chi1 + c320;
  p33 = c332*chi1*chi1 + c331*chi1 + c330;
  p41 = c411*chi1 + c410;
  p42 = c421*chi1 + c420;
    
  p1 = p11*nu + p12*nu*nu + p13*nu*nu*nu; 
  p2 = p21*nu + p22*nu*nu + p23*nu*nu*nu; 
  p3 = p31*nu + p32*nu*nu + p33*nu*nu*nu; 
  p4 = p41*nu + p42*nu*nu;

  Ap[8] = ( 1 + (p1*lambda + p2*lambda*lambda + p3*lambda*lambda*lambda ) ) / ( (1 + lambda*p4*p4 ) );
  if(Ap[8]<0.0){
    Ap[8]=fabs(Ap[8]);
    if(Ap[8]>1.0){Ap[8]=1.0;}
    }
  if(Ap[8]>1.0){Ap[8]=1.0;};

  /* Omega peak */
  double d102,d101,d100,d112,d111,d110,d122,d121,d120,d202,d201,d200,d212,d211,d210,d222,d221,d220,d0;
  double d113,d123,d213,d223,d311,d310,d321,d320,d131,d130,d231,d230,d331,d330,d411,d410,d421,d420;
  double d10,d11,d12,d20,d21,d22,d31,d32,d13,d23,d33,d41,d42;
  double d1,d2,d3,d4;

  /* l=2, m=1 */
  d111 = -15507583.54315245;
  d110 = 29233984.00747841;
  d121 = 188051611.3537233;
  d120 = -303459693.99696875;
  d131 = -561206949.6722764;
  d130 = 781093255.9339185;
  d211 = 65639.32689138538;
  d210 = -28574.823293025787;
  d221 = -715381.2300665657;
  d220 = 202374.5955228343;
  d231 = 1941160.48875847;
  d230 = -200599.97864086388;
  d311 = -80.15396589788257;
  d310 = 41.83106445091412;
  d321 = 869.1588059570613;
  d320 = -413.4870939347111;
  d331 = -2341.85742662733;
  d330 = 992.3696507060471;
  d411 = 35.013368813946784;
  d410 = -20.89552442415687;
  d421 = -174.81894201673433;
  d420 = 230.9684863107324;
    
  d11 = d111*chi1 + d110;
  d12 = d121*chi1 + d120;
  d13 = d131*chi1 + d130;
  d21 = d211*chi1 + d210;
  d22 = d221*chi1 + d220;
  d23 = d231*chi1 + d230;
  d31 = d311*chi1 + d310;
  d32 = d321*chi1 + d320;
  d33 = d331*chi1 + d330;
  d41 = d411*chi1 + d410;
  d42 = d421*chi1 + d420;

  d1 = d11*nu + d12*nu*nu + d13*nu*nu*nu;
  d2 = d21*nu + d22*nu*nu + d23*nu*nu*nu;
  d3 = d31*nu + d32*nu*nu + d33*nu*nu*nu;
  d4 = d41*nu + d42*nu*nu;

  Op[0] = ( 1 + (d1*lambda + d2*lambda*lambda + d3*lambda*lambda*lambda ) ) / ( (1 + lambda*d4*d4 )*(1 + lambda*d4*d4 ) );

  /* l=2, m=2 */
  d102 = 0.1107284045072136; 
  d101 = 0.030700021625216727;
  d100 = -0.06263403886363417;
  d112 = -1.4773209090362833;
  d111 = -0.23044230587158065;
  d110 = 0.9051355077601556;
  d122 = 4.438909358861151;
  d121 = 0.21104611999758255;
  d120 = -2.455588849559558;
  d202 = -0.00013138193195963194;
  d201 = 4.750097996232555e-05;
  d200 = 3.082799336172893e-05;
  d212 = 0.0015827056116190592;
  d211 = -0.0006923427069924896;
  d210 = -0.0003105308152917312;
  d222 = -0.00448167876155956;
  d221 = 0.0021740498780173922;
  d220 = 0.0009061371427948547;
  d0  =  0.05098386692906611;
  
  d10 = d102*chi1*chi1 + d101*chi1 + d100;
  d11 = d112*chi1*chi1 + d111*chi1 + d110;
  d12 = d122*chi1*chi1 + d121*chi1 + d120;
  d20 = d202*chi1*chi1 + d201*chi1 + d200;
  d21 = d212*chi1*chi1 + d211*chi1 + d210;
  d22 = d222*chi1*chi1 + d221*chi1 + d220;

  d1 = d10*nu + d11*nu*nu + d12*nu*nu*nu;
  d2 = d20*nu + d21*nu*nu + d22*nu*nu*nu;

  Op[1] = ( ( 1 + (d1*lambda + d2*lambda*lambda) ) / ( (1 + lambda*d0*nu*nu)*(1 + lambda*d0*nu*nu) ) );
  if(Op[1]>1.0){Op[1]=1.0;};

  /* l=3, m=2 */
  d111 = -0.1855265658329381;
  d110 = 0.169287918459171;
  d121 = 0.7992866688623452;
  d120 = -0.7404661978046917;
  d211 = 3.50257349367862e-05;
  d210 = -4.535640993857267e-05;
  d221 = -0.0001494414775185912;
  d220 = 0.000199499873753364;
  d311 = 1.0541708807321233;
  d310 = -1.197499895272247;
  d321 = -4.434903355213771;
  d320 = 4.265315447317003;

  d11 = d111*chi1 + d110;
  d12 = d121*chi1 + d120;
  d21 = d211*chi1 + d210;
  d22 = d221*chi1 + d220;
  d31 = d311*chi1 + d310;
  d32 = d321*chi1 + d320;

  d1 = d11*nu + d12*nu*nu;
  d2 = d21*nu + d22*nu*nu;
  d3 = d31*nu + d32*nu*nu;

  Op[3] = ( 1 + (d1*lambda + d2*lambda*lambda) ) / ( (1 + lambda*d3*d3 ) );
  if(Op[3]>1.0){Op[3]=1.0;};

  /* l=3, m=3 */
  d113 = -0.17571363206227059;
  d112 = -0.06305375338478839;
  d111 = 0.03980474337705321;
  d110 = -0.013248652179259336;
  d123 = 0.8640147802958323;
  d122 = 0.3102264177830439;
  d121 = -0.19421480848909106;
  d120 = -0.0002218569623679201;
  d213 = 0.00025611940982543836;
  d212 = 7.313903139948755e-05;
  d211 = -4.717653263343681e-05;
  d210 = 7.759871362428855e-05;
  d223 = -0.0012519413471327576;
  d222 = -0.0003946916912993047;
  d221 = 0.00022067454550573677;
  d220 = -0.00028520053777508926;
  d31 = 3.225951474170114e-05;
  d32 = -0.00011929010217509204;

  d11 = d113*chi1*chi1*chi1 + d112*chi1*chi1 + d111*chi1 + d110;
  d12 = d123*chi1*chi1*chi1 + d122*chi1*chi1 + d121*chi1 + d120;
  d21 = d213*chi1*chi1*chi1 + d212*chi1*chi1 + d211*chi1 + d210;
  d22 = d223*chi1*chi1*chi1 + d222*chi1*chi1 + d221*chi1 + d220;

  d1 = d11*nu + d12*nu*nu; 
  d2 = d21*nu + d22*nu*nu;
  d3 = d31*nu + d32*nu*nu;

  Op[4] = ( ( 1 + (d1*lambda + d2*lambda*lambda) ) / ( (1 + lambda*lambda*d3)*(1 + lambda*lambda*d3) ) );
  if(Op[4]>1.0){Op[4]=1.0;};

  /* l=4, m=4 */
  d113 = -0.47548874764193955;
  d112 = -0.027086853803733532;
  d111 = 0.7404314951799879;
  d110 = 0.3719543086750542;
  d123 = 2.163479553205327;
  d122 = 0.5275643836925328;
  d121 = -2.7852826426215676;
  d120 = -1.4670455339066886;
  d213 = 7.505291321097175e-05;
  d212 = 9.023386532929272e-05;
  d211 = -7.837977899026634e-05;
  d210 = -3.119139005074053e-05;
  d223 = -0.0005908587238009375;
  d222 = -0.0008832619825909595;
  d221 = 8.362591371834694e-05;
  d220 = 8.203966372546923e-05;
  d311 = 2.2996565572437153;
  d310 = 2.4505796489855025;
  d321 = -6.258030033905745;
  d320 = -7.485076293892947;   

  d11 = d113*chi1*chi1*chi1 + d112*chi1*chi1 + d111*chi1 + d110;
  d12 = d123*chi1*chi1*chi1 + d122*chi1*chi1 + d121*chi1 + d120;
  d21 = d213*chi1*chi1*chi1 + d212*chi1*chi1 + d211*chi1 + d210;
  d22 = d223*chi1*chi1*chi1 + d222*chi1*chi1 + d221*chi1 + d220;
  d31 = d311*chi1 + d310;
  d32 = d321*chi1 + d320;
    
  d1 = d11*nu + d12*nu*nu; 
  d2 = d21*nu + d22*nu*nu;
  d3 = d31*nu + d32*nu*nu;

  Op[8] = ( 1 + (d1*lambda + d2*lambda*lambda  ) ) / ( (1 + lambda*d3*d3 ) );

  // Peak values for all modes
  for (int k=0; k<KMAX; k++) {
      Apeak[k] = Ap[k]*ap_bbh[k];
      Opeak[k] = Op[k]*op_bbh[k];    
  }

  if((X12==0)&&(chi1==0)&&(chi2==0)){ // TODO: problematic cases where X12=0 and messes up the fits
      Apeak[0] = 0.01;
      Apeak[4] = 0.01;
      Apeak[13] = 0.001;
    }

  /* Apply deviations from merger amplitude, frequency */
  if (VERBOSE & (EOBPars->delta_Alm_mrg_size + EOBPars->delta_Omglm_mrg_size > 0))
    printf("Applying deviations from NR amplitude, frequency at merger.\n");
  if (EOBPars->delta_Alm_mrg_size > 0)
    apply_mode_deviations(Apeak, EOBPars->delta_Alm_mrg, EOBPars->delta_Alm_mrg_k, EOBPars->delta_Alm_mrg_size);
  if (EOBPars->delta_Omglm_mrg_size > 0)
    apply_mode_deviations(Opeak, EOBPars->delta_Omglm_mrg, EOBPars->delta_Omglm_mrg_k, EOBPars->delta_Omglm_mrg_size);
  
  if(VERBOSE) PRFORMd("A22_peak",Apeak[1]);
  if(VERBOSE) PRFORMd("omega22_peak",Opeak[1]);
}

/** 
 *  Function: QNM_bhns_td
 *  -----------------------
 *   \details QNMs for tidally disrupted BHNS 
 *   From Gonzalez+2022 and updated for Gonzalez+ in prep.
 *
 *   @param[in] af                               : spin of the remnant BH
 *   @param[in] alpha1                           : real part of the complex QNM frequency for the fundamental mode                                                                    
 *   @param[in] alpha2                           : real part of the complex QNM frequency for the first overtone                                                                    
 *   @param[in] omega1                           : imaginary part of the complex QNM frequency for the fundamental mode
 *   @param[in] omega2                           : imaginary part of the complex QNM frequency for the first overtone
 *   @param[in] alpha21                          : final spin of the remnant BH
 *   @param[out] lambda                          : Quadrupolar tidal parameter of the NS
 *   @param[out] nu                              : Symmetric mass ratio
 *   @param[out] chi1                            : dimensionless spin of the BH
 */
void QNM_bhns_td(double af, double *alpha1, double *alpha2, double *omega1, double *omega2, double *alpha21, double lambda, double nu, double chi1)
{
  double alpha_bbh[KMAX],omega_bbh[KMAX],Ap[KMAX],Op[KMAX];

  /** QNM fits for the tidal disruption cases in BHNS */
  for (int k=0; k<KMAX; k++) {
    alpha21[k] = alpha1[k] = alpha2[k] = alpha_bbh[k] = omega1[k] = omega2[k]= omega_bbh[k] = 0.;
    Ap[k] = Op[k] = 1.;
  }

  // BBH fits
  double af2 = pow(af,2);
  double af3 = pow(af,3);

    /* l = 2, m = 1*/
  omega_bbh[0]  = +0.373672*(1 - 0.79546*af - 0.1908*af2 + 0.11460*af3)/(1 - 0.96337*af - 0.1495*af2 + 0.19522*af3);
  alpha_bbh[0]  = +0.0889623*(1 - 1.31253*af - 0.21033*af2 + 0.52502*af3)/(1 - 1.30041*af - 0.1566*af2 + 0.46204*af3);
  alpha21[0] = +0.184952*(1 - 1.1329*af - 0.3520*af2 + 0.4924*af3)/(1 - 1.10334*af - 0.3037*af2 + 0.4262*af3);
  
    /* l = 2, m = 2*/
  /* double alpha1_c = 0.1211263886*af3 + 0.7015835813*af2 - 1.8226060896*af + 1;
  double alpha1_d = 0.0811633377*af3 + 0.7201166020*af2 - 1.8002031358*af + 1;
  alpha_bbh[1] = 0.0889623157 * (alpha1_c/alpha1_d);
  
  double omega1_c = -0.0598837831*af3 + 0.8082136788*af2 - 1.7408467418*af + 1;
  double omega1_d = -0.2358960279*af3 + 1.3152369374*af2 - 2.0764065380*af + 1;
  omega_bbh[1] =  0.3736716844*(omega1_c/omega1_d); */
  omega_bbh[1]  = +0.373672*(1 - 1.5367*af + 0.5503*af2)/(1 - 1.8700*af + 0.9848*af2 - 0.10943*af3);
  alpha_bbh[1]  = +0.08896*(1 - 1.90036*af + 0.86200*af2 + 0.0384893*af3)/(1 - 1.87933*af + 0.88062*af2);
  alpha21[1] = +0.184953*(1 - 1.89397*af + 0.88126*af2 + 0.0130256*af3)/(1 - 1.83901*af + 0.84162*af2);
  
  /* l = 3, m = 2*/
  omega_bbh[3]  = +0.599443*(1 - 0.251*af - 0.891*af2 + 0.2706*af3)/(1 - 0.475*af - 0.911*af2 + 0.4609*af3);
  alpha_bbh[3]  = +0.0927030*(1 - 1.58277*af + 0.2783*af2 + 0.30503*af3)/(1 - 1.56797*af + 0.3290*af2 + 0.24155*af3);
  alpha21[3] = +0.188595*(1 - 1.5212*af + 0.1563*af2 + 0.3652*af3)/(1 - 1.4968*af + 0.1968*af2 + 0.3021*af3);
  
  /* l = 3, m = 3*/
  omega_bbh[4]  = +0.599443*(1 - 1.84922*af + 0.9294*af2 - 0.07613*af3)/(1 - 2.18719*af + 1.4903*af2 - 0.3014*af3);
  alpha_bbh[4]  = +0.0927030*(1 - 1.8310*af + 0.7568*af2 + 0.0745*af3)/(1 - 1.8098*af + 0.7926*af2 + 0.0196*af3);
  alpha21[4] = +0.188595*(1 - 1.8011*af + 0.7046*af2 + 0.0968*af3)/(1 - 1.7653*af + 0.7176*af2 + 0.0504*af3);
  
  /* l = 4, m = 4*/
  omega_bbh[8]  = +0.809178*(1 - 1.83156*af + 0.9016*af2 - 0.06579*af3)/(1 - 2.17745*af + 1.4753*af2 - 0.2961*af3);
  alpha_bbh[8]  = +0.0941640*(1 - 1.8662*af + 0.8248*af2 + 0.0417*af3)/(1 - 1.8514*af + 0.8736*af2 - 0.0198*af3);
  alpha21[8] = +0.190170*(1 - 1.8546*af + 0.8041*af2 + 0.0507*af3)/(1 - 1.8315*af + 0.8391*af2 - 0.0051*af3);

  /* l = 5, m = 5*/
  omega_bbh[13]  = +1.012295*(1 - 1.5659*af + 0.5783*af2)/(1 - 1.9149*af + 1.0668*af2 - 0.14663*af3);
  alpha_bbh[13]  = +0.0948705*(1 - 1.8845*af + 0.8585*af2 + 0.0263*af3)/(1 - 1.8740*af + 0.9147*af2 - 0.0384*af3);
  alpha21[13] = +0.190947*(1 - 1.8780*af + 0.8467*af2 + 0.0315*af3)/(1 - 1.8619*af + 0.8936*af2 - 0.0293*af3);
  //


  // Inverse damping time
  double c101,c100,c111,c110,c201,c200,c211,c210,c301,c300,c311,c310,c401,c400;
  double c113,c112,c123,c122,c121,c120,c213,c212,c223,c222,c221,c220;
  double c132,c131,c130,c231,c232,c230,c322,c321,c320,c332,c331,c330,c411,c410,c421,c420;
  double c10,c11,c20,c21,c30,c31,c40,c12,c22,c13,c23,c32,c33,c41,c42;
  double c1,c2,c3,c4;

  /* l = 2, m = 1*/
  c113 = 1.017251288502842;
  c112 = -0.06623359086717126;
  c111 = -0.18161663646183168;
  c110 = -0.04062047724559591;
  c123 = -4.795790111971402;
  c122 = 0.6400132907804421;
  c121 = 1.154845779515909;
  c120 = 0.2470137810794439;
  c213 = -0.0006508038132935816;
  c212 = 0.00020189004480464088;
  c211 = 4.007524560345456e-05;
  c210 = -1.4580017533437743e-05;
  c223 = 0.0035685097544596012;
  c222 = -0.0003690285105227662;
  c221 = 5.645560621097093e-05;
  c220 = 0.00011406588975419034;
  c311 = 0.832457585187816;
  c310 = 0.07742292789237941;
  c321 = -2.3647908436141116;
  c320 = 0.4852248611895286;
   
  c11 = c113*chi1*chi1*chi1 + c112*chi1*chi1 + c111*chi1 + c110;
  c12 = c123*chi1*chi1*chi1 + c122*chi1*chi1 + c121*chi1 + c120;
  c21 = c213*chi1*chi1*chi1 + c212*chi1*chi1 + c211*chi1 + c210;
  c22 = c223*chi1*chi1*chi1 + c222*chi1*chi1 + c221*chi1 + c220;
  c31 = c311*chi1 + c310;
  c32 = c321*chi1 + c320;

  c1 = c11*nu + c12*nu*nu; 
  c2 = c21*nu + c22*nu*nu;
  c3 = c31*nu + c32*nu*nu;

  Ap[0] = ( 1 + (c1*lambda + c2*lambda*lambda  ) ) / ( (1 + lambda*c3*c3 )*(1 + lambda*c3*c3 ) );
  //if(Ap[0]>1.0){Ap[0]=1.0;}
  if(Ap[0]<0.0){Ap[0]=1.0;}

  /* l = 2, m = 2*/
  c101 = 0.070454839545944;
  c100 = 0.006934183212582287;
  c111 = -0.28991238536521063;
  c110 = 0.07310357931603813;
  c201 = -5.195489438172421e-05;
  c200 = 4.9944435533574586e-05;
  c211 = 0.00017505988793666536;
  c210 = -0.0002454447042581111;
  c301 = 6.38200707246774e-09;
  c300 = -1.1341504662808922e-08;
  c311 = -1.1554624161065344e-08;
  c310 = 5.2615121720619686e-08;
  c401 = -0.0025816920977441917;
  c400 = 0.05188115116342483;

  c10 = c101*chi1 + c100;
  c11 = c111*chi1 + c110;
  c20 = c201*chi1 + c200;
  c21 = c211*chi1 + c210;
  c30 = c301*chi1 + c300;
  c31 = c311*chi1 + c310;
  c40 = c401*chi1 + c400;

  c1 = c10*nu + c11*nu*nu;
  c2 = c20*nu + c21*nu*nu;
  c3 = c30*nu + c31*nu*nu;

  Ap[1] = (1 + lambda*c1 + lambda*lambda*c2 + lambda*lambda*lambda*c3 ) / ( 1 + lambda*c40*nu );
  if(Ap[1]>1.0){Ap[1]=1.0;}

  /* l = 3, m = 2*/
  c122 = 13.76187812933075;
  c121 = 11.217068236820518;
  c120 = -3.9975070427233783;
  c132 = -61.37639466623299;
  c131 = -50.96198257968971;
  c130 = 18.43384336942182;
  c222 = -0.06616749911912588;
  c221 = -0.051956237187389;
  c220 = 0.020803793852828833;
  c232 = 0.29667798696415937;
  c231 = 0.23894376648286575;
  c230 = -0.09583888462419352;
  c322 = 5.450459709653773e-05;
  c321 = 4.248357623382005e-05;
  c320 = -1.6118209299630243e-05;
  c332 = -0.0002452746197305893;
  c331 = -0.00019597450119635202;
  c330 = 7.435022952695207e-05;
  c411 = 0.24054582587634482;
  c410 = 0.7681438117400963;
  c421 = -0.38845884374894185;
  c420 = -2.9576614516373314;

  c12 = c122*chi1*chi1 + c121*chi1 + c120;
  c13 = c132*chi1*chi1 + c131*chi1 + c130;
  c22 = c222*chi1*chi1 + c221*chi1 + c220;
  c23 = c232*chi1*chi1 + c231*chi1 + c230;
  c32 = c322*chi1*chi1 + c321*chi1 + c320;
  c33 = c332*chi1*chi1 + c331*chi1 + c330;
  c41 = c411*chi1 + c410;
  c42 = c421*chi1 + c420;

  c1 = c12*nu*nu + c13*nu*nu*nu;
  c2 = c22*nu*nu + c23*nu*nu*nu;
  c3 = c32*nu*nu + c33*nu*nu*nu;
  c4 = c41*nu + c42*nu*nu;

  Ap[3] = ( 1 + (c1*lambda + c2*lambda*lambda + c3*lambda*lambda*lambda) ) / ( (1 + lambda*c4*c4 )*(1 + lambda*c4*c4 ) );
  if(Ap[3]>1.0){Ap[3]=1.0;}
  if(Ap[3]<0.0){Ap[3]=1.0;}

  /* l = 3, m = 3*/
  c113 = 1.5710948456754454;
  c112 = 0.9146276153180174;
  c111 = 0.11132698735711116;
  c110 = -0.033742317650047565;
  c123 = -7.69280560477062;
  c122 = -1.676964030989382;
  c121 = 1.0988217833664355;
  c120 = 0.36291501890859346;
  c213 = -0.0016040433566614856;
  c212 = 9.097724819433372e-05;
  c211 = 0.0003795563942833713;
  c210 = 9.331364459422936e-05;
  c223 = 0.007077732496151121;
  c222 = -0.0014471026733462629;
  c221 = -0.002095780509407605;
  c220 = -0.0004263783108491086;
  c311 = 2.0781424448527246;
  c310 = 0.6641727408438288;
    
  c11 = c113*chi1*chi1*chi1 + c112*chi1*chi1 + c111*chi1 + c110;
  c12 = c123*chi1*chi1*chi1 + c122*chi1*chi1 + c121*chi1 + c120;
  c21 = c213*chi1*chi1*chi1 + c212*chi1*chi1 + c211*chi1 + c210;
  c22 = c223*chi1*chi1*chi1 + c222*chi1*chi1 + c221*chi1 + c220;
  c31 = c311*chi1 + c310;

  c1 = c11*nu + c12*nu*nu; 
  c2 = c21*nu + c22*nu*nu;

  Ap[4] = ( ( 1 + (c1*lambda + c2*lambda*lambda) ) / ( (1 + lambda*c31*c31*nu*nu) ) );
  if(Ap[4]>1.0){Ap[4]=1.0;}

  /* l = 4, m = 4*/
  c112 = 11.242961865950528;
  c111 = -17.389634170001713;
  c110 = 2.7194273570640384;
  c122 = -98.4553492936875;
  c121 = 177.19502092590486;
  c120 = -26.020314732655766;
  c132 = 215.09730936215993;
  c131 = -448.25588490848196;
  c130 = 61.99122487260588;
  c212 = 0.08766877382876546;
  c211 = -0.01207862504552977;
  c210 = 0.01896684410949216;
  c222 = -0.8929260852043949;
  c221 = 0.07249439442967934;
  c220 = -0.1715098874768106;
  c232 = 2.2569730153627905;
  c231 = -0.0749977671495978;
  c230 = 0.38758152889538017;
  c311 = -3.439564725617774;
  c310 = 3.694711197855806;
  c321 = 17.522114041121114;
  c320 = -17.079544591334244;

  c11 = c112*chi1*chi1 + c111*chi1 + c110;
  c12 = c122*chi1*chi1 + c121*chi1 + c120;
  c13 = c132*chi1*chi1 + c131*chi1 + c130;
  c21 = c212*chi1*chi1 + c211*chi1 + c210;
  c22 = c222*chi1*chi1 + c221*chi1 + c220;
  c23 = c232*chi1*chi1 + c231*chi1 + c230;
  c31 = c311*chi1 + c310;
  c32 = c321*chi1 + c320;

  c1 = c11*nu + c12*nu*nu + c13*nu*nu*nu;
  c2 = c21*nu + c22*nu*nu + c23*nu*nu*nu;
  c3 = c31*nu + c32*nu*nu;

  Ap[8] = ( 1 + (c1*lambda + c2*lambda*lambda  ) ) / ( (1 + lambda*c3*c3 )*(1 + lambda*c3*c3 ) );
  if(Ap[8]>1.0){Ap[8]=1.0;}
  if(Ap[8]<0.0){
    Ap[8]=fabs(Ap[8]);
    if(Ap[8]>1.0){Ap[8]=1.0;}
    }

  // Frequency
  double d101,d100,d111,d110,d201,d200,d211,d210,d301,d300,d311,d310,d401,d400;
  double d113,d112,d123,d122,d121,d120,d213,d212,d223,d222,d221,d220,d321,d320,d411,d410,d131,d130,d231,d230;
  double d10,d11,d20,d21,d30,d31,d40,d12,d22,d32,d41,d13,d23;
  double d1,d2,d3,d4;

  /* l = 2, m = 1*/
  d113 = 3.333420017262863;
  d112 = 1.6733100922615411;
  d111 = 0.6428123126575698;
  d110 = 0.15345063221041327;
  d123 = -16.73476556088638;
  d122 = -9.6165271440306;
  d121 = -3.9666698955903628;
  d120 = -0.9188466914621823;
  d213 = -0.004333646061173743;
  d212 = -0.006690849241886476;
  d211 = -0.004053501377389173;
  d210 = -0.0006315097800444771;
  d223 = 0.024848664430444187;
  d222 = 0.03887710063247387;
  d221 = 0.023247889326064287;
  d220 = 0.0037972065617504825;
  d311 = -0.036893725055918904;
  d310 = -0.007460440219362085;
  d321 = 0.25839558187126116;
  d320 = 0.09182425362740411;

  d11 = d113*chi1*chi1*chi1 + d112*chi1*chi1 + d111*chi1 + d110;
  d12 = d123*chi1*chi1*chi1 + d122*chi1*chi1 + d121*chi1 + d120;
  d21 = d213*chi1*chi1*chi1 + d212*chi1*chi1 + d211*chi1 + d210;
  d22 = d223*chi1*chi1*chi1 + d222*chi1*chi1 + d221*chi1 + d220;
  d31 = d311*chi1 + d310;
  d32 = d321*chi1 + d320;

  d1 = d11*nu + d12*nu*nu; 
  d2 = d21*nu + d22*nu*nu;
  d3 = d31*nu + d32*nu*nu; 

  Op[0] =  ( 1 + (d1*lambda + d2*lambda*lambda  ) ) / ( (1 + lambda*lambda*d3*d3 )*(1 + lambda*lambda*d3*d3 ) );
  if(Op[0]>1.0){Op[0]=1.0;}
  if(Op[0]<0.0){Op[0]=1.0;}

  /* l = 2, m = 2*/
  d101 = -1.87204534166561;
  d100 = 2.2235911158585906;	
  d111 = 6.587073531924587;
  d110 = -7.455030280325616;
  d201 = -8.856884879548845e-05;
  d200 = -0.00019454696765165838;
  d211 = 0.0004981829301131264;
  d210 = 0.002728475621165506;
  d301 = 0.08686135769230043;
  d300 = 0.7715191055994437;
  d311 = -0.2756808384095736;
  d310 = -0.4740847695179648;

  d10 = d101*chi1 + d100;
  d11 = d111*chi1 + d110;
  d20 = d201*chi1 + d200;
  d21 = d211*chi1 + d210;
  d30 = d301*chi1 + d300;
  d31 = d311*chi1 + d310;

  d1 = d10*nu + d11*nu*nu;
  d2 = d20*nu + d21*nu*nu;
  d3 = d30*nu + d31*nu*nu;

  Op[1] = ( 1 + lambda*d1 + lambda*lambda*d2 ) / (( 1 + lambda*d3*d3 )*( 1 + lambda*d3*d3 ));
  if(Op[1]>1.0){Op[1]=1.0;}
  if(Op[1]<0.0){Op[1]=0.0;}

  /* l = 3, m = 2 */
  d111 = 0.026527628808116717;
  d110 = -0.02926417906773356;
  d121 = -0.15885016325333506;
  d120 = 0.10134101940210988;
  d211 = 4.410583198188113e-06;
  d210 = 5.7505152559658577e-05;
  d221 = 2.665524286657608e-05;
  d220 = -0.00023554822250791658;
  d311 = -2.8273536377662663e-08;
  d310 = -2.9510126501533417e-08;
  d321 = 1.2052971451081137e-07;
  d320 = 1.302175520779655e-07;
  d411 = 0.5314015395157192;
  d410 = 0.6120130890466131;

  d11 = d111*chi1 + d110;
  d12 = d121*chi1 + d120;
  d21 = d211*chi1 + d210;
  d22 = d221*chi1 + d220;
  d31 = d311*chi1 + d310;
  d32 = d321*chi1 + d320;
  d41 = d411*chi1 + d410;

  d1 = d11*nu + d12*nu*nu;
  d2 = d21*nu + d22*nu*nu;
  d3 = d31*nu + d32*nu*nu;
  d4 = d41*nu*nu;

  Op[3] = ( 1 + (d1*lambda + d2*lambda*lambda + d3*lambda*lambda*lambda) ) / ( (1 + lambda*d4*d4 ) );
  if(Op[3]>1.0){Op[3]=1.0;}
  if(Op[3]<0.0){Op[3]=0.0;}

  /* l = 3, m = 3 */
  d113 = -0.04221771673205502;
  d112 = -0.02242685682240473;
  d111 = 0.028995664124488466;
  d110 = 0.003520493408631722;
  d123 = -0.7437282727260565;
  d122 = -0.20713956840265466;
  d121 = 0.5106082146377563;
  d120 = 0.25791618347153406;
  d213 = 1.1142401781174221e-05;
  d212 = 4.331182069595073e-06;
  d211 = -9.533567551600753e-06;
  d210 = -1.4195826636478271e-06;
  d223 = -0.00020041720650518658;
  d222 = 0.001972740713879979;
  d221 = 0.0025947326664176266;
  d220 = 0.0007966369219430112;
  d311 = 0.41309805786920845;
  d310 = 0.24603136128672057;
    
  d11 = d113*chi1*chi1*chi1 + d112*chi1*chi1 + d111*chi1 + d110;
  d12 = d123*chi1*chi1*chi1 + d122*chi1*chi1 + d121*chi1 + d120;
  d21 = d213*chi1*chi1*chi1 + d212*chi1*chi1 + d211*chi1 + d210;
  d22 = d223*chi1*chi1*chi1 + d222*chi1*chi1 + d221*chi1 + d220;
  d31 = d311*chi1 + d310;

  d1 = d11*nu + d12*nu*nu; 
  d2 = d21*nu + d22*nu*nu;

  Op[4] = ( ( 1 + (d1*lambda + d2*lambda*lambda) ) / ( (1 + lambda*lambda*(d31*nu*nu)*(d31*nu*nu)) ) );
  if(Op[4]>1.0){Op[4]=1.0;}

  /* l = 4, m = 4 */
  d111 = -0.8256636431354363;
  d110 = -0.369302975353016;
  d121 = 8.411035539159677;
  d120 = 3.826674376367923;
  d131 = -21.3219318080355;
  d130 = -9.781112688287969;
  d211 = 0.0008933783780522222;
  d210 = 0.0004863231203924392;
  d221 = -0.009256763786840302;
  d220 = -0.004995621509254127;
  d231 = 0.023730168002922283;
  d230 = 0.012722108040021144;
  d321 = -0.06501172490507952;
  d320 = 0.6538422755200813;

  d11 = d111*chi1 + d110;
  d12 = d121*chi1 + d120;
  d13 = d131*chi1 + d130;
  d21 = d211*chi1 + d210;
  d22 = d221*chi1 + d220;
  d23 = d231*chi1 + d230;
  d32 = d321*chi1 + d320;

  d1 = d11*nu + d12*nu*nu + d13*nu*nu*nu;
  d2 = d21*nu + d22*nu*nu + d23*nu*nu*nu;
  d3 = d32*nu*nu;

  Op[8] = ( 1 + (d1*lambda + d2*lambda*lambda  ) ) / ( (1 + lambda*d3*d3 )*(1 + lambda*d3*d3 ) );
  if(Op[8]>1.0){Op[8]=1.0;}
  if(Op[8]<0.0){Op[8]=0.0;}

  for (int k=0; k<KMAX; k++) {
      alpha2[k] = Ap[k]*alpha_bbh[k];
      omega1[k] = Op[k]*omega_bbh[k]; 
      alpha21[k] = alpha21[k]*Ap[k];
    
    
    if(k==1){alpha21[k] = alpha2[k] - alpha1[k];}

    alpha1[k] = alpha2[k];
  } 

  /* Apply user-input deviations from fitted values of alpha, omega */
  if (EOBPars->delta_alphalm0_size > 0 || EOBPars->delta_omglm0_size > 0) {
    if (VERBOSE) printf("Applying user-input deviations from QNM frequencies.\n");
    QNM_deviations(alpha1, omega1, alpha21);
  }
}

/** 
 *  Function: eob_nqc_point_BHNS_HM
 *  -------------------------------
 *   \details Fits for the NR point used to determine NQC corrections for BHNSs
 *   From Gonzalez et al 2022 and updated for Gonzalez+ in prep.
 *   
 *   @param[in]   dyn      : EOB dynamics
 *   @param[out]  A_tmp    : RWZ-normalized Amplitude at the NQC point 
 *   @param[out]  dA_tmp   : first time derivative of the RWZ-normalized amplitude at the NQC point 
 *   @param[out]  omg_tmp  : omega at the NQC point 
 *   @param[out]  domg_tmp : first time derivative of the frequency at the NQC point 
 *   @param[in]   abh      : spin of the final BH
 *
 */
void eob_nqc_point_BHNS_HM(Dynamics *dyn, double *A_tmp, double *dA_tmp, double *omg_tmp, double *domg_tmp)
{
  /** BHNS fits just for omega^NQC and domega^NQC, and amplitude. With BBH from arxiv:2001.09082 **/
  const double nu   = EOBPars->nu;
  const double X1   = EOBPars->X1;
  const double X2   = EOBPars->X2;
  const double chi1 = EOBPars->chi1;
  const double chi2 = EOBPars->chi2;
  const double aK   = EOBPars->a1 + EOBPars->a2;
  const double lambda = EOBPars->LambdaBl2;

  const double nu2  = SQ(nu);
  const double nu3  = nu2*nu;
  const double X12  = X1 - X2;  
  const double aK2  = SQ(aK);
  const double aK3  = aK2*aK;
  const double aK4  = aK2*aK2;
  double a12  = X1*chi1 - X2*chi2;
  const double Shat = 0.5*(aK + a12*X12);
  const double Shat2 = Shat*Shat;

  double pA[5], pdA[5], pomg[5], pdomg[5];
  double c_p1,     c_p2,     c_p3,   c_p4;
  double c_pdA1,   c_pdA2,   c_pdA3, c_pdA4;
  double c_pdomg1, c_pdomg2;
  double n0, d1;

  double a0_omg_tmp, a1_omg_tmp, a2_omg_tmp, b0_omg_tmp, b1_omg_tmp, b2_omg_tmp, b3_omg_tmp, c11_omg_tmp, c12_omg_tmp, c21_omg_tmp, c22_omg_tmp, c31_omg_tmp, c32_omg_tmp;
  double a0_domg_tmp, a1_domg_tmp, a2_domg_tmp, b1_domg_tmp, b2_domg_tmp, b3_domg_tmp, b4_domg_tmp, c1_domg_tmp, c2_domg_tmp;
  double scale, a0_A_tmp, a1_A_tmp , a2_A_tmp, b0_A_tmp, b1_A_tmp, b2_A_tmp, b3_A_tmp, b4_A_tmp, b5_A_tmp, b6_A_tmp, c11_A_tmp, c12_A_tmp, c21_A_tmp, c22_A_tmp, c31_A_tmp, c32_A_tmp;
  double a0_dA_tmp, a1_dA_tmp, a2_dA_tmp, b1_dA_tmp, b2_dA_tmp, b3_dA_tmp, c1_dA_tmp, c2_dA_tmp;
  double omg_tmp_nu, omg_tmp_equal, omg_tmp_den, domg_tmp_nu, domg_tmp_equal, A_tmp_scale_nu, A_tmp_scale_equal, A_tmp_scale_den, dA_tmp_scale_nu, dA_tmp_scale_equal;

  double p1[2], p2[2], p3[2], p4[2]; 
  double pn0[2], pd1[2], ppdomg1[2], ppdomg2[2], pdA1[2],pdA2[2],pdA3[2],pdA4[2];


  /* l=2, m=2 */	
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
    A_tmp[1]   =  c_p1*aK3 + c_p2*aK2 + c_p3*aK + c_p4;
	      
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
    dA_tmp[1]  =  c_pdA1*aK3   + c_pdA2*aK2 + c_pdA3*aK+ c_pdA4;
	    
    pn0[0]     =  0.46908067;
    pn0[1]     =  0.27022141;
    pd1[0]     =  0.64131115;
    pd1[1]     = -0.37878384;
    n0         =  pn0[0]*nu + pn0[1];
    d1         =  pd1[0]*nu + pd1[1];
    omg_tmp[1] =  n0/(1 + d1*aK);
	    
    ppdomg1[0] =  0.00061175;
    ppdomg1[1] =  0.00074001;
    ppdomg2[0] =  0.02504442;
    ppdomg2[1] =  0.00548217;
    c_pdomg1   =  ppdomg1[0]*nu + ppdomg1[1];
    c_pdomg2   =  ppdomg2[0]*nu + ppdomg2[1];
    domg_tmp[1] =  c_pdomg1*aK   + c_pdomg2;

  /* l=2, m=1 */
  /* BBH Pompili fits */
  const double delta = sqrt(1-4.*nu);
  const double chi_S = (chi1+chi2)/2;
  const double chi_A = (chi1-chi2)/2;
  const double chi   = chi_S + chi_A*delta/(1.-2.*nu);
  const double chisq = chi*chi;
  const double chi3  = chisq*chi;
  const double chi4  = chi3*chi;
  const double chi_21A  = delta*chi_S/(1-1.3*nu) + chi_A;
  const double chi_21A2 = SQ(chi_21A);
  const double chi_21A3 = chi_21A2*chi_21A;
  const double chi_21D  = delta*chi_S/(1-2.*nu)  + chi_A;
  A_tmp[0]    = fabs(-0.033175*chi_21A3*delta + 0.086356*chi_21A2*delta*nu - 0.049897*chi_21A2*delta + 0.012706*chi_21A*delta + 0.168668*chi_21A*nu - 0.285597*chi_21A + 1.067921*delta*nu2 - 0.189346*delta*nu + 0.431426*delta);
  dA_tmp[0]   = chi_21D*delta*(0.023534*nu - 0.008064) + delta*(0.006743 - 0.0297*nu) + 0.008256*fabs( (chi_21D - delta*(5.471011*nu2 + 1.235589*nu + 0.815482)) );
  omg_tmp[0]  = 0.01009*chi3 - 0.077343*chisq*nu + 0.02411*chisq + 0.168854*chi*nu2 - 0.159382*chi*nu + 0.047635*chi + 1.965157*nu3 - 0.53085*nu2 + 0.237904*nu + 0.176526;
  domg_tmp[0] = 0.00149*chi3 + 0.008965*chisq*nu - 0.002739*chisq - 0.033831*chi*nu2 + 0.005752*chi*nu - 0.002003*chi + 0.204368*nu3 - 0.120705*nu2 + 0.035144*nu + 0.006579;
  
  int l       = LINDEX[0];
  double fact = sqrt((l+2)*(l+1)*l*(l-1));
  A_tmp[0]    = A_tmp[0]/fact;
  dA_tmp[0]   = dA_tmp[0]/fact;

  /* l=3, m=1 */
  A_tmp[2]    = 0.00520201*X12*(1 - 4.9441*nu + 8.9339*nu2);
  dA_tmp[2]   = -0.00043382*(1 - 9.0479*nu + 23.054*nu2)/(1 + 88.626*nu2);
  omg_tmp[2]  = 0.485186*(1 - 0.4421*nu - 6.8184*nu2);
  domg_tmp[2] = 0.0673118*(1 + 13.318*nu)/(1 + 70.552*nu);
	  
  /* l=3, m=3 */
  /* BBH Pompili fits */
  const double chi_33   = chi_S*delta+chi_A;
  A_tmp[4]    = fabs(- 0.088371*SQ(chi_33)*delta*nu + 0.036258*SQ(chi_33)*delta + 1.057731*chi_33*nu2 - 0.466709*chi_33*nu + 0.099543*chi_33 + 1.96267*delta*nu2 + 0.027833*delta*nu + 0.558808*delta);
  dA_tmp[4]   = SQ(chi_33)*delta*(0.004941*nu - 0.002094) + 0.001781*sqrt(fabs( (SQ(chi_33) + chi_33*delta*(39.247538*nu - 2.986889) + SQ(delta)*(85.173306*nu + 4.637906)) ));
  omg_tmp[4]  = 0.045141*chi3 - 0.346675*chisq*nu + 0.119419*chisq + 0.745924*chi*nu2 - 0.478915*chi*nu + 0.17467*chi - 8.887163*nu3 + 4.226831*nu2 + 0.427167;
  domg_tmp[4] = 0.001697*chi3 + 0.016231*chisq*nu - 0.003985*chisq - 0.154378*chi*nu2 + 0.050618*chi*nu - 0.002721*chi - 0.255402*nu3 + 0.08663*nu2 + 0.027405*nu + 0.009736;

  l       = LINDEX[4];
  fact = sqrt((l+2)*(l+1)*l*(l-1));
  A_tmp[4]    = A_tmp[4]/fact;
  dA_tmp[4]   = dA_tmp[4]/fact;
	  
  /* l=4, m=1 */
  A_tmp[5]    = 0.00043987*X12*(1-2.*nu)*(1 - 8.4975*nu + 27.31*nu2)/(1 - 1.2002*nu);
  dA_tmp[5]   = -0.00001827*(1 - 2.8242*nu - 3.1871*nu2);
  omg_tmp[5]  = 0.735051*(1 - 8.3628*nu + 20.529*nu2)/(1 - 7.4883*nu + 18.695*nu2);
  domg_tmp[5] = 0.0739078*(1 + 0.99186*nu - 19.435*nu2);
	
  /* l=4, m=4 */
  /* BBH Pompili fits */
  const double chi_44A  = (1 - 5.*nu)*chi_S + chi_A*delta;
  const double chi_44D  = (1 - 7.*nu)*chi_S + chi_A*delta;
  A_tmp[8]    = fabs( 0.031483*SQ(chi_44A) - 0.180165*chi_44A*nu + 0.063931*chi_44A + 6.239418*nu3 - 1.947473*nu2 - 0.615307*nu + 0.262533 );
  dA_tmp[8]   = -0.001251*SQ(chi_44D)*chi_44D + 0.006387*SQ(chi_44D)*nu - 0.001223*SQ(chi_44D) - 0.034308*chi_44D*nu2 + 0.014373*chi_44D*nu - 0.000681*chi_44D + 1.134679*nu3 - 0.417056*nu2 + 0.024004*nu + 0.003498;
  omg_tmp[8]  = 0.042529*chi3 - 0.415864*chisq*nu + 0.155222*chisq + 0.768712*chi*nu2 - 0.592568*chi*nu + 0.244508*chi - 13.651335*nu3 + 5.490329*nu2 + 0.574041;
  domg_tmp[8] = 0.001812*chi3 + 0.024687*chisq*nu - 0.00568*chisq  - 0.162693*chi*nu2 + 0.061205*chi*nu - 0.003623*chi - 0.536664*nu3  + 0.094797*nu2 + 0.045406*nu + 0.013038;
	
  l       = LINDEX[8];
  fact = sqrt((l+2)*(l+1)*l*(l-1));
  A_tmp[8]    = A_tmp[8]/fact;
  dA_tmp[8]   = dA_tmp[8]/fact;

  /* l=5, m=5 */
  b1_A_tmp = -0.0829447;
  b2_A_tmp = +8.434479;
  b3_A_tmp = +1.197758;
  b4_A_tmp = -0.260704;
	
  A_tmp[13] = 0.00516272*X12*(1-2*nu)*(1 - 0.68569*nu + 7.7267*nu2) + 0.034503277967117711167*(0.04402930 + b1_A_tmp*nu)/(1 + b2_A_tmp*nu)*a12/(1 + (-0.5749409+b3_A_tmp*nu)/(1 + b4_A_tmp*nu)*a12);
	
  b1_dA_tmp = +0.00171106;
  b2_dA_tmp = +0.00475246;
  b3_dA_tmp = +6.508129;
	  
  dA_tmp[13] = -0.0000658001*(1 - 3.3212*nu) + 0.034503277967117711167*((0.0002250231 + b1_dA_tmp*nu)*a12 + (0.0002322527 + b2_dA_tmp*nu)*SQ(a12))/(1 + (-0.8223640 + b3_dA_tmp*nu)*a12);
	
  b1_omg_tmp = -1.711526;
  b2_omg_tmp = -1.546737;
  b3_omg_tmp = -2.023890;
	
  omg_tmp[13] =  0.855016*(1 + 1.1639*nu)*(1 + (-0.3577869 + b1_omg_tmp*nu)*Shat + (-0.0627735 + b2_omg_tmp*nu)*Shat2)/(1 + (-0.7053108 + b3_omg_tmp*nu)*Shat);
	
  b1_domg_tmp = -0.01040098;
  b2_domg_tmp = +0.1384369;
  b3_domg_tmp = -0.03244646;
  b4_domg_tmp = +0.2250689;
	  
  domg_tmp[13] =  0.0178326*(1 + 3.1304*nu) + (b1_domg_tmp + b2_domg_tmp*nu)*Shat + (b3_domg_tmp + b4_domg_tmp*nu)*Shat2;

  /** BHNS fits **/
  double Ap[KMAX],Op[KMAX],dAp[KMAX],dOp[KMAX];

  for (int k=0; k<KMAX; k++) {
    Ap[k] = Op[k] = dAp[k] = dOp[k] = 1.;
  }

  /* Omega */
  double f102,f101,f100,f112,f111,f110,f122,f121,f120,f202,f201,f200,f212,f211,f210,f222,f221,f220;
  double f113,f123,f213,f223,f133,f132,f131,f130,f233,f231,f232,f230,f313,f312,f311,f310,f323,f322,f321;
  double f320,f333,f332,f331,f330,f411,f410;
  double f10,f11,f12,f20,f21,f22,f13,f23,f31,f32,f33,f41;
  double f0,f1,f2,f3,f4;

  /* l = 2, m = 1 */
  f113 = 353.17736170387593;
  f112 = 368.77721952946587;
  f111 = 39.99806622388038;
  f110 = 2.520352205338529;
  f123 = -3848.97211049468;
  f122 = -3889.276983501113;
  f121 = -373.68496396133315;
  f120 = -22.818001440139962;
  f133 = 10345.905013323933;
  f132 = 10144.999043641952;
  f131 = 858.125344517891;
  f130 = 50.73245981496382;
  f213 = -0.42947400824171;
  f212 = -1.5050812191894574;
  f211 = -0.4600316680456819;
  f210 = -0.01467825143630204;
  f223 = 5.06471627790655;
  f222 = 15.88763644495495;
  f221 = 4.598576386393477;
  f220 = 0.13865459323013118;
  f233 = -14.50183880723711;
  f232 = -41.37414593942178;
  f231 = -11.449763631646169;
  f230 = -0.3140104217877563;
  f313 = -0.0001845625359841723;
  f312 = 0.0014207429350569227;
  f311 = 0.0006607364734238603;
  f310 = 1.3681280159083329e-05;
  f323 = 0.0013536760570583132;
  f322 = -0.014939518176558928;
  f321 = -0.006692583959047551;
  f320 = -0.00012986782954713027;
  f333 = -0.0022679726827504647;
  f332 = 0.038853528883517646;
  f331 = 0.0168387350523842;
  f330 = 0.00029842358066044174;
  f411 = -0.7140362944053328;
  f410 = 0.38012978267558917;   
    
  f11 = f113*chi1*chi1*chi1 + f112*chi1*chi1 + f111*chi1 + f110;
  f12 = f123*chi1*chi1*chi1 + f122*chi1*chi1 + f121*chi1 + f120;
  f13 = f133*chi1*chi1*chi1 + f132*chi1*chi1 + f131*chi1 + f130;
  f21 = f213*chi1*chi1*chi1 + f212*chi1*chi1 + f211*chi1 + f210;
  f22 = f223*chi1*chi1*chi1 + f222*chi1*chi1 + f221*chi1 + f220;
  f23 = f233*chi1*chi1*chi1 + f232*chi1*chi1 + f231*chi1 + f230;
  f31 = f313*chi1*chi1*chi1 + f312*chi1*chi1 + f311*chi1 + f310;
  f32 = f323*chi1*chi1*chi1 + f322*chi1*chi1 + f321*chi1 + f320;
  f33 = f333*chi1*chi1*chi1 + f332*chi1*chi1 + f331*chi1 + f330;
  f41 = f411*chi1 + f410;    

  f1 = f11*nu + f12*nu*nu + f13*nu*nu*nu;
  f2 = f21*nu + f22*nu*nu + f23*nu*nu*nu;
  f3 = f31*nu + f32*nu*nu + f33*nu*nu*nu;
  f4 = f41*nu; 

  Op[0] = ( 1 + (f1*lambda + f2*lambda*lambda + f3*lambda*lambda*lambda ) ) / ( (1 + lambda*f4*f4 )*(1 + lambda*f4*f4 ) );
  if (Op[0]<0) {Op[0]=1.;}

  /* l = 2, m = 2 */
  f102 = 0.0027730975740547707; 
  f101 = -0.0013160726188961565;
  f100 = 0.04353588700832432;
  f112 = -0.8370275073634688;
  f111 = 0.6536801867953376;
  f110 = 0.3361764212985733;
  f122 = 3.7510346000111587;
  f121 = -3.519057364463967;
  f120 = -0.5788117577224763;
  f202 = -7.458373443450189e-07;
  f201 = 7.853120442799991e-07;
  f200 = 6.100942695272458e-07;
  f212 = 0.00021701656160855214;
  f211 = -0.0012948176199725767;
  f210 = 0.001076170139911975;
  f222 = -0.0006750422874570453;
  f221 = 0.00556807427720544;
  f220 = -0.004130159066577092;
  f0  =  0.023306095434056148;
  
  f10 = f102*chi1*chi1 + f101*chi1 + f100;
  f11 = f112*chi1*chi1 + f111*chi1 + f110;
  f12 = f122*chi1*chi1 + f121*chi1 + f120;
  f20 = f202*chi1*chi1 + f201*chi1 + f200;
  f21 = f212*chi1*chi1 + f211*chi1 + f210;
  f22 = f222*chi1*chi1 + f221*chi1 + f220;

  f1 = f10*nu + f11*nu*nu + f12*nu*nu*nu;
  f2 = f20*nu + f21*nu*nu + f22*nu*nu*nu;

  Op[1] = -1*( ( 1 + (f1*lambda + f2*lambda*lambda) ) / ( (1 + lambda*f0*nu)*(1 + lambda*f0*nu) ) );

  /* l = 3, m = 3 */
  f113 = -19.588729051102252;
  f112 = -43.22678085367127;
  f111 = -27.92748291741765;
  f110 = -5.466513479762885;
  f123 = 128.1454115145983;
  f122 = 265.3361448594152;
  f121 = 167.97216621019794;
  f120 = 33.22919845056078;
  f213 = -0.02463250915949673;
  f212 = -0.038803444817938;
  f211 = -0.026023558342637802;
  f210 = -0.006942576670153899;
  f223 = 0.17994082791820154;
  f222 = 0.31139003107962354;
  f221 = 0.20872506335227003;
  f220 = 0.051899816703641084;
  f321 = 6.299581311689821;
  f320 = 4.7821277096672485;
    
  f11 = f113*chi1*chi1*chi1 + f112*chi1*chi1 + f111*chi1 + f110;
  f12 = f123*chi1*chi1*chi1 + f122*chi1*chi1 + f121*chi1 + f120;
  f21 = f213*chi1*chi1*chi1 + f212*chi1*chi1 + f211*chi1 + f210;
  f22 = f223*chi1*chi1*chi1 + f222*chi1*chi1 + f221*chi1 + f220;
  f32 = f321*chi1 + f320;
    
  f1 = f11*nu + f12*nu*nu; 
  f2 = f21*nu + f22*nu*nu;
  f3 = f32*nu*nu;

  Op[4] = ( 1 + (f1*lambda + f2*lambda*lambda ) ) / ( (1 + lambda*f3*f3 )*(1 + lambda*f3*f3 ) );
  if (Op[4]<0) {Op[4]=1.;}

  /* l = 4, m = 4 */
  f113 = -0.10881642042836596;
  f112 = -0.15513315145142748;
  f111 = 0.10299713694482752;
  f110 = -0.068436178916252;
  f123 = 0.7197962146049899;
  f122 = 0.9224408403513145;
  f121 = -0.6218941042566677;
  f120 = 0.6079384879824817;
  f213 = 0.0010700402803381132;
  f212 = 0.0008737110930575237;
  f211 = -0.0001569042738746961;
  f210 = -0.0003180098977010195;
  f223 = -0.0060384494590680245;
  f222 = -0.0048973961553262875;
  f221 = 0.0006454895547057534;
  f220 = 0.0019251960072054005;
  f321 = -0.12849210212491663;
  f320 = 1.7682780795364919;
    
  f11 = f113*chi1*chi1*chi1 + f112*chi1*chi1 + f111*chi1 + f110;
  f12 = f123*chi1*chi1*chi1 + f122*chi1*chi1 + f121*chi1 + f120;
  f21 = f213*chi1*chi1*chi1 + f212*chi1*chi1 + f211*chi1 + f210;
  f22 = f223*chi1*chi1*chi1 + f222*chi1*chi1 + f221*chi1 + f220;
  f32 = f321*chi1 + f320;
    
  f1 = f11*nu + f12*nu*nu; 
  f2 = f21*nu + f22*nu*nu;
  f3 = f32*nu*nu; 

  Op[8] = ( 1 + (f1*lambda + f2*lambda*lambda ) ) / ( (1 + lambda*f3*f3 )*(1 + lambda*f3*f3 ) );
  if (Op[8]<0) {Op[8]=1.;}

  
  /* derivative omega */
  double e102,e101,e100,e112,e111,e110,e122,e121,e120,e202,e201,e200,e212,e211,e210,e222,e221,e220;
  double e113,e123,e213,e223,e311,e310,e321,e320,e313,e312,e323,e322,e421,e420;
  double e10,e11,e12,e20,e21,e22,e31,e32,e42;
  double e0,e1,e2,e3,e4;

  /* l = 2, m = 1 */
  e113 = 109.90379753602357;
  e112 = 282.7345236349762;
  e111 = 177.70606252810788;
  e110 = 30.091199861511928;
  e123 = 749.4113683426934;
  e122 =  268.3051253328002;
  e121 = -175.71115696673235;
  e120 = -47.81664204077215;
  e213 = -5.600074527093804;
  e212 = -3.5152273527724467;
  e211 = 0.4334580249920891;
  e210 = 0.3203529758149393;
  e223 = 26.65854564108471;
  e222 = 16.965103396559954;
  e221 = -1.8949788120917779;
  e220 = -1.5067162366006845;
  e313 = 181.68764838355855;
  e312 = 46.48101845933214;
  e311 = -77.46914070876984;
  e310 =  -20.618000090591245;
  e323 = -827.1664409469339;
  e322 = -198.70098787921546;
  e321 = 379.24240378523564;
  e320 = 103.76437073875854;

  e11 = e113*chi1*chi1*chi1 + e112*chi1*chi1 + e111*chi1 + e110;
  e12 = e123*chi1*chi1*chi1 + e122*chi1*chi1 + e121*chi1 + e120;
  e21 = e213*chi1*chi1*chi1 + e212*chi1*chi1 + e211*chi1 + e210;
  e22 = e223*chi1*chi1*chi1 + e222*chi1*chi1 + e221*chi1 + e220;
  e31 = e313*chi1*chi1*chi1 + e312*chi1*chi1 + e311*chi1 + e310;
  e32 = e323*chi1*chi1*chi1 + e322*chi1*chi1 + e321*chi1 + e320;

  e1 = e11*nu + e12*nu*nu;
  e2 = e21*nu + e22*nu*nu;
  e3 = e31*nu + e32*nu*nu;

  dOp[0] = ( 1 + (e1*lambda + e2*lambda*lambda  ) ) / ( (1 + lambda*e3*e3 )*(1 + lambda*e3*e3 ) );

  /* l = 2, m = 2 */
  e102 = -0.0026209606150563616; 
  e101 = 5.121199315158528e-05;
  e100 = 0.03580195978164114;
  e112 = -0.31662958435109995;
  e111 = -0.5067440136512827;
  e110 = 0.48032840592888776;
  e122 = 2.007792466450656;
  e121 = 1.115530157663557;
  e120 = -2.3169394297754167;
  e202 = 3.276112836329881e-07;
  e201 = 5.226090580513704e-07;
  e200 = -3.2932886267689516e-07;
  e212 = 0.0011459688540299546;
  e211 = -0.001315350180849681;
  e210 = 0.00024114701127108795;
  e222 = -0.0052276257814133755;
  e221 = 0.006309356574634241;
  e220 = -0.0012963692153788548;
  e0  =  0.01733398003688234;
  
  e10 = e102*chi1*chi1 + e101*chi1 + e100;
  e11 = e112*chi1*chi1 + e111*chi1 + e110;
  e12 = e122*chi1*chi1 + e121*chi1 + e120;
  e20 = e202*chi1*chi1 + e201*chi1 + e200;
  e21 = e212*chi1*chi1 + e211*chi1 + e210;
  e22 = e222*chi1*chi1 + e221*chi1 + e220;

  e1 = e10*nu + e11*nu*nu + e12*nu*nu*nu;
  e2 = e20*nu + e21*nu*nu + e22*nu*nu*nu;

  dOp[1] = ( ( 1 + (e1*lambda + e2*lambda*lambda) ) / ( (1 + lambda*e0*nu)*(1 + lambda*e0*nu) ) );

  if (dOp[1]<0.){
    dOp[1]=1.;
    Op[1] =1.;
  }

  /* l = 3, m = 3 */
  e113 = -1592054853.5290756;
  e112 = 1193131216.0498104;
  e111 = 39838533.1798002;
  e110 = -158420220.57876483;
  e123 = 8991032546.305927;
  e122 = -6770511898.425056;
  e121 = -185433414.46684813;
  e120 = 888834576.8629049;
  e213 = 7920499.5760512855;
  e212 = -6043560.756636857;
  e211 = 95758.20206124203;
  e210 = 665869.5430042462;
  e223 = -44608471.17271389;
  e222 = 33975409.50367452;
  e221 = -530037.7676815926;
  e220 = -3770237.4332308765;
  e313 = -6958.454707644981;
  e312 = 5179.7378134802575;
  e311 = -431.5691613088688;
  e310 = -398.967911171984;
  e323 = 38203.77095038259;
  e322 = -27483.22126180374;
  e321 = 1539.8252232858072;
  e320 = 2427.971914618684;
  e421 = -11.410927119363691;
  e420 = 9.884174594022015;
    
  e11 = e113*chi1*chi1*chi1 + e112*chi1*chi1 + e111*chi1 + e110;
  e12 = e123*chi1*chi1*chi1 + e122*chi1*chi1 + e121*chi1 + e120;
  e21 = e213*chi1*chi1*chi1 + e212*chi1*chi1 + e211*chi1 + e210;
  e22 = e223*chi1*chi1*chi1 + e222*chi1*chi1 + e221*chi1 + e220;
  e31 = e313*chi1*chi1*chi1 + e312*chi1*chi1 + e311*chi1 + e310;
  e32 = e323*chi1*chi1*chi1 + e322*chi1*chi1 + e321*chi1 + e320;
  e42 = e421*chi1 + e420;
    
  e1 = e11*nu + e12*nu*nu; 
  e2 = e21*nu + e22*nu*nu;
  e3 = e31*nu + e32*nu*nu;
  e4 = e42*nu*nu;

  dOp[4] = ( 1 + (e1*lambda + e2*lambda*lambda + e3*lambda*lambda*lambda) ) / ( (1 + lambda*lambda*e4*e4 )*(1 + lambda*lambda*e4*e4 ) );

  /* l = 4, m = 4 */
  e112 = 179210.69776361363;
  e111 = -135573.3404546459;
  e110 = 20893.938274313627;
  e122 = -564011.9694244653;
  e121 = 239102.42671224466;
  e120 = 36894.6182860516;
  e212 = -286.55834017745804;
  e211 = 129.8284659110896;
  e210 = 13.42318091379192;
  e222 = 1422.0780987349697;
  e221 = -479.76042078095696;
  e220 = -157.05526942865345;
  e321 = -3268.575047986779;
  e320 = 2250.3024799305063;
    
  e11 = e112*chi1*chi1 + e111*chi1 + e110;
  e12 = e122*chi1*chi1 + e121*chi1 + e120;
  e21 = e212*chi1*chi1 + e211*chi1 + e210;
  e22 = e222*chi1*chi1 + e221*chi1 + e220;
  e32 = e321*chi1 + e320;
    
  e1 = e11*nu + e12*nu*nu; 
  e2 = e21*nu + e22*nu*nu;
  e3 = e32*nu*nu; 

  dOp[8] = ( 1 + (e1*lambda + e2*lambda*lambda ) ) / ( (1 + lambda*e3*e3 ) );

  /* Amplitude */
  double c102,c101,c100,c113,c112,c111,c110,c122,c123,c121,c120,c202,c201,c200,c213,c212,c211,c210,c223,c222,c221,c220;
  double c132,c131,c130,c232,c231,c230,c311,c310,c133,c233,c321,c320,c312,c322,c411,c410,c421,c420;
  double c10,c11,c12,c20,c21,c22,c13,c23,c30,c31,c32,c41,c42;
  double c1,c2,c3,c4;

  /* l = 2, m = 1 */
  c113 = -93.42972864623421;
  c112 = -1229.384079825144;
  c111 = 1070.5956514039542;
  c110 = -178.81408417308748;
  c123 = 1866.8412412088521;
  c122 = 8972.415771021775;
  c121 = -9543.800468423478;
  c120 = 2089.7290556140033;
  c213 = 0.8669677161029467;
  c212 = 1.3680430387869333;
  c211 = -1.3923206119416818;
  c210 = 0.1540730801885621;
  c223 = -7.058523132126191;
  c222 = -6.568053356330434;
  c221 = 8.635081951736485;
  c220 = -1.28381717345134;
  c311 = -296.5534943063397;
  c310 = 175.54999284090823;

  c11 = c113*chi1*chi1*chi1 + c112*chi1*chi1 + c111*chi1 + c110;
  c12 = c123*chi1*chi1*chi1 + c122*chi1*chi1 + c121*chi1 + c120;
  c21 = c213*chi1*chi1*chi1 + c212*chi1*chi1 + c211*chi1 + c210;
  c22 = c223*chi1*chi1*chi1 + c222*chi1*chi1 + c221*chi1 + c220;
  c31 = c311*chi1 + c310;

  c1 = c11*nu + c12*nu*nu; 
  c2 = c21*nu + c22*nu*nu;
  c3 = c31*nu*nu; 

  Ap[0] = ( 1 + (c1*lambda + c2*lambda*lambda ) ) / ( (1 + lambda*c3*c3 ) );
  if (Ap[0]<0.){Ap[0]=1.;}
  if (Ap[0]>1.){Ap[0]=1.;}

  /* l = 2, m = 2 */
  c102 = 0.0009785096527311836; 
  c101 = 0.0010143191550383841;
  c100 = 0.06298925167477198;
  c112 = 0.20384994126329886;
  c111 = -2.5851799086634304;
  c110 = 1.7404500991629586;
  c122 = -0.1315604591952073;
  c121 = 7.110295022425293;
  c120 = -3.495497940497986;
  c202 = -5.034425005770115e-07;
  c201 = 3.2828002996716215e-07;
  c200 = 3.34808127112953e-07;
  c212 = 0.00023990151371606222;
  c211 = -0.00011455874202441466;
  c210 = -1.738883897518922e-05;
  c222 = -0.001185391692557779;
  c221 = 0.0011437351997991667;
  c220 = -0.0006164680289672201;
  c31  = 0.06456023857357258;
  c32  = 1.1422176143066296;
  
  c10 = c102*chi1*chi1 + c101*chi1 + c100;
  c11 = c112*chi1*chi1 + c111*chi1 + c110;
  c12 = c122*chi1*chi1 + c121*chi1 + c120;
  c20 = c202*chi1*chi1 + c201*chi1 + c200;
  c21 = c212*chi1*chi1 + c211*chi1 + c210;
  c22 = c222*chi1*chi1 + c221*chi1 + c220;

  c1 = c10*nu + c11*nu*nu + c12*nu*nu*nu;
  c2 = c20*nu + c21*nu*nu + c22*nu*nu*nu;
  c3 = c31*nu + c32*nu*nu;

  Ap[1] =  ( 1 + (c1*lambda + c2*lambda*lambda) ) / ( 1 + lambda*c3  );
  if (Ap[1]<0.){Ap[1]=1.;}
  if (Ap[1]>1.){Ap[1]=1.;}

  /* l = 3, m = 3 */
  c113 = -86.20386287795816;
  c112 = -10.699331077074865;
  c111 = 17.246970292171067;
  c110 = 5.102406475415944;
  c123 = 849.0739042449288;
  c122 = 88.2894277432953;
  c121 = -167.72590918207533;
  c120 = -44.32732463839548;
  c133 = -2086.682310126456;
  c132 = -177.3505956893603;
  c131 = 409.4698953574064;
  c130 = 96.41556256518713;
  c213 = 0.03733441603389203;
  c212 = -0.023131489516298686;
  c211 = -0.03212617167840329;
  c210 = 0.016757348842393934;
  c223 = -0.32630630833371593;
  c222 = 0.274104930253967;
  c221 = 0.3217506163567745;
  c220 = -0.15060952414887815;
  c233 = 0.7201771670241224;
  c232 = -0.767024732086749;
  c231 = -0.7987937160274846;
  c230 = 0.33851257061927315;
  c311 = -1.0964469876953258;
  c310 = 2.414952851569057;
  c321 = 6.91754931598061;
  c320 = -9.844527125356167;
    
  c11 = c113*chi1*chi1*chi1 + c112*chi1*chi1 + c111*chi1 + c110;
  c12 = c123*chi1*chi1*chi1 + c122*chi1*chi1 + c121*chi1 + c120;
  c13 = c133*chi1*chi1*chi1 + c132*chi1*chi1 + c131*chi1 + c130;
  c21 = c213*chi1*chi1*chi1 + c212*chi1*chi1 + c211*chi1 + c210;
  c22 = c223*chi1*chi1*chi1 + c222*chi1*chi1 + c221*chi1 + c220;
  c23 = c233*chi1*chi1*chi1 + c232*chi1*chi1 + c231*chi1 + c230;
  c31 = c311*chi1 + c310;
  c32 = c321*chi1 + c320;
    
  c1 = c11*nu + c12*nu*nu + c13*nu*nu*nu;
  c2 = c21*nu + c22*nu*nu + c23*nu*nu*nu;
  c3 = c31*nu + c32*nu*nu;

  Ap[4] = ( 1 + (c1*lambda + c2*lambda*lambda ) ) / ( (1 + lambda*c3*c3 )*(1 + lambda*c3*c3 ) );
  if (Ap[4]<0.){Ap[4]=1.;}

  /* l = 4, m = 4 */
  c112 = -383999414.8996635;
  c111 = -440539181.38602424;
  c110 = -110989340.22848652;
  c122 = 2235094349.7057714;
  c121 = 2590824770.6686583;
  c120 = 664021933.0407408;
  c212 = 2141767.8243446765;
  c211 = 2393114.3223998966;
  c210 = 627456.6514837451;
  c222 = -12511651.64787277;
  c221 = -14180990.752931684;
  c220 = -3801120.272961312;
  c312 = -2873.1756871481566;
  c311 = -3474.4177866366244;
  c310 = -1038.9430866381745;
  c322 = 16977.480734824097;
  c321 = 20957.257721739938;
  c320 = 6426.048679948954;
  c411 = -3.046924825109896;
  c410 = -1.9714971745008207;
  c421 = 22.014330706283033;
  c420 = 14.81026562279348;
    
  c11 = c112*chi1*chi1 + c111*chi1 + c110;
  c12 = c122*chi1*chi1 + c121*chi1 + c120;
  c21 = c212*chi1*chi1 + c211*chi1 + c210;
  c22 = c222*chi1*chi1 + c221*chi1 + c220;
  c31 = c312*chi1*chi1 + c311*chi1 + c310;
  c32 = c322*chi1*chi1 + c321*chi1 + c320;
  c41 = c411*chi1 + c410;
  c42 = c421*chi1 + c420;
    
  c1 = c11*nu + c12*nu*nu; 
  c2 = c21*nu + c22*nu*nu;
  c3 = c31*nu + c32*nu*nu; 
  c4 = c41*nu + c42*nu*nu;

  Ap[8] = ( 1 + (c1*lambda + c2*lambda*lambda + c3*lambda*lambda*lambda) ) / ( (1 + lambda*lambda*c4 )*(1 + lambda*lambda*c4 ) );
  if (Ap[8]<0.){Ap[8]=1.;}
  if (Ap[8]>1.){Ap[8]=1.;}
  
  /* derivative amplitude */
  double b101,b100,b111,b110,b201,b200,b211,b210,b301,b300,b311,b310,b113,b123,b213,b223;
  double b112,b122,b121,b120,b132,b131,b130,b212,b222,b221,b220,b232,b231,b230,b133,b233,b321,b320,b312,b322;
  double b10,b11,b20,b21,b30,b31,b12,b13,b22,b23,b32;
  double b1,b2,b3;

  /* l = 2, m = 1 */
  b112 = -368789.825574789;
  b111 = 333038.7324440759;
  b110 = -133958.7873246825;
  b122 = 3690938.226536074;
  b121 = -3389666.6770304046;
  b120 = 1371584.2975333077;
  b132 = -9062348.686131006;
  b131 = 8510705.839035401;
  b130 = -3462442.3106556754;
  b212 = -742.8322268026368;
  b211 = -168.19680873964367;
  b210 = 420.00912816902934;
  b222 = 8927.666668277603;
  b221 = 988.1724219968736;
  b220 = -4203.188483558627;
  b232 = -25826.98653518626;
  b231 = -1133.0371849440526;
  b230 = 10545.571404324717;
  b31  = -9.334649682992511;
  b30 = 5.349037446237097;

  b11 = b112*chi1*chi1 + b111*chi1 + b110;
  b12 = b122*chi1*chi1 + b121*chi1 + b120;
  b13 = b132*chi1*chi1 + b131*chi1 + b130;
  b21 = b212*chi1*chi1 + b211*chi1 + b210;
  b22 = b222*chi1*chi1 + b221*chi1 + b220;
  b23 = b232*chi1*chi1 + b231*chi1 + b230;

  b1 = b11*nu + b12*nu*nu + b13*nu*nu*nu;
  b2 = b21*nu + b22*nu*nu + b23*nu*nu*nu;
  b3 = b31*chi1 + b30;

  dAp[0] = ( ( 1 + (b1*lambda + b2*lambda*lambda) ) / ( (1 + lambda*lambda*(b3*nu)*(b3*nu) ) ) );

  /* l = 2, m = 2 */
  b101 = 0.006692383187572945;
  b100 = -0.0029099423327064626;	
  b111 = 0.09298127163261365;
  b110 = 0.15371784703359265;
  b201 = -1.8017891913004304e-06;
  b200 = 8.63563033494384e-07;
  b211 = 1.1895204934755798e-05;
  b210 = -2.1308452612548277e-05;
  b301 = 0.6202485959739941;
  b300 = 0.18408910774027512;
  b311 = -1.6022272014694021;
  b310 = 0.44432056462055053;

  b10 = b101*chi1 + b100;
  b11 = b111*chi1 + b110;
  b20 = b201*chi1 + b200;
  b21 = b211*chi1 + b210;
  b30 = b301*chi1 + b300;
  b31 = b311*chi1 + b310;

  b1 = b10*nu + b11*nu*nu;
  b2 = b20*nu + b21*nu*nu;
  b3 = b30*nu + b31*nu*nu;

  dAp[1] = -1*( 1 + lambda*b1 + lambda*lambda*b2 ) / (( 1 + lambda*b3*b3 )*( 1 + lambda*b3*b3 ));

  /* l = 3, m = 3 */
  b113 = -104532308.71274745;
  b112 = -28052057.502293073;
  b111 = -10334156.319118189;
  b110 = 7773885.670298806;
  b123 = 1129395588.3281617;
  b122 = 285802735.60785395;
  b121 = 98375791.70591767;
  b120 = -78213782.01803419;
  b133 = -3022996822.7874002;
  b132 = -717768537.4658437;
  b131 = -229117863.72980285;
  b130 = 194180576.92894575;
  b213 = 119189.0360537401;
  b212 = 71501.19416486716;
  b211 = 32974.62770801159;
  b210 = -20045.02496446706;
  b223 = -1283346.8875412862;
  b222 = -739058.2125105071;
  b221 = -329248.2986502753;
  b220 = 203540.475443693;
  b233 = 3416520.679349472;
  b232 = 1886771.5391024328;
  b231 = 813712.6677610069;
  b230 = -511038.39130576607;
  b321 = -27.251430551104228;
  b320 = 48.99759965335871;
    
  b11 = b113*chi1*chi1*chi1 + b112*chi1*chi1 + b111*chi1 + b110;
  b12 = b123*chi1*chi1*chi1 + b122*chi1*chi1 + b121*chi1 + b120;
  b13 = b133*chi1*chi1*chi1 + b132*chi1*chi1 + b131*chi1 + b130;
  b21 = b213*chi1*chi1*chi1 + b212*chi1*chi1 + b211*chi1 + b210;
  b22 = b223*chi1*chi1*chi1 + b222*chi1*chi1 + b221*chi1 + b220;
  b23 = b233*chi1*chi1*chi1 + b232*chi1*chi1 + b231*chi1 + b230;
  b32 = b321*chi1 + b320;
    
  b1 = b11*nu + b12*nu*nu + b13*nu*nu*nu;
  b2 = b21*nu + b22*nu*nu + b23*nu*nu*nu;
  b3 = b32*nu*nu;

  dAp[4] = ( 1 + (b1*lambda + b2*lambda*lambda ) ) / ( (1 + lambda*b3*b3 )*(1 + lambda*b3*b3 ) );
  if (dAp[4]<0.){dAp[4]=Ap[4];}

  /* l = 4, m = 4 */
  b112 = 13904.899442924849;
  b111 = 2188.159558986493;
  b110 = -3428.942380746902;
  b122 = -203682.21284236634;
  b121 = -23100.79333975182;
  b120 = 56591.65825035908;
  b132 = 651883.3845677918;
  b131 = 59884.13132225108;
  b130 = -191074.31755052076;
  b212 = 170.04957281900852;
  b211 = 3.7413046389109774;
  b210 = -59.036266119590664;
  b222 = -1518.7747992805148;
  b221 = -41.139279970982884;
  b220 = 523.3990858776148;
  b232 = 3378.001071520739;
  b231 = 113.2252224885313;
  b230 = -1152.993306319554;
  b312 = 6.2737580723088335;
  b311 = 0.3047913280091431;
  b310 = -2.2403560162621146;
  b322 = -30.425718021789745;
  b321 = -1.3870049232772124;
  b320 = 10.95952129988658;
     
  b11 = b112*chi1*chi1 + b111*chi1 + b110;
  b12 = b122*chi1*chi1 + b121*chi1 + b120;
  b13 = b132*chi1*chi1 + b131*chi1 + b130;
  b21 = b212*chi1*chi1 + b211*chi1 + b210;
  b22 = b222*chi1*chi1 + b221*chi1 + b220;
  b23 = b232*chi1*chi1 + b231*chi1 + b230;
  b31 = b312*chi1*chi1 + b311*chi1 + b310;
  b32 = b322*chi1*chi1 + b321*chi1 + b320;
    
  b1 = b11*nu + b12*nu*nu + b13*nu*nu*nu;
  b2 = b21*nu + b22*nu*nu + b23*nu*nu*nu;
  b3 = b31*nu + b32*nu*nu;

  dAp[8] = ( 1 + (b1*lambda + b2*lambda*lambda ) ) / ( (1 + lambda*lambda*b3*b3 )*(1 + lambda*lambda*b3*b3 ) );

  // Peak values for all modes
  for (int k=0; k<KMAX; k++) {
    omg_tmp[k] = Op[k]*omg_tmp[k];
    domg_tmp[k] = dOp[k]*domg_tmp[k];  
    A_tmp[k] = Ap[k]*A_tmp[k];
    dA_tmp[k] = dAp[k]*dA_tmp[k];
  }
	  
}

void eob_nqc_point_test(double Mbh, double c1A, double c2A, double c3A, double c4A, 
			    double c1phi, double c2phi, double c3phi, double c4phi,
			    double alpha1, double omega1,
			    double *A_tmp, double *dA_tmp, double *d2A_tmp, double *omg_tmp, double *domg_tmp, double *d2omg_tmp)
{
  
  double tau, Mbh2, x, x2, dA_tmp1, dA_tmp2, omg_tmp1, omg_tmp2, domg_tmp_n1,domg_tmp_n2,domg_tmp_d1,domg_tmp_d2;

  /* the time variable in the post-peak template is given in units of Mbh*/
  tau  = 2./Mbh;
  Mbh2 = SQ(Mbh);
  
  *A_tmp = exp(-alpha1*tau)*(c1A*tanh(c2A*tau + c3A) + c4A);
  
  double fact = cosh(c3A+ c2A*tau);
  dA_tmp1     = c1A*c2A*exp(-alpha1*tau)/(fact*fact);
  dA_tmp2     = -alpha1*exp(-alpha1*tau)*(c4A + c1A*tanh(c3A+ c2A*tau));
  *dA_tmp     = (dA_tmp1 + dA_tmp2)/Mbh;

  double fact2 = 1.0/cosh(c2A*tau+c3A);
  double d2A_tmp1 = -c1A*c2A*exp(-alpha1*tau)*(fact2*fact2)*(alpha1+2*c2A*tanh(c2A*tau+c3A));
  double d2A_tmp2 = alpha1*alpha1*exp(-alpha1*tau)*(c1A*tanh(c2A*tau+c3A)+c4A);
  double d2A_tmp3 = -alpha1*exp(-alpha1*tau)*c1A*c2A * (fact*fact);
  *d2A_tmp = (d2A_tmp1 + d2A_tmp2 + d2A_tmp3)/Mbh2;
  
  x        = exp(-c2phi*tau);
  x2       = x*x;
  omg_tmp1 = c1phi*(-2.*c2phi*c4phi*exp(-2.*c2phi*tau) - c2phi*c3phi*x);
  omg_tmp2 = 1 + c4phi*exp(-2.*c2phi*tau) + c3phi*x;
  *omg_tmp = (omg_tmp1/omg_tmp2)/Mbh + omega1/Mbh2;
	  
  domg_tmp_n1 = c1phi*SQ(-2.*c2phi*c4phi*x2 - c2phi*c3phi*x);
  domg_tmp_d1 = SQ(1 + c4phi*x2 + c3phi*x);
  domg_tmp_n2 = c1phi*(4.*SQ(c2phi)*c4phi*x2 + SQ(c2phi)*c3phi*x);
  domg_tmp_d2 = 1 + c4phi*x2 + c3phi*x;
  *domg_tmp   = -(domg_tmp_n1/domg_tmp_d1 - domg_tmp_n2/domg_tmp_d2)/Mbh2;

  double d2omg_tmp1 = (-c2phi*c3phi*x-8*c2phi*c4phi*x2) / (1+c3phi*x+c4phi*x2);
  double d2omg_tmp2 = (-1)*(c3phi*x+4*c4phi*x2)*(-c2phi*c3phi*x-2*c4phi*c2phi*x2) / ((1+c3phi*x+c4phi*x2)*(1+c3phi*x+c4phi*x2));
  double d2omg_tmp3 = 2*(c3phi*x+2*c4phi*x2) / (1+c3phi*x+c4phi*x2);
  double d2omg_tmp4 = (-c2phi*c3phi*x-4*c2phi*c4phi*x2) / (1+c3phi*x+c4phi*x2);
  double d2omg_tmp5 = (-1)*(c3phi*x+2*c4phi*x2)*(-c2phi*c3phi*x-2*c4phi*c2phi*x2);
  double d2omg_tmp6 = d2omg_tmp4 + d2omg_tmp5;
  *d2omg_tmp = c1phi*c2phi*c2phi*( (d2omg_tmp1+d2omg_tmp2) - d2omg_tmp3*d2omg_tmp6 ) / (Mbh2*Mbh);
}

/** p_phi LSO */
double pph_lso_orbital(const double nu)
{
  return 3.46343 - 0.774482*nu -0.69200*nu*nu; 
} 

/** p_phi LSO for the spinning case; calibrated in the regime |chi|<0.5, arbitrary nu. Residuals up to 1% */
double pph_lso_spin(const double nu, const double a0)
{
  /** Schwarzschild limit */
  double Y_0 = sqrt(12);

  double b_nu_1  =  9.767215568147924 ; 
  double b_nu_2  = -6.519453488314084 ; 
  double c_nu_1  =  9.999832399778665 ; 
  double c_nu_2  = -3.924142989771926 ; 
  double b_chi_1 = -2.5829502974100746; 
  double b_chi_2 =  1.1894184903673561; 
  double c_chi_1 = -2.298897170468296 ; 
  double c_chi_2 =  0.6259081961126761;

  double nu2 = nu * nu ;
  double a02 = a0 * a0 ;

  double expansion_num_nu  = 1. +  b_nu_1  * nu + b_nu_2  * nu2 ;
  double expansion_den_nu  = 1. +  c_nu_1  * nu + c_nu_2  * nu2 ; 
  double expansion_num_chi = 1. +  b_chi_1 * a0 + b_chi_2 * a02 ; 
  double expansion_den_chi = 1. +  c_chi_1 * a0 + c_chi_2 * a02 ;

  double nu_pade  = expansion_num_nu  / expansion_den_nu ;
  double chi_pade = expansion_num_chi / expansion_den_chi ;

  double fit_model = Y_0 * nu_pade * chi_pade ;
  
  return fit_model; 
} 

/** Horizon radius, Nagar fit 13/07/20 */
double horizon_radius(const double nu)
{
  return 2.02156 -6.85222*nu + 13.2842*nu*nu - 72.9496*nu*nu*nu; 
}

/** 
 *  Function: eob_nqc_timeshift_bhns
 *  --------------------------------
 *  \details This function computes the time-shift for (2,2) tNQC with respect
 *   to the peak of the EOB pure orbital frequency, from
 *   https://arxiv.org/abs/1806.01772 and adapted for BHNS cases.
 *  From Gonzalez+ in prep.
 *   
 *   @param[in] nu   : symmmetric mass ratio
 *   @param[in] chi1 : dimensionless spin of body 1
 *
 *   @return DeltaT_nqc
 */
double eob_nqc_timeshift_bhns(double nu, double chi1)
{

  double DeltaT_nqc = 4.0;    

  if (EOBPars->d_delta_t_nqc) DeltaT_nqc += EOBPars->d_delta_t_nqc;
    
  return DeltaT_nqc;  
}

/** 
 *  Function: eob_nqc_deltat_lm_bhns
 *  ---------------------------
 *  \details This function computes the time-shift of the peak of the modes w.r.t
 *   the (2,2) mode, as defined in Eq. 3.50 of https://arxiv.org/1904.09550
 *   From Gonzalez+ in prep.
 *   
 *   @param[out] Dt_lm : KMAX dimensional array of time-shifts 
 *
 */
void eob_nqc_deltat_lm_bhns(double *Dt_lm)
{
  
  const double nu    = EOBPars->nu;
  const double nu2   = SQ(nu);
  const double nu3   = nu2*nu;
  const double nu4   = SQ(nu2);
  const double X1    = EOBPars->X1;
  const double X2    = EOBPars->X2;
  const double X12   = X1 - X2;
  const double X12_2 = SQ(X12);
  const double X12_3 = X12_2*X12;

  const double chi1  = EOBPars->chi1;
  const double chi2  = EOBPars->chi2;
  const double aK    = X1*chi1 + X2*chi2;
  const double aK2   = SQ(aK);
  const double a12   = X1*chi1 - X2*chi2;
  const double a12_2 = SQ(a12);
  const double Shat  = 0.5*(aK + a12*X12);
  const double Shat2 = SQ(Shat);
  const double lambda = EOBPars->LambdaBl2;
  
  const int usespins = EOBPars->use_spins;  

  double Dt_lm_TP[KMAX], dt_bhns[KMAX];
  for (int k=0; k<KMAX; k++) {
    Dt_lm[k] = Dt_lm_TP[k] = 0;
    dt_bhns[k] = 1.;
  }

  Dt_lm_TP[0]  = 11.75925;
  Dt_lm_TP[2]  = 12.9338;
  Dt_lm_TP[3]  = 9.16665;
  Dt_lm_TP[4]  = 3.42593;
  Dt_lm_TP[5]  = 13.1116;
  Dt_lm_TP[6]  = 11.66665;
  Dt_lm_TP[7]  = 9.53705;
  Dt_lm_TP[8]  = 5.27778;
  Dt_lm_TP[13] = 6.561811;

    
    double Dt_21_orb, Dt_32_orb, Dt_33_orb, Dt_42_orb, Dt_43_orb, Dt_44_orb;
    double Dt_32_S, Dt_33_S1, Dt_33_S2, Dt_42_S, Dt_43_S, Dt_44_S;
    double c11Dtlm, c12Dtlm, c21Dtlm, c22Dtlm, b1Dtlm, b2Dtlm;
    double a1Dt, a2Dt, a3Dt, a4Dt;
    double b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15;
    double n1, n2, d1, d2;

    /* (l = 2, m = 1) */
    Dt_21_orb	= (Dt_lm_TP[0]*X12_2 + 4*6.6264*nu)*(1-2.0728*nu*X12);

    c11Dtlm   = -1976.13; 
    c12Dtlm   = 3719.88;
    c21Dtlm   = -2545.41; 
    c22Dtlm   = 5277.62;
    b1Dtlm    = 0.0472289;
    b2Dtlm    = 0.115583;
    b1        = (b1Dtlm + c11Dtlm*X12)/(1 + c12Dtlm*X12);
    b2        = (b2Dtlm + c21Dtlm*X12)/(1 + c22Dtlm*X12);
    Dt_lm[0]  = (1 + b1*aK + b2*aK2)*Dt_21_orb;

    /* (l = 3, m = 2) */
    a1Dt      = -11.3497;
    a2Dt      = 32.9144; 
    a3Dt      = -8.36579;
    a4Dt      = 20.1017;
    Dt_32_orb = Dt_lm_TP[3]*(1 + a1Dt*nu+ a2Dt*nu2)/(1 + a3Dt*nu+ a4Dt*nu2);

   
      
      b1 = +2.497188;
      b2 = -7.532596;
      b3 = +4.645986;
      b4 = -3.652524;
      b5 = +3.398687;
      b6 = +7.054185;
      b7 = -12.260185;
      b8 = +5.724802;
      b9 = -3.242611;
      b10 = +2.714232;
      b11 = +2.614565;
      b12 = -9.507583;
      b13 = +7.321586;
      b14 = -3.937568;
      b15 = +4.584970;
      Dt_32_S  = (1 + (+0.15477 + b1*X12 + b2*X12_2 + b3*X12_3)/(1 + b4*X12 + b5*X12_2)*Shat + (-0.755639 + b6*X12 + b7*X12_2 + b8*X12_3)/(1 + b9*X12 + b10*X12_2)*Shat2)/(1 + (+0.21816 + b11*X12 + b12*X12_2 + b13*X12_3)/(1 + b14*X12 + b15*X12_2)*Shat);

  

    Dt_lm[3] = Dt_32_orb*Dt_32_S/Dt_lm_TP[3];

    /* (l = 3, m = 3) */
    Dt_33_orb = Dt_lm_TP[4]*(1 + 0.183349*nu+ 4.22361*nu2);
    Dt_33_S1  = 1 + (-0.49791 - 1.9478*nu)/(1 + 13.9828*nu)*Shat + (-0.18754 + 1.25084*nu)/(1 - 3.41811*nu)*Shat2;
    Dt_33_S2  = 1 + (-1.07291 - 1043.15*nu)/(1 + 1033.85*nu)*Shat;
    Dt_lm[4]  = Dt_33_orb*Dt_33_S1/Dt_33_S2; 

    /* (l = 4, m = 2) */
    Dt_42_orb  = Dt_lm_TP[6]*(1 - 9.8446172795*nu + 23.3229430582*nu2)/(1 - 5.7604819848*nu + 7.1217930024*nu2);
    
    if (nu < 6./25.) {
      b1       = 0.;
      b2       = 0.;
      b3       = 0.;
      b4       = 0.;
    } else {
      b1       = +24.604717;
      b2       = -0.808279;
      b3       = +62.471781;
      b4       = +48.340961;
    }
    Dt_42_S  = (1 + (-1.3002045 + b1*nu)/(1 + b2*nu)*Shat)/(1 + (-0.9494348+b3*nu)/(1+b4*nu)*Shat);
    Dt_lm[6] = Dt_42_orb*Dt_42_S;

    /* (l = 4, m = 3) */
    a1Dt      = -11.2377;
    a2Dt      = 38.3177; 
    a3Dt      = -7.29734;
    a4Dt      = 21.4267;
    Dt_43_orb = Dt_lm_TP[7]*(1 + a1Dt*nu + a2Dt*nu2)/(1 + a3Dt*nu + a4Dt*nu2);

    b1      = +3.215984;
    b2      = +42.133767;
    b3      = -9.440398;
    b4      = +35.160776;
    b5      = +1.133942;
    b6      = -10.356311;
    b7      = -6.701429;
    b8      = +10.726960;
    b9      = -6.036207;
    b10     = +67.730599;
    b11     = -3.082275;
    b12     = +11.547917;
    Dt_43_S = (1 + (-1.371832 + b1*nu + b2*nu2)/(1 + b3*nu + b4*nu2)*Shat +((+0.362375+b5*nu+b6*nu2)/(1+b7*nu+b8*nu2))*Shat2)/(1+((-1.0808402+b9*nu+b10*nu2)/(1+b11*nu+b12*nu2))*Shat);
    Dt_lm[7] 		= Dt_43_orb*Dt_43_S;

    /* (l = 4, m = 4) */
    n1        = -8.35574;
    n2        = 17.5288;
    d1        = -6.50259;
    d2        = 10.1575;
    Dt_44_orb = Dt_lm_TP[8]*(1 + n1*nu + n2*nu2)/(1 + d1*nu + d2*nu2);
    b1        = 0.00159701 - 2.28656*X12 + 1.66532*X12_2;
    b2        = -1.14134 - 0.589331*X12 + 0.708784*X12_2;
    Dt_44_S   = (1 + b1*Shat)/(1 + b2*Shat);
    Dt_lm[8]  = Dt_44_orb*Dt_44_S;

    // These fits are taken from the testparticle limit and approixmate the spinning case
    /* (l = 3, m = 1) */
    Dt_lm[2]  = Dt_lm_TP[2]*(1 - 25.615*nu2)/(1 + 0.88803*nu + 16.292*nu2);

    /* (l = 4, m = 1) */	
    Dt_lm[5]  = Dt_lm_TP[5]*(1 - 9.6225*nu + 38.451*nu2)/(1 - 7.7998*nu + 32.405*nu2);	

    /* (l = 5, m = 5) */
    Dt_lm[13] = Dt_lm_TP[13]*(1 - 91.2039*nu + 2556.5123*nu2 - 11325.217*nu3 + 27767.2164*nu4)/(1 - 66.0362*nu + 1762.4169*nu2);


    // BHNS section
    double d113,d112,d111,d110,d123,d122,d121,d120,d213,d212,d211,d210,d223,d222,d221,d220,d311,d310,d321,d320,d411,d410;
    double d312,d322,d132,d131,d130,d232,d231,d230,d332,d331,d330,d133,d233;
    double d11,d12,d21,d22,d31,d30,d32,d41,d13,d23,d33;
    double d3,d4;

    /* (l = 3, m = 2) */
    d113 = -1.5072839537387006;
    d112 = -1.471007307099315;
    d111 = 0.12255034780592934;
    d110 = 0.2751384354900841;
    d123 = 14.314520385981604;
    d122 = 13.520854973779135;
    d121 = -1.9252417774544774;
    d120 = -2.9092668321280564;
    d133 = -32.98083488076057;
    d132 = -29.470273418624345;
    d131 = 7.050737575463469;
    d130 = 7.6697745157484265;
    d213 = 8.290211271634893e-06;
    d212 = 1.972060714610703e-05;
    d211 = 1.3734942172468023e-05;
    d210 = 2.936349065235699e-06;
    d223 = -9.92079106433095e-05;
    d222 = -0.0002326678556004386;
    d221 = -0.00016130469554140212;
    d220 = -3.449512968210697e-05;
    d233 = 0.00029936622381147843;
    d232 = 0.0006889011279026787;
    d231 = 0.00047429001025871097;
    d230 = 0.00010134761352693667;
    d311 = -0.14220847798371816;
    d310 = -0.13330070151638176;
    d321 = 1.0876547994820773;
    d320 = 0.8943849484295967;

    d11 = d113*chi1*chi1*chi1 + d112*chi1*chi1 + d111*chi1 + d110;
    d12 = d123*chi1*chi1*chi1 + d122*chi1*chi1 + d121*chi1 + d120;
    d13 = d133*chi1*chi1*chi1 + d132*chi1*chi1 + d131*chi1 + d130;
    d21 = d213*chi1*chi1*chi1 + d212*chi1*chi1 + d211*chi1 + d210;
    d22 = d223*chi1*chi1*chi1 + d222*chi1*chi1 + d221*chi1 + d220;
    d23 = d233*chi1*chi1*chi1 + d232*chi1*chi1 + d231*chi1 + d230;
    d31 = d311*chi1 + d310;
    d32 = d321*chi1 + d320;

    d1 = d11*nu + d12*nu*nu + d13*nu*nu*nu;
    d2 = d21*nu + d22*nu*nu + d23*nu*nu*nu;
    d3 = d31*nu + d32*nu*nu;

    dt_bhns[3] = ( 1 + (d1*lambda*lambda + d2*lambda*lambda*lambda*lambda ) ) / ( (1 + lambda*lambda*d3*d3 )*(1 + lambda*lambda*d3*d3 ) ); ///Dt_lm_TP[3];
    
    for (int k=0; k<KMAX; k++) {
      Dt_lm[k] = Dt_lm[k]*dt_bhns[k];
    } 
}
