/** \file TEOBResumSWaveform.c
 *  \brief EOB IMR Waveform
 * 
 *  This file contains functions for the calculation of the EOB-factorized 
 *  waveform, including NQCs and ringdown attachment. 
 */

#include "TEOBResumS.h"

/** hlmNewt coefficients for amplitude */
static const double ChlmNewt_ampli[35] = {2.1137745587232057, 6.341323676169617, 0.1412325034218127, 1.7864655618418102, 4.9229202032627635, 0.023872650234580958, 0.2250735048768909, 1.7053495827316825, 4.763908164911493, 0.001122165903318321, 0.06333806741197714, 0.2945348827200268, 1.755276012972272, 5.0817902739730565, 0.00014954736544380072, 0.005296595280255638, 0.10342548105284892, 0.3713362832603404, 1.8983258440274462, 5.727111757630886, 5.184622059790144e-06, 0.0012191691413436815, 0.011783593824950922, 0.14639129995388936, 0.4653654097044924, 2.125638973894669, 6.685178460621457, 5.54485779151375621e-7, 0.0000763473331250837455, 0.00353250998285463003, 
0.0204988821800401766, 0.19579402814926015, 0.584571015778149663, 2.44207899966355693, 7.9955401278763745};

/* 
   {8./3.*sqrt(Pi/5.), 8.*sqrt(Pi/5.),
   1./3.*sqrt(2.*Pi/35.), 8./3.*sqrt(Pi/7.), 3.*sqrt(6.*Pi/7.), 
   1./105.*sqrt(2.*Pi), 8./63.*sqrt(Pi), 9./5*sqrt(2*Pi/7.), 64./9.*sqrt(Pi/7.),
   1./180.*sqrt(Pi/77.), 16./135.*sqrt(Pi/11.), 9./20.*sqrt(3*Pi/22.), 256./45.*sqrt(Pi/33.), 125./12.*sqrt(5.*Pi/66.), 
   1./2079.*sqrt(2.*Pi/65.), 16./1485.*sqrt(Pi/13.), 81./385.*sqrt(Pi/13.), 1024./495.*sqrt(2.*Pi/195.), 625./63*sqrt(5.*Pi/429.), 432./5*sqrt(Pi/715.), 
   1./108108.*sqrt(Pi/10.), 8./3003.*sqrt(Pi/15.), 243./20020.*sqrt(3.*Pi/10.), 1024./1365.*sqrt(2.*Pi/165.), 3125./3276.*sqrt(5.*Pi/66.), 648./35.*sqrt(3.*Pi/715.), 16807./180.*sqrt(7.*Pi/4290.), 
   131072./315.*sqrt(2.*Pi/17017.), 0., 0., 0., 0., 0., 0., 0.
   };
*/

/** hlmNewt additive coefficients for phase */
static const double ChlmNewt_phase[35] = {4.71238898038469, 3.141592653589793, 4.71238898038469, 0.0, 1.5707963267948966, 1.5707963267948966, 0.0, 4.71238898038469, 3.141592653589793, 1.5707963267948966, 3.141592653589793, 4.71238898038469, 0.0, 1.5707963267948966, 4.71238898038469, 3.141592653589793, 1.5707963267948966, 0.0, 4.71238898038469, 3.141592653589793, 4.71238898038469, 0.0, 1.5707963267948966, 3.141592653589793, 4.71238898038469, 0.0, 1.5707963267948966, 4.7123889803846898577, 0.0, 1.5707963267948966192, 3.1415926535897932385, 
   4.7123889803846898577, 0.0, 1.5707963267948966192, 3.1415926535897932385};

/*{
  3.*Pi/2., Pi, 
  3.*Pi/2., 0., Pi/2., 
  Pi/2., 0., 3.*Pi/2., Pi, 
  Pi/2., Pi, 3.*Pi/2., 0., Pi/2., 
  3.*Pi/2.,  Pi, Pi/2., 0., 3.*Pi/2., Pi, 
  3.*Pi/2., 0., Pi/2., Pi, 3.*Pi/2., 0., Pi/2., 
  Pi, 0, 0, 0, 0, 0, 0, 0
  };
*/

/* Negative amplitude modes */
int nNegAmp[35];	
	/* 	
	= {0, 0,	
	   0, 1, 0,	
	   1, 1, 1, 1,	
	   1, 0, 1, 0, 1,	
	   0, 0, 0, 0, 0, 0,	
	   0, 1, 0, 1, 0, 1, 0,	
	   1, 1, 1, 1, 1, 1, 1, 1};	
	*/

/**
 *  Function: eob_wav_hlmNewt_v1
 *  ----------------------------
 *   Leading-order (Newtonian) prefactor  of the multipolar resummed waveform. 
 *   Reference: Damour, Iyer & Nagar, PRD 79, 064004 (2009)
 *   
 *   @param[in] r: relative separation
 *   @param[in] Omega: orbital frequency
 *   @param[in] phi: orbital phase
 *   @param[in] nu: symmetric mass ratio
 *   @param[out] hlmNewt: Newtonia prefactor of multipolar waveform
 */
void eob_wav_hlmNewt_v1(double r,
			double Omega,
			double phi,
			double nu,
			Waveform_lm_t *hlmNewt)
{
  /* Shorthands */
  double nu2   = nu*nu;
  double nu3   = nu*nu2;
  
  double vphi  = r*Omega;
  double vphi2 = vphi*vphi;
  double vphi3 = vphi*vphi2;
  double vphi4 = vphi*vphi3;
  double vphi5 = vphi*vphi4;
  double vphi6 = vphi*vphi5;
  double vphi7 = vphi*vphi6;
  double vphi8 = vphi*vphi7;
  double vphi9 = vphi*vphi8;
  
  /* Polynomials in nu */
  const double p1 = 1.;
  const double p2 = sqrt(1.-4.*nu); 
  const double p3 = (3.*nu-1.);
  const double p4 = (2.*nu-1.)*sqrt(1.-4.*nu);
  const double p5 = 1.-5.*nu+5.*nu2;
  const double p6 = (1.-4.*nu+3.*nu2)*sqrt(1.-4.*nu);
  const double p7 = 7.*nu3 - 14.*nu2 + 7.*nu -1.;
  const double p8 = (4.*nu3 - 10.*nu2 + 6.*nu -1.)*sqrt(1.-4.*nu);
    
  const double phix2 = 2. * phi;
  const double phix3 = 3. * phi;
  const double phix4 = 4. * phi;
  const double phix5 = 5. * phi;
  const double phix6 = 6. * phi;
  const double phix7 = 7. * phi;
  
  const double pv23 = p2 * vphi3;
  const double pv34 = p3 * vphi4;
  const double pv45 = p4 * vphi5;
  const double pv56 = p5 * vphi6;
  const double pv67 = p6 * vphi7;
  const double pv78 = p7 * vphi8;
  const double pv89 = p8 * vphi9;
  
  double phim[35] = {
    phi,phix2,
    phi,phix2,phix3,
    phi,phix2,phix3,phix4,
    phi,phix2,phix3,phix4,phix5,
    phi,phix2,phix3,phix4,phix5,phix6,
    phi,phix2,phix3,phix4,phix5,phix6,phix7,
    phi,phix2,phix3,phix4,phix5,phix6,phix7,8.*phi
  };

  double Alm[35] = {
    pv23, p1 * vphi2, 
    pv23, pv34, pv23, 
    pv45, pv34, pv45, pv34, 
    pv45, pv56, pv45, pv56, pv45, 
    pv67, pv56, pv67, pv56, pv67, pv56, 
    pv67, pv78, pv67, pv78, pv67, pv78, pv67, 
    pv89, pv78, pv89, pv78, pv89, pv78, pv89, pv78
  };
    
  /* Compute hlmNewt (without phase factor) in complex Polar coords */
#pragma omp simd
  for (int k = 0; k < KMAX; k++) {
    hlmNewt->phase[k] = - phim[k] + ChlmNewt_phase[k];
    hlmNewt->ampli[k] = ChlmNewt_ampli[k] * Alm[k];
  }
  
}

/**
 *  Function: eob_wav_hlmNewt_HM
 *  ----------------------------
 *   Leading-order (Newtonian) prefactor  of the multipolar resummed waveform.
 *    New terms vphi and vOmg for higher modes. 
 *    Reference: arXiv:2001.09082
 *   
 *   @param[in]  r      : relative separation
 *   @param[in]  Omega  : orbital frequency
 *   @param[in]  phi    : orbital phase
 *   @param[in]  nu     : symmetric mass ratio
 *   @param[out] hlmNewt: Newtonia prefactor of multipolar waveform
 */
void eob_wav_hlmNewt_HM(double r,
			double Omega,
			double phi,
			double nu,
			Waveform_lm_t *hlmNewt)
{
  /** Shorthands */
  double nu2   = nu*nu;
  double nu3   = nu*nu2;

  double vOmg  = pow(fabs(Omega),1./3.);
  double vOmg2 = vOmg*vOmg;
  double vOmg3 = vOmg*vOmg2;
  double vOmg4 = vOmg*vOmg3;
  double vOmg5 = vOmg*vOmg4;

  double vphi  = r*Omega;
  double vphi2 = vphi*vphi;
  double vphi3 = vphi*vphi2;
  double vphi4 = vphi*vphi3;
  double vphi5 = vphi*vphi4;
  double vphi6 = vphi*vphi5;
  double vphi7 = vphi*vphi6;
  double vphi8 = vphi*vphi7;
  double vphi9 = vphi*vphi8;

  /** Polynomials in nu */
  const double p1 = 1.;
  const double p2 = 1.;//sqrt(1.-4.*nu); 
  const double p3 = (3.*nu-1.);
  const double p4 = (2.*nu-1.);//*sqrt(1.-4.*nu);
  const double p5 = 1.-5.*nu+5.*nu2;
  const double p6 = (1.-4.*nu+3.*nu2)*sqrt(1.-4.*nu);
  const double p7 = 7.*nu3 - 14.*nu2 + 7.*nu -1.;
  const double p8 = (4.*nu3 - 10.*nu2 + 6.*nu -1.)*sqrt(1.-4.*nu);
    
  const double phix2 = 2. * phi;
  const double phix3 = 3. * phi;
  const double phix4 = 4. * phi;
  const double phix5 = 5. * phi;
  const double phix6 = 6. * phi;
  const double phix7 = 7. * phi;
  
  const double pv45 = p4 * vphi5;
  const double pv56 = p5 * vphi6;
  const double pv67 = p6 * vphi7;
  const double pv78 = p7 * vphi8;
  const double pv89 = p8 * vphi9;

  const double fact22 = p1*vphi2;
  const double fact21 = p2*vphi3;
  const double fact33 = p2*vphi*vOmg2;
  const double fact32 = p3*vphi2*vOmg2;
  const double fact31 = p2*vOmg3;
  const double fact44 = p3*vphi2*vOmg2;
  const double fact43 = p4*vphi2*vOmg3;
  const double fact42 = p3*vphi*vOmg3;
  const double fact41 = p4*vOmg5;
  const double fact55 = p4*vphi*vOmg4;

  double phim[35] = {
    phi,phix2,
    phi,phix2,phix3,
    phi,phix2,phix3,phix4,
    phi,phix2,phix3,phix4,phix5,
    phi,phix2,phix3,phix4,phix5,phix6,
    phi,phix2,phix3,phix4,phix5,phix6,phix7,
    phi,phix2,phix3,phix4,phix5,phix6,phix7,8.*phi
  };

  double Alm[35] = {
    fact21, fact22, 
    fact31, fact32, fact33, 
    fact41, fact42, fact43, fact44, 
    pv45,   pv56,   pv45,   pv56,   fact55, 
    pv67,   pv56,   pv67,   pv56,   pv67,   pv56, 
    pv67,   pv78,   pv67,   pv78,   pv67,   pv78, pv67, 
    pv89,   pv78,   pv89,   pv78,   pv89,   pv78, pv89, pv78
  };
    
  /* Compute hlmNewt (without phase factor) in complex Polar coords */
#pragma omp simd
  for (int k = 0; k < KMAX; k++) {
    hlmNewt->phase[k] = - phim[k] + ChlmNewt_phase[k];
    hlmNewt->ampli[k] = ChlmNewt_ampli[k] * Alm[k];
  }

  /* Correcting phase in case of negative freuqency */
  if (Omega < 0.) {
    hlmNewt->phase[2]  += Pi; 
    hlmNewt->phase[3]  += 2./3.*Pi; 
    hlmNewt->phase[4]  += 2./3.*Pi; 
    hlmNewt->phase[5]  += 5./3.*Pi; 
    hlmNewt->phase[6]  += Pi; 
    hlmNewt->phase[7]  += Pi; 
    hlmNewt->phase[8]  += 2./3.*Pi; 
    hlmNewt->phase[13] += 4./3.*Pi; 
  }
  
}

/**
 *  Function: eob_wav_hlmNewt_ecc
 *  ----------------------------
 *   Leading-order (Newtonian) prefactor  of the multipolar resummed waveform.
 *   Eccentric version 
 *   See: https://arxiv.org/abs/2001.11736 
 *   
 *   @param[in]  dyn    : Dynamics structure
 *   @param[out] hlmNewt: Newtonia prefactor of multipolar waveform
 */
void eob_wav_hlmNewt_ecc(Dynamics *dyn, Waveform_lm_t *hlmNewt)
{
  /** Extracting variables */
  const double nu = EOBPars->nu;
  const double phi = dyn->phi; 

  const double rw     = dyn->r_omega;

  const double r     = dyn->r;
  const double rdot  = dyn->rdot;
  const double r2dot = dyn->r2dot;
  const double r3dot = dyn->r3dot;
  const double r4dot = dyn->r4dot;
  const double r5dot = dyn->r5dot;
  const double Omega     = dyn->Omg;
  const double Omegadot  = dyn->Omegadot;
  const double Omega2dot = dyn->Omega2dot;
  const double Omega3dot = dyn->Omega3dot;
  const double Omega4dot = dyn->Omega4dot;
  
  /* Shorthands */  
  double nu2   = nu*nu;
  double nu3   = nu*nu2;
  
  double r2 = r*r;
  double r3 = r2*r;
  double r4 = r3*r;
  double r5 = r4*r;
  double Omega2 = Omega*Omega;
  double Omega3 = Omega2*Omega;
  double Omega4 = Omega3*Omega;
  double rdot2 = rdot*rdot;
  double rdot3 = rdot2*rdot;
  double rdot4 = rdot3*rdot;
  double rdot5 = rdot4*rdot;
  double Omegadot2 = Omegadot*Omegadot;
  double r2dot2 = r2dot*r2dot;

  double vOmg  = pow(Omega,1./3.);
  double vOmg2 = vOmg*vOmg;
  double vOmg3 = vOmg*vOmg2;
  double vOmg4 = vOmg*vOmg3;
  double vOmg5 = vOmg*vOmg4;

  double vphi  = rw*Omega;
  double vphi2 = vphi*vphi;
  double vphi3 = vphi*vphi2;
  double vphi4 = vphi*vphi3;
  double vphi5 = vphi*vphi4;
  double vphi6 = vphi*vphi5;
  double vphi7 = vphi*vphi6;
  double vphi8 = vphi*vphi7;
  double vphi9 = vphi*vphi8;

  /* Eccentric variables */  
  // Non-factorized form
  /*
    double Re_vphi22_ecc = vphi2 + sigmoid*(-0.5*(rdot2 + r*r2dot));
    double Im_vphi22_ecc = sigmoid*(2.*r*Omega*rdot + 0.5*Omegadot*r2);
  */
  double Re_vphi22_ecc = vphi2*(1 + -0.5/r2/Omega2*(rdot2 + r*r2dot));
  double Im_vphi22_ecc = vphi2*(2.*rdot/r/Omega + 0.5*Omegadot/Omega2);  
  
  double Re_vphi21_ecc = vphi3 - 6.*r*Omega*rdot2 - 3.*r2*Omega*r2dot - 6.*r2*rdot*Omegadot - r3*Omega2dot;
  double Im_vphi21_ecc = 3.*r2*(2.*rdot*Omega2 + r*Omega*Omegadot);
  
  double Re_vphi33_ecc = vphi3 - r*Omega*(r*r2dot + 2.*rdot2) - 1./9*r3*Omega2dot - r2*rdot*Omegadot;
  double Im_vphi33_ecc = Omega*r3*Omegadot + 3.*r2*Omega2*rdot  - 2./9.*rdot2*rdot - 2./3.*r*rdot*r2dot - 1./9.*r2*r3dot;

  double Re_vphi32_ecc = vphi4 - 3.*r2*Omega2*(r*r2dot + 3.*rdot2) - r4*Omega*Omega2dot - 9.*Omega*r3*rdot*Omegadot - 0.75*r4*Omegadot2;
  double Im_vphi32_ecc = 6.*r3*Omega3*rdot -0.5*r3*Omega*r3dot - 4.5*r2*Omega*rdot*r2dot - 3.*r*Omega*rdot*rdot2 - 1.5*r3*rdot*Omega2dot
    - 0.125*r4*Omega3dot - 1.5*r3*r2dot*Omegadot - 4.5*r2*rdot2*Omegadot + 3.*r4*Omega2*Omegadot;
  
  double Re_vphi31_ecc = vphi3 - 9.*r2*Omega*r2dot - 18.*r*Omega*rdot2 - 9.*r2*rdot*Omegadot - r3*Omega2dot;
  double Im_vphi31_ecc = 9.*r2*Omega2*rdot - 3.*r2*r3dot - 18.*r*rdot*r2dot - 6.*rdot3 + 3.*Omega*r3*Omegadot;

  double Re_vphi44_ecc = vphi4 + 0.015625*r3*r4dot - 0.1875*r4*Omegadot2 + 0.1875*r2*rdot*r3dot + 0.09375*rdot4
    + 0.140625*r2*r2dot2 -4.5*r2*Omega2*rdot2 - 1.5*r3*Omega2*r2dot - 0.25*r4*Omega*Omega2dot - 3.*Omega*r3*rdot*Omegadot + 0.5625*r*rdot2*r2dot;
  double Im_vphi44_ecc = 4.*r3*Omega3*rdot - 1.5*r*Omega*rdot2*rdot - 0.015625*r4*Omega3dot - 0.25*r3*Omega*r3dot + 1.5*r4*Omega2*Omegadot
    - 0.25*r3*rdot*Omega2dot - 2.25*r2*Omega*rdot*r2dot - 0.375*r3*r2dot*Omegadot - 1.125*r2*rdot2*Omegadot;

  double Re_vphi43_ecc = vphi5 + 5./81.*r4*Omega*r4dot + 20./27.*r*Omega*r2*r2dot2 + 80./81.*r3*Omega*rdot*r3dot + 40./9.*r2*Omega*rdot2*r2dot
    + 40./27.*r*Omega*rdot4 - 10./3.*r4*Omega3*r2dot - 40./3.*r3*Omega3*rdot2 + 10./27.*r4*r2dot*Omega2dot + 40./27.*r3*rdot2*Omega2dot
    - 10./9.*r5*Omega2*Omega2dot + 20./81.*r4*rdot*Omega3dot + 1./81.*r5*Omega4dot + 20./81.*r4*r3dot*Omegadot + 80./27.*r3*rdot*r2dot*Omegadot
    + 80./27.*r2*rdot2*rdot*Omegadot - 40./3.*r4*Omega2*rdot*Omegadot - 5./3.*r*Omega*r4*Omegadot2;
  double Im_vphi43_ecc = 20./3.*r4*Omega4*rdot - 20./27.*r4*Omega2*r3dot - 80./9.*r3*Omega2*rdot*r2dot - 80./9.*r2*Omega2*rdot2*rdot
    - 80./27.*r4*Omega*rdot*Omega2dot - 5./27.*r5*Omega*Omega3dot - 10./3.*r4*Omega*Omegadot*r2dot - 40./3.*r3*Omega*rdot2*Omegadot
    + 10./3.*r5*Omega3*Omegadot - 10./27.*r5*Omegadot*Omega2dot - 20./9.*r4*Omegadot2*rdot;
  
  double Re_vphi42_ecc = vphi4 + 2.25*r2*r2dot2 + 0.25*r3*r4dot + 3.*r2*rdot*r3dot + 9.*r*rdot2*r2dot + 1.5*rdot4 - 6.*r3*Omega2*r2dot
    - 18.*r2*Omega2*rdot2 - r4*Omega*Omega2dot - 12*r3*Omega*rdot*Omegadot - 0.75*r4*Omegadot2;
  double Im_vphi42_ecc =  8.*r3*Omega3*rdot - 2.*r3*Omega*r3dot - 18.*r2*Omega*rdot*r2dot - 12.*r*Omega*rdot2*rdot
    - 2.*r3*rdot*Omega2dot - 0.125*r4*Omega3dot - 3.*r3*r2dot*Omegadot - 9.*r2*rdot2*Omegadot + 3.*r4*Omega2*Omegadot;
  
  double Re_vphi41_ecc = vphi5 + 60.*r3*Omega*r2dot2 + 5.*r4*Omega*r4dot + 80.*r3*Omega*rdot*r3dot + 360.*r2*Omega*rdot2*r2dot
    + 120.*r*Omega*rdot4 - 30.*r4*Omega3*r2dot - 120.*r3*Omega3*rdot2 + 30.*r4*r2dot*Omega2dot + 120.*r3*rdot2*Omega2dot
    - 10.*r5*Omega2*Omega2dot + 20.*r4*rdot*Omega3dot + r5*Omega4dot + 20.*r4*r3dot*Omegadot + 240.*r3*rdot*r2dot*Omegadot
    + 240.*r2*rdot3*Omegadot - 120.*r4*Omega2*rdot*Omegadot - 15.*r5*Omega*Omegadot2;
  double Im_vphi41_ecc = 20.*r4*Omega4*rdot - 20.*r4*Omega2*r3dot - 240.*r3*Omega2*rdot*r2dot - 240.*r2*Omega2*rdot3
    - 80.*r4*Omega*rdot*Omega2dot - 5.*r5*Omega*Omega3dot - 90.*r4*Omega*Omegadot*r2dot - 360.*r3*Omega*rdot2*Omegadot
    + 10.*r5*Omega3*Omegadot - 10.*r5*Omegadot*Omega2dot - 60.*r4*Omegadot2*rdot;
    
  double Re_vphi55_ecc = vphi5 + 0.48*r3*Omega*r2dot2 + 0.04*r4*Omega*r4dot + 0.64*r3*Omega*rdot*r3dot + 2.88*r2*Omega*rdot2*r2dot
    + 0.96*r*Omega*rdot4 - 2.*r4*Omega3*r2dot- 8.*r3*Omega3*rdot2 + 0.08*r4*r2dot*Omega2dot + 0.32*r3*rdot2*Omega2dot - 0.4*r5*Omega2*Omega2dot
    + 0.04*r4*rdot*Omega3dot + 0.0016*r5*Omega4dot + 0.08*r4*r3dot*Omegadot + 0.96*r3*rdot*r2dot*Omegadot + 0.96*r2*rdot3*Omegadot
    - 6.*r4*Omega2*rdot*Omegadot - 0.6*r5*Omega*Omegadot2;
  double Im_vphi55_ecc = 5.*r4*Omega4*rdot + 0.0016*r4*r5dot + 0.064*r3*r2dot*r3dot + 0.288*r2*r2dot2*rdot + 0.032*r3*rdot*r4dot
    + 0.192*r2*rdot2*r3dot + 0.384*r*rdot2*rdot*r2dot + 0.0384*rdot5 - 0.4*r4*Omega2*r3dot - 4.8*r3*Omega2*rdot*r2dot - 4.8*r2*Omega2*rdot2*rdot
    - 0.8*r4*Omega*rdot*Omega2dot - 0.04*r5*Omega*Omega3dot - 1.2*r4*Omega*r2dot*Omegadot - 4.8*r3*Omega*rdot2*Omegadot
    + 2.*r5*Omega3*Omegadot - 0.08*r5*Omegadot*Omega2dot - 0.6*r4*Omegadot2*rdot;

  double A22_ecc   = sqrt(SQ(Re_vphi22_ecc) + SQ(Im_vphi22_ecc));
  double A21_ecc   = sqrt(SQ(Re_vphi21_ecc) + SQ(Im_vphi21_ecc));
  double A33_ecc   = sqrt(SQ(Re_vphi33_ecc) + SQ(Im_vphi33_ecc));
  double A32_ecc   = sqrt(SQ(Re_vphi32_ecc) + SQ(Im_vphi32_ecc));
  double A31_ecc   = sqrt(SQ(Re_vphi31_ecc) + SQ(Im_vphi31_ecc));
  double A44_ecc   = sqrt(SQ(Re_vphi44_ecc) + SQ(Im_vphi44_ecc));
  double A43_ecc   = sqrt(SQ(Re_vphi43_ecc) + SQ(Im_vphi43_ecc));
  double A42_ecc   = sqrt(SQ(Re_vphi42_ecc) + SQ(Im_vphi42_ecc));
  double A41_ecc   = sqrt(SQ(Re_vphi41_ecc) + SQ(Im_vphi41_ecc));
  double A55_ecc   = sqrt(SQ(Re_vphi55_ecc) + SQ(Im_vphi55_ecc));

  double phi22_ecc = atan2(Im_vphi22_ecc,Re_vphi22_ecc);
  double phi21_ecc = atan2(Im_vphi21_ecc,Re_vphi21_ecc);
  double phi33_ecc = atan2(Im_vphi33_ecc,Re_vphi33_ecc);
  double phi32_ecc = atan2(Im_vphi32_ecc,Re_vphi32_ecc);
  double phi31_ecc = atan2(Im_vphi31_ecc,Re_vphi31_ecc);
  double phi44_ecc = atan2(Im_vphi44_ecc,Re_vphi44_ecc);
  double phi43_ecc = atan2(Im_vphi43_ecc,Re_vphi43_ecc);
  double phi42_ecc = atan2(Im_vphi42_ecc,Re_vphi42_ecc);
  double phi41_ecc = atan2(Im_vphi41_ecc,Re_vphi41_ecc);
  double phi55_ecc = atan2(Im_vphi55_ecc,Re_vphi55_ecc);
  
  
  /* Polynomials in nu */
  const double p1 = 1.;
  const double p2 = 1.;//sqrt(1.-4.*nu); 
  const double p3 = (3.*nu-1.);
  const double p4 = (2.*nu-1.);//*sqrt(1.-4.*nu);
  const double p5 = 1.-5.*nu+5.*nu2;
  const double p6 = (1.-4.*nu+3.*nu2)*sqrt(1.-4.*nu);
  const double p7 = 7.*nu3 - 14.*nu2 + 7.*nu -1.;
  const double p8 = (4.*nu3 - 10.*nu2 + 6.*nu -1.)*sqrt(1.-4.*nu);
    
  const double phix2 = 2. * phi;
  const double phix3 = 3. * phi;
  const double phix4 = 4. * phi;
  const double phix5 = 5. * phi;
  const double phix6 = 6. * phi;
  const double phix7 = 7. * phi;
  
  const double pv45 = p4 * vphi5;
  const double pv56 = p5 * vphi6;
  const double pv67 = p6 * vphi7;
  const double pv78 = p7 * vphi8;
  const double pv89 = p8 * vphi9;

  const double fact22 = p1*A22_ecc;
  const double fact21 = p2*A21_ecc;
  const double fact33 = p2*A33_ecc;
  const double fact32 = p3*A32_ecc;
  const double fact31 = p2*A31_ecc;
  const double fact44 = p3*A44_ecc;
  const double fact43 = p4*A43_ecc;
  const double fact42 = p3*A42_ecc;
  const double fact41 = p4*A41_ecc;
  const double fact55 = p4*A55_ecc;

  double phiecc[35] = {
    phi21_ecc, phi22_ecc,
    phi31_ecc, phi32_ecc, phi33_ecc,
    phi41_ecc, phi42_ecc, phi43_ecc, phi44_ecc,
    0.,        0.,        0.,        0.,        phi55_ecc,
    0.,        0.,        0.,        0.,        0.,        0.,
    0.,        0.,        0.,        0.,        0.,        0., 0., 
    0.,        0.,        0.,        0.,        0.,        0., 0., 0. 
  };

  double phim[35] = {
    phi, phix2,
    phi, phix2, phix3,
    phi, phix2, phix3, phix4,
    phi, phix2, phix3, phix4, phix5,
    phi, phix2, phix3, phix4, phix5, phix6,
    phi, phix2, phix3, phix4, phix5, phix6, phix7,
    phi, phix2, phix3, phix4, phix5, phix6, phix7, 8.*phi
  };

  double Alm[35] = {
    fact21, fact22, 
    fact31, fact32, fact33, 
    fact41, fact42, fact43, fact44, 
    pv45,   pv56,   pv45,   pv56,   fact55, 
    pv67,   pv56,   pv67,   pv56,   pv67,   pv56, 
    pv67,   pv78,   pv67,   pv78,   pv67,   pv78, pv67, 
    pv89,   pv78,   pv89,   pv78,   pv89,   pv78, pv89, pv78
  };
    
  /* Compute hlmNewt (without phase factor) in complex Polar coords */
#pragma omp simd
  for (int k = 0; k < KMAX; k++) {
    hlmNewt->phase[k] = - phim[k] + phiecc[k] + ChlmNewt_phase[k];
    hlmNewt->ampli[k] = ChlmNewt_ampli[k] * Alm[k];
  }
}

/**
 *  Function: eob_wav_hlmNewt_ecc
 *  ----------------------------
 *   Leading-order (Newtonian) prefactor  of the multipolar resummed waveform.
 *    Eccentric version with sigmoid function
 *   See: https://arxiv.org/abs/2001.11736, 
 *        https://arxiv.org/abs/2101.08624
 *   
 *   @param[in]  dyn    : Dynamics structure
 *   @param[out] hlmNewt: Newtonia prefactor of multipolar waveform
 */
void eob_wav_hlmNewt_ecc_sigmoid(Dynamics *dyn, Waveform_lm_t *hlmNewt)
{
  /** Extracting variables */
  const double nu  = EOBPars->nu;
  const double tOmg_pk = dyn->tOmg_pk;

  const double t   = dyn->t;
  const double phi = dyn->phi;
  const double rw  = dyn->r_omega;

  const double r     = dyn->r;
  const double rdot  = dyn->rdot;
  const double r2dot = dyn->r2dot;
  const double r3dot = dyn->r3dot;
  const double r4dot = dyn->r4dot;
  const double r5dot = dyn->r5dot;
  const double Omega     = dyn->Omg;
  const double Omegadot  = dyn->Omegadot;
  const double Omega2dot = dyn->Omega2dot;
  const double Omega3dot = dyn->Omega3dot;
  const double Omega4dot = dyn->Omega4dot;
  
  /* Shorthands */
  double nu2   = nu*nu;
  double nu3   = nu*nu2;
  
  double r2 = r*r;
  double r3 = r2*r;
  double r4 = r3*r;
  double r5 = r4*r;
  double Omega2 = Omega*Omega;
  double Omega3 = Omega2*Omega;
  double Omega4 = Omega3*Omega;
  double rdot2 = rdot*rdot;
  double rdot3 = rdot2*rdot;
  double rdot4 = rdot3*rdot;
  double rdot5 = rdot4*rdot;
  double Omegadot2 = Omegadot*Omegadot;
  double r2dot2 = r2dot*r2dot;

  double vOmg  = pow(Omega,1./3.);
  double vOmg2 = vOmg*vOmg;
  double vOmg3 = vOmg*vOmg2;
  double vOmg4 = vOmg*vOmg3;
  double vOmg5 = vOmg*vOmg4;

  double vphi  = rw*Omega;
  double vphi2 = vphi*vphi;
  double vphi3 = vphi*vphi2;
  double vphi4 = vphi*vphi3;
  double vphi5 = vphi*vphi4;
  double vphi6 = vphi*vphi5;
  double vphi7 = vphi*vphi6;
  double vphi8 = vphi*vphi7;
  double vphi9 = vphi*vphi8;

  /* Defining sigmoid function to switch off Newtonian factor near the end of the evolution */
  double t0      = tOmg_pk - EOBPars->delta_t0_sigmoid_Newt;
  double alpha   = EOBPars->alpha_sigmoid_Newt;
  double sigmoid = 1./(1. + exp(alpha*(t - t0)));
    
  /* Eccentric variables */  
  // Non-factorized form
  /*
    double Re_vphi22_ecc = vphi2 + sigmoid*(-0.5*(rdot2 + r*r2dot));
    double Im_vphi22_ecc = sigmoid*(2.*r*Omega*rdot + 0.5*Omegadot*r2);
  */
  double Re_vphi22_ecc = vphi2*(1 + sigmoid*(-0.5/r2/Omega2*(rdot2 + r*r2dot)));
  double Im_vphi22_ecc = vphi2*sigmoid*(2.*rdot/r/Omega + 0.5*Omegadot/Omega2);  
  
  double Re_vphi21_ecc = vphi3 + sigmoid*(- 6.*r*Omega*rdot2 - 3.*r2*Omega*r2dot - 6.*r2*rdot*Omegadot - r3*Omega2dot);
  double Im_vphi21_ecc = sigmoid*(3.*r2*(2.*rdot*Omega2 + r*Omega*Omegadot));
  
  double Re_vphi33_ecc = vphi3 + sigmoid*(- r*Omega*(r*r2dot + 2.*rdot2) - 1./9*r3*Omega2dot - r2*rdot*Omegadot);
  double Im_vphi33_ecc = sigmoid*(Omega*r3*Omegadot + 3.*r2*Omega2*rdot  - 2./9.*rdot2*rdot - 2./3.*r*rdot*r2dot - 1./9.*r2*r3dot);

  double Re_vphi32_ecc = vphi4 + sigmoid*(- 3.*r2*Omega2*(r*r2dot + 3.*rdot2) - r4*Omega*Omega2dot
					  - 9.*Omega*r3*rdot*Omegadot - 0.75*r4*Omegadot2);
  double Im_vphi32_ecc = sigmoid*(6.*r3*Omega3*rdot -0.5*r3*Omega*r3dot - 4.5*r2*Omega*rdot*r2dot
				  - 3.*r*Omega*rdot*rdot2 - 1.5*r3*rdot*Omega2dot - 0.125*r4*Omega3dot
				  - 1.5*r3*r2dot*Omegadot - 4.5*r2*rdot2*Omegadot + 3.*r4*Omega2*Omegadot);
				  
  double Re_vphi31_ecc = vphi3 + sigmoid*(- 9.*r2*Omega*r2dot - 18.*r*Omega*rdot2 - 9.*r2*rdot*Omegadot - r3*Omega2dot);
  double Im_vphi31_ecc = sigmoid*(9.*r2*Omega2*rdot - 3.*r2*r3dot - 18.*r*rdot*r2dot - 6.*rdot3 + 3.*Omega*r3*Omegadot);

  double Re_vphi44_ecc = vphi4 + sigmoid*(0.015625*r3*r4dot - 0.1875*r4*Omegadot2 + 0.1875*r2*rdot*r3dot + 0.09375*rdot4
					  + 0.140625*r2*r2dot2 -4.5*r2*Omega2*rdot2 - 1.5*r3*Omega2*r2dot
					  - 0.25*r4*Omega*Omega2dot - 3.*Omega*r3*rdot*Omegadot + 0.5625*r*rdot2*r2dot);
  double Im_vphi44_ecc = sigmoid*(4.*r3*Omega3*rdot - 1.5*r*Omega*rdot2*rdot - 0.015625*r4*Omega3dot
				  - 0.25*r3*Omega*r3dot + 1.5*r4*Omega2*Omegadot - 0.25*r3*rdot*Omega2dot
				  - 2.25*r2*Omega*rdot*r2dot - 0.375*r3*r2dot*Omegadot - 1.125*r2*rdot2*Omegadot);

  double Re_vphi43_ecc = vphi5 + sigmoid*(5./81.*r4*Omega*r4dot + 20./27.*r*Omega*r2*r2dot2 + 80./81.*r3*Omega*rdot*r3dot
					  + 40./9.*r2*Omega*rdot2*r2dot + 40./27.*r*Omega*rdot4 - 10./3.*r4*Omega3*r2dot
					  - 40./3.*r3*Omega3*rdot2 + 10./27.*r4*r2dot*Omega2dot + 40./27.*r3*rdot2*Omega2dot
					  - 10./9.*r5*Omega2*Omega2dot + 20./81.*r4*rdot*Omega3dot + 1./81.*r5*Omega4dot
					  + 20./81.*r4*r3dot*Omegadot + 80./27.*r3*rdot*r2dot*Omegadot
					  + 80./27.*r2*rdot2*rdot*Omegadot - 40./3.*r4*Omega2*rdot*Omegadot
					  - 5./3.*r*Omega*r4*Omegadot2);
  double Im_vphi43_ecc = sigmoid*(20./3.*r4*Omega4*rdot - 20./27.*r4*Omega2*r3dot - 80./9.*r3*Omega2*rdot*r2dot
				  - 80./9.*r2*Omega2*rdot2*rdot - 80./27.*r4*Omega*rdot*Omega2dot - 5./27.*r5*Omega*Omega3dot
				  - 10./3.*r4*Omega*Omegadot*r2dot - 40./3.*r3*Omega*rdot2*Omegadot
				  + 10./3.*r5*Omega3*Omegadot - 10./27.*r5*Omegadot*Omega2dot - 20./9.*r4*Omegadot2*rdot);

  double Re_vphi42_ecc = vphi4 + sigmoid*(2.25*r2*r2dot2 + 0.25*r3*r4dot + 3.*r2*rdot*r3dot + 9.*r*rdot2*r2dot + 1.5*rdot4
					  - 6.*r3*Omega2*r2dot - 18.*r2*Omega2*rdot2 - r4*Omega*Omega2dot
					  - 12*r3*Omega*rdot*Omegadot - 0.75*r4*Omegadot2);
  double Im_vphi42_ecc =  sigmoid*(8.*r3*Omega3*rdot - 2.*r3*Omega*r3dot - 18.*r2*Omega*rdot*r2dot - 12.*r*Omega*rdot2*rdot
				   - 2.*r3*rdot*Omega2dot - 0.125*r4*Omega3dot - 3.*r3*r2dot*Omegadot
				   - 9.*r2*rdot2*Omegadot + 3.*r4*Omega2*Omegadot);

  double Re_vphi41_ecc = vphi5 + sigmoid*(60.*r3*Omega*r2dot2 + 5.*r4*Omega*r4dot + 80.*r3*Omega*rdot*r3dot
					  + 360.*r2*Omega*rdot2*r2dot + 120.*r*Omega*rdot4 - 30.*r4*Omega3*r2dot
					  - 120.*r3*Omega3*rdot2 + 30.*r4*r2dot*Omega2dot
					  + 120.*r3*rdot2*Omega2dot - 10.*r5*Omega2*Omega2dot
					  + 20.*r4*rdot*Omega3dot + r5*Omega4dot + 20.*r4*r3dot*Omegadot
					  + 240.*r3*rdot*r2dot*Omegadot + 240.*r2*rdot3*Omegadot
					  - 120.*r4*Omega2*rdot*Omegadot - 15.*r5*Omega*Omegadot2);
  double Im_vphi41_ecc = sigmoid*(20.*r4*Omega4*rdot - 20.*r4*Omega2*r3dot - 240.*r3*Omega2*rdot*r2dot
				  - 240.*r2*Omega2*rdot3 - 80.*r4*Omega*rdot*Omega2dot - 5.*r5*Omega*Omega3dot
				  - 90.*r4*Omega*Omegadot*r2dot - 360.*r3*Omega*rdot2*Omegadot
				  + 10.*r5*Omega3*Omegadot - 10.*r5*Omegadot*Omega2dot - 60.*r4*Omegadot2*rdot);

  double Re_vphi55_ecc = vphi5 + sigmoid*(0.48*r3*Omega*r2dot2 + 0.04*r4*Omega*r4dot + 0.64*r3*Omega*rdot*r3dot + 2.88*r2*Omega*rdot2*r2dot
					  + 0.96*r*Omega*rdot4 - 2.*r4*Omega3*r2dot- 8.*r3*Omega3*rdot2 + 0.08*r4*r2dot*Omega2dot
					  + 0.32*r3*rdot2*Omega2dot - 0.4*r5*Omega2*Omega2dot + 0.04*r4*rdot*Omega3dot + 0.0016*r5*Omega4dot
					  + 0.08*r4*r3dot*Omegadot + 0.96*r3*rdot*r2dot*Omegadot + 0.96*r2*rdot3*Omegadot
					  - 6.*r4*Omega2*rdot*Omegadot - 0.6*r5*Omega*Omegadot2);
  double Im_vphi55_ecc = sigmoid*(5.*r4*Omega4*rdot + 0.0016*r4*r5dot + 0.064*r3*r2dot*r3dot + 0.288*r2*r2dot2*rdot + 0.032*r3*rdot*r4dot
				  + 0.192*r2*rdot2*r3dot + 0.384*r*rdot2*rdot*r2dot + 0.0384*rdot5 - 0.4*r4*Omega2*r3dot
				  - 4.8*r3*Omega2*rdot*r2dot - 4.8*r2*Omega2*rdot2*rdot - 0.8*r4*Omega*rdot*Omega2dot
				  - 0.04*r5*Omega*Omega3dot - 1.2*r4*Omega*r2dot*Omegadot - 4.8*r3*Omega*rdot2*Omegadot
				  + 2.*r5*Omega3*Omegadot - 0.08*r5*Omegadot*Omega2dot - 0.6*r4*Omegadot2*rdot);

  double A22_ecc   = sqrt(SQ(Re_vphi22_ecc) + SQ(Im_vphi22_ecc));
  double A21_ecc   = sqrt(SQ(Re_vphi21_ecc) + SQ(Im_vphi21_ecc));
  double A33_ecc   = sqrt(SQ(Re_vphi33_ecc) + SQ(Im_vphi33_ecc));
  double A32_ecc   = sqrt(SQ(Re_vphi32_ecc) + SQ(Im_vphi32_ecc));
  double A31_ecc   = sqrt(SQ(Re_vphi31_ecc) + SQ(Im_vphi31_ecc));
  double A44_ecc   = sqrt(SQ(Re_vphi44_ecc) + SQ(Im_vphi44_ecc));
  double A43_ecc   = sqrt(SQ(Re_vphi43_ecc) + SQ(Im_vphi43_ecc));
  double A42_ecc   = sqrt(SQ(Re_vphi42_ecc) + SQ(Im_vphi42_ecc));
  double A41_ecc   = sqrt(SQ(Re_vphi41_ecc) + SQ(Im_vphi41_ecc));
  double A55_ecc   = sqrt(SQ(Re_vphi55_ecc) + SQ(Im_vphi55_ecc));

  double phi22_ecc = atan2(Im_vphi22_ecc,Re_vphi22_ecc);
  double phi21_ecc = atan2(Im_vphi21_ecc,Re_vphi21_ecc);
  double phi33_ecc = atan2(Im_vphi33_ecc,Re_vphi33_ecc);
  double phi32_ecc = atan2(Im_vphi32_ecc,Re_vphi32_ecc);
  double phi31_ecc = atan2(Im_vphi31_ecc,Re_vphi31_ecc);
  double phi44_ecc = atan2(Im_vphi44_ecc,Re_vphi44_ecc);
  double phi43_ecc = atan2(Im_vphi43_ecc,Re_vphi43_ecc);
  double phi42_ecc = atan2(Im_vphi42_ecc,Re_vphi42_ecc);
  double phi41_ecc = atan2(Im_vphi41_ecc,Re_vphi41_ecc);
  double phi55_ecc = atan2(Im_vphi55_ecc,Re_vphi55_ecc);
  
  /** Polynomials in nu */
  const double p1 = 1.;
  const double p2 = 1.;//sqrt(1.-4.*nu); 
  const double p3 = (3.*nu-1.);
  const double p4 = (2.*nu-1.);//*sqrt(1.-4.*nu);
  const double p5 = 1.-5.*nu+5.*nu2;
  const double p6 = (1.-4.*nu+3.*nu2)*sqrt(1.-4.*nu);
  const double p7 = 7.*nu3 - 14.*nu2 + 7.*nu -1.;
  const double p8 = (4.*nu3 - 10.*nu2 + 6.*nu -1.)*sqrt(1.-4.*nu);
    
  const double phix2 = 2. * phi;
  const double phix3 = 3. * phi;
  const double phix4 = 4. * phi;
  const double phix5 = 5. * phi;
  const double phix6 = 6. * phi;
  const double phix7 = 7. * phi;
  
  const double pv45 = p4 * vphi5;
  const double pv56 = p5 * vphi6;
  const double pv67 = p6 * vphi7;
  const double pv78 = p7 * vphi8;
  const double pv89 = p8 * vphi9;

  const double fact22 = p1*A22_ecc;
  const double fact21 = p2*A21_ecc;
  const double fact33 = p2*A33_ecc;
  const double fact32 = p3*A32_ecc;
  const double fact31 = p2*A31_ecc;
  const double fact44 = p3*A44_ecc;
  const double fact43 = p4*A43_ecc;
  const double fact42 = p3*A42_ecc;
  const double fact41 = p4*A41_ecc;
  const double fact55 = p4*A55_ecc;

  double phiecc[35] = {
    phi21_ecc, phi22_ecc,
    phi31_ecc, phi32_ecc, phi33_ecc,
    phi41_ecc, phi42_ecc, phi43_ecc, phi44_ecc,
    0.,        0.,        0.,        0.,        phi55_ecc,
    0.,        0.,        0.,        0.,        0.,        0.,
    0.,        0.,        0.,        0.,        0.,        0., 0., 
    0.,        0.,        0.,        0.,        0.,        0., 0., 0. 
  };

  double phim[35] = {
    phi, phix2,
    phi, phix2, phix3,
    phi, phix2, phix3, phix4,
    phi, phix2, phix3, phix4, phix5,
    phi, phix2, phix3, phix4, phix5, phix6,
    phi, phix2, phix3, phix4, phix5, phix6, phix7,
    phi, phix2, phix3, phix4, phix5, phix6, phix7, 8.*phi
  };

  double Alm[35] = {
    fact21, fact22, 
    fact31, fact32, fact33, 
    fact41, fact42, fact43, fact44, 
    pv45,   pv56,   pv45,   pv56,   fact55, 
    pv67,   pv56,   pv67,   pv56,   pv67,   pv56, 
    pv67,   pv78,   pv67,   pv78,   pv67,   pv78, pv67, 
    pv89,   pv78,   pv89,   pv78,   pv89,   pv78, pv89, pv78
  };
    
  /* Compute hlmNewt (without phase factor) in complex Polar coords */
#pragma omp simd
  for (int k = 0; k < KMAX; k++) {
    hlmNewt->phase[k] = - phim[k] + phiecc[k] + ChlmNewt_phase[k];
    hlmNewt->ampli[k] = ChlmNewt_ampli[k] * Alm[k];
  }
}

/**
 * Function: eob_wav_hhatlmTail 
 * ----------------------------
 *   Tail contribution to the resummed wave.   
 *   Ref. Damour, Iyer & Nagar, PRD 79, 064004 (2009)
 * 
 *   @param[in]  Omega: orbital frequency
 *   @param[in]  Hreal: 
 *   @param[in]  bphys:
 *   @param[out] tlm  : tail contribution to the resummed wave
*/
void eob_wav_hhatlmTail(double Omega, double Hreal, double bphys, Waveform_lm_t *tlm)
{
  double k;
  double hhatk;
  
  gsl_sf_result num_rad;
  gsl_sf_result num_phase;
  gsl_sf_result denom_rad;
  gsl_sf_result denom_phase;
  
  double ratio_rad;
  double ratio_ang;
  double tlm_rad;
  double tlm_phase;
  
  for (int i = 0; i < KMAX; i++) {
    k     = MINDEX[i] * Omega;
    hhatk = k * Hreal;
    
    gsl_sf_lngamma_complex_e(LINDEX[i] + 1., -2.*hhatk, &num_rad, &num_phase);
    gsl_sf_lngamma_complex_e(LINDEX[i] + 1., 0., &denom_rad, &denom_phase);
    
    ratio_rad     = num_rad.val-denom_rad.val;
    ratio_ang     = num_phase.val-0.;
    
    if (Omega > 0.) {
      tlm_rad       = ratio_rad + Pi * hhatk;
      tlm_phase     = ratio_ang + 2.*hhatk*log(2.*k*bphys);
    } else {
      tlm_rad       = ratio_rad - Pi * hhatk;
      tlm_phase     = ratio_ang + 2.*hhatk*log(-2.*k*bphys);
    }
    
    tlm->ampli[i] = exp(tlm_rad);
    tlm->phase[i] = tlm_phase;
  }

}

/**
 * Function: eob_wav_speedyTail 
 * ----------------------------
 *   Alternative implementation of the phase of the tail factor   
 *   Ref. Damour, Iyer & Nagar, PRD 79, 064004 (2009); Nagar et al (2018)
 * 
 *   @param[in]  Omega: orbital frequency
 *   @param[in]  Hreal: 
 *   @param[in]  bphys:
 *   @param[out] tlm  : tail contribution to the resummed wave
*/
void eob_wav_speedyTail(double Omega, double Hreal, double bphys, Waveform_lm_t *tlm)
{
  double x;
  double x2;
  double x3;
  double x4;
  double x5;
  double tlm_ang;
  double num_ang;
  
  /** Fit coefficients*/
  const double b1[] = {
    0.1113090643348557, 0.1112593821157397, 0.0424759238428813, 0.0424489015884926, 0.0424717446800903, 0.0215953972500844, 0.0215873812155663, 0.0215776183122621, 0.0216017621863542, 0.0128123696874894, 0.0128097056242375, 0.0128038943888768, 0.0128025242617949, 0.0128202485907368, 0.0083762045692408, 0.0083751913886140, 0.0083724067460769, 0.0083694435961860, 0.0083710364141552, 0.0083834483913443, 0.0058540393221396, 0.0058536069384738, 0.0058522594457692, 0.0058502436535615, 0.0058491157293566, 0.0058514875071582, 0.0058602498033381, 0.0042956812356573, 0.0042954784390887, 0.0042947951664056, 0.0042935886137697, 0.0042923691461384, 0.0042922256848799, 0.0042945927126022, 0.0043009106861259};
  
  const double b2[] = {
    0.0004643273300862, 0.0009375605440004, 0.0000597134489198, 0.0002551406918111, 0.0001741036904709, 0.0000124649041611, 0.0000685496215625, 0.0001131160409390, 0.0000419907542591, 0.0000035218982282, 0.0000219211271097, 0.0000473186962874, 0.0000524142634057, 0.0000106823372552, 0.0000012237574387, 0.0000081742188269, 0.0000201940563214, 0.0000295722761753, 0.0000260539631956, 0.0000018994753518, 0.0000004932942990, 0.0000034477210351, 0.0000092294406360, 0.0000155143073237, 0.0000183386499818, 0.0000137922469695, -0.0000007075155453, 0.0000002223410995, 0.0000016045317657, 0.0000045260028113, 0.0000082655700107, 0.0000112393599417, 0.0000115758243113, 0.0000076838709956, -0.0000014020591745};
  
  const double b3[] = {
    -0.0221835462237291, -0.0235386333304348, -0.0042911639711832, -0.0047431560217121, -0.0046577314472149, -0.0013089557502947, -0.0014343968205390, -0.0014978542575474, -0.0014329302934532, -0.0005167994164556, -0.0005573939123058, -0.0005921030407223, -0.0005978284714483, -0.0005673965369076, -0.0002409269302708, -0.0002561516055118, -0.0002723768586352, -0.0002815958312453, -0.0002792078156272, -0.0002646630240693, -0.0001261183503407, -0.0001325622938779, -0.0001403198638518, -0.0001464084186977, -0.0001485971591029, -0.0001459023931717, -0.0001384829633836, -0.0000719062974278, -0.0000749128468013, -0.0000788187384314, -0.0000824202283094, -0.0000846673495936, -0.0000849054394951, -0.0000829269749240, -0.0000788883333858};
  
  const double b4[] = {
    0.0058366730167965, 0.0070452306758401, 0.0006914295465364, 0.0010322294603561, 0.0010057563135650, 0.0001394203795507, 0.0002309706405978, 0.0002596611624417, 0.0002409588083156, 0.0000386949167221, 0.0000679154947896, 0.0000830199015202, 0.0000850120755064, 0.0000780125513602, 0.0000133034384660, 0.0000241813441339, 0.0000311573885555, 0.0000340233089866, 0.0000335167900637, 0.0000307571022927, 0.0000053305073331, 0.0000099143129290, 0.0000132296989826, 0.0000150959309402, 0.0000156304390748, 0.0000151274875147, 0.0000139320508803, 0.0000023959090314, 0.0000045285807761, 0.0000061918979830, 0.0000072894226381, 0.0000078251853305, 0.0000078772667984, 0.0000075606242809, 0.0000069956215270
  };
  
  double Tlm_real[KMAX];
  eob_flx_Tlm(Omega*Hreal, Tlm_real);
  
  /** Pre-computed psi */
  const double psi[] = {0.9227843350984671394, 0.9227843350984671394,
			1.256117668431800473, 1.256117668431800473, 1.256117668431800473,
			1.506117668431800473, 1.506117668431800473, 1.506117668431800473, 1.506117668431800473,
			1.706117668431800473, 1.706117668431800473, 1.706117668431800473, 1.706117668431800473, 1.706117668431800473,
			1.872784335098467139, 1.872784335098467139, 1.872784335098467139, 1.872784335098467139, 1.872784335098467139, 1.872784335098467139,
			2.015641477955609997, 2.015641477955609997, 2.015641477955609997, 2.015641477955609997, 2.015641477955609997, 2.015641477955609997, 2.015641477955609997,
			2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997};
  
  double k;
  for (int i=0; i<KMAX; i++) {
    k  = MINDEX[i] * Omega;
    x  = k * Hreal; /* hathatk */
    x2 = x * x;
    x3 = x2 * x;
    x4 = x3 * x;
    x5 = x4 * x;      
    num_ang   = 1. + b1[i]*x2 + b2[i]*x3 + b3[i]*x4 + b4[i]*x5; 
    tlm_ang   = (- 2. * psi[i] * x * num_ang) + 2.*x* log(2. * k * bphys);
    tlm->ampli[i] = Tlm_real[i];
    tlm->phase[i] = tlm_ang;
  }
  
}

/**
 * Function: eob_wav_deltalm_v1
 * ----------------------------
 *  Residual phase corrections delta_{lm} up to l=m=5.
 * Reference(s)
 * Damour, Iyer & Nagar, PRD 79, 064004 (2008)
 * Fujita & Iyer, PRD 82 044051 (2010)
 * Faye et al., Class. Q. Grav. 29 175004 (2012)
 * Damour, Nagar & Bernuzzi, PRD 87, 084035 (2013)
 * @note: this routine can be optimized: precompute coefficients c(nu)
 * 
 * @param[in]  Hreal : 
 * @param[in]  Omega : orbital frequency
 * @param[in]  nu    : symmetric mass ratio
 * @param[out] dlm   : residual phase corrections
*/
void eob_wav_deltalm_v1(double Hreal,double Omega,double nu, double *dlm)
{
    
  /* Useful shorthands*/
  const double Pi2 = SQ(Pi);
  double nu2    = SQ(nu);
  double y      = cbrt(Hreal*Omega*Hreal*Omega);
  double sqrt_y = sqrt(y);
  double y3     = y*y*y;
  double y32    = Hreal*Omega;
  
  /* Leading order contributions*/
  double delta22LO = 7./3.   * y32;
  double delta21LO = 2./3.   * y32;
  double delta33LO = 13./10. * y32;
  double delta31LO = 13./30. * y32;
  
  /* Init phase */
#pragma omp simd
  for (int k = 0; k < KMAX; k++) {
    dlm[k] = 0.;
  }
  
  /* Residual phases in Pade-resummed form when possible */
  double num;
  double den;

  /* l=2 */
  /* Pade(1,2) approximant */
  num        = 69020.*nu + 5992.*Pi*sqrt_y;
  den        = 5992.*Pi*sqrt_y + 2456.*nu*(28.+493.*nu* y);
  dlm[0] = delta21LO*num/den;
  /* Pade(2,2) approximant */
  num        = (808920.*nu*Pi*sqrt(y) + 137388.*Pi2*y + 35.*nu2*(136080. + (154975. - 1359276.*nu)*y));
  den        = (808920.*nu*Pi*sqrt(y) + 137388.*Pi2*y + 35.*nu2*(136080. + (154975. + 40404.*nu)*y));
  dlm[1] = delta22LO*num/den;
  
  /* l=3 */
  /* Pade(1,2) approximant */
  num        = 4641.*nu + 1690.*Pi*sqrt_y;
  den        = num + 18207.*nu2*y;
  dlm[2] = delta31LO*num/den;
  /* Taylor-expanded form */
  num        = 1.  + 94770.*Pi/(566279.*nu)*sqrt_y;
  den        = num + 80897.* nu/3159.*y;
  dlm[3] = (10.+33.*nu)/(15.*(1.-3.*nu)) * y32 + 52./21.*Pi*y3;
  /* Pade(1,2) approximant */
  dlm[4] = delta33LO*num/den;
  
  /* l=4 */
  dlm[5] =   (2.+507.*nu)/(10.*(1.-2.*nu))*y32   + 1571./3465.*Pi*y3;
  dlm[6] =  7.*(1.+6.*nu)/(15.*(1.-3.*nu))*y32   + 6284./3465.*Pi*y3;
  dlm[7] = (486.+4961.*nu)/(810.*(1.-2.*nu))*y32 + 1571./385.*Pi*y3;
  dlm[8] =  (112.+219.*nu)/(120.*(1.-3.*nu))*y32 + 25136./3465.*Pi*y3;
  
  /* l=5 */
  dlm[13] = (96875. + 857528.*nu)/(131250.*(1.-2.*nu))*y32;
  
}

/**
 * Function: eob_wav_deltalm_HM
 * ----------------------------
 *  Residual phase corrections delta_{lm} up to l=m=5 for higher modes.
 *  Ref: arXiv:2001.09082
 *  
 *  @param[in]  Hreal : 
 *  @param[in]  Omega : orbital frequency
 *  @param[in]  nu    : symmetric mass ratio
 *  @param[out] dlm   : residual phase corrections
*/
void eob_wav_deltalm_HM(double Hreal,double Omega,double nu, double *dlm)
{
    
  /* Useful shorthands*/
  const double Pi2 = SQ(Pi);
  double nu2    = SQ(nu);
  double y      = cbrt(Hreal*Omega*Hreal*Omega);
  double sqrt_y = sqrt(y);
  double y3     = y*y*y;
  double y32    = Hreal*Omega;
  
  /* Leading order contributions*/
  double delta22LO = 7./3.   * y32;
  double delta21LO = 2./3.   * y32;
  double delta33LO = 13./10. * y32;
  double delta32LO = (10. + 33.*nu)/(15.*(1. - 3.*nu))*y32;        
  double delta31LO = 13./30. * y32;
  double delta42LO = 7.*(1. + 6.*nu)/(15.*(1. - 3.*nu))*y32;
  double delta43LO = (486. + 4961.*nu)/(810.*(1. - 2.*nu))*y32;
  double delta44LO = (112. + 219.*nu)/(120.*(1. - 3.*nu))*y32;
	  
  /* Init phase */
#pragma omp simd
  for (int k = 0; k < KMAX; k++) {
    dlm[k] = 0.;
  }
  
  /* Residual phases in Pade-resummed form when possible */
  double num;
  double den;

  /* l=2 */
  /* Pade(2,1) approximant */
  num        = 856*Pi*sqrt_y + 2625*nu*(4.-75.*nu*y);
  den        = 856*Pi*sqrt_y + 10500*nu;
  dlm[0] = delta21LO*num/den;
  /* Pade(2,2) approximant */
  num        = (808920.*nu*Pi*sqrt(y) + 137388.*Pi2*y + 35.*nu2*(136080. + (154975. - 1359276.*nu)*y));
  den        = (808920.*nu*Pi*sqrt(y) + 137388.*Pi2*y + 35.*nu2*(136080. + (154975. + 40404.*nu)*y));
  dlm[1] = delta22LO*num/den;
	  
  /* l=3 */
  /* Pade(1,2) approximant */
  num        = 4641.*nu + 1690.*Pi*sqrt_y;
  den        = num + 18207.*nu2*y;
  dlm[2] = delta31LO*num/den;
  /* Pade(0,2) approximant */
  den   = 1. + (260./7.)/(10. + 33.*nu)*(-1. + 3.*nu)*Pi*y32
    +(8./1323.)/SQ((10. + 33.*nu))*(558110. + 167433.*nu -5525289.*nu2
				    +146250.*Pi2 - 1393470.*nu*Pi2 + 2864160.*nu2*Pi2)*y3;
  dlm[3] = delta32LO/den; 
  /* Pade(1,2) approximant */
  num        = 1.  + 94770.*Pi/(566279.*nu)*sqrt_y;
  den        = num + 80897.* nu/3159.*y;
	
	
  dlm[4] = delta33LO*num/den;
	  
  /* l=4 */
  dlm[5] =   (2.+507.*nu)/(10.*(1.-2.*nu))*y32   + 1571./3465.*Pi*y3;
  /* Pade(0,2) approximant */
  den = 1. + (6284./1617.)*(-1. + 3.*nu)/(1. + 6.*nu)*Pi*y32 +
    1./21789075.*(-1. + 3.*nu)/SQ(1. + 6.*nu)*
    (-858240537. -5149443222.*nu - 216169600.*Pi2 + 1664631600.*nu*Pi2)*y3;
  dlm[6] =  delta42LO/den;
  /* Pade(0,2) approximant */
  den = 1. + 254502./77./(486. + 4961.*nu)*(-1. + 2.*nu)*Pi*y32
    + 27./59290./SQ(486. + 4961.*nu)*(53627484834. + 440164685191.*nu - 1094839309718.*nu2 + 16934563080.*Pi2 -153861917640.*nu*Pi2 + 239985582960.*nu2*Pi2)*y3;
  dlm[7] =  delta43LO/den;
  /* Pade(0,2) approximant */
  den = 1. + 201088./231.*(-1. + 3.*nu)/(112. + 219.*nu)*Pi*y32
    + 64./444675.*(-1. + 3.*nu)/SQ(112. + 219.*nu)*
    (-13731848592. -26850668229.*nu -3458713600.*Pi2 + 19327698800.*nu*Pi2)*y3;
  dlm[8] =  delta44LO/den;
	  
  /* l=5 */
  dlm[13] = (96875. + 857528.*nu)/(131250.*(1.-2.*nu))*y32;
  
}

/**
 * Function: eob_wav_hlm_v1
 * ------------------------
 * Resummed amplitudes in the general nu-dependent case.
 * Refs:
 *  * Damour, Iyer & Nagar, PRD 79, 064004 (2009)     [theory]
 *  * Fujita & Iyer, PRD 82, 044051 (2010)            [test-mass 5.5PN]
 *  * Damour, Nagar & Bernuzzi, PRD 87, 084035 (2013) [complete information]
 * 
 * @param[in]  x    : frequency
 * @param[in]  nu   : symmetric mass ratio
 * @param[out] rholm: residual amplitudes
 * @param[out] flm  : residual amplitudes 
 */
void eob_wav_flm_v1(double x,double nu, double *rholm, double *flm)
{

  /* Coefficients */
  static double clm[KMAX][6];
  
  const double nu2 = nu*nu;
  const double nu3 = nu*nu2;
  const double nu4 = nu*nu3;
  const double Pi2 = SQ(Pi);
  
  //static int firstcall = 1;
  if (EOBPars->firstcall[FIRSTCALL_EOBWAVFLMV1]) {
    if (0) printf("Precompute some rholm coefs\n");
    EOBPars->firstcall[FIRSTCALL_EOBWAVFLMV1] = 0;
    
    for (int k=0; k<KMAX; k++) clm[k][0] = 1.;
    for (int k=0; k<KMAX; k++) for (int n=1; n<6; n++) clm[k][n] = 0.;

    /* (2,1) */
    clm[0][1] = (-1.0535714285714286 + 0.27380952380952384 *nu);
    clm[0][2] = (-0.8327841553287982 - 0.7789824263038548  *nu + 0.13116496598639457*nu2);
    /* clm[0][3] = (2.9192806270460925  - 1.019047619047619   *el1); */
    /* clm[0][4] = (-1.28235780892213   + 1.073639455782313   *el1); */
    /* clm[0][5] = (-3.8466571723355227 + 0.8486467106683944  *el1)*PMTERMS_eps; */

    /* (2,2) */
    clm[1][1] = (-1.0238095238095237 + 0.6547619047619048*nu);
    clm[1][2] = (-1.94208238851096   - 1.5601379440665155*nu + 0.4625614134542706*nu2);
    /* clm[1][3] = (12.736034731834051  - 2.902228713904598 *nu - 1.9301558466099282*nu2 + 0.2715020968103451*nu3 - 4.076190476190476*el2); */
    /* clm[1][4] = (-2.4172313935587004 + 4.173242630385488 *el2); */
    /* clm[1][5] = (-30.14143102836864  + 7.916297736025627 *el2); */

    /* (3,1) */
    clm[2][1] = (-0.7222222222222222 - 0.2222222222222222*nu);
    clm[2][2] = (0.014169472502805836 - 0.9455667789001122*nu - 0.46520763187429853*nu2);
    /* clm[2][3] = (1.9098284139598072 - 0.4126984126984127*el1); */
    /* clm[2][4] = (0.5368150316615179 + 0.2980599647266314*el1); */
    /* clm[2][5] = (1.4497991763035063 - 0.0058477188106817735*el1)*PMTERMS_eps; */

    /* (3,2) */
    clm[3][1] = (0.003703703703703704*(328. - 1115.*nu + 320.*nu2))/(-1. + 3.*nu);
    clm[3][2] = (6.235191420376606e-7*(-1.444528e6 + 8.050045e6*nu - 4.725605e6*nu2 - 2.033896e7*nu3 + 3.08564e6*nu4))/((-1. + 3.*nu)*(-1. + 3.*nu));
    /* clm[3][3] = (6.220997955214429 - 1.6507936507936507*el2); */
    /* clm[3][4] = (-3.4527288879001268 + 2.005408583186361*el2)*PMTERMS_eps; */

    /* (3,3) */
    clm[4][1] =  (-1.1666666666666667 + 0.6666666666666666*nu);
    clm[4][2] = (-1.6967171717171716 - 1.8797979797979798*nu + 0.45151515151515154*nu2);
    /* clm[4][3] = (14.10891386831863 - 3.7142857142857144*el3); */
    /* clm[4][4] = (-6.723375314944128 + 4.333333333333333*el3); */
    /* clm[4][5] = (-29.568699895427518 + 6.302092352092352*el3)*PMTERMS_eps; */

    /* (4,1) */
    clm[5][1] = (0.001893939393939394*(602. - 1385.*nu + 288.*nu2))/(-1. + 2.*nu);
    clm[5][2] = (- 0.36778992787515513);
    /* clm[5][3] = (0.6981550175535535 - 0.2266955266955267*el1); */
    /* clm[5][4] = (-0.7931524512893319 + 0.2584672482399755*el1)*PMTERMS_eps; */

    /* (4,2) */
    clm[6][1] = (0.0007575757575757576*(1146. - 3530.*nu + 285.*nu2))/(-1. + 3.*nu);
    clm[6][2] = - (3.1534122443213353e-9*(1.14859044e8 - 2.95834536e8*nu - 1.204388696e9*nu2 + 3.04798116e9*nu3 + 3.79526805e8*nu4))/((-1. + 3.*nu)*(-1. + 3.*nu));
    /* clm[6][3] = 4.550378418934105e-12*(8.48238724511e11 - 1.9927619712e11*el2); */
    /* clm[6][4] = (-0.6621921297263365 + 0.787251738160829*el2)*PMTERMS_eps; */
    
    /* (4,3) */
    clm[7][1] = (0.005681818181818182*(222. - 547.*nu + 160.*nu2))/(-1. + 2.*nu);
    clm[7][2] = (- 0.9783218202252293);
    /* clm[7][3] = (8.519456157072423 - 2.0402597402597404*el3)*PMTERMS_eps; */
    /* clm[7][4] = (-5.353216984886716 + 2.5735094451003544*el3)*PMTERMS_eps; */

    /* (4,4) */
    clm[8][1] = (0.0007575757575757576*(1614. - 5870.*nu + 2625.*nu2))/(-1. + 3.*nu);
    clm[8][2] = (3.1534122443213353e-9*(-5.11573572e8 + 2.338945704e9*nu - 3.13857376e8*nu2 - 6.733146e9*nu3 + 1.252563795e9*nu4))/((-1. + 3.*nu)*(-1. + 3.*nu));
    /* clm[8][3] = (15.108111214795123 - 3.627128427128427*el4); */
    /* clm[8][4] = (-8.857121657199649 + 4.434988849534304*el4)*PMTERMS_eps; */

    /* (5,1) */
    clm[9][1] = (0.002564102564102564*(319. - 626.*nu + 8.*nu2))/(-1. + 2.*nu);
    clm[9][2] = (- 0.1047896120973044);
    /* clm[9][3] = (0.642701885362399 - 0.14414918414918415*el1)*PMTERMS_eps; */
    /* clm[9][4] = (-0.07651588046467575 + 0.11790664036817883*el1)*PMTERMS_eps; */

    /* (5,2) */
    clm[10][1] = (0.00007326007326007326*(-15828. + 84679.*nu - 104930.*nu2 + 21980.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[10][2] = (- 0.4629337197600934)*PMTERMS_eps; 
    /* clm[10][3] = (2.354458371550237 - 0.5765967365967366*el2)*PMTERMS_eps; */

    /* (5,3) */
    clm[11][1] = (0.002564102564102564*(375. - 850.*nu + 176.*nu2))/(-1. + 2.*nu);
    clm[11][2] = (- 0.5788010707241477);
    /* clm[11][3] = (5.733973288504755 - 1.2973426573426574*el3)*PMTERMS_eps; */
    /* clm[11][4] = (-1.9573287625526001 + 1.2474448628294783*el3)*PMTERMS_eps; */

    /* (5,4) */
    clm[12][1] = (0.00007326007326007326*(-17448. + 96019.*nu - 127610.*nu2 + 33320.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[12][2] = (- 1.0442142414362194)*PMTERMS_eps;
    /* clm[12][3] = (10.252052781721588 - 2.3063869463869464*el4)*PMTERMS_eps; */

    /* (5,5) */
    clm[13][1] = (0.002564102564102564*(487. - 1298.*nu + 512.*nu2))/(-1. + 2.*nu);
    clm[13][2] = (- 1.5749727622804546);
    /* clm[13][3] = (15.939827047208668 - 3.6037296037296036*el5)*PMTERMS_eps; */
    /* clm[13][4] = (-10.272578060123237 + 4.500041838503377*el5)*PMTERMS_eps; */
    
    /* (6,1) */
    clm[14][1] = (0.006944444444444444*(-161. + 694.*nu - 670.*nu2 + 124.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[14][2] = (- 0.29175486850885135)*PMTERMS_eps;
    /* clm[14][3] = (0.21653486654395454 - 0.10001110001110002*el1)*PMTERMS_eps; */

    /* (6,2) */
    clm[15][1] = (0.011904761904761904*(-74. + 378.*nu - 413.*nu2 + 49.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[15][2] = ( - 0.24797525070634313)*PMTERMS_eps;
    /* clm[15][3] = (1.7942694138754138 - 0.40004440004440006*el2)*PMTERMS_eps; */

    /* (6,3) */
    clm[16][1] = (0.006944444444444444*(-169. + 742.*nu - 750.*nu2 + 156.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[16][2] = (- 0.5605554442947213)*PMTERMS_eps;
    /* clm[16][3] = (4.002558222882566 - 0.9000999000999002*el3)*PMTERMS_eps; */

    /* (6,4) */
    clm[17][1] = (0.011904761904761904*(-86. + 462.*nu - 581.*nu2 + 133.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[17][2] = (- 0.7228451986855349)*PMTERMS_eps;
    /* clm[17][3] = (7.359388663371044 - 1.6001776001776002*el4)*PMTERMS_eps; */

    /* (6,5) */
    clm[18][1] = (0.006944444444444444*(-185. + 838.*nu - 910.*nu2 + 220.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[18][2] = (- 1.0973940686333457)*PMTERMS_eps;
    /* clm[18][3] = (11.623366217471297 - 2.5002775002775004*el5)*PMTERMS_eps; */

    /* (6,6) */
    clm[19][1] = (0.011904761904761904*(-106. + 602.*nu - 861.*nu2 + 273.*nu3))/(1. - 5.*nu + 5.*nu2); 
    clm[19][2] = (- 1.5543111183867486)*PMTERMS_eps;
    /* clm[19][3] = (16.645950799433503 - 3.6003996003996006*el6)*PMTERMS_eps; */

    /* (7,1) */
    clm[20][1] = (0.0014005602240896359*(-618. + 2518.*nu - 2083.*nu2 + 228.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[20][2] = ( - 0.1508235111143767)*PMTERMS_eps;
    /* clm[20][3] = (0.2581280702019663 - 0.07355557607658449*el1)*PMTERMS_eps; */

    /* (7,2) */
    clm[21][1] = (0.00006669334400426837*(16832. - 123489.*nu + 273924.*nu2 - 190239.*nu3 + 32760.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[21][2] = (- 0.351319484450667)*PMTERMS_eps;
    
    /* (7,3) */
    clm[22][1] = (0.0014005602240896359*(-666. + 2806.*nu - 2563.*nu2 + 420.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[22][2] = (- 0.37187416047628863)*PMTERMS_eps;
    /* clm[22][3] = (3.0835293524055283 - 0.6620001846892604*el3)*PMTERMS_eps; */

    /* (7,4) */
    clm[23][1] = (0.00006669334400426837*(17756. - 131805.*nu + 298872.*nu2 - 217959.*nu3 + 41076.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[23][2] = (- 0.6473746896670599)*PMTERMS_eps;
    
    /* (7,5) */
    clm[24][1] = (0.0014005602240896359*(-762. + 3382.*nu - 3523.*nu2 + 804.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[24][2] = (- 0.8269193364414116)*PMTERMS_eps;
    /* clm[24][3] = (8.750589067052443 - 1.838889401914612*el5)*PMTERMS_eps; */

    /* (7,6) */
    clm[25][1] = (0.0006002400960384153*(2144. - 16185.*nu + 37828.*nu2 - 29351.*nu3 + 6104.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[25][2] = (- 1.1403265020692532)*PMTERMS_eps;
    
    /* (7,7) */
    clm[26][1] = (0.0014005602240896359*(-906. + 4246.*nu - 4963.*nu2 + 1380.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[26][2] = (- 1.5418467934923434)*PMTERMS_eps;
    /* clm[26][3] = (17.255875091408523 - 3.6042232277526396*el7)*PMTERMS_eps; */

    /* (8,1) */
    clm[27][1] = (0.00005482456140350877*(20022. - 126451.*nu + 236922.*nu2 - 138430.*nu3 + 21640.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
    clm[27][2] = (- 0.26842133517043704)*PMTERMS_eps;

    /* (8,2) */
    clm[28][1] = (0.0003654970760233918*(2462. - 17598.*nu + 37119.*nu2 - 22845.*nu3 + 3063.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[28][2] = (- 0.2261796441029474)*PMTERMS_eps;

    /* (8,3) */
    clm[29][1] = (0.00005482456140350877*(20598. - 131059.*nu + 249018.*nu2 - 149950.*nu3 + 24520.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
    clm[29][2] = (- 0.4196774909106648)*PMTERMS_eps;

    /* (8,4) */
    clm[30][1] = (0.0003654970760233918*(2666. - 19434.*nu + 42627.*nu2 - 28965.*nu3 + 4899.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[30][2] = (- 0.47652059150068155)*PMTERMS_eps;

    /* (8,5) */
    clm[31][1] = (0.00027412280701754384*(4350. - 28055.*nu + 54642.*nu2 - 34598.*nu3 + 6056.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
    clm[31][2] = (- 0.7220789990670207)*PMTERMS_eps;

    /* (8,6) */
    clm[32][1] = (0.0010964912280701754*(1002. - 7498.*nu + 17269.*nu2 - 13055.*nu3 + 2653.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[32][2] = (- 0.9061610303170207)*PMTERMS_eps;

    /* (8,7) */
    clm[33][1] = (0.00005482456140350877*(23478. - 154099.*nu + 309498.*nu2 - 207550.*nu3 + 38920.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
    clm[33][2] = (- 1.175404252991305)*PMTERMS_eps;

    /* (8,8) */
    clm[34][1] = (0.0003654970760233918*(3482. - 26778.*nu + 64659.*nu2 - 53445.*nu3 + 12243.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[34][2] = (- 1.5337092502821381)*PMTERMS_eps;

  }

  /* Compute EulerLogs */
  const double el1 = Eulerlog(x,1);
  const double el2 = Eulerlog(x,2);
  const double el3 = Eulerlog(x,3);
  const double el4 = Eulerlog(x,4);
  const double el5 = Eulerlog(x,5);
  const double el6 = Eulerlog(x,6);
  const double el7 = Eulerlog(x,7);

  /* Coefs with Eulerlogs */
  clm[0][3] = (2.9192806270460925  - 1.019047619047619   *el1);
  clm[0][4] = (-1.28235780892213   + 1.073639455782313   *el1);
  clm[0][5] = (-3.8466571723355227 + 0.8486467106683944  *el1)*PMTERMS_eps;
  
  clm[1][3] = (12.736034731834051  - 2.902228713904598 *nu - 1.9301558466099282*nu2 + 0.2715020968103451*nu3 - 4.076190476190476*el2);
  clm[1][4] = (-2.4172313935587004 + 4.173242630385488 *el2);
  clm[1][5] = (-30.14143102836864  + 7.916297736025627 *el2);

  clm[2][3] = (1.9098284139598072 - 0.4126984126984127*el1+ (-4.646868015386534 + (0.21354166666666666)*Pi2)*nu + 2.3020866307903347*nu2 - 0.5813492634480288*nu3);  
  clm[2][4] = (0.5368150316615179 + 0.2980599647266314*el1);
  clm[2][5] = (1.4497991763035063 - 0.0058477188106817735*el1)*PMTERMS_eps;
    
  clm[3][3] = (6.220997955214429 - 1.6507936507936507*el2);
  clm[3][4] = (-3.4527288879001268 + 2.005408583186361*el2)*PMTERMS_eps;
  
  clm[4][3] = (14.10891386831863 - 3.7142857142857144*el3 + (-5.031429681429682 + (0.21354166666666666)*Pi2)*nu - 1.7781727531727531*nu2 + 0.25923767590434255*nu3);
  clm[4][4] = (-6.723375314944128 + 4.333333333333333*el3);
  clm[4][5] = (-29.568699895427518 + 6.302092352092352*el3)*PMTERMS_eps;
  
  clm[5][3] = (0.6981550175535535 - 0.2266955266955267*el1);
  clm[5][4] = (-0.7931524512893319 + 0.2584672482399755*el1)*PMTERMS_eps;
  
  clm[6][3] = 4.550378418934105e-12*(8.48238724511e11 - 1.9927619712e11*el2);
  clm[6][4] = (-0.6621921297263365 + 0.787251738160829*el2)*PMTERMS_eps;
  
  clm[7][3] = (8.519456157072423 - 2.0402597402597404*el3)*PMTERMS_eps;
  clm[7][4] = (-5.353216984886716 + 2.5735094451003544*el3)*PMTERMS_eps;
  
  clm[8][3] = (15.108111214795123 - 3.627128427128427*el4);
  clm[8][4] = (-8.857121657199649 + 4.434988849534304*el4)*PMTERMS_eps;
  
  clm[9][3] = (0.642701885362399 - 0.14414918414918415*el1)*PMTERMS_eps;
  clm[9][4] = (-0.07651588046467575 + 0.11790664036817883*el1)*PMTERMS_eps;
  
  clm[10][3] = (2.354458371550237 - 0.5765967365967366*el2)*PMTERMS_eps;
  
  clm[11][3] = (5.733973288504755 - 1.2973426573426574*el3)*PMTERMS_eps;
  clm[11][4] = (-1.9573287625526001 + 1.2474448628294783*el3)*PMTERMS_eps;
  
  clm[12][3] = (10.252052781721588 - 2.3063869463869464*el4)*PMTERMS_eps;
  
  clm[13][3] = (15.939827047208668 - 3.6037296037296036*el5)*PMTERMS_eps;
  clm[13][4] = (-10.272578060123237 + 4.500041838503377*el5)*PMTERMS_eps;
  
  clm[14][3] = (0.21653486654395454 - 0.10001110001110002*el1)*PMTERMS_eps;
  
  clm[15][3] = (1.7942694138754138 - 0.40004440004440006*el2)*PMTERMS_eps;
  
  clm[16][3] = (4.002558222882566 - 0.9000999000999002*el3)*PMTERMS_eps;
  
  clm[17][3] = (7.359388663371044 - 1.6001776001776002*el4)*PMTERMS_eps;
  
  clm[18][3] = (11.623366217471297 - 2.5002775002775004*el5)*PMTERMS_eps;
  
  clm[19][3] = (16.645950799433503 - 3.6003996003996006*el6)*PMTERMS_eps;
  
  clm[20][3] = (0.2581280702019663 - 0.07355557607658449*el1)*PMTERMS_eps;
  
  clm[22][3] = (3.0835293524055283 - 0.6620001846892604*el3)*PMTERMS_eps;
  
  clm[24][3] = (8.750589067052443 - 1.838889401914612*el5)*PMTERMS_eps;
  
  clm[26][3] = (17.255875091408523 - 3.6042232277526396*el7)*PMTERMS_eps;
    
  /* rho_lm */
  const double x2  = x*x;
  const double x3  = x*x2;
  const double x4  = x*x3;
  const double x5  = x*x4;
  const double xn[] = {1.,x,x2,x3,x4,x5};

  for (int k=0; k<KMAX; k++) {
    /* Note: the two sums give different result */
#if (1)
    rholm[k] = clm[k][0];
    for (int n=1; n<6; n++) {
      rholm[k] += clm[k][n]*xn[n];
    }	     
#else  
    rholm[k] = x5*clm[k][5];
    for (int n=5; n-- >1; ) { // 4,3,2,1 // 
      rholm[k] += clm[k][n]*xn[n];
    } 
    rholm[k] += clm[k][0]; 
#endif
  }

  /* Amplitudes */
#pragma omp simd
  for (int k = 0; k < KMAX; k++) {
      flm[k] = gsl_pow_int(rholm[k], LINDEX[k]);
  }

}

/**
 * Function: eob_wav_flm_HM
 * -------------------------
 *   Resummed amplitudes in the general nu-dependent case.
 *   Function introduced for higher modes
 *   (2,2) is taken at 3^{+2} PN order
 *   (4,4), (4,2) and (5,5) and l >= 7 modes are Talyor expanded
 *   (3,2) is resummed with a Padé 32
 *   (2,1) is resummed with a Padé 51
 *   All other modes are resummed with a Padé 42
 *   Ref: arXiv:2001.09082
 *   
 *   @param[in] x    : frequency
 *   @param[in] nu   : symmetric mass ratio
 *   @param[out] rholm : resummed amplitudes
 *   @param[out] flm   : resummed amplitudes
*/
void eob_wav_flm_HM(double x,double nu, double *rholm, double *flm)
{
  /** Coefficients */
  static double clm[KMAX][7];
  
  const double nu2 = nu*nu;
  const double nu3 = nu*nu2;
  const double nu4 = nu*nu3;
  const double Pi2 = SQ(Pi);
  
  //static int firstcall = 1;
  if (EOBPars->firstcall[FIRSTCALL_EOBWAVFLMHM]) {
    if (0) printf("Precompute some rholm coefs\n");
    EOBPars->firstcall[FIRSTCALL_EOBWAVFLMHM] = 0;
    
    for (int k=0; k<KMAX; k++) clm[k][0] = 1.;
    for (int k=0; k<KMAX; k++) for (int n=1; n<7; n++) clm[k][n] = 0.;

    /** (2,1) */
    clm[0][1] = (-1.0535714285714286 + 0.27380952380952384 *nu);
    clm[0][2] = (-0.8327841553287982 - 0.7789824263038548  *nu + 0.13116496598639457*nu2);

    /** (2,2) */
    clm[1][1] = (-1.0238095238095237 + 0.6547619047619048*nu);
    clm[1][2] = (-1.94208238851096   - 1.5601379440665155*nu + 0.4625614134542706*nu2);

    /** (3,1) */
    clm[2][1] = (-0.7222222222222222 - 0.2222222222222222*nu);
    clm[2][2] = (0.014169472502805836 - 0.9455667789001122*nu - 0.46520763187429853*nu2);

    /** (3,2) */
    clm[3][1] = (0.003703703703703704*(328. - 1115.*nu + 320.*nu2))/(-1. + 3.*nu);
    clm[3][2] = (6.235191420376606e-7*(-1.444528e6 + 8.050045e6*nu - 4.725605e6*nu2 - 2.033896e7*nu3 + 3.08564e6*nu4))/((-1. + 3.*nu)*(-1. + 3.*nu));

    /** (3,3) */
    clm[4][1] =  (-1.1666666666666667 + 0.6666666666666666*nu);
    clm[4][2] = (-1.6967171717171716 - 1.8797979797979798*nu + 0.45151515151515154*nu2);
    
    /** (4,1) */
    clm[5][1] = (0.001893939393939394*(602. - 1385.*nu + 288.*nu2))/(-1. + 2.*nu);
    clm[5][2] = (- 0.36778992787515513);
    
    /** (4,2) */
    clm[6][1] = (0.0007575757575757576*(1146. - 3530.*nu + 285.*nu2))/(-1. + 3.*nu);
    clm[6][2] = - (3.1534122443213353e-9*(1.14859044e8 - 2.95834536e8*nu - 1.204388696e9*nu2 + 3.04798116e9*nu3 + 3.79526805e8*nu4))/((-1. + 3.*nu)*(-1. + 3.*nu));
        
    /** (4,3) */
    clm[7][1] = (0.005681818181818182*(222. - 547.*nu + 160.*nu2))/(-1. + 2.*nu);
    clm[7][2] = (- 0.9783218202252293);
    
    /** (4,4) */
    clm[8][1] = (0.0007575757575757576*(1614. - 5870.*nu + 2625.*nu2))/(-1. + 3.*nu);
    clm[8][2] = (3.1534122443213353e-9*(-5.11573572e8 + 2.338945704e9*nu - 3.13857376e8*nu2 - 6.733146e9*nu3 + 1.252563795e9*nu4))/((-1. + 3.*nu)*(-1. + 3.*nu));
    
    /** (5,1) */
    clm[9][1] = (0.002564102564102564*(319. - 626.*nu + 8.*nu2))/(-1. + 2.*nu);
    clm[9][2] = (- 0.1047896120973044);
    
    /** (5,2) */
    clm[10][1] = (0.00007326007326007326*(-15828. + 84679.*nu - 104930.*nu2 + 21980.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[10][2] = (- 0.4629337197600934)*PMTERMS_eps; 
    
    /** (5,3) */
    clm[11][1] = (0.002564102564102564*(375. - 850.*nu + 176.*nu2))/(-1. + 2.*nu);
    clm[11][2] = (- 0.5788010707241477);
    
    /** (5,4) */
    clm[12][1] = (0.00007326007326007326*(-17448. + 96019.*nu - 127610.*nu2 + 33320.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[12][2] = (- 1.0442142414362194)*PMTERMS_eps;
    
    /** (5,5) */
    clm[13][1] = (0.002564102564102564*(487. - 1298.*nu + 512.*nu2))/(-1. + 2.*nu);
    clm[13][2] = (- 1.5749727622804546);
        
    /** (6,1) */
    clm[14][1] = (0.006944444444444444*(-161. + 694.*nu - 670.*nu2 + 124.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[14][2] = (- 0.29175486850885135)*PMTERMS_eps;
    
    /** (6,2) */
    clm[15][1] = (0.011904761904761904*(-74. + 378.*nu - 413.*nu2 + 49.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[15][2] = ( - 0.24797525070634313)*PMTERMS_eps;
    
    /** (6,3) */
    clm[16][1] = (0.006944444444444444*(-169. + 742.*nu - 750.*nu2 + 156.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[16][2] = (- 0.5605554442947213)*PMTERMS_eps;
    
    /** (6,4) */
    clm[17][1] = (0.011904761904761904*(-86. + 462.*nu - 581.*nu2 + 133.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[17][2] = (- 0.7228451986855349)*PMTERMS_eps;
    
    /** (6,5) */
    clm[18][1] = (0.006944444444444444*(-185. + 838.*nu - 910.*nu2 + 220.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[18][2] = (- 1.0973940686333457)*PMTERMS_eps;
  
    /** (6,6) */
    clm[19][1] = (0.011904761904761904*(-106. + 602.*nu - 861.*nu2 + 273.*nu3))/(1. - 5.*nu + 5.*nu2); 
    clm[19][2] = (- 1.5543111183867486)*PMTERMS_eps;
  
    /** (7,1) */
    clm[20][1] = (0.0014005602240896359*(-618. + 2518.*nu - 2083.*nu2 + 228.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[20][2] = ( - 0.1508235111143767)*PMTERMS_eps;
  
    /** (7,2) */
    clm[21][1] = (0.00006669334400426837*(16832. - 123489.*nu + 273924.*nu2 - 190239.*nu3 + 32760.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[21][2] = (- 0.351319484450667)*PMTERMS_eps;
    
    /** (7,3) */
    clm[22][1] = (0.0014005602240896359*(-666. + 2806.*nu - 2563.*nu2 + 420.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[22][2] = (- 0.37187416047628863)*PMTERMS_eps;
  
    /** (7,4) */
    clm[23][1] = (0.00006669334400426837*(17756. - 131805.*nu + 298872.*nu2 - 217959.*nu3 + 41076.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[23][2] = (- 0.6473746896670599)*PMTERMS_eps;
    
    /** (7,5) */
    clm[24][1] = (0.0014005602240896359*(-762. + 3382.*nu - 3523.*nu2 + 804.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[24][2] = (- 0.8269193364414116)*PMTERMS_eps;
  
    /** (7,6) */
    clm[25][1] = (0.0006002400960384153*(2144. - 16185.*nu + 37828.*nu2 - 29351.*nu3 + 6104.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[25][2] = (- 1.1403265020692532)*PMTERMS_eps;
    
    /** (7,7) */
    clm[26][1] = (0.0014005602240896359*(-906. + 4246.*nu - 4963.*nu2 + 1380.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[26][2] = (- 1.5418467934923434)*PMTERMS_eps;
  
    /** (8,1) */
    clm[27][1] = (0.00005482456140350877*(20022. - 126451.*nu + 236922.*nu2 - 138430.*nu3 + 21640.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
    clm[27][2] = (- 0.26842133517043704)*PMTERMS_eps;

    /** (8,2) */
    clm[28][1] = (0.0003654970760233918*(2462. - 17598.*nu + 37119.*nu2 - 22845.*nu3 + 3063.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[28][2] = (- 0.2261796441029474)*PMTERMS_eps;

    /** (8,3) */
    clm[29][1] = (0.00005482456140350877*(20598. - 131059.*nu + 249018.*nu2 - 149950.*nu3 + 24520.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
    clm[29][2] = (- 0.4196774909106648)*PMTERMS_eps;

    /** (8,4) */
    clm[30][1] = (0.0003654970760233918*(2666. - 19434.*nu + 42627.*nu2 - 28965.*nu3 + 4899.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[30][2] = (- 0.47652059150068155)*PMTERMS_eps;

    /** (8,5) */
    clm[31][1] = (0.00027412280701754384*(4350. - 28055.*nu + 54642.*nu2 - 34598.*nu3 + 6056.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
    clm[31][2] = (- 0.7220789990670207)*PMTERMS_eps;

    /** (8,6) */
    clm[32][1] = (0.0010964912280701754*(1002. - 7498.*nu + 17269.*nu2 - 13055.*nu3 + 2653.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[32][2] = (- 0.9061610303170207)*PMTERMS_eps;

    /** (8,7) */
    clm[33][1] = (0.00005482456140350877*(23478. - 154099.*nu + 309498.*nu2 - 207550.*nu3 + 38920.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
    clm[33][2] = (- 1.175404252991305)*PMTERMS_eps;

    /** (8,8) */
    clm[34][1] = (0.0003654970760233918*(3482. - 26778.*nu + 64659.*nu2 - 53445.*nu3 + 12243.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[34][2] = (- 1.5337092502821381)*PMTERMS_eps;

  }

  /** Compute EulerLogs */
  const double el1 = Eulerlog(x,1);
  const double el2 = Eulerlog(x,2);
  const double el3 = Eulerlog(x,3);
  const double el4 = Eulerlog(x,4);
  const double el5 = Eulerlog(x,5);
  const double el6 = Eulerlog(x,6);
  const double el7 = Eulerlog(x,7);

  const double logx  = log(x);
  const double log2x = logx*logx;
  
  /** Coefs with Eulerlogs */
  /* l = 2 */
  clm[0][3] = (2.9192806270460925  - 1.019047619047619   *el1);
  clm[0][4] = (-1.28235780892213   + 1.073639455782313   *el1);
  clm[0][5] = (-3.8466571723355227 + 0.8486467106683944  *el1)*PMTERMS_eps;
  clm[0][6] = -0.53614531347435562831 - 1.5600367264240624649*logx + 0.1298072562358276644*log2x;
  
  clm[1][3] = (12.736034731834051  - 2.902228713904598 *nu - 1.9301558466099282*nu2 + 0.2715020968103451*nu3 - 4.076190476190476*el2);
  clm[1][4] = (-2.4172313935587004 + 4.173242630385488 *el2);
  clm[1][5] = (-30.14143102836864  + 7.916297736025627 *el2);

  /* l = 3 */
  clm[2][3] = (1.9098284139598072 - 0.4126984126984127*el1+ (-4.646868015386534 + (0.21354166666666666)*Pi2)*nu + 2.3020866307903347*nu2 - 0.5813492634480288*nu3);  
  clm[2][4] = (0.5368150316615179 + 0.2980599647266314*el1);
  clm[2][5] = (1.4497991763035063 - 0.0058477188106817735*el1)*PMTERMS_eps;
    
  clm[3][3] = (6.220997955214429 - 1.6507936507936507*el2);
  clm[3][4] = (-3.4527288879001268 + 2.005408583186361*el2)*PMTERMS_eps;
  clm[3][5] = -4.5270419165221037846 + 0.74342737585535939445*logx;
  clm[3][6] = 2.9819373321032360597 - 4.1746594065141637002*logx + 0.34063995968757873520*log2x;
    
  clm[4][3] = (14.10891386831863 - 3.7142857142857144*el3 + (-5.031429681429682 + (0.21354166666666666)*Pi2)*nu - 1.7781727531727531*nu2 + 0.25923767590434255*nu3);
  clm[4][4] = (-6.723375314944128 + 4.333333333333333*el3);
  clm[4][5] = (-29.568699895427518 + 6.302092352092352*el3)*PMTERMS_eps;
  clm[4][6] = 4.2434215365016546735 - 18.544715779075533606*logx + 1.7244897959183673469*log2x;

  /* l = 4 */
  clm[5][3] = (0.6981550175535535 - 0.2266955266955267*el1);
  clm[5][4] = (-0.7931524512893319 + 0.2584672482399755*el1)*PMTERMS_eps;
  clm[5][5] = -0.73874769521337814809 + 0.041688165706484035882*logx;
  clm[5][6] = -0.60643478766899956315 - 0.078148517222175236907*logx + 0.0064238577279702820828*log2x;
  
  clm[6][3] = 4.550378418934105e-12*(8.48238724511e11 - 1.9927619712e11*el2);
  clm[6][4] = (-0.6621921297263365 + 0.787251738160829*el2)*PMTERMS_eps;
  clm[6][5] = 0.066452311537142230164 + 0.16421729454462630378*logx;
  clm[6][6] = 4.7020949396983490140 - 1.4492579183627638456*logx + 0.10278172364752451333*log2x;
  
  clm[7][3] = (8.519456157072423 - 2.0402597402597404*el3)*PMTERMS_eps;
  clm[7][4] = (-5.353216984886716 + 2.5735094451003544*el3)*PMTERMS_eps;
  clm[7][5] = -6.0091714045534954041 + 0.99801531141158133484*logx;
  clm[7][6] = 4.9466129199935928986 - 6.3245279024451771064*logx + 0.52033247596559284871*log2x;
  
  clm[8][3] = (15.108111214795123 - 3.627128427128427*el4);
  clm[8][4] = (-8.857121657199649 + 4.434988849534304*el4)*PMTERMS_eps;
  clm[8][5] = -14.633690582678763747 + 2.9256460798810267523*logx;
  clm[8][6] = 5.3708202812535269509 - 18.028080626090983076*logx + 1.6445075783603922132*log2x;

  /* l = 5 */
  clm[9][3] = (0.642701885362399 - 0.14414918414918415*el1)*PMTERMS_eps;
  clm[9][4] = (-0.07651588046467575 + 0.11790664036817883*el1)*PMTERMS_eps;
  clm[9][5] =  0.2112256289378054518 + 0.007552668545567954*logx;
  clm[9][6] =  0.6954698480021733513 - 0.045868829656707944*logx + 0.0025973734113594253*log2x;
  
  clm[10][3] = (2.354458371550237 - 0.5765967365967366*el2)*PMTERMS_eps;
  clm[10][4] = -0.36261168338728289637 + 0.33429938266861343784*logx;
  clm[10][5] = -1.2210109278350054619 + 0.13346303603712902107*logx;
  clm[10][6] = 0.33075256130927413448 + 0.041557974581750805527*log2x - 0.55630499019814556741*logx;
  
  clm[11][3] = (5.733973288504755 - 1.2973426573426574*el3)*PMTERMS_eps;
  clm[11][4] = (-1.9573287625526001 + 1.2474448628294783*el3)*PMTERMS_eps;
  clm[11][5] = -1.2893903635541462354 + 0.37545165958302052977*logx;
  clm[11][6] = 5.5345161119771993654 + 0.21038724632011345298*log2x - 2.758182643809987602*logx;
  
  clm[12][3] = (10.252052781721588 - 2.3063869463869464*el4)*PMTERMS_eps;
  clm[12][4] = 0.97104909920221130334 + 1.4740600527677450754*logx;
  clm[12][5] = -7.1507394017332213549 + 1.2041810478399218889*logx;
  clm[12][6] = 6.1586857187986271095 + 0.66492759330801288843*log2x - 8.0193253436649521323*logx;
  
  clm[13][3] = (15.939827047208668 - 3.6037296037296036*el5)*PMTERMS_eps;
  clm[13][4] = (-10.272578060123237 + 4.500041838503377*el5)*PMTERMS_eps;
  clm[13][5] = -14.676578889547799889 + 2.8378879842489309945*logx;
  clm[13][6] = 5.9456506579632350247 + 1.6233583820996408409*log2x - 17.987106103978598244*logx;

  /* l = 6 */
  clm[14][3] = (0.21653486654395454 - 0.10001110001110002*el1)*PMTERMS_eps;
  clm[14][4] = -0.41910058747759839184 + 0.05590898299231632565*logx;
  clm[14][5] = -0.54703676004200114476 + 0.014589362666582032864*logx;
  clm[14][6] = -0.67046877793201698851 - 0.010612705664941704748*logx + 0.001250277515678781080*log2x;
  
  clm[15][3] = (1.7942694138754138 - 0.40004440004440006*el2)*PMTERMS_eps;
  clm[15][4] = 0.11894436114207465266 + 0.17621003335289049575*logx;
  clm[15][5] = -0.020168519020176921191 + 0.049600555197339365552*logx;
  clm[15][6] = 1.3092135879268909484 + 0.020004440250860497281*log2x - 0.29998579925017448295*logx;
  
  clm[16][3] = (4.002558222882566 - 0.9000999000999002*el3)*PMTERMS_eps;
  clm[16][4] = -0.17301414399537528118 + 0.52818362193362193362*logx;
  clm[16][5] = -1.8598787273050426274 + 0.25227794970506685875*logx;
  clm[16][6] = 1.5038409846301927305 + 0.10127247876998126748*log2x - 1.3388799755686528906*logx;
  
  clm[17][3] = (7.359388663371044 - 1.6001776001776002*el4)*PMTERMS_eps;
  clm[17][4] = 1.1251294602654682041 + 0.81913853342424770996*logx;
  clm[17][5] = -2.5262527575555087974 + 0.57834034766625984931*logx;
  clm[17][6] = 6.2193866198899176839 + 0.32007104401376795649*log2x - 4.058212416583905376*logx;

  clm[18][3] = (11.623366217471297 - 2.5002775002775004*el5)*PMTERMS_eps;
  clm[18][4] = 1.1674823841175551524 + 1.6060810331643664977*logx;
  clm[18][5] = -8.0397483562766587382 + 1.3718948493709684487*logx;
  clm[18][6] = 6.9514240915914444298 + 0.78142344729923817503*log2x - 9.3656857835417176404*logx;
  
  clm[19][3] = (16.645950799433503 - 3.6003996003996006*el6)*PMTERMS_eps;
  clm[19][4] = 2.5867656736437215138 + 2.2716807002521288236*logx;
  clm[19][5] = -14.719871717781358093 + 2.7980705647681528943*logx;
  clm[19][6] = 6.3131827696196368972 + 1.6203596603197002797*log2x - 18.073907939714613453*logx;

  /* l = 7 */
  clm[20][3] = (0.2581280702019663 - 0.07355557607658449*el1)*PMTERMS_eps;
  
  clm[22][3] = (3.0835293524055283 - 0.6620001846892604*el3)*PMTERMS_eps;
  
  clm[24][3] = (8.750589067052443 - 1.838889401914612*el5)*PMTERMS_eps;
  
  clm[26][3] = (17.255875091408523 - 3.6042232277526396*el7)*PMTERMS_eps;

  /** rho_lm */
  const double x2  = x*x;
  const double x3  = x*x2;
  const double x4  = x*x3;
  const double x5  = x*x4;
  const double x6  = x*x5;
  const double xn[] = {1.,x,x2,x3,x4,x5,x6};

  /** Initializing Padé approximants */
  double cden,n1,n2,n3,n4,n5,d1,d2 = 0.;
  double k = 0.;

  // Padé (3,2) : Used for (3,1) multipole
  const int kmaxPade32 = 1;
  int kPade32[] = {2};
  
  for (int i=0; i<kmaxPade32; i++) {
    int k = kPade32[i];
    rholm[k] = Pade32(x,clm[k]);
  }

  // Padé (4,2) - Used for (3,2), (3,3), (4,1), (4,3), (5,2), (5,3), (5,4) and l=6 except (6,1) 
  const int kmaxPade42 = 12;
  int kPade42[] = {3,4,5,7,10,11,12,15,16,17,18,19};
  
  for (int i=0; i<kmaxPade42; i++) {
    int k = kPade42[i];
    rholm[k] = Pade42(x,clm[k]);
  }

  // Padé (5,1) - Used for (2,1)
  const int kmaxPade51 = 1;
  int kPade51[] = {0};
  
  for (int i=0; i<kmaxPade51; i++) {
    int k = kPade51[i];
    rholm[k] = Pade51(x,clm[k]);
  }

  
  // Taylor series : (2,2) at 5PN;
  //                 (4,4), (4,2), (5,5), (5,1), (6,1) and l>6  at 6PN
  const int kmaxTaylor = 21;
  int kTaylor[kmaxTaylor];
  kTaylor[0] = 1;
  kTaylor[1] = 6;
  kTaylor[2] = 8;
  kTaylor[3] = 9;
  kTaylor[4] = 13;
  kTaylor[5] = 14;
  for (int i=6; i<kmaxTaylor; i++) {
    kTaylor[i] = 14 + i;
  }
  
  for (int i=0; i<kmaxTaylor; i++) {
    int k = kTaylor[i];
    rholm[k] = Taylorseries(x,clm[k],6);
  }
  
  if (kmaxTaylor+kmaxPade32+kmaxPade42+kmaxPade51 != KMAX) {
    errorexit("Wrong function: not all multipoles are written.\n");
  }
  
  /** Amplitudes */
#pragma omp simd
  for (int k = 0; k < KMAX; k++) {
      flm[k] = gsl_pow_int(rholm[k], LINDEX[k]);
  }

}

/**
 * Function: eob_wav_flm_Kerr
 * --------------------------
 *   Computes the resummed amplitudes for a particle orbiting around Kerr
 * 
 *   @param[in] x     : orbital frequency
 *   @param[in] nu    : symmetric mass ratio  
 *   @param[out] rholm: resummed amplitudes
 *   @param[out] flm  : resummed amplitudes
*/
void eob_wav_flm_Kerr(double x,double nu, double *rholm, double *flm)
{
  /** Coefficients */
  static double clm[KMAX][9];
  
  const double Pi2 = SQ(Pi);
  
  //static int firstcall = 1;
  if (EOBPars->firstcall[FIRSTCALL_EOBWAVFLMKERR]) {
    if (0) printf("Precompute some rholm coefs\n");
    EOBPars->firstcall[FIRSTCALL_EOBWAVFLMKERR] = 0;
    
    for (int k=0; k<KMAX; k++) clm[k][0] = 1.;
    for (int k=0; k<KMAX; k++) for (int n=1; n<7; n++) clm[k][n] = 0.;

    /** (2,1) */
    clm[0][1] = -1.0535714285714286;
    clm[0][2] = -0.8327841553287982;

    /** (2,2) */
    clm[1][1] = -1.0238095238095237;
    clm[1][2] = -1.94208238851096;

    /** (3,1) */
    clm[2][1] = -0.7222222222222222;
    clm[2][2] = 0.014169472502805836;

    /** (3,2) */
    clm[3][1] = -1.2148148148148148;
    clm[3][2] = -0.9006908592093777;

    /** (3,3) */
    clm[4][1] = -1.1666666666666667;
    clm[4][2] = -1.6967171717171716;
    
    /** (4,1) */
    clm[5][1] = -1.1401515151515151;
    clm[5][2] = -0.36778992787515513;
    
    /** (4,2) */
    clm[6][1] = -0.8681818181818182;
    clm[6][2] = -0.362197915720643;
        
    /** (4,3) */
    clm[7][1] = -1.2613636363636365;
    clm[7][2] = -0.9783218202252293;
    
    /** (4,4) */
    clm[8][1] = -1.2227272727272727;
    clm[8][2] = -1.6132023658160024;
    
    /** (5,1) */
    clm[9][1] = -0.8179487179487179;
    clm[9][2] = -0.1047896120973044;
    
    /** (5,2) */
    clm[10][1] = -1.1595604395604395;
    clm[10][2] = -0.4629337197600934; 
    
    /** (5,3) */
    clm[11][1] = -0.9615384615384616;
    clm[11][2] = -0.5788010707241477;
    
    /** (5,4) */
    clm[12][1] = -1.2782417582417582;
    clm[12][2] = -1.0442142414362194;
    
    /** (5,5) */
    clm[13][1] = -1.2487179487179487;
    clm[13][2] = -1.5749727622804546;
        
    /** (6,1) */
    clm[14][1] = -1.1180555555555556;
    clm[14][2] = -0.29175486850885135;
    
    /** (6,2) */
    clm[15][1] = -0.8809523809523809;
    clm[15][2] = -0.24797525070634313;
    
    /** (6,3) */
    clm[16][1] = -1.173611111111111;
    clm[16][2] = -0.5605554442947213;
    
    /** (6,4) */
    clm[17][1] = -1.0238095238095237;
    clm[17][2] = -0.7228451986855349;
    
    /** (6,5) */
    clm[18][1] = -1.284722222222222;
    clm[18][2] = -1.0973940686333457;
  
    /** (6,6) */
    clm[19][1] = -1.2619047619047619; 
    clm[19][2] = -1.5543111183867486;
  
    /** (7,1) */
    clm[20][1] = -0.865546218487395;
    clm[20][2] = -0.1508235111143767;
  
    /** (7,2) */
    clm[21][1] = -1.122582366279845;
    clm[21][2] = -0.351319484450667;
    
    /** (7,3) */
    clm[22][1] = -0.9327731092436975;
    clm[22][2] = -0.37187416047628863;
  
    /** (7,4) */
    clm[23][1] = -1.184207016139789;
    clm[23][2] = -0.6473746896670599;
    
    /** (7,5) */
    clm[24][1] = -1.0672268907563025;
    clm[24][2] = -0.8269193364414116;
  
    /** (7,6) */
    clm[25][1] = -1.2869147659063624;
    clm[25][2] = -1.1403265020692532;
    
    /** (7,7) */
    clm[26][1] = -1.26890756302521;
    clm[26][2] = -1.5418467934923434;
  
    /** (8,1) */
    clm[27][1] = -1.0976973684210527;
    clm[27][2] = -0.26842133517043704;

    /** (8,2) */
    clm[28][1] = -0.8998538011695906;
    clm[28][2] = -0.2261796441029474;

    /** (8,3) */
    clm[29][1] = -1.1292763157894736;
    clm[29][2] = -0.4196774909106648;

    /** (8,4) */
    clm[30][1] = -0.9744152046783625;
    clm[30][2] = -0.47652059150068155;

    /** (8,5) */
    clm[31][1] = -1.1924342105263157;
    clm[31][2] = -0.7220789990670207;

    /** (8,6) */
    clm[32][1] = -1.0986842105263157;
    clm[32][2] = -0.9061610303170207;

    /** (8,7) */
    clm[33][1] = -1.287171052631579;
    clm[33][2] = -1.175404252991305;

    /** (8,8) */
    clm[34][1] = -1.2726608187134503;
    clm[34][2] = -1.5337092502821381;

  }

  /** Compute EulerLogs */
  const double el1 = Eulerlog(x,1);
  const double el2 = Eulerlog(x,2);
  const double el3 = Eulerlog(x,3);
  const double el4 = Eulerlog(x,4);
  const double el5 = Eulerlog(x,5);
  const double el6 = Eulerlog(x,6);
  const double el7 = Eulerlog(x,7);

  const double logx  = log(x);
  const double log2x = logx*logx;
  
  /** Coefs with Eulerlogs */
  /* l = 2 */
  clm[0][3] = (2.9192806270460925  - 1.019047619047619   *el1);
  clm[0][4] = (-1.28235780892213   + 1.073639455782313   *el1);
  clm[0][5] = (-3.8466571723355227 + 0.8486467106683944  *el1);
  clm[0][6] = -0.53614531347435562831 - 1.5600367264240624649*logx + 0.1298072562358276644*log2x;
  
  clm[1][3] = (12.736034731834051 - 4.076190476190476*el2);
  clm[1][4] = (-2.4172313935587004 + 4.173242630385488 *el2);
  clm[1][5] = (-30.14143102836864  + 7.916297736025627 *el2);
  clm[1][6] = 0.49742318178380476016 + 2.0769160997732426304*log2x - 21.360317428096729531*logx;
  
  /* l = 3 */
  clm[2][3] = (1.9098284139598072 - 0.4126984126984127*el1);  
  clm[2][4] = (0.5368150316615179 + 0.2980599647266314*el1);
  clm[2][5] = (1.4497991763035063 - 0.0058477188106817735*el1);
  clm[2][6] = 4.0291559090570805988 + 0.0212899974804736709499*log2x - 0.39311169961214638857*logx;
  
  clm[3][3] = (6.220997955214429 - 1.6507936507936507*el2);
  clm[3][4] = (-3.4527288879001268 + 2.005408583186361*el2);
  clm[3][5] = -4.5270419165221037846 + 0.74342737585535939445*logx;
  clm[3][6] = 2.9819373321032360597 - 4.1746594065141637002*logx + 0.34063995968757873520*log2x;
    
  clm[4][3] = (14.10891386831863 - 3.7142857142857144*el3);
  clm[4][4] = (-6.723375314944128 + 4.333333333333333*el3);
  clm[4][5] = (-29.568699895427518 + 6.302092352092352*el3);
  clm[4][6] = 4.2434215365016546735 - 18.544715779075533606*logx + 1.7244897959183673469*log2x;

  /* l = 4 */
  clm[5][3] = (0.6981550175535535 - 0.2266955266955267*el1);
  clm[5][4] = (-0.7931524512893319 + 0.2584672482399755*el1);
  clm[5][5] = -0.73874769521337814809 + 0.041688165706484035882*logx;
  clm[5][6] = -0.60643478766899956315 - 0.078148517222175236907*logx + 0.0064238577279702820828*log2x;
  
  clm[6][3] = 4.550378418934105e-12*(8.48238724511e11 - 1.9927619712e11*el2);
  clm[6][4] = (-0.6621921297263365 + 0.787251738160829*el2);
  clm[6][5] = 0.066452311537142230164 + 0.16421729454462630378*logx;
  clm[6][6] = 4.7020949396983490140 - 1.4492579183627638456*logx + 0.10278172364752451333*log2x;
  
  clm[7][3] = (8.519456157072423 - 2.0402597402597404*el3);
  clm[7][4] = (-5.353216984886716 + 2.5735094451003544*el3);
  clm[7][5] = -6.0091714045534954041 + 0.99801531141158133484*logx;
  clm[7][6] = 4.9466129199935928986 - 6.3245279024451771064*logx + 0.52033247596559284871*log2x;
  
  clm[8][3] = (15.108111214795123 - 3.627128427128427*el4);
  clm[8][4] = (-8.857121657199649 + 4.434988849534304*el4);
  clm[8][5] = -14.633690582678763747 + 2.9256460798810267523*logx;
  clm[8][6] = 5.3708202812535269509 - 18.028080626090983076*logx + 1.6445075783603922132*log2x;

  /* l = 5 */
  clm[9][3] = (0.642701885362399 - 0.14414918414918415*el1);
  clm[9][4] = (-0.07651588046467575 + 0.11790664036817883*el1);
  clm[9][5] =  0.2112256289378054518 + 0.007552668545567954*logx;
  clm[9][6] =  0.6954698480021733513 - 0.045868829656707944*logx + 0.0025973734113594253*log2x;
  
  clm[10][3] = (2.354458371550237 - 0.5765967365967366*el2);
  clm[10][4] = -0.36261168338728289637 + 0.33429938266861343784*logx;
  clm[10][5] = -1.2210109278350054619 + 0.13346303603712902107*logx;
  clm[10][6] = 0.33075256130927413448 + 0.041557974581750805527*log2x - 0.55630499019814556741*logx;
  
  clm[11][3] = (5.733973288504755 - 1.2973426573426574*el3);
  clm[11][4] = (-1.9573287625526001 + 1.2474448628294783*el3);
  clm[11][5] = -1.2893903635541462354 + 0.37545165958302052977*logx;
  clm[11][6] = 5.5345161119771993654 + 0.21038724632011345298*log2x - 2.758182643809987602*logx;
  
  clm[12][3] = (10.252052781721588 - 2.3063869463869464*el4);
  clm[12][4] = 0.97104909920221130334 + 1.4740600527677450754*logx;
  clm[12][5] = -7.1507394017332213549 + 1.2041810478399218889*logx;
  clm[12][6] = 6.1586857187986271095 + 0.66492759330801288843*log2x - 8.0193253436649521323*logx;
  
  clm[13][3] = (15.939827047208668 - 3.6037296037296036*el5);
  clm[13][4] = (-10.272578060123237 + 4.500041838503377*el5);
  clm[13][5] = -14.676578889547799889 + 2.8378879842489309945*logx;
  clm[13][6] = 5.9456506579632350247 + 1.6233583820996408409*log2x - 17.987106103978598244*logx;

  /* l = 6 */
  clm[14][3] = (0.21653486654395454 - 0.10001110001110002*el1);
  clm[14][4] = -0.41910058747759839184 + 0.05590898299231632565*logx;
  clm[14][5] = -0.54703676004200114476 + 0.014589362666582032864*logx;
  clm[14][6] = -0.67046877793201698851 - 0.010612705664941704748*logx + 0.001250277515678781080*log2x;
  
  clm[15][3] = (1.7942694138754138 - 0.40004440004440006*el2);
  clm[15][4] = 0.11894436114207465266 + 0.17621003335289049575*logx;
  clm[15][5] = -0.020168519020176921191 + 0.049600555197339365552*logx;
  clm[15][6] = 1.3092135879268909484 + 0.020004440250860497281*log2x - 0.29998579925017448295*logx;
  
  clm[16][3] = (4.002558222882566 - 0.9000999000999002*el3);
  clm[16][4] = -0.17301414399537528118 + 0.52818362193362193362*logx;
  clm[16][5] = -1.8598787273050426274 + 0.25227794970506685875*logx;
  clm[16][6] = 1.5038409846301927305 + 0.10127247876998126748*log2x - 1.3388799755686528906*logx;
  
  clm[17][3] = (7.359388663371044 - 1.6001776001776002*el4);
  clm[17][4] = 1.1251294602654682041 + 0.81913853342424770996*logx;
  clm[17][5] = -2.5262527575555087974 + 0.57834034766625984931*logx;
  clm[17][6] = 6.2193866198899176839 + 0.32007104401376795649*log2x - 4.058212416583905376*logx;

  clm[18][3] = (11.623366217471297 - 2.5002775002775004*el5);
  clm[18][4] = 1.1674823841175551524 + 1.6060810331643664977*logx;
  clm[18][5] = -8.0397483562766587382 + 1.3718948493709684487*logx;
  clm[18][6] = 6.9514240915914444298 + 0.78142344729923817503*log2x - 9.3656857835417176404*logx;
  
  clm[19][3] = (16.645950799433503 - 3.6003996003996006*el6);
  clm[19][4] = 2.5867656736437215138 + 2.2716807002521288236*logx;
  clm[19][5] = -14.719871717781358093 + 2.7980705647681528943*logx;
  clm[19][6] = 6.3131827696196368972 + 1.6203596603197002797*log2x - 18.073907939714613453*logx;

  /* l = 7 */
  clm[20][3] = 0.16468579927775824743 - 0.036777788038292239973*logx;
  clm[20][4] = -0.11461516127647669863 + 0.031832875360874795943*logx;
  clm[20][5] = -0.080774412348589137719 + 0.0055469551229555596135*logx;
  clm[20][6] = 0.030761646076864233744 + 0.0006763028464947758845*log2x - 0.0093836267188368350366*logx;

  clm[21][3] = 0.60518781438139242832 - 0.14711115215316895989*logx;
  clm[21][4] = -0.42604721720872623266 + 0.1651443852902587657*logx;
  clm[21][5] = -0.73175760961036224033 + 0.051683014131394955905*logx;
  clm[21][6] = -0.3216279886799197889 + 0.010820845543916414152*log2x - 0.14225943344402801078*logx;
  clm[21][7] = -0.76144280538773684905 - 0.012147290395838407563*log2x + 0.12243085882528399631*logx;
  clm[21][8] = -1.7107570672811356954 - 0.0038015738778090124244*log2x + 0.12635028549960758173*logx;
 
  clm[22][3] = 1.5152673760874757137 - 0.33100009234463015976*logx;
  clm[22][4] = 0.23890581023665524063 + 0.30874798529625166162*logx;
  clm[22][5] = -0.40951574151546501449 + 0.12309038145823335302*logx;
  clm[22][6] = 2.138249656492112722 + 0.054780530566076846645*log2x - 0.77102827271845867626*logx;
  clm[22][7] = 3.8828696441423173767 - 0.051097805822138907374*log2x + 0.17228084445624686514*logx;
  clm[22][8] = 4.1659219244282881166 - 0.020371463814705495917*log2x + 0.23576040054503462879*logx;

  clm[23][3] = 2.3765941182906920513 - 0.58844460861267583957*logx;
  clm[23][4] = 0.034831715172197790382 + 0.69684023412876298568*logx;
  clm[23][5] = -2.5496448149313856708 + 0.38094414588688557689*logx;
  clm[23][6] = 2.6151793026413180884 + 0.17313352870266262643*log2x - 2.2501669045899251545*logx;
 
  clm[24][3] = 3.4549539737325957224 - 0.91944470095730599932*logx;
  clm[24][4] = 1.2445994163412265058 + 0.9812561094250240497*logx;
  clm[24][5] = -3.6121737425863076541 + 0.7603066020101875292*logx;
  clm[24][6] = 6.7427941691635346059 + 0.42268927905923492781*log2x - 5.2559186857515106152*logx;
  
  clm[25][3] = 4.6386119841850458377 - 1.324000369378520639*logx;
  clm[25][4] = 1.3353529096239084184 + 1.7038756254186964286*logx;
  clm[25][5] = -8.7464194960250207018 + 1.5097927099518075214*logx;
  clm[25][6] = 7.4931076636745106598 + 0.87648848905722954631*log2x - 10.453118081371813659*logx;

  clm[26][3] = 5.6637092577775648284 - 1.8021116138763197587*logx;
  clm[26][4] = 2.5475814336400496379 + 2.2867130562632292736*logx;
  clm[26][5] = -14.753691230416999582 + 2.7785800133705155933*logx;
  clm[26][6] = 6.5798142868167883268 + 1.6238031344339568987*log2x - 18.194396598558345462*logx;

  /* l = 8 */
  clm[27][3] = -0.0402600243878044614 - 0.028207774337805297558*logx;
  clm[27][4] = -0.38616527816334081443 + 0.030963599659623775641*logx;
  clm[27][5] = -0.49517911331245971841 + 0.0075715684499400889988*logx;
  clm[27][6] = -0.67742559772045169648 + 0.00039783926654627354656*log2x - 0.00082612129478740295921*logx;

  clm[28][3] = 0.55112429012141994351 - 0.11283109735122119023*logx;
  clm[28][4] = -0.10872490919342399656 + 0.1015314918416325184*logx;
  clm[28][5] = -0.16034999914759608591 + 0.025520097442644222935*logx;
  clm[28][6] = 0.32507998024315148283 + 0.006365428264740376745*log2x - 0.093572230072162669854*logx;

  clm[29][3] = 1.0926607342566844069 - 0.25386996904024767802*logx;
  clm[29][4] = -0.35309773473557329475 + 0.28668934332735864429*logx;
  clm[29][5] = -1.011356410663755538 + 0.10654351162437930413*logx;
  clm[29][6] = 0.27959450210499743463 + 0.032224980590248157272*log2x - 0.43629687195904222137*logx;

  clm[30][3] = 1.9484233096763321564 - 0.45132438940488476092*logx;
  clm[30][4] = 0.37783359676055286399 + 0.43977734727829779701*logx;
  clm[30][5] = -0.90667499119882721697 + 0.21506536499789960964*logx;
  clm[30][6] = 2.9741287296867279309 + 0.10184685223584602792*log2x - 1.3815833067955163488*logx;

  clm[31][3] = 2.7774354272726199981 - 0.70519435844513243894*logx;
  clm[31][4] = 0.23547928350848787334 + 0.84089787808013325367*logx;
  clm[31][5] = -3.2311975154762676935 + 0.50920603649377099228*logx;
  clm[31][6] = 3.5778742849855167012 + 0.2486495415914209666*log2x - 3.1847361552290032316*logx;

  clm[32][3] = 3.7287103355726041884 - 1.0154798761609907121*logx;
  clm[32][4] = 1.3530481267263253384 + 1.115691706045299006*logx;
  clm[32][5] = -4.553227506130683765 + 0.92018829084824385117*logx;
  clm[32][6] = 7.1368206877669610918 + 0.51559968944397051635*log2x - 6.3288803127163295107*logx;

  clm[33][3] = 4.8013970579549365794 - 1.3821809425524595803*logx;
  clm[33][4] = 1.4789424352383168596 + 1.779103298752557348*logx;
  clm[33][5] = -9.3195431866151982863 + 1.6246213582796916264*logx;
  clm[33][6] = 7.8768628245826477659 + 0.9552120789776027853*log2x - 11.346602024237823561*logx;

  clm[34][3] = 5.6960230130299035028 - 1.8052975576195390437*logx;
  clm[34][4] = 2.5376043424259770642 + 2.2975314677014747625*logx;
  clm[34][5] = -14.777949692855076906 + 2.7688015636328382236*logx;
  clm[34][6] = 6.7885227283912752537 + 1.6295496357735364467*log2x - 18.318413973625189326*logx;
 
  /** rho_lm */
  // Taylor - Used at 4PN for lm = 53
  //          Used at 6PN for lm = 42, 51, 61, 77, 88
  //          Used at 8PN for lm = 72
  const int kmaxT4PN = 1;
  int kT4PN[] = {11};

  for (int i=0; i<kmaxT4PN; i++) {
    int k = kT4PN[i];
    rholm[k] = Taylorseries(x,clm[k],4);
  }

  const int kmaxT6PN = 6;
  int kT6PN[] = {2,6,9,14,26,34};

  for (int i=0; i<kmaxT6PN; i++) {
    int k = kT6PN[i];
    rholm[k] = Taylorseries(x,clm[k],6);
  }

  const int kmaxT8PN = 1;
  int kT8PN[] = {21};
  for (int i=0; i<kmaxT8PN; i++) {
    int k = kT8PN[i];
    rholm[k] = Taylorseries(x,clm[k],8);
  }

  // Padé (2,3) : Used for lm = 71
  const int kmaxP23 = 1;
  int kP23[] = {20};

  for (int i=0; i<kmaxP23; i++) {
    int k = kP23[i];
    rholm[k] = Pade23(x,clm[k]);
  }

  // Padé (5,1) : Used for lm = 21, 82
  const int kmaxP51 = 2;
  int kP51[] = {0,28};

  for (int i=0; i<kmaxP51; i++) {
    int k = kP51[i];
    rholm[k] = Pade51(x,clm[k]);
  }

  // Padé (4,2) : Used for lm = 22, 33, 32, 44, 43, 41, 55, 54, 52,
  //                            66, 65, 64, 63, 62, 76, 75, 74, 87, 86, 85, 84
  const int kmaxP42 = 21;
  int kP42[] = {1,3,4,5,7,8,10,12,13,15,16,17,18,19,23,24,25,30,31,32,33};

  for (int i=0; i<kmaxP42; i++) {
    int k = kP42[i];
    rholm[k] = Pade42(x,clm[k]);
  }

  // Padé (3,3) : Used for lm = 81
  const int kmaxP33 = 1;
  int kP33[] = {27};

  for (int i=0; i<kmaxP33; i++) {
    int k = kP33[i];
    rholm[k] = Pade33(x,clm[k]);
  }

  // Padé (1,5) : Used for lm = 83
  const int kmaxP15 = 1;
  int kP15[] = {29};

  for (int i=0; i<kmaxP15; i++) {
    int k = kP15[i];
    rholm[k] = Pade15(x,clm[k]);
  }

  // Padé (6,2) : Used for lm = 73
  const int kmaxP62 = 1;
  int kP62[] = {22};

  for (int i=0; i<kmaxP62; i++) {
    int k = kP62[i];
    rholm[k] = Pade62(x,clm[k]);
  }
  
  if (kmaxT4PN+kmaxT6PN+kmaxT8PN+kmaxP23+kmaxP15+kmaxP42+kmaxP33+kmaxP51+kmaxP62 != KMAX) {
    errorexit("Wrong function: not all multipoles are written.\n");
  }
  
  /** Amplitudes */
#pragma omp simd
  for (int k = 0; k < KMAX; k++) {
      flm[k] = gsl_pow_int(rholm[k], LINDEX[k]);
  }
  
}

/**
 * Function: eob_wav_flm_HM_4PN22
 * ------------------------------
 *   Resummed amplitudes in the general nu-dependent case.
 *   Function introduced for higher modes
 *   (2,2) is taken at 4 PN order and resummed with a Pade' 22
 *   (4,4), (4,2) and (5,5) and l >= 7 modes are Talyor expanded
 *   (3,2) is resummed with a Padé 32
 *   (2,1) is resummed with a Padé 51
 *   All other modes are resummed with a Padé 42
 *   Nagar et al in preparation
 *   
 *   @param[in] x    : frequency
 *   @param[in] nu   : symmetric mass ratio
 *   @param[out] rholm : resummed amplitudes
 *   @param[out] flm   : resummed amplitudes
*/
void eob_wav_flm_HM_4PN22(double x,double nu, double *rholm, double *flm)
{
  /** Coefficients */
  static double clm[KMAX][7];
  
  const double nu2 = nu*nu;
  const double nu3 = nu*nu2;
  const double nu4 = nu*nu3;
  const double Pi2 = SQ(Pi);
  
  //static int firstcall = 1;
  if (EOBPars->firstcall[FIRSTCALL_EOBWAVFLMHM4PN22]) {

    if (0) printf("Precompute some rholm coefs\n");
    EOBPars->firstcall[FIRSTCALL_EOBWAVFLMHM4PN22] = 0;
    
    for (int k=0; k<KMAX; k++) clm[k][0] = 1.;
    for (int k=0; k<KMAX; k++) for (int n=1; n<7; n++) clm[k][n] = 0.;

    /** (2,1) */
    clm[0][1] = (-1.0535714285714286 + 0.27380952380952384 *nu);
    clm[0][2] = (-0.8327841553287982 - 0.7789824263038548  *nu + 0.13116496598639457*nu2);

    /** (2,2) */
    clm[1][1] = (-1.0238095238095237 + 0.6547619047619048*nu);
    clm[1][2] = (-1.94208238851096   - 1.5601379440665155*nu + 0.4625614134542706*nu2);

    /** (3,1) */
    clm[2][1] = (-0.7222222222222222 - 0.2222222222222222*nu);
    clm[2][2] = (0.014169472502805836 - 0.9455667789001122*nu - 0.46520763187429853*nu2);

    /** (3,2) */
    clm[3][1] = (0.003703703703703704*(328. - 1115.*nu + 320.*nu2))/(-1. + 3.*nu);
    clm[3][2] = (6.235191420376606e-7*(-1.444528e6 + 8.050045e6*nu - 4.725605e6*nu2 - 2.033896e7*nu3 + 3.08564e6*nu4))/((-1. + 3.*nu)*(-1. + 3.*nu));

    /** (3,3) */
    clm[4][1] =  (-1.1666666666666667 + 0.6666666666666666*nu);
    clm[4][2] = (-1.6967171717171716 - 1.8797979797979798*nu + 0.45151515151515154*nu2);
    
    /** (4,1) */
    clm[5][1] = (0.001893939393939394*(602. - 1385.*nu + 288.*nu2))/(-1. + 2.*nu);
    clm[5][2] = (- 0.36778992787515513);
    
    /** (4,2) */
    clm[6][1] = (0.0007575757575757576*(1146. - 3530.*nu + 285.*nu2))/(-1. + 3.*nu);
    clm[6][2] = - (3.1534122443213353e-9*(1.14859044e8 - 2.95834536e8*nu - 1.204388696e9*nu2 + 3.04798116e9*nu3 + 3.79526805e8*nu4))/((-1. + 3.*nu)*(-1. + 3.*nu));
        
    /** (4,3) */
    clm[7][1] = (0.005681818181818182*(222. - 547.*nu + 160.*nu2))/(-1. + 2.*nu);
    clm[7][2] = (- 0.9783218202252293);
    
    /** (4,4) */
    clm[8][1] = (0.0007575757575757576*(1614. - 5870.*nu + 2625.*nu2))/(-1. + 3.*nu);
    clm[8][2] = (3.1534122443213353e-9*(-5.11573572e8 + 2.338945704e9*nu - 3.13857376e8*nu2 - 6.733146e9*nu3 + 1.252563795e9*nu4))/((-1. + 3.*nu)*(-1. + 3.*nu));
    
    /** (5,1) */
    clm[9][1] = (0.002564102564102564*(319. - 626.*nu + 8.*nu2))/(-1. + 2.*nu);
    clm[9][2] = (- 0.1047896120973044);
    
    /** (5,2) */
    clm[10][1] = (0.00007326007326007326*(-15828. + 84679.*nu - 104930.*nu2 + 21980.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[10][2] = (- 0.4629337197600934)*PMTERMS_eps; 
    
    /** (5,3) */
    clm[11][1] = (0.002564102564102564*(375. - 850.*nu + 176.*nu2))/(-1. + 2.*nu);
    clm[11][2] = (- 0.5788010707241477);
    
    /** (5,4) */
    clm[12][1] = (0.00007326007326007326*(-17448. + 96019.*nu - 127610.*nu2 + 33320.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[12][2] = (- 1.0442142414362194)*PMTERMS_eps;
    
    /** (5,5) */
    clm[13][1] = (0.002564102564102564*(487. - 1298.*nu + 512.*nu2))/(-1. + 2.*nu);
    clm[13][2] = (- 1.5749727622804546);
        
    /** (6,1) */
    clm[14][1] = (0.006944444444444444*(-161. + 694.*nu - 670.*nu2 + 124.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[14][2] = (- 0.29175486850885135)*PMTERMS_eps;
    
    /** (6,2) */
    clm[15][1] = (0.011904761904761904*(-74. + 378.*nu - 413.*nu2 + 49.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[15][2] = ( - 0.24797525070634313)*PMTERMS_eps;
    
    /** (6,3) */
    clm[16][1] = (0.006944444444444444*(-169. + 742.*nu - 750.*nu2 + 156.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[16][2] = (- 0.5605554442947213)*PMTERMS_eps;
    
    /** (6,4) */
    clm[17][1] = (0.011904761904761904*(-86. + 462.*nu - 581.*nu2 + 133.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[17][2] = (- 0.7228451986855349)*PMTERMS_eps;
    
    /** (6,5) */
    clm[18][1] = (0.006944444444444444*(-185. + 838.*nu - 910.*nu2 + 220.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[18][2] = (- 1.0973940686333457)*PMTERMS_eps;
  
    /** (6,6) */
    clm[19][1] = (0.011904761904761904*(-106. + 602.*nu - 861.*nu2 + 273.*nu3))/(1. - 5.*nu + 5.*nu2); 
    clm[19][2] = (- 1.5543111183867486)*PMTERMS_eps;
  
    /** (7,1) */
    clm[20][1] = (0.0014005602240896359*(-618. + 2518.*nu - 2083.*nu2 + 228.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[20][2] = ( - 0.1508235111143767)*PMTERMS_eps;
  
    /** (7,2) */
    clm[21][1] = (0.00006669334400426837*(16832. - 123489.*nu + 273924.*nu2 - 190239.*nu3 + 32760.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[21][2] = (- 0.351319484450667)*PMTERMS_eps;
    
    /** (7,3) */
    clm[22][1] = (0.0014005602240896359*(-666. + 2806.*nu - 2563.*nu2 + 420.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[22][2] = (- 0.37187416047628863)*PMTERMS_eps;
  
    /** (7,4) */
    clm[23][1] = (0.00006669334400426837*(17756. - 131805.*nu + 298872.*nu2 - 217959.*nu3 + 41076.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[23][2] = (- 0.6473746896670599)*PMTERMS_eps;
    
    /** (7,5) */
    clm[24][1] = (0.0014005602240896359*(-762. + 3382.*nu - 3523.*nu2 + 804.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[24][2] = (- 0.8269193364414116)*PMTERMS_eps;
  
    /** (7,6) */
    clm[25][1] = (0.0006002400960384153*(2144. - 16185.*nu + 37828.*nu2 - 29351.*nu3 + 6104.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[25][2] = (- 1.1403265020692532)*PMTERMS_eps;
    
    /** (7,7) */
    clm[26][1] = (0.0014005602240896359*(-906. + 4246.*nu - 4963.*nu2 + 1380.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[26][2] = (- 1.5418467934923434)*PMTERMS_eps;
  
    /** (8,1) */
    clm[27][1] = (0.00005482456140350877*(20022. - 126451.*nu + 236922.*nu2 - 138430.*nu3 + 21640.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
    clm[27][2] = (- 0.26842133517043704)*PMTERMS_eps;

    /** (8,2) */
    clm[28][1] = (0.0003654970760233918*(2462. - 17598.*nu + 37119.*nu2 - 22845.*nu3 + 3063.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[28][2] = (- 0.2261796441029474)*PMTERMS_eps;

    /** (8,3) */
    clm[29][1] = (0.00005482456140350877*(20598. - 131059.*nu + 249018.*nu2 - 149950.*nu3 + 24520.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
    clm[29][2] = (- 0.4196774909106648)*PMTERMS_eps;

    /** (8,4) */
    clm[30][1] = (0.0003654970760233918*(2666. - 19434.*nu + 42627.*nu2 - 28965.*nu3 + 4899.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[30][2] = (- 0.47652059150068155)*PMTERMS_eps;

    /** (8,5) */
    clm[31][1] = (0.00027412280701754384*(4350. - 28055.*nu + 54642.*nu2 - 34598.*nu3 + 6056.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
    clm[31][2] = (- 0.7220789990670207)*PMTERMS_eps;

    /** (8,6) */
    clm[32][1] = (0.0010964912280701754*(1002. - 7498.*nu + 17269.*nu2 - 13055.*nu3 + 2653.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[32][2] = (- 0.9061610303170207)*PMTERMS_eps;

    /** (8,7) */
    clm[33][1] = (0.00005482456140350877*(23478. - 154099.*nu + 309498.*nu2 - 207550.*nu3 + 38920.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
    clm[33][2] = (- 1.175404252991305)*PMTERMS_eps;

    /** (8,8) */
    clm[34][1] = (0.0003654970760233918*(3482. - 26778.*nu + 64659.*nu2 - 53445.*nu3 + 12243.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[34][2] = (- 1.5337092502821381)*PMTERMS_eps;

  }

  /** Compute EulerLogs */
  const double el1 = Eulerlog(x,1);
  const double el2 = Eulerlog(x,2);
  const double el3 = Eulerlog(x,3);
  const double el4 = Eulerlog(x,4);
  const double el5 = Eulerlog(x,5);
  const double el6 = Eulerlog(x,6);
  const double el7 = Eulerlog(x,7);

  const double logx  = log(x);
  const double log2x = logx*logx;
  
  /** Coefs with Eulerlogs */
  /* l = 2 */
  clm[0][3] = (2.9192806270460925  - 1.019047619047619   *el1);
  clm[0][4] = (-1.28235780892213   + 1.073639455782313   *el1);
  clm[0][5] = (-3.8466571723355227 + 0.8486467106683944  *el1)*PMTERMS_eps;
  clm[0][6] = -0.53614531347435562831 - 1.5600367264240624649*logx + 0.1298072562358276644*log2x;
  
  clm[1][3] = (12.736034731834051  - 2.902228713904598 *nu - 1.9301558466099282*nu2 + 0.2715020968103451*nu3 - 4.076190476190476*el2);
  clm[1][4] = -2.417231393558700   - 50.70249714211636 *nu + 10.48174467376330 *nu2 - 0.1792973174305593*nu3 + 0.1368191977995206*nu4 + (4.173242630385488+19.99773242630385*nu)*el2;

  /* l = 3 */
  clm[2][3] = (1.9098284139598072 - 0.4126984126984127*el1+ (-4.646868015386534 + (0.21354166666666666)*Pi2)*nu + 2.3020866307903347*nu2 - 0.5813492634480288*nu3);  
  clm[2][4] = (0.5368150316615179 + 0.2980599647266314*el1);
  clm[2][5] = (1.4497991763035063 - 0.0058477188106817735*el1)*PMTERMS_eps;
    
  clm[3][3] = (6.220997955214429 - 1.6507936507936507*el2);
  clm[3][4] = (-3.4527288879001268 + 2.005408583186361*el2)*PMTERMS_eps;
  clm[3][5] = -4.5270419165221037846 + 0.74342737585535939445*logx;
  clm[3][6] = 2.9819373321032360597 - 4.1746594065141637002*logx + 0.34063995968757873520*log2x;
    
  clm[4][3] = (14.10891386831863 - 3.7142857142857144*el3 + (-5.031429681429682 + (0.21354166666666666)*Pi2)*nu - 1.7781727531727531*nu2 + 0.25923767590434255*nu3);
  clm[4][4] = (-6.723375314944128 + 4.333333333333333*el3);
  clm[4][5] = (-29.568699895427518 + 6.302092352092352*el3)*PMTERMS_eps;
  clm[4][6] = 4.2434215365016546735 - 18.544715779075533606*logx + 1.7244897959183673469*log2x;

  /* l = 4 */
  clm[5][3] = (0.6981550175535535 - 0.2266955266955267*el1);
  clm[5][4] = (-0.7931524512893319 + 0.2584672482399755*el1)*PMTERMS_eps;
  clm[5][5] = -0.73874769521337814809 + 0.041688165706484035882*logx;
  clm[5][6] = -0.60643478766899956315 - 0.078148517222175236907*logx + 0.0064238577279702820828*log2x;
  
  clm[6][3] = 4.550378418934105e-12*(8.48238724511e11 - 1.9927619712e11*el2);
  clm[6][4] = (-0.6621921297263365 + 0.787251738160829*el2)*PMTERMS_eps;
  clm[6][5] = 0.066452311537142230164 + 0.16421729454462630378*logx;
  clm[6][6] = 4.7020949396983490140 - 1.4492579183627638456*logx + 0.10278172364752451333*log2x;
  
  clm[7][3] = (8.519456157072423 - 2.0402597402597404*el3)*PMTERMS_eps;
  clm[7][4] = (-5.353216984886716 + 2.5735094451003544*el3)*PMTERMS_eps;
  clm[7][5] = -6.0091714045534954041 + 0.99801531141158133484*logx;
  clm[7][6] = 4.9466129199935928986 - 6.3245279024451771064*logx + 0.52033247596559284871*log2x;
  
  clm[8][3] = (15.108111214795123 - 3.627128427128427*el4);
  clm[8][4] = (-8.857121657199649 + 4.434988849534304*el4)*PMTERMS_eps;
  clm[8][5] = -14.633690582678763747 + 2.9256460798810267523*logx;
  clm[8][6] = 5.3708202812535269509 - 18.028080626090983076*logx + 1.6445075783603922132*log2x;

  /* l = 5 */
  clm[9][3] = (0.642701885362399 - 0.14414918414918415*el1)*PMTERMS_eps;
  clm[9][4] = (-0.07651588046467575 + 0.11790664036817883*el1)*PMTERMS_eps;
  clm[9][5] =  0.2112256289378054518 + 0.007552668545567954*logx;
  clm[9][6] =  0.6954698480021733513 - 0.045868829656707944*logx + 0.0025973734113594253*log2x;
  
  clm[10][3] = (2.354458371550237 - 0.5765967365967366*el2)*PMTERMS_eps;
  clm[10][4] = -0.36261168338728289637 + 0.33429938266861343784*logx;
  clm[10][5] = -1.2210109278350054619 + 0.13346303603712902107*logx;
  clm[10][6] = 0.33075256130927413448 + 0.041557974581750805527*log2x - 0.55630499019814556741*logx;
  
  clm[11][3] = (5.733973288504755 - 1.2973426573426574*el3)*PMTERMS_eps;
  clm[11][4] = (-1.9573287625526001 + 1.2474448628294783*el3)*PMTERMS_eps;
  clm[11][5] = -1.2893903635541462354 + 0.37545165958302052977*logx;
  clm[11][6] = 5.5345161119771993654 + 0.21038724632011345298*log2x - 2.758182643809987602*logx;
  
  clm[12][3] = (10.252052781721588 - 2.3063869463869464*el4)*PMTERMS_eps;
  clm[12][4] = 0.97104909920221130334 + 1.4740600527677450754*logx;
  clm[12][5] = -7.1507394017332213549 + 1.2041810478399218889*logx;
  clm[12][6] = 6.1586857187986271095 + 0.66492759330801288843*log2x - 8.0193253436649521323*logx;
  
  clm[13][3] = (15.939827047208668 - 3.6037296037296036*el5)*PMTERMS_eps;
  clm[13][4] = (-10.272578060123237 + 4.500041838503377*el5)*PMTERMS_eps;
  clm[13][5] = -14.676578889547799889 + 2.8378879842489309945*logx;
  clm[13][6] = 5.9456506579632350247 + 1.6233583820996408409*log2x - 17.987106103978598244*logx;

  /* l = 6 */
  clm[14][3] = (0.21653486654395454 - 0.10001110001110002*el1)*PMTERMS_eps;
  clm[14][4] = -0.41910058747759839184 + 0.05590898299231632565*logx;
  clm[14][5] = -0.54703676004200114476 + 0.014589362666582032864*logx;
  clm[14][6] = -0.67046877793201698851 - 0.010612705664941704748*logx + 0.001250277515678781080*log2x;
  
  clm[15][3] = (1.7942694138754138 - 0.40004440004440006*el2)*PMTERMS_eps;
  clm[15][4] = 0.11894436114207465266 + 0.17621003335289049575*logx;
  clm[15][5] = -0.020168519020176921191 + 0.049600555197339365552*logx;
  clm[15][6] = 1.3092135879268909484 + 0.020004440250860497281*log2x - 0.29998579925017448295*logx;
  
  clm[16][3] = (4.002558222882566 - 0.9000999000999002*el3)*PMTERMS_eps;
  clm[16][4] = -0.17301414399537528118 + 0.52818362193362193362*logx;
  clm[16][5] = -1.8598787273050426274 + 0.25227794970506685875*logx;
  clm[16][6] = 1.5038409846301927305 + 0.10127247876998126748*log2x - 1.3388799755686528906*logx;
  
  clm[17][3] = (7.359388663371044 - 1.6001776001776002*el4)*PMTERMS_eps;
  clm[17][4] = 1.1251294602654682041 + 0.81913853342424770996*logx;
  clm[17][5] = -2.5262527575555087974 + 0.57834034766625984931*logx;
  clm[17][6] = 6.2193866198899176839 + 0.32007104401376795649*log2x - 4.058212416583905376*logx;

  clm[18][3] = (11.623366217471297 - 2.5002775002775004*el5)*PMTERMS_eps;
  clm[18][4] = 1.1674823841175551524 + 1.6060810331643664977*logx;
  clm[18][5] = -8.0397483562766587382 + 1.3718948493709684487*logx;
  clm[18][6] = 6.9514240915914444298 + 0.78142344729923817503*log2x - 9.3656857835417176404*logx;
  
  clm[19][3] = (16.645950799433503 - 3.6003996003996006*el6)*PMTERMS_eps;
  clm[19][4] = 2.5867656736437215138 + 2.2716807002521288236*logx;
  clm[19][5] = -14.719871717781358093 + 2.7980705647681528943*logx;
  clm[19][6] = 6.3131827696196368972 + 1.6203596603197002797*log2x - 18.073907939714613453*logx;

  /* l = 7 */
  clm[20][3] = (0.2581280702019663 - 0.07355557607658449*el1)*PMTERMS_eps;
  
  clm[22][3] = (3.0835293524055283 - 0.6620001846892604*el3)*PMTERMS_eps;
  
  clm[24][3] = (8.750589067052443 - 1.838889401914612*el5)*PMTERMS_eps;
  
  clm[26][3] = (17.255875091408523 - 3.6042232277526396*el7)*PMTERMS_eps;

  /** rho_lm */
  const double x2  = x*x;
  const double x3  = x*x2;
  const double x4  = x*x3;
  const double x5  = x*x4;
  const double x6  = x*x5;
  const double xn[] = {1.,x,x2,x3,x4,x5,x6};

  /** Initializing Padé approximants */
  double cden,n1,n2,n3,n4,n5,d1,d2 = 0.;
  double k = 0.;

  // Padé (2,2) : Used for (2,2) multipole
  const int kmaxPade22 = 1;
  int kPade22[] = {1};

  for (int i=0; i<kmaxPade22; i++) {
    int k = kPade22[i];
    rholm[k] = Pade22(x,clm[k]);
  }

  // Padé (3,2) : Used for (3,1) multipole
  const int kmaxPade32 = 1;
  int kPade32[] = {2};
  
  for (int i=0; i<kmaxPade32; i++) {
    int k = kPade32[i];
    rholm[k] = Pade32(x,clm[k]);
  }

  // Padé (4,2) - Used for (3,2), (3,3), (4,1), (4,3), (5,2), (5,3), (5,4) and l=6 except (6,1) 
  const int kmaxPade42 = 12;
  int kPade42[] = {3,4,5,7,10,11,12,15,16,17,18,19};
  
  for (int i=0; i<kmaxPade42; i++) {
    int k = kPade42[i];
    rholm[k] = Pade42(x,clm[k]);
  }

  // Padé (5,1) - Used for (2,1)
  const int kmaxPade51 = 1;
  int kPade51[] = {0};
  
  for (int i=0; i<kmaxPade51; i++) {
    int k = kPade51[i];
    rholm[k] = Pade51(x,clm[k]);
  }

  // Taylor series : (4,4), (4,2), (5,5), (5,1), (6,1) and l>6  at 6PN
  const int kmaxTaylor = 20;
  int kTaylor[kmaxTaylor];
  kTaylor[0] = 6;
  kTaylor[1] = 8;
  kTaylor[2] = 9;
  kTaylor[3] = 13;
  kTaylor[4] = 14;
  for (int i=5; i<kmaxTaylor; i++) {
    kTaylor[i] = 15+i;
  }
  
  for (int i=0; i<kmaxTaylor; i++) {
    int k = kTaylor[i];
    rholm[k] = Taylorseries(x,clm[k],6);
  }
  
  if (kmaxTaylor+kmaxPade32+kmaxPade22+kmaxPade42+kmaxPade51 != KMAX) {
    errorexit("Wrong function: not all multipoles are written.\n");
  }
  
  /** Amplitudes */
#pragma omp simd
  for (int k = 0; k < KMAX; k++) {
      flm[k] = gsl_pow_int(rholm[k], LINDEX[k]);
  }

}

/** 
 * Function: eob_wav_flm_s_SSLO
 *   Resummed amplitudes for the spin case. 
 *   This function computes the residual amplitude corrections flm's as 
 *   introduced in Damour, Iyer & Nagar, PRD 79, 064004 (2008).
 *   The orbital part is taken at the usual 3^{+2} PN order, i.e. 3PN terms
 *   are integrated by the 4PN and 5PN test-particle terms, with the higher
 *   modes obtained by Fujita & Iyer.
 *   It only includes spin-spin interaction at LO for the (2,2) mode.
 *   Note that the variables called here (a1,a2)
 *   are what we usually call tilde{a}_1 and tilde{a}_2 and are defined as
 *   a1 = X1*chi1, a2=X2*chi2 and are passed here as parameters. Special 
 *   combinations of these quantities are used here to write the spin-dependent
 *   part of the waveform in particularly compact form, so that the (spinning)
 *   test-particle limit is recovered just by visual inspection of the equations
 * 
 *   @param[in] x: x = (M omega)^{2/3}
 *   @param[in] nu: symmetric mass ratio
 *   @param[in] X1: mass fraction of body 1
 *   @param[in] X2: mass fraction of body 2
 *   @param[in] chi1: dimensionless spin of body 1
 *   @param[in] chi2: dimensionless spin of body 2
 *   @param[in] a1  : a1 = X1*chi1
 *   @param[in] a2  : a2 = X2*chi2
 *   @param[in] C_Q1: spin-induced quadrupole of body 1
 *   @param[in] C_Q2: spin-induced quadrupole of body 2
 *   @param[in] usetidal: flag for tidal effects
 *   @param[out] rholm: residual amplitude corrections
 *   @param[out] flm: residual amplitude corrections
 */
void eob_wav_flm_s_SSLO(double x, double nu, double X1, double X2, double chi1, double chi2, double a1, double a2, double C_Q1, double C_Q2, int usetidal,
			double *rholm, double *flm)
{

  /** Orbital part */
  //double rholm_orb[KMAX], flm_orb[KMAX];
  eob_wav_flm(x,nu, rholm, flm);

  /** Spin corrections */
  double rho22S;
  double rho32S;
  double rho44S;
  double rho42S;
  double f21S;
  double f33S;
  double f31S;
  double f43S;
  double f41S;
      
  const double a0      = a1+a2;
  const double a12     = a1-a2;
  const double X12     = X1-X2;
  const double a0X12   = a0*X12;
  const double a12X12  = a12*X12;
  
  const double v  = sqrt(x);
  const double v2 = x;
  const double v3 = v*v2;
  const double v4 = v3*v;
  const double v5 = v4*v;
     
  /** l=m=2 multipole */
  /* spin-orbit */
  const double cSO_lo    = (-0.5*a0 - a12X12/6.);
  const double cSO_nlo   = (-52./63.-19./504.*nu)*a0 - (50./63.+209./504.*nu)*a12X12;
  
  /* SPIN-SPIN contribution */
  double cSS_lo;
  if (usetidal) {
#if (EXCLUDESPINSPINTIDES)
    /* Switch off spin-spin-tidal couplings */
    /* See also: eob_dyn_s_get_rc() */
    cSS_lo = 0.;
    /* Above code switch off everything, 
       Alt. one can set C_Q1=C_Q2=0, but keep the term: */
    /*
      cSS_lo = a1*a2;
    */
#else
    cSS_lo = 0.5*(C_Q1*a1*a1 + 2.*a1*a2 + C_Q2*a2*a2);
#endif
  } else {
    cSS_lo = 0.5*a0*a0; 
  }

  /* rho_22^S: Eq. (80) of Damour & Nagar, PRD 90, 044018 (2014) */
  rho22S = cSO_lo*v3 + cSS_lo*v4 + cSO_nlo*v5 ;
    
  /** l>=3, m=even: multipoles rewritten in compact and self-explanatory form */
  rho32S = (a0-a12X12)/(3.*(1.-3.*nu))*v;
  rho44S = (-19./30.*a0 -  (1.-21.*nu)/(30.-90.*nu)*a12X12)*v3;
  rho42S = ( -1./30.*a0 - (19.-39.*nu)/(30.-90.*nu)*a12X12)*v3;
  
  /** l>=2, m=odd: multipoles rewritten in compact and self-explanatory form */
  f21S = -1.5*a12*v + ((110./21. + 79./84.*nu)*a12 - 13./84.*a0X12)*v3;
  f33S = ((-0.25 + 2.5*nu)*a12 - 1.75*a0X12)*v3;
  f31S = ((-2.25 + 6.5*nu)*a12 + 0.25*a0X12)*v3;
  f43S = (( 5. -10.*nu)*a12 - 5.*a0X12)/(-4.+8.*nu)*v;
  f41S = f43S;
    
  /** Amplitudes (correct with spin terms) */
  flm[0] = gsl_pow_int(rholm[0], 2);
  flm[0] = (X12*flm[0] + f21S);
  
  flm[1] = gsl_pow_int(rholm[1]+ rho22S, 2);
  
  flm[2] = gsl_pow_int(rholm[2], 3);
  flm[2] = (X12*flm[2] + f31S);
  
  flm[3] = gsl_pow_int(rholm[3]+ rho32S, 3);
  
  flm[4] = gsl_pow_int(rholm[4], 3);
  flm[4] = (X12*flm[4] + f33S);
  
  flm[5] = gsl_pow_int(rholm[5], 4);
  flm[5] = (X12*flm[5] + f41S);
  
  flm[6] = gsl_pow_int(rholm[6] + rho42S, 4);
  
  flm[7] = gsl_pow_int(rholm[7], 4);
  flm[7] = (X12*flm[7] + f43S);
  
  flm[8] = gsl_pow_int(rholm[8] + rho44S, 4);
      
}


/**
 * Function: eob_wav_flm_s_SSLO
 * ----------------------------
 *   Resummed amplitudes for the spin case. 
 *   This function computes the residual amplitude corrections flm's as 
 *   introduced in Damour, Iyer & Nagar, PRD 79, 064004 (2008).
 *   The orbital part is taken at the usual 3^{+2} PN order, i.e. 3PN terms
 *   are integrated by the 4PN and 5PN test-particle terms, with the higher
 *   modes obtained by Fujita & Iyer.
 *   The function includes spin-spin interaction at NLO for the (2,2) mode
 *   and at LO for the (2,1),(3,1) and (3,3) modes. 
 *   Note that the variables called here (a1,a2)
 *   are what we usually cal tilde{a}_1 and tilde{a}_2 and are defined as
 *   a1 = X1*chi1, a2=X2*chi2 and are passed here as parameters.
 * 
 *   @param[in] x: x = (M omega)^{2/3}
 *   @param[in] nu: symmetric mass ratio
 *   @param[in] X1: mass fraction of body 1
 *   @param[in] X2: mass fraction of body 2
 *   @param[in] chi1: dimensionless spin of body 1
 *   @param[in] chi2: dimensionless spin of body 2
 *   @param[in] a1  : a1 = X1*chi1
 *   @param[in] a2  : a2 = X2*chi2
 *   @param[in] C_Q1: spin-induced quadrupole of body 1
 *   @param[in] C_Q2: spin-induced quadrupole of body 2
 *   @param[in] usetidal: flag for tidal effects
 *   @param[out] rholm: residual amplitude corrections
 *   @param[out] flm: residual amplitude corrections
*/
void eob_wav_flm_s_SSNLO(double x, double nu, double X1, double X2, double chi1, double chi2, double a1, double a2, double C_Q1, double C_Q2, int usetidal,
			 double *rholm, double *flm)
{

  /** Orbital part */
  //double rholm_orb[KMAX], flm_orb[KMAX];
  eob_wav_flm(x,nu, rholm, flm);

  /** Spin corrections */
  double rho22S;
  double rho32S;
  double rho44S;
  double rho42S;
  double f21S;
  double f33S;
  double f31S;
  double f43S;
  double f41S;
      
  const double a0      = a1+a2;
  const double a12     = a1-a2;
  const double X12     = X1-X2;
  const double a0X12   = a0*X12;
  const double a12X12  = a12*X12;
  
  const double v  = sqrt(x);
  const double v2 = x;
  const double v3 = v2*v;
  const double v4 = v3*v;
  const double v5 = v4*v;
  const double v6 = v5*v;
  const double v7 = v6*v;
     
  /** l=m=2 multipole */
  /* spin-orbit */
  const double cSO_lo    = (-0.5*a0 - a12X12/6.);
  const double cSO_nlo   = (-52./63.-19./504.*nu)*a0 - (50./63.+209./504.*nu)*a12X12;
  const double cSO_nnlo  = (32873./21168 + 477563./42336.*nu + 147421./84672.*nu*nu)*a0 - (23687./63504 - 171791./127008.*nu + 50803./254016.*nu*nu)*a12X12;
  
  /* SPIN-SPIN contribution */
  double cSS_lo = 0.;
  double cSS_nlo = 0.;
  if (usetidal) {
#if (EXCLUDESPINSPINTIDES)
    /* Switch off spin-spin-tidal couplings */
    /* See also: eob_dyn_s_get_rc() */
    cSS_lo  = 0.;
    cSS_nlo = 0.;
    /* Above code switch off everything, 
       Alt. one can set C_Q1=C_Q2=0, but keep the term: */
    /*
      cSS_lo = a1*a2;
    */
#else
    cSS_lo  = 0.5*(C_Q1*a1*a1 + 2.*a1*a2 + C_Q2*a2*a2);
    cSS_nlo = (-85./63. + 383./252.*nu)*a1*a2 + (-2./3. - 5./18.*nu)*(a1*a1 + a2*a2) + (1./7. + 27./56.*nu)*(C_Q1*a1*a1 + C_Q2*a2*a2) + 2./9.*X12*(a1*a1 - a2*a2) + 55./84.*X12*(C_Q1*a1*a1 - C_Q2*a2*a2);
#endif
  } else {
    cSS_lo  = 0.5*a0*a0;
    cSS_nlo = 1./504.*(2.*(19. - 70.*nu)*a12*a12 + (-302. + 243.*nu)*a0*a0 + 442.*X12*a0*a12);
  }
    
  /* rho_22^S: Eq. (80) of Damour & Nagar, PRD 90, 044018 (2014) */
  rho22S = cSO_lo*v3 + cSS_lo*v4 + cSO_nlo*v5;

  // Adding NLO SS term w.r.t. eob_wav_flm_s_SSLO
  rho22S += cSS_nlo*v6;

  /** l>=3, m=even: multipoles rewritten in compact and self-explanatory form */
  rho32S = (a0-a12X12)/(3.*(1.-3.*nu))*v;
  rho44S = (-19./30.*a0 -  (1.-21.*nu)/(30.-90.*nu)*a12X12)*v3;
  rho42S = ( -1./30.*a0 - (19.-39.*nu)/(30.-90.*nu)*a12X12)*v3;
  
  /** l>=2, m=odd*/
  /* spin-orbit */
  f21S = -1.5*a12*v + ((110./21. + 79./84.*nu)*a12 - 13./84.*a0X12)*v3;
  f33S = ((-0.25 + 2.5*nu)*a12 - 1.75*a0X12)*v3;
  f31S = ((-2.25 + 6.5*nu)*a12 + 0.25*a0X12)*v3;
  f43S = (( 5. -10.*nu)*a12 - 5.*a0X12)/(-4.+8.*nu)*v;
  f41S = f43S;
  
  /* SPIN-SPIN contribution */
  double c21SS_lo;
  double c33SS_lo;
  double c31SS_lo;
  if (usetidal) {
#if (EXCLUDESPINSPINTIDES)
    /* Switch off spin-spin-tidal couplings */
    /* See also: eob_dyn_s_get_rc() */
    c21SS_lo  = 0.;
    c33SS_lo  = 0.;
    c31SS_lo  = 0.;
    /* Above code switch off everything, 
       Alt. one can set C_Q1=C_Q2=0, but keep the term: */
#else
    c21SS_lo  = -19./8.*(a1*a1 - a2*a2) - (C_Q1*a1*a1 - C_Q2*a2*a2) + 1./8.*(-9.*a1*a1 + 10*a1*a2 -9.*a2*a2 + 12.*(C_Q1*a1*a1 + C_Q2*a2*a2))*X12;
    c33SS_lo  = 3.*(a1*a2 + 0.5*(C_Q1*a1*a1 + C_Q2*a2*a2))*X12;
    c31SS_lo  = -4.*(C_Q1*a1*a1 - C_Q2*a2*a2) + 3.*(a1*a2 + 0.5*(C_Q1*a1*a1 + C_Q2*a2*a2))*X12;
#endif
  } else {
    c21SS_lo  = 1./8.*(-27.*(a1*a1 - a2*a2) + (3.*a1*a1 + 10.*a1*a2 + 3.*a2*a2)*X12);
    c33SS_lo  = 3./2.*a0*a0*X12;
    c31SS_lo  = -4.*(a1*a1 - a2*a2) + 3./2.*a0*a0*X12;
  }

  // Adding LO SS term w.r.t. eob_wav_flm_s_SSLO
  f21S += c21SS_lo*v4;
  f33S += c33SS_lo*v4;
  f31S += c31SS_lo*v4;
  
  /** Amplitudes (correct with spin terms) */
  flm[0] = gsl_pow_int(rholm[0], 2);
  flm[0] = (X12*flm[0] + f21S);
  
  flm[1] = gsl_pow_int(rholm[1]+ rho22S, 2);
  
  flm[2] = gsl_pow_int(rholm[2], 3);
  flm[2] = (X12*flm[2] + f31S);
  
  flm[3] = gsl_pow_int(rholm[3]+ rho32S, 3);
  
  flm[4] = gsl_pow_int(rholm[4], 3);
  flm[4] = (X12*flm[4] + f33S);
  
  flm[5] = gsl_pow_int(rholm[5], 4);
  flm[5] = (X12*flm[5] + f41S);
  
  flm[6] = gsl_pow_int(rholm[6] + rho42S, 4);
  
  flm[7] = gsl_pow_int(rholm[7], 4);
  flm[7] = (X12*flm[7] + f43S);
  
  flm[8] = gsl_pow_int(rholm[8] + rho44S, 4);
      
}

/**
 * Function: eob_wav_flm_s_HM
 * --------------------------
 *    Resummed amplitudes for the spin case. 
 *    This function computes the residual amplitude corrections flm's as 
 *    introduced in Damour, Iyer & Nagar, PRD 79, 064004 (2008).
 *    The m=even modes have the usual structure (rho_lm^orb + rho_lm^spin)^l.
 *    The m=odd modes are factorized as f_lm^orb * f_lm^spin.
 *    Some of the f_lm^spin are inverse-resummed
 *    Function introduced for higher modes.
 *    Ref: arXiv:2001.09082
 *    Note that the variables called here (a1,a2)
 *    are what we usually cal tilde{a}_1 and tilde{a}_2 and are defined as
 *    a1 = X1*chi1, a2=X2*chi2 and are passed here as parameters.
 * 
 *   @param[in] x: x = (M omega)^{2/3}
 *   @param[in] nu: symmetric mass ratio
 *   @param[in] X1: mass fraction of body 1
 *   @param[in] X2: mass fraction of body 2
 *   @param[in] chi1: dimensionless spin of body 1
 *   @param[in] chi2: dimensionless spin of body 2
 *   @param[in] a1  : a1 = X1*chi1
 *   @param[in] a2  : a2 = X2*chi2
 *   @param[in] C_Q1: spin-induced quadrupole of body 1
 *   @param[in] C_Q2: spin-induced quadrupole of body 2
 *   @param[in] usetidal: flag for tidal effects
 *   @param[out] rholm: residual amplitude corrections
 *   @param[out] flm: residual amplitude corrections
 */
void eob_wav_flm_s_HM(double x, double nu, double X1, double X2, double chi1, double chi2, double a1, double a2,
		      double C_Q1, double C_Q2, int usetidal, double *rholm, double *flm)
{
  /** Orbital part */
  //double rholm_orb[KMAX], flm_orb[KMAX];
  eob_wav_flm(x,nu, rholm, flm);

  /** Spin corrections */
  double rho22S;
  double rho32S;
  double rho44S;
  double rho42S;
  double f21S;
  double f33S;
  double f31S;
  double f43S;
  double f41S;
  double f55S;

  double nu2 = nu*nu;

  const double el1 = Eulerlog(x,1);
  const double el2 = Eulerlog(x,2);
  const double el3 = Eulerlog(x,3);
  const double el4 = Eulerlog(x,4);
  const double el5 = Eulerlog(x,5);
  const double el6 = Eulerlog(x,6);
  const double el7 = Eulerlog(x,7);
  
  const double a0      = a1+a2;
  const double a12     = a1-a2;
  const double X12     = X1-X2;
  const double a0X12   = a0*X12;
  const double a12X12  = a12*X12;
  
  const double v  = sqrt(x);
  const double v2 = x;
  const double v3 = v2*v;
  const double v4 = v3*v;
  const double v5 = v4*v;
  const double v6 = v5*v;
  const double v7 = v6*v;
  const double v8 = v7*v;
  const double v9 = v8*v;

  double logx = log(x);
     
  /** l=m=2 multipole */
  /* spin-orbit */
  const double cSO_lo    = (-0.5*a0 - a12X12/6.);
  const double cSO_nlo   = (-52./63.-19./504.*nu)*a0 - (50./63.+209./504.*nu)*a12X12;
  const double cSO_nnlo  = (32873./21168 + 477563./42336.*nu + 147421./84672.*nu*nu)*a0 - (23687./63504 - 171791./127008.*nu + 50803./254016.*nu*nu)*a12X12; // Not used for the moment

  
  /* SPIN-SPIN contribution */
  double cSS_lo = 0.;
  double cSS_nlo = 0.;
  if (usetidal) {
#if (EXCLUDESPINSPINTIDES)
    /* Switch off spin-spin-tidal couplings */
    /* See also: eob_dyn_s_get_rc() */
    cSS_lo  = 0.;
    cSS_nlo = 0.;
    /* Above code switch off everything, 
       Alt. one can set C_Q1=C_Q2=0, but keep the term: */
    /*
      cSS_lo = a1*a2;
    */
#else
    cSS_lo  = 0.5*(C_Q1*a1*a1 + 2.*a1*a2 + C_Q2*a2*a2);
    cSS_nlo = (-85./63. + 383./252.*nu)*a1*a2 + (-2./3. - 5./18.*nu)*(a1*a1 + a2*a2) + (1./7. + 27./56.*nu)*(C_Q1*a1*a1 + C_Q2*a2*a2) + 2./9.*X12*(a1*a1 - a2*a2) + 55./84.*X12*(C_Q1*a1*a1 - C_Q2*a2*a2);
#endif
  } else {
    cSS_lo  = 0.5*a0*a0;
    cSS_nlo = 1./504.*(2.*(19. - 70.*nu)*a12*a12 + (-302. + 243.*nu)*a0*a0 + 442.*X12*a0*a12);
  }
  
  /* New cubic-in-spin term */
  double cS3_lo = (7./12.*a0 - 0.25*X12*a12)*a0*a0;
  
  /* rho_22^S: Eq. (80) of Damour & Nagar, PRD 90, 044018 (2014) */
  rho22S = cSO_lo*v3 + cSS_lo*v4 + cSO_nlo*v5 + cSS_nlo*v6 + cS3_lo*v7; //+ cSO_nnlo*v7;
  
  /* l=3, m=2 */
  /* spin-orbit coefficients */
  double c32_SO_lo   = 0.;
  double c32_SO_nlo  = 0.;
  double c32_SO_nnlo = 0.;
  c32_SO_lo   = (a0-a12X12)/(3.*(1.-3.*nu));
  c32_SO_nlo  = (-(0.884567901234568 - 3.41358024691358*nu + 2.4598765432098766*nu2)*a0 + (1.10679012345679 - 2.635802469135802*nu - 1.8734567901234567*nu2)*X12*a12)/SQ(1-3*nu);
  c32_SO_nnlo = -1.335993162073409*a0 - 0.9586570436879079*a12*X12;
  
  rho32S = c32_SO_lo*v + c32_SO_nlo*v3 + c32_SO_nnlo*v5;
  
  /* l=4, m=4 */
  /* spin-orbit coefficients */
  double c44_SO_lo   = 0.;
  double c44_SO_nlo  = 0.;
  double c44_SO_nnlo = 0.;
  c44_SO_lo   = -19./30.*a0 - (1. - 21.*nu)/(30. - 90.*nu)*a12*X12;
  c44_SO_nlo  = -199./550.*a0 - 491./550.*a12*X12;
  c44_SO_nnlo = 1.9942241584173401*a0 + 0.012143034995307722*a12*X12;
  
  rho44S = c44_SO_lo*v3 + c44_SO_nlo*v5 + c44_SO_nnlo*v7;
  
  /* l=4, m=2 */
  /* spin-orbit coefficients */
  double c42_SO_lo   = 0.;
  double c42_SO_nlo  = 0.;
  double c42_SO_nnlo = 0.;
  double c42_SO_n3lo = 0.;
  c42_SO_lo   = -1./30.*a0 - (19.-39.*nu)/(30.-90.*nu)*a12X12;
  c42_SO_nlo  = -219./550.*a0 + 92./275.*a12*X12;
  c42_SO_nnlo = -1.245162901492446947*a0 + 0.6414503261889625526*a12*X12;
  c42_SO_n3lo = -(2.6105076622404808654 + 0.16912938912938912939*logx)*a0 + -(1.3519220653268411505 - 0.10290524290524290524*logx)*X12*a12;
  
  rho42S = c42_SO_lo*v3 + c42_SO_nlo*v5 + c42_SO_nnlo*v7 + c42_SO_n3lo*v9;
  
  /** l>=2, m=odd*/
  double if210s = 1. + 13./84.*a0*v3 - 1./8.*(3.*a1+a2)*(a1+3.*a2)*v4 + a0*(14705./7056. - 12743./7056.*nu)*v5;
  double if211s = 1. - 9./4.*a0*v3 + (349./252. + 74./63.*nu)*v2 + (65969./31752. + 89477./31752.*nu + 46967./31752.*nu2 - 0.5*a0*a0)*v4;
  f21S = X12/if210s - 1.5*v*a12/if211s;

  
  double if330s = 1. + 7./4.*a0*v3 - 1.5*a0*a0*v4 + + 1./60.*a0*(211. - 127.*nu)*v5;
  double f331s = (10.*nu -1. + (-169. + 671.*nu + 182.*nu2)/15.*x);
  f33S = X12/if330s + 0.25*a12*v3*f331s;

  double if310s = 1. - 0.25*a0*v3 - 1.5*a0*a0*v4 + 1./36.*a0*(13. - 449.*nu)*v5;
  double f311s  = 26.*nu - 9. - 16.*a0*v + (9. - 95.*nu + 66.*nu2)/9.*v2;
  f31S = X12/if310s + 0.25*a12*v3*f311s;

  double f430s = 1. - 1.25/(2.*nu - 1.)*a0*v;
  double f431s = 1.;            
  f43S = X12*f430s - 1.25*a12*v*f431s;

  double f410s = 1. - 1.25/(2.*nu - 1.)*a0*v;
  double f411s = 1.;
  f41S   = X12*f410s - 1.25*a12*v*f411s;

  double if550s = 1. + 10./3.*a0*v3 - 2.5*a0*a0*v4;
  double f551s  = 1.;
  f55S   = X12/if550s + 10.*nu*(1. - 3.*nu)/(3. - 6.*nu)*a12*v3*f551s;
	    
  /** Amplitudes (correct with spin terms) */
  flm[0] = gsl_pow_int(rholm[0], 2);
  flm[0] = flm[0]*f21S;
  
  flm[1] = gsl_pow_int(rholm[1]+ rho22S, 2);
  
  flm[2] = gsl_pow_int(rholm[2], 3);
  flm[2] = flm[2]*f31S;
  
  flm[3] = gsl_pow_int(rholm[3]+ rho32S, 3);
  
  flm[4] = gsl_pow_int(rholm[4], 3);
  flm[4] = flm[4]*f33S;
  
  flm[5] = gsl_pow_int(rholm[5], 4);
  flm[5] = flm[5]*f41S;
  
  flm[6] = gsl_pow_int(rholm[6] + rho42S, 4);
  
  flm[7] = gsl_pow_int(rholm[7], 4);
  flm[7] = flm[7]*f43S;
  
  flm[8] = gsl_pow_int(rholm[8] + rho44S, 4);

  flm[13] = gsl_pow_int(rholm[13], 5);
  flm[13] = flm[13]*f55S;
}

/**
 * Function: eob_wav_flm_s_HM
 * --------------------------
 *    Resummed amplitudes for the spin case. 
 *    This function computes the residual amplitude corrections flm's as 
 *    introduced in Damour, Iyer & Nagar, PRD 79, 064004 (2008).
 *    The m=even modes have the usual structure (rho_lm^orb + rho_lm^spin)^l.
 *    The m=odd modes are factorized as f_lm^orb * f_lm^spin.
 *    Some of the f_lm^spin are inverse-resummed
 *    Function introduced for higher modes when using the 4PN22 orbital part.
 *    Ref: Nagar et al in preparation
 *    Note that the variables called here (a1,a2)
 *    are what we usually cal tilde{a}_1 and tilde{a}_2 and are defined as
 *    a1 = X1*chi1, a2=X2*chi2 and are passed here as parameters.
 * 
 *   @param[in] x: x = (M omega)^{2/3}
 *   @param[in] nu: symmetric mass ratio
 *   @param[in] X1: mass fraction of body 1
 *   @param[in] X2: mass fraction of body 2
 *   @param[in] chi1: dimensionless spin of body 1
 *   @param[in] chi2: dimensionless spin of body 2
 *   @param[in] a1  : a1 = X1*chi1
 *   @param[in] a2  : a2 = X2*chi2
 *   @param[in] C_Q1: spin-induced quadrupole of body 1
 *   @param[in] C_Q2: spin-induced quadrupole of body 2
 *   @param[in] usetidal: flag for tidal effects
 *   @param[out] rholm: residual amplitude corrections
 *   @param[out] flm: residual amplitude corrections
 */
void eob_wav_flm_s_HM_4PN22(double x, double nu, double X1, double X2, double chi1, double chi2, double a1, double a2,
		      double C_Q1, double C_Q2, int usetidal, double *rholm, double *flm)
{

/** Orbital part */
  //double rholm_orb[KMAX], flm_orb[KMAX];
  eob_wav_flm(x,nu, rholm, flm);

  /** Spin corrections */
  double rho22S;
  double rho32S;
  double rho44S;
  double rho42S;
  double f21S;
  double f33S;
  double f31S;
  double f43S;
  double f41S;
  double f55S;

  double nu2 = nu*nu;

  const double el1 = Eulerlog(x,1);
  const double el2 = Eulerlog(x,2);
  const double el3 = Eulerlog(x,3);
  const double el4 = Eulerlog(x,4);
  const double el5 = Eulerlog(x,5);
  const double el6 = Eulerlog(x,6);
  const double el7 = Eulerlog(x,7);
  
  const double a0      = a1+a2;
  const double a12     = a1-a2;
  const double X12     = X1-X2;
  const double a0X12   = a0*X12;
  const double a12X12  = a12*X12;
  
  const double v  = sqrt(x);
  const double v2 = x;
  const double v3 = v2*v;
  const double v4 = v3*v;
  const double v5 = v4*v;
  const double v6 = v5*v;
  const double v7 = v6*v;
  const double v8 = v7*v;
  const double v9 = v8*v;
  const double v11 = v9*v2;

  double logx = log(x);
     
  /** l=m=2 multipole */
  /* spin-orbit */
  const double c_a_n3lo     = -8494939./467775.   + 2536./315.*el2;
  const double c_sigma_n3lo = -14661629./8731800. + 214./315.*el2;
  const double c_a_n4lo     = -890245226581./26698351680. + 328./6615.*el2;
  const double c_sigma_n4lo = -90273995723./88994505600.  + 428./6615.*el2;

  const double cSO_lo    = (-0.5*a0 - a12X12/6.);
  const double cSO_nlo   = (-52./63.-19./504.*nu)*a0 - (50./63.+209./504.*nu)*a12X12;
  const double cSO_nnlo  = (32873./21168 + 477563./42336.*nu + 147421./84672.*nu*nu)*a0 - (23687./63504 - 171791./127008.*nu + 50803./254016.*nu*nu)*a12X12;
  const double cSO_n3lo  = a0*0.5*(c_a_n3lo + c_sigma_n3lo) + a12*X12*0.5*(c_a_n3lo - c_sigma_n3lo);      
  const double cSO_n4lo  = a0*0.5*(c_a_n4lo + c_sigma_n4lo) + a12*X12*0.5*(c_a_n4lo - c_sigma_n4lo);      
  
  /* SPIN-SPIN contribution */
  double cSS_lo = 0.;
  double cSS_nlo = 0.;
  if (usetidal) {
#if (EXCLUDESPINSPINTIDES)
    /* Switch off spin-spin-tidal couplings */
    /* See also: eob_dyn_s_get_rc() */
    cSS_lo  = 0.;
    cSS_nlo = 0.;
    /* Above code switch off everything, 
       Alt. one can set C_Q1=C_Q2=0, but keep the term: */
    /*
      cSS_lo = a1*a2;
    */
#else
    cSS_lo  = 0.5*(C_Q1*a1*a1 + 2.*a1*a2 + C_Q2*a2*a2);
    cSS_nlo = (-85./63. + 383./252.*nu)*a1*a2 + (-2./3. - 5./18.*nu)*(a1*a1 + a2*a2) + (1./7. + 27./56.*nu)*(C_Q1*a1*a1 + C_Q2*a2*a2) + 2./9.*X12*(a1*a1 - a2*a2) + 55./84.*X12*(C_Q1*a1*a1 - C_Q2*a2*a2);
#endif
  } else {
    cSS_lo  = 0.5*a0*a0;
    cSS_nlo = 1./504.*(2.*(19. - 70.*nu)*a12*a12 + (-302. + 243.*nu)*a0*a0 + 442.*X12*a0*a12);
  }
  
  /* New cubic-in-spin term */
  double cS3_lo = (7./12.*a0 - 0.25*X12*a12)*a0*a0;
  
  /* rho_22^S: Eq. (80) of Damour & Nagar, PRD 90, 044018 (2014) */
  rho22S = cSO_lo*v3 + cSS_lo*v4 + cSO_nlo*v5 + cSS_nlo*v6 + cS3_lo*v7 + cSO_nnlo*v7 + cSO_n3lo*v9 + cSO_n4lo*v11;
  
  /* l=3, m=2 */
  /* spin-orbit coefficients */
  double c32_SO_lo   = 0.;
  double c32_SO_nlo  = 0.;
  double c32_SO_nnlo = 0.;
  c32_SO_lo   = (a0-a12X12)/(3.*(1.-3.*nu));
  c32_SO_nlo  = (-(0.884567901234568 - 3.41358024691358*nu + 2.4598765432098766*nu2)*a0 + (1.10679012345679 - 2.635802469135802*nu - 1.8734567901234567*nu2)*X12*a12)/SQ(1-3*nu);
  c32_SO_nnlo = -1.335993162073409*a0 - 0.9586570436879079*a12*X12;
  
  rho32S = c32_SO_lo*v + c32_SO_nlo*v3 + c32_SO_nnlo*v5;
  
  /* l=4, m=4 */
  /* spin-orbit coefficients */
  double c44_SO_lo   = 0.;
  double c44_SO_nlo  = 0.;
  double c44_SO_nnlo = 0.;
  c44_SO_lo   = -19./30.*a0 - (1. - 21.*nu)/(30. - 90.*nu)*a12*X12;
  c44_SO_nlo  = -199./550.*a0 - 491./550.*a12*X12;
  c44_SO_nnlo = 1.9942241584173401*a0 + 0.012143034995307722*a12*X12;
  
  rho44S = c44_SO_lo*v3 + c44_SO_nlo*v5 + c44_SO_nnlo*v7;
  
  /* l=4, m=2 */
  /* spin-orbit coefficients */
  double c42_SO_lo   = 0.;
  double c42_SO_nlo  = 0.;
  double c42_SO_nnlo = 0.;
  double c42_SO_n3lo = 0.;
  c42_SO_lo   = -1./30.*a0 - (19.-39.*nu)/(30.-90.*nu)*a12X12;
  c42_SO_nlo  = -219./550.*a0 + 92./275.*a12*X12;
  c42_SO_nnlo = -1.245162901492446947*a0 + 0.6414503261889625526*a12*X12;
  c42_SO_n3lo = -(2.6105076622404808654 + 0.16912938912938912939*logx)*a0 + -(1.3519220653268411505 - 0.10290524290524290524*logx)*X12*a12;
  
  rho42S = c42_SO_lo*v3 + c42_SO_nlo*v5 + c42_SO_nnlo*v7 + c42_SO_n3lo*v9;
  
  /** l>=2, m=odd*/
  double if210s = 1. + 13./84.*a0*v3 - 1./8.*(3.*a1+a2)*(a1+3.*a2)*v4 + a0*(14705./7056. - 12743./7056.*nu)*v5;
  double if211s = 1. - 9./4.*a0*v3 + (349./252. + 74./63.*nu)*v2 + (65969./31752. + 89477./31752.*nu + 46967./31752.*nu2 - 0.5*a0*a0)*v4;
  f21S = X12/if210s - 1.5*v*a12/if211s;

  
  double if330s = 1. + 7./4.*a0*v3 - 1.5*a0*a0*v4 + + 1./60.*a0*(211. - 127.*nu)*v5;
  double f331s = (10.*nu -1. + (-169. + 671.*nu + 182.*nu2)/15.*x);
  f33S = X12/if330s + 0.25*a12*v3*f331s;

  double if310s = 1. - 0.25*a0*v3 - 1.5*a0*a0*v4 + 1./36.*a0*(13. - 449.*nu)*v5;
  double f311s  = 26.*nu - 9. - 16.*a0*v + (9. - 95.*nu + 66.*nu2)/9.*v2;
  f31S = X12/if310s + 0.25*a12*v3*f311s;

  double f430s = 1. - 1.25/(2.*nu - 1.)*a0*v;
  double f431s = 1.;            
  f43S = X12*f430s - 1.25*a12*v*f431s;

  double f410s = 1. - 1.25/(2.*nu - 1.)*a0*v;
  double f411s = 1.;
  f41S   = X12*f410s - 1.25*a12*v*f411s;

  double if550s = 1. + 10./3.*a0*v3 - 2.5*a0*a0*v4;
  double f551s  = 1.;
  f55S   = X12/if550s + 10.*nu*(1. - 3.*nu)/(3. - 6.*nu)*a12*v3*f551s;
	    
  /** Amplitudes (correct with spin terms) */
  flm[0] = gsl_pow_int(rholm[0], 2);
  flm[0] = flm[0]*f21S;
  
  flm[1] = gsl_pow_int(rholm[1]+ rho22S, 2);
  
  flm[2] = gsl_pow_int(rholm[2], 3);
  flm[2] = flm[2]*f31S;
  
  flm[3] = gsl_pow_int(rholm[3]+ rho32S, 3);
  
  flm[4] = gsl_pow_int(rholm[4], 3);
  flm[4] = flm[4]*f33S;
  
  flm[5] = gsl_pow_int(rholm[5], 4);
  flm[5] = flm[5]*f41S;
  
  flm[6] = gsl_pow_int(rholm[6] + rho42S, 4);
  
  flm[7] = gsl_pow_int(rholm[7], 4);
  flm[7] = flm[7]*f43S;
  
  flm[8] = gsl_pow_int(rholm[8] + rho44S, 4);

  flm[13] = gsl_pow_int(rholm[13], 5);
  flm[13] = flm[13]*f55S;

}


/**
 * Function: eob_wav_flm_s_Kerr
 * ----------------------------
 *   Resummed amplitudes for a particle orbiting around Kerr 
 *   Adding inverse-resummed factorized spinning term
 * 
 *   @param[in] x: x = (M omega)^{2/3}
 *   @param[in] nu: symmetric mass ratio
 *   @param[in] X1: mass fraction of body 1
 *   @param[in] X2: mass fraction of body 2
 *   @param[in] chi1: dimensionless spin of body 1
 *   @param[in] chi2: dimensionless spin of body 2
 *   @param[in] a1  : a1 = X1*chi1
 *   @param[in] a2  : a2 = X2*chi2
 *   @param[in] C_Q1: spin-induced quadrupole of body 1
 *   @param[in] C_Q2: spin-induced quadrupole of body 2
 *   @param[in] usetidal: flag for tidal effects
 *   @param[out] rholm: residual amplitude corrections
 *   @param[out] flm: residual amplitude corrections
*/
void eob_wav_flm_s_Kerr(double x, double nu, double X1, double X2, double chi1, double chi2, double a1, double a2,
		      double C_Q1, double C_Q2, int usetidal, double *rholm, double *flm)
{
  /** Orbital part */
  eob_wav_flm(x,nu, rholm, flm);

  /** Spin corrections */
  //choosing spin variable as \tilde{a} = \tilde{a}_1
  //a different choice could be \tilde{a}_0
  const double ta   = a1;
  const double ta2  = ta*ta;
  const double ta3  = ta2*ta;
  const double ta4  = ta3*ta;
  const double ta5  = ta4*ta;
  const double ta6  = ta5*ta;
  const double ta7  = ta6*ta;
  const double ta8  = ta7*ta;
  const double ta9  = ta8*ta;
  const double ta10 = ta9*ta;
  const double ta11 = ta10*ta;
  const double ta12 = ta11*ta;

  /** Coefficients of the inverse-resummed series*/
  static double clms[KMAX][17];
  
  //static int firstcall = 1;
  if (EOBPars->firstcall[FIRSTCALL_EOBWAVFLMKERRS]) {
    if (0) printf("Precompute some rholm spin coefs\n");
    EOBPars->firstcall[FIRSTCALL_EOBWAVFLMKERRS] = 0;
    
    for (int k=0; k<KMAX; k++) clms[k][0] = 1.;
    for (int k=0; k<KMAX; k++) for (int n=1; n<17; n++) clms[k][n] = 0.;
    
    /* These are series in v = sqrt(x) */
    /** (2,1) */
    clms[0][1] = 0.75*ta;
    clms[0][2] = 0.84375*ta2;
    clms[0][3] = -0.96130952380952380952*ta + 1.0546875*ta3;
    clms[0][4] = -0.66294642857142857143*ta2 + 1.38427734375*ta4;
    clms[0][5] = 0.92230017006802721088*ta - 1.0555245535714285714*ta3 + 1.8687744140625*ta5;
    clms[0][6] = -1.6378569657029478458*ta2 - 1.61279296875*ta4 + 2.5695648193359375*ta6;
    clms[0][7] = 1.3114143975591679673*ta + 0.21708253613945578231*ta3 - 2.484283447265625*ta5 + 3.579036712646484375*ta7;
    clms[0][8] = -0.1940010182915157405*ta2 - 0.050144740513392857143*ta4 - 3.84198760986328125*ta6 + 5.0330203771591186523*ta8;
    
    /** (2,2) */
    clms[1][1] = 0.;
    clms[1][2] = 0.;
    clms[1][3] = 0.66666666666666666667*ta;
    clms[1][4] = -0.5*ta2;
    clms[1][5] = 2.3015873015873015873*ta;
    clms[1][6] = -0.42063492063492063492*ta2;
    clms[1][7] = 2.4711514235323759133*ta - ta3;
    clms[1][8] = 0.3450491307634164777*ta2 + 0.375*ta4;

    /** (3,1) */
    clms[2][1] = 0.;
    clms[2][2] = 0.;
    clms[2][3] = 0.66666666666666666667*ta;
    clms[2][4] = 0.83333333333333333333*ta2;
    clms[2][5] = 0.037037037037037037037*ta;
    clms[2][6] = 1.5*ta2;
    clms[2][7] = 0.49971941638608305275*ta + ta3;
    clms[2][8] = 0.28072390572390572391*ta2 + 1.9305555555555555556*ta4;
    
    /** (3,2) */
    clms[3][1] = 0.;
    clms[3][2] = 0.;
    clms[3][3] = -0.22222222222222222222*ta;
    clms[3][4] = -0.33333333333333333333*ta2;
    clms[3][5] = 2.024691358024691358*ta;
    clms[3][6] = -1.5604938271604938272*ta2;
    
    /** (3,3) */
    clms[4][1] = 0.;
    clms[4][2] = 0.;
    clms[4][3] = 0.66666666666666666667*ta;
    clms[4][4] = -0.5*ta2;
    clms[4][5] = 2.1111111111111111111*ta;
    clms[4][6] = -0.27777777777777777778*ta2;
    clms[4][7] = 1.8106060606060606061*ta - ta3;
    clms[4][8] = 1.4739898989898989899*ta2 + 0.375*ta4;
    
    /** (4,1) */
    clms[5][1] = 0.;
    clms[5][2] = 0.;
    clms[5][3] = 0.;
    clms[5][4] = -0.375*ta2;
    clms[5][5] = 0.3613455988455988456*ta + 0.83333333333333333333*ta3;
    clms[5][6] = -2.2035984848484848485*ta2;
    
    /** (4,2) */
    clms[6][1] = 0.;
    clms[6][2] = 0.;
    clms[6][3] = 0.66666666666666666667*ta;
    clms[6][4] = -0.5*ta2;
    clms[6][5] = 0.64242424242424242424*ta;
    clms[6][6] = -0.57626262626262626263*ta2;
    clms[6][7] = 1.4029188992825356462*ta - ta3;
    clms[6][8] = 0.73493203765931038658*ta2 + 0.375*ta4;
 
    /** (4,3) */
    clms[7][1] = 0.;
    clms[7][2] = 0.;
    clms[7][3] = 0.;
    clms[7][4] = -0.375*ta2;
    clms[7][5] = 1.9663961038961038961*ta;
    clms[7][6] = -1.3323863636363636364*ta2;
 
    /** (4,4) */
    clms[8][1] = 0.;
    clms[8][2] = 0.;
    clms[8][3] = 0.66666666666666666667*ta;
    clms[8][4] = -0.5*ta2;
    clms[8][5] = 2.069696969696969697*ta;
    clms[8][6] = -0.22171717171717171717*ta2;
    clms[8][7] = 1.5997759815941634123*ta - ta3;
    clms[8][8] = 1.6391290527654164018*ta2 + 0.375*ta4;
 
    /** (5,1) */
    clms[9][1] = 0.;
    clms[9][2] = 0.;
    clms[9][3] = 0.66666666666666666667*ta;
    clms[9][4] = -0.5*ta2;
    clms[9][5] = 0.40273504273504273504*ta;
    clms[9][6] = -0.62649572649572649573*ta2;
    clms[9][7] = 0.99677780908550139319*ta - ta3;
    clms[9][8] = -2.0113928806236498544*ta2 + 1.975*ta4;
    
    /** (5,2) */
    clms[10][1] = 0.;
    clms[10][2] = 0.;
    clms[10][3] = 0.13333333333333333333*ta;
    clms[10][4] = -0.4*ta2;
    clms[10][5] = 0.74344322344322344322*ta;
    clms[10][6] = -1.7020268620268620269*ta2;
    
    /** (5,3) */
    clms[11][1] = 0.;
    clms[11][2] = 0.;
    clms[11][3] = 0.66666666666666666667*ta;
    clms[11][4] = -0.5*ta2;
    clms[11][5] = 0.95794871794871794872*ta;
    clms[11][6] = -0.48290598290598290598*ta2;
    clms[11][7] = 1.6612980182210951442*ta - ta3;
    clms[11][8] = 0.77681036911806142575*ta2 + 0.375*ta4;
    
    /** (5,4) */
    clms[12][1] = 0.;
    clms[12][2] = 0.;
    clms[12][3] = 0.13333333333333333333*ta;
    clms[12][4] = -0.4*ta2;
    clms[12][5] = 1.9909157509157509158*ta;
    clms[12][6] = -1.1547741147741147741*ta2;

    /** (5,5) */
    clms[13][1] = 0.;
    clms[13][2] = 0.;
    clms[13][3] = 0.66666666666666666667*ta;
    clms[13][4] = -0.5*ta2;
    clms[13][5] = 2.0683760683760683761*ta;
    clms[13][6] = -0.1957264957264957265*ta2;
    clms[13][7] = 1.5248364171441094518*ta - ta3;
    clms[13][8] = 1.6970320278012585705*ta2 + 0.375*ta4;
 
    /** (6,1) */
    clms[14][1] = 0.;
    clms[14][2] = 0.;
    clms[14][3] = 0.22222222222222222222*ta;
    clms[14][4] = -0.41666666666666666667*ta2;
    clms[14][5] = 0.4749037999037999038*ta;
    clms[14][6] = -1.6037808641975308642*ta2;

    /** (6,2) */
    clms[15][1] = 0.;
    clms[15][2] = 0.;
    clms[15][3] = 0.66666666666666666667*ta;
    clms[15][4] = -0.5*ta2;
    clms[15][5] = 0.60634920634920634921*ta;
    clms[15][6] = -0.56349206349206349206*ta2;
    
    /** (6,3) */
    clms[16][1] = 0.;
    clms[16][2] = 0.;
    clms[16][3] = 0.22222222222222222222*ta;
    clms[16][4] = -0.41666666666666666667*ta2;
    clms[16][5] = 0.99130591630591630592*ta;
    clms[16][6] = -1.4093364197530864198*ta2;
    
    /** (6,4) */
    clms[17][1] = 0.;
    clms[17][2] = 0.;
    clms[17][3] = 0.66666666666666666667*ta;
    clms[17][4] = -0.5*ta2;
    clms[17][5] = 1.1587301587301587302*ta;
    clms[17][6] = -0.42063492063492063492*ta2;
    clms[17][7] = 1.7486994797919167667*ta - ta3;
    clms[17][8] = 0.92409691149186947506*ta2 + 0.375*ta4;
    
    /** (6,5) */
    clms[18][1] = 0.;
    clms[18][2] = 0.;
    clms[18][3] = 0.22222222222222222222*ta;
    clms[18][4] = -0.41666666666666666667*ta2;
    clms[18][5] = 2.0241101491101491101*ta;
    clms[18][6] = -1.0204475308641975309*ta2;

    /** (6,6) */
    clms[19][1] = 0.;
    clms[19][2] = 0.;
    clms[19][3] = 0.66666666666666666667*ta;
    clms[19][4] = -0.5*ta2;
    clms[19][5] = 2.0793650793650793651*ta;
    clms[19][6] = -0.18253968253968253968*ta2;
    clms[19][7] = 1.5019098548510313216*ta - ta3;
    clms[19][8] = 1.7235985303212193968*ta2 + 0.375*ta4;

    /** (7,1) */
    clms[20][1] = 0.;
    clms[20][2] = 0.;
    clms[20][3] = 0.66666666666666666667*ta;
    clms[20][4] = -0.5*ta2;
    clms[20][5] = 0.50920368147258903561*ta;
    clms[20][6] = -0.57889822595704948646*ta2;
    clms[20][7] = 1.0465146866514689186*ta - ta3;
    clms[20][8] = -1.6079243474544551632*ta2 + 0.375*ta4;
    
    /** (7,2) */
    clms[21][1] = 0.;
    clms[21][2] = 0.;
    clms[21][3] = 0.28571428571428571429*ta;
    clms[21][4] = -0.4285714285714285714*ta2;
    clms[21][5] = 0.6373148171467499199*ta;
    clms[21][6] = -1.393517406962785114*ta2;
    
    /** (7,3) */
    clms[22][1] = 0.;
    clms[22][2] = 0.;
    clms[22][3] = 0.66666666666666666667*ta;
    clms[22][4] = -0.5*ta2;
    clms[22][5] = 0.77330932372949179672*ta;
    clms[22][6] = -0.51167133520074696545*ta2;
    clms[22][7] = 1.4026251780742896074*ta - ta3;
    clms[22][8] = -0.50048216124441506235*ta2 + 0.375*ta4;

    /** (7,4) */
    clms[23][1] = 0.;
    clms[23][2] = 0.;
    clms[23][3] = 0.28571428571428571429*ta;
    clms[23][4] = -0.42857142857142857143*ta2;
    clms[23][5] = 1.1687540595103620314*ta;
    clms[23][6] = -1.2148859543817527011*ta2;
 
    /** (7,5) */
    clms[24][1] = 0.;
    clms[24][2] = 0.;
    clms[24][3] = 0.66666666666666666667*ta;
    clms[24][4] = -0.5*ta2;
    clms[24][5] = 1.3015206082432973189*ta;
    clms[24][6] = -0.37721755368814192344*ta2;
    clms[24][7] = 1.7784672337040497927*ta - ta3;
    clms[24][8] = 1.0544815246347949212*ta2 + 0.375*ta4;
    
    /** (7,6) */
    clms[25][1] = 0.;
    clms[25][2] = 0.;
    clms[25][3] = 0.28571428571428571429*ta;
    clms[25][4] = -0.42857142857142857143*ta2;
    clms[25][5] = 2.0544861301163822172*ta;
    clms[25][6] = -0.91716686674669867947*ta2;
 
    /** (7,7) */
    clms[26][1] = 0.;
    clms[26][2] = 0.;
    clms[26][3] = 0.66666666666666666667*ta;
    clms[26][4] = -0.5*ta2;
    clms[26][5] = 2.0938375350140056022*ta;
    clms[26][6] = -0.17553688141923436041*ta2;
    clms[26][7] = 1.5012829991089870893*ta - ta3;
    clms[26][8] = 1.7371253371014343515*ta2 + 0.375*ta4;

    /** (8,1) */
    clms[27][1] = 0.;
    clms[27][2] = 0.;
    clms[27][3] = 0.33333333333333333333*ta;
    clms[27][4] = -0.4375*ta2;
    clms[27][5] = 0.52597539521881627145*ta;
    clms[27][6] = -1.3146786236424394319*ta2;
 
    /** (8,2) */
    clms[28][1] = 0.;
    clms[28][2] = 0.;
    clms[28][3] = 0.66666666666666666667*ta;
    clms[28][4] = -0.5*ta2;
    clms[28][5] = 0.61170286828181565024*ta;
    clms[28][6] = -0.54459064327485380117*ta2;
    clms[28][7] = 1.1615662118509148364*ta - ta3;
    clms[28][8] = -1.1838784065700039846*ta2 + 0.375*ta4;

    /** (8,3) */
    clms[29][1] = 0.;
    clms[29][2] = 0.;
    clms[29][3] = 0.33333333333333333333*ta;
    clms[29][4] = -0.4375*ta2;
    clms[29][5] = 0.78506060825139772508*ta;
    clms[29][6] = -1.234885390559732665*ta2;
    
    /** (8,4) */
    clms[30][1] = 0.;
    clms[30][2] = 0.;
    clms[30][3] = 0.66666666666666666667*ta;
    clms[30][4] = -0.5*ta2;
    clms[30][5] = 0.91109718741297688666*ta;
    clms[30][6] = -0.47002923976608187135*ta2;
    clms[30][7] = 1.5136986515872326891*ta - ta3;
    clms[30][8] = -0.096799224338467182419*ta2 + 0.375*ta4;
 
    /** (8,5) */
    clms[31][1] = 0.;
    clms[31][2] = 0.;
    clms[31][3] = 0.33333333333333333333*ta;
    clms[31][4] = -0.4375*ta2;
    clms[31][5] = 1.3032310343165606324*ta;
    clms[31][6] = -1.0752989243943191312*ta2;
 
    /** (8,6) */
    clms[32][1] = 0.;
    clms[32][2] = 0.;
    clms[32][3] = 0.66666666666666666667*ta;
    clms[32][4] = -0.5*ta2;
    clms[32][5] = 1.410087719298245614*ta;
    clms[32][6] = -0.34576023391812865497*ta2;
    clms[32][7] = 1.7869038599648017931*ta - ta3;
    clms[32][8] = 1.1578169829554871106*ta2 + 0.375*ta4;
 
    /** (8,7) */
    clms[33][1] = 0.;
    clms[33][2] = 0.;
    clms[33][3] = 0.33333333333333333333*ta;
    clms[33][4] = -0.4375*ta2;
    clms[33][5] = 2.0804866734143049933*ta;
    clms[33][6] = -0.83591922514619883041*ta2;

    /** (8,8) */
    clms[34][1] = 0.;
    clms[34][2] = 0.;
    clms[34][3] = 0.66666666666666666667*ta;
    clms[34][4] = -0.5*ta2;
    clms[34][5] = 2.1086744639376218324*ta;
    clms[34][6] = -0.17178362573099415205*ta2;
    clms[34][7] = 1.5106585502091811726*ta - ta3;
    clms[34][8] = 1.7441965706739483286*ta2 + 0.375*ta4;
  }

  /** Coefs with logs */
  double logx  = log(x);
  double log2x = logx*logx;

  /** (2,1) */
  clms[0][9] = 2.1432276762345683892*ta - 2.9637925001030204049*ta3 + 0.16360146658761160714*ta5 - 5.9490709304809570313*ta7 + 7.1301122009754180908*ta9 - 0.16349206349206349206*ta*logx;
  clms[0][10] = 10.1604098863899707794*ta10 + 0.53420903707856251647*ta2 + 1.70142648592442816773*ta4 + 0.60125255584716796875*ta6 - 9.2106116499219621931*ta8 - 0.117857142857142857143*ta2*logx;
  clms[0][11] = 3.2244044759352393094*ta + 14.547859610058367252*ta11 - 3.4105249991378461458*ta3 + 0.034083771572926933558*ta5 + 1.4600903647286551339*ta7 - 14.249239635254655566*ta9 - 0.0065381708238851095994*ta*logx - 0.12723214285714285714*ta3*logx;
  clms[0][12] = -0.0089759790102565521099*ta - 22.020712419013891901*ta10 + 20.912548189458902925*ta12 + 1.1671393292614780147*ta2 - 3.4104369378224115654*ta4 + 0.60554137456984747024*ta6 + 3.0681325337108300657*ta8 + 0.033838813303099017385*ta2*logx - 0.15234375*ta4*logx;
    
  /** (2,2) */
  clms[1][9] = 6.197427210971866706*ta - 3.9744268077601410935*ta3 - 2.6666666666666666667*ta*logx;
  clms[1][10] = 2.8533526606088964366*ta2 - 0.011904761904761904762*ta4;
  clms[1][11] = 27.840967036880135309*ta - 9.8558411018728479046*ta3 + 1.25*ta5 + 0.084656084656084656085*ta*logx;
  clms[1][12] = 7.627431430434046489*ta2 - 2.8063009154278995549*ta4 - 0.3125*ta6 + 0.50793650793650793651*ta2*logx;
    
  /** (3,1) */
  clms[2][9] = 3.0338047663005778393*ta - 1.7345679012345679012*ta3 + 0.44814814814814814815*ta*logx;
  clms[2][10] = -1.1299102480134886064*ta2 + 6.7314814814814814815*ta4 - 0.71111111111111111111*ta2*logx;
  clms[2][11] = 7.0978695345383967876*ta - 7.5534200024940765682*ta3 + 3.8796296296296296296*ta5 + 0.97833894500561167228*ta*logx;
  clms[2][12] = 1.3404128655316451151*ta - 2.7950944724296102379*ta2 + 3.1447109988776655443*ta4 + 4.3603395061728395062*ta6 - 1.1911335578002244669*ta2*logx;
    
  /** (3,2) */
  clms[3][7] = 4.2843103702525058965*ta + 0.2962962962962962963*ta3 + 0.71111111111111111111*ta*logx;
  clms[3][8] = -2.2875071704701334331*ta2 + 0.22222222222222222222*ta4;
  clms[3][9] = -1.1981029652358133921*ta - 6.576314586191129401*ta3 - 0.69099887766554433221*ta*logx;
  clms[3][10] = 6.7020643276582255754*ta + 10.16886776640903892*ta2 + 3.37448559670781893*ta4 - 0.98765432098765432099*ta2*logx;
  clms[3][11] = 3.5324513835421456474*ta - 8.9360857702109674338*ta2 - 19.63998330979916144*ta3 - 0.34567901234567901235*ta5 - 1.3954719359702524686*ta*logx - 0.94814814814814814815*ta3*logx;
  clms[3][12] = -2.7417535885874559172*ta + 40.020756001133736386*ta2 - 2.6307151995492187536*ta4 - 0.17283950617283950617*ta6 + 3.4807731637361266991*ta2*logx;
    
  /** (3,3) */
  clms[4][9] = 3.0630805494318437582*ta - 3.7839506172839506173*ta3 - 2.1*ta*logx;
  clms[4][10] = 5.711309369642702976*ta2 - 0.083333333333333333333*ta4;
  clms[4][11] = 15.699495272936482402*ta - 7.3698559670781893004*ta3 + 1.25*ta5 + 0.074545454545454545455*ta*logx;
  clms[4][12] = 13.389478502644973202*ta2 - 5.3829732510288065844*ta4 - 0.3125*ta6 + 0.31818181818181818182*ta2*logx;
    
  /** (4,1) */
  clms[5][7] = 0.24996650483782792129*ta + 1.7171717171717171717*ta3 + 0.051767676767676767677*ta*logx;
  clms[5][8] = -2.2993855313457586185*ta2 + 1.9800347222222222222*ta4;
  clms[5][9] = -0.14878438410793557782*ta - 1.4993503918219827311*ta3 - 1.9791666666666666667*ta5 + 0.012560335287608014881*ta*logx;
  clms[5][10] = 0.11926416555294585442*ta - 2.5504108026385059532*ta2 + 7.4710554202741702742*ta4 + 1.7361111111111111111*ta6 + 0.36624579124579124579*ta2*logx;
  clms[5][11] = -0.78412366076737686749*ta - 0.14908020694118231803*ta2 - 1.4161189856680101673*ta3 - 6.5873842592592592593*ta5 + 0.072625781386641738411*ta*logx - 0.74289772727272727273*ta3*logx;
  clms[5][12] = 0.046702435207935493752*ta - 3.9207187992494699478*ta2 + 7.3714562254921597414*ta4 + 0.63828059238215488215*ta6 + 0.90173875121791788458*ta2*logx + 0.21569865319865319865*ta4*logx;
    
  /** (4,2) */
  clms[6][9] = 4.0268821778988419428*ta - 1.8152637485970819304*ta3 + 0.36848484848484848485*ta*logx;
  clms[6][10] = 10.306452573767249751*ta2 + 0.065909090909090909091*ta4 + 1.1151515151515151515*ta2*logx;
  clms[6][11] = 7.2092959679688420526*ta - 0.95702406459982217558*ta3 + 1.25*ta5 + 0.91744070071342798616*ta*logx;
  clms[6][12] = 1.8735316188680948768*ta + 9.2704964875929844336*ta2 - 3.4381464658611123258*ta4 - 0.3125*ta6 + 1.3407360316451225542*ta2*logx;
    
  /** (4,3) */
  clms[7][7] = 2.8428435865102365785*ta + 0.46590909090909090909*ta*logx;
  clms[7][8] = -0.5103470866254957164*ta2 + 0.2578125*ta4;
  clms[7][9] = 0.86096433988270195774*ta - 5.3356128246753246753*ta3 - 0.84624920534011443102*ta*logx;
  clms[7][10] = 3.2201324699295380694*ta + 5.0215868579408468653*ta2 + 2.770478219696969697*ta4 - 0.39469696969696969697*ta2*logx;
  clms[7][11] = 7.8154393981818957165*ta - 4.0251655874119225868*ta2 - 11.16678948744902902*ta3 - 0.9288359748960511834*ta*logx - 0.640625*ta3*logx;
  clms[7][12] = -1.531838124189827296*ta + 21.753915537816650071*ta2 - 4.222309464328285351*ta4 - 0.2041015625*ta6 + 2.3050859651711924439*ta2*logx;
    
  /** (4,4) */
  clms[8][9] = 2.4433734619111654609*ta - 3.7258698092031425365*ta3 - 2.0557575757575757576*ta*logx;
  clms[8][10] = 6.2831750535965411998*ta2 - 0.11136363636363636364*ta4;
  clms[8][11] = 13.449423970553870329*ta - 7.0318198580319792441*ta3 + 1.25*ta5 + 0.060613971886699159426*ta*logx;
  clms[8][12] = 15.321564778193561136*ta2 - 5.1337335050845151855*ta4 - 0.3125*ta6 + 0.46418874055237691601*ta2*logx;
    
  /** (5,1) */
  clms[9][9] = 2.4310490379465778321*ta - 2.8422412155745489079*ta3 + 0.068114108114108114108*ta*logx;
  clms[9][10] = -3.9266821912029810108*ta2 + 5.8346153846153846154*ta4 + 0.11764235764235764236*ta2*logx;
  clms[9][11] = 5.0329509017016084144*ta - 8.1023818776981170144*ta3 + 5.5166666666666666667*ta5 + 0.10948841471918394995*ta*logx;
  clms[9][12] = 0.101527299581346527388*ta - 8.5582844933153709418*ta2 + 13.014627403302616978*ta4 - 4.3125*ta6 + 0.26055494363186670879*ta2*logx;
    
  /** (5,2) */
  clms[10][7] = 0.86612970193593326616*ta - 0.18666666666666666667*ta3 + 0.089297369297369297369*ta*logx;
  clms[10][8] = -1.1544619678246051872*ta2 + 0.28*ta4;
  clms[10][9] = 0.7917293162706524283*ta - 3.4368503595170261837*ta3 - 0.022952032094889237746*ta*logx;
  clms[10][10] = 0.28053595936951014147*ta + 0.31120798949137087903*ta2 + 3.4084297924297924298*ta4 + 0.43448817848817848818*ta2*logx;
  clms[10][11] = 1.2370784757232070562*ta - 0.33664315124341216976*ta2 + 1.0618687980017224638*ta3 + 0.224*ta5 + 0.067212945017654593793*ta*logx + 0.47875324675324675325*ta3*logx;
  clms[10][12] = 0.088448031208912618076*ta + 1.0940182676326825551*ta2 - 0.096850502240490030478*ta4 - 0.224*ta6 + 0.99685163351097417031*ta2*logx;
    
  /** (5,3) */
  clms[11][9] = 4.2440104520863026147*ta - 2.2774548907882241216*ta3 + 0.1594005994005994006*ta*logx;
  clms[11][10] = 6.220542420977369858*ta2 + 0.019230769230769230769*ta4 + 0.70585414585414585415*ta2*logx;
  clms[11][11] = 7.90386884653014077*ta - 3.3995520908341421162*ta3 + 1.25*ta5 + 0.77096775751720806666*ta*logx;
  clms[11][12] = 1.827491392464237493*ta + 6.3772557882025475026*ta2 - 2.4145952259200122448*ta4 - 0.3125*ta6 + 0.51487384410461333538*ta2*logx;
    
  /** (5,4) */
  clms[12][7] = 2.5706537471616748009*ta - 0.18666666666666666667*ta3 + 0.35718947718947718948*ta*logx;
  clms[12][8] = 0.00022756242536462316682*ta2 + 0.28*ta4;
  clms[12][9] = 1.5873882331599828046*ta - 5.0400151946818613485*ta3 - 1.010173572459286745*ta*logx;
  clms[12][10] = 2.2442876749560811317*ta + 4.6400216555128308984*ta2 + 2.3209572649572649573*ta4 - 0.26665601065601065601*ta2*logx;
  clms[12][11] = 9.4305496551807111492*ta - 2.6931452099472973581*ta2 - 9.787214674612952795*ta3 + 0.224*ta5 - 0.75239512164096151537*ta*logx - 0.50006526806526806527*ta3*logx;
  clms[12][12] = -1.2093808162357422587*ta + 19.126696720272226*ta2 - 4.6959888564381848875*ta4 - 0.224*ta6 + 1.9263805641651795498*ta2*logx;
    
  /** (5,5) */
  clms[13][9] = 2.2417987079623586638*ta - 3.7078822412155745489*ta3 - 2.0773670773670773671*ta*logx;
  clms[13][10] = 6.5563587035650417216*ta2 - 0.12435897435897435897*ta4;
  clms[13][11] = 12.848844156281746158*ta - 6.9532102591931651761*ta3 + 1.25*ta5 + 0.0491165326330161495*ta*logx;
  clms[13][12] = 16.266814378442378433*ta2 - 5.1197728972515297302*ta4 - 0.3125*ta6 + 0.56911749219441527134*ta2*logx;
    
  /** (6,1) */
  clms[14][7] = 0.56892806352368092085*ta - 0.31481481481481481481*ta3 + 0.011732711732711732712*ta*logx;
  clms[14][8] = -1.6033912503693957942*ta2 + 0.29513888888888888889*ta4;
  clms[14][9] = 0.75121500050859794385*ta - 4.4250539269829393286*ta3 + 0.93333333333333333333*ta5 + 0.0056634269869563987211*ta*logx;
  clms[14][10] = 0.014217197523238524845*ta - 2.3705054864491033002*ta2 + 2.498488940329218107*ta4 + 0.073620515287181953849*ta2*logx;
    
  /** (6,2) */
    
  /** (6,3) */
  clms[16][7] = 1.2236359044183690609*ta - 0.31481481481481481481*ta3 + 0.10559440559440559441*ta*logx;
  clms[16][8] = -0.5567932520230314348*ta2 + 0.29513888888888888889*ta4;
  clms[16][9] = 1.3645438762626977901*ta - 3.6855156874446997904*ta3 - 0.10443257722669487375*ta*logx;
  clms[16][10] = 0.38386433312744017083*ta + 1.7466464752462693341*ta2 + 2.7049704218106995885*ta4 + 0.41609779109779109779*ta2*logx;
  clms[16][11] = 2.6988788534801879986*ta - 0.44784172198201353263*ta2 - 3.4373123153115704585*ta3 + 0.38040123456790123457*ta5 + 0.029486979827706358961*ta*logx + 0.18762626262626262626*ta3*logx;
  clms[16][12] = 0.082016610407290104641*ta + 4.1653548210535425076*ta2 - 1.8015926128654852153*ta4 - 0.2377507716049382716*ta6 + 0.9524074203920460605*ta2*logx;
    
  /** (6,4) */
  clms[17][9] = 4.2286874622205514798*ta - 2.545855379188712522*ta3 - 0.059132401989544846688*ta*logx;
  clms[17][10] = 5.5594354239492184991*ta2 - 0.011904761904761904762*ta4 + 0.53753653753653753654*ta2*logx;
  clms[17][11] = 8.539660694380415602*ta - 4.4828595860341892088*ta3 + 1.25*ta5 + 0.64650615046773510159*ta*logx;
  clms[17][12] = 1.6368582968553984161*ta + 7.1480353202353490562*ta2 - 2.7647974453530009086*ta4 - 0.3125*ta6 + 0.19564342701597603558*ta2*logx;
    
  /** (6,5) */
  clms[18][7] = 2.4537414136867447362*ta - 0.31481481481481481481*ta3 + 0.29331779331779331779*ta*logx;
  clms[18][8] = 0.28831684364261978641*ta2 + 0.29513888888888888889*ta4;
  clms[18][9] = 1.9767789053923672042*ta - 4.87310587503488738057*ta3 - 1.153442826236943884*ta*logx;
  clms[18][10] = 1.7771496904048156057*ta + 4.7621066674198963425*ta2 + 1.9846000514403292181*ta4 - 0.21354417187750521084*ta2*logx;
  clms[18][11] = 10.508208333282895882*ta - 2.0733413054722848733*ta2 - 9.2768927367069239217*ta3 + 0.38040123456790123457*ta5 - 0.65101418461487422227*ta*logx - 0.41553354053354053354*ta3*logx;
  clms[18][12] = -1.0515111247481700306*ta + 18.449083514387999662*ta2 - 4.9140202178393864855*ta4 - 0.2377507716049382716*ta6 + 1.7366316869157778249*ta2*logx;
    
  /** (6,6) */
  clms[19][9] = 2.1812186536338070473*ta - 3.7045855379188712522*ta3 - 2.1095475952618809762*ta*logx;
  clms[19][10] = 6.7199790687764824788*ta2 - 0.13095238095238095238*ta4;
    
  /** (7,1) */
  clms[20][9] = 2.2077220984432431257*ta - 2.6558527114549523513*ta3 + 0.022948037023667275768*ta*logx;
  clms[20][10] = -2.9774969573145188973*ta2 + 0.067226890756302521008*ta4 + 0.035527497712371661951*ta2*logx;
    
  /** (7,2) */
  clms[21][7] = 0.7995896565508844514*ta - 0.408163265306122449*ta3 + 0.027813083275268149218*ta*logx;
  clms[21][8] = -1.1592930085331132437*ta2 + 0.30612244897959183673*ta4;
  clms[21][9] = 1.0801827493164306069*ta - 3.999371896610792169*ta3 + 0.0019182597873516965612*ta*logx;
  clms[21][10] = 0.05533900612446692301*ta - 1.1831882849548764438*ta2 + 2.2425732197640961146*ta4 + 0.16408401522447140694*ta2*logx;
  clms[21][11] = 1.648649058114679313*ta - 0.06324457842796219772*ta2 - 6.197318869119520098*ta3 + 0.4956268221574344023*ta5 + 0.03350773508625026294*ta*logx + 0.05094360532535802644*ta3*logx;
  clms[21][12] = 0.031231792889999899199*ta - 1.6802001225328727291*ta2 - 2.5234866368730122104*ta4 - 0.24781341107871720117*ta6 + 0.33274715813475757533*ta2*logx;
  clms[21][13] = 2.338345281742978114*ta + 0.7602866485968883665*ta2 - 12.674159707537393383*ta3 + 11.066333125705446062*ta5 - 0.0023690507789108150684*ta*log2x + 0.11354373498668644652*ta*logx + 0.3128433282189350443*ta3*logx;
  clms[21][14] = 0.06797482667513325822*ta - 2.0465393072612297536*ta2 + 0.6860613932397195402*ta3 - 0.9519707914718048409*ta4 - 2.6978751591339483639*ta6 + 0.003094270405108003355*ta2*log2x + 0.5110544196066189615*ta2*logx + 0.21785220417816313317*ta4*logx;
  clms[21][15] = 3.0140459871385259285*ta + 1.5672672212875105207*ta2 - 25.931081853031804871*ta3 - 0.17686688254677189431*ta4 + 21.690156819063400939*ta5 - 0.5664306538942107455*ta7 - 0.0013389785622021850709*ta*log2x + 0.17727378890275170836*ta*logx + 0.8927146823431971037*ta3*logx + 0.419563601727896017*ta5*logx;
  clms[21][16] = 2.1890700049699915766e-40*(1.2916368170372250336e39*ta + 1.210921307977269819e39*ta2 + 1.0782437960169702961e40*ta3 - 1.8558447699769338357e40*ta4 + 1.0577498758382079791e41*ta6 + 9.7032755794961810547e38*ta8 - 1.6977875492072041457e38*ta2*log2x - 4.306524615965977082e37*ta*logx + 4.263846084512403994e39*ta2*logx + 4.444809433726965622e39*ta4*logx);
    
  /** (7,3) */
  clms[22][9] = 2.9714082846020944976*ta - 2.5771012108547122553*ta3 + 0.1065667525751559365*ta*logx;
  clms[22][10] = 0.45747586385332266568*ta2 + 0.033613445378151260504*ta4 + 0.26645623284278746464*ta2*logx;
  clms[22][11] = 5.2017019395102489322*ta - 5.3642030059854542589*ta3 + 1.25*ta5 + 0.27013779991077721315*ta*logx;
  clms[22][12] = 0.51490240660923170844*ta - 0.30506446371916066344*ta2 - 0.26027312231193293422*ta4 - 0.3125*ta6 + 0.41369243252519782946*ta2*logx;
  clms[22][13] = 9.3749738139594926204*ta + 2.4942904310002277676*ta2 - 8.5078328282082442166*ta3 + 4.6374899959983993597*ta5 + 0.22313250959270240249*ta*logx + 0.62237849685328676925*ta3*logx;
  clms[22][14] = 0.75553666098639847862*ta - 1.2019295900876114704*ta2 + 1.9530493665835735312*ta4 + 0.80812324929971988796*ta6 + 0.39707215581169703034*ta2*logx + 0.16933766513598446372*ta4*logx;
  clms[22][15] = 19.000719055695188081*ta + 2.4995219379755869047*ta2 - 15.349237125204803645*ta3 + 15.593673716342448653*ta5 - 1.4583333333333333333*ta7 - 0.034266627337817146525*ta*log2x + 0.3683151456627067275*ta*logx + 1.7864422209261106784*ta3*logx;
  clms[22][16] = 0.65388518968373147792*ta + 6.6248152375530218914*ta2 + 6.0914596749314341185*ta3 + 1.2839976586806555446*ta4 + 4.6794216956562767565*ta6 + 0.2734375*ta8 - 0.18005699453524984395*ta2*log2x + 1.4777258263268110824*ta2*logx + 0.49269059873866812122*ta4*logx;
    
  /** (7,4) */
  clms[23][7] = 1.4496675860582969489*ta - 0.40816326530612244898*ta3 + 0.11125233310107259687*ta*logx;
  clms[23][8] = -0.17610885417862078495*ta2 + 0.30612244897959183673*ta4;
  clms[23][9] = 1.7628403801577378217*ta - 3.7607791701429156411*ta3 - 0.207705700047528419744*ta*logx;
  clms[23][10] = 0.44271204899573538406*ta + 2.622186613722073201*ta2 + 2.2440214180910459422*ta4 + 0.37965587833535012407*ta2*logx;
  clms[23][11] = 3.9202302700514431948*ta - 0.50595662742369758179*ta2 - 4.8525125253466300192*ta3 + 0.49562682215743440233*ta5 - 0.009889324658014566267*ta*logx + 0.067759549152106174915*ta3*logx;
  clms[23][12] = 0.050748696060563808514*ta + 6.2949499758488010221*ta2 - 2.7279912119428394758*ta4 - 0.24781341107871720117*ta6 + 0.89361634401065439329*ta2*logx;
    
  /** (7,5) */
  clms[24][9] = 4.1411343707612865777*ta - 2.7243601144161368251*ta3 - 0.25934557972373098424*ta*logx;
  clms[24][10] = 5.4513763110008233384*ta2 - 0.033613445378151260504*ta4 + 0.44409372140464577439*ta2*logx;
  clms[24][11] = 9.0893186520691012541*ta - 5.0782146105532813157*ta3 + 1.25*ta5 + 0.54775152531476423262*ta*logx;
  clms[24][12] = 1.4302844628034214123*ta + 8.1811523906474301347*ta2 - 3.1180728855738810701*ta4 - 0.3125*ta6 + 0.057358137902275061783*ta2*logx;
    
  /** (7,6) */
  clms[25][7] = 2.3826585899891806148*ta - 0.40816326530612244898*ta3 + 0.25031774947741334296*ta*logx;
  clms[25][8] = 0.48630272791427036407*ta2 + 0.30612244897959183673*ta4;
  clms[25][9] = 2.2227922945116004366*ta - 4.7539031596654645874*ta3 - 1.2750080970954459669*ta*logx;
  clms[25][10] = 1.4941531653606069212*ta + 4.9580393284749195172*ta2 + 1.7264048476533470531*ta4 - 0.18332495835497035977*ta2*logx;
  clms[25][11] = 11.328548064369656554*ta - 1.7076036175549793385*ta2 - 8.9963953810820538925*ta3 + 0.49562682215743440233*ta5 - 0.58066111462207537202*ta*logx - 0.35759678496773334709*ta3*logx;
  clms[25][12] = -0.94869241550492118732*ta + 18.239931407916544305*ta2 - 5.0366249456958292461*ta4 - 0.24781341107871720117*ta6 + 1.6217540601204584763*ta2*logx;
    
  /** (7,7) */
  clms[26][9] = 2.1800107879019123863*ta - 3.7071532316630355846*ta3 - 2.14108848667672197084*ta*logx;
  clms[26][10] = 6.8277992433864045393*ta2 - 0.13445378151260504202*ta4;
  clms[26][11] = 12.893365173034811523*ta - 6.949257549546930005*ta3 + 1.25*ta5 + 0.033398701890324545579*ta*logx;
  clms[26][12] = 17.205804172171744272*ta2 - 5.1308309505440073354*ta4 - 0.3125*ta6 + 0.71355129707433477239*ta2*logx;
    
  /** (8,1) */
  clms[27][7] = 0.69899771528581122717*ta - 0.47916666666666666667*ta3 + 0.004469388970162963971*ta*logx;
  clms[27][8] = -1.383789167108242964*ta2 + 0.314453125*ta4;
  clms[27][9] = 1.0512485555498172281*ta - 3.958518359508628149*ta3 + 0.0027096906576395740482*ta*logx;
  clms[27][10] = 0.0037888411496807416939*ta - 2.1207386143913336474*ta2 + 1.9663921359388053467*ta4 + 0.030696401981231702594*ta2*logx;
  clms[27][11] = 1.6896813024244264682*ta - 0.0042624462933908344057*ta2 - 6.874142838316291662*ta3 + 0.583984375*ta5 + 0.0069181393624218524795*ta*logx + 0.006897689433244773802*ta3*logx;
  clms[27][12] = 0.0027317932958053683653*ta - 3.2553476998186743811*ta2 - 2.7117120206382908029*ta4 - 0.2554931640625*ta6 + 0.06170483253101548895*ta2*logx;
    
  /** (8,2) */
  clms[28][9] = 2.3106291009701235825*ta - 2.6029947554070361088*ta3 + 0.051117235018163810733*ta*logx;
  clms[28][10] = -1.8242168958625063865*ta2 + 0.050073099415204678363*ta4 + 0.091261679496973614621*ta2*logx;
  clms[28][11] = 4.2690352889287602691*ta - 6.2778898227616422959*ta3 + 1.25*ta5 + 0.080395141759776256234*ta*logx;
  clms[28][12] = 0.09865845189756911705*ta - 3.5286405321886623232*ta2 - 0.77876901052349267185*ta4 - 0.3125*ta6 + 0.19422599380233403372*ta2*logx;
    
  /** (8,3) */
  clms[29][7] = 1.0014844909878653913*ta - 0.47916666666666666667*ta3 + 0.040224500731466675739*ta*logx;
  clms[29][8] = -0.78844627961990137989*ta2 + 0.314453125*ta4;
  clms[29][9] = 1.3783122542484726548*ta - 3.7872575753366597665*ta3 - 0.021862061119994556527*ta*logx;
  clms[29][10] = 0.10229871104138002574*ta - 0.12744904307415542009*ta2 + 1.9877325964651211362*ta4 + 0.22523614156741091726*ta2*logx;
  clms[29][11] = 2.2512065045690092825*ta - 0.11508604992155252895*ta2 - 5.8121610732455340336*ta3 + 0.583984375*ta5 + 0.036718237405704139414*ta*logx + 0.042095550782421912453*ta3*logx;
  clms[29][12] = 0.050033504632424358351*ta + 0.29877711388960533985*ta2 - 2.4184188444303750323*ta4 - 0.2554931640625*ta6 + 0.44554470979302590744*ta2*logx;
    
  /** (8,4) */
  clms[30][9] = 3.1856593383140305417*ta - 2.6057620904112132182*ta3 + 0.055100042716141787349*ta*logx;
  clms[30][10] = 1.5470217536448190474*ta2 + 0.01279239766081871345*ta4 + 0.29203737439031556679*ta2*logx;
  clms[30][11] = 5.6478412888967329656*ta - 5.2045835137972462277*ta3 + 1.25*ta5 + 0.29387531006594144315*ta*logx;
  clms[30][12] = 0.63141409214444234912*ta + 1.8891664098961039713*ta2 - 1.0485049955295607804*ta4 - 0.3125*ta6 + 0.36205783847954289538*ta2*logx;
    
  /** (8,5) */
  clms[31][7] = 1.6001562643525075384*ta - 0.47916666666666666667*ta3 + 0.11173472425407409928*ta*logx;
  clms[31][8] = 0.091359664199942321232*ta2 + 0.314453125*ta4;
  clms[31][9] = 2.0547913312529041779*ta - 3.7907677530244690332*ta3 - 0.3176683755486000068*ta*logx;
  clms[31][10] = 0.47360514371009271174*ta + 3.2224895631384367854*ta2 + 1.9161278032320384294*ta4 + 0.34214774733350584744*ta2*logx;
  clms[31][11] = 4.9775933141905786522*ta - 0.53280578667385430071*ta2 - 5.5686406682653805182*ta3 + 0.583984375*ta5 - 0.04313414105981932841*ta*logx + 0.0059117848579439136715*ta3*logx;
  clms[31][12] = 0.011961462610107330221*ta + 7.8823385068976805763*ta2 - 3.2750758880263531356*ta4 - 0.2554931640625*ta6 + 0.84275839456194558729*ta2*logx;
    
  /** (8,6) */
  clms[32][9] = 4.0365010802691321658*ta - 2.8537605588044184535*ta3 - 0.43615826897560643691*ta*logx;
  clms[32][10] = 5.4921410744235233458*ta2 - 0.049342105263157894737*ta4 + 0.38329905388728918141*ta2*logx;
  clms[32][11] = 9.567125696371277875*ta - 5.4522876136593818958*ta3 + 1.25*ta5 + 0.46930982493146416888*ta*logx;
  clms[32][12] = 1.2430964939093708748*ta + 9.1332400782678827316*ta2 - 3.401695991206706042*ta4 - 0.3125*ta6 + 0.0033927074500644648926*ta2*logx;
    
  /** (8,7) */
  clms[33][7] = 2.3318576973637191899*ta - 0.47916666666666666667*ta3 + 0.21900005953798523458*ta*logx;
  clms[33][8] = 0.6338690020376092052*ta2 + 0.314453125*ta4;
  clms[33][9] = 2.3933726267215332998*ta - 4.6611123846355480127*ta3 - 1.378034874374851155*ta*logx;
  clms[33][10] = 1.299572514340494401*ta + 5.1485472634072771991*ta2 + 1.523006327668128655*ta4 - 0.16290452753301050515*ta2*logx;
  clms[33][11] = 11.98738156679547853*ta - 1.4620190786330562011*ta2 - 8.8064030699260740948*ta3 + 0.583984375*ta5 - 0.52696407147237839226*ta*logx - 0.31481258558585377471*ta3*logx;
  clms[33][12] = -0.87136059365697232105*ta + 18.174389380077077065*ta2 - 5.1125308457054834997*ta4 - 0.2554931640625*ta6 + 1.5450254569475920835*ta2*logx;
    
  /** (8,8) */
  clms[34][9] = 2.2058573983497587983*ta - 3.71206952566601689409*ta3 - 2.16950218683964813996*ta*logx;
  clms[34][10] = 6.9028605005619598973*ta2 - 0.1363304093567251462*ta4;
  clms[34][11] = 13.097578254912827461*ta - 6.9647129334411318612*ta3 + 1.25*ta5 + 0.028097856152823292034*ta*logx;
  clms[34][12] = 17.456381660169843359*ta2 - 5.1382329729999205222*ta4 - 0.3125*ta6 + 0.76617783510402053664*ta2*logx;

  /** Amplitudes (correct with spin terms) */
  // Taylor - Used at 3PN for lm = 62
  //          Used at 5PN for lm = 61, 66, 71
  //          Used at 8PN for lm = 72, 73
  //          Used at 6PN for all the other ones
  double rholms[KMAX];
  int nTaylor[KMAX];

  for (int k=0; k<KMAX; k++) nTaylor[k] = 6;
  nTaylor[15] = 3;
  nTaylor[14] = 5;
  nTaylor[19] = 5;
  nTaylor[20] = 5;
  nTaylor[21] = 8;
  nTaylor[22] = 8;

  double v = sqrt(x);
  for (int k=0; k<KMAX; k++) {
    rholms[k] = Taylorseries(v,clms[k],nTaylor[k]*2);
    flm[k]    = gsl_pow_int(rholm[k]/rholms[k],LINDEX[k]);
  }
  
}

/** Calculate tidal correction to multipolar waveform amplitude
    Ref. Damour, Nagar & Villain, Phys.Rev. D85 (2012) 123007 */
#define use_fmode_22amplitude_correction (1)
/**
 * Function: eob_wav_hlmTidal
 * --------------------------
 *   Calculate tidal correction to multipolar waveform amplitude
 *   Ref. Damour, Nagar & Villain, Phys.Rev. D85 (2012) 123007
 *   @param[in] x: x = (M omega)^{2/3}
 *   @param[in] dyn: dynamics structure
 *   @param[out] hTidallm: tidal correction to multipolar waveform amplitude
*/
void eob_wav_hlmTidal(double x, Dynamics *dyn, double *hTidallm)
{
  const double nu       = EOBPars->nu;  
  const double XA       = EOBPars->X1;
  const double XB       = EOBPars->X2;

  double khatA_2 = 0.5*EOBPars->kapA2; 
  double khatB_2 = 0.5*EOBPars->kapB2; 
  const double kapA2j   = EOBPars->japA2;
  const double kapB2j   = EOBPars->japB2;
  const double kapT2j   = EOBPars->japT2;
  
#if (use_fmode_22amplitude_correction)
  if (EOBPars->use_tidal_fmode_model) {
    const double fact22A  = eob_wav_hlmTidal_fmode_fact22A(x, dyn->dress_tides_fmode_A[2], EOBPars->bomgfA[2], XB); 
    const double fact22B  = eob_wav_hlmTidal_fmode_fact22A(x, dyn->dress_tides_fmode_B[2], EOBPars->bomgfB[2], XA);         
    khatA_2 *= fact22A;
    khatB_2 *= fact22B;
  }
#endif
  
  const double x5 = gsl_pow_int(x,5);
  const double x6 = gsl_pow_int(x,6);
  
  double hA[KMAX], hB[KMAX], betaA1[KMAX],betaB1[KMAX];

  memset(hTidallm, 0., KMAX*sizeof(double));
  memset(hA, 0., KMAX*sizeof(double));
  memset(hB, 0., KMAX*sizeof(double));
  memset(betaA1, 0., KMAX*sizeof(double));
  memset(betaB1, 0., KMAX*sizeof(double));

  /** l=2 */
  hA[1]     = 2 * khatA_2 *(XA/XB+3);
  hB[1]     = 2 * khatB_2 *(XB/XA+3);

  betaA1[1] = (-202. + 560*XA - 340*XA*XA + 45*XA*XA*XA)/(42*(3-2*XA));
  betaB1[1] = (-202. + 560*XB - 340*XB*XB + 45*XB*XB*XB)/(42*(3-2*XB));
    
  hA[0]     = 3 * khatA_2 * (3-4*XA);
  hB[0]     = 3 * khatB_2 * (3-4*XB);
    
  /** l=3 */  
  hA[2] = 12 * khatA_2 * XB;
  hB[2] = 12 * khatB_2 * XA;

  betaA1[2] = (-6. -5.*XA +131.*XA*XA -130.*XA*XA*XA)/(36.*(1.-XA));
  betaB1[2] = (-6. -5.*XB +131.*XB*XB -130.*XB*XB*XB)/(36.*(1.-XB));

  hA[4] = hA[2];
  hB[4] = hB[2];

  betaA1[4] = ( (XA-3.)*(10.*XA*XA - 25.*XA+ 14.) )/(12.*(1.-XA));
  betaB1[4] = ( (XB-3.)*(10.*XB*XB - 25.*XB+ 14.) )/(12.*(1.-XB));
  
  /** l=2 */
  /* (2,1) */
  hTidallm[0] = ( -hA[0] + hB[0] )*x5;
  /* (2,2) */
  hTidallm[1] = ( hA[1]*(1. + betaA1[1]*x) + hB[1]*(1. + betaB1[1]*x) )*x5;
  
  /** l=3 */
  /* (3,1) */
  hTidallm[2] = ( -hA[2]*(1. + betaA1[2]*x) + hB[2]*(1. + betaB1[2]*x) )*x5;
  /* (3,2) */
  hTidallm[3] = 8.*( khatA_2*(1. -2.*XB + 3.*XB*XB) +khatB_2*(1. -2.*XA + 3.*XA*XA) )*x5/(1.-3.*nu);
  /* (3,3) */
  hTidallm[4] = ( -hA[4]*(1. + betaA1[4]*x) + hB[4]*(1. + betaB1[4]*x) )*x5;
  
  if ( (EOBPars->use_tidal_gravitomagnetic==TIDES_GM_GSF) || (EOBPars->use_tidal_gravitomagnetic==TIDES_GM_PN) ) {
    const double fourtnine= 1.5555555555555555556;  // 14/9 = 112/(3*24)
    const double fourthird = 1.3333333333333333333; // 32/24 = 4/3
    hTidallm[0] += 0.5*( -1.*kapA2j/XB + kapB2j/XA )*x5;
    hTidallm[1] += fourtnine*kapT2j*x6;
    hTidallm[2] += 0.5*( kapA2j*(4. - 17.*XB) - kapB2j*(4. - 17.*XA) )*x6;
    hTidallm[3] += fourthird*kapT2j*x5/(1.-3.*nu);
    hTidallm[4] += 0.5*( kapA2j*(4. - 9.*XB) - kapB2j*(4. - 9.*XA) )*x6;
  }


  /* OLD STUFF 
   // l=2 
  hA[1]     = 2 * khatA_2 *(XA/XB+3);
  hB[1]     = 2 * khatB_2 *(XB/XA+3);

  betaA1[1] = (-202. + 560*XA - 340*XA*XA + 45*XA*XA*XA)/(42*(3-2*XA));
  betaB1[1] = (-202. + 560*XB - 340*XB*XB + 45*XB*XB*XB)/(42*(3-2*XB));
    
  hA[0]     = 3 * khatA_2 * XB * (3-4*XA)/XA;
  hB[0]     = 3 * khatB_2 * XA * (3-4*XB)/XB;
    
  // l=3 
  hA[2] = hA[4];
  hB[2] = hB[4];
  
  hA[4] = 12 * khatA_2 * XB*XB/XA;
  hB[4] = 12 * khatB_2 * XA*XA/XB;
  
  // l=2 
  // (2,1) 
  hTidallm[0] = ( -hA[0] + hB[0] )*x5;
  // (2,2) 
  hTidallm[1] = ( hA[1]*(1. + betaA1[1]*x) + hB[1]*(1. + betaB1[1]*x) )*x5;
  
  // l=3 
  // (3,1) 
  hTidallm[2] = ( -hA[2] + hB[2] )*x5;
  // (3,3) 
  hTidallm[4] = ( -hA[4] + hB[4] )*x5; */
  
}

/**
 * Function: eob_wav_hlmNQC_find_a1a2a3
 * ------------------------------------
 *   Computes the factors and the coefficients that build the
 *   NQC corrections to the waveform in the spinning case
 * 
 *   @param[in]  dyn : dynamics structure
 *   @param[in]  h   : multipolar waveform (modified by this function)
 *   @param[out] hnqc: NQC part of the waveform
*/
void eob_wav_hlmNQC_find_a1a2a3_ecc(Dynamics *dyn, Waveform_lm *h, Waveform_lm *hnqc)
{
  double A_tmp, dA_tmp, omg_tmp, domg_tmp;
  double alpha1[KMAX], omega1[KMAX];
  double c1A[KMAX], c2A[KMAX], c3A[KMAX], c4A[KMAX];
  double c1phi[KMAX], c2phi[KMAX], c3phi[KMAX], c4phi[KMAX];
	    
  const double nu   = EOBPars->nu;
  const double chi1 = EOBPars->chi1;
  const double chi2 = EOBPars->chi2;
  const double X1   = EOBPars->X1;
  const double X2   = EOBPars->X2;
  const double aK   = EOBPars->a1+EOBPars->a2;
  const double Mbh  = EOBPars->Mbhf;
  const double abh  = EOBPars->abhf;
    
  double *t       = h->time;
  double *r       = dyn->data[EOB_RAD];
  double *w       = dyn->data[EOB_MOMG]; /* Omega */
  double *pph     = dyn->data[EOB_PPHI];
  double *pr_star = dyn->data[EOB_PRSTAR];
  double *Omg_orb = dyn->data[EOB_OMGORB]; /* Omega orbital */
  double *ddotr   = dyn->data[EOB_DDOTR];
  
  double P[2], M[4];
  double max_A[KMAX],max_dA[KMAX],d2max[KMAX],d3max[KMAX],max_omg[KMAX],max_domg[KMAX],maxd2omg[KMAX], DeltaT[KMAX];
  double ai[KMAX][2] = {0.};
  double bi[KMAX][2] = {0.};
  
  const int size = h->size;
  for (int i = 0; i < size; i++) {
    hnqc->time[i] = t[i];
  }
  
  double *omg[KMAX], *domg[KMAX];
  double *n1[KMAX],*n2[KMAX],*n3[KMAX],*n4[KMAX],*n5[KMAX],*n6[KMAX],
    *d_n4[KMAX],*d_n5[KMAX],*d_n6[KMAX],*d2_n4[KMAX],*d2_n5[KMAX],*d2_n6[KMAX];  double *m11[KMAX], *m12[KMAX], *m13[KMAX], *m21[KMAX], *m22[KMAX];
  double *p1tmp[KMAX], *p2tmp[KMAX]; /* RWZ amplitude and derivative */

  for (int k=0; k<KMAX; k++) {
    omg[k]  = (double*) calloc (size,sizeof(double));
    domg[k] = (double*) calloc (size,sizeof(double));
    m11[k] = (double*) calloc (size,sizeof(double));
    m12[k] = (double*) calloc (size,sizeof(double));
    m13[k] = (double*) calloc (size,sizeof(double));
    m21[k] = (double*) calloc (size,sizeof(double));
    m22[k] = (double*) calloc (size,sizeof(double));
    p1tmp[k] = (double*) calloc (size,sizeof(double));
    p2tmp[k] = (double*) calloc (size,sizeof(double));
    n1[k] = (double*) calloc (size,sizeof(double));
    n2[k] = (double*) calloc (size,sizeof(double));
    n4[k] = (double*) calloc (size,sizeof(double));
    n5[k] = (double*) calloc (size,sizeof(double));
    d_n4[k] = (double*) calloc (size,sizeof(double));
    d_n5[k] = (double*) calloc (size,sizeof(double));
    d2_n4[k] = (double*) calloc (size,sizeof(double));
    d2_n5[k] = (double*) calloc (size,sizeof(double));
  }

  /** omega derivatives */
  const double dt = t[1]-t[0];
  for (int k=0; k<KMAX; k++) {
    if(h->kmask[k]){
      D0(h->phase[k], dt, size, omg[k]);
      D0(omg[k], dt, size, domg[k]);
    }
  }
  /** NR fits */
  for (int k=0; k<KMAX; k++) {   
    max_A[k]    = 0.;
    max_dA[k]   = 0.;
    max_omg[k]  = 0.;
    max_domg[k] = 0.;
  }
  
  if (EOBPars->use_flm == USEFLM_HM || EOBPars->use_flm == USEFLM_HM_4PN22) {
    
    /* Higher modes */
    /* Choosing modes using kpostpeak array */
    int kpostpeak_size = EOBPars->kpostpeak_size;  
    int *kpostpeak     = EOBPars->kpostpeak;

    QNMHybridFitCab_HM(nu, X1, X2, chi1, chi2, aK,  Mbh, abh,  
            c1A, c2A, c3A, c4A, c1phi, c2phi, c3phi, c4phi,
            alpha1, omega1);
    
    /* 22, 31, 33, 41 and 55 fitted directly + 44 dA */
    eob_nqc_point_HM(dyn, max_A, max_dA, max_omg, max_domg);
    /* 21 fitted directly at tpeak_22*/
    eob_nqc_point_HM_peak22(dyn, max_A, max_dA, max_omg, max_domg);
    
    /* Over-writing fits using postpeak quantities for modes in kpostpeak */
    for (int j=0; j<kpostpeak_size; j++) {
      int k = kpostpeak[j];
      
      /* Normalizing c1A and c4A */
      int l = LINDEX[k];
      c1A[k] /= sqrt((l+2)*(l+1)*l*(l-1));
      c4A[k] /= sqrt((l+2)*(l+1)*l*(l-1));
      
      eob_nqc_point_postpeak(Mbh,c1A[k],c2A[k],c3A[k],c4A[k],
            c1phi[k],c2phi[k],c3phi[k],c4phi[k],alpha1[k],omega1[k],
            &A_tmp,&dA_tmp,&omg_tmp,&domg_tmp);
      
      max_A[k]    = A_tmp;
      max_dA[k]   = dA_tmp;
      max_omg[k]  = omg_tmp;
      max_domg[k] = domg_tmp;
    }
  } else {
    
    eob_nqc_point(dyn, &A_tmp, &dA_tmp, &omg_tmp, &domg_tmp);
	
    /** Switch on the 22 values (only) */
    max_A[1]    = A_tmp;
    max_dA[1]   = dA_tmp;
    max_omg[1]  = omg_tmp;
    max_domg[1] = domg_tmp;
  
  }

  if (VERBOSE) {
    printf("NR values for NQC determination:\n");
    PRFORMd("A22_mrg",max_A[1]);
    PRFORMd("dA22_mrg",max_dA[1]);
    PRFORMd("omg22_mrg",max_omg[1]);
    PRFORMd("domg22_mrg",max_domg[1]);
  }
  
  /** NQC corrections to AMPLITUDE (n1,n2,n3) and PHASE (n4,n5,n6)
   * NQC basis for (2,2) waveform : AMPLITUDE
   * note: n3 and n6 are not used
   */
  double pr_star2, r2, w2;
  for (int k=0; k<KMAX; k++) {   
    if(h->kmask[k]) {
      for (int j=0; j<size; j++) {
        pr_star2 = SQ(pr_star[j]);
        r2       = SQ(r[j]);
        w2       = SQ(w[j]);                  //CHECKME: Omg or Omg_orbital ?
        n1[k][j]  = pr_star2/(r2*w2);         /* [pr*\/(r Omg)]^2 */
        n2[k][j]  = ddotr[j]/(r[j]*w2);       /* [ddot{r}/(r Omg^2)] */
        n4[k][j]  = pr_star[j]/(r[j]*w[j]);   /* pr*\/(r Omg) */
        n5[k][j]  = n4[k][j]*r2*w2;           /* (pr*)*(r Omg) */
      }
    }
  }

  if (EOBPars->use_flm == USEFLM_HM || EOBPars->use_flm == USEFLM_HM_4PN22) {
    for (int j=0; j<size; j++) {
      /* l=2,m=1 */
      if(h->kmask[0]) {
        n2[0][j] = cbrt(SQ(w[j]))*n1[0][j];
        n5[0][j] = cbrt(SQ(w[j]))*n4[0][j];
      }
      /* l=2,m=2 */
      if(h->kmask[1] && (EOBPars->use_flm == USEFLM_HM_4PN22)) {
        n2[1][j] = SQ(pr_star[j])*n1[1][j];
      }
      /* l=3, l=4 & l=5 */
      // FIXME: condition on spins bad for non-spinning limit
      if ((chi1 <= 0.) || (chi2 <= 0.)) {
	      for (int k=2; k<14; k++) {
	        if(h->kmask[k]){
	        n5[k][j]  = cbrt(SQ(w[j]))*n4[k][j];
	      }
	    }
      }
      /* l=3, m=2 */
      if(h->kmask[3]){
        n2[3][j] = cbrt(SQ(w[j]))*n1[3][j];
        n5[3][j] = cbrt(SQ(w[j]))*n4[3][j];
      }
      /* l=4, m=2 */
      if(h->kmask[6]){
        n2[6][j] = cbrt(SQ(w[j]))*n1[6][j];
        n5[6][j] = cbrt(SQ(w[j]))*n4[6][j];
      }
      /* l=4, m=3 */
      if(h->kmask[7]){
        n2[7][j] = cbrt(SQ(w[j]))*n1[7][j];
        n5[7][j] = cbrt(SQ(w[j]))*n4[7][j];
      }
      /* l=5, m=5 */
      if(h->kmask[13]){
        n5[13][j] = cbrt(SQ(w[j]))*n4[13][j];
      }
    }
  }  

#if (DEBUG)
  FILE* fp_dbg;
  fp_dbg = fopen("nqc_nfunc.txt", "w");
  for (int j=0; j<size; j++) {
    fprintf(fp_dbg, "%20.12f\t%.16e\t%.16e\t%.16e\t%.16e\n", t[j], n1[1][j], n2[1][j], n4[1][j], n5[1][j]);
  }
  fclose(fp_dbg);
#endif    
  
  /** Derivatives for the phase */
  for (int k=0; k<KMAX; k++) {
    if(h->kmask[k]) {
      D0(n4[k],dt,size,d_n4[k]);
      D0(n5[k],dt,size,d_n5[k]);
      D0(d_n4[k],dt,size,d2_n4[k]);
      D0(d_n5[k],dt,size,d2_n5[k]);
    }
  }

#if (DEBUG)
  fp_dbg = fopen("nqc_dfunc.txt", "w");
  for (int j=0; j<size; j++) {
    fprintf(fp_dbg, "%f\t%.16e\t%.16e\t%.16e\t%.16e\n", t[j], d_n4[1][j], d_n5[1][j], d2_n4[1][j], d2_n5[1][j]);  
  }
  fclose(fp_dbg);
#endif    

  /** Find max Omg */
  int Omgmax_index = 0;
  double Omg_max   = Omg_orb[0];
  for (int j=0; j<size; j++) {
    if (Omg_orb[j] > Omg_max) {
      Omg_max = Omg_orb[j];
      Omgmax_index = j;
    }
  }
  //TODO:
  //Test the search backwards
  /*
    int Omgmax_index = size-1;
    double Omg_max   = Omg_orb[Omgmax_index];
    for (int j=(size-2); j--; ) {
    if (Omg_orb[j] < Omg_max) 
    break;
    Omg_max = Omg_orb[j];
    Omgmax_index = j;
    }
  */

  /** Time */
  double tOmgOrb_pk = t[Omgmax_index];
  double DeltaT_nqc = eob_nqc_timeshift(nu, chi1);
  double tNQC = tOmgOrb_pk - DeltaT_nqc;

  if (VERBOSE) {
    printf("NQC info:\n");
    PRFORMd("DeltaT_tNQC",DeltaT_nqc);
    PRFORMd("tNQC[bare]",tNQC);
  }

  /** Find jmax: t[jmax] <= tNQC */
  double tmrg[KMAX];
  tmrg[1] = tNQC - 2.;
	
  int jmax = 0;
  for (int j=0; j<size; j++) {
    if(t[j] > tNQC) {
      jmax = j-2;
      break;
    }
  }
	     
  double dtmrg[KMAX];
  double t_NQC[KMAX];
  int    j_NQC[KMAX];

  if (EOBPars->use_flm == USEFLM_HM || EOBPars->use_flm == USEFLM_HM_4PN22) {
    int modesatpeak22[KMAX]; 
    set_multipolar_idx_mask (modesatpeak22, KMAX, EOBPars->knqcpeak22, EOBPars->knqcpeak22_size, 0);      
    eob_nqc_deltat_lm(dyn, dtmrg);

    for (int k=0; k<KMAX; k++) {   

      if(h->kmask_nqc[k]){
        tmrg[k]  = tmrg[1] + dtmrg[k];
        t_NQC[k] = tmrg[k] + 2.;
        if(modesatpeak22[k]) t_NQC[k] = tmrg[1];

        j_NQC[k] = size-1;
        for (int j=size-2; j>=0; j--) {
          if(t[j] < t_NQC[k]) {
            break;
          }
        j_NQC[k] = j;
        } 
      }
    }
  }
  
  /** Solve the linear systems */
  
  /* Regge-Wheeler-Zerilli normalized amplitude. 
     The ringdown coefficient refer to this normalization.
     Nagar & Rezzolla, CQG 22 (2005) R167 */      
  for (int k=0; k<KMAX; k++) {
    if(h->kmask[k]){  
      double nlm = 1./(sqrt( (LINDEX[k]+2)*(LINDEX[k]+1)*LINDEX[k]*(LINDEX[k]-1) ) );
      if (h->ampli[k][0] > 0.) {	
	nNegAmp[k] = 0;	
      } else {		
	nNegAmp[k] = 1;	
      }
      for (int j=0; j<size; j++) {
	p1tmp[k][j] = fabs(h->ampli[k][j] * nlm);
      }
    }
  }
  
  /* Matrix elements: waveform amplitude at all points */
  for (int k=0; k<KMAX; k++) {
    if(h->kmask[k]) {
      for (int j=0; j<size; j++) {
	m11[k][j] = n1[k][j] * p1tmp[k][j];
	m12[k][j] = n2[k][j] * p1tmp[k][j];
      }
    }
  }
  
  /* Take FD derivatives */
  for (int k=0; k<KMAX; k++) {
    if(h->kmask[k]) {
      D0(m11[k],dt,size, m21[k]);
      D0(m12[k],dt,size, m22[k]);
      D0(p1tmp[k],dt,size, p2tmp[k]);
    }
  }

#if (DEBUG)
  fp_dbg = fopen("nqc_amp_func.txt", "w");
  for (int j=0; j<size; j++) {
    fprintf(fp_dbg, "%e\t%e\t%e\n", t[j], p1tmp[1][j], p2tmp[1][j]);
  }
  fclose(fp_dbg);  
#endif
  
  double detM = 1.;
  double oodetM = 1.;
  for (int k=0; k<KMAX; k++) {
    if(h->kmask_nqc[k]) {  
      
      if (EOBPars->use_flm == USEFLM_HM || EOBPars->use_flm == USEFLM_HM_4PN22) {
	      jmax = j_NQC[k];
      }
      
      
      /* Computation of ai coefficients at Omega peak */
      P[0]     = max_A[k]  - p1tmp[k][jmax];
      P[1]     = max_dA[k] - p2tmp[k][jmax];
      
      M[0]     = m11[k][jmax];
      M[1]     = m12[k][jmax];
      M[2]     = m21[k][jmax];
      M[3]     = m22[k][jmax];
      
      /* detM     = M[0]*M[3]-M[1]*M[2]; 
	 ai[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
	 ai[k][1] = (M[0]*P[1] - M[2]*P[0])/detM; */
      /* safe version (amplitude can be zero) */
      oodetM   = 1.0/(M[0]*M[3]-M[1]*M[2]);
      if (isfinite(oodetM)) {
        ai[k][0] = (M[3]*P[0] - M[1]*P[1])*oodetM;
        ai[k][1] = (M[0]*P[1] - M[2]*P[0])*oodetM;
      }
      
      /* Computation of bi coefficients at Omega peak */
      P[0]     = omg[k][jmax]   - max_omg[k];
      P[1]     = domg[k][jmax]  - max_domg[k];
      
      M[0]     = d_n4[k][jmax];
      M[1]     = d_n5[k][jmax];
      M[2]     = d2_n4[k][jmax];
      M[3]     = d2_n5[k][jmax];
      
      /* detM     =  M[0]*M[3] - M[1]*M[2];
	 bi[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
	 bi[k][1] = (M[0]*P[1] - M[2]*P[0])/detM; */
      /* safe version (phase can be zero) */
      oodetM   = 1.0/(M[0]*M[3]-M[1]*M[2]);
      if (isfinite(oodetM)) {
        bi[k][0] = (M[3]*P[0] - M[1]*P[1])*oodetM;
        bi[k][1] = (M[0]*P[1] - M[2]*P[0])*oodetM;
      }
    }
  }
  
  if (VERBOSE){
    printf("NQC coefficients for 22 mode:\n");
    PRFORMd("a1",ai[1][0]);
    PRFORMd("a2",ai[1][1]);
    PRFORMd("b1",bi[1][0]);
    PRFORMd("b2",bi[1][1]);
  }

  /* Defining sigmoid function to switch on NQCs near the end of the evolution */
  /*
    Old configuration used for arXiv:2001.11736
    double t0 = tNQC - 30.;
    double alpha = 0.09;
  */
  double t0 = tNQC - EOBPars->delta_t0_sigmoid_NQC;
  double alpha = EOBPars->alpha_sigmoid_NQC;
  double *sigmoid;
  sigmoid = (double*) calloc (size, sizeof(double));
  double fact[KMAX];
  for (int k=0; k<KMAX; k++) {
    if(h->kmask[k]){
      // Correction needed to have sigmoid = 1 when attaching ringdown
      fact[k] = 1. + exp(-alpha*(t[j_NQC[k]] - t0));
    }
  }   

  for (int j=0; j<size; j++) {
    sigmoid[j] = 1./(1. + exp(-alpha*(t[j] - t0)));
    for (int k=0; k<KMAX; k++) {
      if(h->kmask[k]){
	n1[k][j] = n1[k][j]*sigmoid[j]*fact[k];
	n2[k][j] = n2[k][j]*sigmoid[j]*fact[k];
	n4[k][j] = n4[k][j]*sigmoid[j]*fact[k];
	n5[k][j] = n5[k][j]*sigmoid[j]*fact[k];
      }
    }
  }
    
  free(sigmoid);
  
  /** Set amplitude and phase */
  for (int k=0; k<KMAX; k++) {
    if(h->kmask_nqc[k]) {  
      for (int j=0; j<size; j++) {
	hnqc->ampli[k][j] = 1. + ai[k][0]*n1[k][j] + ai[k][1]*n2[k][j];
	hnqc->phase[k][j] =      bi[k][0]*n4[k][j] + bi[k][1]*n5[k][j];
      }
    }
  }
  
  /** Multiply waveform to NQC */
  for (int k=0; k<KMAX; k++) {
    if(h->kmask_nqc[k]){    
      for (int j=0; j<size; j++) {
	h->ampli[k][j] *= hnqc->ampli[k][j];
	h->phase[k][j] -= hnqc->phase[k][j];
      }
    }
  }
  
  if (EOBPars->output_nqc_coefs) {
    /** Output the NQC coefficients */
    FILE* fp;
    char fname[STRLEN];
    strcpy(fname, EOBPars->output_dir);
    strcat(fname, "/nqc_coefs.txt");
    fp = fopen(fname, "w");
    fprintf(fp, "# q=%e chizA=%e chizB=%e f0=%e\n",EOBPars->q,EOBPars->chi1,EOBPars->chi2,EOBPars->initial_frequency);
    fprintf(fp, "# M=%e LambdaA=[%e,%e,%e] LambdaBl2=[%e,%e,%e]\n",EOBPars->M,
	    EOBPars->LambdaAl2, EOBPars->LambdaAl3,EOBPars->LambdaAl4,
	    EOBPars->LambdaBl2,EOBPars->LambdaBl3,EOBPars->LambdaBl4);
    for (int k=0; k<KMAX; k++) {
      if(h->kmask_nqc[k]){ 
        fprintf(fp, "%d %d %d %e %e %e %e\n", k, LINDEX[k], MINDEX[k], 
	      ai[k][0], ai[k][1], 
	      bi[k][0], bi[k][1]);
      }
    }  
    fclose(fp);  
  }
  
  /** Free mem */
  for (int k=0; k<KMAX; k++) {
    free(omg[k]);
    free(domg[k]);
    free(m11[k]);
    free(m12[k]);
    free(m13[k]);
    free(m21[k]);
    free(m22[k]);
    free(p1tmp[k]);
    free(p2tmp[k]);
    free(n1[k]);
    free(n2[k]);
    free(n4[k]);
    free(n5[k]);
    free(d_n4[k]);
    free(d_n5[k]);
    free(d2_n4[k]);
    free(d2_n5[k]);
  }
  
}

/**
 * Function: eob_wav_hlmNQC_find_a1a2a3
 * ------------------------------------
 *   Computes the factors and the coefficients that build the
 *   NQC corrections to the waveform in the spinning quasi-circ case
 * 
 *   @param[in]  dyn : dynamics structure
 *   @param[in]  h   : multipolar waveform (modified by this function)
 *   @param[out] hnqc: NQC part of the waveform
*/
void eob_wav_hlmNQC_find_a1a2a3_circ(Dynamics *dyn, Waveform_lm *h, Waveform_lm *hnqc)
{  
  double A_tmp, dA_tmp, omg_tmp, domg_tmp;
  double alpha1[KMAX], omega1[KMAX];
  double c1A[KMAX], c2A[KMAX], c3A[KMAX], c4A[KMAX];
  double c1phi[KMAX], c2phi[KMAX], c3phi[KMAX], c4phi[KMAX];
	    
  const double nu   = EOBPars->nu;
  const double chi1 = EOBPars->chi1;
  const double chi2 = EOBPars->chi2;
  const double X1   = EOBPars->X1;
  const double X2   = EOBPars->X2;
  const double aK   = EOBPars->a1+EOBPars->a2;
  const double Mbh  = EOBPars->Mbhf;
  const double abh  = EOBPars->abhf;
    
  double *t       = h->time;
  double *r       = dyn->data[EOB_RAD];
  double *w       = dyn->data[EOB_MOMG]; /* Omega */
  double *pph     = dyn->data[EOB_PPHI];
  double *pr_star = dyn->data[EOB_PRSTAR];
  double *Omg_orb = dyn->data[EOB_OMGORB]; /* Omega orbital */
  double *ddotr   = dyn->data[EOB_DDOTR];
  
  double P[2], M[4];
  double max_A[KMAX],max_dA[KMAX],d2max[KMAX],d3max[KMAX],max_omg[KMAX],max_domg[KMAX],maxd2omg[KMAX], DeltaT[KMAX];
  double ai[KMAX][2] = {0.};
  double bi[KMAX][2] = {0.};
  
  const int size = h->size;
  for (int i = 0; i < size; i++) {
    hnqc->time[i] = t[i];
  }
  
  double *omg[KMAX], *domg[KMAX];
  double *n1[KMAX],*n2[KMAX],*n3[KMAX],*n4[KMAX],*n5[KMAX],*n6[KMAX],
    *d_n4[KMAX],*d_n5[KMAX],*d_n6[KMAX],*d2_n4[KMAX],*d2_n5[KMAX],*d2_n6[KMAX];  double *m11[KMAX], *m12[KMAX], *m13[KMAX], *m21[KMAX], *m22[KMAX];
  double *p1tmp[KMAX], *p2tmp[KMAX]; /* RWZ amplitude and derivative */

  for (int k=0; k<KMAX; k++) {
    omg[k]  = (double*) calloc (size,sizeof(double));
    domg[k] = (double*) calloc (size,sizeof(double));
    m11[k] = (double*) calloc (size,sizeof(double));
    m12[k] = (double*) calloc (size,sizeof(double));
    m13[k] = (double*) calloc (size,sizeof(double));
    m21[k] = (double*) calloc (size,sizeof(double));
    m22[k] = (double*) calloc (size,sizeof(double));
    p1tmp[k] = (double*) calloc (size,sizeof(double));
    p2tmp[k] = (double*) calloc (size,sizeof(double));
    n1[k] = (double*) calloc (size,sizeof(double));
    n2[k] = (double*) calloc (size,sizeof(double));
    n4[k] = (double*) calloc (size,sizeof(double));
    n5[k] = (double*) calloc (size,sizeof(double));
    d_n4[k] = (double*) calloc (size,sizeof(double));
    d_n5[k] = (double*) calloc (size,sizeof(double));
    d2_n4[k] = (double*) calloc (size,sizeof(double));
    d2_n5[k] = (double*) calloc (size,sizeof(double));
  }

  /** omega derivatives */
  const double dt = t[1]-t[0];
  for (int k=0; k<KMAX; k++) {
    if(h->kmask_nqc[k]){
      D0(h->phase[k], dt, size, omg[k]);
      D0(omg[k], dt, size, domg[k]);
    }
  }
  /** NR fits */
  for (int k=0; k<KMAX; k++) {   
    max_A[k]    = 0.;
    max_dA[k]   = 0.;
    max_omg[k]  = 0.;
    max_domg[k] = 0.;
  }
  
  if (EOBPars->use_flm == USEFLM_HM || EOBPars->use_flm == USEFLM_HM_4PN22) {
    /* Higher modes */
    /* Choosing modes using kpostpeak array */
    int kpostpeak_size = EOBPars->kpostpeak_size;  
    int *kpostpeak     = EOBPars->kpostpeak;

    QNMHybridFitCab_HM(nu, X1, X2, chi1, chi2, aK,  Mbh, abh,  
            c1A, c2A, c3A, c4A, c1phi, c2phi, c3phi, c4phi,
            alpha1, omega1);
    
    /* 22, 31, 33, 41 and 55 fitted directly + 44 dA */
    eob_nqc_point_HM(dyn, max_A, max_dA, max_omg, max_domg);

    /* 21 fitted directly at tpeak_22*/
    eob_nqc_point_HM_peak22(dyn, max_A, max_dA, max_omg, max_domg);

    /* Over-writing fits using postpeak quantities for modes in kpostpeak */
    for (int j=0; j<kpostpeak_size; j++) {
      int k = kpostpeak[j];
      
      /* Normalizing c1A and c4A */
      int l = LINDEX[k];
      c1A[k] /= sqrt((l+2)*(l+1)*l*(l-1));
      c4A[k] /= sqrt((l+2)*(l+1)*l*(l-1));
      
      eob_nqc_point_postpeak(Mbh,c1A[k],c2A[k],c3A[k],c4A[k],
            c1phi[k],c2phi[k],c3phi[k],c4phi[k],alpha1[k],omega1[k],
            &A_tmp,&dA_tmp,&omg_tmp,&domg_tmp);
      
      max_A[k]    = A_tmp;
      max_dA[k]   = dA_tmp;
      max_omg[k]  = omg_tmp;
      max_domg[k] = domg_tmp;
    }
  } else {
    
    eob_nqc_point(dyn, &A_tmp, &dA_tmp, &omg_tmp, &domg_tmp);
	
    /** Switch on the 22 values (only) */
    max_A[1]    = A_tmp;
    max_dA[1]   = dA_tmp;
    max_omg[1]  = omg_tmp;
    max_domg[1] = domg_tmp;
  
  }

  if (VERBOSE) {
    printf("NR values for NQC determination:\n");
    PRFORMd("A22_mrg",max_A[1]);
    PRFORMd("dA22_mrg",max_dA[1]);
    PRFORMd("omg22_mrg",max_omg[1]);
    PRFORMd("domg22_mrg",max_domg[1]);
  }
  
  /** NQC corrections to AMPLITUDE (n1,n2,n3) and PHASE (n4,n5,n6)
   * NQC basis for (2,2) waveform : AMPLITUDE
   * note: n3 and n6 are not used
   */
  double pr_star2, r2, w2;
  for (int k=0; k<KMAX; k++) {   
    if(h->kmask_nqc[k]) {
      for (int j=0; j<size; j++) {
        pr_star2 = SQ(pr_star[j]);
        r2       = SQ(r[j]);
        w2       = SQ(w[j]); //CHECKME: Omg or Omg_orbital ?
        n1[k][j]  = pr_star2/(r2*w2);         /* [pr*\/(r Omg)]^2 */
        n2[k][j]  = ddotr[j]/(r[j]*w2);       /* [ddot{r}/(r Omg^2)] */
        n4[k][j]  = pr_star[j]/(r[j]*w[j]);   /* pr*\/(r Omg) */
        n5[k][j]  = n4[k][j]*r2*w2;              /* (pr*)*(r Omg) */
      }
    }
  }

  if (EOBPars->use_flm == USEFLM_HM || EOBPars->use_flm == USEFLM_HM_4PN22) {
    for (int j=0; j<size; j++) {
      /* l=2,m=1 */
      if(h->kmask_nqc[0]) {
      n2[0][j] = cbrt(SQ(w[j]))*n1[0][j];
      n5[0][j] = cbrt(SQ(w[j]))*n4[0][j];
      }
      /* l=2,m=2 */
      if(h->kmask_nqc[1] && EOBPars->use_flm == USEFLM_HM_4PN22) {
        n2[1][j] = SQ(pr_star[j])*n1[1][j];
      }
      /* l=3, l=4 & l=5 */
      for (int k=2; k<14; k++) {   
	      if(h->kmask_nqc[k]) {
	        n5[k][j]  = cbrt(SQ(w[j]))*n4[k][j];
	      }
      }
    }
  }  

#if (DEBUG)
  FILE* fp_dbg;
  fp_dbg = fopen("nqc_nfunc.txt", "w");
  for (int j=0; j<size; j++) {
    fprintf(fp_dbg, "%20.12f\t%.16e\t%.16e\t%.16e\t%.16e\n", t[j], n1[1][j], n2[1][j], n4[1][j], n5[1][j]);
  }
  fclose(fp_dbg);
#endif    
  
  /** Derivatives for the phase */
  for (int k=0; k<KMAX; k++) {
    if(h->kmask_nqc[k]) {
      D0(n4[k],dt,size,d_n4[k]);
      D0(n5[k],dt,size,d_n5[k]);
      D0(d_n4[k],dt,size,d2_n4[k]);
      D0(d_n5[k],dt,size,d2_n5[k]);
    }
  }

#if (DEBUG)
  fp_dbg = fopen("nqc_dfunc.txt", "w");
  for (int j=0; j<size; j++) {
    fprintf(fp_dbg, "%f\t%.16e\t%.16e\t%.16e\t%.16e\n", t[j], d_n4[1][j], d_n5[1][j], d2_n4[1][j], d2_n5[1][j]);  
  }
  fclose(fp_dbg);
#endif    

  /** Find max Omg */
  int Omgmax_index = 0;
  double Omg_max   = Omg_orb[0];
  for (int j=0; j<size; j++) {
    if (Omg_orb[j] > Omg_max) {
      Omg_max = Omg_orb[j];
      Omgmax_index = j;
    }
  }
  //TODO:
  //Test the search backwards
  /*
    int Omgmax_index = size-1;
    double Omg_max   = Omg_orb[Omgmax_index];
    for (int j=(size-2); j--; ) {
    if (Omg_orb[j] < Omg_max) 
    break;
    Omg_max = Omg_orb[j];
    Omgmax_index = j;
    }
  */

  /** Time */
  double tOmgOrb_pk = t[Omgmax_index];
  double DeltaT_nqc = eob_nqc_timeshift(nu, chi1);
  double tNQC = tOmgOrb_pk - DeltaT_nqc;

  if (VERBOSE) {
    printf("NQC info:\n");
    PRFORMd("DeltaT_tNQC",DeltaT_nqc);
    PRFORMd("tNQC[bare]",tNQC);
  }

  /** Find jmax: t[jmax] <= tNQC */
  double tmrg[KMAX];
  tmrg[1] = tNQC - 2.;
	
  int jmax = 0;
  for (int j=0; j<size; j++) {
    if(t[j] > tNQC) {
      jmax = j-2;
      break;
    }
  }
	     
  double dtmrg[KMAX];
  double t_NQC[KMAX];
  int    j_NQC[KMAX];

  if (EOBPars->use_flm == USEFLM_HM || EOBPars->use_flm == USEFLM_HM_4PN22) {
    int modesatpeak22[KMAX]; 
    set_multipolar_idx_mask (modesatpeak22, KMAX, EOBPars->knqcpeak22, EOBPars->knqcpeak22_size, 0);
    eob_nqc_deltat_lm(dyn, dtmrg);
    
    for (int k=0; k<KMAX; k++) {   
      if(h->kmask_nqc[k]) {
	      tmrg[k]  = tmrg[1] + dtmrg[k];
	      t_NQC[k] = tmrg[k] + 2.;
	      if(modesatpeak22[k]) t_NQC[k] = tmrg[1];
	      j_NQC[k] = size-1;
	      for (int j=size-2; j>=0; j--) {
	        if(t[j] < t_NQC[k]) {
	          break;
	        }
	      j_NQC[k] = j;
	      }
      }
    }
  }
  
  /** Solve the linear systems */
  
  /* Regge-Wheeler-Zerilli normalized amplitude. 
     The ringdown coefficient refer to this normalization.
     Nagar & Rezzolla, CQG 22 (2005) R167 */      
  for (int k=0; k<KMAX; k++) {
    if(h->kmask_nqc[k]){  
      double nlm = 1./(sqrt( (LINDEX[k]+2)*(LINDEX[k]+1)*LINDEX[k]*(LINDEX[k]-1) ) );
      if (h->ampli[k][0] > 0.) {	
	nNegAmp[k] = 0;	
      } else {		
	nNegAmp[k] = 1;	
      }
      for (int j=0; j<size; j++) {
	p1tmp[k][j] = fabs(h->ampli[k][j] * nlm);
      }
    }
  }
  
  /* Matrix elements: waveform amplitude at all points */
  for (int k=0; k<KMAX; k++) {
    if(h->kmask_nqc[k]) {
      for (int j=0; j<size; j++) {
	m11[k][j] = n1[k][j] * p1tmp[k][j];
	m12[k][j] = n2[k][j] * p1tmp[k][j];
      }
    }
  }
  
  /* Take FD derivatives */
  for (int k=0; k<KMAX; k++) {
    if(h->kmask_nqc[k]) {
      D0(m11[k],dt,size, m21[k]);
      D0(m12[k],dt,size, m22[k]);
      D0(p1tmp[k],dt,size, p2tmp[k]);
    }
  }

#if (DEBUG)
  fp_dbg = fopen("nqc_amp_func.txt", "w");
  for (int j=0; j<size; j++) {
    fprintf(fp_dbg, "%e\t%e\t%e\n", t[j], p1tmp[1][j], p2tmp[1][j]);
  }
  fclose(fp_dbg);  
#endif
  
  double detM = 1.;
  double oodetM = 1.;
  for (int k=0; k<KMAX; k++) {
    if(h->kmask_nqc[k]) {  
      
      if (EOBPars->use_flm == USEFLM_HM || EOBPars->use_flm == USEFLM_HM_4PN22) {
	jmax = j_NQC[k];
      }
            
      /* Computation of ai coefficients at Omega peak */
      P[0]     = max_A[k]  - p1tmp[k][jmax];
      P[1]     = max_dA[k] - p2tmp[k][jmax];
      
      M[0]     = m11[k][jmax];
      M[1]     = m12[k][jmax];
      M[2]     = m21[k][jmax];
      M[3]     = m22[k][jmax];
      
      /* detM     = M[0]*M[3]-M[1]*M[2]; 
	 ai[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
	 ai[k][1] = (M[0]*P[1] - M[2]*P[0])/detM; */
      /* safe version (amplitude can be zero) */
      oodetM   = 1.0/(M[0]*M[3]-M[1]*M[2]);
      if (isfinite(oodetM)) {
	ai[k][0] = (M[3]*P[0] - M[1]*P[1])*oodetM;
	ai[k][1] = (M[0]*P[1] - M[2]*P[0])*oodetM;
      }
      
      /* Computation of bi coefficients at Omega peak */
      P[0]     = omg[k][jmax]   - max_omg[k];
      P[1]     = domg[k][jmax]  - max_domg[k];
      
      M[0]     = d_n4[k][jmax];
      M[1]     = d_n5[k][jmax];
      M[2]     = d2_n4[k][jmax];
      M[3]     = d2_n5[k][jmax];
      
      /* detM     =  M[0]*M[3] - M[1]*M[2];
	 bi[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
	 bi[k][1] = (M[0]*P[1] - M[2]*P[0])/detM; */
      /* safe version (phase can be zero) */
      oodetM   = 1.0/(M[0]*M[3]-M[1]*M[2]);
      if (isfinite(oodetM)) {
	bi[k][0] = (M[3]*P[0] - M[1]*P[1])*oodetM;
	bi[k][1] = (M[0]*P[1] - M[2]*P[0])*oodetM;
      }
    }
  }
  
  if (VERBOSE){
    printf("NQC coefficients for 22 mode:\n");
    PRFORMd("a1",ai[1][0]);
    PRFORMd("a2",ai[1][1]);
    PRFORMd("b1",bi[1][0]);
    PRFORMd("b2",bi[1][1]);
  }
  
  /** Set amplitude and phase */
  for (int k=0; k<KMAX; k++) {
    if(h->kmask_nqc[k]) {  
      for (int j=0; j<size; j++) {
	hnqc->ampli[k][j] = 1. + ai[k][0]*n1[k][j] + ai[k][1]*n2[k][j];
	hnqc->phase[k][j] =      bi[k][0]*n4[k][j] + bi[k][1]*n5[k][j];
      }
    }
  }
  
  /** Multiply waveform to NQC */
  for (int k=0; k<KMAX; k++) {
    if(h->kmask_nqc[k]){    
      for (int j=0; j<size; j++) {
	h->ampli[k][j] *= hnqc->ampli[k][j];
	h->phase[k][j] -= hnqc->phase[k][j];
      }
    }
  }
  
  if (EOBPars->output_nqc_coefs) {
    /** Output the NQC coefficients */
    FILE* fp;
    char fname[STRLEN];
    strcpy(fname, EOBPars->output_dir);
    strcat(fname, "/nqc_coefs.txt");
    fp = fopen(fname, "w");
    fprintf(fp, "# q=%e chizA=%e chizB=%e f0=%e\n",EOBPars->q,EOBPars->chi1,EOBPars->chi2,EOBPars->initial_frequency);
    fprintf(fp, "# M=%e LambdaA=[%e,%e,%e] LambdaBl2=[%e,%e,%e]\n",EOBPars->M,
	    EOBPars->LambdaAl2, EOBPars->LambdaAl3,EOBPars->LambdaAl4,
	    EOBPars->LambdaBl2,EOBPars->LambdaBl3,EOBPars->LambdaBl4);
    for (int k=0; k<KMAX; k++) {
      if(h->kmask_nqc[k]){
        fprintf(fp, "%d %d %d %e %e %e %e\n", k, LINDEX[k], MINDEX[k], 
          ai[k][0], ai[k][1], 
          bi[k][0], bi[k][1]);
      }
    }  
    fclose(fp);  
  }
  
  /** Free mem */
  for (int k=0; k<KMAX; k++) {
    free(omg[k]);
    free(domg[k]);
    free(m11[k]);
    free(m12[k]);
    free(m13[k]);
    free(m21[k]);
    free(m22[k]);
    free(p1tmp[k]);
    free(p2tmp[k]);
    free(n1[k]);
    free(n2[k]);
    free(n4[k]);
    free(n5[k]);
    free(d_n4[k]);
    free(d_n5[k]);
    free(d2_n4[k]);
    free(d2_n5[k]);
  }
  
}

/**
 * Function: eob_wav_hlmNQC_find_a1a2a3_mrg_HM
 * -------------------------------------------
 *   Computes the factors and the coefficients that build the
 *   NQC corrections to the waveform in the eccentric case.
 *   This routine works around merger with dyn_mrg and hlm_mrg and
 *   then add everything also to hlm.
 *   With respect to the quasi-circ case, here we also consider the sigmoid
 * 
 *   @param[in] dyn_mrg : dynamics around merger
 *   @param[in] hlm_mrg : multipolar waveform around merger (modified by this function)
 *   @param[in] hnqc    : NQC part of the waveform
 *   @param[in] dyn     : full dynamics
 *   @param[in] hlm     : full multipolar waveform (modified by this function)
*/
void eob_wav_hlmNQC_find_a1a2a3_mrg_ecc(Dynamics *dyn_mrg, Waveform_lm *hlm_mrg, Waveform_lm *hnqc,
				       Dynamics *dyn, Waveform_lm *hlm)
{  
  double A_tmp, dA_tmp, omg_tmp, domg_tmp;
  double alpha1[KMAX], omega1[KMAX];
  double c1A[KMAX], c2A[KMAX], c3A[KMAX], c4A[KMAX];
  double c1phi[KMAX], c2phi[KMAX], c3phi[KMAX], c4phi[KMAX];
	  
  const double nu   = EOBPars->nu;
  const double chi1 = EOBPars->chi1;
  const double chi2 = EOBPars->chi2;
  const double X1   = EOBPars->X1;
  const double X2   = EOBPars->X2;
  const double aK   = EOBPars->a1+EOBPars->a2;
  const double Mbh  = EOBPars->Mbhf;
  const double abh  = EOBPars->abhf;
    
  double *t       = hlm_mrg->time;
  double *r       = dyn_mrg->data[EOB_RAD];
  double *w       = dyn_mrg->data[EOB_MOMG]; /* Omega */
  double *pph     = dyn_mrg->data[EOB_PPHI];
  double *pr_star = dyn_mrg->data[EOB_PRSTAR];
  double *Omg_orb = dyn_mrg->data[EOB_OMGORB]; /* Omega orbital */
  double *ddotr   = dyn_mrg->data[EOB_DDOTR];
  
  double P[2], M[4];   
  double max_A[KMAX],max_dA[KMAX],d2max[KMAX],d3max[KMAX],max_omg[KMAX],max_domg[KMAX],maxd2omg[KMAX], DeltaT[KMAX];
  double ai[KMAX][2] = {0.};
  double bi[KMAX][2] = {0.};
  
  const int size = hlm_mrg->size;
  for (int i = 0; i < size; i++) {
    hnqc->time[i] = t[i];
  }
  
  double *omg[KMAX], *domg[KMAX];
  double *n1[KMAX],*n2[KMAX],*n3[KMAX],*n4[KMAX],*n5[KMAX],*n6[KMAX],
    *d_n4[KMAX],*d_n5[KMAX],*d_n6[KMAX],*d2_n4[KMAX],*d2_n5[KMAX],*d2_n6[KMAX];
  double *m11[KMAX], *m12[KMAX], *m13[KMAX], *m21[KMAX], *m22[KMAX];
  double *p1tmp[KMAX], *p2tmp[KMAX]; /* RWZ amplitude and derivative */
  
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask[k]){
      omg[k]  = (double*) calloc (size,sizeof(double));
      domg[k] = (double*) calloc (size,sizeof(double));
      m11[k] = (double*) calloc (size,sizeof(double));
      m12[k] = (double*) calloc (size,sizeof(double));
      m13[k] = (double*) calloc (size,sizeof(double));
      m21[k] = (double*) calloc (size,sizeof(double));
      m22[k] = (double*) calloc (size,sizeof(double));
      p1tmp[k] = (double*) calloc (size,sizeof(double));
      p2tmp[k] = (double*) calloc (size,sizeof(double));
      n1[k] = (double*) calloc (size,sizeof(double));
      n2[k] = (double*) calloc (size,sizeof(double));
      n4[k] = (double*) calloc (size,sizeof(double));
      n5[k] = (double*) calloc (size,sizeof(double));
      d_n4[k] = (double*) calloc (size,sizeof(double));
      d_n5[k] = (double*) calloc (size,sizeof(double));
      d2_n4[k] = (double*) calloc (size,sizeof(double));
      d2_n5[k] = (double*) calloc (size,sizeof(double));
    }
  }

  /** omega derivatives */
  const double dt = t[1]-t[0];
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){
      D0(hlm_mrg->phase[k], dt, size, omg[k]);
      D0(omg[k], dt, size, domg[k]);
    }
  }
  
  /** NR fits */
  for (int k=0; k<KMAX; k++) {   
    max_A[k]    = 0.;
    max_dA[k]   = 0.;
    max_omg[k]  = 0.;
    max_domg[k] = 0.;
  }

  /* Higher modes */
  /* Choosing modes using kpostpeak array */
  int kpostpeak_size = EOBPars->kpostpeak_size;  
  int *kpostpeak     = EOBPars->kpostpeak;

  QNMHybridFitCab_HM(nu, X1, X2, chi1, chi2, aK,  Mbh, abh,  
          c1A, c2A, c3A, c4A, c1phi, c2phi, c3phi, c4phi,
          alpha1, omega1);
  
  /* 22, 31, 33, 41 and 55 fitted directly + 44 dA */
  eob_nqc_point_HM(dyn, max_A, max_dA, max_omg, max_domg);

  /* 21 fitted directly at tpeak_22*/
  eob_nqc_point_HM_peak22(dyn, max_A, max_dA, max_omg, max_domg);

  /* Over-writing fits using postpeak quantities for modes in kpostpeak */
  for (int j=0; j<kpostpeak_size; j++) {
    int k = kpostpeak[j];
    
    /* Normalizing c1A and c4A */
    int l = LINDEX[k];
    c1A[k] /= sqrt((l+2)*(l+1)*l*(l-1));
    c4A[k] /= sqrt((l+2)*(l+1)*l*(l-1));
    
    eob_nqc_point_postpeak(Mbh,c1A[k],c2A[k],c3A[k],c4A[k],
          c1phi[k],c2phi[k],c3phi[k],c4phi[k],alpha1[k],omega1[k],
          &A_tmp,&dA_tmp,&omg_tmp,&domg_tmp);
    
    max_A[k]    = A_tmp;
    max_dA[k]   = dA_tmp;
    max_omg[k]  = omg_tmp;
    max_domg[k] = domg_tmp;
    
  }

  if (VERBOSE) {
    printf("NR values for NQC determination:\n");
    PRFORMd("A22_mrg",max_A[1]);
    PRFORMd("dA22_mrg",max_dA[1]);
    PRFORMd("omg22_mrg",max_omg[1]);
    PRFORMd("domg22_mrg",max_domg[1]);
  }

  /** NQC corrections to AMPLITUDE (n1,n2,n3) and PHASE (n4,n5,n6)
   * NQC basis for (2,2) waveform : AMPLITUDE
   * note: n3 and n6 are not used
   */
  double pr_star2, r2, w2;
  for (int k=0; k<KMAX; k++) {   
    if(hlm_mrg->kmask_nqc[k]){ 
      for (int j=0; j<size; j++) {
        pr_star2 = SQ(pr_star[j]);
        r2       = SQ(r[j]);
        w2       = SQ(w[j]); //CHECKME: Omg or Omg_orbital ?
        n1[k][j]  = pr_star2/(r2*w2);         /* [pr*\/(r Omg)]^2 */
        n2[k][j]  = ddotr[j]/(r[j]*w2);       /* [ddot{r}/(r Omg^2)] */
        //n3[k][j]  = n1[k][j]*pr_star2;
        n4[k][j]  = pr_star[j]/(r[j]*w[j]);   /* pr*\/(r Omg) */
        n5[k][j]  = n4[k][j]*r2*w2;              /* (pr*)*(r Omg) */
        //n6[k][j]  = n5[k][j]*pr_star2;
      }
    }
  }

  for (int j=0; j<size; j++) {
    /* l=2,m=1 */
    if(hlm_mrg->kmask_nqc[0]){
      n2[0][j] = cbrt(SQ(w[j]))*n1[0][j];
      n5[0][j] = cbrt(SQ(w[j]))*n4[0][j];
    }
    /* l=2,m=2 */
    if(hlm_mrg->kmask_nqc[1] && EOBPars->use_flm == USEFLM_HM_4PN22) {
      n2[1][j] = SQ(pr_star[j])*n1[1][j];
    }

    /* l=3, l=4 & l=5 */
    // FIXME: condition on spins bad for non-spinning limit
    if ((chi1 <= 0.) || (chi2 <= 0.)) {
      for (int k=2; k<14; k++) {
        if(hlm_mrg->kmask_nqc[k]){
          n5[k][j]  = cbrt(SQ(w[j]))*n4[k][j];
        }
      }
    }
    /* l=3, m=2 */
    if(hlm_mrg->kmask_nqc[3]){
      n2[3][j] = cbrt(SQ(w[j]))*n1[3][j];
      n5[3][j] = cbrt(SQ(w[j]))*n4[3][j];
    }
    /* l=4, m=2 */
    if(hlm_mrg->kmask_nqc[6]){
      n2[6][j] = cbrt(SQ(w[j]))*n1[6][j];
      n5[6][j] = cbrt(SQ(w[j]))*n4[6][j];
    }
    /* l=4, m=3 */
    if(hlm_mrg->kmask_nqc[7]){
      n2[7][j] = cbrt(SQ(w[j]))*n1[7][j];
      n5[7][j] = cbrt(SQ(w[j]))*n4[7][j];
    }
    /* l=5, m=5 */
    if(hlm_mrg->kmask_nqc[13]){
      n5[13][j] = cbrt(SQ(w[j]))*n4[13][j];
    }      
  }
      
#if (DEBUG)
  FILE* fp_dbg;
  fp_dbg = fopen("nqc_nfunc.txt", "w");
  for (int j=0; j<size; j++) {
    fprintf(fp_dbg, "%20.12f\t%.16e\t%.16e\t%.16e\t%.16e\n", t[j], n1[1][j], n2[1][j], n4[1][j], n5[1][j]);
  }
  fclose(fp_dbg);
#endif    
      
  /** Derivatives for the phase */
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){   
      D0(n4[k],dt,size, d_n4[k]);
      D0(n5[k],dt,size, d_n5[k]);
      D0(d_n4[k],dt,size, d2_n4[k]);
      D0(d_n5[k],dt,size, d2_n5[k]);
    }
  }
    
#if (DEBUG)
  fp_dbg = fopen("nqc_dfunc.txt", "w");
  for (int j=0; j<size; j++) {
    fprintf(fp_dbg, "%f\t%.16e\t%.16e\t%.16e\t%.16e\n", t[j], d_n4[1][j], d_n5[1][j], d2_n4[1][j], d2_n5[1][j]);  
  }
  fclose(fp_dbg);
#endif    

  /** Find max Omg */
  int Omgmax_index = size-1;
  double Omg_max = Omg_orb[Omgmax_index];
  for (int j = size-2; j-- ; ) {
    if (Omg_orb[j] < Omg_max) 
      break;
    Omgmax_index = j;
    Omg_max = Omg_orb[j]; 
  }
  
  /** Time */
  double tOmgOrb_pk = t[Omgmax_index];
  double DeltaT_nqc = eob_nqc_timeshift(nu, chi1);
  double tNQC = tOmgOrb_pk - DeltaT_nqc;

  if (VERBOSE) {
    printf("NQC info:\n");
    PRFORMd("DeltaT_tNQC",DeltaT_nqc);
    PRFORMd("tNQC[bare]",tNQC);
  }

  /** Find jmax: t[jmax] <= tNQC */
  double tmrg[KMAX];
  tmrg[1] = tNQC - 2.;
	
  int jmax = 0;
  for (int j=0; j<size; j++) {
    if(t[j] > tNQC) {
      jmax = j-2;
      break;
    }
  }
  
  double dtmrg[KMAX];
  double t_NQC[KMAX];
  int    j_NQC[KMAX];

  int modesatpeak22[KMAX]; 
  set_multipolar_idx_mask (modesatpeak22, KMAX, EOBPars->knqcpeak22, EOBPars->knqcpeak22_size, 0);  
  eob_nqc_deltat_lm(dyn, dtmrg);
  
  for (int k=0; k<KMAX; k++) {   
    if(hlm_mrg->kmask_nqc[k]){
      tmrg[k]  = tmrg[1] + dtmrg[k];
      t_NQC[k] = tmrg[k] + 2.;
      if(modesatpeak22[k]) t_NQC[k] = tmrg[1];

      j_NQC[k] = size-1;
      for (int j=size-2; j>=0; j--) {
	      if(t[j] < t_NQC[k]) {
	        break;
	      }
	    j_NQC[k] = j;
      } 
    }
  }
  
  /** Solve the linear systems */
  
  /* Regge-Wheeler-Zerilli normalized amplitude. 
     The ringdown coefficient refer to this normalization.
     Nagar & Rezzolla, CQG 22 (2005) R167 */      
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){
      double nlm = 1./(sqrt( (LINDEX[k]+2)*(LINDEX[k]+1)*LINDEX[k]*(LINDEX[k]-1) ) );

      if (hlm->ampli[k][0] > 0.) {	
        nNegAmp[k] = 0;
      } else {	
        nNegAmp[k] = 1;	
      }
      
      for (int j=0; j<size; j++) {
        p1tmp[k][j] = fabs(hlm_mrg->ampli[k][j] * nlm);      
      }
    }
  }
  
  /* Matrix elements: waveform amplitude at all points */
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){
      for (int j=0; j<size; j++) {
        m11[k][j] = n1[k][j] * p1tmp[k][j];
        m12[k][j] = n2[k][j] * p1tmp[k][j];
      }
    }
  }

  /* Take FD derivatives */
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){
      D0(m11[k],dt,size, m21[k]);
      D0(m12[k],dt,size, m22[k]);
      D0(p1tmp[k],dt,size, p2tmp[k]);
    }
  }

#if (DEBUG)
  fp_dbg = fopen("nqc_amp_func.txt", "w");
  for (int j=0; j<size; j++) {
    fprintf(fp_dbg, "%e\t%e\t%e\n", t[j], p1tmp[1][j], p2tmp[1][j]);
  }
  fclose(fp_dbg);  
#endif

  double detM = 1.;
  double oodetM = 1.;
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){

      jmax = j_NQC[k];
      

      /* Computation of ai coefficients at Omega peak */
      P[0]     = max_A[k]  - p1tmp[k][jmax];
      P[1]     = max_dA[k] - p2tmp[k][jmax];
      
      M[0]     = m11[k][jmax];
      M[1]     = m12[k][jmax];
      M[2]     = m21[k][jmax];
      M[3]     = m22[k][jmax];
      
      /* detM     = M[0]*M[3]-M[1]*M[2]; 
        ai[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
        ai[k][1] = (M[0]*P[1] - M[2]*P[0])/detM; */
      /* safe version (amplitude can be zero) */
      oodetM   = 1.0/(M[0]*M[3]-M[1]*M[2]);
      if (isfinite(oodetM)) {
        ai[k][0] = (M[3]*P[0] - M[1]*P[1])*oodetM;
        ai[k][1] = (M[0]*P[1] - M[2]*P[0])*oodetM;
      }

      /* Computation of bi coefficients at Omega peak */
      P[0]     = omg[k][jmax]   - max_omg[k];
      P[1]     = domg[k][jmax]  - max_domg[k];
      
      M[0]     = d_n4[k][jmax];
      M[1]     = d_n5[k][jmax];
      M[2]     = d2_n4[k][jmax];
      M[3]     = d2_n5[k][jmax];

      /* detM     =  M[0]*M[3] - M[1]*M[2];
        bi[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
        bi[k][1] = (M[0]*P[1] - M[2]*P[0])/detM; */
      /* safe version (phase can be zero) */
      oodetM   = 1.0/(M[0]*M[3]-M[1]*M[2]);
      if (isfinite(oodetM)) {
        bi[k][0] = (M[3]*P[0] - M[1]*P[1])*oodetM;
        bi[k][1] = (M[0]*P[1] - M[2]*P[0])*oodetM;
      }
    }

  }
  
  if (VERBOSE){
    printf("NQC coefficients for 22 mode:\n");
    PRFORMd("a1",ai[1][0]);
    PRFORMd("a2",ai[1][1]);
    PRFORMd("b1",bi[1][0]);
    PRFORMd("b2",bi[1][1]);
  }

  /* Defining sigmoid function to switch on NQCs near the end of the evolution */
  /*
    Old configuration used for arXiv:2001.11736
    double t0 = tNQC - 30.;
    double alpha = 0.09;
  */
  double t0 = tNQC - EOBPars->delta_t0_sigmoid_NQC;
  double alpha =EOBPars->alpha_sigmoid_NQC;
  double *sigmoid;
  sigmoid = (double*) calloc (size, sizeof(double));
  double fact[KMAX];
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){
      // Correction needed to have sigmoid = 1 when attaching ringdown
      fact[k] = 1. + exp(-alpha*(t[j_NQC[k]] - t0));
    }
  }
    
  for (int j=0; j<size; j++) {
    sigmoid[j] = 1./(1. + exp(-alpha*(t[j] - t0)));
    for (int k=0; k<KMAX; k++) {
      if(hlm_mrg->kmask_nqc[k]){
        n1[k][j] = n1[k][j]*sigmoid[j]*fact[k];
        n2[k][j] = n2[k][j]*sigmoid[j]*fact[k];
        n4[k][j] = n4[k][j]*sigmoid[j]*fact[k];
        n5[k][j] = n5[k][j]*sigmoid[j]*fact[k];
      }
    }
  }
    
  free(sigmoid);
  
  /** Set amplitude and phase */
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){
      for (int j=0; j<size; j++) {
        hnqc->ampli[k][j] = 1. + ai[k][0]*n1[k][j] + ai[k][1]*n2[k][j];
        hnqc->phase[k][j] = bi[k][0]*n4[k][j] + bi[k][1]*n5[k][j];
      }
    }
  }
  
  /** Multiply merger waveform to NQC */
  for (int k=0; k<KMAX; k++) {
    for (int j=0; j<size; j++) {
      if(hlm_mrg->kmask_nqc[k]){
        hlm_mrg->ampli[k][j] *= hnqc->ampli[k][j];
        hlm_mrg->phase[k][j] -= hnqc->phase[k][j];
      }
    }
  }

  /** Multiply full waveform to NQC */
  r       = dyn->data[EOB_RAD];
  w       = dyn->data[EOB_MOMG]; /* Omega */
  pph     = dyn->data[EOB_PPHI];
  pr_star = dyn->data[EOB_PRSTAR];
  Omg_orb = dyn->data[EOB_OMGORB]; /* Omega orbital */
  ddotr   = dyn->data[EOB_DDOTR];

  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask[k]){
      free(n1[k]);
      free(n2[k]);
      free(n4[k]);
      free(n5[k]);
    }
  }

  const int fullsize = hlm->size;
  
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask[k]){
      n1[k] = (double*) calloc (fullsize,sizeof(double));
      n2[k] = (double*) calloc (fullsize,sizeof(double));
      n4[k] = (double*) calloc (fullsize,sizeof(double));
      n5[k] = (double*) calloc (fullsize,sizeof(double));
    }
  }
	
  for (int k=0; k<KMAX; k++) {   
    if(hlm_mrg->kmask_nqc[k]){   
      for (int j=0; j<fullsize; j++) {
        pr_star2 = SQ(pr_star[j]);
        r2       = SQ(r[j]);
        w2       = SQ(w[j]);                  //CHECKME: Omg or Omg_orbital ?
        n1[k][j]  = pr_star2/(r2*w2);         /* [pr*\/(r Omg)]^2 */
        n2[k][j]  = ddotr[j]/(r[j]*w2);       /* [ddot{r}/(r Omg^2)] */
        n4[k][j]  = pr_star[j]/(r[j]*w[j]);   /* pr*\/(r Omg) */
        n5[k][j]  = n4[k][j]*r2*w2;           /* (pr*)*(r Omg) */
      }
    }
  }

  for (int j=0; j<fullsize; j++) {
    /* l=2,m=1 */
    if(hlm->kmask_nqc[0]){
      n2[0][j] = cbrt(SQ(w[j]))*n1[0][j];
      n5[0][j] = cbrt(SQ(w[j]))*n4[0][j];
    }
    /* l=2,m=2 */
    if(hlm->kmask_nqc[1] && EOBPars->use_flm == USEFLM_HM_4PN22) {
      n2[1][j] = SQ(pr_star[j])*n1[1][j];
    }

    /* l=3, l=4 & l=5 */
    // FIXME: condition on spins bad for non-spinning limit
    if ((chi1 <= 0.) || (chi2 <= 0.)) {
      for (int k=2; k<14; k++) {
	      if(hlm->kmask_nqc[k]){
	        n5[k][j]  = cbrt(SQ(w[j]))*n4[k][j];
	      }
      }
    }
    /* l=3, m=2 */
    if(hlm->kmask_nqc[3]){
      n2[3][j] = cbrt(SQ(w[j]))*n1[3][j];
      n5[3][j] = cbrt(SQ(w[j]))*n4[3][j];
    }
    /* l=4, m=2 */
    if(hlm->kmask_nqc[6]){
      n2[6][j] = cbrt(SQ(w[j]))*n1[6][j];
      n5[6][j] = cbrt(SQ(w[j]))*n4[6][j];
    }
    /* l=4, m=3 */
    if(hlm->kmask_nqc[7]){
      n2[7][j] = cbrt(SQ(w[j]))*n1[7][j];
      n5[7][j] = cbrt(SQ(w[j]))*n4[7][j];
    }
    /* l=5, m=5 */
    if(hlm->kmask_nqc[13]){
      n5[13][j] = cbrt(SQ(w[j]))*n4[13][j];
    }      
  }

  /* Defining sigmoid function to switch on NQCs near the end of the evolution */
  sigmoid = (double*) calloc (fullsize, sizeof(double));
    
  for (int j=0; j<fullsize; j++) {
    sigmoid[j] = 1./(1. + exp(-alpha*(hlm->time[j] - t0)));
    for (int k=0; k<KMAX; k++) {
      if(hlm->kmask_nqc[k]){
        n1[k][j] = n1[k][j]*sigmoid[j]*fact[k];
        n2[k][j] = n2[k][j]*sigmoid[j]*fact[k];
        n4[k][j] = n4[k][j]*sigmoid[j]*fact[k];
        n5[k][j] = n5[k][j]*sigmoid[j]*fact[k];
      }
    }
  }
  free(sigmoid);
  
  for (int k=0; k<KMAX; k++) {
    if(hlm->kmask_nqc[k]){
      for (int j=0; j<fullsize; j++) {
	      hlm->ampli[k][j] *= (1. + ai[k][0]*n1[k][j] + ai[k][1]*n2[k][j]);
        hlm->phase[k][j] -= (bi[k][0]*n4[k][j] + bi[k][1]*n5[k][j]);
      }
    }
  }

  if (EOBPars->output_nqc_coefs) {
    /** Output the NQC coefficients */
    FILE* fp;
    char fname[STRLEN];
    strcpy(fname, EOBPars->output_dir);
    strcat(fname, "/nqc_coefs.txt");
    fp = fopen(fname, "w");
    fprintf(fp, "# q=%e chizA=%e chizB=%e f0=%e\n",EOBPars->q,EOBPars->chi1,EOBPars->chi2,EOBPars->initial_frequency);
    fprintf(fp, "# M=%e LambdaA=[%e,%e,%e] LambdaBl2=[%e,%e,%e]\n",EOBPars->M,
	    EOBPars->LambdaAl2,EOBPars->LambdaAl3,EOBPars->LambdaAl4,
	    EOBPars->LambdaBl2,EOBPars->LambdaBl3,EOBPars->LambdaBl4);
    for (int k=0; k<KMAX; k++) {
      if(hlm->kmask_nqc[k]){
        fprintf(fp, "%d %d %d %e %e %e %e\n", k, LINDEX[k], MINDEX[k], 
		ai[k][0], ai[k][1], 
		bi[k][0], bi[k][1]);
      }
    }  
    fclose(fp);  
  }
  
  /** Free mem */
  for (int k=0; k<KMAX; k++) {
    if(hlm->kmask[k]){
      free(omg[k]);
      free(domg[k]);
      free(m11[k]);
      free(m12[k]);
      free(m13[k]);
      free(m21[k]);
      free(m22[k]);
      free(p1tmp[k]);
      free(p2tmp[k]);
      free(n1[k]);
      free(n2[k]);
      free(n4[k]);
      free(n5[k]);
      free(d_n4[k]);
      free(d_n5[k]);
      free(d2_n4[k]);
      free(d2_n5[k]);
    }
  }

}

/**
 * Function: eob_wav_hlmNQC_find_a1a2a3_mrg_HM
 * -------------------------------------------
 *   Computes the factors and the coefficients that build the
 *   NQC corrections to the waveform in the spinning case.
 *   This routine works around merger with dyn_mrg and hlm_mrg and
 *   then add everything also to hlm
 *   
 * 
 *   @param[in] dyn_mrg : dynamics around merger
 *   @param[in] hlm_mrg : multipolar waveform around merger (modified by this function)
 *   @param[in] hnqc    : NQC part of the waveform
 *   @param[in] dyn     : full dynamics
 *   @param[in] hlm     : full multipolar waveform (modified by this function)
*/
void eob_wav_hlmNQC_find_a1a2a3_mrg_HM(Dynamics *dyn_mrg, Waveform_lm *hlm_mrg, Waveform_lm *hnqc,
				       Dynamics *dyn, Waveform_lm *hlm)
{
  double A_tmp, dA_tmp, omg_tmp, domg_tmp;
  double alpha1[KMAX], omega1[KMAX];
  double c1A[KMAX], c2A[KMAX], c3A[KMAX], c4A[KMAX];
  double c1phi[KMAX], c2phi[KMAX], c3phi[KMAX], c4phi[KMAX];
	  
  const double nu   = EOBPars->nu;
  const double chi1 = EOBPars->chi1;
  const double chi2 = EOBPars->chi2;
  const double X1   = EOBPars->X1;
  const double X2   = EOBPars->X2;
  const double aK   = EOBPars->a1+EOBPars->a2;
  const double Mbh  = EOBPars->Mbhf;
  const double abh  = EOBPars->abhf;
    
  double *t       = hlm_mrg->time;
  double *r       = dyn_mrg->data[EOB_RAD];
  double *w       = dyn_mrg->data[EOB_MOMG]; /* Omega */
  double *pph     = dyn_mrg->data[EOB_PPHI];
  double *pr_star = dyn_mrg->data[EOB_PRSTAR];
  double *Omg_orb = dyn_mrg->data[EOB_OMGORB]; /* Omega orbital */
  double *ddotr   = dyn_mrg->data[EOB_DDOTR];
  
  double P[2], M[4];   
  double max_A[KMAX],max_dA[KMAX],d2max[KMAX],d3max[KMAX],max_omg[KMAX],max_domg[KMAX],maxd2omg[KMAX], DeltaT[KMAX];
  double ai[KMAX][2] = {0.};
  double bi[KMAX][2] = {0.};
  
  const int size = hlm_mrg->size;
  for (int i = 0; i < size; i++) {
    hnqc->time[i] = t[i];
  }
  
  double *omg[KMAX], *domg[KMAX];
  double *n1[KMAX],*n2[KMAX],*n3[KMAX],*n4[KMAX],*n5[KMAX],*n6[KMAX],
    *d_n4[KMAX],*d_n5[KMAX],*d_n6[KMAX],*d2_n4[KMAX],*d2_n5[KMAX],*d2_n6[KMAX];
  double *m11[KMAX], *m12[KMAX], *m13[KMAX], *m21[KMAX], *m22[KMAX];
  double *p1tmp[KMAX], *p2tmp[KMAX]; /* RWZ amplitude and derivative */
  
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask[k]){
      omg[k]  = (double*) calloc (size,sizeof(double));
      domg[k] = (double*) calloc (size,sizeof(double));
      m11[k] = (double*) calloc (size,sizeof(double));
      m12[k] = (double*) calloc (size,sizeof(double));
      m13[k] = (double*) calloc (size,sizeof(double));
      m21[k] = (double*) calloc (size,sizeof(double));
      m22[k] = (double*) calloc (size,sizeof(double));
      p1tmp[k] = (double*) calloc (size,sizeof(double));
      p2tmp[k] = (double*) calloc (size,sizeof(double));
      n1[k] = (double*) calloc (size,sizeof(double));
      n2[k] = (double*) calloc (size,sizeof(double));
      n4[k] = (double*) calloc (size,sizeof(double));
      n5[k] = (double*) calloc (size,sizeof(double));
      d_n4[k] = (double*) calloc (size,sizeof(double));
      d_n5[k] = (double*) calloc (size,sizeof(double));
      d2_n4[k] = (double*) calloc (size,sizeof(double));
      d2_n5[k] = (double*) calloc (size,sizeof(double));
    }
  }

  /** omega derivatives */
  const double dt = t[1]-t[0];
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){
      D0(hlm_mrg->phase[k], dt, size, omg[k]);
      D0(omg[k], dt, size, domg[k]);
    }
  }
  
  /** NR fits */
  for (int k=0; k<KMAX; k++) {   
    max_A[k]    = 0.;
    max_dA[k]   = 0.;
    max_omg[k]  = 0.;
    max_domg[k] = 0.;
  }
  
  /* 31, 33, 41 and 55 fitted directly + 44 dA */
  eob_nqc_point_HM(dyn, max_A, max_dA, max_omg, max_domg);
  
  /* 21 fitted directly at tpeak_22*/
  eob_nqc_point_HM_peak22(dyn, max_A, max_dA, max_omg, max_domg);

  /* Over-writing fits using postpeak quantities for modes in kpostpeak */
  int kpostpeak_size = EOBPars->kpostpeak_size;  
  int *kpostpeak     = EOBPars->kpostpeak;

  QNMHybridFitCab_HM(nu, X1, X2, chi1, chi2, aK,  Mbh, abh,  
		     c1A, c2A, c3A, c4A, c1phi, c2phi, c3phi, c4phi,
		     alpha1, omega1);

  for (int j=0; j<kpostpeak_size; j++) {
    int k = kpostpeak[j];
    
    /* Normalizing c1A and c4A */
    int l = LINDEX[k];
    c1A[k] /= sqrt((l+2)*(l+1)*l*(l-1));
    c4A[k] /= sqrt((l+2)*(l+1)*l*(l-1));
    
    eob_nqc_point_postpeak(Mbh,c1A[k],c2A[k],c3A[k],c4A[k],
			   c1phi[k],c2phi[k],c3phi[k],c4phi[k],alpha1[k],omega1[k],
			   &A_tmp,&dA_tmp,&omg_tmp,&domg_tmp);
    
    max_A[k]    = A_tmp;
    max_dA[k]   = dA_tmp;
    max_omg[k]  = omg_tmp;
    max_domg[k] = domg_tmp;
  }
  
  if (VERBOSE) {
    printf("NR values for NQC determination:\n");
    PRFORMd("A22_mrg",max_A[1]);
    PRFORMd("dA22_mrg",max_dA[1]);
    PRFORMd("omg22_mrg",max_omg[1]);
    PRFORMd("domg22_mrg",max_domg[1]);
  }

  /** NQC corrections to AMPLITUDE (n1,n2,n3) and PHASE (n4,n5,n6)
   * NQC basis for (2,2) waveform : AMPLITUDE
   * note: n3 and n6 are not used
   */
  double pr_star2, r2, w2;
  for (int k=0; k<KMAX; k++) {   
    if(hlm_mrg->kmask_nqc[k]){ 
      for (int j=0; j<size; j++) {
	pr_star2 = SQ(pr_star[j]);
	r2       = SQ(r[j]);
	w2       = SQ(w[j]); //CHECKME: Omg or Omg_orbital ?
	n1[k][j]  = pr_star2/(r2*w2);         /* [pr*\/(r Omg)]^2 */
	n2[k][j]  = ddotr[j]/(r[j]*w2);       /* [ddot{r}/(r Omg^2)] */
	//n3[k][j]  = n1[k][j]*pr_star2;
	n4[k][j]  = pr_star[j]/(r[j]*w[j]);   /* pr*\/(r Omg) */
	n5[k][j]  = n4[k][j]*r2*w2;              /* (pr*)*(r Omg) */
	//n6[k][j]  = n5[k][j]*pr_star2;
      }
    }
  }

  for (int j=0; j<size; j++) {
    /* l=2,m=1 */
    if(hlm_mrg->kmask_nqc[0]){
      n2[0][j] = cbrt(SQ(w[j]))*n1[0][j];
      n5[0][j] = cbrt(SQ(w[j]))*n4[0][j];
    }
    /* l=2,m=2 */
    // FIXME: add a function specific for this RR!
    if(hlm_mrg->kmask_nqc[1] && EOBPars->use_flm == USEFLM_HM_4PN22) {
      n2[1][j] = SQ(pr_star[j])*n1[1][j];
    }
    /* l=3, l=4 & l=5 */
    for (int k=2; k<14; k++) {
      if(hlm_mrg->kmask_nqc[k]){
	n5[k][j]  = cbrt(SQ(w[j]))*n4[k][j];
      }
    }
  }
      
#if (DEBUG)
  FILE* fp_dbg;
  fp_dbg = fopen("nqc_nfunc.txt", "w");
  for (int j=0; j<size; j++) {
    fprintf(fp_dbg, "%20.12f\t%.16e\t%.16e\t%.16e\t%.16e\n", t[j], n1[1][j], n2[1][j], n4[1][j], n5[1][j]);
  }
  fclose(fp_dbg);
#endif    
      
  /** Derivatives for the phase */
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){   
      D0(n4[k],dt,size, d_n4[k]);
      D0(n5[k],dt,size, d_n5[k]);
      D0(d_n4[k],dt,size, d2_n4[k]);
      D0(d_n5[k],dt,size, d2_n5[k]);
    }
  }
    
#if (DEBUG)
  fp_dbg = fopen("nqc_dfunc.txt", "w");
  for (int j=0; j<size; j++) {
    fprintf(fp_dbg, "%f\t%.16e\t%.16e\t%.16e\t%.16e\n", t[j], d_n4[1][j], d_n5[1][j], d2_n4[1][j], d2_n5[1][j]);  
  }
  fclose(fp_dbg);
#endif    

  /** Find max Omg */
  int Omgmax_index = 0;
  double Omg_max   = Omg_orb[0];
  for (int j=0; j<size; j++) {
    if (Omg_orb[j] > Omg_max) {
      Omg_max = Omg_orb[j];
      Omgmax_index = j;
    }
  }
  
  /** Time */
  double tOmgOrb_pk = t[Omgmax_index];
  double DeltaT_nqc = eob_nqc_timeshift(nu, chi1);
  double tNQC = tOmgOrb_pk - DeltaT_nqc;

  if (VERBOSE) {
    printf("NQC info:\n");
    PRFORMd("DeltaT_tNQC",DeltaT_nqc);
    PRFORMd("tNQC[bare]",tNQC);
  }

  /** Find jmax: t[jmax] <= tNQC */
  double tmrg[KMAX];
  tmrg[1] = tNQC - 2.;
	
  int jmax = 0;
  for (int j=0; j<size; j++) {
    if(t[j] > tNQC) {
      jmax = j-2;
      break;
    }
  }
  
  double dtmrg[KMAX];
  double t_NQC[KMAX];
  int    j_NQC[KMAX];

  /* Usually, the NQC match point is at tmrg_lm + 2
   * except for modes in knqcpeak22 where it is at tmrg_22
  */
  int modesatpeak22[KMAX]; 
  set_multipolar_idx_mask (modesatpeak22, KMAX, EOBPars->knqcpeak22, EOBPars->knqcpeak22_size, 0);
  eob_nqc_deltat_lm(dyn, dtmrg);

  for (int k=0; k<KMAX; k++) {   
    if(hlm_mrg->kmask_nqc[k]){
      tmrg[k]  = tmrg[1] + dtmrg[k];
      t_NQC[k] = tmrg[k] + 2.;

      if(modesatpeak22[k]) t_NQC[k] = tmrg[1];

      j_NQC[k] = size-1;
      for (int j=size-2; j>=0; j--) {
	      if(t[j] < t_NQC[k]) {
	        break;
	      }
	      j_NQC[k] = j;
      }
    }
  }

  /** Solve the linear systems */
  
  /* Regge-Wheeler-Zerilli normalized amplitude. 
     The ringdown coefficient refer to this normalization.
     Nagar & Rezzolla, CQG 22 (2005) R167 */      
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){
      double nlm = 1./(sqrt( (LINDEX[k]+2)*(LINDEX[k]+1)*LINDEX[k]*(LINDEX[k]-1) ) );

      if (hlm->ampli[k][0] > 0.) {	
        nNegAmp[k] = 0;
      } else {	
        nNegAmp[k] = 1;	
      }
      
      for (int j=0; j<size; j++) {
        p1tmp[k][j] = fabs(hlm_mrg->ampli[k][j] * nlm);      
      }
    }
  }
  
  /* Matrix elements: waveform amplitude at all points */
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){
      for (int j=0; j<size; j++) {
        m11[k][j] = n1[k][j] * p1tmp[k][j];
        m12[k][j] = n2[k][j] * p1tmp[k][j];
      }
    }
  }

  /* Take FD derivatives */
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){
      D0(m11[k],dt,size, m21[k]);
      D0(m12[k],dt,size, m22[k]);
      D0(p1tmp[k],dt,size, p2tmp[k]);
    }
  }

#if (DEBUG)
  fp_dbg = fopen("nqc_amp_func.txt", "w");
  for (int j=0; j<size; j++) {
    fprintf(fp_dbg, "%e\t%e\t%e\n", t[j], p1tmp[1][j], p2tmp[1][j]);
  }
  fclose(fp_dbg);  
#endif

  double detM = 1.;
  double oodetM = 1.;
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){

      jmax = j_NQC[k];

      /* Computation of ai coefficients at Omega peak */
      P[0]     = max_A[k]  - p1tmp[k][jmax];
      P[1]     = max_dA[k] - p2tmp[k][jmax];
      
      M[0]     = m11[k][jmax];
      M[1]     = m12[k][jmax];
      M[2]     = m21[k][jmax];
      M[3]     = m22[k][jmax];
      
      /* detM     = M[0]*M[3]-M[1]*M[2]; 
        ai[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
        ai[k][1] = (M[0]*P[1] - M[2]*P[0])/detM; */
      /* safe version (amplitude can be zero) */
      oodetM   = 1.0/(M[0]*M[3]-M[1]*M[2]);
      if (isfinite(oodetM)) {
        ai[k][0] = (M[3]*P[0] - M[1]*P[1])*oodetM;
        ai[k][1] = (M[0]*P[1] - M[2]*P[0])*oodetM;
      }

      /* Computation of bi coefficients at Omega peak */
      P[0]     = omg[k][jmax]   - max_omg[k];
      P[1]     = domg[k][jmax]  - max_domg[k];
      
      M[0]     = d_n4[k][jmax];
      M[1]     = d_n5[k][jmax];
      M[2]     = d2_n4[k][jmax];
      M[3]     = d2_n5[k][jmax];

      /* detM     =  M[0]*M[3] - M[1]*M[2];
        bi[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
        bi[k][1] = (M[0]*P[1] - M[2]*P[0])/detM; */
      /* safe version (phase can be zero) */
      oodetM   = 1.0/(M[0]*M[3]-M[1]*M[2]);
      if (isfinite(oodetM)) {
        bi[k][0] = (M[3]*P[0] - M[1]*P[1])*oodetM;
        bi[k][1] = (M[0]*P[1] - M[2]*P[0])*oodetM;
      }
    }

  }
  
  if (VERBOSE){
    printf("NQC coefficients for 22 mode:\n");
    PRFORMd("a1",ai[1][0]);
    PRFORMd("a2",ai[1][1]);
    PRFORMd("b1",bi[1][0]);
    PRFORMd("b2",bi[1][1]);
  }

  /** Set amplitude and phase */
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){
      for (int j=0; j<size; j++) {
        hnqc->ampli[k][j] = 1. + ai[k][0]*n1[k][j] + ai[k][1]*n2[k][j];
        hnqc->phase[k][j] =      bi[k][0]*n4[k][j] + bi[k][1]*n5[k][j];
      }
    }
  }
  
  /** Multiply merger waveform to NQC */
  for (int k=0; k<KMAX; k++) {
    for (int j=0; j<size; j++) {
      if(hlm_mrg->kmask_nqc[k]){
        hlm_mrg->ampli[k][j] *= hnqc->ampli[k][j];
        hlm_mrg->phase[k][j] -= hnqc->phase[k][j];
      }
    }
  }

  /** Multiply full waveform to NQC */
  r       = dyn->data[EOB_RAD];
  w       = dyn->data[EOB_MOMG]; /* Omega */
  pph     = dyn->data[EOB_PPHI];
  pr_star = dyn->data[EOB_PRSTAR];
  Omg_orb = dyn->data[EOB_OMGORB]; /* Omega orbital */
  ddotr   = dyn->data[EOB_DDOTR];

  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask[k]){
      free(n1[k]);
      free(n2[k]);
      free(n4[k]);
      free(n5[k]);
    }
  }

  const int fullsize = hlm->size;
  
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask[k]){
      n1[k] = (double*) calloc (fullsize,sizeof(double));
      n2[k] = (double*) calloc (fullsize,sizeof(double));
      n4[k] = (double*) calloc (fullsize,sizeof(double));
      n5[k] = (double*) calloc (fullsize,sizeof(double));
    }
  }
	
  for (int k=0; k<KMAX; k++) {   
    if(hlm_mrg->kmask_nqc[k]){   
      for (int j=0; j<fullsize; j++) {
        pr_star2 = SQ(pr_star[j]);
        r2       = SQ(r[j]);
        w2       = SQ(w[j]); //CHECKME: Omg or Omg_orbital ?
        n1[k][j]  = pr_star2/(r2*w2);         /* [pr*\/(r Omg)]^2 */
        n2[k][j]  = ddotr[j]/(r[j]*w2);       /* [ddot{r}/(r Omg^2)] */
        n4[k][j]  = pr_star[j]/(r[j]*w[j]);   /* pr*\/(r Omg) */
        n5[k][j]  = n4[k][j]*r2*w2;           /* (pr*)*(r Omg) */
      }
    }
  }

  for (int j=0; j<fullsize; j++) {
    /* l=2,m=1 */
    if(hlm_mrg->kmask_nqc[0]){
      n2[0][j] = cbrt(SQ(w[j]))*n1[0][j];
      n5[0][j] = cbrt(SQ(w[j]))*n4[0][j];
    }
    /* l=3, l=4 & l=5 */
    for (int k=2; k<14; k++) {   
      if(hlm_mrg->kmask_nqc[k]){   
	n5[k][j]  = cbrt(SQ(w[j]))*n4[k][j];
      }
    }
  }
  
  for (int k=0; k<KMAX; k++) {
    if(hlm->kmask_nqc[k]){
      for (int j=0; j<fullsize; j++) {
        hlm->ampli[k][j] *= (1. + ai[k][0]*n1[k][j] + ai[k][1]*n2[k][j]);
        hlm->phase[k][j] -= (bi[k][0]*n4[k][j] + bi[k][1]*n5[k][j]);
      }
    }
  }

  if (EOBPars->output_nqc_coefs) {
    /** Output the NQC coefficients */
    FILE* fp;
    char fname[STRLEN];
    strcpy(fname, EOBPars->output_dir);
    strcat(fname, "/nqc_coefs.txt");
    fp = fopen(fname, "w");
    fprintf(fp, "# q=%e chizA=%e chizB=%e f0=%e\n",EOBPars->q,EOBPars->chi1,EOBPars->chi2,EOBPars->initial_frequency);
    fprintf(fp, "# M=%e LambdaA=[%e,%e,%e] LambdaBl2=[%e,%e,%e]\n",EOBPars->M,
	    EOBPars->LambdaAl2,EOBPars->LambdaAl3,EOBPars->LambdaAl4,
	    EOBPars->LambdaBl2,EOBPars->LambdaBl3,EOBPars->LambdaBl4);
    for (int k=0; k<KMAX; k++) {
      if(hlm->kmask_nqc[k]){
        fprintf(fp, "%d %d %d %e %e %e %e\n", k, LINDEX[k], MINDEX[k], 
        ai[k][0], ai[k][1], 
        bi[k][0], bi[k][1]);
      }
    }  
    fclose(fp);  
  }
  
  /** Free mem */
  for (int k=0; k<KMAX; k++) {
    if(hlm->kmask[k]){
      free(omg[k]);
      free(domg[k]);
      free(m11[k]);
      free(m12[k]);
      free(m13[k]);
      free(m21[k]);
      free(m22[k]);
      free(p1tmp[k]);
      free(p2tmp[k]);
      free(n1[k]);
      free(n2[k]);
      free(n4[k]);
      free(n5[k]);
      free(d_n4[k]);
      free(d_n5[k]);
      free(d2_n4[k]);
      free(d2_n5[k]);
    }
  }

}

/**
 * Function: eob_wav_hlmNQC_find_a1a2a3_mrg_22
 * -------------------------------------------
 *   Computes the factors and the coefficients that build the
 *   NQC corrections to the waveform in the spinning case.
 *   This routine works around merger with dyn_mrg and hlm_mrg and
 *   then add everything also to hlm
 *   Specific for the 22 mode
 * 
 *   @param[in] dyn_mrg : dynamics around merger
 *   @param[in] hlm_mrg : multipolar waveform around merger (modified by this function)
 *   @param[in] hnqc    : NQC part of the waveform
 *   @param[in] dyn     : full dynamics
 *   @param[in] hlm     : full multipolar waveform (modified by this function)
*/
void eob_wav_hlmNQC_find_a1a2a3_mrg_22(Dynamics *dyn_mrg, Waveform_lm *hlm_mrg, Waveform_lm *hnqc,
				       Dynamics *dyn, Waveform_lm *hlm)
{
#define KMAX_22 (2) // 0, 1 -> (2,1), (2,2) // k=0 unused.

  double A_tmp, dA_tmp, omg_tmp, domg_tmp;
  double alpha1[KMAX_22], omega1[KMAX_22];
  double c1A[KMAX_22], c2A[KMAX_22], c3A[KMAX_22], c4A[KMAX_22];
  double c1phi[KMAX_22], c2phi[KMAX_22], c3phi[KMAX_22], c4phi[KMAX_22];
	  
  const double nu   = EOBPars->nu;
  const double chi1 = EOBPars->chi1;
  const double chi2 = EOBPars->chi2;
  const double X1   = EOBPars->X1;
  const double X2   = EOBPars->X2;
  const double aK   = EOBPars->a1+EOBPars->a2;
  const double Mbh  = EOBPars->Mbhf;
  const double abh  = EOBPars->abhf;
    
  double *t       = hlm_mrg->time;
  double *r       = dyn_mrg->data[EOB_RAD];
  double *w       = dyn_mrg->data[EOB_MOMG]; /* Omega */
  double *pph     = dyn_mrg->data[EOB_PPHI];
  double *pr_star = dyn_mrg->data[EOB_PRSTAR];
  double *Omg_orb = dyn_mrg->data[EOB_OMGORB]; /* Omega orbital */
  double *ddotr   = dyn_mrg->data[EOB_DDOTR];
  
  double P[2], M[4];   
  double max_A[KMAX_22],max_dA[KMAX_22],d2max[KMAX_22],d3max[KMAX_22],max_omg[KMAX_22],max_domg[KMAX_22],maxd2omg[KMAX_22], DeltaT[KMAX_22];
  double ai[KMAX_22][2];
  double bi[KMAX_22][2];
  
  for (int k=0; k<KMAX_22; k++) {  
    max_A[k]    = 0.;
    max_dA[k]   = 0.;
    max_omg[k]  = 0.;
    max_domg[k] = 0.;
  }

  const int size = hlm_mrg->size;
  for (int i = 0; i < size; i++) {
    hnqc->time[i] = t[i];
  }
  
  double *omg[KMAX_22], *domg[KMAX_22];
  double *n1[KMAX_22],*n2[KMAX_22],*n3[KMAX_22],*n4[KMAX_22],*n5[KMAX_22],*n6[KMAX_22],
    *d_n4[KMAX_22],*d_n5[KMAX_22],*d_n6[KMAX_22],*d2_n4[KMAX_22],*d2_n5[KMAX_22],*d2_n6[KMAX_22];
  double *m11[KMAX_22], *m12[KMAX_22], *m13[KMAX_22], *m21[KMAX_22], *m22[KMAX_22];
  double *p1tmp[KMAX_22], *p2tmp[KMAX_22]; /* RWZ amplitude and derivative */
  
  const int k22 = 1; // (2,2) mode
  
  omg[k22]  = (double*) calloc (size,sizeof(double));
  domg[k22] = (double*) calloc (size,sizeof(double));
  m11[k22] = (double*) calloc (size,sizeof(double));
  m12[k22] = (double*) calloc (size,sizeof(double));
  m13[k22] = (double*) calloc (size,sizeof(double));
  m21[k22] = (double*) calloc (size,sizeof(double));
  m22[k22] = (double*) calloc (size,sizeof(double));
  p1tmp[k22] = (double*) calloc (size,sizeof(double));
  p2tmp[k22] = (double*) calloc (size,sizeof(double));
  n1[k22] = (double*) calloc (size,sizeof(double));
  n2[k22] = (double*) calloc (size,sizeof(double));
  n4[k22] = (double*) calloc (size,sizeof(double));
  n5[k22] = (double*) calloc (size,sizeof(double));
  d_n4[k22] = (double*) calloc (size,sizeof(double));
  d_n5[k22] = (double*) calloc (size,sizeof(double));
  d2_n4[k22] = (double*) calloc (size,sizeof(double));
  d2_n5[k22] = (double*) calloc (size,sizeof(double));

  /** omega derivatives */
  const double dt = t[1]-t[0];
  D0(hlm_mrg->phase[k22], dt, size, omg[k22]);
  D0(omg[k22], dt, size, domg[k22]);
  
  /** NR fits */
  eob_nqc_point(dyn, &A_tmp, &dA_tmp, &omg_tmp, &domg_tmp);
	
  /** Switch on the 22 values (only) */
  max_A[k22]    = A_tmp;
  max_dA[k22]   = dA_tmp;
  max_omg[k22]  = omg_tmp;
  max_domg[k22] = domg_tmp;

  if (VERBOSE) {
    printf("NR values for NQC determination:\n");
    PRFORMd("A22_mrg",max_A[1]);
    PRFORMd("dA22_mrg",max_dA[1]);
    PRFORMd("omg22_mrg",max_omg[1]);
    PRFORMd("domg22_mrg",max_domg[1]);
  }
  
  /** NQC corrections to AMPLITUDE (n1,n2,n3) and PHASE (n4,n5,n6)
   * NQC basis for (2,2) waveform : AMPLITUDE
   * note: n3 and n6 are not used
   */
  double pr_star2, r2, w2;
  for (int j=0; j<size; j++) {
    pr_star2 = SQ(pr_star[j]);
    r2       = SQ(r[j]);
    w2       = SQ(w[j]); 
    n1[k22][j]  = pr_star2/(r2*w2);         /* [pr*\/(r Omg)]^2 */
    n2[k22][j]  = ddotr[j]/(r[j]*w2);       /* [ddot{r}/(r Omg^2)] */
    //n3[k22][j]  = n1[k22][j]*pr_star2;
    n4[k22][j]  = pr_star[j]/(r[j]*w[j]);   /* pr*\/(r Omg) */
    n5[k22][j]  = n4[k22][j]*r2*w2;         /* (pr*)*(r Omg) */
    //n6[k22][j]  = n5[k22][j]*pr_star2;
  }
  
#if (DEBUG)
  FILE* fp_dbg;
  fp_dbg = fopen("nqc_nfunc.txt", "w");
  for (int j=0; j<size; j++) {
    fprintf(fp_dbg, "%20.12f\t%.16e\t%.16e\t%.16e\t%.16e\n", t[j], n1[k22][j], n2[k22][j], n4[k22][j], n5[k22][j]);
  }
  fclose(fp_dbg);
#endif    

  /** Derivatives for the phase */
  D0(n4[k22],dt,size, d_n4[k22]);
  D0(n5[k22],dt,size, d_n5[k22]);
  D0(d_n4[k22],dt,size, d2_n4[k22]);
  D0(d_n5[k22],dt,size, d2_n5[k22]);
  
#if (DEBUG)
  fp_dbg = fopen("nqc_dfunc.txt", "w");
  for (int j=0; j<size; j++) {
    fprintf(fp_dbg, "%f\t%.16e\t%.16e\t%.16e\t%.16e\n", t[j], d_n4[k22][j], d_n5[k22][j], d2_n4[k22][j], d2_n5[k22][j]);  
  }
  fclose(fp_dbg);
#endif    

  /** Find max Omg */
  int Omgmax_index = 0;
  double Omg_max   = Omg_orb[0];
  for (int j=0; j<size; j++) {
    if (Omg_orb[j] > Omg_max) {
      Omg_max = Omg_orb[j];
      Omgmax_index = j;
    }
  }
  
  /** Time */
  double tOmgOrb_pk = t[Omgmax_index];
  double DeltaT_nqc = eob_nqc_timeshift(nu, chi1);
  double tNQC = tOmgOrb_pk - DeltaT_nqc;

  if (VERBOSE) {
    printf("NQC info:\n");
    PRFORMd("DeltaT_tNQC",DeltaT_nqc);
    PRFORMd("tNQC[bare]",tNQC);
  }

  /** Find jmax: t[jmax] <= tNQC */
  double tmrg[KMAX_22];
  tmrg[k22] = tNQC - 2.;
	
  int jmax = 0;
  for (int j=0; j<size; j++) {
    if(t[j] > tNQC) {
      jmax = j-2;
      break;
    }
  }
  
  double dtmrg[KMAX_22];
  double t_NQC[KMAX_22];
  int    j_NQC[KMAX_22];

  /** Solve the linear systems */
  
  /* Regge-Wheeler-Zerilli normalized amplitude. 
     The ringdown coefficient refer to this normalization.
     Nagar & Rezzolla, CQG 22 (2005) R167 */      
  const double nlm = 1./(sqrt( (LINDEX[k22]+2)*(LINDEX[k22]+1)*LINDEX[k22]*(LINDEX[k22]-1) ) );

  if (hlm->ampli[k22][0] > 0.) {	
    nNegAmp[k22] = 0;
  } else {		
    nNegAmp[k22] = 1;
  }
      
  for (int j=0; j<size; j++) {
    p1tmp[k22][j] = fabs(hlm_mrg->ampli[k22][j] * nlm);      
  }

  /* Matrix elements: waveform amplitude at all points */
  for (int j=0; j<size; j++) {
    m11[k22][j] = n1[k22][j] * p1tmp[k22][j];
    m12[k22][j] = n2[k22][j] * p1tmp[k22][j];
  }

  /* Take FD derivatives */
  D0(m11[k22],dt,size, m21[k22]);
  D0(m12[k22],dt,size, m22[k22]);
  D0(p1tmp[k22],dt,size, p2tmp[k22]);

#if (DEBUG)
  fp_dbg = fopen("nqc_amp_func.txt", "w");
  for (int j=0; j<size; j++) {
    fprintf(fp_dbg, "%e\t%e\t%e\n", t[j], p1tmp[k22][j], p2tmp[k22][j]);
  }
  fclose(fp_dbg);  
#endif

  double detM = 1.;
  double oodetM = 1.;
  ai[k22][0] = ai[k22][1] = 0.;
  bi[k22][0] = bi[k22][1] = 0.;

  /* Computation of ai coefficients at Omega peak */
  P[0]     = max_A[k22]  - p1tmp[k22][jmax];
  P[1]     = max_dA[k22] - p2tmp[k22][jmax];
  
  M[0]     = m11[k22][jmax];
  M[1]     = m12[k22][jmax];
  M[2]     = m21[k22][jmax];
  M[3]     = m22[k22][jmax];
      
  /* detM     = M[0]*M[3]-M[1]*M[2]; 
     ai[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
     ai[k][1] = (M[0]*P[1] - M[2]*P[0])/detM; */
  /* safe version (amplitude can be zero) */
  oodetM   = 1.0/(M[0]*M[3]-M[1]*M[2]);
  if (isfinite(oodetM)) {
    ai[k22][0] = (M[3]*P[0] - M[1]*P[1])*oodetM;
    ai[k22][1] = (M[0]*P[1] - M[2]*P[0])*oodetM;
  }

  /* Computation of bi coefficients at Omega peak */
  P[0]     = omg[k22][jmax]   - max_omg[k22];
  P[1]     = domg[k22][jmax]  - max_domg[k22];
  
  M[0]     = d_n4[k22][jmax];
  M[1]     = d_n5[k22][jmax];
  M[2]     = d2_n4[k22][jmax];
  M[3]     = d2_n5[k22][jmax];
  
  /* detM     =  M[0]*M[3] - M[1]*M[2];
     bi[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
     bi[k][1] = (M[0]*P[1] - M[2]*P[0])/detM; */
  /* safe version (phase can be zero) */
  oodetM   = 1.0/(M[0]*M[3]-M[1]*M[2]);
  if (isfinite(oodetM)) {
    bi[k22][0] = (M[3]*P[0] - M[1]*P[1])*oodetM;
    bi[k22][1] = (M[0]*P[1] - M[2]*P[0])*oodetM;
  }
  
  if (VERBOSE){
    printf("NQC coefficients for 22 mode:\n");
    PRFORMd("a1",ai[k22][0]);
    PRFORMd("a2",ai[k22][1]);
    PRFORMd("b1",bi[k22][0]);
    PRFORMd("b2",bi[k22][1]);
  }

  /** Set amplitude and phase */
  for (int j=0; j<size; j++) {
    hnqc->ampli[k22][j] = 1. + ai[k22][0]*n1[k22][j] + ai[k22][1]*n2[k22][j];
    hnqc->phase[k22][j] =      bi[k22][0]*n4[k22][j] + bi[k22][1]*n5[k22][j];
  }
  
  /** Multiply merger waveform to NQC */
  for (int j=0; j<size; j++) {
    hlm_mrg->ampli[k22][j] *= hnqc->ampli[k22][j];
    hlm_mrg->phase[k22][j] -= hnqc->phase[k22][j];
  }

  /** Multiply full waveform to NQC */
  r       = dyn->data[EOB_RAD];
  w       = dyn->data[EOB_MOMG]; /* Omega */
  pph     = dyn->data[EOB_PPHI];
  pr_star = dyn->data[EOB_PRSTAR];
  Omg_orb = dyn->data[EOB_OMGORB]; /* Omega orbital */
  ddotr   = dyn->data[EOB_DDOTR];

  free(n1[k22]);
  free(n2[k22]);
  free(n4[k22]);
  free(n5[k22]);

  const int fullsize = hlm->size;
  
  n1[k22] = (double*) calloc (fullsize,sizeof(double));
  n2[k22] = (double*) calloc (fullsize,sizeof(double));
  n4[k22] = (double*) calloc (fullsize,sizeof(double));
  n5[k22] = (double*) calloc (fullsize,sizeof(double));
	
  for (int j=0; j<fullsize; j++) {
    pr_star2 = SQ(pr_star[j]);
    r2       = SQ(r[j]);
    w2       = SQ(w[j]); 
    n1[k22][j]  = pr_star2/(r2*w2);         /* [pr*\/(r Omg)]^2 */
    n2[k22][j]  = ddotr[j]/(r[j]*w2);       /* [ddot{r}/(r Omg^2)] */
    n4[k22][j]  = pr_star[j]/(r[j]*w[j]);   /* pr*\/(r Omg) */
    n5[k22][j]  = n4[k22][j]*r2*w2;              /* (pr*)*(r Omg) */
  }
	  	  
  for (int j=0; j<fullsize; j++) {
    hlm->ampli[k22][j] *= (1. + ai[k22][0]*n1[k22][j] + ai[k22][1]*n2[k22][j]);
    hlm->phase[k22][j] -= (bi[k22][0]*n4[k22][j] + bi[k22][1]*n5[k22][j]);
  }

  if (EOBPars->output_nqc_coefs) {
    /** Output the NQC coefficients */
    FILE* fp;
    char fname[STRLEN];
    strcpy(fname, EOBPars->output_dir);
    strcat(fname, "/nqc_coefs.txt");
    fp = fopen(fname, "w");
    fprintf(fp, "# q=%e chizA=%e chizB=%e f0=%e\n",EOBPars->q,EOBPars->chi1,EOBPars->chi2,EOBPars->initial_frequency);
    fprintf(fp, "# M=%e LambdaA=[%e,%e,%e] LambdaBl2=[%e,%e,%e]\n",EOBPars->M,
	    EOBPars->LambdaAl2,EOBPars->LambdaAl3,EOBPars->LambdaAl4,
	    EOBPars->LambdaBl2,EOBPars->LambdaBl3,EOBPars->LambdaBl4);
    fprintf(fp, "%d %d %d %e %e %e %e\n", k22, LINDEX[k22], MINDEX[k22], 
	    ai[k22][0], ai[k22][1], 
	    bi[k22][0], bi[k22][1]);
    fclose(fp);  
  }
  
  /** Free mem */
  free(omg[k22]);
  free(domg[k22]);
  free(m11[k22]);
  free(m12[k22]);
  free(m13[k22]);
  free(m21[k22]);
  free(m22[k22]);
  free(p1tmp[k22]);
  free(p2tmp[k22]);
  free(n1[k22]);
  free(n2[k22]);
  free(n4[k22]);
  free(n5[k22]);
  free(d_n4[k22]);
  free(d_n5[k22]);
  free(d2_n4[k22]);
  free(d2_n5[k22]);
}

/** 
 * Function: eob_wav_hlmNQC_nospin201602
 *   NQC corrections to the RWZ multipolar waveform
 *   Nagar, Damour, Reisswig, Pollney http://arxiv.org/abs/1506.08457
 *   Nonspinning case, Current fits: 9/02/2016
 * 
 *   @param[in] nu      : symmetric mass ratio
 *   @param[in] r       : radial separation
 *   @param[in] prstar  : (tortoise) radial momentum
 *   @param[in] Omega   : orbital frequency
 *   @param[in] ddotr   : radial acceleration
 *   @param[out] hlmnqc : NQC part of the waveform
*/
void eob_wav_hlmNQC_nospin201602(double  nu, double  r, double  prstar, double  Omega, double  ddotr,
				 Waveform_lm_t *hlmnqc)
{       
  const double xnu  = 1-4*nu;
  const double xnu2 = SQ(xnu);

  double a1[KMAX], a2[KMAX], a3[KMAX];
  double b1[KMAX], b2[KMAX], b3[KMAX];
  double n[KMAX][6];    

  const int k21 = 0;
  const int k22 = 1;
  const int k33 = 4;
  int k;

  /** NR fits */  
  for (int k = 0; k < KMAX; k++) {
    a1[k] = 0.;
    a2[k] = 0.;
    a3[k] = 0.;
    b1[k] = 0.;
    b2[k] = 0.;
    b3[k] = 0.;
  }
  
  /* (2,1) */
  a1[k21] =  0.0162387198*(7.32653082*xnu2 + 1.19616248*xnu + 0.73496656);
  a2[k21] =                -1.80492460*xnu2 + 1.78172686*xnu + 0.30865284;
  a3[k21] =                                                           0.0;
  
  b1[k21] =  -0.0647955017*(3.59934444*xnu2 - 4.08628784*xnu + 1.37890907);
  b2[k21] =   1.3410693180*(0.38491989*xnu2 + 0.10969453*xnu + 0.97513971);
  b3[k21] =                                                            0.0;
  
  /* (2,2) */
  a1[k22]   = -0.0805236959*( 1 - 2.00332326*xnu2)/( 1 + 3.08595088*xnu2);
  a2[k22]   =  1.5299534255*( 1 + 1.16438929*xnu2)/( 1 + 1.92033923*xnu2);
  a3[k22]   =  0.0;
  
  b1[k22]   = 0.146768094955*( 0.07417121*xnu + 1.01691256);
  b2[k22]   = 0.896911234248*(-0.61072011*xnu + 0.94295129);
  b3[k22]   = 0.0;
  
  /* (3,3) */
  a1[k33]   = -0.0377680000*(1 - 14.61548907*xnu2)/( 1 + 2.44559263*xnu2);
  a2[k33]   =  1.9898000000*(1 + 2.09750346 *xnu2)/( 1 + 2.57489466*xnu2);
  a3[k33]   =  0.0;
  
  b1[k33]   = 0.1418400000*(1.07430512 - 1.23906804*xnu + 4.44910652*xnu2);
  b2[k33]   = 0.6191300000*(0.80672432 + 4.07432829*xnu - 7.47270977*xnu2);
  b3[k33]   = 0.0;
  
  /** NQC corrections to the modulus and phase */
  for (int k = 0; k < KMAX; k++) {
    for (int j = 0; j < 6; j++) {
      n[k][j] = 0.;
    }    
  }

  k = k21;
  n[k][0] = (prstar/(r*Omega))*(prstar/(r*Omega));
  n[k][1] = ddotr/(r*Omega*Omega);
  n[k][2] = n[k][0]*prstar*prstar;
  n[k][3] = prstar/(r*Omega);
  n[k][4] = n[k][3]*cbrt(Omega*Omega);
  n[k][5] = n[k][4]*prstar*prstar;

  k = k22;
  n[k][0] = (prstar/(r*Omega))*(prstar/(r*Omega));
  n[k][1] = ddotr/(r*Omega*Omega);
  n[k][2] = n[k][0]*prstar*prstar;
  n[k][3] = prstar/(r*Omega);
  /* n[k][4] = n[k][3]*cbrt(Omega*Omega); */
  n[k][4] = n[k][3]*(r*Omega)*(r*Omega);
  n[k][5] = n[k][4]*prstar*prstar;

  k = k33;
  n[k][0] = (prstar/(r*Omega))*(prstar/(r*Omega));
  n[k][1] = ddotr/(r*Omega*Omega);
  n[k][2] = n[k][0]*prstar*prstar;
  n[k][3] = prstar/(r*Omega);
  n[k][4] = n[k][3]*cbrt(Omega*Omega);
  n[k][5] = n[k][4]*prstar*prstar;

  /** NQC factor */
  for (int k = 0; k < KMAX; k++) {
    hlmnqc->ampli[k] = 1.; 
    hlmnqc->phase[k] = 0.; 
  }

  k = k21; /* (2,1) */
  hlmnqc->ampli[k] = 1. + a1[k]*n[k][0] + a2[k]*n[k][1] + a3[k]*n[k][2];
  hlmnqc->phase[k] =      b1[k]*n[k][3] + b2[k]*n[k][4] + b3[k]*n[k][5];
  
  k = k22; /* (2,2) */
  hlmnqc->ampli[k] = 1. + a1[k]*n[k][0] + a2[k]*n[k][1] + a3[k]*n[k][2];
  hlmnqc->phase[k] =      b1[k]*n[k][3] + b2[k]*n[k][4] + b3[k]*n[k][5];
  
  k = k33; /* (3,3) */
  hlmnqc->ampli[k] = 1. + a1[k]*n[k][0] + a2[k]*n[k][1] + a3[k]*n[k][2];
  hlmnqc->phase[k] =      b1[k]*n[k][3] + b2[k]*n[k][4] + b3[k]*n[k][5];
  
}

 /**
 * Function: eob_wav_hlmNQC
 * ------------------------
 *   Generic routine for NQC
 * 
 *   @param[in] r       : radial separation
 *   @param[in] prstar  : (tortoise) radial momentum
 *   @param[in] Omega   : orbital frequency
 *   @param[in] ddotr   : radial acceleration
 *   @param[in] nqc     : NQC coefficients
 *   @param[out] hlmnqc : NQC part of the waveform
 */
void eob_wav_hlmNQC(double  nu, double  r, double  prstar, double  Omega, double  ddotr, NQCcoefs *nqc, 
		    Waveform_lm_t *hlmnqc)
{      
  const int maxk = MIN(KMAX, nqc->maxk+1);
  
  /* Multipoles with special treatment */
  const int k22 = 1;

  /* Shorthand */
  const double n0 = (prstar/(r*Omega))*(prstar/(r*Omega));
  const double n1 = ddotr/(r*Omega*Omega);
  const double n2 = n0*SQ(prstar);
  const double n3 = prstar/(r*Omega);
  const double n4 = n3*cbrt(Omega*Omega);
  const double n5 = n4*SQ(prstar);
  
  const double n4_k = n3*SQ((r*Omega)); 
  const double n5_k = n4*SQ(prstar);
  
  /** n functions */
  for (int k = 0; k < maxk; k++) {
    if (nqc->activemode[k]) {
      nqc->n[k][0] = n0;
      nqc->n[k][1] = n1;
      nqc->n[k][2] = n2;
      nqc->n[k][3] = n3;
      nqc->n[k][4] = n4;
      nqc->n[k][5] = n5;
    }
  }

  /** Change special multipoles */
  int k = k22;
  nqc->n[k][4] = n4_k;
  nqc->n[k][5] = n5_k;
  
  /** NQC wave factor */
  for (int k = 0; k < KMAX; k++) {
    hlmnqc->ampli[k] = 1.;
    hlmnqc->phase[k] = 0.;
  }
  
  for (int k = 0; k < maxk; k++) {
    if (nqc->activemode[k]) {
      hlmnqc->ampli[k] += nqc->a1[k]*nqc->n[k][0] + nqc->a2[k]*nqc->n[k][1] + nqc->a3[k]*nqc->n[k][2]; 
      hlmnqc->phase[k] += nqc->b1[k]*nqc->n[k][3] + nqc->b2[k]*nqc->n[k][4] + nqc->b3[k]*nqc->n[k][5]; 
    }
  }
  
}


 /**
 * Function: eob_wav_hlmNQC
 * ------------------------
 *   Generic routine for NQC with sigmoid (ecc systems)
 * 
 *   @param[in] r       : radial separation
 *   @param[in] prstar  : (tortoise) radial momentum
 *   @param[in] Omega   : orbital frequency
 *   @param[in] ddotr   : radial acceleration
 *   @param[in] nqc     : NQC coefficients
 *   @param[out] hlmnqc : NQC part of the waveform
 */
void eob_wav_hlmNQC_ecc_sigmoid(double  nu, double  r, double  prstar, double  Omega, double  ddotr, double t, double tOmg_pk, NQCcoefs *nqc, 
		    Waveform_lm_t *hlmnqc)
{      
  const int maxk = MIN(KMAX, nqc->maxk+1);
  
  /* Multipoles with special treatment */
  const int k22 = 1;

  /* Shorthand */
  const double n0 = (prstar/(r*Omega))*(prstar/(r*Omega));
  const double n1 = ddotr/(r*Omega*Omega);
  const double n2 = n0*SQ(prstar);
  const double n3 = prstar/(r*Omega);
  const double n4 = n3*cbrt(Omega*Omega);
  const double n5 = n4*SQ(prstar);
  
  const double n4_k = n3*SQ((r*Omega)); 
  const double n5_k = n4*SQ(prstar);

  /* time shift */
  double chi1 = EOBPars->chi1;
  double DeltaT_nqc = eob_nqc_timeshift(nu, chi1);
  double tNQC = tOmg_pk - DeltaT_nqc;

  /* Defining sigmoid function to switch on NQCs near the end of the evolution */
  /*
    Old configuration used for arXiv:2001.11736
    double t0 = tNQC - 30.;
    double alpha = 0.09;
  */
  double t0 = tNQC - EOBPars->delta_t0_sigmoid_NQC;
  double alpha = EOBPars->alpha_sigmoid_NQC;
  double sigmoid = 1./(1. + exp(-alpha*(t - t0)));

  /** n functions */
  for (int k = 0; k < maxk; k++) {
    if (nqc->activemode[k]) {
      nqc->n[k][0] = n0*sigmoid;
      nqc->n[k][1] = n1*sigmoid;
      nqc->n[k][2] = n2*sigmoid;
      nqc->n[k][3] = n3*sigmoid;
      nqc->n[k][4] = n4*sigmoid;
      nqc->n[k][5] = n5*sigmoid;
    }
  }

  /** Change special multipoles */
  int k = k22;
  nqc->n[k][4] = n4_k*sigmoid;
  nqc->n[k][5] = n5_k*sigmoid;
  
  /** NQC wave factor */
  for (int k = 0; k < KMAX; k++) {
    hlmnqc->ampli[k] = 1.;
    hlmnqc->phase[k] = 0.;
  }
  
  for (int k = 0; k < maxk; k++) {
    if (nqc->activemode[k]) {
      hlmnqc->ampli[k] += nqc->a1[k]*nqc->n[k][0] + nqc->a2[k]*nqc->n[k][1] + nqc->a3[k]*nqc->n[k][2]; 
      hlmnqc->phase[k] += nqc->b1[k]*nqc->n[k][3] + nqc->b2[k]*nqc->n[k][4] + nqc->b3[k]*nqc->n[k][5]; 
    }
  }
  
}

/** 
 * Function: eob_wav_ringdown_template
 * -----------------------------------
 *   Ringdown waveform template
 * 
 *   @param[in] x             : time
 *   @param[in] a1,a2,a3,a4   : amplitude coefficients
 *   @param[in] b1,b2,b3,b4   : phase coefficients
 *   @param[in] sigmar,sigmai : Complex QNM frequency
 *   @param[out] psi          : waveform
 */
void eob_wav_ringdown_template(double x, double a1, double a2, double a3, double a4, double b1, double b2, double b3, double b4, double sigmar, double sigmai, double *psi)
{  
  double amp   = ( a1 * tanh(a2*x +a3) + a4 ) ;
  double phase = -b1*log((1. + b3*exp(-b2*x) + b4*exp(-2.*b2*x))/(1.+b3+b4));   
  psi[0] = amp * exp(-sigmar*x); /* amplitude */
  psi[1] = - (phase - sigmai*x); /* phase, minus sign in front by convention */
}

/** 
 * Function: eob_wav_ringdown_v1
 * -----------------------------
 *   Ringdown calculation and match to the dynamics
 * 
 *   @param[in]  dyn : dynamics
 *   @param[out] hlm : waveform
 *   @return         : 0 if successful
 */
int eob_wav_ringdown_v1(Dynamics *dyn, Waveform_lm *hlm)
{

  const double Mbh   = EOBPars->Mbhf;
  const double abh   = EOBPars->abhf;
  const double nu    = EOBPars->nu;  
  const double q     = EOBPars->q;
  const double chi1  = EOBPars->chi1;
  const double chi2  = EOBPars->chi2;
  const double X1    = EOBPars->X1;
  const double X2    = EOBPars->X2;
  const double aK    = EOBPars->a1+EOBPars->a2;

  const double xnu   = (1.-4.*nu);
  const double ooMbh = 1./Mbh;
  /* const double dt = par_get_d("dt"); */	
  const double dt = dyn->dt;
  
  /* double *Omega = dyn->data[EOB_MOMG]; */
  double *Omega = dyn->data[EOB_OMGORB]; /* use this for spin */
  
  /* Note:
     dynsize < size , since the wf has been extended 
     but the two time arrays agree up to dynsize */
  const int dynsize = dyn->size; 
  const int size = hlm->size; 
  double *t = hlm->time;

  /*
    for (int j = 0; j< size ; j++) {
    printf(" %d %.8e", j, hlm->time[j]);
    if (j<dynsize) printf(" %.8e (%.8e)", dyn->time[j], dyn->time[j]-hlm->time[j]);
    printf("\n");
    }
  */
    
  if (VERBOSE) {
    PRFORMi("ringdown_dynamics_size",dynsize);
    PRFORMi("ringdown_waveform_size",size);
  }
    
  const int k21 = 0;
  const int k22 = 1;
  const int k33 = 4;
  
  /** Find peak of Omega */
  /* Assume a monotonically increasing function, 
     start from after the peak */
  int index_pk = dynsize-1;
  double Omega_pk = Omega[index_pk];
  for (int j = dynsize-2; j-- ; ) {
    if (Omega[j] < Omega_pk) 
      break;
    index_pk = j;
    Omega_pk = Omega[j]; 
  }
  
  if (VERBOSE) PRFORMi("ringdown_index_pk",index_pk);
  if (index_pk >= dynsize-2) {
    if (VERBOSE) printf("No omega-maximum found.\n");
  }
  
#if (1)
  
  /* This is a hard-fix that always guarantees the 7 points */
  /* Make sure to comment the following line in main:
     dt_merger_interp = MIN(dt_merger_interp, (dyn->time[size-1] - dyn->tMOmgpeak)/4 ); 
     and uncomment:
     dt_merger_interp = MIN(dt_merger_interp, dyn->dt);
  */
  double *Omega_ptr = &Omega[index_pk-3];  
  double tOmg_pk; /* New interpolated value of the Omega peak */
  const int n = 7; /* USE 7, it seems we need at least 7 points to determine t_Omega_peak properly */
  double tmax = dyn->time[index_pk];   
  
  if ( (index_pk + (n-1)/2) > (dynsize-1) ) { 
    /* Here there are not enough points after the Omega peak 
       We always need 3; we compute what we need by linear extrapolation */
    double Omega_pk_grid[7]; /* Temporary buffer for the 7-point interp */
    const int ni = (index_pk + (n-1)/2) - (dynsize-1) ; /* Pts to extrap, 0 <  ni <= 3 */
    if (VERBOSE) {
      printf("q=%.1f\tchi1=%.3f\tchi2=%.3f\n", q, chi1, chi2);
      printf("%d vs. %d\n", (index_pk + (n-1)/2), (dynsize-1));
      printf("Extrapolate ni=%d points\n",ni);
    } 
    /* errorexit("Not enough points to interpolate.\n"); */

    /* Copy the pts we have */
    for (int j = 0; j < (7-ni); j++) 
      Omega_pk_grid[j] = Omega_ptr[j];
    /* Extrapolate the others */
    if (ni==1) {
      Omega_pk_grid[6] = 2.*Omega_pk_grid[5]-Omega_pk_grid[4];
      //Omega_pk_grid[6] =3.*Omega_pk_grid[5]-3.*Omega_pk_grid[4]+Omega_pk_grid[3];//quadratic, PLEASE CHECK
    } else if (ni==2) {
      Omega_pk_grid[5] = 2.*Omega_pk_grid[4]-Omega_pk_grid[3];
      Omega_pk_grid[6] = 2.*Omega_pk_grid[5]-Omega_pk_grid[4];
    } else if (ni==3) {
      Omega_pk_grid[4] = 2.*Omega_pk_grid[3]-Omega_pk_grid[2];
      Omega_pk_grid[5] = 2.*Omega_pk_grid[4]-Omega_pk_grid[3];
      Omega_pk_grid[6] = 2.*Omega_pk_grid[5]-Omega_pk_grid[4];
    } else errorexit("Wrong counting (ni)\n");
    /* Now we have 7 */
    tOmg_pk = find_max(n, dt, tmax, Omega_pk_grid, NULL);
  } else {    
    /* Here everything is good */
    tOmg_pk = find_max(n, dt, tmax, Omega_ptr, NULL);
  }

  /* Scale peak value by BH mass */
  tOmg_pk *= ooMbh;

#else

  const int n = 7; /* USE 7, it seems we need at least 7 points to determine t_Omega_peak properly */
  if ( (index_pk + (n-1)/2) > (dynsize-1) ) { 
    printf("%.1f\t%.3f\t%.3f\n", q, chi1, chi2);
    printf("%d vs. %d\n", (index_pk + (n-1)/2), (dynsize-1));
    errorexit("Not enough points to interpolate.\n");
  }
  
  double tmax = dyn->time[index_pk];  
  double *Omega_ptr = &Omega[index_pk-3];
  double tOmg_pk = find_max(n, dt, tmax, Omega_ptr, NULL);
  tOmg_pk *= ooMbh;

#endif

  if (VERBOSE) PRFORMd("ringdown_Omega_pk",Omega_pk);
  if (VERBOSE) PRFORMd("ringdown_Omega_ptr",Omega_ptr[0]);
  if (VERBOSE) PRFORMd("ringdown_tmax",tmax);
  if (VERBOSE) PRFORMd("ringdown_tOmg_pk",tOmg_pk/ooMbh);
  if (VERBOSE) PRFORMd("ringdown_tOmg_pk",tOmg_pk);
  
  /** Merger time t_max(A22) */
  double DeltaT_nqc = eob_nqc_timeshift(nu, chi1);
  double tmrg[KMAX], tmatch[KMAX], dtmrg[KMAX];

  /* nonspinning case */
  /* tmrg[k22]  = tOmg_pk-3./Mbh; */ /* OLD */       
  double tmrgA22 = tOmg_pk-(DeltaT_nqc + 2.)/Mbh;
  if (VERBOSE) PRFORMd("ringdown_tmrgA22",tmrgA22);

  for (int k=0; k<KMAX; k++) {
    tmrg[k] = tmrgA22;
  }
  
  /* The following values are the difference between the time of the peak of
     the 22 waveform and the other modes. */
  /* These specific values refer to the 21 and 33 in the
     nonspinning case. They are different in the spinning case, which
     is however not implemented. These are here only as placeholder */
  dtmrg[k21] = 5.70364338 + 1.85804796*xnu  + 4.0332262*xnu*xnu; //k21
  dtmrg[k33] = 4.29550934 - 0.85938*xnu;                         //k33
  tmrg[k21]  = tmrg[k22] + dtmrg[k21]/Mbh;     // t_max(A21) => peak of 21 mode
  tmrg[k33]  = tmrg[k22] + dtmrg[k33]/Mbh;     // t_max(A33) => peak of 33 mode
	  
  // for (int k=0; k<KMAX; k++) {
  //   tmrg[k] = tmrgA22 + dtmrg[k]/Mbh;
  // }	  
	
  /** Postmerger-Ringdown matching time */
  for (int k=0; k<KMAX; k++) {
    tmatch[k] = 2.*ooMbh + tmrg[k];
  }
  
  /** Compute QNM */
  double sigma[2][KMAX];
  double a1[KMAX], a2[KMAX], a3[KMAX], a4[KMAX];
  double b1[KMAX], b2[KMAX], b3[KMAX], b4[KMAX];

  QNMHybridFitCab(nu, X1, X2, chi1, chi2, aK,  Mbh, abh,  
		  a1, a2, a3, a4, b1, b2, b3, b4, 
		  sigma[0],sigma[1]);
  
  /** Define a time vector for each multipole, scale by mass
      Ringdown of each multipole has its own starting time */
  double *t_lm[KMAX];
  for (int k=0; k<KMAX; k++) {
    t_lm[k] =  malloc ( size * sizeof(double) );
    for (int j = 0; j < size; j++ ) {  
      t_lm[k][j] = t[j] * ooMbh;
    }
  }  
  
  /** Find attachment index */
  int idx[KMAX];
  for (int k = 0; k < KMAX; k++) {
    for (int j = size-1; j-- ; ) {  
      if (t_lm[k][j] < tmatch[k]) {
	      idx[k] = j - 1;
	      break;
      }
    }
  }
  
  /** Compute Ringdown waveform for t>=tmatch */
  double t0, tm, psi[2];
  double Deltaphi[KMAX];
  for (int k = 0; k < KMAX; k++) {
    if(hlm->kmask[k]){

      /* Calculate Deltaphi */
      t0 = t_lm[k][idx[k]] - tmrg[k]; 
      eob_wav_ringdown_template(t0, a1[k], a2[k], a3[k], a4[k], b1[k], b2[k], b3[k], b4[k], sigma[0][k], sigma[1][k], psi);
      Deltaphi[k] = psi[1] - hlm->phase[k][idx[k]];
      /* Compute and attach ringdown */
      for (int j = idx[k]; j < size ; j++ ) {   
        tm = t_lm[k][j] - tmrg[k];
        eob_wav_ringdown_template(tm, a1[k], a2[k], a3[k], a4[k], b1[k], b2[k], b3[k], b4[k], sigma[0][k], sigma[1][k], psi);
        hlm->phase[k][j] = psi[1] - Deltaphi[k];
        hlm->ampli[k][j] = psi[0];
      }
    }
  }

  /** Free mem. */
  for (int k=0; k<KMAX; k++) {
    free(t_lm[k]);
  }
  
  return 0;

}
	
/** 
 * Function: eob_wav_ringdown_HM
 * -----------------------------
 *   Ringdown calculation and match to the dynamics
 *   This refers to the higher modes paper: arXiv:2001.09082
 * 
 *   @param[in]  dyn : dynamics
 *   @param[out] hlm : waveform
 *   @return         : 0 if successful
*/
int eob_wav_ringdown_HM(Dynamics *dyn, Waveform_lm *hlm)
{
  
  const double Mbh   = EOBPars->Mbhf;
  const double abh   = EOBPars->abhf;
  const double nu    = EOBPars->nu;
  const double q    = EOBPars->q;
  const double chi1  = EOBPars->chi1;
  const double chi2  = EOBPars->chi2;
  const double X1    = EOBPars->X1;
  const double X2    = EOBPars->X2;
  const double aK    = EOBPars->a1+EOBPars->a2;
	
  const double xnu   = (1.-4.*nu);
  const double ooMbh = 1./Mbh;
  /* const double dt = par_get_d("dt"); */	
  const double dt = 0.5;//dyn->dt;
	  
  /* double *Omega = dyn->data[EOB_MOMG]; */
  double *Omega = dyn->data[EOB_OMGORB]; /* use this for spin */
	  
  /* Note:
     dynsize < size , since the wf has been extended 
     but the two time arrays agree up to dynsize */
  const int dynsize = dyn->size; 
  const int size = hlm->size; 
  double *t = hlm->time;
	    
  if (VERBOSE) {
    PRFORMi("ringdown_dynamics_size",dynsize);
    PRFORMi("ringdown_waveform_size",size);
  }
  
  /** Find peak of Omega */
  /* Assume a monotonically increasing function, 
     start from after the peak */
  int index_pk = dynsize-1;
  double Omega_pk = Omega[index_pk];
  for (int j = dynsize-2; j-- ; ) {
    if (Omega[j] < Omega_pk) 
      break;
    index_pk = j;
    Omega_pk = Omega[j]; 
  }
  
  if (VERBOSE) PRFORMi("ringdown_index_pk",index_pk);
  if (index_pk >= dynsize-2) {
    if (VERBOSE) printf("No omega-maximum found.\n");
  }
  
  double tOmg_pk = dyn->time[index_pk]*ooMbh;

  if (VERBOSE) PRFORMd("ringdown_Omega_pk",Omega_pk);
  if (VERBOSE) PRFORMd("ringdown_tOmg_pk",tOmg_pk/ooMbh);
  if (VERBOSE) PRFORMd("ringdown_tOmg_pk",tOmg_pk);
	  
  /** Merger time t_max(A22) */
  double DeltaT_nqc = eob_nqc_timeshift(nu, chi1);
  double tmrg[KMAX], tmatch[KMAX], dtmrg[KMAX];
 
  /* nonspinning case */
  double tmrgA22 = tOmg_pk-(DeltaT_nqc + 2.)/Mbh;
  if (VERBOSE) PRFORMd("ringdown_tmrgA22",tmrgA22);
  if (tmrgA22 < 0.0) {
    printf("ERROR: the A22 peak time is negative, dynamics is too short\n");
    return 1;
  }
  
  /* The following values are the difference between the time of the peak of
     the 22 waveform and the other modes. 
     For modes in knqcpeak22 we impose tmrg[k] = tmrg[1], to attach the ringdown there
  */
  int modesatpeak22[KMAX]; 
  set_multipolar_idx_mask (modesatpeak22, KMAX, EOBPars->knqcpeak22, EOBPars->knqcpeak22_size, 0);
  eob_nqc_deltat_lm(dyn, dtmrg);
  for (int k=0; k<KMAX; k++) {
    tmrg[k] = tmrgA22 + dtmrg[k]/Mbh;
    if (modesatpeak22[k]) tmrg[k] = tmrgA22;
  }	  

  /** Postmerger-Ringdown matching time */
  int idx[KMAX];
  for (int k = 0; k < KMAX; k++) {
    if(hlm->kmask[k]){
      int j  = size-1;
      idx[k] = size-1;
      for (j = size-1; j>0; j--) {  
	      if ( (t[j] < tmrg[k]*Mbh) || (fabs(t[j] - tmrg[k]*Mbh)<1.e-8) ) {
	        break;
	      }
      }
      idx[k]    = j;
      tmatch[k] = (t[idx[k]])*ooMbh;	    
    }
  }
  
  /** Compute QNM */
  double sigma[2][KMAX];
  double a1[KMAX], a2[KMAX], a3[KMAX], a4[KMAX];
  double b1[KMAX], b2[KMAX], b3[KMAX], b4[KMAX];

  QNMHybridFitCab_HM(nu, X1, X2, chi1, chi2, aK,  Mbh, abh,  
		     a1, a2, a3, a4, b1, b2, b3, b4, 
		     sigma[0],sigma[1]);
  
  /* Overwrite the modes attached at the peak */
  QNMHybridFitCab_HM_Pompili23(nu, X1, X2, chi1, chi2, aK,  Mbh, abh,
                                a1, a2, a3, a4, b1, b2, b3, b4,
                                sigma[0], sigma[1]
                              );
  /** Define a time vector for each multipole, scale by mass
      Ringdown of each multipole has its own starting time */
  double *t_lm[KMAX];
  for (int k=0; k<KMAX; k++) {
    t_lm[k] =  malloc ( size * sizeof(double) );
    for (int j = 0; j < size; j++ ) {  
      t_lm[k][j] = t[j] * ooMbh;   // CHECKME
    }
  }
	  
  /** Compute Ringdown waveform for t>=tmatch */
  double t0, tm, psi[2];
  double Deltaphi[KMAX];
  int n0 = 2./dt*ooMbh; // CHECKME: I think this is to get the t + 2, we need to remove it for the modes in modesatpeak22
  int index_rng;
	  
  for (int k = 0; k < KMAX; k++) {
    double fact = 1.; // this is set to 1/Mbh below for the modes attached at the peak of the (2,2)
    if(hlm->kmask[k]){

      /* Ringdown attachment index */      
      index_rng = idx[k]+n0;
      if (modesatpeak22[k]){
          index_rng = idx[k];
          fact      = ooMbh;
      }
      if (index_rng > dynsize -1) index_rng = dynsize - 1;
      
      /* Calculate Deltaphi */
      t0  = t_lm[k][index_rng] - tmatch[k];
      t0 /= fact;
      eob_wav_ringdown_template(t0, a1[k], a2[k], a3[k], a4[k], b1[k], b2[k], b3[k], b4[k], sigma[0][k], sigma[1][k], psi);
      
      Deltaphi[k] = psi[1] - hlm->phase[k][index_rng];
      
      /* Compute and attach ringdown */
      for (int j = index_rng-1; j < size ; j++ ) {
        tm  = t_lm[k][j] - tmatch[k];
        tm /= fact;

        //printf("j =%d, k=%d, tlm_k_rd = %.8e, tm = %.8e\n", j,k, t_lm[k][j], tm);

        eob_wav_ringdown_template(tm, a1[k], a2[k], a3[k], a4[k], b1[k], b2[k], b3[k], b4[k], sigma[0][k], sigma[1][k], psi);
        hlm->phase[k][j] = psi[1] - Deltaphi[k];
        hlm->ampli[k][j] = psi[0];

        if(nNegAmp[k]==1) {
          hlm->ampli[k][j] = -hlm->ampli[k][j];
        }
      }
    }
  } 
	  
  /** Free mem. */
  for (int k=0; k<KMAX; k++) {
    free(t_lm[k]);
  }
  
  return 0;

}

/** 
 * Function: eob_wav_hlm
 * ---------------------
 *   Main routine for factorized EOB waveform,
 *   quasi-circular case
 * 
 *   @param[in]  dyn : dynamics
 *   @param[out] hlm : waveform
*/
void eob_wav_hlm_circ(Dynamics *dyn, Waveform_lm_t *hlm)
{
  
  const double t   = dyn->t;
  
  const double nu = EOBPars->nu;  
  const double chi1 = EOBPars->chi1;  
  const double chi2 = EOBPars->chi2;  
  const double a1 = EOBPars->a1;  
  const double a2 = EOBPars->a2;  
  const double X1 = EOBPars->X1;  
  const double X2 = EOBPars->X2;  
  const double C_Q1 = EOBPars->C_Q1;  
  const double C_Q2 = EOBPars->C_Q2;  
  const int usetidal = EOBPars->use_tidal;
  const int usespins = EOBPars->use_spins;
  const int usespeedytail = EOBPars->use_speedytail;
  const double X12 = X1-X2; /* sqrt(1-4nu) */

  const double phi = dyn->phi; 
  const double r   = dyn->r;
  const double pph = dyn->pphi;
  const double prstar = dyn->prstar;
  const double Omega  = dyn->Omg;
  const double ddotr  = dyn->ddotr;
  const double H      = dyn->H;
  const double Heff   = dyn->Heff;
  const double jhat   = dyn->jhat;
  const double rw     = dyn->r_omega;

  hlm->time = t;

  /** Source term */
  double source[] = {
    jhat,Heff,
    Heff,jhat,Heff,
    jhat,Heff,jhat,Heff,
    Heff,jhat,Heff,jhat,Heff,
    jhat,Heff,jhat,Heff,jhat,Heff,
    Heff,jhat,Heff,jhat,Heff,jhat,Heff,
    jhat,Heff,jhat,Heff,jhat,Heff,jhat,Heff
  };
  
  /** Newtonian waveform */
  Waveform_lm_t hNewt;
  eob_wav_hlmNewt(rw,Omega,phi,nu, &hNewt);

  if (usetidal) {
    /* Need to correct some of the m=odd modes. 
       The Newtonian factor has a different normalization when entering the point-mass 
       and the tidal term. The factor X12 = sqrt*1-4nu) is re-introduced in the point-mass term 
       in eob_wav_hlm() */
    double vphi3 = gsl_pow_int(rw*Omega,3); 
    hNewt.ampli[0] = ChlmNewt_ampli[0] * vphi3;
    hNewt.ampli[2] = ChlmNewt_ampli[2] * vphi3;
    hNewt.ampli[4] = ChlmNewt_ampli[4] * vphi3; 
    double p4_vphi5 = (2.*nu-1) * gsl_pow_int(rw*Omega,5); 
    hNewt.ampli[5]  = ChlmNewt_ampli[5]  * p4_vphi5;
    hNewt.ampli[7]  = ChlmNewt_ampli[7]  * p4_vphi5; 
    hNewt.ampli[9]  = ChlmNewt_ampli[9]  * p4_vphi5; 
    hNewt.ampli[11] = ChlmNewt_ampli[11] * p4_vphi5;
    hNewt.ampli[13] = ChlmNewt_ampli[13] * p4_vphi5; 
  }

  if (usespins) {
    /* Special treatment when spin is on because of the singularity in the sqrt(1-4*nu) 
       for m=odd mode and nu=1/4. See discussion in 
       Damour & Nagar, PRD 90, 044018, Sec. 4, Eq.(89). 
       This is not done for multipoles l>4 because no spinning information is included there */ 
    
    double vphi3 = gsl_pow_int(rw*Omega,3);
    hNewt.ampli[0] = ChlmNewt_ampli[0] * vphi3; /* (2,1) */
    
    double p4_vphi5 = (2.*nu-1) * gsl_pow_int(rw*Omega,5);
    hNewt.ampli[9]  = ChlmNewt_ampli[9]  * p4_vphi5 * X12; /* (5,1) */
    hNewt.ampli[11] = ChlmNewt_ampli[11] * p4_vphi5 * X12; /* (5,3) */

    if (!(EOBPars->use_flm == USEFLM_HM || EOBPars->use_flm == USEFLM_HM_4PN22)) {
      hNewt.ampli[2] = ChlmNewt_ampli[2] * vphi3; /* (3,1) */
      hNewt.ampli[4] = ChlmNewt_ampli[4] * vphi3; /* (3,3) */
      hNewt.ampli[5]  = ChlmNewt_ampli[5]  * p4_vphi5; /* (4,1) */
      hNewt.ampli[7]  = ChlmNewt_ampli[7]  * p4_vphi5; /* (4,3) */
      hNewt.ampli[13] = ChlmNewt_ampli[13] * p4_vphi5 * X12; /* (5,5) */
    }

  }

  /** Compute corrections */
  double rholm[KMAX], flm[KMAX];
  double x = SQ(rw*Omega);
  if (usespins){
    /* eob_wav_flm_s_old(x, nu, X1,X2,chi1,chi2,a1,a2,C_Q1,C_Q2,usetidal,rholm,flm); */ 
    eob_wav_flm_s(x, nu, X1,X2,chi1,chi2,a1,a2,C_Q1,C_Q2,usetidal,rholm,flm); 
  } else {
    /* eob_wav_flm_old(x, nu, rholm,flm); */
    eob_wav_flm(x, nu, rholm,flm);
  }
  
  /** Computing the tail */
#define RTAIL (1.213061319425267e+00)
  const double Hreal = H * nu;
  Waveform_lm_t tlm;
  if (usespeedytail) {
    eob_wav_speedyTail(Omega,Hreal, RTAIL, &tlm); 
  } else {
    eob_wav_hhatlmTail(Omega,Hreal, RTAIL, &tlm); 
  }
  
  /** Residual phase corrections delta_{lm} */
  double dlm[KMAX];
  eob_wav_deltalm(Hreal, Omega, nu, dlm); 

  /** Point-mass h_lm */
  for (int k = 0; k < KMAX; k++) {
    hlm->ampli[k] =  hNewt.ampli[k] * flm[k] * source[k] * tlm.ampli[k];
    hlm->phase[k] = -( hNewt.phase[k] + tlm.phase[k] + dlm[k]); /* Minus sign by convention */
  }
  
  /** NQC */
  if (!(EOBPars->nqc_coefs_hlm == NQC_HLM_NONE) &&
      !(EOBPars->nqc_coefs_hlm == NQC_HLM_COMPUTE)) {

    /* Add NQC correction */    
    Waveform_lm_t hNQC; 
    eob_wav_hlmNQC(nu,r,prstar,Omega,ddotr, NQC->hlm, &hNQC); 
    const int maxk = MIN(KMAX, NQC->hlm->maxk+1);
    for (int k = 0; k < maxk; k++) {
      if (NQC->hlm->activemode[k]) {
	hlm->ampli[k] *= hNQC.ampli[k];
	hlm->phase[k] -= hNQC.phase[k];
      }
    }
    
  }
  
  if (usetidal) {   
    /** Tidal contribution */
    double hlmtidal[KMAX];
    eob_wav_hlmTidal(x, dyn, hlmtidal);
    if( !(usespins) ) { 
      /* Correct normalization of point-mass wave for some of the m=odd modes */
      hlm->ampli[0] *= X12;
      hlm->ampli[2] *= X12;
      hlm->ampli[4] *= X12;
      hlm->ampli[5] *= X12;
      hlm->ampli[7] *= X12;
      hlm->ampli[9] *= X12;
      hlm->ampli[11] *= X12;
      hlm->ampli[13] *= X12;
    }
    /* Add tidal contribution to waveform */
    for (int k = 0; k < KMAX; k++) {
      hlm->ampli[k] += (hNewt.ampli[k] * tlm.ampli[k] * hlmtidal[k]);
    }
  }
  
}

/** 
 * Function: eob_wav_hlm_ecc
 * -------------------------
 *   Main routine for factorized EOB waveform,
 *   generic case
 * 
 *   @param[in]  dyn : dynamics
 *   @param[out] hlm : waveform
*/
void eob_wav_hlm_ecc(Dynamics *dyn, Waveform_lm_t *hlm)
{
  
  const double nu = EOBPars -> nu;
  const double chi1 = EOBPars -> chi1;
  const double chi2 = EOBPars -> chi2;
  const double a1 = EOBPars -> a1;
  const double a2 = EOBPars -> a2;
  const double X1 = EOBPars -> X1;
  const double X2 = EOBPars -> X2;
  const double C_Q1 = EOBPars -> C_Q1;
  const double C_Q2 = EOBPars -> C_Q2;
  const double ecc = EOBPars -> ecc;
  const double r_hyp = EOBPars -> r_hyp;
  const int usetidal = EOBPars -> use_tidal;
  const int usespins = EOBPars -> use_spins;
  const int usespeedytail = EOBPars -> use_speedytail;
  const double X12 = X1 - X2; /* sqrt (1 - 4 nu)*/

  const double t   = dyn->t;
  const double phi = dyn->phi; 
  const double r   = dyn->r;
  const double pph = dyn->pphi;
  const double prstar = dyn->prstar;
  const double Omega  = dyn->Omg;
  const double ddotr  = dyn->ddotr;
  const double H      = dyn->H;
  const double Heff   = dyn->Heff;
  const double jhat   = dyn->jhat;
  const double rw     = dyn->r_omega;
  hlm->time = t;

  /** Source term */
  double source[] = {
    jhat,Heff,
    Heff,jhat,Heff,
    jhat,Heff,jhat,Heff,
    Heff,jhat,Heff,jhat,Heff,
    jhat,Heff,jhat,Heff,jhat,Heff,
    Heff,jhat,Heff,jhat,Heff,jhat,Heff,
    jhat,Heff,jhat,Heff,jhat,Heff,jhat,Heff
  };
  
  /** Newtonian waveform */
  Waveform_lm_t hNewt;
  eob_wav_hlmNewt_ecc(dyn, &hNewt);

  /** Compute corrections */
  double rholm[KMAX], flm[KMAX];
  double x = SQ(rw*Omega);
  if (usespins){
    /* eob_wav_flm_s_old(x, nu, X1,X2,chi1,chi2,a1,a2,C_Q1,C_Q2,usetidal,rholm,flm); */ 
    eob_wav_flm_s(x, nu, X1,X2,chi1,chi2,a1,a2,C_Q1,C_Q2,usetidal,rholm,flm); 
  } else {
    /* eob_wav_flm_old(x, nu, rholm,flm); */
    eob_wav_flm(x, nu, rholm,flm);
  }
  
  /** Computing the tail */
#define RTAIL (1.213061319425267e+00)
  const double Hreal = H * nu;
  Waveform_lm_t tlm;
  if (usespeedytail) {
    eob_wav_speedyTail(Omega,Hreal, RTAIL, &tlm); 
  } else {
    eob_wav_hhatlmTail(Omega,Hreal, RTAIL, &tlm); 
  }
  
  /** Residual phase corrections delta_{lm} */
  double dlm[KMAX];
  eob_wav_deltalm(Hreal, Omega, nu, dlm); 

  /** Point-mass h_lm */
  for (int k = 0; k < KMAX; k++) {
    hlm->ampli[k] =  hNewt.ampli[k] * flm[k] * source[k] * tlm.ampli[k];
    hlm->phase[k] = -( hNewt.phase[k] + tlm.phase[k] + dlm[k]); /* Minus sign by convention */    
  }
  
  /** NQC */
  if (!(EOBPars->nqc_coefs_hlm == NQC_HLM_NONE) &&
      !(EOBPars->nqc_coefs_hlm == NQC_HLM_COMPUTE)) {

    /* Add NQC correction */    
    Waveform_lm_t hNQC; 
    eob_wav_hlmNQC(nu,r,prstar,Omega,ddotr, NQC->hlm, &hNQC); 
    const int maxk = MIN(KMAX, NQC->hlm->maxk+1);
    
    for (int k = 0; k < maxk; k++) {
      if (NQC->hlm->activemode[k]) {
	hlm->ampli[k] *= hNQC.ampli[k];
	hlm->phase[k] -= hNQC.phase[k];
      }
    }
    
  }
  
  if (usetidal) {   
    /** Tidal contribution */
    double hlmtidal[KMAX];
    eob_wav_hlmTidal(x, dyn, hlmtidal);
    if( !(usespins) ) { 
      /* Correct normalization of point-mass wave for some of the m=odd modes */
      hlm->ampli[0] *= X12;
      hlm->ampli[2] *= X12;
      hlm->ampli[4] *= X12;
      hlm->ampli[5] *= X12;
      hlm->ampli[7] *= X12;
      hlm->ampli[9] *= X12;
      hlm->ampli[11] *= X12;
      hlm->ampli[13] *= X12;
    }
    /* Add tidal contribution to waveform */
    for (int k = 0; k < KMAX; k++) {
      hlm->ampli[k] += (hNewt.ampli[k] * tlm.ampli[k] * hlmtidal[k]);
    }
  }
  
}


/** 
 * Function: eob_wav_hlm_ecc
 * -------------------------
 *   Main routine for factorized EOB waveform,
 *   generic case with sigmoid
 * 
 *   @param[in]  dyn : dynamics
 *   @param[out] hlm : waveform
*/
void eob_wav_hlm_ecc_sigmoid(Dynamics *dyn, Waveform_lm_t *hlm)
{
  
  const double nu = EOBPars -> nu;
  const double chi1 = EOBPars -> chi1;
  const double chi2 = EOBPars -> chi2;
  const double a1 = EOBPars -> a1;
  const double a2 = EOBPars -> a2;
  const double X1 = EOBPars -> X1;
  const double X2 = EOBPars -> X2;
  const double C_Q1 = EOBPars -> C_Q1;
  const double C_Q2 = EOBPars -> C_Q2;
  const double ecc = EOBPars -> ecc;
  const double r_hyp = EOBPars -> r_hyp;
  const int usetidal = EOBPars -> use_tidal;
  const int usespins = EOBPars -> use_spins;
  const int usespeedytail = EOBPars -> use_speedytail;
  const double X12 = X1 - X2; /* sqrt (1 - 4 nu)*/

  const double t   = dyn->t;
  const double phi = dyn->phi; 
  const double r   = dyn->r;
  const double pph = dyn->pphi;
  const double prstar = dyn->prstar;
  const double Omega  = dyn->Omg;
  const double ddotr  = dyn->ddotr;
  const double H      = dyn->H;
  const double Heff   = dyn->Heff;
  const double jhat   = dyn->jhat;
  const double rw     = dyn->r_omega;
  hlm->time = t;

  /** Source term */
  double source[] = {
    jhat,Heff,
    Heff,jhat,Heff,
    jhat,Heff,jhat,Heff,
    Heff,jhat,Heff,jhat,Heff,
    jhat,Heff,jhat,Heff,jhat,Heff,
    Heff,jhat,Heff,jhat,Heff,jhat,Heff,
    jhat,Heff,jhat,Heff,jhat,Heff,jhat,Heff
  };
  
  /** Newtonian waveform */
  Waveform_lm_t hNewt;
  eob_wav_hlmNewt_ecc_sigmoid(dyn, &hNewt);

  /** Compute corrections */
  double rholm[KMAX], flm[KMAX];
  double x = SQ(rw*Omega);
  if (usespins){
    /* eob_wav_flm_s_old(x, nu, X1,X2,chi1,chi2,a1,a2,C_Q1,C_Q2,usetidal,rholm,flm); */ 
    eob_wav_flm_s(x, nu, X1,X2,chi1,chi2,a1,a2,C_Q1,C_Q2,usetidal,rholm,flm); 
  } else {
    /* eob_wav_flm_old(x, nu, rholm,flm); */
    eob_wav_flm(x, nu, rholm,flm);
  }
  
  /** Computing the tail */
#define RTAIL (1.213061319425267e+00)
  const double Hreal = H * nu;
  Waveform_lm_t tlm;
  if (usespeedytail) {
    eob_wav_speedyTail(Omega,Hreal, RTAIL, &tlm); 
  } else {
    eob_wav_hhatlmTail(Omega,Hreal, RTAIL, &tlm); 
  }
  
  /** Residual phase corrections delta_{lm} */
  double dlm[KMAX];
  eob_wav_deltalm(Hreal, Omega, nu, dlm); 

  /** Point-mass h_lm */
  for (int k = 0; k < KMAX; k++) {
    hlm->ampli[k] =  hNewt.ampli[k] * flm[k] * source[k] * tlm.ampli[k];
    hlm->phase[k] = -( hNewt.phase[k] + tlm.phase[k] + dlm[k]); /* Minus sign by convention */    
  }
  
  /** NQC */
  if (!(EOBPars->nqc_coefs_hlm == NQC_HLM_NONE) &&
      !(EOBPars->nqc_coefs_hlm == NQC_HLM_COMPUTE)) {

    /* Add NQC correction */
    Waveform_lm_t hNQC;
    double tOmg_pk = dyn->tOmg_pk;
    eob_wav_hlmNQC_ecc_sigmoid(nu,r,prstar,Omega,ddotr, t, tOmg_pk, NQC->hlm, &hNQC); 
    const int maxk = MIN(KMAX, NQC->hlm->maxk+1);
  
    for (int k = 0; k < maxk; k++) {
      if (NQC->hlm->activemode[k]) {
        hlm->ampli[k] *= hNQC.ampli[k];
        hlm->phase[k] -= hNQC.phase[k];
      }
    }
    
  }
  
  if (usetidal) {   
    /** Tidal contribution */
    double hlmtidal[KMAX];
    eob_wav_hlmTidal(x, dyn, hlmtidal);
    if( !(usespins) ) { 
      /* Correct normalization of point-mass wave for some of the m=odd modes */
      hlm->ampli[0] *= X12;
      hlm->ampli[2] *= X12;
      hlm->ampli[4] *= X12;
      hlm->ampli[5] *= X12;
      hlm->ampli[7] *= X12;
      hlm->ampli[9] *= X12;
      hlm->ampli[11] *= X12;
      hlm->ampli[13] *= X12;
    }
    /* Add tidal contribution to waveform */
    for (int k = 0; k < KMAX; k++) {
      hlm->ampli[k] += (hNewt.ampli[k] * tlm.ampli[k] * hlmtidal[k]);
    }
  }
  
}


/**
  * Function: eob_wav_flm_old
  * -------------------------
  *   Computes the resummed amplitudes in the general nu-dependent case.
  *   Refs:
  *   * Damour, Iyer & Nagar, PRD 79, 064004 (2009)     [theory]
  *   * Fujita & Iyer, PRD 82, 044051 (2010)            [test-mass 5.5PN]
  *   * Damour, Nagar & Bernuzzi, PRD 87, 084035 (2013) [complete information]
  *   @note Old version: (i) coefficients computed all the times, (i) (3,1) and (3,3) mode do not contain nu-corrections.
  * 
  *   @param[in]  x     : x = (M omega)^2
  *   @param[in]  nu    : symmetric mass ratio
  *   @param[out] rholm : resummed amplitudes
  *   @param[out] flm   : resummed amplitudes
*/
void eob_wav_flm_old(double x,double nu, double *rholm, double *flm)
{
  
  /** Shorthands */
  const double x2  = x*x;
  const double x3  = x*x2;
  const double x4  = x*x3;
  const double x5  = x*x4;
  const double nu2 = nu*nu;
  const double nu3 = nu*nu2;
  const double nu4 = nu*nu3;
  
  /** Compute EulerLogs */
  const double el1 = Eulerlog(x,1);
  const double el2 = Eulerlog(x,2);
  const double el3 = Eulerlog(x,3);
  const double el4 = Eulerlog(x,4);
  const double el5 = Eulerlog(x,5);
  const double el6 = Eulerlog(x,6);
  const double el7 = Eulerlog(x,7);
  
  /** l=2 ------------------------------------------------------------------
   *  (2,1) */
  rholm[0] = 1.                                                                           +
    (-1.0535714285714286 + 0.27380952380952384 *nu                          )*x  +
    (-0.8327841553287982 - 0.7789824263038548  *nu + 0.13116496598639457*nu2)*x2 +
    (2.9192806270460925  - 1.019047619047619   *el1                         )*x3 +
    (-1.28235780892213   + 1.073639455782313   *el1                         )*x4 +
    (-3.8466571723355227 + 0.8486467106683944  *el1                     )*PMTERMS_eps*x5 ;
  
  /** (2,2) */
  rholm[1] = 1.                                                                                                                         +
    (-1.0238095238095237 + 0.6547619047619048*nu                                                                          )*x  +
    (-1.94208238851096   - 1.5601379440665155*nu + 0.4625614134542706*nu2                                                 )*x2 +
    (12.736034731834051  - 2.902228713904598 *nu - 1.9301558466099282*nu2 + 0.2715020968103451*nu3 - 4.076190476190476*el2)*x3 +
    (-2.4172313935587004 + 4.173242630385488 *el2                                                                         )*x4 +
    (-30.14143102836864  + 7.916297736025627 *el2                                                                         )*x5 ;
  
  /** l=3 ------------------------------------------------------------------
   *  (3,1) */
  rholm[2] = 1. + (-0.7222222222222222 - 0.2222222222222222*nu)*x + (0.014169472502805836 - 0.9455667789001122*nu - 0.46520763187429853*nu2)*x2 + x3*(1.9098284139598072 - 0.4126984126984127*el1) + x4*(0.5368150316615179 + 0.2980599647266314*el1) + PMTERMS_eps*x5*(1.4497991763035063 - 0.0058477188106817735*el1);
  /** (3,2) */
  rholm[3] = 1. + (0.003703703703703704*(328. - 1115.*nu + 320.*nu2)*x)/(-1. + 3.*nu) + (6.235191420376606e-7*(-1.444528e6 + 8.050045e6*nu - 4.725605e6*nu2 - 2.033896e7*nu3 + 3.08564e6*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + x3*(6.220997955214429 - 1.6507936507936507*el2) + PMTERMS_eps*x4*(-3.4527288879001268 + 2.005408583186361*el2);
  /** (3,3) */
  rholm[4] = 1. + (-1.1666666666666667 + 0.6666666666666666*nu)*x + (-1.6967171717171716 - 1.8797979797979798*nu + 0.45151515151515154*nu2)*x2 + x3*(14.10891386831863 - 3.7142857142857144*el3) + x4*(-6.723375314944128 + 4.333333333333333*el3) + PMTERMS_eps*x5*(-29.568699895427518 + 6.302092352092352*el3);
  
  /** l=4 ------------------------------------------------------------------
   *  (4,1) */
  rholm[5] = 1. + (0.001893939393939394*(602. - 1385.*nu + 288.*nu2)*x)/(-1. + 2.*nu) - 0.36778992787515513*x2 + x3*(0.6981550175535535 - 0.2266955266955267*el1) + PMTERMS_eps*x4*(-0.7931524512893319 + 0.2584672482399755*el1);
  /** (4,2) */
  rholm[6] = 1. + (0.0007575757575757576*(1146. - 3530.*nu + 285.*nu2)*x)/(-1. + 3.*nu) - (3.1534122443213353e-9*(1.14859044e8 - 2.95834536e8*nu - 1.204388696e9*nu2 + 3.04798116e9*nu3 + 3.79526805e8*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + 4.550378418934105e-12*x3*(8.48238724511e11 - 1.9927619712e11*el2) + PMTERMS_eps*x4*(-0.6621921297263365 + 0.787251738160829*el2);
  /** (4,3) */
  rholm[7] = 1. + (0.005681818181818182*(222. - 547.*nu + 160.*nu2)*x)/(-1. + 2.*nu) - 0.9783218202252293*x2 + PMTERMS_eps*(x3*(8.519456157072423 - 2.0402597402597404*el3) +      x4*(-5.353216984886716 + 2.5735094451003544*el3));
  /** (4,4) */
  rholm[8] = 1. + (0.0007575757575757576*(1614. - 5870.*nu + 2625.*nu2)*x)/(-1. + 3.*nu) + (3.1534122443213353e-9*(-5.11573572e8 + 2.338945704e9*nu - 3.13857376e8*nu2 - 6.733146e9*nu3 + 1.252563795e9*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + x3*(15.108111214795123 - 3.627128427128427*el4) + PMTERMS_eps*x4*(-8.857121657199649 + 4.434988849534304*el4);
  
  /** l=5 ------------------------------------------------------------------
   *  (5,1) */
  rholm[9] = 1. + (0.002564102564102564*(319. - 626.*nu + 8.*nu2)*x)/(-1. + 2.*nu) - 0.1047896120973044*x2 + PMTERMS_eps*(x3*(0.642701885362399 - 0.14414918414918415*el1) +      x4*(-0.07651588046467575 + 0.11790664036817883*el1));
  /** (5,2) */
  rholm[10] = 1. + (0.00007326007326007326*(-15828. + 84679.*nu - 104930.*nu2 + 21980.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + PMTERMS_eps*(- 0.4629337197600934*x2 +      x3*(2.354458371550237 - 0.5765967365967366*el2));
  /** (5,3) */
  rholm[11] = 1. + (0.002564102564102564*(375. - 850.*nu + 176.*nu2)*x)/(-1. + 2.*nu) - 0.5788010707241477*x2 + PMTERMS_eps*(x3*(5.733973288504755 - 1.2973426573426574*el3) +      x4*(-1.9573287625526001 + 1.2474448628294783*el3));
  /** (5,4) */
  rholm[12] = 1. + (0.00007326007326007326*(-17448. + 96019.*nu - 127610.*nu2 + 33320.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + PMTERMS_eps*(- 1.0442142414362194*x2 +   x3*(10.252052781721588 - 2.3063869463869464*el4));
    /**  (5,5) */
  rholm[13] = 1. + (0.002564102564102564*(487. - 1298.*nu + 512.*nu2)*x)/(-1. + 2.*nu) - 1.5749727622804546*x2 + PMTERMS_eps*(x3*(15.939827047208668 - 3.6037296037296036*el5) +      x4*(-10.272578060123237 + 4.500041838503377*el5));
  
  /** l=6 ------------------------------------------------------------------
   *  (6,1) */
  rholm[14] = 1. + (0.006944444444444444*(-161. + 694.*nu - 670.*nu2 + 124.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 0.29175486850885135*x2 + x3*(0.21653486654395454 - 0.10001110001110002*el1));
  /** (6,2) */
  rholm[15] = 1. + (0.011904761904761904*(-74. + 378.*nu - 413.*nu2 + 49.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + PMTERMS_eps*( - 0.24797525070634313*x2 + x3*(1.7942694138754138 - 0.40004440004440006*el2));
  /** (6,3) */
  rholm[16] = 1. + (0.006944444444444444*(-169. + 742.*nu - 750.*nu2 + 156.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 0.5605554442947213*x2 + x3*(4.002558222882566 - 0.9000999000999002*el3));
  /** 6,4)  */
  rholm[17] = 1. + (0.011904761904761904*(-86. + 462.*nu - 581.*nu2 + 133.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + PMTERMS_eps*(- 0.7228451986855349*x2 + x3*(7.359388663371044 - 1.6001776001776002*el4));
  /** 6,5)  */
  rholm[18] = 1. + (0.006944444444444444*(-185. + 838.*nu - 910.*nu2 + 220.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 1.0973940686333457*x2 + x3*(11.623366217471297 - 2.5002775002775004*el5));
  /** (6,6) */
  rholm[19] = 1. + (0.011904761904761904*(-106. + 602.*nu - 861.*nu2 + 273.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + PMTERMS_eps*(- 1.5543111183867486*x2 + x3*(16.645950799433503 - 3.6003996003996006*el6));
  
  /** l=7 ------------------------------------------------------------------
   *  (7,1) */
  rholm[20] = 1. + (0.0014005602240896359*(-618. + 2518.*nu - 2083.*nu2 + 228.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*( - 0.1508235111143767*x2 + x3*(0.2581280702019663 - 0.07355557607658449*el1));
  /** (7,2) */
  rholm[21] = 1. + (0.00006669334400426837*(16832. - 123489.*nu + 273924.*nu2 - 190239.*nu3 + 32760.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 0.351319484450667*PMTERMS_eps*x2;
  /** (7,3) */
  rholm[22] = 1. + (0.0014005602240896359*(-666. + 2806.*nu - 2563.*nu2 + 420.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 0.37187416047628863*x2 + x3*(3.0835293524055283 - 0.6620001846892604*el3));
  /** (7,4) */
  rholm[23] = 1. + (0.00006669334400426837*(17756. - 131805.*nu + 298872.*nu2 - 217959.*nu3 + 41076.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 0.6473746896670599*PMTERMS_eps*x2;
  /** (7,5) */
  rholm[24] = 1. + (0.0014005602240896359*(-762. + 3382.*nu - 3523.*nu2 + 804.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 0.8269193364414116*x2 + x3*(8.750589067052443 - 1.838889401914612*el5));
  /** (7,6) */
  rholm[25] = 1. + (0.0006002400960384153*(2144. - 16185.*nu + 37828.*nu2 - 29351.*nu3 + 6104.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - 1.1403265020692532*PMTERMS_eps*x2;
  /** (7,7) */
  rholm[26] = 1. + (0.0014005602240896359*(-906. + 4246.*nu - 4963.*nu2 + 1380.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 1.5418467934923434*x2 + x3*(17.255875091408523 - 3.6042232277526396*el7));
  
  /** l=8 ------------------------------------------------------------------
   *  (8,1) */
  rholm[27] = 1. + (0.00005482456140350877*(20022. - 126451.*nu + 236922.*nu2 - 138430.*nu3 + 21640.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - 0.26842133517043704*PMTERMS_eps*x2;
  /** (8,2) */
  rholm[28] = 1. + (0.0003654970760233918*(2462. - 17598.*nu + 37119.*nu2 - 22845.*nu3 + 3063.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - 0.2261796441029474*PMTERMS_eps*x2;
  /** (8,3) */
  rholm[29] = 1. + (0.00005482456140350877*(20598. - 131059.*nu + 249018.*nu2 - 149950.*nu3 + 24520.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - 0.4196774909106648*PMTERMS_eps*x2;
  /** (8,4) */
  rholm[30] = 1. + (0.0003654970760233918*(2666. - 19434.*nu + 42627.*nu2 - 28965.*nu3 + 4899.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - 0.47652059150068155*PMTERMS_eps*x2;
  /** (8,5) */
  rholm[31] = 1. + (0.00027412280701754384*(4350. - 28055.*nu + 54642.*nu2 - 34598.*nu3 + 6056.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3)- 0.7220789990670207*PMTERMS_eps*x2;
  /** (8,6) */
  rholm[32] = 1. + (0.0010964912280701754*(1002. - 7498.*nu + 17269.*nu2 - 13055.*nu3 + 2653.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 0.9061610303170207*PMTERMS_eps*x2;
  /** (8,7) */
  rholm[33] = 1. + (0.00005482456140350877*(23478. - 154099.*nu + 309498.*nu2 - 207550.*nu3 + 38920.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3)- 1.175404252991305*PMTERMS_eps*x2;
  /** (8,8) */
  rholm[34] = 1. + (0.0003654970760233918*(3482. - 26778.*nu + 64659.*nu2 - 53445.*nu3 + 12243.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 1.5337092502821381*PMTERMS_eps*x2;
  
  /** Amplitudes */
  int k;
  for (k = 0; k < KMAX; k++) {
      flm[k] = gsl_pow_int(rholm[k], LINDEX[k]);
  }
  
}

/**
  * Function: eob_wav_flm_s_old
  * ---------------------------
  *   Computes the resummed amplitudes for the spin case. 
  *   This function computes the residual amplitude corrections flm's as 
  *   introduced in Damour, Iyer & Nagar, PRD 79, 064004 (2008).
  *   The orbital part is taken at the usual 3^{+2} PN order, i.e. 3PN terms
  *   are integrated by the 4PN and 5PN test-particle terms, with the higher
  *   modes obtained by Fujita & Iyer.
  *   Note that the variables called here (a1,a2)
  *   are what we usually cal tilde{a}_1 and tilde{a}_2 and are defined as
  *   a1 = X1*chi1, a2=X2*chi2 and are passed here as parameters. Special 
  *   combinations of these quantities are used here to write the spin-dependent
  *   part of the waveform in particularly compact form, so that the (spinning)
  *   test-particle limit is recovered just by visual inspection of the equations 
  * 
  *   @param[in]  x       : x = (M omega)^2
  *   @param[in]  nu      : symmetric mass ratio
  *   @param[in]  X1      : X1 = M1/M
  *   @param[in]  X2      : X2 = M2/M
  *   @param[in]  chi1    : chi1 = S1/M1^2
  *   @param[in]  chi2    : chi2 = S2/M2^2
  *   @param[in]  a1      : a1 = X1*chi1
  *   @param[in]  a2      : a2 = X2*chi2
  *   @param[in]  C_Q1    : spin-induced quadrupole of body 1
  *   @param[in]  C_Q2    : spin-induced quadrupole of body 2
  *   @param[in]  usetidal: flag for use tidal corrections
  *   @param[out] rholm   : resummed amplitudes
  *   @param[out] flm     : resummed amplitudes
*/
void eob_wav_flm_s_old(double x, double nu, double X1, double X2, double chi1, double chi2, double a1, double a2, double C_Q1, double C_Q2, int usetidal,
		       double *rholm, double *flm)
{

  /* declaring the spin-dependent terms */
  double rho22S;
  double rho32S;
  double rho44S;
  double rho42S;
  double f21S;
  double f33S;
  double f31S;
  double f43S;
  double f41S;
      
  const double a0      = a1+a2;
  const double a12     = a1-a2;
  const double X12     = X1-X2;
  const double a0X12   = a0*X12;
  const double a12X12  = a12*X12;
  
  /* Shorthands */
  const double nu2 = nu*nu;
  const double nu3 = nu*nu2;
  const double nu4 = nu*nu3;
  const double x2 = x*x;
  const double x3 = x*x2;
  const double x4 = x*x3;
  const double x5 = x*x4;
  const double v  = sqrt(x);
  const double v2 = x;
  const double v3 = v*v2;
  const double v4 = v3*v;
  const double v5 = v4*v;
  
  /* Compute EulerLogs */
  const double el1 = Eulerlog(x,1);
  const double el2 = Eulerlog(x,2);
  const double el3 = Eulerlog(x,3);
  const double el4 = Eulerlog(x,4);
  const double el5 = Eulerlog(x,5);
  const double el6 = Eulerlog(x,6);
  const double el7 = Eulerlog(x,7);
    
  /** l=m=2 multipole */
  /* spin-orbit */
  const double cSO_lo    = (-0.5*a0 - a12X12/6.);
  const double cSO_nlo   = (-52./63.-19./504.*nu)*a0 - (50./63.+209./504.*nu)*a12X12;
  
  /* SPIN-SPIN contribution */
  double cSS_lo;
  if (usetidal) {
#if (EXCLUDESPINSPINTIDES)
    /* Switch off spin-spin-tidal couplings */
    /* See also: eob_dyn_s_get_rc() */
    cSS_lo = 0.;
    /* Above code switch off everything, 
       Alt. one can set C_Q1=C_Q2=0, but keep the term: */
    /*
      cSS_lo = a1*a2;
    */

#else
    cSS_lo = 0.5*(C_Q1*a1*a1 + 2.*a1*a2 + C_Q2*a2*a2);
#endif
  } else {
    cSS_lo = 0.5*a0*a0; 
  }

  /* rho_22^S: Eq. (80) of Damour & Nagar, PRD 90, 044018 (2014) */
  rho22S = cSO_lo*v3 + cSS_lo*v4 + cSO_nlo*v5 ;
    
  /** l>=3, m=even: multipoles rewritten in compact and self-explanatory form */
  rho32S = (a0-a12X12)/(3.*(1.-3.*nu))*v;
  rho44S = (-19./30.*a0 -  (1.-21.*nu)/(30.-90.*nu)*a12X12)*v3;
  rho42S = ( -1./30.*a0 - (19.-39.*nu)/(30.-90.*nu)*a12X12)*v3;
  
  /** l>=2, m=odd: multipoles rewritten in compact and self-explanatory form */
  f21S = -1.5*a12*v + ((110./21. + 79./84.*nu)*a12 - 13./84.*a0X12)*v3;
  f33S = ((-0.25 + 2.5*nu)*a12 - 1.75*a0X12)*v3;
  f31S = ((-2.25 + 6.5*nu)*a12 + 0.25*a0X12)*v3;
  f43S = (( 5. -10.*nu)*a12 - 5.*a0X12)/(-4.+8.*nu)*v;
  f41S = f43S;
 
   
  // l=2 ------------------------------------------------------------------
  //(2,2)
  rholm[1] = 1. + (-1.0238095238095237 + 0.6547619047619048*nu)*x + (-1.94208238851096 - 1.5601379440665155*nu + 0.4625614134542706*nu2)*x2 + x3*(12.736034731834051 - 2.902228713904598*nu - 1.9301558466099282*nu2 + 0.2715020968103451*nu3 - 4.076190476190476*el2) + x4*(-2.4172313935587004 + 4.173242630385488*el2) + x5*(-30.14143102836864 + 7.916297736025627*el2);
    
    //(2,1)
    rholm[0] = 1. + (-1.0535714285714286 + 0.27380952380952384*nu)*x + (-0.8327841553287982 - 0.7789824263038548*nu + 0.13116496598639457*nu2)*x2 + x3*(2.9192806270460925 - 1.019047619047619*el1) + x4*(-1.28235780892213 + 1.073639455782313*el1) + PMTERMS_eps*x5*(-3.8466571723355227 + 0.8486467106683944*el1);
    
     
    // l=3 ------------------------------------------------------------------
    //(3,3)
    rholm[4] = 1. + (-1.1666666666666667 + 0.6666666666666666*nu)*x + (-1.6967171717171716 - 1.8797979797979798*nu + 0.45151515151515154*nu2)*x2 + x3*(14.10891386831863 - 3.7142857142857144*el3) + x4*(-6.723375314944128 + 4.333333333333333*el3) + PMTERMS_eps*x5*(-29.568699895427518 + 6.302092352092352*el3);
    
    //(3,2)
    rholm[3] = 1. + (0.003703703703703704*(328. - 1115.*nu + 320.*nu2)*x)/(-1. + 3.*nu) + (6.235191420376606e-7*(-1.444528e6 + 8.050045e6*nu - 4.725605e6*nu2 - 2.033896e7*nu3 + 3.08564e6*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + x3*(6.220997955214429 - 1.6507936507936507*el2) + PMTERMS_eps*x4*(-3.4527288879001268 + 2.005408583186361*el2);
    

    
    //(3,1)
    rholm[2] = 1. + (-0.7222222222222222 - 0.2222222222222222*nu)*x + (0.014169472502805836 - 0.9455667789001122*nu - 0.46520763187429853*nu2)*x2 + x3*(1.9098284139598072 - 0.4126984126984127*el1) + x4*(0.5368150316615179 + 0.2980599647266314*el1) + PMTERMS_eps*x5*(1.4497991763035063 - 0.0058477188106817735*el1);
    
    // l=4 ------------------------------------------------------------------
    //(4,4)
    rholm[8] = 1. + (0.0007575757575757576*(1614. - 5870.*nu + 2625.*nu2)*x)/(-1. + 3.*nu) + (3.1534122443213353e-9*(-5.11573572e8 + 2.338945704e9*nu - 3.13857376e8*nu2 - 6.733146e9*nu3 + 1.252563795e9*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + x3*(15.108111214795123 - 3.627128427128427*el4) + PMTERMS_eps*x4*(-8.857121657199649 + 4.434988849534304*el4);
    
    
    //(4,3)
    rholm[7] = 1. + (0.005681818181818182*(222. - 547.*nu + 160.*nu2)*x)/(-1. + 2.*nu) - 0.9783218202252293*x2 + PMTERMS_eps*(x3*(8.519456157072423 - 2.0402597402597404*el3) +      x4*(-5.353216984886716 + 2.5735094451003544*el3));
    
    //(4,2)
    rholm[6] = 1. + (0.0007575757575757576*(1146. - 3530.*nu + 285.*nu2)*x)/(-1. + 3.*nu) - (3.1534122443213353e-9*(1.14859044e8 - 2.95834536e8*nu - 1.204388696e9*nu2 + 3.04798116e9*nu3 + 3.79526805e8*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + 4.550378418934105e-12*x3*(8.48238724511e11 - 1.9927619712e11*el2) + PMTERMS_eps*x4*(-0.6621921297263365 + 0.787251738160829*el2);
    
    
    //(4,1)
    rholm[5] = 1. + (0.001893939393939394*(602. - 1385.*nu + 288.*nu2)*x)/(-1. + 2.*nu) - 0.36778992787515513*x2 + x3*(0.6981550175535535 - 0.2266955266955267*el1) + PMTERMS_eps*x4*(-0.7931524512893319 + 0.2584672482399755*el1);
    
    // l=5 ------------------------------------------------------------------
    //(5,5)
    rholm[13] = 1. + (0.002564102564102564*(487. - 1298.*nu + 512.*nu2)*x)/(-1. + 2.*nu) - 1.5749727622804546*x2 + PMTERMS_eps*(x3*(15.939827047208668 - 3.6037296037296036*el5) +      x4*(-10.272578060123237 + 4.500041838503377*el5));
    //(5,4)
    rholm[12] = 1. + (0.00007326007326007326*(-17448. + 96019.*nu - 127610.*nu2 + 33320.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + PMTERMS_eps*(- 1.0442142414362194*x2 +   x3*(10.252052781721588 - 2.3063869463869464*el4));
    //(5,3)
    rholm[11] = 1. + (0.002564102564102564*(375. - 850.*nu + 176.*nu2)*x)/(-1. + 2.*nu) - 0.5788010707241477*x2 + PMTERMS_eps*(x3*(5.733973288504755 - 1.2973426573426574*el3) +      x4*(-1.9573287625526001 + 1.2474448628294783*el3));
    //(5,2)
    rholm[10] = 1. + (0.00007326007326007326*(-15828. + 84679.*nu - 104930.*nu2 + 21980.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + PMTERMS_eps*(- 0.4629337197600934*x2 +      x3*(2.354458371550237 - 0.5765967365967366*el2));
    //(5,1)
    rholm[9] = 1. + (0.002564102564102564*(319. - 626.*nu + 8.*nu2)*x)/(-1. + 2.*nu) - 0.1047896120973044*x2 + PMTERMS_eps*(x3*(0.642701885362399 - 0.14414918414918415*el1) +      x4*(-0.07651588046467575 + 0.11790664036817883*el1));
    
    // l=6 ------------------------------------------------------------------
    //(6,6)
    rholm[19] = 1. + (0.011904761904761904*(-106. + 602.*nu - 861.*nu2 + 273.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + PMTERMS_eps*(- 1.5543111183867486*x2 + x3*(16.645950799433503 - 3.6003996003996006*el6));
    //(6,5)
    rholm[18] = 1. + (0.006944444444444444*(-185. + 838.*nu - 910.*nu2 + 220.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 1.0973940686333457*x2 + x3*(11.623366217471297 - 2.5002775002775004*el5));
    //(6,4)
    rholm[17] = 1. + (0.011904761904761904*(-86. + 462.*nu - 581.*nu2 + 133.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + PMTERMS_eps*(- 0.7228451986855349*x2 + x3*(7.359388663371044 - 1.6001776001776002*el4));
    //(6,3)
    rholm[16] = 1. + (0.006944444444444444*(-169. + 742.*nu - 750.*nu2 + 156.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 0.5605554442947213*x2 + x3*(4.002558222882566 - 0.9000999000999002*el3));
    //(6,2)
    rholm[15] = 1. + (0.011904761904761904*(-74. + 378.*nu - 413.*nu2 + 49.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + PMTERMS_eps*( - 0.24797525070634313*x2 + x3*(1.7942694138754138 - 0.40004440004440006*el2));
    //(6,1)
    rholm[14] = 1. + (0.006944444444444444*(-161. + 694.*nu - 670.*nu2 + 124.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 0.29175486850885135*x2 + x3*(0.21653486654395454 - 0.10001110001110002*el1));
    
    // l=7 ------------------------------------------------------------------
    //(7,7)
    rholm[26] = 1. + (0.0014005602240896359*(-906. + 4246.*nu - 4963.*nu2 + 1380.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 1.5418467934923434*x2 + x3*(17.255875091408523 - 3.6042232277526396*el7));
    //(7,6)
    rholm[25] = 1. + (0.0006002400960384153*(2144. - 16185.*nu + 37828.*nu2 - 29351.*nu3 + 6104.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - 1.1403265020692532*PMTERMS_eps*x2;
    //(7,5)
    rholm[24] = 1. + (0.0014005602240896359*(-762. + 3382.*nu - 3523.*nu2 + 804.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 0.8269193364414116*x2 + x3*(8.750589067052443 - 1.838889401914612*el5));
    //(7,4)
    rholm[23] = 1. + (0.00006669334400426837*(17756. - 131805.*nu + 298872.*nu2 - 217959.*nu3 + 41076.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 0.6473746896670599*PMTERMS_eps*x2;
    //(7,3)
    rholm[22] = 1. + (0.0014005602240896359*(-666. + 2806.*nu - 2563.*nu2 + 420.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 0.37187416047628863*x2 + x3*(3.0835293524055283 - 0.6620001846892604*el3));
    //(7,2)
    rholm[21] = 1. + (0.00006669334400426837*(16832. - 123489.*nu + 273924.*nu2 - 190239.*nu3 + 32760.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 0.351319484450667*PMTERMS_eps*x2;
    //(7,1)
    rholm[20] = 1. + (0.0014005602240896359*(-618. + 2518.*nu - 2083.*nu2 + 228.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*( - 0.1508235111143767*x2 + x3*(0.2581280702019663 - 0.07355557607658449*el1));
    
        
    // l=8 ------------------------------------------------------------------
    //(8,8)
    rholm[34] = 1. + (0.0003654970760233918*(3482. - 26778.*nu + 64659.*nu2 - 53445.*nu3 + 12243.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 1.5337092502821381*PMTERMS_eps*x2;
    //(8,7)
    rholm[33] = 1. + (0.00005482456140350877*(23478. - 154099.*nu + 309498.*nu2 - 207550.*nu3 + 38920.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3)- 1.175404252991305*PMTERMS_eps*x2;
    //(8,6)
    rholm[32] = 1. + (0.0010964912280701754*(1002. - 7498.*nu + 17269.*nu2 - 13055.*nu3 + 2653.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 0.9061610303170207*PMTERMS_eps*x2;
    //(8,5)
    rholm[31] = 1. + (0.00027412280701754384*(4350. - 28055.*nu + 54642.*nu2 - 34598.*nu3 + 6056.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3)- 0.7220789990670207*PMTERMS_eps*x2;
    //(8,4)
    rholm[30] = 1. + (0.0003654970760233918*(2666. - 19434.*nu + 42627.*nu2 - 28965.*nu3 + 4899.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - 0.47652059150068155*PMTERMS_eps*x2;
    //(8,3)
    rholm[29] = 1. + (0.00005482456140350877*(20598. - 131059.*nu + 249018.*nu2 - 149950.*nu3 + 24520.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - 0.4196774909106648*PMTERMS_eps*x2;
    //(8,2)
    rholm[28] = 1. + (0.0003654970760233918*(2462. - 17598.*nu + 37119.*nu2 - 22845.*nu3 + 3063.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - 0.2261796441029474*PMTERMS_eps*x2;
    //(8,1)
    rholm[27] = 1. + (0.00005482456140350877*(20022. - 126451.*nu + 236922.*nu2 - 138430.*nu3 + 21640.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - 0.26842133517043704*PMTERMS_eps*x2;
    
    /** Amplitudes */
    flm[0] = gsl_pow_int(rholm[0], 2);
    flm[0] = (X12*flm[0] + f21S);
    
    flm[1] = gsl_pow_int(rholm[1]+ rho22S, 2);

    flm[2] = gsl_pow_int(rholm[2], 3);
    flm[2] = (X12*flm[2] + f31S);

    flm[3] = gsl_pow_int(rholm[3]+ rho32S, 3);

    flm[4] = gsl_pow_int(rholm[4], 3);
    flm[4] = (X12*flm[4] + f33S);
    
    flm[5] = gsl_pow_int(rholm[5], 4);
    flm[5] = (X12*flm[5] + f41S);
    
    flm[6]=gsl_pow_int(rholm[6] + rho42S, 4);

    flm[7] = gsl_pow_int(rholm[7], 4);
    flm[7] = (X12*flm[7] + f43S);

    flm[8] = gsl_pow_int(rholm[8] + rho44S, 4);

    /** Amplitudes */
    int k;
    for (k = 9; k < KMAX; k++) {
      flm[k] = gsl_pow_int(rholm[k], LINDEX[k]);
    }
    
}


/**
  * Function: eob_wav_hlmTidal_fmode_fact22A
  * ----------------------------------------
  *   Computes the f-mode resonance correction for 22 amplitude of star A
  * 
  *   @param[in]  x       : x = (M omega)^2
  *   @param[in]  alpha   : alpha = M omega_22
  *   @param[in]  bomgf   : omega_22
  *   @param[in]  XB      : X1-X2
  *   
  *   @return f-mode resonance correction for 22 amplitude of star A
*/
double eob_wav_hlmTidal_fmode_fact22A(double x, double alpha, double bomgf, double XB)
{
  const double Omega2 = gsl_pow_int(x,3);
  const double y2 = SQ(bomgf)/Omega2;
  return ( y2 * (-1. + alpha*( 1 + 6.*XB/y2 ) )/( 3.*(1. + 2*XB) ) ); 
  /* return ( (-1. + alpha)*SQ(bomgf) + 6.*alpha*XB*Omega2 )/( (1. + 2*XB)*(3.*Omega2) ); */
}


/**
 * Function: prolong_euler_angles_TD
 * ---------------------------------
 *   Time domain routine to interpolate and prolong the Euler angles computed from the
 *   dynamics beyond merger.
 *   The prolongation is done by:
 *   (1) identifying the end of the dynamics (merger time) by finding the maximum 
         of the co-precessing A_{22};
 *   (2) interpolating the Euler angles on the waveform time grid via spline;
 *   (3) fixing the values to the last value beyond merger
 *  
 *   @param[out]  alpha   : alpha euler angle
 *   @param[out]  beta    : beta euler angle
 *   @param[out]  gamma   : gamma euler angle
 *   @param[in]   dyn     : EOB dynamics
 *   @param[in]   spin    : spin dynamics
 *   @param[in]   hlm     : multipolar waveform
 *
*/
void prolong_euler_angles_TD(double *alpha, double *beta, double *gamma, Dynamics *dyn, DynamicsSpin *spin, Waveform_lm *hlm){

  /* First, unwrap alpha and gamma */
  unwrap_euler(spin->data[EOB_EVOLVE_SPIN_alp], spin->size);
  unwrap_euler(spin->data[EOB_EVOLVE_SPIN_gam], spin->size);

  /* Determine merger time (as maximum of the co-precessing A_{22}) */
  int jmax = 0;
  for(int i=1; i<hlm->size; i++){
    if(hlm->ampli[1][i]>hlm->ampli[1][jmax]) jmax = i;
  }

  double tmax_wav = hlm->time[jmax];

  /* Find the corresponding time in the spin dynamics */
  const int tmax_dyn_idx = find_point_bisection(tmax_wav, spin->size, spin->time, 1);
  const int tmax_wav_idx = find_point_bisection(spin->time[tmax_dyn_idx], hlm->size, hlm->time, 1);

  /* Interpolation */
  interp_spline_omp(spin->time, spin->data[EOB_EVOLVE_SPIN_alp], spin->size, hlm->time, tmax_wav_idx, alpha);
  interp_spline_omp(spin->time, spin->data[EOB_EVOLVE_SPIN_bet], spin->size, hlm->time, tmax_wav_idx, beta);  
  interp_spline_omp(spin->time, spin->data[EOB_EVOLVE_SPIN_gam], spin->size, hlm->time, tmax_wav_idx, gamma); 

  /* Now, prolong the angles based on user request 
     for t > tM_idx, fix the values to the last  */
  for(int j=tmax_wav_idx; j < hlm->size; j++){
    alpha[j] = alpha[tmax_wav_idx-1];
    beta[j]  = beta[tmax_wav_idx-1];
    gamma[j] = gamma[tmax_wav_idx-1];
  }
}


/**
 * Function: prolong_euler_angles
 * ------------------------------
 *   Time domain routine to prolong the Euler angles computed from the 
 *   dynamics beyond merger.  
 *   The prolongation is done by (i) interpolating the Euler angles
 *   on the EOB Orbital frequency grid via spline, (ii) evaluating the spline
 *   at 1/2 of the frequency of the co-precessing 22 mode up to the peak of omg22/2, 
 *   (iii) estimating their values beyond merger
 * 
 *   @param[out]  alpha   : alpha euler angle
 *   @param[out]  beta    : beta euler angle
 *   @param[out]  gamma   : gamma euler angle
 *   @param[in]   dyn     : EOB dynamics
 *   @param[in]   spin    : spin dynamics
 *   @param[in]   hlm     : multipolar waveform
 * 
*/
void prolong_euler_angles(double *alpha, double *beta, double *gamma, Dynamics *dyn, DynamicsSpin *spin, Waveform_lm *hlm){
  
  /* choose whether to use MOmega (from the dynamics) or MOmega_22 for the interpolation */
  int map_from_22 = 1;
  spin->data[EOB_EVOLVE_SPIN_alp][0] = spin->data[EOB_EVOLVE_SPIN_alp][1]; //alpha_initial_condition(EOBPars);
  spin->data[EOB_EVOLVE_SPIN_gam][0] = spin->data[EOB_EVOLVE_SPIN_alp][0];
  /* First, unwrap alpha and gamma */
  unwrap_euler(spin->data[EOB_EVOLVE_SPIN_alp], spin->size);
  unwrap_euler(spin->data[EOB_EVOLVE_SPIN_gam], spin->size);
  double *omega;
  int size_omega;
  if(map_from_22){
    double *omg22_eob;
    omg22_eob = malloc ( hlm->size * sizeof(double) );
    //D0_x_4(hlm->phase[1], hlm->time, hlm->size, omg22_eob);
    D0(hlm->phase[1], hlm->time[1]-hlm->time[0], hlm->size, omg22_eob);
    for(int i =0; i < hlm->size; i++) omg22_eob[i] = omg22_eob[i]/2;
    omega = omg22_eob;
    size_omega = hlm->size;
  } else {
    omega = dyn->data[EOB_MOMG];
    size_omega = dyn->size;
  }
  

  /* Then, do the interpolations.
  find the max of omega*/
  int omg_jmax =0;
  for(int i=0; i < size_omega; i++) {
    omg_jmax = i;
    if(i >5 && omega[i+1] <= omega[i])
      break;
    if(omega[i+1] > spin->data[EOB_EVOLVE_SPIN_Momg][spin->size -1])
      break;
  }

  /* allocate temporary angles*/
  int size_tmp = omg_jmax +1;
  double *alpha_tmp, *beta_tmp, *gamma_tmp;
  int tM_idx;
  alpha_tmp = malloc ( size_tmp * sizeof(double) );
  beta_tmp  = malloc ( size_tmp * sizeof(double) );
  gamma_tmp = malloc ( size_tmp * sizeof(double) );


  /* first interpolation: angles(omega_PN)->angles(omega_EOB) */
  if (VERBOSE) printf("First interpolation: angles(omega_PN)->angles(omega_EOB)\n");
  interp_spline_omp(spin->data[EOB_EVOLVE_SPIN_Momg], spin->data[EOB_EVOLVE_SPIN_alp], spin->size, omega, size_tmp, alpha_tmp);
  interp_spline_omp(spin->data[EOB_EVOLVE_SPIN_Momg], spin->data[EOB_EVOLVE_SPIN_bet], spin->size, omega, size_tmp, beta_tmp);
  interp_spline_omp(spin->data[EOB_EVOLVE_SPIN_Momg], spin->data[EOB_EVOLVE_SPIN_gam], spin->size, omega, size_tmp, gamma_tmp);

  if (VERBOSE) printf("Second interpolation: omega_EOB->t_EOB up to omega peak\n");
  if(!map_from_22){
    /* map omega_dyn_EOB->t_EOB up to the peak of the EOB dynamics */
    tM_idx = find_point_bisection(dyn->time[omg_jmax], hlm->size, hlm->time, 1);
    interp_spline_omp(dyn->time, alpha_tmp, size_tmp, hlm->time, tM_idx+1, alpha);
    interp_spline_omp(dyn->time, beta_tmp,  size_tmp, hlm->time, tM_idx+1, beta);  
    interp_spline_omp(dyn->time, gamma_tmp, size_tmp, hlm->time, tM_idx+1, gamma); 
  } else {
    /* map omega22_EOB->t_EOB up to the peak of the EOB dynamics */
    tM_idx = find_point_bisection(hlm->time[omg_jmax], hlm->size, hlm->time, 1);
    interp_spline_omp(hlm->time, alpha_tmp, size_tmp, hlm->time, tM_idx+1, alpha);
    interp_spline_omp(hlm->time, beta_tmp,  size_tmp, hlm->time, tM_idx+1, beta);  
    interp_spline_omp(hlm->time, gamma_tmp, size_tmp, hlm->time, tM_idx+1, gamma); 
  }

  if (VERBOSE) printf("Correct for backward/forward integration\n");
  for(int i =0; i<tM_idx+1;i++){
    if(omega[i] <  spin->omg_backward){
      alpha[i] =  alpha[i] - Pi;
      beta[i]  = -beta[i]; 
      gamma[i] =  gamma[i] - Pi;
    }
  }
  /* Now, prolong the angles based on user request */
  if (EOBPars->ringdown_eulerangles == RD_EULERANGLES_CONSTANT) {
    //for t > tM_idx, fix the values to the last
    for(int j=tM_idx+1; j < hlm->size; j++){
      alpha[j] = alpha[tM_idx];
      beta[j]  = beta[tM_idx];
      gamma[j] = gamma[tM_idx];
    }
  } else if (EOBPars->ringdown_eulerangles == RD_EULERANGLES_QNMs){
    /* use QNM for alpha_dot, and fix beta constant */
    /* Table VIII or https://arxiv.org/pdf/gr-qc/0512160.pdf */
  
    //final spin (assume merger ~ max omega)
    double SAmrg[3], SBmrg[3], Lmrg[3], Jmrg[3];
    //FIXME: use the max of omega!
    omg_jmax = spin->size-1;

    SAmrg[0] = spin->data[EOB_EVOLVE_SPIN_SxA][omg_jmax];
    SAmrg[1] = spin->data[EOB_EVOLVE_SPIN_SyA][omg_jmax];
    SAmrg[2] = spin->data[EOB_EVOLVE_SPIN_SzA][omg_jmax];

    SBmrg[0] = spin->data[EOB_EVOLVE_SPIN_SxB][omg_jmax];
    SBmrg[1] = spin->data[EOB_EVOLVE_SPIN_SyB][omg_jmax];
    SBmrg[2] = spin->data[EOB_EVOLVE_SPIN_SzB][omg_jmax];

    //final L
    double nu    = EOBPars->nu;
    double nu2   = nu*nu;
    double v2mrg = pow(spin->data[EOB_EVOLVE_SPIN_Momg][omg_jmax], 0.6666666666666);
    double v4mrg = v2mrg*v2mrg;
    double vmrg  = sqrt(v2mrg);
    const double L2PN = nu/vmrg*(1 + v2mrg*(1.5+0.1666666666666667*nu) + v4mrg*(3.375 - 2.375*nu + 0.04166666666666666*nu2));
    Lmrg[0]  = L2PN*spin->data[EOB_EVOLVE_SPIN_Lx][omg_jmax];
    Lmrg[1]  = L2PN*spin->data[EOB_EVOLVE_SPIN_Ly][omg_jmax];
    Lmrg[2]  = L2PN*spin->data[EOB_EVOLVE_SPIN_Lz][omg_jmax];

    for(int i=0; i<3;i++)
      Jmrg[i] = SAmrg[i]+SBmrg[i]+Lmrg[i];

    double adot, JdotL;
    vect_dot3(Jmrg, Lmrg, &JdotL);

    if(JdotL>0){
      /** (l,m,n)=(2,2,0) */
      double f10 = 1.5251; 
      double f20 = -1.1568;
      double f30 = 0.1292;
      double omega220  = (f10 + f20*pow(1. - EOBPars->abhf, f30));  
      /** (l,m,n)=(2,1,0) */
      f10 = 0.6; 
      f20 = -0.2339;
      f30 = 0.4175;
      double omega210  = (f10 + f20*pow(1. - EOBPars->abhf, f30));  
      adot = omega220-omega210;
    } else {
      /** (l,m,n)=(2,-2,0) */
      double f10 = 0.2938; 
      double f20 = 0.0782;
      double f30 = 1.3546;
      double omega2m20  = (f10 + f20*pow(1. - EOBPars->abhf, f30));  
      /** (l,m,n)=(2,-1,0) */
      f10 = 0.3441; 
      f20 = 0.0293;
      f30 = 2.0010;
      double omega2m10  = (f10 + f20*pow(1. - EOBPars->abhf, f30)); 
      adot = omega2m10 - omega2m20;
    }

    for(int j=tM_idx+1; j < hlm->size; j++){
      double dt= hlm->time[j]-hlm->time[tM_idx];
      beta[j]  = beta[tM_idx];
      alpha[j] = alpha[tM_idx] + dt*adot;
      gamma[j] = gamma[tM_idx] + dt*adot*cos(beta[j]);
    }
  } else {
    errorexit("Need to specify angles for ringdown!\n");
  }

  /* free */
  free(alpha_tmp);
  free(beta_tmp);
  free(gamma_tmp);
  if (map_from_22) free(omega);
}

/**
  * Function: prolong_euler_angles_FD
  * ---------------------------------
  *   Frequency domain routine to prolong the Euler angles computed from the 
  *   dynamics until the end of the waveform.
  *   The prolongation is done by assuming that angles remain constant between
  *   the final frequency of the dynamics and the final frequency of the waveform.
  * 
  *   @param[in]   spin    : spin dynamics
  *   @param[in]   hlm     : multipolar waveform (frequency domain)
  *
*/
void prolong_euler_angles_FD(DynamicsSpin *spin, WaveformFD_lm *hlm)
{
  int spinsize = spin->size;
  /* First, unwrap alpha and gamma */
  spin->data[EOB_EVOLVE_SPIN_alp][0] = spin->data[EOB_EVOLVE_SPIN_alp][1];
  spin->data[EOB_EVOLVE_SPIN_gam][0] = spin->data[EOB_EVOLVE_SPIN_gam][1];
  unwrap_HM(spin->data[EOB_EVOLVE_SPIN_alp], spinsize);
  unwrap_HM(spin->data[EOB_EVOLVE_SPIN_gam], spinsize);
  
  /*prolong up to final hlm frequency (constant angles by default), add just N=2 points*/

  double omg_fin = hlm->freq[hlm->size-1];
  int N = 2;
  double domg = (omg_fin - spin->data[EOB_EVOLVE_SPIN_Momg][spinsize-1])/N;
  if (spin->data[EOB_EVOLVE_SPIN_Momg][spinsize-1] < omg_fin){
    DynamicsSpin_push(&spin, spinsize + N);

    /* fill the new points with constant angles */
    for(int i = 0; i<N;i++){
      spin->data[EOB_EVOLVE_SPIN_Momg][spinsize+i] = spin->data[EOB_EVOLVE_SPIN_Momg][spinsize-1] + (1+i)*domg;
      spin->data[EOB_EVOLVE_SPIN_alp][spinsize+i]  = spin->data[EOB_EVOLVE_SPIN_alp][spinsize-1];
      spin->data[EOB_EVOLVE_SPIN_bet][spinsize+i]  = spin->data[EOB_EVOLVE_SPIN_bet][spinsize-1];
      spin->data[EOB_EVOLVE_SPIN_gam][spinsize+i]  = spin->data[EOB_EVOLVE_SPIN_gam][spinsize-1];
    }
  }

}


/** Twist TD multipoles */
/**
  * Function: twist_hlm_TD
  * ----------------------
  *   Time domain routine to compute the twisted inertial multipoles
  *   from the co-precessing multipoles.
  * 
  *   @param[in]   dyn     : EOB dynamics
  *   @param[in]   hlm     : multipolar waveform
  *   @param[in]   spin    : spin dynamics
  *   @param[out]  hTlm    : twisted inertial multipoles
  *   @param[out]  hTlm_neg: twisted inertial multipoles (m<0)
  *   @param[out]  hTl0    : twisted inertial multipoles (m=0)
  *
*/
void twist_hlm_TD(Dynamics *dyn, Waveform_lm *hlm, DynamicsSpin *spin, int interp_spin_abc,
		  Waveform_lm *hTlm, Waveform_lm *hTlm_neg, Waveform_lm *hTl0)
{  

  int activemode[KMAX], activemode_inertial[KMAX];
  set_multipolar_idx_mask (activemode_inertial, KMAX, EOBPars->use_mode_lm_inertial, EOBPars->use_mode_lm_inertial_size, 1);
  set_multipolar_idx_mask (activemode,          KMAX, EOBPars->use_mode_lm,          EOBPars->use_mode_lm_size, 1);

  const int size  = hlm->size;
  double *alpha, *beta, *gamma;
  alpha = malloc ( size * sizeof(double) );
  beta  = malloc ( size * sizeof(double) );
  gamma = malloc ( size * sizeof(double) );
  
  /* Euler angles */
  if (EOBPars->model == MODEL_GIOTTO)
    prolong_euler_angles(alpha, beta, gamma, dyn, spin, hlm);
  else
    prolong_euler_angles_TD(alpha, beta, gamma, dyn, spin, hlm);
  
#if (DEBUG)    
    /*output angles */
    char fname[STRLEN*2];
    if (EOBPars->output_hpc + EOBPars->output_multipoles + EOBPars->output_dynamics){
      sprintf(fname,"%s/anglesint.txt",EOBPars->output_dir);
      FILE* fp;
      if ((fp = fopen(fname, "w+")) == NULL) errorexits("error opening file",fname);
      for (int i = 0; i < hlm->size; i++) {
        fprintf(fp, "%.9e  %.16e  %.16e  %.16e\n", hlm->time[i], alpha[i], beta[i], gamma[i]);
      }
      fclose(fp);
    }
#endif

  int emm_sign[3] = {1,-1,0}; //to loop over m>0, m<0, m=0
  int zero_flag = 1;
  /* Loop over modes */
  for (int k = 0; k < KMAX; k++ ) {
    /* compute twisted inertial modes */
    if (!activemode_inertial[k]) continue;  
    
    int ell = LINDEX[k];

    for (int q=0; q <=2; q++){
      int emm = MINDEX[k]*emm_sign[q];

      /* if you already computed the m=0 for this l, skip */
      if(q==2 && !zero_flag) continue;
      double eps = pow(-1., ell);
      
      // for each time ...
      for (int i = 0; i < size; i++) {
        // ... do the twist (sum up on m')
        double sumr = 0;
        double sumi = 0;
        for (int n = 1; n <= ell; n++) {
        // Uncomment below for previous summation
        // for (int n = -ell; n <= ell; n++) {
        //   if (n==0) continue; // skip m=0 modes
          int j = KINDEX[ell][abs(n)-1]; // map to linear index (ell,n) -> j
          if (!activemode[j]) continue;  // skip if co-precessing mode is not available
          double cosng = cos( n * gamma[i] );
          double sinng = sin( n * gamma[i] );
    
          // d^l_{m,s}(angle)
          //Checked by SA on April 2021
          double dl_mn = wigner_d_function(ell, n,emm, -beta[i]);
          double dl_mnn= wigner_d_function(ell,-n,emm, -beta[i]);

          // hlm modes are given as phase/amplitude
          // but here we need real/imag
          double real=0, imag=0;
          rmap_twist(&real,&imag, &(hlm->phase[j][i]), &(hlm->ampli[j][i]), 0);

          // Here we need to deal with m<0 modes
          // H_{l-m} = (-)^l H^{*}_{lm}
          double hln_real_n, hln_imag_n;
          hln_real_n =   eps * real;
          hln_imag_n = - eps * imag;	  
          double hln_real_p, hln_imag_p;
          hln_real_p = real;
          hln_imag_p = imag;

          sumr += dl_mn*(cosng * hln_real_p - sinng * hln_imag_p) + dl_mnn*(cosng * hln_real_n + sinng * hln_imag_n);
          sumi += dl_mn*(sinng * hln_real_p + cosng * hln_imag_p) + dl_mnn*(-sinng* hln_real_n + cosng * hln_imag_n);

          // Previous summation, less efficient
          // double hln_real, hln_imag;
          // if (n<0) {
          //   hln_real =   eps * real;
          //   hln_imag = - eps * imag;	  
          // } else {
          //   hln_real = real;
          //   hln_imag = imag;
          // }
          // sumr += dl_mn*(cosng * hln_real - sinng * hln_imag);
          // sumi += dl_mn*(sinng * hln_real + cosng * hln_imag);
    
        } // n (m')

        double hTlm_real =   sumr * cos( emm * alpha[i] ) + sumi * sin( emm * alpha[i] );
        double hTlm_imag = - sumr * sin( emm * alpha[i] ) + sumi * cos( emm * alpha[i] );
        
        // Re-map back into phase/ampli
        if(q==0)  // m>0 
          rmap_twist(&hTlm_real,&hTlm_imag, &(hTlm->phase[k][i]), &(hTlm->ampli[k][i]), 1);
        if(q==1)  // m<0
          rmap_twist(&hTlm_real,&hTlm_imag, &(hTlm_neg->phase[k][i]), &(hTlm_neg->ampli[k][i]), 1);
        if(q==2){ // m=0
          rmap_twist(&hTlm_real,&hTlm_imag, &(hTl0->phase[k][i]), &(hTl0->ampli[k][i]), 1);
          zero_flag = 0; //avoid re-computation of m=0 for fixed l
        }
      }// i (times)
    } // q (m>0,<0,=0)
    if(LINDEX[k+1]-LINDEX[k]) zero_flag = 1;
  }// k (active modes)
  
  hTlm->size = size;
  for(int i = 0; i < size; i++){
    hTlm->time[i]     = hlm->time[i];
    hTlm_neg->time[i] = hlm->time[i];
    hTl0->time[i]     = hlm->time[i];
  }

  free(alpha);
  free(beta);
  free(gamma);
  
}

/**
  * Function: compute_hpc_old
  * -------------------------
  *   Time domain routine to compute the (h+, hx) polarizations
  *   from the inertial multipoles.
  *   Note: deprecated
  * 
  *   @param[in]   hlm     : multipolar waveform
  *   @param[in]   nu      : symmetric mass ratio
  *   @param[in]   M       : total mass
  *   @param[in]   distance: distance to source
  *   @param[in]   amplitude_prefactor: amplitude prefactor
  *   @param[in]   phi     : reference phase
  *   @param[in]   iota    : inclination angle
  *   @param[out]  hpc     : (h+, hx) polarizations
  *
*/
void compute_hpc_old(Waveform_lm *hlm, double nu, double M, double distance, double amplitude_prefactor, double phi, double iota, Waveform *hpc)
{  
#ifdef _OPENMP
  if (USETIMERS) openmp_timer_start("compute_hpc");
#endif
#pragma omp parallel 
  {
    double Y_real[KMAX], Y_imag[KMAX];
    static const int mneg = 1; /* m>0 modes only, add m<0 modes afterwards */
    double Y_real_mneg[KMAX], Y_imag_mneg[KMAX];
    double Aki, cosPhi, sinPhi;
    double sumr, sumi;
    int activemode[KMAX];
    double Msun = M;
    set_multipolar_idx_mask (activemode, KMAX, EOBPars->use_mode_lm, EOBPars->use_mode_lm_size, 1);
    if (!(EOBPars->use_geometric_units)) Msun = M/MSUN_S;
#if (DEBUG)
    printf("h+,x: nu = %e M = %e D = %e Mpc phi = %e iota = %e prefactor = %e\n",
	   nu,Msun,distance,phi,iota,amplitude_prefactor);
#endif
    
    /* Precompute Ylm */
    for (int k = 0; k < KMAX; k++ ) {
      if (!activemode[k]) continue;
      spinsphericalharm(&Y_real[k], &Y_imag[k], -2, LINDEX[k], MINDEX[k], phi, iota);
      /* add m<0 modes */
      if ( (mneg) && (MINDEX[k]!=0) )
	spinsphericalharm(&Y_real_mneg[k], &Y_imag_mneg[k], -2, LINDEX[k], -MINDEX[k], phi, iota); 
      }
      
    /* Sum up  hlm * Ylm 
     * Note because EOB code defines phase>0, 
     * but the convention is h_lm = A_lm Exp[-I phi_lm] we have
     * h_{l,m>=0} = A_lm ( cos(phi_lm) - I*sin(phi_lm) ) for m>0 and
     * h_{l,m<0}  = (-)^l A_l|m| ( cos(phi_l|m|) + I*sin(phi_l|m|) ) for m<0 below
     * We now agree with, e.g., LALSimSphHarmMode.c: 64-74
     */
#pragma omp for
    for (int i = 0; i < hlm->size; i++) {
        hpc->time[i] = hlm->time[i]*M; 
        sumr = sumi = 0.;
	
        /* Loop over modes */
        for (int k = 0; k < KMAX; k++ ) {
	  if (!activemode[k]) continue;
	  Aki  = amplitude_prefactor * hlm->ampli[k][i];
	  cosPhi = cos( hlm->phase[k][i] );
	  sinPhi = sin( hlm->phase[k][i] );
	  sumr += Aki*(cosPhi*Y_real[k] + sinPhi*Y_imag[k]);
	  sumi += Aki*(cosPhi*Y_imag[k] - sinPhi*Y_real[k]); 
          
	  /* add m<0 modes */
	  if ( (mneg) && (MINDEX[k]!=0) ) { 
	    /* H_{l-m} = (-)^l H^{*}_{lm} */
	    if (LINDEX[k] % 2) {
	      sumr -= Aki*(cosPhi*Y_real_mneg[k] - sinPhi*Y_imag_mneg[k]); 
	      sumi -= Aki*(cosPhi*Y_imag_mneg[k] + sinPhi*Y_real_mneg[k]); 
	    }
	    else { 
	      sumr += Aki*(cosPhi*Y_real_mneg[k] - sinPhi*Y_imag_mneg[k]);
	      sumi += Aki*(cosPhi*Y_imag_mneg[k] + sinPhi*Y_real_mneg[k]); 
	    }
	  }    
        }
        /* h = h+ - i hx */
        hpc->real[i] = sumr;
        hpc->imag[i] = -sumi;
    }
  }
#ifdef _OPENMP
  if (USETIMERS) openmp_timer_stop("compute_hpc");
#endif
}

/**
  * Function: compute_hpc
  * ---------------------
  *   Time domain routine to compute the (h+, hx) polarizations
  *   from the inertial multipoles.
  * 
  *   @param[in]   hlm     : multipolar waveform
  *   @param[in]   hlm_neg : multipolar waveform (m<0)
  *   @param[in]   hl0     : multipolar waveform (m=0)
  *   @param[in]   nu      : symmetric mass ratio
  *   @param[in]   M       : total mass
  *   @param[in]   distance: distance to source
  *   @param[in]   amplitude_prefactor: amplitude prefactor
  *   @param[in]   phi     : reference phase
  *   @param[in]   iota    : inclination angle
  *   @param[out]  hpc     : (h+, hx) polarizations
  *
*/
void compute_hpc(Waveform_lm *hlm, Waveform_lm *hlm_neg, Waveform_lm *hl0, double nu, double M, double distance, double amplitude_prefactor, double phi, double iota, Waveform *hpc)
{  
#ifdef _OPENMP
  if (USETIMERS) openmp_timer_start("compute_hpc");
#endif
#pragma omp parallel 
  {
    double Y_real[KMAX], Y_imag[KMAX];
    /* m<0 */
    int mneg = 0; 
    if (hlm_neg == NULL) mneg =1;
    double Y_real_mneg[KMAX], Y_imag_mneg[KMAX];
    /* m=0*/
    int mzero     = 1;
    int zero_flag = 1;
    if (hl0 == NULL){
      mzero =0;
      zero_flag = 0;
    }
    double Y_real_m0[KMAX], Y_imag_m0[KMAX];

    double Aki, cosPhi, sinPhi;
    double sumr, sumi;
    int activemode[KMAX];
    double Msun = M;
    set_multipolar_idx_mask (activemode, KMAX, EOBPars->use_mode_lm_inertial, EOBPars->use_mode_lm_inertial_size, 1);
    if (!(EOBPars->use_geometric_units)) Msun = M/MSUN_S;
#if (DEBUG)
    printf("h+,x: nu = %e M = %e D = %e Mpc phi = %e iota = %e prefactor = %e\n",
	   nu,Msun,distance,phi,iota,amplitude_prefactor);
#endif
    
    /* Precompute Ylm */
    for (int k = 0; k < KMAX; k++ ) {
      if (!activemode[k]) continue;
      spinsphericalharm(&Y_real[k], &Y_imag[k], -2, LINDEX[k], MINDEX[k], phi, iota);
      /* add m<0 modes */
      if ( (MINDEX[k]!=0) )
	      spinsphericalharm(&Y_real_mneg[k], &Y_imag_mneg[k], -2, LINDEX[k], -MINDEX[k], phi, iota);
      if ( (mzero !=0) )
  	    spinsphericalharm(&Y_real_m0[k], &Y_imag_m0[k], -2, LINDEX[k], 0, phi, iota);
      }
      
    /* Sum up  hlm * Ylm 
     * Note because EOB code defines phase>0, 
     * but the convention is h_lm = A_lm Exp[-I phi_lm] we have
     * h_{l,m>=0} = A_lm ( cos(phi_lm) - I*sin(phi_lm) ) for m>0 and
     * h_{l,m<0}  = (-)^l A_l|m| ( cos(phi_l|m|) + I*sin(phi_l|m|) ) for m<0 below
     * We now agree with, e.g., LALSimSphHarmMode.c: 64-74
     */
#pragma omp for
    for (int i = 0; i < hlm->size; i++) {
        hpc->time[i] = hlm->time[i]*M; 
        sumr = sumi = 0.;
	
        /* Loop over modes */
        for (int k = 0; k < KMAX; k++ ) {
	        if (!activemode[k]) continue;
          Aki  = amplitude_prefactor * hlm->ampli[k][i];
          cosPhi = cos( hlm->phase[k][i] );
          sinPhi = sin( hlm->phase[k][i] );
          sumr += Aki*(cosPhi*Y_real[k] + sinPhi*Y_imag[k]);
          sumi += Aki*(cosPhi*Y_imag[k] - sinPhi*Y_real[k]); 
          
          /* precessing wf, there is no symmetry between +m and -m*/
          if (EOBPars->use_spins==MODE_SPINS_GENERIC && !mneg){
            if (!activemode[k]) continue;
            Aki  = amplitude_prefactor * hlm_neg->ampli[k][i];
            cosPhi = cos( hlm_neg->phase[k][i] );
            sinPhi = sin( hlm_neg->phase[k][i] );
            sumr += Aki*(cosPhi*Y_real_mneg[k] + sinPhi*Y_imag_mneg[k]);
            sumi += Aki*(cosPhi*Y_imag_mneg[k] - sinPhi*Y_real_mneg[k]); 
          }

          /* add m<0 modes */
          if ( (mneg) && (MINDEX[k]!=0) ) {
            /* H_{l-m} = (-)^l H^{*}_{lm} */
            if (LINDEX[k] % 2) {
              sumr -= Aki*(cosPhi*Y_real_mneg[k] - sinPhi*Y_imag_mneg[k]); 
              sumi -= Aki*(cosPhi*Y_imag_mneg[k] + sinPhi*Y_real_mneg[k]); 
            }
            else { 
              sumr += Aki*(cosPhi*Y_real_mneg[k] - sinPhi*Y_imag_mneg[k]);
              sumi += Aki*(cosPhi*Y_imag_mneg[k] + sinPhi*Y_real_mneg[k]); 
            }
          }   

          /* add m=0 modes */
          if (EOBPars->use_spins==MODE_SPINS_GENERIC){
            if ( (mzero) && zero_flag){
              Aki  = amplitude_prefactor * hl0->ampli[k][i];
              cosPhi = cos( hl0->phase[k][i] );
              sinPhi = sin( hl0->phase[k][i] );
              sumr += Aki*(cosPhi*Y_real_m0[k] + sinPhi*Y_imag_m0[k]); 
              sumi += Aki*(cosPhi*Y_imag_m0[k] - sinPhi*Y_real_m0[k]);
              zero_flag = 0;
            }
            if(LINDEX[k+1]-LINDEX[k]) zero_flag = 1;
          }
        }
        /* h = h+ - i hx */
        hpc->real[i] = sumr;
        hpc->imag[i] = -sumi;
    }
  }
#ifdef _OPENMP
  if (USETIMERS) openmp_timer_stop("compute_hpc");
#endif
}

/**
  * Function: SPA
  * -------------
  *   Routine to compute the mode-by-mode stationary phase approximation
  *   of the multipolar waveform.
  *   See https://arxiv.org/abs/2012.00027
  * 
  *   @param[in]   TDlm    : multipolar waveform (time domain)
  *   @param[out]  FDlm    : multipolar waveform (frequency domain)
  *
*/
void SPA(Waveform_lm *TDlm, WaveformFD_lm *FDlm)
{
  const int size = TDlm->size;
  double tmpf0 = EOBPars->initial_frequency;  
  double tmpdf = EOBPars->df;
  double tmpsrate = EOBPars->srate_interp/2.;
  /* Determine the size of output frequency array */
  const int interp_size = get_uniform_size(tmpsrate,tmpf0,tmpdf);

  /* if necessary, transform f0, df and srate_interp to geom units */
  if (!(EOBPars->use_geometric_units)){
    double Msun = EOBPars->M;
    double conv = time_units_factor(Msun);
    tmpsrate = tmpsrate/conv; //FIXME: this must be transformed in geom units 
    tmpf0    = tmpf0/conv;
    tmpdf    = tmpdf/conv;
    if (EOBPars->interp_freqs){
      for (int i=0; i < EOBPars->freqs_size; i++)
        EOBPars->freqs[i] = EOBPars->freqs[i]/conv;
    }
  }
  const double half_srate_interp = tmpsrate;
  const double f0 = tmpf0;
  const double df = tmpdf;
  int nmin[KMAX]; // Fmin > f0 index

  const double Pio4 = Pi/4.;
  int *activemode = TDlm->kmask;
  
  /* For each active mode... */
  for (int k = 0; k < KMAX; k++ ) {
    if (!activemode[k]) continue;
    
    /* Compute frequencies */
    D0_x_4(TDlm->phase[k], TDlm->time, size, FDlm->F[k]);
    D0_x_4(FDlm->F[k], TDlm->time, size, FDlm->Fdot[k]);
    
    for (int i=0; i < size; i++){
      FDlm->F[k][i]     = FDlm->F[k][i]/(TwoPi);
      FDlm->Fdot[k][i]  = FDlm->Fdot[k][i]/(TwoPi);
      FDlm->phase[k][i] = (TwoPi * FDlm->F[k][i] * TDlm->time[i] - TDlm->phase[k][i]- Pio4);
      FDlm->ampli[k][i] = TDlm->ampli[k][i]/sqrt(fabs(FDlm->Fdot[k][i])); 
    }

    /* Track multipoles w\ Fmin > f0 */
    nmin[k] = (FDlm->F[k][0]>f0)?(1):(0);

    /* Get last index until Fdot is monotonically increasing (for attachment) */
    int n = 1;
    while(FDlm->Fdot[k][n] > FDlm->Fdot[k][n-1] && n<size) n++;
    /* Prolong the waveform, if necessary  
       - Fill the points between Fmax and half_srate_interp such that 
         FDlm->F[k][n] > half_srate_interp (strictly larger)
	 for later interp;
       - The amplitude behaves as ~ f^{-10/3}
       - For the phase, we express it as phi(f) = (a + b*f);
    */
    if (FDlm->F[k][n] < half_srate_interp) { 
      const int n1 = n-1;
      double Fn1 = FDlm->F[k][n-1];
      double An1 = FDlm->ampli[k][n-1];
      double pn1 =  FDlm->phase[k][n-1];
      double b = TwoPi * TDlm->time[n-1];
      double dfk = (half_srate_interp - FDlm->F[k][n1])/(size-n-1); // one point more
      for (int i=n; i < size; i++){
        FDlm->F[k][i] = Fn1 + (i-n1)*dfk;
        FDlm->ampli[k][i] = An1/pow((FDlm->F[k][i]),10./3)*pow(Fn1,10./3);
        FDlm->phase[k][i] = (pn1 + b*(FDlm->F[k][i] - Fn1));
      }
    } else{
      /*update values with index >= n to assure F is formally increasing (for interpolation).
        Note that we do not care about frequencies higher than srate_interp/2, so
        we can fill F however we want. We also fill A, to avoid inf which may appear
        and would mess up the spline 
      */
      double Fn1 = FDlm->F[k][n-1];
      const int n1 = n-1;
      for (int i=n; i < size; i++){
        FDlm->F[k][i] = Fn1 + (i-n1);
        FDlm->ampli[k][i] = FDlm->ampli[k][n1];
      }
    }
  }

  
  /* Interpolate each mode */
  if (EOBPars->interp_freqs)
    WaveformFD_lm_interp_ap_freqs(FDlm, "");
  else
    WaveformFD_lm_interp_ap (FDlm, interp_size, f0, df, "");

  /* Correct Fmin > f0 */
  for (int k = 0; k < KMAX; k++ ) {
    if (!activemode[k]) continue;
    if (nmin[k]) {
      /* f0 < Fmin , the interpolation returned 0s at those points
	 -> fill points with linear extrapolation starting from innermost */
      int i0 = 0;
      while(FDlm->ampli[k][i0]==0. && i0<size) i0++;
      for (int i = i0-1; i>=0; i--){
        FDlm->ampli[k][i] = 2*FDlm->ampli[k][i+1] - FDlm->ampli[k][i+2]; //TODO: check sign!
        FDlm->phase[k][i] = 2*FDlm->phase[k][i+1] - FDlm->phase[k][i+2]; 
      }
    }
  }
  
}

/**
  * Function: twist_hlm_FD
  * ----------------------
  *   Frequency domain routine to compute the twisted plus and cross
  *   polarizations from the co-precessing frequency domain multipoles.
  *   E16-E19 of https://arxiv.org/pdf/2004.06503.pdf.
  *   Note the change: -m'<-->m' due to EOB phase > 0
  * 
  *   @param[in]   hlm     : multipolar waveform (frequency domain)
  *   @param[in]   spin    : spin dynamics
  *   @param[in]   M       : total mass
  *   @param[in]   amplitude_prefactor: amplitude prefactor
  *   @param[in]   phi     : reference phase
  *   @param[in]   iota    : inclination angle
  *   @param[out]  hpc     : (h+, hx) polarizations
  *
*/
void twist_hlm_FD(WaveformFD_lm *hlm, DynamicsSpin *spin, double M, double amplitude_prefactor, double phi, double iota, WaveformFD *hpc)
{
  
  const int size = hlm->size;
  int *activemode = hlm->kmask;
  int activemode_inertial[KMAX];
  set_multipolar_idx_mask (activemode_inertial, KMAX, EOBPars->use_mode_lm_inertial, EOBPars->use_mode_lm_inertial_size, 1);
  /* initialize the splines for angles */
  gsl_spline *alpha;
  gsl_spline *beta;
  gsl_spline *gamma;

  double *frequencies;

  prolong_euler_angles_FD(spin, hlm);
  frequencies = malloc(spin->size*sizeof(double));

  alpha = gsl_spline_alloc(gsl_interp_cspline, spin->size);
  beta  = gsl_spline_alloc(gsl_interp_cspline, spin->size);
  gamma = gsl_spline_alloc(gsl_interp_cspline, spin->size);
  gsl_interp_accel *acc_al = gsl_interp_accel_alloc ();
  gsl_interp_accel *acc_bt = gsl_interp_accel_alloc ();
  gsl_interp_accel *acc_gm = gsl_interp_accel_alloc ();


  double conv = 1.;
  double Msun = M;
  if (!(EOBPars->use_geometric_units)) {
    Msun = M/MSUN_S;
    conv = time_units_factor(Msun);
  }

  for(int i=0;i<spin->size;i++)
    frequencies[i] = spin->data[EOB_EVOLVE_SPIN_Momg][i];
  
  /* spline for  the angles */
  gsl_spline_init (alpha, frequencies, spin->data[EOB_EVOLVE_SPIN_alp], spin->size);  
  gsl_spline_init (beta,  frequencies, spin->data[EOB_EVOLVE_SPIN_bet], spin->size);  
  gsl_spline_init (gamma, frequencies, spin->data[EOB_EVOLVE_SPIN_gam], spin->size);  
  /* precompute Ylms */
  double Y_real[KMAX] = {0};
  double Y_imag[KMAX] = {0};
  double Y_real_mneg[KMAX] = {0};
  double Y_imag_mneg[KMAX] = {0};
  
  static const int mneg = 1; /* m>0 modes only, add m<0 modes afterwards */
  for (int k = 0; k < KMAX; k++ ) {
    if (!activemode_inertial[k]) continue;
    spinsphericalharm(&Y_real[k], &Y_imag[k], -2, LINDEX[k], MINDEX[k], phi, iota);
    /* add m<0 modes */
    if ( (mneg) && (MINDEX[k]!=0) ) 
      spinsphericalharm(&Y_real_mneg[k], &Y_imag_mneg[k], -2, LINDEX[k], -MINDEX[k], phi, iota); 
  }

  double f022 = hlm->freq[0];
  double omg0 = frequencies[0];
  double omgM = frequencies[spin->size-1];
  /* Loop over frequencies */
  for(int i=0; i<size;i++){
    double hpr = 0.;
    double hpi = 0.;
    double hcr = 0.;
    double hci = 0.;
    double f   = hlm->freq[i];

    //printf("loop over modes:\n");
    /* Loop over modes */
    for(int k=0;k<KMAX;k++){
      if(!activemode[k])
        continue;
      int emm = MINDEX[k];
      int ell = LINDEX[k];
      double f0lm = f022/2.*emm;
      double omg  = 2.*Pi*f/emm;

      if(f < f0lm || omg > omgM){ 
        continue;
      }    
      
      /* avoid interpolation error (due to small numerical differences) */
      if(k == 1 && f == f0lm) 
        omg = omg0;
      double alph = gsl_spline_eval(alpha, omg, acc_al); //evaluate spline
      double bet  = gsl_spline_eval(beta , omg, acc_bt); //evaluate spline
      double gam  = gsl_spline_eval(gamma, omg, acc_gm); //evaluate spline
      double cosmg = cos(-emm * gam);
      double sinmg = sin(-emm * gam);
      double eps   = pow (-1., ell);
      double sumpr = 0.;
      double sumpi = 0.;
      double sumcr = 0.;
      double sumci = 0.;
      /* loop over n */
      for(int n=-ell; n<=ell; n++){
        int j = 0;
        if(n!=0){
          j = KINDEX[ell][abs(n)-1]; // map to linear index (ell,n) -> j
          if(!activemode_inertial[j]) continue;
        }
        double cosma = cos(n * alph);
        double sinma = sin(n * alph);
        double dl_mn = wigner_d_function(ell, emm, n, -bet); //CHECKME
        double dl_mnn= wigner_d_function(ell,-emm, n, -bet);
        
        double rAlm=0;double iAlm=0;
        if (n==0) {
          double Y_real_0, Y_imag_0;
          spinsphericalharm(&Y_real_0, &Y_imag_0, -2, ell, 0, phi, iota);
          rAlm  = cosma *Y_real_0 + sinma*Y_imag_0;
          iAlm  = cosma *Y_imag_0 - sinma*Y_real_0; 
        } else if(n>0) {
          rAlm  = cosma *Y_real[j] + sinma*Y_imag[j];
          iAlm  = cosma *Y_imag[j] - sinma*Y_real[j]; 
        } else {
          rAlm = cosma *Y_real_mneg[j] + sinma*Y_imag_mneg[j];
          iAlm = cosma *Y_imag_mneg[j] - sinma*Y_real_mneg[j]; 
        }
        sumpr += rAlm*(dl_mn + eps*dl_mnn);
        sumpi += iAlm*(dl_mn - eps*dl_mnn);
        sumcr += rAlm*(dl_mn - eps*dl_mnn); 
        sumci += iAlm*(dl_mn + eps*dl_mnn); 
      }
      /* Add stuff to compute h+, hx */
      double Amplm = hlm->ampli[k][i];
      double clm   = cos(hlm->phase[k][i]);
      double slm   = sin(hlm->phase[k][i]);
      double tmpr  = cosmg*clm + sinmg*slm;
      double tmpi  = cosmg*slm - sinmg*clm;
      hpr +=  Amplm*(sumpr*tmpr - sumpi*tmpi); 
      hpi +=  Amplm*(sumpi*tmpr + sumpr*tmpi);
      hcr += -Amplm*(sumci*tmpr + sumcr*tmpi); 
      hci +=  Amplm*(sumcr*tmpr - sumci*tmpi); 
    }
    hpc->preal[i] = amplitude_prefactor*0.5*hpr; 
    hpc->pimag[i] = amplitude_prefactor*0.5*hpi;
    hpc->creal[i] = amplitude_prefactor*0.5*hcr;
    hpc->cimag[i] = amplitude_prefactor*0.5*hci;
    hpc->freq[i]  = f*conv; 
  }
  /* free */
  gsl_spline_free (alpha);
  gsl_spline_free (beta);
  gsl_spline_free (gamma);
  gsl_interp_accel_free (acc_al);
  gsl_interp_accel_free (acc_bt);
  gsl_interp_accel_free (acc_gm);

  free(frequencies);
}

/**
  * Function: compute_hpc_FD
  * ------------------------
  *   Frequency domain routine to compute the (h+, hx) polarizations
  *   from the frequency domain multipoles.
  *   @note : we assume that all modes have a common frequency array
  * 
  *   @param[in]   hflm    : multipolar waveform (frequency domain)
  *   @param[in]   nu      : symmetric mass ratio
  *   @param[in]   M       : total mass
  *   @param[in]   distance: distance to source
  *   @param[in]   amplitude_prefactor: amplitude prefactor
  *   @param[in]   phi     : reference phase
  *   @param[in]   iota    : inclination angle
  *   @param[out]  hpc     : (h+, hx) polarizations
  *
*/
void compute_hpc_FD(WaveformFD_lm *hflm, double nu, double M, double distance, double amplitude_prefactor, double phi, double iota, WaveformFD *hpc)
{  
#ifdef _OPENMP
  if (USETIMERS) openmp_timer_start("compute_hpc_FD");
#endif
  
  double Y_real[KMAX] = {0};
  double Y_imag[KMAX] = {0};
  double Y_real_mneg[KMAX] = {0};
  double Y_imag_mneg[KMAX] = {0};
  static const int mneg = 1; /* m>0 modes only, add m<0 modes afterwards */
  double conv=1.;
  double sumr=0., sumi=0.;
  double sumpr=0., sumpi=0., sumcr=0., sumci=0.;
  
  int *activemode = hflm->kmask;
  
  double Msun = M;
  if (!(EOBPars->use_geometric_units)) {
    Msun = M/MSUN_S;
    conv = time_units_factor(Msun);
    for (int i = 0; i < hflm->size; i++)
	    hflm->freq[i] = hflm->freq[i]*conv; //CHECKME!
  }
  
    /* Precompute Ylm */
  for (int k = 0; k < KMAX; k++ ) {
    if (!activemode[k]) continue;
    spinsphericalharm(&Y_real[k], &Y_imag[k], -2, LINDEX[k], MINDEX[k], phi, iota);
    /* add m<0 modes */
    if ( (mneg) && (MINDEX[k]!=0) ) 
      spinsphericalharm(&Y_real_mneg[k], &Y_imag_mneg[k], -2, LINDEX[k], -MINDEX[k], phi, iota); 
  }
  
  /* Sum up  hlm * Ylm 
   * Note because EOB code defines phase>0, 
   * but the convention is h_lm = A_lm Exp[-I phi_lm] we have
   * h_{l,m>=0} = A_lm ( cos(phi_lm) - I*sin(phi_lm) ) for m>0 and
   * h_{l,m<0}  = (-)^l A_l|m| ( cos(phi_l|m|) + I*sin(phi_l|m|) ) for m<0 below
   * We now agree with, e.g., LALSimSphHarmMode.c: 64-74
   */
  
  const double pm = Pi/2.;
  for (int i = 0; i < hflm->size; i++) {
    hpc->freq[i] = hflm->freq[i]; 
    sumpr = sumpi = sumcr = sumci = 0.;
    
    for (int k=0; k< KMAX; k++){
      if (!activemode[k]) continue;
      double Aki  = 0.5*amplitude_prefactor * hflm->ampli[k][i];
      double cosPhi = cos( hflm->phase[k][i] );
      double sinPhi = sin( hflm->phase[k][i] );
      double cosPhipm = cos( hflm->phase[k][i] + pm );
      double sinPhipm = sin( hflm->phase[k][i] + pm );
      
      /* H_{l-m} = (-)^l H^{*}_{lm} */
      if (LINDEX[k] % 2) {
        sumpr +=  Aki * (cosPhi*(Y_real[k] - Y_real_mneg[k]) - sinPhi* (Y_imag[k] + Y_imag_mneg[k]));
        sumpi +=  Aki * (cosPhi*(Y_imag[k] + Y_imag_mneg[k]) + sinPhi* (Y_real[k] - Y_real_mneg[k]));
        sumcr +=  Aki * (cosPhipm* (Y_real[k] + Y_real_mneg[k]) - sinPhipm* (Y_imag[k] - Y_imag_mneg[k]));
        sumci +=  Aki * (cosPhipm* (Y_imag[k] - Y_imag_mneg[k]) + sinPhipm* (Y_real[k] + Y_real_mneg[k]));
      } else {
        sumpr +=  Aki * (cosPhi* (Y_real[k] + Y_real_mneg[k]) - sinPhi* (Y_imag[k] - Y_imag_mneg[k]));
        sumpi +=  Aki * (cosPhi* (Y_imag[k] - Y_imag_mneg[k]) + sinPhi* (Y_real[k] + Y_real_mneg[k]));
        sumcr +=  Aki * (cosPhipm* (Y_real[k] - Y_real_mneg[k]) - sinPhipm* (Y_imag[k] + Y_imag_mneg[k]));
        sumci +=  Aki * (cosPhipm* (Y_imag[k] + Y_imag_mneg[k]) + sinPhipm* (Y_real[k] - Y_real_mneg[k]));
      }
      
      hpc->preal[i] = sumpr; 
      hpc->pimag[i] = sumpi;
      hpc->creal[i] = sumcr;
      hpc->cimag[i] = sumci;	  
    }
  }
        
#ifdef _OPENMP
    if (USETIMERS) openmp_timer_stop("compute_hpc_FD");
#endif
}


/**
 * Function: eob_wav_hlmNQC_find_a1a2a3_mrg_BHNS_HM
 * ------------------------------------------------
 *   Routine to compute the NQC coefficients a1, a2, a3 for the BHNS case
 *   from the dynamics and waveform multipoles.
 *   Based on NQC fits from arxiv:2001.0908
 * 
 *   @param[in]      dyn_mrg : dynamics (merger)
 *   @param[out]     hlm_mrg : multipolar waveform (merger)
 *   @param[in]      hnqc    : NQC waveform
 *   @param[in]      dyn     : dynamics
 *   @param[in,out]  hlm     : multipolar waveform
 *
 */
void eob_wav_hlmNQC_find_a1a2a3_mrg_BHNS_HM(Dynamics *dyn_mrg, Waveform_lm *hlm_mrg, Waveform_lm *hnqc,
				       Dynamics *dyn, Waveform_lm *hlm)
{
  if (VERBOSE) PRSECTN("BHNS NQC");
  double A_tmp, dA_tmp, omg_tmp, domg_tmp;
  double alpha1[KMAX], omega1[KMAX];
  double c1A[KMAX], c2A[KMAX], c3A[KMAX], c4A[KMAX];
  double c1phi[KMAX], c2phi[KMAX], c3phi[KMAX], c4phi[KMAX];
	  
  const double nu   = EOBPars->nu;
  const double chi1 = EOBPars->chi1;
  const double chi2 = EOBPars->chi2;
  const double X1   = EOBPars->X1;
  const double X2   = EOBPars->X2;
  const double aK   = EOBPars->a1+EOBPars->a2;
  const double Mbh  = EOBPars->Mbhf;
  const double abh  = EOBPars->abhf;
  const double kt2 = EOBPars->kapT2;
  const double lambda = EOBPars->LambdaBl2;

    
  double *t       = hlm_mrg->time;
  double *r       = dyn_mrg->data[EOB_RAD];
  double *w       = dyn_mrg->data[EOB_MOMG]; /* Omega */
  double *pph     = dyn_mrg->data[EOB_PPHI];
  double *pr_star = dyn_mrg->data[EOB_PRSTAR];
  double *Omg_orb = dyn_mrg->data[EOB_OMGORB]; /* Omega orbital */
  double *ddotr   = dyn_mrg->data[EOB_DDOTR];
  
  double P[2], M[4];   
  double max_A[KMAX],max_dA[KMAX],d2max[KMAX],d3max[KMAX],max_omg[KMAX],max_domg[KMAX],maxd2omg[KMAX], DeltaT[KMAX], max_d2A[KMAX], max_d2omg[KMAX];
  double ai[KMAX][2] = {0.};
  double bi[KMAX][2] = {0.};
  
  const int size = hlm_mrg->size;
  for (int i = 0; i < size; i++) {
    hnqc->time[i] = t[i];
  }
  
  double *omg[KMAX], *domg[KMAX];
  double *n1[KMAX],*n2[KMAX],*n3[KMAX],*n4[KMAX],*n5[KMAX],*n6[KMAX],
    *d_n4[KMAX],*d_n5[KMAX],*d_n6[KMAX],*d2_n4[KMAX],*d2_n5[KMAX],*d2_n6[KMAX];
  double *m11[KMAX], *m12[KMAX], *m13[KMAX], *m21[KMAX], *m22[KMAX];
  double *p1tmp[KMAX], *p2tmp[KMAX]; /* RWZ amplitude and derivative */
  
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask[k]){
      omg[k]  = (double*) calloc (size,sizeof(double));
      domg[k] = (double*) calloc (size,sizeof(double));
      m11[k] = (double*) calloc (size,sizeof(double));
      m12[k] = (double*) calloc (size,sizeof(double));
      m13[k] = (double*) calloc (size,sizeof(double));
      m21[k] = (double*) calloc (size,sizeof(double));
      m22[k] = (double*) calloc (size,sizeof(double));
      p1tmp[k] = (double*) calloc (size,sizeof(double));
      p2tmp[k] = (double*) calloc (size,sizeof(double));
      n1[k] = (double*) calloc (size,sizeof(double));
      n2[k] = (double*) calloc (size,sizeof(double));
      n4[k] = (double*) calloc (size,sizeof(double));
      n5[k] = (double*) calloc (size,sizeof(double));
      d_n4[k] = (double*) calloc (size,sizeof(double));
      d_n5[k] = (double*) calloc (size,sizeof(double));
      d2_n4[k] = (double*) calloc (size,sizeof(double));
      d2_n5[k] = (double*) calloc (size,sizeof(double));
    }
  }

  /** omega derivatives */
  const double dt = t[1]-t[0];
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){
      D0(hlm_mrg->phase[k], dt, size, omg[k]);
      D0(omg[k], dt, size, domg[k]);
    }
  }
  
  /** NR fits */
  for (int k=0; k<KMAX; k++) {   
    max_A[k]    = 0.;
    max_dA[k]   = 0.;
    max_omg[k]  = 0.;
    max_domg[k] = 0.;
  }

  /* Higher modes */
  /* Choosing modes using kpostpeak array */
  int kpostpeak_size = EOBPars->kpostpeak_size;  
  int *kpostpeak     = EOBPars->kpostpeak;
  /* old option: 21, 32, 42, 43 and 44 extracted from postpeak */
  // int kpostpeak_size = 5;
  // int kpostpeak[kpostpeak_size] = {0,3,6,7,8};

  QNMHybridFitCab_BHNS_HM(nu, X1, X2, chi1, chi2, aK,  Mbh, abh,  
			  c1A, c2A, c3A, c4A, c1phi, c2phi, c3phi, c4phi,
			  alpha1, omega1);

  /* 22, 31, 33, 41 and 55 fitted directly + 44 dA */
  eob_nqc_point_BHNS_HM(dyn, max_A, max_dA, max_omg, max_domg, abh, kt2);
  
  // Over-writing fits using postpeak quantities for modes in kpostpeak
  for (int j=0; j<kpostpeak_size; j++) {
    int k = kpostpeak[j];
    
    /* Normalizing c1A and c4A */
    int l = LINDEX[k];
    c1A[k] /= sqrt((l+2)*(l+1)*l*(l-1));
    c4A[k] /= sqrt((l+2)*(l+1)*l*(l-1));
    
    eob_nqc_point_postpeak(Mbh,c1A[k],c2A[k],c3A[k],c4A[k],
			   c1phi[k],c2phi[k],c3phi[k],c4phi[k],alpha1[k],omega1[k],
			   &A_tmp,&dA_tmp,&omg_tmp,&domg_tmp);
    
    max_A[k]    = A_tmp;
    max_dA[k]   = dA_tmp;
    max_omg[k]  = omg_tmp;
    max_domg[k] = domg_tmp;
  }

  if (VERBOSE) {
    printf("NR values for NQC determination:\n");
    PRFORMd("A22_mrg",max_A[1]);    
    PRFORMd("dA22_mrg",max_dA[1]);
    PRFORMd("omg22_mrg",max_omg[1]);
    PRFORMd("domg22_mrg",max_domg[1]);
  }

  /** NQC corrections to AMPLITUDE (n1,n2,n3) and PHASE (n4,n5,n6)
   * NQC basis for (2,2) waveform : AMPLITUDE
   * note: n3 and n6 are not used
   */
  double pr_star2, r2, w2;
  for (int k=0; k<KMAX; k++) {   
    if(hlm_mrg->kmask_nqc[k]){ 
      for (int j=0; j<size; j++) {
	pr_star2 = SQ(pr_star[j]);
	r2       = SQ(r[j]);
	w2       = SQ(w[j]); //CHECKME: Omg or Omg_orbital ?
	n1[k][j]  = pr_star2/(r2*w2);         /* [pr*\/(r Omg)]^2 */
	n2[k][j]  = ddotr[j]/(r[j]*w2);       /* [ddot{r}/(r Omg^2)] */
	//n3[k][j]  = n1[k][j]*pr_star2;
	n4[k][j]  = pr_star[j]/(r[j]*w[j]);   /* pr*\/(r Omg) */
	n5[k][j]  = n4[k][j]*r2*w2;              /* (pr*)*(r Omg) */
	//n6[k][j]  = n5[k][j]*pr_star2;
      }
    }
  }

  for (int j=0; j<size; j++) {
    /* l=2,m=1 */
    if(hlm_mrg->kmask_nqc[0]){
      n2[0][j] = cbrt(SQ(w[j]))*n1[0][j];
      n5[0][j] = cbrt(SQ(w[j]))*n4[0][j];
    }
    /* l=3 & l=4 */
    for (int k=2; k<14; k++) {
      if(hlm_mrg->kmask_nqc[k]){
	n5[k][j]  = cbrt(SQ(w[j]))*n4[k][j];
      }
    }
  }
      
#if (DEBUG)
  FILE* fp_dbg;
  fp_dbg = fopen("nqc_nfunc.txt", "w");
  for (int j=0; j<size; j++) {
    fprintf(fp_dbg, "%20.12f\t%.16e\t%.16e\t%.16e\t%.16e\n", t[j], n1[1][j], n2[1][j], n4[1][j], n5[1][j]);
  }
  fclose(fp_dbg);
#endif    
      
  /** Derivatives for the phase */
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){   
      D0(n4[k],dt,size, d_n4[k]);
      D0(n5[k],dt,size, d_n5[k]);
      D0(d_n4[k],dt,size, d2_n4[k]);
      D0(d_n5[k],dt,size, d2_n5[k]);
    }
  }
    
#if (DEBUG)
  fp_dbg = fopen("nqc_dfunc.txt", "w");
  for (int j=0; j<size; j++) {
    fprintf(fp_dbg, "%f\t%.16e\t%.16e\t%.16e\t%.16e\n", t[j], d_n4[1][j], d_n5[1][j], d2_n4[1][j], d2_n5[1][j]);  
  }
  fclose(fp_dbg);
#endif    

  /** Find max Omg */
  int Omgmax_index = 0;
  double Omg_max   = Omg_orb[0];
  for (int j=0; j<size; j++) {
    if (Omg_orb[j] > Omg_max) {
      Omg_max = Omg_orb[j];
      Omgmax_index = j;
    }
  }
  
  /** Time */
  double tOmgOrb_pk = t[Omgmax_index];
  double DeltaT_nqc = eob_nqc_timeshift(nu, chi1);
  double tNQC = tOmgOrb_pk - DeltaT_nqc;

  if (VERBOSE) {
    printf("NQC info:\n");
    PRFORMd("DeltaT_tNQC",DeltaT_nqc);
    PRFORMd("tNQC[bare]",tNQC);
  }

  /** Find jmax: t[jmax] <= tNQC */
  double tmrg[KMAX];
  tmrg[1] = tNQC - 2.;
	
  int jmax = 0;
  for (int j=0; j<size; j++) {
    if(t[j] > tNQC) {
      jmax = j-2;
      break;
    }
  }
  
  double dtmrg[KMAX];
  double t_NQC[KMAX];
  int    j_NQC[KMAX];

  eob_nqc_deltat_lm(dyn, dtmrg);
  
  for (int k=0; k<KMAX; k++) {   
    if(hlm_mrg->kmask_nqc[k]){
      tmrg[k]  = tmrg[1] + dtmrg[k];
      t_NQC[k] = tmrg[k] + 2.;
      
      j_NQC[k] = size-1;
      for (int j=size-2; j>=0; j--) {
	if(t[j] < t_NQC[k]) {
	  break;
	}
	j_NQC[k] = j;
      } 
    }
  }
  
  /** Solve the linear systems */
  
  /* Regge-Wheeler-Zerilli normalized amplitude. 
     The ringdown coefficient refer to this normalization.
     Nagar & Rezzolla, CQG 22 (2005) R167 */      
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){
      double nlm = 1./(sqrt( (LINDEX[k]+2)*(LINDEX[k]+1)*LINDEX[k]*(LINDEX[k]-1) ) );

      if (hlm->ampli[k][0] > 0.) {	
        nNegAmp[k] = 0;
      } else {	
        nNegAmp[k] = 1;	
      }
      
      for (int j=0; j<size; j++) {
        p1tmp[k][j] = fabs(hlm_mrg->ampli[k][j] * nlm);      
      }
    }
  }
  
  /* Matrix elements: waveform amplitude at all points */
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){
      for (int j=0; j<size; j++) {
        m11[k][j] = n1[k][j] * p1tmp[k][j];
        m12[k][j] = n2[k][j] * p1tmp[k][j];
      }
    }
  }

  /* Take FD derivatives */
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){
      D0(m11[k],dt,size, m21[k]);
      D0(m12[k],dt,size, m22[k]);
      D0(p1tmp[k],dt,size, p2tmp[k]);
    }
  }

#if (DEBUG)
  fp_dbg = fopen("nqc_amp_func.txt", "w");
  for (int j=0; j<size; j++) {
    fprintf(fp_dbg, "%e\t%e\t%e\n", t[j], p1tmp[1][j], p2tmp[1][j]);
  }
  fclose(fp_dbg);  
#endif

  double detM = 1.;
  double oodetM = 1.;
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){

      jmax = j_NQC[k];

      /* Computation of ai coefficients at Omega peak */
      P[0]     = max_A[k]  - p1tmp[k][jmax];
      P[1]     = max_dA[k] - p2tmp[k][jmax];
      
      M[0]     = m11[k][jmax];
      M[1]     = m12[k][jmax];
      M[2]     = m21[k][jmax];
      M[3]     = m22[k][jmax];
      
      /* detM     = M[0]*M[3]-M[1]*M[2]; 
        ai[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
        ai[k][1] = (M[0]*P[1] - M[2]*P[0])/detM; */
      /* safe version (amplitude can be zero) */
      oodetM   = 1.0/(M[0]*M[3]-M[1]*M[2]);
      if (isfinite(oodetM)) {
        ai[k][0] = (M[3]*P[0] - M[1]*P[1])*oodetM;
        ai[k][1] = (M[0]*P[1] - M[2]*P[0])*oodetM;
      }

      /* Computation of bi coefficients at Omega peak */
      P[0]     = omg[k][jmax]   - max_omg[k];
      P[1]     = domg[k][jmax]  - max_domg[k];
      
      M[0]     = d_n4[k][jmax];
      M[1]     = d_n5[k][jmax];
      M[2]     = d2_n4[k][jmax];
      M[3]     = d2_n5[k][jmax];

      /* detM     =  M[0]*M[3] - M[1]*M[2];
        bi[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
        bi[k][1] = (M[0]*P[1] - M[2]*P[0])/detM; */
      /* safe version (phase can be zero) */
      oodetM   = 1.0/(M[0]*M[3]-M[1]*M[2]);
      if (isfinite(oodetM)) {
        bi[k][0] = (M[3]*P[0] - M[1]*P[1])*oodetM;
        bi[k][1] = (M[0]*P[1] - M[2]*P[0])*oodetM;
      }
    }

  }

   if((nu==0.25)&&(chi1==0)&&(chi2==0)){// TODO: problematic cases where X12=0 and messes up the fits
     ai[0][0] = ai[1][0];
     ai[0][1] = ai[1][1];
     ai[4][0] = ai[1][0];
     ai[4][1] = ai[1][1];
     ai[13][0] = ai[1][0];
     ai[13][1] = ai[1][1];
   }
  
  if (VERBOSE){
    printf("NQC coefficients for 22 mode:\n");
    PRFORMd("a1",ai[1][0]);
    PRFORMd("a2",ai[1][1]);
    PRFORMd("b1",bi[1][0]);
    PRFORMd("b2",bi[1][1]);
  }

  /** Set amplitude and phase */
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask_nqc[k]){
      for (int j=0; j<size; j++) {
        hnqc->ampli[k][j] = 1. + ai[k][0]*n1[k][j] + ai[k][1]*n2[k][j];
        hnqc->phase[k][j] =      bi[k][0]*n4[k][j] + bi[k][1]*n5[k][j];
      }
    }
  }
  
  /** Multiply merger waveform to NQC */
  for (int k=0; k<KMAX; k++) {
    for (int j=0; j<size; j++) {
      if(hlm_mrg->kmask_nqc[k]){
        hlm_mrg->ampli[k][j] *= hnqc->ampli[k][j];
        hlm_mrg->phase[k][j] -= hnqc->phase[k][j];
      }
    }
  }

  /** Multiply full waveform to NQC */
  r       = dyn->data[EOB_RAD];
  w       = dyn->data[EOB_MOMG]; /* Omega */
  pph     = dyn->data[EOB_PPHI];
  pr_star = dyn->data[EOB_PRSTAR];
  Omg_orb = dyn->data[EOB_OMGORB]; /* Omega orbital */
  ddotr   = dyn->data[EOB_DDOTR];

  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask[k]){
      free(n1[k]);
      free(n2[k]);
      free(n4[k]);
      free(n5[k]);
    }
  }

  const int fullsize = hlm->size;
  
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask[k]){
      n1[k] = (double*) calloc (fullsize,sizeof(double));
      n2[k] = (double*) calloc (fullsize,sizeof(double));
      n4[k] = (double*) calloc (fullsize,sizeof(double));
      n5[k] = (double*) calloc (fullsize,sizeof(double));
    }
  }
	
  for (int k=0; k<KMAX; k++) {   
    if(hlm_mrg->kmask_nqc[k]){   
      for (int j=0; j<fullsize; j++) {
	pr_star2 = SQ(pr_star[j]);
	r2       = SQ(r[j]);
	w2       = SQ(w[j]); //CHECKME: Omg or Omg_orbital ?
	n1[k][j]  = pr_star2/(r2*w2);         /* [pr*\/(r Omg)]^2 */
	n2[k][j]  = ddotr[j]/(r[j]*w2);       /* [ddot{r}/(r Omg^2)] */
	n4[k][j]  = pr_star[j]/(r[j]*w[j]);   /* pr*\/(r Omg) */
	n5[k][j]  = n4[k][j]*r2*w2;              /* (pr*)*(r Omg) */
      }
    }
  }

  for (int j=0; j<fullsize; j++) {
    /* l=2,m=1 */
    if(hlm_mrg->kmask_nqc[0]){
      n2[0][j] = cbrt(SQ(w[j]))*n1[0][j];
      n5[0][j] = cbrt(SQ(w[j]))*n4[0][j];
    }
    /* l=3 & l=4 */
    for (int k=2; k<14; k++) {   
      if(hlm_mrg->kmask_nqc[k]){   
	n5[k][j]  = cbrt(SQ(w[j]))*n4[k][j];
      }
    }
  }
  
  for (int k=0; k<KMAX; k++) {
    if(hlm->kmask_nqc[k]){
      for (int j=0; j<fullsize; j++) {
        hlm->ampli[k][j] *= (1. + ai[k][0]*n1[k][j] + ai[k][1]*n2[k][j]);
        hlm->phase[k][j] -= (bi[k][0]*n4[k][j] + bi[k][1]*n5[k][j]);
      }
    }
  }

  if (EOBPars->output_nqc_coefs) {
    /** Output the NQC coefficients */
    FILE* fp;
    char fname[STRLEN];
    strcpy(fname, EOBPars->output_dir);
    strcat(fname, "/nqc_coefs.txt");
    fp = fopen(fname, "w");
    fprintf(fp, "# q=%e chizA=%e chizB=%e f0=%e\n",EOBPars->q,EOBPars->chi1,EOBPars->chi2,EOBPars->initial_frequency);
    fprintf(fp, "# M=%e LambdaA=[%e,%e,%e] LambdaBl2=[%e,%e,%e]\n",EOBPars->M,
	    EOBPars->LambdaAl2,EOBPars->LambdaAl3,EOBPars->LambdaAl4,
	    EOBPars->LambdaBl2,EOBPars->LambdaBl3,EOBPars->LambdaBl4);
    for (int k=0; k<KMAX; k++) {
      if(hlm->kmask_nqc[k]){
        fprintf(fp, "%d %d %d %e %e %e %e\n", k, LINDEX[k], MINDEX[k], 
		ai[k][0], ai[k][1], 
		bi[k][0], bi[k][1]);
      }
    }  
    fclose(fp);  
  }
  
  /** Free mem */
  for (int k=0; k<KMAX; k++) {
    if(hlm->kmask[k]){
      free(omg[k]);
      free(domg[k]);
      free(m11[k]);
      free(m12[k]);
      free(m13[k]);
      free(m21[k]);
      free(m22[k]);
      free(p1tmp[k]);
      free(p2tmp[k]);
      free(n1[k]);
      free(n2[k]);
      free(n4[k]);
      free(n5[k]);
      free(d_n4[k]);
      free(d_n5[k]);
      free(d2_n4[k]);
      free(d2_n5[k]);
    }
  }

}

/**
 * Function: eob_wav_ringdown_template_td
 * --------------------------------------
 *   Ringdown waveform template for tidal disruption cases
 *   Based on arxiv:2001.0908 and Gonzalez et al. 2022
 * 
 *   @param[in]   x           : time
 *   @param[in]   a1,a2,a3,a4 : amplitude coefficients
 *   @param[in]   b1,b2,b3,b4 : phase coefficients
 *   @param[in]   sigmai      : imaginary part of the QNM complex frequency
 *   @param[out]  psi         : waveform
 *   @param[in]   alpha2      : inverse damping time of the first overtone
 *   @param[in]   Amrg        : Amrg
 *
 */
void eob_wav_ringdown_template_td(double x, double a1, double a2, double a3, double a4, double b1, double b2, double b3, double b4, double sigmai, double *psi, double alpha2, double Amrg)
{  
  double amp   = ( a1 * tanh(a2*x +a3) + a4 ) ;
  double phase = -b1*log((1. + b3*exp(-b2*x) + b4*exp(-2.*b2*x))/(1.+b3+b4));
  
  double sigmar = -(alpha2)*log( (Amrg/exp(1))/amp );

  psi[0] = amp * exp(-sigmar*x); /* amplitude */
  psi[1] = - (phase - sigmai*x); /* phase, minus sign in front by convention */
  
}

/**
 * Function: eob_wav_ringdown_bhns
 * -------------------------------
 *   Ringdown calculation and match to the dynamics
 *   for BHNS systems
 * 
 *   @param[in]  dyn : dynamics
 *   @param[out] hlm : waveform
 *   @return         : 0 if successful
 */
int eob_wav_ringdown_bhns(Dynamics *dyn, Waveform_lm *hlm)
{
  if (VERBOSE) PRSECTN("entered BHNS ringdown model");
  const double Mbh   = EOBPars->Mbhf;
  const double abh   = EOBPars->abhf;
  const double nu    = EOBPars->nu;
  const double q    = EOBPars->q;
  const double chi1  = EOBPars->chi1;
  const double chi2  = EOBPars->chi2;
  const double X1    = EOBPars->X1;
  const double X2    = EOBPars->X2;
  const double aK    = EOBPars->a1+EOBPars->a2;
  const double M = EOBPars->M;
  const int binary = EOBPars->binary;
	
  const double xnu   = (1.-4.*nu);
  const double ooMbh = 1./Mbh;
  /* const double dt = par_get_d("dt"); */	
  const double dt = 0.5;//dyn->dt;
	  
  /* double *Omega = dyn->data[EOB_MOMG]; */
  double *Omega = dyn->data[EOB_OMGORB]; /* use this for spin */
	  
  /* Note:
     dynsize < size , since the wf has been extended 
     but the two time arrays agree up to dynsize */
  const int dynsize = dyn->size; 
  const int size = hlm->size; 
  double *t = hlm->time;
	    
  if (VERBOSE) {
    PRFORMi("ringdown_dynamics_size",dynsize);
    PRFORMi("ringdown_waveform_size",size);
  }
  
  /** Find peak of Omega */
  /* Assume a monotonically increasing function, 
     start from after the peak */
  int index_pk = dynsize-1;
  double Omega_pk = Omega[index_pk];
  for (int j = dynsize-2; j-- ; ) {
    if (Omega[j] < Omega_pk) 
      break;
    index_pk = j;
    Omega_pk = Omega[j]; 
  }
  
  if (VERBOSE) PRFORMi("ringdown_index_pk",index_pk);
  if (index_pk >= dynsize-2) {
    if (VERBOSE) printf("No omega-maximum found.\n");
  }
  
  double tOmg_pk = dyn->time[index_pk]*ooMbh;  
	
  if (VERBOSE) PRFORMd("ringdown_Omega_pk",Omega_pk);
  if (VERBOSE) PRFORMd("ringdown_tOmg_pk",tOmg_pk/ooMbh);
  if (VERBOSE) PRFORMd("ringdown_tOmg_pk",tOmg_pk);
	  
  /** Merger time t_max(A22) */
  double DeltaT_nqc = eob_nqc_timeshift(nu, chi1);
  double tmrg[KMAX], tmatch[KMAX], dtmrg[KMAX];
 
  /* nonspinning case */
  double tmrgA22 = tOmg_pk-(DeltaT_nqc + 2.)/Mbh;
  if (VERBOSE) PRFORMd("ringdown_tmrgA22",tmrgA22);
  
  /* The following values are the difference between the time of the peak of
     the 22 waveform and the other modes. */
  eob_nqc_deltat_lm(dyn, dtmrg);	  
  for (int k=0; k<KMAX; k++) {
    tmrg[k] = tmrgA22 + dtmrg[k]/Mbh;
  }	  
    
  /** Postmerger-Ringdown matching time */
  int idx[KMAX];
  for (int k = 0; k < KMAX; k++) {
    if(hlm->kmask[k]){
      int j  = size-1;
      idx[k] = size-1;
      for (j = size-1; j-- ; ) {  
	      if (t[j] * ooMbh < tmrg[k]) {
	        break;
	      }
      }
	    idx[k]    = j;
      tmatch[k] = (t[idx[k]])*ooMbh;	    
    }
  }
  
  /** Compute QNM */
  double sigma[2][KMAX];
  double a1[KMAX], a2[KMAX], a3[KMAX], a4[KMAX];
  double b1[KMAX], b2[KMAX], b3[KMAX], b4[KMAX];

  /** 
   * BHNS only part
   **/
  double Apeak[KMAX], alpha2[KMAX];
  
  
  /** Postpeak coefficients calculation */   
   if(binary == BINARY_BHNS){
     if (VERBOSE) PRSECTN("No tidal disruption cases");
    QNMHybridFitCab_BHNS_HM(nu, X1, X2, chi1, chi2, aK,  Mbh, abh,  
	     a1, a2, a3, a4, b1, b2, b3, b4, 
	     sigma[0],sigma[1]);
   }else{
     if (VERBOSE) PRSECTN("Tidal disruption cases");
     postpeak_coef(a1, a2, a3, a4, b1, b2, b3, b4, sigma[0],sigma[1], nu, chi1, chi2, X1, X2, aK, Mbh, abh, Apeak, alpha2);
   }
  
  if (VERBOSE) PRFORMd("a1", a1[1] );
  if (VERBOSE) PRFORMd("a2", a2[1] );
  if (VERBOSE) PRFORMd("a3", a3[1] );
  if (VERBOSE) PRFORMd("a4", a4[1] );
  if (VERBOSE) PRFORMd("b1", b1[1] );
  if (VERBOSE) PRFORMd("b2", b2[1] );
  if (VERBOSE) PRFORMd("b3", b3[1] );
  if (VERBOSE) PRFORMd("b4", b4[1] );
 /**
  * end of BHNS only part
  **/

  
 /** Define a time vector for each multipole, scale by mass
      Ringdown of each multipole has its own starting time */
  double *t_lm[KMAX];
  for (int k=0; k<KMAX; k++) {
    t_lm[k] =  malloc ( size * sizeof(double) );
    for (int j = 0; j < size; j++ ) {  
      t_lm[k][j] = t[j] * ooMbh;
    }
  }
	  
  /** Compute Ringdown waveform for t>=tmatch */
  double t0, tm, psi[2];
  double Deltaphi[KMAX];
  int n0 = 2/dt*ooMbh;
  int index_rng;
	  
  for (int k = 0; k < KMAX; k++) {
    if(hlm->kmask[k]){

      /* Ringdown attachment index */      
      index_rng = idx[k]+n0;
      if (index_rng > dynsize -1) index_rng = dynsize - 1;
      
      /* Calculate Deltaphi */
      t0 = t_lm[k][index_rng] - tmatch[k];
      if(binary == BINARY_BHNS_TD){
        eob_wav_ringdown_template_td(t0, a1[k], a2[k], a3[k], a4[k], b1[k], b2[k], b3[k], b4[k], sigma[1][k], psi, alpha2[k], Apeak[k]);
      }else{
        eob_wav_ringdown_template(t0, a1[k], a2[k], a3[k], a4[k], b1[k], b2[k], b3[k], b4[k], sigma[0][k], sigma[1][k], psi);
      }
      Deltaphi[k] = psi[1] - hlm->phase[k][index_rng];
     
      /* Compute and attach ringdown */
      for (int j = index_rng-1; j < size ; j++ ) {
	
	      tm = t_lm[k][j] - tmatch[k];
	      if(binary == BINARY_BHNS_TD){
          eob_wav_ringdown_template_td(tm, a1[k], a2[k], a3[k], a4[k], b1[k], b2[k], b3[k], b4[k], sigma[1][k], psi, alpha2[k], Apeak[k]);
        }else{
          eob_wav_ringdown_template(tm, a1[k], a2[k], a3[k], a4[k], b1[k], b2[k], b3[k], b4[k], sigma[0][k], sigma[1][k], psi);
        }
	      hlm->phase[k][j] = psi[1] - Deltaphi[k];
	      hlm->ampli[k][j] = psi[0];
	
	      if(nNegAmp[k]==1) {
	        hlm->ampli[k][j] = -hlm->ampli[k][j];
	      }
      }
    }
  } 
	  
  /** Free mem. */
  for (int k=0; k<KMAX; k++) {
    free(t_lm[k]);
  }
  
  return 0;
  
}
