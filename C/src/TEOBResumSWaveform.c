/** \file TEOBResumSWaveform.c
 *  \brief EOB IMR Waveform
 * 
 *  This file contains functions for the calculation of the EOB-factorized 
 *  waveform, including NQCs and ringdown attachment. 
 */

#include "TEOBResumS.h"

/** hlmNewt coefficients for amplitude */
static const double ChlmNewt_ampli[54] = 
{2.1137745587232057, 6.341323676169617, 
0.1412325034218127, 1.7864655618418102, 4.9229202032627635, 
0.023872650234580958, 0.2250735048768909, 1.7053495827316825, 4.763908164911493, 
0.001122165903318321, 0.06333806741197714, 0.2945348827200268, 1.755276012972272, 5.0817902739730565, 
0.00014954736544380072, 0.005296595280255638, 0.10342548105284892, 0.3713362832603404, 1.8983258440274462, 5.727111757630886, 
5.184622059790144e-06, 0.0012191691413436815, 0.011783593824950922, 0.14639129995388936, 0.4653654097044924, 2.125638973894669, 6.685178460621457, 
5.54485779151375621e-7, 0.0000763473331250837455, 0.00353250998285463003, 0.0204988821800401766, 0.19579402814926015, 0.584571015778149663, 2.44207899966355693, 7.9955401278763745,
0., 0., 0., 0., 0., 0., 0., 0., 0.,
0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

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
static const double ChlmNewt_phase[54] = {4.71238898038469, 3.141592653589793, 4.71238898038469, 0.0, 1.5707963267948966, 1.5707963267948966, 0.0, 4.71238898038469, 3.141592653589793, 1.5707963267948966, 3.141592653589793, 4.71238898038469, 0.0, 1.5707963267948966, 4.71238898038469, 3.141592653589793, 1.5707963267948966, 0.0, 4.71238898038469, 3.141592653589793, 4.71238898038469, 0.0, 1.5707963267948966, 3.141592653589793, 4.71238898038469, 0.0, 1.5707963267948966, 4.7123889803846898577, 0.0, 1.5707963267948966192, 3.1415926535897932385, 
   4.7123889803846898577, 0.0, 1.5707963267948966192, 3.1415926535897932385,
   0., 0., 0., 0., 0., 0., 0., 0., 0.,
  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

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
int nNegAmp[54];	
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
  
  double phim[KMAX] = {
    phi,phix2,
    phi,phix2,phix3,
    phi,phix2,phix3,phix4,
    phi,phix2,phix3,phix4,phix5,
    phi,phix2,phix3,phix4,phix5,phix6,
    phi,phix2,phix3,phix4,phix5,phix6,phix7,
    phi,phix2,phix3,phix4,phix5,phix6,phix7,8.*phi,
    0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0
  };

  double Alm[KMAX] = {
    pv23, p1 * vphi2, 
    pv23, pv34, pv23, 
    pv45, pv34, pv45, pv34, 
    pv45, pv56, pv45, pv56, pv45, 
    pv67, pv56, pv67, pv56, pv67, pv56, 
    pv67, pv78, pv67, pv78, pv67, pv78, pv67, 
    pv89, pv78, pv89, pv78, pv89, pv78, pv89, pv78,
    0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0
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

  double phim[KMAX] = {
    phi,phix2,
    phi,phix2,phix3,
    phi,phix2,phix3,phix4,
    phi,phix2,phix3,phix4,phix5,
    phi,phix2,phix3,phix4,phix5,phix6,
    phi,phix2,phix3,phix4,phix5,phix6,phix7,
    phi,phix2,phix3,phix4,phix5,phix6,phix7,8.*phi,
    0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0
  };

  double Alm[KMAX] = {
    fact21, fact22, 
    fact31, fact32, fact33, 
    fact41, fact42, fact43, fact44, 
    pv45,   pv56,   pv45,   pv56,   fact55, 
    pv67,   pv56,   pv67,   pv56,   pv67,   pv56, 
    pv67,   pv78,   pv67,   pv78,   pv67,   pv78, pv67, 
    pv89,   pv78,   pv89,   pv78,   pv89,   pv78, pv89, pv78,
    0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0
  };
    
  /* Compute hlmNewt (without phase factor) in complex Polar coords */
#pragma omp simd
  for (int k = 0; k < KMAX; k++) {
    hlmNewt->phase[k] = - phim[k] + ChlmNewt_phase[k];
    hlmNewt->ampli[k] = ChlmNewt_ampli[k] * Alm[k];
  }

  /* Correcting phase in case of negative frequency */
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
 *   @param[out] hlmNewt: Newtonian prefactor of multipolar waveform
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

  // Set all the unused multipoles to zero (ell = 9, 10, only used in the 22PN flux)
  for (int k=35; k<KMAX; k++) {
    flm[k] = 0.;
  }

  /* Amplitudes */
#pragma omp simd
  for (int k = 0; k < 35; k++) {
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
    for (int k=35; k<KMAX; k++) clm[k][0] = 0.;
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

  // Set all the unused multipoles to zero (ell = 9, 10, only used in the 22PN flux)
  for (int k=35; k<KMAX; k++) {
    flm[k] = 0.;
  }
  
  if (kmaxTaylor+kmaxPade32+kmaxPade42+kmaxPade51+19 != KMAX) {
    errorexit("Wrong function: not all multipoles are written.\n");
  }
  
  /** Amplitudes */
#pragma omp simd
  for (int k = 0; k < 35; k++) {
      flm[k] = gsl_pow_int(rholm[k], LINDEX[k]);
  }

}

/** Taylor expanded amplitudes up to 22PN
 *  Ref: arXiv:2001.09082 */
void eob_wav_flm_22PN(double x,double nu, double *rholm, double *flm)
{
  /** Coefficients */
  static double clm[KMAX][45]; // up to v^44 = x^22
  
  const double nu2 = nu*nu;
  const double nu3 = nu*nu2;
  const double nu4 = nu*nu3;
  const double nu5 = nu*nu4;
  const double Pi2 = SQ(Pi);
  
  //static int firstcall = 1;
  if (EOBPars->firstcall[FIRSTCALL_EOBWAVFLMHM]) {
    if (0) printf("Precompute some rholm coefs\n");
    EOBPars->firstcall[FIRSTCALL_EOBWAVFLMHM] = 0;
    
    for (int k=0; k<KMAX; k++) clm[k][0] = 1.;
    for (int k=0; k<KMAX; k++) for (int n=1; n<45; n++) clm[k][n] = 0.;
  }

  /** Compute Logs */
  const double v     = sqrt(x);
  const double logv  = log(v);
  //printf("x, v = %.5f %.5f\n", x, v);
  const double logv2 = logv*logv;
  const double logv3 = logv2*logv;
  const double logv4 = logv3*logv;
  const double logv5 = logv4*logv;
  const double logv6 = logv5*logv;
  const double logv7 = logv6*logv;

  
  /** Coefs with Eulerlogs */
  // FIXME: separate constant parts (in firstcall) from log terms
  /* l = 2 */
  clm[0][2] = -1.053571428571429 + 0.2738095238095238*nu;
  clm[0][4] = -0.8327841553287982 - 0.7789824263038549*nu + 0.1311649659863946*nu2;
  clm[0][6] = 1.624720394052015 - 1.019047619047619*logv;
  clm[0][8] = 0.08155386512520203 + 1.073639455782313*logv;
  clm[0][10] = -2.768567922179298 + 0.8486467106683943*logv;
  clm[0][12] = -0.5361453134743556 - 3.120073452848125*logv + 0.5192290249433107*logv2;
  clm[0][14] = -0.996335447573493 + 1.459749053216591*logv - 0.5470448655652737*logv2;
  clm[0][16] = -10.54257790186801 + 10.44083666609426*logv - 0.4324057049596105*logv2;
  clm[0][18] = -10.32906936731376 + 6.990263491277985*logv + 2.335901437369426*logv2 - 0.1763730338696325*logv3;
  clm[0][20] = -14.27817030672405 + 40.39518474564478*logv - 1.5298989309567*logv2 + 0.1858215892555057*logv3;
  clm[0][21] = -1.340412865531645;
  clm[0][22] = 33.8967533555577 + 148.0787709813931*logv - 35.289807970345*logv2 + 0.1468806680338994*logv3;
  clm[0][23] = 1.412220697613698;
  clm[0][24] = 169.0923080456532 + 194.7477688870285*logv - 60.84906786630315*logv2 - 1.046918966496203*logv3 + 0.04493313005726352*logv4;
  clm[0][25] = 1.116274596013625;
  clm[0][26] = 490.0493367192367 + 455.1102828785716*logv - 182.6484009976958*logv2 + 0.7867118827954247*logv3 - 0.04734026202461692*logv4;
  clm[0][27] = -47.67875956208643 + 17.75727900909065*logv;
  clm[0][28] = 1886.815243083545 + 566.8832844338378*logv - 720.5426076920055*logv2 + 101.9211796099349*logv3 - 0.03741959876101723*logv4;
  clm[0][29] = 47.82919920597807 - 18.70856181314908*logv;
  clm[0][30] = 3937.08550451999 + 20.22302286404711*logv - 1129.925694731496*logv2 + 201.1840217341791*logv3 + 0.3312855024570407*logv4 - 0.009157799840242279*logv5;
  clm[0][31] = 41.60350546951885 - 14.78798060052336*logv;
  clm[0][32] = 8861.798719512112 - 1130.67916432647*logv - 2457.095443275486*logv2 + 555.47068230031*logv3 - 0.2684537804891634*logv4 + 0.009648396260255259*logv5;
  clm[0][33] = -630.6932132154036 + 657.1478557614215*logv - 117.6208338173576*logv2;
  clm[0][34] = 21634.83300598225 - 9718.483093371418*logv - 4115.425545354976*logv2 + 2340.073782376239*logv3 - 231.7409137915696*logv4 + 0.00762647060462637*logv5;
  clm[0][35] = 584.9755434907948 - 660.507848518708*logv + 123.921949914716*logv2;
  clm[0][36] = 37978.98805214537 - 22417.79371776483*logv - 3756.710608294051*logv2 + 3716.14747559116*logv3 - 444.2211049422257*logv4 - 0.0806792095846884*logv5 + 0.001555372353818927*logv6;
  clm[0][37] = 604.3690277804423 - 580.9768637070492*logv + 97.95276673965709*logv2;
  clm[0][38] = 78421.52101798659 - 60066.08093788073*logv - 1598.572457760974*logv2 + 8560.348330512021*logv3 - 1295.847863431751*logv4 + 0.06857850995120722*logv5 - 0.001638695872773512*logv6;
  clm[0][39] = -4075.092794994711 + 9337.670032826439*logv - 4521.847354477074*logv2 + 519.3986661585537*logv3;
  clm[0][40] = 171251.1547708798 - 165476.1064332664*logv + 27809.38514563435*logv2 + 17163.82382104097*logv3 - 5713.157561507821*logv4 + 424.5939856553635*logv5 - 0.00129528945189686*logv6;
  clm[0][41] = 3398.417180779158 - 8804.297163710187*logv + 4553.15822412416*logv2 - 547.2235947027619*logv3;
  clm[0][42] = 305251.5713460777 - 300727.6504066538*logv + 63476.80049104249*logv2 + 20909.54353919977*logv3 - 8648.396843668505*logv4 + 746.191691438565*logv5 + 0.01593778225006525*logv6 - 0.0002264283562702383*logv7;
  clm[0][43] = 4648.916019909193 - 9340.477962861954*logv + 4080.833005195625*logv2 - 432.5469794757556*logv3;

  clm[1][2] = -1.023809523809524 + 0.6547619047619048*nu;
  clm[1][4] = -1.94208238851096 - 1.560137944066515*nu + 0.4625614134542706*nu2;
  clm[1][6] = 4.732393863861011 - 4.076190476190476*logv - 2.902228713904598*nu - 1.930155846609928*nu2 + 0.2715020968103451*nu3;
  clm[1][8] = 5.776972352223222 + 4.173242630385488*logv;
  clm[1][10] = -14.59770105471162 + 7.916297736025627*logv;
  clm[1][12] = 0.4974231817838048 - 42.72063485619346*logv + 8.307664399092971*logv2;
  clm[1][14] = 104.3177089938809 + 0.4403253291725122*logv - 8.505465932404708*logv2;
  clm[1][16] = -53.86659117081395 + 87.9402970965808*logv - 16.13416871913794*logv2;
  clm[1][18] = -77.99672500530377 - 509.2438045585568*logv + 134.8223049123425*logv2 - 11.28787416765648*logv3;
  clm[1][20] = 1761.547589750222 - 446.8144489684867*logv - 49.78799745830767*logv2 + 11.55663307641021*logv3;
  clm[1][21] = -171.5728467880506;
  clm[1][22] = -890.4577156155182 + 75.3511402926456*logv + 41.07933637032255*logv2 + 21.92198162473346*logv3;
  clm[1][23] = 175.6579145687184;
  clm[1][24] = -5105.190117862853 - 10780.80762822826*logv + 3915.743408482643*logv2 - 248.0713644073432*logv3 + 11.50288129465946*logv4;
  clm[1][25] = 333.2086040937622;
  clm[1][26] = 30315.24787592984 - 18414.46653373547*logv + 3827.787820767194*logv2 + 134.0775490250378*logv3 - 11.77675942072278*logv4;
  clm[1][27] = -22169.39592882982 + 9091.726852654413*logv;
  clm[1][28] = -57187.42910591078 + 12641.78660543444*logv + 19036.5562040386*logv2 - 3757.976962800744*logv3 - 22.33954317949029*logv4;
  clm[1][29] = 20874.78465970963 - 9308.196539622376*logv;
  clm[1][30] = -310058.6208055035 - 46187.2434962891*logv + 147523.4009279655*logv2 - 35921.34249029682*logv3 + 318.9166560524071*logv4 - 9.377587036408094*logv5;
  clm[1][31] = 43982.48817042576 - 17656.88260169232*logv;
  clm[1][32] = 372873.6925490746 - 372410.7305056553*logv + 268448.7053331478*logv2 - 40060.20726791388*logv3 - 204.3258185389913*logv4 + 9.600862918227334*logv5;
  clm[1][33] = -979411.1649772588 + 1.227027326959803e6*logv - 240887.4676579484*logv2;
  clm[1][34] = -1.603324113110717e6 + 1.736487139592903e6*logv + 181175.8770271784*logv2 - 294532.9747054492*logv3 + 38810.95056144209*logv4 + 18.21204663013685*logv5;
  clm[1][35] = 757885.3256284707 - 1.159669556466518e6*logv + 246622.8835545662*logv2;
  clm[1][36] = -8.527047667952658e6 + 1.021477798965746e7*logv - 572312.8506590727*logv2 - 2.159195066061053e6*logv3 + 403911.9306563247*logv4 - 313.8966501876038*logv5 + 6.370805161242324*logv6;
  clm[1][37] = 2.026653230520771e6 - 2.429218973112111e6*logv + 467823.308551505*logv2;
  clm[1][38] = 2.84387511493005e6 - 1.064660497630082e6*logv + 7.870004501707526e6*logv2 - 3.055724647658831e6*logv3 + 266198.4574566851*logv4 + 221.7612541520891*logv5 - 6.522490998414761*logv6;
  clm[1][39] = -1.742163470414021e7 + 5.942993402117034e7*logv - 3.38950350333873e7*logv2 + 4.254913873170872e6*logv3;
  clm[1][40] = -2.194754424599069e7 + 5.500040992298716e7*logv - 2.203812256448159e7*logv2 - 4.104521016010827e6*logv3 + 2.800874096853552e6*logv4 - 289141.126342549*logv5 - 12.37262850428344*logv6;
  clm[1][41] = 5.980215256835689e6 - 4.72652574267284e7*logv + 3.214334281556101e7*logv2 - 4.356221346341607e6*logv3;
  clm[1][42] = -9.810136771804681e7 + 4.010620602238436e8*logv - 3.147005611883822e8*logv2 + 4.497372331093769e6*logv3 + 2.922575995738459e7*logv4 - 4.092585653923977e6*logv5 + 249.8706753499774*logv6 - 3.709802189131585*logv7;
  clm[1][43] = 3.989155944454545e7 - 1.217206957769425e8*logv + 6.692652098858528e7*logv2 - 8.263393297716107e6*logv3;
  clm[1][44] = -5.571949408177452e7 + 1.959978515295102e8*logv + 1.375489272496405e8*logv2 - 1.248898386226567e8*logv3 + 2.280845694693602e7*logv4 - 772286.445029612*logv5 - 188.1490154881398*logv6 + 3.798130812682337*logv7;

  /* l=3 */
  clm[2][2] = -0.7222222222222222 - 0.2222222222222222*nu;
  clm[2][4] = 0.01416947250280584 - 0.9455667789001122*nu - 0.4652076318742985*nu2;
  clm[2][6] = 1.385551684086816 - 0.4126984126984127*logv - 2.539296242237244*nu + 2.302086630790334*nu2 - 0.5813492634480289*nu3;
  clm[2][8] = 0.9154593365697892 + 0.2980599647266314*logv;
  clm[2][10] = 1.44237045159571 - 0.005847718810681774*logv;
  clm[2][12] = 4.029155909057081 - 0.7862233992242928*logv + 0.08515998992189468*logv2;
  clm[2][14] = 7.491041633637398 - 0.222958090627257*logv - 0.06150443716581283*logv2;
  clm[2][16] = 13.68246817639354 - 0.5983020500868071*logv + 0.001206672135537509*logv2;
  clm[2][18] = 28.35559136100461 - 2.223770443753643*logv + 0.2064795814331186*logv2 - 0.01171513088872625*logv3;
  clm[2][20] = 57.67389261138288 - 3.503599428041022*logv + 0.01405394222395017*logv2 + 0.008460927864080071*logv3;
  clm[2][22] = 115.8032318299523 - 8.510698185489489*logv + 0.1240860532635022*logv2 - 0.0001659972249945779*logv3;
  clm[2][24] = 237.5974620260455 - 23.58729920976238*logv + 0.5976091406934322*logv2 - 0.03449093808995427*logv3 + 0.001208703980582867*logv4;
  clm[2][26] = 482.5470042237087 - 62.53866244448888*logv + 1.797590299401554*logv2 + 0.002462343151621072*logv3 - 0.0008729528748654042*logv4;
  clm[2][27] = -0.002026323303902714;
  clm[2][28] = 961.9534590697911 - 169.0624275846552*logv + 8.664166135014963*logv2 - 0.01715627929469155*logv3 + 0.0000171266978169009*logv4;
  clm[2][29] = 0.001463455719485293;
  clm[2][30] = 1896.51955760205 - 440.3043522014814*logv + 32.68324697096954*logv2 - 0.1044581640452629*logv3 + 0.004186544522114774*logv4 - 0.00009976604284176048*logv5;
  clm[2][31] = -0.00002871193233644418;
  clm[2][32] = 3645.322742174219 - 1086.759686947645*logv + 113.6340253823079*logv2 - 2.192280031732195*logv3 - 0.0007075748181221594*logv4 + 0.00007205325316349368*logv5;
  clm[2][33] = -0.07371901748467261 + 0.01588894781155461*logv;
  clm[2][34] = 6804.160152181997 - 2554.825930810783*logv + 355.8393392683031*logv2 - 13.43875643280251*logv3 + 0.001778990108803667*logv4 - 1.413632200760074e-6*logv5;
  clm[2][35] = 0.04935880250363309 - 0.01147535119723389*logv;
  clm[2][36] = 12377.07890085131 - 5732.016661606501*logv + 1003.094207027645*logv2 - 53.54389853600818*logv3 + 0.01339903512557408*logv4 - 0.0003973873175734808*logv5 + 6.862214586999398e-6*logv6;
  clm[2][37] = -0.003927486583909853 + 0.00022513800911434*logv;
  clm[2][38] = 21775.27385750876 - 12224.79171654166*logv + 2627.137842998088*logv2 - 196.2961034299442*logv3 + 2.840165939554994*logv4 + 0.00009583668992412362*logv5 - 4.956043868388454e-6*logv6;
  clm[2][39] = -1.041893988130982 + 0.5863054632544544*logv - 0.06229476364212681*logv2;
  clm[2][40] = 36980.34520642485 - 24915.34139770045*logv + 6425.56535318099*logv2 - 628.1448261273665*logv3 + 17.80177642189373*logv4 - 0.0001475717004604658*logv5 + 9.723396089884108e-8*logv6;
  clm[2][41] = 0.5955272931628665 - 0.3929974567734887*logv + 0.04499066263042492*logv2;
  clm[2][42] = 60666.20927184655 - 48777.04851682075*logv + 14820.83352536496*logv2 - 1791.733708534699*logv3 + 71.23435070161489*logv4 - 0.001349266387954836*logv5 + 0.00003089863244229468*logv6 - 4.045750096643636e-7*logv7;
  clm[2][43] = -0.1432609043745549 + 0.0309134476224812*logv - 0.0008826839404959046*logv2;

  clm[3][2]  = 1.214814814814815/(-1. + 3.*nu) - (4.12962962962963*nu)/(-1. + 3.*nu) + (1.185185185185185*nu2)/(-1. + 3.*nu);
  clm[3][4]  =  -0.9006908592093777/pow(1. - 3.*nu,2) + (5.019357151764559*nu)/pow(1. - 3.*nu,2) - (2.946505175208879*nu2)/pow(1. - 3.*nu,2) - (12.6817308891383*nu3)/pow(1. - 3.*nu,2) + (1.923955605437087*nu4)/pow(1. - 3.*nu,2);
  clm[3][6]  = 2.979648070988587 - 1.650793650793651*logv;
  clm[3][8]  = 0.4849109714557107 + 2.005408583186361*logv;
  clm[3][10] = -4.527041916522104 + 1.486854751710719*logv;
  clm[3][12] = 2.981937332103236 - 8.349318813028327*logv + 1.362559838750315*logv2;
  clm[3][14] = 10.72131796946142 + 3.366976418364392*logv - 1.655257878185568*logv2;
  clm[3][16] = -13.55256602723969 + 10.56306329541586*logv - 1.227245191888212*logv2;
  clm[3][18] = -7.134753974410764 - 32.03204531982684*logv + 9.723053689860247*logv2 - 0.7497683768784802*logv3;
  clm[3][20] = 37.21233167087015 + 3.320431531427986*logv - 6.218903508882837*logv2 + 0.9108297319116352*logv3;
  clm[3][22] = -115.6445165603524 + 63.79673021578208*logv - 11.26907231721168*logv2 + 0.6753095235786989*logv3;
  clm[3][24] = -240.9868820205025 - 123.1411723671516*logv + 54.69951313519433*logv2 - 6.908354699846434*logv3 + 0.309428219029214*logv4;
  clm[3][26] = -192.9976740013522 + 253.8314662318033*logv - 48.70365820833436*logv2 + 5.31484866968053*logv3 - 0.3758979845984526*logv4;
  clm[3][27] = -1.037477531598189;
  clm[3][28] = -1719.605317063841 + 1187.385350436742*logv - 197.9662270509499*logv2 + 7.604340080887957*logv3 - 0.2786991684610503*logv4;
  clm[3][29] = 1.260343075422986;
  clm[3][30] = -3810.009993282014 + 1471.851602132424*logv - 111.4553947338178*logv2 - 48.88779279270105*logv3 + 3.494093626665713*logv4 - 0.1021604278699627*logv5;
  clm[3][31] = 0.9344465293455976;
  clm[3][32] = -2654.449385345373 + 8495.956344885118*logv - 2776.600625953676*logv2 + 209.5777258410875*logv3 - 2.974587859040989*logv4 + 0.124106001264251*logv5;
  clm[3][33] = -141.1951925702934 + 32.54056511806385*logv;
  clm[3][34] = -9355.994380794742 + 23453.1345375199*logv - 9766.41932909049*logv2 + 948.7655987485625*logv3 - 3.717467268768494*logv4 + 0.09201496355539439*logv5;
  clm[3][35] = 167.2675486638445 - 39.53076058787015*logv;
  clm[3][36] = -5233.017463425225 + 38067.81334336466*logv - 25873.81319031021*logv2 + 2965.001586431733*logv3 + 29.26621721591043*logv4 - 1.365906363086054*logv5 + 0.02810763094834953*logv6;
  clm[3][37] = 129.0856017867594 - 29.30898955534763*logv;
  clm[3][38] = 106564.2227202865 + 75397.0274845208*logv - 93250.19174981053*logv2 + 19355.69977381382*logv3 - 987.3695397633907*logv4 + 1.239992365938875*logv5 - 0.0341455664854024*logv6;
  clm[3][39] = -7090.632913410583 + 4496.221371632673*logv - 510.3187037563029*logv2;
  clm[3][40] = 282345.2302749777 + 55597.94439354783*logv - 241470.3381572856*logv2 + 68260.16547283211*logv3 - 4800.490699004713*logv4 + 1.418571706284772*logv5 - 0.02531628626920904*logv6;
  clm[3][41] = 8032.05899234506 - 5328.509301372881*logv + 619.9427216002494*logv2;
  clm[3][42] = 849091.8542268969 - 395227.6223481429*logv - 534715.5436373495*logv2 + 206352.2123465735*logv3 - 17402.3329642525*logv4 - 13.10492493626098*logv5 + 0.4342157406851593*logv6 - 0.006628556958340933*logv7;

  clm[4][2] = -1.166666666666667 + 0.6666666666666667*nu;
  clm[4][4] = -1.696717171717172 - 1.87979797979798*nu + 0.4515151515151515*nu2;
  clm[4][6] = 5.309863370123018 - 3.714285714285714*logv - 2.923857908280391*nu - 1.778172753172753*nu2 + 0.2592376759043426*nu3;
  clm[4][8] = 3.542183599617419 + 4.333333333333333*logv;
  clm[4][10] = -14.63919982033249 + 6.302092352092352*logv;
  clm[4][12] = 4.243421536501655 - 37.08943155815107*logv + 6.897959183673469*logv2;
  clm[4][14] = 71.24767224144459 + 7.104913605873704*logv - 8.047619047619048*logv2;
  clm[4][16] = -69.78772913353108 + 83.84119683967105*logv - 11.70388579674294*logv2;
  clm[4][18] = -69.19476239290907 - 300.3395244395172*logv + 101.1335249894268*logv2 - 8.540330417881438*logv3;
  clm[4][20] = 853.9704787572587 - 105.387098873188*logv - 50.82351699876943*logv2 + 9.963718820861678*logv3;
  clm[4][22] = -339.9628065465908 + 818.5013258658652*logv - 210.4295568479133*logv2 + 14.49052527215793*logv3;
  clm[4][24] = -273.1545381160376 - 3894.754610602063*logv + 1237.083623323762*logv2 - 165.1454097244646*logv3 + 7.930306816604193*logv4;
  clm[4][26] = 15127.8666721651 - 2401.744979704086*logv - 194.4660930462477*logv2 + 109.5122407399099*logv3 - 9.252024619371558*logv4;
  clm[4][27] = -39.88412159071711;
  clm[4][28] = 2230.873079022809 + 7051.051973552006*logv - 2339.204430867415*logv2 + 328.2859466723239*logv3 - 13.45548775271807*logv4;
  clm[4][29] = 46.53147518916997;
  clm[4][30] = 15881.67452722046 - 83075.05985007992*logv + 21665.9966368085*logv2 - 2559.392470684648*logv3 + 190.4294634800389*logv4 - 5.891085063763115*logv5;
  clm[4][31] = 67.67207398182533;
  clm[4][32] = 322567.8560350076 - 84022.99225946555*logv + 16649.47699761884*logv2 - 370.1075881760622*logv3 - 144.9501180217917*logv4 + 6.8729325743903*logv5;
  clm[4][33] = -11797.03112981659 + 2814.679437973465*logv;
  clm[4][34] = 76346.77384481159 - 15123.20026477306*logv + 42844.01194159982*logv2 - 3425.015075408338*logv3 - 367.751485266987*logv4 + 9.995505187733426*logv5;
  clm[4][35] = 13374.85032773748 - 3283.792677635709*logv;
  clm[4][36] = 146249.7481142762 - 2.023970405372989e6*logv + 745748.4273048159*logv2 - 95630.42747908866*logv3 + 3504.309692162421*logv4 - 169.0071447481617*logv5 + 3.646862182329547*logv6;
  clm[4][37] = 20240.76745205043 - 4775.714935288816*logv;
  clm[4][38] = 6.106268229886903e6 - 1.757995500611342e6*logv + 1.398014176101096e6*logv2 - 306035.8125222888*logv3 + 16543.87526775864*logv4 + 139.8133643872955*logv5 - 4.254672546051138*logv6;
  clm[4][39] = -1.252973962146134e6 + 845694.0850483589*logv - 99317.97445420655*logv2;
  clm[4][40] = -3.175101201696982e6 + 2.817262924442127e6*logv + 3.865477583005913e6*logv2 - 1.341224020704408e6*logv3 + 101128.0727170837*logv4 + 319.9233291818797*logv5 - 6.187693687644502*logv6;
  clm[4][41] = 1.345957974322137e6 - 959236.4974256593*logv + 115870.9701965743*logv2;
  clm[4][42] = -1.848205187253995e7 - 1.95289969739077e7*logv + 2.353200091128772e7*logv2 - 8.169751118581583e6*logv3 + 761414.6003627911*logv4 - 3569.73930725743*logv5 + 121.6752966864953*logv6 - 1.935069729399352*logv7;
  clm[4][43] = 2.193662084745242e6 - 1.450749937042828e6*logv + 168514.5127166196*logv2;

  /* l=4 */
  clm[5][2] = 1.140151515151515/(-1. + 2.*nu) - (2.623106060606061*nu)/(-1. + 2.*nu) + (0.5454545454545455*nu2)/(-1. + 2.*nu);
  clm[5][4] = -0.3677899278751551;
  clm[5][6] = 0.4101694432072357 - 0.2266955266955267*logv;
  clm[5][8] = -0.4648052623565983 + 0.2584672482399755*logv;
  clm[5][10] = -0.7387476952133781 + 0.08337633141296807*logv;
  clm[5][12] = -0.6064347876689996 - 0.1562970344443505*logv + 0.02569543091188113*logv2;
  clm[5][14] = -1.18107965637858 + 0.1775562070982992*logv - 0.02929668448665235*logv2;
  clm[5][16] = -2.263124072737848 + 0.1907568494545724*logv - 0.009450520681801792*logv2;
  clm[5][18] = -3.706658167287949 + 0.07024476996126146*logv + 0.02489235795421463*logv2 - 0.001941679748079137*logv3;
  clm[5][20] = -6.679770676424708 + 0.3442188781720788*logv - 0.02830782637984619*logv2 + 0.00221380910671144*logv3;
  clm[5][22] = -12.58797150037285 + 0.5749885358417176*logv - 0.02426128409421681*logv2 + 0.0007141302545026751*logv3;
  clm[5][24] = -23.74242314635865 + 0.8979545383624628*logv + 0.001662724950989621*logv2 - 0.002423284247983942*logv3 + 0.0001100425282912094*logv4;
  clm[5][26] = -45.95886035247543 + 2.160837440985962*logv - 0.04996973130507833*logv2 + 0.00275737732303317*logv3 - 0.0001254651553623258*logv4;
  clm[5][28] = -90.74167234200507 + 5.28078513282558*logv - 0.07293254092756606*logv2 + 0.002032756568656165*logv3 - 0.00004047253354342361*logv4;
  clm[5][30] = -180.4543300058243 + 14.22513647575558*logv - 0.2364800507805363*logv2 - 0.001004399796227395*logv3 + 0.0001680705387904571*logv4 - 4.989229781976621e-6*logv5;
  clm[5][32] = -359.9645949855758 + 40.19701045533469*logv - 1.249123253602897*logv2 + 0.004776340604061357*logv3 - 0.0001913122527701756*logv4 + 5.688477895359708e-6*logv5;
  clm[5][33] = -3.446128067861758e-6;
  clm[5][34] = -716.860508561584 + 110.3150628017081*logv - 4.938513947158461*logv2 + 0.006153164760398559*logv3 - 0.0001265077189415723*logv4 + 1.834988461665758e-6*logv5;
  clm[5][35] = 3.929108137978747e-6;
  clm[5][36] = -1414.210912216234 + 292.8885291800591*logv - 18.31500909039385*logv2 + 0.196796281756955*logv3 + 0.0001153092767516145*logv4 - 9.013602427849329e-6*logv5 + 1.88506012205033e-7*logv6;
  clm[5][37] = 1.267451193527424e-6;
  clm[5][38] = -2750.874802817697 + 746.9357664165309*logv - 61.27039064970969*logv2 + 1.361448460374703*logv3 - 0.000337175103632049*logv4 + 0.00001026265291181843*logv5 - 2.149254154307384e-7*logv6;
  clm[5][39] = -0.0001221164372933856 + 0.00001953054543510397*logv;
  clm[5][40] = -5263.141538774922 + 1827.798453531415*logv - 185.0286379636443*logv2 + 5.919914322619015*logv3 - 0.0003882655601914555*logv4 + 6.248237993546056e-6*logv5 - 6.93306126329222e-8*logv6;
  clm[5][41] = 0.0001392214193683998 - 0.0000222677809695693*logv;
  
  clm[6][2] = 0.8681818181818182/(-1. + 3.*nu) - (2.674242424242424*nu)/(-1. + 3.*nu) + (0.2159090909090909*nu2)/(-1. + 3.*nu);
  clm[6][4] = -0.362197915720643/pow(1. - 3.*nu,2) + (0.9328882481155208*nu)/pow(1. - 3.*nu,2) + (3.797934060888606*nu2)/pow(1. - 3.*nu,2) - (9.611541110404747*nu3)/pow(1. - 3.*nu,2) - (1.196804473935156*nu4)/pow(1. - 3.*nu,2);
  clm[6][6] = 2.079331428035551 - 0.9067821067821068*logv;
  clm[6][8] = 0.8835845511552439 + 0.7872517381608291*logv;
  clm[6][10] = 0.06645231153714223 + 0.3284345890892526*logv;
  clm[6][12] = 4.702094939698349 - 2.898515836725528*logv + 0.4111268945900981*logv2;
  clm[6][14] = 9.706789086659544 + 0.07826280740706208*logv - 0.356932894848676*logv2;
  clm[6][16] = 12.40009352668919 + 0.306654264541846*logv - 0.148909304317234*logv2;
  clm[6][18] = 31.48745492752259 - 9.010940574157729*logv + 1.773453224135425*logv2 - 0.1242675038770648*logv3;
  clm[6][20] = 72.0721168798459 - 7.404358430862183*logv - 0.4342326860074991*logv2 + 0.1078867874569062*logv3;
  clm[6][22] = 130.8174956967989 - 10.35501858434278*logv - 0.3053889325355716*logv2 + 0.04500943089607977*logv3;
  clm[6][24] = 289.036016052724 - 48.34738677061541*logv + 6.750901990075693*logv2 - 0.6748711622829651*logv3 + 0.0281708872425496*logv4;
  clm[6][26] = 644.4887797012424 - 71.37025143987651*logv + 2.277953590346271*logv2 + 0.2517776382439659*logv3 - 0.0244574521060317*logv4;
  clm[6][28] = 1328.160464664332 - 138.0347664499968*logv + 4.105889871861052*logv2 + 0.142589541257063*logv3 - 0.01020343664325272*logv4;
  clm[6][30] = 2954.41324368615 - 416.1993746900165*logv + 36.32168667345375*logv2 - 3.001274714800263*logv3 + 0.1844614943807857*logv4 - 0.00510897129674406*logv5;
  clm[6][32] = 6639.334448384189 - 885.1527465070923*logv + 47.22987155819863*logv2 - 0.2277151181539897*logv3 - 0.08439960600244658*logv4 + 0.004435515989445979*logv5;
  clm[6][33] = -0.00705767028298088;
  clm[6][34] = 14501.17822571699 - 2229.808444602314*logv + 143.3290972020545*logv2 - 1.006849504621694*logv3 - 0.04372322173383112*logv4 + 0.001850458755157289*logv5;
  clm[6][35] = 0.006127341018406128;
  clm[6][36] = 32406.07034333956 - 6484.206944006967*logv + 655.5802403238486*logv2 - 25.29128822100557*logv3 + 0.9333297444416572*logv4 - 0.03916078447813328*logv5 + 0.0007721206259918151*logv6;
  clm[6][37] = 0.002556273466339195;
  clm[6][38] = 71693.25511453399 - 17083.87730696127*logv + 2202.693552138683*logv2 - 90.05777647988499*logv3 - 0.0833677118667691*logv4 + 0.02026156516314422*logv5 - 0.0006703410889292577*logv6;
  clm[6][39] = -0.9679858717548558 + 0.1599942282043718*logv;
  clm[6][40] = 153844.6556775055 - 46065.49018132709*logv + 7907.076803347628*logv2 - 441.8677051616003*logv3 + 0.1624220380056738*logv4 + 0.009996734517837534*logv5 - 0.0002796604814191536*logv6;
  clm[6][41] = 0.8214109129353904 - 0.1389040799410682*logv;
  clm[6][42] = 323698.3617247479 - 122669.3423728269*logv + 27907.01072888079*logv2 - 2289.755889964479*logv3 + 42.66921953473507*logv4 - 0.2215163751225839*logv5 + 0.006780960793692491*logv6 - 0.0001000207382752539*logv7;
  
  clm[7][2] = 1.261363636363636/(-1. + 2.*nu) - (3.107954545454545*nu)/(-1. + 2.*nu) + (0.9090909090909091*nu2)/(-1. + 2.*nu);
  clm[7][4] = -0.9783218202252293;
  clm[7][6] = 3.686131565231406 - 2.04025974025974*logv;
  clm[7][8] = 0.743362898003657 + 2.573509445100354*logv;
  clm[7][10] = -6.009171404553495 + 1.996030622823163*logv;
  clm[7][12] = 4.946612919993593 - 12.64905580489035*logv + 2.081329903862371*logv2;  clm[7][14] = 18.98430683931334 + 4.952111234417816*logv - 2.625313856008218*logv2;   clm[7][16] = -17.80822924068248 + 17.27748630424444*logv - 2.036210460035837*logv2;
  clm[7][18] = 0.8156837250491665 - 59.0763213821536*logv + 18.13530345472942*logv2 - 1.415484536349691*logv3;
  clm[7][20] = 117.1264487738335 - 4.603320766125338*logv - 11.65077660990744*logv2 + 1.785440721986542*logv3;
  clm[7][22] = -94.27688407033645 + 96.57873336773973*logv - 22.74349157770715*logv2 + 1.384799408102294*logv3;
  clm[7][24] = -96.74223430545005 - 375.8273568121428*logv + 123.3856985935259*logv2 - 15.89153364327114*logv3 + 0.7219890281186247*logv4;
  clm[7][26] = 634.6812648900818 - 52.80306699347598*logv - 46.70780933174426*logv2 + 12.41141457348714*logv3 - 0.9106907059223561*logv4;
  clm[7][28] = -1506.978891775617 + 723.7357314789912*logv - 172.8460210467084*logv2 + 18.94837050975608*logv3 - 0.7063376201716572*logv4;
  clm[7][30] = -3244.997410951023 - 2903.35556880409*logv + 957.1638435186035*logv2 - 135.783413156892*logv3 + 9.920503252350164*logv4 - 0.2946090293959375*logv5;
  clm[7][32] = -1055.851965233171 + 1052.428703342151*logv - 401.8455512188769*logv2 + 78.00509933852762*logv3 - 8.61973644340931*logv4 + 0.3716091166244212*logv5;
  clm[7][33] = -0.6104712488375068;
  clm[7][34] = -34329.892952971 + 11413.58869216156*logv - 1933.158095119156*logv2 + 176.845845280216*logv3 - 11.44034722413615*logv4 + 0.2882224418934217*logv5;
  clm[7][35] = 0.770026234329128;
  clm[7][36] = -85226.10176892666 - 10523.18351008881*logv + 5569.284873847756*logv2 - 1078.288370412259*logv3 + 100.2773061936844*logv4 - 4.788608936556825*logv5 + 0.1001798236322549*logv6;
  clm[7][37] = 0.5972373433578786;
  clm[7][38] = -140939.3950913998 + 76856.5986609465*logv - 19620.20154238238*logv2 + 1673.655858395749*logv3 - 69.12683529341212*logv4 + 4.451375664751876*logv5 - 0.1263631866270488*logv6;
  clm[7][39] = -183.7430042796191 + 31.13799778973127*logv;
  clm[7][40] = -684077.5886217394 + 317771.5278751064*logv - 84419.82065231645*logv2 + 7194.236343175255*logv3 - 124.9102234895702*logv4 + 5.392730925256984*logv5 - 0.0980081074057501*logv6;
  clm[7][41] = 228.4745243494189 - 39.27633812113831*logv;
  
  clm[8][2] = 1.222727272727273/(-1. + 3.*nu) - (4.446969696969697*nu)/(-1. + 3.*nu) + (1.988636363636364*nu2)/(-1. + 3.*nu);
  clm[8][4] = -1.613202365816002/pow(1. - 3.*nu,2) + (7.375660021796385*nu)/pow(1. - 3.*nu,2) - (0.9897216924489652*nu2)/pow(1. - 3.*nu,2) - (21.23238503920322*nu3)/pow(1. - 3.*nu,2) + (3.949850007946599*nu4)/pow(1. - 3.*nu,2);
  clm[8][6] = 5.472074339668243 - 3.627128427128427*logv;
  clm[8][8] = 2.925123431023674 + 4.434988849534304*logv;
  clm[8][10] = -14.63369058267876 + 5.851292159762054*logv;
  clm[8][12] = 5.370820281253527 - 36.05616125218197*logv + 6.578030313441569*logv2;
  clm[8][14] = 65.30459993664586 + 9.208464683118369*logv - 8.043137065071736*logv2;
  clm[8][16] = -75.2652996565951 + 79.22545405916301*logv - 10.61169406405332*logv2;
  clm[8][18] = -76.72203517730976 - 277.1826851682075*logv + 94.78485656717226*logv2 - 7.953120248132144*logv3;
  clm[8][20] = 692.9942182342524 - 77.62699564148358*logv - 52.64183462217874*logv2 + 9.724497030670667*logv3;
  clm[8][22] = -558.7630444601686 + 782.8295659920613*logv - 191.1000363197621*logv2 + 12.82999239990593*logv3;
  clm[8][24] = -984.4198929205594 - 3117.53161173697*logv + 1101.401378974578*logv2 - 150.1383912429906*logv3 + 7.211747134092697*logv4;
  clm[8][26] = 9952.132062418888 - 1160.583776822272*logv - 308.6221853472857*logv2 + 107.101276600929*logv3 - 8.817999904867889*logv4;
  clm[8][28] = -5109.050120775663 + 9515.870316179451*logv - 2556.225889141014*logv2 + 288.3804366413404*logv3 - 11.63400753838511*logv4;
  clm[8][30] = -6507.50925115469 - 49350.46662256401*logv + 15297.43930610327*logv2 - 2214.322925487458*logv3 + 168.3693368147121*logv4 - 5.231586607865918*logv5;
  clm[8][32] = 193665.4491065994 - 18984.41564463201*logv - 3084.430315321398*logv2 + 1110.669696134406*logv3 - 136.5217784871756*logv4 + 6.396803625072417*logv5;
  clm[8][33] = -14.45410873954484;
  clm[8][34] = 9456.909625346922 + 127656.8505895217*logv - 35781.09417112569*logv2 + 4720.876951991524*logv3 - 313.4861344828937*logv4 + 8.439607892780612*logv5;
  clm[8][35] = 17.6734329588071;
  clm[8][36] = 195540.1331416446 - 987604.7063141745*logv + 265910.3309328462*logv2 - 35599.40827113254*logv3 + 2952.316199788473*logv4 - 145.517394246187*logv5 + 3.162606084062475*logv6;
  clm[8][37] = 23.31740241439549;
  clm[8][38] = 4.587389694283606e6 - 598258.4961439706*logv + 24010.72109518641*logv2 + 9733.502023347284*logv3 - 1851.998967982397*logv4 + 127.6212650602657*logv5 - 3.86700471187639*logv6;
  clm[8][39] = -7597.898177090714 + 1310.672717450213*logv;
  clm[8][40] = 2.883410123761246e6 + 1.055915032646601e6*logv - 216716.8789916695*logv2 + 45694.41139821931*logv3 - 5991.438969803028*logv4 + 265.1242625543783*logv5 - 5.101923616953667*logv6;
  clm[8][41] = 9151.166925397304 - 1602.595277245943*logv;
  clm[8][42] = 1.10603185406636e7 - 2.38940791967202e7*logv + 6.790833682275566e6*logv2 - 891625.6733438705*logv3 + 61408.13786525783*logv4 - 2931.250771037528*logv5 + 102.1008523425845*logv6 - 1.63873977590176*logv7;
  
  /* l=5 */
  clm[9][2] = 0.8179487179487179/(-1. + 2.*nu) - (1.605128205128205*nu)/(-1. + 2.*nu) + (0.02051282051282051*nu2)/(-1. + 2.*nu);
  clm[9][4] = -0.1047896120973044;
  clm[9][6] = 0.4595801176156909 - 0.1441491841491841*logv;
  clm[9][8] = 0.07326833469224709 + 0.1179066403681788*logv;
  clm[9][10] = 0.2112256289378055 + 0.01510533709113591*logv;
  clm[9][12] = 0.6954698480021734 - 0.09173765931341589*logv + 0.0103894936454377*logv2;
  clm[9][14] = 1.246204046232756 + 0.01028758250508941*logv - 0.00849807300742212*logv2;
  clm[9][16] = 2.359916265689464 - 0.02777696094558434*logv - 0.001088711008992826*logv2;
  clm[9][18] = 4.716558813975265 - 0.1224199022140409*logv + 0.008449104034246973*logv2 - 0.0004992123442376592*logv3;
  clm[9][20] = 9.275629264210046 - 0.1729560416399149*logv - 0.002244167522672963*logv2 + 0.0004083300969533674*logv3;
  clm[9][22] = 18.25689736149785 - 0.3444685825917005*logv + 0.001809498928846822*logv2 + 0.00005231226790685031*logv3;
  clm[9][24] = 36.11660196112405 - 0.7081691462754341*logv + 0.01074600830371054*logv2 - 0.0004942516927241667*logv3 + 0.00001799026303476506*logv4;
  clm[9][26] = 71.51957020653513 - 1.364902839552711*logv + 0.01171832540320107*logv2 + 0.0001800356873037547*logv3 - 0.0000147151125848463*logv4;
  clm[9][28] = 141.879663689403 - 2.703326654984462*logv + 0.02510248826739793*logv2 - 0.00007769567646348809*logv3 - 1.885192684941505e-6*logv4;
  clm[9][30] = 282.1002994594472 - 5.440879738997143*logv + 0.05349511327302148*logv2 - 0.0006243346668801238*logv3 + 0.00002099267035636245*logv4 - 5.186563478181219e-7*logv5;
  clm[9][32] = 562.0751345159587 - 11.15148922968072*logv + 0.1001863843710637*logv2 - 0.0005143817942264499*logv3 - 9.0900380186022e-6*logv4 + 4.242342947537971e-7*logv5;
  clm[9][34] = 1122.265895951212 - 23.71219456411683*logv + 0.2142430854937493*logv2 - 0.001217748385169895*logv3 + 2.466587916618658e-6*logv4 + 5.434979749966559e-8*logv5;
  clm[9][36] = 2245.018302911014 - 52.72548636931074*logv + 0.5550006756125668*logv2 - 0.002711082641094472*logv3 + 0.00002695358874891025*logv4 - 6.969280399125082e-7*logv5 + 1.246064823196296e-8*logv6;
  clm[9][38] = 4497.621984475868 - 122.4964061078011*logv + 1.722594190051092*logv2 - 0.004889387143783602*logv3 + 0.00001632250107761989*logv4 + 3.370806627871926e-7*logv5 - 1.019217124614406e-8*logv6;
  clm[9][39] = -4.961918099547508e-9;
  clm[9][40] = 9018.505162575807 - 294.7705061468165*logv + 5.939811219012561*logv2 - 0.02555665117428381*logv3 + 0.00004424274815447412*logv4 - 6.150078063558686e-8*logv5 - 1.305746494708361e-9*logv6;
  clm[9][41] = 4.058594548091423e-9;
  
  clm[10][2] = -1.15956043956044/(1. - 5.*nu + 5.*nu2) + (6.203589743589744*nu)/(1. - 5.*nu + 5.*nu2) - (7.687179487179487*nu2)/(1. - 5.*nu + 5.*nu2) + (1.61025641025641*nu3)/(1. - 5.*nu + 5.*nu2);
  clm[10][4] = -0.4629337197600934;
  clm[10][6] = 1.222304898271311 - 0.5765967365967366*logv;
  clm[10][8] = -0.3626116833872829 + 0.6685987653372269*logv;
  clm[10][10] = -1.221010927835005 + 0.266926072074258*logv;
  clm[10][12] = 0.3307525613092741 - 1.112609980396291*logv + 0.1662318983270032*logv2;
  clm[10][14] = 0.1920326995102395 + 0.6819876853707688*logv - 0.1927559330930262*logv2;
  clm[10][16] = -2.180265106537973 + 0.8928305478329907*logv - 0.07695435103530125*logv2;
  clm[10][18] = -1.631924541211857 - 1.358262815770904*logv + 0.438341220227426*logv2 - 0.03194959003121019*logv3;
  clm[10][20] = -1.950854512964425 + 0.8129701352685247*logv - 0.3329542452989289*logv2 + 0.03704748066036593*logv3;
  clm[10][22] = -8.712267452265219 + 2.06483071358303*logv - 0.3118322158023505*logv2 + 0.01479054255795813*logv3;
  clm[10][24] = -14.56651593073065 - 1.483558757642448*logv + 0.8113521445117501*logv2 - 0.1068469883530284*logv3 + 0.004605507336899856*logv4;
  clm[10][26] = -28.18101968315925 + 2.995833040120419*logv - 0.5971114134054209*logv2 + 0.09019751832532053*logv3 - 0.005340364111974427*logv4;
  clm[10][28] = -73.26734973030815 + 7.458221607031352*logv - 0.8666444003054761*logv2 + 0.0703953197147043*logv3 - 0.002132044642853452*logv4;
  clm[10][30] = -158.2133907307724 + 3.693154825280175*logv + 1.501207635668794*logv2 - 0.2526038011992649*logv3 + 0.01865943020754285*logv4 - 0.0005311041001657568*logv5;
  clm[10][32] = -357.7952498579033 + 24.37671036635808*logv - 1.749650616880006*logv2 + 0.2030155574634444*logv3 - 0.01677919464907281*logv4 + 0.0006158473038405567*logv5;
  clm[10][34] = -868.0019392800233 + 64.08284996009644*logv - 3.35263051026827*logv2 + 0.2261221800764066*logv3 - 0.0116554456108678*logv4 + 0.000245865996669571*logv5;
  clm[10][36] = -2032.671467239365 + 143.6220665393683*logv - 1.09011012955591*logv2 - 0.5781429700923455*logv3 + 0.05265057700026299*logv4 - 0.002527448855324798*logv5 + 0.00005103881515812028*logv6;
  clm[10][38] = -4810.6304165439 + 481.1544786762609*logv - 24.96702313294012*logv2 + 0.5862784147607627*logv3 - 0.04465761468465213*logv4 + 0.002370561081390453*logv5 - 0.00005918259093939398*logv6;
  clm[10][39] = -0.00004064803307149319;
  clm[10][40] = -11495.2791108801 + 1478.331450832749*logv - 102.3971068701585*logv2 + 2.024929174512621*logv3 - 0.04224663510523071*logv4 + 0.001518001998171691*logv5 - 0.00002362758855329646*logv6;
  
  clm[11][2] = 0.9615384615384615/(-1. + 2.*nu) - (2.179487179487179*nu)/(-1. + 2.*nu) + (0.4512820512820513*nu2)/(-1. + 2.*nu);
  clm[11][4] = -0.5788010707241476;
  clm[11][6] = 2.660600792814397 - 1.297342657342657*logv;
  clm[11][8] = 0.9978370986881283 + 1.247444862829478*logv;
  clm[11][10] = -1.289390363554146 + 0.7509033191660411*logv;
  clm[11][12] = 5.534516111977199 - 5.516365287619975*logv + 0.8415489852804538*logv2;
  clm[11][14] = 13.37736964717371 + 0.6907079676987834*logv - 0.809181716615821*logv2;
  clm[11][16] = 10.03228055710087 + 2.867805289285182*logv - 0.4870894537471466*logv2;
  clm[11][18] = 34.93324756403466 - 20.29434794429616*logv + 4.917590103684364*logv2 - 0.3639257989492536*logv3;
  clm[11][20] = 98.7255401790005 - 12.08737197853274*logv - 1.735813708138272*logv2 + 0.3499286528358208*logv3;
  clm[11][22] = 143.9681657421122 - 5.94213849841859*logv - 2.635440982663179*logv2 + 0.2106406420959689*logv3;
  clm[11][24] = 346.2859166359856 - 114.7452247833317*logv + 23.80255573766739*logv2 - 2.705769071812915*logv3 + 0.1180341157710936*logv4;
  clm[11][26] = 878.176755281876 - 127.0341816926767*logv + 2.37432944539667*logv2 + 1.3075419160646*logv3 - 0.11349434208759*logv4;
  clm[11][28] = 1601.058897856996 - 175.9181770318077*logv - 1.967321718028918*logv2 + 1.474913794752041*logv3 - 0.06831827259028697*logv4;
  clm[11][30] = 3705.446467044735 - 883.7661782978658*logv + 140.8712920049335*logv2 - 15.81553619065005*logv3 + 1.065422658593875*logv4 - 0.03062613868231228*logv5;
  clm[11][32] = 8996.76228588485 - 1278.242573747789*logv + 68.72179787953713*logv2 + 2.223426814456598*logv3 - 0.604702907867135*logv4 + 0.02944821027145412*logv5;
  clm[11][34] = 18454.24168338591 - 2314.637922818737*logv + 96.36464949410922*logv2 + 3.901885292229647*logv3 - 0.5870921767963577*logv4 + 0.01772644186146858*logv5;
  clm[11][36] = 43268.55602537584 - 8211.262969595608*logv + 1027.842729484291*logv2 - 100.3745861977236*logv3 + 7.219528195045492*logv4 - 0.3251835835585632*logv5 + 0.006622099357042627*logv6;
  clm[11][38] = 105097.4665137779 - 14671.58415804216*logv + 938.7690834320644*logv2 - 16.76855490254335*logv3 - 2.06274048053198*logv4 + 0.2037666936904899*logv5 - 0.006367403227925603*logv6;
  clm[11][39] = -0.007910900150224882;
  clm[11][40] = 234899.7885251351 - 32459.42742768366*logv + 2104.788064654028*logv2 - 38.46692006034082*logv3 - 2.428138502360837*logv4 + 0.1805426691650564*logv5 - 0.003832878198297962*logv6;
  clm[11][41] = 0.007606634759831617;
  
  clm[12][2] = -1.278241758241758/(1. - 5.*nu + 5.*nu2) + (7.034358974358974*nu)/(1. - 5.*nu + 5.*nu2) - (9.348717948717949*nu2)/(1. - 5.*nu + 5.*nu2) + (2.441025641025641*nu3)/(1. - 5.*nu + 5.*nu2);
  clm[12][4] = -1.044214241436219;
  clm[12][6] = 4.124773279437509 - 2.306386946386946*logv;
  clm[12][8] = 0.9710490992022113 + 2.94812010553549*logv;
  clm[12][10] = -7.150739401733221 + 2.408362095679844*logv;
  clm[12][12] = 6.158685718798627 - 16.0386506873299*logv + 2.659710373232052*logv2;
  clm[12][14] = 25.61702283571343 + 6.101331051811898*logv - 3.399752863893981*logv2;
  clm[12][16] = -22.65932978080194 + 23.30621185483266*logv - 2.777307449824551*logv2;
  clm[12][18] = 2.833263325555434 - 83.93937581391293*logv + 26.02063230436739*logv2 - 2.044773761997452*logv3;
  clm[12][20] = 183.9336990904313 - 10.68542153262521*logv - 16.65473965555006*logv2 + 2.613715208742238*logv3;
  clm[12][22] = -104.2622355399484 + 143.6349003349161*logv - 34.73424702864574*logv2 + 2.135181882792854*logv3;
  clm[12][24] = -52.14895383387937 - 600.1338892889957*logv + 198.5063094652831*logv2 - 25.78970925377271*logv3 + 1.179009878246363*logv4;
  clm[12][26] = 1411.657184951861 - 151.3318607121659*logv - 70.70284658817904*logv2 + 20.19892492834681*logv3 - 1.507059659754033*logv4;
  clm[12][28] = -1282.659592596175 + 1088.159670246273*logv - 293.2416727508439*logv2 + 32.74448481595133*logv3 - 1.231138905658836*logv4;
  clm[12][30] = -1917.512346050212 - 5601.583697582023*logv + 1743.833218426178*logv2 - 247.1710513174765*logv3 + 18.20596677734418*logv4 - 0.543850598569735*logv5;
  clm[12][32] = 10086.77831130398 - 937.0574639759549*logv - 476.6626520515402*logv2 + 139.107232781099*logv3 - 15.91047108867021*logv4 + 0.6951725453366107*logv5;
  clm[12][34] = -25667.90208269209 + 11463.07778929654*logv - 2680.776693544287*logv2 + 340.6347347905221*logv3 - 22.36355331808892*logv4 + 0.5678965402401297*logv5;
  clm[12][36] = -56476.49106930208 - 57618.4646658878*logv + 17672.59010353309*logv2 - 2477.873186486669*logv3 + 206.478677809795*logv4 - 9.936685926158626*logv5 + 0.2090549868876607*logv6;
  clm[12][38] = 11554.03938990409 + 6524.026827937141*logv - 5384.965309911169*logv2 + 1184.957408448934*logv3 - 141.139464735497*logv4 + 9.305952115904399*logv5 - 0.2672228140084911*logv6;
  clm[12][39] = -0.3329886869216722;
  clm[12][40] = -601813.8258609105 + 172390.4802636982*logv - 31364.12679920941*logv2 + 3617.340243171601*logv3 - 272.6828887934495*logv4 + 11.92251838595499*logv5 - 0.2182981945513574*logv6;
  
  clm[13][2] = 1.248717948717949/(-1. + 2.*nu) - (3.328205128205128*nu)/(-1. + 2.*nu) + (1.312820512820513*nu2)/(-1. + 2.*nu);
  clm[13][4] = -1.574972762280455;
  clm[13][6] = 5.561803803137422 - 3.603729603729604*logv;
  clm[13][8] = 2.686645836960601 + 4.500041838503377*logv;
  clm[13][10] = -14.6765788895478 + 5.675775968497862*logv;
  clm[13][12] = 5.945650657963235 - 35.9742122079572*logv + 6.493433528398563*logv2;
  clm[13][14] = 63.89993464362914 + 10.21134952600708*logv - 8.108466995718206*logv2;
  clm[13][16] = -77.85134521438472 + 77.98127383041856*logv - 10.2269809409064*logv2;
  clm[13][18] = -78.69138322820005 - 273.3751616382571*logv + 93.52613021099908*logv2 - 7.800192878713426*logv3;
  clm[13][20] = 651.8569833837604 - 69.108781893588*logv - 54.24449873755437*logv2 + 9.74024085111138*logv3;
  clm[13][22] = -633.5605022216484 + 771.6292678778662*logv - 185.7220355956478*logv2 + 12.28509132450762*logv3;
  clm[13][24] = -1187.565950340686 - 2965.471497866674*logv + 1073.460460880897*logv2 - 146.8298708728066*logv3 + 7.027446497955103*logv4;
  clm[13][26] = 8548.252255242965 - 1001.519985714124*logv - 324.3404067792052*logv2 + 108.2194308325415*logv3 - 8.775298575651628*logv4;
  clm[13][28] = -7846.017515101346 + 9457.658894567254*logv - 2475.081570931726*logv2 + 277.4059258014163*logv3 - 11.06803682266246*logv4;
  clm[13][30] = -15364.57911108881 - 44135.98819743921*logv + 14505.74059781701*logv2 - 2139.69502281804*logv3 + 163.349957905819*logv4 - 5.065003396661347*logv5;
  clm[13][32] = 145241.9015543615 - 12472.94981795313*logv - 3532.096902201742*logv2 + 1119.156591893042*logv3 - 136.2912754893975*logv4 + 6.3247606517284*logv5;
  clm[13][34] = -95297.78456055093 + 140879.3279866681*logv - 35920.84788885175*logv2 + 4516.256877819965*logv3 - 298.8523580470891*logv4 + 7.977242390599606*logv5;
  clm[13][36] = -143692.2827580565 - 777061.1821991202*logv + 231211.2932077125*logv2 - 33078.73634942783*logv3 + 2831.034395721566*logv4 - 140.1246309796823*logv5 + 3.042150447256582*logv6;
  clm[13][38] = 2.990907482988476e6 - 205607.6947329314*logv - 56064.30652650905*logv2 + 14988.88842910431*logv3 - 1837.046412967574*logv4 + 126.19119362526*logv5 - 3.798787866189629*logv6;
  clm[13][39] = -6.057028930111704;
  clm[13][40] = -466789.1721522239 + 2.250386582634537e6*logv - 575370.0979179772*logv2 + 73826.34354105101*logv3 - 5675.351907457648*logv4 + 250.6615419779415*logv5 - 4.791304093188419*logv6;
  clm[13][41] = 7.563520740934359;
  
  /* l = 6 */
  clm[14][2] = -1.118055555555556/(1. - 4.*nu + 3.*nu2) + (4.819444444444444*nu)/(1. - 4.*nu + 3.*nu2) - (4.652777777777778*nu2)/(1. - 4.*nu + 3.*nu2) + (0.8611111111111111*nu3)/(1. - 4.*nu + 3.*nu2);
  clm[14][4] = -0.2917548685088514;
  clm[14][6] = 0.08948448095612106 - 0.1000111000111*logv;
  clm[14][8] = -0.4191005874775984 + 0.1118179659846327*logv;
  clm[14][10] = -0.5470367600420011 + 0.02917872533316407*logv;
  clm[14][12] = -0.670468777932017 - 0.02122541132988341*logv + 0.005001110062715124*logv2;
  clm[14][14] = -1.121914900357191 + 0.05563992717755679*logv - 0.005591518889563438*logv2;
  clm[14][16] = -1.863290198464164 + 0.05829132211850571*logv - 0.001459098208745744*logv2;
  clm[14][18] = -3.091320153387928 + 0.06247185416367579*logv + 0.001675254997123208*logv2 - 0.000166722172882907*logv3;
  clm[14][20] = -5.309288914701494 + 0.1212440670484604*logv - 0.003468642156240916*logv2 + 0.0001864046516260279*logv3;
  clm[14][22] = -9.285007340596687 + 0.1940815706474527*logv - 0.00309398820085199*logv2 + 0.00004864200562696251*logv3;
  clm[14][24] = -16.44921385210954 + 0.3157647205847183*logv - 0.002819440983415458*logv2 - 0.0000763125239785403*logv3 + 4.168516976565079e-6*logv4;
  clm[14][26] = -29.53186430873094 + 0.5476969423235019*logv - 0.006599177789782477*logv2 + 0.0001385146784948799*logv3 - 4.660633564076234e-6*logv4;
  clm[14][28] = -53.61120992757774 + 0.9537686079689937*logv - 0.01011377565844615*logv2 + 0.0001091150030537719*logv3 - 1.216185122374659e-6*logv4;
  clm[14][30] = -98.26909244745144 + 1.684686605591762*logv - 0.01607043314172434*logv2 + 0.00008132878718029996*logv3 + 2.419693963237519e-6*logv4 - 8.337959364824367e-8*logv5;
  clm[14][32] = -181.7343126088659 + 3.025949183091317*logv - 0.02830283263322823*logv2 + 0.0002406841961416637*logv3 - 4.035325816722677e-6*logv4 + 9.322301789838355e-8*logv5;
  clm[14][34] = -338.8348336079449 + 5.50736888971166*logv - 0.04900580111179493*logv2 + 0.0003515183235135234*logv3 - 2.87745982069401e-6*logv4 + 2.432640238116479e-8*logv5;
  clm[14][36] = -636.5355314577565 + 10.19825208976693*logv - 0.08624744472371682*logv2 + 0.0005431308573453331*logv3 - 1.654029582205626e-6*logv4 - 5.863376882341427e-8*logv5 + 1.389814146539896e-9*logv6;
  clm[14][38] = -1204.338220421203 + 19.34520069401661*logv - 0.1576358522788692*logv2 + 0.0009770256502296743*logv3 - 6.605229760367747e-6*logv4 + 9.215823428800077e-8*logv5 - 1.553889427728634e-9*logv6;
  
  clm[15][2] = -0.880952380952381/(1. - 5.*nu + 5.*nu2) + (4.5*nu)/(1. - 5.*nu + 5.*nu2) - (4.916666666666667*nu2)/(1. - 5.*nu + 5.*nu2) + (0.5833333333333333*nu3)/(1. - 5.*nu + 5.*nu2);
  clm[15][4] = -0.2479752507063431;
  clm[15][6] = 1.008778223534509 - 0.4000444000444*logv;
  clm[15][8] = 0.1189443611420747 + 0.352420066705781*logv;
  clm[15][10] = -0.02016851902017692 + 0.09920111039467873*logv;
  clm[15][12] = 1.309213587926891 - 0.599971598500349*logv + 0.08001776100344199*logv2;
  clm[15][14] = 2.208502648633813 + 0.1254496937815717*logv - 0.0704918370744608*logv2;
  clm[15][16] = 3.305438010536408 + 0.05677449072942549*logv - 0.01984242434578877*logv2;
  clm[15][18] = 7.444305496545044 - 0.9448565135395048*logv + 0.1592951033694906*logv2 - 0.01067021906450605*logv3;
  clm[15][20] = 14.98272155565497 - 0.7104328195204352*logv - 0.05970310894717618*logv2 + 0.009399954890160089*logv3;
  clm[15][22] = 28.46390694380826 - 1.263068707952848*logv - 0.02109847734692948*logv2 + 0.002645950247612489*logv3;
  clm[15][24] = 57.66267256223212 - 3.798181902220405*logv + 0.2925137452144732*logv2 - 0.02648061471242105*logv3 + 0.00106714034600066*logv4;
  clm[15][26] = 115.9750597030864 - 6.15625461023423*logv + 0.09049207822823792*logv2 + 0.0125765283947161*logv3 - 0.0009400998286196293*logv4;
  clm[15][28] = 231.0493861030565 - 11.94774113585906*logv + 0.2360064891802927*logv2 + 0.004112562597975413*logv3 - 0.0002646243948383675*logv4;
  clm[15][30] = 466.9561333507499 - 25.76079845161901*logv + 0.9870190049774344*logv2 - 0.05538277881152276*logv3 + 0.003172304560978949*logv4 - 0.00008538070389580152*logv5;
  clm[15][32] = 945.3390872660703 - 49.35114713098581*logv + 1.227604062148186*logv2 - 0.002918812647663349*logv3 - 0.001719366694336515*logv4 + 0.00007521633438439658*logv5;
  clm[15][34] = 1915.219013032986 - 98.69337392412439*logv + 2.490759613329867*logv2 - 0.02861489204338667*logv3 - 0.0005412283322622002*logv4 + 0.00002117230145404543*logv5;
  clm[15][36] = 3903.680910271429 - 202.6421809991913*logv + 5.852115730285185*logv2 - 0.1716247053038245*logv3 + 0.007433986910436765*logv4 - 0.0002957331200137012*logv5 + 5.692678744227415e-6*logv6;
  clm[15][38] = 7985.949671139861 - 407.8011252344942*logv + 10.4652771800987*logv2 - 0.1572427684693882*logv3 - 0.0008496272996041578*logv4 + 0.0001744946427427303*logv5 - 5.014978893724151e-6*logv6;
  clm[15][40] = 16397.77857210496 - 837.6577792910523*logv + 21.46935532696616*logv2 - 0.3437858898261829*logv3 + 0.002512370018473757*logv4 + 0.00005369834027909613*logv5 - 1.411643438790464e-6*logv6;
  
  clm[16][2] = -1.173611111111111/(1. - 4.*nu + 3.*nu2) + (5.152777777777778*nu)/(1. - 4.*nu + 3.*nu2) - (5.208333333333333*nu2)/(1. - 4.*nu + 3.*nu2) + (1.083333333333333*nu3)/(1. - 4.*nu + 3.*nu2);
  clm[16][4] = -0.5605554442947213;
  clm[16][6] = 1.870243941313376 - 0.9000999000999001*logv;
  clm[16][8] = -0.1730141439953753 + 1.056367243867244*logv;
  clm[16][10] = -1.859878727305043 + 0.5045558994101337*logv;
  clm[16][12] = 1.503840984630193 - 2.677759951137306*logv + 0.4050899150799251*logv2;
  clm[16][14] = 2.718574580050657 + 1.322714407624572*logv - 0.4754180253368565*logv2;
  clm[16][16] = -2.456893385077829 + 2.231466961843858*logv - 0.2270753573269383*logv2;
  clm[16][18] = 1.952876911228635 - 5.753095535781691*logv + 1.652634505144822*logv2 - 0.1215404640316392*logv3;
  clm[16][20] = 10.02035978518508 + 0.7057915188343507*logv - 1.120488821264864*logv2 + 0.1426412390371321*logv3;
  clm[16][22] = -1.762050398586885 + 5.484527878345057*logv - 1.255125073728925*logv2 + 0.06813016881504209*logv3;
  clm[16][24] = 7.815410253003804 - 16.66700799382937*logv + 5.063539407769625*logv2 - 0.6301129182514564*logv3 + 0.02734963988324348*logv4;
  clm[16][26] = 34.28904557513508 - 1.115171955670535*logv - 2.316747568367597*logv2 + 0.4937618283691619*logv3 - 0.03209784125186214*logv4;
  clm[16][28] = -13.6142375402125 + 13.60932209782288*logv - 4.218481781597067*logv2 + 0.4518437149585279*logv3 - 0.01533098953605218*logv4;
  clm[16][30] = -19.40411340424668 - 63.48635554377662*logv + 17.66093690842844*logv2 - 2.409947864150596*logv3 + 0.1720046921603512*logv4 - 0.004923481625335141*logv5;
  clm[16][32] = -29.09988473416285 - 5.739491493308159*logv - 6.104825873941489*logv2 + 1.468978412451647*logv3 - 0.1465676992850818*logv4 + 0.005778252740844714*logv5;
  clm[16][34] = -534.5054728156719 + 55.91993831271635*logv - 13.87109045200235*logv2 + 1.873940059795112*logv3 - 0.1186124897588802*logv4 + 0.002759884429966636*logv5;
  clm[16][36] = -1419.16509952046 - 222.7858684597667*logv + 71.28500893263598*logv2 - 9.546183621046467*logv3 + 0.7761094747366284*logv4 - 0.03640332363757386*logv5 + 0.000738604219851309*logv6;
  clm[16][38] = -3713.937491013746 + 137.4369005380124*logv - 21.92112478133935*logv2 + 4.921126060721236*logv3 - 0.543870511033832*logv4 + 0.03276843488137251*logv5 - 0.0008668341191310501*logv6;
  
  clm[17][2] = -1.023809523809524/(1. - 5.*nu + 5.*nu2) + (5.5*nu)/(1. - 5.*nu + 5.*nu2) - (6.916666666666667*nu2)/(1. - 5.*nu + 5.*nu2) + (1.583333333333333*nu3)/(1. - 5.*nu + 5.*nu2);
  clm[17][4] = -0.7228451986855348;
  clm[17][6] = 3.108265310049143 - 1.6001776001776*logv;
  clm[17][8] = 1.125129460265468 + 1.638277066848495*logv;
  clm[17][10] = -2.526252757555509 + 1.15668069533252*logv;
  clm[17][12] = 6.219386619889918 - 8.116424833167811*logv + 1.280284176055072*logv2;
  clm[17][14] = 17.32479759750863 + 1.417066308730376*logv - 1.310767132627812*logv2;
  clm[17][16] = 6.58361927465768 + 6.314101316069022*logv - 0.9254472696144746*logv2;
  clm[17][18] = 36.4981018188328 - 33.99059591405658*logv + 9.008258320123517*logv2 - 0.682894020128387*logv3;
  clm[17][20] = 130.4604702662439 - 16.64858243689519*logv - 3.708043209187966*logv2 + 0.6991534015600153*logv3;
  clm[17][22] = 146.7546663938798 + 7.719373056904413*logv - 6.869362061011746*logv2 + 0.4936266636608675*logv3;
  clm[17][24] = 389.7428528273225 - 208.8451110283208*logv + 51.36652711997797*logv2 - 6.146098693635042*logv3 + 0.273187928576169*logv4;
  clm[17][26] = 1180.214989111305 - 190.6606220249511*logv - 0.5956689896147897*logv2 + 3.350935932059407*logv3 - 0.2796924030660778*logv4;
  clm[17][28] = 1820.396832249497 - 156.9178049433848*logv - 24.35076718283097*logv2 + 4.633518198036662*logv3 - 0.1974725825101306*logv4;
  clm[17][30] = 4450.059044161869 - 1708.511285925715*logv + 347.6568075210313*logv2 - 42.92512594444447*logv3 + 2.995236298839372*logv4 - 0.08742984078930076*logv5;
  clm[17][32] = 12252.72850909396 - 2016.182676825665*logv + 97.70613400707267*logv2 + 10.43709214273402*logv3 - 1.889821468247219*logv4 + 0.08951150366523649*logv5;
  clm[17][34] = 22172.60482522548 - 2807.184689401682*logv + 18.31666743523555*logv2 + 24.58664678581961*logv3 - 2.24143675753904*logv4 + 0.06319824063638678*logv5;
  clm[17][36] = 53597.73189894204 - 16096.57378469419*logv + 2778.833617437393*logv2 - 319.0288045060139*logv3 + 24.43698651309845*logv4 - 1.13028872286328*logv5 + 0.02331721213635549*logv6;
  clm[17][38] = 141317.9432714796 - 22488.8882575078*logv + 1469.127660000635*logv2 + 3.958024165503231*logv3 - 9.302283655116637*logv4 + 0.7806049678505766*logv5 - 0.02387238385388777*logv6;
  clm[17][40] = 285451.8256225642 - 39012.51261392419*logv + 1817.663288372833*logv2 + 75.11080082598259*logv3 - 15.23717856304559*logv4 + 0.8414567538105182*logv5 - 0.01685473483949665*logv6;
  
  clm[18][2] = -1.284722222222222/(1. - 4.*nu + 3.*nu2) + (5.819444444444444*nu)/(1. - 4.*nu + 3.*nu2) - (6.319444444444444*nu2)/(1. - 4.*nu + 3.*nu2) + (1.527777777777778*nu3)/(1. - 4.*nu + 3.*nu2);
  clm[18][4] = -1.097394068633346;
  clm[18][6] = 4.423065177222889 - 2.5002775002775*logv;
  clm[18][8] = 1.167482384117555 + 3.212162066328733*logv;
  clm[18][10] = -8.039748356276659 + 2.743789698741937*logv;
  clm[18][12] = 6.951424091591444 - 18.73137156708344*logv + 3.125693789196953*logv2;
  clm[18][14] = 31.02879410169328 + 6.937977188779185*logv - 4.015648270843307*logv2;
  clm[18][16] = -27.10883401396521 + 28.52133730804786*logv - 3.430117824628823*logv2;
  clm[18][18] = 2.622411276949867 - 105.7532025288746*logv + 33.00847952485751*logv2 - 2.605033951295422*logv3;
  clm[18][20] = 241.7990501442648 - 16.60187880615444*logv - 20.99606069957623*logv2 + 3.34674500687259*logv3;
  clm[18][22] = -124.0502254255965 + 189.2030790361593*logv - 46.18146644521696*logv2 + 2.858748806740083*logv3;
  clm[18][24] = -48.32695128923596 - 816.4905248777904*logv + 272.1178138086722*logv2 - 35.50406185945683*logv3 + 1.628326943970734*logv4;
  clm[18][26] = 2099.913275702211 - 243.4993155127599*logv - 93.29105600836908*logv2 + 27.76865470488834*logv3 - 2.091947809962401*logv4;
  clm[18][28] = -1371.373890322351 + 1521.869818857354*logv - 420.6264182060486*logv2 + 47.26133202741203*logv3 - 1.786916330109346*logv4;
  clm[18][30] = -1568.19436585955 - 8262.723823348972*logv + 2581.107798722944*logv2 - 367.9266147267831*logv3 + 27.18927027754992*logv4 - 0.8142538442211296*logv5;
  clm[18][32] = 19603.05596289914 - 2069.408115549888*logv - 649.1950642255883*logv2 + 204.3148556855846*logv3 - 23.77679519809236*logv4 + 1.046090008200757*logv5;
  clm[18][34] = -24355.73217397576 + 16263.06557219583*logv - 4112.335345673659*logv2 + 530.9113134090461*logv3 - 35.02503541000433*logv4 + 0.8935573390101678*logv5;
  clm[18][36] = -43922.23705287004 - 94294.37745984936*logv + 28301.43888427573*logv2 - 3983.493814032424*logv3 + 333.5331679016582*logv4 - 16.09480572503187*logv5 + 0.3393100943700918*logv6;
  clm[18][38] = 165695.2655007736 - 11439.47954899095*logv - 6978.304541533387*logv2 + 1851.086062877051*logv3 - 226.5211258195106*logv4 + 15.09980327045615*logv5 - 0.4359192184615763*logv6;
  
  clm[19][2] = -1.261904761904762/(1. - 5.*nu + 5.*nu2) + (7.166666666666667*nu)/(1. - 5.*nu + 5.*nu2) - (10.25*nu2)/(1. - 5.*nu + 5.*nu2) + (3.25*nu3)/(1. - 5.*nu + 5.*nu2);
  clm[19][4] = -1.554311118386749;
  clm[19][6] = 5.62108684125061 - 3.6003996003996*logv;
  clm[19][8] = 2.586765673643722 + 4.543361400504258*logv;
  clm[19][10] = -14.71987171778136 + 5.596141129536306*logv;
  clm[19][12] = 6.313182769619637 - 36.14781587942923*logv + 6.481438641278801*logv2;
  clm[19][14] = 63.84368318281148 + 10.76308190957363*logv - 8.178958285423249*logv2;
  clm[19][16] = -79.38448856754049 + 77.72597711240944*logv - 10.07417214328114*logv2;
  clm[19][18] = -79.55126473742501 - 274.7070443739797*logv + 93.71385238874981*logv2 - 7.778589698024908*logv3;
  clm[19][20] = 642.0084467426511 - 65.8978622922516*logv - 55.51735879938505*logv2 + 9.815839380840956*logv3;
  clm[19][22] = -672.8817769180182 + 772.653377984057*logv - 184.4386315898924*logv2 + 12.09034845300874*logv3;
  clm[19][24] = -1286.093641566693 - 2945.655689664447*logv + 1074.695219959818*logv2 - 146.8415942486533*logv3 + 7.001507810110332*logv4;
  clm[19][26] = 8083.122884630402 - 948.8640219302947*logv - 336.2451704476271*logv2 + 110.0030326283783*logv3 - 8.835236046091609*logv4;
  clm[19][28] = -9081.412142283361 + 9525.667068906208*logv - 2464.046651280587*logv2 + 274.7764664871577*logv3 - 10.88252143472615*logv4;
  clm[19][30] = -19062.57529419995 - 42911.95215755722*logv + 14418.61603258085*logv2 - 2137.943046511198*logv3 + 163.1107779448557*logv4 - 5.041645184343184*logv5;
  clm[19][32] = 128836.7968306846 - 10867.89812129537*logv - 3680.851610291534*logv2 + 1141.116112064817*logv3 - 138.0553782616331*logv4 + 6.362076065956875*logv5;
  clm[19][34] = -135552.6890371023 + 143842.5429340854*logv - 35817.45731464893*logv2 + 4481.146483210583*logv3 - 295.4145944479028*logv4 + 7.83628516498562*logv5;
  clm[19][36] = -272813.488447405 - 731373.9813019432*logv + 226668.9273396297*logv2 - 32873.06408530208*logv3 + 2824.509750965002*logv4 - 139.7311135823901*logv5 + 3.025322884510962*logv6;
  clm[19][38] = 2.427901541557683e6 - 146358.7112684716*logv - 59688.29008420252*logv2 + 15255.57033671866*logv3 - 1863.531166887379*logv4 + 127.5240208569374*logv5 - 3.817669354263832*logv6;
  clm[19][40] = -1.875720037781619e6 + 2.402167477600229e6*logv - 582484.3379698158*logv2 + 73494.40743192533*logv3 - 5617.693317383513*logv4 + 247.3495543605843*logv5 - 4.702292996105257*logv6;
  
  clm[20][2] = -0.865546218487395/(1. - 4.*nu + 3.*nu2) + (3.526610644257703*nu)/(1. - 4.*nu + 3.*nu2) - (2.917366946778711*nu2)/(1. - 4.*nu + 3.*nu2) + (0.319327731092437*nu3)/(1. - 4.*nu + 3.*nu2);
  clm[20][4] = -0.1508235111143767;
  clm[20][6] = 0.1646857992777582 - 0.07355557607658448*logv;
  clm[20][8] = -0.1146151612764767 + 0.06366575072174959*logv;
  clm[20][10] = -0.08077441234858914 + 0.01109391024591112*logv;
  clm[20][12] = 0.03076164607686423 - 0.01876725343767367*logv + 0.002705211385979104*logv2;
  clm[20][14] = 0.1220226348136702 + 0.01418966441494504*logv - 0.002341485485343258*logv2;
  clm[20][16] = 0.3657824651154519 + 0.006944942014027768*logv - 0.0004080094795399577*logv2;
  clm[20][18] = 0.9329610214277329 - 0.00474750392129144*logv + 0.0009349262386680864*logv2 - 0.0000663277939682095*logv3;
  clm[20][20] = 2.118621613645649 - 0.007010548868877552*logv - 0.0007336707010837378*logv2 + 0.00005740977124979477*logv3;
  clm[20][22] = 4.594625401072714 - 0.02615839112912772*logv - 0.0002923273506791129*logv2 + 0.00001000379077075633*logv3;
  clm[20][24] = 9.706347949454449 - 0.06944429002147967*logv + 0.0002881244565539782*logv2 - 0.00002892289615518018*logv3 + 1.219694773805164e-6*logv4;
  clm[20][26] = 20.14054003763214 - 0.1561547078556184*logv + 0.0001664083592034193*logv2 + 0.00002318170013297335*logv3 - 1.055702199175898e-6*logv4;
  clm[20][28] = 41.30448293064695 - 0.3402236501026327*logv + 0.0009312380363316936*logv2 + 8.07235905099173e-6*logv3 - 1.839586482731503e-7*logv4;
  clm[20][30] = 84.03043794259975 - 0.720395092040001*logv + 0.00259704571272463*logv2 - 0.00001039051084383189*logv3 + 6.421912932661668e-7*logv4 - 1.794307034496764e-8*logv5;
  clm[20][32] = 169.9834202684805 - 1.495562634305508*logv + 0.005744193462736853*logv2 - 1.36873828449541e-6*logv3 - 5.21782598169245e-7*logv4 + 1.553055668514006e-8*logv5;
  clm[20][34] = 342.4694540049825 - 3.069204284903071*logv + 0.01259273166794561*logv2 - 0.00002199531668875542*logv3 - 1.650822972313233e-7*logv4 + 2.70623686960027e-9*logv5;
  clm[20][36] = 687.9921921065454 - 6.24722669590638*logv + 0.02673709152587413*logv2 - 0.0000650960508574273*logv3 + 2.622140934002791e-7*logv4 - 1.107044541157542e-8*logv5 + 2.199688126344624e-10*logv6;
  clm[20][38] = 1379.28005637384 - 12.6453297853608*logv + 0.05552113451388867*logv2 - 0.0001405669694692686*logv3 - 3.332758401424974e-8*logv4 + 9.080867923423045e-9*logv5 - 1.903931739609212e-10*logv6;
  
  clm[21][2] = 1.122582366279845/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - (8.235894357743097*nu)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) + (18.26890756302521*nu2)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - (12.68767507002801*nu3)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) + (2.184873949579832*nu4)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
  clm[21][4] = -0.351319484450667;
  clm[21][6] = 0.6051878143813924 - 0.2942223043063379*logv;
  clm[21][8] = -0.4260472172087262 + 0.3302887705805175*logv;
  clm[21][10] = -0.7317576096103622 + 0.1033660282627899*logv;
  clm[21][12] = -0.3216279886799198 - 0.284518866888056*logv + 0.04328338217566566*logv2;
  clm[21][14] = -0.7614428053877368 + 0.244861717650568*logv - 0.04858916158335363*logv2;
  clm[21][16] = -1.710757067281136 + 0.2527005709992152*logv - 0.01520629551123605*logv2;
  clm[21][18] = -2.393904400030632 - 0.05869646614192913*logv + 0.05751722117650645*logv2 - 0.004244978813965408*logv3;
  clm[21][20] = -4.027167600748205 + 0.3691863313683247*logv - 0.05360301427623388*logv2 + 0.004765338361789099*logv3;
  clm[21][22] = -7.462463859414905 + 0.6124770461029678*logv - 0.04267720002219957*logv2 + 0.001491343768426331*logv3;
  clm[21][24] = -12.84141989575447 + 0.5858690456748547*logv + 0.03685772810987156*logv2 - 0.007176919950664992*logv3 + 0.0003122418620941219*logv4;
  clm[21][26] = -22.95591997353326 + 1.414855635890697*logv - 0.0820264508692257*logv2 + 0.006981320483659549*logv3 - 0.0003505172084012445*logv4;
  clm[21][28] = -42.40967077658477 + 2.477550670868785*logv - 0.108147903105183*logv2 + 0.004725144291060049*logv3 - 0.0001096966500148232*logv4;
  clm[21][30] = -78.18131317003929 + 3.8737355915617*logv - 0.05586590932106153*logv2 - 0.006938480121676191*logv3 + 0.0006408816518758884*logv4 - 0.00001837370403324687*logv5;
  clm[21][32] = -146.7794290151865 + 7.436487293603099*logv - 0.2550111852672496*logv2 + 0.01138668689108439*logv3 - 0.0006403434744540262*logv4 + 0.00002062599615096781*logv5;
  clm[21][34] = -280.4304967830408 + 13.54177915878974*logv - 0.4134216150015225*logv2 + 0.01257157470365798*logv3 - 0.0003872524942897806*logv4 + 6.455040228409431e-6*logv5;
  clm[21][36] = -541.3226461134253 + 24.35594849464146*logv - 0.5655329130763367*logv2 + 0.001039583069926688*logv3 + 0.0007957194551910756*logv4 - 0.00004436053365415123*logv5 + 9.00992256550758e-7*logv6;
  clm[21][38] = -1059.625759302543 + 46.34509633657497*logv - 1.218974266496503*logv2 + 0.03111122611640746*logv3 - 0.001129268869655173*logv4 + 0.00004514381678205216*logv5 - 1.011438019358567e-6*logv6;
  
  clm[22][2] = -0.9327731092436975/(1. - 4.*nu + 3.*nu2) + (3.929971988795518*nu)/(1. - 4.*nu + 3.*nu2) - (3.589635854341737*nu2)/(1. - 4.*nu + 3.*nu2) + (0.5882352941176471*nu3)/(1. - 4.*nu + 3.*nu2);
  clm[22][4] = -0.3718741604762886;
  clm[22][6] = 1.515267376087476 - 0.6620001846892603*logv;
  clm[22][8] = 0.2389058102366552 + 0.6174959705925033*logv;
  clm[22][10] = -0.409515741515465 + 0.2461807629164667*logv;
  clm[22][12] = 2.138249656492113 - 1.542056545436917*logv + 0.2191221222643074*logv2;
  clm[22][14] = 3.882869644142317 + 0.3445616889124937*logv - 0.2043912232885556*logv2;
  clm[22][16] = 4.165921924428288 + 0.4715208010900693*logv - 0.08148585525882198*logv2;
  clm[22][18] = 11.33111399454784 - 3.244787204260044*logv + 0.6888131146343659*logv2 - 0.04835296180282472*logv3;
  clm[22][20] = 24.72860928391453 - 1.754680073876463*logv - 0.2804494498572441*logv2 + 0.04510234252196256*logv3;
  clm[22][22] = 42.61620399211888 - 2.160568154980297*logv - 0.2224128990249436*logv2 + 0.0179812170769675*logv3;
  clm[22][24] = 90.51336407767374 - 12.72296217032742*logv + 1.824745079671782*logv2 - 0.1913633717736554*logv3 + 0.008002417410935679*logv4;
  clm[22][26] = 189.344126829522 - 16.34037667206031*logv + 0.1753063156080104*logv2 + 0.09860469555917146*logv3 - 0.007464439769864373*logv4;
  clm[22][28] = 368.0909833342251 - 29.09978281078297*logv + 0.4634424072065594*logv2 + 0.06371804061231687*logv3 - 0.002975892256472541*logv4;
  clm[22][30] = 762.2208703333631 - 79.41901262758027*logv + 6.705519812669277*logv2 - 0.6003676252252755*logv3 + 0.0381855951414314*logv4 - 0.001059520360799994*logv5;
  clm[22][32] = 1579.682178504302 - 136.3214383467089*logv + 4.924385766738155*logv2 + 0.08068816078444449*logv3 - 0.0223960502289371*logv4 + 0.0009882921012504149*logv5;
  clm[22][34] = 3201.80281685076 - 267.0780127886743*logv + 9.663483991903678*logv2 - 0.03480534343304703*logv3 - 0.01296807958447866*logv4 + 0.0003940082446800324*logv5;
  clm[22][36] = 6646.515997217026 - 605.153280224857*logv + 35.75833381537822*logv2 - 2.225156904724671*logv3 + 0.1373853697973921*logv4 - 0.005918353599769061*logv5 + 0.0001169004457552713*logv6;
  clm[22][38] = 13844.47114405644 - 1162.393543644032*logv + 47.69703212341268*logv2 - 0.8507706953905084*logv3 - 0.03805743432025452*logv4 + 0.003769828739531367*logv5 - 0.0001090415922591186*logv6;
  
  clm[23][2] = 1.184207016139789/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - (8.790516206482593*nu)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) + (19.9327731092437*nu2)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - (14.53641456582633*nu3)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) + (2.739495798319328*nu4)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
  clm[23][4] = -0.64737468966706;
  clm[23][6] = 2.376594118290692 - 1.176889217225352*logv;
  clm[23][8] = 0.03483171517219779 + 1.393680468257526*logv;
  clm[23][10] = -2.549644814931386 + 0.7618882917737712*logv;
  clm[23][12] = 2.615179302641318 - 4.50033380917985*logv + 0.6925341148106505*logv2;
  clm[23][14] = 5.85160030915473 + 1.976120998160804*logv - 0.8201037576749307*logv2;
  clm[23][16] = -3.279694702450431 + 4.103352460539171*logv - 0.4483290576593969*logv2;
  clm[23][18] = 5.800055775131681 - 12.81546346640257*logv + 3.650521829963816*logv2 - 0.2716786440937861*logv3;
  clm[23][20] = 28.35800120088367 - 0.2084397571230767*logv - 2.349797645702164*logv2 + 0.3217237564712062*logv3;
  clm[23][22] = 5.083486956949533 + 11.8860189630802*logv - 3.063475250293594*logv2 + 0.1758778779093824*logv3;
  clm[23][24] = 37.9623606825984 - 50.12342562588101*logv + 14.7219733711686*logv2 - 1.825294955655846*logv3 + 0.07993391669609521*logv4;
  clm[23][26] = 150.7240012967134 - 13.55421130032578*logv - 5.525043743584685*logv2 + 1.387457272518919*logv3 - 0.09465830497904939*logv4;
  clm[23][28] = 70.65308894510197 + 30.50229728484663*logv - 12.65636933187836*logv2 + 1.456343471532978*logv3 - 0.05174719451500726*logv4;
  clm[23][30] = 221.936250976296 - 260.7704727268617*logv + 68.11192171846053*logv2 - 9.16147443621748*logv3 + 0.6527331584317712*logv4 - 0.01881467293004479*logv5;
  clm[23][32] = 750.5050549154053 - 104.5476013567672*logv - 15.1129699016347*logv2 + 5.056961211443847*logv3 - 0.5452225796024562*logv4 + 0.02228046769013441*logv5;
  clm[23][34] = 144.4914143766073 + 76.31729461279377*logv - 49.3909227016876*logv2 + 7.554693984493226*logv3 - 0.5033839439941136*logv4 + 0.01218014304927498*logv5;
  clm[23][36] = 250.7998386331813 - 1484.699683762236*logv + 366.7497037589568*logv2 - 47.84743288790288*logv3 + 3.859217894464863*logv4 - 0.1808699437151728*logv5 + 0.003690464282831905*logv6;
  clm[23][38] = 1144.587140358962 - 568.2850209976019*logv - 42.75659064071929*logv2 + 20.44418647408209*logv3 - 2.536314544481956*logv4 + 0.1605804805626198*logv5 - 0.004370273696542837*logv6;
  
  clm[24][2] = -1.067226890756303/(1. - 4.*nu + 3.*nu2) + (4.736694677871148*nu)/(1. - 4.*nu + 3.*nu2) - (4.934173669467787*nu2)/(1. - 4.*nu + 3.*nu2) + (1.126050420168067*nu3)/(1. - 4.*nu + 3.*nu2);
  clm[24][4] = -0.8269193364414115;
  clm[24][6] = 3.454953973732596 - 1.838889401914612*logv;
  clm[24][8] = 1.244599416341227 + 1.962512218850048*logv;
  clm[24][10] = -3.612173742586308 + 1.520613204020375*logv;
  clm[24][12] = 6.742794169163535 - 10.51183737150302*logv + 1.69075711623694*logv2;
  clm[24][14] = 21.12760355640137 + 2.149445448771333*logv - 1.804421460185642*logv2;
  clm[24][16] = 2.611325742762683 + 10.08118096539906*logv - 1.398119752642245*logv2;
  clm[24][18] = 36.3517946288593 - 48.4707004497653*logv + 13.48861831974738*logv2 - 1.036371780753273*logv3;
  clm[24][20] = 163.6610433329257 - 20.86416903191162*logv - 6.056907775782471*logv2 + 1.106043833240888*logv3;
  clm[24][22] = 139.7261526797444 + 28.16680507719338*logv - 12.43086837569363*logv2 + 0.8569958652471009*logv3;
  clm[24][24] = 413.1898311939644 - 320.7987728292861*logv + 86.37863713410732*logv2 - 10.61173023623844*logv3 + 0.4764432710176421*logv4;
  clm[24][26] = 1518.159404179245 - 254.7553849998199*logv - 6.582851697449791*logv2 + 6.213925615349316*logv3 - 0.5084730707499205*logv4;
  clm[24][28] = 1937.434951244078 - 64.59883531946086*logv - 63.36096883016284*logv2 + 9.55771858933207*logv3 - 0.3939801535218842*logv4;
  clm[24][30] = 5032.244022274792 - 2800.255396265202*logv + 648.5823358884759*logv2 - 83.87662026365699*logv3 + 5.955902898311719*logv4 - 0.1752252963375746*logv5;
  clm[24][32] = 16184.74536176648 - 2815.108723993503*logv + 106.193436089763*logv2 + 25.4853263650941*logv3 - 4.006567658459951*logv4 + 0.1870051481922015*logv5;
  clm[24][34] = 25206.94712416333 - 2645.221641501253*logv - 223.3215766212016*logv2 + 66.18428364075122*logv3 - 5.284863849001053*logv4 + 0.1448971857752169*logv5;
  clm[24][36] = 63634.65820795422 - 27882.58086391009*logv + 5684.054280937648*logv2 - 699.9384358750966*logv3 + 55.21563987564169*logv4 - 2.586712846164492*logv5 + 0.05370332339708555*logv6;
  clm[24][38] = 189918.9170462098 - 32721.06590084349*logv + 1983.664830327566*logv2 + 78.62968233714428*logv3 - 24.17776173814722*logv4 + 1.896430026686453*logv5 - 0.0573136308523518*logv6;
  
  clm[25][2] = 1.286914765906363/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - (9.714885954381753*nu)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) + (22.70588235294118*nu2)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - (17.61764705882353*nu3)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) + (3.663865546218487*nu4)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
  clm[25][4] = -1.140326502069253;
  clm[25][6] = 4.638611984185046 - 2.648000738757041*logv;
  clm[25][8] = 1.335352909623908 + 3.407751250837393*logv;
  clm[25][10] = -8.746419496025021 + 3.019585419903615*logv;
  clm[25][12] = 7.493107663674511 - 20.90623616274363*logv + 3.505953956228918*logv2;
  clm[25][14] = 35.49881475159081 + 7.561292734918962*logv - 4.511863914858824*logv2;
  clm[25][16] = -31.01700952327795 + 32.99377532581307*logv - 3.997932211322382*logv2;
  clm[25][18] = 1.35768889941637 - 124.6486300755616*logv + 39.09696876620986*logv2 - 3.094589555380782*logv3;
  clm[25][20] = 292.0224392766461 - 22.11447849184849*logv - 24.7039945644814*logv2 + 3.982472993239134*logv3;
  clm[25][22] = -146.487087638166 + 231.4565117349629*logv - 56.70299740221018*logv2 + 3.528842483027413*logv3;
  clm[25][24] = -62.85557979033369 - 1017.778838432662*logv + 340.9788614491865*logv2 - 44.58710098918076*logv3 + 2.048618857199534*logv4;
  clm[25][26] = 2724.534975230418 - 330.7812965067045*logv - 113.7451433197564*logv2 + 34.77428251193308*logv3 - 2.636397857044298*logv4;
  clm[25][28] = -1568.43735420117 + 1962.213715891363*logv - 546.5497896760764*logv2 + 61.54150029570663*logv3 - 2.336094375503455*logv4;
  clm[25][30] = -1612.1672192226 - 10904.24701463748*logv + 3421.626389654395*logv2 - 489.0880119794015*logv3 + 36.18797612545074*logv4 - 1.084948849459194*logv5;
  clm[25][32] = 28383.14869696013 - 3166.794222809963*logv - 817.4461481414688*logv2 + 268.8753055405059*logv3 - 31.60598497767222*logv4 + 1.396236694622156*logv5;
  clm[25][34] = -26211.90727983269 + 21726.99308853211*logv - 5644.821433891322*logv2 + 731.7732875275*logv3 - 48.34795277823778*logv4 + 1.237195926427864*logv5;
  clm[25][36] = -41718.27851668892 - 132262.5498072458*logv + 39661.71814973421*logv2 - 5603.254099467025*logv3 + 470.0350558159861*logv4 - 22.69828269474664*logv5 + 0.4788242258135914*logv6;
  clm[25][38] = 301034.5415156125 - 25831.38168556064*logv - 9004.47555232293*logv2 + 2558.143892327178*logv3 - 317.4850796258341*logv4 + 21.28536529323461*logv5 - 0.6162059664731932*logv6;
  
  clm[26][2] = -1.26890756302521/(1. - 4.*nu + 3.*nu2) + (5.946778711484594*nu)/(1. - 4.*nu + 3.*nu2) - (6.950980392156863*nu2)/(1. - 4.*nu + 3.*nu2) + (1.932773109243697*nu3)/(1. - 4.*nu + 3.*nu2);
  clm[26][4] = -1.541846793492343;
  clm[26][6] = 5.663709257777565 - 3.60422322775264*logv;
  clm[26][8] = 2.54758143364005 + 4.573426112526459*logv;
  clm[26][10] = -14.753691230417 + 5.557160026741031*logv;
  clm[26][12] = 6.579814286816788 - 36.38879319711669*logv + 6.495212537735828*logv2;
  clm[26][14] = 64.26967218822669 + 11.08940690618708*logv - 8.241824312589159*logv2;
  clm[26][16] = -80.33253536689753 + 77.80740204738514*logv - 10.01462262435925*logv2;
  clm[26][18] = -79.84695569995812 - 277.6154981592823*logv + 94.36633828972789*logv2 - 7.803398632565879*logv3;
  clm[26][20] = 643.4947569445058 - 64.97672508780842*logv - 56.51578082196189*logv2 + 9.901791542163426*logv3;
  clm[26][22] = -696.8935584805936 + 777.2023904145858*logv - 184.6068854925862*logv2 + 12.03164515996424*logv3;
  clm[26][24] = -1344.740560326935 - 2970.238327371462*logv + 1085.459565975261*logv2 - 147.9605836522042*logv3 + 7.031297651726782*logv4;
  clm[26][26] = 7939.986495987045 - 936.0675776024468*logv - 345.1781162301413*logv2 + 111.7876417243546*logv3 - 8.922066768157513*logv4;
  clm[26][28] = -9824.44705569599 + 9640.38542990259*logv - 2476.186976404769*logv2 + 275.1177453613578*logv3 - 10.84118373840518*logv4;
  clm[26][30] = -21209.45476310728 - 42900.07559088104*logv + 14557.47749643779*logv2 - 2160.411829968671*logv3 + 164.4865822622791*logv4 - 5.068473263519251*logv5;
  clm[26][32] = 122232.8799556034 - 10241.75573213954*logv - 3793.500562012632*logv2 + 1164.614751610858*logv3 - 140.273472743515*logv4 + 6.431424057070646*logv5;
  clm[26][34] = -157836.5567483436 + 147008.9309078772*logv - 36127.6492695139*logv2 + 4503.50270267017*logv3 - 295.9493912318696*logv4 + 7.81480924925883*logv5;
  clm[26][36] = -340371.1912133246 - 720832.9244287086*logv + 227917.817720976*logv2 - 33232.52979339402*logv3 + 2856.395682928409*logv4 - 141.0350003744899*logv5 + 3.044651510936552*logv6;
  clm[26][38] = 2.197980950148844e6 - 122753.3396536351*logv - 61901.24536404951*logv2 + 15593.40564981101*logv3 - 1900.38621828469*logv4 + 129.6223142706865*logv5 - 3.863381329003523*logv6;
  
  clm[27][2] = 1.097697368421053/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - (6.932620614035088*nu)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) + (12.98914473684211*nu2)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - (7.589364035087719*nu3)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) + (1.18640350877193*nu4)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
  clm[27][4] = -0.2684213351704371;
  clm[27][6] = -0.04026002438780446 - 0.0564155486756106*logv;
  clm[27][8] = -0.3861652781633408 + 0.06192719931924755*logv;
  clm[27][10] = -0.4951791133124597 + 0.01514313689988018*logv;
  clm[27][12] = -0.6774255977204517 - 0.001652242589574806*logv + 0.001591357066185094*logv2;
  clm[27][14] = -1.090375566041843 + 0.02609257894448544*logv - 0.001746828463769625*logv2;
  clm[27][16] = -1.779434173367453 + 0.02898896159304093*logv - 0.0004271541884383126*logv2;
  clm[27][18] = -2.963169996457922 + 0.03774630622071774*logv + 0.0001572802465303912*logv2 - 0.00002992576067588068*logv3;
  clm[27][20] = -5.056398961141743 + 0.06371971141246999*logv - 0.0008575003133942658*logv2 + 0.00003284942874191245*logv3;
  clm[27][22] = -8.772561781866856 + 0.1024994421960989*logv - 0.0008474213928105546*logv2 + 8.032712636610854e-6*logv3;
  clm[27][24] = -15.42564526731568 + 0.1697179324318313*logv - 0.001043755504823867*logv2 - 5.038931629240938e-6*logv3 + 4.22069552016205e-7*logv4;
  clm[27][26] = -27.43979368225512 + 0.2899278030292494*logv - 0.00186805466076116*logv2 + 0.00001841003048699726*logv3 - 4.633046365388409e-7*logv4;
  clm[27][28] = -49.28847708187 + 0.5022380302613717*logv - 0.002952913179784615*logv2 + 0.00001649456559908705*logv3 - 1.13292472686978e-7*logv4;
  clm[27][30] = -89.27899508322265 + 0.882270588499165*logv - 0.004857110777861639*logv2 + 0.00001908866467419949*logv3 + 1.004222076866794e-7*logv4 - 4.762257071252675e-9*logv5;
  clm[27][32] = -162.9056481575626 + 1.568644108180782*logv - 0.008315581181605996*logv2 + 0.00003661683386755949*logv3 - 2.918744550333514e-7*logv4 + 5.22751705485861e-9*logv5;
  clm[27][34] = -299.1791406910273 + 2.816250526821153*logv - 0.01437822971535366*logv2 + 0.00005672792682087644*logv3 - 2.405166473049253e-7*logv4 + 1.278291401490498e-9*logv5;
  clm[27][36] = -552.6241579851332 + 5.099193509179734*logv - 0.02523065910099707*logv2 + 0.00009259978807561437*logv3 - 2.595760009141831e-7*logv4 - 1.464275632922441e-9*logv5 + 4.477755760150434e-11*logv6;
  
  clm[28][2] = 0.8998538011695906/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - (6.432017543859649*nu)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) + (13.56688596491228*nu2)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - (8.349780701754386*nu3)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) + (1.119517543859649*nu4)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
  clm[28][4] = -0.2261796441029474;
  clm[28][6] = 0.5511242901214199 - 0.2256621947024424*logv;
  clm[28][8] = -0.108724909193424 + 0.203062983683265*logv;
  clm[28][10] = -0.1603499991475961 + 0.05104019488528845*logv;
  clm[28][12] = 0.3250799802431515 - 0.1871444601443253*logv + 0.02546171305896151*logv2;
  clm[28][14] = 0.5117442997880728 + 0.081024812723737*logv - 0.02291181928039592*logv2;
  clm[28][16] = 0.8626870722037884 + 0.05038370894564951*logv - 0.005758921197927283*logv2;
  clm[28][18] = 2.1111394887332 - 0.1482114362433531*logv + 0.02819886106767676*logv2 - 0.001915248683256364*logv3;
  clm[28][20] = 4.409840322979726 - 0.07243190621947299*logv - 0.01551591462434826*logv2 + 0.001723443807813292*logv3;
  clm[28][22] = 8.891308677765211 - 0.1755046655000746*logv - 0.007286912670623534*logv2 + 0.000433190265547563*logv3;
  clm[28][24] = 18.32246401659411 - 0.5533410533575773*logv + 0.02713907199771176*logv2 - 0.002653938402837774*logv3 + 0.0001080498053161485*logv4;
  clm[28][26] = 37.27492394959283 - 0.9919788763085359*logv + 0.001542142006400079*logv2 + 0.001646560053888139*logv3 - 0.00009722902802937046*logv4;
  clm[28][28] = 75.32727036174022 - 2.046365972068056*logv + 0.01719360226025485*logv2 + 0.0006686352905651749*logv3 - 0.00002443866651179922*logv4;
  clm[28][30] = 152.2819576067931 - 4.333276671846018*logv + 0.07472799086498939*logv2 - 0.00297315449162011*logv3 + 0.0001797815641122823*logv4 - 4.876551240962739e-6*logv5;
  clm[28][32] = 307.1320049770177 - 8.674606326283328*logv + 0.1090817833218316*logv2 + 0.0005161196664337046*logv3 - 0.0001199395500661346*logv4 + 4.388183170778605e-6*logv5;
  clm[28][34] = 618.5695620544179 - 17.57820331523107*logv + 0.2345119461559303*logv2 - 0.001063556866518657*logv3 - 0.0000445199736447001*logv4 + 1.102976624130739e-6*logv5;
  clm[28][36] = 1245.505282276923 - 35.666505079728*logv + 0.5163349212771814*logv2 - 0.006858861047099605*logv3 + 0.0002286582334073508*logv4 - 9.470579120924788e-6*logv5 + 1.834088759357618e-7*logv6;
  clm[28][38] = 2507.079643736949 - 71.80609074433956*logv + 1.00651172678062*logv2 - 0.007766990658408049*logv3 - 0.00007230303287080817*logv4 + 6.633904877174163e-6*logv5 - 1.650411741790371e-7*logv6;
  clm[29][2] = 1.129276315789474/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - (7.185252192982456*nu)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) + (13.65230263157895*nu2)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - (8.22094298245614*nu3)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) + (1.344298245614035*nu4)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
  
  clm[29][4] = -0.4196774909106648;
  clm[29][6] = 1.092660734256684 - 0.5077399380804954*logv;
  clm[29][8] = -0.3530977347355733 + 0.5733786866547173*logv;
  clm[29][10] = -1.011356410663756 + 0.2130870232487586*logv;
  clm[29][12] = 0.2795945021049974 - 0.8725937439180844*logv + 0.1288999223609926*logv2;
  clm[29][14] = 0.1352618160911201 + 0.5381728935164843*logv - 0.145563629429371*logv2;
  clm[29][16] = -1.511424810568664 + 0.6468821710760524*logv - 0.05409639599504089*logv2;
  clm[29][18] = -0.8230627422015665 - 0.9477509719676868*logv + 0.302206809693804*logv2 - 0.02181587953271702*logv3;
  clm[29][20] = -0.550065313500036 + 0.5613518567051008*logv - 0.227737601037403*logv2 + 0.02463615606441366*logv3;
  clm[29][22] = -3.630341576999292 + 1.281263058275853*logv - 0.198084150673681*logv2 + 0.009155633584300005*logv3;
  clm[29][24] = -4.117793015940577 - 1.052164968263592*logv + 0.4956717165937958*logv2 - 0.06480255594638107*logv3 + 0.002769198330778321*logv4;
  clm[29][26] = -5.332895990598469 + 1.392184533988381*logv - 0.3594849074257944*logv2 + 0.0539641689025325*logv3 - 0.003127190088671703*logv4;
  clm[29][28] = -15.02300260683365 + 3.156746673847175*logv - 0.4769197689786912*logv2 + 0.03925580238905036*logv3 - 0.001162170207295047*logv4;
  clm[29][30] = -24.72639442761109 - 0.7008356688777928*logv + 0.8412240840431906*logv2 - 0.1356068999633292*logv3 + 0.009959017154087578*logv4 - 0.0002812065178003992*logv5;
  clm[29][32] = -45.8680673578273 + 5.348336043580834*logv - 0.8183601049010202*logv2 + 0.1072157784055385*logv3 - 0.008807322039835296*logv4 + 0.0003175598603976218*logv5;
  clm[29][34] = -110.0455796481566 + 11.1886704601853*logv - 1.247087182154438*logv2 + 0.1099695668249651*logv3 - 0.005710364063361521*logv4 + 0.0001180160458181967*logv5;
  clm[29][36] = -232.2019871048044 + 7.037012898243219*logv + 1.444829533828377*logv2 - 0.2788930255904016*logv3 + 0.02486276160812705*logv4 - 0.001187331858235367*logv5 + 0.00002379662998930107*logv6;
  
  clm[30][2] = 0.9744152046783626/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - (7.103070175438596*nu)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) + (15.58004385964912*nu2)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - (10.58662280701754*nu3)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) + (1.790570175438596*nu4)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
  clm[30][4] = -0.4765205915006815;
  clm[30][6] = 1.948423309676332 - 0.9026487788097695*logv;
  clm[30][8] = 0.3778335967605529 + 0.8795546945565956*logv;
  clm[30][10] = -0.9066749911988272 + 0.4301307299957992*logv;
  clm[30][12] = 2.974128729686728 - 2.763166613591033*logv + 0.4073874089433841*logv2;
  clm[30][14] = 6.006986383195681 + 0.6376756576241908*logv - 0.3969644854689554*logv2;
  clm[30][16] = 4.529369465400148 + 1.297038122174402*logv - 0.1941284890796314*logv2;
  clm[30][18] = 15.4023429159389 - 7.217989931679307*logv + 1.700405845751938*logv2 - 0.1225759157284073*logv3;
  clm[30][20] = 37.91423656792668 - 3.291267968276383*logv - 0.7295218036272335*logv2 + 0.1194398360131337*logv3;
  clm[30][22] = 57.19376701702843 - 1.950058392358427*logv - 0.8014019016321848*logv2 + 0.05840994786663833*logv3;
  clm[30][24] = 129.2053829732866 - 30.7019551394593*logv + 5.808121681699477*logv2 - 0.6480197443694149*logv3 + 0.02766075016093401*logv4;
  clm[30][26] = 290.2321728118671 - 32.66324004057749*logv + 0.03215976559432336*logv2 + 0.3524076321375527*logv3 - 0.02695305552962357*logv4;
  clm[30][28] = 532.8771893910477 - 49.64984177763758*logv - 0.3253636816647978*logv2 + 0.3061239652649384*logv3 - 0.01318091702804085*logv4;
  clm[30][30] = 1141.927666428424 - 191.7528758091341*logv + 25.00896498993805*logv2 - 2.666732279653826*logv3 + 0.1770131268409914*logv4 - 0.004993588470745845*logv5;
  clm[30][32] = 2470.643572387351 - 284.887176084416*logv + 11.70669917359238*logv2 + 0.5754802999056636*logv3 - 0.1095171598382017*logv4 + 0.004865828531801324*logv5;
  clm[30][34] = 4898.923978564774 - 525.5329041759924*logv + 19.82587356455908*logv2 + 0.5329267557812712*logv3 - 0.08374770433402771*logv4 + 0.002379547731790794*logv5;
  clm[30][36] = 10421.09022166571 - 1446.217523180099*logv + 139.140773902715*logv2 - 12.12356685729754*logv3 + 0.8434520867623728*logv4 - 0.03751276465024688*logv5 + 0.0007512427558328802*logv6;
  clm[30][38] = 22350.38377543279 - 2537.268627176915*logv + 131.0852122608999*logv2 - 1.743249421054461*logv3 - 0.2952856379632813*logv4 + 0.02518556901982776*logv5 - 0.0007320223636880331*logv6;
  
  clm[31][2] = 1.192434210526316/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - (7.690515350877193*nu)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) + (14.97861842105263*nu2)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - (9.484100877192982*nu3)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) + (1.660087719298246*nu4)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
  clm[31][4] = -0.7220789990670206;
  clm[31][6] = 2.77743542727262 - 1.410388716890265*logv;
  clm[31][8] = 0.2354792835084879 + 1.681795756160267*logv;
  clm[31][10] = -3.231197515476268 + 1.018412072987542*logv;
  clm[31][12] = 3.577874284985517 - 6.369472310458006*logv + 0.9945981663656839*logv2;
  clm[31][14] = 9.197565237498466 + 2.591980246872178*logv - 1.185992879301186*logv2;
  clm[31][16] = -4.573907678257224 + 6.327932937314457*logv - 0.7181784484432271*logv2;
  clm[31][18] = 9.186630156314592 - 21.68504699988174*logv + 6.220999696050501*logv2 - 0.4675900105606357*logv3;
  clm[31][20] = 51.05732026410876 - 1.830340663310938*logv - 3.88990695818561*logv2 + 0.5575703250928633*logv3;
  clm[31][22] = 9.779882933919075 + 21.47115242939105*logv - 5.711102091968533*logv2 + 0.3376369267993614*logv3;
  clm[31][24] = 69.38467937262716 - 101.4310829222686*logv + 30.03246221851754*logv2 - 3.737663359158753*logv3 + 0.1648709187563301*logv4;
  clm[31][26] = 317.2182838926112 - 34.39095364268541*logv - 10.15163848483552*logv2 + 2.79819432212023*logv3 - 0.1965977238459528*logv4;
  clm[31][28] = 157.1584599266376 + 67.01880356798942*logv - 27.90451972812282*logv2 + 3.271999135567879*logv3 - 0.1190498279908309*logv4;
  clm[31][30] = 524.4163011993937 - 625.7812500416758*logv + 166.3713764303666*logv2 - 22.46251614449526*logv3 + 1.604546632870285*logv4 - 0.0465064167114519*logv5;
  clm[31][32] = 2059.402905851572 - 288.1758997929052*logv - 30.30693416211035*logv2 + 11.83752021798056*logv3 - 1.328455128468301*logv4 + 0.05545584229572801*logv5;
  clm[31][34] = 1026.302958589135 + 203.0894944083293*logv - 131.3670721334237*logv2 + 20.13976447833545*logv3 - 1.360686719806631*logv4 + 0.03358130682919895*logv5;
  clm[31][36] = 2861.496354608876 - 4300.123126981623*logv + 1068.079892530676*logv2 - 140.5183752736862*logv3 + 11.36047140429968*logv4 - 0.5334664743386564*logv5 + 0.01093202089880477*logv6;
  
  clm[32][2] = 1.098684210526316/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - (8.221491228070175*nu)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) + (18.93530701754386*nu2)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - (14.31469298245614*nu3)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) + (2.908991228070175*nu4)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
  clm[32][4] = -0.9061610303170206;
  clm[32][6] = 3.728710335572604 - 2.030959752321981*logv;
  clm[32][8] = 1.353048126726325 + 2.231383412090598*logv;
  clm[32][10] = -4.553227506130684 + 1.840376581696488*logv;
  clm[32][12] = 7.136820687766961 - 12.65776062543266*logv + 2.062398757775882*logv2;
  clm[32][14] = 24.64233266583133 + 2.838713060161271*logv - 2.265924950977449*logv2;
  clm[32][16] = -1.497762392798969 + 13.85516003644761*logv - 1.868865383270737*logv2;
  clm[32][18] = 34.98343517953965 - 62.80097180152478*logv + 18.01731482058861*logv2 - 1.396216290093889*logv3;
  clm[32][20] = 196.3890417878832 - 24.68543105127737*logv - 8.555836749349667*logv2 + 1.534000792405786*logv3;
  clm[32][22] = 125.5022100368004 + 52.78703620618591*logv - 18.74870164349015*logv2 + 1.265196791976887*logv3;
  clm[32][24] = 418.4725469252681 - 442.0913871982069*logv + 125.7553800492672*logv2 - 15.69317756698094*logv3 + 0.7089147726792502*logv4;
  clm[32][26] = 1870.651434753064 - 316.3049194859453*logv - 14.89110311968426*logv2 + 9.632853960591411*logv3 - 0.7788734673515447*logv4;
  clm[32][28] = 1952.502092643075 + 90.52942444745752*logv - 115.8343067527977*logv2 + 15.86028401482912*logv3 - 0.6423909408179859*logv4;
  clm[32][30] = 5404.121991444188 - 4087.755495054164*logv + 1023.147958605832*logv2 - 135.8480145911727*logv3 + 9.742958061550423*logv4 - 0.2879554742276088*logv5;
  clm[32][32] = 20530.20958650449 - 3617.139914733715*logv + 90.9136770956221*logv2 + 46.74840714512519*logv3 - 6.841044833046276*logv4 + 0.3163721328684912*logv5;
  clm[32][34] = 27100.22344633354 - 1711.999550933126*logv - 636.7816887447184*logv2 + 128.4780973062579*logv3 - 9.661249417349364*logv4 + 0.2609340292115162*logv5;
  clm[32][36] = 71763.79312103799 - 42945.75178089501*logv + 9664.054079031757*logv2 - 1236.960744325847*logv3 + 99.16845271777349*logv4 - 4.678463284658176*logv5 + 0.09747099643617717*logv6;
  clm[32][38] = 247239.7774001666 - 43587.31568091222*logv + 2276.787411696883*logv2 + 217.5784180316715*logv3 - 47.23585351926251*logv4 + 3.570876083682873*logv5 - 0.1070898447686947*logv6;
  
  clm[33][2] = 1.287171052631579/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - (8.448410087719298*nu)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) + (16.96809210526316*nu2)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - (11.37883771929825*nu3)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) + (2.133771929824561*nu4)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
  clm[33][4] = -1.175404252991305;
  clm[33][6] = 4.801397057954937 - 2.764361885104919*logv;
  clm[33][8] = 1.478942435238317 + 3.558206597505115*logv;
  clm[33][10] = -9.319543186615198 + 3.249242716559383*logv;
  clm[33][12] = 7.876862824582648 - 22.69320404847565*logv + 3.820848315910411*logv2;
  clm[33][14] = 39.23851838996583 + 8.037340555565651*logv - 4.918085348736*logv2;
  clm[33][16] = -34.41379864487549 + 36.83537410437634*logv - 4.491041360555763*logv2;
  clm[33][18] = -0.4028514552805957 - 141.0058999527918*logv + 44.3868184596104*logv2 - 3.520735817756687*logv3;
  clm[33][20] = 335.7165218872334 - 27.15066302187201*logv - 27.86893260242245*logv2 + 4.531789228579577*logv3;
  clm[33][22] = -169.0604182604 + 269.9063269538522*logv - 66.21774330808568*logv2 + 4.13828785385003*logv3;
  clm[33][24] = -86.75861587297136 - 1201.767317789548*logv + 404.0670527753738*logv2 - 52.89838927680245*logv3 + 2.433146975532571*logv4;
  clm[33][26] = 3290.365931064441 - 411.663654060078*logv - 131.9845427031609*logv2 + 41.12339033441915*logv3 - 3.131876353703602*logv4;
  clm[33][28] = -1810.595573344593 + 2390.210850948554*logv - 666.7839138324942*logv2 + 75.11907810814349*logv3 - 2.859931303193915*logv4;
  clm[33][30] = -1857.748874434925 - 13444.11674899475*logv + 4233.81942339611*logv2 - 605.9577574289288*logv3 + 44.84926212344715*logv4 - 1.34521975200411*logv5;
  clm[33][32] = 36603.16043162465 - 4217.209436541402*logv - 976.4142967906078*logv2 + 330.4431472233343*logv3 - 39.09280571334907*logv4 + 1.731527924207922*logv5;
  clm[33][34] = -29378.99501164079 + 27390.77955624367*logv - 7194.369772428028*logv2 + 933.1227957185185*logv3 - 61.66016601617108*logv4 + 1.58117701771354*logv5;
  clm[33][36] = -44072.37169435979 - 170539.7246016607*logv + 51216.03072766674*logv2 - 7249.954195669449*logv3 + 608.4865680290047*logv4 - 29.3801641573904*logv5 + 0.6197790349217422*logv6;
  
  clm[34][2] = 1.27266081871345/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - (9.787280701754386*nu)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) + (23.63267543859649*nu2)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - (19.53399122807018*nu3)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) + (4.474780701754386*nu4)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
  clm[34][4] = -1.533709250282138;
  clm[34][6] = 5.696023013029904 - 3.610595115239078*logv;
  clm[34][8] = 2.537604342425977 + 4.59506293540295*logv;
  clm[34][10] = -14.77794969285508 + 5.537603127265676*logv;
  clm[34][12] = 6.788522728391275 - 36.63682794725038*logv + 6.518198543094146*logv2;
  clm[34][14] = 64.87096758962416 + 11.29040938086284*logv - 8.295455894391014*logv2;
  clm[34][16] = -80.91544069742709 + 78.00512004803278*logv - 9.997021400719047*logv2;
  clm[34][18] = -79.78326372889528 - 280.9367396718201*logv + 95.15294311923451*logv2 - 7.844858606618066*logv3;
  clm[34][20] = 649.648007602446 - 65.15951210945412*logv - 57.3057058873432*logv2 + 9.983844176989805*logv3;
  clm[34][22] = -712.1950342951454 + 782.5342680005674*logv - 185.3192952511785*logv2 + 12.03173221212557*logv3;
  clm[34][24] = -1381.728091857961 - 3010.965124099805*logv + 1099.235108245833*logv2 - 149.4371282348339*logv3 + 7.081152041199107*logv4;
  clm[34][26] = 7936.178800222597 - 942.1373539860641*logv - 351.9225743392615*logv2 + 113.4074245090202*logv3 - 9.011904754186876*logv4;
  clm[34][28] = -10319.3594156685 + 9760.482381220705*logv - 2496.047390458268*logv2 + 276.5910081786453*logv3 - 10.86042838824131*logv4;
  clm[34][30] = -22625.58401119238 - 43367.72192283886*logv + 14776.83372771198*logv2 - 2190.946466053268*logv3 + 166.4075200848711*logv4 - 5.113434594043745*logv5;
  clm[34][32] = 119691.2153113268 - 10069.56038068149*logv - 3877.533716296967*logv2 + 1186.468996337164*logv3 - 142.4791517110795*logv4 + 6.507667856893392*logv5;
  clm[34][34] = -172583.1493681104 + 150059.8640829464*logv - 36571.38603887291*logv2 + 4546.272981154833*logv3 - 298.0044114721901*logv4 + 7.842521937597581*logv5;
  clm[34][36] = -384036.8514093059 - 723849.7467513309*logv + 231324.8448910075*logv2 - 33792.97535923978*logv3 + 2901.426365551152*logv4 - 142.925984508547*logv5 + 3.077090327891477*logv6;
  clm[34][38] = 2.090929497068402e6 - 111462.2686650418*logv - 63604.89997801167*logv2 + 15925.23264182338*logv3 - 1938.151976404402*logv4 + 131.8526007554969*logv5 - 3.916092295949607*logv6;
  
  clm[35][2] = 0.8936507936507937/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - (5.482539682539683*nu)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) + (9.711111111111111*nu2)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - (4.984126984126984*nu3)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) + (0.6349206349206349*nu4)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
  clm[35][4] = -0.1748231756470704;
  clm[35][6] = 0.02956057389859245 - 0.04466285539556953*logv;
  clm[35][8] = -0.1851774385668289 + 0.03991299617096135*logv;
  clm[35][10] = -0.1980665035872493 + 0.00780810221371936*logv;
  clm[35][12] = -0.2144330670952057 - 0.00378530344983222*logv + 0.0009973853260427773*logv2;
  clm[35][14] = -0.2864192102621554 + 0.01047344152055795*logv - 0.0008913141881937835*logv2;
  clm[35][16] = -0.3590464122400539 + 0.009277162395814435*logv - 0.0001743660700425871*logv2;
  clm[35][18] = -0.4088846507330296 + 0.009190968327656362*logv + 0.0001395791779723477*logv2 - 0.00001484869219690385*logv3;
  clm[35][20] = -0.3788153780233501 + 0.01352880167091797*logv - 0.0002830805441858906*logv2 + 0.00001326954556643947*logv3;
  clm[35][22] = -0.1042783406568641 + 0.01657906260049806*logv - 0.000216795938306954*logv2 + 2.595895524068605e-6*logv3;
  clm[35][24] = 0.7832966735810512 + 0.01872833674131438*logv - 0.0001935842743257717*logv2 - 2.8975343895137e-6*logv3 + 1.657962481009095e-7*logv4;
  clm[35][26] = 3.091521773383798 + 0.01773035165529451*logv - 0.0003212797034510843*logv2 + 4.946771311591826e-6*logv3 - 1.481639486997017e-7*logv4;
  clm[35][28] = 8.558754698668745 + 0.005613754735453107*logv - 0.0003828918039394679*logv2 + 3.370848548333141e-6*logv3 - 2.898502660337057e-8*logv4;
  clm[35][30] = 20.85751820993533 - 0.03392110111024586*logv - 0.000427783195627532*logv2 + 2.663140119074563e-6*logv3 + 4.150370906614798e-8*logv4 - 1.480986770811778e-9*logv5;
  clm[35][32] = 47.65342266148895 - 0.1370335597364572*logv - 0.000415316132336022*logv2 + 5.108791363713679e-6*logv3 - 6.341173573908688e-8*logv4 + 1.323485003122272e-9*logv5;
  clm[35][34] = 104.8043342753109 - 0.3818765618257362*logv - 0.0001474574389248753*logv2 + 5.896697173971054e-6*logv3 - 3.923767936844576e-8*logv4 + 2.589108103646152e-10*logv5;
  clm[35][36] = 224.888171171151 - 0.9333502428528438*logv + 0.0007333155678176945*logv2 + 6.495556633885453e-6*logv3 - 2.678698916690853e-8*logv4 - 4.52473834379934e-10*logv5 + 1.102418299791966e-11*logv6;
  
  clm[36][2] = -1.098925765592432/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (10.25485008818342*nu)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) - (32.88359788359788*nu2)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (42.33333333333333*nu3)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) - (19.87301587301587*nu4)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (2.804232804232804*nu5)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4);
  clm[36][4] = -0.3075612113861835;
  clm[36][6] = 0.315060944434847 - 0.1786514215822781*logv;
  clm[36][8] = -0.4190213897713871 + 0.1963246502364814*logv;
  clm[36][10] = -0.5891806155404224 + 0.05494624763770924*logv;
  clm[36][12] = -0.5203134603563333 - 0.09572678660656996*logv + 0.01595816521668444*logv2;
  clm[36][14] = -0.9536947261879051 + 0.1182011694960216*logv - 0.01753683892819547*logv2;
  clm[36][16] = -1.673374629632685 + 0.1173883843119538*logv - 0.004908112625544323*logv2;
  clm[36][18] = -2.625083644393215 + 0.06047369051836322*logv + 0.01207393190615118*logv2 - 0.0009503163006018461*logv3;
  clm[36][20] = -4.4350032256442 + 0.208944179123203*logv - 0.01442999439585555*logv2 + 0.001044327068193852*logv3;
  clm[36][22] = -7.696012982394446 + 0.3283565389823094*logv - 0.01156936012931545*logv2 + 0.0002922804326131403*logv3;
  clm[36][24] = -13.3279706138773 + 0.4696245036070181*logv - 0.001722678166985301*logv2 - 0.0009288087739710433*logv3 + 0.00004244383951383284*logv4;
  clm[36][26] = -23.4818983634427 + 0.8532316629489346*logv - 0.02296368399447732*logv2 + 0.001089868077621656*logv3 - 0.00004664262883242108*logv4;
  clm[36][28] = -41.89464693283499 + 1.45688715963291*logv - 0.03219658999284587*logv2 + 0.0007534873454764211*logv3 - 0.00001305407869675519*logv4;
  clm[36][30] = -75.32631042439981 + 2.480128866646128*logv - 0.04097154084154163*logv2 - 0.0001628278782948303*logv3 + 0.00005085353722845016*logv4 - 1.516530453311261e-6*logv5;
  clm[36][32] = -136.7056462526658 + 4.406483087254939*logv - 0.08285166790294592*logv2 + 0.00167444070296138*logv3 - 0.00005897386826131119*logv4 + 1.666554389449316e-6*logv5;
  clm[36][34] = -250.1605083184883 + 7.83288897566124*logv - 0.1381622446285008*logv2 + 0.002102233860818804*logv3 - 0.00003653483264098694*logv4 + 4.664259433244495e-7*logv5;
  clm[36][36] = -461.0834666424791 + 14.03001000964885*logv - 0.2298084208090586*logv2 + 0.002309712920490579*logv3 + 0.00002119518188562349*logv4 - 2.151814300366368e-6*logv5 + 4.515505355947891e-8*logv6;
  
  clm[37][2] = 0.9317460317460317/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - (5.787301587301587*nu)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) + (10.51111111111111*nu2)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - (5.746031746031746*nu3)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) + (0.8253968253968254*nu4)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
  clm[37][4] = -0.3033763003078555;
  clm[37][6] = 0.9457044588574083 - 0.4019656985601258*logv;
  clm[37][8] = -0.04439442470593296 + 0.3745299445314188*logv;
  clm[37][10] = -0.3136652453710446 + 0.1219468664798337*logv;
  clm[37][12] = 0.8177189144602828 - 0.5798093022395803*logv + 0.08078821140946496*logv2;
  clm[37][14] = 1.284913980605547 + 0.2038854139512792*logv - 0.07527409539262846*logv2;
  clm[37][16] = 1.531339926541965 + 0.1866573750934606*logv - 0.02450922868589236*logv2;
  clm[37][18] = 4.099688921517976 - 0.7459593750377647*logv + 0.156661679453346*logv2 - 0.0108246966115429*logv3;
  clm[37][20] = 8.534185314589889 - 0.294781922601334*logv - 0.07836839667660186*logv2 + 0.01008586811265982*logv3;
  clm[37][22] = 15.82754001445737 - 0.4836146720553714*logv - 0.04968940801600465*logv2 + 0.003283956409964866*logv3;
  clm[37][24] = 32.8312668904386 - 2.37460122169945*logv + 0.2537217974402495*logv2 - 0.02636782878387597*logv3 + 0.001087789183790067*logv4;
  clm[37][26] = 66.71394954726967 - 3.35321378269502*logv - 0.003886887015073011*logv2 + 0.0155104255682308*logv3 - 0.00101354325537265*logv4;
  clm[37][28] = 132.4610477208316 - 6.490854980872872*logv + 0.06463486397837842*logv2 + 0.008289053242732853*logv3 - 0.0003300094580931326*logv4;
  clm[37][30] = 268.4994664778055 - 15.14400260923173*logv + 0.6877648347085564*logv2 - 0.05057429501713011*logv3 + 0.003190078545072457*logv4 - 0.00008745078782966471*logv5;
  clm[37][32] = 542.9903848998359 - 28.22759096547058*logv + 0.6150319024085699*logv2 + 0.01146849552505824*logv3 - 0.002062122424494344*logv4 + 0.00008148192453335426*logv5;
  clm[37][34] = 1093.602530922594 - 56.45992445811537*logv + 1.310332375639502*logv2 - 0.003486923997638345*logv3 - 0.0009969044719572358*logv4 + 0.00002653049647077092*logv5;
  clm[37][36] = 2214.208734923214 - 117.4711217411974*logv + 3.589874294196811*logv2 - 0.1303226907833458*logv3 + 0.007016679743743772*logv4 - 0.0002998998877126532*logv5 + 5.85870283659742e-6*logv6;
  
  clm[38][2] = -1.13644380311047/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (10.66754850088183*nu)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) - (34.53439153439153*nu2)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (45.22222222222222*nu3)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) - (21.93650793650794*nu4)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (3.216931216931217*nu5)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4);
  clm[38][4] = -0.4870081031822398;
  clm[38][6] = 1.520717488141046 - 0.7146056863291125*logv;
  clm[38][8] = -0.2390154717046289 + 0.8121092038962241*logv;
  clm[38][10] = -1.360795771264292 + 0.3480187598223838*logv;
  clm[38][12] = 0.99138555780063 - 1.717764580297387*logv + 0.255330643466951*logv2;
  clm[38][14] = 1.507908055789311 + 0.8879560590370946*logv - 0.2901689275122252*logv2;
  clm[38][16] = -1.421395256036264 + 1.279759451779282*logv - 0.1243480923591406*logv2;
  clm[38][18] = 1.414739738709451 - 2.951609134558672*logv + 0.8392385620767554*logv2 - 0.06082024323851815*logv3;
  clm[38][20] = 5.67378758670618 + 0.5319069067040282*logv - 0.5735104748093973*logv2 + 0.0691187885320854*logv3;
  clm[38][22] = 1.708383796535126 + 2.499548124071645*logv - 0.5670705214708855*logv2 + 0.02961995129467317*logv3;
  clm[38][24] = 10.11527957660128 - 7.058135318101678*logv + 2.055218377569162*logv2 - 0.2536171205604553*logv3 + 0.01086562291554121*logv4;
  clm[38][26] = 28.45962057104517 - 0.7562076386101605*logv - 0.9913993631081915*logv2 + 0.1976484336706299*logv3 - 0.01234816982930192*logv4;
  clm[38][28] = 30.26221594377128 + 3.112510419974246*logv - 1.520235093470874*logv2 + 0.1612339446932344*logv3 - 0.005291646405991202*logv4;
  clm[38][30] = 72.84766569160553 - 25.55703494218354*logv + 5.795311332889403*logv2 - 0.7753307904378088*logv3 + 0.05490423123856927*logv4 - 0.001552927184190731*logv5;
  clm[38][32] = 160.7643979126525 - 12.65154071841114*logv - 1.876366353561502*logv2 + 0.4809361548124924*logv3 - 0.04621454753618396*logv4 + 0.001764814475155348*logv5;
  clm[38][34] = 214.8168868233563 - 10.85568428884155*logv - 3.325829227901758*logv2 + 0.5346070166657341*logv3 - 0.03347759999714602*logv4 + 0.0007562881223528649*logv5;
  clm[38][36] = 410.3141906975307 - 117.750734726572*logv + 20.12282549141787*logv2 - 2.467132477205625*logv3 + 0.1980411246969817*logv4 - 0.009218328079950404*logv5 + 0.0001849550993796256*logv6;
  
  clm[39][2] = 1.007936507936508/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - (6.396825396825397*nu)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) + (12.11111111111111*nu2)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - (7.26984126984127*nu3)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) + (1.206349206349206*nu4)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
  clm[39][4] = -0.5649969034605498;
  clm[39][6] = 2.313746482428831 - 1.116571384889238*logv;
  clm[39][8] = 0.5160368238423473 + 1.125433062547089*logv;
  clm[39][10] = -1.456038053273834 + 0.6308593749550774*logv;
  clm[39][12] = 3.736613377188061 - 4.1241154969115*logv + 0.6233658287767358*logv2;
  clm[39][14] = 8.370456006831638 + 0.976687831554885*logv - 0.6283131766241702*logv2;
  clm[39][16] = 4.333344214577687 + 2.496234251154498*logv - 0.352199762981975*logv2;
  clm[39][18] = 19.10406514874335 - 12.63307195248469*logv + 3.162558858231137*logv2 - 0.2320108155766226*logv3;
  clm[39][20] = 53.82146465856748 - 5.206188960923941*logv - 1.412221407057504*logv2 + 0.2338521712558021*logv3;
  clm[39][22] = 70.0345864832841 + 0.171096685566136*logv - 1.87957936701002*logv2 + 0.1310853923701484*logv3;
  clm[39][24] = 169.1342341955406 - 59.09736824909026*logv + 12.96325443292851*logv2 - 1.497204257851507*logv3 + 0.06476415941441779*logv4;
  clm[39][26] = 417.4573003184734 - 54.67301997805689*logv - 0.6009985964371365*logv2 + 0.8482860683078297*logv3 - 0.06527816067961158*logv4;
  clm[39][28] = 706.4567082034331 - 66.06107894128951*logv - 3.562841454317388*logv2 + 0.8804339802904781*logv3 - 0.03659154952437146*logv4;
  clm[39][30] = 1575.917913532102 - 385.8215461862765*logv + 64.36691263675563*logv2 - 7.466308304325709*logv3 + 0.5072958618729175*logv4 - 0.01446276143370877*logv5;
  clm[39][32] = 3635.358740515833 - 498.4991258447533*logv + 20.50509685884359*logv2 + 1.9743818066417*logv3 - 0.3268642137319187*logv4 + 0.01457754525461122*logv5;
  clm[39][34] = 6885.421757213782 - 830.5263940004251*logv + 24.11474661008281*logv2 + 2.866128264190659*logv3 - 0.2962561026935995*logv4 + 0.008171415425534119*logv5;
  clm[39][36] = 15090.31907481183 - 2971.287900452708*logv + 391.1013036841294*logv2 - 40.07280717262211*logv3 + 2.944827062648679*logv4 - 0.1332422199116431*logv5 + 0.002691450927226479*logv6;
  
  clm[40][2] = -1.198973865640532/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (11.35537918871252*nu)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) - (37.28571428571429*nu2)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (50.03703703703704*nu3)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) - (25.37566137566138*nu4)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (3.904761904761905*nu5)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4);
  clm[40][4] = -0.7860146347141633;
  clm[40][6] = 3.100564753186921 - 1.607862794240503*logv;
  clm[40][8] = 0.4205826127965186 + 1.927785469830124*logv;
  clm[40][10] = -3.877004904107114 + 1.263803686885443*logv;
  clm[40][12] = 4.385474596517428 - 8.179979488639322*logv + 1.292611382551439*logv2;
  clm[40][14] = 12.53313661145835 + 3.154118813872999*logv - 1.549807266108652*logv2;
  clm[40][16] = -6.19911229231497 + 8.744770361629299*logv - 1.016011463683539*logv2;
  clm[40][18] = 11.85566722370085 - 31.57657261635518*logv + 9.144459385126414*logv2 - 0.6927805831387458*logv3;
  clm[40][20] = 76.28128825993325 - 3.966239692274667*logv - 5.615040162044386*logv2 + 0.8306258138065642*logv3;
  clm[40][22] = 11.60303604307496 + 33.84468507384949*logv - 9.048930239384213*logv2 + 0.5445356769928663*logv3;
  clm[40][24] = 97.77936684877117 - 167.3269330823687*logv + 50.20512565893009*logv2 - 6.277512480507936*logv3 + 0.2784740310502573*logv4;
  clm[40][26] = 523.727680092908 - 62.22049127960053*logv - 15.94172191322331*logv2 + 4.659792816234262*logv3 - 0.3338830854888286*logv4;
  clm[40][28] = 229.1979910517696 + 127.7256246743062*logv - 50.41563387985028*logv2 + 5.931762270642291*logv3 - 0.2188846637933485*logv4;
  clm[40][30] = 841.6886697179521 - 1162.7190293375*logv + 316.2426587199725*logv2 - 42.94488808961563*logv3 + 3.076650664145605*logv4 - 0.08954960673757666*logv5;
  clm[40][32] = 3887.887542103975 - 553.6995123392211*logv - 52.41350362213934*logv2 + 22.07619442231817*logv3 - 2.536476277889114*logv4 + 0.1073676381567417*logv5;
  clm[40][34] = 1925.1291450582 + 523.2611253171925*logv - 276.8599054684644*logv2 + 41.60802142694813*logv3 - 2.819271558398994*logv4 + 0.0703873014286333*logv5;
  clm[40][36] = 6060.262795390097 - 9040.818202032667*logv + 2301.254982784303*logv2 - 305.8337084756566*logv3 + 24.80816647467751*logv4 - 1.167294444860693*logv5 + 0.02399724681870303*logv6;
  
  clm[41][2] = 1.122222222222222/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - (7.311111111111111*nu)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) + (14.51111111111111*nu2)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - (9.555555555555556*nu3)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) + (1.777777777777778*nu4)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
  clm[41][4] = -0.9687136927674012;
  clm[41][6] = 3.949345025494101 - 2.188479914382907*logv;
  clm[41][8] = 1.450819522709787 + 2.455960792807485*logv;
  clm[41][10] = -5.367092940655841 + 2.120010459409152*logv;
  clm[41][12] = 7.434506107411189 - 14.56163245680869*logv + 2.394722167828708*logv2;
  clm[41][14] = 27.83675774273276 + 3.466861610135333*logv - 2.68741043278555*logv2;
  clm[41][16] = -5.525024154870556 + 17.47917508734539*logv - 2.319800154349304*logv2;
  clm[41][18] = 32.82509532964537 - 76.51083723549713*logv + 22.41025607139729*logv2 - 1.74693378827354*logv3;
  clm[41][20] = 227.6748843968245 - 28.13501223362784*logv - 11.06146667227816*logv2 + 1.960447917951418*logv3;
  clm[41][22] = 106.5431856853826 + 79.56961093870061*logv - 25.40012715685198*logv2 + 1.692278681058607*logv3;
  clm[41][24] = 409.3339481118024 - 566.6425955401473*logv + 167.1535905082001*logv2 - 21.07257551092244*logv3 + 0.9557823768483714*logv4;
  clm[41][26] = 2223.644518888002 - 373.9917304948807*logv - 24.74381827573214*logv2 + 13.37114164040265*logv3 - 1.072600222907617*logv4;
  clm[41][28] = 1881.628490166671 + 293.0545382761048*logv - 177.8723626517579*logv2 + 23.10585605135314*logv3 - 0.9258794757587897*logv4;
  clm[41][30] = 5572.002975881062 - 5502.377802288338*logv + 1448.84525951395*logv2 - 195.5776673553648*logv3 + 14.11406462483662*logv4 - 0.4183421068507631*logv5;
  clm[41][32] = 25094.16943004614 - 4391.654378469716*logv + 53.80322834450295*logv2 + 72.89791413281601*logv3 - 10.21638088139456*logv4 + 0.4694728087991897*logv5;
  clm[41][34] = 27796.91265086389 - 65.65453103967191*logv - 1201.856582433087*logv2 + 208.3748331639634*logv3 - 15.14564300847432*logv4 + 0.4052537271674974*logv5;
  clm[41][36] = 77423.94103751671 - 60615.2994867837*logv + 14539.48454014588*logv2 - 1905.714498786844*logv3 + 154.2849849046803*logv4 - 7.30904240577787*logv5 + 0.1525888830305872*logv6;
  
  clm[42][2] = -1.28651595318262/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (12.31834215167549*nu)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) - (41.13756613756614*nu2)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (56.77777777777778*nu3)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) - (30.19047619047619*nu4)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (4.867724867724868*nu5)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4);
  clm[42][4] = -1.204473373789698;
  clm[42][6] = 4.928545277855381 - 2.85842274531645*logv;
  clm[42][8] = 1.60242913440268 + 3.677406462789674*logv;
  clm[42][10] = -9.792728053970096 + 3.442894087768516*logv;
  clm[42][12] = 8.15706823027834 - 24.18468537909053*logv + 4.085290295471216*logv2;
  clm[42][14] = 42.40602198077982 + 8.40929942043115*logv - 5.255791138505858*logv2;
  clm[42][16] = -37.36595074158928 + 40.15310680233446*logv - 4.920623385096529*logv2;
  clm[42][18] = -2.368917546960077 - 155.2235416979344*logv + 48.9955165815015*logv2 - 3.892495567265162*logv3;
  clm[42][20] = 373.9021206042986 - 31.71348789824722*logv - 30.58372092749245*logv2 + 5.007757644979262*logv3;
  clm[42][22] = -190.7232240497488 + 304.6240265434721*logv - 74.76841689384072*logv2 + 4.688407268365315*logv3;
  clm[42][24] = -115.3183195716514 - 1368.681618115125*logv + 461.3433920762529*logv2 - 60.43277918398773*logv3 + 2.781599466378549*logv4;
  clm[42][26] = 3801.916543435435 - 485.8330102894277*logv - 148.1677802312328*logv2 + 46.82932591903706*logv3 - 3.578572088860266*logv4;
  clm[42][28] = -2069.546593506033 + 2796.724854727269*logv - 779.4889689191681*logv2 + 87.80079646957263*logv3 - 3.350362493800596*logv4;
  clm[42][30] = -2216.443926556076 - 15842.48946028316*logv + 5002.569335354607*logv2 - 716.3325687132011*logv3 + 53.01106388909733*logv4 - 1.590197436611309*logv5;
  clm[42][32] = 44270.21739349133 - 5207.070322004987*logv - 1124.136930453927*logv2 + 388.0197048967238*logv3 - 46.10510879332555*logv4 + 2.045814370910557*logv5;
  clm[42][34] = -33175.91359577166 + 33036.80315925257*logv - 8713.214534643884*logv2 + 1129.078216704689*logv3 - 74.57744885207107*logv4 + 1.915350471466953*logv5;
  clm[42][36] = -48961.24587129322 - 208089.4555246855*logv + 62601.41663172181*logv2 - 8870.329900817125*logv3 + 744.3768764497379*logv4 - 35.92266767673773*logv5 + 0.7575760870589466*logv6;
  
  clm[43][2] = 1.274603174603175/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - (8.53015873015873*nu)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) + (17.71111111111111*nu2)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - (12.6031746031746*nu3)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) + (2.53968253968254*nu4)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
  clm[43][4] = -1.528069729721782;
  clm[43][6] = 5.721457463135437 - 3.617691287041132*logv;
  clm[43][8] = 2.542380301532624 + 4.611120799196872*logv;
  clm[43][10] = -14.79452734489813 + 5.528084547205788*logv;
  clm[43][12] = 6.959896607928711 - 36.87161926664798*logv + 6.543845124166662*logv2;
  clm[43][14] = 65.52303989923325 + 11.41680439484552*logv - 8.340805769374333*logv2;
  clm[43][16] = -81.25980682950376 + 78.23573736926671*logv - 9.999451650226551*logv2;
  clm[43][18] = -79.48272268422933 - 284.243322067472*logv + 95.94980423678454*logv2 - 7.891203829814776*logv3;
  clm[43][20] = 657.7427384197156 - 65.90005310521496*logv - 57.93941672504719*logv2 + 10.05815345292264*logv3;
  clm[43][22] = -721.9521179784645 + 787.717350627965*logv - 186.2196497857725*logv2 + 12.05830970340456*logv3;
  clm[43][24] = -1405.107122952243 - 3056.847818145596*logv + 1113.54553429601*logv2 - 150.9837918677505*logv3 + 7.136984834846633*logv4;
  clm[43][26] = 7999.381016107217 - 957.7484356080324*logv - 357.0661016732984*logv2 + 114.8346821173742*logv3 - 9.096823527590232*logv4;
  clm[43][28] = -10662.67405590583 + 9872.286019490392*logv - 2517.81586218447*logv2 + 278.4692868307404*logv3 - 10.90581048761255*logv4;
  clm[43][30] = -23609.11833651756 - 44031.62926773829*logv + 15020.46177038967*logv2 - 2223.341662725677*logv3 + 168.4595977255917*logv4 - 5.163881570553871*logv5;
  clm[43][32] = 119152.9690223181 - 10139.70072486662*logv - 3939.741651944966*logv2 + 1205.939002722184*logv3 - 144.5271195342046*logv4 + 6.581899843102792*logv5;
  clm[43][34] = -183087.3284117445 + 152819.0271141667*logv - 37041.03599757172*logv2 + 4595.217861733413*logv3 - 300.6091984589552*logv4 + 7.890771115831545*logv5;
  clm[43][36] = -414873.2600612241 - 732871.8211149285*logv + 235532.9637272103*logv2 - 34413.71151840571*logv3 + 2949.843561372946*logv4 - 144.9724725399336*logv5 + 3.11355489418417*logv6;
  
  clm[44][2] = -1.082536231884058/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) + (9.00536231884058*nu)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) - (25.44130434782609*nu2)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) + (28.58202898550725*nu3)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) - (11.78876811594203*nu4)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) + (1.515217391304348*nu5)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4);
  clm[44][4] = -0.2588238147191746;
  clm[44][6] = -0.1062194752115815 - 0.03624753583563652*logv;
  clm[44][8] = -0.3662630016509777 + 0.03923927085859232*logv;
  clm[44][10] = -0.473906218960086 + 0.00938172549914943*logv;
  clm[44][12] = -0.6759559087350896 + 0.002222842295700946*logv + 0.0006569419270778769*logv2;
  clm[44][14] = -1.073936641574548 + 0.01503779871303578*logv - 0.0007111634383055365*logv2;
  clm[44][16] = -1.747206558084181 + 0.01759913009109304*logv - 0.000170032215615262*logv2;
  clm[44][18] = -2.913072326076835 + 0.02450641490502319*logv - 0.00001079252903169146*logv2 - 7.937508681229153e-6*logv3;
  clm[44][20] = -4.963912039797517 + 0.03970565403994341*logv - 0.0003044696256179324*logv2 + 8.592640738324806e-6*logv3;
  clm[44][22] = -8.597514807339407 + 0.06414667289564332*logv - 0.0003265962339150949*logv2 + 2.054416276242294e-6*logv3;
  clm[44][24] = -15.0906288942477 + 0.1066864918462409*logv - 0.0004429092380028714*logv2 - 2.259577118870745e-7*logv3 + 7.192878259213243e-8*logv4;
  clm[44][26] = -26.78521126157571 + 0.1817638384618356*logv - 0.0007351515372726995*logv2 + 4.064528955115578e-6*logv3 - 7.786551327129466e-8*logv4;
  clm[44][28] = -47.99363862223175 + 0.314567721040896*logv - 0.00117768829948323*logv2 + 4.038336982768583e-6*logv3 - 1.861688189860188e-8*logv4;
  clm[44][30] = -86.69483957106254 + 0.5518511632811061*logv - 0.001953270053794121*logv2 + 5.320482936680519e-6*logv3 + 5.276882599503137e-9*logv4 - 5.214482249244057e-10*logv5;
  clm[44][32] = -157.711237369104 + 0.9791689817926402*logv - 0.003328270562433531*logv2 + 9.087504226615597e-6*logv3 - 4.032810237941305e-8*logv4 + 5.644865965322969e-10*logv5;
  clm[44][34] = -288.6802995366307 + 1.753948049794207*logv - 0.005754935148357305*logv2 + 0.00001441612715372376*logv3 - 3.743075570221365e-8*logv4 + 1.349632187534769e-10*logv5;
  
  clm[45][2] = -0.9150197628458498/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (8.398550724637681*nu)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) - (26.21086956521739*nu2)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (32.10652173913043*nu3)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) - (13.62608695652174*nu4)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (1.715217391304348*nu5)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4);
  clm[45][4] = -0.2214372499152394;
  clm[45][6] = 0.311820260455624 - 0.1449901433425461*logv;
  clm[45][8] = -0.1985965335658797 + 0.1326688465762823*logv;
  clm[45][10] = -0.2395753904501332 + 0.03210621860658977*logv;
  clm[45][12] = -0.06291143334331852 - 0.07124849527614473*logv + 0.01051107083324603*logv2;
  clm[45][14] = -0.1003441553085346 + 0.0526194868259877*logv - 0.009617837541092712*logv2;
  clm[45][16] = -0.0884999241195518 + 0.04050177160910765*logv - 0.002327542618978285*logv2;
  clm[45][18] = 0.1912827813225116 - 0.009760888845881007*logv + 0.007052764698138652*logv2 - 0.0005080005555986658*logv3;
  clm[45][20] = 0.7758318815481701 + 0.02956858905331916*logv - 0.005541844705963088*logv2 + 0.0004648305479094511*logv3;
  clm[45][22] = 2.104767390672864 + 0.02145300240113368*logv - 0.003354163772307657*logv2 + 0.0001124902459871822*logv3;
  clm[45][24] = 5.150051640066042 - 0.03536168483374984*logv + 0.002415478319459188*logv2 - 0.0004320882494929972*logv3 + 0.0000184137683435859*logv4;
  clm[45][26] = 11.69817610355992 - 0.1023304702500356*logv - 0.00354260203004252*logv2 + 0.0003513128544541412*logv3 - 0.00001684896194284639*logv4;
  clm[45][28] = 25.55653757918318 - 0.2989789993366779*logv - 0.002255307524993711*logv2 + 0.0001823081273424137*logv3 - 4.077494222579956e-6*logv4;
  clm[45][30] = 54.58685085349815 - 0.7566158047767067*logv + 0.003502461989539207*logv2 - 0.0002156641106327894*logv3 + 0.000018968917063251*logv4 - 5.339629823225915e-7*logv5;
  clm[45][32] = 114.7058632981158 - 1.710624210304025*logv + 0.006358328994274426*logv2 + 0.0002538197288106933*logv3 - 0.00001575999685415456*logv4 + 4.885866814532803e-7*logv5;
  clm[45][34] = 238.34338900527 - 3.757047632749863*logv + 0.02108237943291789*logv2 + 0.0001464613441294199*logv3 - 7.340465258438271e-6*logv4 + 1.182392943620543e-7*logv5;
  clm[45][36] = 491.3203555195273 - 8.054139178880908*logv + 0.05584713086256594*logv2 - 0.0002362471220676284*logv3 + 0.00001199688011907642*logv4 - 6.459513819580346e-7*logv5 + 1.290322822442766e-8*logv6;
  
  clm[46][2] = -1.102826086956522/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) + (9.208260869565217*nu)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) - (26.17173913043478*nu2)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) + (29.71826086956522*nu3)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) - (12.49891304347826*nu4)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) + (1.63695652173913*nu5)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4);
  clm[46][4] = -0.3558536021969772;
  clm[46][6] = 0.6848219014317448 - 0.3262278225207287*logv;
  clm[46][8] = -0.3974674308039463 + 0.3597725529668819*logv;
  clm[46][10] = -0.735442421011556 + 0.1160893457808775*logv;
  clm[46][12] = -0.1930546490489573 - 0.3552234647349253*logv + 0.05321229609330803*logv2;
  clm[46][14] = -0.5598370252747916 + 0.2750345142771023*logv - 0.05868390827855471*logv2;
  clm[46][16] = -1.485905998801564 + 0.2868288025931592*logv - 0.0189357872459758*logv2;
  clm[46][18] = -1.836176494480783 - 0.149891810784853*logv + 0.07944283161857534*logv2 - 0.005786443828616053*logv3;
  clm[46][20] = -2.915008139386642 + 0.3702348756898984*logv - 0.06857375609527935*logv2 + 0.006381441204906355*logv3;
  clm[46][22] = -5.503958917347518 + 0.6253147661027195*logv - 0.05443695583974804*logv2 + 0.002059126880323491*logv3;
  clm[46][24] = -8.826071157005258 + 0.3890479296455576*logv + 0.06785942413840212*logv2 - 0.01097688925428395*logv3 + 0.0004719247425869809*logv4;
  clm[46][26] = -14.87565959465179 + 1.240506557989916*logv - 0.1005716861978733*logv2 + 0.01003537208548958*logv3 - 0.0005204509172051639*logv4;
  clm[46][28] = -26.50214842280045 + 2.135470813724743*logv - 0.1280181410701982*logv2 + 0.006751626654946652*logv3 - 0.0001679361196154584*logv4;
  clm[46][30] = -45.79983300869119 + 2.824304652179728*logv - 0.007100815169141102*logv2 - 0.01304444487416082*logv3 + 0.001085927426641712*logv4 - 0.00003079099623356123*logv5;
  clm[46][32] = -80.77635143065299 + 5.605237152917323*logv - 0.2708490635512019*logv2 + 0.01634769986622234*logv3 - 0.001028747622888569*logv4 + 0.00003395711388975133*logv5;
  clm[46][34] = -145.9246037756114 + 9.74202120514073*logv - 0.4159118482444811*logv2 + 0.01708675651577554*logv3 - 0.0006184983289642443*logv4 + 0.00001095708692494632*logv5;
  
  clm[47][2] = -0.9600790513833992/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (8.894202898550725*nu)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) - (28.19347826086957*nu2)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (35.57608695652174*nu3)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) - (16.10434782608696*nu4)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (2.210869565217391*nu5)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4);
  clm[47][4] = -0.376741451928498;
  clm[47][6] = 1.310864679084881 - 0.5799605733701844*logv;
  clm[47][8] = 0.04851076625407583 + 0.5568079971210189*logv;
  clm[47][10] = -0.5424798099895575 + 0.2184951884727674*logv;
  clm[47][12] = 1.404035580337795 - 1.17685192714197*logv + 0.1681771333319365*logv2;
  clm[47][14] = 2.387460975849373 + 0.3718366135599075*logv - 0.161463342633705*logv2;
  clm[47][16] = 2.145722689898741 + 0.4715681802606644*logv - 0.06335929739264635*logv2;
  clm[47][18] = 6.608614441722657 - 2.049249658165187*logv + 0.4620702545224794*logv2 - 0.03251203555831461*logv3;
  clm[47][20] = 14.48327883582346 - 0.7434872798275493*logv - 0.2238089772045778*logv2 + 0.03121412425737004*logv3;
  clm[47][22] = 24.4136798359547 - 0.7589159501565517*logv - 0.182258252872612*logv2 + 0.0122486314813904*logv3;
  clm[47][24] = 52.18822007821059 - 6.833461618519658*logv + 1.039135981030489*logv2 - 0.1126818253458089*logv3 + 0.004713924695957991*logv4;
  clm[47][26] = 108.4935890233932 - 7.975356063175518*logv - 0.05363740089820516*logv2 + 0.06568878325118752*logv3 - 0.004525740350388126*logv4;
  clm[47][28] = 209.3490334414997 - 14.10925766958941*logv + 0.04658672283961687*logv2 + 0.04403273898269992*logv3 - 0.001775930834236817*logv4;
  clm[47][30] = 430.1087623596167 - 39.82384340519678*logv + 3.252454163417193*logv2 - 0.303669021089441*logv3 + 0.01972389955071901*logv4 - 0.0005467780938983337*logv5;
  clm[47][32] = 880.8409368998292 - 66.97428890393544*logv + 1.918323310877034*logv2 + 0.07852392102826274*logv3 - 0.0127751934991378*logv4 + 0.0005249501937071354*logv5;
  clm[47][34] = 1765.851422895049 - 130.7907003200528*logv + 3.921798075076343*logv2 + 0.03085203405543926*logv3 - 0.007660014524699314*logv4 + 0.0002059939729779548*logv5;
  clm[47][36] = 3612.67228438352 - 293.6807293646586*logv + 15.6849230337134*logv2 - 0.963973426417805*logv3 + 0.06136390891199014*logv4 - 0.00268058300048556*logv5 + 0.00005285162280725568*logv6;
  
  clm[48][2] = -1.143405797101449/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) + (9.614057971014493*nu)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) - (27.63260869565217*nu2)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) + (31.99072463768116*nu3)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) - (13.91920289855072*nu4)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) + (1.880434782608696*nu5)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4);
  clm[48][4] = -0.5498980542792365;
  clm[48][6] = 1.889978256288099 - 0.9061883958909131*logv;
  clm[48][8] = -0.1075228707431728 + 1.036141065127733*logv;
  clm[48][10] = -1.752643058019641 + 0.4983112357108356*logv;
  clm[48][12] = 1.719501225226593 - 2.729771325880267*logv + 0.4105887044236731*logv2;
  clm[48][14] = 3.222925006514556 + 1.260388252994542*logv - 0.4694695048624013*logv2;
  clm[48][16] = -1.531395138822354 + 2.14752334168744*logv - 0.2257819296716104*logv2;
  clm[48][18] = 3.978871611881167 - 6.108254247850642*logv + 1.697683375310207*logv2 - 0.1240235731442055*logv3;
  clm[48][20] = 14.62341011639791 + 0.1933967012275318*logv - 1.098001532986521*logv2 + 0.1418092725103202*logv3;
  clm[48][22] = 7.669106287038343 + 4.615345193536262*logv - 1.226445289632606*logv2 + 0.06820032155675718*logv3;
  clm[48][24] = 28.77261597365961 - 19.33796352859551*logv + 5.3464624186339*logv2 - 0.652009559046661*logv3 + 0.02809718070005173*logv4;
  clm[48][26] = 81.81817471529102 - 5.936923164773147*logv - 2.089966899164935*logv2 + 0.490830438619264*logv3 - 0.0321264792946461*logv4;
  clm[48][28] = 90.42792581554792 + 4.171913828660569*logv - 3.838028575109466*logv2 + 0.4470107175790772*logv3 - 0.01545058499769056*logv4;
  clm[48][30] = 216.6751246546928 - 87.56748670507984*logv + 19.54042934338128*logv2 - 2.550180690895954*logv3 + 0.17924681194398*logv4 - 0.0050922678215274*logv5;
  clm[48][32] = 513.9110190336289 - 55.48254837002383*logv - 3.735828096670766*logv2 + 1.414778575538806*logv3 - 0.1472545860993883*logv4 + 0.005822528547527597*logv5;
  clm[48][34] = 707.5390827893906 - 50.18311014125886*logv - 9.305637929802868*logv2 + 1.772897451392269*logv3 - 0.1186105321350536*logv4 + 0.002800228166926683*logv5;
  
  clm[49][2] = -1.035177865612648/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (9.720289855072464*nu)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) - (31.49782608695652*nu2)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (41.35869565217391*nu3)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) - (20.23478260869565*nu4)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (3.03695652173913*nu5)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4);
  clm[49][4] = -0.6403519511464899;
  clm[49][6] = 2.622545148321783 - 1.304911290082915*logv;
  clm[49][8] = 0.6466669574448678 + 1.350815284081879*logv;
  clm[49][10] = -2.019766194588097 + 0.8356024906776778*logv;
  clm[49][12] = 4.403298643937682 - 5.531236885058741*logv + 0.8513967374929285*logv2;
  clm[49][14] = 10.82719175407577 + 1.339396909631842*logv - 0.8813470575075019*logv2;
  clm[49][16] = 3.649239296959992 + 3.986148784394269*logv - 0.5451935620533526*logv2;
  clm[49][18] = 22.18681751909363 - 19.12348160325213*logv + 4.98494707622713*logv2 - 0.3703324050314274*logv3;
  clm[49][20] = 71.66761428626726 - 7.369859993656796*logv - 2.298364287093795*logv2 + 0.3833599086076317*logv3;
  clm[49][22] = 80.04969686440051 + 4.522368004165657*logv - 3.481948204135555*logv2 + 0.2371430781346467*logv3;
  clm[49][24] = 206.8844497185923 - 97.58792916800516*logv + 23.4295349200379*logv2 - 2.766850134076546*logv3 + 0.1208127341022671*logv4;
  clm[49][26] = 568.0995036414298 - 81.29356708126762*logv - 1.903501243028765*logv2 + 1.619321618278488*logv3 - 0.1250626682268133*logv4;
  clm[49][28] = 873.0177009995099 - 71.47865547594167*logv - 10.43240911317616*logv2 + 1.897824325308716*logv3 - 0.07736267000572882*logv4;
  clm[49][30] = 2031.740549223675 - 675.1907846405018*logv + 130.8350894778517*logv2 - 15.92247896891551*logv3 + 1.097885709416174*logv4 - 0.0315299801431667*logv5;
  clm[49][32] = 5078.337095011453 - 772.4798311170062*logv + 29.37075688790038*logv2 + 4.74891505507083*logv3 - 0.730398888383148*logv4 + 0.03263913754741249*logv5;
  clm[49][34] = 9012.949020662097 - 1116.076567221628*logv + 10.03025988396237*logv2 + 8.411607504680028*logv3 - 0.7441596374912714*logv4 + 0.02019028430428688*logv5;
  clm[49][36] = 20415.69223036054 - 5411.873330029131*logv + 870.6707803699643*logv2 - 97.42217555091364*logv3 + 7.379663367212416*logv4 - 0.3374886652466778*logv5 + 0.006857304510818058*logv6;
  
  clm[50][2] = -1.204275362318841/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) + (10.22275362318841*nu)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) - (29.82391304347826*nu2)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) + (35.39942028985507*nu3)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) - (16.04963768115942*nu4)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) + (2.245652173913043*nu5)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4);
  clm[50][4] = -0.8409269252192606;
  clm[50][6] = 3.36568477946288 - 1.77612925594619*logv;
  clm[50][8] = 0.5883308594508383 + 2.13894870322969*logv;
  clm[50][10] = -4.476425680169599 + 1.493594913994802*logv;
  clm[50][12] = 5.056138071228948 - 9.885163207371356*logv + 1.577317566913983*logv2;
  clm[50][14] = 15.74344258064181 + 3.66047975697938*logv - 1.899524684387208*logv2;
  clm[50][16] = -8.027404008201136 + 11.23644084516793*logv - 1.326408811639301*logv2;
  clm[50][18] = 13.78910916740288 - 41.91692583873047*logv + 12.24857384507332*logv2 - 0.9338399588379286*logv3;
  clm[50][20] = 102.66299407927 - 6.433455445065282*logv - 7.429469786939283*logv2 + 1.124600454777358*logv3;
  clm[50][22] = 10.58743930669787 + 48.38683706705819*logv - 12.89662645549511*logv2 + 0.7852911652324602*logv3;
  clm[50][24] = 121.0285091898161 - 243.9123716302552*logv + 74.10813015869774*logv2 - 9.30601970696849*logv3 + 0.4146551178159076*logv4;
  clm[50][26] = 759.1623857216791 - 95.24171308763651*logv - 22.59846031882036*logv2 + 6.872552755161909*logv3 - 0.4993589422451137*logv4;
  clm[50][28] = 277.7976842041892 + 213.0348413043547*logv - 79.69863243427252*logv2 + 9.362905188851921*logv3 - 0.3486946532513614*logv4;
  clm[50][30] = 1138.52380329595 - 1854.8004954335*logv + 514.4706618390914*logv2 - 70.23093837018632*logv3 + 5.044365146399711*logv4 - 0.1472962171761295*logv5;
  clm[50][32] = 6174.131791542555 - 889.2048314338737*logv - 81.38090750991232*logv2 + 35.58353302958194*logv3 - 4.15016547110566*logv4 + 0.177385205307978*logv5;
  clm[50][34] = 2687.305429213063 + 1092.661493286257*logv - 494.2586714548255*logv2 + 72.72997524073174*logv3 - 4.9245190029218*logv4 + 0.123865355006351*logv5;
  
  clm[51][2] = -1.140316205533597/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (10.8768115942029*nu)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) - (36.12391304347826*nu2)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (49.45434782608696*nu3)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) - (26.01739130434783*nu4)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (4.193478260869565*nu5)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4);
  clm[51][4] = -1.019423991363057;
  clm[51][6] = 4.13051673830906 - 2.319842293480737*logv;
  clm[51][8] = 1.538995432800233 + 2.645353761538311*logv;
  clm[51][10] = -6.072921432276953 + 2.364902890152961*logv;
  clm[51][12] = 7.661841006938047 - 16.24778096344635*logv + 2.690834133310984*logv2;
  clm[51][14] = 30.72242499215001 + 4.030703251311645*logv - 3.068401768617466*logv2;
  clm[51][16] = -9.358901755970383 + 20.88332673187858*logv - 2.743100872275835*logv2;
  clm[51][18] = 30.19610168111017 - 89.39351962646247*logv + 26.57775402655966*logv2 - 2.080770275732135*logv3;
  clm[51][20] = 257.0967138771569 - 31.25672962593244*logv - 13.49177731645165*logv2 + 2.372736065409964*logv3;
  clm[51][22] = 84.74350888847101 + 107.1341313746844*logv - 32.10480030230214*logv2 + 2.121187139596461*logv3;
  clm[51][24] = 389.5016116768593 - 690.5004457053954*logv + 208.9774595317963*logv2 - 26.53077070202677*logv3 + 1.206764722165246*logv4;
  clm[51][26] = 2568.47911984676 - 427.4379662854392*logv - 35.49077327572523*logv2 + 17.25054579199341*logv3 - 1.376093368951278*logv4;
  clm[51][28] = 1743.755833504844 + 528.3274512772239*logv - 245.9978824676647*logv2 + 30.92085958094758*logv3 - 1.230204909705825*logv4;
  clm[51][30] = 5562.333787999707 - 6988.092623416122*logv + 1906.226342762123*logv2 - 260.1941461229521*logv3 + 18.8542140230323*logv4 - 0.5599007681518937*logv5;
  clm[51][32] = 29731.64581483237 - 5123.202366860865*logv - 1.126908911273982*logv2 + 102.5133631244191*logv3 - 13.95858370325641*logv4 + 0.6384639194143135*logv5;
  clm[51][34] = 27399.64027212418 + 2176.838288713641*logv - 1888.153222710886*logv2 + 301.8698388968886*logv3 - 21.46764348723607*logv4 + 0.5707762758366449*logv5;
  clm[51][36] = 80553.5263014818 - 80191.35293607062*logv + 20100.49281167585*logv2 - 2676.313234573568*logv3 + 218.0470425391673*logv4 - 10.35653090158111*logv5 + 0.2164802470185193*logv6;
  
  clm[52][2] = -1.285434782608696/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) + (11.03434782608696*nu)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) - (32.74565217391304*nu2)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) + (39.94434782608696*nu3)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) - (18.89021739130435*nu4)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4) + (2.732608695652174*nu5)/(1. - 8.*nu + 21.*nu2 - 20.*nu3 + 5.*nu4);
  clm[52][4] = -1.228894846397012;
  clm[52][6] = 5.030529583961025 - 2.936050402686558*logv;
  clm[52][8] = 1.70937994890574 + 3.774101311105569*logv;
  clm[52][10] = -10.18950383822775 + 3.608097208623382*logv;
  clm[52][12] = 8.366940465798152 - 25.44694447903912*logv + 4.310195983557951*logv2;
  clm[52][14] = 45.11947729663597 + 8.705833558757519*logv - 5.540475837125687*logv2;
  clm[52][16] = -39.94163736011378 + 43.03787602445976*logv - 5.296777631155464*logv2;
  clm[52][18] = -4.385335946370249 - 167.6534714482438*logv + 53.03094317686705*logv2 - 4.218317551061102*logv3;
  clm[52][20] = 407.460809520858 - 35.83658676679717*logv - 32.92852871606018*logv2 + 5.422372104222673*logv3;
  clm[52][22] = -211.0523913581381 + 335.8988428275194*logv - 82.44261471441168*logv2 + 5.183868698965051*logv3;
  clm[52][24] = -145.8569369648549 - 1519.742053244957*logv + 513.175652310128*logv2 - 67.2405754170668*logv3 + 3.096298236113181*logv4;
  clm[52][26] = 4264.387684807067 - 553.5459732151909*logv - 162.5252429611796*logv2 + 51.9452637182827*logv3 - 3.980089450029835*logv4;
  clm[52][28] = -2330.592381148574 + 3177.940260693186*logv - 884.1063629917584*logv2 + 99.53655462669031*logv3 - 3.805024945267646*logv4;
  clm[52][30] = -2637.525635010891 - 18084.68523200147*logv + 5721.957459333214*logv2 - 819.3846731638803*logv3 + 60.61523288420706*logv4 - 1.818177536595557*logv5;
  clm[52][32] = 51392.42610271649 - 6131.531094545562*logv - 1260.224356764995*logv2 + 441.3167344144748*logv3 - 52.60222085945583*logv4 + 2.337148646497724*logv5;
  clm[52][34] = -37257.95714008222 + 38537.15674191327*logv - 10174.42166855137*logv2 + 1316.477509737049*logv3 - 86.89819951085121*logv4 + 2.234349004557094*logv5;
  clm[53][2] = -1.275494071146245/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (12.36376811594203*nu)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) - (42.07173913043478*nu2)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (59.86304347826087*nu3)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) - (33.45217391304348*nu4)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4) + (5.680434782608696*nu5)/(1. - 9.*nu + 27.*nu2 - 30.*nu3 + 9.*nu4);
  
  clm[53][4] = -1.523974913889576;
  clm[53][6] = 5.742049707288021 - 3.624753583563652*logv;
  clm[53][8] = 2.554702907152556 + 4.623351705201544*logv;
  clm[53][10] = -14.80533661374933 + 5.52403353038235*logv;
  clm[53][12] = 7.105103270853597 - 37.08703463822647*logv + 6.569419270778769*logv2;
  clm[53][14] = 66.17195105673082 + 11.49660897426948*logv - 8.379255330752209*logv2;
  clm[53][16] = -81.44452822895783 + 78.46613224958377*logv - 10.0116301674896*logv2;
  clm[53][18] = -79.02768150067391 - 287.3751948136805*logv + 96.70942970783317*logv2 - 7.937508681229153*logv3;
  clm[53][20] = 666.533234814899 - 66.92150627132467*logv - 58.45528908767589*logv2 + 10.12424526257964*logv3;
  clm[53][22] = -727.9839798353665 + 792.4497047336987*logv - 187.1579303993269*logv2 + 12.09656410897396*logv3;
  clm[53][24] = -1419.344064986892 - 3103.211108495551*logv + 1127.401581462422*logv2 - 152.4851412491087*logv3 + 7.192878259213243*logv4;
  clm[53][26] = 8096.001136693161 - 978.2162601183887*logv - 361.0302312222648*logv2 + 116.0819975511436*logv3 - 9.174473574103215*logv4;
  clm[53][28] = -10904.37341710349 + 9971.9569491493*logv - 2539.186836107474*logv2 + 280.4419455035808*logv3 - 10.96176602570271*logv4;
  clm[53][30] = -24306.54461732997 - 44766.82210352143*logv + 15264.70895789946*logv2 - 2255.061598661878*logv3 + 170.4730658933476*logv4 - 5.214482249244057*logv5;
  clm[53][32] = 119694.6180865245 - 10346.7661306268*logv - 3985.490471408743*logv2 + 1223.040621359429*logv3 - 146.3814345367805*logv4 + 6.651041193008132*logv5;
  clm[53][34] = -190835.8420153851 + 155248.761011927*logv - 37494.76176359543*logv2 + 4644.458871973679*logv3 - 303.3466543654318*logv4 + 7.946740136770437*logv5;
  clm[53][36] = -437586.5697006182 - 744578.5209195252*logv + 239955.8140059951*logv2 - 35035.50979299194*logv3 + 2997.605896497844*logv4 - 146.995260056266*logv5 + 3.150202203229408*logv6;
  
  //printf("Computing flm 22PN\n");

  /** rho_lm */
  // Taylor series to 22 PN
  for (int k=0; k<KMAX; k++) {
    rholm[k] = Taylorseries(v,clm[k],44);
  }
    
  /** Amplitudes */
#pragma omp simd 
  for (int k = 0; k < KMAX; k++) {
      flm[k] = gsl_pow_int(rholm[k], LINDEX[k]);
      // printf("flm k =%.d %.5f\n", k, flm[k]);
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
  
  if (kmaxT4PN+kmaxT6PN+kmaxT8PN+kmaxP23+kmaxP15+kmaxP42+kmaxP33+kmaxP51+kmaxP62+19 != KMAX) {
    errorexit("Wrong function: not all multipoles are written.\n");
  }

  // Set all the unused multipoles to zero (ell = 9, 10, only used in the 22PN flux)
  for (int k=35; k<KMAX; k++) {
    flm[k] = 0.;
  }
  
  /** Amplitudes */
#pragma omp simd
  for (int k = 0; k < 35; k++) {
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
    for (int k=35; k<KMAX; k++) clm[k][0] = 0.;
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
  
  if (kmaxTaylor+kmaxPade32+kmaxPade22+kmaxPade42+kmaxPade51+19 != KMAX) {
    errorexit("Wrong function: not all multipoles are written.\n");
  }

  // Set all the unused multipoles to zero (ell = 9, 10, only used in the 22PN flux)
  for (int k=35; k<KMAX; k++) {
    flm[k] = 0.;
  }
  
  /** Amplitudes */
#pragma omp simd
  for (int k = 0; k < 35; k++) {
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
  eob_wav_flm(x, nu, rholm, flm);

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
  eob_wav_flm(x, nu, rholm, flm);

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

  // f-odd corrections flag
  int fodd_flag;
  if (fabs(chi1) < 1e-14 && fabs(chi2) < 1e-14) { 
    fodd_flag = 0; // if zero spins, spin corrections for f_modd modes set to 1
  } else {
    fodd_flag = 1;
  } 
  
  /** l>=2, m=odd*/

  if (fodd_flag) {

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

  } else {

    f21S = 1.;
    f33S = 1.; 
    f31S = 1.;
    f43S = 1.;
    f41S = 1.;
    f55S = 1.;

  }
	    
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
 */
void eob_wav_ringdown_v1(Dynamics *dyn, Waveform_lm *hlm)
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
  
}
	
/** 
 * Function: eob_wav_ringdown_HM
 * -----------------------------
 *   Ringdown calculation and match to the dynamics
 *   This refers to the higher modes paper: arXiv:2001.09082
 * 
 *   @param[in]  dyn : dynamics
 *   @param[out] hlm : waveform
*/
void eob_wav_ringdown_HM(Dynamics *dyn, Waveform_lm *hlm)
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
        tm = t_lm[k][j] - tmatch[k];
        tm /= fact;

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
    
    double vphi3    = gsl_pow_int(rw*Omega,3);
    hNewt.ampli[0]  = ChlmNewt_ampli[0] * vphi3; /* (2,1) */
    
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
 */
void eob_wav_ringdown_bhns(Dynamics *dyn, Waveform_lm *hlm)
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
  
  
}
