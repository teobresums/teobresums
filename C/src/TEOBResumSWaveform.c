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

/** Leading-order (Newtonian) prefactor  of the multipolar resummed waveform. 
    Reference: Damour, Iyer & Nagar, PRD 79, 064004 (2009) */
void eob_wav_hlmNewt_v1(double r,
			double Omega,
			double phi,
			double nu,
			Waveform_lm_t *hlmNewt)
{
  /** Shorthands */
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
  
  /** Polynomials in nu */
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
    
  /** Compute hlmNewt (without phase factor) in complex Polar coords */
#pragma omp simd
  for (int k = 0; k < KMAX; k++) {
    hlmNewt->phase[k] = - phim[k] + ChlmNewt_phase[k];
    hlmNewt->ampli[k] = ChlmNewt_ampli[k] * Alm[k];
  }
  
}

/** Leading-order (Newtonian) prefactor  of the multipolar resummed waveform.
    New terms vphi and vOmg for higher modes. 
    Reference: arXiv:2001.09082 */
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
    
  /** Compute hlmNewt (without phase factor) in complex Polar coords */
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

/** Leading-order (Newtonian) prefactor  of the multipolar resummed waveform.
    Eccentric version */
void eob_wav_hlmNewt_ecc(double r,
			 double rdot,
			 double r2dot,
			 double r3dot,
			 double r4dot,
			 double r5dot,
			 double Omega,
			 double Omegadot,
			 double Omega2dot,
			 double Omega3dot,
			 double Omega4dot,
			 double phi,
			 double nu,
			 Waveform_lm_t *hlmNewt)
{  
  /** Shorthands */
  double nu2   = nu*nu;
  double nu3   = nu*nu2;

  double vOmg  = pow(Omega,1./3.);
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

  /* Eccentric variables */
  double SQrdot = SQ(rdot);
  
  double Re_vphi22_ecc = vphi2 - 0.5*(SQrdot + r*r2dot);
  double Im_vphi22_ecc = 2.*vphi*rdot + 0.5*Omegadot*SQ(r);

  double Re_vphi21_ecc = vphi3 - 6.*vphi*SQrdot - 3.*vphi*r*r2dot
    - 6.*SQ(r)*rdot*Omegadot - SQ(r)*r*Omega2dot;
  double Im_vphi21_ecc = 3.*SQ(r)*(2.*rdot*SQ(Omega) + vphi*Omegadot);

  double Re_vphi33_ecc = vphi3 - vphi*(r*r2dot + 2.*SQrdot)
    - 1./9*r*SQ(r)*Omega2dot - SQ(r)*rdot*Omegadot;
  double Im_vphi33_ecc = vphi*SQ(r)*Omegadot + 3.*vphi2*rdot
    - 2./9.*SQrdot*rdot -2./3.*r*rdot*r2dot - 1./9.*SQ(r)*r3dot; 

  double Re_vphi32_ecc = vphi4 - 3.*vphi2*(r*r2dot + 3.*SQrdot)
    - SQ(SQ(r))*Omega*Omega2dot - 9.*vphi*SQ(r)*rdot*Omegadot
    - 0.75*SQ(SQ(r)*Omegadot);
  double Im_vphi32_ecc = 6.*vphi3*rdot -0.5*vphi*SQ(r)*r3dot
    - 4.5*Omega*SQ(r)*rdot*r2dot - 3.*vphi*rdot*SQrdot
    - 1.5*SQ(r)*r*rdot*Omega2dot - 0.125*SQ(SQ(r))*Omega3dot
    - 1.5*SQ(r)*r*r2dot*Omegadot - 4.5*SQ(r)*SQrdot*Omegadot
    + 3.*SQ(r)*vphi2*Omegadot;

  double Re_vphi31_ecc = vphi3 - 9.*vphi*r*r2dot - 18.*vphi*SQrdot
    -9.*SQ(r)*rdot*Omegadot - SQ(r)*r*Omega2dot;
  double Im_vphi31_ecc = 9.*vphi2*rdot - 3.*SQ(r)*r3dot - 18.*r*rdot*r2dot
    - 6.*SQrdot*rdot + 3.*vphi*SQ(r)*Omegadot;

  double Re_vphi44_ecc = vphi4 + 0.015625*SQ(r)*r*r4dot
    - 0.1875*SQ(SQ(r)*Omegadot) + 0.1875*SQ(r)*rdot*r3dot
    + 0.09375*SQ(SQrdot) + 0.140625*SQ(r*r2dot) -4.5*vphi2*SQrdot
    - 1.5*vphi2*r*r2dot - 0.25*vphi*SQ(r)*r*Omega2dot
    - 3.*vphi*SQ(r)*rdot*Omegadot + 0.5625*r*SQrdot*r2dot;
  double Im_vphi44_ecc = 4.*vphi3*rdot - 1.5*vphi*SQrdot*rdot
    - 0.015625*SQ(SQ(r))*Omega3dot - 0.25*vphi*SQ(r)*r3dot
    + 1.5*vphi2*SQ(r)*Omegadot - 0.25*SQ(r)*r*rdot*Omega2dot
    - 2.25*vphi*r*rdot*r2dot - 0.375*SQ(r)*r*r2dot*Omegadot
    - 1.125*SQ(r*rdot)*Omegadot;

  double Re_vphi43_ecc = vphi5 + 5./81.*SQ(SQ(r))*Omega*r4dot
    + 20./27.*vphi*SQ(r*r2dot) + 80./81.*vphi*SQ(r)*rdot*r3dot
    + 40./9.*vphi*r*SQrdot*r2dot + 40./27.*vphi*SQ(SQrdot)
    - 10./3.*vphi3*r*r2dot - 40./3.*vphi3*SQrdot
    + 10./27.*SQ(SQ(r))*r2dot*Omega2dot + 40./27.*SQ(r)*r*SQrdot*Omega2dot
    - 10./9.*vphi2*SQ(r)*r*Omega2dot + 20./81.*SQ(SQ(r))*rdot*Omega3dot
    + 1./81.*SQ(SQ(r))*r*Omega4dot + 20./81.*SQ(SQ(r))*r3dot*Omegadot
    + 80./27.*SQ(r)*r*rdot*r2dot*Omegadot
    + 80./27.*SQ(r)*SQrdot*rdot*Omegadot - 40./3.*vphi2*SQ(r)*rdot*Omegadot
    - 5./3.*vphi*SQ(SQ(r)*Omegadot);
  double Im_vphi43_ecc = 20./3.*vphi4*rdot - 20./27.*vphi2*SQ(r)*r3dot
    - 80./9.*vphi2*r*rdot*r2dot - 80./9.*vphi2*SQrdot*rdot 
    - 80./27.*SQ(SQ(r))*Omega*rdot*Omega2dot - 5./27.*vphi*SQ(SQ(r))*Omega3dot
    - 10./3.*SQ(SQ(r))*Omega*Omegadot*r2dot
    - 40./3.*vphi*SQ(r*rdot)*Omegadot + 10./3.*vphi3*SQ(r)*Omegadot
    - 10./27.*SQ(SQ(r))*r*Omegadot*Omega2dot
    - 20./9.*SQ(SQ(r)*Omegadot)*rdot;

  double Re_vphi42_ecc = vphi4 + 2.25*SQ(r*r2dot) + 0.25*SQ(r)*r*r4dot
    + 3.*SQ(r)*rdot*r3dot + 9.*r*SQrdot*r2dot + 1.5*SQ(SQrdot)
    - 6.*vphi2*r*r2dot - 18.*vphi2*SQrdot - SQ(SQ(r))*Omega*Omega2dot
    - 12*vphi*SQ(r)*rdot*Omegadot - 0.75*SQ(SQ(r)*Omegadot);
  double Im_vphi42_ecc =  8.*vphi3*rdot - 2.*vphi*SQ(r)*r3dot
    - 18.*vphi*r*rdot*r2dot - 12.*vphi*SQrdot*rdot
    - 2.*SQ(r)*r*rdot*Omega2dot - 0.125*SQ(SQ(r))*Omega3dot
    - 3.*SQ(r)*r*r2dot*Omegadot - 9.*SQ(r*rdot)*Omegadot
    + 3.*vphi2*SQ(r)*Omegadot;

  double Re_vphi41_ecc = vphi5 + 60.*vphi*SQ(r*r2dot)
    + 5.*SQ(SQ(r))*Omega*r4dot + 80.*vphi*SQ(r)*rdot*r3dot
    + 360.*vphi*r*SQrdot*r2dot + 120.*vphi*SQ(SQrdot) - 30.*vphi3*r*r2dot
    - 120.*vphi3*SQrdot + 30.*SQ(SQ(r))*r2dot*Omega2dot
    + 120.*SQ(r*rdot)*r*Omega2dot - 10.*vphi2*SQ(r)*r*Omega2dot
    + 20.*SQ(SQ(r))*rdot*Omega3dot + SQ(SQ(r))*r*Omega4dot
    + 20.*SQ(SQ(r))*r3dot*Omegadot + 240.*SQ(r)*r*rdot*r2dot*Omegadot
    + 240.*SQ(r*rdot)*rdot*Omegadot - 120.*vphi2*SQ(r)*rdot*Omegadot
    - 15.*vphi*SQ(SQ(r)*Omegadot);
  double Im_vphi41_ecc = 20.*vphi4*rdot - 20.*vphi2*SQ(r)*r3dot
    - 240.*vphi2*r*rdot*r2dot - 240.*vphi2*SQrdot*rdot
    - 80.*SQ(SQ(r))*Omega*rdot*Omega2dot - 5.*vphi*SQ(SQ(r))*Omega3dot
    - 90.*SQ(SQ(r))*Omega*Omegadot*r2dot - 360.*vphi*SQ(r*rdot)*Omegadot
    + 10.*vphi3*SQ(r)*Omegadot - 10.*SQ(SQ(r))*r*Omegadot*Omega2dot
    - 60.*SQ(SQ(r)*Omegadot)*rdot;

  double Re_vphi55_ecc = vphi5 + 0.48*vphi*SQ(r*r2dot)
    + 0.04*SQ(SQ(r))*Omega*r4dot + 0.64*vphi*SQ(r)*rdot*r3dot
    + 2.88*vphi*r*SQrdot*r2dot + 0.96*vphi*SQ(SQrdot) - 2.*vphi3*r*r2dot
    - 8.*vphi3*SQrdot + 0.08*SQ(SQ(r))*r2dot*Omega2dot
    + 0.32*SQ(r*rdot)*r*Omega2dot - 0.4*vphi2*SQ(r)*r*Omega2dot
    + 0.04*SQ(SQ(r))*rdot*Omega3dot + 0.0016*SQ(SQ(r))*r*Omega4dot
    + 0.08*SQ(SQ(r))*r3dot*Omegadot + 0.96*SQ(r)*r*rdot*r2dot*Omegadot
    + 0.96*SQ(r*rdot)*rdot*Omegadot - 6.*vphi2*SQ(r)*rdot*Omegadot
    - 0.6*vphi*SQ(SQ(r)*Omegadot);
  double Im_vphi55_ecc = 5.*vphi4*rdot + 0.0016*SQ(SQ(r))*r5dot
    + 0.064*SQ(r)*r*r2dot*r3dot + 0.288*SQ(r*r2dot)*rdot
    + 0.032*SQ(r)*r*rdot*r4dot + 0.192*SQ(r*rdot)*r3dot
    + 0.384*r*SQrdot*rdot*r2dot + 0.0384*SQ(SQrdot)*rdot
    - 0.4*vphi2*SQ(r)*r3dot - 4.8*vphi2*r*rdot*r2dot
    - 4.8*vphi2*SQrdot*rdot - 0.8*SQ(SQ(r))*Omega*rdot*Omega2dot
    - 0.04*vphi*SQ(SQ(r))*Omega3dot - 1.2*SQ(SQ(r))*Omega*r2dot*Omegadot
    - 4.8*vphi*SQ(r*rdot)*Omegadot + 2.*vphi3*SQ(r)*Omegadot
    - 0.08*SQ(SQ(r))*r*Omegadot*Omega2dot - 0.6*SQ(SQ(r)*Omegadot)*rdot;

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

  double phi22_ecc = atan(Im_vphi22_ecc/Re_vphi22_ecc);
  double phi21_ecc = atan(Im_vphi21_ecc/Re_vphi21_ecc);
  double phi33_ecc = atan(Im_vphi33_ecc/Re_vphi33_ecc);
  double phi32_ecc = atan(Im_vphi32_ecc/Re_vphi32_ecc);
  double phi31_ecc = atan(Im_vphi31_ecc/Re_vphi31_ecc);
  double phi44_ecc = atan(Im_vphi44_ecc/Re_vphi44_ecc);
  double phi43_ecc = atan(Im_vphi43_ecc/Re_vphi43_ecc);
  double phi42_ecc = atan(Im_vphi42_ecc/Re_vphi42_ecc);
  double phi41_ecc = atan(Im_vphi41_ecc/Re_vphi41_ecc);
  double phi55_ecc = atan(Im_vphi55_ecc/Re_vphi55_ecc);
  
  if (Re_vphi22_ecc < 0.)    phi22_ecc = phi22_ecc + Pi;
  if (Re_vphi21_ecc < 0.)    phi21_ecc = phi21_ecc + Pi;
  if (Re_vphi33_ecc < 0.)    phi33_ecc = phi33_ecc + Pi;
  if (Re_vphi32_ecc < 0.)    phi32_ecc = phi32_ecc + Pi;
  if (Re_vphi31_ecc < 0.)    phi31_ecc = phi31_ecc + Pi;
  if (Re_vphi44_ecc < 0.)    phi44_ecc = phi44_ecc + Pi;
  if (Re_vphi43_ecc < 0.)    phi43_ecc = phi43_ecc + Pi;
  if (Re_vphi42_ecc < 0.)    phi42_ecc = phi42_ecc + Pi;
  if (Re_vphi41_ecc < 0.)    phi41_ecc = phi41_ecc + Pi;
  if (Re_vphi55_ecc < 0.)    phi55_ecc = phi55_ecc + Pi;
  
  
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
    
  /** Compute hlmNewt (without phase factor) in complex Polar coords */
#pragma omp simd
  for (int k = 0; k < KMAX; k++) {
    hlmNewt->phase[k] = - phim[k] + phiecc[k] + ChlmNewt_phase[k];
    hlmNewt->ampli[k] = ChlmNewt_ampli[k] * Alm[k];
  }
}

/** Tail contribution to the resummed wave.   
    Ref. Damour, Iyer & Nagar, PRD 79, 064004 (2009) */
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

/** Alternative implementation of the phase of the tail factor */
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

/** Residual phase corrections delta_{lm} up to l=m=5.
    Reference(s)
    Damour, Iyer & Nagar, PRD 79, 064004 (2008)
    Fujita & Iyer, PRD 82 044051 (2010)
    Faye et al., Class. Q. Grav. 29 175004 (2012)
    Damour, Nagar & Bernuzzi, PRD 87, 084035 (2013) */   
//TODO: this routine can be optimized: precompute coefficients c(nu)
void eob_wav_deltalm_v1(double Hreal,double Omega,double nu, double *dlm)
{
    
  /** Useful shorthands*/
  const double Pi2 = SQ(Pi);
  double nu2    = SQ(nu);
  double y      = cbrt(Hreal*Omega*Hreal*Omega);
  double sqrt_y = sqrt(y);
  double y3     = y*y*y;
  double y32    = Hreal*Omega;
  
  /** Leading order contributions*/
  double delta22LO = 7./3.   * y32;
  double delta21LO = 2./3.   * y32;
  double delta33LO = 13./10. * y32;
  double delta31LO = 13./30. * y32;
  
  /** Init phase */
#pragma omp simd
  for (int k = 0; k < KMAX; k++) {
    dlm[k] = 0.;
  }
  
  /** Residual phases in Pade-resummed form when possible */
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

/** Residual phase corrections delta_{lm} up to l=m=5 for higher modes
 *  Ref: arXiv:2001.09082 */
void eob_wav_deltalm_HM(double Hreal,double Omega,double nu, double *dlm)
{
    
  /** Useful shorthands*/
  const double Pi2 = SQ(Pi);
  double nu2    = SQ(nu);
  double y      = cbrt(Hreal*Omega*Hreal*Omega);
  double sqrt_y = sqrt(y);
  double y3     = y*y*y;
  double y32    = Hreal*Omega;
  
  /** Leading order contributions*/
  double delta22LO = 7./3.   * y32;
  double delta21LO = 2./3.   * y32;
  double delta33LO = 13./10. * y32;
  double delta32LO = (10. + 33.*nu)/(15.*(1. - 3.*nu))*y32;        
  double delta31LO = 13./30. * y32;
  double delta42LO = 7.*(1. + 6.*nu)/(15.*(1. - 3.*nu))*y32;
  double delta43LO = (486. + 4961.*nu)/(810.*(1. - 2.*nu))*y32;
  double delta44LO = (112. + 219.*nu)/(120.*(1. - 3.*nu))*y32;
	  
  /** Init phase */
#pragma omp simd
  for (int k = 0; k < KMAX; k++) {
    dlm[k] = 0.;
  }
  
  /** Residual phases in Pade-resummed form when possible */
  double num;
  double den;

  /* l=2 */
  /* Pade(2,1) approximant */
  num        = 5992.*Pi*sqrt_y + 2465.*nu*(28.-493.*nu* y);
  den        = 69020.*nu + 5992.*Pi*sqrt_y;
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

/** Resummed amplitudes in the general nu-dependent case.
 *  Refs:
 *  . Damour, Iyer & Nagar, PRD 79, 064004 (2009)     [theory]
 *  . Fujita & Iyer, PRD 82, 044051 (2010)            [test-mass 5.5PN]
 *  . Damour, Nagar & Bernuzzi, PRD 87, 084035 (2013) [complete information]
 */
void eob_wav_flm_v1(double x,double nu, double *rholm, double *flm)
{

  /** Coefficients */
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

    /** (2,1) */
    clm[0][1] = (-1.0535714285714286 + 0.27380952380952384 *nu);
    clm[0][2] = (-0.8327841553287982 - 0.7789824263038548  *nu + 0.13116496598639457*nu2);
    /* clm[0][3] = (2.9192806270460925  - 1.019047619047619   *el1); */
    /* clm[0][4] = (-1.28235780892213   + 1.073639455782313   *el1); */
    /* clm[0][5] = (-3.8466571723355227 + 0.8486467106683944  *el1)*PMTERMS_eps; */

    /** (2,2) */
    clm[1][1] = (-1.0238095238095237 + 0.6547619047619048*nu);
    clm[1][2] = (-1.94208238851096   - 1.5601379440665155*nu + 0.4625614134542706*nu2);
    /* clm[1][3] = (12.736034731834051  - 2.902228713904598 *nu - 1.9301558466099282*nu2 + 0.2715020968103451*nu3 - 4.076190476190476*el2); */
    /* clm[1][4] = (-2.4172313935587004 + 4.173242630385488 *el2); */
    /* clm[1][5] = (-30.14143102836864  + 7.916297736025627 *el2); */

    /** (3,1) */
    clm[2][1] = (-0.7222222222222222 - 0.2222222222222222*nu);
    clm[2][2] = (0.014169472502805836 - 0.9455667789001122*nu - 0.46520763187429853*nu2);
    /* clm[2][3] = (1.9098284139598072 - 0.4126984126984127*el1); */
    /* clm[2][4] = (0.5368150316615179 + 0.2980599647266314*el1); */
    /* clm[2][5] = (1.4497991763035063 - 0.0058477188106817735*el1)*PMTERMS_eps; */

    /** (3,2) */
    clm[3][1] = (0.003703703703703704*(328. - 1115.*nu + 320.*nu2))/(-1. + 3.*nu);
    clm[3][2] = (6.235191420376606e-7*(-1.444528e6 + 8.050045e6*nu - 4.725605e6*nu2 - 2.033896e7*nu3 + 3.08564e6*nu4))/((-1. + 3.*nu)*(-1. + 3.*nu));
    /* clm[3][3] = (6.220997955214429 - 1.6507936507936507*el2); */
    /* clm[3][4] = (-3.4527288879001268 + 2.005408583186361*el2)*PMTERMS_eps; */

    /** (3,3) */
    clm[4][1] =  (-1.1666666666666667 + 0.6666666666666666*nu);
    clm[4][2] = (-1.6967171717171716 - 1.8797979797979798*nu + 0.45151515151515154*nu2);
    /* clm[4][3] = (14.10891386831863 - 3.7142857142857144*el3); */
    /* clm[4][4] = (-6.723375314944128 + 4.333333333333333*el3); */
    /* clm[4][5] = (-29.568699895427518 + 6.302092352092352*el3)*PMTERMS_eps; */

    /** (4,1) */
    clm[5][1] = (0.001893939393939394*(602. - 1385.*nu + 288.*nu2))/(-1. + 2.*nu);
    clm[5][2] = (- 0.36778992787515513);
    /* clm[5][3] = (0.6981550175535535 - 0.2266955266955267*el1); */
    /* clm[5][4] = (-0.7931524512893319 + 0.2584672482399755*el1)*PMTERMS_eps; */

    /** (4,2) */
    clm[6][1] = (0.0007575757575757576*(1146. - 3530.*nu + 285.*nu2))/(-1. + 3.*nu);
    clm[6][2] = - (3.1534122443213353e-9*(1.14859044e8 - 2.95834536e8*nu - 1.204388696e9*nu2 + 3.04798116e9*nu3 + 3.79526805e8*nu4))/((-1. + 3.*nu)*(-1. + 3.*nu));
    /* clm[6][3] = 4.550378418934105e-12*(8.48238724511e11 - 1.9927619712e11*el2); */
    /* clm[6][4] = (-0.6621921297263365 + 0.787251738160829*el2)*PMTERMS_eps; */
    
    /** (4,3) */
    clm[7][1] = (0.005681818181818182*(222. - 547.*nu + 160.*nu2))/(-1. + 2.*nu);
    clm[7][2] = (- 0.9783218202252293);
    /* clm[7][3] = (8.519456157072423 - 2.0402597402597404*el3)*PMTERMS_eps; */
    /* clm[7][4] = (-5.353216984886716 + 2.5735094451003544*el3)*PMTERMS_eps; */

    /** (4,4) */
    clm[8][1] = (0.0007575757575757576*(1614. - 5870.*nu + 2625.*nu2))/(-1. + 3.*nu);
    clm[8][2] = (3.1534122443213353e-9*(-5.11573572e8 + 2.338945704e9*nu - 3.13857376e8*nu2 - 6.733146e9*nu3 + 1.252563795e9*nu4))/((-1. + 3.*nu)*(-1. + 3.*nu));
    /* clm[8][3] = (15.108111214795123 - 3.627128427128427*el4); */
    /* clm[8][4] = (-8.857121657199649 + 4.434988849534304*el4)*PMTERMS_eps; */

    /** (5,1) */
    clm[9][1] = (0.002564102564102564*(319. - 626.*nu + 8.*nu2))/(-1. + 2.*nu);
    clm[9][2] = (- 0.1047896120973044);
    /* clm[9][3] = (0.642701885362399 - 0.14414918414918415*el1)*PMTERMS_eps; */
    /* clm[9][4] = (-0.07651588046467575 + 0.11790664036817883*el1)*PMTERMS_eps; */

    /** (5,2) */
    clm[10][1] = (0.00007326007326007326*(-15828. + 84679.*nu - 104930.*nu2 + 21980.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[10][2] = (- 0.4629337197600934)*PMTERMS_eps; 
    /* clm[10][3] = (2.354458371550237 - 0.5765967365967366*el2)*PMTERMS_eps; */

    /** (5,3) */
    clm[11][1] = (0.002564102564102564*(375. - 850.*nu + 176.*nu2))/(-1. + 2.*nu);
    clm[11][2] = (- 0.5788010707241477);
    /* clm[11][3] = (5.733973288504755 - 1.2973426573426574*el3)*PMTERMS_eps; */
    /* clm[11][4] = (-1.9573287625526001 + 1.2474448628294783*el3)*PMTERMS_eps; */

    /** (5,4) */
    clm[12][1] = (0.00007326007326007326*(-17448. + 96019.*nu - 127610.*nu2 + 33320.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[12][2] = (- 1.0442142414362194)*PMTERMS_eps;
    /* clm[12][3] = (10.252052781721588 - 2.3063869463869464*el4)*PMTERMS_eps; */

    /** (5,5) */
    clm[13][1] = (0.002564102564102564*(487. - 1298.*nu + 512.*nu2))/(-1. + 2.*nu);
    clm[13][2] = (- 1.5749727622804546);
    /* clm[13][3] = (15.939827047208668 - 3.6037296037296036*el5)*PMTERMS_eps; */
    /* clm[13][4] = (-10.272578060123237 + 4.500041838503377*el5)*PMTERMS_eps; */
    
    /** (6,1) */
    clm[14][1] = (0.006944444444444444*(-161. + 694.*nu - 670.*nu2 + 124.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[14][2] = (- 0.29175486850885135)*PMTERMS_eps;
    /* clm[14][3] = (0.21653486654395454 - 0.10001110001110002*el1)*PMTERMS_eps; */

    /** (6,2) */
    clm[15][1] = (0.011904761904761904*(-74. + 378.*nu - 413.*nu2 + 49.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[15][2] = ( - 0.24797525070634313)*PMTERMS_eps;
    /* clm[15][3] = (1.7942694138754138 - 0.40004440004440006*el2)*PMTERMS_eps; */

    /** (6,3) */
    clm[16][1] = (0.006944444444444444*(-169. + 742.*nu - 750.*nu2 + 156.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[16][2] = (- 0.5605554442947213)*PMTERMS_eps;
    /* clm[16][3] = (4.002558222882566 - 0.9000999000999002*el3)*PMTERMS_eps; */

    /** (6,4) */
    clm[17][1] = (0.011904761904761904*(-86. + 462.*nu - 581.*nu2 + 133.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[17][2] = (- 0.7228451986855349)*PMTERMS_eps;
    /* clm[17][3] = (7.359388663371044 - 1.6001776001776002*el4)*PMTERMS_eps; */

    /** (6,5) */
    clm[18][1] = (0.006944444444444444*(-185. + 838.*nu - 910.*nu2 + 220.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[18][2] = (- 1.0973940686333457)*PMTERMS_eps;
    /* clm[18][3] = (11.623366217471297 - 2.5002775002775004*el5)*PMTERMS_eps; */

    /** (6,6) */
    clm[19][1] = (0.011904761904761904*(-106. + 602.*nu - 861.*nu2 + 273.*nu3))/(1. - 5.*nu + 5.*nu2); 
    clm[19][2] = (- 1.5543111183867486)*PMTERMS_eps;
    /* clm[19][3] = (16.645950799433503 - 3.6003996003996006*el6)*PMTERMS_eps; */

    /** (7,1) */
    clm[20][1] = (0.0014005602240896359*(-618. + 2518.*nu - 2083.*nu2 + 228.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[20][2] = ( - 0.1508235111143767)*PMTERMS_eps;
    /* clm[20][3] = (0.2581280702019663 - 0.07355557607658449*el1)*PMTERMS_eps; */

    /** (7,2) */
    clm[21][1] = (0.00006669334400426837*(16832. - 123489.*nu + 273924.*nu2 - 190239.*nu3 + 32760.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[21][2] = (- 0.351319484450667)*PMTERMS_eps;
    
    /** (7,3) */
    clm[22][1] = (0.0014005602240896359*(-666. + 2806.*nu - 2563.*nu2 + 420.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[22][2] = (- 0.37187416047628863)*PMTERMS_eps;
    /* clm[22][3] = (3.0835293524055283 - 0.6620001846892604*el3)*PMTERMS_eps; */

    /** (7,4) */
    clm[23][1] = (0.00006669334400426837*(17756. - 131805.*nu + 298872.*nu2 - 217959.*nu3 + 41076.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[23][2] = (- 0.6473746896670599)*PMTERMS_eps;
    
    /** (7,5) */
    clm[24][1] = (0.0014005602240896359*(-762. + 3382.*nu - 3523.*nu2 + 804.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[24][2] = (- 0.8269193364414116)*PMTERMS_eps;
    /* clm[24][3] = (8.750589067052443 - 1.838889401914612*el5)*PMTERMS_eps; */

    /** (7,6) */
    clm[25][1] = (0.0006002400960384153*(2144. - 16185.*nu + 37828.*nu2 - 29351.*nu3 + 6104.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[25][2] = (- 1.1403265020692532)*PMTERMS_eps;
    
    /** (7,7) */
    clm[26][1] = (0.0014005602240896359*(-906. + 4246.*nu - 4963.*nu2 + 1380.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[26][2] = (- 1.5418467934923434)*PMTERMS_eps;
    /* clm[26][3] = (17.255875091408523 - 3.6042232277526396*el7)*PMTERMS_eps; */

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

  /** Coefs with Eulerlogs */
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
    
  /** rho_lm */
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

  /** Amplitudes */
#pragma omp simd
  for (int k = 0; k < KMAX; k++) {
      flm[k] = gsl_pow_int(rholm[k], LINDEX[k]);
  }

}

/** Resummed amplitudes in the general nu-dependent case.
 *  Function introduced for higher modes
 *  (2,2) is taken at 3^{+2} PN order
 *  (4,4), (4,2) and (5,5) and l >= 7 modes are Talyor expanded 
 *  (3,2) is resummed with a Padé 32
 *  (2,1) is resummed with a Padé 51
 *  All other modes are resummed with a Padé 42
 *  Ref: arXiv:2001.09082 */
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
    
    cden = SQ(clm[k][3]) - clm[k][2]*clm[k][4];

    n1 = clm[k][1]*SQ(clm[k][3]) - clm[k][1]*clm[k][2]*clm[k][4] - clm[k][3]*clm[k][4] + clm[k][2]*clm[k][5];
    n2 = clm[k][2]*SQ(clm[k][3]) - SQ(clm[k][2])*clm[k][4] - clm[k][1]*clm[k][3]*clm[k][4] + SQ(clm[k][4]) + clm[k][1]*clm[k][2]*clm[k][5] - clm[k][3]*clm[k][5];
    n3 = SQ(clm[k][3])*clm[k][3] - 2.*clm[k][2]*clm[k][3]*clm[k][4] + clm[k][1]*SQ(clm[k][4]) + SQ(clm[k][2])*clm[k][5] - clm[k][1]*clm[k][3]*clm[k][5];
    
    d1 = - clm[k][3]*clm[k][4] + clm[k][2]*clm[k][5];
    d2 = SQ(clm[k][4]) - clm[k][3]*clm[k][5];

    rholm[k] = (cden + n1*x + n2*x2 + n3*x3)/(cden + d1*x + d2*x2);
  }

  // Padé (4,2) - Used for (3,2), (3,3), (4,1), (4,3), (5,2), (5,3), (5,4) and l=6 except (6,1) 
  const int kmaxPade42 = 12;
  int kPade42[] = {3,4,5,7,10,11,12,15,16,17,18,19};
  
  for (int i=0; i<kmaxPade42; i++) {
    int k = kPade42[i];
    
    cden = SQ(clm[k][4]) - clm[k][3]*clm[k][5];

    n1 = clm[k][1]*SQ(clm[k][4]) - clm[k][1]*clm[k][3]*clm[k][5] - clm[k][4]*clm[k][5] + clm[k][3]*clm[k][6];
    n2 = clm[k][2]*SQ(clm[k][4]) - clm[k][2]*clm[k][3]*clm[k][5] - clm[k][1]*clm[k][4]*clm[k][5] + SQ(clm[k][5])
	  + clm[k][1]*clm[k][3]*clm[k][6] - clm[k][4]*clm[k][6];
    n3 = clm[k][3]*SQ(clm[k][4]) - SQ(clm[k][3])*clm[k][5] - clm[k][2]*clm[k][4]*clm[k][5] + clm[k][1]*SQ(clm[k][5])
	  + clm[k][2]*clm[k][3]*clm[k][6] - clm[k][1]*clm[k][4]*clm[k][6];
    n4 = SQ(clm[k][4])*clm[k][4] - 2*clm[k][3]*clm[k][4]*clm[k][5] + clm[k][2]*SQ(clm[k][5]) + SQ(clm[k][3])*clm[k][6]
	  - clm[k][2]*clm[k][4]*clm[k][6];
    d1 = - clm[k][4]*clm[k][5] + clm[k][3]*clm[k][6];
    d2 = SQ(clm[k][5]) - clm[k][4]*clm[k][6];

    rholm[k] = (cden + n1*x + n2*x2 + n3*x3 + n4*x4)/(cden + d1*x + d2*x2);
  }

  // Padé (5,1) - Used for (2,1)
  const int kmaxPade51 = 1;
  int kPade51[] = {0};
  
  for (int i=0; i<kmaxPade51; i++) {
    int k = kPade51[i];
    
    n1 = clm[k][1]*clm[k][5] -           clm[k][6];
    n2 = clm[k][2]*clm[k][5] - clm[k][1]*clm[k][6];
    n3 = clm[k][3]*clm[k][5] - clm[k][2]*clm[k][6];
    n4 = clm[k][4]*clm[k][5] - clm[k][3]*clm[k][6];
    n5 = clm[k][5]*clm[k][5] - clm[k][4]*clm[k][6];

    d1 = - clm[k][6];
      
    rholm[k] = (clm[k][5] + n1*x + n2*x2 + n3*x3 + n4*x4 + n5*x5)/(clm[k][5] + d1*x);
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
    
    /* Note: the two sums give different result */
#if (1)
    rholm[k] = clm[k][0];
    for (int n=1; n<7; n++) {	
      rholm[k] += clm[k][n]*xn[n];
    }
#else  
    rholm[k] = x6*clm[k][6];
    for (int n=6; n-- >1; ) { // 4,3,2,1 //
      rholm[k] += clm[k][n]*xn[n];
    } 
    rholm[k] += clm[k][0];
#endif
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
  
/** Resummed amplitudes for the spin case. 
    This function computes the residual amplitude corrections flm's as 
    introduced in Damour, Iyer & Nagar, PRD 79, 064004 (2008).
    The orbital part is taken at the usual 3^{+2} PN order, i.e. 3PN terms
    are integrated by the 4PN and 5PN test-particle terms, with the higher
    modes obtained by Fujita & Iyer.
    It only includes spin-spin interaction at LO for the (2,2) mode.
    Note that the variables called here (a1,a2)
    are what we usually cal tilde{a}_1 and tilde{a}_2 and are defined as
    a1 = X1*chi1, a2=X2*chi2 and are passed here as parameters. Special 
    combinations of these quantities are used here to write the spin-dependent
    part of the waveform in particularly compact form, so that the (spinning)
    test-particle limit is recovered just by visual inspection of the equations */
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
    Resummed amplitudes for the spin case. 
    This function computes the residual amplitude corrections flm's as 
    introduced in Damour, Iyer & Nagar, PRD 79, 064004 (2008).
    The orbital part is taken at the usual 3^{+2} PN order, i.e. 3PN terms
    are integrated by the 4PN and 5PN test-particle terms, with the higher
    modes obtained by Fujita & Iyer.
    The function includes spin-spin interaction at NLO for the (2,2) mode
    and at LO for the (2,1),(3,1) and (3,3) modes. 
    Note that the variables called here (a1,a2)
    are what we usually cal tilde{a}_1 and tilde{a}_2 and are defined as
    a1 = X1*chi1, a2=X2*chi2 and are passed here as parameters. */
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
    Resummed amplitudes for the spin case. 
    This function computes the residual amplitude corrections flm's as 
    introduced in Damour, Iyer & Nagar, PRD 79, 064004 (2008).
    The m=even modes have the usual structure (rho_lm^orb + rho_lm^spin)^l.
    The m=odd modes are factorized as f_lm^orb * f_lm^spin.
    Some of the f_lm^spin are inverse-resummed
    Function introduced for higher modes.
    Ref: arXiv:2001.09082
    Note that the variables called here (a1,a2)
    are what we usually cal tilde{a}_1 and tilde{a}_2 and are defined as
    a1 = X1*chi1, a2=X2*chi2 and are passed here as parameters. */
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

/** Calculate tidal correction to multipolar waveform amplitude
    Ref. Damour, Nagar & Villain, Phys.Rev. D85 (2012) 123007 */
void eob_wav_hlmTidal(double x, Dynamics *dyn, double *hTidallm)
{
  const double nu       = dyn->nu;
  const double XA       = dyn->X1;
  const double XB       = dyn->X2;
  const double khatA_2  = dyn->khatA2;
  const double khatB_2  = dyn->khatB2;
  const double kapA2j   = dyn->kapA2j;
  const double kapB2j   = dyn->kapB2j;
  const double kapT2j   = dyn->kapT2j;
  
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

#if(USEGRAVITOMAGNETICTERMS)
  
  if ( (dyn->use_tidal_gravitomagnetic==TIDES_GM_GSF) || (dyn->use_tidal_gravitomagnetic==TIDES_GM_PN) ) {
    const double fourtnine= 1.5555555555555555556;  // 14/9 = 112/(3*24)
    const double fourthird = 1.3333333333333333333; // 32/24 = 4/3
    hTidallm[0] += 0.5*( -1.*kapA2j/XB + kapB2j/XA )*x5;
    hTidallm[1] += fourtnine*kapT2j*x6;
    hTidallm[2] += 0.5*( kapA2j*(4. - 17.*XB) - kapB2j*(4. - 17.*XA) )*x6;
    hTidallm[3] += fourthird*kapT2j*x5/(1.-3.*nu);
    hTidallm[4] += 0.5*( kapA2j*(4. - 9.*XB) - kapB2j*(4. - 9.*XA) )*x6;
  }

#endif

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


/** Computes the factors and the coefficients that build the  
    NQC corrections to the waveform in the spinning case */
void eob_wav_hlmNQC_find_a1a2a3(Dynamics *dyn, Waveform_lm *h, Waveform_lm *hnqc)
{
  double A_tmp, dA_tmp, omg_tmp, domg_tmp;
  double alpha1[KMAX], omega1[KMAX];
  double c1A[KMAX], c2A[KMAX], c3A[KMAX], c4A[KMAX];
  double c1phi[KMAX], c2phi[KMAX], c3phi[KMAX], c4phi[KMAX];
	    
  const double nu   = dyn->nu;
  const double chi1 = dyn->chi1;
  const double chi2 = dyn->chi2;
  const double X1   = dyn->X1;
  const double X2   = dyn->X2;
  const double aK   = dyn->a1+dyn->a2;
  const double Mbh  = dyn->Mbhf;
  const double abh  = dyn->abhf;
  const double ecc  = dyn->ecc;

  FILE* fp;
    
  double *t       = h->time;
  double *r       = dyn->data[EOB_RAD];
  double *w       = dyn->data[EOB_MOMG]; /* Omega */
  double *pph     = dyn->data[EOB_PPHI];
  double *pr_star = dyn->data[EOB_PRSTAR];
  double *Omg_orb = dyn->data[EOB_OMGORB]; /* Omega orbital */
  double *ddotr   = dyn->data[EOB_DDOTR];
  
  double P[2], M[4];
  double max_A[KMAX],max_dA[KMAX],d2max[KMAX],d3max[KMAX],max_omg[KMAX],max_domg[KMAX],maxd2omg[KMAX], DeltaT[KMAX];
  double ai[KMAX][2];
  double bi[KMAX][2];
  
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
  
  if (EOBPars->use_flm == USEFLM_HM) {
    
    /* Higher modes */
    /* 21, 32, 42, 43 and 44 extracted from postpeak */
    int K_HM[5] = {0,3,6,7,8};
	
    QNMHybridFitCab_HM(nu, X1, X2, chi1, chi2, aK,  Mbh, abh,  
		       c1A, c2A, c3A, c4A, c1phi, c2phi, c3phi, c4phi,
		       alpha1, omega1);
	    
    for (int j=0; j<5; j++) {
      int k = K_HM[j];
      
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
    
    /* 22, 31, 33, 41 and 55 fitted directly + 44 dA */
    eob_nqc_point_HM(dyn, max_A, max_dA, max_omg, max_domg);
	    
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
	w2       = SQ(w[j]); //CHECKME: Omg or Omg_orbital ?
	n1[k][j]  = pr_star2/(r2*w2);         /* [pr*\/(r Omg)]^2 */
	n2[k][j]  = ddotr[j]/(r[j]*w2);       /* [ddot{r}/(r Omg^2)] */
	n4[k][j]  = pr_star[j]/(r[j]*w[j]);   /* pr*\/(r Omg) */
	n5[k][j]  = n4[k][j]*r2*w2;              /* (pr*)*(r Omg) */
      }
    }
  }

  if (EOBPars->use_flm == USEFLM_HM) {
    for (int j=0; j<size; j++) {
      /* l=2,m=1 */
      n2[0][j] = cbrt(SQ(w[j]))*n1[0][j];
      n5[0][j] = cbrt(SQ(w[j]))*n4[0][j];
      /* l=3 & l=4 */
      for (int k=2; k<14; k++) {   
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

  if (EOBPars->use_flm == USEFLM_HM) {
    
    eob_nqc_deltat_lm(dyn, dtmrg);
    
    for (int k=0; k<KMAX; k++) {   
      if(h->kmask[k]) {
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
    if(h->kmask[k]) {  
      
      if (EOBPars->use_flm == USEFLM_HM) {
	jmax = j_NQC[k];
      }
      
      ai[k][0] = ai[k][1] = 0.;
      bi[k][0] = bi[k][1] = 0.;
      
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

  
  if (ecc != 0.) {
    double t0 = tNQC - 30.;
    double alpha = 0.09;
    double *smooth_theta;
    smooth_theta = (double*) calloc (size, sizeof(double));
   
    for (int j=0; j<size; j++) {
      smooth_theta[j] = 1./(1. + exp(-alpha*(t[j] - t0)));

      for (int k=0; k<KMAX; k++) {
	if(hnqc->kmask[k]){
	  n1[k][j] = n1[k][j]*smooth_theta[j];
	  n2[k][j] = n2[k][j]*smooth_theta[j];
	  n4[k][j] = n4[k][j]*smooth_theta[j];
	  n5[k][j] = n5[k][j]*smooth_theta[j];
	}
      }
    }
    
    free(smooth_theta);
  }

  /** Set amplitude and phase */
  for (int k=0; k<KMAX; k++) {
    if(h->kmask[k]) {  
      for (int j=0; j<size; j++) {
	hnqc->ampli[k][j] = 1. + ai[k][0]*n1[k][j] + ai[k][1]*n2[k][j];
	hnqc->phase[k][j] =      bi[k][0]*n4[k][j] + bi[k][1]*n5[k][j];
      }
    }
  }
  
  /** Multiply waveform to NQC */
  for (int k=0; k<KMAX; k++) {
    if(h->kmask[k]){    
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
      fprintf(fp, "%d %d %d %e %e %e %e\n", k, LINDEX[k], MINDEX[k], 
	      ai[k][0], ai[k][1], 
	      bi[k][0], bi[k][1]);
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

/** Computes the factors and the coefficients that build the  
    NQC corrections to the waveform in the spinning case. 
    This routine works around merger with dyn and h and 
    then add everything also to hlm */
void eob_wav_hlmNQC_find_a1a2a3_mrg_HM(Dynamics *dyn_mrg, Waveform_lm *hlm_mrg, Waveform_lm *hnqc,
				       Dynamics *dyn, Waveform_lm *hlm)
{

  double A_tmp, dA_tmp, omg_tmp, domg_tmp;
  double alpha1[KMAX], omega1[KMAX];
  double c1A[KMAX], c2A[KMAX], c3A[KMAX], c4A[KMAX];
  double c1phi[KMAX], c2phi[KMAX], c3phi[KMAX], c4phi[KMAX];
	  
  const double nu   = dyn->nu;
  const double chi1 = dyn->chi1;
  const double chi2 = dyn->chi2;
  const double X1   = dyn->X1;
  const double X2   = dyn->X2;
  const double aK   = dyn->a1+dyn->a2;
  const double Mbh  = dyn->Mbhf;
  const double abh  = dyn->abhf;
  const double ecc  = dyn->ecc;

  FILE* fp;
    
  double *t       = hlm_mrg->time;
  double *r       = dyn_mrg->data[EOB_RAD];
  double *w       = dyn_mrg->data[EOB_MOMG]; /* Omega */
  double *pph     = dyn_mrg->data[EOB_PPHI];
  double *pr_star = dyn_mrg->data[EOB_PRSTAR];
  double *Omg_orb = dyn_mrg->data[EOB_OMGORB]; /* Omega orbital */
  double *ddotr   = dyn_mrg->data[EOB_DDOTR];
  
  double P[2], M[4];   
  double max_A[KMAX],max_dA[KMAX],d2max[KMAX],d3max[KMAX],max_omg[KMAX],max_domg[KMAX],maxd2omg[KMAX], DeltaT[KMAX];
  double ai[KMAX][2];
  double bi[KMAX][2];
  
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
    if(hlm_mrg->kmask[k]){
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
  /* 21, 32, 42, 43 and 44 extracted from postpeak */
  int K_HM[5] = {0,3,6,7,8};

  QNMHybridFitCab_HM(nu, X1, X2, chi1, chi2, aK,  Mbh, abh,  
		     c1A, c2A, c3A, c4A, c1phi, c2phi, c3phi, c4phi,
		     alpha1, omega1);
  
  for (int j=0; j<5; j++) {
    int k = K_HM[j];
    
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
  
  /* 22, 31, 33, 41 and 55 fitted directly + 44 dA */
  eob_nqc_point_HM(dyn, max_A, max_dA, max_omg, max_domg);
	    
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
    if(hlm_mrg->kmask[k]){ 
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
    if(hlm_mrg->kmask[0]){
      n2[0][j] = cbrt(SQ(w[j]))*n1[0][j];
      n5[0][j] = cbrt(SQ(w[j]))*n4[0][j];
    }
    /* l=3 & l=4 */
    for (int k=2; k<14; k++) {
      if(hlm_mrg->kmask[k]){
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
    if(hlm_mrg->kmask[k]){   
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
    if(hlm_mrg->kmask[k]){
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
    if(hlm_mrg->kmask[k]){
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
    if(hlm_mrg->kmask[k]){
      for (int j=0; j<size; j++) {
        m11[k][j] = n1[k][j] * p1tmp[k][j];
        m12[k][j] = n2[k][j] * p1tmp[k][j];
      }
    }
  }

  /* Take FD derivatives */
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask[k]){
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
    if(hlm_mrg->kmask[k]){

      jmax = j_NQC[k];
      
      ai[k][0] = ai[k][1] = 0.;
      bi[k][0] = bi[k][1] = 0.;

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

  if (ecc != 0.) {
    double t0 = tNQC - 30.;
    double alpha = 0.09;
    double *smooth_theta;
    smooth_theta = (double*) calloc (size, sizeof(double));

    for (int j=0; j<size; j++) {
      smooth_theta[j] = 1./(1. + exp(-alpha*(t[j] - t0)));

      for (int k=0; k<KMAX; k++) {
	if(hnqc->kmask[k]){
	  n1[k][j] = n1[k][j]*smooth_theta[j];
	  n2[k][j] = n2[k][j]*smooth_theta[j];
	  n4[k][j] = n4[k][j]*smooth_theta[j];
	  n5[k][j] = n5[k][j]*smooth_theta[j];
	}
      }
    }
    
    free(smooth_theta);
  }
  
  /** Set amplitude and phase */
  for (int k=0; k<KMAX; k++) {
    if(hlm_mrg->kmask[k]){
      for (int j=0; j<size; j++) {
        hnqc->ampli[k][j] = 1. + ai[k][0]*n1[k][j] + ai[k][1]*n2[k][j];
        hnqc->phase[k][j] =      bi[k][0]*n4[k][j] + bi[k][1]*n5[k][j];
      }
    }
  }
  
  /** Multiply merger waveform to NQC */
  for (int k=0; k<KMAX; k++) {
    for (int j=0; j<size; j++) {
      if(hlm_mrg->kmask[k]){
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
    if(hlm_mrg->kmask[k]){   
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
    if(hlm_mrg->kmask[0]){
      n2[0][j] = cbrt(SQ(w[j]))*n1[0][j];
      n5[0][j] = cbrt(SQ(w[j]))*n4[0][j];
    }
    /* l=3 & l=4 */
    for (int k=2; k<14; k++) {   
      if(hlm_mrg->kmask[k]){   
	n5[k][j]  = cbrt(SQ(w[j]))*n4[k][j];
      }
    }
  }

  if (ecc != 0.) {
    double t0 = tNQC - 30.;
    double alpha = 0.09;
    double *smooth_theta;
    smooth_theta = (double*) calloc (fullsize, sizeof(double));
    
    for (int j=0; j<fullsize; j++) {
      smooth_theta[j] = 1./(1. + exp(-alpha*(hlm->time[j] - t0)));

      for (int k=0; k<KMAX; k++) {
	if(hlm->kmask[k]){
	  n1[k][j] = n1[k][j]*smooth_theta[j];
	  n2[k][j] = n2[k][j]*smooth_theta[j];
	  n4[k][j] = n4[k][j]*smooth_theta[j];
	  n5[k][j] = n5[k][j]*smooth_theta[j];
	}
      }
    }

    free(smooth_theta);
  }

  for (int k=0; k<KMAX; k++) {
    if(hlm->kmask[k]){
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
      if(hlm->kmask[k]){
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

/** Computes the factors and the coefficients that build the  
    NQC corrections to the waveform in the spinning case. 
    This routine works around merger with dyn and h and 
    then add everything also to hlm 
    This is specific for the 22 mode */
void eob_wav_hlmNQC_find_a1a2a3_mrg_22(Dynamics *dyn_mrg, Waveform_lm *hlm_mrg, Waveform_lm *hnqc,
				       Dynamics *dyn, Waveform_lm *hlm)
{
#define KMAX_22 (2) // 0, 1 -> (2,1), (2,2) // k=0 unused.

  double A_tmp, dA_tmp, omg_tmp, domg_tmp;
  double alpha1[KMAX_22], omega1[KMAX_22];
  double c1A[KMAX_22], c2A[KMAX_22], c3A[KMAX_22], c4A[KMAX_22];
  double c1phi[KMAX_22], c2phi[KMAX_22], c3phi[KMAX_22], c4phi[KMAX_22];
	  
  const double nu   = dyn->nu;
  const double chi1 = dyn->chi1;
  const double chi2 = dyn->chi2;
  const double X1   = dyn->X1;
  const double X2   = dyn->X2;
  const double aK   = dyn->a1+dyn->a2;
  const double Mbh  = dyn->Mbhf;
  const double abh  = dyn->abhf;
  const double ecc  = dyn->ecc;

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
    n5[k22][j]  = n4[k22][j]*r2*w2;              /* (pr*)*(r Omg) */
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

  if (ecc != 0.) {
    double t0 = tNQC - 30.;
    double alpha = 0.09;
    double *smooth_theta;
    smooth_theta = (double*) calloc (size, sizeof(double));

    for (int j=0; j<size; j++) {
      smooth_theta[j] = 1./(1. + exp(-alpha*(t[j] - t0)));
      
      n1[k22][j] = n1[k22][j]*smooth_theta[j];
      n2[k22][j] = n2[k22][j]*smooth_theta[j];
      n4[k22][j] = n4[k22][j]*smooth_theta[j];
      n5[k22][j] = n5[k22][j]*smooth_theta[j];
    }
    
    free(smooth_theta);
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

  if (ecc != 0.) {
    double t0 = tNQC - 30.;
    double alpha = 0.09;
    double *smooth_theta;
    smooth_theta = (double*) calloc (fullsize, sizeof(double));
    
    for (int j=0; j<fullsize; j++) {
      smooth_theta[j] = 1./(1. + exp(-alpha*(hlm->time[j] - t0)));

      n1[k22][j] = n1[k22][j]*smooth_theta[j];
      n2[k22][j] = n2[k22][j]*smooth_theta[j];
      n4[k22][j] = n4[k22][j]*smooth_theta[j];
      n5[k22][j] = n5[k22][j]*smooth_theta[j];
    }

    free(smooth_theta);
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

/** NQC corrections to the RWZ multipolar waveform
    Nagar, Damour, Reisswig, Pollney http://arxiv.org/abs/1506.08457 
    Nonspinning case, Current fits: 9/02/2016 */
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

/** Generic routine for NQC */
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

/** Ringdown waveform template */
void eob_wav_ringdown_template(double x, double a1, double a2, double a3, double a4, double b1, double b2, double b3, double b4, double sigmar, double sigmai, double *psi)
{  
  double amp   = ( a1 * tanh(a2*x +a3) + a4 ) ;
  double phase = -b1*log((1. + b3*exp(-b2*x) + b4*exp(-2.*b2*x))/(1.+b3+b4));   
  psi[0] = amp * exp(-sigmar*x); /* amplitude */
  psi[1] = - (phase - sigmai*x); /* phase, minus sign in front by convention */
}

/** Ringdown calculation and match to the dynamics */ 
void eob_wav_ringdown_v1(Dynamics *dyn, Waveform_lm *hlm)
{

  const double Mbh   = dyn->Mbhf;
  const double abh   = dyn->abhf;
  const double nu    = dyn->nu;
  const double q     = dyn->q;
  const double chi1  = dyn->chi1;
  const double chi2  = dyn->chi2;
  const double X1    = dyn->X1;
  const double X2    = dyn->X2;
  const double aK    = dyn->a1+dyn->a2;

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
  
  /* This is a hard-fix that always guarantee the 7 points */
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
	
/** Ringdown calculation and match to the dynamics 
    This refers to the higher modes paper: arXiv:2001.09082 
*/ 
void eob_wav_ringdown_HM(Dynamics *dyn, Waveform_lm *hlm)
{
  
  const double Mbh   = dyn->Mbhf;
  const double abh   = dyn->abhf;
  const double nu    = dyn->nu;
  const double q    = dyn->q;
  const double chi1  = dyn->chi1;
  const double chi2  = dyn->chi2;
  const double X1    = dyn->X1;
  const double X2    = dyn->X2;
  const double aK    = dyn->a1+dyn->a2;
	
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
      idx[k] = size-1;
      for (int j = size-1; j-- ; ) {  
	if (t[j] * ooMbh < tmrg[k]) {
	  break;
	}
	idx[k] = j;
      }
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
	
      eob_wav_ringdown_template(t0, a1[k], a2[k], a3[k], a4[k], b1[k], b2[k], b3[k], b4[k], sigma[0][k], sigma[1][k], psi);
      Deltaphi[k] = psi[1] - hlm->phase[k][index_rng];
      
      /* Compute and attach ringdown */
      for (int j = index_rng-1; j < size ; j++ ) {
	
	tm = t_lm[k][j] - tmatch[k];
	
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

/** Main routine for factorized EOB waveform */
void eob_wav_hlm(Dynamics *dyn, Waveform_lm_t *hlm)
{

  const double nu = dyn->nu;  
  const double chi1 = dyn->chi1;  
  const double chi2 = dyn->chi2;  
  const double a1 = dyn->a1;  
  const double a2 = dyn->a2;  
  const double X1 = dyn->X1;  
  const double X2 = dyn->X2;  
  const double C_Q1 = dyn->C_Q1;  
  const double C_Q2 = dyn->C_Q2;  
  const double ecc = dyn->ecc;  
  const double r_hyp = dyn->r_hyp;  
  const int usetidal = dyn->use_tidal;
  const int usespins = dyn->use_spins;
  const int usespeedytail = EOBPars->use_speedytail;
  const double X12 = X1-X2; /* sqrt(1-4nu) */

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

    if ((STREQUAL(use_flm_opt[EOBPars->use_flm],"HM"))) {
    } else {
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
  if (!(STREQUAL(nqc_hlm_opt[EOBPars->nqc_coefs_hlm],"none")) &&
      !(STREQUAL(nqc_hlm_opt[EOBPars->nqc_coefs_hlm],"compute"))) {
    /** Add NQC correction */    
    Waveform_lm_t hNQC; 
    /* eob_wav_hlmNQC_nospin2016(nu,r,prstar,Omega,ddotr, &hNQC); */ 
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

/** Routine for factorized EOB waveform and eccentric systems */
void eob_wav_hlm_ecc(Dynamics *dyn, Waveform_lm_t *hlm)
{
  
  const double nu = dyn->nu;  
  const double chi1 = dyn->chi1;  
  const double chi2 = dyn->chi2;  
  const double a1 = dyn->a1;  
  const double a2 = dyn->a2;  
  const double X1 = dyn->X1;  
  const double X2 = dyn->X2;  
  const double C_Q1 = dyn->C_Q1;  
  const double C_Q2 = dyn->C_Q2;  
  const double ecc = dyn->ecc;  
  const double r_hyp = dyn->r_hyp;  
  const int usetidal = dyn->use_tidal;
  const int usespins = dyn->use_spins;
  const int usespeedytail = EOBPars->use_speedytail;
  const double X12 = X1-X2; /* sqrt(1-4nu) */

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
  const double rdot  = dyn->rdot;
  const double r2dot = dyn->r2dot;
  const double r3dot = dyn->r3dot;
  const double r4dot = dyn->r4dot;
  const double r5dot = dyn->r5dot;
  const double Omegadot  = dyn->Omegadot;
  const double Omega2dot = dyn->Omega2dot;
  const double Omega3dot = dyn->Omega3dot;
  const double Omega4dot = dyn->Omega4dot;
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
  eob_wav_hlmNewt_ecc(rw, rdot, r2dot, r3dot, r4dot, r5dot,
		      Omega, Omegadot, Omega2dot, Omega3dot, Omega4dot,
		      phi, nu, &hNewt);
  
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

    if (!(EOBPars->use_flm == USEFLM_HM)) {
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

/** Old flux routines */

/** Resummed amplitudes in the general nu-dependent case.
 *  Refs:
 *  . Damour, Iyer & Nagar, PRD 79, 064004 (2009)     [theory]
 *  . Fujita & Iyer, PRD 82, 044051 (2010)            [test-mass 5.5PN]
 *  . Damour, Nagar & Bernuzzi, PRD 87, 084035 (2013) [complete information]
 */
/* Old version
   - coefficients computed all the times.
   - (3,1) and (3,3) mode do not contain nu-corrections.
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

/** Resummed amplitudes for the spin case. 
    This function computes the residual amplitude corrections flm's as 
    introduced in Damour, Iyer & Nagar, PRD 79, 064004 (2008).
    The orbital part is taken at the usual 3^{+2} PN order, i.e. 3PN terms
    are integrated by the 4PN and 5PN test-particle terms, with the higher
    modes obtained by Fujita & Iyer.
    Note that the variables called here (a1,a2)
    are what we usually cal tilde{a}_1 and tilde{a}_2 and are defined as
    a1 = X1*chi1, a2=X2*chi2 and are passed here as parameters. Special 
    combinations of these quantities are used here to write the spin-dependent
    part of the waveform in particularly compact form, so that the (spinning)
    test-particle limit is recovered just by visual inspection of the equations 
*/
/* Old version
   - repeats the rholm computation for the orbital part.
   - (3,1) and (3,3) mode do not contain nu-corrections.
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
