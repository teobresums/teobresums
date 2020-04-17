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

/** Coefficients for Newtonian flux */
/*
static const double CNlm[35] = {
  8./45, 32./5,
  1./1260, 32./63, 243./28,
  1./44100, 32./3969, 729./700, 8192./567,
  1./19958400, 256./400950, 2187./70400, 131072./66825, 1953125./76032,
  1./1123782660, 128./28667925, 59049./15415400, 4194304./47779875, 48828125./13621608, 839808./17875,
  1./9.3498717312e11, 32./1.35270135e8, 1594323./3.2064032e10, 4194304./3.07432125e8, 1220703125./5.666588928e9, 5668704./875875, 96889010407./1111968000,
  1./8.174459284992e13, 32./3.4493884425e10, 177147./3.96428032e10, 4194304./1.5679038375e10, 30517578125./8.00296713216e11, 51018336./1.04229125e8, 4747561509943./4.083146496e11, 274877906944./1688511825.
};
*/
static const double CNlm[35] = {
  0.17777777777777778, 6.4,
  0.0007936507936507937, 0.5079365079365079, 8.678571428571429,
  2.2675736961451248e-05, 0.008062484252960444, 1.0414285714285714, 14.447971781305114,
  5.010421677088344e-08, 0.0006384836014465644, 0.03106534090909091, 1.9614216236438458, 25.688197074915823,
  8.898517797026696e-10, 4.464920289836115e-06, 0.003830520129221428, 0.08778390483440988, 3.584607999290539, 46.98226573426573,
  1.0695333890657086e-12, 2.3656367312710969e-07, 4.972309783123969e-05, 0.013643024456211269, 0.21542115380351798, 6.472046810332524, 87.1329124642076,
  1.2233225038333268e-14, 9.27700678929842e-10, 4.468579053461083e-06, 0.0002675102834551229, 0.03813282951314997, 0.4894825318738884, 11.627213264559293, 162.79300083906728
};

/** Newtonian partial fluxes */
void eob_flx_FlmNewt(double x, double nu, double *Nlm)
{
  
  /** Shorthands*/
  const double nu2 = nu*nu;
  const double nu3 = nu2*nu;
  const double x5  = x*x*x*x*x;
  const double x6  = x*x5;
  const double x7  = x*x6;
  const double x8  = x*x7;
  const double x9  = x*x8;
  const double x10 = x*x9;
  const double x11 = x*x10;
  const double x12 = x*x11;
  
  const double sp2 = 1.-4.*nu;
  const double sp4 = (1-4*nu)*SQ((1-2*nu));
  const double sp3 = (1.-3.*nu)*(1.-3.*nu);
  const double sp5 = (1.-5.*nu+5.*nu2)*(1.-5.*nu+5.*nu2);
  const double sp6 = (1-4*nu)*(3*nu2-4*nu +1)*(3*nu2-4*nu +1);
  const double sp7 = (1 - 7*nu + 14*nu2 - 7*nu3)*(1 - 7*nu + 14*nu2 - 7*nu3);
  const double sp8 = (1 - 4*nu)*(1 - 6*nu + 10*nu2 - 4*nu3)*(1 - 6*nu + 10*nu2 - 4*nu3);

  double spx[] = {
    sp2 * x6, x5,
    sp2 * x6, sp3 * x7, sp2 * x6,
    sp4 * x8, sp3 * x7, sp4 * x8, sp3 * x7,
    sp4 * x8, sp5 * x9, sp4 * x8, sp5 * x9, sp4 * x8,
    sp6 * x10, sp5 * x9, sp6 * x10, sp5 * x9, sp6 * x10, sp5 * x9,
    sp6 * x10, sp7 * x11, sp6 * x10, sp7 * x11, sp6 * x10, sp7 * x11, sp6 * x10,
    sp8 * x12, sp7 * x11, sp8 * x12, sp7 * x11, sp8 * x12, sp7 * x11, sp8 * x12, (7*nu3-14*nu2+7*nu-1)*(7*nu3-14*nu2+7*nu-1) * x11
  };

  /** Newtonian partial fluxes*/
  for (int k = 0; k < KMAX; k++) {
    Nlm[k] = CNlm[k] * spx[k];
  }

}

/** Tail term (modulus) */
void eob_flx_Tlm(const double w, double *MTlm)
{
  double hhatk, x2, y, prod, fl;
  for (int k = 0; k < KMAX; k++) {
    hhatk = MINDEX[k] * w;
    x2    = 4.*hhatk*hhatk;
    prod  = 1.;
    for (int j=1; j <= LINDEX[k]; j++) {
      prod *= ( j*j + x2 );
    }
    y  = 4.*Pi*hhatk;
    y /= ( 1. - exp(-y) );
    MTlm[k] = sqrt( 1./( SQ(fact(LINDEX[k])) ) * y * prod );
  }
}

/** Compute horizon-absorbed fluxes. no spin case.
 * Nagar & Akcay, PRD 85, 044025 (2012)
 * Bernuzzi, Nagar & Zenginoglu, PRD 86, 104038 (2012)
 */
double eob_flx_HorizonFlux(double x, double Heff, double jhat, double nu)
{
  double rhoHlm[2]; /* only 21,22 multipoles -> k=0,1 */
  double FlmHLO[2];
  double FlmH[2];
  
  /** Shorthands */
  double nu2 = nu*nu;
  double nu3 = nu*nu2;
  double x2  = x*x;
  double x3  = x*x2;
  double x4  = x*x3;
  double x5  = x*x4;
  double x9  = x4*x5;
  double x10 = x*x9;
    
  /** The Newtonian asymptotic contribution */
  const double FNewt22 = 32./5.*x5;
  
  /** Compute leading-order part (nu-dependent) */
  FlmHLO[1] = 32./5.*(1-4*nu+2*nu2)*x9;
  FlmHLO[0] = 32./5.*(1-4*nu+2*nu2)*x10;
    
  /** Compute rho_lm */
  double c1[2];
  double c2[2];
  double c3[2];
  double c4[2];
    
  c1[0] = 0.58121;
  c2[0] = 1.01059;
  c3[0] = 7.955729;
  c4[0] = 1.650228;
  
  c1[1] = (4.-21.*nu + 27.*nu2 - 8.*nu3)/(4.*(1.-4.*nu+2.*nu2));
  c2[1] =  4.78752;
  c3[1] = 26.760136;
  c4[1] = 43.861478;
    
  rhoHlm[1] = 1. + c1[1]*x + c2[1]*x2 + c3[1]*x3 + c4[1]*x4;
  rhoHlm[0] = 1. + c1[0]*x + c2[0]*x2 + c3[0]*x3 + c4[0]*x4;
    
  /** Compute horizon multipolar flux (only l=2) */
  const double Heff2 = Heff*Heff;
  const double jhat2 = jhat*jhat;
  
  FlmH[1] = FlmHLO[1] * Heff2 * gsl_pow_int(rhoHlm[1],4);
  FlmH[0] = FlmHLO[0] * jhat2 * gsl_pow_int(rhoHlm[0],4);
    
  /** Sum over multipoles and normalize to the 22 Newtonian multipole */
  double hatFH = (FlmH[0]+FlmH[1])/FNewt22;
  
  return hatFH;
}

/** Compute horizon-absorbed fluxes. spin case. */
double eob_flx_HorizonFlux_s(double x, double Heff, double jhat, double nu, double X1, double X2, double chi1, double chi2)
{
    
  double x2 = x*x;
  double x3 = x2*x;
  double x4 = x3*x;
  double x5 = x4*x;
  double v5 = sqrt(x5);
  
  double cv5[2];
  double cv8[2];
  
  /** Coefficients of the v^5 term (Alvi leading order) */
  cv5[0] = -1./4.*chi1*(1.+3.*chi1*chi1)*X1*X1*X1;
  cv5[1] = -1./4.*chi2*(1.+3.*chi2*chi2)*X2*X2*X2;
  
  /** Coefficients of the v^8=x^4 term */
  cv8[0] = 0.5*(1.+sqrt(1.-chi1*chi1))*(1.+3.*chi1*chi1)*X1*X1*X1*X1;
  cv8[1] = 0.5*(1.+sqrt(1.-chi2*chi2))*(1.+3.*chi2*chi2)*X2*X2*X2*X2;
  
  double FH22_S = (cv5[0]+cv5[1])*v5;
  double FH22   = (cv8[0]+cv8[1])*x4;
  double FH21   =  0.0;
  
  /** Newton-normalized horizon flux: use only l=2 fluxes */
  double hatFH  = FH22_S + FH22 + FH21;
    
  return hatFH;
}

/** Flux calculation for Newton-Normalized energy flux
    Use the DIN resummation procedure.
    Add non-QC and non-K corrections to (2,2) partial flux. */
double eob_flx_Flux(double x, double Omega, double r_omega, double E, double Heff, double jhat, double r, double pr_star, double ddotr, Dynamics *dyn)
{
  return eob_flx_Flux_s(x, Omega, r_omega, E, Heff, jhat, r, pr_star, ddotr,dyn);
}

/** Flux calculation for spinning systems */
//FIXME: NQC are not applied in spin case!
double eob_flx_Flux_s(double x, double Omega, double r_omega, double E, double Heff, double jhat, double r, double pr_star, double ddotr, Dynamics *dyn)
{
  const double nu = dyn->nu;
  const double chi1 = dyn->chi1;
  const double chi2 = dyn->chi2;
  const double X1 = dyn->X1;
  const double X2 = dyn->X2;
  const double a1 = dyn->a1;
  const double a2 = dyn->a2;
  const double C_Q1 = dyn->C_Q1;
  const double C_Q2 = dyn->C_Q2;
  const double X12  = X1-X2; /* sqrt(1-4nu) */
  const double X12sq = SQ(X12); /* (1-4nu) */

  const int usetidal = dyn->use_tidal;
  const int usespins = dyn->use_spins;
  
  double prefact[] = {
    jhat, Heff,
    Heff, jhat, Heff,
    jhat, Heff, jhat, Heff,
    Heff, jhat, Heff, jhat, Heff,
    jhat, Heff, jhat, Heff, jhat, Heff,
    Heff, jhat, Heff, jhat, Heff, jhat, Heff,
    jhat, Heff, jhat, Heff, jhat, Heff, jhat, Heff};
  
  double FNewt22, sum_k=0.;
  double rholm[KMAX], flm[KMAX], FNewtlm[KMAX], MTlm[KMAX], hlmTidal[KMAX], hlmNQC[KMAX];
  double Modhhatlm[KMAX];

  /** Newtonian flux */
  eob_flx_FlmNewt(x, nu, FNewtlm);

  /* Correct amplitudes for specific multipoles and cases */
  if (usespins) {
    /* Correct (2,1), (3,1) and (3,3) ( sp2 = 1 ) */
    double x6 = gsl_pow_int(x, 6);
    FNewtlm[0] = CNlm[0] * x6; /* (2,1) */
    FNewtlm[2] = CNlm[2] * x6; /* (3,1) */
    FNewtlm[4] = CNlm[4] * x6; /* (3,3) */
    /* Correct (4,1), (4,3)  ( sp4 = (1-2nu)^2 ) */
    double sp4x8 = SQ((1-2*nu)) * gsl_pow_int(x, 8);
    FNewtlm[5]  = CNlm[5] * sp4x8; /* (4,1) */
    FNewtlm[7]  = CNlm[7] * sp4x8; /* (4,3) */
    
    /* Correcting (5,5) for Higher Modes */
    if (EOBPars->use_flm == USEFLM_HM ){
      FNewtlm[13] = CNlm[13] * sp4x8;
    }
    
  } else {
    if (usetidal) {
      /* Correct (2,1), (3,1) and (3,3) ( sp2 = 1 ) */
      double x6 = gsl_pow_int(x, 6);
      FNewtlm[0] = CNlm[0] * x6; /* (2,1) */
      FNewtlm[2] = CNlm[2] * x6; /* (3,1) */
      FNewtlm[4] = CNlm[4] * x6; /* (3,3) */
    }
  }

  /** Tail term */
  eob_flx_Tlm(E*Omega, MTlm);

  /** Amplitudes */
  if (usespins) {
    /* eob_wav_flm_s_old(x,nu, X1,X2,chi1,chi2,a1,a2,C_Q1,C_Q2, usetidal, rholm, flm); */
    eob_wav_flm_s(x,nu, X1,X2,chi1,chi2,a1,a2,C_Q1,C_Q2, usetidal, rholm, flm);
  } else {
    /* eob_wav_flm_old(x,nu, rholm, flm); */
    eob_wav_flm(x,nu, rholm, flm);
  }
  
  FNewt22 = FNewtlm[1];

  /** NQC correction to the modulus of the (l,m) waveform */
  for (int k = 0; k < KMAX; k++) hlmNQC[k] = 1.; /* no NQC */
  if (!(STREQUAL(nqc_flx_opt[EOBPars->nqc_coefs_flx],"none"))) {
    Waveform_lm_t hNQC;
    /* eob_wav_hlmNQC_nospin201602(nu,r,pr_star,Omega,ddotr, &hNQC); */
    eob_wav_hlmNQC(nu,r,pr_star,Omega,ddotr, NQC->flx, &hNQC);
    const int maxk = MIN(KMAX, NQC->hlm->maxk+1);
    /*
      for (int k = 0; k < maxk; k++) {
      if (NQC->hlm->activemode[k]) {
	hlmNQC[k] = hNQC.ampli[k];
      }
    }
    */
    /* Use only the 22:  */
    hlmNQC[1] = hNQC.ampli[1];
  }
    
  /** Compute modulus of hhat_lm (with NQC) */
  for (int k = 0; k < KMAX; k++) {
    Modhhatlm[k] = prefact[k] * MTlm[k] * flm[k] * hlmNQC[k];
  }
  
  if (usetidal) {
    /** Tidal amplitudes */
    eob_wav_hlmTidal(x,dyn, hlmTidal);
    if (!(usespins)) {
      /* Correct normalization of (2,1) (3,1), (3,3) point-mass amplitudes */
      Modhhatlm[0] *= X12;
      Modhhatlm[2] *= X12;
      Modhhatlm[4] *= X12;
    }
    /* Add tidal amplitudes */
    for (int k = 0; k < KMAX; k++) {
      Modhhatlm[k] += MTlm[k] * hlmTidal[k];
    }
  }

  /** Total multipolar flux */
  for (int k = KMAX; k--;) sum_k += SQ(Modhhatlm[k]) * FNewtlm[k];
  
  /** Normalize to the 22 Newtonian multipole */
  double hatf = sum_k/(FNewt22);
  
  /** Horizon flux */
  if (!(usetidal)) {
    double hatFH;
    if (usespins) {
      hatFH = eob_flx_HorizonFlux_s(x, Heff, jhat, nu, X1, X2, chi1, chi2);
    } else {
      hatFH = eob_flx_HorizonFlux(x,Heff,jhat,nu);
    }
    hatf += hatFH;
  }

  /* return Fphi */
  return (-32./5. * nu * gsl_pow_int(r_omega,4) * gsl_pow_int(Omega,5) * hatf);
}

/** Radial flux calculation for eccentric systems */
double eob_flx_Fr_ecc(double r, double prstar, double pphi, Dynamics *dyn)
{
  const double nu = dyn->nu;
  double nu2 = nu*nu;

  double u  = 1/r;
  double u2 = u*u;
  double u3 = u2*u;
  double u4 = u3*u;
  double u5 = u4*u;
  double pphi2 = pphi*pphi;
  double pphi4 = pphi2*pphi2;
  double pphi6 = pphi4*pphi2;
  double prstar2 = prstar*prstar;
  double prstar4 = prstar2*prstar2;
  
  double a1, a2, F0PN, F1PN, F2PN;
  
  F0PN = -8./15. + 56./5.*pphi2*u;

  F1PN = prstar2*(-1228./105. + 556./105.*nu) + u*(-1984./105. + 16./21.*nu)
    + prstar2*pphi2*u*(-124./105. - 436./105.*nu) + pphi2*u2*(-1696./35. - 1268./105.*nu)
    + pphi4*u3*(1252./105. - 2588./105.*nu);

  F2PN = prstar4*(323./315. + 1061./315.*nu - 1273./315.*nu2)
    + pphi6*u5*(-3229./315. - 718./63.*nu + 3277./105.*nu2)
    + pphi4*u4*(-35209./315. + 1606./15.*nu + 25217./315.*nu2)
    + u2*(59554./2835. + 9686./105.*nu - 3548./315.*nu2)
    + prstar2*pphi2*u2*(-1774./21. + 10292./315.*nu - 8804./315.*nu2)
    + prstar2*u*(20666./315. + 17590./189.*nu - 218./189.*nu2)
    + prstar4*pphi2*u*(-461./315. - 983./315.*nu + 131./63.*nu2)
    + pphi2*u3*(-29438./315. + 58424./315.*nu - 1752./35.*nu2)
    + prstar2*u3*pphi4*(-628./105. - 1052./105.*nu + 194./7.*nu2);

  a1 = F1PN/F0PN;
  a2 = F2PN/F0PN;
  
  /* return Fphi */
  return nu*u4*prstar*F0PN/(1 - a1 + (a1*a1 - a2));
}

/** Non-circular flux for eccentric systems */
double eob_flx_Fphi_ecc(double x, double Omg, double r_omega, double jhat, double rdot, Dynamics *dyn)
{  
  const double nu = dyn->nu;
  const double chi1 = dyn->chi1;
  const double chi2 = dyn->chi2;
  const double aK2  = dyn->aK2;
  const double X1 = dyn->X1;
  const double X2 = dyn->X2;
  const double a1 = dyn->a1;
  const double a2 = dyn->a2;
  const double S1 = dyn->S1;
  const double S2 = dyn->S2;
  const double S  = S1 + S2;
  const double Sstar = X2*a1 + X1*a2;
  const double C_Q1 = dyn->C_Q1;
  const double C_Q2 = dyn->C_Q2;
  const double C_Oct1 = dyn->C_Oct1;
  const double C_Oct2 = dyn->C_Oct2;
  const double C_Hex1 = dyn->C_Hex1;
  const double C_Hex2 = dyn->C_Hex2;
  const double X12   = X1-X2; /* sqrt(1-4nu) */
  const double X12sq = SQ(X12); /* (1-4nu) */
  const double c3    = dyn->cN3LO;

  const int usetidal = dyn->use_tidal;
  const int usespins = dyn->use_spins;

  double r      = dyn->r;
  double pphi   = dyn->pphi;
  double prstar = dyn->prstar;

  double pphi2 = pphi*pphi;
  double prstar2 = prstar*prstar;
  double prstar3 = prstar2*prstar;
  double prstar4 = prstar3*prstar;

  double Frdot = 0.;
  double pphi2dot = 0.;
  
  double A, dA, d2A, B, dB, d2B, Q, sqrtAbyB, dsqrtAbyB_dr, d2sqrtAbyB_d2r, fact;
  double ggm[26], G, dG_dr, dG_dprstar,d2G_dr2, d2G_dr_dprstar, d2G_dprstar2,
    d3G_dr2_dprstar, d3G_dr_dprstar2, d3G_dprstar3;
  double  rc, drc_dr, d2rc_dr2, uc, uc2, uc3, uc4;
  double H, Heff, Heff_orb, E, dHeff_dr, dHeff_dpphi, EHeff_orb;
  double Adot, prstardot, sqrtAbyBdot, dAbyrc2, d2Abyrc2, Omgdot_0, Heffdot, HSOdot, Edot,
    Heff_orbdot, EHeff_orbdot, Omgdot, Omg2dot, r2dot, r3dot, EHeff_orb2dot,
    HSO2dot, Heff_orb2dot, prstar2dot, Heff2dot, E2dot;
  double Fr, Fphi, Fphi_Newt;
  
  /* Computing metric, centrifugal radius and ggm functions*/
  if(usespins) {
    eob_metric_s(r, dyn, &A, &B, &dA, &d2A, &dB, &d2B);
    eob_dyn_s_get_rc(r, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc, &drc_dr, &d2rc_dr2);
    eob_dyn_s_GS(r, rc, drc_dr, d2rc_dr2, aK2, prstar, 0.0, nu, chi1, chi2, X1, X2, c3, ggm);
    G = ggm[2]*S + ggm[3]*Sstar;    // tildeG = GS*S+GSs*Ss
    dG_dr      = ggm[6]*S   + ggm[7]*Sstar;
    dG_dprstar = ggm[4]  *S + ggm[5]  *Sstar;
    d2G_dr2        = ggm[14] *S + ggm[15] *Sstar;
    d2G_dprstar2   = ggm[16] *S + ggm[17] *Sstar;
    d2G_dr_dprstar = ggm[18] *S + ggm[19] *Sstar;
    d3G_dprstar3    = ggm[20] *S + ggm[21] *Sstar;
    d3G_dr2_dprstar = ggm[22] *S + ggm[23] *Sstar;
    d3G_dr_dprstar2 = ggm[24] *S + ggm[25] *Sstar;
  } else {
    eob_metric(r, dyn, &A, &B, &dA, &d2A, &dB, &d2B);
    rc       = r;   //Nonspinning case: rc = r; G = 0;
    drc_dr   = 1.0;
    d2rc_dr2 = 0.0;
    G        = 0.0;
    dG_dr    = 0.0;
    dG_dprstar = 0.0;
    d2G_dr2        = 0.0;
    d2G_dprstar2   = 0.0;
    d2G_dr_dprstar = 0.0;
    d3G_dprstar3    = 0.0;
    d3G_dr2_dprstar = 0.0;
    d3G_dr_dprstar2 = 0.0;
  }
  Q   = 2.*nu*(4. - 3.*nu);
  uc  = 1./rc;
  uc2 = uc*uc;
  uc3 = uc2*uc;
  uc4 = uc3*uc;

  sqrtAbyB = sqrt(A/B);
  fact = (dA*B - A*dB)/SQ(B);
  dsqrtAbyB_dr   = 0.5/sqrtAbyB*fact;
  d2sqrtAbyB_d2r = -0.25*sqrtAbyB*SQ(dA*B - A*dB)/(SQ(A)*SQ(B)) + 0.5/sqrtAbyB*(d2A*SQ(B) - A*B*d2B - 2.*dA*dB*B + 2.*A*SQ(dB))/(SQ(B)*B);

  /** Circular Hamiltonians, ref: arXiv: 1406.6913 */
  if(usespins) {
    eob_ham_s(nu, r, rc, drc_dr, pphi, prstar, S, Sstar, chi1, chi2, X1, X2, aK2, c3, A, dA, &H, &Heff, &Heff_orb, &dHeff_dr, NULL, &dHeff_dpphi, NULL);
    E = nu*H;
  } else {
    eob_ham(nu, r, pphi, prstar, A, dA,
	    &H, &Heff_orb, &dHeff_dr, NULL, &dHeff_dpphi);
    Heff = Heff_orb;
    E = nu*H;
  }

  /* Fluxes */
  Fphi = eob_flx_Flux_s(x, Omg, r_omega, E, Heff, jhat, r, prstar, rdot, dyn);
  Fr   = eob_flx_Fr_ecc(r, prstar, pphi, dyn);

  Edot = nu*(rdot*Fr + Omg*Fphi);
  
  Adot = dA*rdot;
  sqrtAbyBdot = dsqrtAbyB_dr*rdot;

  prstardot = -sqrtAbyB/E*(pphi*dG_dr + 1./(2.*Heff_orb)
			    *(dA*(1. + pphi2*uc2 + Q*prstar4*uc2)
	       - 2.*A*uc3*drc_dr*(pphi2 + Q*prstar4))) + sqrtAbyB*Fr;
    
  dAbyrc2  = (dA*uc2 - 2.*A*uc3*drc_dr);
  d2Abyrc2 = d2A*uc2 - 4.*dA*uc3*drc_dr  + 6.*A*uc4*SQ(drc_dr) - 2.*A*uc3*d2rc_dr2;

  EHeff_orb = E*Heff_orb;
  Omgdot_0 = 1./EHeff_orb*dAbyrc2*rdot*pphi + 1./E*(dG_dr*rdot + dG_dprstar*prstardot);
  
  /* Begin iteration */
  int iter = 2; // Hard-fixed to 2: seems to be enough
    
  for(int n = 1; n <= iter; n++){
    Heffdot  = 1./nu*E*Edot;
    HSOdot   = Fphi*G + pphi*(dG_dr*rdot + dG_dprstar*prstardot);
    Heff_orbdot  = Heffdot - HSOdot;
    EHeff_orbdot = Edot*Heff_orb + E*Heff_orbdot;

    Omgdot = Omgdot_0 + 1./EHeff_orb*A*uc2*(Fphi - pphi/EHeff_orb*EHeff_orbdot)
      - G*Edot/SQ(E);

    r2dot = dsqrtAbyB_dr/sqrtAbyB*SQ(rdot) + sqrtAbyB*
      (1./EHeff_orb*((prstardot - prstar/EHeff_orb*EHeff_orbdot)*(1. + 2.*Q*A*uc2*prstar2)
       + 2.*Q*SQ(prstar)*(dAbyrc2*prstar*rdot + 2.*A*uc2*prstardot))
       + 1./E*(dG_dprstar*(Fphi - pphi*Edot/E)
       + pphi/E*(d2G_dr_dprstar*rdot + d2G_dprstar2*prstardot)));
    
    E2dot    = nu*(r2dot*Fr + Omgdot*Fphi);
    Heff2dot = 1./nu*(SQ(Edot) + E*E2dot);
    
    prstar2dot = dsqrtAbyB_dr/sqrtAbyB*rdot*prstardot + sqrtAbyB
      *(Frdot -0.5/EHeff_orb*(rdot*(d2A + d2Abyrc2*(pphi2 + Q*prstar4))
        - dA*EHeff_orbdot/EHeff_orb + dAbyrc2*(pphi*(2.*Fphi
	- pphi*EHeff_orbdot/EHeff_orb)
	+ Q*prstar3*(4.*prstardot - prstar*EHeff_orbdot/EHeff_orb)))
	- 1./E*(dG_dr*(Fphi - pphi*Edot/E) + pphi*(d2G_dr2*rdot + d2G_dr_dprstar*prstardot)));
    
    HSO2dot = pphi2dot*G + 2.*Fphi*(dG_dr*rdot + dG_dprstar*prstardot)
      + pphi*(d2G_dr2*rdot + dG_dr*r2dot + 2.*d2G_dr_dprstar*rdot*prstardot
	      + d2G_dprstar2*SQ(prstardot) + dG_dprstar*prstar2dot);
    
    Heff_orb2dot  = Heff2dot - HSO2dot;
    EHeff_orb2dot = E2dot*Heff_orb + 2.*Edot*Heff_orbdot + E*Heff_orb2dot;
    
    Omg2dot = 1./EHeff_orb*(d2Abyrc2*pphi*SQ(rdot) + dAbyrc2*pphi*r2dot
	      + 2.*dAbyrc2*rdot*(Fphi - pphi*EHeff_orbdot/EHeff_orb)
	      + A*uc2*(pphi2dot - 2.*Fphi*EHeff_orbdot/EHeff_orb
              + pphi*(2*SQ(EHeff_orbdot/EHeff_orb) - EHeff_orb2dot/EHeff_orb)))
              + 1./E*(d2G_dr2*SQ(rdot) + dG_dr*(r2dot - 2.*rdot*Edot/E)
              + dG_dprstar*(prstar2dot - 2.*prstardot*Edot/E)
	      + d2G_dprstar2*SQ(prstardot) + 2*d2G_dr_dprstar*rdot*prstardot
              + G*(2.*SQ(Edot/E) - E2dot/E));

    r3dot = rdot/sqrtAbyB*(dsqrtAbyB_dr*(3.*r2dot
	    - 2.*dsqrtAbyB_dr/sqrtAbyB*SQ(rdot)) + d2sqrtAbyB_d2r*SQ(rdot))
            + 1./EHeff_orb*sqrtAbyB*(prstar2dot - 2.*prstardot*EHeff_orbdot/EHeff_orb
	    + prstar*(2.*SQ(EHeff_orbdot/EHeff_orb) - EHeff_orb2dot/EHeff_orb)
	    + 2.*Q*prstar*(A*uc2*(6.*SQ(prstardot)
	    + prstar*(3.*prstar2dot - 6.*prstardot*EHeff_orbdot/EHeff_orb
	    + prstar*(2.*SQ(EHeff_orbdot/EHeff_orb) - EHeff_orb2dot/EHeff_orb)))
	    + dAbyrc2*(prstar2*(r2dot - 2.*rdot*EHeff_orbdot/EHeff_orb)
		       + 6.*prstar*prstardot*rdot) + d2Abyrc2*SQ(rdot)*prstar2))
            + sqrtAbyB/E*(dG_dprstar*(pphi2dot - 2.*Fphi*Edot/E
	    - pphi*E2dot/E + 2.*pphi*SQ(Edot/E)) + pphi/E*(d3G_dr2_dprstar*SQ(rdot)
	    + 2*d3G_dr_dprstar2*rdot*prstardot + d3G_dprstar3*SQ(prstardot)
	    + d2G_dr_dprstar*r2dot + d2G_dprstar2*prstar2dot)
	    + 2.*(Fphi - pphi*Edot/E)*(d2G_dr_dprstar*rdot + d2G_dprstar2*prstardot));
    
    Fphi_Newt = Fphi_NewtPref(r, Omg, rdot, r2dot, r3dot, Omgdot, Omg2dot);
    Fphi      = Fphi_Newt*Fphi;
    Edot      = nu*(rdot*Fr + Omg*Fphi);
  }
  
  /* return F_NC */  
  return Fphi_Newt;  
}

/* Generic Newtonian prefactor */
double Fphi_NewtPref(double r, double Omg, double rdot, double r2dot, double r3dot, double Omgdot, double Omg2dot)
{
  double u, u2, u3, u4;
  double Omg2, Omg3, Omg4, Omg5;
  double FphiNewtNC;
  
  u  = 1./r;
  u2 = u*u;
  u3 = u2*u;
  u4 = u3*u;
  Omg2 = Omg*Omg;
  Omg3 = Omg2*Omg;
  Omg4 = Omg3*Omg;
  Omg5 = Omg4*Omg;

  FphiNewtNC  =  1. - r3dot*(Omgdot*u/Omg5/8. + 0.5*rdot*u2/Omg4)
    + r2dot*(0.375*rdot*Omgdot*u2/Omg5 - 2.*u/Omg2 + Omg2dot/8.*u/Omg5)
    + SQ(rdot)*(Omg2dot*u2/Omg5/8. + 4.*u2/Omg2) + SQ(rdot)*rdot*0.75*Omgdot*u3/Omg5 
    + 0.75*SQ(Omgdot)/Omg4 + SQ(r2dot)*0.75*u2/Omg4 + SQ(rdot)*SQ(rdot)*u4/Omg4*0.75
    - Omg2dot/Omg3/4. + rdot*Omgdot*u/Omg3*3.;

  return FphiNewtNC;
}
