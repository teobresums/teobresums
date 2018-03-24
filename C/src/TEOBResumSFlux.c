/**
 * Copyright (C) 2017 Sebastiano Bernuzzi, Gregorio Carullo, Walter Del Pozzo, Alessandro Nagar, Ka Wa Tsang
 * This file is part of TEOBResumS
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
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 */

#include "TEOBResumS.h"

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
  const double X12 = X1-X2;
  //const double sqrt_one_4nu = sqrt(1.-4.*nu);

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
  
  double FNewt22, Modhhatlm;  
  double rholm[KMAX], flm[KMAX], FNewtlm[KMAX], MTlm[KMAX], hlmTidal[KMAX], hlmNQC[KMAX];
  double sum_k=0.; /* sum */

  eob_flx_FlmNewt(x, nu, usetidal, usespins, FNewtlm);
  eob_flx_Tlm(E*Omega, MTlm);
  if (usespins) {
    eob_wav_flm_s(x,nu, X1,X2,chi1,chi2,a1,a2,C_Q1,C_Q2, usetidal, rholm, flm);
  } else {
    eob_wav_flm(x,nu, rholm, flm);
  }
  
  FNewt22 = FNewtlm[1];


  /** Tidal amplitude */
  if (usetidal) {
    eob_wav_hlmTidal(x,dyn, hlmTidal);
  }

  /*
  printf("x %.12e\n",x);
  printf("FNewt22 %.12e\n",FNewt22);
  printf("f22 %.12e\n",flm[1]);
  printf("T22 %.12e\n",MTlm[1]);
  printf("hT22 %.12e\n",hlmTidal[1]);
  DBGSTOP
  */
  
  /** NQC correction to the modulus of the (l,m) waveform */  
  int k;
  Waveform_lm_t NQC;  
  if ( (!(usetidal)) && (!(usespins)) ) {
    eob_wav_hlmNQC(nu,r,pr_star,Omega,ddotr, &NQC);
    for (k = 0; k < KMAX; k++) {
      //hlmNQC[k] = NQC.ampli[k];
      hlmNQC[k] = 1.;
    }
    // Set NQC only in 22:
    k=1;
    hlmNQC[k] = NQC.ampli[k];
  } else {
    for (k = 0; k < KMAX; k++) {
      hlmNQC[k] = 1.;
    }
    //memset(hlmNQC, 1., KMAX*sizeof(double));//FIXME: does not work?!
  }

  /** Sum up */
  /* for (k = 0; k < KMAX; k++) { */
  for (k = KMAX; k--;) { 
    /* Compute modulus of hhat_lm (with NQC) */
    Modhhatlm = prefact[k] * MTlm[k] * flm[k] * hlmNQC[k]; 
    if (usetidal) {
      /* Re-introduce the nu-dependent prefactor in the tidal case */
      if (k==0) Modhhatlm *= X12;
      if (k==2) Modhhatlm *= X12;
      if (k==4) Modhhatlm *= X12;
      Modhhatlm += MTlm[k] * hlmTidal[k];
    }  	
    /* Total flux multipoles */
    //printf("%d %.12e %.12e %.12e %.12e %.12e %.12e\n",k,Modhhatlm, FNewtlm[k],prefact[k],MTlm[k],flm[k],hlmTidal[k]);
    sum_k += SQ(Modhhatlm) * FNewtlm[k];     
  }
  //DBGSTOP

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

/** Coefficients for Newtonian flux */
static const double CNlm[35] = {
  8./45, 32./5,
  1./1260, 32./63, 243./28, 
  1./44100, 32./3969, 729./700, 8192./567, 
  1./19958400, 256./400950, 2187./70400, 131072./66825, 1953125./76032, 
  1./1123782660, 128./28667925, 59049./15415400, 4194304./47779875, 48828125./13621608, 839808./17875,
  1./9.3498717312e11, 32./1.35270135e8, 1594323./3.2064032e10, 4194304./3.07432125e8, 1220703125./5.666588928e9, 5668704./875875, 96889010407./1111968000, 
  1./8.174459284992e13, 32./3.4493884425e10, 177147./3.96428032e10, 4194304./1.5679038375e10, 30517578125./8.00296713216e11, 51018336./1.04229125e8, 4747561509943./4.083146496e11, 274877906944./1688511825.
};

/** Newtonian partial fluxes */
void eob_flx_FlmNewt(double x, double nu, int usetidal, int usespins, double *Nlm)
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
  
  double sp2 = 1.-4.*nu;
  double sp4 = (1-4*nu)*SQ((1-2*nu));
  const double sp3 = (1.-3.*nu)*(1.-3.*nu);
  const double sp5 = (1.-5.*nu+5.*nu2)*(1.-5.*nu+5.*nu2);
  const double sp6 = (1-4*nu)*(3*nu2-4*nu +1)*(3*nu2-4*nu +1);
  const double sp7 = (1 - 7*nu + 14*nu2 - 7*nu3)*(1 - 7*nu + 14*nu2 - 7*nu3);
  const double sp8 = (1 - 4*nu)*(1 - 6*nu + 10*nu2 - 4*nu3)*(1 - 6*nu + 10*nu2 - 4*nu3);
  
  if (usespins) {
    sp2 = 1.;
    sp4 = SQ((2*nu-1));
  } else {
    //TODO: check case with tides and spins 
    if (usetidal) {
      sp2 = 1.;
    } 
  }

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
  int k;
  for (k = 0; k < KMAX; k++) {
    Nlm[k] = CNlm[k] * spx[k];
  }
  
}

/* Factorials evaluated for the tail term */
static const double f14[] = {1.,         1.,          2.,
			     6.,         24.,         120.,
			     720.,       5040.,       40320.,
			     362880.,    3628800.,    39916800.,
			     479001600., 6227020800., 87178291200.};

/** Tail term (modulus) */
void eob_flx_Tlm(const double w, double *MTlm)
{
  double hhatk, x2, y, prod;
  int k, j;    
  for (k = 0; k < KMAX; k++ ) {
    hhatk = MINDEX[k] * w;
    x2    = 4.*hhatk*hhatk;
    prod  = 1.;
    for (j=1; j <= LINDEX[k]; j++ ) {
      prod *= ( j*j + x2 );
    }
    y  = 4.*Pi*hhatk;
    y /= ( 1. - exp(-y) ); 
    MTlm[k] = sqrt( 1./(f14[LINDEX[k]]*f14[LINDEX[k]]) * y * prod );
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
