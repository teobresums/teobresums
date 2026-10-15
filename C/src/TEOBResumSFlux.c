/** \file TEOBResumSFlux.c      
 *  \brief Compute the flux for the TEOBResumS model
 * 
 *  This file contains the routines required for the computation of the radiation reaction
 *  for the TEOBResumS rhs
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

/**
 * Function: eob_flx_FlmNewt
 * -------------------------
 *   Compute the Newtonian partial fluxes
 * 
 *   @param[in] x       :  frequency parameter (v_phi^2)
 *   @param[in] nu     :  symmetric mass ratio
 *   @param[in,out] Nlm:  Newtonian partial fluxes
 */void eob_flx_FlmNewt(double x, double nu, double *Nlm)
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

  /* Newtonian partial fluxes*/
  for (int k = 0; k < KMAX; k++) {
    Nlm[k] = CNlm[k] * spx[k];
  }
  
}

/**
  * Function: eob_flx_Tlm
  * ---------------------
  *   Compute the modulus of the tail term
  * 
  *   @param[in] w    :  E*Omega
  *   @param[out] MTlm:  modulus of the tail term
  */
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

/**
  * Function: eob_flx_HorizonFlux
  * -----------------------------
  *   Compute the horizon-absorbed fluxes
  *   Nagar & Akcay, PRD 85, 044025 (2012)
  *   Bernuzzi, Nagar & Zenginoglu, PRD 86, 104038 (2012)
  * 
  *   @param[in] x   :  frequency parameter, v_phi^2
  *   @param[in] Heff:  effective Hamiltonian
  *   @param[in] jhat:  angular momentum
  *   @param[in] nu  :  symmetric mass ratio
  *   
  *   @return hatFH:  horizon-absorbed fluxes
  */
double eob_flx_HorizonFlux(double x, double Heff, double jhat, double nu)
{
  double rhoHlm[2]; /* only 21,22 multipoles -> k=0,1 */
  double FlmHLO[2];
  double FlmH[2];
  
  /* Shorthands */
  double nu2 = nu*nu;
  double nu3 = nu*nu2;
  double x2  = x*x;
  double x3  = x*x2;
  double x4  = x*x3;
  double x5  = x*x4;
  double x9  = x4*x5;
  double x10 = x*x9;
    
  /* The Newtonian asymptotic contribution */
  const double FNewt22 = 32./5.*x5;
  
  /* Compute leading-order part (nu-dependent) */
  FlmHLO[1] = 32./5.*(1-4*nu+2*nu2)*x9;
  FlmHLO[0] = 32./5.*(1-4*nu+2*nu2)*x10;
    
  /* Compute rho_lm */
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
    
  /* Compute horizon multipolar flux (only l=2) */
  const double Heff2 = Heff*Heff;
  const double jhat2 = jhat*jhat;
  
  FlmH[1] = FlmHLO[1] * Heff2 * gsl_pow_int(rhoHlm[1],4);
  FlmH[0] = FlmHLO[0] * jhat2 * gsl_pow_int(rhoHlm[0],4);
    
  /* Sum over multipoles and normalize to the 22 Newtonian multipole */
  double hatFH = (FlmH[0]+FlmH[1])/FNewt22;
  
  return hatFH;
}

/**
  * Function: eob_flx_HorizonFlux_s
  * -------------------------------
  *   Compute the horizon-absorbed fluxes, spin case
  * 
  *   @param[in] x   :  frequency parameter, v_phi^2
  *   @param[in] Heff:  effective Hamiltonian
  *   @param[in] jhat:  angular momentum
  *   @param[in] nu  :  symmetric mass ratio
  *   @param[in] X1  :  mass fraction of body 1
  *   @param[in] X2  :  mass fraction of body 2
  *   @param[in] chi1:  dimensionless spin 1
  *   @param[in] chi2:  dimensionless spin 2
  * 
  *   @return hatFH:  horizon-absorbed fluxes
  */
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

/**
 * Function: eob_flx_Flux
 * ----------------------
 *   Flux calculation for Newton-Normalized energy flux 
 *   Use the DIN resummation procedure. 
 *   Add non-QC and non-K corrections to (2,2) partial flux.
 * 
 *   @param[in] x      :   frequency parameter, v_phi^2
 *   @param[in] Omega  :  orbital frequency
 *   @param[in] r_omega:  r*psi^(1./3) (from generalized Kepler's law) 
 *   @param[in] E      :  energy
 *   @param[in] Heff   :  effective Hamiltonian
 *   @param[in] jhat   :  angular momentum
 *   @param[in] r      :  radial separation
 *   @param[in] pr_star:  (tortoise) radial momentum
 *   @param[in] ddotr  :  radial acceleration
 *   @param[in] dyn    :  dynamics structure
 * 
 *   @return[out] Fphi:  energy flux
 */
double eob_flx_Flux(double x, double Omega, double r_omega, double E, double Heff, double jhat, double r, double pr_star, double ddotr, Dynamics *dyn)
{
  return eob_flx_Flux_s(x, Omega, r_omega, E, Heff, jhat, r, pr_star, ddotr,dyn);
}

/**
 * Function: eob_flx_Flux_s
 * ------------------------
 *   Flux calculation for spinning systems
 * 
 *   @param[in] x      :  frequency parameter
 *   @param[in] Omega  :  orbital frequency
 *   @param[in] r_omega: r*psi^(1./3) (from generalized Kepler's law) 
 *   @param[in] E      :  energy 
 *   @param[in] Heff   :  effective Hamiltonian
 *   @param[in] jhat   :  angular momentum
 *   @param[in] r      :  radial separation
 *   @param[in] pr_star:  (tortoise) radial momentum
 *   @param[in] ddotr  :  radial acceleration
 *   @param[in] dyn    :  dynamics structure
 * 
 *   @return[out] Fphi:  energy flux
*/
double eob_flx_Flux_s(double x, double Omega, double r_omega, double E, double Heff, double jhat, double r, double pr_star, double ddotr, Dynamics *dyn)
{
  const double nu = EOBPars->nu;
  const double chi1 = EOBPars->chi1;
  const double chi2 = EOBPars->chi2;
  const double X1 = EOBPars->X1;
  const double X2 = EOBPars->X2;
  const double a1 = EOBPars->a1;
  const double a2 = EOBPars->a2;
  const double C_Q1 = EOBPars->C_Q1;
  const double C_Q2 = EOBPars->C_Q2;
  const double X12  = X1-X2; /* sqrt(1-4nu) */
  const double X12sq = SQ(X12); /* (1-4nu) */

  const int usetidal = EOBPars->use_tidal;
  const int usespins = EOBPars->use_spins;
  
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
    if (EOBPars->use_flm == USEFLM_HM || EOBPars->use_flm == USEFLM_HM_4PN22 || EOBPars->use_flm == USEFLM_HM_6PN3p3){
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

  /* Tail term */
  eob_flx_Tlm(E*Omega, MTlm);

  /* Amplitudes */
  if (usespins) {
    /* eob_wav_flm_s_old(x,nu, X1,X2,chi1,chi2,a1,a2,C_Q1,C_Q2, usetidal, rholm, flm); */
    eob_wav_flm_s(x,nu, X1,X2,chi1,chi2,a1,a2,C_Q1,C_Q2, usetidal, rholm, flm);
  } else {
    /* eob_wav_flm_old(x,nu, rholm, flm); */
    eob_wav_flm(x,nu, rholm, flm);
  }
  
  FNewt22 = FNewtlm[1];

  /* NQC correction to the modulus of the (l,m) waveform */  
  for (int k = 0; k < KMAX; k++) hlmNQC[k] = 1.; /* no NQC */
  
  if (!(EOBPars->nqc_coefs_flx == NQC_FLX_NONE)) {
    
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
  
  /* Compute modulus of hhat_lm (with NQC) */  
  for (int k = 0; k < KMAX; k++) { 
    Modhhatlm[k] = prefact[k] * MTlm[k] * flm[k] * hlmNQC[k];
  }

  if (usetidal) {
    /* Tidal amplitudes */
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

  /* Total multipolar flux */
  for (int k = KMAX; k--;) sum_k += SQ(Modhhatlm[k]) * FNewtlm[k];
  
  /* Normalize to the 22 Newtonian multipole */
  double hatf = sum_k/(FNewt22);
    
  /* Horizon flux */ 
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

/**
 * Function: eob_flx_Flux_ecc
 * ------------------------
 *   Flux calculation for eccentric systems
 *   See https://arxiv.org/abs/2001.11736
 * 
 *   @param[in] x        :  frequency parameter
 *   @param[in] Omega    :  orbital frequency
 *   @param[in] r_omega  : r*psi^(1./3) (from generalized Kepler's law) 
 *   @param[in] E        :  energy 
 *   @param[in] Heff     :  effective Hamiltonian
 *   @param[in] jhat     :  angular momentum
 *   @param[in] r        :  radial separation
 *   @param[in] pphi     :  orbital angular momentum
 *   @param[in] pr_star  :  (tortoise) radial momentum
 *   @param[in] ddotr    :  radial acceleration
 *   @param[in,out] Fphi :  angular momentum flux
 *   @param[in,out] Fr   :  radial flux
 *   @param[in] dyn      :  dynamics structure
 * 
*/
void eob_flx_Flux_ecc(double x, double Omega, double r_omega, double E, double Heff, double jhat, double r, double pr_star, double pphi, double rdot, double ddotr, double prsdot, double *Fphi, double *Fr, Dynamics *dyn)
{
  const double nu = EOBPars -> nu;
  const double chi1 = EOBPars -> chi1;
  const double chi2 = EOBPars -> chi2;
  const double X1 = EOBPars -> X1;
  const double X2 = EOBPars -> X2;
  const double a1 = EOBPars -> a1;
  const double a2 = EOBPars -> a2;
  const double C_Q1 = EOBPars -> C_Q1;
  const double C_Q2 = EOBPars -> C_Q2;
  const double X12 = X1 - X2; /* sqrt (1 - 4 nu)*/
  const double X12sq = SQ (X12); /* (1 - 4 nu)*/

  const int usetidal = EOBPars -> use_tidal;
  const int usespins = EOBPars -> use_spins;
  
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
  double Modhhatlm[KMAX], Flm[KMAX];

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
    if (EOBPars->use_flm == USEFLM_HM || EOBPars->use_flm == USEFLM_HM_4PN22 || EOBPars->use_flm == USEFLM_HM_6PN3p3){
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

  /* Amplitudes */
  if (usespins) {
    /* eob_wav_flm_s_old(x,nu, X1,X2,chi1,chi2,a1,a2,C_Q1,C_Q2, usetidal, rholm, flm); */
    eob_wav_flm_s(x,nu, X1,X2,chi1,chi2,a1,a2,C_Q1,C_Q2, usetidal, rholm, flm);
  } else {
    /* eob_wav_flm_old(x,nu, rholm, flm); */
    eob_wav_flm(x,nu, rholm, flm);
  }
  
  FNewt22 = FNewtlm[1];

  /* NQC correction to the modulus of the (l,m) waveform */
  for (int k = 0; k < KMAX; k++) hlmNQC[k] = 1.; /* no NQC */
  
  if (!(EOBPars->nqc_coefs_flx == NQC_FLX_NONE)) {
    
    Waveform_lm_t hNQC;
    /* eob_wav_hlmNQC_nospin201602(nu,r,pr_star,Omega,ddotr, &hNQC); */
    eob_wav_hlmNQC(nu,r,pr_star,Omega,ddotr, NQC->flx, &hNQC);
    const int maxk = MIN(KMAX, NQC->hlm->maxk+1);
    /*
      for (int k = 0; k < maxk; k++) {
      if (NQC->hlm->activemode[k]) {
	hlmNQC[k] = hNQC.ampli[k];
      }
    */
    /* Use only the 22:  */
    hlmNQC[1] = hNQC.ampli[1];
  }
    
  /* Compute modulus of hhat_lm (with NQC) */
  for (int k = 0; k < KMAX; k++) {
    Modhhatlm[k] = prefact[k] * MTlm[k] * flm[k] * hlmNQC[k];
  }
  
  if (usetidal) {
    /* Tidal amplitudes */
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

  /* Total multipolar flux */
  for (int k = 0; k < KMAX; k++) {
    Flm[k] = SQ(Modhhatlm[k]) * FNewtlm[k];
    sum_k += Flm[k];    
  }

  /* Normalize to the 22 Newtonian multipole */
  double hatf = sum_k/(FNewt22);

  /* Horizon flux */
  double hatFH = 0.;
  if (!(usetidal)) {
    if (usespins) {
      hatFH = eob_flx_HorizonFlux_s(x, Heff, jhat, nu, X1, X2, chi1, chi2);
    } else {
      hatFH = eob_flx_HorizonFlux(x,Heff,jhat,nu);
    }
  }

  /* Compute circular, asymptotic and horizon Fphi separately, then sum */
  double Fphi_lo  = -32./5. * nu * gsl_pow_int(r_omega,4) * gsl_pow_int(Omega,5);
  double Fphi_inf = Fphi_lo * hatf;
  double Fphi_H   = Fphi_lo * hatFH;
  *Fphi           = Fphi_inf + Fphi_H;

  /* Compute Fr using the infinity Fphi */
  *Fr = eob_flx_Fr(r, pr_star, pphi, prsdot, dyn, Fphi_inf);
  
  /* Compute non-circular Fphi and Fr */
  double Fphi_NC[KMAX];
  for (int k = 0; k < KMAX; k++) Fphi_NC[k] = 1.;
  eob_flx_Fphi_ecc(r, pr_star, pphi, Omega, rdot, *Fphi, Fphi_lo, FNewtlm, Flm, Fphi_H, Fr, dyn, Fphi_NC);

  /* Adding non-circular corrections and re-compute flux */
  sum_k = 0.;
  for (int k = 0; k < KMAX; k++) sum_k += Flm[k] * Fphi_NC[k];

  hatf = sum_k/(FNewt22);

  /* Compute non-circular Fphi, infinity */
  *Fphi = Fphi_lo * hatf;

  /* Add horizon Fphi */
  *Fphi += Fphi_H;
}

/**
  * Function: eob_flx_Fr
  * --------------------
  *   Radial flux calculation for eccentric systems
  *   See https://arxiv.org/abs/2001.11736
  *       https://arxiv.org/abs/1210.2834
  * 
  *   @param[in] r        :  radial separation
  *   @param[in] pr_star  :  (tortoise) radial momentum
  *   @param[in] pphi     :  orbital angular momentum
  *   @param[in] dyn      :  dynamics structure
  *   @param[in] Fphi     :  angular momentum flux
  * 
  *   @return[out] Fr     :  radial flux
  */
double eob_flx_Fr_ecc(double r, double prstar, double pphi, double prsdot, Dynamics *dyn, double Fphi)
{
  const double nu = EOBPars->nu;
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

  /*
    Old configuration used for arXiv:2001.11736
    F2PN = prstar4*(323./315. + 1061./315.*nu - 1273./315.*nu2)
    + pphi6*u5*(-3229./315. - 718./63.*nu + 3277./105.*nu2)
    + pphi4*u4*(-35209./315. + 1606./15.*nu + 25217./315.*nu2)
    + u2*(59554./2835. + 9686./105.*nu - 3548./315.*nu2)
    + prstar2*pphi2*u2*(-1774./21. + 10292./315.*nu - 8804./315.*nu2)
    + prstar2*u*(20666./315. + 17590./189.*nu - 218./189.*nu2)
    + prstar4*pphi2*u*(-461./315. - 983./315.*nu + 131./63.*nu2)
    + pphi2*u3*(-29438./315. + 58424./315.*nu - 1752./35.*nu2)
    + prstar2*u3*pphi4*(-628./105. - 1052./105.*nu + 194./7.*nu2);
  */
  F2PN = prstar4*(323./315. + 1061./315.*nu - 1273./315.*nu2)
    + pphi6*u5*(-3229./315. - 718./63.*nu + 3277./105.*nu2)
    + pphi4*u4*(-35209./315. + 1606./15.*nu + 25217./315.*nu2)
    + u2*(59554./2835. + 9854./105.*nu - 3548./315.*nu2)
    + prstar2*pphi2*u2*(-1774./21. + 10292./315.*nu - 8804./315.*nu2)
    + prstar2*u*(20666./315. + 17590./189.*nu - 218./189.*nu2)
    + prstar4*pphi2*u*(-461./315. - 983./315.*nu + 131./63.*nu2)
    + pphi2*u3*(-29438./315. + 9568./315.*nu - 1752./35.*nu2)
    + prstar2*u3*pphi4*(-628./105. - 1052./105.*nu + 194./7.*nu2);
  
  a1 = F1PN/F0PN;
  a2 = F2PN/F0PN;
  
  /* return Fr */
  return nu*u4*prstar*F0PN/(1 - a1 + (a1*a1 - a2));
}

/**
  * Function: eob_flx_Fr_ecc_BD
  * ---------------------------
  *   Radial flux calculation for eccentric systems
  *   Circular expression from Bini-Damour inverse-resummed
  *   See https://arxiv.org/abs/1210.2834
  *       https://arxiv.org/abs/2001.11736
  * 
  *   @param[in] r        :  radial separation
  *   @param[in] pr_star  :  (tortoise) radial momentum
  *   @param[in] pphi     :  orbital angular momentum
  *   @param[in] dyn      :  dynamics structure
  *   @param[in] Fphi     :  angular momentum flux
  * 
  *   @return[out] Fr     :  radial flux
*/
double eob_flx_Fr_ecc_BD(double r, double prstar, double pphi, double prsdot, Dynamics *dyn, double Fphi)
{
  const double nu = EOBPars->nu;
  double nu2 = nu*nu;

  double u  = 1/r;
  double u2 = u*u;
  double u4 = u2*u2;

  double c1 = -573./280. - 118./35.*nu;
  double c2 = (-237433. + 175311.*nu + 70794.*nu2)/15120.;

  double if2 = 1./(1. - c1*u + (c1*c1 - c2)*u2);
  
  /* return Fr */
  return nu*32./3.*u4*prstar*if2;
}


/** Radial flux calculation for eccentric systems 
    Eq. 6 of arXiv:2304.09662
*/

/**
  * Function: eob_flx_Fr_ecc_next
  * -----------------------------
  *   Radial flux calculation for eccentric systems
  *   See Eq. 6 of arXiv:2304.09662
  * 
  *   @param[in] r        :  radial separation
  *   @param[in] pr_star  :  (tortoise) radial momentum
  *   @param[in] pphi     :  orbital angular momentum
  *   @param[in] dyn      :  dynamics structure
  *   @param[in] Fphi     :  angular momentum flux
  * 
  *   @return[out] Fr     :  radial flux
*/
double eob_flx_Fr_ecc_next(double r, double prstar, double pphi, double prsdot, Dynamics *dyn, double Fphi)
{
  const double nu = EOBPars->nu;
  const double nu2 = nu*nu;

  const double u  = 1./r;
  const double u2 = u*u;

  //double c1 = 5317./1680 - 227./140*nu;
  //double c2 = 1296935./1016064 - 274793./70560*nu + 753./560*nu2;
  // Correct coefficients
  double c1 = 1957./1680 - 227./140*nu;
  double c2 = -25672541./5080320 + 165703./70560.*nu + 753./560.*nu2;
  
  double a[]  = {1., c1, c2};
  double hatf_prstar = Pade02(u, a);
  const double Frstar = -5./3.*prstar/pphi*Fphi*hatf_prstar;
  
  /* Eq. 6 of the paper above gives Fr*. 
     Compute the conversion from Fr to Fr*
  */
  double A, B, pl_hold;
  eob_metric_s(r, prstar, dyn, &A, &B, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
  const double sqrtAbyB = sqrt(A/B);
  
  /* return Fr */
  return Frstar/sqrtAbyB;
}

/**
  * Function: eob_flx_Fr_ecc_impqc_full
  * -----------------------------
  *   Radial flux calculation for eccentric systems
  *   Fixing a typo in the previous version of the code, plus full expression in BD gauge 
  * 
  *   @param[in] r        :  radial separation
  *   @param[in] pr_star  :  (tortoise) radial momentum
  *   @param[in] pphi     :  orbital angular momentum
  *   @param[in] dyn      :  dynamics structure
  *   @param[in] Fphi     :  angular momentum flux
  * 
  *   @return[out] Fr     :  radial flux
*/
double eob_flx_Fr_ecc_impqc_full(double r, double prstar, double pphi, double prsdot, Dynamics *dyn, double Fphi)
{
  double nu = EOBPars->nu;
  double nu2 = nu*nu;
  double nu3 = nu2*nu;

  double r2 = r * r;
  double r3 = r2 * r;
  double r4 = r3 * r;
  double r5 = r4 * r;
  double r6 = r5 * r;
  double r7 = r6 * r;
  double r8 = r7 * r;
  double r9 = r8 * r;
  double r10 = r9 * r;

  double logr = log(r);
  double sqrt_r = sqrt(r);
  
  double prstar2 = prstar * prstar;
  double prstar4 = prstar2 * prstar2;
  double prstar6 = prstar4 * prstar2;

  double prsdot2 = prsdot * prsdot;
  double prsdot3 = prsdot2 * prsdot;
  double prsdot4 = prsdot3 * prsdot;
  double prsdot5 = prsdot4 * prsdot;
  double prsdot6 = prsdot5 * prsdot;

  /* Newtonian contribution */
  double f0 = 1.0 + 0.25 * prstar2 * r + 0.0625 * prstar4 * r2 
    + prsdot * (0.55 * r2 + 0.0125 * prstar2 * r3 - 0.028125 * prstar4 * r4) 
    + prsdot2 * (-0.275 * r4 - 0.075 * prstar2 * r5) 
    + prsdot3 * (0.1375 * r6 + 0.071875 * prstar2 * r7) 
    - 0.06875 * prsdot4 * r8 + 0.034375 * prsdot5 * r10;

  /* 1PN */
  double f1;
  if (EOBPars->flux_nc_order >= FLUX_NC_ORDER_1PN)
    f1 = (-1.5178571428571428 + 2.288095238095238 * nu) * prstar2 
      + (1.1648809523809525 - 1.6214285714285714 * nu) / r 
      + (-1.3454241071428572 + 0.8674107142857143 * nu) * prstar4 * r 
      + prsdot * ((2.9647321428571427 - 1.3625 * nu) * r 
		  + (0.6417410714285714 + 1.2059523809523809 * nu) * prstar2 * r2 
		  + (0.14171316964285716 - 0.009933035714285714 * nu) * prstar4 * r3) 
      + prsdot2 * ((-1.4866071428571428 + 0.5238095238095238 * nu) * r3 
		   + (-0.86796875 - 0.7486607142857142 * nu) * prstar2 * r4) 
      + prsdot3 * ((1.1324776785714286 - 0.6760416666666667 * nu) * r5 
		   + (0.9021205357142857 + 0.24010416666666667 * nu) * prstar2 * r6) 
      + (-0.7608258928571429 + 0.5450892857142857 * nu) * prsdot4 * r7 
      + (0.4777064732142857 - 0.37607886904761906 * nu) * prsdot5 * r9;
  else
    f1 = 0.0;

  /* 1.5PN */
  double f15;
  if (EOBPars->flux_nc_order >= FLUX_NC_ORDER_15PN)
    f15 = (1.3625 * prstar2 - 1.325 / r + 0.9334375 * prstar4 * r  
	   + prsdot * (-1.59375 * r + 1.3263020833333334 * prstar2 * r2 + 0.46173177083333333 * prstar4 * r3) 
	   + prsdot2 * (0.371875 * r3 - 0.4298177083333333 * prstar2 * r4) 
	   + prsdot3 * (-0.18359375 * r5 - 0.21380208333333334 * prstar2 * r6)
	   + 0.2498046875 * prsdot4 * r7 - 0.35216796875 * prsdot5 * r9) * M_PI / sqrt_r;
  else
    f15 = 0.0;

  /* 2PN */
  double f2;
  if (EOBPars->flux_nc_order >= FLUX_NC_ORDER_2PN)
    f2 = (-6.483223497732427 - 4.3551321570294785 * nu + 0.5247236394557823 * nu2) * prstar4 
      + (-5.053331483056185 + 2.3483985260770974 * nu + 1.3446428571428573 * nu2) / r2 
      + (-16.105473523321365 + 7.908641581632653 * nu - 6.336968537414966 * nu2) * prstar2 / r 
      + prsdot * (7.97368782871945 - 12.197115929705216 * nu + 7.585480442176871 * nu2 
		  + (-6.689099222588341 - 15.166263640873016 * nu - 7.027237457482993 * nu2) * prstar2 * r 
		  + (-0.3731170908033352 - 6.749221185870182 * nu - 0.5030514810090703 * nu2) * prstar4 * r2) 
      + prsdot2 * ((-6.9162761154809145 + 8.336807327097505 * nu + 0.7025687358276644 * nu2) * r2 
		   + (0.055171278580876795 - 0.6395532171201814 * nu - 1.645280612244898 * nu2) * prstar2 * r3) 
      + prsdot3 * ((4.7562702004991815 - 5.866097115929705 * nu - 0.22435870181405895 * nu2) * r4 
		   + (4.31563323643983 + 0.15075600552721088 * nu - 0.7230513038548753 * nu2) * prstar2 * r5)
      + (-3.9656034954491055 + 5.564095273526077 * nu - 0.370999858276644 * nu2) * prsdot4 * r6 
      + (3.0519113738160195 - 4.683648667800454 * nu + 0.7389247626133787 * nu2) * prsdot5 * r8;
  else
    f2 = 0.0;

  /* 2.5PN: logarithmic contribution */
  double f25log = (0.7071428571428572 - 1.6714285714285715 * nu) / r3
    + (-3.367689732142857 - 2.9773995535714284 * nu) * prstar6
    + (-87.37232142857142 + 23.879464285714285 * nu) * prstar2 / r2
    + (-39.06495535714286 + 2.0792410714285716 * nu) * prstar4 / r
    + prsdot * ((24.9375 + 2.6625 * nu) / r + (-66.6 + 46.6875 * nu) * prstar2
		+ (-10.54296875 + 14.14921875 * nu) * prstar4 * r)
    + prsdot2 * ((20.8875 - 12.45 * nu) * r
		 + (11.669866071428572 + 16.637276785714285 * nu) * prstar2 * r2
		 + (7.798549107142857 - 0.19720982142857144 * nu) * prstar4 * r3)
    + prsdot3 * ((-7.680803571428571 - 4.238839285714286 * nu) * r3
		 + (-3.033705357142857 - 1.3051339285714285 * nu) * prstar2 * r4)
    + prsdot4 * ((5.496428571428571 + 1.8763392857142858 * nu) * r5
		 + (3.693861607142857 - 0.3073660714285714 * nu) * prstar2 * r6)
    + (-4.1574776785714285 - 1.7728794642857142 * nu) * prsdot5 * r7
    + (2.0787388392857142 + 0.8864397321428571 * nu) * prsdot6 * r9;

  /* 2.5PN: non-logarithmic contribution */
  double f25pr = (2.149579400518661 + 20.609865314183814 * nu) / r3
    + (-29.28641329644961 - 2168.7665133429837 * nu - 1.2233461364146352 * nu2) * prstar6
    + (34.828639867545974 - 1136.3033436464975 * nu - 0.19577740673815822 * nu2) * prstar4 / r
    + (140.56916399950222 - 105.99383331122496 * nu + 0.5751618008314436 * nu2) * prstar2 / r2
    + prsdot * ((-24.322917997640975 + 105.60631685946906 * nu + 4.848053665910808 * nu2) / r
		+ (127.00745294851276 + 932.4908477285135 * nu - 3.7621974600025196 * nu2) * prstar2
		+ (-5.088853048131568 + 1501.23148964328 * nu - 8.199563499444917 * nu2) * prstar4 * r)
    + prsdot2 * ((-9.366276906784378 - 11.016367724321071 * nu - 12.652521494708994 * nu2) * r
		 + (24.39905112233646 - 899.7446138177474 * nu + 16.790513638904635 * nu2) * prstar2 * r2
		 + (57.257611514331074 - 4764.619801516164 * nu + 21.349330117246552 * nu2) * prstar4 * r3)
    + prsdot3 * ((-12.894068315919867 + 679.302859129737 * nu + 26.923550288170823 * nu2) * r3
		 + (-32.859126558265416 + 4114.924674009206 * nu - 43.08160969683012 * nu2) * prstar2 * r4)
    + prsdot4 * ((-7.4432960753838415 - 228.58301619711165 * nu - 45.58718385357615 * nu2) * r5
		 + (-4.211558913302724 - 4281.606130723667 * nu + 83.45447463816255 * nu2) * prstar2 * r6)
    + (19.926421890064375 + 1573.4648855982523 * nu + 67.87506565777353 * nu2) * prsdot5 * r7
    + (-9.963210945032188 - 786.7324427991225 * nu - 33.937532828886766 * nu2) * prsdot6 * r9;

  /* 2.5PN: pi-dependent contribution */
  double f25pi = (16.62626488095238 - 7.095535714285714 * nu) / r2 + (-5.7873418898809526 - 4.04765625 * nu) * prstar2 / r
    + (-7.86878412155878 + 4.307093331473214 * nu) * prstar4
    + prsdot * (4.939304315476191 - 2.193080357142857 * nu + (1.29674556671627 - 0.4616908482142857 * nu) * prstar2 * r
		+ (-3.3829658726283482 + 6.583151002914186 * nu) * prstar4 * r2)
    + prsdot2 * ((-1.0277994791666667 + 0.2753348214285714 * nu) * r2 + (12.161341979011656 + 3.8923490978422617 * nu) * prstar2 * r3)
    + prsdot3 * ((-11.173858351934523 - 0.20172061011904763 * nu) * r4 + (-11.537350173223587 - 3.5896308051215278 * nu) * prstar2 * r5)
    + (15.658773949032739 + 0.1655180431547619 * nu) * prsdot4 * r6 + (-6.002374965122768 + 0.2285600353422619 * nu) * prsdot5 * r8;

  /* Combine 2.5PN contributions */
  double f25;
  if ((prstar != 0.0) && (EOBPars->flux_nc_order >= FLUX_NC_ORDER_2PN))
    f25 = (f25log * logr + f25pr) / prstar + f25pi * M_PI / sqrt_r;
  else
    f25 = 0.0;
  
  /* 3PN: rational contribution */
  double f3r = (-66.80850356556769 - 76.28681518464018 * nu + 4.444192080241187 * nu2 + 1.9269074675324676 * nu3) / r3
    + (-101.79818896054712 + 449.3471327189581 * nu - 75.15908409517306 * nu2 + 5.782592614065829 * nu3) * prstar2 / r2
    + (42.53014839394018 + 96.74411737041655 * nu - 1.62864164892824 * nu2 - 2.2118504231147726 * nu3) * prstar4 / r
    + prsdot * ((4.488164528620473 - 186.59576170142213 * nu + 39.4806924646888 * nu2 + 0.07876645923520924 * nu3) / r
		+ (-129.38871169101913 + 208.68131897822232 * nu - 36.37685489836718 * nu2 + 20.210584966986435 * nu3) * prstar2
		+ (-48.087255090940744 + 33.80748979530793 * nu - 40.79924193246698 * nu2 - 3.7895778979205317 * nu3) * prstar4 * r)
    + prsdot2 * ((36.63532451381616 + 14.408244301665782 * nu - 10.995149000482227 * nu2 - 9.813822444562241 * nu3) * r
		 + (-8.996568015776575 - 36.04718274806616 * nu - 0.7903821578347248 * nu2 - 7.296047045712224 * nu3) * prstar2 * r2)
    + prsdot3 * ((11.722444008802995 - 36.10913117199574 * nu + 16.73579216645623 * nu2 + 0.288995204412934 * nu3) * r3
		 + (7.0377299000974745 - 4.345356662718373 * nu + 8.2768650094812 * nu2 + 1.8558167218682207 * nu3) * prstar2 * r4)
    + (-17.494416541885585 + 35.17464631974073 * nu - 14.140388573443222 * nu2 + 0.9663049595965731 * nu3) * prsdot4 * r5
    + (16.889053044840693 - 31.376208024853966 * nu + 14.37360057568499 * nu2 - 0.3818776391948825 * nu3) * prsdot5 * r7;

  /* 3PN: logarithmic contribution */
  double f3log = -55.42879556318837 / r3 + (194.61085952715075 + 0.35059523809523807 * nu) * prstar2 / r2
    + (45.50803439031006 - 0.24479166666666666 * nu) * prstar4 / r
    + prsdot * ((-67.60033706363402 + 0.8589285714285714 * nu) / r + (137.15218914998695 - 1.0254464285714286 * nu) * prstar2
		+ (83.17001417058054 + 0.38508184523809524 * nu) * prstar4 * r)
    + prsdot2 * ((-31.97699932276737 - 1.9764880952380952 * nu) * r + (-27.658114099875093 + 0.8965773809523809 * nu) * prstar2 * r2)
    + prsdot3 * ((7.891377012711018 + 3.0302083333333334 * nu) * r3 + (25.529714173288085 - 2.141108630952381 * nu) * prstar2 * r4)
    + (3.05692109442316 - 3.4256696428571427 * nu) * prsdot4 * r5 + (-2.7616699582431465 + 2.4600818452380953 * nu) * prsdot5 * r7;

  /* 3PN: pi-dependent contribution */
  double f3pi = -6.57335179965218 / r3 + (7.827986181029701 - 0.06883473875661375 * nu) * prstar6
    + (-83.38301929297154 + 0.01904761904761905 * nu) * prstar2 / r2 + (8.821572533219296 + 0.2010168650793651 * nu) * prstar4 / r
    + prsdot * ((-0.6501326039101754 + 0.1523809523809524 * nu) / r + (-119.40570425488659 + 0.04851190476190476 * nu) * prstar2
		+ (-12.023064924404025 - 1.0518642526455027 * nu) * prstar4 * r)
    + prsdot2 * ((22.150198283387226 - 0.6154761904761905 * nu) * r + (-7.35819630936021 + 0.5338789682539683 * nu) * prstar2 * r2
		 + (-16.620460820849985 + 2.4524408895502647 * nu) * prstar4 * r3)
    + prsdot3 * ((7.199227046125088 + 1.8268849206349207 * nu) * r3 + (6.685184066998772 - 2.4753513558201057 * nu)* prstar2 * r4)
    + prsdot4 * ((0.022242534105316736 - 3.8092509920634923 * nu) * r5 + (-2.0848379231174476 + 6.253129133597883 * nu) * prstar2 * r6)
    + (-2.581031847774284 + 6.565525793650794 * nu) * prsdot5 * r7 + (1.290515923887142 - 3.282762896825397 * nu) * prsdot6 * r9;

  /* 3PN: pi-logarithmic contribution */
  double f3pilog = 1.35 / r3 + 31.10625 * prstar2 / r2 + 1.37578125 * prstar4 / r + 1.7853515625 * prstar6
    + prsdot * (-0.45 / r + 34.734375 * prstar2 + 6.59484375 * prstar4 * r)
    + prsdot2 * (-7.2 * r + 1.81171875 * prstar2 * r2 + 8.6750390625 * prstar4 * r3)
    + prsdot3 * (-0.815625 * r3 + 3.416015625 * prstar2 * r4)
    + prsdot4 * (-2.687109375 * r5 - 3.82705078125 * prstar2 * r6) + 2.31515625 * prsdot5 * r7
    - 1.157578125 * prsdot6 * r9;

  /* 3PN: pi-squared contribution */
  double f3pi2 = (-3.2 + 5.765625 * nu) / r3 + (-7.0875 - 12.171875 * nu) * prstar2 / r2
    + (-1.6191015625 - 2.1220703125 * nu) * prstar4 / r
    + prsdot * ((0.48125 + 4.1640625 * nu) / r + (-8.4375 - 11.8595703125 * nu) * prstar2
		+ (-0.0853515625 - 2.4984375 * nu) * prstar4 * r)
    + prsdot2 * ((1.2640625 + 2.178125 * nu) * r + (3.1233723958333335 + 3.7716796875 * nu) * prstar2 * r2)
    + prsdot3 * ((-0.3697916666666667 - 1.7056640625 * nu) * r3 + (0.27386067708333334 - 1.115087890625 * nu) * prstar2 * r4)
    + (-1.0970703125 + 1.1611328125 * nu) * prsdot4 * r5 + (0.17977213541666667 - 0.734716796875 * nu) * prsdot5 * r7;

  /* 3PN: squared-logarithmic contribution */
  double f3log2 = 15.1875 / r3 - 32.7375 * prstar2 / r2 - 4.640625 * prstar4 / r
    + prsdot * (8.60625 / r - 17.0015625 * prstar2 - 23.923828125 * prstar4 * r)
    + prsdot2 * (10.040625 * r + 11.93203125 * prstar2 * r2)
    + prsdot3 * (-5.4421875 * r3 - 9.6328125 * prstar2 * r4)
    - 2.63671875 * prsdot4 * r5 + 1.782421875 * prsdot5 * r7;


  /* Combine 3PN contributions */
  double f3;
  if ((prstar != 0.0) && (EOBPars->flux_nc_order >= FLUX_NC_ORDER_3PN))
    f3 = f3r + f3log * logr + (f3pi + f3pilog * logr) * M_PI / sqrt_r / prstar
      + f3pi2 * M_PI * M_PI + f3log2 * logr * logr;
  else
    f3 = 0.0;

  double a2 = f1/f0;
  double a3 = f15/f0;
  double a4 = f2/f0;
  double a5 = f25/f0;
  double a6 = f3/f0;
  
  double hatf_prstar = (1.0 - a2 - a3);
  if (EOBPars->flux_nc_order >= FLUX_NC_ORDER_2PN)
    hatf_prstar += (a2 * a2 - a4);
  if (EOBPars->flux_nc_order >= FLUX_NC_ORDER_3PN)
    hatf_prstar += (- a2*a2*a2 + 2.*a2*a3 + a3*a3 + 2.*a2*a4 - a5 - a6);
  hatf_prstar = f0/hatf_prstar;
  
  const double Frstar = -5./3.*prstar/pphi*Fphi*hatf_prstar;
  
  /* We need to output Fr, not Fr*; compute the conversion from Fr to Fr* */
  double A, B, pl_hold;
  eob_metric_s(r, prstar, dyn, &A, &B, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
  const double sqrtAbyB = sqrt(A/B);
  
  /* return Fr */
  return Frstar/sqrtAbyB;    
}

/** Non-circular flux for eccentric systems */

/**
  * Function: eob_flx_Fphi_ecc
  * --------------------------
  *   Non-circular flux calculation for eccentric systems
  *   obtained via an iterative procedure (two iterations)
  *   See https://arxiv.org/abs/2001.11736
  *   and App. A of https://arxiv.org/pdf/2407.04762
  * 
  *   @param[in] r        :  radial separation
  *   @param[in] pr_star  :  (tortoise) radial momentum
  *   @param[in] pphi     :  orbital angular momentum
  *   @param[in] Omg      :  orbital frequency
  *   @param[in] rdot     :  radial velocity
  *   @param[in] Fphi     :  angular momentum flux
  *   @param[in] Fphi_lo  :  leading order angular momentum flux
  *   @param[in] FlmNewt  :  Newtonian multipolar flux
  *   @param[in] Flm      :  multipolar flux
  *   @param[in] Fphi_H   :  horizon angular momentum flux
  *   @param[in] Fr       :  pointer to radial flux
  *   @param[in] dyn      :  dynamics structure
  *   @param[in] hatflm_NC :  non-circular multipolar flux (empty, to be filled)
  * 
  *   @return[out] Fr     :  radial flux after iterative procedure
  *   @return[out] hatflm_NC :  non-circular corrections to the multipolar flux
*/
void eob_flx_Fphi_ecc(double r, double prstar, double pphi, double Omg, double rdot, double Fphi, double Fphi_lo, double *FlmNewt, double *Flm, double Fphi_H, double *Fr, Dynamics *dyn, double *hatflm_NC)
{  
  const double nu     = EOBPars -> nu;
  const double chi1   = EOBPars -> chi1;
  const double chi2   = EOBPars -> chi2;
  const double aK2    = EOBPars -> aK2;
  const double X1     = EOBPars -> X1;
  const double X2     = EOBPars -> X2;
  const double a1     = EOBPars -> a1;
  const double a2     = EOBPars -> a2;
  const double S1     = EOBPars -> S1;
  const double S2     = EOBPars -> S2;
  const double S      = S1 + S2;
  const double Sstar  = X2*a1 + X1*a2;
  const double C_Q1   = EOBPars -> C_Q1;
  const double C_Q2   = EOBPars -> C_Q2;
  const double C_Oct1 = EOBPars -> C_Oct1;
  const double C_Oct2 = EOBPars -> C_Oct2;
  const double C_Hex1 = EOBPars -> C_Hex1;
  const double C_Hex2 = EOBPars -> C_Hex2;
  const double X12    = X1 - X2; /* sqrt (1 - 4 nu)*/
  const double X12sq  = SQ (X12); /* (1 - 4 nu)*/
  const double c3     = EOBPars -> cN3LO;

  const int usetidal  = EOBPars -> use_tidal;
  const int usespins  = EOBPars -> use_spins;

  double pphi2   = pphi*pphi;
  double prstar2 = prstar*prstar;
  double prstar3 = prstar2*prstar;
  double prstar4 = prstar3*prstar;

  double Frdot = 0.;
  double pphi2dot = 0.;
  
  double A, dA, d2A, B, dB, d2B, Q, dQ, dQ_dprstar, d2Q, d2Q_dprstar2, ddQ_drdprstar, d3Q_dr2dprstar, d3Q_drdprstar2, d3Q_dprstar3;
  double sqrtAbyB, oosqrtAbyB, dsqrtAbyB_dr, d2sqrtAbyB_d2r, fact;
  double sqA,sqB;
  double ggm[26], G, dG_dr, dG_dprstar,d2G_dr2, d2G_dr_dprstar, d2G_dprstar2,
    d3G_dr2_dprstar, d3G_dr_dprstar2, d3G_dprstar3;
  double  rc, drc_dr, d2rc_dr2, uc, uc2, uc3, uc4;
  double H, Heff, Heff_orb, E, dHeff_dr, dHeff_dprstar, dHeff_dpphi, d2Heff_dr2, EHeff_orb;
  double Adot, prstardot, sqrtAbyBdot, dAbyrc2, d2Abyrc2, Omgdot_0, Heffdot, HSOdot, Edot,
    Heff_orbdot, EHeff_orbdot, Omgdot, Omg2dot, r2dot, r3dot, EHeff_orb2dot,
    HSO2dot, Heff_orb2dot, prstar2dot, Heff2dot, E2dot;
  double Fphi_Newt, sum_k;
  double oneby_EHeff_orb, oneby_E;
  double rdot2;
  
  /* Computing metric, centrifugal radius and ggm functions*/
  if(usespins) {
    eob_metric_s(r, prstar, dyn, &A, &B, &dA, &d2A, &dB, &d2B, &Q, &dQ, &dQ_dprstar, &d2Q, &ddQ_drdprstar, &d2Q_dprstar2, &d3Q_dr2dprstar, &d3Q_drdprstar2, &d3Q_dprstar3);
    eob_dyn_s_get_rc(r, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc, &drc_dr, &d2rc_dr2);
    eob_dyn_s_GS(r, rc, drc_dr, d2rc_dr2, aK2, prstar, 0.0, nu, chi1, chi2, X1, X2, c3, ggm);
    G = ggm[2]*S + ggm[3]*Sstar;    // tildeG = GS*S+GSs*Ss
    dG_dr           = ggm[6]*S   + ggm[7]*Sstar;
    dG_dprstar      = ggm[4]  *S + ggm[5]  *Sstar;
    d2G_dr2         = ggm[14] *S + ggm[15] *Sstar;
    d2G_dprstar2    = ggm[16] *S + ggm[17] *Sstar;
    d2G_dr_dprstar  = ggm[18] *S + ggm[19] *Sstar;
    d3G_dprstar3    = ggm[20] *S + ggm[21] *Sstar;
    d3G_dr2_dprstar = ggm[22] *S + ggm[23] *Sstar;
    d3G_dr_dprstar2 = ggm[24] *S + ggm[25] *Sstar;
  } else {
    eob_metric(r, prstar, dyn, &A, &B, &dA, &d2A, &dB, &d2B, &Q, &dQ, &dQ_dprstar, &d2Q, &ddQ_drdprstar, &d2Q_dprstar2, &d3Q_dr2dprstar, &d3Q_drdprstar2, &d3Q_dprstar3);
    rc              = r;   //Nonspinning case: rc = r; G = 0;
    drc_dr          = 1.0;
    d2rc_dr2        = 0.0;
    G               = 0.0;
    dG_dr           = 0.0;
    dG_dprstar      = 0.0;
    d2G_dr2         = 0.0;
    d2G_dprstar2    = 0.0;
    d2G_dr_dprstar  = 0.0;
    d3G_dprstar3    = 0.0;
    d3G_dr2_dprstar = 0.0;
    d3G_dr_dprstar2 = 0.0;
  }
  // Q   = 2.*nu*(4. - 3.*nu); // this is z3 !!
  uc  = 1./rc;
  uc2 = uc*uc;
  uc3 = uc2*uc;
  uc4 = uc3*uc;

  sqrtAbyB = sqrt(A/B);
  oosqrtAbyB = 1./sqrtAbyB;
  sqB      = SQ(B);
  sqA      = SQ(A);
  fact = (dA*B - A*dB)/sqB;
  dsqrtAbyB_dr   = 0.5*oosqrtAbyB*fact;
  /*d2sqrtAbyB_d2r = -0.25*sqrtAbyB*SQ(dA*B - A*dB)/(SQ(A)*SQ(B)) + 0.5/sqrtAbyB*(d2A*SQ(B) - A*B*d2B - 2.*dA*dB*B + 2.*A*SQ(dB))/(SQ(B)*B);*/
  d2sqrtAbyB_d2r = -0.25*sqrtAbyB*SQ(dA*B - A*dB)/(sqA*sqB) + 0.5*oosqrtAbyB*(d2A*sqB - A*B*d2B - 2.*dA*dB*B + 2.*A*SQ(dB))/(sqB*B);

  /* Circular Hamiltonians, ref: arXiv: 1406.6913 */
  if(usespins) {
    eob_ham_s(nu, r, rc, drc_dr, d2rc_dr2, pphi, prstar, S, Sstar, chi1, chi2, X1, X2, aK2, c3, A, dA, d2A, Q, dQ, dQ_dprstar, d2Q, d2Q_dprstar2, &H, &Heff, &Heff_orb, &dHeff_dr, NULL, &dHeff_dpphi, NULL, NULL);
    E = nu*H;
  } else {
    eob_ham(nu, r, pphi, prstar, A, dA, Q, dQ, dQ_dprstar,
	    &H, &Heff_orb, &dHeff_dr, NULL, &dHeff_dpphi);
    Heff = Heff_orb;
    E = nu*H;
  }
  
  Adot = dA*rdot;
  sqrtAbyBdot = dsqrtAbyB_dr*rdot;
    
  dAbyrc2  = (dA*uc2 - 2.*A*uc3*drc_dr);
  d2Abyrc2 = d2A*uc2 - 4.*dA*uc3*drc_dr  + 6.*A*uc4*SQ(drc_dr) - 2.*A*uc3*d2rc_dr2;

  /* shorthands */
  EHeff_orb       = E*Heff_orb;
  oneby_EHeff_orb = 1.0/EHeff_orb;
  oneby_E         = 1.0/E;
  rdot2           = SQ(rdot);
  
  /* Begin iteration */
  int iter = 2; // Hard-fixed to 2: seems to be enough
    
  for(int n = 1; n <= iter; n++){
    
    prstardot = - sqrtAbyB/E*(pphi*dG_dr + 1./(2.*Heff_orb)
			        *(dA*(1. + pphi2*uc2 + Q)
	            + A*(-2.*uc3*pphi2*drc_dr + dQ))) + sqrtAbyB* (*Fr);

    Edot = nu*(rdot* (*Fr) + Omg*Fphi);
    Heffdot  = 1./nu*E*Edot;
    HSOdot   = Fphi*G + pphi*(dG_dr*rdot + dG_dprstar*prstardot);
    Heff_orbdot  = Heffdot - HSOdot;
    EHeff_orbdot = Edot*Heff_orb + E*Heff_orbdot;

    Omgdot = oneby_EHeff_orb*dAbyrc2*rdot*pphi + oneby_E*(dG_dr*rdot + dG_dprstar*prstardot) 
           + oneby_EHeff_orb*A*uc2*(Fphi - pphi*oneby_EHeff_orb*EHeff_orbdot)
           - G*Edot*SQ(oneby_E);

    r2dot = dsqrtAbyB_dr*oosqrtAbyB*rdot2 
          + sqrtAbyB*(oneby_EHeff_orb*(-oneby_EHeff_orb*EHeff_orbdot)*(prstar + 0.5*A*dQ_dprstar)
          + 0.5*oneby_EHeff_orb*(rdot*(dA*dQ_dprstar + A*ddQ_drdprstar) + prstardot*(2. + A*d2Q_dprstar2))
          + oneby_E*(pphi*(rdot*d2G_dr_dprstar + prstardot*d2G_dprstar2) + dG_dprstar*(Fphi - pphi*Edot*oneby_E)));
    
    
    E2dot    = nu*(r2dot* (*Fr) + Omgdot*Fphi);
    Heff2dot = 1./nu*(SQ(Edot) + E*E2dot);
    
    // derivatives of A*(1 + pphi2*u2 + Q) wrt dr2, drdprstar, drdpphi
    double der1 = d2A*(1. + pphi2*uc2 + Q) + 2.*dA*(-2.*uc3*pphi2*drc_dr + dQ) + A*(d2Q + 6.*pphi2*uc4*SQ(drc_dr) - 2.*pphi2*uc3*d2rc_dr2);
    double der2 = dA*dQ_dprstar + A*ddQ_drdprstar;
    double der3 = 2.*dA*pphi*uc2 - 4.*A*pphi*uc3*drc_dr;

    prstar2dot = dsqrtAbyB_dr*oosqrtAbyB*rdot*prstardot 
             + 0.5*sqrtAbyB*oneby_EHeff_orb*(oneby_EHeff_orb*EHeff_orbdot*(dA*(1. + pphi2*uc2 + Q) + A*(dQ - 2.*uc3*pphi2*drc_dr))
             - (der1*rdot + der2*prstardot + der3*Fphi)) 
             - sqrtAbyB*oneby_E*(pphi*(rdot*d2G_dr2 + prstardot*d2G_dr_dprstar) + dG_dr*(Fphi - pphi*Edot*oneby_E))
             + sqrtAbyB*Frdot;
    
    HSO2dot = pphi2dot*G + 2.*Fphi*(dG_dr*rdot + dG_dprstar*prstardot)
            + pphi*(d2G_dr2*rdot2 + dG_dr*r2dot + 2.*d2G_dr_dprstar*rdot*prstardot
	          + d2G_dprstar2*SQ(prstardot) + dG_dprstar*prstar2dot);
    
    Heff_orb2dot  = Heff2dot - HSO2dot;
    EHeff_orb2dot = E2dot*Heff_orb + 2.*Edot*Heff_orbdot + E*Heff_orb2dot;
    
    Omg2dot = oneby_EHeff_orb*(d2Abyrc2*pphi*SQ(rdot) + dAbyrc2*pphi*r2dot
	      + 2.*dAbyrc2*rdot*(Fphi - pphi*EHeff_orbdot*oneby_EHeff_orb)
	      + A*uc2*(pphi2dot - 2.*Fphi*EHeff_orbdot*oneby_EHeff_orb
              + pphi*(2.0*SQ(EHeff_orbdot*oneby_EHeff_orb) - EHeff_orb2dot*oneby_EHeff_orb)))
              + oneby_E*(d2G_dr2*SQ(rdot) + dG_dr*(r2dot - 2.*rdot*Edot*oneby_E)
              + dG_dprstar*(prstar2dot - 2.*prstardot*Edot*oneby_E)
	      + d2G_dprstar2*SQ(prstardot) + 2*d2G_dr_dprstar*rdot*prstardot
              + G*(2.*SQ(Edot*oneby_E) - E2dot*oneby_E));

    
    double comb1 = sqrtAbyB*oneby_EHeff_orb*(-EHeff_orbdot*oneby_EHeff_orb)*(prstar + 0.5*A*dQ_dprstar); 
    double comb2 = sqrtAbyB*0.5*oneby_EHeff_orb*(rdot*(dA*dQ_dprstar + A*ddQ_drdprstar) + prstardot*(2. + A*d2Q_dprstar2)); 
    double comb3 = sqrtAbyB*oneby_E*(pphi*(rdot*d2G_dr_dprstar + prstardot*d2G_dprstar2) + dG_dprstar*(Fphi - pphi*Edot*oneby_E)); 
    double comb4 = -EHeff_orbdot*oneby_EHeff_orb;
    double comb5 = prstar + 0.5*A*dQ_dprstar;

    double D1 = 2.*rdot*r2dot*oosqrtAbyB*dsqrtAbyB_dr + rdot2*rdot*(-SQ(dsqrtAbyB_dr)*(B/A) + oosqrtAbyB*d2sqrtAbyB_d2r);
    double D2 = comb1*rdot*oosqrtAbyB*dsqrtAbyB_dr + sqrtAbyB*(comb5*oneby_EHeff_orb*(2.*SQ(comb4) - EHeff_orb2dot*oneby_EHeff_orb) 
                + oneby_EHeff_orb*comb4*(prstardot*(1. + 0.5*A*d2Q_dprstar2) + 0.5*rdot*(dA*dQ_dprstar + A*ddQ_drdprstar)));

    double D3 = rdot*dsqrtAbyB_dr*comb2*oosqrtAbyB - EHeff_orbdot*oneby_EHeff_orb*comb2 
                + 0.5*sqrtAbyB*oneby_EHeff_orb*(r2dot*(dA*dQ_dprstar + A*ddQ_drdprstar) 
                + rdot2*(d2A*dQ_dprstar + 2.*dA*ddQ_drdprstar + A*d3Q_dr2dprstar)
                + 2.*rdot*prstardot*(d2Q_dprstar2*dA + A*d3Q_drdprstar2)
                + SQ(prstardot)*A*d3Q_dprstar3 + prstar2dot*(2. + A*d2Q_dprstar2));
    
    double D4 = rdot*dsqrtAbyB_dr*comb3*oosqrtAbyB - Edot*comb3*oneby_E 
                + sqrtAbyB*oneby_E*(d2G_dr_dprstar*(2.*Fphi*rdot + pphi*r2dot - rdot*pphi*Edot*oneby_E)
                                  + d2G_dprstar2*(2.*Fphi*prstardot + pphi*prstar2dot - prstardot*pphi*Edot*oneby_E)
                                  + pphi*(2.*rdot*prstardot*d3G_dr_dprstar2 + rdot2*d3G_dr2_dprstar + SQ(prstardot)*d3G_dprstar3)
                                  + dG_dprstar*(pphi2dot - (Fphi*Edot + pphi*E2dot)*oneby_E + pphi*SQ(Edot*oneby_E)));

    r3dot = D1 + D2 + D3 + D4; 
    
    // NC corrections to the modes, Newtonian (FlmNewt_nc) times PN corrections (hatflm_NC)
    for (int k=0; k < KMAX; k++){
      hatflm_NC[k] = eob_flx_FlmNewt_nc[k](r, Omg, rdot, r2dot, r3dot, Omgdot, Omg2dot);
      hatflm_NC[k] = hatflm_NC[k] * eob_flx_hatflm_nc[k](r, prstar, prstardot);
    }

    sum_k = 0.;
    for (int k = 0; k < KMAX; k++) sum_k += Flm[k] * hatflm_NC[k];
    sum_k = sum_k/FlmNewt[1];

    Fphi  = Fphi_lo * sum_k;
    *Fr    = eob_flx_Fr(r, prstar, pphi, prstardot, dyn, Fphi);
    Fphi  = Fphi + Fphi_H;
  } // end iteration
  
  /* Saving useful variables */
  dyn->rdot  = rdot;
  dyn->r2dot = r2dot;
  dyn->r3dot = r3dot;
  dyn->r4dot = 0.;
  dyn->r5dot = 0.;
  dyn->Omegadot  = Omgdot;
  dyn->Omega2dot = Omg2dot;
  dyn->Omega3dot = 0.;
  dyn->Omega4dot = 0.;
  dyn->prsdot = prstardot;
}

/* Generic Newtonian prefactors */
double eob_flx_FlmNewt_nc_22_old(double r, double Omg, double rdot, double r2dot, double r3dot, double Omgdot, double Omg2dot)
{
  double u, u2, u3, u4;
  /*double Omg2, Omg3, Omg4, Omg5;*/
  double invOmg,invOmg2,invOmg3,invOmg4,invOmg5;
  double rdot2, rdot3, rdot4;
  double FphiNewtNC;
  
  u  = 1./r;
  u2 = u*u;
  u3 = u2*u;
  u4 = u3*u;
  /*Omg2 = Omg*Omg;
  Omg3 = Omg2*Omg;
  Omg4 = Omg3*Omg;
  Omg5 = Omg4*Omg;

  FphiNewtNC  =  1. - r3dot*(Omgdot*u/Omg5/8. + 0.5*rdot*u2/Omg4)
    + r2dot*(0.375*rdot*Omgdot*u2/Omg5 - 2.*u/Omg2 + Omg2dot/8.*u/Omg5)
    + SQ(rdot)*(Omg2dot*u2/Omg5/8. + 4.*u2/Omg2) + SQ(rdot)*rdot*0.75*Omgdot*u3/Omg5 
    + 0.75*SQ(Omgdot)/Omg4 + SQ(r2dot)*0.75*u2/Omg4 + SQ(rdot)*SQ(rdot)*u4/Omg4*0.75
    - Omg2dot/Omg3/4. + rdot*Omgdot*u/Omg3*3.;*/

  invOmg  = 1./Omg;
  invOmg2 = invOmg*invOmg;
  invOmg3 = invOmg*invOmg2;
  invOmg4 = invOmg*invOmg3;
  invOmg5 = invOmg*invOmg4;

  rdot2 = SQ(rdot);
  rdot3 = rdot2*rdot;
  rdot4 = rdot2*rdot2;

  FphiNewtNC  =  1. - r3dot*(0.125*Omgdot*u*invOmg5 + 0.5*rdot*u2*invOmg4)
              + r2dot*(0.375*rdot*Omgdot*u2*invOmg5 - 2.0*u*invOmg2             + 0.125*Omg2dot*u*invOmg5)
              + rdot2*(0.125*Omg2dot*u2*invOmg5     + 4.*u2*invOmg2)            + rdot3*0.75*Omgdot*u3*invOmg5 
              + 0.75*SQ(Omgdot)*invOmg4             + SQ(r2dot)*0.75*u2*invOmg4 + rdot4*u4*invOmg4*0.75
              - 0.25*Omg2dot*invOmg3                + 3.0*rdot*Omgdot*u*invOmg3;
  
  return FphiNewtNC;
}

/*========================================================================== */
/*(2,2) Mode (this should in principle be identical to the previous one... is it?)                                                            */
/*========================================================================== */
double eob_flx_FlmNewt_nc_22(double r, double Omg, double rdot, double r2dot, double r3dot, double Omgdot, double Omg2dot)
{
    const double invOmg  = 1.0 / Omg;
    const double invOmg2 = invOmg  * invOmg;
    const double invOmg3 = invOmg2 * invOmg;
    const double invOmg4 = invOmg3 * invOmg;
    const double invOmg5 = invOmg4 * invOmg;

    const double u   = 1.0 / r;
    const double u2  = u   * u;
    const double u3  = u2  * u;
    const double u4  = u3  * u;

    const double rdot2 = rdot  * rdot;
    const double rdot3 = rdot2 * rdot;
    const double rdot4 = rdot3 * rdot;

    const double r2dot2 = r2dot * r2dot;
    const double Omgdot2 = Omgdot * Omgdot;

    double FphiNewtNC = 1.0 + invOmg2 * (-2.0 * u * r2dot + 4.0 * u2 * rdot2) 
        + invOmg3 * (-0.25 * Omg2dot + 3.0 * u * Omgdot * rdot) 
        + invOmg4 * (0.75 * Omgdot2 + 0.75 * u4 * rdot4 + u2 * (0.75 * r2dot2 - 0.5 * r3dot * rdot)) 
        + invOmg5 * (u * (0.125 * Omg2dot * r2dot - 0.125 * Omgdot * r3dot) 
        + 0.75 * u3 * Omgdot * rdot3 + u2 * (0.375 * Omgdot * r2dot * rdot + 0.125 * Omg2dot * rdot2));

    return FphiNewtNC;
}

/*========================================================================== */
/* (2,1) Mode                                                                 */
/*========================================================================== */
double eob_flx_FlmNewt_nc_21(double r, double Omg, double rdot, double r2dot, double r3dot, double Omgdot, double Omg2dot)
{
    const double invOmg  = 1.0 / Omg;
    const double invOmg2 = invOmg  * invOmg;
    const double invOmg3 = invOmg2 * invOmg;
    const double invOmg4 = invOmg3 * invOmg;
    const double invOmg5 = invOmg4 * invOmg;
    const double invOmg6 = invOmg5 * invOmg;
    const double invOmg7 = invOmg6 * invOmg;

    const double u   = 1.0 / r;
    const double u2  = u   * u;
    const double u3  = u2  * u;
    const double u4  = u3  * u;

    const double rdot2 = rdot  * rdot;
    const double rdot3 = rdot2 * rdot;
    const double rdot4 = rdot3 * rdot;

    const double r2dot2 = r2dot * r2dot;
    const double Omgdot2 = Omgdot * Omgdot;
    const double Omg2dot2 = Omg2dot * Omg2dot;
    const double Omgdot3 = Omgdot2 * Omgdot;

    double FphiNewtNC = 1.0 + invOmg2 * u * (-12.0 * r2dot + 30.0 * u * rdot2) 
        + invOmg3 * (-5.0 * Omg2dot + 30.0 * u * Omgdot * rdot) 
        + invOmg4 * (15.0 * Omgdot2 - 36.0 * u3 * r2dot * rdot2 
        + 96.0 * u4 * rdot4 + u2 * (27.0 * r2dot2 - 6.0 * r3dot * rdot)) 
        + invOmg5 * u * (21.0 * Omg2dot * r2dot - 3.0 * Omgdot * r3dot 
        + 33.0 * u * Omgdot * r2dot * rdot + 132.0 * u2 * Omgdot * rdot3) 
        + invOmg6 * (4.0 * Omg2dot2 + 114.0 * u2 * Omgdot2 * rdot2 
        + u * (-6.0 * Omgdot2 * r2dot + 30.0 * Omg2dot * Omgdot * rdot)) 
        + invOmg7 * (3.0 * Omg2dot * Omgdot2 + 18.0 * u * Omgdot3 * rdot);

    return FphiNewtNC;
}

/*========================================================================== */
/*(3,3) Mode                                                                 */
/*========================================================================== */
double eob_flx_FlmNewt_nc_33(double r, double Omg, double rdot, double r2dot, double r3dot, double Omgdot, double Omg2dot)
{
    const double invOmg  = 1.0 / Omg;
    const double invOmg2 = invOmg  * invOmg;
    const double invOmg3 = invOmg2 * invOmg;
    const double invOmg4 = invOmg3 * invOmg;
    const double invOmg5 = invOmg4 * invOmg;
    const double invOmg6 = invOmg5 * invOmg;
    const double invOmg7 = invOmg6 * invOmg;

    const double u   = 1.0 / r;
    const double u2  = u   * u;
    const double u3  = u2  * u;
    const double u4  = u3  * u;
    const double u5  = u4  * u;
    const double u6  = u5  * u;
    const double u7  = u6  * u;

    const double rdot2 = rdot  * rdot;
    const double rdot3 = rdot2 * rdot;
    const double rdot4 = rdot3 * rdot;
    const double rdot6 = rdot3 * rdot3;

    const double r2dot2 = r2dot  * r2dot;
    const double r2dot3 = r2dot2 * r2dot;
    const double r3dot2 = r3dot  * r3dot;
    const double rdot5   = rdot4   * rdot;
    
    const double Omgdot2 = Omgdot  * Omgdot;
    const double Omg2dot2 = Omg2dot * Omg2dot;
    const double Omgdot3 = Omgdot2 * Omgdot;

    double FphiNewtNC = 1.0 + 0.037037037037037035 * ( 
        invOmg2 * (-80.0 * u * r2dot + 162.0 * u2 * rdot2) 
        + invOmg3 * (-13.0 * Omg2dot + 139.0 * u * Omgdot * rdot) 
        + invOmg4 * (47.0 * Omgdot2 - 6.0 * u5 * rdot2 - 80.0 * u3 * r2dot * rdot2 
        + u2 * (59.0 * r2dot2 - 34.0 * r3dot * rdot) + u4 * (2.0 * r2dot + 132.0 * rdot4)) 
        + invOmg5 * u * (15.88888888888889 * Omg2dot * r2dot - 16.0 * Omgdot * r3dot 
        + 164.0 * u2 * Omgdot * rdot3 + u * rdot * (Omgdot * r2dot + 2.0 * Omg2dot * rdot)) 
        + invOmg6 * (1.1111111111111112 * Omg2dot2 + 12.0 * u7 * rdot4 
        + u * (-7.0 * Omgdot2 * r2dot + 11.555555555555555 * Omg2dot * Omgdot * rdot) 
        + u3 * (-6.0 * r2dot3 + 7.333333333333333 * r2dot * r3dot * rdot) 
        + u2 * (1.1111111111111112 * r3dot2 + 78.0 * Omgdot2 * rdot2) 
        + u4 * (28.0 * r2dot2 * rdot2 - 11.555555555555555 * r3dot * rdot3) 
        + u5 * (-2.0 * r2dot2 + 6.666666666666667 * r2dot * rdot4) 
        + u6 * (2.0 * r2dot * rdot2 + 4.444444444444445 * rdot6)) 
        + invOmg7 * (0.7777777777777778 * Omg2dot * Omgdot2 + 7.0 * u * Omgdot3 * rdot 
        + u2 * (-0.6666666666666666 * Omg2dot * r2dot2 
        + 1.5555555555555556 * Omgdot * r2dot * r3dot + 0.2222222222222222 * Omg2dot * r3dot * rdot) 
        + u3 * (3.3333333333333335 * Omgdot * r2dot2 * rdot 
        + (5.333333333333333 * Omg2dot * r2dot - 3.5555555555555554 * Omgdot * r3dot) * rdot2) 
        + u4 * (-0.2222222222222222 * Omg2dot * r2dot + 17.77777777777778 * Omgdot * r2dot * rdot3 
        + 2.2222222222222223 * Omg2dot * rdot4) 
        + u5 * (-2.0 * Omgdot * r2dot * rdot + 0.6666666666666666 * Omg2dot * rdot2 
        + 8.88888888888889 * Omgdot * rdot5) + 6.0 * u6 * Omgdot * rdot3));

    return FphiNewtNC;
}

/*========================================================================== */
/*(3,2) Mode                                                                 */
/*========================================================================== */
double eob_flx_FlmNewt_nc_32(double r, double Omg, double rdot, double r2dot, double r3dot, double Omgdot, double Omg2dot)
{
    const double invOmg  = 1.0 / Omg;
    const double invOmg2 = invOmg  * invOmg;
    const double invOmg3 = invOmg2 * invOmg;
    const double invOmg4 = invOmg3 * invOmg;
    const double invOmg5 = invOmg4 * invOmg;
    const double invOmg6 = invOmg5 * invOmg;
    const double invOmg7 = invOmg6 * invOmg;
    const double invOmg8 = invOmg7 * invOmg;
    const double invOmg9 = invOmg8 * invOmg;

    const double u   = 1.0 / r;
    const double u2  = u   * u;
    const double u3  = u2  * u;
    const double u4  = u3  * u;
    const double u5  = u4  * u;
    const double u6  = u5  * u;
    const double u7  = u6  * u;
    const double u8 = u4 * u4;
    const double u10 = u4  * u6; 
    const double u9  = u3  * u6;

    const double rdot2 = rdot  * rdot;
    const double rdot3 = rdot2 * rdot;
    const double rdot4 = rdot3 * rdot;
    const double rdot5 = rdot4 * rdot;
    const double rdot6 = rdot3 * rdot3;

    const double r2dot2 = r2dot  * r2dot;
    const double r3dot2 = r3dot  * r3dot;

    const double Omgdot2 = Omgdot  * Omgdot;
    const double Omgdot3 = Omgdot2 * Omgdot;
    const double Omgdot4 = Omgdot2 * Omgdot2;

    const double Omg2dot2 = Omg2dot * Omg2dot;

    double FphiNewtNC = 1.0 + 0.0625 * (
          invOmg2 * (-142.0 * u * r2dot + 336.0 * u2 * rdot2) 
        + invOmg3 * (-52.0 * Omg2dot + 344.0 * u * Omgdot * rdot) 
        + invOmg4 * (172.0 * Omgdot2 - 12.0 * u5 * rdot2 - 846.0 * u3 * r2dot * rdot2 
        + u2 * (324.0 * r2dot2 - 76.0 * r3dot * rdot) + u4 * (4.0 * r2dot + 1340.0 * rdot4)) 
        + invOmg5 * u * (220.0 * Omg2dot * r2dot - 46.0 * Omgdot * r3dot 
        + 1988.0 * u2 * Omgdot * rdot3 + u * rdot * (-136.0 * Omgdot * r2dot - 146.0 * Omg2dot * rdot)) 
        + invOmg6 * (36.0 * Omg2dot2 + 108.0 * u7 * rdot4 + 560.0 * u6 * rdot6 
        + u * Omgdot * (-153.5 * Omgdot * r2dot + 156.0 * Omg2dot * rdot) 
        + u3 * r2dot * (-126.0 * r2dot2 + 66.0 * r3dot * rdot) + u2 * (3.0 * r3dot2 + 1712.0 * Omgdot2 * rdot2) 
        + u4 * rdot2 * (837.0 * r2dot2 - 184.0 * r3dot * rdot) 
        + u5 * r2dot * (-12.0 * r2dot - 108.0 * rdot4)) 
        + invOmg7 * (23.0 * Omg2dot * Omgdot2 + 498.0 * u * Omgdot3 * rdot 
        + 108.0 * u6 * Omgdot * rdot3 + u2 * (-96.0 * Omg2dot * r2dot2 + 8.0 * Omgdot * r2dot * r3dot 
        + 14.0 * Omg2dot * r3dot * rdot) 
        + u3 * (198.0 * Omgdot * r2dot2 * rdot + (486.0 * Omg2dot * r2dot - 164.0 * Omgdot * r3dot) * rdot2) 
        + u4 * (-4.0 * Omg2dot * r2dot + 1420.0 * Omgdot * r2dot * rdot3 + 30.0 * Omg2dot * rdot4) 
        + u5 * rdot * (-36.0 * Omgdot * r2dot + 12.0 * Omg2dot * rdot + 860.0 * Omgdot * rdot4)) 
        + invOmg8 * (42.0 * Omgdot4 + u * (-18.0 * Omg2dot2 * r2dot - 5.0 * Omg2dot * Omgdot * r3dot) 
        + 9.0 * u5 * Omgdot2 * rdot2 + u2 * (72.5 * Omgdot2 * r2dot2 + 
        (75.0 * Omg2dot * Omgdot * r2dot - 87.0 * Omgdot2 * r3dot) * rdot + 85.0 * Omg2dot2 * rdot2) 
        + u3 * Omgdot * (389.0 * Omgdot * r2dot * rdot2 + 590.0 * Omg2dot * rdot3) 
        + u4 * Omgdot2 * (-3.0 * r2dot + 1355.0 * rdot4)) 
        + invOmg9 * (375.0 * u3 * Omgdot3 * rdot3 
        + u * (6.5 * Omg2dot * Omgdot2 * r2dot - 7.5 * Omgdot3 * r3dot + 25.0 * Omg2dot2 * Omgdot * rdot) 
        + u2 * (16.5 * Omgdot3 * r2dot * rdot + 182.5 * Omg2dot * Omgdot2 * rdot2))
    );

    return FphiNewtNC;
}

/*========================================================================== */
/*(3,1) Mode                                                                 */
/*========================================================================== */
double eob_flx_FlmNewt_nc_31(double r, double Omg, double rdot, double r2dot, double r3dot, double Omgdot, double Omg2dot)
{
    const double invOmg  = 1.0 / Omg;
    const double invOmg2 = invOmg  * invOmg;
    const double invOmg3 = invOmg2 * invOmg;
    const double invOmg4 = invOmg3 * invOmg;
    const double invOmg5 = invOmg4 * invOmg;
    const double invOmg6 = invOmg5 * invOmg;
    const double invOmg7 = invOmg6 * invOmg;

    const double u   = 1.0 / r;
    const double u2  = u   * u;
    const double u3  = u2  * u;
    const double u4  = u3  * u;
    const double u5  = u4  * u;
    const double u6  = u5  * u;
    const double u7  = u6  * u;

    const double rdot2 = rdot  * rdot;
    const double rdot3 = rdot2 * rdot;
    const double rdot4 = rdot2 * rdot2;

    const double r2dot2 = r2dot  * r2dot;
    const double r3dot2 = r3dot  * r3dot;

    const double Omgdot2 = Omgdot  * Omgdot;
    const double Omgdot3 = Omgdot2 * Omgdot;
    const double Omg2dot2 = Omg2dot * Omg2dot;

    double FphiNewtNC = 1.0 
        + invOmg2 * u * (-24.0 * r2dot + 54.0 * u * rdot2) 
        + invOmg3 * (Omg2dot + 57.0 * u * Omgdot * rdot) 
        + invOmg4 * (21.0 * Omgdot2 - 18.0 * u5 * rdot2 - 288.0 * u3 * r2dot * rdot2 
        + u2 * (153.0 * r2dot2 - 102.0 * r3dot * rdot) + u4 * (6.0 * r2dot + 396.0 * rdot4)) 
        + invOmg5 * (u * (-3.0 * Omg2dot * r2dot - 48.0 * Omgdot * r3dot) 
        + 300.0 * u3 * Omgdot * rdot3 + u2 * rdot * (-117.0 * Omgdot * r2dot - 90.0 * Omg2dot * rdot)) 
        + invOmg6 * (-2.0 * Omg2dot2 + 324.0 * u7 * rdot4 + u * Omgdot * (-69.0 * Omgdot * r2dot - 24.0 * Omg2dot * rdot) 
        + u3 * r2dot * (-162.0 * r2dot2 + 198.0 * r3dot * rdot) + u2 * (30.0 * r3dot2 + 42.0 * Omgdot2 * rdot2) 
        + u4 * rdot2 * (756.0 * r2dot2 - 312.0 * r3dot * rdot) + u5 * r2dot * (-54.0 * r2dot + 180.0 * rdot4) 
        + u6 * rdot2 * (54.0 * r2dot + 120.0 * rdot4)) 
        + invOmg7 * (-3.0 * Omg2dot * Omgdot2 - 27.0 * u * Omgdot3 * rdot + 162.0 * u6 * Omgdot * rdot3 
        + u2 * (-18.0 * Omg2dot * r2dot2 + 42.0 * Omgdot * r2dot * r3dot + 6.0 * Omg2dot * r3dot * rdot) 
        + u3 * rdot * (90.0 * Omgdot * r2dot2 + (144.0 * Omg2dot * r2dot - 96.0 * Omgdot * r3dot) * rdot) 
        + u4 * (-6.0 * Omg2dot * r2dot + 480.0 * Omgdot * r2dot * rdot3 + 60.0 * Omg2dot * rdot4) 
        + u5 * rdot * (-54.0 * Omgdot * r2dot + 18.0 * Omg2dot * rdot + 240.0 * Omgdot * rdot4));

    return FphiNewtNC;
}

/*========================================================================== */
/*(4,4) Mode                                                                 */
/*========================================================================== */
double eob_flx_FlmNewt_nc_44(double r, double Omg, double rdot, double r2dot, double r3dot, double Omgdot, double Omg2dot)
{
    const double invOmg  = 1.0 / Omg;
    const double invOmg2 = invOmg  * invOmg;
    const double invOmg3 = invOmg2 * invOmg;
    const double invOmg4 = invOmg3 * invOmg;
    const double invOmg5 = invOmg4 * invOmg;
    const double invOmg6 = invOmg5 * invOmg;
    const double invOmg7 = invOmg6 * invOmg;
    const double invOmg8 = invOmg7 * invOmg;
    const double invOmg9 = invOmg8 * invOmg;

    const double u   = 1.0 / r;
    const double u2  = u   * u;
    const double u3  = u2  * u;
    const double u4  = u3  * u;
    const double u5  = u4  * u;
    const double u6  = u5  * u;
    const double u7  = u6  * u;
    const double u8  = u7  * u;
    const double u9  = u8  * u;
    const double u10 = u5  * u5;

    const double rdot2 = rdot  * rdot;
    const double rdot3 = rdot2 * rdot;
    const double rdot4 = rdot3 * rdot;
    const double rdot5 = rdot4 * rdot;
    const double rdot6 = rdot3 * rdot3;

    const double r2dot2 = r2dot  * r2dot;
    const double r2dot3 = r2dot2 * r2dot;
    const double r2dot4 = r2dot3 * r2dot;
    const double r3dot2 = r3dot  * r3dot;

    const double Omgdot2 = Omgdot  * Omgdot;
    const double Omgdot3 = Omgdot2 * Omgdot;
    const double Omgdot4 = Omgdot2 * Omgdot2;
    const double Omg2dot2 = Omg2dot * Omg2dot;

    double FphiNewtNC = 1.0 + 0.000244140625 * ( 
        invOmg2 * u * (-16032.0 * r2dot + 32768.0 * u * rdot2) 
        + invOmg3 * (-2880.0 * Omg2dot + 30080.0 * u * Omgdot * rdot) 
        + invOmg4 * (11456.0 * Omgdot2 - 2112.0 * u5 * rdot2 - 40560.0 * u3 * r2dot * rdot2 
        + u2 * (18324.5 * r2dot2 - 8192.0 * r3dot * rdot) + u4 * (704.0 * r2dot + 53312.0 * rdot4)) 
        + invOmg5 * u * (5666.0 * Omg2dot * r2dot - 4648.0 * Omgdot * r3dot + 74816.0 * u2 * Omgdot * rdot3 
        + u * rdot * (-14900.0 * Omgdot * r2dot - 3264.0 * Omg2dot * rdot)) 
        + invOmg6 * (434.0 * Omg2dot2 + 6432.0 * u7 * rdot4 
        + u * Omgdot * (-6448.0 * Omgdot * r2dot + 1928.0 * Omg2dot * rdot) 
        + u3 * r2dot * (-5817.0 * r2dot2 + 4997.0 * r3dot * rdot) 
        + u2 * (430.5 * r3dot2 + 46280.0 * Omgdot2 * rdot2) 
        + u4 * rdot2 * (20842.5 * r2dot2 - 11247.0 * r3dot * rdot) 
        + u5 * (-1166.0 * r2dot2 + 128.0 * r3dot * rdot - 5127.0 * r2dot * rdot4) 
        + u6 * rdot2 * (714.0 * r2dot + 13224.0 * rdot4)) 
        + invOmg7 * (208.0 * Omg2dot * Omgdot2 + 9312.0 * u * Omgdot3 * rdot + 6264.0 * u6 * Omgdot * rdot3 
        + u2 * (-1737.5 * Omg2dot * r2dot2 + 1934.5 * Omgdot * r2dot * r3dot + 184.5 * Omg2dot * r3dot * rdot) 
        + u3 * (4732.5 * Omgdot * r2dot2 * rdot + (6150.0 * Omg2dot * r2dot - 9097.0 * Omgdot * r3dot) * rdot2) 
        + u4 * (-188.0 * Omg2dot * r2dot + 48.0 * Omgdot * r3dot + 15126.0 * Omgdot * r2dot * rdot3 + 1936.0 * Omg2dot * rdot4) 
        + u5 * (-2328.0 * Omgdot * r2dot * rdot + 564.0 * Omg2dot * rdot2 + 23076.0 * Omgdot * rdot5)) 
        + invOmg8 * (510.0 * Omgdot4 - 119.0 * u * Omg2dot2 * r2dot + 162.0 * u10 * rdot4 
        + u2 * (1086.0 * Omgdot2 * r2dot2 + (1253.0 * Omg2dot * r2dot - 2026.0 * Omgdot * r3dot) * Omgdot * rdot + 553.0 * Omg2dot2 * rdot2) 
        + u3 * (-105.0 * r2dot * r3dot2 + 4320.0 * Omgdot2 * r2dot * rdot2 + 5296.0 * Omg2dot * Omgdot * rdot3) 
        + u5 * (-7.0 * r3dot2 + (576.0 * Omgdot2 + 108.0 * r2dot3) * rdot2 + 1170.0 * r2dot * r3dot * rdot3) 
        + u4 * (-192.0 * Omgdot2 * r2dot + 418.5 * r2dot4 - 405.0 * r2dot2 * r3dot * rdot + 426.0 * r3dot2 * rdot2 + 17724.0 * Omgdot2 * rdot4) 
        + u6 * (174.0 * r2dot3 + 97.0 * r2dot * r3dot * rdot + 3330.0 * r2dot2 * rdot4 - 882.0 * r3dot * rdot5) 
        + u9 * rdot2 * (-108.0 * r2dot + 1266.0 * rdot4) 
        + u7 * rdot2 * (-246.0 * r2dot2 - 457.0 * r3dot * rdot + 366.0 * r2dot * rdot4) 
        + u8 * r2dot2 * (18.0 + 105.0 * r2dot * rdot2 + 90.0 * rdot6)) 
        + invOmg9 * (804.0 * u8 * Omgdot * rdot5 + u * Omgdot * (85.0 * Omg2dot * Omgdot * r2dot - 85.0 * Omgdot2 * r3dot + 119.0 * Omg2dot2 * rdot) 
        + u2 * Omgdot2 * rdot * (255.0 * Omgdot * r2dot + 1135.0 * Omg2dot * rdot) 
        + u3 * (76.5 * Omg2dot * r2dot3 - 73.5 * Omgdot * r2dot2 * r3dot + (-3.0 * Omg2dot * r2dot + 102.0 * Omgdot * r3dot) * r3dot * rdot + 3126.0 * Omgdot3 * rdot3) 
        + u4 * (238.5 * Omgdot * r2dot3 * rdot + (193.5 * Omg2dot * r2dot + 354.0 * Omgdot * r3dot) * r2dot * rdot2 + 60.0 * Omg2dot * r3dot * rdot3) 
        + u5 * (17.0 * Omg2dot * r2dot2 + 7.0 * Omgdot * r2dot * r3dot - 7.0 * Omg2dot * r3dot * rdot + 1419.0 * Omgdot * r2dot2 * rdot3 + (651.0 * Omg2dot * r2dot - 189.0 * Omgdot * r3dot) * rdot4) 
        + u6 * (93.0 * Omgdot * r2dot2 * rdot + (-48.0 * Omg2dot * r2dot - 89.0 * Omgdot * r3dot) * rdot2 + 2019.0 * Omgdot * r2dot * rdot5 + 135.0 * Omg2dot * rdot6) 
        + u7 * rdot3 * (-207.0 * Omgdot * r2dot + 96.0 * Omg2dot * rdot + 450.0 * Omgdot * rdot4)));

    return FphiNewtNC;
}

/*========================================================================== */
/*(4,2) Mode                                                                 */
/*========================================================================== */
double eob_flx_FlmNewt_nc_42(double r, double Omg, double rdot, double r2dot, double r3dot, double Omgdot, double Omg2dot)
{
    const double invOmg  = 1.0 / Omg;
    const double invOmg2 = invOmg  * invOmg;
    const double invOmg3 = invOmg2 * invOmg;
    const double invOmg4 = invOmg3 * invOmg;
    const double invOmg5 = invOmg4 * invOmg;
    const double invOmg6 = invOmg5 * invOmg;
    const double invOmg7 = invOmg6 * invOmg;
    const double invOmg8 = invOmg7 * invOmg;
    const double invOmg9 = invOmg8 * invOmg;

    const double u   = 1.0 / r;
    const double u2  = u   * u;
    const double u3  = u2  * u;
    const double u4  = u3  * u;
    const double u5  = u4  * u;
    const double u6  = u5  * u;
    const double u7  = u6  * u;
    const double u8  = u7  * u;
    const double u9  = u8  * u;
    const double u10 = u5  * u5;

    const double rdot2 = rdot  * rdot;
    const double rdot3 = rdot2 * rdot;
    const double rdot4 = rdot3 * rdot;
    const double rdot5 = rdot4 * rdot;
    const double rdot6 = rdot3 * rdot3;
    const double rdot8 = rdot4 * rdot4;

    const double r2dot2 = r2dot  * r2dot;
    const double r2dot3 = r2dot2 * r2dot;
    const double r2dot4 = r2dot3 * r2dot;
    const double r3dot2 = r3dot  * r3dot;

    const double Omgdot2 = Omgdot  * Omgdot;
    const double Omgdot3 = Omgdot2 * Omgdot;
    const double Omgdot4 = Omgdot2 * Omgdot2;
    const double Omg2dot2 = Omg2dot * Omg2dot;

    double FphiNewtNC = 1.0 + 0.0625 * (
        invOmg2 * u * (-234.0 * r2dot + 512.0 * u * rdot2) 
        + invOmg3 * (-12.0 * Omg2dot + 536.0 * u * Omgdot * rdot) 
        + invOmg4 * (212.0 * Omgdot2 - 132.0 * u5 * rdot2 - 2544.0 * u3 * r2dot * rdot2 
        + u2 * (1038.5 * r2dot2 - 548.0 * r3dot * rdot) + u4 * (44.0 * r2dot + 3236.0 * rdot4)) 
        + invOmg5 * (u * (128.0 * Omg2dot * r2dot - 304.0 * Omgdot * r3dot) + 4124.0 * u3 * Omgdot * rdot3 
        + u2 * rdot * (-1496.0 * Omgdot * r2dot - 486.0 * Omg2dot * rdot)) 
        + invOmg6 * (2.0 * Omg2dot2 + 1608.0 * u7 * rdot4 
        + u * Omgdot * (-638.5 * Omgdot * r2dot - 184.0 * Omg2dot * rdot) 
        + u3 * r2dot * (-1389.0 * r2dot2 + 1343.0 * r3dot * rdot) 
        + u2 * (115.5 * r3dot2 + 2024.0 * Omgdot2 * rdot2) 
        + u4 * rdot2 * (4723.5 * r2dot2 - 2739.0 * r3dot * rdot) 
        + u5 * (-278.0 * r2dot2 + 32.0 * r3dot * rdot - 1491.0 * r2dot * rdot4) 
        + u6 * rdot2 * (138.0 * r2dot + 3444.0 * rdot4)) 
        + invOmg7 * (-11.0 * Omg2dot * Omgdot2 + 222.0 * u * Omgdot3 * rdot + 1404.0 * u6 * Omgdot * rdot3 
        + u2 * (-297.5 * Omg2dot * r2dot2 + 506.5 * Omgdot * r2dot * r3dot + 130.5 * Omg2dot * r3dot * rdot) 
        + u3 * (1456.5 * Omgdot * r2dot2 * rdot + (1128.0 * Omg2dot * r2dot - 2089.0 * Omgdot * r3dot) * rdot2) 
        + u4 * (-20.0 * Omg2dot * r2dot + 12.0 * Omgdot * r3dot + 2694.0 * Omgdot * r2dot * rdot3 + 16.0 * Omg2dot * rdot4) 
        + u5 * rdot * (-528.0 * Omgdot * r2dot + 60.0 * Omg2dot * rdot + 4836.0 * Omgdot * rdot4)) 
        + invOmg8 * (6.0 * Omgdot4 + u * Omg2dot * (-17.0 * Omg2dot * r2dot - 3.0 * Omgdot * r3dot) 
        + 162.0 * u10 * rdot4 
        + u2 * (409.5 * Omgdot2 * r2dot2 + (170.0 * Omg2dot * r2dot - 415.0 * Omgdot * r3dot) * Omgdot * rdot + 88.0 * Omg2dot2 * rdot2) 
        + u3 * (-105.0 * r2dot * r3dot2 + 615.0 * Omgdot2 * r2dot * rdot2 + 766.0 * Omg2dot * Omgdot * rdot3) 
        + u5 * (-7.0 * r3dot2 + (63.0 * Omgdot2 + 108.0 * r2dot3) * rdot2 + 1170.0 * r2dot * r3dot * rdot3) 
        + u4 * (-21.0 * Omgdot2 * r2dot + 418.5 * r2dot4 - 405.0 * r2dot2 * r3dot * rdot + 426.0 * r3dot2 * rdot2 + 3129.0 * Omgdot2 * rdot4) 
        + u6 * (174.0 * r2dot3 + 97.0 * r2dot * r3dot * rdot + 3330.0 * r2dot2 * rdot4 - 882.0 * r3dot * rdot5) 
        + u9 * rdot2 * (-108.0 * r2dot + 1266.0 * rdot4) 
        + u7 * rdot2 * (-246.0 * r2dot2 - 457.0 * r3dot * rdot + 366.0 * r2dot * rdot4) 
        + u8 * (18.0 * r2dot2 + 105.0 * r2dot * rdot4 + 90.0 * rdot8)) 
        + invOmg9 * (804.0 * u8 * Omgdot * rdot5 
        + u * Omgdot * (11.5 * Omg2dot * Omgdot * r2dot - 8.5 * Omgdot2 * r3dot + 14.0 * Omg2dot2 * rdot) 
        + u2 * Omgdot2 * rdot * (43.5 * Omgdot * r2dot + 137.5 * Omg2dot * rdot) 
        + u3 * (76.5 * Omg2dot * r2dot3 - 73.5 * Omgdot * r2dot2 * r3dot + (-3.0 * Omg2dot * r2dot + 102.0 * Omgdot * r3dot) * r3dot * rdot + 381.0 * Omgdot3 * rdot3) 
        + u4 * (238.5 * Omgdot * r2dot3 * rdot + (193.5 * Omg2dot * r2dot + 354.0 * Omgdot * r3dot) * r2dot * rdot2 + 60.0 * Omg2dot * r3dot * rdot3) 
        + u5 * (17.0 * Omg2dot * r2dot2 + 7.0 * Omgdot * r2dot * r3dot - 7.0 * Omg2dot * r3dot * rdot + 1419.0 * Omgdot * r2dot2 * rdot3 + (651.0 * Omg2dot * r2dot - 189.0 * Omgdot * r3dot) * rdot4) 
        + u6 * (93.0 * Omgdot * r2dot2 * rdot + (-48.0 * Omg2dot * r2dot - 89.0 * Omgdot * r3dot) * rdot2 + 2019.0 * Omgdot * r2dot * rdot5 + 135.0 * Omg2dot * rdot6) 
        + u7 * rdot3 * (-207.0 * Omgdot * r2dot + 96.0 * Omg2dot * rdot + 450.0 * Omgdot * rdot4)) 
    );

    return FphiNewtNC;
}


/**
 *  Function: eob_flx_flm_nc_22
 *  ---------------------------
 *  Computes the non-circular correction factor for the (2,2) mode of the gravitational
 *  wave flux in the Effective One Body (EOB) formalismm up to 2PN order.
 *  Depending on `use_pade_22`, it either resums the series using a Padé approximant or returns the Taylor-expanded form.
 * 
 *  @param[in] r The radial separation between the two bodies.
 *  @param[in] prstar The radial momentum conjugate to the tortoise coordinate.
 *  @param[in] prsdot The time derivative of the radial momentum.
 *  @return The non-circular correction factor for the (2,2) mode of the gravitational wave flux.
 */
#define use_pade_22 (1)
double eob_flx_flm_nc_22(double r, double prstar, double prsdot) {

  const double nu = EOBPars -> nu;
  const double nu2 = nu * nu;
  const double nu3 = nu2 * nu;
  
  const double r2  = r * r;
  const double r3  = r2 * r;
  const double r4  = r3 * r;
  const double r5  = r4 * r;
  const double r6  = r5 * r;
  const double r7  = r6 * r;
  const double r8  = r7 * r;
  const double r9  = r8 * r;
  const double r10 = r9 * r;
  const double r11 = r10 * r;

  const double logr = log(r);
  const double sqrt_r = sqrt(r);
  
  const double prstar2 = prstar  * prstar;
  const double prstar3 = prstar2 * prstar;
  const double prstar4 = prstar3 * prstar;
  const double prstar6 = prstar4 * prstar2;

  const double prsdot2 = prsdot  * prsdot;
  const double prsdot3 = prsdot2 * prsdot;
  const double prsdot4 = prsdot3 * prsdot;
  const double prsdot5 = prsdot4 * prsdot;
  const double prsdot6 = prsdot5 * prsdot;

  double f22_15PN = 0.0, f22_2PN = 0.0, f22_25PN = 0.0, f22_3PN = 0.0;

  const double f22_1PN = (2.6964285714285716 - 0.4642857142857143 * nu) * prstar2 
    + (0.6622023809523809 - 0.5803571428571429 * nu) * prstar4 * r 
    + (0.16555059523809523 - 0.14508928571428573 * nu) * prstar6 * r2 
    + prsdot * r * (2.2976190476190474 - 0.30952380952380953 * nu 
		    + (1.1607142857142858 - 0.8571428571428571 * nu) * prstar2 * r 
		    + (-0.04092261904761905 + 0.07589285714285714 * nu) * prstar4 * r2) 
    + prsdot2 * r3 * (0.2976190476190476 - 0.14285714285714285 * nu 
		      + (-0.5059523809523809 + 0.39285714285714285 * nu) * prstar2 * r 
		      + (-0.10602678571428571 + 0.060267857142857144 * nu) * prstar4 * r2) 
    + prsdot3 * r5 * (-0.1488095238095238 + 0.07142857142857142 * nu 
		      + (0.21577380952380953 - 0.17857142857142858 * nu) * prstar2 * r) 
    + prsdot4 * r7 * (0.0744047619047619 - 0.03571428571428571 * nu 
		      + (-0.08928571428571429 + 0.08035714285714286 * nu) * prstar2 * r) 
    + prsdot5 * r9 * (-0.03720238095238095 + 0.017857142857142856 * nu) 
    + prsdot6 * r11 * (0.018601190476190476 - 0.008928571428571428 * nu);

  if (EOBPars->flux_nc_order >= FLUX_NC_ORDER_15PN)
    f22_15PN = (-1.875 * prstar2 - 0.609375 * prstar4 * r - 0.24231770833333333 * prstar6 * r2
		+ prsdot * r * (-1.25 - 0.34375 * prstar2 * r + 0.04296875 * prstar4 * r2) 
		+ prsdot2 * r3 * (0.125 + 0.234375 * prstar2 * r + 0.37890625 * prstar4 * r2) 
		+ prsdot3 * r5 * (-0.041666666666666667 + 0.13671875 * prstar2 * r) 
		+ prsdot4 * r7 * (0.0625 - 0.3323567708333333 * prstar2 * r) 
		- 0.074609375 * prsdot5 * r9 + 0.07254774305555555 * prsdot6 * r11) / sqrt_r * M_PI;

  if (EOBPars->flux_nc_order >= FLUX_NC_ORDER_2PN)
    f22_2PN = (0.6313303099017384 - 5.15239984882842 * nu + 1.028533635676493 * nu2) * prstar2 / r 
      + (5.323932350718065 + 0.8621267951625095 * nu + 0.31949168556311414 * nu2) * prstar4 
      + (0.34901738473167043 + 1.2075243291761149 * nu - 0.08792753212396069 * nu2) * prstar6 * r 
      + prsdot * (19.462915721844293 + 6.640495086923658 * nu - 2.512660619803477 * nu2 
		  + (-1.7667824074074074 - 0.25462962962962965 * nu + 1.460978835978836 * nu2) * prstar2 * r 
		  + (2.817982922335601 - 0.3347328514739229 * nu + 0.061685090702947844 * nu2) * prstar4 * r2) 
      + prsdot2 * r2 * (14.228670634920634 + 0.8015873015873016 * nu - 2.3333333333333333 * nu2 
			+ (-1.3015518707482994 + 2.5677437641723357 * nu + 0.2128684807256236 * nu2) * prstar2 * r 
			+ (-0.9038362740929705 - 0.37875566893424034 * nu + 0.40072278911564624 * nu2) * prstar4 * r2) 
      + prsdot3 * r4 * (-0.5165934429327287 + 0.7197656840513983 * nu - 0.002928949357520786 * nu2 
			+ (1.6302585270219199 - 2.389644746787604 * nu + 0.24104780801209372 * nu2) * prstar2 * r) 
      + prsdot4 * r6 * (0.518713388133031 - 0.596489984882842 * nu + 0.055035903250188964 * nu2 
			+ (-1.1746622260015118 + 1.6294052343159486 * nu - 0.26747921390778534 * nu2) * prstar2 * r) 
      + prsdot5 * r8 * (-0.38956502739984883 + 0.41654856386999245 * nu - 0.05430366591080877 * nu2) 
      + prsdot6 * r10 * (0.25988668036659107 - 0.26742606764928195 * nu + 0.040544690098261525 * nu2);

  /* 2.5PN: pi-dependent contribution */
  const double f22_25pi = (-0.4494047619047619 - 0.15178571428571427 * nu) * prstar4
    + (9.732142857142858 - 5.571428571428571 * nu) * prstar2 / r + (2.202188740079365 - 0.24211309523809524 * nu) * prstar6 * r
    + (-7.4466362847222225 + 0.4365234375 * nu) * prsdot5 * r8 + (9.099389260912698 - 0.4668526785714286 * nu) * prsdot6 * r10
    + prsdot * (-5.446428571428571 + 5.339285714285714 * nu + (-3.4174107142857144 - 3.122767857142857 * nu) * prstar2 * r
		+ (-0.7711123511904762 - 0.6710379464285714 * nu) * prstar4 * r2)
    + prsdot2 * ((-1.8988095238095237 + 0.6547619047619048 * nu) * r2 + (-11.482142857142858 + 0.7693452380952381 * nu) * prstar2 * r3
		 + (-14.419642857142858 + 0.9197358630952381 * nu) * prstar4 * r4)
    + prsdot3 * ((-2.7038690476190474 + 0.05952380952380952 * nu) * r4 + (13.877139136904763 - 0.44131324404761907 * nu) * prstar2 * r5)
    + prsdot4 * ((5.29265873015873 - 0.30853174603174605 * nu) * r6 + (-10.659365699404763 + 0.36904761904761907 * nu) * prstar2 * r7);

  /* 2.5PN: rational contribution */
  const double f22_25pr = (37.99357326333484 + 1883.174823500558 * nu + 1.1719635770975056 * nu2) * prstar4
    + (-91.37281999328036 + 157.01212938094432 * nu) / r2
    + (9.120427201923917 + 875.6339729375509 * nu + 0.9207057823129252 * nu2) * prstar2 / r
    + (49.657475269341376 + 3717.314263188837 * nu - 151.835689484127 * nu2) * prsdot4 * r6
    + (-68.54681856877869 - 4509.771126570733 * nu + 266.64274199263036 * nu2) * prsdot5 * r8
    + prsdot * (-48.63778403743601 - 194.37156309306687 * nu + 8.080498866213151 * nu2
		+ (27.841346826311565 - 234.7320469828701 * nu + 6.187322845804989 * nu2) * prstar2 * r
		+ (-31.541446593764704 + 1416.2744536327664 * nu - 0.4717616213151927 * nu2) * prstar4 * r2)
    + prsdot2 * ((14.273003592661553 + 1177.5635757053863 * nu - 29.32405045351474 * nu2) * r2
		 + (-77.8468363773136 + 4033.0575850398454 * nu - 19.536706349206348 * nu2) * prstar2 * r3)
    + prsdot3 * ((-29.406835329275054 - 1476.0992116718098 * nu + 74.7610189909297 * nu2) * r4
		 + (127.10547177464468 - 1122.2895707174903 * nu + 42.06435126133787 * nu2) * prstar2 * r5);
    
  /* 2.5PN: logarithmic contribution */
  const double f22_25log = (-13.133928571428571 - 2.4575892857142856 * nu) * prstar4
    + (53.67857142857143 - 14.785714285714286 * nu) / r2 + (9.107142857142858 + 4.741071428571429 * nu) * prstar2 / r
    + (-9.207589285714286 - 2.986607142857143 * nu) * prsdot4 * r6 + (10.322544642857142 + 3.0870535714285716 * nu) * prsdot5 * r8
    + prsdot * (24.696428571428573 - 27.964285714285715 * nu + (-8.620535714285714 - 7.888392857142857 * nu) * prstar2 * r
		+ (4.675223214285714 + 6.716517857142857 * nu) * prstar4 * r2)
    + prsdot2 * ((-8.116071428571429 - 8.089285714285714 * nu) * r2 + (9.441964285714286 + 5.564732142857143 * nu) * prstar2 * r3)
    + prsdot3 * ((8.004464285714286 + 4.205357142857143 * nu) * r4 + (-16.295758928571427 - 5.815848214285714 * nu) * prstar2 * r5);

  /* Combine 2.5PN contributions */
  if (EOBPars->flux_nc_order >= FLUX_NC_ORDER_25PN)
    f22_25PN = f22_25pr * prstar + f22_25pi * M_PI / sqrt_r + f22_25log * prstar * logr;
  else
    f22_25PN = 0.0;

  /* 3PN: rational contribution */
  const double f22_3 = (4.86195068131201 + 7.9255356977402585 * nu + 3.6513035378894223 * nu2 - 1.16366773756761 * nu3) * prstar6
    + (23.98279534644493 - 376.908949620374 * nu + 46.647747314413984 * nu2 + 7.312542946471518 * nu3) * prstar2 / r2
    + (-40.84569365987045 - 57.82256732262956 * nu + 34.22599987565974 * nu2 - 3.1379291886307192 * nu3) * prstar4 / r
    + (1.198430168995401 + 7.346149588306773 * nu - 1.6216392213150292 * nu2 - 0.10431066516313114 * nu3) * prsdot5 * r7
    + (-0.06060499500017613 - 5.5881126443824485 * nu + 1.2855057524753528 * nu2 + 0.021472226005601856 * nu3) * prsdot6 * r9
    + prsdot * ((83.39078369629715 - 30.770945222107354 * nu - 33.4151099082009 * nu2 + 16.3665120389185 * nu3) * prstar2
		+ (-49.49757563333576 + 163.0639346166557 * nu - 22.259952026746245 * nu2 - 2.8804850824408645 * nu3) / r
		+ (39.039045613026246 - 94.6921930189791 * nu + 34.875503962783405 * nu2 - 2.4664518397245865 * nu3) * prstar4 * r)
    + prsdot2 * ((17.2614688265503 + 31.543513595074028 * nu - 0.8937698576304018 * nu2 - 1.4995425178248307 * nu3) * r
		 + (64.63727837255283 + 67.80545927955623 * nu - 22.011939834856502 * nu2 + 10.495045408735885 * nu3) * prstar2 * r2
		 + (8.7628044616431 + 62.244551162052176 * nu - 20.776682387120314 * nu2 + 2.2773965511529286 * nu3) * prstar4 * r3)
    + prsdot3 * ((32.26528655677794 + 34.104897728401056 * nu - 4.446300805675806 * nu2 - 0.3066528379028379 * nu3) * r3
		 + (-5.649879188742489 - 59.2944662467264 * nu + 13.435651141619466 * nu2 + 1.8244519263333039 * nu3) * prstar2 * r4)
    + prsdot4 * ((-1.1651089297229191 - 9.691648393085579 * nu + 1.8869445896444195 * nu2 + 0.25099661377297433 * nu3) * r5
		 + (-14.785881293937564 + 44.92075609343021 * nu - 9.884160942515534 * nu2 - 0.8082175977052678 * nu3) * prstar2 * r6);

  /* 3PN: pi-dependent contribution */
  const double f22_3pi = (-5.622593047730334 + 0.11908895502645503 * nu) * prstar4 + 57.205329415187485 / r2
    + (-33.16618179338275 - 0.12301587301587301 * nu) * prstar2 / r + (-3.948689405893674 - 10.447627314814815 * nu) * prsdot4 * r6
    + (4.21026488248026 + 21.076554232804234 * nu) * prsdot5 * r8
    + prsdot * (55.66614779681274 + 0.25396825396825395 * nu + (9.486688022974704 + 1.072089947089947 * nu) * prstar2 * r
		+ (1.9690664978697896 - 0.1371527777777778 * nu) * prstar4 * r2)
    + prsdot2 * ((0.47608117216259416 - 1.2797619047619047 * nu) * r2 + (10.831421781564131 - 2.9187334656084656 * nu) * prstar2 * r3)
    + prsdot3 * ((2.568091360917606 + 4.266534391534392 * nu) * r4 + (-14.013567832065746 + 6.029017857142857 * nu) * prstar2 * r5);

  /* 3PN: pi-squared contribution */
  const double f22_3pi2 = (0.7419704861111112 - 0.150146484375 * nu) * prstar6 + (4.25 + 12.97265625 * nu) * prstar2 / r2
    + (2.015625 - 0.6005859375 * nu) * prstar4 / r + (-0.03958333333333333 - 0.06005859375 * nu) * prsdot5 * r7
    + (0.058203125 + 0.030029296875 * nu) * prsdot6 * r9
    + prsdot * ((4.5 + 4.083984375 * nu) * prstar2 + (2.4166666666666665 - 6.40625 * nu) / r
		+ (1.21875 + 1.3212890625 * nu) * prstar4 * r)
    + prsdot2 * ((1.125 - 2.08203125 * nu) * r + (-0.4375 - 1.921875 * nu) * prstar2 * r2
		 + (-1.0983072916666667 - 1.14111328125 * nu) * prstar4 * r3)
    + prsdot3 * ((-0.4166666666666667 - 1.0944010416666667 * nu) * r3 + (-0.890625 + 0.90087890625 * nu) * prstar2 * r4)
    + prsdot4 * ((0.09895833333333333 + 0.1201171875 * nu) * r5 + (0.580078125 - 0.42041015625 * nu) * prstar2 * r6);

  /* 3PN: squared-logarithmic contribution */
  const double f22_3log2 = 10.828125 * prstar6 + 33.75 * prstar2 / r2 - 1.6875 * prstar4 / r
    + 0.421875 * prsdot5 * r7 - 0.2109375 * prsdot6 * r9
    + prsdot * (1.6875 * prstar2 - 20.25 / r + 23.765625 * prstar4 * r)
    + prsdot2 * (-3.375 * r + 9.5625 * prstar2 * r2 - 18.2109375 * prstar4 * r3)
    + prsdot3 * (1.6875 * r3 - 7.734375 * prstar2 * r4) + prsdot4 * (-0.84375 * r5 + 8.71875 * prstar2 * r6);

  /* 3PN: pi-logarithmic contribution */
  const double f22_3pilog = -1.44140625 * prstar4 - 20.25 / r2 + 7.78125 * prstar2 / r
    + 3.908203125 * prsdot4 * r6 - 3.928125 * prsdot5 * r8
    + prsdot * (-18.0 - 7.5 * prstar2 * r + 1.94296875 * prstar4 * r2)
    + prsdot2 * (-1.125 * r2 - 3.64453125 * prstar2 * r3)
    + prsdot3 * (-2.34375 * r4 + 4.234375 * prstar2 * r5);

  /* 3PN: logarithmic contribution */
  const double f22_3log = (-30.493746451684274 + 0.036458333333333336 * nu) * prstar6
    + (-180.39350437457028 - 0.30952380952380953 * nu) * prstar2 / r2 + (0.34609485849068733 + 0.37797619047619047 * nu) * prstar4 / r
    + (-2.0988707940123277 - 2.3206845238095237 * nu) * prsdot5 * r7 + (1.1425101077766158 + 3.255580357142857 * nu) * prsdot6 * r9
    + prsdot * ((-52.91481994667629 + 1.5863095238095237 * nu) * prstar2 + 53.04005759058326 / r
		+ (-71.99938367959112 - 0.9159226190476191 * nu) * prstar4 * r)
    + prsdot2 * ((-11.034657040893933 + 0.30952380952380953 * nu) * r + (-16.58582343067974 - 3.462797619047619 * nu) * prstar2 * r2
		 + (23.02829445619136 + 2.2239583333333335 * nu) * prstar4 * r3)
    + prsdot3 * ((-13.564812622097634 - 0.8125 * nu) * r3 + (7.687301731901243 + 6.086309523809524 * nu) * prstar2 * r4)
    + prsdot4 * ((3.173409719231131 + 1.4985119047619047 * nu) * r5 + (1.3151836437173188 - 9.618675595238095 * nu) * prstar2 * r6);

  /* Combine 3PN contributions */
  if (EOBPars->flux_nc_order >= FLUX_NC_ORDER_25PN)
    f22_3PN = f22_3 + f22_3pi * M_PI * prstar / sqrt_r + (f22_3log + f22_3pilog * M_PI * prstar / sqrt_r) * logr
      + f22_3pi2 * M_PI * M_PI + f22_3log2 * logr * logr;
  else
    f22_25PN = 0.0;

  double a2 = f22_1PN;
  double a3 = f22_15PN;
  double a4 = f22_2PN;
  double a5 = f22_25PN;
  double a6 = f22_3PN;
  
#if use_pade_22
  double den = 1.0 - a2 - a3;
  if (EOBPars->flux_nc_order >= FLUX_NC_ORDER_2PN)
    den += (a2*a2-a4);
  if (EOBPars->flux_nc_order >= FLUX_NC_ORDER_3PN)
    den += (-a2*a2*a2 + 2.*a2*a3 + a3*a3 + 2.*a2*a4 - a5 -a6);
  return 1.0/den;
#else
  return 1.0 + a2 + a3 + a4;
#endif
}

/**
 * Function: eob_flx_flm_nc_21
 * ---------------------------
 * Computes the non-circular correction factor for the (2,1) mode of the gravitational
 * wave flux in the Effective One Body (EOB) formalismm up to 1PN order.
 * 
 *  @param[in] r The radial separation between the two bodies.
 *  @param[in] prstar The radial momentum conjugate to the tortoise coordinate.
 *  @param[in] prsdot The time derivative of the radial momentum.
 *  @return The non-circular correction factor for the (2,1) mode of the gravitational wave flux.
 */
double eob_flx_flm_nc_21(double r, double prstar, double prsdot) {
    const double nu = EOBPars -> nu;
    const double prstar2 = prstar * prstar;

    const double f21_1PN = (0.21428571428571427 - 0.42857142857142855 * nu) * prstar2 
        + prsdot * r * (4.642857142857143 + 0.38095238095238093 * nu);

    return 1.0 + f21_1PN;
}

/**
 * Function: eob_flx_flm_nc_33
 * ---------------------------
 * Computes the non-circular correction factor for the (3,3) mode of the gravitational
 * wave flux in the Effective One Body (EOB) formalismm up to 1PN order.
 * 
 *  @param[in] r The radial separation between the two bodies.
 *  @param[in] prstar The radial momentum conjugate to the tortoise coordinate.
 *  @param[in] prsdot The time derivative of the radial momentum.
 *  @return The non-circular correction factor for the (3,3) mode of the gravitational wave flux.
 */
double eob_flx_flm_nc_33(double r, double prstar, double prsdot) {
    const double nu = EOBPars -> nu;
    const double r2 = r * r;
    const double r3 = r2 * r;
    const double r5 = r3 * r2;
    const double r7 = r5 * r2;
    const double r9 = r7 * r2;
    const double r11 = r9 * r2;

    const double prstar2 = prstar  * prstar;
    const double prstar4 = prstar2 * prstar2;
    const double prstar6 = prstar4 * prstar2;

    const double prsdot2 = prsdot  * prsdot;
    const double prsdot3 = prsdot2 * prsdot;
    const double prsdot4 = prsdot3 * prsdot;
    const double prsdot5 = prsdot4 * prsdot;
    const double prsdot6 = prsdot5 * prsdot;

    const double f33_1PN = (4.091449474165524 - 0.7288523090992227 * nu) * prstar2 
        + (0.38950325624105026 - 0.3167689357802654 * nu) * prstar4 * r 
        + (0.25195551286395695 - 0.2050987654398423 * nu) * prstar6 * r2 
        + prsdot * r * (3.9670781893004117 - 1.4650205761316872 * nu 
        + (0.9247743399549526 - 0.584700841673864 * nu) * prstar2 * r 
        + (0.5468741277645747 - 0.4043296120045938 * nu) * prstar4 * r2) 
        + prsdot2 * r3 * (0.3682365493065082 - 0.16491388507849414 * nu 
        + (0.0666113453798258 - 0.07074350680508279 * nu) * prstar2 * r 
        + (-0.3150000050720085 + 0.26821867957530765 * nu) * prstar4 * r2) 
        + prsdot3 * r5 * (-0.12453499071392685 + 0.07311441909826302 * nu 
        + (-0.16334396294667833 + 0.155088514175098253 * nu) * prstar2 * r) 
        + prsdot4 * r7 * (0.06405247451948778 - 0.04731203568327539 * nu 
        + (0.1418970418297452 - 0.14408840072128107 * nu) * prstar2 * r) 
        + prsdot5 * r9 * (-0.045222577270558346 + 0.036976041295534064 * nu) 
        + prsdot6 * r11 * (0.03644723802353617 - 0.030758819492607913 * nu);

    return 1.0 + f33_1PN;
}

/**
 * Function: eob_flx_flm_nc_31
 * ---------------------------
 * Computes the non-circular correction factor for the (3,1) mode of the gravitational
 * wave flux in the Effective One Body (EOB) formalismm up to 1PN order.
 * 
 *  @param[in] r The radial separation between the two bodies.
 *  @param[in] prstar The radial momentum conjugate to the tortoise coordinate.
 *  @param[in] prsdot The time derivative of the radial momentum.
 *  @return The non-circular correction factor for the (3,1) mode of the gravitational wave flux.
 */
double eob_flx_flm_nc_31(double r, double prstar, double prsdot) {
    const double nu = EOBPars -> nu;
    const double r2 = r * r;
    const double r3 = r2 * r;
    const double r5 = r3 * r2;
    const double r7 = r5 * r2;
    const double r9 = r7 * r2;
    const double r11 = r9 * r2;

    const double prstar2 = prstar  * prstar;
    const double prstar4 = prstar2 * prstar2;
    const double prstar6 = prstar4 * prstar2;

    const double prsdot2 = prsdot  * prsdot;
    const double prsdot3 = prsdot2 * prsdot;
    const double prsdot4 = prsdot3 * prsdot;
    const double prsdot5 = prsdot4 * prsdot;
    const double prsdot6 = prsdot5 * prsdot;

    const double f31_1PN = (68.0 + 354.0 * nu) * prstar2 
        + (-3194.0 - 17240.0 * nu) * prstar4 * r 
        + (154936.0 + 835960.0 * nu) * prstar6 * r2 
        + prsdot * r * (-5.333333333333333 - 41.333333333333333 * nu 
        + (1575.0 + 8370.0 * nu) * prstar2 * r + (-132734.0 - 711080.0 * nu) * prstar4 * r2) 
        + prsdot2 * r3 * (-89.33333333333333 - 423.3333333333333 * nu 
        + (27224.0 + 140900.0 * nu) * prstar2 * r + (-3420720.0 - 18075420.0 * nu) * prstar4 * r2) 
        + prsdot3 * r5 * (-897.6666666666666 - 4226.666666666667 * nu 
        + (405698.0 + 2063540.0 * nu) * prstar2 * r) 
        + prsdot4 * r7 * (-9364.333333333334 - 43913.333333333333 * nu 
        + (5587272.0 + 28066140.0 * nu) * prstar2 * r) 
        + prsdot5 * r9 * (-99947.66666666667 - 467566.6666666667 * nu) 
        + prsdot6 * r11 * (-1081064.3333333333 - 5050353.333333333 * nu);

    return 1.0 + f31_1PN;
}
