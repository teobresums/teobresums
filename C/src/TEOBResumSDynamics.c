/** \file TEOBResumSDynamics.c
 *  \brief Implementation of the EOB and spins dynamics
 *  
 *  This file contains the implementation of the EOB and spins dynamics,
 *  as well as routines to compute the EOB hamiltonian and the gyro-gravitomagnetic
 *  functions.  
 */

#include "TEOBResumS.h"

/**
 * Function: eob_dyn_rhs
 * ---------------------
 *   Computes the r.h.s. of the EOB Hamiltonian dynamics, no
 *   spins version
 * 
 *  @param[in]  t     : time
 *  @param[in]  y     : EOB variables
 *  @param[out] dy    : r.h.s. of EOB variables
 *  @param[in]  d     : dummy pointer to Dynamics
 * 
 *  @return GSL_SUCCESS
*/
int eob_dyn_rhs(double t, const double y[], double dy[], void *d)
{
  
  (void)(t); /* avoid unused parameter warning */
  Dynamics *dyn = d;

  const double nu = EOBPars->nu;
  const double z3 = 2.0*nu*(4.0-3.0*nu);

  /* Unpack y */ 
  const double phi    = y[EOB_EVOLVE_PHI];
  const double r      = y[EOB_EVOLVE_RAD];
  const double pphi   = y[EOB_EVOLVE_PPHI];
  const double prstar = y[EOB_EVOLVE_PRSTAR];
  
  if ((EOBPars->use_tidal)&&(EOBPars->use_tidal_fmode_model))
    /* Update the dressing factors for the f-mode resonances */
    fmode_resonance_dressing_factors(r, dyn);
  
  /* Compute EOB Metric */
  double A, B, Q, dA, d2A, dB, dQ, dQ_dprstar, ddQ_drdprstar, d2Q_dprstar2, pl_hold;
  eob_metric(r, prstar, dyn, &A, &B, &dA, &d2A, &dB, &pl_hold, &Q, &dQ, &dQ_dprstar, &pl_hold, &ddQ_drdprstar, &d2Q_dprstar2, &pl_hold, &pl_hold, &pl_hold);

  /* Failsafe to avoid event horizon crossing */
  if (A < 1e-6){
    A = fabs(A);
    B = fabs(B);
    dyn->ode_stop = true;
  }

  /* Compute Hamiltonian */
  double H, Heff, dHeff_dr,dHeff_dprstar;
  eob_ham(nu, r, pphi, prstar, A, dA, Q, dQ, dQ_dprstar, &H, &Heff, &dHeff_dr, &dHeff_dprstar, NULL);
  double E = nu*H;

  /* Shorthands */
  const double u  = 1./r;
  const double u2 = u*u;
  const double u3 = u2*u;
  const double pphi2    = SQ(pphi);
  const double prstar2  = prstar*prstar;
  const double prstar3  = prstar2*prstar;
  const double prstar4  = prstar3*prstar;
  const double sqrtAbyB = sqrt(A/B);
  const double divHE    = 1./(Heff*E);
  const double Omega    = A*pphi*u2*divHE;

  /* d\phi/dt */
  dy[EOB_EVOLVE_PHI] = Omega;
  
  /* dr/dt (conservative part of) */
  dy[EOB_EVOLVE_RAD] = sqrtAbyB*dHeff_dprstar/E;
  
  /* dp_{r*}/dt (conservative part of) */
  dy[EOB_EVOLVE_PRSTAR] = - sqrtAbyB*dHeff_dr/E;       
  
  /* Compute flux */
  const double sqrtW = sqrt(A*(1. + pphi2*u2));
  const double psi   = 2.*(1.0 + 2.0*nu*(sqrtW - 1.0))/(SQ(r)*dA);
  /*const double psi = 2.*(1.0 + 2.0*nu*(Heff - 1.0))/(r2*dA); */
  const double r_omega = r*cbrt(psi);
  const double v_phi   = r_omega*Omega;
  const double x       = v_phi * v_phi;
  const double jhat    = pphi/(r_omega*v_phi);  
  const double tmpE    = 1./Heff+nu/(E*E);
  const double dprstar_dt    = dy[EOB_EVOLVE_PRSTAR];
  const double dr_dt         = dy[EOB_EVOLVE_RAD];
  const double ddotr_dr      = sqrtAbyB*( (prstar + 0.5*A*dQ_dprstar)*(0.5*(dA/A - dB/B) - dHeff_dr*tmpE) + 0.5*dA*dQ_dprstar + 0.5*A*ddQ_drdprstar )*divHE;
  const double ddotr_dprstar = sqrtAbyB*( 1. + 0.5*A*d2Q_dprstar2 - (prstar + 0.5*A*dQ_dprstar)*tmpE*dHeff_dprstar )*divHE;
  
  /* Approximate ddot(r) without Flux */
  const double ddotr = dprstar_dt*ddotr_dprstar + dr_dt*ddotr_dr;
  
  /* Compute flux and dp_{\phi}/dt */
  if (dyn->noflx) dy[EOB_EVOLVE_PPHI] = 0.;
  else            dy[EOB_EVOLVE_PPHI] = eob_flx_Flux(x,Omega,r_omega,E,Heff,jhat,r, prstar,ddotr,dyn);

  if(dyn->store) {
    /* Store values */
    dyn->t = t;
    dyn->r = r;
    dyn->phi = y[EOB_EVOLVE_PHI];
    dyn->pphi = pphi;
    dyn->prstar = prstar;
    dyn->Omg = Omega;
    dyn->Omg_orb = Omega;
    dyn->H = H;
    dyn->E = E;
    dyn->Heff = Heff;
    dyn->A = A;
    dyn->dA = dA;
    dyn->d2A = d2A;
    dyn->B = B;
    dyn->dB = dB;
    dyn->psi = psi;
    dyn->r_omega = r_omega;
    dyn->v_phi = v_phi;
    dyn->jhat = jhat;
    dyn->ddotr = ddotr;
  }
  
  return GSL_SUCCESS;

}

/**
  * Function: eob_ham
  * -----------------
  *   Computes the EOB nonspinning Hamiltonian
  * 
  *  @param[in]  nu            : symmetric mass ratio
  *  @param[in]  r             : radial separation
  *  @param[in]  pphi          : orbital angular momentum
  *  @param[in]  prstar        : (tortoise) radial momentum
  *  @param[in]  A             : A radial potential
  *  @param[in]  dA            : drvt A,r
  *  @param[in]  Q             : Q potential
  *  @param[in]  dQ            : drvt Q,r
  *  @param[in]  dQ_dprstar    : drvt Q,prstar
  *  @param[out] H             : real EOB Hamiltonian divided by mu=m1m2/(m1+m2)
  *  @param[out] Heff          : effective EOB Hamiltonian (divided by mu)
  *  @param[out] dHeff_dr      : drvt Heff,r
  *  @param[out] dHeff_dprstar : drvt Heff,prstar
  *  @param[out] dHeff_dpphi   : drvt Heff,pphi
  * 
  */
 void eob_ham(double nu, double r, double pphi, double prstar, 
             double A, double dA,
             double Q, double dQ, double dQ_dprstar,
	           double *H, 
	           double *Heff, 
	           double *dHeff_dr, 
	           double *dHeff_dprstar,
	           double *dHeff_dpphi
	          )
{

  const double pphi2 = SQ(pphi);
  const double u = 1./r;
  const double u2 = SQ(u);
  const double u3 = u2*u;
  const double prstar2 = SQ(prstar);
  const double prstar3 = prstar2*prstar;
  const double prstar4 = prstar2*prstar2;

  *Heff = sqrt(A*(1.0 + pphi2*u2) + prstar2 + Q);
  *H    = sqrt( 1.0 + 2.0*nu*(*Heff - 1) )/nu;  
   
  if (dHeff_dr != NULL)      *dHeff_dr      = 0.5*(A*(-2.*pphi2*u3 + dQ) + dA*(1. + pphi2*u2 + Q))/(*Heff);
  if (dHeff_dprstar != NULL) *dHeff_dprstar = 0.5*(2.*prstar + A*dQ_dprstar)/(*Heff);
  if (dHeff_dpphi != NULL)   *dHeff_dpphi   = A*pphi*u2/(*Heff);
}

/**
 * Function: eob_dyn_rhs_s
 * -----------------------
 *   Computes the r.h.s. of the EOB Hamiltonian dynamics, spins version
 * 
 *  @param[in]  t     : time
 *  @param[in]  y     : EOB variables
 *  @param[out] dy    : r.h.s. of EOB variables
 *  @param[in]  d     : dummy pointer to Dynamics
 * 
 *  @return GSL_SUCCESS
 */
int eob_dyn_rhs_s(double t, const double y[], double dy[], void *d)
{
      
  (void)(t); /* avoid unused parameter warning */
  Dynamics *dyn = d;
  
  const int usetidal = EOBPars->use_tidal;
  const int usespins = EOBPars->use_spins;
  

  /* Shorthands */
  const double r      = y[EOB_EVOLVE_RAD];
  const double prstar = y[EOB_EVOLVE_PRSTAR];
  const double pphi   = y[EOB_EVOLVE_PPHI];
  const double pphi2  = pphi*pphi;
  
  /* Unpack values */
  const double nu    = EOBPars->nu;
  const double S     = EOBPars->S;
  const double Sstar = EOBPars->Sstar;
  const double chi1  = EOBPars->chi1;
  const double chi2  = EOBPars->chi2;
  const double X1    = EOBPars->X1;
  const double X2    = EOBPars->X2;
  const double c3    = EOBPars->cN3LO;
  const double aK2   = EOBPars->aK2;
  const double a1    = EOBPars->a1;
  const double a2    = EOBPars->a2;
  double C_Q1  = EOBPars->C_Q1;
  double C_Q2  = EOBPars->C_Q2;
  double C_Oct1 = EOBPars->C_Oct1;
  double C_Oct2 = EOBPars->C_Oct2;
  double C_Hex1 = EOBPars->C_Hex1;
  double C_Hex2 = EOBPars->C_Hex2;

  if ((EOBPars->use_tidal)&&(EOBPars->use_tidal_fmode_model)) {
    /* Update the dressing factors for the f-mode resonances */
    fmode_resonance_dressing_factors(r, dyn);
    /* Update the QOH */
    fmode_resonance_dress_QOH(dyn);
    C_Q1 = dyn->dressed_C_Q1;
    C_Q2 = dyn->dressed_C_Q2;
    C_Oct1 = dyn->dressed_C_Oct1;
    C_Oct2 = dyn->dressed_C_Oct2;
    C_Hex1 = dyn->dressed_C_Hex1;
    C_Hex2 = dyn->dressed_C_Hex2;
  }
  
  /* Compute Metric */
  double A, B, dA, d2A, dB, Q, dQ, dQ_dprstar, d2Q_dprstar2, pl_hold;
  eob_metric_s(r, prstar, d, &A, &B, &dA, &d2A, &dB, &pl_hold, &Q, &dQ, &dQ_dprstar, &pl_hold, &pl_hold, &d2Q_dprstar2, &pl_hold, &pl_hold, &pl_hold);
  
  /* Failsafe to avoid event horizon crossing */
  if (A < 1e-6){
    A = fabs(A);
    B = fabs(B);
    dyn->ode_stop = true;
  }

  /* Compute centrifugal radius */
  double rc, drc_dr, d2rc_dr2;
  eob_dyn_s_get_rc(r, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, 
       &rc, &drc_dr, &d2rc_dr2);
  
  if ((EOBPars->use_tidal)&&(EOBPars->use_tidal_fmode_model)) {
    /* Add derivative terms to rc' and rc'' from f-mode resonances u-dependent terms */
    eob_dyn_s_rc_add_QOH_drvts(dyn, rc, r, a1, a2,
			       &drc_dr, &d2rc_dr2);
  }
  
  const double uc     = 1./rc;
  const double uc2    = uc*uc;
  const double uc3    = uc2*uc;
  
  /* Compute Hamiltonian */
  double Heff_orb, Heff, H, dHeff_dr, dHeff_dprstar, d2Heff_dprstar20, dHeff_dpphi, d2Heff_dr2;
  eob_ham_s(nu, r, rc, drc_dr, d2rc_dr2, pphi, prstar, S, Sstar, chi1, chi2, X1, X2, aK2, c3, A, dA, d2A, Q, dQ, dQ_dprstar, 0., d2Q_dprstar2,
	    &H, &Heff, &Heff_orb, &dHeff_dr, &dHeff_dprstar, &dHeff_dpphi, &d2Heff_dprstar20, &pl_hold);
  // 0. is d2Q which only enters d2Heff_dr2, not needed here
  
  /* H follows the same convention of Heff, i.e. it is the energy per unit mass,
     while E is the real energy.*/
  double E = nu*H;
  const double ooH = 1./E;
    
  const double sqrtAbyB       = sqrt(A/B);
  const double dp_rstar_dt_0  = - sqrtAbyB*dHeff_dr*ooH;
  const double ddotr_dp_rstar = sqrtAbyB*d2Heff_dprstar20*ooH;
  const double Omg            = dHeff_dpphi*ooH;
  const double ddotr          = dp_rstar_dt_0*ddotr_dp_rstar; /* approximate ddot(r)_0 without Fphi, order pr_star^2 neglected */
  
  /* r evol eqn rhs */
  dy[EOB_EVOLVE_RAD] = sqrtAbyB*dHeff_dprstar*ooH;
  
  /* phi evol eqn rhs */
  dy[EOB_EVOLVE_PHI] = Omg;
  
  /* dp_{r*}/dt */
  dy[EOB_EVOLVE_PRSTAR] = -sqrtAbyB*dHeff_dr*ooH;
    
  /* Compute here the new r_omg radius
     Compute same quantities with prstar=0. This to obtain psi.
     Procedure consistent with the nonspinning case. */
  double ggm0[26];
  eob_dyn_s_GS(r, rc, drc_dr, d2rc_dr2, aK2, 0., pphi, nu, chi1, chi2, X1, X2, c3, ggm0);
  
  const double GS_0       = ggm0[2];
  const double GSs_0      = ggm0[3];
  const double dGS_dr_0   = ggm0[6];
  const double dGSs_dr_0  = ggm0[7];
  const double Heff_orb_0 = sqrt(A*(1.0 + pphi2*uc2));    /* effective Hamiltonian H_0^eff */
  const double Heff_0     = Heff_orb_0 + (GS_0*S + GSs_0*Sstar)*pphi;
  const double H0         = sqrt(1.0 + 2.0*nu*(Heff_0 - 1.0) );
  const double ooH0       = 1./H0;
  const double Gtilde     = GS_0*S     + GSs_0*Sstar;
  const double dGtilde_dr = dGS_dr_0*S + dGSs_dr_0*Sstar;
  const double duc_dr     = -uc2*drc_dr;
  const double psic       = fabs((duc_dr + dGtilde_dr*rc*sqrt(A/pphi2 + A*uc2)/A)/(-0.5*dA));
  //NOTE Different from Matlab code.
  //     Added absolute value to avoid NaN
  const double r_omg      = pow( ((1./sqrt(rc*rc*rc*psic))+Gtilde)*ooH0, -2./3. );
  const double v_phi      = r_omg*Omg;
  const double x          = v_phi*v_phi;
  const double jhat       = pphi/(r_omg*v_phi);

  /* Compute flux and dp_{\phi}/dt */
  if (dyn->noflx) dy[EOB_EVOLVE_PPHI] = 0.;
  else            dy[EOB_EVOLVE_PPHI] = eob_flx_Flux_s(x,Omg,r_omg,E,Heff,jhat,r,prstar,ddotr,dyn);

  if (dyn->store) {
    /* Store values */
    dyn->t = t;
    dyn->r = r;
    dyn->phi = y[EOB_EVOLVE_PHI];
    dyn->pphi = pphi;
    dyn->prstar = prstar;
    dyn->Omg = Omg;
    dyn->Omg_orb = ooH*pphi*A*uc2/Heff_orb;
    dyn->H = H;
    dyn->E = E;
    dyn->Heff = Heff;
    dyn->A = A;
    dyn->dA = dA;
    dyn->d2A = d2A;
    dyn->B = B;
    dyn->dB = dB;
    //    dyn->psi = psi;
    dyn->r_omega = r_omg;
    dyn->v_phi = v_phi;
    dyn->jhat = jhat;
    dyn->ddotr = ddotr;
  }
    
  return GSL_SUCCESS;
}

/**
 * Function: eob_dyn_rhs_ecc
 * -------------------------
 *   Computes the r.h.s. of the EOB Hamiltonian dynamics, eccentric version
 * 
 *  @param[in]  t     : time
 *  @param[in]  y     : EOB variables
 *  @param[out] dy    : r.h.s. of EOB variables
 *  @param[in]  d     : dummy pointer to Dynamics
 * 
 *  @return GSL_SUCCESS
*/
int eob_dyn_rhs_ecc(double t, const double y[], double dy[], void *d)
{
  
  (void)(t); /* avoid unused parameter warning */
  Dynamics *dyn = d;
  
  const int usetidal = EOBPars->use_tidal;
  const int usespins = EOBPars->use_spins;
  
  /* Shorthands */
  const double r      = y[EOB_EVOLVE_RAD];
  const double prstar = y[EOB_EVOLVE_PRSTAR];
  const double pphi   = y[EOB_EVOLVE_PPHI];
  const double pphi2  = pphi*pphi;
  
  /* Unpack values */
  const double nu    = EOBPars->nu;
  const double S     = EOBPars->S;
  const double Sstar = EOBPars->Sstar;
  const double chi1  = EOBPars->chi1;
  const double chi2  = EOBPars->chi2;
  const double X1    = EOBPars->X1;
  const double X2    = EOBPars->X2;
  const double c3    = EOBPars->cN3LO;
  const double aK2   = EOBPars->aK2;
  const double a1    = EOBPars->a1;
  const double a2    = EOBPars->a2;
  double C_Q1  = EOBPars->C_Q1;
  double C_Q2  = EOBPars->C_Q2;
  double C_Oct1 = EOBPars->C_Oct1;
  double C_Oct2 = EOBPars->C_Oct2;
  double C_Hex1 = EOBPars->C_Hex1;
  double C_Hex2 = EOBPars->C_Hex2;

  if ((EOBPars->use_tidal)&&(EOBPars->use_tidal_fmode_model)) {
    /* Update the dressing factors for the f-mode resonances */
    fmode_resonance_dressing_factors(r, dyn);
    /* Update the QOH */
    fmode_resonance_dress_QOH(dyn);
    C_Q1 = dyn->dressed_C_Q1;
    C_Q2 = dyn->dressed_C_Q2;
    C_Oct1 = dyn->dressed_C_Oct1;
    C_Oct2 = dyn->dressed_C_Oct2;
    C_Hex1 = dyn->dressed_C_Hex1;
    C_Hex2 = dyn->dressed_C_Hex2;
  }
  
  /* Compute Metric */
  double A, B, dA, d2A, dB, Q, dQ, dQ_dprstar, d2Q_dprstar2, pl_hold;
  eob_metric_s(r, prstar, d, &A, &B, &dA, &d2A, &dB, &pl_hold, &Q, &dQ, &dQ_dprstar, &pl_hold, &pl_hold, &d2Q_dprstar2, &pl_hold, &pl_hold, &pl_hold);

  /* Failsafe to avoid event horizon crossing */
  if (A < 1e-6){
    A = fabs(A);
    B = fabs(B);
    dyn->ode_stop = true;
  }
  
  /* Compute centrifugal radius */
  double rc, drc_dr, d2rc_dr2;
  eob_dyn_s_get_rc(r, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, 
       &rc, &drc_dr, &d2rc_dr2);
  
  if ((EOBPars->use_tidal)&&(EOBPars->use_tidal_fmode_model)) {
    /* Add derivative terms to rc' and rc'' from f-mode resonances u-dependent terms */
    eob_dyn_s_rc_add_QOH_drvts(dyn, rc, r, a1, a2,
			       &drc_dr, &d2rc_dr2);
  }
  
  const double uc     = 1./rc;
  const double uc2    = uc*uc;
  const double uc3    = uc2*uc;
  
  /* Compute Hamiltonian */
  double Heff_orb, Heff, H, dHeff_dr, dHeff_dprstar, d2Heff_dprstar20, dHeff_dpphi, d2Heff_dr2;
  eob_ham_s(nu, r, rc, drc_dr, d2rc_dr2, pphi, prstar, S, Sstar, chi1, chi2, X1, X2, aK2, c3, A, dA, d2A, Q, dQ, dQ_dprstar, 0., d2Q_dprstar2, 
	    &H, &Heff, &Heff_orb, &dHeff_dr, &dHeff_dprstar, &dHeff_dpphi, &d2Heff_dprstar20, &pl_hold);
  // 0. is d2Q which only enters d2Heff_dr2, not needed here

  /* H follows the same convention of Heff, i.e. it is the energy per unit mass,
     while E is the real energy.*/
  double E = nu*H;
  const double ooH = 1./E;
    
  const double sqrtAbyB       = sqrt(A/B);
  const double dp_rstar_dt_0  = - sqrtAbyB*dHeff_dr*ooH;
  const double ddotr_dp_rstar = sqrtAbyB*d2Heff_dprstar20*ooH;
  const double Omg            = dHeff_dpphi*ooH;
  const double ddotr          = dp_rstar_dt_0*ddotr_dp_rstar; /* approximate ddot(r)_0 without Fphi, order pr_star^2 neglected */
  
  /* r evol eqn rhs */
  dy[EOB_EVOLVE_RAD] = sqrtAbyB*dHeff_dprstar*ooH;
  
  /* phi evol eqn rhs */
  dy[EOB_EVOLVE_PHI] = Omg;
  
  /* Compute here the new r_omg radius */
  double r_omg       = eob_dyn_get_romg(r, prstar, pphi, dyn);
  const double v_phi = r_omg*Omg;
  const double x     = v_phi*v_phi;
  const double jhat  = pphi/(r_omg*v_phi);
  
  /* Compute fluxes and dp_{\phi}/dt */
  double Fphi = 0.;
  double Fr   = 0.;
  eob_flx_Flux_ecc(x, Omg, r_omg, E, Heff, jhat, r, prstar, pphi, dy[EOB_EVOLVE_RAD], ddotr, &Fphi, &Fr, dyn);
  
  if (dyn->noflx) dy[EOB_EVOLVE_PPHI] = 0.;
  else            dy[EOB_EVOLVE_PPHI] = Fphi;
  
  /* dp_{r*}/dt */
  double Frstar = sqrtAbyB*Fr;
  dy[EOB_EVOLVE_PRSTAR] = -sqrtAbyB*dHeff_dr*ooH + Frstar;
  
  if (dyn->store) {
    /* Store values */
    dyn->t = t;
    dyn->r = r;
    dyn->phi = y[EOB_EVOLVE_PHI];
    dyn->pphi = pphi;
    dyn->prstar = prstar;
    dyn->Omg = Omg;
    dyn->Omg_orb = ooH*pphi*A*uc2/Heff_orb;
    dyn->H = H;
    dyn->E = E;
    dyn->Heff = Heff;
    dyn->A = A;
    dyn->dA = dA;
    dyn->d2A = d2A;
    dyn->B = B;
    dyn->dB = dB;
    //    dyn->psi = psi;
    dyn->r_omega = r_omg;
    dyn->v_phi = v_phi;
    dyn->jhat = jhat;
    dyn->ddotr = ddotr;
  }
    
  return GSL_SUCCESS;
}

/**
 * Function: eob_ham_s
 * ------------------
 *   Computes the EOB spinning Hamiltonian
 * 
 *  @param[in]  nu            : symmetric mass ratio
 *  @param[in]  r             : radial separation
 *  @param[in]  rc            : centrifugal radius
 *  @param[in]  drc_dr        : drvt rc,r
 *  @param[in]  d2rc_dr2      : drvt drc,r
 *  @param[in]  pphi          : orbital angular momentum
 *  @param[in]  prstar        : (tortoise) radial momentum
 *  @param[in]  S             : S1 + S2
 *  @param[in]  Sstar         : X2*a1 + X1*a2
 *  @param[in]  chi1          : dimensionless spin of body 1
 *  @param[in]  chi2          : dimensionless spin of body 2
 *  @param[in]  X1            : mass fraction of body 1
 *  @param[in]  X2            : mass fraction of body 2
 *  @param[in]  aK2           : Kerr spin parameter squared (a1 + a2)^2
 *  @param[in]  c3            : CN3LO parameter
 *  @param[in]  A             : A radial potential
 *  @param[in]  dA            : drvt A,r
 *  @param[in]  d2A           : drvt dA,r
 *  @param[in]  Q             : Q potential
 *  @param[in]  dQ            : drvt Q,r
 *  @param[in]  dQ_dprstar    : drvt Q,prstar
 *  @param[in]  d2Q           : 2nd derivative of Q wrt to r
 *  @param[in]  d2Q_dprstar2  : 2nd derivative of Q wrt to pr_star
 *  @param[out] H             : real EOB Hamiltonian divided by mu=m1m2/(m1+m2)
 *  @param[out] Heff          : effective EOB Hamiltonian (divided by mu)
 *  @param[out] Heff_orb      : effective EOB Hamiltonian for the orbital part (divided by mu)
 *  @param[out] dHeff_dr      : drvt Heff,r
 *  @param[out] dHeff_dprstar : drvt Heff,prstar
 *  @param[out] dHeff_dpphi   : drvt Heff,pphi
 *  @param[out] d2Heff_dprstar20 : 2nd derivative of Heff wrt to pr_star neglecting all pr_star^2 terms
 *  @param[out] d2Heff_dr2    : drvt Heff, d2r
 * 
 */
void eob_ham_s(double nu,
               double r,
               double rc,
               double drc_dr,
	             double d2rc_dr2,
               double pphi,
               double prstar,
               double S,
               double Sstar,
               double chi1,
               double chi2,
               double X1,
               double X2,
               double aK2,
               double c3,
               double A,
               double dA,
	             double d2A,
               double Q,
               double dQ,
               double dQ_dprstar,
               double d2Q,
               double d2Q_dprstar2,
               double *H,             /* real EOB Hamiltonian divided by mu=m1m2/(m1+m2) */
               double *Heff,          /* effective EOB Hamiltonian (divided by mu) */
               double *Heff_orb,
               double *dHeff_dr,      /* drvt Heff,r */
               double *dHeff_dprstar, /* drvt Heff,prstar */
               double *dHeff_dpphi,    /* drvt Heff,pphi */
               double *d2Heff_dprstar20,
	       double *d2Heff_dr2
               )
{
    /* Shorthands */
    const double pphi2    = SQ(pphi);
    const double prstar2 = SQ(prstar);
    const double prstar3 = prstar2*prstar;
    const double prstar4 = prstar2*prstar2;
    const double uc  = 1./rc;
    const double uc2 = uc*uc;
    const double uc3 = uc2*uc;
    const double uc4 = uc3*uc;
    
    /* Compute spin-related functions*/
    double ggm[26];
    eob_dyn_s_GS(r, rc, drc_dr, d2rc_dr2, aK2, prstar, pphi, nu, chi1, chi2, X1, X2, c3, ggm);
    const double GS              = ggm[2];
    const double GSs             = ggm[3];
    const double dGS_dprstar     = ggm[4];
    const double dGSs_dprstar    = ggm[5];
    const double dGS_dr          = ggm[6];
    const double dGSs_dr         = ggm[7];
    const double dGSs_dpphi      = ggm[9];
    const double d2GS_dprstar20  = ggm[12];
    const double d2GSs_dprstar20 = ggm[13];
    const double d2GS_dr2        = ggm[14];
    const double d2GSs_dr2       = ggm[15];
    
    /* Compute Hamiltonian and its derivatives */
    *Heff_orb         = sqrt( prstar2+A*(1. + pphi2*uc2 +  Q) );
    *Heff             = *Heff_orb + (GS*S + GSs*Sstar)*pphi;
    *H                = sqrt( 1. + 2.*nu*(*Heff - 1.) )/nu;
    
    double ooHeff_orb          = 1./(*Heff_orb);
    double dHefforb_dr         = 0.5*ooHeff_orb*(A*(-2.*pphi2*uc3*drc_dr + dQ) + dA*(1. + pphi2*uc2 + Q));
    double dHefforb_dprstar    = 0.5*ooHeff_orb*(2.*prstar + A*dQ_dprstar);
    double d2Hefforb_dprstar20 = 0.5*ooHeff_orb*( 2. + A*d2Q_dprstar2 );
    double d2Hefforb_dr2       = ooHeff_orb*(-SQ(dHefforb_dr) + A*(3.*pphi2*uc4*SQ(drc_dr) - pphi2*uc3*d2rc_dr2 + 0.5*d2Q) + 
                                 dA*(-2.*pphi2*uc3*drc_dr + dQ) + 0.5*d2A*(1. + pphi2*uc2 + Q));

    if (dHeff_dr != NULL)         *dHeff_dr         = dHefforb_dr + pphi*(dGS_dr*S + dGSs_dr*Sstar);
    if (dHeff_dprstar != NULL)    *dHeff_dprstar    = dHefforb_dprstar + pphi*(dGS_dprstar*S + dGSs_dprstar*Sstar);
    if (d2Heff_dprstar20 != NULL) *d2Heff_dprstar20 = d2Hefforb_dprstar20 + pphi*(d2GS_dprstar20*S + d2GSs_dprstar20*Sstar); /* second derivative of Heff wrt to pr_star neglecting all pr_star^2 terms */
    if (dHeff_dpphi != NULL)      *dHeff_dpphi      = GS*S + (GSs + pphi*dGSs_dpphi)*Sstar + pphi*A*uc2*ooHeff_orb;
    if (d2Heff_dr2 != NULL)       *d2Heff_dr2       = d2Hefforb_dr2 + pphi*(d2GS_dr2*S + d2GSs_dr2*Sstar);

}


/**
 * Function: eob_dyn_s_GS
 * ----------------------
 *   Computes the gyro-gravitomagnetic functions GS and GS*, that are called GS and GSs.
 * 
 *  @param[in]  r             : radial separation
 *  @param[in]  rc            : centrifugal radius
 *  @param[in]  drc_dr        : drvt rc,r
 *  @param[in]  d2rc_dr2      : drvt drc,r
 *  @param[in]  aK2           : squared Kerr paramete
 *  @param[in]  prstar        : (tortoise) radial momentum
 *  @param[in]  pphi          : orbital angular momentum
 *  @param[in]  nu            : symmetric mass ratio
 *  @param[in]  chi1          : dimensionless spin of body 1
 *  @param[in]  chi2          : dimensionless spin of body 2
 *  @param[in]  X1            : mass fraction of body 1
 *  @param[in]  X2            : mass fraction of body 2
 *  @param[in]  cN3LO         : CN3LO parameter
 *  @param[out] ggm           : gyrp-gravitomagnetic functions, output structure
 *
 */
void eob_dyn_s_GS(double r, double rc, double drc_dr, double d2rc_dr2, double aK2, double prstar, double pph, double nu, double chi1, double chi2, double X1, double X2, double cN3LO,
	  double *ggm)
{
  
  static double c10,c20,c30,c02,c12,c04;
  static double cs10,cs20,cs30,cs40,cs02,cs12,cs04;

  /* Compute the nu-dep. coefficient at first call only */
  //static int firstcall = 1;  
  if (EOBPars->firstcall[FIRSTCALL_EOBDYNSGS]) {
    EOBPars->firstcall[FIRSTCALL_EOBDYNSGS] = 0;      
    double nu2   = nu*nu;
    /* coefficients of hat{GS} */
    c10 =  5./16.*nu;
    c20 =  51./8.*nu + 41./256.*nu2;
    c30 =  nu*cN3LO;
    c02 =  27./16.*nu;
    c12 =  12.*nu - 49./128.*nu2;
    c04 = -5./16.*nu + 169./256.*nu2;
    /* coefficients of hat{GS*} */
    cs10 = 3./4.   + nu/2.;
    cs20 = 27./16. + 29./4.*nu + 3./8.*nu2;
    cs02 = 5./4.   + 3./2.*nu;
    cs12 = 4.   + 11.*nu     - 7./8.*nu2;
    cs04 = 5./48.  + 25./12.*nu + 3./8.*nu2;
    cs30 = nu*cN3LO + 135./32.;
    cs40 = 2835./256.;
  }

  double u  = 1./r;
  double u2 = u*u;
  double u3 = u2*u;
  
  double uc      = 1./rc;
  double uc2     = uc*uc;
  double uc3     = uc2*uc;
  double uc4     = uc3*uc;
  double prstar2 = prstar*prstar;
  double prstar3 = prstar2*prstar;
  double prstar4 = prstar2*prstar2;
  
  double GS0        = 2.*u*uc2;
  double dGS0_duc   = 2.*u2/drc_dr + 4.*u*uc;
  double d2GS0_duc2 = 4.*u3/(uc2*SQ(drc_dr))
    + 2.*u2*d2rc_dr2/(uc2*SQ(drc_dr)*drc_dr) + 4.*u2/(uc*drc_dr) + 4.*u;
  
  double GSs0          =  3./2.*uc3;
  double dGSs0_duc     =  9./2.*uc2;
  double d2GSs0_duc2   =  9.*uc;
  double dGSs0_dprstar =  0.0;
  double dGSs0_dpph    =  0.0;
  
  double hGS  = 1./(1.  + c10*uc + c20*uc2 + c30*uc3 + c02*prstar2 + c12*uc*prstar2 + c04*prstar4);   
  double hGSs = 1./(1.  + cs10*uc + cs20*uc2  + cs30*uc3 + cs40*uc4 + cs02*prstar2 + cs12*uc*prstar2 + cs04*prstar4); 
  
  /* complete gyro-gravitomagnetic functions */
  double GS  =  GS0*hGS; 
  double GSs = GSs0*hGSs; 
  
  /* Get derivatives of gyro-gravitomagnetic functions */
  double dhGS_dprstar  = -2.*prstar*hGS*hGS *( c02 +  c12*uc +  2.*c04*prstar2);
  double dhGSs_dprstar = -2.*prstar*hGSs*hGSs*(cs02 + cs12*uc + 2.*cs04*prstar2);
  
  double dGS_dprstar  = GS0 *dhGS_dprstar; 
  double dGSs_dprstar = GSs0*dhGSs_dprstar + dGSs0_dprstar*hGSs; 
  
  /* derivatives of hat{G} with respect to uc */
  double dhGS_duc  = -hGS*hGS*(c10 + 2.*c20*uc  + 3.*c30*uc2 + c12*prstar2);
  double dhGSs_duc = -hGSs*hGSs*(cs10 + 2.*cs20*uc + 3.*cs30*uc2 + 4.*cs40*uc3 + cs12*prstar2);

  /* second derivative of  hat{G} with respect to uc^2 */
  double d2hGS_duc2  = hGS*hGS*(2.*hGS*SQ(c10 + c12*prstar2 + 2.*c20*uc + 3.*c30*uc2)
			       - (2.*c20 + 6.*c30*uc));
  double d2hGSs_duc2 = hGSs*hGSs*(2.*hGSs*SQ(cs10 + cs12*prstar2 + 2.*cs20*uc + 3.*cs30*uc2 + 4.*cs40*uc3)
			       - (2.*cs20 + 6.*cs30*uc + 12.*cs40*uc2));

  /* second derivative of  hat{G} with respect to prstar^2 */
  double d2hGS_dprstar2  = hGS*hGS*(2.*hGS*SQ(2.*c02*prstar + 4.*c04*prstar3 + 2.*c12*prstar*uc)
				   - (2.*c02 + 12.*c04*prstar2 + 2.*c12*uc));
  double d2hGSs_dprstar2 = hGSs*hGSs*(2.*hGSs*SQ(2.*cs02*prstar + 4.*cs04*prstar3 + 2.*cs12*prstar*uc)
				   - (2.*cs02 + 12.*cs04*prstar2 + 2.*cs12*uc));
  
  /* second derivative of  hat{G} with respect to uc/prstar */
  double d2hGS_duc_dprstar  = hGS*hGS*(-2.*c12*prstar
			      + 2.*hGS*(2.*c02*prstar + 4.*c04*prstar3 + 2.*c12*prstar*uc)
			      *(c10 + c12*prstar2 + 2.*c20*uc + 3.*c30*uc2));
  double d2hGSs_duc_dprstar = hGSs*hGSs*(-2.*cs12*prstar
			      + 2.*hGSs*(2.*cs02*prstar + 4.*cs04*prstar3 + 2.*cs12*prstar*uc)
			      *(cs10 + cs12*prstar2 + 2.*cs20*uc + 3.*cs30*uc2 + 4.*cs40*uc3));

  /* third derivative of  hat{G} with respect to prstar^3 */
  double d3hGS_dprstar3  = hGS*hGS*(-24.*c04*prstar + 6.*hGS*(2.*c02 + 12.*c04*prstar2 + 2.*c12*uc)
				    *(2.*c02*prstar + 4.*c04*prstar3 + 2.*c12*prstar*uc)
                                      - 6.*hGS*hGS*SQ(2.*c02*prstar + 4.*c04*prstar3 + 2.*c12*prstar*uc)
				    *(2.*c02*prstar + 4.*c04*prstar3 + 2.*c12*prstar*uc));
  double d3hGSs_dprstar3 = hGSs*hGSs*(-24.*cs04*prstar + 6.*hGSs*(2.*cs02 + 12.*cs04*prstar2 + 2.*cs12*uc)
				    *(2.*cs02*prstar + 4.*cs04*prstar3 + 2.*cs12*prstar*uc)
                                      - 6.*hGSs*hGSs*SQ(2.*cs02*prstar + 4.*cs04*prstar3 + 2.*cs12*prstar*uc)
				    *(2.*cs02*prstar + 4.*cs04*prstar3 + 2.*cs12*prstar*uc));

  /* third derivative of  hat{G} with respect to prstar^2/uc */
  double d3hGS_duc_dprstar2  = hGS*hGS*(-2.*c12 + 2.*hGS*(2.*c02 + 12.*c04*prstar2 + 2.*c12*uc)
					*(c10 + c12*prstar2 + 2.*c20*uc + 3.*c30*uc2)
					+ hGS*8.*c12*prstar*(2.*c02*prstar + 4.*c04*prstar3 + 2.*c12*prstar*uc)
					- hGS*hGS*6.*SQ(2.*c02*prstar + 4.*c04*prstar3 + 2.*c12*prstar*uc)
					*(c10 + c12*prstar2 + 2.*c20*uc + 3.*c30*uc2));
  double d3hGSs_duc_dprstar2 = hGSs*hGSs*(-2.*cs12 + 2.*hGSs*(2.*cs02 + 12.*cs04*prstar2 + 2.*cs12*uc)
					*(cs10 + cs12*prstar2 + 2.*cs20*uc + 3.*cs30*uc2 + 4.*cs40*uc3)
					+ hGSs*8.*cs12*prstar*(2.*cs02*prstar + 4.*cs04*prstar3 + 2.*cs12*prstar*uc)
					- hGSs*hGSs*6.*SQ(2.*cs02*prstar + 4.*cs04*prstar3 + 2.*cs12*prstar*uc)
					*(cs10 + cs12*prstar2 + 2.*cs20*uc + 3.*cs30*uc2 + 4.*cs40*uc3));

  /* third derivative of  hat{G} with respect to prstar/uc^2 */
  double d3hGS_duc2_dprstar  = SQ(hGS)*hGS*(-6.*hGS*(2*c02*prstar + 4*c04*prstar3 + 2*c12*prstar*uc)
					   *SQ(c10 + c12*prstar2 + 2.*c20*uc + 3.*c30*uc2)
					   + 2.*(2.*c02*prstar + 4.*c04*prstar3 + 2.*c12*prstar*uc)*(2.*c20 + 6.*c30*uc)
					   + 8.*c12*prstar*(c10 + c12*prstar2 + 2.*c20*uc + 3.*c30*uc2));
  double d3hGSs_duc2_dprstar = SQ(hGSs)*hGSs*(-6.*hGSs*(2*cs02*prstar + 4*cs04*prstar3 + 2*cs12*prstar*uc)
					   *SQ(cs10 + cs12*prstar2 + 2.*cs20*uc + 3.*cs30*uc2 + 4.*cs40*uc3)
					   + 2.*(2.*cs02*prstar + 4.*cs04*prstar3 + 2.*cs12*prstar*uc)*(2.*cs20 + 6.*cs30*uc + 12.*cs40*uc2)
					   + 8.*cs12*prstar*(cs10 + cs12*prstar2 + 2.*cs20*uc + 3.*cs30*uc2 + 4.*cs40*uc3));
  
  /* derivatives of G with respect to uc */
  double dGS_duc  =  dGS0_duc*hGS  +  GS0*dhGS_duc;
  double dGSs_duc = dGSs0_duc*hGSs + GSs0*dhGSs_duc;
  
  /* derivatives of (G,G*) with respect to r */
  double dGS_dr  = -drc_dr*uc2*dGS_duc; 
  double dGSs_dr = -drc_dr*uc2*dGSs_duc; 
    
  /* derivatives of (G,G*) with respect to pph */
  double dGS_dpph  = 0.; 
  double dGSs_dpph = dGSs0_dpph*hGSs;    
  
  /* For initial data: compute the two ratios of ggm.dG_dprstar/prstar for GS and GSs */
  const double dGS_dprstarbyprstar  = -2.*GS0*hGS*hGS *( c02  +  c12*uc +  2.*c04*prstar2);
  const double dGSs_dprstarbyprstar = -2.*GSs0*hGSs*hGSs*(cs02 + cs12*uc + 2.*cs04*prstar2);
  
  /* For NQC: Second derivatives neglecting all pr_star^2 terms */
  const double d2GS_dprstar20  =  GS0*(-2.*hGS*hGS *( c02 +  c12*uc +  2.*c04*prstar2));
  const double d2GSs_dprstar20 =  GSs0*(-2.*hGSs*hGSs*(cs02 + cs12*uc + 2.*cs04*prstar2));

  /* Second derivatives with respect to uc^2 */
  double d2GS_duc2  = d2GS0_duc2*hGS + 2.*dGS0_duc*dhGS_duc + GS0*d2hGS_duc2;
  double d2GSs_duc2 = d2GSs0_duc2*hGSs + 2.*dGSs0_duc*dhGSs_duc + GSs0*d2hGSs_duc2;

  /* Second derivatives with respect to r^2 */
  double d2GS_dr2  =  2*uc3*SQ(drc_dr)*dGS_duc - uc2*d2rc_dr2*dGS_duc + uc4*SQ(drc_dr)*d2GS_duc2;
  double d2GSs_dr2 =  2*uc3*SQ(drc_dr)*dGSs_duc - uc2*d2rc_dr2*dGSs_duc + uc4*SQ(drc_dr)*d2GSs_duc2;
	    
  /* Second derivatives with respect to prstar^2 */
  double d2GS_dprstar2  = GS0*d2hGS_dprstar2;
  double d2GSs_dprstar2 = GSs0*d2hGSs_dprstar2;

  /* Second derivatives with respect to prstar/uc */
  double d2GS_duc_dprstar  = dGS0_duc*dhGS_dprstar + GS0*d2hGS_duc_dprstar;
  double d2GSs_duc_dprstar = dGSs0_duc*dhGSs_dprstar + GSs0*d2hGSs_duc_dprstar;

  /* Second derivatives with respect to prstar/r */
  double d2GS_dr_dprstar  = -uc2*drc_dr*d2GS_duc_dprstar;  
  double d2GSs_dr_dprstar = -uc2*drc_dr*d2GSs_duc_dprstar;  
  
  /* Third derivatives with respect to prstar^3 */
  double d3GS_dprstar3  = GS0*d3hGS_dprstar3;
  double d3GSs_dprstar3 = GSs0*d3hGSs_dprstar3;
  
  /* Third derivatives with respect to prstar^2/uc */
  double d3GS_duc_dprstar2  = dGS0_duc*d2hGS_dprstar2 + GS0*d3hGS_duc_dprstar2;
  double d3GSs_duc_dprstar2  = dGSs0_duc*d2hGSs_dprstar2 + GSs0*d3hGSs_duc_dprstar2;

  /* Third derivatives with respect to prstar^2/r */
  double d3GS_dr_dprstar2  = -uc2*drc_dr*d3GS_duc_dprstar2;
  double d3GSs_dr_dprstar2 = -uc2*drc_dr*d3GSs_duc_dprstar2;

  /* Third derivatives with respect to prstar/uc^2 */
  double d3GS_duc2_dprstar  = d2GS0_duc2*dhGS_dprstar + 2.*dGS0_duc*d2hGS_duc_dprstar + GS0*d3hGS_duc2_dprstar;
  double d3GSs_duc2_dprstar = d2GSs0_duc2*dhGSs_dprstar + 2.*dGSs0_duc*d2hGSs_duc_dprstar + GSs0*d3hGSs_duc2_dprstar;

  /* Third derivatives with respect to prstar/r^2 */
  double d3GS_dr2_dprstar  = 2*uc3*SQ(drc_dr)*d2GS_duc_dprstar - uc2*d2rc_dr2*d2GS_duc_dprstar
    + uc4*SQ(drc_dr)*d3GS_duc2_dprstar;
  double d3GSs_dr2_dprstar = 2*uc3*SQ(drc_dr)*d2GSs_duc_dprstar - uc2*d2rc_dr2*d2GSs_duc_dprstar
    + uc4*SQ(drc_dr)*d3GSs_duc2_dprstar;
  
  ggm[0]=hGS;
  ggm[1]=hGSs;
  ggm[2]=GS;
  ggm[3]=GSs;
  ggm[4]=dGS_dprstar;
  ggm[5]=dGSs_dprstar;
  ggm[6]=dGS_dr;
  ggm[7]=dGSs_dr;
  ggm[8]=dGS_dpph;
  ggm[9]=dGSs_dpph;
  ggm[10]=dGS_dprstarbyprstar;
  ggm[11]=dGSs_dprstarbyprstar;
  ggm[12]=d2GS_dprstar20;
  ggm[13]=d2GSs_dprstar20;
  ggm[14]=d2GS_dr2;
  ggm[15]=d2GSs_dr2;
  ggm[16]=d2GS_dprstar2;
  ggm[17]=d2GSs_dprstar2;
  ggm[18]=d2GS_dr_dprstar;
  ggm[19]=d2GSs_dr_dprstar;
  ggm[20]=d3GS_dprstar3;
  ggm[21]=d3GSs_dprstar3;
  ggm[22]=d3GS_dr2_dprstar;
  ggm[23]=d3GSs_dr2_dprstar;
  ggm[24]=d3GS_dr_dprstar2;
  ggm[25]=d3GSs_dr_dprstar2;
    
}


/**
 * Function: eob_dyn_s_get_rc_LO
 * -----------------------------
 *   Computes the radius rc that includes of LO spin-square coupling.
 *   The S1*S2 term coincides with the BBH one, no effect of structure.
 *   The self-spin couplings, S1*S1 and S2*S2 get a EOS-dependent coefficient, CQ, that describe the quadrupole
 *   deformation due to spin. Notation of Levi-Steinhoff, JCAP 1412 (2014), no.12, 003. Notation analogous to
 *   the parameter a of Poisson, PRD 57, (1998) 5287-5290 or C_ES^2 in Porto & Rothstein, PRD 78 (2008), 044013
 *    
 *  @param[in]  r             : radial separation
 *  @param[in]  nu            : symmetric mass ratio
 *  @param[in]  at1           : spin of body 1
 *  @param[in]  at2           : spin of body 2
 *  @param[in]  aK2           : kerr parameter squared
 *  @param[in]  C_Q1          : spin-induced quadrupole of body 1
 *  @param[in]  C_Q2          : spin-induced quadrupole of body 2
 *  @param[in]  C_Oct1        : spin-induced octupole of body 1
 *  @param[in]  C_Oct2        : spin-induced octupole of body 2
 *  @param[in]  C_Hex1        : spin-induced hexadecapole of body 1
 *  @param[in]  C_Hex2        : soin-induced hexadecapole of body 2
 *  @param[in]  usetidal      : tidal flag
 *  @param[out] rc            : centrifugal radius
 *  @param[out] drc_dr        : drvt rc,r
 *  @param[out] d2rc_dr2      : drvt drc,r
 * 
 */
void eob_dyn_s_get_rc_LO(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, double C_Oct1, double C_Oct2, double C_Hex1, double C_Hex2, int usetidal, 
		      double *rc, double *drc_dr, double *d2rc_dr2)
{

  double u   = 1./r;
  double u2  = u*u;
  double u3  = u*u2;
  double r2  = r*r;
    
  if (usetidal) {    
#if (EXCLUDESPINSPINTIDES)
    /* Switch off spin-spin-tidal couplings */
    /* See also: eob_wav_flm_s() */
    double rc2 = r2;
    *rc = r;
    *drc_dr = 1.;
    *d2rc_dr2 = 0.;
    /* Above code switch off everything, 
       Alt. one can set C_Q1=C_Q2=0, but keep centrifugal radius */
    /* 
       double a02  = 2.*at1*at2; 
       double rc2  = r2 + a02*(1.+2.*u); 
       *rc         = sqrt(rc2); 
       *drc_dr     = r/(*rc)*(1.-a02*u3); 
       *d2rc_dr2   = 1./(*rc)*(1.-(*drc_dr)*r/(*rc)*(1.-a02*u3)+2.*a02*u3); 
       */
#else
    /* BNS effective spin parameter */
    double a02  = C_Q1*at1*at1 + 2.*at1*at2 + C_Q2*at2*at2;
    double rc2  = r2 + a02*(1.+2.*u); /* tidally-modified centrifugal radius */
    *rc         = sqrt(rc2);
    *drc_dr     = r/(*rc)*(1.-a02*u3);
    *d2rc_dr2   = 1./(*rc)*(1.-(*drc_dr)*r/(*rc)*(1.-a02*u3)+2.*a02*u3);
#endif
  } else {
    /*
    double X12      = sqrt(1.-4.*nu);   
    double alphanu2 = 1. + 0.5/aK2*(- at2*at2*(5./4. + 5./4.*X12 + nu/2.) - at1*at1*(5./4. - 5./4.*X12 +nu/2.) + at1*at2*(-2.+nu));
    double rc2 = r2 + aK2*(1. + 2.*alphanu2/r);
    *rc         = sqrt(rc2);
    *drc_dr     = r/(*rc)*(1.+aK2*(-alphanu2*u3 ));
    *d2rc_dr2   = 1./(*rc)*(1.-(*drc_dr)*r/(*rc)*(1.-alphanu2*aK2*u3)+ 2.*alphanu2*aK2*u3);
    */
    /* Following implementation is regular (avoids 1/aK2) */
    double X12 = sqrt(1.-4.*nu);   
    double c_ss_nlo = (- at2*at2*(1.25 + 1.25*X12 + 0.5*nu) - at1*at1*(1.25 - 1.25*X12 + 0.5*nu) + at1*at2*(-2.+nu));
    double rc2   = r2 + aK2*(1. + 2.*u) + u*c_ss_nlo;
    *rc          = sqrt(rc2);
    double divrc = 1.0/(*rc);
    *drc_dr      = r*divrc*(1-(aK2 + 0.5*c_ss_nlo)*u3);	
    *d2rc_dr2    = divrc*(1.-(*drc_dr)*r*divrc*(1.-(aK2+0.5*c_ss_nlo)*u3)+ (2.*aK2 + c_ss_nlo)*u3);
  }
  
}

/**
  * Function: eob_dyn_s_get_rc_NLO
  * ------------------------------
  *   Computes the radius rc that includes of NLO spin-square coupling.
  *   The S1*S2 term coincides with the BBH one, no effect of structure.
  *   The self-spin couplings, S1*S1 and S2*S2 get a EOS-dependent coefficient, CQ, that describe the quadrupole
  *   deformation due to spin. Notation of Levi-Steinhoff, JCAP 1412 (2014), no.12, 003. Notation analogous to
  *   the parameter a of Poisson, PRD 57, (1998) 5287-5290 or C_ES^2 in Porto & Rothstein, PRD 78 (2008), 044013
  *    
  *  @param[in]  r             : radial separation
  *  @param[in]  nu            : symmetric mass ratio
  *  @param[in]  at1           : spin of body 1
  *  @param[in]  at2           : spin of body 2
  *  @param[in]  aK2           : kerr parameter squared
  *  @param[in]  C_Q1          : spin-induced quadrupole of body 1
  *  @param[in]  C_Q2          : spin-induced quadrupole of body 2
  *  @param[in]  C_Oct1        : spin-induced octupole of body 1
  *  @param[in]  C_Oct2        : spin-induced octupole of body 2
  *  @param[in]  C_Hex1        : spin-induced hexadecapole of body 1
  *  @param[in]  C_Hex2        : soin-induced hexadecapole of body 2
  *  @param[in]  usetidal      : tidal flag
  *  @param[out] rc            : centrifugal radius
  *  @param[out] drc_dr        : drvt rc,r
  *  @param[out] d2rc_dr2      : drvt drc,r
  *
*/
void eob_dyn_s_get_rc_NLO(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, double C_Oct1, double C_Oct2, double C_Hex1, double C_Hex2, int usetidal, 
		      double *rc, double *drc_dr, double *d2rc_dr2)
{

  double u   = 1./r;
  double u2  = u*u;
  double u3  = u*u2;
  double r2  = r*r;
  double X12 = sqrt(1.-4.*nu);   
    
  if (usetidal) {

    /* BNS effective spin parameter */
    double a02      = C_Q1*at1*at1 + 2.*at1*at2 + C_Q2*at2*at2;

    double delta_a2 = X12*(at1*at1*(C_Q1+0.25) - at2*at2*(C_Q2+0.25))
      + at1*at1*(-17./4.+3.*C_Q1-0.5*nu)
      + at2*at2*(-17./4.+3.*C_Q2-0.5*nu)
      + at1*at2*(nu-2.0);

    double rc2 = r2 + a02*(1. + 2.*u) + delta_a2*u;
    *rc         = sqrt(rc2);
    double divrc = 1.0/(*rc);
    *drc_dr     = divrc*(r - (a02 + 0.5*delta_a2)*u2);
    *d2rc_dr2   = divrc*(1 + (2.*a02 + delta_a2)*u3 - (*drc_dr)*(*drc_dr));

  } else {
    
    double c_ss_nlo = (- at2*at2*(1.25 + 1.25*X12 + 0.5*nu) - at1*at1*(1.25 - 1.25*X12 + 0.5*nu) + at1*at2*(-2.+nu));
    double rc2   = r2 + aK2*(1. + 2.*u) + u*c_ss_nlo;
    *rc          = sqrt(rc2);
    double divrc = 1.0/(*rc);
    *drc_dr      = r*divrc*(1-(aK2 + 0.5*c_ss_nlo)*u3);	
    *d2rc_dr2    = divrc*(1.-(*drc_dr)*r*divrc*(1.-(aK2+0.5*c_ss_nlo)*u3)+ (2.*aK2 + c_ss_nlo)*u3);
    
  }
  
}

/**
 * Function: eob_dyn_s_get_rc_NNLO
 * --------------------------------
 *   Computes the radius rc that includes of NNLO spin-square coupling.
 *   The S1*S2 term coincides with the BBH one, no effect of structure.
 *   The self-spin couplings, S1*S1 and S2*S2 get a EOS-dependent coefficient, CQ, that describe the quadrupole
 *   deformation due to spin. Notation of Levi-Steinhoff, JCAP 1412 (2014), no.12, 003. Notation analogous to
 *   the parameter a of Poisson, PRD 57, (1998) 5287-5290 or C_ES^2 in Porto & Rothstein, PRD 78 (2008), 044013
 *    
 *  @param[in]  r             : radial separation
 *  @param[in]  nu            : symmetric mass ratio
 *  @param[in]  at1           : spin of body 1
 *  @param[in]  at2           : spin of body 2
 *  @param[in]  aK2           : kerr parameter squared
 *  @param[in]  C_Q1          : spin-induced quadrupole of body 1
 *  @param[in]  C_Q2          : spin-induced quadrupole of body 2
 *  @param[in]  C_Oct1        : spin-induced octupole of body 1
 *  @param[in]  C_Oct2        : spin-induced octupole of body 2
 *  @param[in]  C_Hex1        : spin-induced hexadecapole of body 1
 *  @param[in]  C_Hex2        : soin-induced hexadecapole of body 2
 *  @param[in]  usetidal      : tidal flag
 *  @param[out] rc            : centrifugal radius
 *  @param[out] drc_dr        : drvt rc,r
 *  @param[out] d2rc_dr2      : drvt drc,r
 * 
 */
void eob_dyn_s_get_rc_NNLO(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, double C_Oct1, double C_Oct2, double C_Hex1, double C_Hex2, int usetidal, 
			   double *rc, double *drc_dr, double *d2rc_dr2)
{

  double u   = 1./r;
  double u2  = u*u;
  double u3  = u*u2;
  double u4  = u*u3;
  double u5  = u*u4;
  double r2  = r*r;
  double X12 = sqrt(1.-4.*nu);   
    
  if (usetidal) {

    /* BNS effective spin parameter */
    double a02      = C_Q1*at1*at1 + 2.*at1*at2 + C_Q2*at2*at2;
    
    double delta_a2 = X12*(at1*at1*(C_Q1+0.25) - at2*at2*(C_Q2+0.25))
      + at1*at1*(-17./4.+3.*C_Q1-0.5*nu)
      + at2*at2*(-17./4.+3.*C_Q2-0.5*nu)
      + at1*at2*(nu-2.0);

    double delta_a2_nnlo  =
        (  387./28.  - 207./28.*nu              )     *a02
      + (-2171./112. - 269./28.*nu + 0.375*nu*nu)     *(at1*at1+at2*at2)
      + (- 281./7    - 187./56.*nu - 0.75 *nu*nu)     *at1*at2
      +    163./28.                               *X12*(C_Q1*at1*at1-C_Q2*at2*at2)
      + (  -29./112. - 2.625   *nu              ) *X12*(at1*at1-at2*at2);
    
    double alphanu2 = 1. + 0.5/a02*delta_a2;
        
    double rc2   =  r2 + a02*(1. + 2.*u) + delta_a2*u + delta_a2_nnlo*u2;
    *rc          = sqrt(rc2);
    double divrc = 1.0/(*rc);
    *drc_dr      = divrc*(r - (a02 + 0.5*delta_a2)*u2 - delta_a2_nnlo*u3);
    *d2rc_dr2    = divrc*(1 + (2.*a02 + delta_a2)*u3
			  + 3*delta_a2_nnlo*u4 - (*drc_dr)*(*drc_dr));

  } else {

    double a0  = at1 + at2;
    double a12 = at1 - at2;
    
    double c_ss_nlo = -1.125*a0*a0 -(0.125+0.5+nu)*a12*a12 + 1.25*X12*a0*a12;

    double c_ss_nnlo = - (189./32. + 417./32.*nu              )    *a0 *a0
                       + ( 11./32. - 127./32.*nu + 0.375*nu*nu)    *a12*a12
                       + ( 89./16. -  2.625*nu              )*X12*a0 *a12;

    
    double rc2   = r2 + aK2*(1. + 2.*u) + u*c_ss_nlo + u2*c_ss_nnlo;
    *rc          = sqrt(rc2);
    double divrc = 1.0/(*rc);
    *drc_dr      = r*divrc*(1-(aK2 + 0.5*c_ss_nlo)*u3 - 0.5*u4*c_ss_nnlo);	
    *d2rc_dr2    = 1./r*(*drc_dr) + r*divrc*((3.*aK2+c_ss_nlo)*u4 + 2.*c_ss_nnlo*u5);
    
  }
  
}

/**
 * Function: eob_dyn_s_get_rc_NNLO_S4
 * -----------------------------------
 *   Computes the radius rc that includes NNLO spin-square coupling,
 *   with the addition of the LO spin^4 coefficient that depends on C_Q, C_Oct and C_Hex.
 *   The S1*S2 term coincides with the BBH one, no effect of structure.
 *   The self-spin couplings, S1*S1 and S2*S2 get a EOS-dependent coefficient, CQ, that describe the quadrupole
 *   deformation due to spin. Notation of Levi-Steinhoff, JCAP 1412 (2014), no.12, 003. Notation analogous to
 *   the parameter a of Poisson, PRD 57, (1998) 5287-5290 or C_ES^2 in Porto & Rothstein, PRD 78 (2008), 044013
 *    
 *  @param[in]  r             : radial separation
 *  @param[in]  nu            : symmetric mass ratio
 *  @param[in]  at1           : spin of body 1, chi1* X1
 *  @param[in]  at2           : spin of body 2, chi2* X2
 *  @param[in]  aK2           : kerr parameter squared
 *  @param[in]  C_Q1          : spin-induced quadrupole of body 1
 *  @param[in]  C_Q2          : spin-induced quadrupole of body 2
 *  @param[in]  C_Oct1        : spin-induced octupole of body 1
 *  @param[in]  C_Oct2        : spin-induced octupole of body 2
 *  @param[in]  C_Hex1        : spin-induced hexadecapole of body 1
 *  @param[in]  C_Hex2        : soin-induced hexadecapole of body 2
 *  @param[in]  usetidal      : tidal flag
 *  @param[out] rc            : centrifugal radius
 *  @param[out] drc_dr        : drvt rc,r
 *  @param[out] d2rc_dr2      : drvt drc,r
 * 
 */
void eob_dyn_s_get_rc_NNLO_S4(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, double C_Oct1, double C_Oct2, double C_Hex1, double C_Hex2, int usetidal, double *rc, double *drc_dr, double *d2rc_dr2)
{

  double u   = 1./r;
  double u2  = u*u;
  double u3  = u*u2;
  double u4  = u*u3;
  double u5  = u*u4;
  double r2  = r*r;
  double X12 = sqrt(1.-4.*nu);   
    
  if (usetidal) {

    /* BNS effective spin parameter */
    double a02      = C_Q1*at1*at1 + 2.*at1*at2 + C_Q2*at2*at2;
    
    double delta_a2 = X12*(at1*at1*(C_Q1+0.25) - at2*at2*(C_Q2+0.25))
      + at1*at1*(-17./4.+3.*C_Q1-0.5*nu)
      + at2*at2*(-17./4.+3.*C_Q2-0.5*nu)
      + at1*at2*(nu-2.0);

    double delta_a2_nnlo  =
        (  387./28.  - 207./28.*nu              )     *a02
      + (-2171./112. - 269./28.*nu + 0.375*nu*nu)     *(at1*at1+at2*at2)
      + (- 281./7    - 187./56.*nu - 0.75 *nu*nu)     *at1*at2
      +    163./28.                               *X12*(C_Q1*at1*at1-C_Q2*at2*at2)
      + (  -29./112. - 2.625   *nu              ) *X12*(at1*at1-at2*at2);

    double delta_a4_lo = 0.75*(C_Hex1 - C_Q1*C_Q1)*at1*at1*at1*at1
                         + 3.*(C_Oct1 - C_Q1)     *at1*at1*at1*at2
                         + 3.*(C_Q1*C_Q2 - 1)     *at1*at1*at2*at2
                         + 3.*(C_Oct2 - C_Q2)     *at1*at2*at2*at2
                       + 0.75*(C_Hex2 - C_Q2*C_Q2)*at2*at2*at2*at2;
    
    double rc2   =  r2 + a02*(1. + 2.*u) + delta_a2*u + (delta_a2_nnlo+delta_a4_lo)*u2;
    *rc          = sqrt(rc2);
    double divrc = 1.0/(*rc);
    *drc_dr      = divrc*(r - (a02 + 0.5*delta_a2)*u2 - (delta_a2_nnlo+delta_a4_lo)*u3);
    *d2rc_dr2    = divrc*(1 + (2.*a02 + delta_a2)*u3
			  + 3*(delta_a2_nnlo+delta_a4_lo)*u4 - (*drc_dr)*(*drc_dr));

  } else {

    double a0  = at1 + at2;
    double a12 = at1 - at2;
    
    double c_ss_nlo = -1.125*a0*a0 -(0.125+0.5+nu)*a12*a12 + 1.25*X12*a0*a12;

    double c_ss_nnlo = - (189./32. + 417./32.*nu              )    *a0 *a0
                       + ( 11./32. - 127./32.*nu + 0.375*nu*nu)    *a12*a12
                       + ( 89./16. -  2.625*nu              )*X12*a0 *a12;

    
    double rc2   = r2 + aK2*(1. + 2.*u) + u*c_ss_nlo + u2*c_ss_nnlo;
    *rc          = sqrt(rc2);
    double divrc = 1.0/(*rc);
    *drc_dr      = r*divrc*(1-(aK2 + 0.5*c_ss_nlo)*u3 - 0.5*u4*c_ss_nnlo);	
    *d2rc_dr2    = 1./r*(*drc_dr) + r*divrc*((3.*aK2+c_ss_nlo)*u4 + 2.*c_ss_nnlo*u5);
    
  }
  
}

/**
  * Function: eob_dyn_s_get_rc_NOSPIN
  * ---------------------------------
  *   Non-spinning case -- rc = r
  *    
  *  @param[in]  r             : radial separation
  *  @param[in]  nu            : symmetric mass ratio
  *  @param[in]  at1           : spin of body 1,  chi1* X1
  *  @param[in]  at2           : spin of body 2, chi2* X2
  *  @param[in]  aK2           : kerr parameter squared
  *  @param[in]  C_Q1          : spin-induced quadrupole of body 1
  *  @param[in]  C_Q2          : spin-induced quadrupole of body 2
  *  @param[in]  C_Oct1        : spin-induced octupole of body 1
  *  @param[in]  C_Oct2        : spin-induced octupole of body 2
  *  @param[in]  C_Hex1        : spin-induced hexadecapole of body 1
  *  @param[in]  C_Hex2        : soin-induced hexadecapole of body 2
  *  @param[in]  usetidal      : tidal flag
  *  @param[out] rc            : centrifugal radius
  *  @param[out] drc_dr        : drvt rc,r
  *  @param[out] d2rc_dr2      : drvt drc,r
  *
*/
void eob_dyn_s_get_rc_NOSPIN(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, double C_Oct1, double C_Oct2, double C_Hex1, double C_Hex2, int usetidal, 
		      double *rc, double *drc_dr, double *d2rc_dr2)
{
    *rc = r;
    *drc_dr = 1.;
    *d2rc_dr2 = 0.;
}

/**
  * Function: eob_dyn_s_get_rc_NOTIDES
  * -------------------------------------
  *   LO case with C_Q1 = 0 for tidal part
  *    
  *  @param[in]  r             : radial separation
  *  @param[in]  nu            : symmetric mass ratio
  *  @param[in]  at1           : spin of body 1
  *  @param[in]  at2           : spin of body 2
  *  @param[in]  aK2           : kerr parameter squared
  *  @param[in]  C_Q1          : spin-induced quadrupole of body 1
  *  @param[in]  C_Q2          : spin-induced quadrupole of body 2
  *  @param[in]  C_Oct1        : spin-induced octupole of body 1
  *  @param[in]  C_Oct2        : spin-induced octupole of body 2
  *  @param[in]  C_Hex1        : spin-induced hexadecapole of body 1
  *  @param[in]  C_Hex2        : soin-induced hexadecapole of body 2
  *  @param[in]  usetidal      : tidal flag
  *  @param[out] rc            : centrifugal radius
  *  @param[out] drc_dr        : drvt rc,r
  *  @param[out] d2rc_dr2      : drvt drc,r
  * 
  */
void eob_dyn_s_get_rc_NOTIDES(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, double C_Oct1, double C_Oct2, double C_Hex1, double C_Hex2, int usetidal, 
		      double *rc, double *drc_dr, double *d2rc_dr2)
{

  double u   = 1./r;
  double u2  = u*u;
  double u3  = u*u2;
  double r2  = r*r;
    
  if (usetidal) {    
    /*  We set C_Q1=C_Q2=0, but keep centrifugal radius */
     
    double a02  = 2.*at1*at2; 
    double rc2  = r2 + a02*(1.+2.*u); 
    *rc         = sqrt(rc2); 
    *drc_dr     = r/(*rc)*(1.-a02*u3); 
    *d2rc_dr2   = 1./(*rc)*(1.-(*drc_dr)*r/(*rc)*(1.-a02*u3)+2.*a02*u3); 
      
  } else {

    double X12 = sqrt(1.-4.*nu);   
    double c_ss_nlo = (- at2*at2*(1.25 + 1.25*X12 + 0.5*nu) - at1*at1*(1.25 - 1.25*X12 + 0.5*nu) + at1*at2*(-2.+nu));
    double rc2   = r2 + aK2*(1. + 2.*u) + u*c_ss_nlo;
    *rc          = sqrt(rc2);
    double divrc = 1.0/(*rc);
    *drc_dr      = r*divrc*(1-(aK2 + 0.5*c_ss_nlo)*u3);	
    *d2rc_dr2    = divrc*(1.-(*drc_dr)*r*divrc*(1.-(aK2+0.5*c_ss_nlo)*u3)+ (2.*aK2 + c_ss_nlo)*u3);
  }
  
}

/**
 * Function: eob_dyn_s_rc_add_QOH_drvts
 * ------------------------------------
 *   This routine implements the additional terms in the derivatives of rc
 *   that arise if the quadrupole, octupole and hexapole coefficients varies with u.
 *   In particular, the coefficients
 *   a02, delta_a2, delta_a2_nnlo, delta_a4_lo
 *   in the definition of rc
 *   r_c^2   =  r^2 + a02 (1 + 2 u) + delta_a2 u + (delta_a2_nnlo + delta_a4_lo) u^2 
 *   all become functions of u(r), and terms proportional top them need to be added to 
 *   drc_dr, d2rc_dr2  
 *  
 *   This is needed only if the f-mode resonance model is included. 
 *   More or less terms are needed depending on the option for rc:
 *   CENTRAD_LO, CENTRAD_NLO, CENTRAD_NNLO, CENTRAD_NNLOS4, CENTRAD_NOSPIN, CENTRAD_NOTIDES, ...
 *
 *   @note 2nd drvts are not used in the dynamics (but only LSO computation)
 *       Not yet implemented for the dressing factors, and also not below here.
 *       1st drvts of Oct_i and Hex_i are also omitted for simplicity.
 * 
 *   @param[in]  dyn           : dynamics structure
 *   @param[in]  rc            : centrifugal radius
 *   @param[in]  r             : radial separation
 *   @param[in]  at1           : spin of body 1
 *   @param[in]  at2           : spin of body 2
 *   @param[out] drc_dr        : drvt rc,r
 *   @param[out] d2rc_dr2      : drvt drc,r
*/
void eob_dyn_s_rc_add_QOH_drvts(Dynamics *dyn, double rc, double r,
				double at1, double at2, 
				double *drc_dr, double *d2rc_dr2)
{
  if (!(EOBPars->use_tidal_fmode_model)) return;
  if (EOBPars->centrifugal_radius == CENTRAD_NOSPIN) return;
  if (EOBPars->centrifugal_radius == CENTRAD_NOTIDES) return;

  const double u = 1./r;
  const double u2 = SQ(u);
  const double u3 = u*u2;
  const double u4 = u2*u2;
  const double divrc = 1./rc;
  const double nu = EOBPars->nu; 
  const double X12 = sqrt(1.-4.*nu);
  const double at12 = SQ(at1);
  const double at22 = SQ(at2);
  
  const double coef_a02 = 0.5*(1+2*u)*divrc;
  const double coef_delta_a2 = 0.5*r*divrc;
  const double coef_delta_a2_nnlo = u2;
  
  /* Compute r-drvts of dressed quadrupoles */
  const double dCQ1_dr = - u2 * dyn->dressed_C_Q1_u;
  const double dCQ2_dr = - u2 * dyn->dressed_C_Q2_u;
  /* const double d2CQ1_dr2 = 2 * u3 * dyn->dressed_C_Q1_u + u4 * dyn->dressed_C_Q1_uu; */
  /* const double d2CQ2_dr2 = 2 * u3 * dyn->dressed_C_Q2_u + u4 * dyn->dressed_C_Q2_uu; */

  /* Add drvts terms from a02 */
  const double da02_dr = dCQ1_dr * at12 + dCQ2_dr * at22;
  /* double d2a02_dr2 = d2CQ1_dr2 * at12 + d2CQ2_dr2 * at22; */

  const double drc_dr_a02_term = coef_a02 * da02_dr;
  /* double d2rc_dr2_a02_term = 0; */
 
  *drc_dr += drc_dr_a02_term;
  /* *d2rc_dr2 += 0; */
  
  if (EOBPars->centrifugal_radius == CENTRAD_LO) 
    return;
  
  /* Drvts terms from delta_a2 */
  const double ddelta_a2_dr = X12*(at12*dCQ1_dr - at22*dCQ2_dr) + 3.*da02_dr;
  /* double d2delta_a2_dr2 = 0; */
  
  const double drc_dr_delta_a2_term = coef_delta_a2 * ddelta_a2_dr;
  /* double d2rc_dr2_delta_a2_term = 0; */

  *drc_dr += drc_dr_delta_a2_term;
  /* *d2rc_dr2 += 0; */

  if (EOBPars->centrifugal_radius == CENTRAD_NLO) 
    return;
  
  /* Drvts terms from delta_a2_nnlo */
  const double ddelta_a2_nnlo_dr = ( (387 - 207*nu)*da02_dr + 163.*X12*(at12*dCQ1_dr - at22*dCQ2_dr) )/28;

  const double drc_dr_delta_a2_nnlo_term = coef_delta_a2_nnlo * ddelta_a2_nnlo_dr;
  /* double d2rc_dr2_delta_a2_nnlo_term = 0; */

  *drc_dr += drc_dr_delta_a2_nnlo_term;
  /* *d2rc_dr2 += 0; */
  
  if (EOBPars->centrifugal_radius == CENTRAD_NNLO) 
    return;
  
  //TODO Correction for the following case not implement, we stop at NNLO.
  // if (EOBPars->centrifugal_radius == CENTRAD_NNLOS4)  
  /* Compute r-drvts of octupole and hexapole */
  /* Drvts terms from delta_a4_lo */
  /* double drc_dr_delta_a4_lo_term = 0; */
  /* double d2rc_dr2_delta_a4_lo_term = 0; */
  /* *drc_dr += drc_dr_delta_a4_lo_term; */
  /* *d2rc_dr2 += d2rc_dr2_delta_a4_lo_term; */
  //return;

}

// function to evaluate r_omega
/**
  * Function: eob_dyn_get_romg
  * ---------------------------
  *   Compute r_omega
  * 
  *  @param[in]  r             : radial separation
  *  @param[in]  prstar        : radial momentum
  *  @param[in]  pphi          : orbital angular momentum
  *  @param[in]  dyn           : dynamics structure
  *  
  *  @return r_omg         : r_omega
  *
*/
double eob_dyn_get_romg(double r, double prstar, double pphi, Dynamics *dyn)
{
  /* Unpack values */
  const double nu = EOBPars -> nu;
  const double S = EOBPars -> S;
  const double Sstar = EOBPars -> Sstar;
  const double chi1 = EOBPars -> chi1;
  const double chi2 = EOBPars -> chi2;
  const double X1 = EOBPars -> X1;
  const double X2 = EOBPars -> X2;
  const double c3 = EOBPars -> cN3LO;
  const double aK2 = EOBPars -> aK2;
  const double a1 = EOBPars -> a1;
  const double a2 = EOBPars -> a2;
  const double C_Q1 = EOBPars -> C_Q1;
  const double C_Q2 = EOBPars -> C_Q2;
  const double C_Oct1 = EOBPars -> C_Oct1;
  const double C_Oct2 = EOBPars -> C_Oct2;
  const double C_Hex1 = EOBPars -> C_Hex1;
  const double C_Hex2 = EOBPars -> C_Hex2;
  const int usetidal = EOBPars -> use_tidal;
  const int usespins = EOBPars -> use_spins;
  
  const double pphi2  = pphi*pphi;
  
  /* Compute Metric, centrifugal radius and Hamiltonan */
  double A, dA, rc, drc_dr, pl_hold;

  if (usespins) {
    eob_metric_s(r, 0., dyn, &A, &pl_hold, &dA, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);

    eob_dyn_s_get_rc(r, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc, &drc_dr, &pl_hold);
  } else {
    eob_metric(r, 0., dyn, &A, &pl_hold, &dA, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);

    rc = r;
    drc_dr = 1.;
  }
  
  double uc     = 1./rc;
  double uc2    = uc*uc;
    
  /* Compute here the new r_omg radius
     Compute same quantities with prstar=0. This to obtain psi.
     Procedure consistent with the nonspinning case. */
  double psic, r_omg;
  if (usespins) {
    double ggm0[26];
    eob_dyn_s_GS(r, rc, drc_dr, 0., aK2, 0., pphi, nu, chi1, chi2, X1, X2, c3, ggm0);
  
    const double GS_0       = ggm0[2];
    const double GSs_0      = ggm0[3];
    const double dGS_dr_0   = ggm0[6];
    const double dGSs_dr_0  = ggm0[7];
    const double Heff_orb_0 = sqrt(A*(1.0 + pphi2*uc2));    /* effective Hamiltonian H_0^eff */
    const double Heff_0     = Heff_orb_0 + (GS_0*S + GSs_0*Sstar)*pphi;
    const double H0         = sqrt(1.0 + 2.0*nu*(Heff_0 - 1.0) );
    const double ooH0       = 1./H0;
    const double Gtilde     = GS_0*S     + GSs_0*Sstar;
    const double dGtilde_dr = dGS_dr_0*S + dGSs_dr_0*Sstar;
    const double duc_dr     = -uc2*drc_dr;
    psic       = fabs((duc_dr + dGtilde_dr*rc*sqrt(A/pphi2 + A*uc2)/A)/(-0.5*dA));
    r_omg      = pow( ((1./sqrt(rc*rc*rc*psic))+Gtilde)*ooH0, -2./3. );

  } else {
    double V = A*(1. + SQ(pphi/r));
    psic  = 2.*(1. + 2.*nu*(sqrt(V) - 1.))/(SQ(r)*dA);
    r_omg = r*pow(psic,1./3.);
  }

  return r_omg;
}

/**
 * Function: eob_dyn_fLR
 * ---------------------
 *  Root function to compute light-ring
 *  @note This is for the non-spinning case
 * 
 *  @param[in]  r             : radial separation
 *  @param[in]  params        : dummy parameter
 * 
 *  @return     fLR(r) = 0
*/
double eob_dyn_fLR(double r, void  *params)
{
  Dynamics *dyn = params;     
  double A, dA, pl_hold;
  //if (EOBPars->use_spins) eob_metric_s(r, dyn, &A,&B,&dA,&d2A,&dB);
  //else
  eob_metric(r, 0., dyn, &A, &pl_hold, &dA, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
  double u = 1./r;
  double dA_u = (-dA)*SQ(r);
  return A + 0.5 * u * dA_u;
}

/**
  * Function: eob_dyn_adiabLR
  * -------------------------
  *  Root finder for adiabatic light-ring
  *  Look for r in [1.8, 5.6] for BHNS, [2.1, 5.9] for BNS
  *  and [1.8, 3.1] for BBH
  * 
  *  @note This is for the non-spinning case
  * 
  *  @param[in]  dyn           : dynamics structure
  *  @param[out] rLR           : light-ring radius
  * 
  *  @return     status        : GSL status
*/
int eob_dyn_adiabLR(Dynamics *dyn, double *rLR)
{
  int status;
  int iter = 0, max_iter = 200;
  const double epsabs = 0.; /* if converges, precision is |r-r*| = epsabs + epsrel r*  */
  const double epsrel = 1e-10; 
  const gsl_root_fsolver_type *T;
  double x, x_lo, x_hi;
    
  /* Set interval to search root */
  if (EOBPars->use_tidal) {
    /* Tides are always temporarily set as = NNLO to compute LR, 
       But we may want to define different searches intervals */
    const int tides = EOBPars->use_tidal;
    if (tides == TIDES_TEOBRESUM_BHNS) {
      /* BHNS */
      x_lo = 1.8; 
      x_hi = 5.6; // nu~1/4 kappaT2 ~ 600
    } else {
      /* BNS */
      x_lo = 2.1; // nu~1/4 kappaT2 ~ 12  
      x_hi = 5.9; // nu~1/4 kappaT2 ~ 600 
    }
  } else {
    /* BBH */
    x_lo = 1.8; // 1.818461553848201e+00 nu = 1/4
    x_hi = 3.1; // 3. nu = 0 
    /* x_lo = 0.9*eob_approxLR(EOBPars->nu); 
       x_hi = 1.1*eob_approxLR(EOBPars->nu); */
  }  
  
  gsl_root_fsolver *s;
  gsl_function F;
  F.function = &eob_dyn_fLR;
  F.params = dyn;
  //T = gsl_root_fsolver_bisection;
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      x      = gsl_root_fsolver_root (s);
      x_lo   = gsl_root_fsolver_x_lower (s);
      x_hi   = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi, epsabs, epsrel);
    }
  while (status == GSL_CONTINUE && iter < max_iter);
  gsl_root_fsolver_free (s);

  *rLR = 0.;
  if (isfinite(x)) *rLR = x;

  //if (status == ???) {
  //  return ROOT_ERRORS_BRACKET;
  //}
  if (status == GSL_SUCCESS) {
    return ROOT_ERRORS_NO;
  } 
  if (iter >= max_iter) {
    return ROOT_ERRORS_MAXITS;
  }
  if (status != GSL_SUCCESS) {
    return ROOT_ERRORS_NOSUCC;
  }
    
  return status;
}

/**
 * Function: eob_dyn_fLSO
 * ----------------------
 *  Root function to compute LSO
 *  @note This is for the non-spinning case
 * 
 *  @param[in]  r             : radial separation
 *  @param[in]  params        : dummy parameter
 * 
 *  @return     fLSO(r) = 0
*/
double eob_dyn_fLSO(double r, void  *params)
{
  Dynamics *dyn = params;    
  double A,B,dA,d2A,dB,pl_hold;
  //if (EOBPars->use_spins) eob_metric_s(r, dyn, &A,&B,&dA,&d2A,&dB);
  //else                
  eob_metric (r, 0., dyn, &A, &B, &dA, &d2A, &dB, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
  double u = 1./r;
  double u2  = SQ(u);
  double dA_u = (-dA)*SQ(r);
  double d2A_u = d2A*SQ(r)*SQ(r) + 2*dA*SQ(r)*r;
  dB = u2*dA_u + 2.*A*u;
  double d2B = d2A_u*u2 + 4.*u*dA_u + 2*A;
  return ( dA_u*d2B - d2A_u*(dB) );
}

/**
  * Function: eob_dyn_adiabLSO
  * --------------------------
  *  Root finder for adiabatic LSO
  *  Look for r in [4.5, 6.2] for BBH, [4.5, 36] for if
  *  tides are included
  * 
  *  @note This is for the non-spinning case
  * 
  *  @param[in]  dyn           : dynamics structure
  *  @param[out] rLSO          : LSO radius
  * 
  *  @return     status        : GSL status
*/
int eob_dyn_adiabLSO(Dynamics *dyn, double *rLSO)
{
  int status;
  int iter = 0, max_iter = 200;
  const double epsabs = 0.; /* if converges, precision is |r-r*| = epsabs + epsrel r*  */
  const double epsrel = 1e-10; 
  const gsl_root_fsolver_type *T;
  double x;
  double x_lo = 4.5; // 4.532648e+00 nu= 1/4
  double x_hi = 6.2; // 6 nu=0 
  if (EOBPars->use_tidal) x_hi = 36.; 
  
  gsl_root_fsolver *s;
  gsl_function F;
  F.function = &eob_dyn_fLSO;
  F.params = dyn;
  //T = gsl_root_fsolver_bisection;
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      x      = gsl_root_fsolver_root (s);
      x_lo   = gsl_root_fsolver_x_lower (s);
      x_hi   = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi, epsabs, epsrel);
    }
  while (status == GSL_CONTINUE && iter < max_iter);
  gsl_root_fsolver_free (s);

  *rLSO = 0.;
  if (isfinite(x)) *rLSO = x;

  //if (status == ???) {
  //  return ROOT_ERRORS_BRACKET;
  //}
  if (status == GSL_SUCCESS) {
    return ROOT_ERRORS_NO;
  } 
  if (iter >= max_iter) {
    return ROOT_ERRORS_MAXITS;
  }
  if (status != GSL_SUCCESS) {
    return ROOT_ERRORS_NOSUCC;
  }
  
  return status;
}

/**
 * Function: eob_dyn_fLSO_s
 * ------------------------
 *  Root function to compute LSO in spinning case
 * 
 *  @param[in]  x             : vector of variables (r, pphi)
 *  @param[in]  params        : dynamics structure
 *  @param[out] f             : vector of functions (dHeff_dr, d2Heff_dr2)
 * 
 *  @return     GSL_SUCCESS
*/
int eob_dyn_fLSO_s (const gsl_vector *x, void * params, gsl_vector *f) {

  Dynamics *dyn = params;
  double r    = gsl_vector_get(x,0);
  double pphi = gsl_vector_get(x,1);

  const double nu    = EOBPars->nu;
  const double S     = EOBPars->S;
  const double Sstar = EOBPars->Sstar;
  const double chi1  = EOBPars->chi1;
  const double chi2  = EOBPars->chi2;
  const double X1    = EOBPars->X1;
  const double X2    = EOBPars->X2;
  const double c3    = EOBPars->cN3LO;
  const double aK2   = EOBPars->aK2;
  const double a1    = EOBPars->a1;
  const double a2    = EOBPars->a2;
  const double C_Q1  = EOBPars->C_Q1;
  const double C_Q2  = EOBPars->C_Q2;
  const double C_Oct1 = EOBPars->C_Oct1;
  const double C_Oct2 = EOBPars->C_Oct2;
  const double C_Hex1 = EOBPars->C_Hex1;
  const double C_Hex2 = EOBPars->C_Hex2;
  const int usetidal = EOBPars->use_tidal;
  const int usespins = EOBPars->use_spins;
  
  double A, B, dA, d2A, dB, pl_hold;
  eob_metric_s(r, 0., dyn, &A, &B, &dA, &d2A, &dB, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold); 

  /* Compute centrifugal radius */
  double rc, drc_dr, d2rc_dr2;
  eob_dyn_s_get_rc(r, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc, &drc_dr, &d2rc_dr2, &pl_hold, &pl_hold, &pl_hold);
   
  /* Compute Hamiltonian */
  double Heff_orb, Heff, H, dHeff_dr, dHeff_dprstar, d2Heff_dprstar20, dHeff_dpphi, d2Heff_dr2;
  eob_ham_s(nu, r, rc, drc_dr, d2rc_dr2, pphi, 0., S, Sstar, chi1, chi2, X1, X2, aK2, c3, A, dA, d2A, 0., 0., 0., 0., 0., 
	    &H, &Heff, &Heff_orb, &dHeff_dr, &dHeff_dprstar, &dHeff_dpphi, &d2Heff_dprstar20, &d2Heff_dr2);
   
  gsl_vector_set (f, 0, dHeff_dr);
  gsl_vector_set (f, 1, d2Heff_dr2);
  
  return GSL_SUCCESS;
}

/**
 * Function: eob_dyn_LSO_s
 * -----------------------
 *  Root finder for LSO in spinning case
 *   
 *  @param[in]  dyn           : dynamics structure
 *  @param[out] rLSO          : LSO radius
 *  @param[out] pphiLSO       : LSO pphi
 * 
 *  @return     status        : GSL status
*/
int eob_dyn_LSO_s(Dynamics *dyn, double *rLSO, double *pphiLSO)
{
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  int status;
  size_t i, iter = 0;
  
  const size_t n = 2;

  gsl_multiroot_function f = {&eob_dyn_fLSO_s, n, dyn};

  double x_init[2] = {3.0, 2.0};
  gsl_vector *x = gsl_vector_alloc(n);

  gsl_vector_set (x, 0, x_init[0]);
  gsl_vector_set (x, 1, x_init[1]);

  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, 2);
  gsl_multiroot_fsolver_set (s, &f, x);

  int max_iter = 200;
  const double epsrel = 1e-10;
  
  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate(s);
      
      if (status) break;

      status = gsl_multiroot_test_residual (s->f, epsrel);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  double x0 = gsl_vector_get(s->x,0);
  double x1 = gsl_vector_get(s->x,1);
  
  *rLSO = 0.;
  if (isfinite(x0)) *rLSO = x0;
  
  *pphiLSO = 0.;
  if (isfinite(x1)) *pphiLSO = x1;
  
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);
  
  //if (status == ???) {
  //  return ROOT_ERRORS_BRACKET;
  //}
  if (status == GSL_SUCCESS) {
    return ROOT_ERRORS_NO;
  } 
  if (iter >= max_iter) {
    return ROOT_ERRORS_MAXITS;
  }
  if (status != GSL_SUCCESS) {
    return ROOT_ERRORS_NOSUCC;
  }
  
  return status;
}

/**
 * Function: eob_spin_dyn_alpha
 * ----------------------------
 *   Compute alpha from Lhat
 * 
 *   @param[in]  Lhx           : x-component of hatL
 *   @param[in]  Lhy           : y-component of hatL
 *   @param[in]  Lhz           : z-component of hatL
 * 
 *   @return     alpha         : alpha euler angle
*/
double eob_spin_dyn_alpha(double Lhx, double Lhy, double Lhz)
{
  return atan2(Lhy,Lhx);
}

/**
 * Function: alpha_initial_condition
 * ---------------------------------
 *   Compute the alpha initial condition by
 *   evaluating the lim t-->0 of alpha(t) at NLO.
 *   See Eq. A5 of https://arxiv.org/2111.03675
 *
 *   @param[in]  eobp         : eob parameters
 * 
 *   @return    alpha0        : initial alpha   
*/
double alpha_initial_condition(EOBParameters *eobp)
{
  double q 	   = eobp->q;
  double nu    = eobp->nu;
  double f0 	 = eobp->initial_frequency;
  double chi1x = eobp->chi1x;
  double chi1y = eobp->chi1y;
  double chi1z = eobp->chi1z;
  double chi2x = eobp->chi2x;
  double chi2y = eobp->chi2y;
  double chi2z = eobp->chi2z;

  /*Convert f0 to v used here */
  double v = cbrt(Pi*f0);
  if (!(eobp->use_geometric_units)) v = cbrt(Pi*f0/time_units_factor(eobp->M));

  /* precompute some quantities and coefficients*/
  double alpha_x_NLO, alpha_y_NLO;

  alpha_x_NLO =  -3.*q*(chi1y + chi2y*q)*(chi1z + chi2z*q)*v + q*(chi1y*(4. + 3.*q) + chi2y*q*(3. + 4.*q));
  alpha_y_NLO =  +3.*q*(chi1x + chi2x*q)*(chi1z + chi2z*q)*v - q*(chi1x*(4. + 3.*q) + chi2x*q*(3. + 4.*q));

  return atan2(alpha_y_NLO, alpha_x_NLO);

}

/**
 * Function: eob_spin_dyn_beta
 * ---------------------------
 *   Compute beta from Lhat
 * 
 *   @param[in]  Lhx           : x-component of hatL
 *   @param[in]  Lhy           : y-component of hatL
 *   @param[in]  Lhz           : z-component of hatL
 * 
 *   @return     beta          : beta euler angle
*/
double eob_spin_dyn_beta(double Lhx, double Lhy, double Lhz)
{
  /* possible numerical errors give Lhz = 1.000....X > 1., get around the issue */
  if (Lhz > 1)
    Lhz = 1.;
  if (Lhz < -1)
    Lhz = -1.;
  return acos(Lhz);
}

/**
  * Function: eob_spin_dyn_rhs_PN
  * -----------------------------
  *   Compute the r.h.s. of the PN precessing equations
  *   https://arxiv.org/abs/1307.4418
  *   https://arxiv.org/abs/1703.03967
  * 
  *   @param[in]  t             : time
  *   @param[in]  y             : vector of variables
  *   @param[out] dy            : vector of derivatives
  *   @param[in]  d             : dynamics structure
  * 
  *   @return     status        : GSL status
*/
int eob_spin_dyn_rhs_PN(double t, const double y[], double dy[], void *d)
{
  (void)(t); /* avoid unused parameter warning */
  Dynamics *dyn = d;  

  const double q   = EOBPars->q; // Assume q = MA/MB >=1
  const double nu  = EOBPars->nu; 
  const double nu2 = nu*nu;
  const double nu3 = nu2*nu;
  const double nu4 = nu3*nu;
  const double nu5 = nu4*nu;

  const double MA = nu_to_X1(nu); 
  const double MA4= SQ(MA)*SQ(MA);
  const double MB = 1. - MA;
  const double MB4= SQ(MB)*SQ(MB);
  const double dm = MA - MB; 
  const double ma_o_mb = MA/MB;
  const double mb_o_ma = MB/MA;

  const double Pi2 = SQ(Pi);
  const double Pi3 = Pi2*Pi;
  const double oothree = 0.3333333333333333; // 1/3
  const double eleven_o_three = 3.6666666666666665; // 11/3

  /* Unpack y */ 
  double SA[IN3], SB[IN3], Lh[IN3]; 
  SA[Ix] = y[EOB_EVOLVE_SPIN_SxA];
  SA[Iy] = y[EOB_EVOLVE_SPIN_SyA];
  SA[Iz] = y[EOB_EVOLVE_SPIN_SzA];

  SB[Ix] = y[EOB_EVOLVE_SPIN_SxB];
  SB[Iy] = y[EOB_EVOLVE_SPIN_SyB];
  SB[Iz] = y[EOB_EVOLVE_SPIN_SzB];

  Lh[Ix] = y[EOB_EVOLVE_SPIN_Lx];
  Lh[Iy] = y[EOB_EVOLVE_SPIN_Ly];
  Lh[Iz] = y[EOB_EVOLVE_SPIN_Lz];

  const double omg = y[EOB_EVOLVE_SPIN_Momg]; // M omega
  const double lnomg = log(omg);
  const double v  = pow(omg, oothree); // division by M !
  const double v2 = v*v;
  const double v3 = v2*v;
  const double v4 = v3*v;
  const double v5 = v4*v;
  const double v6 = v5*v;
  const double v7 = v6*v;
  const double v8 = v7*v;
  const double v9 = v8*v;
  const double v10= v9*v;
  /* rhs */

  /* spins and Lhat */
  
  // Note: alpha and beta are not actually evolved,
  //       their rhs is left to 0 
  for(int v=0; v<EOB_EVOLVE_SPIN_NVARS; v++)
    dy[v] = 0.; 

  double qSAB[IN3], SABq[IN3];
  double OmgA[IN3], OmgB[IN3];
  double Omg_x_SA[IN3], Omg_x_SB[IN3];
  
  /* NLO */
  for(int a=Ix; a<IN3; a++) qSAB[a] = 1./q*SA[a] + SB[a];
  for(int a=Ix; a<IN3; a++) SABq[a] = SA[a] + SB[a]*q;

  double qSABLh, SABqLh;
  vect_dot3(qSAB, Lh, &qSABLh);
  vect_dot3(SABq, Lh, &SABqLh);
  
  /* compute the mass-dep. coefficients at first call only */
  static double v5_cA, v5_cB, v7_cA, v7_cB, v9_cA, v9_cB;
  static double csA, csB, csAL, csBL, csANL, csBNL;
  static double L2PN_v2, L2PN_v4;

  if (EOBPars->firstcall[FIRSTCALL_SPINDYN]){
    EOBPars->firstcall[FIRSTCALL_SPINDYN] = 0;

    v5_cA = nu*(2+1.5/q);
    v5_cB = nu*(2+1.5*q);
    v7_cA = 0.5625 + 1.25*nu - 0.04166666666666666*nu2 + dm*(-0.5625+0.625*nu);
    v7_cB = 0.5625 + 1.25*nu - 0.04166666666666666*nu2 - dm*(-0.5625+0.625*nu);
    v9_cA = 0.84375 + 0.1875*nu - 3.28125*nu2 - 0.02083333333333*nu3
			      + dm*(-0.84375+4.875*nu-0.15625*nu2);
    v9_cB = 0.84375 + 0.1875*nu - 3.28125*nu2 - 0.02083333333333*nu3
			      - dm*(-0.84375+4.875*nu-0.15625*nu2);

    csA  = -0.25*(3.+1./MA);
    csB  = -0.25*(3.+1./MB);
    csAL = -0.08333333333333333*(1.+27./MA);
    csBL = -0.08333333333333333*(1.+27./MB);
    csANL= 0.0625*(-31.+7/MA) + 0.0208333333333333*nu*(9. + 22/MA);
    csBNL= 0.0625*(-31.+7/MB) + 0.0208333333333333*nu*(9. + 22/MB);

    L2PN_v2 = 1.5+0.1666666666666667*nu;
    L2PN_v4 = 3.375 - 2.375*nu + 0.04166666666666666*nu2;
  }

  /* NLO */
  /* 
  for(int a=Ix; a<IN3; a++) 
  OmgA[a] = (v5*v5_cA - 1.5*v6*qSABLh)*Lh[a] + 0.5*v6*SB[a];
  for(int a=Ix; a<IN3; a++) 
  OmgB[a] = (v5*v5_cB - 1.5*v6*SABqLh)*Lh[a] + 0.5*v6*SA[a];
  
  vect_cross3(OmgA, SA, Omg_x_SA);
  vect_cross3(OmgB, SB, Omg_x_SB);
  
  dy[EOB_EVOLVE_SPIN_SxA] += Omg_x_SA[Ix];
  dy[EOB_EVOLVE_SPIN_SyA] += Omg_x_SA[Iy];
  dy[EOB_EVOLVE_SPIN_SzA] += Omg_x_SA[Iz];
  
  dy[EOB_EVOLVE_SPIN_SxB] += Omg_x_SB[Ix];
  dy[EOB_EVOLVE_SPIN_SyB] += Omg_x_SB[Iy];
  dy[EOB_EVOLVE_SPIN_SzB] += Omg_x_SB[Iz];
  
  const double v_o_nu = v/nu;
  dy[EOB_EVOLVE_SPIN_Lx] += -(v_o_nu) * ( Omg_x_SA[Ix] + Omg_x_SB[Ix] );
  dy[EOB_EVOLVE_SPIN_Ly] += -(v_o_nu) * ( Omg_x_SA[Iy] + Omg_x_SB[Iy] );
  dy[EOB_EVOLVE_SPIN_Lz] += -(v_o_nu) * ( Omg_x_SA[Iz] + Omg_x_SB[Iz] );  
  */

  
  /* N4LO */
  double OmgANLO[IN3],  OmgBNLO[IN3],  OmgANNLO[IN3],  OmgBNNLO[IN3],  OmgAN4LO[IN3],  OmgBN4LO[IN3];
  double SdotANLO[IN3], SdotBNLO[IN3], SdotANNLO[IN3], SdotBNNLO[IN3], SdotAN4LO[IN3], SdotBN4LO[IN3];
  double LNdotN4LO[IN3];

  for(int a=Ix; a <IN3; a++){
    OmgANLO[a] = (v5*v5_cA  - 1.5*v6*qSABLh)*Lh[a] + 0.5*v6*SB[a];
    OmgBNLO[a] = (v5*v5_cB  - 1.5*v6*SABqLh)*Lh[a] + 0.5*v6*SA[a];
    OmgANNLO[a]= OmgANLO[a]  + v7*v7_cA*Lh[a];
    OmgBNNLO[a]= OmgBNLO[a]  + v7*v7_cB*Lh[a];
    OmgAN4LO[a]= OmgANNLO[a] + v9*v9_cA*Lh[a];
    OmgBN4LO[a]= OmgBNNLO[a] + v9*v9_cB*Lh[a];
  }  

  vect_cross3(OmgANLO,  SA, SdotANLO);
  vect_cross3(OmgANNLO, SA, SdotANNLO);
  vect_cross3(OmgAN4LO, SA, SdotAN4LO);
  vect_cross3(OmgBNLO,  SB, SdotBNLO);
  vect_cross3(OmgBNNLO, SB, SdotBNNLO);
  vect_cross3(OmgBN4LO, SB, SdotBN4LO);

  dy[EOB_EVOLVE_SPIN_SxA] = SdotAN4LO[Ix];
  dy[EOB_EVOLVE_SPIN_SyA] = SdotAN4LO[Iy];
  dy[EOB_EVOLVE_SPIN_SzA] = SdotAN4LO[Iz];
  dy[EOB_EVOLVE_SPIN_SxB] = SdotBN4LO[Ix];
  dy[EOB_EVOLVE_SPIN_SyB] = SdotBN4LO[Iy];
  dy[EOB_EVOLVE_SPIN_SzB] = SdotBN4LO[Iz];

  // Lhdot
  const double L2PN = 1. + v2*L2PN_v2 + v4*L2PN_v4;
  const double v_o_nu = v/nu;

  double SALh, SBLh;
  double dSBNLOSA,  dSANLOSB;
  double dSANNLOLh, dSBNNLOLh;
  vect_dot3(SA, Lh, &SALh);
  vect_dot3(SB, Lh, &SBLh);  
  vect_dot3(SdotBNLO,  SA, &dSBNLOSA);
  vect_dot3(SdotANLO,  SB, &dSANLOSB);
  vect_dot3(SdotANNLO, Lh, &dSANNLOLh);
  vect_dot3(SdotBNNLO, Lh, &dSBNNLOLh);

  // Eq. (4c) of https://arxiv.org/abs/2005.05338

  LNdotN4LO[Ix] = (v_o_nu*(-SdotAN4LO[Ix] -SdotBN4LO[Ix]) -v3*(csA*SdotANNLO[Ix] + csB*SdotBNNLO[Ix])
                            -v3*(csAL*(-v_o_nu*(SdotANLO[Ix]+SdotBNLO[Ix])*SALh + Lh[Ix]*(-v_o_nu*dSBNLOSA + dSANNLOLh)) 
                                +csBL*(-v_o_nu*(SdotBNLO[Ix]+SdotANLO[Ix])*SBLh + Lh[Ix]*(-v_o_nu*dSANLOSB + dSBNNLOLh))))/L2PN; 
  LNdotN4LO[Iy] = (v_o_nu*(-SdotAN4LO[Iy] -SdotBN4LO[Iy]) -v3*(csA*SdotANNLO[Iy] + csB*SdotBNNLO[Iy])
                            -v3*(csAL*(-v_o_nu*(SdotANLO[Iy]+SdotBNLO[Iy])*SALh + Lh[Iy]*(-v_o_nu*dSBNLOSA + dSANNLOLh)) 
                                +csBL*(-v_o_nu*(SdotBNLO[Iy]+SdotANLO[Iy])*SBLh + Lh[Iy]*(-v_o_nu*dSANLOSB + dSBNNLOLh))))/L2PN;
  LNdotN4LO[Iz]= (v_o_nu*(-SdotAN4LO[Iz] -SdotBN4LO[Iz]) -v3*(csA*SdotANNLO[Iz] + csB*SdotBNNLO[Iz])
                            -v3*(csAL*(-v_o_nu*(SdotANLO[Iz]+SdotBNLO[Iz])*SALh + Lh[Iz]*(-v_o_nu*dSBNLOSA + dSANNLOLh)) 
                                +csBL*(-v_o_nu*(SdotBNLO[Iz]+SdotANLO[Iz])*SBLh + Lh[Iz]*(-v_o_nu*dSANLOSB + dSBNNLOLh))))/L2PN; 

  // Eq. (7) of https://arxiv.org/abs/2005.05338
  double LNdotN4LOLh;
  double LNdotN4LOperp[IN3];
  vect_dot3(LNdotN4LO, Lh, &LNdotN4LOLh);
  for(int a=Ix; a<IN3; a++) 
    LNdotN4LOperp[a] = LNdotN4LO[a] - LNdotN4LOLh*Lh[a];

  dy[EOB_EVOLVE_SPIN_Lx] = LNdotN4LOperp[Ix];
  dy[EOB_EVOLVE_SPIN_Ly] = LNdotN4LOperp[Iy];
  dy[EOB_EVOLVE_SPIN_Lz] = LNdotN4LOperp[Iz];
  
  
  /* dot{gamma} = + dot{alpha} cos{beta}, where cos{beta} = Lh[Iz] */ 
  const double div = SQ(Lh[Ix]) + SQ(Lh[Iy]);
  if (div == 0.)
    /* Lh is along z
       alpha is undefined, beta = 0 */
    dy[EOB_EVOLVE_SPIN_gam] = 0.;
  else
    dy[EOB_EVOLVE_SPIN_gam] = Lh[Iz] * (Lh[Ix] * dy[EOB_EVOLVE_SPIN_Ly] - Lh[Iy] * dy[EOB_EVOLVE_SPIN_Lx] )/( SQ(Lh[Ix]) + SQ(Lh[Iy]) );

  /* dot omg (Rad.React.) */
  if(EOBPars->spin_flx == SPIN_FLX_PN) {
    
    static double beta3A, beta3B, beta5A, beta5B, beta6A, beta6B, beta7A, beta7B, beta8A, beta8B;
    static double b6, b8, b9, b10, b11;
    static double a0, a2, a4_nosigma, a5_nobeta, a6_nobeta, a7_nobeta, a8_nobeta, a9, a10, a11, a10_tidal, a12_tidal;
    static double sigma4_SASB, sigma4_SALh_SBLh, sigma4_SA2, sigma4_SALh2, sigma4_SB2, sigma4_SBLh2;
    
    if (EOBPars->firstcall[FIRSTCALL_PNMOMG]){
      EOBPars->firstcall[FIRSTCALL_PNMOMG] = 0;

      sigma4_SASB      = 247./(48*nu);
      sigma4_SALh_SBLh = - 721./(48.*nu);
      sigma4_SA2       = 233./(96*SQ(MA));
      sigma4_SALh2     = -719./(96*SQ(MA));
      sigma4_SB2       = 233./(96*SQ(MB));
      sigma4_SBLh2     = - 719./(96*SQ(MB));

      beta3A = 113./12 + 25./4*mb_o_ma;
      beta3B = 113./12 + 25./4*ma_o_mb; 
      beta5A = (31319./1008-1159./24*nu) + mb_o_ma*( 809./84 - 281./8*nu);
      beta5B = (31319./1008-1159./24*nu) + ma_o_mb*( 809./84 - 281./8*nu); 
      beta6A = Pi * (75./2+151./6*mb_o_ma);
      beta6B = Pi * (75./2+151./6*ma_o_mb);
      beta7A = (130325./756 - 796069./2016*nu+100019./864*nu2)+ mb_o_ma*(1195759./18144-257023./1008*nu+2903./32*nu2);
      beta7B = (130325./756 - 796069./2016*nu+100019./864*nu2)+ ma_o_mb*(1195759./18144-257023./1008*nu+2903./32*nu2);
      beta8A = Pi*(76927./504 -220055./672*nu) + mb_o_ma*(1665./28-50483./224*nu);
      beta8B = Pi*(76927./504 -220055./672*nu) + ma_o_mb*(1665./28-50483./224*nu);

      b6  = -1712./315;  
      b8  = - 856./315*nu + 124741./4410;  
      b9  = - 6848./105* Pi;  
      b10 = 3090781./26460*nu -  2354./945*nu2 - 11821184./1964655;  
      b11 = 311233./5880*Pi - 3424./315*Pi*nu;

      a0 = 96./5*nu;
      a2 = -743./336-11./4*nu;
      a4_nosigma = 34103./18144 + 13661./2016*nu + 59./18*nu2;
      a5_nobeta  = -4159./672*Pi - 189./8*Pi*nu;
      a6_nobeta  = 16447322263./139708800 + 16./3*Pi2 - 856./105*log(16.) - 1712./105*EulerGamma
                   + nu *( 451./48*Pi2- 56198689./217728) + nu2*541./896 - nu3*5605./2592;
      a7_nobeta  = - 4415./4032*Pi + 358675./6048*Pi*nu +91495./1512*Pi*nu2;
      a8_nobeta  = 3971984677513./25427001600 + 127751./1470*Log2  - 47385./1568*Log3 + 124741./4410*EulerGamma -361./126*Pi2 + 82651980013./838252800*nu - 1712./315*nu*Log2
                   - 856./315*EulerGamma*nu  - 31495./8064*Pi2*nu + 54732199./93312*nu2- 3157./144*Pi2*nu2  - 18927373./435456*nu3 -95./3888*nu4;
      a9 =  343801320119./745113600*Pi- 13696./105*Pi*Log2 -  6848./105*Pi*EulerGamma - 51438847./48384*Pi*nu + 205./6*Pi3*nu + 42680611./145152*Pi*nu2  +  9731./1344*Pi*nu3;
      a10= 29619150939541789./36248733480960  -107638990./392931*Log2 + 616005./3136*Log3 - 11821184./1964655*EulerGamma - 21512./1701*Pi2 - 884576519037433./228843014400*nu 
          + 2105111./8820*nu*Log2 - 15795./3136*nu*Log3 + 3090781./26460*EulerGamma*nu+ 14555455./217728*Pi2*nu  + 1175999369413./914457600*nu2 - 4708./945*nu2*Log2
          - 126809./3024*Pi2*nu2 - 2354./945*EulerGamma*nu2 - 9007327699./11757312*nu3 + 9799./384*Pi2*nu3 + 51439207./1741824*nu4 - 34613./186624*nu5;
      a11= 91347297344213./81366405120*Pi+ 5069891./17640*Pi*Log2- 142155./784*Pi*Log3  + 311233./5880*Pi*EulerGamma - 1903651780081./4470681600*Pi*nu- 6848./315*Pi*nu*Log2
      - 3424./315*Pi*EulerGamma*nu - 26035./16128*Pi3*nu + 1760705531./290304*Pi*nu2 - 112955./576*Pi3*nu2 - 7030123./13608*Pi*nu3 + 49187./6048*Pi*nu4;  

      /* Tidal terms from  Eq. A21 of 
        https://arxiv.org/pdf/1402.5156.pdf */

      a10_tidal = 6.*MA4*(12. - 11.*MA)*EOBPars->LambdaAl2 + 6.*MB4*(12. - 11.*MB)*EOBPars->LambdaBl2;
      a12_tidal = MA4*(4421./56. - 12263./56.*MA + 1893./4.*MA*MA - 661./2.*MA*MA*MA)*EOBPars->LambdaAl2 + MB4*(4421./56. - 12263./56.*MB + 1893./4.*MB*MB - 661./2.*MB*MB*MB)*EOBPars->LambdaBl2;

    }

    /* PN omega_dot, taken from: 
       https://arxiv.org/abs/1307.4418 , App.A */

    const double vlo = 0.33333333333*v9;
    double SAdotLh, SBdotLh, SAdotSB, SA2, SB2;
    vect_dot3(SA, Lh, &SAdotLh);
    vect_dot3(SB, Lh, &SBdotLh);    
    vect_dot3(SA, SB, &SAdotSB);
    vect_dot3(SA, SA, &SA2);
    vect_dot3(SB, SB, &SB2);

    double a[12], b[12], beta[9];
    for (int i=0; i<12; i++) a[i] = 0;
    for (int i=0; i<12; i++) b[i] = 0;
    for (int i=0; i<9; i++) beta[i]= 0;
    
    // precomputed spin independent coefs as first call
        
    double sigma4 = sigma4_SASB*SAdotSB + sigma4_SALh_SBLh*SAdotLh*SBdotLh 
                  + sigma4_SA2*SA2 + sigma4_SALh2*SQ(SAdotLh)
                  + sigma4_SB2*SB2 + sigma4_SBLh2*SQ(SBdotLh);
    
    beta[3] = beta3A*SAdotLh + beta3B*SBdotLh; 
    beta[5] = beta5A*SAdotLh + beta5B*SBdotLh; 
    beta[6] = beta6A*SAdotLh + beta6B*SBdotLh; 
    beta[7] = beta7A*SAdotLh + beta7B*SBdotLh;
    beta[8] = beta8A*SAdotLh + beta8B*SBdotLh;
    
    b[6]  = b6;  
    b[8]  = b8;  
    b[9]  = b9;  
    b[10] = b10;  
    b[11] = b11;

    a[0] = a0;
    a[2] = a2;
    a[3] = 4*Pi-beta[3];
    a[4] = a4_nosigma - sigma4;
    a[5] = a5_nobeta  - beta[5];
    a[6] = a6_nobeta  - beta[6];
    a[7] = a7_nobeta  - beta[7];
    a[8] = a8_nobeta  - beta[8];
    a[9] = a9;
    a[10]= a10;
    a[11]= a11;
    
    // Eq.(A1) of https://arxiv.org/abs/1307.4418 for Momega
    dy[EOB_EVOLVE_SPIN_Momg] = 0.;
    for (int i=2; i<8; i++)
      dy[EOB_EVOLVE_SPIN_Momg] += (a[i] + b[i]*lnomg)*pow(omg,(double)i*oothree);
    dy[EOB_EVOLVE_SPIN_Momg] += 1.;  
    
    // Add LO 5PN tidal term if BNS evolution
    if(EOBPars->use_tidal){
      dy[EOB_EVOLVE_SPIN_Momg]+= a10_tidal*pow(omg, 10.*oothree);
    }
    dy[EOB_EVOLVE_SPIN_Momg] *= a[0]*pow(omg, eleven_o_three); // LO

  } else if (EOBPars->spin_flx == SPIN_FLX_EOB_HYBRIDv1){

    /*hybrid v1: PN expressions for j(v), u(j) + EOB flux */

    const double X1    = EOBPars->X1;
    const double X2    = EOBPars->X2;
    const double X12   = X1-X2;

    /* spin variables */
    double SAdotLh, SBdotLh;
    vect_dot3(SA, Lh, &SAdotLh);
    vect_dot3(SB, Lh, &SBdotLh);    
    const double c1    = SAdotLh/SQ(X1);
    const double c2    = SBdotLh/SQ(X2);
    const double a1    = X1*c1;
    const double a2    = X2*c2;
    const double a0    = a1+a2;
    const double aK2   = a0*a0;
    const double aAB   = a1-a2;
    const double aCq   = EOBPars->C_Q1*a1*a1 + 2*a1*a2+EOBPars->C_Q2*a2*a2; 
    const double Sl    = X1*a1 + X2*a2;
    const double Sigmal= X2*c2 - X1*c1;
    const double S     = SAdotLh + SBdotLh;
    const double Sstar = X2*a1+X1*a2;

    /* j(v), app G of https://arxiv.org/pdf/2004.06503.pdf */
    double jhat_orb   = 1.+v2*(nu/6.+1.5)
                       + v4*(nu2/24.-2.375*nu+3.375)
                       + v6*(0.005401234567901234*nu3+1.2916666666666667*nu2 + (1.7083333333333333*Pi2-47.84027777777778)*nu + 8.4375)
                       + v8*(-55./31104.*nu4 -215./1728.*nu3 + (356035./3456. - 2255.*Pi2/576.)*nu2 + nu*(-64./3.*log(16*v2) -6455./1536.*Pi2 - 128./3.*EulerGamma + 98869./5760.) + 2835./128.);
    double jhat_so    = v3*(-35./6.*Sl - 2.5*X12*Sigmal)
                       + v5*((-77./8.+427./72.*nu)*Sl+X12*(-21./8.+35./12.*nu)*Sigmal)
                       + v7*((-405./16. + 1101./16*nu - 29./16.*nu2)*Sl + X12*( -81./6. + 117./4.*nu - 15./16.*nu2)*Sigmal);
    double jhat_ss    = v4*((0.5+X12*0.5-nu)*c1*c1+2*nu*c1*c2+(0.5+X12*0.5-nu)*c2*c2);
    double jhat       = jhat_orb + jhat_so + jhat_ss; 

    /* derivative of jhat*/
    double djhatdx_orb= (nu/6.+1.5)
                       + 2.*v2*(nu2/24.-2.375*nu+3.375)
                       + 3.*v4*(0.005401234567901234*nu3+1.2916666666666667*nu2 + (1.7083333333333333*Pi2-47.84027777777778)*nu + 8.4375)
                       + 4.*v6*(-55./31104.*nu4 -215./1728.*nu3 + (356035./3456. - 2255.*Pi2/576.)*nu2 + nu*(-64./3.*log(16*v2) -6455./1536.*Pi2 - 128./3.*EulerGamma + 98869./5760.-16./3.) + 2835./128.);
    double djhatdx_so = 1.5*v*(-36./5.*Sl - 2.5*X12*Sigmal)
                       + 2.5*v3*((-77./8.+427./72.*nu)*Sl+X12*(-21./8.+35./12.*nu)*Sigmal)
                       + 3.5*v5*((-405./16. + 1101./16*nu - 29./16.*nu2)*Sl + X12*( -81./6. + 117./4.*nu - 15./16.*nu2)*Sigmal);
    double djhatdx_ss = 2*v2*((0.5+X12*0.5-nu)*c1*c1+2*nu*c1*c2+(0.5+X12*0.5-nu)*c2*c2);
    
    double djhatdx    = djhatdx_orb+djhatdx_so+djhatdx_ss;
    double djhatdomg  = 2./3.*djhatdx/v;
    double djdomg     = -0.33333333333/v4*jhat+djhatdomg/v;

    /* one over j powers*/
    double ooj       = v/jhat;
    double ooj2      = ooj*ooj;   double ooj3      = ooj2*ooj;
    double ooj4      = ooj2*ooj2; double ooj5      = ooj4*ooj;
    double ooj6      = ooj3*ooj3; double ooj7      = ooj6*ooj;
    double ooj8      = ooj4*ooj4; double ooj9      = ooj8*ooj;
    double ooj10     = ooj5*ooj5;

    /* u(j), https://arxiv.org/pdf/1812.07923.pdf 
        Eq. (12),(17) and (18)
    */
    double delta_a2_nlo  = -33./8.*a0*a0+3*aCq-0.125*(1+4*nu)*aAB*aAB+X12*(0.25*a0*aAB+EOBPars->C_Q1*a1*a1-EOBPars->C_Q2*a2*a2);
    double delta_a2_nnlo = -(4419./224+1263/224.*nu)*a0*a0+(387/28-207/28*nu)*aCq
                           + (11./32 -127/32*nu+3./8.*nu2)*aAB*aAB
                           + X12*(-(29/112+21/8*nu)*a0*aAB + 163/28*(EOBPars->C_Q1*a1*a1-EOBPars->C_Q2*a2*a2));

    double u_orb    = ooj2 + 3.*ooj4 + ooj6*(18. - 3*nu) + ooj8*(135.+(-311./3+41./16.*Pi2)*nu) 
                    + ooj10*(1134. - (163063./120.+64.*EulerGamma - 31921./1024.*Pi2+128.*Log2+64*log(ooj))*nu 
                    + (1321./12. - 205./64.*Pi2)*nu2);
    double u_so     = - 3./4.*(7.*a0+X12*aAB)*ooj5+
                    + ((-465./8+11./4*nu)*a0-(87./8+nu/4.)*X12*aAB)*ooj7
                    + ((-1269./2.+1273./8.*nu*25./32.*nu2)*a0 + (-531./4.+103./8.*nu+5./32*nu2)*X12*aAB)*ooj9;
    double u_ss     = 2*aCq*ooj6+(441./8*a0*a0+22.*aCq+(9./8-9./2*nu)*aAB*aAB+63./4*X12*a0*aAB+5./2*delta_a2_nlo)*ooj8
                    + ((9009./8.-1155./16*nu)*a0*a0+(234.-45./2.*nu)*aCq+(261./8.-2073./16.*nu-15./4.*nu2)*aAB*aAB
                    + (1557./4.-15./4.*nu)*X12*a0*aAB+29*delta_a2_nlo+3*delta_a2_nnlo)*ooj10;

    const double u  = u_orb+u_so+u_ss;
    const double r  = 1./u;

    /* Compute Metric */
    double A, B, dA, d2A, dB, pl_hold;
    eob_metric_s(r, 0., NULL, &A, &B, &dA, &d2A, &dB, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
  
    /* Compute centrifugal radius */
    double rc, drc_dr, d2rc_dr;
    eob_dyn_s_get_rc(r, nu, a1, a2, aK2, EOBPars->C_Q1, EOBPars->C_Q2, EOBPars->C_Oct1, EOBPars->C_Oct2, EOBPars->C_Hex1, EOBPars->C_Hex2, EOBPars->use_tidal, &rc, &drc_dr, &d2rc_dr);
    const double uc     = 1./rc;
    const double uc2    = uc*uc;
    const double uc3    = uc2*uc;

    /* Compute Hamiltonian */
    double Heff_orb, Heff, H, dHeff_dr, dHeff_dprstar, d2Heff_dprstar20, dHeff_dpphi;
    eob_ham_s(nu, r, rc, drc_dr, d2rc_dr, 1./ooj, 0., S, Sstar, c1, c2, X1, X2, aK2, EOBPars->cN3LO, A, dA, d2A, 0., 0., 0., 0., 0.,
	    &H, &Heff, &Heff_orb, &dHeff_dr, &dHeff_dprstar, &dHeff_dpphi, &d2Heff_dprstar20, NULL);

    /* EOB Flux */
    //double sqrtW = sqrt(A*(1.+jhat*jhat*SQ(u)));
    double psi   = 2.*(1.0 + 2.0*nu*(Heff- 1.0))/(SQ(r)*dA);
    double r_omg = 1./v2;//r*pow(psi, 1./3);
    double v_phi = r_omg*omg;
    //double flx   = eob_flx_Flux_s(v_phi*v_phi, omg, r_omg, nu*H, Heff, jhat, r, 0, 0);
    double flx   = eob_flx_Flux_s(v2, omg, r_omg, nu*H, Heff, jhat, r, 0, 0, dyn);

    if (VERBOSE)
      printf("F=%.10f, djdomg=%.10f, omg=%.10f\n",flx, djdomg, omg);
    
    if(dy[EOB_EVOLVE_SPIN_Momg]<0.)
       dyn->spins->omg_stop = 1e-10;

    dy[EOB_EVOLVE_SPIN_Momg] = flx/djdomg; //flx is dj/dt. The energy flux is omg*flx (modulo nu terms)

    if (DEBUG){
      /* print to file omega, 3.5PN nonspinning flux divided by newtonian prefactor and EOB flux */
      double en_flx_pn_LO = 32./5.*nu2*v10;
      double hat_en_flx   = 1. + (-1247./336. - 35./12.*nu)*v2 + 4*Pi*v3
                              + (-44711./9072. + 9271./504*nu+65./18.*nu2)*v4 + (-8191./672.-583./24*nu)*Pi*v5
                              + (6643739519./69854400. + 16./3.*Pi2 - 1712/105*EulerGamma - 856./105.*log(16*v2) + (-134543./7776. + 41./48.*Pi2)*nu - 94403./3024.*nu2 - 775./3024*nu3)*v6
                              + (-16285./504. + 214745./1728.*nu + 193385./3024.*nu2)*Pi*v7;
      double newt         = -32./5.*nu*pow(omg, 7./3.); 
      double en_flx       = en_flx_pn_LO*hat_en_flx;
      FILE *f;
      f = fopen("fluxes.txt", "a");
      fprintf(f, "%.10f %.10f %.10f\n", omg, en_flx/(nu*omg*newt), flx/newt); 
      fclose(f);
    }
  
  } else if (EOBPars->spin_flx == SPIN_FLX_EOB_HYBRIDv2){

    /* hybrid v2: use eob_dyn_r0_eob for r(omega), use j_circ for j(r) and EOB flux*/

    const double X1    = EOBPars->X1;
    const double X2    = EOBPars->X2;
    const double X12   = X1-X2;

    /* spin variables */
    double SAdotLh, SBdotLh;
    vect_dot3(SA, Lh, &SAdotLh);
    vect_dot3(SB, Lh, &SBdotLh);    
    const double c1    = SAdotLh/SQ(X1);
    const double c2    = SBdotLh/SQ(X2);
    const double a1    = X1*c1;
    const double a2    = X2*c2;
    const double a0    = a1+a2;
    const double aK2   = a0*a0;
    const double aAB   = a1-a2;
    const double aCq   = EOBPars->C_Q1*a1*a1 + 2*a1*a2+EOBPars->C_Q2*a2*a2; 
    const double Sl    = X1*a1 + X2*a2;
    const double Sigmal= X2*c2 - X1*c1;
    const double S     = SAdotLh + SBdotLh;
    const double Sstar = X2*a1+X1*a2;

    /*Compute r(omega) numerically by inverting hamilton Eq.*/
    double r = 1./v2;//eob_dyn_r0_eob(omg/Pi, NULL);
    double u = 1./r;
    double u2= u*u;
    double u3= u2*u;

    /* Compute j(u) on circular orbits */
    double ggm[26]; 
    double A, B, dA, d2A, dB, pl_hold, rc, drc_dr, d2rc_dr;
    eob_metric_s(r, 0., NULL, &A, &B, &dA, &d2A, &dB, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    eob_dyn_s_get_rc(r, nu, a1, a2, aK2, EOBPars->C_Q1, EOBPars->C_Q2, EOBPars->C_Oct1, EOBPars->C_Oct2, EOBPars->C_Hex1, EOBPars->C_Hex2, EOBPars->use_tidal, &rc, &drc_dr, &d2rc_dr);
    eob_dyn_s_GS(r, rc, drc_dr, d2rc_dr, aK2, 0.0, 0.0, nu, c1, c2, X1, X2, EOBPars->cN3LO, ggm);

    double GS     = ggm[2];
    double GSs    = ggm[3];  
    double dGS_dr = ggm[6];
    double dGSs_dr= ggm[7];
    double G     = GS*S + GSs*Sstar;    // tildeG = GS*S+GSs*Ss
    double dG_dr = dGS_dr*S + dGSs_dr*Sstar;

    double uc       = 1./rc;
    double uc2      = uc*uc;
    double uc3      = uc2*uc;
    double uc4      = uc3*uc;
    double dAuc2_dr = uc2*(dA-2*A*uc*drc_dr);
    
    // Quadratic equation a*x^2+b*x+c=0 
    double a_coeff = SQ(dAuc2_dr)  - 4*A*uc2*SQ(dG_dr);
    double b_coeff = 2*dA*dAuc2_dr - 4*A*SQ(dG_dr);
    double c_coeff = SQ(dA);
      
    double Delta = SQ(b_coeff) - 4*a_coeff*c_coeff;
      
    if (S==0 && Sstar==0)
  	  Delta=0;             // dG_dr=0 -> Set Delta=0 to avoid num. errors          
      
    double sol_p   = (-b_coeff + sqrt(Delta))/(2*a_coeff); 
    double sol_m   = (-b_coeff - sqrt(Delta))/(2*a_coeff);
    
    double j02=0;
    if (dG_dr > 0)
	    j02 = sol_p;
    else
	    j02 = sol_m;
    
    double j = sqrt(j02);

    //Compute hamiltonian and derivatives
    double Heff_orb, Heff, H, dHeff_dr, dHeff_dprstar, d2Heff_dprstar20, dHeff_dpphi;
    eob_ham_s(nu, r, rc, drc_dr, d2rc_dr, j, 0., S, Sstar, c1, c2, X1, X2, aK2, EOBPars->cN3LO, A, dA, d2A, 0., 0., 0., 0., 0.,
	    &H, &Heff, &Heff_orb, &dHeff_dr, &dHeff_dprstar, &dHeff_dpphi, &d2Heff_dprstar20, NULL);

    /* Compute dj_circ/dr */
    double c10,c20,c30;
    double cs10,cs20,cs30,cs40;
    c10 =  5./16.*nu;
    c20 =  51./8.*nu + 41./256.*nu2;
    c30 =  nu*EOBPars->cN3LO;

    cs10 = 3./4.   + nu/2.;
    cs20 = 27./16. + 29./4.*nu + 3./8.*nu2;
    cs30 = nu*EOBPars->cN3LO + 135./32.;
    cs40 = 2835./256.;

    // Auxiliary quantities 
    double duc_dr    = -uc2*drc_dr;
    double d2uc_dr   = -2*uc*duc_dr*drc_dr - uc2*d2rc_dr;
    double d2Auc2_dr =  2*duc_dr*dAuc2_dr/uc+uc2*(d2A-2*dA*uc*drc_dr-2*A*duc_dr*drc_dr-2*A*uc*d2rc_dr);

    // LO of Giro-Gravitomagn. coeffs and derivatives
    double GS0      = 2*u*uc2;
    double dGS0_dr  = -2*u2*uc2 + 4*u*uc*duc_dr;
    double d2GS0_dr = 4*(u3*uc2-u2*uc*duc_dr) + 4*(-u2*uc*duc_dr + u*duc_dr*duc_dr+u*uc*d2uc_dr);

    double GSs0     = 1.5*uc3;
    double dGSs0_dr = 4.5*uc2*duc_dr;
    double d2GSs0_dr= 9.*uc*duc_dr*duc_dr+4.5*uc2*d2uc_dr;

    //hatGS, hatGSs and derivatives with pr_star=0
    double hGS      = 1./(1.+c10*uc + c20*uc2 + c30*uc3);  
    double dhGS_dr  = -hGS*hGS*(c10 + 2.*c20*uc  + 3.*c30*uc2)*duc_dr;
    double d2hGS_dr = -2*hGS*dhGS_dr*(c10 + 2.*c20*uc  + 3.*c30*uc2)*duc_dr
                      -hGS*hGS*(2.*c20+6.*c30*uc)*duc_dr*duc_dr +
                      dhGS_dr/duc_dr*d2uc_dr;
           dGS_dr   = dGS0_dr*hGS+GS0*dhGS_dr;
    double d2GS_dr  = d2GS0_dr*hGS+2*dGS0_dr*dhGS_dr+GS0*d2hGS_dr;

    double hGSs     = 1./(1.+ cs10*uc + cs20*uc2  + cs30*uc3 + cs40*uc4); 
    double dhGSs_dr = -hGSs*hGSs*(cs10 + 2.*cs20*uc + 3.*cs30*uc2 + 4.*cs40*uc3)*duc_dr;
    double d2hGSs_dr= -2.*hGSs*dhGSs_dr*(cs10 + 2.*cs20*uc + 3.*cs30*uc2 + 4.*cs40*uc3)*duc_dr
                      -hGS*hGS*(2.*cs20+ 6.*cs30*uc + 12.*cs40*uc2)*duc_dr*duc_dr
                      +dhGSs_dr/duc_dr*d2uc_dr;

           dGSs_dr   = dGSs0_dr*hGSs+GSs0*dhGSs_dr;

    double d2GSs_dr  = d2GSs0_dr*hGSs+2*dGSs0_dr*dhGSs_dr+GSs0*d2hGSs_dr;

    // second derivative of G = GS*S+GSs*Ss
    double d2G_dr    = d2GS_dr*S + d2GSs_dr*Sstar;    
    
    // derivatives of a,b,c coeffs
    double da_dr = 2.*dAuc2_dr*d2Auc2_dr-4*(dA*uc2+2.*A*uc*duc_dr)*dG_dr*dG_dr - 8.*A*uc2*dG_dr*d2G_dr;
    double db_dr = 2.*d2A*dAuc2_dr+2*dA*d2Auc2_dr-4.*(dA*dG_dr+2.*A*d2G_dr)*dG_dr;
    double dc_dr = 2*dA*d2A;

    // dj_dr
    double djdr_p_num = a_coeff*(-db_dr + 0.5/sqrt(Delta)*(2.*b_coeff*db_dr-4.*(da_dr*c_coeff+a_coeff*dc_dr)) )-da_dr*(-b_coeff + sqrt(Delta));
    double djdr_p     = djdr_p_num/(4.*SQ(a_coeff)*j);

    double djdr_n_num = a_coeff*(-db_dr - 0.5/sqrt(Delta)*(2.*b_coeff*db_dr-4.*(da_dr*c_coeff+a_coeff*dc_dr)) )-da_dr*(-b_coeff - sqrt(Delta));
    double djdr_n     = djdr_n_num/(4.*SQ(a_coeff)*j);

    double djdr_zero  = (b_coeff*da_dr - db_dr*a_coeff)/(4.*SQ(a_coeff)*j);

    double djdr=0;

    if (Delta==0)
      djdr = djdr_zero;
    else if (dG_dr > 0)
      djdr    = djdr_p;
    else 
      djdr = djdr_n;

    /* Compute domg/dr via Hamilton Eqs. */
    double E     = nu*H;
    double ooH   = 1./E; 
    double dHeff_orb_dr = 0.5/sqrt(Heff_orb)*(dA*(1.+j02*uc2)+A*(2.*j*djdr*uc2+2.*j02*uc*duc_dr));
    double domg_dr = -1.5*pow(r, -5./2);//-omg/Heff_orb*dHeff_orb_dr + ooH/Heff_orb*(dA*j*uc2+A*djdr*uc2+2*A*j*uc*duc_dr + dHeff_orb_dr*G + Heff_orb*dG_dr);  

    /* Compute flx */
    const double psic       = fabs((duc_dr + dG_dr*rc*sqrt(A/j02 + A*uc2)/A)/(-0.5*dA));
    const double r_omg      = pow( ((1./sqrt(rc*rc*rc*psic))+G)/E, -2./3. ); // 1./v2;
    const double v_phi      = r_omg*omg;
    const double x          = v_phi*v_phi;
    double jhat             = j/(r_omg*v_phi);
    double flx              = eob_flx_Flux_s(v_phi*v_phi, omg, r_omg, nu*H, Heff, jhat, r, 0, 0, dyn);
    
    double djdomg= djdr/domg_dr;
        
    if (DEBUG)
      printf("F=%.10f, djdomg=%.10f, omg=%.10f\n",flx, djdomg, omg);

    dy[EOB_EVOLVE_SPIN_Momg] = flx/djdomg; //flx is dj/dt

    if(dy[EOB_EVOLVE_SPIN_Momg]<0.)
       dyn->spins->omg_stop = 1e-10;

    if (DEBUG){
      /* print to file omega, 3.5PN nonspinning flux divided by newtonian prefactor and EOB flux */
      double en_flx_pn_LO = 32./5.*nu2*v10;
      double hat_en_flx   = 1. + (-1247./336. - 35./12.*nu)*v2 + 4*Pi*v3
                              + (-44711./9072. + 9271./504*nu+65./18.*nu2)*v4 + (-8191./672.-583./24*nu)*Pi*v5
                              + (6643739519./69854400. + 16./3.*Pi2 - 1712/105*EulerGamma - 856./105.*log(16*v2) + (-134543./7776. + 41./48.*Pi2)*nu - 94403./3024.*nu2 - 775./3024*nu3)*v6
                              + (-16285./504. + 214745./1728.*nu + 193385./3024.*nu2)*Pi*v7;
      double newt         = -32./5.*nu*pow(omg, 7./3.); 
      double en_flx       = en_flx_pn_LO*hat_en_flx;
      FILE *f;
      f = fopen("fluxes.txt", "a");
      fprintf(f, "%.10f %.10f %.10f\n", omg, en_flx/(nu*omg*newt), flx/newt); 
      fclose(f);
    }

  } else if (EOBPars->spin_flx == SPIN_FLX_EOB){
      /* use the exact omega22 from the EOB dynamics */
      /* do not explicitly integrate omgdot here*/
      dy[EOB_EVOLVE_SPIN_Momg] = 0.;  
  } else {
    errorexit("specify an omega dot for the spin-precession evolution");
  }

  return GSL_SUCCESS;
}

/**
 * Function: eob_spin_dyn_integrate
 * --------------------------------
 *   Integrate the spin-precession dynamics
 *   The initial data are those stored in dyn->y
 *   integration starts at dyn->t
 * 
 *   @param[out] dyn: the spin-precession dynamics
 *   @param[in]  eobdyn: the EOB dynamics
 *   @param[in]  hlm: the waveform
 * 
 *   @return GSL_SUCCESS if all goes well
*/
int eob_spin_dyn_integrate(DynamicsSpin *dyn, Dynamics *eobdyn, Waveform_lm *hlm)
{
  const int chunk = dyn->size;
  int size = chunk; // can change

  /* index of closest element to dyn->t in dyn->time 
     initial data refer to this time */
  const int i0 = find_point_bisection(dyn->t, dyn->size, dyn->time, 1);//CHECKME: 0 or 1?
  if (DEBUG) printf("initial data at index %d\n",i0);
  
  /* Set initial data */
  for (int v=0; v<EOB_EVOLVE_SPIN_NVARS; v++)
    dyn->data[v][i0] = dyn->y[v];
  
  /* GSL integrator memory */
  if (DEBUG) printf("omg_stop = %e\n", dyn->omg_stop);
  
  dyn->dt = EOBPars->spin_odes_dt; 
  const double ode_abstol = EOBPars->ode_abstol*100;
  const double ode_reltol = EOBPars->ode_reltol*100;
  double tstop;
  if (dyn->dt < 0)
    tstop = MIN(dyn->t_stop, EOBPars->ode_tmax);
  else
    tstop = MAX(dyn->t_stop, EOBPars->ode_tmax);

  gsl_odeiv2_system sys          = {p_eob_spin_dyn_rhs, NULL, EOB_EVOLVE_SPIN_NVARS, dyn};
#if (USERK45)
  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_driver * d          = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, dyn->dt, ode_abstol, ode_reltol);    
#else
  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_driver * d          = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, dyn->dt, ode_abstol, ode_reltol);    
#endif
  gsl_odeiv2_step * s            = gsl_odeiv2_step_alloc (T, EOB_EVOLVE_SPIN_NVARS);
  gsl_odeiv2_control * c         = gsl_odeiv2_control_y_new (ode_abstol, ode_reltol);
  gsl_odeiv2_evolve * e          = gsl_odeiv2_evolve_alloc (EOB_EVOLVE_SPIN_NVARS);
  
  /** Solve ODE */
  if (VERBOSE) PRSECTN("ODE Precession evolution");
  if (DEBUG) printf("t = %e\ndt = %e\ntstop = %e\n",dyn->t,dyn->dt,dyn->t_stop);
  int iter = 0;
  int status;
  int GSLSTATUS = OK;
  int breakit = 0; // flag to break while allowing for one more step
  
  /* quantities to handle backward integration */
  int eps     = 1;
  int deltap  = 0 ;
  if (dyn->dt<0){
    eps    = -1;
    deltap =  1;
  }

  /* PN/ Exact-EOB flux auxiliary flags/variables */
  double omg0eob=-1;
  double omgeobmax=100;
  double tshift=0;
  int omg_jmax =0;
  int dN       =1;
  int spin_flx = SPIN_FLX_PN; //default it to PN
  gsl_spline *omg_sp    = NULL; 
  gsl_interp_accel *acc = NULL;

  double *omega_eob = NULL;
  double *time_eob  = NULL;
  int size_eob;
  int map_from_22 = 0; //use the 22 phase rather than the orbital phase

  if(EOBPars->spin_flx==SPIN_FLX_EOB){
    /* set spin_flx*/
    spin_flx = SPIN_FLX_EOB;

    /* choose which omega_eob to use: orbital or omega22*/
    if(map_from_22){
      double *omg22_eob;
      omg22_eob = malloc ( hlm->size * sizeof(double) );
      D0_x_4(hlm->phase[1], hlm->time, hlm->size, omg22_eob);
      //D0(hlm->phase[1], hlm->time[1]-hlm->time[0], hlm->size, omg22_eob);
      for(int i =0; i < hlm->size; i++) omg22_eob[i] = omg22_eob[i]/2;
      omega_eob  = omg22_eob;
      time_eob   = hlm->time;
      size_eob   = hlm->size;
    } else {
      omega_eob  = eobdyn->data[EOB_MOMG];
      time_eob   = eobdyn->time;
      size_eob   = eobdyn->size;
    }

    /* alloc */
    omg_sp = gsl_spline_alloc (gsl_interp_cspline, size_eob);
    acc    = gsl_interp_accel_alloc ();

    /* compute spline for Momega_eob, compute Momega(t=0) and find the maximum value */
    gsl_spline_init (omg_sp, time_eob, omega_eob, size_eob);    
    
    
  }

  while (1) {
    if (VERBOSE) printf("iter %09d | t=%.f | momg= %.17f | alpha=%.6f | beta=%.6f | Lhz=%.20f\n", iter, dyn->t, dyn->y[EOB_EVOLVE_SPIN_Momg], dyn->y[EOB_EVOLVE_SPIN_alp], dyn->y[EOB_EVOLVE_SPIN_bet], dyn->y[EOB_EVOLVE_SPIN_Lz]); 
    iter++;

    //GSLSTATUS = gsl_odeiv2_evolve_apply_fixed_step (e, c, s, &sys, &dyn->t, dyn->dt, dyn->y);//uniform
    GSLSTATUS = gsl_odeiv2_evolve_apply (e, c, s, &sys, &dyn->t, tstop, &dyn->dt, dyn->y);

    /** Check for failures ... */
    if (GSLSTATUS != GSL_SUCCESS) {
      printf("GSL Error = %d", GSLSTATUS);
      return ERROR_ODEINT;
    }
    if (dyn->y[EOB_EVOLVE_SPIN_Momg] < 0.){
      if (VERBOSE) printf("Stop: Momg < 0.\n");
      break;
    }
    if (EOBPars->spin_flx != SPIN_FLX_EOB && dyn->dt > 0 && (dyn->y[EOB_EVOLVE_SPIN_Momg] - dyn->data[EOB_EVOLVE_SPIN_Momg][i0+iter-1] < 1e-15) ){
      iter--; //don't count this iteration
      if (VERBOSE) printf("Stop: Forward integration, Momg reached a plateau.\n");
      break;
    }
    if (EOBPars->spin_flx == SPIN_FLX_EOB && dyn->dt > 0 && (dyn->y[EOB_EVOLVE_SPIN_Momg] > omgeobmax) ){
      iter--; //don't count this iteration
      if (VERBOSE) printf("Stop: Forward integration, EOB flux, Momg > max(Momg_eob).\n");
      break;
    }
    /* Update alpha and beta angles */
    dyn->y[EOB_EVOLVE_SPIN_alp]  = eob_spin_dyn_alpha(dyn->y[EOB_EVOLVE_SPIN_Lx],
			          dyn->y[EOB_EVOLVE_SPIN_Ly],
			          dyn->y[EOB_EVOLVE_SPIN_Lz]) + deltap*Pi;
    dyn->y[EOB_EVOLVE_SPIN_bet]  = eps*eob_spin_dyn_beta(dyn->y[EOB_EVOLVE_SPIN_Lx],
						    dyn->y[EOB_EVOLVE_SPIN_Ly],
						    dyn->y[EOB_EVOLVE_SPIN_Lz]);
  

    if(EOBPars->spin_flx==SPIN_FLX_EOB && dyn->t-tshift > time_eob[size_eob-1]){
      iter--;
      if (VERBOSE) printf("Stop: Forward integration, EOB flux, t_spindyn > max(t_eob).\n");
      break;
    }

    if(spin_flx==SPIN_FLX_EOB && (dyn->y[EOB_EVOLVE_SPIN_Momg] >= omega_eob[0])){
      EOBPars->spin_flx = spin_flx;
      dyn->y[EOB_EVOLVE_SPIN_Momg] = gsl_spline_eval(omg_sp, dyn->t, acc); //evaluate spline
    }

    /* Update size and push arrays (if needed) */
    if (iter==size) {
      size += chunk;
      dyn->size = size;
      DynamicsSpin_push (&dyn, size);
    }
    
    dyn->time[iter]  = dyn->t; 
    for (int v=0; v<EOB_EVOLVE_SPIN_NVARS; v++)
      dyn->data[v][i0+iter]  = dyn->y[v];   

    /* Stop integration */
    if (breakit)
      break;
    if ((dyn->omg_stop>0.) && (dyn->y[EOB_EVOLVE_SPIN_Momg] > dyn->omg_stop) && dyn->dt > 0.){
      if (VERBOSE) printf("Stop: Forward integration, Momg > Momg_stop.\n");
      break;
    }
    if (iter>1 && dyn->dt > 0 && EOBPars->model == MODEL_GIOTTO && dyn->y[EOB_EVOLVE_SPIN_Momg] < dyn->data[EOB_EVOLVE_SPIN_Momg][i0+iter-1]){
      iter--;
      if (VERBOSE) printf("Stop: Forward integration, Reached maximum of Momg.\n");
      break;
    }
    if ((dyn->t_stop>0.) && (dyn->t > dyn->t_stop)){
      if (VERBOSE) printf("Stop: t_spindyn > t_stop.\n");
      break;
    }
    if (dyn->dt < 0. && (dyn->y[EOB_EVOLVE_SPIN_Momg] < dyn->omg_stop)){
      breakit = 1; //break at next iteration
      if (VERBOSE) printf("Backward integration. Momg < Momg_stop. Stop at next iteration.\n");
    }
  } /* end time iteration */

  /* Resize to actual size */
  size = i0 + iter + 1;
  DynamicsSpin_push (&dyn, size);
  
  gsl_spline_free(omg_sp);
  gsl_interp_accel_free(acc);
  /* Free ODE system solver */
  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
  gsl_odeiv2_driver_free (d);

  if(map_from_22)
    free(omega_eob);

  return OK;
}

/**
  * Function: eob_spin_dyn_integrate_backwards
  * ------------------------------------------
  *   Integrate the spin-precession dynamics backwards in time
  *   to omg0  
  *   The initial data are those stored in dyn->y
  *   integration starts at dyn->t
  * 
  *   @param[out] dyn: the spin-precession dynamics
  *   @param[in]  eobdyn: the EOB dynamics
  *   @param[in]  hlm: the waveform
  *   @param[in]  omg0: the frequency to which we want to integrate
  * 
  *   @return GSL_SUCCESS if all goes well
*/
int eob_spin_dyn_integrate_backwards(DynamicsSpin *dyn, Dynamics *eobdyn, Waveform_lm *hlm, double omg0)
{
  // Alloc new dynamics
  DynamicsSpin *spindyn_tmp = NULL;
  DynamicsSpin_alloc(&spindyn_tmp, 10);
  EOBPars->spin_odes_dt = -EOBPars->spin_odes_dt/10;
  spindyn_tmp->omg_stop = 0.99*omg0; // slightly below omg0, to avoid interpolation issues at the edges
  spindyn_tmp->t_stop   = -1e15;
  //eob_spin_dyn(spindyn_tmp, EOBPars->initial_frequency/time_units_factor(EOBPars->M));
  int Nint = 0;
  dyn->t = 0.;
  double m1 = nu_to_X1(EOBPars->nu);
  double m2 = 1 - m1;
  const double M12 = SQ(m1);
  const double M22 = SQ(m2);  
  spindyn_tmp->y[EOB_EVOLVE_SPIN_SxA] = dyn->data[EOB_EVOLVE_SPIN_SxA][Nint];
  spindyn_tmp->y[EOB_EVOLVE_SPIN_SyA] = dyn->data[EOB_EVOLVE_SPIN_SyA][Nint];
  spindyn_tmp->y[EOB_EVOLVE_SPIN_SzA] = dyn->data[EOB_EVOLVE_SPIN_SzA][Nint];
  spindyn_tmp->y[EOB_EVOLVE_SPIN_SxB] = dyn->data[EOB_EVOLVE_SPIN_SxB][Nint];
  spindyn_tmp->y[EOB_EVOLVE_SPIN_SyB] = dyn->data[EOB_EVOLVE_SPIN_SyB][Nint];
  spindyn_tmp->y[EOB_EVOLVE_SPIN_SzB] = dyn->data[EOB_EVOLVE_SPIN_SzB][Nint];
  spindyn_tmp->y[EOB_EVOLVE_SPIN_Lx]  = dyn->data[EOB_EVOLVE_SPIN_Lx][Nint]; 
  spindyn_tmp->y[EOB_EVOLVE_SPIN_Ly]  = dyn->data[EOB_EVOLVE_SPIN_Ly][Nint];
  spindyn_tmp->y[EOB_EVOLVE_SPIN_Lz]  = dyn->data[EOB_EVOLVE_SPIN_Lz][Nint];
  spindyn_tmp->y[EOB_EVOLVE_SPIN_alp] = dyn->data[EOB_EVOLVE_SPIN_alp][Nint];
  spindyn_tmp->y[EOB_EVOLVE_SPIN_bet] = dyn->data[EOB_EVOLVE_SPIN_bet][Nint];
  spindyn_tmp->y[EOB_EVOLVE_SPIN_gam] = dyn->data[EOB_EVOLVE_SPIN_gam][Nint]; 
  spindyn_tmp->y[EOB_EVOLVE_SPIN_Momg] = dyn->data[EOB_EVOLVE_SPIN_Momg][Nint];
  
  for (int v=0; v<EOB_EVOLVE_SPIN_NVARS; v++)
    spindyn_tmp->data[v][0] = spindyn_tmp->y[v];

  eob_spin_dyn_integrate(spindyn_tmp, eobdyn, hlm);

  int size = spindyn_tmp->size;

  /*shift by t0*/
  for(int i=0; i<size;i++)
    spindyn_tmp->time[i] = spindyn_tmp->time[i] + dyn->time[Nint];

  /* rearrange variables*/
  double tmp;
  for(int v=0; v<EOB_EVOLVE_SPIN_NVARS; v++){
    for(int i=0; i<floor(size/2);i++){
      tmp = spindyn_tmp->data[v][i];
      spindyn_tmp->data[v][i] = spindyn_tmp->data[v][size-1-i];
      spindyn_tmp->data[v][size-1-i] = tmp;
    }
  }
  for(int i=0; i<floor(size/2);i++){
    tmp = spindyn_tmp->time[i];
    spindyn_tmp->time[i] = spindyn_tmp->time[size-1-i];
    spindyn_tmp->time[size-1-i] = tmp;
  }
  
  /*join at t0*/
  DynamicsSpin_join(spindyn_tmp, dyn, dyn->time[Nint]);

  /*copy in dyn*/
  DynamicsSpin_push(&dyn, spindyn_tmp->size);
  for (int i = 0; i < spindyn_tmp->size; i++) 
    dyn->time[i] = spindyn_tmp->time[i] - spindyn_tmp->time[0]; //start dynamics at t=0
  
  dyn->time_backward = - spindyn_tmp->time[0];             // keep track of when forward integration starts (time)
  dyn->omg_backward  = dyn->data[EOB_EVOLVE_SPIN_Momg][0]; // keep track of when forward integration starts (omega)

  for (int v = 0; v < EOB_EVOLVE_SPIN_NVARS; v++) {
    for (int i = 0; i < spindyn_tmp->size; i++) {
      dyn->data[v][i] = spindyn_tmp->data[v][i];
    }
  }
  if (EOBPars->spin_interp_domain == 0){
    for(int v=0; v < EOB_EVOLVE_SPIN_NVARS; v++)
      gsl_spline_init (dyn->spline[v], dyn->time, dyn->data[v], dyn->size);  
  } else {
    for(int v=0; v < EOB_EVOLVE_SPIN_NVARS; v++)
      gsl_spline_init (dyn->spline[v], dyn->data[EOB_EVOLVE_SPIN_Momg], dyn->data[v], dyn->size); 
  }
  /*free*/ 
  DynamicsSpin_free(spindyn_tmp);
  return OK;
}

/**
  * Function: eob_spin_dyn
  * ----------------------
  *   Main driver routine for the spin-precession dynamics
  *   The initial data are those stored in dyn->y
  *   integration starts at dyn->t
  * 
  *   @param[out] dyn: the spin-precession dynamics
  *   @param[in]  eobdyn: the EOB dynamics
  *   @param[in]  hlm: the waveform
  *   @param[in]  omg0: the initial frequency
  * 
  *   @return GSL_SUCCESS if all goes well
*/
int eob_spin_dyn(DynamicsSpin *dyn, Dynamics *eobdyn, Waveform_lm *hlm, double omg0)
{
  const int chunk = dyn->size;

  //FIXME add option in EOBpars and add PN_abc and EOB rhs (when coded)
  p_eob_spin_dyn_rhs = eob_spin_dyn_rhs_PN;
  
  /* Set the stopping frequency as the NR merger (if not set) */
  if (EOBPars->spin_odes_omg_stop<0) {
    const double fact = 1.1; // need to go slightly above for ringdown attachment
    dyn->omg_stop = fact * eob_mrg_momg(EOBPars->nu, EOBPars->X1, EOBPars->X2, EOBPars->chi1, EOBPars->chi2);
    EOBPars->spin_odes_omg_stop = dyn->omg_stop;
  }
  dyn->t_stop = EOBPars->spin_odes_t_stop;
  
  /** Initial data */
  dyn->t = 0.;
  double m1 = nu_to_X1(EOBPars->nu);
  double m2 = 1 - m1;
  const double M12 = SQ(m1);
  const double M22 = SQ(m2);  
  dyn->y[EOB_EVOLVE_SPIN_SxA] = EOBPars->chi1x*M12; 
  dyn->y[EOB_EVOLVE_SPIN_SyA] = EOBPars->chi1y*M12;
  dyn->y[EOB_EVOLVE_SPIN_SzA] = EOBPars->chi1z*M12;
  dyn->y[EOB_EVOLVE_SPIN_SxB] = EOBPars->chi2x*M22;
  dyn->y[EOB_EVOLVE_SPIN_SyB] = EOBPars->chi2y*M22;
  dyn->y[EOB_EVOLVE_SPIN_SzB] = EOBPars->chi2z*M22;
  dyn->y[EOB_EVOLVE_SPIN_Lx]  = 0; 
  dyn->y[EOB_EVOLVE_SPIN_Ly]  = 0;
  dyn->y[EOB_EVOLVE_SPIN_Lz]  = 1.;
  dyn->y[EOB_EVOLVE_SPIN_alp] = alpha_initial_condition(EOBPars); 
  dyn->y[EOB_EVOLVE_SPIN_bet] = eob_spin_dyn_beta(dyn->y[EOB_EVOLVE_SPIN_Lx],
						  dyn->y[EOB_EVOLVE_SPIN_Ly],
						  dyn->y[EOB_EVOLVE_SPIN_Lz]);
  dyn->y[EOB_EVOLVE_SPIN_gam] = alpha_initial_condition(EOBPars); 
  
  double time_unit_fact = 1;
  if(!EOBPars->use_geometric_units)
    time_unit_fact = time_units_factor(EOBPars->M);

  dyn->y[EOB_EVOLVE_SPIN_Momg] = omg0;
  
  for (int v=0; v<EOB_EVOLVE_SPIN_NVARS; v++)
    dyn->data[v][0] = dyn->y[v];
  
  /** Integrate ODEs */
  return eob_spin_dyn_integrate(dyn, eobdyn, hlm);
  
}

/**
  * Function: eob_spin_dyn_abc_interp
  * ---------------------------------
  *   Helper routine to interpolate Euler angles at given time
  * 
  *   @param[in]  dyn: the spin-precession dynamics
  *   @param[in]  eobdyn: the EOB dynamics
  *   @param[in]  hlm: the waveform
  *   @param[in]  time: the time at which we want to interpolate
  *   @param[out] alpha_p: the alpha angle
  *   @param[out] beta_p: the beta angle
  *   @param[out] gamma_p: the gamma angle
  *   @param[in]  continue_integration: flag to continue integration if time > tmax
  * 
  *   @return GSL_SUCCESS if all goes well
*/
void eob_spin_dyn_abc_interp(DynamicsSpin *dyn, Dynamics *eobdyn, Waveform_lm *hlm, double time,
			     double *alpha_p, double *beta_p, double *gamma_p,
			     int continue_integration)
{
  const int smax = dyn->size-1;
  double alpha, beta, gamma;
  int interp = 1;
  
  if (time >= dyn->time[smax]) {
    
    /* The precessing dynamics is too short ! 
       We have two options 
       1. continue the integration to time we need
       2. return the last angle
     */
    if (continue_integration) {      

      dyn->omg_stop = -1; // use the tstop not Momg_stop
      dyn->t = dyn->time[smax];
      dyn->t_stop = time + dyn->dt;
      eob_spin_dyn_integrate(dyn, eobdyn, hlm);      

    } else {

      alpha = dyn->data[EOB_EVOLVE_SPIN_alp][smax];
      beta = dyn->data[EOB_EVOLVE_SPIN_bet][smax];
      gamma = dyn->data[EOB_EVOLVE_SPIN_gam][smax];

      interp = 0;// skip interpolation below      
      if (VERBOSE) errorexit("Spin dynamics too short to interp!\n");
    }
  }
  
  /* Interp */
  if (interp) {
    alpha = interp_spline_pt(dyn->time, dyn->data[EOB_EVOLVE_SPIN_alp], dyn->size, time);
    beta  = interp_spline_pt(dyn->time, dyn->data[EOB_EVOLVE_SPIN_bet], dyn->size, time);
    gamma = interp_spline_pt(dyn->time, dyn->data[EOB_EVOLVE_SPIN_gam], dyn->size, time);
  }
  
  *alpha_p = alpha;
  *beta_p = beta;
  *gamma_p = gamma;
}

/**
  * Function: eob_spin_dyn_Sproj_interp
  * -----------------------------------
  *   Helper routine to interpolate spins parallel/perpendicular to orbital ang. mom. 
  *   at given time or frequency
  * 
  *   @param[in]   dyn: the spin-precession dynamics
  *   @param[in]   var: the time or the frequency at which we want to interpolate
  *   @param[out]  SApara: the spin of body A parallel to orbital ang. mom.
  *   @param[out]  SBpara: the spin of body B parallel to orbital ang. mom.
  *   @param[out]  Spara: the total spin parallel to orbital ang. mom.
  *   @param[out]  SAperp: the spin of body A perpendicular to orbital ang. mom.
  *   @param[out]  SBperp: the spin of body B perpendicular to orbital ang. mom.
  *   @param[out]  Sperp: the total spin perpendicular to orbital ang. mom.
  *   @param[in]   interp_domain: flag to interpolate in time or in frequency
  * 
  *   @return GSL_SUCCESS if all goes well
*/
void eob_spin_dyn_Sproj_interp(DynamicsSpin *dyn, double var,
			       double *SApara, double *SBpara, double *Spara,
			       double *SAperp, double *SBperp, double *Sperp, // these are 3-vectors
			       int interp_domain)
{
  const int smax = dyn->size-1;
  double SA[IN3], SB[IN3], Lh[IN3];
  int interp = 1;
  //printf("spindyn_t0=%.10f, spindyn tmax = %.10f\n", dyn->time[0], dyn->time[smax]);
  if(interp_domain == 0){ //interpolate in time
    double time = var;
    if (time >= dyn->time[smax]){

      SA[Ix] = dyn->data[EOB_EVOLVE_SPIN_SxA][smax];
      SA[Iy] = dyn->data[EOB_EVOLVE_SPIN_SyA][smax];
      SA[Iz] = dyn->data[EOB_EVOLVE_SPIN_SzA][smax];
      
      SB[Ix] = dyn->data[EOB_EVOLVE_SPIN_SxB][smax];
      SB[Iy] = dyn->data[EOB_EVOLVE_SPIN_SyB][smax];
      SB[Iz] = dyn->data[EOB_EVOLVE_SPIN_SzB][smax];
      
      Lh[Ix] = dyn->data[EOB_EVOLVE_SPIN_Lx][smax];
      Lh[Iy] = dyn->data[EOB_EVOLVE_SPIN_Ly][smax];
      Lh[Iz] = dyn->data[EOB_EVOLVE_SPIN_Lz][smax];

      interp = 0;
    }
    if(interp){
      SA[Ix] = gsl_spline_eval(dyn->spline[EOB_EVOLVE_SPIN_SxA], time, dyn->accel[EOB_EVOLVE_SPIN_SxA]);
      SA[Iy] = gsl_spline_eval(dyn->spline[EOB_EVOLVE_SPIN_SyA], time, dyn->accel[EOB_EVOLVE_SPIN_SyA]);
      SA[Iz] = gsl_spline_eval(dyn->spline[EOB_EVOLVE_SPIN_SzA], time, dyn->accel[EOB_EVOLVE_SPIN_SzA]);

      SB[Ix] = gsl_spline_eval(dyn->spline[EOB_EVOLVE_SPIN_SxB], time, dyn->accel[EOB_EVOLVE_SPIN_SxB]);
      SB[Iy] = gsl_spline_eval(dyn->spline[EOB_EVOLVE_SPIN_SyB], time, dyn->accel[EOB_EVOLVE_SPIN_SyB]);
      SB[Iz] = gsl_spline_eval(dyn->spline[EOB_EVOLVE_SPIN_SzB], time, dyn->accel[EOB_EVOLVE_SPIN_SzB]);

      Lh[Ix] = gsl_spline_eval(dyn->spline[EOB_EVOLVE_SPIN_Lx], time, dyn->accel[EOB_EVOLVE_SPIN_Lx]);
      Lh[Iy] = gsl_spline_eval(dyn->spline[EOB_EVOLVE_SPIN_Ly], time, dyn->accel[EOB_EVOLVE_SPIN_Ly]);
      Lh[Iz] = gsl_spline_eval(dyn->spline[EOB_EVOLVE_SPIN_Lz], time, dyn->accel[EOB_EVOLVE_SPIN_Lz]);
    }

  } else {// interpolate in omega
    double omg = var;
    if (omg >= dyn->data[EOB_EVOLVE_SPIN_Momg][smax]){
      SA[Ix] = dyn->data[EOB_EVOLVE_SPIN_SxA][smax];
      SA[Iy] = dyn->data[EOB_EVOLVE_SPIN_SyA][smax];
      SA[Iz] = dyn->data[EOB_EVOLVE_SPIN_SzA][smax];
      
      SB[Ix] = dyn->data[EOB_EVOLVE_SPIN_SxB][smax];
      SB[Iy] = dyn->data[EOB_EVOLVE_SPIN_SyB][smax];
      SB[Iz] = dyn->data[EOB_EVOLVE_SPIN_SzB][smax];
      
      Lh[Ix] = dyn->data[EOB_EVOLVE_SPIN_Lx][smax];
      Lh[Iy] = dyn->data[EOB_EVOLVE_SPIN_Ly][smax];
      Lh[Iz] = dyn->data[EOB_EVOLVE_SPIN_Lz][smax];

      interp = 0;
    } else if (omg <= dyn->data[EOB_EVOLVE_SPIN_Momg][0]){
      //this can happen if the initial EOB and PN omg0 are very close
      SA[Ix] = dyn->data[EOB_EVOLVE_SPIN_SxA][0];
      SA[Iy] = dyn->data[EOB_EVOLVE_SPIN_SyA][0];
      SA[Iz] = dyn->data[EOB_EVOLVE_SPIN_SzA][0];
      
      SB[Ix] = dyn->data[EOB_EVOLVE_SPIN_SxB][0];
      SB[Iy] = dyn->data[EOB_EVOLVE_SPIN_SyB][0];
      SB[Iz] = dyn->data[EOB_EVOLVE_SPIN_SzB][0];
      
      Lh[Ix] = dyn->data[EOB_EVOLVE_SPIN_Lx][0];
      Lh[Iy] = dyn->data[EOB_EVOLVE_SPIN_Ly][0];
      Lh[Iz] = dyn->data[EOB_EVOLVE_SPIN_Lz][0];  
      interp = 0;  
    }
    if(interp){
      SA[Ix] = gsl_spline_eval(dyn->spline[EOB_EVOLVE_SPIN_SxA], omg, dyn->accel[EOB_EVOLVE_SPIN_SxA]);
      SA[Iy] = gsl_spline_eval(dyn->spline[EOB_EVOLVE_SPIN_SyA], omg, dyn->accel[EOB_EVOLVE_SPIN_SyA]);
      SA[Iz] = gsl_spline_eval(dyn->spline[EOB_EVOLVE_SPIN_SzA], omg, dyn->accel[EOB_EVOLVE_SPIN_SzA]);

      SB[Ix] = gsl_spline_eval(dyn->spline[EOB_EVOLVE_SPIN_SxB], omg, dyn->accel[EOB_EVOLVE_SPIN_SxB]);
      SB[Iy] = gsl_spline_eval(dyn->spline[EOB_EVOLVE_SPIN_SyB], omg, dyn->accel[EOB_EVOLVE_SPIN_SyB]);
      SB[Iz] = gsl_spline_eval(dyn->spline[EOB_EVOLVE_SPIN_SzB], omg, dyn->accel[EOB_EVOLVE_SPIN_SzB]);

      Lh[Ix] = gsl_spline_eval(dyn->spline[EOB_EVOLVE_SPIN_Lx], omg, dyn->accel[EOB_EVOLVE_SPIN_Lx]);
      Lh[Iy] = gsl_spline_eval(dyn->spline[EOB_EVOLVE_SPIN_Ly], omg, dyn->accel[EOB_EVOLVE_SPIN_Ly]);
      Lh[Iz] = gsl_spline_eval(dyn->spline[EOB_EVOLVE_SPIN_Lz], omg, dyn->accel[EOB_EVOLVE_SPIN_Lz]);
    }
  }

  /* Total spin */
  //double SAB[IN3], S2;
  //for (int i=0; i<IN3; i++) SAB[i] = SA[i]+SB[i];
  //vect_dot3(SAB, SAB, &S2);

  /* direction of Lh */
  //FIXME: this should not be needed
  double n[IN3], normL2;
  vect_dot3(Lh, Lh, &normL2); 
  const double oonormL = 1./sqrt(normL2); 
  for (int i=0; i<IN3; i++) n[i] = Lh[i]*oonormL;
  
  /* Projections */
  double LSA, LSB, LS;    
  if (SApara) {
    vect_dot3(SA, n, &LSA);
    *SApara = LSA;
  }
  if (SBpara) {
    vect_dot3(SB, n, &LSB);
    *SBpara = LSB;
  }
  if (Spara) {
    *Spara = *SApara + *SBpara;
  }
  if (SAperp) {
    for (int i=0; i<IN3; i++) SAperp[i] = SA[i] - LS * n[i];
  }
  if (SBperp) {
    for (int i=0; i<IN3; i++) SBperp[i] = SB[i] - LS * n[i];
  }
  if (Sperp) {
    for (int i=0; i<IN3; i++) Sperp[i] = SAperp[i] - SBperp[i];
  }

}


