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

/** r.h.s. of EOB Hamiltonian dynamics, no spins version */ 
int eob_dyn_rhs(double t, const double y[], double dy[], void *d)
{
  
  (void)(t); /* avoid unused parameter warning */
  Dynamics *dyn = d;  

  const double nu = dyn->nu;  
  const double z3 = 2.0*nu*(4.0-3.0*nu);

  /** Unpack y */ 
  const double phi    = y[EOB_EVOLVE_PHI];
  const double r      = y[EOB_EVOLVE_RAD];
  const double pphi   = y[EOB_EVOLVE_PPHI];
  const double prstar = y[EOB_EVOLVE_PRSTAR];
 
  /** Compute EOB Metric */
  double A, B, dA, d2A, dB, pl_hold;
  eob_metric(r, d, &A, &B, &dA, &d2A, &dB, &pl_hold);

  /** Compute Hamiltonian */
  double H, Heff, dHeff_dr,dHeff_dprstar;
  eob_ham(nu, r,pphi,prstar,A,dA, &H,&Heff,&dHeff_dr,&dHeff_dprstar,NULL);
  double E = nu*H;

  /** Shorthands */
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

  /** d\phi/dt */
  dy[EOB_EVOLVE_PHI] = Omega;
  
  /** dr/dt (conservative part of) */
  dy[EOB_EVOLVE_RAD] = sqrtAbyB*(prstar+4.0*nu*(4.0-3.0*nu)*A*u2*prstar3)*divHE;
  
  /** dp_{r*}/dt (conservative part of) */
  dy[EOB_EVOLVE_PRSTAR] = - 0.5*sqrtAbyB*( pphi2*u2*(dA-2.0*A*u) + dA + 2.0*nu*(4.0-3.0*nu)*(dA*u2 - 2.0*A*u3)*prstar4 )*divHE;        
  
  /** Compute flux */
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
  const double ddotr_dr      = sqrtAbyB*( (prstar + z3*2.*A*u2*prstar3)*(0.5*(dA/A-dB/B)-dHeff_dr*tmpE)+ 2.0*z3*(dA*u2 - 2.*A*u3)*prstar3)*divHE;
  const double ddotr_dprstar = sqrtAbyB*( 1.+z3*6.*A*u2*prstar2-(prstar + z3*2.*A*u2*prstar3)*dHeff_dprstar*tmpE)*divHE;
  
  /* Approximate ddot(r) without Flux */
  const double ddotr = dprstar_dt*ddotr_dprstar + dr_dt*ddotr_dr;
  
  /** Compute flux and dp_{\phi}/dt */
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

/* EOB nonspinning Hamiltonian */
void eob_ham(double nu, double r, double pphi, double prstar, double A, double dA,
	     double *H, /* real EOB Hamiltonian divided by mu=m1m2/(m1+m2) */
	     double *Heff, /* effective EOB Hamiltonian (divided by mu) */
	     double *dHeff_dr, /* drvt Heff,r */
	     double *dHeff_dprstar, /* drvt Heff,prstar */
	     double *dHeff_dpphi /* drvt Heff,pphi */
	     )
{
  const double z3 = 2.0*nu*(4.0-3.0*nu);
  const double pphi2 = SQ(pphi);
  const double u = 1./r;
  const double u2 = SQ(u);
  const double u3 = u2*u;
  const double prstar2 = SQ(prstar);
  const double prstar3 = prstar2*prstar;
  const double prstar4 = prstar2*prstar2;

  *Heff = sqrt(A*(1.0 + pphi2*u2) + prstar2 + z3*A*u2*prstar4);
  *H    = sqrt( 1.0 + 2.0*nu*(*Heff - 1) )/nu;  
   
  if (dHeff_dr != NULL)      *dHeff_dr      = 0.5*(dA + (pphi2 + z3*prstar4)*(dA*u2 - 2*A*u3))/(*Heff);
  if (dHeff_dprstar != NULL) *dHeff_dprstar = (prstar + z3*2.0*A*u2*prstar3)/(*Heff);
  if (dHeff_dpphi != NULL)   *dHeff_dpphi   = A*pphi*u2/(*Heff);
}

/** r.h.s. of EOB Hamiltonian dynamics, spins version */
int eob_dyn_rhs_s(double t, const double y[], double dy[], void *d)
{
  
  (void)(t); /* avoid unused parameter warning */
  Dynamics *dyn = d;
  
  /* Unpack values */
  const double nu    = dyn->nu;
  const double S     = dyn->S;
  const double Sstar = dyn->Sstar;
  const double chi1  = dyn->chi1;
  const double chi2  = dyn->chi2;
  const double X1    = dyn->X1;
  const double X2    = dyn->X2;
  const double c3    = dyn->cN3LO;
  const double aK2   = dyn->aK2;
  const double a1    = dyn->a1;
  const double a2    = dyn->a2;
  const double C_Q1  = dyn->C_Q1;
  const double C_Q2  = dyn->C_Q2;
  const double C_Oct1 = dyn->C_Oct1;
  const double C_Oct2 = dyn->C_Oct2;
  const double C_Hex1 = dyn->C_Hex1;
  const double C_Hex2 = dyn->C_Hex2;
  const int usetidal = dyn->use_tidal;
  const int usespins = dyn->use_spins;
  
  /* Shorthands */
  const double r      = y[EOB_EVOLVE_RAD];
  const double prstar = y[EOB_EVOLVE_PRSTAR];
  const double pphi   = y[EOB_EVOLVE_PPHI];
  const double pphi2  = pphi*pphi;
  
  /** Compute Metric */
  double A, B, dA, d2A, dB, pl_hold;
  eob_metric_s(r, d, &A, &B, &dA, &d2A, &dB, &pl_hold);
  
  /* Compute centrifugal radius */
  double rc, drc_dr, d2rc_dr;
  eob_dyn_s_get_rc(r, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc, &drc_dr, &d2rc_dr);
  const double uc     = 1./rc;
  const double uc2    = uc*uc;
  const double uc3    = uc2*uc;
  
  /* Compute Hamiltonian */
  double Heff_orb, Heff, H, dHeff_dr, dHeff_dprstar, d2Heff_dprstar20, dHeff_dpphi;
  eob_ham_s(nu, r, rc, drc_dr, pphi, prstar, S, Sstar, chi1, chi2, X1, X2, aK2, c3, A, dA, 
	    &H, &Heff, &Heff_orb, &dHeff_dr, &dHeff_dprstar, &dHeff_dpphi, &d2Heff_dprstar20);
  
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
  const double psic       = fabs((duc_dr + dGtilde_dr*rc*sqrt(A/pphi2 + A*uc2)/A)/(-0.5*dA));
  // FIXME: Different from Matlab code.
  //        Added absolute value to avoid NaN
  const double r_omg      = pow( ((1./sqrt(rc*rc*rc*psic))+Gtilde)*ooH0, -2./3. );
  const double v_phi      = r_omg*Omg;
  const double x          = v_phi*v_phi;
  const double jhat       = pphi/(r_omg*v_phi);
      
  /** Compute flux and dp_{\phi}/dt */
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

/** r.h.s. of EOB Hamiltonian dynamics, eccentric version */
int eob_dyn_rhs_ecc(double t, const double y[], double dy[], void *d)
{

  (void)(t); /* avoid unused parameter warning */
  Dynamics *dyn = d;
  
  /* Unpack values */
  const double nu    = dyn->nu;
  const double S     = dyn->S;
  const double Sstar = dyn->Sstar;
  const double chi1  = dyn->chi1;
  const double chi2  = dyn->chi2;
  const double X1    = dyn->X1;
  const double X2    = dyn->X2;
  const double c3    = dyn->cN3LO;
  const double aK2   = dyn->aK2;
  const double a1    = dyn->a1;
  const double a2    = dyn->a2;
  const double C_Q1  = dyn->C_Q1;
  const double C_Q2  = dyn->C_Q2;
  const double C_Oct1 = dyn->C_Oct1;
  const double C_Oct2 = dyn->C_Oct2;
  const double C_Hex1 = dyn->C_Hex1;
  const double C_Hex2 = dyn->C_Hex2;
  const int usetidal = dyn->use_tidal;
  const int usespins = dyn->use_spins;
  
  /* Shorthands */
  const double r      = y[EOB_EVOLVE_RAD];
  const double prstar = y[EOB_EVOLVE_PRSTAR];
  const double pphi   = y[EOB_EVOLVE_PPHI];
  
  const double pphi2  = pphi*pphi;
  
  /** Compute Metric, centrifugal radius and Hamiltonan */
  double A, B, dA, d2A, dB, pl_hold;
  double rc, drc_dr, d2rc_dr2;
  double Heff_orb, Heff, H, dHeff_dr, dHeff_dprstar, d2Heff_dprstar20, dHeff_dpphi;

  if (usespins) {
    eob_metric_s(r, d, &A, &B, &dA, &d2A, &dB, &pl_hold);

    eob_dyn_s_get_rc(r, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc, &drc_dr, &d2rc_dr2);

    eob_ham_s(nu, r, rc, drc_dr, pphi, prstar, S, Sstar, chi1, chi2, X1, X2, aK2, c3, A, dA, 
	      &H, &Heff, &Heff_orb, &dHeff_dr, &dHeff_dprstar, &dHeff_dpphi, &d2Heff_dprstar20);
  } else {
    eob_metric(r, d, &A, &B, &dA, &d2A, &dB, &pl_hold);

    eob_ham(nu, r, pphi, prstar, A, dA, 
	    &H, &Heff_orb, &dHeff_dr, &dHeff_dprstar, &dHeff_dpphi);
    
    Heff = Heff_orb;
    d2Heff_dprstar20 = 1/Heff_orb;

    rc = r;
    drc_dr = 1.;
    d2rc_dr2 = 0.;
  }
  
  double uc     = 1./rc;
  double uc2    = uc*uc;
  double uc3    = uc2*uc;
  double E      = nu*H;
  double ooH    = 1./E;
    
  const double sqrtAbyB       = sqrt(A/B);
  const double dp_rstar_dt_0  = - sqrtAbyB*dHeff_dr*ooH;
  const double ddotr_dp_rstar = sqrtAbyB*d2Heff_dprstar20*ooH;
  const double Omg            = dHeff_dpphi*ooH;
  const double ddotr          = dp_rstar_dt_0*ddotr_dp_rstar; /* approximate ddot(r)_0 without Fphi, order pr_star^2 neglected */
  
  /* r evol eqn rhs */
  dy[EOB_EVOLVE_RAD] = sqrtAbyB*dHeff_dprstar*ooH;
  
  /* phi evol eqn rhs */
  dy[EOB_EVOLVE_PHI] = Omg;
    
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
  psic       = (duc_dr + dGtilde_dr*rc*sqrt(A/pphi2 + A*uc2)/A)/(-0.5*dA);
  r_omg      = pow( ((1./sqrt(rc*rc*rc*psic))+Gtilde)*ooH0, -2./3. );

  } else {
    double V = A*(1. + SQ(pphi/r));
    psic  = 2.*(1. + 2.*nu*(sqrt(V) - 1.))/(SQ(r)*dA);
    r_omg = r*pow(psic,1./3.);
  }
  const double v_phi      = r_omg*Omg;
  const double x          = v_phi*v_phi;
  const double jhat       = pphi/(r_omg*v_phi);
  
  /** Compute flux and dp_{\phi}/dt */
  double Fphi = eob_flx_Flux_s(x,Omg,r_omg,E,Heff,jhat,r,prstar,ddotr,dyn);
  double Fphi_NC = eob_flx_Fphi_ecc(x, Omg, r_omg, jhat, dy[EOB_EVOLVE_RAD], dyn);
  if (dyn->noflx) dy[EOB_EVOLVE_PPHI] = 0.;
  else            dy[EOB_EVOLVE_PPHI] = Fphi*Fphi_NC;

  double Fr = eob_flx_Fr_ecc(r, prstar, pphi, dyn);
  double Frstar = sqrtAbyB*Fr;
  
  /* dp_{r*}/dt */
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

/* EOB spinning Hamiltonian */
void eob_ham_s(double nu,
               double r,
               double rc,
               double drc_dr,
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
               double *H,             /* real EOB Hamiltonian divided by mu=m1m2/(m1+m2) */
               double *Heff,          /* effective EOB Hamiltonian (divided by mu) */
               double *Heff_orb,
               double *dHeff_dr,      /* drvt Heff,r */
               double *dHeff_dprstar, /* drvt Heff,prstar */
               double *dHeff_dpphi,    /* drvt Heff,pphi */
               double *d2Heff_dprstar20
               )
{
    /* Shorthands */
    const double z3      = 2.0*nu*(4.0-3.0*nu);
    const double pphi2    = SQ(pphi);
    const double prstar2 = SQ(prstar);
    const double prstar3 = prstar2*prstar;
    const double prstar4 = prstar2*prstar2;
    const double uc  = 1./rc;
    const double uc2 = uc*uc;
    const double uc3 = uc2*uc;
    
    /* Compute spin-related functions*/
    double ggm[26];
    eob_dyn_s_GS(r, rc, drc_dr, 0., aK2, prstar, pphi, nu, chi1, chi2, X1, X2, c3, ggm);
    const double GS              = ggm[2];
    const double GSs             = ggm[3];
    const double dGS_dprstar     = ggm[4];
    const double dGSs_dprstar    = ggm[5];
    const double dGS_dr          = ggm[6];
    const double dGSs_dr         = ggm[7];
    const double dGSs_dpphi      = ggm[9];
    const double d2GS_dprstar20  = ggm[12];
    const double d2GSs_dprstar20 = ggm[13];
    
    /* Compute Hamiltonian and its derivatives */
    *Heff_orb         = sqrt( prstar2+A*(1. + pphi2*uc2 +  z3*prstar4*uc2) );
    *Heff             = *Heff_orb + (GS*S + GSs*Sstar)*pphi;
    *H                = sqrt( 1. + 2.*nu*(*Heff - 1.) )/nu;
    if (dHeff_dr != NULL)         *dHeff_dr         = pphi*(dGS_dr*S + dGSs_dr*Sstar) + 1./(2.*(*Heff_orb))*( dA*(1. + pphi2*uc2 + z3*prstar4*uc2) - 2.*A*uc3*drc_dr*(pphi2 + z3*prstar4) );
    if (dHeff_dprstar != NULL)    *dHeff_dprstar    = pphi*(dGS_dprstar*S + dGSs_dprstar*Sstar) + (prstar/(*Heff_orb))*(1. + 2.*A*uc2*z3*prstar2);
    if (d2Heff_dprstar20 != NULL) *d2Heff_dprstar20 = pphi*(d2GS_dprstar20*S + d2GSs_dprstar20*Sstar) +  (1./(*Heff_orb))*(1. + 2.*A*uc2*z3*prstar2); /* second derivative of Heff wrt to pr_star neglecting all pr_star^2 terms */
    if (dHeff_dpphi != NULL)      *dHeff_dpphi      = GS*S + (GSs + pphi*dGSs_dpphi)*Sstar + pphi*A*uc2/(*Heff_orb);
}


/** Computes the gyro-gravitomagnetic functions GS and GS*, that are called GS and GSs.
    r      => BL radius
    aK2    => squared Kerr parameter
    prstar => r* conjugate momentum
    nu     => symmetric mass ratio
    the CN3LO parameter is hard-coded in this routine 
    ggm is the output structure. */

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


/** Define radius rc that includes of LO spin-square coupling.  */
/* 
   The S1*S2 term coincides with the BBH one, no effect of structure.
   The self-spin couplings, S1*S1 and S2*S2 get a EOS-dependent coefficient, CQ, that describe the quadrupole
   deformation due to spin. Notation of Levi-Steinhoff, JCAP 1412 (2014), no.12, 003. Notation analogous to
   the parameter a of Poisson, PRD 57, (1998) 5287-5290 or C_ES^2 in Porto & Rothstein, PRD 78 (2008), 044013
   
   The implementation uses the I-Love-Q fits of Table I of Yunes-Yagi
   paper, PRD 88, 023009, the bar{Q}(bar{\lambda)^{tid}) relation, line 3 of the table. 
   The dimensionless bar{\lambda} love number is related to our apsidal constant as lambda = 2/3 k2/(C^5) so that both quantities have to appear here.  
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

// tidal rc with NLO coefficient that depends on C_Qi
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

// tidal rc with NNLO coefficient that depends on C_Qi
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
      + (-2171./212. - 269./28.*nu + 0.375*nu*nu)     *(at1*at1+at2*at2)
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

    double c_ss_nnlo = - (189./32. + 417.32*nu              )    *a0 *a0
                       + ( 11./32. - 127.32*nu + 0.375*nu*nu)    *a12*a12
                       + ( 87.16   -  2.625*nu              )*X12*a0 *a12;

    
    double rc2   = r2 + aK2*(1. + 2.*u) + u*c_ss_nlo + u2*c_ss_nnlo;
    *rc          = sqrt(rc2);
    double divrc = 1.0/(*rc);
    *drc_dr      = r*divrc*(1-(aK2 + 0.5*c_ss_nlo)*u3 - 0.5*u4*c_ss_nnlo);	
    *d2rc_dr2    = 1./r*(*drc_dr) + r*divrc*((3.*aK2+c_ss_nlo)*u4 + 2.*c_ss_nnlo*u5);
    
  }
  
}

// tidal rc @ NNLO with the addition of the LO spin^4 coefficient that depends on C_Q, C_Oct and C_Hex
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
      + (-2171./212. - 269./28.*nu + 0.375*nu*nu)     *(at1*at1+at2*at2)
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

    double c_ss_nnlo = - (189./32. + 417.32*nu              )    *a0 *a0
                       + ( 11./32. - 127.32*nu + 0.375*nu*nu)    *a12*a12
                       + ( 87.16   -  2.625*nu              )*X12*a0 *a12;

    
    double rc2   = r2 + aK2*(1. + 2.*u) + u*c_ss_nlo + u2*c_ss_nnlo;
    *rc          = sqrt(rc2);
    double divrc = 1.0/(*rc);
    *drc_dr      = r*divrc*(1-(aK2 + 0.5*c_ss_nlo)*u3 - 0.5*u4*c_ss_nnlo);	
    *d2rc_dr2    = 1./r*(*drc_dr) + r*divrc*((3.*aK2+c_ss_nlo)*u4 + 2.*c_ss_nnlo*u5);
    
  }
  
}

/* Non-spinning case -- rc = r */
void eob_dyn_s_get_rc_NOSPIN(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, double C_Oct1, double C_Oct2, double C_Hex1, double C_Hex2, int usetidal, 
		      double *rc, double *drc_dr, double *d2rc_dr2)
{
    *rc = r;
    *drc_dr = 1.;
    *d2rc_dr2 = 0.;
}

/* LO case with C_Q1 = 0 for tidal part*/
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

/** Root function to compute light-ring */
//TODO: THIS IS FOR NOSPIN
double eob_dyn_fLR(double r, void  *params)
{
  Dynamics *dyn = params;     
  double A,B,dA,d2A,dB,pl_hold;
  //if (dyn->use_spins) eob_metric_s(r, dyn, &A,&B,&dA,&d2A,&dB);
  //else
  eob_metric (r, dyn, &A,&B,&dA,&d2A,&dB,&pl_hold);
  double u = 1./r;
  double dA_u = (-dA)*SQ(r);
  return A + 0.5 * u * dA_u;
}

/** Root finder for adiabatic light-ring */
int eob_dyn_adiabLR(Dynamics *dyn, double *rLR)
{
  int status;
  int iter = 0, max_iter = 200;
  const double epsabs = 0.; /* if converges, precision is |r-r*| = epsabs + epsrel r*  */
  const double epsrel = 1e-10; 
  const gsl_root_fsolver_type *T;
  double x, x_lo, x_hi;
    
  /* Set interval to search root */
  if (dyn->use_tidal) {
    /* Tides are always temporarily set as = NNLO to compute LR, 
       But we may want to define different searches intervals */
    const int tides = EOBPars->use_tidal;
    if (tides == TIDES_TEOBRESUM_BHNS) {
      /* BHNS */
      //FIXME best interval
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
    /* x_lo = 0.9*eob_approxLR(dyn->nu); 
       x_hi = 1.1*eob_approxLR(dyn->nu); */
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

/** Root function to compute LSO */
//TODO: THIS IS FOR NOSPIN
double eob_dyn_fLSO(double r, void  *params)
{
  Dynamics *dyn = params;    
  double A,B,dA,d2A,dB,pl_hold;
  //if (dyn->use_spins) eob_metric_s(r, dyn, &A,&B,&dA,&d2A,&dB);
  //else                
  eob_metric  (r, dyn, &A,&B,&dA,&d2A,&dB,&pl_hold);
  double u = 1./r;
  double u2  = SQ(u);
  double dA_u = (-dA)*SQ(r);
  double d2A_u = d2A*SQ(r)*SQ(r) + 2*dA*SQ(r)*r;
  dB = u2*dA_u + 2.*A*u;
  double d2B = d2A_u*u2 + 4.*u*dA_u + 2*A;
  return ( dA_u*d2B - d2A_u*(dB) );
}

/** Root finder for adiabatic LSO */
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
  if (dyn->use_tidal) x_hi = 36.; 
  
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
