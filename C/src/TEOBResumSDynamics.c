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

  const double nu = EOBPars->nu; 
  const double z3 = 2.0*nu*(4.0-3.0*nu);

  /** Unpack y */ 
  const double phi    = y[EOB_EVOLVE_PHI];
  const double r      = y[EOB_EVOLVE_RAD];
  const double pphi   = y[EOB_EVOLVE_PPHI];
  const double prstar = y[EOB_EVOLVE_PRSTAR];
 
  /** Compute EOB Metric */
  double A, B, dA, d2A, dB;
  eob_metric(r, d, &A, &B, &dA, &d2A, &dB);

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
  const int usetidal = EOBPars->use_tidal;
  const int usespins = EOBPars->use_spins;

  /* Updated spins parallel to L, if required */
  if (usespins == MODE_SPINS_GENERIC) {
    
    double SA, SB; // projections of the spin parallel to hatL
    eob_spin_dyn_Sproj_interp(dyn->spins, t, &SA, &SB, NULL, 
			      NULL, NULL, NULL, 
			      EOBPars->spin_interp_integrate);
    
    const double M2 = SQ(EOBPars->M);
    const double XA = EOBPars->X1;
    const double XB = EOBPars->X2;
    
    EOBPars->chi1 = SA / (SQ(XA)*M2);
    EOBPars->chi2 = SB / (SQ(XB)*M2);
    
    set_spin_vars(XA,XB, EOBPars->chi1,EOBPars->chi2, 
		  &EOBPars->S1, &EOBPars->S2,
		  &EOBPars->a1, &EOBPars->a2,
		  &EOBPars->aK, &EOBPars->aK2,
		  &EOBPars->S, &EOBPars->Sstar);
    
  }
   
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
  const double C_Q1  = EOBPars->C_Q1;
  const double C_Q2  = EOBPars->C_Q2;
  const double C_Oct1 = EOBPars->C_Oct1;
  const double C_Oct2 = EOBPars->C_Oct2;
  const double C_Hex1 = EOBPars->C_Hex1;
  const double C_Hex2 = EOBPars->C_Hex2;
   
  /* Shorthands */
  const double r      = y[EOB_EVOLVE_RAD];
  const double prstar = y[EOB_EVOLVE_PRSTAR];
  const double pphi   = y[EOB_EVOLVE_PPHI];
  const double pphi2  = pphi*pphi;
  
  /** Compute Metric */
  double A, B, dA, d2A, dB;
  eob_metric_s(r, d, &A, &B, &dA, &d2A, &dB);
  
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
  double ggm0[14];
  eob_dyn_s_GS(r, rc, drc_dr, aK2, 0., pphi, nu, chi1, chi2, X1, X2, c3, ggm0);
  
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
    double ggm[14];
    eob_dyn_s_GS(r, rc, drc_dr, aK2, prstar, pphi, nu, chi1, chi2, X1, X2, c3, ggm);
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

void eob_dyn_s_GS(double r, double rc, double drc_dr, double aK2, double prstar, double pph, double nu, double chi1, double chi2, double X1, double X2, double cN3LO,
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

  double u   = 1./r;
  double u2  = u*u;
  
  double uc      = 1./rc;
  double uc2     = uc*uc;
  double uc3     = uc2*uc;
  double uc4     = uc3*uc;
  double prstar2 = prstar*prstar;
  double prstar4 = prstar2*prstar2;
  
  double GS0       = 2.*u*uc2;
  double dGS0_duc  = 2.*u2/drc_dr + 4.*u*uc;
  
  double GSs0          =  3./2.*uc3;
  double dGSs0_duc     =  9./2.*uc2;
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
    *drc_dr = 1;
    *d2rc_dr2 = 0;
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
    *drc_dr = 1;
    *d2rc_dr2 = 0;
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
  double A,B,dA,d2A,dB;
  //if (EOBPars->use_spins) eob_metric_s(r, dyn, &A,&B,&dA,&d2A,&dB);
  //else
  eob_metric (r, dyn, &A,&B,&dA,&d2A,&dB);
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
  if (EOBPars->use_tidal) {
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

/** Root function to compute LSO */
//TODO: THIS IS FOR NOSPIN
double eob_dyn_fLSO(double r, void  *params)
{
  Dynamics *dyn = params;    
  double A,B,dA,d2A,dB;
  //if (EOBPars->use_spins) eob_metric_s(r, dyn, &A,&B,&dA,&d2A,&dB);
  //else                
  eob_metric  (r, dyn, &A,&B,&dA,&d2A,&dB);
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



/** spin dynamics */


/** Compute alpha from Lhat */
double eob_spin_dyn_alpha(double Lhx, double Lhy, double Lhz)
{
  return atan2(Lhy,Lhx);
}

/** Compute beta from Lhat */
double eob_spin_dyn_beta(double Lhx, double Lhy, double Lhz)
{
  return acos(Lhz);
}

/** r.h.s. of the PN precessing equations 
    https://arxiv.org/abs/1307.4418
    https://arxiv.org/abs/1703.03967
*/
int eob_spin_dyn_rhs_PN(double t, const double y[], double dy[], void *d)
{
  (void)(t); /* avoid unused parameter warning */
  DynamicsSpin *dyn = d;  

  const double q = EOBPars->q; // Assume q = MA/MB >=1 //CHECKME paper convention!
  const double nu = EOBPars->nu; 
  const double nu2 = nu*nu;
  const double nu3 = nu2*nu;
  const double nu4 = nu3*nu;
  const double nu5 = nu4*nu;

  const double M = EOBPars->M; 
  const double M2 = SQ(M);
  const double MA = M*nu_to_X1(nu); 
  const double MB = M - MA;
  const double dm = MA - MB; 
  const double ma_o_mb = MA/MB;
  const double mb_o_ma = MB/MA;

  const double Pi2 = SQ(Pi);
  const double Pi3 = Pi2*Pi;
  const double oothree = 0.3333333333333333; // 1/3
  const double eleven_o_three = 3.6666666666666665; // 11/3

  /** Unpack y */ 
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
  const double v = pow(omg/M, oothree); // division by M !
  const double v2 = v*v;
  const double v3 = v2*v;
  const double v4 = v3*v;
  const double v5 = v4*v;
  const double v6 = v5*v;
  const double v7 = v6*v;
  const double v8 = v7*v;
  const double v9 = v8*v;
  
  /** rhs */

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
  
  for(int a=Ix; a<IN3; a++) 
    OmgA[a] = v5*(nu*(2+1.5/q) - 1.5*v*qSABLh)*Lh[a] + 0.5*v6*SB[a];
  for(int a=Ix; a<IN3; a++) 
    OmgB[a] = v5*(nu*(2+1.5*q) - 1.5*v*SABqLh)*Lh[a] + 0.5*v6*SA[a];
  
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
  
  
  /* N4LO */
  /*
  double OmgANLO[IN3],  OmgBNLO[IN3],  OmgANNLO[IN3],  OmgBNNLO[IN3],  OmgAN4LO[IN3],  OmgBN4LO[IN3];
  double SdotANLO[IN3], SdotBNLO[IN3], SdotANNLO[IN3], SdotBNNLO[IN3], SdotAN4LO[IN3], SdotBN4LO[IN3];
  double LNdotN4LO[IN3];

  const double cv7A = v7*( 0.5625 + 1.25*nu - 0.04166666666666666*nu2 + dm*(-0.5625+0.625*nu) );
  const double cv7B = v7*( 0.5625 + 1.25*nu - 0.04166666666666666*nu2 - dm*(-0.5625+0.625*nu) );
  const double cv9A = v9*( 0.84375 + 0.1875*nu - 3.28125*nu2 - 0.02083333333333*nu3
			  + dm*(-0.84375+4.875*nu-0.15625*nu2) );
  const double cv9B = v9*( 0.84375 + 0.1875*nu - 3.28125*nu2 - 0.02083333333333*nu3
			  - dm*(-0.84375+4.875*nu-0.15625*nu2) );

  for(int a=Ix; a <IN3; a++){
    OmgANLO[a] = v5*(nu*(2+1.5/q) - 1.5*v*qSABLh)*Lh[a] + 0.5*v6*SB[a];
    OmgBNLO[a] = v5*(nu*(2+1.5*q) - 1.5*v*SABqLh)*Lh[a] + 0.5*v6*SA[a];
    OmgANNLO[a]= OmgANLO[a]  + cv7A*Lh[a];
    OmgBNNLO[a]= OmgBNLO[a]  + cv7B*Lh[a];
    OmgAN4LO[a]= OmgANNLO[a] + cv9A*Lh[a];
    OmgBN4LO[a]= OmgBNNLO[a] + cv9B*Lh[a];
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
  const double csA = -0.25*(3+M/MA);
  const double csB = -0.25*(3+M/MB);
  const double csAL = -0.08333333333333333*(1+27.*M/MA);
  const double csBL = -0.08333333333333333*(1+27.*M/MB);
  //const double csANL= 0.0625*(-31.+7*M/MA) + 0.0208333333333333*nu*(9. + 22* M/MA)
  //const double csBNL= 0.0625*(-31.+7*M/MB) + 0.0208333333333333*nu*(9. + 22* M/MB)
  const double L2PN = 1 + v2*(1.5+0.1666666666666667*nu) + v4*(3.375 - 2.375*nu + 0.04166666666666666*nu2);
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
  
  */
  
  /* dot gamma = dot alpha(t) * cos(beta(t)) = dot alpha(t) * Lhz */
 
  const double div = SQ(Lh[Ix]) + SQ(Lh[Iy]);
  if (div == 0.)
    /* Lh is along z
       alpha is undefined, beta = 0 */
    dy[EOB_EVOLVE_SPIN_gam] = 0.;
  else
    dy[EOB_EVOLVE_SPIN_gam] = Lh[Iz] * ( Lh[Iy] * dy[EOB_EVOLVE_SPIN_Lx] - Lh[Ix] * dy[EOB_EVOLVE_SPIN_Ly] )/( SQ(Lh[Ix]) + SQ(Lh[Iy]) );
    
  /* dot omg (Rad.React.) */

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
  
  //TODO: precompted spin independent coefs as first call       
  // sgima, beta, a_3-8 are spin-dependend, precompute the fractions
  
  // https://arxiv.org/abs/1307.4418 , App.A 
  
  double sigma4 = ( 247./48*SAdotSB - 721./48*SAdotLh*SBdotLh )/(nu*SQ(M2))
    + (233./96*SA2 - 719./96*SQ(SAdotLh))/(M2*SQ(MA))
    + (233./96*SB2 - 719./96*SQ(SBdotLh))/(M2*SQ(MB)); 
  
  
  beta[3] = (113./12 + 25./4*mb_o_ma)*SAdotLh/M2
    + (113./12 + 25./4*ma_o_mb)*SBdotLh/M2; 

  beta[5] = ((31319./1008-1159./24*nu) + mb_o_ma*( 809./84 - 281./8*nu))*SAdotLh/M2
    + ((31319./1008-1159./24*nu) + ma_o_mb*( 809./84 - 281./8*nu))*SBdotLh/M2; 

  beta[6] = Pi/M2 * (75./2+151./6*mb_o_ma)*SAdotLh
    + Pi/M2 * (75./2+151./6*ma_o_mb)*SBdotLh; 

  beta[7] = ((130325./756 - 796069./2016*nu+100019./864*nu2)+ mb_o_ma*(1195759./18144-257023./1008*nu+2903./32*nu2))*SAdotLh/M2 
    + ((130325./756 - 796069./2016*nu+100019./864*nu2)+ ma_o_mb*(1195759./18144-257023./1008*nu+2903./32*nu2))*SBdotLh/M2;
  
  beta[8] = Pi/M2*((76927./504 -220055./672*nu) + mb_o_ma*(1665./28-50483./224*nu))*SAdotLh
    + Pi/M2*((76927./504 -220055./672*nu) + ma_o_mb*(1665./28-50483./224*nu))*SBdotLh;

  
  b[6] = -1712./315;  

  b[8] = - 856./315*nu + 124741./4410;  

  b[9] = - 6848./105* Pi;  

  b[10] = 3090781./26460*nu -  2354./945*nu2 - 11821184./1964655;  

  b[11] = 311233./5880*Pi - 3424./315*Pi*nu;

  
  a[0] = 96./5*nu;

  a[2] = -743./336-11./4*nu;

  a[3] = 4*Pi-beta[3];

  a[4] = 34103./18144 + 13661./2016*nu + 59./18*nu2 - sigma4;

  a[5] = -4159./672*Pi - 189./8*Pi*nu - beta[5];

  a[6] = 16447322263./139708800 + 16./3*Pi2 - 856./105*log(16.) - 1712./105*EulerGamma - beta[6]
    + nu *( 451./48*Pi2- 56198689./217728) + nu2*541./896 - nu3*5605./2592;

  a[7] = - 4415./4032*Pi + 358675./6048*Pi*nu +91495./1512*Pi*nu2 - beta[7];

  a[8] = 3971984677513./25427001600 + 127751./1470*Log2  - 47385./1568*Log3 + 124741./4410*EulerGamma -361./126*Pi2 + 82651980013./838252800*nu - 1712./315*nu*Log2
    - 856./315*EulerGamma*nu  - 31495./8064*Pi2*nu + 54732199./93312*nu2- 3157./144*Pi2*nu2  - 18927373./435456*nu3 -95./3888*nu4 -beta[8];

  a[9] = 343801320119./745113600*Pi- 13696./105*Pi*Log2 -  6848./105*Pi*EulerGamma - 51438847./48384*Pi*nu + 205./6*Pi3*nu + 42680611./145152*Pi*nu2  +  9731./1344*Pi*nu3;

  a[10] = 29619150939541789./36248733480960  -107638990./392931*Log2 + 616005./3136*Log3 - 11821184./1964655*EulerGamma - 21512./1701*Pi2 - 884576519037433./228843014400*nu 
    + 2105111./8820*nu*Log2 - 15795./3136*nu*Log3 + 3090781./26460*EulerGamma*nu+ 14555455./217728*Pi2*nu  + 1175999369413./914457600*nu2 - 4708./945*nu2*Log2
    - 126809./3024*Pi2*nu2 - 2354./945*EulerGamma*nu2 - 9007327699./11757312*nu3 + 9799./384*Pi2*nu3 + 51439207./1741824*nu4 - 34613./186624*nu5;

  a[11] =  91347297344213./81366405120*Pi+ 5069891./17640*Pi*Log2- 142155./784*Pi*Log3  + 311233./5880*Pi*EulerGamma - 1903651780081./4470681600*Pi*nu- 6848./315*Pi*nu*Log2
    - 3424./315*Pi*EulerGamma*nu - 26035./16128*Pi3*nu + 1760705531./290304*Pi*nu2 - 112955./576*Pi3*nu2 - 7030123./13608*Pi*nu3 + 49187./6048*Pi*nu4;  
  
  // Eq.(A1) of https://arxiv.org/abs/1307.4418 for Momega
  dy[EOB_EVOLVE_SPIN_Momg] = 0.;
  for (int i=11; i<=2; i--)
    dy[EOB_EVOLVE_SPIN_Momg] += (a[i] + b[i]*lnomg)*pow(omg,(double)i*oothree);
  dy[EOB_EVOLVE_SPIN_Momg] += 1.;  
  dy[EOB_EVOLVE_SPIN_Momg] *= a[0]*pow(omg, eleven_o_three); // LO
      
  return GSL_SUCCESS;
}

/** r.h.s. of the PN precessing equations for the angle */ 
int eob_spin_dyn_rhs_PN_abc(double t, const double y[], double dy[], void *d)
{

  //TODO
  // SB: This is a better approach, just evolve alpha, beta, gamma and omega
  //     Not sure why we did not do this from the beginning ...
  
  return GSL_SUCCESS;
}

/** r.h.s. of the EOB precessing equations */ 
int eob_spin_dyn_rhs_EOB(double t, const double y[], double dy[], void *d)
{

  //TODO
  
  return GSL_SUCCESS;
}


/** Precessing dynamics ODE integration 
    The initial data are those stored in dyn->y 
    integration starts at dyn->t */
int eob_spin_dyn_integrate(DynamicsSpin *dyn)
{
  const int chunk = dyn->size;
  int size = chunk; // can change

  /* index of closest element to dyn->t in dyn->time 
     initial data refer to this time */
  const int i0 = find_point_bisection(dyn->t, dyn->size, dyn->time, 0);//CHECKME: 0 or 1?
  if (DEBUG) printf("initial data at index %d\n",i0);
  
  /* Set initial data */
  for (int v=0; v<EOB_EVOLVE_SPIN_NVARS; v++)
    dyn->data[v][i0] = dyn->y[v];
  
  /* GSL integrator memory */
  //dyn->omg_stop = EOBPars->spin_odes_omg_stop; 
  //dyn->t_stop = EOBPars->spin_odes_t_stop;
  if (DEBUG) printf("omg_stop = %e\n", dyn->omg_stop);
  
  dyn->dt = EOBPars->spin_odes_dt; 
  const double ode_abstol = EOBPars->ode_abstol;
  const double ode_reltol = EOBPars->ode_reltol;
  const double tstop = MAX(dyn->t_stop, EOBPars->ode_tmax);
  
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
  while (1) {
    if (VERBOSE) printf("iter %09d\n", iter); 
    iter++;

    //GSLSTATUS = gsl_odeiv2_evolve_apply_fixed_step (e, c, s, &sys, &dyn->t, dyn->dt, dyn->y);//uniform
    GSLSTATUS = gsl_odeiv2_evolve_apply (e, c, s, &sys, &dyn->t, tstop, &dyn->dt, dyn->y);
    
    /** Check for failures ... */
    if (GSLSTATUS != GSL_SUCCESS) {
      printf("GSL Error = %d", GSLSTATUS);
      return ERROR_ODEINT;
    }
    
    /** Update alpha and beta angles */
    dyn->y[EOB_EVOLVE_SPIN_alp] = eob_spin_dyn_alpha(dyn->y[EOB_EVOLVE_SPIN_Lx],
						     dyn->y[EOB_EVOLVE_SPIN_Ly],
						     dyn->y[EOB_EVOLVE_SPIN_Lz]);
    dyn->y[EOB_EVOLVE_SPIN_bet] = eob_spin_dyn_beta(dyn->y[EOB_EVOLVE_SPIN_Lx],
						    dyn->y[EOB_EVOLVE_SPIN_Ly],
						    dyn->y[EOB_EVOLVE_SPIN_Lz]);
    
    /** Update size and push arrays (if needed) */
    if (iter==size) {
      /* if (DEBUG)  printf("Push memory\n"); */ 
      size += chunk;
      dyn->size = size;
      DynamicsSpin_push (&dyn, size);
    }
    
    dyn->time[iter]  = dyn->t; 
    for (int v=0; v<EOB_EVOLVE_SPIN_NVARS; v++)
      dyn->data[v][i0+iter]  = dyn->y[v];   

    /** Stop integration */
    if ((dyn->omg_stop>0) && (dyn->y[EOB_EVOLVE_SPIN_Momg] > dyn->omg_stop))
      break;
    if ((dyn->t_stop>0) && (dyn->t > dyn->t_stop))
      break;
    
  } /* end time iteration */

  /** Resize to actual size */
  size = i0 + iter + 1;
  DynamicsSpin_push (&dyn, size);
  
  /** Free ODE system solver */
  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
  gsl_odeiv2_driver_free (d);
  
  return OK;
}

/** Precessing dynamics main driver routine */
int eob_spin_dyn(DynamicsSpin *dyn)
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
  double m1 = EOBPars->M*nu_to_X1(EOBPars->nu);
  double m2 = EOBPars->M - m1;
  const double M12 = SQ(m1);
  const double M22 = SQ(m2);  
  dyn->y[EOB_EVOLVE_SPIN_SxA] = EOBPars->chi1x *M12; 
  dyn->y[EOB_EVOLVE_SPIN_SyA] = EOBPars->chi1y *M12;
  dyn->y[EOB_EVOLVE_SPIN_SzA] = EOBPars->chi1z *M12;
  dyn->y[EOB_EVOLVE_SPIN_SxB] = EOBPars->chi2x *M22;
  dyn->y[EOB_EVOLVE_SPIN_SyB] = EOBPars->chi2y *M22;
  dyn->y[EOB_EVOLVE_SPIN_SzB] = EOBPars->chi2z *M22;
  dyn->y[EOB_EVOLVE_SPIN_Lx] = 0; //FIXME Lh t=0 ?
  dyn->y[EOB_EVOLVE_SPIN_Ly] = 0;
  dyn->y[EOB_EVOLVE_SPIN_Lz] = 1.;
  dyn->y[EOB_EVOLVE_SPIN_alp] =  eob_spin_dyn_alpha(dyn->y[EOB_EVOLVE_SPIN_Lx],
						    dyn->y[EOB_EVOLVE_SPIN_Ly],
						    dyn->y[EOB_EVOLVE_SPIN_Lz]);
  dyn->y[EOB_EVOLVE_SPIN_bet] = eob_spin_dyn_beta(dyn->y[EOB_EVOLVE_SPIN_Lx],
						  dyn->y[EOB_EVOLVE_SPIN_Ly],
						  dyn->y[EOB_EVOLVE_SPIN_Lz]);
  dyn->y[EOB_EVOLVE_SPIN_gam] = Pi/2.; // P.8 https://arxiv.org/abs/2004.09442
  dyn->y[EOB_EVOLVE_SPIN_Momg] = Pi * EOBPars->initial_frequency; 
  
  for (int v=0; v<EOB_EVOLVE_SPIN_NVARS; v++)
    dyn->data[v][0] = dyn->y[v];
  
  /** Integrate ODEs */
  return eob_spin_dyn_integrate(dyn);
  
}

/** Helper routine to interpolate Euler angles at given time */
void eob_spin_dyn_abc_interp(DynamicsSpin *dyn, double time,
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
      eob_spin_dyn_integrate(dyn);      

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
    beta = interp_spline_pt(dyn->time, dyn->data[EOB_EVOLVE_SPIN_bet], dyn->size, time);
    gamma = interp_spline_pt(dyn->time, dyn->data[EOB_EVOLVE_SPIN_gam], dyn->size, time);
  }
  
  *alpha_p = alpha;
  *beta_p = beta;
  *gamma_p = gamma;
}

/** Helper routine to interpolate spins parallel/perpendicular to orbital ang. mom. at given time */
void eob_spin_dyn_Sproj_interp(DynamicsSpin *dyn, double time,
			       double *SApara, double *SBpara, double *Spara,
			       double *SAperp, double *SBperp, double *Sperp, // these are 3-vectors
			       int continue_integration)
{
  const int smax = dyn->size-1;
  double SA[IN3], SB[IN3], Lh[IN3];
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
      eob_spin_dyn_integrate(dyn);      
      
    } else {
      
      SA[Ix] = dyn->data[EOB_EVOLVE_SPIN_SxA][smax];
      SA[Iy] = dyn->data[EOB_EVOLVE_SPIN_SyA][smax];
      SA[Iz] = dyn->data[EOB_EVOLVE_SPIN_SzA][smax];
      
      SB[Ix] = dyn->data[EOB_EVOLVE_SPIN_SxB][smax];
      SB[Iy] = dyn->data[EOB_EVOLVE_SPIN_SyB][smax];
      SB[Iz] = dyn->data[EOB_EVOLVE_SPIN_SzB][smax];
      
      Lh[Ix] = dyn->data[EOB_EVOLVE_SPIN_Lx][smax];
      Lh[Iy] = dyn->data[EOB_EVOLVE_SPIN_Ly][smax];
      Lh[Iz] = dyn->data[EOB_EVOLVE_SPIN_Lz][smax];
      
      interp = 0;// skip interpolation below
      if (VERBOSE) errorexit("Spin dynamics too short to interp!\n");
    }
  }

  
  /* Interp */
  if (interp) {
    SA[Ix] = interp_spline_pt(dyn->time, dyn->data[EOB_EVOLVE_SPIN_SxA], dyn->size, time);
    SA[Iy] = interp_spline_pt(dyn->time, dyn->data[EOB_EVOLVE_SPIN_SyA], dyn->size, time);
    SA[Iz] = interp_spline_pt(dyn->time, dyn->data[EOB_EVOLVE_SPIN_SzA], dyn->size, time);
    
    SB[Ix] = interp_spline_pt(dyn->time, dyn->data[EOB_EVOLVE_SPIN_SxB], dyn->size, time);
    SB[Iy] = interp_spline_pt(dyn->time, dyn->data[EOB_EVOLVE_SPIN_SyB], dyn->size, time);
    SB[Iz] = interp_spline_pt(dyn->time, dyn->data[EOB_EVOLVE_SPIN_SzB], dyn->size, time);
    
    Lh[Ix] = interp_spline_pt(dyn->time, dyn->data[EOB_EVOLVE_SPIN_Lx], dyn->size, time);
    Lh[Iy] = interp_spline_pt(dyn->time, dyn->data[EOB_EVOLVE_SPIN_Ly], dyn->size, time);
    Lh[Iz] = interp_spline_pt(dyn->time, dyn->data[EOB_EVOLVE_SPIN_Lz], dyn->size, time);
  }

  /* Total spin */
  //double SAB[IN3], S2;
  //for (int i=0; i<IN3; i++) SAB[i] = SA[i]+SB[i];
  //vect_dot3(SAB, SAB, &S2);

  /* direction of Lh */
  //FIXME: this should not be needed
  double n[IN3], normL;
  vect_dot3(Lh, Lh, &normL); 
  const double oonormL = 1./normL; 
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


