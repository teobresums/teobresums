/**
 * Copyright (C) 2017 Alessandro Nagar, Gregorio Carullo, Ka Wa Tsang, Philipp Fleig, Sebastiano Bernuzzi, Walter Del Pozzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#include "TEOBResumS.h"

/** r.h.s. of EOB Hamiltonian dynamics, no spins version */ 
int eob_dyn_rhs(double t, const double y[], double dy[], void *dyn)
{
  
  (void)(t); /* avoid unused parameter warning */
  Dynamics *d = dyn;  

  const double nu = d->nu;  
  const double z3 = 2.0*nu*(4.0-3.0*nu);
  
  const double r      = y[EOB_EVOLVE_RAD];
  const double pphi   = y[EOB_EVOLVE_PPHI];
  const double prstar = y[EOB_EVOLVE_PRSTAR];
  
  /** Compute EOB Metric */
  double A, B, dA, d2A, dB;
  Metric(r, d, &A, &B, &dA, &d2A, &dB);

  const double ooA = 1./A;
  const double ooB = 1./B;
  const double sqrAB = sqrt(A/B);  
  
  /** r.h.s eqns */  
  const double prstar2 = prstar*prstar;
  const double prstar3 = prstar2*prstar;
  const double prstar4 = prstar3*prstar;
  const double pphi2   = pphi*pphi;
  
  const double r2 = r*r;
  const double u  = 1./r;
  const double u2 = u*u;
  const double u3 = u2*u;
  
  const double Heff  = sqrt(prstar2+A*(1. + pphi2*u2 +  z3*prstar4*u2));
  const double H     = sqrt( 1. + 2.*nu*(Heff - 1.) )/nu;
  const double E     = H*nu;
  const double ooHE  = 1./(Heff*E);
  const double Omega = A*pphi*u2*ooHE;

  const double sqrW = sqrt(A*(1. + pphi2*u2));  

  /* r evol eqn rhs */
  dy[EOB_EVOLVE_RAD]  = (sqrAB)*(prstar+2.0*z3*A*prstar3*u2)*ooHE;
  
  /* phi evol eqn rhs */
  dy[EOB_EVOLVE_PHI] = Omega;
  
  /* prstar evol eqn rhs */
  dy[EOB_EVOLVE_PRSTAR]  = (-0.5*sqrAB)*(dA + ( pphi2 + z3*prstar4 )*(dA*u2-2.0*A*u3))*ooHE;
  
  /* pphi evol eqn rhs */
  const double psi = 2.*(1.0 + 2.0*nu*(sqrW - 1.0))/(r2*dA);
  /* psi = 2.*(1.0 + 2.0*nu*(Heff - 1.0))/(r2*dA); */
   
  /** Compute flux */
  const double r_omega = r*cbrt(psi);
  const double v_phi   = r_omega*Omega;
  const double x       = v_phi * v_phi;
  const double jhat    = pphi/(r_omega*v_phi);  
  const double tmpE    = 1./Heff+nu/(E*E);
  
  const double dHeff_dr      = 0.5*(dA + (pphi2 + z3*prstar4)*(dA*u2 - 2.*A*u3))/Heff;
  const double dHeff_dprstar = (prstar + z3*2.0*A*u2*prstar3)/Heff;
  const double dprstar_dt    = dy[EOB_EVOLVE_PRSTAR];
  const double dr_dt         = dy[EOB_EVOLVE_RAD];
  const double ddotr_dr      = sqrAB*( (prstar + z3*2.*A*u2*prstar3)*(0.5*(dA*ooA-dB*ooB)-dHeff_dr*tmpE)+ 2.0*z3*(dA*u2 - 2.*A*u3)*prstar3)/E;
  const double ddotr_dprstar = sqrAB*( 1.+z3*6.*A*u2*prstar2-(prstar + z3*2.*A*u2*prstar3)*dHeff_dprstar*tmpE)*ooHE;
  
  /* Approximate ddot(r) without Flux */
  const double ddotr = dprstar_dt*ddotr_dprstar + dr_dt*ddotr_dr;
  
  dy[EOB_EVOLVE_PPHI] = eob_flx_Flux(x,Omega,r_omega,E,Heff,jhat,r, prstar,ddotr,dyn);

  if(d->store) {
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

/** r.h.s. of EOB Hamiltonian dynamics, spins version */ 
int eob_dyn_rhs_s(double t, const double y[], double dy[], void *dyn)
{
    
  (void)(t); /* avoid unused parameter warning */
  Dynamics *d = dyn;

  const double nu    = d->nu;
  const double S     = d->S;
  const double Sstar = d->Sstar;
  const double chi1  = d->chi1;
  const double chi2  = d->chi2;
  const double X1    = d->X1;
  const double X2    = d->X2;
  const double c3    = d->cN3LO;
  const double aK2   = d->aK2;
  const double a1    = d->a1;
  const double a2    = d->a2;
  const double C_Q1  = d->C_Q1;
  const double C_Q2  = d->C_Q2;
  const int usetidal = d->use_tidal;

  const double z3     = 2.*nu*(4.-3.*nu);    
  const double r      = y[EOB_EVOLVE_RAD];
  const double prstar = y[EOB_EVOLVE_PRSTAR];
  const double pph    = y[EOB_EVOLVE_PPHI];
  const double pphi2  = pph*pph;
  
  /** Compute Metric */
  double A, B, dA, d2A;
  if (usetidal) {
    eob_metric(r, d, &A, &B, &dA, &d2A);
  } else {
    eob_metric_s(r, d, &A, &B, &dA, &d2A);
  }
 
  /* shorthands */
  const double prstar2 = prstar*prstar;
  const double prstar3 = prstar2*prstar;
  const double prstar4 = prstar3*prstar;
  
  double rc, drc_dr, d2rc_dr;
  eob_dyn_s_get_rc(r, nu, a1,a2,aK2, C_Q1,C_Q2, usetidal, &rc,&drc_dr,&d2rc_dr);  
  
  const double uc     = 1./rc;
  const double uc2    = uc*uc;
  const double uc3    = uc2*uc;
  
  double ggm[14];
  eob_dyn_s_GS(r, rc, drc_dr, aK2, prstar, pph, nu, chi1, chi2, X1, X2, c3, ggm);
  
  const double GS              = ggm[2];
  const double GSs             = ggm[3];
  const double dGS_dprstar     = ggm[4];
  const double dGSs_dprstar    = ggm[5];
  const double dGS_dr          = ggm[6];
  const double dGSs_dr         = ggm[7];
  const double dGSs_dpph       = ggm[9];
  const double d2GS_dprstar20  = ggm[12];
  const double d2GSs_dprstar20 = ggm[13];

  const double Heff_orb = sqrt( prstar2+A*(1. + pphi2*uc2 +  z3*prstar4*uc2) );  
  const double Heff = Heff_orb + (GS*S + GSs*Sstar)*pph;
  const double H = sqrt( 1. + 2.*nu*(Heff - 1.) );
  const double ooH = 1./H;
  const double sqrtAbyB = sqrt(A/B);
  const double dHeff_dr = pph*(dGS_dr*S + dGSs_dr*Sstar) + 1./(2.*Heff_orb)*( dA*(1. + pphi2*uc2 + z3*prstar4*uc2) - 2.*A*uc3*drc_dr*(pphi2 + z3*prstar4) );
  const double dp_rstar_dt_0 = - sqrtAbyB*ooH*dHeff_dr;
  const double dHeff_dprstar = pph*(dGS_dprstar*S + dGSs_dprstar*Sstar) + (prstar/Heff_orb)*(1. + 2.*A*uc2*z3*prstar2);
  
  /* second derivative of Heff wrt to pr_star neglecting all pr_star^2 terms */
  const double d2Heff_dprstar20 = pph*(d2GS_dprstar20*S + d2GSs_dprstar20*Sstar) +  (1./Heff_orb)*(1. + 2.*A*uc2*z3*prstar2);
  const double ddotr_dp_rstar = sqrtAbyB*ooH*d2Heff_dprstar20;
  
  const double dHeff_dpph = GS*S + (GSs + pph*dGSs_dpph)*Sstar + pph*A*uc2/Heff_orb;
  const double Omg = ooH*dHeff_dpph; 
  
  /* approximate ddot(r)_0 without Fphi, order pr_star^2 neglected */
  const double ddotr = dp_rstar_dt_0*ddotr_dp_rstar; 
  
  /* r evol eqn rhs */
  dy[EOB_EVOLVE_RAD] = sqrtAbyB*ooH*dHeff_dprstar;
  
  /* phi evol eqn rhs */
  dy[EOB_EVOLVE_PHI] = Omg;
  
  /* dp_{r*}/dt */ 
  dy[EOB_EVOLVE_PRSTAR] = -sqrtAbyB*ooH*dHeff_dr;
  
  /** Compute flux */
  
  /* Compute here the new r_omg radius 
     Compute same quantities with prstar=0. This to obtain psi.
     Procedure consistent with the nonspinning case
  */
  double gmm0[14];
  eob_dyn_s_GS(r, rc, drc_dr, aK2, 0., pph, nu, chi1, chi2, X1, X2, c3, gmm0);
  
  const double GS_0       = ggm0[2];
  const double GSs_0      = ggm0[3];
  const double dGS_dr_0   = ggm0[6];
  const double dGSs_dr_0  = ggm0[7];
  const double Heff_orb_0 = sqrt(A*(1.0 + pphi2*uc2));    /* effective Hamiltonian H_0^eff */
  const double Heff_0     = Heff_orb_0 + (GS_0*S + GSs_0*Sstar)*pph;
  const double H0         = sqrt(1.0 + 2.0*nu*(Heff_0 - 1.0) );
  const double ooH0       = 1./H0;
  const double Gtilde     = GS_0*S     + GSs_0*Sstar;
  const double dGtilde_dr = dGS_dr_0*S + dGSs_dr_0*Sstar;
  const double duc_dr     = -uc2*drc_dr;
  const double psic       = (duc_dr + dGtilde_dr*rc*sqrt(A/pphi2 + A*uc2)/A)/(-0.5*dA);
  const double r_omg      = pow( ((1./sqrt(rc*rc*rc*psic))+Gtilde)*ooH0, -2./3. );
  const double v_phi      = r_omg*Omg;
  const double x          = v_phi*v_phi;
  const double jhat       = pph/(r_omg*v_phi);
  
  const double Fphi = eob_flx_Flux_s(x,Omg,r_omg,H,Heff,jhat,r,prstar,ddotr,dyn);
  
  dy[EOB_EVOLVE_PPHI] = Fphi;

  if (d->store) {
    /* Store values */
    d->t = t;
    d->r = r;
    d->phi = y[EOB_EVOLVE_PHI];
    d->pphi = pphi;
    d->prstar = prstar;
    d->Omg = Omg;
    d->Omg_orb = Omg;
    d->H = H;
    d->E = E;
    d->Heff = Heff;
    d->A = A;
    d->dA = dA;
    d->d2A = d2A;
    d->B = B;
    d->dB = dB;
    d->psi = psi;
    d->r_omega = r_omg;
    d->v_phi = v_phi;
    d->jhat = jhat;
    d->ddotr = ddotr;
  }
  
  return GSL_SUCCESS;
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
  static double cs10,cs20,cs30,cs02,cs12,cs04;

  /* Compute the nu-dep. coefficient at first call only */
  static int firstcall = 1;  
  if (firstcall) {
    firstcall = 0;      
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
  
  double hGS  =  1./(1.  + c10*uc + c20*uc2 + c30*uc3 + c02*prstar2 + c12*uc*prstar2 + c04*prstar4);   
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
  double dhGS_duc  = -hGS*hGS*(  c10 + 2.*c20*uc  + 3.*c30*uc2);
  double dhGSs_duc = -hGSs*hGSs*(cs10 + 2.*cs20*uc + 3.*cs30*uc2 + 4.*cs40*uc3);
  
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
void eob_dyn_s_get_rc(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, int usetidal, 
		      double *rc, double *drc_dr, double *d2rc_dr2)
{

  double u   = 1./r;
  double u2  = u*u;
  double u3  = u*u2;
  double r2  = r*r;
    
  if (usetidal) {    
    /* BNS effective spin parameter */
    double a02      = C_Q1*at1*at1 + 2.*at1*at2 + C_Q2*at2*at2;
    /* tidally-modified centrifugal radius */
    double rc2 = r2 + a02*(1.+2.*u);
    *rc         = sqrt(rc2);
    *drc_dr     = r/rc*(1.-a02*u3);
    *d2rc_dr2   = 1./rc*( 1.-drc_dr*r/rc*(1.-a02*u3)+2.*a02*u3);
#if (0)
    /* NO spin-spin-tidal couplings */
    double rc2 = r2;
    rc = r;
    drc_dr = 1;
    d2rc_dr2 = 0;
#endif
  } else {
    double X12      = sqrt(1.-4.*nu);   
    double alphanu2 = 1. + 0.5/aK2*(- at2*at2*(5./4. + 5./4.*X12 + nu/2.) - at1*at1*(5./4. - 5./4.*X12 +nu/2.) + at1*at2*(-2.+nu));
    double rc2 = r2 + aK2*(1. + 2.*alphanu2/r);
    *rc         = sqrt(rc2);
    *drc_dr     = r/rc*(1.+aK2*(-alphanu2*u3 ));
    *d2rc_dr2   = 1./rc*(1.-drc_dr*r/rc*(1.-alphanu2*aK2*u3)+ 2.*alphanu2*aK2*u3);
  }
  
}

/** Root function to compute light ring */
double eob_dyn_fLR(double r, Dynamics *dyn)
{    
  double A,B,dA,d2A,dB;
  Metric(r, dyn, &A,&B,&dA,&d2A,&dB);
  double dA_u = (-dA*)*SQ(r); 
  return (A +(0.5*dA_u)/r);
}

/** Root finder for adiabtic light-ring */
double eob_dyn_adiabLR(void *params)
{
  int status;
  int iter = 0, max_iter = 200;
  const gsl_root_fsolver_type *T;
  double rLR;
  double x_lo = 0.1, x_hi = 15.;
  
  gsl_root_fsolver *s;
  gsl_function F;
  F.function = &eob_dyn_fLR;
  F.params = params; // Dynamics dyn
  T = gsl_root_fsolver_bisection;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      rLR    = gsl_root_fsolver_root (s);
      x_lo   = gsl_root_fsolver_x_lower (s);
      x_hi   = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi, 0, 1e-10);
    }
  while (status == GSL_CONTINUE && iter < max_iter);
  gsl_root_fsolver_free (s);
  
  return rLR;
}













/** Compute orbital frequency from arrays */
// used anywhere?
void get_Omg_orb(double *r, double *pph, double *pr_star, double *A, double *B, int size, void *params, double *Omg_orb)
{

  // fixme:
  const double nu   = d->nu;
  const double aK2  = d->aK2;
  const double S1   = d->S1;
  const double S2   = d->S2;
  const double a1   = d->a1;
  const double a2   = d->a2;
  const double X1   = d->X1;
  const double X2   = d->X2;
  const double chi1 = d->chi1;
  const double chi2 = d->chi2;
  const double c3   = d->cN3LO;
  const double C_Q1 = d->C_Q1;
  const double C_Q2 = d->C_Q2;
  const double aK2  = d->aK2;
  
  const double z3    = 2.*nu*(4.-3.*nu);
  const double S     = S1 + S2;
  const double Sstar = X2*a1 + X1*a2;
  const int usetidal  = d->use_tidal;

  double rc, drc_dr,d2rc_dr, ggm[14];
  int i;
  for ( i=0; i<size; i++) {
        
    double rc, drc_dr, d2rc_dr;
    s_get_rc(r, nu, a1,a2,aK2, C_Q1,C_Q2, usetidal, &rc,&drc_dr,&d2rc_dr);  
    double drc_dr         = rc_vec[1];
    double uc             = 1./rc;
    double uc2            = uc*uc;

    s_GS(r[i], rc, drc_dr, aK2, pr_star[i], pph[i], nu, chi1, chi2, X1, X2, c3, ggm);
        
    double GS             = ggm[2];
    double GSs            = ggm[3];
    
    double prstar2        = pr_star[i]*pr_star[i];
    double prstar4        = prstar2*prstar2;
    double pphi2          = pph[i]*pph[i];
    double Horbeff        = sqrt(prstar2+A[i]*(1. + pphi2*uc2 +  z3*prstar4*uc2));
    double Heff           = Horbeff + (GS*S + GSs*Sstar)*pph[i];
    double H              = sqrt( 1. + 2.*nu*(Heff - 1.) );
    
    Omg_orb[i]            = pph[i]*A[i]*uc2/(H*Horbeff);
    
  }
  
}


