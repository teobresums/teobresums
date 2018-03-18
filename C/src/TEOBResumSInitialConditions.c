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

/** Initial conditions calculation for non-spinning systems */
/* Post-post-circular initial data at separation r0
   r0       => relative separation
   pph0     => post-post-circular angular momentum
   pr       => post-circular radial momentum
   pprstar0 => post-circular r*-conjugate radial momentum
   j0       => circular angular momentum
   
   Three-step procedure
   1. Compute j0                         =>           circular ID, j!=0, pr =0
   2. From j0, compute pr*               =>      post circular ID, j!=0, pr !=0
   3. From pr* and j0, re-compute pph0   => post-post-circular ID, pph0!=j!=0, pr !=0
*/
void eob_dyn_ic(double r0, Dynamics *dyn, double y_init[])
{
  const double nu = dyn->nu;
  const double z3 = 2.0*nu*(4.0-3.0*nu);    

  /** Build a small grid */
#define N (6)
  const double dr = 1e-8;
    
  double r[2*N], dA[2*N], j[2*N], j2[2*N], djdr[2*N]; /** j:angular momentum */
  double E0[2*N], Omega_j[2*N];
  double Fphi[2*N], Ctmp[2*N], prstar[2*N], pr[2*N], pph[2*N], dprstardt[2*N];

  double A, B, d2A, dB;
  double r2, r3, j3;
  double H0eff, H0, psi, r_omega, v_phi, jhat, x;
  
  for (int i = 0; i < 2*N; i++) {
    
    r[i] = r0+(i-N+1)*dr;
    r2   = SQ(r[i]);
    r3   = r2*r[i];
    
    /** Compute metric  */
    eob_metric(r[i], dyn, &A, &B, &dA[i], &d2A, &dB);
    
    /** Angular momentum for circular orbit: circular ID  */
    j2[i]   =  r3*dA[i]/(2.*A-r[i]*dA[i]);
    j[i]    =  sqrt(j2[i]);
    j3      =  j2[i]*j[i];
    djdr[i] = -j3/r3*( 2.0 - 3.0*A/(r[i]*dA[i]) - A*d2A/(dA[i]*dA[i]) );
    
    /** For circular orbit at r0=r(N)  */
    H0eff      = sqrt(A*(1.0 + j2[i]/r2));                     /** effective Hamiltonian H_0^eff  */
    E0[i]      = sqrt(1.0 + 2.0*nu*(H0eff - 1.0) );            /** real Hamiltonian      H_0  */
    H0         = E0[i]/nu;                                     /** H_0/nu  */
    Omega_j[i] = A*j[i]/(nu*r2*H0*H0eff);                      /** Orbital frequency (from Hamilton's equation)  */
    psi        = 2.*(1.0 + 2.0*nu*(H0eff - 1.0))/(r2*dA[i]);   /** correction factor to the radius  */
    r_omega    = r[i]*cbrt(psi);                               /** EOB-corrected radius  */
    v_phi      = Omega_j[i]*r_omega;                           /** "corrected" azimuthal velocity such that Kepler's law is satisfied, r_omg^3 Omg_i^2 = 1  */
    x          = v_phi * v_phi;
    jhat       = j[i]/(r_omega*v_phi);                         /** Newton-normalized angular momentum  */
        
    Fphi[i] = eob_flx_Flux(x,Omega_j[i],r_omega,E0[i],H0eff,jhat,r[i], 0,0,dyn); 

    /** Radial momentum conjugate to r*: post-circular ID  */
    Ctmp[i]   = sqrt(B/A)*nu*H0*H0eff;
    prstar[i] = Ctmp[i]*Fphi[i]/djdr[i];
    
    /** Radial momentum conjugate to r  */
    pr[i] = prstar[i]*sqrt(B/A);
    
  }
    
  /** prstar by finite diff. */
  D0(prstar, dr, 2*N, dprstardt);
  int i = N-1;
  dprstardt[i] *= Fphi[i]/djdr[i];
  pph[i] = j[i]*sqrt(1. + 2.*Ctmp[i]/dA[i]*dprstardt[i] - z3*gsl_pow_int(prstar[i],4)/j2[i]);
  
  y_init[EOB_ID_RAD]    = r[N-1];
  y_init[EOB_ID_PPHI]   = pph[N-1];
  y_init[EOB_ID_PRSTAR] = prstar[N-1];
  y_init[EOB_ID_PR]     = pr[N-1];
  y_init[EOB_ID_J]      = j[N-1];
  y_init[EOB_ID_E0]     = E0[N-1];
  y_init[EOB_ID_OMGJ]   = Omega_j[N-1];
    
}

/** Initial conditions calculation for spinning systems */
void eob_dyn_ic_s(double r0, Dynamics *dyn, double y_init[])
{
  const double nu   = dyn->nu;
  const double chi1 = dyn->chi1;
  const double chi2 = dyn->chi2;
  const double S1   = dyn->S1;
  const double S2   = dyn->S2;
  const double c3   = dyn->cN3LO;
  const double X1   = dyn->X1;
  const double X2   = dyn->X2;
  const double a1   = dyn->a1;
  const double a2   = dyn->a2;
  const double aK2  = dyn->aK2;
  const double C_Q1 = dyn->C_Q1;
  const double C_Q2 = dyn->C_Q2;

  const double S  = S1 + S2;        
  const double Ss = X2*a1 + X1*a2;  
  const double z3 = 2.0*nu*(4.0-3.0*nu);
    
  /** Build a small grid */
#define N (6)
  const double dr = 1.e-8;

  double r[2*N], dA[2*N], j[2*N], j2[2*N], djdr[2*N]; /** j:angular momentum */
  double E0[2*N], Omega_j[2*N];
  double Fphi[2*N], Ctmp[2*N], prstar[2*N], pr[2*N], pph[2*N], dprstardt[2*N];
  double rc[2*N], drc_dr[2*N], d2rc_dr2[2*N], drc[2*N];
  double A[2*N],B[2*N],d2A[2*N],dB;
  double rorb;
  double pphorb;
  double r2, r3, j3;
  double H0eff, H0, psi, r_omega, v_phi, jhat, x;
  double ggm[14];

  //FIXME: redudant variables (ssame meaning, different name)
  double sqrtAbyB, uc, uc2;
  double pph2, Horbeff0, Heff0, one_H0, dHeff_dprstarbyprstar, dHeff_dpph, Heff, H, Horbeff;
  double ggm0[14];
  double GS_0, GSs_0, dGS_dr_0, dGSs_dr_0, dGSs_dpph_0, dGS_dprstarbyprstar_0, dGSs_dprstarbyprstar_0, GS, GSs, dGS_dr, dGSs_dr;
  double C0;
  double Omg;
  double Gtilde, dGtilde_dr, duc_dr, psic, r_omg;

  int i;

  for (i = 0; i < 2*N; i++) {
    r[i] = r0+(i-N+1)*dr;
        
    /** Compute metric  */
    eob_metric(r[i], dyn, &A[i],&B[i],&dA[i],&d2A[i],&dB);
    
    eob_dyn_s_get_rc(r[i], nu, a1, a2, aK2, C_Q1, C_Q2, dyn->use_tidal, &rc[i], &drc_dr[i], &d2rc_dr2[i]);
    
    /* Compute minimum of Heff0 using bisection method */
    rorb   = r[i];
    pphorb = rorb/sqrt(rorb-3.);
    pph[i] = eob_dyn_bisecHeff0_s(nu,chi1,chi2,X1,X2,c3, pphorb,rorb,A[i],dA[i],rc[i],drc_dr[i],aK2,S,Ss);
  }

  /** Post-circular initial conditions */
  
  /* pph by finite diff. */
  double dpph_dr[2*N];
  D0(pph, dr, 2*N, dpph_dr);
  
  for (i = 0; i < 2*N; i++) {
    
    sqrtAbyB = sqrt(A[i]/B[i]);
    uc  = 1./rc[i];
    uc2 = uc*uc;
    
    /* circular angular momentum */
    pph2 = pph[i]*pph[i];

    /* Orbital effective Hamiltonian */
    Horbeff0 = sqrt(A[i]*(1. + pph2*uc2));
    
    /* Compute gyro-gravitomagnetic coupling functions */
    eob_dyn_s_GS(r[i], rc[i], drc_dr[i], aK2, 0, pph[i], nu, chi1, chi2, X1, X2, c3, ggm0);
    GS_0                   = ggm0[2];
    GSs_0                  = ggm0[3];
    dGS_dr_0               = ggm0[6];
    dGSs_dr_0              = ggm0[7];
    dGSs_dpph_0            = ggm0[9];
    dGS_dprstarbyprstar_0  = ggm0[10];
    dGSs_dprstarbyprstar_0 = ggm0[11];

    /* Final effective Hamiltonian */
    Heff0 = (GS_0*S + GSs_0*Ss)*pph[i] + Horbeff0;

    /* Real Hamiltonian: beware that this is NOT divided by nu */
    H0     = sqrt( 1. + 2.*nu*(Heff0 - 1.));
    one_H0 = 1./H0;

    /* Get gyro-gravitomagnetic (derivative) functions */
    dHeff_dprstarbyprstar = pph[i]*(dGS_dprstarbyprstar_0*S + dGSs_dprstarbyprstar_0*Ss) + 1./Horbeff0;

    C0 = sqrtAbyB*one_H0*dHeff_dprstarbyprstar;

    /* Orbital frequency for circular orbit */
    dHeff_dpph = GS_0*S + (GSs_0 + pph[i]*dGSs_dpph_0)*Ss + pph[i]*A[i]*uc2/Horbeff0;
    Omg        = one_H0*dHeff_dpph;

    /* Flux */ 
    Gtilde     =  GS_0*S     + GSs_0*Ss;
    dGtilde_dr =  dGS_dr_0*S + dGSs_dr_0*Ss;
    duc_dr     = -uc2*drc_dr[i];
    psic       = (duc_dr + dGtilde_dr*rc[i]*sqrt(A[i]/pph2 + A[i]*uc2)/A[i])/(-0.5*dA[i]);
    r_omg      =  pow((pow(rc[i]*rc[i]*rc[i]*psic,-1./2)+Gtilde)*one_H0,-2./3.);
    v_phi      =  r_omg*Omg;
    x          =  v_phi*v_phi;
    jhat       =  pph[i]/(r_omg*v_phi);  /* Newton-normalized angular momentum */

    Fphi[i] = eob_flx_Flux_s(x, Omg, r_omg, H0, Heff0, jhat, r[i], 0., 0., dyn);
    prstar[i] = Fphi[i]/(dpph_dr[i]*C0);
    pr[i]     = prstar[i]* sqrt(B[i]/A[i]);

    j[i]       = pph[i];
    E0[i]      = H0;
    Omega_j[i] = Omg;

  }
  
#if (POSTPOSTCIRCULAR)  
  
  /** Post-post-circular initial data */
  
  double dpi1bydj[2*N];
  D0(prstar, dr, 2*N, dpi1bydj);

  double dpi1dt, prstar4, a,b,c;

  for (i = 0; i < 2*N; i++) {
    
    dpi1dt   = dpi1bydj[i]*Fphi[i];
    prstar4  = prstar[i]*prstar[i]*prstar[i]*prstar[i];
    
    /* Circular angular momentum */
    pph2     = pph[i]*pph[i];              

    /* Still circular, no pr* dependence here */
    Horbeff  = sqrt(A[i]*(1. + pph2*uc2)); 
    
    eob_dyn_s_GS(r[i], rc[i], drc_dr[i], aK2, 0, pph[i], nu, chi1, chi2, X1, X2, c3, ggm0);
    GS      = ggm0[2];
    GSs     = ggm0[3];
    dGS_dr  = ggm0[6];
    dGSs_dr = ggm0[7];

    /* Effective EOB energy */
    Heff     = (GS*S + GSs*Ss)*pph[i] + Horbeff;  

    /* Total EOB energy */
    H        = sqrt( 1. + 2.*nu*(Heff - 1.));     
    
    /* Setting up second order equation for the orbital angular momentum */       
    a = -sqrtAbyB*uc2/(2.*H*Horbeff)*(dA[i]  - 2.*A[i]*uc*drc[i]);                       
    b = -sqrtAbyB/H*(dGS_dr*S + dGSs_dr*Ss); 
    c = -dpi1dt - sqrtAbyB/(2.*H*Horbeff)*(dA[i] + z3*prstar4*uc2*(dA[i] - 2.*A[i]*uc*drc[i]));
    
    /* Fill out the array of the post-circular angular momentum */ 
    pph[i] = 0.5*(-b + sqrt(b*b-4*a*c))/a;      
  
  }

#endif

  y_init[EOB_ID_RAD]    = r[N-1];
  y_init[EOB_ID_PPHI]   = pph[N-1];
  y_init[EOB_ID_PRSTAR] = prstar[N-1];
  y_init[EOB_ID_PR]     = pr[N-1];
  y_init[EOB_ID_J]      = j[N-1];
  y_init[EOB_ID_E0]     = E0[N-1];
  y_init[EOB_ID_OMGJ]   = Omega_j[N-1];
  
}

/** Function for root finder: Derivative of the effective Hamiltonian */
struct DHeff0_tmp_params
{
  double rorb, A, dA, rc, drc_dr, ak2, S, Ss, nu, chi1, chi2, X1, X2, c3;
};

double eob_dyn_DHeff0(double x, void *params)
{
  
  struct DHeff0_tmp_params *p
    = (struct DHeff0_tmp_params *) params;
  
  double rorb   = p->rorb;
  double A      = p->A;
  double dA     = p->dA;
  double rc     = p->rc;
  double drc_dr = p->drc_dr;
  double ak2    = p->ak2;
  double S      = p->S;
  double Ss     = p->Ss;
  double nu     = p->nu;
  double chi1   = p->chi1;
  double chi2   = p->chi2;
  double X1     = p->X1;
  double X2     = p->X2;
  double c3     = p->c3;

  double ggm0[14];
  eob_dyn_s_GS(rorb, rc, drc_dr, ak2, 0., x, nu, chi1, chi2, X1, X2, c3, ggm0);
  double dGS_dr  = ggm0[6];
  double dGSs_dr = ggm0[7];
  
  double x2  = SQ(x);    
  double uc  = 1./rc;
  double uc2 = SQ(uc);
  double uc3 = uc2*uc;
  
  /* Orbital circular effective Hamiltonian */
  double Horbeff0 = sqrt(A*(1. + x2*uc2));    
  double dHeff_dr = x*(dGS_dr*S + dGSs_dr*Ss) + 1./(2.*Horbeff0)*( dA*(1. + x2*uc2) - 2.*A*uc3*drc_dr*x2);
  return dHeff_dr;
}

/** Root finder: Compute minimum of Heff0 */
double eob_dyn_bisecHeff0_s(double nu, double chi1, double chi2, double X1, double X2, double c3,
		       double pph, double rorb, double A, double dA, double rc, double drc_dr, double ak2, double S, double Ss)
{
  
#define max_iter (200)
#define tolerance (1e-14)

  int status;
  int iter = 0;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  
  double r;
  double x_lo = 0.5*pph, x_hi = 1.5*pph;

  gsl_function F;
  struct DHeff0_tmp_params p = {rorb,A,dA,rc,drc_dr,ak2,S,Ss,nu,chi1,chi2,X1,X2,c3};
  F.function = &eob_dyn_DHeff0;
  
  F.params = &p;
  T = gsl_root_fsolver_bisection;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r      = gsl_root_fsolver_root (s);
      x_lo   = gsl_root_fsolver_x_lower (s);
      x_hi   = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi, 0, tolerance);
    }
  while (status == GSL_CONTINUE && iter < max_iter);
  gsl_root_fsolver_free (s);
  
  return r;
}
