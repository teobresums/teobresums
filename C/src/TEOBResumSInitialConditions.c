/** \file TEOBResumSInitialConditions.c
 *  \brief Initial conditions for the TEOBResumS model
 * 
 *  This file contains the routines employed to compute the EOB initial coditions
 *  according to the user's choice of initial frequency or initial separation
 *  and the intrinsic parameters of the source.
 */

#include "TEOBResumS.h"

/**
 * Function: eob_dyn_ic_circ
 * -------------------------
 *   Initial conditions calculation for non-spinning, quasi-circ systems
 *   Post-post-circular initial data at separation r0
 *   
 *   Three-step procedure
 *    1. Compute j0                         =>           circular ID, j!=0, pr =0
 *    2. From j0, compute pr*               =>      post circular ID, j!=0, pr !=0
 *    3. From pr* and j0, re-compute pph0   => post-post-circular ID, pph0!=j!=0, pr !=0
 * 
 *   @param[in] r0:  initial separation
 *   @param[in] dyn: Dynamics
 *   @param[out] y_init: initial data
 *   @return status: 0 if successful
*/
int eob_dyn_ic_circ(double r0, Dynamics *dyn, double y_init[])
{

  int status = 0;
  const double nu = EOBPars->nu;
  const double z3 = 2.0*nu*(4.0-3.0*nu);    

  /* Build a small grid */
#define N (6)
  const double dr = 1e-10;
    
  double r[2*N], dA[2*N], j[2*N], j2[2*N], djdr[2*N]; /** j:angular momentum */
  double E0[2*N], Omega_j[2*N];
  double Fphi[2*N], Ctmp[2*N], prstar[2*N], pr[2*N], pph[2*N], dprstardr[2*N];

  double A, B, Q, d2A, dB, pl_hold;
  double r2, r3, j3;
  double H0eff, H0, psi, r_omega, v_phi, jhat, x, dprstardt;
  
  for (int i = 0; i < 2*N; i++) {
    r[i] = r0+(i-N+1)*dr;
    r2   = SQ(r[i]);
    r3   = r2*r[i];
    
    /* Compute metric  */
    eob_metric(r[i], 0., dyn, &A, &B, &dA[i], &d2A, &dB, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    
    /* Angular momentum for circular orbit: circular ID  */
    j2[i]   =  r3*dA[i]/(2.*A-r[i]*dA[i]);
    j[i]    =  sqrt(j2[i]);
    j3      =  j2[i]*j[i];
    djdr[i] = -j3/r3*( 2.0 - 3.0*A/(r[i]*dA[i]) - A*d2A/(dA[i]*dA[i]) );
    
    /* For circular orbit at r0=r(N)  */
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

    /* Radial momentum conjugate to r*: post-circular ID  */
    Ctmp[i]   = sqrt(B/A)*nu*H0*H0eff;
    prstar[i] = Ctmp[i]*Fphi[i]/djdr[i];
    
    /* Radial momentum conjugate to r  */
    pr[i] = prstar[i]*sqrt(B/A);
    
  }
    
  /* prstar by finite diff. */
  D0(prstar, dr, 2*N, dprstardr);
  
  int i = N-1;
  dprstardt = dprstardr[i] * Fphi[i]/djdr[i];
  // compute Q after having evaluated pr_star
  eob_metric(r[i], prstar[i], dyn, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &Q, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
  pph[i] = j[i]*sqrt(1. + 2.*Ctmp[i]/dA[i]*dprstardt - (Q*r2)/j2[i]);
  //printf("%d %.16e %.16e %.16e %.16e %.16e %.16e\n",i,r[i],Fphi[i],djdr[i],dprstardr[i],dprstardt,pph[i]);
  //printf("%d %.16e %.16e %.16e %.16e %.16e %.16e\n",i,r[i],Fphi[i],djdr[i],dprstardr[i],dprstardt,pph[i]);


  y_init[EOB_ID_RAD]    = r[N-1];
  y_init[EOB_ID_PHI]    = 0.;
  y_init[EOB_ID_PPHI]   = pph[N-1];
  y_init[EOB_ID_PRSTAR] = prstar[N-1];
  y_init[EOB_ID_PR]     = pr[N-1];
  y_init[EOB_ID_J]      = j[N-1];
  y_init[EOB_ID_E0]     = E0[N-1];
  y_init[EOB_ID_OMGJ]   = Omega_j[N-1];
  
  return status;
}
/**
 * Function: eob_dyn_ic_circ_s
 * ----------------------
 *   Initial conditions calculation for spinning quasi-circ systems
 *   Post-post-circular initial data at separation r0
 *   
 *   Three-step procedure
 *    1. Compute j0                         =>           circular ID, j!=0, pr =0
 *    2. From j0, compute pr*               =>      post circular ID, j!=0, pr !=0
 *    3. From pr* and j0, re-compute pph0   => post-post-circular ID, pph0!=j!=0, pr !=0
 * 
 *   @param[in] r0:  initial separation
 *   @param[in] dyn: Dynamics
 *   @param[out] y_init: initial data
 *   @return status: 0 if successful
*/
int eob_dyn_ic_circ_s(double r0, Dynamics *dyn, double y_init[])
{
  int status = 0;
  const double nu   = EOBPars->nu;
  const double chi1 = EOBPars->chi1;
  const double chi2 = EOBPars->chi2;
  const double S1   = EOBPars->S1;
  const double S2   = EOBPars->S2;
  const double c3   = EOBPars->cN3LO;
  const double X1   = EOBPars->X1;
  const double X2   = EOBPars->X2;
  const double a1   = EOBPars->a1;
  const double a2   = EOBPars->a2;
  const double aK2  = EOBPars->aK2;
  const double C_Q1 = EOBPars->C_Q1;
  const double C_Q2 = EOBPars->C_Q2;
  const double C_Oct1 = EOBPars->C_Oct1;
  const double C_Oct2 = EOBPars->C_Oct2;
  const double C_Hex1 = EOBPars->C_Hex1;
  const double C_Hex2 = EOBPars->C_Hex2;

  const double S  = S1 + S2;        
  const double Ss = X2*a1 + X1*a2;  
  const double z3 = 2.0*nu*(4.0-3.0*nu);

  /* Build a small grid */
#define N (6)
  const double dr = 1e-4; /* do not change this */

  double r[2*N], dA[2*N], j[2*N]; /* j:angular momentum */
  double E0[2*N], Omega_j[2*N];
  double Fphi[2*N], Ctmp[2*N], prstar[2*N], pr[2*N], pph[2*N];
  double rc[2*N], drc_dr[2*N], d2rc_dr2[2*N]; //, drc[2*N];
  double A[2*N],B[2*N],d2A[2*N],dB, sqrtAbyB, Q, pl_hold;
  double pphorb, uc, uc2, psic, r_omg, v_phi, jhat, x, Omg;
  double H0eff, H0, Horbeff0, Heff0, one_H0, dHeff_dprstarbyprstar, dHeff_dpph, Heff, H, Horbeff;
  double ggm0[26], GS_0, GSs_0, dGS_dr_0, dGSs_dr_0, dGSs_dpph_0, dGS_dprstarbyprstar_0, dGSs_dprstarbyprstar_0, GS, GSs, dGS_dr, dGSs_dr;
  double C0;
  double Gtilde, dGtilde_dr, duc_dr;

  int i;
  for (i = 0; i < 2*N; i++) {
    r[i] = r0+(i-N+1)*dr;

    /* Compute metric  */
    eob_metric_s(r[i], 0., dyn, &A[i], &B[i], &dA[i], &d2A[i], &dB, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    
    /* Compute initial adiabatic angular momentum  */
    pphorb = r[i]/sqrt(r[i]-3.);
    eob_dyn_s_get_rc(r[i], nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, EOBPars->use_tidal, &rc[i], &drc_dr[i], &d2rc_dr2[i]);
    pph[i] = eob_dyn_bisecHeff0_s(nu,chi1,chi2,X1,X2,c3, pphorb,r[i],A[i],dA[i],rc[i],drc_dr[i],aK2,S,Ss);

  }

  /* Post-circular initial conditions */
  
  /* pph by finite diff. */
  double dpph_dr[2*N]; 
  D0(pph, dr, 2*N, dpph_dr); 

  for (i = 0; i < 2*N; i++) {
    
    sqrtAbyB = sqrt(A[i]/B[i]);
    uc  = 1./rc[i];
    uc2 = uc*uc;
        
    /* Orbital effective Hamiltonian */
    Horbeff0 = sqrt(A[i]*(1. + SQ(pph[i])*uc2));
    
    /* Compute gyro-gravitomagnetic coupling functions */
    eob_dyn_s_GS(r[i], rc[i], drc_dr[i], 0., aK2, 0, pph[i], nu, chi1, chi2, X1, X2, c3, ggm0);
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
    psic       = (duc_dr + dGtilde_dr*rc[i]*sqrt(A[i]/(SQ(pph[i])) + A[i]*uc2)/A[i])/(-0.5*dA[i]);
    r_omg      =  pow((pow(rc[i]*rc[i]*rc[i]*psic,-1./2)+Gtilde)*one_H0,-2./3.);
    v_phi      =  r_omg*Omg;
    x          =  v_phi*v_phi;
    jhat       =  pph[i]/(r_omg*v_phi);  /* Newton-normalized angular momentum */

    Fphi[i]    = eob_flx_Flux_s(x, Omg, r_omg, H0, Heff0, jhat, r[i], 0., 0., dyn);
    prstar[i]  = Fphi[i]/(dpph_dr[i]*C0);
    pr[i]      = prstar[i]/sqrtAbyB; 

    j[i]       = pph[i];
    E0[i]      = H0;
    Omega_j[i] = Omg;

  }
  
#if (POSTPOSTCIRCULAR)  
  
  /* Post-post-circular initial data */
  
  double dpi1bydj, dprstardr[2*N],djdr[2*N];
  D0(prstar, dr, 2*N, dprstardr); 
  D0(j, dr, 2*N, djdr); 

  double dpi1dt, prstar4, a,b,c;

  //for (i = 0; i < 2*N; i++) { //No need of the loop here
  i = N-1;
  
  sqrtAbyB = sqrt(A[i]/B[i]);
  uc  = 1./rc[i];
  uc2 = uc*uc;
  
  dpi1bydj = dprstardr[i]/djdr[i];
  dpi1dt   = dpi1bydj*Fphi[i];
  prstar4  = prstar[i]*prstar[i]*prstar[i]*prstar[i];
  
  /* Still circular, no pr* dependence here */
  Horbeff  = sqrt(A[i]*(1. + SQ(pph[i])*uc2)); 
  eob_dyn_s_GS(r[i], rc[i], drc_dr[i], 0., aK2, 0, pph[i], nu, chi1, chi2, X1, X2, c3, ggm0);
  GS      = ggm0[2];
  GSs     = ggm0[3];
  dGS_dr  = ggm0[6];
  dGSs_dr = ggm0[7];
  
  /* Effective EOB energy */
  Heff     = (GS*S + GSs*Ss)*pph[i] + Horbeff;  
  
  /* Total EOB energy */
  H        = sqrt( 1. + 2.*nu*(Heff - 1.));  

  // compute Q with prstar value
  eob_metric_s(r[i], prstar[i], dyn, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &Q, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);     
  
  /* Setting up second order equation for the orbital angular momentum */       
  a = -sqrtAbyB*uc2/(2.*H*Horbeff)*(dA[i]  - 2.*A[i]*uc*drc_dr[i]);                       
  b = -sqrtAbyB/H*(dGS_dr*S + dGSs_dr*Ss); 
  c = -dpi1dt - sqrtAbyB/(2.*H*Horbeff)*(dA[i] + Q*(dA[i] - 2.*A[i]*uc*drc_dr[i]));
  
  /* Fill out the array of the post-circular angular momentum */ 
  pph[i] = 0.5*(-b + sqrt(b*b-4*a*c))/a;      
  
  //  }
  
#endif

  y_init[EOB_ID_RAD]    = r[N-1];
  y_init[EOB_ID_PHI]    = 0.;
  y_init[EOB_ID_PPHI]   = pph[N-1];
  y_init[EOB_ID_PRSTAR] = prstar[N-1];
  y_init[EOB_ID_PR]     = pr[N-1];
  y_init[EOB_ID_J]      = j[N-1];
  y_init[EOB_ID_E0]     = E0[N-1];
  y_init[EOB_ID_OMGJ]   = Omega_j[N-1];
  return status;
}

/** Initial conditions calculation for eccentric systems */

/**
 * Function: eob_dyn_ic_ecc
 * ------------------------
 *   Initial conditions calculation for eccentric systems
 *   Adiabatic level, assumes pr*0 = 0
 * 
 *   @param[in] r0:  initial separation
 *   @param[in] dyn: Dynamics
 *   @param[out] y_init: initial data
 *   @return status: 0 if successful
*/
int eob_dyn_ic_ecc(double r0, Dynamics *dyn, double y_init[])
{

  int status = 0;

  const double nu = EOBPars -> nu;
  const double chi1 = EOBPars -> chi1;
  const double chi2 = EOBPars -> chi2;
  const double S1 = EOBPars -> S1;
  const double S2 = EOBPars -> S2;
  const double c3 = EOBPars -> cN3LO;
  const double X1 = EOBPars -> X1;
  const double X2 = EOBPars -> X2;
  const double a1 = EOBPars -> a1;
  const double a2 = EOBPars -> a2;
  const double aK2 = EOBPars -> aK2;
  const double ecc = EOBPars -> ecc;
  const double C_Q1 = EOBPars -> C_Q1;
  const double C_Q2 = EOBPars -> C_Q2;
  const double C_Oct1 = EOBPars -> C_Oct1;
  const double C_Oct2 = EOBPars -> C_Oct2;
  const double C_Hex1 = EOBPars -> C_Hex1;
  const double C_Hex2 = EOBPars -> C_Hex2;
  const double S = S1 + S2;
  const double Sstar = X2*a1 + X1*a2;

  const int usetidal = EOBPars -> use_tidal;
  const int usespins = EOBPars -> use_spins;

  double r1, A1, rc1, B1, ggm1[26], G1;
  double r2, A2, rc2, B2, ggm2[26], G2;
  double pl_hold, A12, B12, DA, DB, DG;
  double j0, j02, Heff_orb1, Heff1, H1, dHeff1_dj0, omg_orb1;
  
  r1 = r0/(1-ecc);
  r2 = r0/(1+ecc);

  /* Computing metric, centrifugal radius and ggm functions*/
  if(usespins) {
    eob_metric_s(r1, 0., dyn, &A1, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    eob_dyn_s_get_rc(r1, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc1, &pl_hold, &pl_hold, &pl_hold);
    
    eob_dyn_s_GS(r1, rc1, 0.0, 0.0, aK2, 0.0, 0.0, nu, chi1, chi2, X1, X2, c3, ggm1);
    G1     = ggm1[2]*S + ggm1[3]*Sstar;    // tildeG = GS*S+GSs*Ss

    eob_metric_s(r2, 0., dyn, &A2, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    eob_dyn_s_get_rc(r2, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc2, &pl_hold, &pl_hold);
    
    eob_dyn_s_GS(r2, rc2, 0.0, 0.0, aK2, 0.0, 0.0, nu, chi1, chi2, X1, X2, c3, ggm2);
    G2     = ggm2[2]*S + ggm2[3]*Sstar;    
  } else {
    eob_metric(r1, 0., dyn, &A1, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    rc1 = r1;   //Nonspinning case: rc = r; G = 0;  
    G1  = 0.0;

    eob_metric(r2, 0., dyn, &A2, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    rc2 = r2;
    G2  = 0.0;
  }
  /* Adiabatic angular momentum */
  j0  = eob_dyn_j0(r0, dyn);
  j02 = SQ(j0);

  /* Adiabatic energy */
  Heff_orb1 = sqrt(A1*(1. + j02/SQ(rc1)));
  Heff1     = Heff_orb1 + j0*G1;
  H1        = sqrt(1. + 2.*nu*(Heff1 - 1.))/nu;
    
  /* Orbital frequency */
  dHeff1_dj0 = G1 + j0*A1/(Heff_orb1*SQ(rc1));
  omg_orb1    = dHeff1_dj0/nu/H1;
  
  y_init[EOB_ID_RAD]    = r1;
  y_init[EOB_ID_PHI]    = 0.;
  y_init[EOB_ID_PPHI]   = j0;
  y_init[EOB_ID_PRSTAR] = 0.;
  y_init[EOB_ID_PR]     = 0.;
  y_init[EOB_ID_J]      = j0;
  y_init[EOB_ID_E0]     = H1*nu;
  y_init[EOB_ID_OMGJ]   = omg_orb1;

  return status;

}


/**
 * Function: eob_dyn_j0
 * ------------------------
 *    Compute the conservative angular momentum j0 for a 
 *    given initial semilatus rectum r0. Discriminate from 
 *    circular and eccentric orbits att ecc=1e-10, where
 *    the eob_dyn_ecc_j0 function becomes unstable.
 * 
 *    @param[in] r0:  semilatus rectum
 *    @param[in] dyn: Dynamics
 * 
 *    @return j0: angular momentum
*/
double eob_dyn_j0(double r0, Dynamics *dyn)
{
  const double ecc    = EOBPars->ecc;
  double j0;

  if(ecc > 1e-10){
    j0 = eob_dyn_ecc_j0(r0, dyn);
  } else {
    j0 = eob_dyn_circ_j0(r0, dyn);
  }
  return j0;
}

/**
 * Function: eob_dyn_ecc_j0
 * ------------------------
 *    Compute the conservative angular momentum of an elliptic orbit
 *    from Energy conservation at periastron and apastron
 * 
 *    @param[in] r0:  initial separation
 *    @param[in] dyn: Dynamics
 * 
 *    @return j0: angular momentum
*/
double eob_dyn_ecc_j0(double r0, Dynamics *dyn)
{
  const double nu   = EOBPars->nu;
  const double chi1 = EOBPars->chi1;
  const double chi2 = EOBPars->chi2;
  const double S1   = EOBPars->S1;
  const double S2   = EOBPars->S2;
  const double c3   = EOBPars->cN3LO;
  const double X1   = EOBPars->X1;
  const double X2   = EOBPars->X2;
  const double a1   = EOBPars->a1;
  const double a2   = EOBPars->a2;
  const double aK2  = EOBPars->aK2;
  const double ecc    = EOBPars->ecc;
  const double C_Q1   = EOBPars->C_Q1;
  const double C_Q2   = EOBPars->C_Q2;
  const double C_Oct1 = EOBPars->C_Oct1;
  const double C_Oct2 = EOBPars->C_Oct2;
  const double C_Hex1 = EOBPars->C_Hex1;
  const double C_Hex2 = EOBPars->C_Hex2;
  const double S     = S1 + S2;        
  const double Sstar = X2*a1 + X1*a2;  

  const int usespins = EOBPars->use_spins;
  const int usetidal = EOBPars->use_tidal;  

  double r1, A1, rc1, B1, ggm1[26], G1;
  double r2, A2, rc2, B2, ggm2[26], G2;
  double pl_hold, A12, B12, DA, DB, DG;
  double j0, j02, Heff_orb1, Heff1, H1, dHeff1_dj0, omg_orb1;

  r1 = r0/(1-ecc); // apastron
  r2 = r0/(1+ecc); // periastron

  /* Computing metric, centrifugal radius and ggm functions*/
  if(usespins) {
    eob_metric_s(r1, 0., dyn, &A1, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    eob_dyn_s_get_rc(r1, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc1, &pl_hold, &pl_hold, &pl_hold);
    
    eob_dyn_s_GS(r1, rc1, 0.0, 0.0, aK2, 0.0, 0.0, nu, chi1, chi2, X1, X2, c3, ggm1);
    G1     = ggm1[2]*S + ggm1[3]*Sstar;    // tildeG = GS*S+GSs*Ss

    eob_metric_s(r2, 0., dyn, &A2, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    eob_dyn_s_get_rc(r2, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc2, &pl_hold, &pl_hold);
    
    eob_dyn_s_GS(r2, rc2, 0.0, 0.0, aK2, 0.0, 0.0, nu, chi1, chi2, X1, X2, c3, ggm2);
    G2     = ggm2[2]*S + ggm2[3]*Sstar;    
  } else {
    eob_metric(r1, 0., dyn, &A1, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    rc1 = r1;   //Nonspinning case: rc = r; G = 0;  
    G1  = 0.0;

    eob_metric(r2, 0., dyn, &A2, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    rc2 = r2;
    G2  = 0.0;
  }
  B1 = A1/SQ(rc1);
  B2 = A2/SQ(rc2);

  A12 = A1 + A2;
  DA  = A1 - A2;
  B12 = B1 + B2;
  DB  = B1 - B2;
  DG  = G1 - G2;

  /* Angular momentum */
  j02       = (A12*SQ(DG) - DA*DB + DG*sqrt(4.*A1*A2*SQ(DG) + 2.*DA*(B12*DA - A12*DB)))/(SQ(DB) - 2.*B12*SQ(DG) + SQ(SQ(DG)));
  j0        = sqrt(j02);

  return j0;

}

/**
 * Function: eob_dyn_circ_j0
 * -------------------------
 *    Compute the conservative angular momentum of an circular orbit
 *    from analytically solving dH_dr = 0
 * 
 *    @param[in] r0:  initial separation
 *    @param[in] dyn: Dynamics
 * 
 *    @return j0: circular angular momentum
*/
double eob_dyn_circ_j0(double r0, Dynamics *dyn)
{
  const double nu    = EOBPars->nu;
  const double X1    = EOBPars->X1;
  const double X2    = EOBPars->X2;
  const double chi1  = EOBPars->chi1;
  const double chi2  = EOBPars->chi2;
  const double a1    = EOBPars->a1;
  const double a2    = EOBPars->a2;
  const double aK2   = EOBPars->aK2;
  const double S     = EOBPars->S;
  const double Sstar = EOBPars->Sstar;
  const double c3    = EOBPars->cN3LO;
  const double C_Q1  = EOBPars->C_Q1;
  const double C_Q2  = EOBPars->C_Q2;
  const double C_Oct1 = EOBPars->C_Oct1;
  const double C_Oct2 = EOBPars->C_Oct2;
  const double C_Hex1 = EOBPars->C_Hex1;
  const double C_Hex2 = EOBPars->C_Hex2;

  const int usetidal = EOBPars->use_tidal;  
  const int usespins = EOBPars->use_spins;

  double A,B,dA,rc,drc_dr,G,dG_dr,uc,uc2,dAuc2_dr,j02,j0,H,Heff,Heff_orb,dHeff_dj0,omg_orb;
  double pl_hold,a_coeff,b_coeff,c_coeff,Delta,sol_p,sol_m;
  double ggm[26];

  /* Computing metric, centrifugal radius and ggm functions*/
  if(usespins) {
    eob_metric_s(r0, 0., dyn, &A, &B, &dA, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    eob_dyn_s_get_rc(r0, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc, &drc_dr, &pl_hold);
    eob_dyn_s_GS(r0, rc, drc_dr, 0., aK2, 0.0, 0.0, nu, chi1, chi2, X1, X2, c3, ggm);
    G     = ggm[2]*S + ggm[3]*Sstar;    // tildeG = GS*S+GSs*Ss
    dG_dr = ggm[6]*S + ggm[7]*Sstar;
  } else {
    eob_metric(r0, 0., dyn, &A, &B, &dA, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    rc     = r0;   //Nonspinning case: rc = r; G = 0;
    drc_dr = 1.;  
    G      = 0.;
    dG_dr  = 0.;
  }

  /* Auxiliary variables*/
  uc       = 1./rc;
  uc2      = uc*uc;
  dAuc2_dr = uc2*(dA-2*A*uc*drc_dr);

  /* Circular angular momentum */
  if (usespins) {

    // Quadratic equation a*x^2+b*x+c=0 
    a_coeff = SQ(dAuc2_dr)  - 4*A*uc2*SQ(dG_dr);
    b_coeff = 2*dA*dAuc2_dr - 4*A*SQ(dG_dr);
    c_coeff = SQ(dA);
      
    Delta = SQ(b_coeff) - 4*a_coeff*c_coeff;
      
    if (Delta<0.)  
      Delta=0;             // dG_dr=0 -> Set Delta=0 to avoid num. errors
      
    sol_p   = (-b_coeff + sqrt(Delta))/(2*a_coeff); 
    sol_m   = (-b_coeff - sqrt(Delta))/(2*a_coeff);
      
    if (dG_dr > 0)
      j02 = sol_p;
    else
      j02 = sol_m;

  } else {
    // Linear equation a*x+b=0
    a_coeff = dAuc2_dr;
    b_coeff = dA;    
    j02 = -b_coeff/a_coeff;    
  }
  
  j0 = sqrt(j02);
  return j0;

}

/**
 * Function: eob_dyn_ic_ecc_PA
 * ---------------------------
 *   Initial conditions calculation for eccentric systems
 *   Post-adiabatic level, computes pr*0 != 0 (but very small)
 *   The routine employed is similar to the quasi-circular one.
 *   Formally, this is not correct, but it allows a good approximation
 *   in the quasi-circular limit.
 * 
 *   @param[in] r0:  initial separation
 *   @param[in] dyn: Dynamics
 *   @param[out] y_init: initial data
 *   @return status: 0 if successful
*/
int eob_dyn_ic_ecc_PA(double r0, Dynamics *dyn, double y_init[])
{

  int status = 0;
  const double nu   = EOBPars->nu;
  const double chi1 = EOBPars->chi1;
  const double chi2 = EOBPars->chi2;
  const double S1   = EOBPars->S1;
  const double S2   = EOBPars->S2;
  const double c3   = EOBPars->cN3LO;
  const double X1   = EOBPars->X1;
  const double X2   = EOBPars->X2;
  const double a1   = EOBPars->a1;
  const double a2   = EOBPars->a2;
  const double aK2  = EOBPars->aK2;
  const double ecc    = EOBPars->ecc;
  const double C_Q1   = EOBPars->C_Q1;
  const double C_Q2   = EOBPars->C_Q2;
  const double C_Oct1 = EOBPars->C_Oct1;
  const double C_Oct2 = EOBPars->C_Oct2;
  const double C_Hex1 = EOBPars->C_Hex1;
  const double C_Hex2 = EOBPars->C_Hex2;
  const double S      = S1 + S2;        
  const double Sstar  = X2*a1 + X1*a2;  

  const int use_tidal = EOBPars->use_tidal;
  
  /*----------------*/
  /* 1. compute j0  */
  /*----------------*/

  double r1 = r0/(1-ecc);
  double r2 = r0/(1+ecc);
  double j0 = eob_dyn_j0(r0, dyn);
  
  if (DEBUG) printf("0PA:\np0 = %.8f\npph0 = %.8f\nr0 = %.8f\n", r0, j0, r1);

  /*----------------------*/
  /* 2. compute pr* != 0  */
  /*----------------------*/

  /** Build a small grid */ 
#define N (6)
  double dpph_dr[2*N], p[2*N], r[2*N], pph[2*N]; 
  const double dp = 1e-4;   /* do not change this */
  const double dr = dp/(1. - ecc);
  for (int i=0; i< 2*N; i++) {
    p[i]   = r0+(i-N+1)*dp; /* grid of semilatus rectum */
    r[i]   = p[i]/(1.-ecc); /* grid of r = p/(1-e) */
    pph[i] = eob_dyn_j0(p[i], dyn);
  }

  /* dpph_dr by finite diff. */
  D0(pph, dr, 2*N, dpph_dr);
  int idx = N-1;

  double E0, Omega_j;
  double Fphi, Ctmp, prstar, pr;
  double rc, drc_dr, d2rc_dr2;
  double A,B,dA,d2A,dB, sqrtAbyB, pl_hold;
  double pphorb, uc, uc2, psic, r_omg, v_phi, jhat, x, Omg;
  double H0eff, H0, Horbeff0, Heff0, one_H0, dHeff_dprstarbyprstar, dHeff_dpph, Heff, H, Horbeff;
  double ggm0[26], GS_0, GSs_0, dGS_dr_0, dGSs_dr_0, dGSs_dpph_0, dGS_dprstarbyprstar_0, dGSs_dprstarbyprstar_0, GS, GSs, dGS_dr, dGSs_dr;
  double C0;
  double Gtilde, dGtilde_dr, duc_dr;

  eob_metric_s(r[idx], 0., dyn, &A, &B, &dA, &d2A, &dB, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
  eob_dyn_s_get_rc(r[idx], nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, use_tidal, &rc, &drc_dr, &d2rc_dr2);    
  sqrtAbyB = sqrt(A/B);
  uc  = 1./rc;
  uc2 = uc*uc;
  /* Orbital effective Hamiltonian */
  Horbeff0 = sqrt(A*(1. + SQ(pph[idx])*uc2));
  /* Compute gyro-gravitomagnetic coupling functions */
  eob_dyn_s_GS(r[idx], rc, drc_dr, 0.0, aK2, 0, pph[idx], nu, chi1, chi2, X1, X2, c3, ggm0);
  GS_0                   = ggm0[2];
  GSs_0                  = ggm0[3];
  dGS_dr_0               = ggm0[6];
  dGSs_dr_0              = ggm0[7];
  dGSs_dpph_0            = ggm0[9];
  dGS_dprstarbyprstar_0  = ggm0[10];
  dGSs_dprstarbyprstar_0 = ggm0[11];
  /* Final effective Hamiltonian */
  Heff0 = (GS_0*S + GSs_0*Sstar)*pph[idx] + Horbeff0;

  /* Real Hamiltonian: beware that this is NOT divided by nu */
  H0     = sqrt( 1. + 2.*nu*(Heff0 - 1.));
  one_H0 = 1./H0;

  /* Get gyro-gravitomagnetic (derivative) functions */
  dHeff_dprstarbyprstar = pph[idx]*(dGS_dprstarbyprstar_0*S + dGSs_dprstarbyprstar_0*Sstar) + 1./Horbeff0;

  C0         = sqrtAbyB*one_H0*dHeff_dprstarbyprstar;
  if (DEBUG) printf("sqrtAbyB = %.8f\none_H0 = %.8f\ndHeff_dprstarbyprstar = %.8f\nC0 = %.8f\n", sqrtAbyB,one_H0,dHeff_dprstarbyprstar,C0);

  /* Orbital frequency at apastron */
  dHeff_dpph = GS_0*S + (GSs_0 + pph[idx]*dGSs_dpph_0)*Sstar + pph[idx]*A*uc2/Horbeff0;
  Omg        = one_H0*dHeff_dpph;

  /* Flux */ 
  Gtilde     =  GS_0*S     + GSs_0*Sstar;
  dGtilde_dr =  dGS_dr_0*S + dGSs_dr_0*Sstar;
  duc_dr     = -uc2*drc_dr;
  psic       = (duc_dr + dGtilde_dr*rc*sqrt(A/(SQ(pph[idx])) + A*uc2)/A)/(-0.5*dA);
  r_omg      =  pow((pow(rc*rc*rc*psic,-1./2)+Gtilde)*one_H0,-2./3.);
  v_phi      =  r_omg*Omg;
  x          =  v_phi*v_phi;
  jhat       =  pph[idx]/(r_omg*v_phi);  /* Newton-normalized angular momentum */

  //Fphi    = eob_flx_Flux_s(x, Omg, r_omg, H0, Heff0, jhat, r[idx], 0.0, 0.0, dyn);
  eob_flx_Flux_ecc(x, Omg, r_omg, H0, Heff0, jhat, r[idx], 0.0, pph[idx], 0.0, 0.0, &Fphi, &pl_hold, dyn);
  prstar  = Fphi/(dpph_dr[idx]*C0);
  pr      = prstar/sqrtAbyB; 
  E0      = H0;
  Omega_j = Omg;

  if (DEBUG) printf("\n1PA:\np0 = %.8f\npph0 = %.8f\nFphi = %.8f\npr0 = %.8f\npr*0 = %.8f\nE0 = %.8f\nomega0 = %.8f\n", r0, j0, Fphi, pr, prstar, E0, Omg);

  // Update ICs
  y_init[EOB_ID_RAD]    = r1;
  y_init[EOB_ID_PHI]    = 0.;
  y_init[EOB_ID_PPHI]   = pph[idx];
  y_init[EOB_ID_PRSTAR] = prstar;
  y_init[EOB_ID_PR]     = pr;
  y_init[EOB_ID_J]      = j0;
  y_init[EOB_ID_E0]     = E0;
  y_init[EOB_ID_OMGJ]   = Omega_j;
  
  return status;

}

/**
 * Function: eob_dyn_ic_ecc_ma
 * ---------------------------
 *   Initial conditions calculation for eccentric systems
 *   Adiabatic level, computes pr*0 != 0 at a specified anomaly.
 *   Performs a bidimensional root-finding to find the correct
 *   initial data.
 *
 *   @param[in] r0_kepl:  initial estimate of semilatus rectum
 *   @param[in] dyn: Dynamics
 *   @param[out] y_init: initial data
 *   @return status: 0 if successful
*/
int eob_dyn_ic_ecc_ma(double r0_kepl, Dynamics *dyn, double y_init[])
{
  
  int status = 0;

  const double f0 = EOBPars -> f0;
  const double nu = EOBPars -> nu;
  const double chi1 = EOBPars -> chi1;
  const double chi2 = EOBPars -> chi2;
  const double S1 = EOBPars -> S1;
  const double S2 = EOBPars -> S2;
  const double c3 = EOBPars -> cN3LO;
  const double X1 = EOBPars -> X1;
  const double X2 = EOBPars -> X2;
  const double a1 = EOBPars -> a1;
  const double a2 = EOBPars -> a2;
  const double aK2 = EOBPars -> aK2;
  const double C_Q1 = EOBPars -> C_Q1;
  const double C_Q2 = EOBPars -> C_Q2;
  const double C_Oct1 = EOBPars -> C_Oct1;
  const double C_Oct2 = EOBPars -> C_Oct2;
  const double C_Hex1 = EOBPars -> C_Hex1;
  const double C_Hex2 = EOBPars -> C_Hex2;
  const double S = S1 + S2;
  const double Sstar = X2*a1 + X1*a2;
  const double ecc = EOBPars -> ecc;
  const double zeta = EOBPars -> anomaly;

  const int usetidal = EOBPars -> use_tidal;
  const int usespins = EOBPars -> use_spins;

  double j0    = eob_dyn_j0(r0_kepl, dyn); // initial guess for j0
  double pr0PN = 0.01;// = ecc/j0*sin(zeta+1.);
  
  double r0, pr0abs;

  double omg_orb0 = Pi*f0;
  if (EOBPars->ecc_freq == ECCFREQ_ORBAVGD){
    /* Assume that the user gave as input the orbit averaged 
    frequency. Transform it into an average frequency between apastron
    and periastron and overwrite omg_orb0 */
    if (DEBUG) printf("Orbit-averaged ICs:\nomg_bar = %.8f\n", omg_orb0);
    double omg_orb0_p = eob_dyn_omg_from_omgbar(omg_orb0, 0., ecc);
    double omg_orb0_m = eob_dyn_omg_from_omgbar(omg_orb0, Pi, ecc);
    omg_orb0 = 0.5*(omg_orb0_p + omg_orb0_m);
    if (DEBUG) {
      printf("omg_orb0  = %.8f\n",  omg_orb0);
      printf("omg_orb_p = %.8f\n", omg_orb0_p);
      printf("omg_orb_m = %.8f\n", omg_orb0_m);
    }  
  }
  eob_dyn_rootfind_rpr(dyn, &r0, &pr0abs, omg_orb0, r0_kepl, pr0PN);
  j0  = eob_dyn_j0(r0, dyn);      // Update j0

  if (DEBUG) printf("Initial j0: %.3f\nInitial r0: %.3f\n", j0, r0);

  double pr0 = pr0abs;
  if(fmod(zeta,2.*Pi) >= Pi)
    pr0 = -pr0abs;
      
  double rma, rcma, Ama, Bma, Qma, Hma;
  double pl_hold, omg_orb, dHeffma_dj0;
  
  rma = r0/(1.+ecc*cos(zeta));
    
  /* Computing metric, centrifugal radius and ggm functions*/
  if(usespins) {
    eob_metric_s(rma, pr0, dyn, &Ama, &Bma, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &Qma, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    eob_dyn_s_get_rc(rma, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rcma, &pl_hold, &pl_hold);
  } else {
    eob_metric(rma, pr0, dyn, &Ama, &Bma, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &Qma, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
  }

  /* Energy */
  eob_ham_s(nu, rma, rcma, 0, 0, j0, pr0, S, Sstar, chi1, chi2, X1, X2, aK2, c3, Ama, 0., 0., Qma, 0., 0., 0., 0.,
	    &Hma, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &dHeffma_dj0, &pl_hold, &pl_hold);

  /* Orbital frequency */
  omg_orb   = dHeffma_dj0/nu/Hma;

  y_init[EOB_ID_RAD]    = rma;
  y_init[EOB_ID_PHI]    = 0.;
  y_init[EOB_ID_PPHI]   = j0;
  y_init[EOB_ID_PRSTAR] = pr0;
  y_init[EOB_ID_PR]     = pr0*sqrt(Bma/Ama);
  y_init[EOB_ID_J]      = j0;
  y_init[EOB_ID_E0]     = Hma*nu;
  y_init[EOB_ID_OMGJ]   = omg_orb;
  
  return status;
}


/**
 * Function: eob_dyn_ic_ecc_ma_split
 * ---------------------------------
 *   Initial conditions calculation for eccentric systems
 *   Adiabatic level, computes pr*0 != 0 at a specified anomaly.
 *   Performs two one-dimensional root-finding steps to find the correct
 *   initial data. Determines j0 first, then pr*0.
 *
 *   @param[in] r0_kepl:  initial estimate of semilatus rectum
 *   @param[in] dyn: Dynamics
 *   @param[out] y_init: initial data
 *   @return status: 0 if successful
*/
int eob_dyn_ic_ecc_ma_split(double r0_kepl, Dynamics *dyn, double y_init[])
{
  int status = 0;
  const double f0 = EOBPars -> f0;
  const double nu = EOBPars -> nu;
  const double chi1 = EOBPars -> chi1;
  const double chi2 = EOBPars -> chi2;
  const double S1 = EOBPars -> S1;
  const double S2 = EOBPars -> S2;
  const double c3 = EOBPars -> cN3LO;
  const double X1 = EOBPars -> X1;
  const double X2 = EOBPars -> X2;
  const double a1 = EOBPars -> a1;
  const double a2 = EOBPars -> a2;
  const double aK2 = EOBPars -> aK2;
  const double C_Q1 = EOBPars -> C_Q1;
  const double C_Q2 = EOBPars -> C_Q2;
  const double C_Oct1 = EOBPars -> C_Oct1;
  const double C_Oct2 = EOBPars -> C_Oct2;
  const double C_Hex1 = EOBPars -> C_Hex1;
  const double C_Hex2 = EOBPars -> C_Hex2;
  const double S = S1 + S2;
  const double Sstar = X2*a1 + X1*a2;
  const double ecc = EOBPars -> ecc;
  const double zeta = EOBPars -> anomaly;

  const int usetidal = EOBPars -> use_tidal;
  const int usespins = EOBPars -> use_spins;

  double j0    = eob_dyn_j0(r0_kepl, dyn); // initial guess for j0
  double pr0PN = fabs(ecc/j0*sin(zeta));   // consider positive pr for bisection; correct sign later
  
  double r0, pr0abs;

  double omg_orb0 = Pi*f0;
  if (EOBPars->ecc_freq == ECCFREQ_ORBAVGD){
    /* Assume that the user gave as input the orbit averaged 
    frequency. Transform it into an average frequency between apastron
    and periastron and overwrite omg_orb0 */
    if (DEBUG) printf("Orbit-averaged ICs:\nomg_bar = %.8f\n", omg_orb0);
    double omg_orb0_p = eob_dyn_omg_from_omgbar(omg_orb0, 0., ecc);
    double omg_orb0_m = eob_dyn_omg_from_omgbar(omg_orb0, Pi, ecc);
    omg_orb0 = 0.5*(omg_orb0_p + omg_orb0_m);
    if (DEBUG) {
      printf("omg_orb0  = %.8f\n",  omg_orb0);
      printf("omg_orb_p = %.8f\n", omg_orb0_p);
      printf("omg_orb_m = %.8f\n", omg_orb0_m);
    }
  }
  r0 = eob_dyn_bisecOmegaecc0(dyn, omg_orb0, r0_kepl);
  if (r0 < 0.) return 1;    // Bisection failed

  j0 = eob_dyn_j0(r0, dyn); // Update j0
      
  double rma, rap, rper, rcap, Aap, Bap, Qap, Hap;
  double pl_hold, omg_orb, dHeffap_dj0;
  
  rap  = r0/(1.-ecc);
  rper = r0/(1.+ecc);
  if (DEBUG) printf("Periastron and apastron radii: %.3f %.3f\n", rper, rap);
  if (DEBUG) printf("Initial j0: %.3f\n", j0);
  rma  = r0/(1.+ecc*cos(zeta));
    
  /* Computing metric, centrifugal radius and ggm functions at apastron*/
  if(usespins) {
    eob_metric_s(rap, 0., dyn, &Aap, &Bap, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &Qap, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    eob_dyn_s_get_rc(rap, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rcap, &pl_hold, &pl_hold);
  } else {
    eob_metric(rap, 0., dyn, &Aap, &Bap, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &Qap, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
  }

  /* Energy */
  eob_ham_s(nu, rap, rcap, 0, 0, j0, 0., S, Sstar, chi1, chi2, X1, X2, aK2, c3, Aap, 0., 0., Qap, 0., 0., 0., 0.,
	    &Hap, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &dHeffap_dj0, &pl_hold, &pl_hold);

  /* Unless zeta = 0, pi, solve for pr0 */
  if (DEQUAL(zeta, 0., 1e-10) || DEQUAL(zeta, Pi, 1e-10))
      pr0abs = 0.;
  else 
      pr0abs = eob_dyn_bisecHam0(dyn, pr0PN, j0, Hap, rma);
      if (pr0abs < 0.) return 1; // Bisection failed

  double pr0 = pr0abs;
  if(fmod(zeta,2.*Pi) >= Pi)
    pr0 = -pr0abs;

  /* Recomputing metric, r at r(zeta) */
  double rcma, Ama, Bma, Qma, Hma;
  double dHeffma_dj0;

  if(usespins) {
    eob_metric_s(rma, pr0, dyn, &Ama, &Bma, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &Qma, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    eob_dyn_s_get_rc(rma, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rcma, &pl_hold, &pl_hold);
  } else {
    eob_metric(rma, pr0, dyn, &Ama, &Bma, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &Qma, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
  }

  eob_ham_s(nu, rma, rcma, 0, 0, j0, pr0, S, Sstar, chi1, chi2, X1, X2, aK2, c3, Ama, 0., 0., Qma, 0., 0., 0., 0.,
	    &Hma, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &dHeffma_dj0, &pl_hold, &pl_hold);

  /* Check that the pr0 found is truly a solution */
  double check = Hma/Hap;
  if (fabs(Hma/Hap - 1.) > 1e-5){
    if (DEBUG) printf("Hma/Hap - 1 = %.3e\n", fabs(Hma/Hap - 1.));
    return 1;
  }

  /* Orbital frequency */
  omg_orb   = dHeffma_dj0/nu/Hma;

  y_init[EOB_ID_RAD]    = rma;
  y_init[EOB_ID_PHI]    = 0.;
  y_init[EOB_ID_PPHI]   = j0;
  y_init[EOB_ID_PRSTAR] = pr0;
  y_init[EOB_ID_PR]     = pr0*sqrt(Bma/Ama);
  y_init[EOB_ID_J]      = j0;
  y_init[EOB_ID_E0]     = Hma*nu;
  y_init[EOB_ID_OMGJ]   = omg_orb;
  
  return status;

}

/**
 * Function: eob_dyn_ic_hyp
 * ------------------------
 *   Initial conditions calculation for hyperbolic systems
 * 
 *   @param[in] r0:  initial separation
 *   @param[in] dyn: Dynamics
 *   @param[out] y_init: initial data
 *   @return status: 0 if successful
*/
int eob_dyn_ic_hyp(double r0, Dynamics *dyn, double y_init[])
{
  int status = 0;
  const double H_ADM = EOBPars->H_hyp;
  const double j_ADM = EOBPars->j_hyp;  
  
  const double nu   = EOBPars->nu;
  const double X1   = EOBPars->X1;
  const double X2   = EOBPars->X2;
  
  const int usetidal = EOBPars->use_tidal;

  double A, B, Q, dA, dQ, dQ_dprstar, Heff, H0, E0, Omg0, dHeff0_dpphi, pl_hold;
  double pr0, prstar0, u, u2, a, b, c, Delta;

  u  = 1./r0;
  u2 = u*u;
  
  Heff = 1./(2.*nu)*(SQ(H_ADM) - 1.) + 1.;

  /* Computing metric */
  eob_metric(r0, 1., dyn, &A, &B, &dA, &pl_hold, &pl_hold, &pl_hold, &Q, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
  // since prstar = 1., the Q returned here should be equivalent to Q/(prstar4), which equals z3*u2 for Q3PN
  
  // a = z3*A*u2;
  a = A*Q;
  b = 1;
  c = A*(1. + SQ(j_ADM*u)) - SQ(Heff);

  Delta = SQ(b) - 4*a*c;

  if (Delta == 0.){
    prstar0 = - b/(2.*a);
  } else if (Delta > SQ(b)) {
    prstar0 = -sqrt((-b + sqrt(Delta))/(2.*a));
  } else {
    prstar0 = 0.;
  }
  pr0 = prstar0*sqrt(B/A);

  // compute Q & derivatives with new value of prstar0
  eob_metric(r0, prstar0, dyn, &A, &B, &dA, &pl_hold, &pl_hold, &pl_hold, &Q, &dQ, &dQ_dprstar, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);

  eob_ham(nu, r0, j_ADM, prstar0, A, dA, Q, dQ, dQ_dprstar, &H0, &pl_hold, NULL, NULL, &dHeff0_dpphi);
  E0   = nu*H0;
  Omg0 = dHeff0_dpphi/E0;
  
  y_init[EOB_ID_RAD]    = r0;
  y_init[EOB_ID_PHI]    = 0.;
  y_init[EOB_ID_PPHI]   = j_ADM;
  y_init[EOB_ID_PRSTAR] = prstar0;
  y_init[EOB_ID_PR]     = 0.;
  y_init[EOB_ID_J]      = 0.;
  y_init[EOB_ID_E0]     = E0;
  y_init[EOB_ID_OMGJ]   = Omg0;

  return status;
}

/**
 * Struct: DHeff0_tmp_params
 * -------------------------
 *  Auxiliary structure for the root finder
 * 
*/
struct DHeff0_tmp_params {
  double rorb, A, dA, rc, drc_dr, ak2, S, Ss, nu, chi1, chi2, X1, X2, c3;
};

/**
  * Function: eob_dyn_DHeff0
  * ------------------------
  *   Derivative of the effective Hamiltonian
  * 
  *   @param[in] x:  p_phi
  *   @param[in] params:  parameters
  *   @return dHeff_dr
  */
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

  double ggm0[26];
  eob_dyn_s_GS(rorb, rc, drc_dr, 0., ak2, 0., x, nu, chi1, chi2, X1, X2, c3, ggm0);
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

/**
  * Function: eob_dyn_bisecHeff0_s
  * ------------------------------
  *   Root finder: Compute minimum of Heff0
  * 
  *   @param[in] nu:  symmetric mass ratio
  *   @param[in] chi1:  dimensionless spin of body 1
  *   @param[in] chi2:  dimensionless spin of body 2
  *   @param[in] X1:  mass fraction of body 1
  *   @param[in] X2:  mass fraction of body 2
  *   @param[in] c3:  N3LO coefficient
  *   @param[in] pph:  orbital angular momentum
  *   @param[in] rorb:  orbital separation
  *   @param[in] A:  A potential
  *   @param[in] dA:  dA/dr
  *   @param[in] rc:  centrifugal radius
  *   @param[in] drc_dr:  drc/dr
  *   @param[in] ak2:  kerr spin parameter squared
  *   @param[in] S:  S1 + S2 total spin
  *   @param[in] Ss:  X1*a2 + X2*a1, Sstar
  *   @return rorb
*/
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

/**
 *  Function: eob_dyn_omg_from_omgbar()
 *  -----------------------------------
 *  Compute the instantaneous frequency from
 *  the orbit averaged frequency at a specific
 *  point in the orbit.
 *  Eq. (16) of arXiv 2309.15528, newtonian level
 * 
 *  @param[in] omgbar: orbit-averaged frequency
 *  @param[in] zeta  : true anomaly
 *  @param[in] ecc   : eccentricity
 * 
 *  @return omg: instantaneous frequency
 */
double eob_dyn_omg_from_omgbar(double omgbar, double zeta, double ecc)
{
  double fact1 = (ecc*cos(zeta)+1.);
  double den   = pow((1.-ecc*ecc), 3./2);
  return omgbar*fact1*fact1/den;
}

/**
  * Function: eob_dyn_r0_Kepler
  * ----------------------------
  *   Initial radius from initial frequency using Kepler's law
  * 
  *   @param[in] f0:  initial frequency
  *   @return r0
*/
double eob_dyn_r0_Kepler (double f0)
{
  const double omg_orb0 = Pi*f0; // =2*Pi*(f_orb0/2)
  return pow(omg_orb0, -2./3.);
}

/**
  * Function: eob_dyn_r0_circ
  * -------------------------
  *   Initial radius from initial frequency using EOB circular dynamics
  * 
  *   @param[in] f0:  initial frequency
  *   @param[in] dyn:  Dynamics
  *   @return r0
*/
double eob_dyn_r0_circ (double f0, Dynamics *dyn)
{
  const double omg_orb0 = Pi*f0;
  const double r0_kepl  = eob_dyn_r0_Kepler(f0);
  return eob_dyn_bisecOmegaorb0(dyn,omg_orb0,r0_kepl);
}

/**
  * Struct: Omegaorb0_tmp_params
  * ----------------------------
  *  Auxiliary structure for the root finder
  * 
*/
struct Omegaorb0_tmp_params {
  double omg_orb0;
  Dynamics *dyn;
};

/** Initial radius from initial frequency using EOB circular dynamics */

/**
 * Function: eob_dyn_r0_ecc
 * ------------------------
 *   Initial separation from initial frequency, eccentric case.
 *   Depending on the value of ecc_freq, the reference frequency
 *   is assumed to be ccorresponding to the apastron, periastron 
 *   or average separation
 * 
 *   @param[in] f0:  initial frequency
 *   @param[in] dyn:  Dynamics
 *  
 *   @return r0
*/
double eob_dyn_r0_ecc (double f0, Dynamics *dyn)
{
  const double omg_orb0 = Pi*f0;
  const double r0_kepl  = eob_dyn_r0_Kepler(f0);
  double ecc = EOBPars->ecc;

  double r0;

  /* Correct radius so to start at specified separation */
  if (EOBPars->ecc_freq == ECCFREQ_APASTRON) {
    /* apastron frequency */
    r0 = r0_kepl*(1-ecc);
  } else if (EOBPars->ecc_freq == ECCFREQ_PERIASTRON) {
    /* periastron frequency */
    r0 = r0_kepl*(1+ecc);
  } else {
    /* secular (average) frequency */
    /* use it for both ECCFREQ_AVERAGE and ECCFREQ_ORBAVGD*/
    r0 = r0_kepl;
  }

  return eob_dyn_bisecOmegaecc0(dyn,omg_orb0,r0);
}

/**
  * Struct: Omegaecc0_tmp_params
  * ----------------------------
  *  Auxiliary structure for the root finder, 
  *  eccentric case
  * 
*/
struct Omegaecc0_tmp_params {
  double omg_orb0;
  Dynamics *dyn;
};

/**
  * Function: eob_dyn_Omegaorb0
  * ----------------------------
  *   Root finder helper function: compute omg_orb - omg_orb0
  * 
  *   @param[in] r:  separation
  *   @param[in] params:  parameters
  *   @return omg_orb - omg_orb0
*/
double eob_dyn_Omegaorb0(double r, void *params)
{
 
  /* Unpack parameters */  
  struct Omegaorb0_tmp_params *p
    = (struct Omegaorb0_tmp_params *) params;
  double   omg_orb0 = p->omg_orb0;
  Dynamics *dyn     = p->dyn;

  const double nu    = EOBPars->nu;
  const double X1    = EOBPars->X1;
  const double X2    = EOBPars->X2;
  const double chi1  = EOBPars->chi1;
  const double chi2  = EOBPars->chi2;
  const double a1    = EOBPars->a1;
  const double a2    = EOBPars->a2;
  const double aK2   = EOBPars->aK2;
  const double S     = EOBPars->S;
  const double Sstar = EOBPars->Sstar;
  const double c3    = EOBPars->cN3LO;
  const double C_Q1  = EOBPars->C_Q1;
  const double C_Q2  = EOBPars->C_Q2;
  const double C_Oct1 = EOBPars->C_Oct1;
  const double C_Oct2 = EOBPars->C_Oct2;
  const double C_Hex1 = EOBPars->C_Hex1;
  const double C_Hex2 = EOBPars->C_Hex2;

  const int usetidal = EOBPars->use_tidal;  
  const int usespins = EOBPars->use_spins;

  double A,B,dA,rc,drc_dr,G,dG_dr,uc,uc2,dAuc2_dr,j02,j0,H,Heff,Heff_orb,dHeff_dj0,omg_orb;
  double pl_hold,a_coeff,b_coeff,c_coeff,Delta,sol_p,sol_m;
  double ggm[26];

  /* Computing metric, centrifugal radius and ggm functions*/
  if(usespins) {
    eob_metric_s(r, 0., dyn, &A, &B, &dA, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    eob_dyn_s_get_rc(r, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc, &drc_dr, &pl_hold);
    eob_dyn_s_GS(r, rc, drc_dr, 0., aK2, 0.0, 0.0, nu, chi1, chi2, X1, X2, c3, ggm);
    G     = ggm[2]*S + ggm[3]*Sstar;    // tildeG = GS*S+GSs*Ss
  } else {
    eob_metric(r, 0., dyn, &A, &B, &dA, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    rc     = r;   //Nonspinning case: rc = r; G = 0;
    G      = 0.;
  }

  /* Auxiliary variables*/
  uc       = 1./rc;
  uc2      = uc*uc;

  /* Circular angular momentum */
  j0  = eob_dyn_j0(r, dyn);
  j02 = SQ(j0);
    
  /* Circular Hamiltonians */
  Heff_orb = sqrt(A*(1+j02*uc2));
  Heff     = Heff_orb + j0*G;
  H        = sqrt(1+2*nu*(Heff-1))/nu;
    
  /* Circular orbital frequency */
  dHeff_dj0 = G + A*j0*uc2/Heff_orb;
  omg_orb   = dHeff_dj0/nu/H;

  /* Subtraction of initial evolution frequency */
  return (omg_orb - omg_orb0);
}

/**
  * Function: eob_dyn_bisecOmegaorb0
  * ---------------------------------
  *   Root finder: Compute r0 such that omg_orb = omg_orb0
  * 
  *   @param[in] dyn:  Dynamics
  *   @param[in] omg_orb0:  initial orbital frequency
  *   @param[in] r0_kepl:  initial radius from Kepler's law
  *   @return r0
*/
double eob_dyn_bisecOmegaorb0(Dynamics *dyn, double omg_orb0,double r0_kepl)
{
#define max_iter (200)
#define tolerance (1e-14)

  int status;
  int iter = 0;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;  
  double r0;
  double x_lo = 0.5*r0_kepl, x_hi = 1.5*r0_kepl;
  gsl_function F;
  
  struct  Omegaorb0_tmp_params p = {omg_orb0,dyn};
  
  F.function = &eob_dyn_Omegaorb0;
  F.params = &p;
  T = gsl_root_fsolver_bisection;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  do {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    r0     = gsl_root_fsolver_root (s);
    x_lo   = gsl_root_fsolver_x_lower (s);
    x_hi   = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo, x_hi, 0, tolerance);
  }
  while (status == GSL_CONTINUE && iter < max_iter);
  gsl_root_fsolver_free (s);
  
  return r0;
}

/**
  * Function: eob_dyn_Omegaecc0
  * ----------------------------
  *   Root finder helper function: compute omg_orb - omg_orb0
  *   Eccentric version
  * 
  *   @param[in] r:  semilatus rectum ("p" in celestial mechanics)
  *   @param[in] params:  parameters
  *   @return omg_orb - omg_orb0
*/
double eob_dyn_Omegaecc0(double r, void *params)
{
  /* Unpack parameters */  
  struct Omegaecc0_tmp_params *p
    = (struct Omegaecc0_tmp_params *) params;
  double   omg_orb0 = p->omg_orb0;
  Dynamics *dyn     = p->dyn;

  const double nu = EOBPars -> nu;
  const double X1 = EOBPars -> X1;
  const double X2 = EOBPars -> X2;
  const double chi1 = EOBPars -> chi1;
  const double chi2 = EOBPars -> chi2;
  const double a1 = EOBPars -> a1;
  const double a2 = EOBPars -> a2;
  const double aK2 = EOBPars -> aK2;
  const double S = EOBPars -> S;
  const double Sstar = EOBPars -> Sstar;
  const double c3 = EOBPars -> cN3LO;
  const double C_Q1 = EOBPars -> C_Q1;
  const double C_Q2 = EOBPars -> C_Q2;
  const double C_Oct1 = EOBPars -> C_Oct1;
  const double C_Oct2 = EOBPars -> C_Oct2;
  const double C_Hex1 = EOBPars -> C_Hex1;
  const double C_Hex2 = EOBPars -> C_Hex2;
  const double ecc = EOBPars -> ecc;

  const int usetidal = EOBPars -> use_tidal;
  const int usespins = EOBPars -> use_spins;

  double r1, A1, B1, rc1, G1, ggm1[26];
  double r2, A2, B2, rc2, G2, ggm2[26];
  
  double pl_hold, A12, B12, DA, DB, DG, j0, j02, omg_orb;
  double Heff_orb1, Heff_orb2, Heff1, Heff2, H1, H2, dHeff1_dj0, dHeff2_dj0, omg_orb1, omg_orb2;
  
  r1 = r/(1-ecc);
  r2 = r/(1+ecc);

  /* Computing metric, centrifugal radius and ggm functions*/
  if(usespins) {
    eob_metric_s(r1, 0., dyn, &A1, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    eob_dyn_s_get_rc(r1, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc1, &pl_hold, &pl_hold);
    
    eob_dyn_s_GS(r1, rc1, 0.0, 0.0, aK2, 0.0, 0.0, nu, chi1, chi2, X1, X2, c3, ggm1);
    G1     = ggm1[2]*S + ggm1[3]*Sstar;    // tildeG = GS*S+GSs*Ss

    eob_metric_s(r2, 0., dyn, &A2, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    eob_dyn_s_get_rc(r2, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc2, &pl_hold, &pl_hold);
    
    eob_dyn_s_GS(r2, rc2, 0.0, 0.0, aK2, 0.0, 0.0, nu, chi1, chi2, X1, X2, c3, ggm2);
    G2     = ggm2[2]*S + ggm2[3]*Sstar;    
  } else {
    eob_metric(r1, 0., dyn, &A1, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    rc1 = r1;   //Nonspinning case: rc = r; G = 0;  
    G1  = 0.0;

    eob_metric(r2, 0., dyn, &A2, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    rc2 = r2;
    G2  = 0.0;
  }

  /* Adiabatic angular momentum */
  j0  = eob_dyn_j0(r, dyn);
  j02 = SQ(j0);

  /* Energy */
  Heff_orb1 = sqrt(A1*(1. + j02/SQ(rc1)));
  Heff_orb2 = sqrt(A2*(1. + j02/SQ(rc2)));
  Heff1     = Heff_orb1 + j0*G1;
  Heff2     = Heff_orb2 + j0*G2;
  H1        = sqrt(1. + 2.*nu*(Heff1 - 1.))/nu;
  H2        = sqrt(1. + 2.*nu*(Heff2 - 1.))/nu;
  
  /* Orbital frequency */
  dHeff1_dj0 = G1 + j0*A1/(Heff_orb1*SQ(rc1));
  dHeff2_dj0 = G2 + j0*A2/(Heff_orb2*SQ(rc2));
  omg_orb1   = dHeff1_dj0/nu/H1;
  omg_orb2   = dHeff2_dj0/nu/H2;

  if (EOBPars->ecc_freq == ECCFREQ_APASTRON) {
    /* apastron frequency */
    omg_orb = omg_orb1;
  } else if (EOBPars->ecc_freq == ECCFREQ_PERIASTRON) {
    /* periastron frequency */
    omg_orb = omg_orb2;
  } else {
    /* secular (average) frequency */
    omg_orb = 0.5*(omg_orb1+omg_orb2);
  }
  
  /* Subtraction of initial evolution frequency */
  return (omg_orb - omg_orb0);
}

/**
  * Struct: Ham0_tmp_params
  * ----------------------------
  *  Auxiliary structure for the root finder, 
  *  eccentric case
  * 
*/
struct Ham0_tmp_params {
  double rma;
  double j0;
  double Hap;
  Dynamics *dyn;
};

/**
  * Function: eob_dyn_Ham0
  * ----------------------------
  *   Root finder helper function: compute H(pr) - H(apastron)
  *   Eccentric version
  * 
  *   @param[in] pr:  radial momentum
  *   @param[in] params:  parameters
  *   @return H(pr) - H(apastron)
*/
double eob_dyn_Ham0(double pr, void *params)
{
  /* Unpack parameters */  
  struct Ham0_tmp_params *p
    = (struct Ham0_tmp_params *) params;
  double   rma  = p->rma;
  double   j0   = p->j0;
  double   Hap  = p->Hap;
  Dynamics *dyn = p->dyn;

  const double nu = EOBPars -> nu;
  const double X1 = EOBPars -> X1;
  const double X2 = EOBPars -> X2;
  const double chi1 = EOBPars -> chi1;
  const double chi2 = EOBPars -> chi2;
  const double a1 = EOBPars -> a1;
  const double a2 = EOBPars -> a2;
  const double aK2 = EOBPars -> aK2;
  const double S = EOBPars -> S;
  const double Sstar = EOBPars -> Sstar;
  const double c3 = EOBPars -> cN3LO;
  const double C_Q1 = EOBPars -> C_Q1;
  const double C_Q2 = EOBPars -> C_Q2;
  const double C_Oct1 = EOBPars -> C_Oct1;
  const double C_Oct2 = EOBPars -> C_Oct2;
  const double C_Hex1 = EOBPars -> C_Hex1;
  const double C_Hex2 = EOBPars -> C_Hex2;

  const int usetidal = EOBPars -> use_tidal;
  const int usespins = EOBPars -> use_spins;
  
  double pl_hold, Ama, Bma, Qma, rcma, Hma;

  /* Computing metric, centrifugal radius and ggm functions*/
  if(usespins) {
    eob_metric_s(rma, pr, dyn, &Ama, &Bma, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &Qma, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    eob_dyn_s_get_rc(rma, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rcma, &pl_hold, &pl_hold);
  } else {
    eob_metric(rma, pr, dyn, &Ama, &Bma, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &Qma, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
  }

  /* Computing energy at rma = r(zeta) */
  eob_ham_s(nu, rma, rcma, 0, 0, j0, pr, S, Sstar, chi1, chi2, X1, X2, aK2, c3, Ama, 0., 0., Qma, 0., 0., 0., 0.,
	    &Hma, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
  
  /* Subtraction of apastron energy */
  return (Hma - Hap);
}

/**
  * Function: eob_dyn_bisecOmegaecc0
  * ---------------------------------
  *   Root finder: Compute eccentric p such that omg_orb = omg_orb0
  * 
  *   @param[in] dyn:  Dynamics
  *   @param[in] omg_orb0:  initial orbital frequency
  *   @param[in] r0_kepl:  initial radius from Kepler's law
  *   @return r0
*/
double eob_dyn_bisecOmegaecc0(Dynamics *dyn, double omg_orb0,double r0_kepl)
{
#define max_iter (200)
#define tolerance (1e-14)

  int status;
  int iter = 0;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;  
  gsl_function F;

  double r0;
  double x_lo = 0.5*r0_kepl, x_hi = 2.*r0_kepl;
  if (x_lo < 3.)  x_lo = 3.;
  if (x_hi < 10.) x_hi = 10.;
  struct  Omegaecc0_tmp_params p = {omg_orb0,dyn};
  
  /* Check that the bisection points straddle 0*/
  /* if they do not, slowly decrease x_lo */
  double f_xlo = eob_dyn_Omegaecc0(x_lo, &p);
  double f_xhi = eob_dyn_Omegaecc0(x_hi, &p);
  int  iter_r0 = 0;
  if (VERBOSE) PRSECTN("Bisection for eccentric initial conditions\n");
  if (VERBOSE) printf("initial f(x_lo) = %e, f(x_hi) = %e\n",f_xlo,f_xhi);
  while(f_xlo*f_xhi > 0. && iter_r0 < max_iter){
    x_lo  *= 0.99;
    f_xlo  = eob_dyn_Omegaecc0(x_lo, &p);
    iter_r0++;
    if (DEBUG) printf("\t iter %d: f(x_lo) = %e, x_lo = %e\n",iter_r0, f_xlo,x_lo);
  }
  if (VERBOSE) printf("f(x_lo) = %e, f(x_hi) = %e\n",f_xlo,f_xhi);

  // Did we find a bracketing interval? If not, fail and return a negative number
  if (f_xlo*f_xhi > 0. || isnan(f_xlo)){
    printf("ERROR: eob_dyn_bisecOmegaecc0, could not find a bracketing interval\n");
    return -1.;
  }

  F.function = &eob_dyn_Omegaecc0;
  F.params = &p;
  T = gsl_root_fsolver_bisection;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  do {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    r0     = gsl_root_fsolver_root (s);
    x_lo   = gsl_root_fsolver_x_lower (s);
    x_hi   = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo, x_hi, 0, tolerance);
  }
  while (status == GSL_CONTINUE && iter < max_iter);
  gsl_root_fsolver_free (s);
  
  return r0;
}

/**
  * Function: eob_dyn_bisecHam0
  * ---------------------------------
  *   Root finder: Compute pr0 such that energy at r(zeta) = energy at apastron
  * 
  *   @param[in] dyn:   Dynamics
  *   @param[in] pr0PN: initial guess for pr0
  *   @param[in] j0:    angular momentum
  *   @param[in] Hap:   energy at apastron
  *   @param[in] rma:   desired initial separation r(zeta)
  *   @return r0
*/
double eob_dyn_bisecHam0(Dynamics *dyn, double pr0PN, double j0, double Hap, double rma)
{
#define max_iter (200)
#define tolerance (1e-14)

  int status;
  int iter = 0;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;  
  gsl_function F;

  double pr0;
  double x_lo = 0., x_hi = 2.*pr0PN;
  struct  Ham0_tmp_params p = {rma, j0, Hap, dyn};
  
  /* Check that the bisection points straddle 0*/
  /* if they do not, slowly increase x_hi */
  double f_xlo = eob_dyn_Ham0(x_lo, &p);
  double f_xhi = eob_dyn_Ham0(x_hi, &p);
  int  iter_pr0 = 0;
  if (VERBOSE) PRSECTN("Bisection for pr in eccentric initial conditions\n");
  while(f_xlo*f_xhi > 0. && iter_pr0 < max_iter){
    x_hi  *= 1.01;
    f_xhi  = eob_dyn_Ham0(x_hi, &p);
    iter_pr0++;
    if (DEBUG) printf("\t iter %d: f(x_hi) = %e, x_hi = %e\n", iter_pr0, f_xhi, x_hi);
  }
  if (VERBOSE) printf("f(x_lo) = %e, f(x_hi) = %e\n",f_xlo,f_xhi);
  if (f_xlo*f_xhi > 0.){
    // print an error, the bisection will fail. Return a negative number
    printf("ERROR: eob_dyn_bisecHam0, could not find a bracketing interval\n");
    return -1.;
  }
  F.function = &eob_dyn_Ham0;
  F.params = &p;
  T = gsl_root_fsolver_bisection;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  do {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    pr0    = gsl_root_fsolver_root (s);
    x_lo   = gsl_root_fsolver_x_lower (s);
    x_hi   = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo, x_hi, 0, tolerance);
  }
  while (status == GSL_CONTINUE && iter < max_iter);
  gsl_root_fsolver_free (s);
  
  return pr0;
}

/** 
 * Function: eob_dyn_rootfind_rpr
 * ------------------------------
 * 
 *  Root finder for initial conitions (r, pr) in eccentric case
 *  
 *  @param[in]    dyn:  Dynamics
 *  @param[in,out] r0:  semilatus rectum
 *  @param[in,out] pr0:  prstar
 *  @param[in]    omg_orb0:  initial orbital frequency
 *  @param[in]    rKepl:  initial estimate of the semilatus rectum
 *  @param[in]    pr0PN:  initial estimate of prstar
 * 
 *  @return status
 */
int eob_dyn_rootfind_rpr(Dynamics *dyn, double *r0, double *pr0, double omg_orb0, double rKepl, double pr0PN)
{
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  int status;
  size_t i, iter = 0;
  
  const size_t n = 2;

  struct  Omegaorb0_tmp_params params = {omg_orb0,dyn};
  
  gsl_multiroot_function f = {&eob_dyn_rpr, n, &params};

  double x_init[2] = {rKepl, pr0PN};
  gsl_vector *x = gsl_vector_alloc(n);

  gsl_vector_set (x, 0, x_init[0]);
  gsl_vector_set (x, 1, x_init[1]);

  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, 2);
  gsl_multiroot_fsolver_set (s, &f, x);

  const int max_it    = 200;
  const double epsrel = 1e-10;
  
  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate(s);
      
      if (status) break;

      status = gsl_multiroot_test_residual (s->f, epsrel);
    }
  while (status == GSL_CONTINUE && iter < max_it);

  double x0 = gsl_vector_get(s->x,0);
  double x1 = gsl_vector_get(s->x,1);
  
  *r0 = 0.;
  if (isfinite(x0)) *r0 = fabs(x0);
  
  *pr0 = 0.;
  if (isfinite(x1)) *pr0 = fabs(x1);
  
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);
  
  //if (status == ???) {
  //  return ROOT_ERRORS_BRACKET;
  //}
  if (status == GSL_SUCCESS) {
    return ROOT_ERRORS_NO;
  } 
  if (iter >= max_it) {
    return ROOT_ERRORS_MAXITS;
  }
  if (status != GSL_SUCCESS) {
    return ROOT_ERRORS_NOSUCC;
  }
  
  return status;
}

/** 
 * Function: eob_dyn_rpr
 * ---------------------
 * 
 *  Root function to compute r0 and pr0, uses GSL multiroot
 *  solves the system:
 *    - E(apastron, pr=0) - E(r=r0/(1+ecc*cos(anomaly)), pr=pr0) = 0;
 *    - Omega(r=r0/(1+ecc*cos(anomaly)), pr=pr0) - Omg0 = 0;
 * 
 *   @param[in] x     :  (r0, pr0)
 *   @param[in] params:  parameters to be passed to the function
 *   @param[in,out] f :  system of equations  
 *   @return status
 */
int eob_dyn_rpr(const gsl_vector *x, void * params, gsl_vector *f)
{

  /* Unpack parameters */  
  struct Omegaecc0_tmp_params *p
    = (struct Omegaecc0_tmp_params *) params;
  double   omg_orb0 = p->omg_orb0;
  Dynamics *dyn     = p->dyn;
  
  double r0  = fabs(gsl_vector_get(x,0)); //semilatus rectum
  double pr0 = fabs(gsl_vector_get(x,1));

  const double nu     = EOBPars -> nu;
  const double X1     = EOBPars -> X1;
  const double X2     = EOBPars -> X2;
  const double chi1   = EOBPars -> chi1;
  const double chi2   = EOBPars -> chi2;
  const double a1     = EOBPars -> a1;
  const double a2     = EOBPars -> a2;
  const double aK2    = EOBPars -> aK2;
  const double S      = EOBPars -> S;
  const double Sstar  = EOBPars -> Sstar;
  const double c3     = EOBPars -> cN3LO;
  const double C_Q1   = EOBPars -> C_Q1;
  const double C_Q2   = EOBPars -> C_Q2;
  const double C_Oct1 = EOBPars -> C_Oct1;
  const double C_Oct2 = EOBPars -> C_Oct2;
  const double C_Hex1 = EOBPars -> C_Hex1;
  const double C_Hex2 = EOBPars -> C_Hex2;
  const double ecc    = EOBPars -> ecc;
  const double zeta   = EOBPars -> anomaly;

  const int usetidal  = EOBPars -> use_tidal;
  const int usespins  = EOBPars -> use_spins;

  double r1, A1, B1, rc1, G1, ggm1[26];
  double r2, A2, B2, rc2, G2, ggm2[26];
  double rma, rcma, Ama, Bma, Qma;
  
  double pl_hold, A12, B12, DA, DB, DG, j0, j02, omg_orb;
  
  double Heff_orb1, Heff_orb2, Heff1, Heff2, H1, H2, dHeff1_dj0, dHeff2_dj0, omg_orb1, omg_orb2;
  double Hma, dHeffma_dj0;
  
  r1  = r0/(1.-ecc);
  r2  = r0/(1.+ecc);
  rma = r0/(1.+ecc*cos(zeta));

  /* Computing metric, centrifugal radius and ggm functions*/
  if(usespins) {
    eob_metric_s(r1, 0., dyn, &A1, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    eob_dyn_s_get_rc(r1, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc1, &pl_hold, &pl_hold);
    
    eob_dyn_s_GS(r1, rc1, 0.0, 0.0, aK2, 0.0, 0.0, nu, chi1, chi2, X1, X2, c3, ggm1);
    G1     = ggm1[2]*S + ggm1[3]*Sstar;    // tildeG = GS*S+GSs*Ss

    eob_metric_s(r2, 0., dyn, &A2, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    eob_dyn_s_get_rc(r2, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc2, &pl_hold, &pl_hold);
    
    eob_dyn_s_GS(r2, rc2, 0.0, 0.0, aK2, 0.0, 0.0, nu, chi1, chi2, X1, X2, c3, ggm2);
    G2     = ggm2[2]*S + ggm2[3]*Sstar;

    eob_metric_s(rma, pr0, dyn, &Ama, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &Qma, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    eob_dyn_s_get_rc(rma, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rcma, &pl_hold, &pl_hold);
  } else {
    eob_metric(r1, 0., dyn, &A1, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    rc1 = r1;   //Nonspinning case: rc = r; G = 0;  
    G1  = 0.0;

    eob_metric(r2, 0., dyn, &A2, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
    rc2 = r2;
    G2  = 0.0;
    eob_metric(rma, pr0, dyn, &Ama, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &Qma, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
  }

  /* Adiabatic angular momentum */
  j0  = eob_dyn_j0(r0, dyn);
  j02 = SQ(j0);
  

  /* Energy */
  Heff_orb1 = sqrt(A1*(1. + j02/SQ(rc1)));
  Heff_orb2 = sqrt(A2*(1. + j02/SQ(rc2)));
  Heff1     = Heff_orb1 + j0*G1;
  Heff2     = Heff_orb2 + j0*G2;
  H1        = sqrt(1. + 2.*nu*(Heff1 - 1.))/nu;
  H2        = sqrt(1. + 2.*nu*(Heff2 - 1.))/nu;
  
  /* Energy at point specified by anomaly */
  eob_ham_s(nu, rma, rcma, 0, 0, j0, pr0, S, Sstar, chi1, chi2, X1, X2, aK2, c3, Ama, 0., 0., Qma, 0., 0., 0., 0.,
	    &Hma, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &dHeffma_dj0, &pl_hold, &pl_hold);
  
  /* Orbital frequency */
  dHeff1_dj0 = G1 + j0*A1/(Heff_orb1*SQ(rc1));
  dHeff2_dj0 = G2 + j0*A2/(Heff_orb2*SQ(rc2));
  omg_orb1   = dHeff1_dj0/nu/H1;
  omg_orb2   = dHeff2_dj0/nu/H2;

  if (EOBPars->ecc_freq == ECCFREQ_APASTRON) {
    /* apastron frequency */
    omg_orb = omg_orb1;
  } else if (EOBPars->ecc_freq == ECCFREQ_PERIASTRON) {
    /* periastron frequency */
    omg_orb = omg_orb2;
  } else {
    /* secular (average) frequency */
    omg_orb = 0.5*(omg_orb1+omg_orb2);
  }

  /* Orbital frequency at point specified by anomaly */
  //omg_orb   = dHeffma_dj0/nu/Hma;

  double eq1 = Hma - H1;
  double eq2 = omg_orb - omg_orb0;
    
  gsl_vector_set (f, 0, eq1);
  gsl_vector_set (f, 1, eq2);
  
  return GSL_SUCCESS;

}
