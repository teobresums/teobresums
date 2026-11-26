/** \file TEOBResumSPostAdiab.c
 *  \brief Post-adiabatic dynamics
 * 
 *  This file contains the function to evolve the post-adiabatic dynamics
 *  of the EOB system in both the aligned- and precessing- spins case.
 *
 */
#include <string.h>
#include "TEOBResumS.h"

#define nv (25) /* temp arrays */

/** Post-adiabatic dynamics */

/**
 * Function: eob_dyn_Npostadiabatic
 * -------------------------------
 *   Computes the post-adiabatic dynamics of the EOB system.
 *   If the system is spin precessing and EOBPars->project_spins is set to 1,
 *   the spins are projected on the orbital angular momentum at each iteration.
 * 
 *   @param[in]  dyn  : EOB Dynamics structure
 *   @param[in]  r0   : Initial separation
 *   @param[in]  spin : DynamicsSpins structure
*/
int eob_dyn_Npostadiabatic(Dynamics *dyn, const double r0, DynamicsSpin *spin)
{
  /* Unpack values */
  const double nu    = EOBPars->nu;
  double S           = EOBPars->S;
  double Sstar       = EOBPars->Sstar;
  double chi1        = EOBPars->chi1;
  double chi2        = EOBPars->chi2;
  const double X1    = EOBPars->X1;
  const double X2    = EOBPars->X2;
  const double c3    = EOBPars->cN3LO;
  double aK2   = EOBPars->aK2;
  double a1    = EOBPars->a1;
  double a2    = EOBPars->a2;
  const double C_Q1  = EOBPars->C_Q1;
  const double C_Q2  = EOBPars->C_Q2;
  const double C_Oct1 = EOBPars->C_Oct1;
  const double C_Oct2 = EOBPars->C_Oct2;
  const double C_Hex1 = EOBPars->C_Hex1;
  const double C_Hex2 = EOBPars->C_Hex2;
  const int usetidal = EOBPars->use_tidal;
  const int usespins = EOBPars->use_spins;

 /* Parameters for post adiabatic dynamics */
  const int Npa = EOBPars->postadiabatic_dynamics_N;    
  const int size = dyn->size; // PA clause in main.c sets dyn->size= PA size overwriting PA par files
  double rmin = EOBPars->postadiabatic_dynamics_rmin; 
  //if(usetidal) rmin = EOBPars->postadiabatic_dynamics_rmin;   
  const double dr = EOBPars->postadiabatic_dynamics_dr; //(r0 - rmin)/(size-1); /* Uniform grid spacing */

  /* If generic spins, build a grid of chi1z, chi2z by interpolating the PN dynamics in Omega_circ*/
  double *chi1_grid, *chi2_grid;
  double SApara, SBpara;
  if (usespins == MODE_SPINS_GENERIC && EOBPars->project_spins){
    chi1_grid = (double*) malloc(size * sizeof(double));
    chi2_grid = (double*) malloc(size * sizeof(double));
    for(int i=0; i<size;i++){
      const double Omg_i = pow(r0 - i*dr, -1.5);
      eob_spin_dyn_Sproj_interp(spin, Omg_i, &SApara, &SBpara, NULL, NULL, NULL, NULL, 1);
      chi1_grid[i] = SApara/SQ(X1);
      chi2_grid[i] = SBpara/SQ(X2);
    }
  }
  if (VERBOSE) {
    PRFORMd("post_adiabatic_dynamics_r0",r0);
    PRFORMd("post_adiabatic_dynamics_rmin",rmin);
    PRFORMi("post_adiabatic_dynamics_size",size);
    PRFORMd("post_adiabatic_dynamics_dr",dr);
  }

  /* Additional memory */
  double *buffer[nv]; 
  for (int v=0; v < nv; v++)
    buffer[v] = (double*)malloc(size * sizeof (double));

  double *A_vec                  = buffer[0];
  double *dA_vec                 = buffer[1];
  double *B_vec                  = buffer[2];
  double *sqrtAbyB_vec           = buffer[3];
  double *rc_vec                 = buffer[4];
  double *drc_dr_vec             = buffer[5];
  double *uc2_vec                = buffer[6];
  double *duc_dr_vec             = buffer[7];
  double *dAuc2_dr_vec           = buffer[8];
  double *dG_dr_vec              = buffer[9];
  double *dG_dprstar_vec         = buffer[10];
  double *dG_dprstarbyprstar_vec = buffer[11];
  double *G0_vec                 = buffer[12];
  double *dG_dr0_vec             = buffer[13];
  double *E_vec                  = buffer[14];
  double *Heff_vec               = buffer[15];
  double *Heff_orb_vec           = buffer[16];
  double *dpphi_dr_vec           = buffer[17];
  double *dprstar_dr_vec         = buffer[18];
  double *dphi_dr_vec            = buffer[19];
  double *dt_dr_vec              = buffer[20];
  double *Q_vec                  = buffer[21];
  double *dQ_vec                 = buffer[22];
  double *dQdprstar_vec          = buffer[23];
  double *d2Qdprstar2_vec        = buffer[24];
  
  double ggm[26]; 
  double a_coeff, b_coeff, c_coeff, Delta, sol_p, sol_m, j02, uc, u2, prstar2, dHeff_dpphi, dHeff_dprstar, dHeff_dr, dHeff_dprstarbyprstar, d2Heff_dprstar20,
    H, G, pl_hold, x, jhat, psi, r_omg, v_phi, Fphi, dr_dtbyprstar, prstar4, Heff_orb_f, Heff_f, E_f;
  /* 
   * Compute circular dynamics 
   */
  for (int i = 0; i < size; i++) {
    
    /* Current radius */
    dyn->r = r0 - i*dr;
    if (usespins == MODE_SPINS_GENERIC && EOBPars->project_spins){
      /* Current spins */
      chi1 = chi1_grid[i];
      chi2 = chi2_grid[i];
      set_spin_vars(X1,X2,chi1,chi2, &pl_hold,&pl_hold,&a1,&a2,&pl_hold,&aK2,&S,&Sstar);
      EOBPars->chi1 = chi1_grid[i];
      EOBPars->chi2 = chi2_grid[i];
      set_spin_vars(X1,X2,EOBPars->chi1,EOBPars->chi2, 
          &EOBPars->S1, &EOBPars->S2,
          &EOBPars->a1, &EOBPars->a2,
          &EOBPars->aK, &EOBPars->aK2,
          &EOBPars->S, &EOBPars->Sstar);  
    }

    /** Computing metric functions and centrifugal radius */
    if(usespins){ 
      
      eob_metric_s(dyn->r, 0., dyn, &A_vec[i], &B_vec[i], &dA_vec[i], &pl_hold, &pl_hold, &pl_hold, &Q_vec[i], &dQ_vec[i], &dQdprstar_vec[i], &pl_hold, &pl_hold, &d2Qdprstar2_vec[i], &pl_hold, &pl_hold, &pl_hold);
      eob_dyn_s_get_rc(dyn->r, nu, EOBPars->a1, EOBPars->a2, EOBPars->aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc_vec[i], &drc_dr_vec[i], &pl_hold);
      eob_dyn_s_GS(dyn->r, rc_vec[i], drc_dr_vec[i], pl_hold, EOBPars->aK2, 0.0, 0.0, nu, chi1, chi2, X1, X2, c3, ggm);
      
      G                         = ggm[2] *EOBPars->S+ggm[3] *EOBPars->Sstar;    // tildeG = GS*S+GSs*Ss
      dG_dr_vec[i]              = ggm[6] *EOBPars->S+ggm[7] *EOBPars->Sstar;
      dG_dprstar_vec[i]         = ggm[4] *EOBPars->S+ggm[5] *EOBPars->Sstar;
      dG_dprstarbyprstar_vec[i] = ggm[10]*EOBPars->S+ggm[11]*EOBPars->Sstar;
    
    } else {
      
      eob_metric(dyn->r, 0., dyn, &A_vec[i], &B_vec[i], &dA_vec[i], &pl_hold, &pl_hold, &pl_hold, &Q_vec[i], &dQ_vec[i], &dQdprstar_vec[i], &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
      
      rc_vec[i]                 = dyn->r; //Nonspinning case: rc = r
      drc_dr_vec[i]             = 1;  

      G                         = 0.0;
      dG_dr_vec[i]              = 0.0;
      dG_dprstar_vec[i]         = 0.0;
      dG_dprstarbyprstar_vec[i] = 0.0;
    
    }

    /** Defining circular quantities for the flux calculation.
	  Must not be overwritten in successive iterations, thus
	  we define separate quantities with the subscripts 0. */
    G0_vec[i]        = G;
    dG_dr0_vec[i]    = dG_dr_vec[i];
      
    /* Auxiliary variables*/
    sqrtAbyB_vec[i] = sqrt(A_vec[i]/B_vec[i]);
    uc              = 1./rc_vec[i];
    uc2_vec[i]      = uc*uc;
    duc_dr_vec[i]   = -uc2_vec[i]*drc_dr_vec[i];
    dAuc2_dr_vec[i] = uc2_vec[i]*(dA_vec[i]-2*A_vec[i]*uc*drc_dr_vec[i]);
    
    /** Computing the circular angular momentum by solving eq. (A15) of TEOBResumS paper 
	  (which is equivalent to solve eq.(4)=0 of arXiv:1805.03891). 
    */
                                                                                   
    if (usespins) {
      
      a_coeff = SQ(dAuc2_dr_vec[i]) - 4*A_vec[i]*uc2_vec[i]*SQ(dG_dr_vec[i]);  /* First coefficient of the quadratic equation a*x^2+b*x+c=0 */
      b_coeff = 2*dA_vec[i]*dAuc2_dr_vec[i] - 4*A_vec[i]*SQ(dG_dr_vec[i]);     /* Second coefficient of the quadratic equation */
      c_coeff = SQ(dA_vec[i]);                                                 /* Third coefficient of the quadratic equation */

      /* Delta of the quadratic equation */
      Delta = SQ(b_coeff) - 4*a_coeff*c_coeff; 
      if (Delta<0)  
	   /* If the spins are very small, 
	   numerical fluctuations sometimes make Delta negative (e.g. -1e-30). 
	   Setting it to 0 by hand */
          Delta=0.;                                              
      
      sol_p   = (-b_coeff + sqrt(Delta))/(2*a_coeff); /* Plus  solution of the quadratic equation */
      sol_m   = (-b_coeff - sqrt(Delta))/(2*a_coeff); /* Minus solution of the quadratic equation */
      
      /* dGdr sign determines choice of solution */
      if (dG_dr0_vec[i] > 0) j02 = sol_p;
      else                   j02 = sol_m;
    
    } else {
      
      a_coeff = dAuc2_dr_vec[i];
      b_coeff = dA_vec[i];
      
      j02 = -b_coeff/a_coeff;
    }
  
    /** Define momenta in the circular orbit approximation */
    dyn->pphi                = sqrt(j02);
    dyn->prstar              = 0.0;
    dprstar_dr_vec[i]        = 0.0;
    
    /** Circular Hamiltonians, ref: arXiv: 1406.6913 */
    if(usespins) {
      
      eob_ham_s(nu, dyn->r, rc_vec[i], drc_dr_vec[i], 0., dyn->pphi, dyn->prstar, EOBPars->S, EOBPars->Sstar, chi1, chi2, X1, X2, EOBPars->aK2, c3, A_vec[i], dA_vec[i], 0., Q_vec[i], dQ_vec[i], dQdprstar_vec[i], 0., d2Qdprstar2_vec[i],
                &H,               /* real EOB Hamiltonian divided by mu=m1m2/(m1+m2) */
                &Heff_vec[i],     /* effective EOB Hamiltonian (divided by mu)       */
                &Heff_orb_vec[i],
                NULL,             /* drvt Heff,r      */
                NULL,             /* drvt Heff,prstar */
                &dHeff_dpphi,     /* drvt Heff,pphi   */
                &d2Heff_dprstar20,
		NULL);
      
      E_vec[i] = nu*H;
      
    } else {
      
      eob_ham(nu, dyn->r, dyn->pphi, dyn->prstar, A_vec[i], dA_vec[i], Q_vec[i], dQ_vec[i], dQdprstar_vec[i],
              &H,               /* real EOB Hamiltonian divided by mu=m1m2/(m1+m2) */
              &Heff_orb_vec[i], /* effective EOB Hamiltonian (divided by mu). */
              NULL,             /* drvt Heff,r      */
              NULL,             /* drvt Heff,prstar */
              &dHeff_dpphi);    /* drvt Heff,pphi   */
      
      d2Heff_dprstar20 = 1/Heff_orb_vec[i];
      
      Heff_vec[i] = Heff_orb_vec[i]; /* Heff coincides with Heff_orb for the non-spinning case */
      E_vec[i] = nu*H;
      
    }
    
    /* Circular orbital frequency */
    dyn->Omg     = dHeff_dpphi/E_vec[i];
    
    /* Circular real orbital frequency */
    dyn->Omg_orb = (dyn->pphi*A_vec[i]*uc2_vec[i])/(E_vec[i]*Heff_orb_vec[i]);
    
    /** ddotr */
    dyn->ddotr   = 0.;
    
    dyn->data[EOB_RAD][i]    = dyn->r;
    dyn->data[EOB_PPHI][i]   = dyn->pphi;
    dyn->data[EOB_PRSTAR][i] = dyn->prstar;
    dyn->data[EOB_DDOTR][i]  = dyn->ddotr;
    dyn->data[EOB_MOMG][i]   = dyn->Omg;
    dyn->data[EOB_OMGORB][i] = dyn->Omg_orb;
  } // END r-GRID FOR

  /* update the spins grid in Omg if precession & projection*/
  if (usespins == MODE_SPINS_GENERIC && EOBPars->project_spins){
    for(int i=0; i<size;i++){
      eob_spin_dyn_Sproj_interp(spin, dyn->data[EOB_MOMG][i],&SApara, &SBpara, NULL, NULL, NULL, NULL, 1);
      chi1_grid[i] = SApara/SQ(X1);;
      chi2_grid[i] = SBpara/SQ(X2);;
    }
  }

  /** Computing angular momentum derivative */
  D0(dyn->data[EOB_PPHI],-dr, size, dpphi_dr_vec); /* dJ0/dr */

  /*
   * Post-Adiabatic dynamics 
  */
  
  int parity = 1; /* parity of the post-adiab iteration */

  /* For on PA orders */
  for (int n = 1; n <= Npa; n++) {
    
    /* Separating even and odd orders */
    if (n%2==0) parity = 0;
    else        parity = 1;
    
    /* For on r-grid */
    for (int i = 0; i < size; i++) {
      
      /* Setting loop variables to help reader */
      dyn->r       = dyn->data[EOB_RAD][i];
      dyn->phi     = dyn->data[EOB_PHI][i];
      dyn->pphi    = dyn->data[EOB_PPHI][i];
      dyn->Omg     = dyn->data[EOB_MOMG][i];
      dyn->ddotr   = dyn->data[EOB_DDOTR][i];
      dyn->prstar  = dyn->data[EOB_PRSTAR][i];
      dyn->Omg_orb = dyn->data[EOB_OMGORB][i];
      if (usespins == MODE_SPINS_GENERIC && EOBPars->project_spins){

        chi1         = chi1_grid[i];
        chi2         = chi2_grid[i];
        set_spin_vars(X1,X2,chi1,chi2,&pl_hold,&pl_hold,&a1,&a2,&pl_hold,&aK2,&S,&Sstar);
        EOBPars->chi1 = chi1_grid[i];
        EOBPars->chi2 = chi2_grid[i];
        set_spin_vars(X1,X2,EOBPars->chi1,EOBPars->chi2, 
          &EOBPars->S1, &EOBPars->S2,
          &EOBPars->a1, &EOBPars->a2,
          &EOBPars->aK, &EOBPars->aK2,
          &EOBPars->S, &EOBPars->Sstar);

        /* recompute rc and related quantities*/
        eob_dyn_s_get_rc(dyn->r, nu, EOBPars->a1, EOBPars->a2, EOBPars->aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc_vec[i], &drc_dr_vec[i], &pl_hold);
        eob_dyn_s_GS(dyn->r, rc_vec[i], drc_dr_vec[i], pl_hold, EOBPars->aK2, dyn->prstar, dyn->pphi, nu, chi1, chi2, X1, X2, c3, ggm);
      
        G                         = ggm[2] *EOBPars->S+ggm[3] *EOBPars->Sstar;    // tildeG = GS*S+GSs*Ss
        dG_dr_vec[i]              = ggm[6] *EOBPars->S+ggm[7] *EOBPars->Sstar;
        dG_dprstar_vec[i]         = ggm[4] *EOBPars->S+ggm[5] *EOBPars->Sstar;
        dG_dprstarbyprstar_vec[i] = ggm[10]*EOBPars->S+ggm[11]*EOBPars->Sstar;
        G0_vec[i]        = G;
        dG_dr0_vec[i]    = dG_dr_vec[i];
          
        /* Auxiliary variables*/
        sqrtAbyB_vec[i] = sqrt(A_vec[i]/B_vec[i]);
        uc              = 1./rc_vec[i];
        uc2_vec[i]      = uc*uc;
        duc_dr_vec[i]   = -uc2_vec[i]*drc_dr_vec[i];
        dAuc2_dr_vec[i] = uc2_vec[i]*(dA_vec[i]-2*A_vec[i]*uc*drc_dr_vec[i]);
      }

      if (parity)  {
	
      /* ***********************************
      * Odd PA orders : prstar corrections 
      * ********************************** */
      
      /** Calculating the flux Fphi */
      //FIXME USE C-routines, jhat etc. are already present inside dynamics
      //FIXME Non-spinning routine gives 1e-2 difference between PA and full EOB waveform. Tested cases: bbh q 1 f 0.001 and q 5 f 0.006.

        if (usespins) {

          /* Variables for which Kepler's law is still valid */
          Heff_orb_f = sqrt(A_vec[i]*(1.0 + SQ(dyn->pphi)*uc2_vec[i]));
          Heff_f     = G0_vec[i]*dyn->pphi + Heff_orb_f;
          E_f        = sqrt(1 + 2*nu*(Heff_f - 1));
          psi        = fabs((duc_dr_vec[i] + dG_dr0_vec[i]*rc_vec[i]*sqrt(A_vec[i]/(SQ(dyn->pphi)) + A_vec[i]*uc2_vec[i])/A_vec[i])/(-0.5*dA_vec[i]));
          // FIXME: Different from Matlab code.
          //        Added absolute value to avoid NaN
          r_omg      = 1.0/cbrt(SQ(((1./sqrt(rc_vec[i]*rc_vec[i]*rc_vec[i]*psi))+G0_vec[i])/(E_f)));
          v_phi      = r_omg*dyn->Omg;
          x          = SQ(v_phi);
          jhat       = dyn->pphi/(r_omg*v_phi);

          Fphi = eob_flx_Flux_s(x, dyn->Omg, r_omg, E_vec[i], Heff_vec[i], jhat, dyn->r, dyn->prstar, dyn->ddotr, dyn);
        
        } else {
        
          Heff_orb_f = sqrt(A_vec[i]*(1.0 + SQ(dyn->pphi)*uc2_vec[i]));
          Heff_f     = Heff_orb_f;
          psi        = 2.*(1.0 + 2.0*nu*(Heff_orb_f - 1.0))/(SQ(dyn->r)*dA_vec[i]);
          r_omg      = dyn->r*cbrt(psi);
          v_phi      = r_omg*dyn->Omg;
          x          = SQ(v_phi);
          jhat       = dyn->pphi/(r_omg*v_phi);
          
          Fphi = eob_flx_Flux(x, dyn->Omg, r_omg, E_vec[i], Heff_vec[i], jhat, dyn->r, dyn->prstar, dyn->ddotr, dyn);
        }

        dyn->data[EOB_FPHI][i] = Fphi; /* Fr = 0 in PA */

        /* recomputing Q */
        if(usespins) {
          eob_metric_s(dyn->r, dyn->prstar, dyn, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &Q_vec[i], &dQ_vec[i], &dQdprstar_vec[i], &pl_hold, &pl_hold, &d2Qdprstar2_vec[i], &pl_hold, &pl_hold, &pl_hold);
        } else {
           eob_metric(dyn->r, dyn->prstar, dyn, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &Q_vec[i], &dQ_vec[i], &dQdprstar_vec[i], &pl_hold, &pl_hold, &d2Qdprstar2_vec[i], &pl_hold, &pl_hold, &pl_hold); 
        }

        /** Calculating prstar */
        if (n < 2) {
          dHeff_dprstarbyprstar = dyn->pphi*dG_dprstarbyprstar_vec[i] + 1./Heff_orb_vec[i];	
        } else {
          dHeff_dprstarbyprstar = dyn->pphi*dG_dprstarbyprstar_vec[i] + (1. + 0.5*A_vec[i]*dQdprstar_vec[i]/(dyn->prstar))/Heff_orb_vec[i];	
        }
        dr_dtbyprstar         = sqrtAbyB_vec[i]/(E_vec[i])*dHeff_dprstarbyprstar;
        dyn->prstar           = Fphi/dpphi_dr_vec[i]/dr_dtbyprstar; 
    
        /** Note: p_phi does not change at odd orders 
            Computing first PA using the approximation detailed above A19 of TEOBResumS paper and Hamilton's equations.   
        */
 
        /** New GGM functions */
        eob_dyn_s_GS(dyn->r, rc_vec[i], drc_dr_vec[i], pl_hold, EOBPars->aK2, dyn->prstar, 0.0, nu, chi1, chi2, X1, X2, c3, ggm);
        
        dG_dr_vec[i]              = ggm[6] *EOBPars->S+ggm[7] *EOBPars->Sstar;
        dG_dprstar_vec[i]         = ggm[4] *EOBPars->S+ggm[5] *EOBPars->Sstar;
        dG_dprstarbyprstar_vec[i] = ggm[10]*EOBPars->S+ggm[11]*EOBPars->Sstar;
   
      } else {
	
        /* ***********************************
        * Even PA orders : pphi corrections 
        * ********************************** */
    
        prstar4 = SQ(SQ(dyn->prstar));
        a_coeff = dAuc2_dr_vec[i];                   /* coefficients of the quadratic equation a*x^2+b*x+c=0 */
        b_coeff = 2*Heff_orb_vec[i]*(dG_dr_vec[i] + dG_dprstar_vec[i]*dprstar_dr_vec[i]);
        c_coeff = dA_vec[i] + 2.*dprstar_dr_vec[i]*(dyn->prstar + 0.5*A_vec[i]*dQdprstar_vec[i]) + dA_vec[i]*Q_vec[i] + A_vec[i]*dQ_vec[i];
        Delta   = SQ(b_coeff) - 4*a_coeff*c_coeff;   /* Delta of the quadratic equation */
              
        /* sol_p = (-b_coeff + sqrt(Delta))/(2*a_coeff); */  /* Plus solution - Unphysical */
        sol_m = (-b_coeff - sqrt(Delta))/(2*a_coeff);  /* Minus solution of the quadratic equation */
        dyn->pphi = sol_m;                             
        
        /** Note: prstar and G functions do not change at even orders 
          G does not change because of the chosen gauge,     
          which eliminates the dependence of G from pphi).
        */
                     
      } //END IF-ELSE parity
      
      /** New Hamiltonians */
      if(usespins) {
        
        eob_ham_s(nu, dyn->r, rc_vec[i], drc_dr_vec[i], 0., dyn->pphi, dyn->prstar, EOBPars->S, EOBPars->Sstar, chi1, chi2, X1, X2, EOBPars->aK2, c3, A_vec[i], dA_vec[i], 0., Q_vec[i], dQ_vec[i], dQdprstar_vec[i], 0., d2Qdprstar2_vec[i],
                  &H,               /* real EOB Hamiltonian divided by mu=m1m2/(m1+m2) */
                  &Heff_vec[i],     /* effective EOB Hamiltonian (divided by mu). Heff coincides with Heff_orb for the non-spinning case */
                  &Heff_orb_vec[i],
                  &dHeff_dr,        /* drvt Heff,r      */
                  &dHeff_dprstar,   /* drvt Heff,prstar */
                  &dHeff_dpphi,     /* drvt Heff,pphi   */
                  &d2Heff_dprstar20,
		  NULL);
        
        E_vec[i] = nu*H;
        
      } else {
        
        eob_ham(nu, dyn->r, dyn->pphi, dyn->prstar, A_vec[i], dA_vec[i], Q_vec[i], dQ_vec[i], dQdprstar_vec[i],
                &H,               /* real EOB Hamiltonian divided by mu=m1m2/(m1+m2) */
                &Heff_orb_vec[i], /* effective EOB Hamiltonian (divided by mu). Heff coincides with Heff_orb for the non-spinning case */
                &dHeff_dr,        /* drvt Heff,r      */
                &dHeff_dprstar,   /* drvt Heff,prstar */
                &dHeff_dpphi);    /* drvt Heff,pphi   */
        
        u2      = 1./((dyn->r)*(dyn->r));
        prstar2 = (dyn->prstar)*(dyn->prstar);
        d2Heff_dprstar20 = 0.5*(2. + A_vec[i]*d2Qdprstar2_vec[i])/Heff_orb_vec[i];
        
        Heff_vec[i] = Heff_orb_vec[i]; /* Heff coincides with Heff_orb for the non-spinning case */
        E_vec[i] = nu*H;
        
      }
      
      /** Orbital Frequency */
      dyn->Omg = dHeff_dpphi/E_vec[i];
      
      /** Real Orbital Frequency */
      dyn->Omg_orb = (dyn->pphi*A_vec[i]*uc2_vec[i])/(E_vec[i]*Heff_orb_vec[i]);

      /** ddotr */
      dyn->ddotr = -A_vec[i]/B_vec[i]*dHeff_dr*d2Heff_dprstar20;
      
      /** Time and phase radial derivatives */
      dt_dr_vec[i]   = E_vec[i]/(sqrtAbyB_vec[i]*dHeff_dprstar); /* dt_dr = 1/dr_dt */
      dphi_dr_vec[i] = dyn->Omg*dt_dr_vec[i];                    /* d(phi)_dr = d(phi)_dt*dt_dr */
      
      /* Re-assigning quantities to array elements */
      dyn->data[EOB_PHI][i]    = dyn->phi;
      dyn->data[EOB_PPHI][i]   = dyn->pphi;
      dyn->data[EOB_MOMG][i]   = dyn->Omg;
      dyn->data[EOB_DDOTR][i]  = dyn->ddotr;
      dyn->data[EOB_PRSTAR][i] = dyn->prstar;
      dyn->data[EOB_OMGORB][i] = dyn->Omg_orb;
      //printf("%.10f, %.10f, %.10f, %.10f, \n", dyn->pphi, dyn->Omg, dyn->r, dyn->Omg);
      if (usespins == MODE_SPINS_GENERIC && EOBPars->project_spins){
        /* Updating spins */
        eob_spin_dyn_Sproj_interp(spin, dyn->Omg, &SApara, &SBpara, NULL, NULL, NULL, NULL, 1);
        chi1_grid[i] = SApara/SQ(X1);
        chi2_grid[i] = SBpara/SQ(X2);
        EOBPars->chi1 = chi1_grid[i];
        EOBPars->chi2 = chi2_grid[i];
        set_spin_vars(X1,X2,EOBPars->chi1,EOBPars->chi2, 
          &EOBPars->S1, &EOBPars->S2,
          &EOBPars->a1, &EOBPars->a2,
          &EOBPars->aK, &EOBPars->aK2,
          &EOBPars->S, &EOBPars->Sstar);
      }
    } // END R-GRID FOR
    
    /* Computing derivatives of the momenta */
    if (parity) D0(dyn->data[EOB_PRSTAR],-dr, size, dprstar_dr_vec);
    else        D0(dyn->data[EOB_PPHI],-dr, size, dpphi_dr_vec);

  } // END PA-CORRECTIONS FOR

  /*
   * Computing integrals for time and phase
   */
  
  /** Compute time */
  cumint3(dt_dr_vec, dyn->data[EOB_RAD], size, dyn->time);
  
  /* Set last value for evolution */
  dyn->t = dyn->time[size-1];
  
  /** Compute orbital phase */
  cumint3(dphi_dr_vec, dyn->data[EOB_RAD], size, dyn->data[EOB_PHI]);
 
  /** Update spins */
  // if (usespins == MODE_SPINS_GENERIC && EOBPars->project_spins){
  //   EOBPars->chi1 = chi1_grid[size-1];
  //   EOBPars->chi2 = chi2_grid[size-1];
  //   set_spin_vars(X1,X2,EOBPars->chi1,EOBPars->chi2, 
  //       &EOBPars->S1, &EOBPars->S2,
  //       &EOBPars->a1, &EOBPars->a2,
  //       &EOBPars->aK, &EOBPars->aK2,
  //       &EOBPars->S, &EOBPars->Sstar);
  // }

  /* Free memory */
  for (int v=0; v < nv; v++)
    free(buffer[v]);
  if (usespins == MODE_SPINS_GENERIC && EOBPars->project_spins){
    free(chi1_grid);
    free(chi2_grid);
  }

  return OK;
}

/* ===========================================================================
 *  Fold-aware PA->ODE switch
 *  The PA inspiral must stop above the last stable orbit (LSO): the odd PA orders divide by dpphi/dr along
 *  the circular sequence, which vanishes there. With fold_switch the PA leg stops at
 *      r_sw = rLSO + sqrt(-x_L/6) Delta_r_P,   x_L = fold_switch_xL < 0,
 *  a fixed depth into the transition region in its own scale. The functions below compute rLSO, pphiLSO
 *  and Delta_r_P from the quasi-circular r.h.s. itself (finite differences), so every nu-dressed term enters.
 *
 *  Circular orbits are the zeros of the conservative prstardot at prstar = 0; at the LSO this sequence folds,
 *  prstardot = d(prstardot)/dr = 0. Near the fold, with the flux Fphi held at its LSO value,
 *      d2(delta_r)/dt2 = -A_P delta_r^2 - B_P (t - tLSO),   delta_r = r - rLSO,
 *      A_P = -(1/2) (drdot/dprstar) (d2 prstardot/dr2) > 0,   B_P = -(drdot/dprstar) (d prstardot/dpphi) Fphi > 0,
 *  and delta_r = -Delta_r_P X, t - tLSO = t_P T with
 *      Delta_r_P = 6^(3/5) A_P^(-3/5) B_P^(2/5) ~ nu^(2/5),   t_P = (6/(A_P B_P))^(1/5) ~ nu^(-1/5)
 *  give the first Painleve transcendent X'' = 6 X^2 + T (Ori & Thorne 2000; Buonanno & Damour 2000;
 *  Compere & Kuchler 2021). x_L is the Painleve time of the handover: the inspiral follows X = -sqrt(-T/6)
 *  until T = O(1). Quasi-circular r.h.s. only (eob_dyn_rhs, eob_dyn_rhs_s), where prstardot with noflx set
 *  is exactly the conservative part.
 *  Ref.: Eq. (6) and "Fold-aware switch to the ODEs" (End Matter) of arXiv:2611.32840; the default
 *  fold_switch_xL = -8 gives the handover coefficient c_sw = sqrt(-x_L/6) = sqrt(4/3) of that paper.
 */


/* The finite-difference steps FOLD_STEP_* and the Newton parameters FOLD_NEWTON_* are in TEOBResumS.h */

/** Failure reporting for eob_fold_compute: every failure is recoverable (EOBRun falls back to the standard
    PA range), so all honor errors_to_warnings. Returns ERROR_ROOTFINDER. */
static int eob_fold_fail(const char *msg)
{
  if (EOBPars->errors_to_warnings) {
    char buf[STRLEN];
    snprintf(buf, sizeof(buf), "eob_fold_compute: %s", msg);
    PRWARN(buf);
  } else {
    printf("ERROR(eob_fold_compute): %s\n", msg);
  }
  return ERROR_ROOTFINDER;
}

/**
 * Function: eob_fold_rhs
 * ----------------------
 *   Evaluate the active r.h.s. p_eob_dyn_rhs at a probe point, without
 *   disturbing the ongoing evolution: dyn->store is forced off and
 *   dyn->noflx and dyn->ode_stop are saved and restored.
 *
 *   @param[in]  r      : radius
 *   @param[in]  prstar : tortoise radial momentum
 *   @param[in]  pphi   : angular momentum
 *   @param[in]  noflx  : if 1, dpphi/dt is zeroed (conservative system)
 *   @param[in]  dyn    : Dynamics structure (used read-only up to flags)
 *   @param[out] dy     : r.h.s. values, indexed by EOB_EVOLVE_*
 */
static void eob_fold_rhs(double r, double prstar, double pphi, int noflx,
                         Dynamics *dyn, double dy[EOB_EVOLVE_NVARS])
{
  double y[EOB_EVOLVE_NVARS];
  const int store_save  = dyn->store;
  const int noflx_save  = dyn->noflx;
  const bool stop_save  = dyn->ode_stop;

  y[EOB_EVOLVE_RAD]    = r;
  y[EOB_EVOLVE_PHI]    = 0.; /* the r.h.s. does not depend on phi */
  y[EOB_EVOLVE_PRSTAR] = prstar;
  y[EOB_EVOLVE_PPHI]   = pphi;

  dyn->store = 0;
  dyn->noflx = noflx;
  p_eob_dyn_rhs(0., y, dy, dyn);

  dyn->store    = store_save;
  dyn->noflx    = noflx_save;
  dyn->ode_stop = stop_save;
}

/**
 * Function: eob_fold_prstardot
 * ----------------------------
 *   Conservative part of dprstar/dt from the r.h.s. (noflx set);
 *   proportional to -dHeff/dr, so its zeroes at prstar = 0 are the
 *   circular orbits. Quasi-circular r.h.s. only.
 *
 *   @param[in]  r      : radius
 *   @param[in]  prstar : tortoise radial momentum
 *   @param[in]  pphi   : angular momentum
 *   @param[in]  dyn    : Dynamics structure
 *
 *   @return conservative dprstar/dt
 */
double eob_fold_prstardot(double r, double prstar, double pphi, Dynamics *dyn)
{
  double dy[EOB_EVOLVE_NVARS];
  eob_fold_rhs(r, prstar, pphi, 1, dyn, dy);
  return dy[EOB_EVOLVE_PRSTAR];
}

/**
 * Function: eob_fold_rdot
 * -----------------------
 *   dr/dt from the r.h.s.; odd in prstar, vanishes on circular orbits.
 *
 *   @param[in]  r      : radius
 *   @param[in]  prstar : tortoise radial momentum
 *   @param[in]  pphi   : angular momentum
 *   @param[in]  dyn    : Dynamics structure
 *
 *   @return dr/dt
 */
static double eob_fold_rdot(double r, double prstar, double pphi, Dynamics *dyn)
{
  double dy[EOB_EVOLVE_NVARS];
  eob_fold_rhs(r, prstar, pphi, 1, dyn, dy);
  return dy[EOB_EVOLVE_RAD];
}

/**
 * Function: eob_fold_Fphi
 * -----------------------
 *   Angular momentum flux Fphi = dpphi/dt at full strength from the
 *   r.h.s. (< 0).
 *
 *   @param[in]  r      : radius
 *   @param[in]  prstar : tortoise radial momentum
 *   @param[in]  pphi   : angular momentum
 *   @param[in]  dyn    : Dynamics structure
 *
 *   @return Fphi
 */
double eob_fold_Fphi(double r, double prstar, double pphi, Dynamics *dyn)
{
  double dy[EOB_EVOLVE_NVARS];
  eob_fold_rhs(r, prstar, pphi, 0, dyn, dy);
  return dy[EOB_EVOLVE_PPHI];
}

/**
 * Function: eob_fold_dprstardot_dr
 * --------------------------------
 *   d(prstardot)/dr at fixed (prstar, pphi); 4th-order central 5-point
 *   stencil. Proportional to -d2Heff/dr2 at circular points.
 *
 *   @param[in]  r      : radius
 *   @param[in]  prstar : tortoise radial momentum
 *   @param[in]  pphi   : angular momentum
 *   @param[in]  hstep  : finite-difference step in r
 *   @param[in]  dyn    : Dynamics structure
 *
 *   @return d(prstardot)/dr
 */
double eob_fold_dprstardot_dr(double r, double prstar, double pphi,
                              double hstep, Dynamics *dyn)
{
  const double fm2 = eob_fold_prstardot(r - 2.*hstep, prstar, pphi, dyn);
  const double fm1 = eob_fold_prstardot(r -    hstep, prstar, pphi, dyn);
  const double fp1 = eob_fold_prstardot(r +    hstep, prstar, pphi, dyn);
  const double fp2 = eob_fold_prstardot(r + 2.*hstep, prstar, pphi, dyn);
  return (fm2 - 8.*fm1 + 8.*fp1 - fp2)/(12.*hstep);
}

/**
 * Function: eob_fold_d2prstardot_dr2
 * ----------------------------------
 *   d2(prstardot)/dr2 at fixed (prstar, pphi); 4th-order central 5-point
 *   stencil.
 *
 *   @param[in]  r      : radius
 *   @param[in]  prstar : tortoise radial momentum
 *   @param[in]  pphi   : angular momentum
 *   @param[in]  hstep  : finite-difference step in r
 *   @param[in]  dyn    : Dynamics structure
 *
 *   @return d2(prstardot)/dr2
 */
double eob_fold_d2prstardot_dr2(double r, double prstar, double pphi,
                                double hstep, Dynamics *dyn)
{
  const double fm2 = eob_fold_prstardot(r - 2.*hstep, prstar, pphi, dyn);
  const double fm1 = eob_fold_prstardot(r -    hstep, prstar, pphi, dyn);
  const double f0  = eob_fold_prstardot(r,            prstar, pphi, dyn);
  const double fp1 = eob_fold_prstardot(r +    hstep, prstar, pphi, dyn);
  const double fp2 = eob_fold_prstardot(r + 2.*hstep, prstar, pphi, dyn);
  return (-fm2 + 16.*fm1 - 30.*f0 + 16.*fp1 - fp2)/(12.*hstep*hstep);
}

/**
 * Function: eob_fold_dprstardot_dpphi
 * -----------------------------------
 *   d(prstardot)/dpphi at fixed (r, prstar); 4th-order central 5-point
 *   stencil.
 *
 *   @param[in]  r      : radius
 *   @param[in]  prstar : tortoise radial momentum
 *   @param[in]  pphi   : angular momentum
 *   @param[in]  hstep  : finite-difference step in pphi
 *   @param[in]  dyn    : Dynamics structure
 *
 *   @return d(prstardot)/dpphi
 */
double eob_fold_dprstardot_dpphi(double r, double prstar, double pphi,
                                 double hstep, Dynamics *dyn)
{
  const double fm2 = eob_fold_prstardot(r, prstar, pphi - 2.*hstep, dyn);
  const double fm1 = eob_fold_prstardot(r, prstar, pphi -    hstep, dyn);
  const double fp1 = eob_fold_prstardot(r, prstar, pphi +    hstep, dyn);
  const double fp2 = eob_fold_prstardot(r, prstar, pphi + 2.*hstep, dyn);
  return (fm2 - 8.*fm1 + 8.*fp1 - fp2)/(12.*hstep);
}

/**
 * Function: eob_fold_drdot_dprstar
 * --------------------------------
 *   d(rdot)/dprstar at fixed (r, pphi); 4th-order central 5-point
 *   stencil. At prstar = 0 this is the prstar -> 0 limit of rdot/prstar
 *   (cf. dr_dtbyprstar in TEOBResumSPostAdiab.c); rdot is odd in prstar,
 *   so the stencil is free of even-order contamination there.
 *
 *   @param[in]  r      : radius
 *   @param[in]  prstar : tortoise radial momentum
 *   @param[in]  pphi   : angular momentum
 *   @param[in]  hstep  : finite-difference step in prstar
 *   @param[in]  dyn    : Dynamics structure
 *
 *   @return d(rdot)/dprstar
 */
double eob_fold_drdot_dprstar(double r, double prstar, double pphi,
                              double hstep, Dynamics *dyn)
{
  const double fm2 = eob_fold_rdot(r, prstar - 2.*hstep, pphi, dyn);
  const double fm1 = eob_fold_rdot(r, prstar -    hstep, pphi, dyn);
  const double fp1 = eob_fold_rdot(r, prstar +    hstep, pphi, dyn);
  const double fp2 = eob_fold_rdot(r, prstar + 2.*hstep, pphi, dyn);
  return (fm2 - 8.*fm1 + 8.*fp1 - fp2)/(12.*hstep);
}

/**
 * Function: eob_fold_omegar2
 * --------------------------
 *   Radial epicyclic frequency squared at (r, prstar = 0, pphi):
 *   omega_r^2 = -(drdot/dprstar)(d(prstardot)/dr). Positive on the stable
 *   circular branch, vanishing at the fold (where it crosses zero with
 *   slope 2 A_P in r - rLSO).
 *
 *   @param[in]  r    : radius
 *   @param[in]  pphi : angular momentum
 *   @param[in]  dyn  : Dynamics structure
 *
 *   @return omega_r^2
 */
double eob_fold_omegar2(double r, double pphi, Dynamics *dyn)
{
  const double drdot_dprstar =
    eob_fold_drdot_dprstar(r, 0., pphi, FOLD_STEP_DPRSTAR, dyn);
  const double dprstardot_dr =
    eob_fold_dprstardot_dr(r, 0., pphi, FOLD_STEP_DR, dyn);
  return -drdot_dprstar*dprstardot_dr;
}

/**
 * Function: eob_fold_kerr_isco_seed
 * ---------------------------------
 *   Test-mass Kerr ISCO (Bardeen-Press-Teukolsky), seed of the Newton solve for the
 *   nu-dressed fold; a = aK = a1 + a2, prograde.
 *
 *   @param[in]  a         : dimensionless spin
 *   @param[out] r_seed    : ISCO radius
 *   @param[out] pphi_seed : ISCO angular momentum
 */
static void eob_fold_kerr_isco_seed(double a, double *r_seed, double *pphi_seed)
{
  const double Z1 = 1. + cbrt(1. - a*a)*(cbrt(1. + a) + cbrt(1. - a));
  const double Z2 = sqrt(3.*a*a + Z1*Z1);
  const double r  = 3. + Z2 - SIGN(1., a)*sqrt((3. - Z1)*(3. + Z1 + 2.*Z2));
  const double sr = sqrt(r);
  *r_seed    = r;
  *pphi_seed = (r*r - 2.*a*sr + a*a)/(pow(r, 0.75)*sqrt(r*sr - 3.*sr + 2.*a));
}

/**
 * Function: eob_fold_compute
 * --------------------------
 *   Locate the nu-dressed fold (rLSO, pphiLSO) by a damped 2D Newton on
 *   (prstardot, d(prstardot)/dr) = 0 at prstar = 0, seeded from the
 *   test-mass Kerr ISCO of aK, and compute the transition scale Delta_r_P
 *   (see the section header). Fails loudly (nonzero return, zeroed outputs)
 *   if the Newton does not converge, if A_P or B_P has an unphysical
 *   sign, or if the fold is not adiabatic for this configuration.
 *
 *   To be called once per configuration (spins, nu), after eob_set_params.
 *   Cost: O(100) r.h.s. evaluations (~0.3 ms).
 *   Ref.: Delta_r_P is Eq. (6) of arXiv:2611.32840.
 *
 *   @param[in]  dyn       : Dynamics structure
 *   @param[out] rLSO      : fold radius
 *   @param[out] pphiLSO   : fold angular momentum
 *   @param[out] Delta_r_P : radial scale of the transition region (> 0)
 *
 *   @return OK, or ERROR_ROOTFINDER on failure
 */
int eob_fold_compute(Dynamics *dyn, double *rLSO, double *pphiLSO,
                     double *Delta_r_P)
{
  const double hr  = FOLD_STEP_DR;
  const double hpp = FOLD_STEP_DPPHI;

  *rLSO = 0.;
  *pphiLSO = 0.;
  *Delta_r_P = 0.;

  double r, pp;
  eob_fold_kerr_isco_seed(EOBPars->use_spins ? EOBPars->aK : 0., &r, &pp);

  /* Damped Newton on F = (prstardot, d(prstardot)/dr)(r, 0, pphi) with
     finite-difference Jacobian; note dF1/dr = F2 exactly at the current
     point. */
  int converged = 0;
  for (int it = 0; it < FOLD_NEWTON_ITMAX; it++) {
    const double F1 = eob_fold_prstardot(r, 0., pp, dyn);
    const double F2 = eob_fold_dprstardot_dr(r, 0., pp, hr, dyn);

    const double J11 = F2;
    const double J12 = eob_fold_dprstardot_dpphi(r, 0., pp, hpp, dyn);
    const double J21 = eob_fold_d2prstardot_dr2(r, 0., pp, hr, dyn);
    /* d(d(prstardot)/dr)/dpphi by differencing the dr stencil in pphi */
    const double J22 = (  eob_fold_dprstardot_dr(r, 0., pp - 2.*hpp, hr, dyn)
                      - 8.*eob_fold_dprstardot_dr(r, 0., pp -    hpp, hr, dyn)
                      + 8.*eob_fold_dprstardot_dr(r, 0., pp +    hpp, hr, dyn)
                      -    eob_fold_dprstardot_dr(r, 0., pp + 2.*hpp, hr, dyn))/(12.*hpp);

    const double det = J11*J22 - J12*J21;
    if (det == 0. || !isfinite(det)) break;

    double dr  = -( J22*F1 - J12*F2)/det;
    double dpp = -(-J21*F1 + J11*F2)/det;
    if (fabs(dr)  > FOLD_NEWTON_MAXSTEP) dr  = SIGN(FOLD_NEWTON_MAXSTEP, dr);
    if (fabs(dpp) > FOLD_NEWTON_MAXSTEP) dpp = SIGN(FOLD_NEWTON_MAXSTEP, dpp);
    r  += dr;
    pp += dpp;

    if (r < 1.5 || !isfinite(r) || !isfinite(pp)) break; /* left the sane region */

    if (fabs(dr) < FOLD_NEWTON_TOL*(1. + fabs(r)) &&
        fabs(dpp) < FOLD_NEWTON_TOL*(1. + fabs(pp))) {
      converged = 1;
      break;
    }
  }

  if (!converged)
    return eob_fold_fail("Newton for the nu-dressed fold did not converge.");

  /* Transition constants at (rLSO, 0, pphiLSO), flux at full strength */
  const double drdot_dprstar    = eob_fold_drdot_dprstar(r, 0., pp, FOLD_STEP_DPRSTAR, dyn);
  const double Fphi             = eob_fold_Fphi(r, 0., pp, dyn);
  const double d2prstardot_dr2  = eob_fold_d2prstardot_dr2(r, 0., pp, hr, dyn);
  const double dprstardot_dpphi = eob_fold_dprstardot_dpphi(r, 0., pp, hpp, dyn);

  const double A_P = -0.5*drdot_dprstar*d2prstardot_dr2;
  const double B_P = -drdot_dprstar*dprstardot_dpphi*Fphi;

  if (!(A_P > 0.) || !(B_P > 0.)) {
    char msg[STRLEN];
    snprintf(msg, sizeof(msg), "unphysical fold constants "
             "A_P = %.6e, B_P = %.6e (both must be > 0).", A_P, B_P);
    return eob_fold_fail(msg);
  }

  const double t_P = pow(6./(A_P*B_P), 0.2);
  const double DrP = pow(6., 0.6)*pow(A_P, -0.6)*pow(B_P, 0.4);

  /* Adiabaticity of the fold: eps_P = |Fphi| t_P / pphiLSO (~ nu^(4/5)), the fractional angular-momentum drift
     over one transition time, measures the first neglected correction of the near-fold expansion; refuse
     eps_P > 0.1 */
  const double eps_P = fabs(Fphi)*t_P/pp;
  if (eps_P > 0.1) {
    char msg[STRLEN];
    snprintf(msg, sizeof(msg), "fold is not adiabatic "
             "(|Fphi| t_P / pphiLSO = %.3e > 0.1); "
             "the transition scaling does not apply here.", eps_P);
    return eob_fold_fail(msg);
  }

  /* Normalization identities (algebraic in the definitions; violations
     signal an implementation bug, not a model property) */
  const double id1 = t_P*t_P*t_P*B_P - DrP;
  const double id2 = A_P*DrP*t_P*t_P - 6.;
  if (fabs(id1) > 1e-10*DrP || fabs(id2) > 1e-10*6.) {
    char msg[STRLEN];
    snprintf(msg, sizeof(msg), "normalization identities violated "
             "(t_P^3 B_P - Delta_r_P = %.3e, A_P Delta_r_P t_P^2 - 6 = %.3e).",
             id1, id2);
    return eob_fold_fail(msg);
  }

  *rLSO = r;
  *pphiLSO = pp;
  *Delta_r_P = DrP;
  return OK;
}

/**
 * Function: eob_fold_wiggle_qa
 * ----------------------------
 *   Scaled epicyclic-wiggle amplitude of the stored dynamics just after the PA->ODE switch: r(t) over
 *   about 3 epicyclic periods 2 pi/omega_r(r_sw) after the switch sample i_sw is fitted with
 *   {1, t, sqrt(tLSO - t)} (the adiabatic drift r - rLSO = Delta_r_P sqrt(-T/6)), and rms(residual)/Delta_r_P
 *   is returned. An error of the handover data excites an epicyclic oscillation that the transition window
 *   amplifies (~e^6) and that no ODE tolerance removes; this metric certifies the handover.
 *
 *   @param[in]  dyn       : Dynamics with stored time/data arrays
 *   @param[in]  rLSO      : fold radius of this configuration
 *   @param[in]  Delta_r_P : transition scale of this configuration
 *   @param[in]  i_sw      : index of the switch sample (last PA point)
 *   @param[out] npts_out  : number of samples used (optional, NULL ok)
 *
 *   @return scaled wiggle rms, or -1. if the metric cannot be formed
 *           (window not sampled, bad omega_r^2, fewer than 12 points, or
 *           no rLSO crossing beyond the window)
 */
double eob_fold_wiggle_qa(Dynamics *dyn, double rLSO, double Delta_r_P,
                          int i_sw, int *npts_out)
{
  if (npts_out) *npts_out = 0;
  if (i_sw < 0 || i_sw >= dyn->size - 1) return -1.;

  const double t_sw    = dyn->time[i_sw];
  const double r_sw    = dyn->data[EOB_RAD][i_sw];
  const double pphi_sw = dyn->data[EOB_PPHI][i_sw];

  const double w2 = eob_fold_omegar2(r_sw, pphi_sw, dyn);
  if (!(w2 > 0.)) return -1.;
  const double Twin = 3.*TwoPi/sqrt(w2);

  int n = 0;
  while (i_sw + n < dyn->size && dyn->time[i_sw + n] <= t_sw + Twin) n++;
  if (n < 12) return -1.; /* ~4 samples per period, the minimum for the
                             fit to see an oscillation at all */
  if (npts_out) *npts_out = n;

  /* the sqrt-drift basis needs tLSO: first crossing of rLSO after the
     switch, linearly interpolated; it must lie beyond the fit window */
  double tLSO = 0.;
  {
    int j = i_sw;
    while (j < dyn->size && dyn->data[EOB_RAD][j] > rLSO) j++;
    if (j >= dyn->size || j == i_sw) return -1.;
    const double f = (dyn->data[EOB_RAD][j-1] - rLSO)
                    /(dyn->data[EOB_RAD][j-1] - dyn->data[EOB_RAD][j]);
    tLSO = dyn->time[j-1] + f*(dyn->time[j] - dyn->time[j-1]);
    if (tLSO <= t_sw + Twin) return -1.;
  }

  /* least-squares fit of r(t) on {1, s, sqrt((tLSO - t)/Twin)} in
     s = 2(t-t_sw)/Twin - 1; normal equations solved by Gauss elimination
     with partial pivoting (well-conditioned: smooth basis, 3 functions,
     s in [-1,1]) */
  const int m = 3;
  double M[3][3] = {{0.}}, b[3] = {0.};
  for (int i = 0; i < n; i++) {
    const double t = dyn->time[i_sw + i];
    const double s = 2.*(t - t_sw)/Twin - 1.;
    double pk[3];
    pk[0] = 1.; pk[1] = s; pk[2] = sqrt((tLSO - t)/Twin);
    for (int j = 0; j < m; j++) {
      b[j] += pk[j]*dyn->data[EOB_RAD][i_sw + i];
      for (int k = 0; k < m; k++) M[j][k] += pk[j]*pk[k];
    }
  }
  for (int col = 0; col < m; col++) {
    int piv = col;
    for (int row = col + 1; row < m; row++)
      if (fabs(M[row][col]) > fabs(M[piv][col])) piv = row;
    if (M[piv][col] == 0.) return -1.;
    if (piv != col) {
      for (int k = 0; k < m; k++) SWAPTRS(M[col][k], M[piv][k]);
      SWAPTRS(b[col], b[piv]);
    }
    for (int row = col + 1; row < m; row++) {
      const double fac = M[row][col]/M[col][col];
      for (int k = col; k < m; k++) M[row][k] -= fac*M[col][k];
      b[row] -= fac*b[col];
    }
  }
  double coef[3] = {0.};
  for (int row = m - 1; row >= 0; row--) {
    double s = b[row];
    for (int k = row + 1; k < m; k++) s -= M[row][k]*coef[k];
    coef[row] = s/M[row][row];
  }

  double ss = 0.;
  for (int i = 0; i < n; i++) {
    const double t = dyn->time[i_sw + i];
    const double s = 2.*(t - t_sw)/Twin - 1.;
    double pk[3];
    pk[0] = 1.; pk[1] = s; pk[2] = sqrt((tLSO - t)/Twin);
    double fit = 0.;
    for (int j = 0; j < m; j++) fit += coef[j]*pk[j];
    const double res = dyn->data[EOB_RAD][i_sw + i] - fit;
    ss += res*res;
  }
  return sqrt(ss/n)/Delta_r_P;
}
