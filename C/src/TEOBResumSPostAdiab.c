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
    if(EOBPars->spin_interp_domain==0){
      errorexit("Spin projection in t + PA is not yet implemented. Retry with spin_interp_domain = 1");
    } else {
      chi1_grid = (double*) malloc(size * sizeof(double));
      chi2_grid = (double*) malloc(size * sizeof(double));
      for(int i=0; i<size;i++){
        const double Omg_i = pow(r0 - i*dr, -1.5);
        eob_spin_dyn_Sproj_interp(spin, Omg_i, &SApara, &SBpara, NULL, NULL, NULL, NULL, 1);
        chi1_grid[i] = SApara/SQ(X1);
        chi2_grid[i] = SBpara/SQ(X2);
      }
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
        d2Heff_dprstar20 = d2Heff_dprstar20 = 0.5*(2. + A_vec[i]*d2Qdprstar2_vec[i])/Heff_orb_vec[i];
        
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
