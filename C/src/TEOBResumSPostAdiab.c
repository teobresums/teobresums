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
#include <string.h>
#include "TEOBResumS.h"

/** Post-adiabatic dynamics */
int eob_dyn_Npostadiabatic(Dynamics *dyn, double r0)
{
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
  const double z3    = 2.0*nu*(4.0-3.0*nu);
  const int usetidal = dyn->use_tidal;
  const int usespins = dyn->use_spins;
    
  const int size = 100;//dyn->size;
  printf("\nI am using size: %d\n",size);

  const int Npa = par_get_i("postadiabatic_dynamics_N");

  /* Mem for quantities to be used. FIXME: Can we avoid allocating the memory?
     Can we use already-defined C-structures here? */
  double *A_vec                  = (double*)malloc(size * sizeof (double));
  double *dA_vec                 = (double*)malloc(size * sizeof (double));
  double *d2A_vec                = (double*)malloc(size * sizeof (double));
  double *B_vec                  = (double*)malloc(size * sizeof (double));
  double *dB_vec                 = (double*)malloc(size * sizeof (double));
  double *sqrtAbyB_vec           = (double*)malloc(size * sizeof (double));
  double *rc_vec                 = (double*)malloc(size * sizeof (double));
  double *drc_dr_vec             = (double*)malloc(size * sizeof (double));
  double *d2rc_dr_vec            = (double*)malloc(size * sizeof (double));
  double *uc2_vec                = (double*)malloc(size * sizeof (double));
  double *duc_dr_vec             = (double*)malloc(size * sizeof (double));
  double *dAuc2_dr_vec           = (double*)malloc(size * sizeof (double));
  double *G_vec                  = (double*)malloc(size * sizeof (double));
  double *dG_dr_vec              = (double*)malloc(size * sizeof (double));
  double *dG_dprstar_vec         = (double*)malloc(size * sizeof (double));
  double *dG_dprstarbyprstar_vec = (double*)malloc(size * sizeof (double));
  double *dphi_dr_vec            = (double*)malloc(size * sizeof (double));
  double *dpphi_dr_vec           = (double*)malloc(size * sizeof (double));
  double *dprstar_dr_vec         = (double*)malloc(size * sizeof (double));
  double *H_vec                  = (double*)malloc(size * sizeof (double));
  double *Heff_vec               = (double*)malloc(size * sizeof (double));
  double *Heff_orb_vec           = (double*)malloc(size * sizeof (double));
  double *dr_dt_vec              = (double*)malloc(size * sizeof (double));
  double *dt_dr_vec              = (double*)malloc(size * sizeof (double));
  double *G0_vec                 = (double*)malloc(size * sizeof (double));
  double *dG_dr0_vec             = (double*)malloc(size * sizeof (double));
  double *dt_dr_mat_vec          = (double*)malloc(size * sizeof (double));
  double *dr_dt_mat_vec          = (double*)malloc(size * sizeof (double));
  double *r_mat_vec              = (double*)malloc(size * sizeof (double));
  double *prstar_mat_vec         = (double*)malloc(size * sizeof (double));
  double *pphi_mat_vec           = (double*)malloc(size * sizeof (double));
  double *omg_mat_vec            = (double*)malloc(size * sizeof (double));
  double *t_mat_vec              = (double*)malloc(size * sizeof (double));
  double *phi_mat_vec            = (double*)malloc(size * sizeof (double));
  
  double ggm[14];
    
  double a_coeff, b_coeff, c_coeff, Delta, sol_p, sol_m, j02, uc, dHeff_dpphi, r_omg4, Omg5, dr_dtbypr, dHeff_dprstar,
    dHeff_dprstarbyprstar, ddotr_fake, prstar_fake, x, jhat, psi, r_omg, sqrtW, v_phi, Fphi, dr_dtbyprstar, prstar4, Heff_orb_f, Heff_f, E_f;

  /* Compute radius of inflection point of Pr* */
  //TODO implementing a robust stopping condition (which)
  r0        = 22;
  double ri = 12.;

  /* Build a uniform grid and compute circular dynamics */
  const double dr = 0.1;//FIXME(r0 - ri)/(size-1);
  /*Better to set fixed resolution. For the moment the derivative oscillates when there are too many points (higher order derivative or maximum resolution)*/
    
  for (int i = 0; i < size; i++)
    {
      
      //FIXME dyn->data[EOB_RAD] forse non e' quello da usare, in altri punti usa dyn->r. Io non voglio che dyn->data cambi, quindi non mi piace passargli il puntatore.
      dyn->r = r0 - i*dr;
      
      /* Computing metric functions and centrifugal radius */
      if(usespins)
        {
          eob_metric_s(dyn->r,dyn, &A_vec[i], &B_vec[i], &dA_vec[i], &d2A_vec[i], &dB_vec[i]);
            
          eob_dyn_s_get_rc(dyn->r, nu, a1, a2, aK2, C_Q1, C_Q2, usetidal, &rc_vec[i], &drc_dr_vec[i], &d2rc_dr_vec[i]);
            
          eob_dyn_s_GS(dyn->r, rc_vec[i], drc_dr_vec[i], aK2, 0.0, 0.0, nu, chi1, chi2, X1, X2, c3, ggm);
          G_vec[i]                  = ggm[2] *S+ggm[3] *Sstar;    // tildeG = GS*S+GSs*Ss
          dG_dr_vec[i]              = ggm[6] *S+ggm[7] *Sstar;
          dG_dprstar_vec[i]         = ggm[4] *S+ggm[5] *Sstar;
          dG_dprstarbyprstar_vec[i] = ggm[10]*S+ggm[11]*Sstar;
        }
      else
        {
          eob_metric(dyn->r ,dyn, &A_vec[i], &B_vec[i], &dA_vec[i], &d2A_vec[i], &dB_vec[i]);
            
          rc_vec[i]                 = dyn->r; //Nonspinning case: rc = r
          drc_dr_vec[i]             = 1;
          d2rc_dr_vec[i]            = 0;
            
          G_vec[i]                  = 0.0;
          dG_dr_vec[i]              = 0.0;
          dG_dprstar_vec[i]         = 0.0;
          dG_dprstarbyprstar_vec[i] = 0.0;
        }
        
      /* Auxiliary variables*/
      sqrtAbyB_vec[i] = sqrt(A_vec[i]/B_vec[i]);
      uc              = 1./rc_vec[i];
      uc2_vec[i]      = uc*uc;
      duc_dr_vec[i]   = -uc2_vec[i]*drc_dr_vec[i];
      dAuc2_dr_vec[i] = uc2_vec[i]*(dA_vec[i]-2*A_vec[i]*uc*drc_dr_vec[i]);
      
      /**************************
       * Circular approximation *
       **************************/
      // Computing the circular angular momentum by solving eq. (A15) of TEOBResumS paper (which is equivalent to solve eq.(4)=0 of arXiv:1805.03891). The procedure to choose the physical solution of the quadratic equation is effective but not understood.  FIXME

      if (usespins)
	{
	  a_coeff = dAuc2_dr_vec[i]*dAuc2_dr_vec[i] - 4*A_vec[i]*uc2_vec[i]*dG_dr_vec[i]*dG_dr_vec[i];  /* First coefficient of the quadratic equation a*x^2+b*x+c=0 */
	  b_coeff = 2*dA_vec[i]*dAuc2_dr_vec[i] - 4*A_vec[i]*dG_dr_vec[i]*dG_dr_vec[i];                 /* Second coefficient of the quadratic equation */
	  c_coeff = dA_vec[i]*dA_vec[i]; /* Third coefficient of the quadratic equation */
	  Delta   = b_coeff*b_coeff - 4*a_coeff*c_coeff ; /* Delta of the quadratic equation */
	  
	  sol_p   = (-b_coeff + sqrt(Delta))/(2*a_coeff); /* Plus  solution of the quadratic equation */
	  sol_m   = (-b_coeff - sqrt(Delta))/(2*a_coeff); /* Minus solution of the quadratic equation */
	  
	  /* Effective prescription: If the Tilde G function is negative, take the positive solution and viceversa. */
	  if (G_vec[i] < 0)
	    {j02 = sol_p;}
	  else
	    {j02 = sol_m;}
	}
      else
	{
	  a_coeff = dAuc2_dr_vec[i];
	  b_coeff = dA_vec[i];
	  
	  j02 = -b_coeff/a_coeff;
        }
      
      /* Define momenta in the circular orbit approximation */
      dyn->pphi                = sqrt(j02);
      dyn->prstar              = 0.0;
      dprstar_dr_vec[i]        = 0.0;
        
      // Circular Hamiltonians, ref: arXiv: 1406.6913
      if(usespins)
      {
        eob_ham_s(nu,dyn->r,rc_vec[i],drc_dr_vec[i],dyn->pphi,dyn->prstar,S,Sstar,chi1,chi2,X1,X2,aK2,c3,A_vec[i],dA_vec[i],
        &H_vec[i],        /* real EOB Hamiltonian divided by mu=m1m2/(m1+m2) */
        &Heff_vec[i],     /* effective EOB Hamiltonian (divided by mu)       */
        &Heff_orb_vec[i],
        NULL,             /* drvt Heff,r      */
        NULL,             /* drvt Heff,prstar */
        &dHeff_dpphi,      /* drvt Heff,pphi   */
        NULL);
      }
      else
      {
        //NON spinning hamiltonian
        eob_ham(nu, dyn->r, dyn->pphi, dyn->prstar, A_vec[i], dA_vec[i],
        &H_vec[i],        /* real EOB Hamiltonian divided by mu=m1m2/(m1+m2) */
	&Heff_orb_vec[i], /* effective EOB Hamiltonian (divided by mu). */
        NULL,             /* drvt Heff,r      */
        NULL,             /* drvt Heff,prstar */
        &dHeff_dpphi);     /* drvt Heff,pphi   */

	Heff_vec[i] = Heff_orb_vec[i]; /* Heff coincides with Heff_orb for the non-spinning case */
      }
        
      // Circular orbital frequency
      dyn->Omg     = 1./(H_vec[i]*nu)*dHeff_dpphi;

      /* Defining circular quantities for the flux calculation.
	 Must not be overwritten in successive iterations, thus
         we define separate quantities with the subscripts 0. */
      G0_vec[i]        = G_vec[i];
      dG_dr0_vec[i]    = dG_dr_vec[i];
	
      dyn->data[EOB_RAD][i]    = dyn->r;
      dyn->data[EOB_PPHI][i]   = dyn->pphi;
      dyn->data[EOB_PRSTAR][i] = dyn->prstar;
      dyn->data[EOB_MOMG][i]   = dyn->Omg;
    }
  // END r-GRID FOR

  // Computing angular momentum derivative
  D0(dyn->data[EOB_PPHI],-dr, size, dpphi_dr_vec); /* dJ0/dr */

  /***************
   * PA dynamics *
   ***************/

  // FIXME
  // maximize use of routines already coded e.g. flux routines, hamiltonian, etc
  // check code in *Dynamics.c and *InitialConditions.c
  // store final result into dyn->data[<variable_index>][i]
  // note Dynamics structure has 'storage' arrays dyn->data and point-wise vars
  
  bool prstar_bool = 1;
    
  /* For on PA orders */
  for (int n = 1; n <= Npa; n++)
    {
 
      if (n%2==0) /*Separating even and odd orders*/
        {prstar_bool = 0;}
      else
        {prstar_bool = 1;}

      /* For on r-grid */
      for (int i = 0; i < size; i++)
        {
	  /* Setting loop variables to help reader */
	  dyn->r       = dyn->data[EOB_RAD][i];
	  dyn->phi     = dyn->data[EOB_PHI][i];
	  dyn->pphi    = dyn->data[EOB_PPHI][i];
	  dyn->Omg     = dyn->data[EOB_MOMG][i];
	  dyn->ddotr   = dyn->data[EOB_DDOTR][i]; // Not used for the moment
	  dyn->prstar  = dyn->data[EOB_PRSTAR][i];
	  dyn->Omg_orb = dyn->data[EOB_OMGORB][i];
            
	  if (prstar_bool)  // Odd PA orders - prstar corrections
            {
	      /*****************************
	       * Calculating the flux Fphi *
	       *****************************/
	      //FIXME USE C-routines, jhat etc. are already present inside dynamics

	      if(usespins)
                {
		  // Variables for which Kepler's law is still valid
		  Heff_orb_f = sqrt(A_vec[i]*(1.0 + SQ(dyn->pphi)*uc2_vec[i]));
		  Heff_f     = G0_vec[i]*dyn->pphi + Heff_orb_f;
		  E_f        = sqrt(1 + 2*nu*(Heff_f - 1));
		  psi        = (duc_dr_vec[i] + dG_dr0_vec[i]*rc_vec[i]*sqrt(A_vec[i]/(SQ(dyn->pphi)) + A_vec[i]*uc2_vec[i])/A_vec[i])/(-0.5*dA_vec[i]);
		  r_omg      = 1.0/cbrt(SQ(((1./sqrt(rc_vec[i]*rc_vec[i]*rc_vec[i]*psi))+G0_vec[i])/(E_f)));
		  v_phi      = r_omg*dyn->Omg;
		  x          = v_phi*v_phi;
		  jhat       = dyn->pphi/(r_omg*v_phi);
		    
		  ddotr_fake  = 0.0; //FIXME: To be changed when considering NQCs (dyn->ddotr or dyn->data[EOB_DDOTR][i])
		  prstar_fake = 0.0; //FIXME: To be changed to the true prstar value when considering NQCs.

		  Fphi = eob_flx_Flux_s(x,dyn->Omg,r_omg, nu*H_vec[i], Heff_vec[i],jhat,dyn->r,prstar_fake, ddotr_fake, dyn);
                }
	      //END-IF spins
	      else
                {
		  //NON-spinning
		  psi   = 2.*(1.0 + 2.0*nu*(Heff_orb_vec[i] - 1.0))/(SQ(dyn->r)*dA_vec[i]);
		  r_omg = dyn->r*cbrt(psi);
		  v_phi = r_omg*dyn->Omg;
		  x     = v_phi * v_phi;
		  jhat  = dyn->pphi/(r_omg*v_phi);
                  
		  ddotr_fake  = 0.0; //FIXME: To be changed when considering NQCs (dyn->ddotr or dyn->data[EOB_DDOTR][i])
		  prstar_fake = 0.0; //FIXME: To be changed to the true prstar value when considering NQCs.

		  Fphi = eob_flx_Flux(x,dyn->Omg,r_omg, nu*H_vec[i], Heff_vec[i],jhat,dyn->r,prstar_fake, ddotr_fake, dyn);
                }


	      /**********************
	       * Calculating prstar *
	       **********************/

	      dHeff_dprstarbyprstar = dyn->pphi*dG_dprstarbyprstar_vec[i] + 1./Heff_orb_vec[i]*(1+2*z3*A_vec[i]*uc2_vec[i]*SQ(dyn->prstar));
	      dr_dtbyprstar         = sqrtAbyB_vec[i]/(nu*H_vec[i])*dHeff_dprstarbyprstar;

	      dyn->prstar = Fphi/dpphi_dr_vec[i]/dr_dtbyprstar; /* Computing first PA using the approximation detailed above A19 of TEOBResumS paper and Hamilton's equations. */
		
	      /***************************************
	       * p_phi does not change at odd orders *
	       ***************************************/

	      /*********************
	       * New GGM functions *
	       *********************/
	      eob_dyn_s_GS(dyn->r, rc_vec[i], drc_dr_vec[i], aK2, dyn->prstar, 0.0, nu, chi1, chi2, X1, X2, c3, ggm);
                
	      G_vec[i]                  = ggm[2] *S+ggm[3] *Sstar;    // Tilde G
	      dG_dr_vec[i]              = ggm[6] *S+ggm[7] *Sstar;
	      dG_dprstar_vec[i]         = ggm[4] *S+ggm[5] *Sstar;
	      dG_dprstarbyprstar_vec[i] = ggm[10]*S+ggm[11]*Sstar;
            }
	  //END IF prstar

	  else
            {
	      /********************
	       * Calculating pphi *
	       ********************/
	      prstar4 = SQ(SQ(dyn->prstar));
	      a_coeff = dAuc2_dr_vec[i];           /* coefficients of the quadratic equation a*x^2+b*x+c=0 */
	      b_coeff = 2*Heff_orb_vec[i]*(dG_dr_vec[i] + dG_dprstar_vec[i]*dprstar_dr_vec[i]);
	      c_coeff = dA_vec[i] + 2*dyn->prstar*dprstar_dr_vec[i]*(1+2*z3*A_vec[i]*uc2_vec[i]*SQ(dyn->prstar)) + z3*dAuc2_dr_vec[i]*prstar4;
	      Delta   = SQ(b_coeff) - 4*a_coeff*c_coeff;     // Delta of the quadratic equation
                
	      //sol_p = (-b_coeff + sqrt(Delta))/(2*a_coeff) // Plus solution of the quadratic equation FIXME: understand why.
	      sol_m = (-b_coeff - sqrt(Delta))/(2*a_coeff);  // Minus solution of the quadratic equation
                
	      dyn->pphi = sol_m;                // Choosing minus solution - To be understood FIXME
                
	      /*******************************************************
	       * prstar and G functions do not change at even orders *
               * (G does not change because of the chosen gauge,     *
               * which eliminates the dependence of G from pphi).    *
	       *******************************************************/                
	    }
	  //END IF-ELSE prstar-pphi
            
	  /********************
	   * New Hamiltonians *
	   ********************/
	  if(usespins)
            {
	      eob_ham_s(nu,dyn->r,rc_vec[i],drc_dr_vec[i],dyn->pphi,dyn->prstar,S,Sstar,chi1,chi2,X1,X2,aK2,c3,A_vec[i],dA_vec[i],
			&H_vec[i],        /* real EOB Hamiltonian divided by mu=m1m2/(m1+m2) */
			&Heff_vec[i],     /* effective EOB Hamiltonian (divided by mu). Heff coincides with Heff_orb for the non-spinning case */
			&Heff_orb_vec[i],
			NULL,             /* drvt Heff,r      */
			&dHeff_dprstar,   /* drvt Heff,prstar */
			&dHeff_dpphi,      /* drvt Heff,pphi   */
			NULL);
            }
	  else
            {
	      //NON spinning hamiltonian
	      eob_ham(nu, dyn->r, dyn->pphi, dyn->prstar, A_vec[i], dA_vec[i],
		      &H_vec[i],        /* real EOB Hamiltonian divided by mu=m1m2/(m1+m2) */
		      &Heff_orb_vec[i], /* effective EOB Hamiltonian (divided by mu). Heff coincides with Heff_orb for the non-spinning case */
		      NULL,             /* drvt Heff,r      */
		      &dHeff_dprstar,   /* drvt Heff,prstar */
		      &dHeff_dpphi);     /* drvt Heff,pphi   */

	      Heff_vec[i] = Heff_orb_vec[i]; /* Heff coincides with Heff_orb for the non-spinning case */
            }
      
	  /*********************
	   * Orbital Frequency *
	   *********************/
	  dyn->Omg = 1./(H_vec[i]*nu)*dHeff_dpphi;

	  /*********
	   * dr_dt *
	   *********/
	  //    dyn->dy[EOB_EVOLVE_RAD]  FIXME: to be checked

	  dt_dr_vec[i]   = 1.0/((sqrtAbyB_vec[i]*dHeff_dprstar)/(nu*H_vec[i])); /* dt_dr = 1/dr_dt */
	  dphi_dr_vec[i] = dyn->Omg*dt_dr_vec[i];                               /* d(phi)_dr = d(phi)_dt*dt_dr */
      
	  /* Re-assigning quantities to array elements */

	  dyn->data[EOB_PHI][i]    = dyn->phi;
	  dyn->data[EOB_PPHI][i]   = dyn->pphi;
	  dyn->data[EOB_MOMG][i]   = dyn->Omg;
	  dyn->data[EOB_DDOTR][i]  = dyn->ddotr;
	  dyn->data[EOB_PRSTAR][i] = dyn->prstar;
	  dyn->data[EOB_OMGORB][i] = dyn->Omg_orb;
        }
      // END R-GRID FOR

      /*Computing derivatives of the momenta*/
      if (prstar_bool)
	{D0(dyn->data[EOB_PRSTAR],-dr, size, dprstar_dr_vec);}
      else
	{D0(dyn->data[EOB_PPHI],-dr, size, dpphi_dr_vec);}
    }
  // END PA-CORRECTIONS FOR

  printf("\n\n\nI am debugging Post-adiab!!!!!!!!!!!\n\n\n");
  char ro_string[256]; //size of the number

  sprintf(ro_string, "r0_%f", r0);
  printf("I am computing the dynamics with r0: %f\n", r0);
  char outputadiab[256]     = "Post_adiab_";
  strcat(outputadiab, ro_string);
  strcat(outputadiab, ".dat"); 
 
  FILE* Post_adiab_debug = fopen(outputadiab, "w");
  for (int kt = 0; kt < size; kt++)
    {
      fprintf(Post_adiab_debug, "%20.12f\t%20.12f\t%20.12f\t%20.12f\t%20.12f\n", dyn->data[EOB_RAD][kt], dyn->data[EOB_PPHI][kt], dyn->data[EOB_PRSTAR][kt], dt_dr_vec[kt], dphi_dr_vec[kt]);
    }
  fclose(Post_adiab_debug);
  
  

  /***********************
   * Computing integrals *
   ***********************/
  
  /** Try to integrate from from Matlab data, to see if cumint behaves properly */
  FILE *matlab;
  matlab = fopen("/mnt/c/Users/giuli/Repositories/teobresums/Matlab_Dynamics_old/Matlabdynam_q1_chi1_0.0_chi2_0.0.txt", "r");
  fscanf(matlab, "%*[^\n]\n");
  for (int v=0; v<size; v++){
      fscanf(matlab, "\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &r_mat_vec[v], &prstar_mat_vec[v], &pphi_mat_vec[v],
      &dr_dt_mat_vec[v], &omg_mat_vec[v], &t_mat_vec[v], &phi_mat_vec[v]);
      dt_dr_mat_vec[v]=1./dr_dt_mat_vec[v];
  }
  fclose(matlab);
  
  
  /** Compute time */
  cumint3(dt_dr_vec, dyn->data[EOB_RAD], size, dyn->time);

  /* Set last value for evolution */
  dyn->t = dyn->time[size-1];

  /** Compute orbital phase */
  cumint3(dphi_dr_vec, dyn->data[EOB_RAD], size, dyn->data[EOB_PHI]);
  
  
  /** Print Post-adiab dynamics on file for comparison with Matlab code */ 
  int lAL2, lBL2;
  if (usetidal)
    {lAL2 = par_get_i("LambdaAl2");
    lBL2 = par_get_i("LambdaAl2");
    }
  else
    {lAL2 = 0;
    lBL2 = 0;
    }
  char q_string[256]; //size of the number
  sprintf(q_string, "_q_%1.0f", dyn->q);
  char chi1_string[256]; 
  sprintf(chi1_string, "_chi1_%2.1f", dyn->chi1);
  char chi2_string[256]; 
  sprintf(chi2_string, "_chi2_%2.1f", dyn->chi2);
  char lAL2_string[256]; 
  sprintf(lAL2_string, "_lAL2_%3d", lAL2); 
  char lBL2_string[256]; 
  sprintf(lBL2_string, "_lBL2_%3d", lBL2); 
  char post_adiab_dyn[256]     = "Post_adiab_dynamics";
  strcat(post_adiab_dyn, q_string);
  strcat(post_adiab_dyn, chi1_string);
  strcat(post_adiab_dyn, chi2_string);
  strcat(post_adiab_dyn, lAL2_string);
  strcat(post_adiab_dyn, lBL2_string);
  strcat(post_adiab_dyn, ".dat");
  printf("I'm running with chi1 = %f and chi2=%f\n", dyn->chi1, dyn->chi2);
  
  FILE* Post_adiab_dynamics = fopen(post_adiab_dyn, "w");
  fprintf(Post_adiab_dynamics, "#8PA\tr\tp_r*\tp_phi\tdr_dt\tMOmg\tt\tphi\n");
  for (int kt = 0; kt < size; kt++)
    {
      fprintf(Post_adiab_dynamics, "\t%20.14f\t%20.14f\t%20.14f\t%20.14f\t%20.14f\t%20.14f\t%20.14f\n", dyn->data[EOB_RAD][kt], dyn->data[EOB_PRSTAR][kt], dyn->data[EOB_PPHI][kt],
      dr_dt_vec[kt], dyn->data[EOB_MOMG][kt], dyn->time[kt], dyn->data[EOB_PHI][kt]);
    }
  fclose(Post_adiab_dynamics);
    
    
  /* Free memory */
  free(A_vec);
  free(dA_vec);
  free(d2A_vec);
  free(B_vec);
  free(dB_vec);
  free(sqrtAbyB_vec);
  free(rc_vec);
  free(drc_dr_vec);
  free(d2rc_dr_vec);
  free(uc2_vec);
  free(duc_dr_vec);
  free(dAuc2_dr_vec);
  free(G_vec);
  free(dG_dr_vec);
  free(dG_dprstar_vec);
  free(dG_dprstarbyprstar_vec);
  free(dphi_dr_vec);
  free(dpphi_dr_vec);
  free(dprstar_dr_vec);
  free(H_vec);
  free(Heff_vec);
  free(Heff_orb_vec);
  free(G0_vec);
  free(dG_dr0_vec);

  printf("\n\nDAJECHEGIRO\n\n");
  
  
  return OK;
}
