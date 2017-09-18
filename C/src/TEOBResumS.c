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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#include "TEOBResumS.h"

int TEOBResumS(
	       Waveform **hpp,       /** (h+,hx) return array  **/
	       Waveform_lm **hlm      /** hlm return array    **/
	       )
{

  int size = par_get_i("size"); /* can change runtime! */ 

  double q    = par_get_d("q");
  double chi1 = par_get_d("chi1");
  double chi2 = par_get_d("chi2");

  /** Alloc memory for dynamical vars */
  Dynamics *dyn;
  Dynamics_alloc(&dyn, size);
      
  /** Set parameters */
  dyn->t             = 0.;
  dyn->phi           = 0.;
  dyn->r             = 0.;
  dyn->pph           = 0.;
  dyn->prstar        = 0.;
  dyn->Omg           = 0.;
  dyn->Omg_orb       = 0.;
  dyn->ddotr         = 0.;
  dyn->A             = 0.;
  dyn->H             = 0.;
  dyn->Heff          = 0.;
  dyn->jhat          = 0.;
  dyn->r_omega       = 0.;
  dyn->MOmg_prev     = 0.;

  // ...
  const int usespins = par_get_i("use_spins");

  /* Compute light-ring (if needed) */
  if (par_get_i("use_tidal")) {
    dyn->rLR = AdiabLR(dyn);
    dyn->rLSO = 6.0; //FIXME
    par_set_d("rLR", dyn->rLR);
  }
    
  /** Computing the initial conditions */
  gsl_odeiv2_system sys = {rhs, NULL , EOB_EVOLVE_VARS, dyn};
  if (usespins) {
    sys = {s_RHS, NULL, EOB_EVOLVE_VARS, dyn};
    s_initial(dyn->y0, dyn);
  } else {
    sys     = {rhs, NULL, EOB_EVOLVE_VARS, dyn};
    initial(dyn->y0, dyn);
  }
    
  /** Initial conditions: t, r, phi, prstar, pphi */
  dyn->t = 0.0;
  dyn->y[EOB_EVOLVE_RAD]    = dyn->y0[EOB_ID_RAD];
  dyn->y[EOB_EVOLVE_PHI]    = 0.;
  dyn->y[EOB_EVOLVE_PRSTAR] = dyn->y0[EOB_ID_PRSTAR];
  dyn->y[EOB_EVOLVE_PPH]    = dyn->y0[EOB_ID_PPH];

  /** Final BH */
  if (!(par_get_i("use_tidal"))) {
    double BH_final_mass = HealyBBHFitRemnant(chi1, chi2, q);
    par_set_d("BH_final_mass", BH_final_mass);
  }
    
  /** Initialize ODE system solver */
  dyn->dt            = par_get_d("dt");
  dyn->t1            = par_get_d("ode_tmax");
  dyn->t_stop        = par_get_d("ode_tmax");
  dyn->ode_stop          = false;
  dyn->ode_stop_MOmgpeak = false;
  in j;
  for (j=0; j<ODE_TSTEP_NOPT; j++) {
    if (STREQUAL(par_get_s("ode_timestep"),ode_tstep_opt[j])) {
      if (DEBUG) printf("ode_timestep = %s\n",ode_tstep_opt[j]);
      break;
    }
  }
  if (j==ODE_TSTEP_NOPT) {
    if (DEBUG) printf("ode_timestep '%s' undefined, set to default\n",par_get_s("ode_timestep"));
    j = 0;
  }
  dyn->ode_timestep  = j;

  const int ode_tstep = dyn->ode_timestep;
  const double ode_abstol = par_get_d("ode_abstol");
  const double ode_reltol = par_get_d("ode_relstol");

  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_driver * d          = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, dyn->dt, ode_abstol, ode_reltol);    
  gsl_odeiv2_step * s            = gsl_odeiv2_step_alloc (T, EOB_EVOLVE_NVARS);
  gsl_odeiv2_control * c         = gsl_odeiv2_control_y_new (ode_abstol, ode_reltol);
  gsl_odeiv2_evolve * e          = gsl_odeiv2_evolve_alloc (EOB_EVOLVE_NVARS);
    
  /** Sove ODE */
  const int solver_scheme = par_get_i("solver_scheme");
  int STATUS = OK;
  int iter = 0;
  int k;
  while (dyn->ode_stop) {
    
    if (ode_tstep == ODE_TSTEP_UNIFORM) {
      /* Uniform timestepping */
      dyn->ti = dyn->t + dyn->dy;
      STATUS = gsl_odeiv2_driver_apply (d, &dyn->t, dyn->ti, dyn->y);
      if (STATUS != GSL_SUCCESS) {
	printf ("ODE solver failed. Error = %d\n", STATUS);
	return STATUS;
      }
    } 
    
    if (ode_tstep == ODE_TSTEP_ADAPTIVE) {
      /* Adaptive timestepping */
      STATUS = gsl_odeiv2_evolve_apply (e, c, s, &sys, &dyn->t, dyn->t1, &dyn->dt, dyn->y);
      if (STATUS != GSL_SUCCESS) {
	printf ("ODE solver failed. Error = %d\n", STATUS);
	return STATUS;	
      }
    }
    
    if (ode_tstep == ODE_TSTEP_ADAPTIVE_UNIFORM_AFTER_LSO) {
      /* Adaptive timestepping until LSO ... */
      if (y[EOB_EVOLVE_RAD]>dyn->rLSO) {
	STATUS = gsl_odeiv2_evolve_apply (e, c, s, &sys, &dyn->t, dyn->t1, &dyn->dt, dyn->y);
	if (STATUS != GSL_SUCCESS) {
	  printf ("ODE solver failed. Error = %d\n", STATUS);
	  return STATUS;
	}
      } else {
	/* ... uniform afterwards */
	dyn->ti = dyn->t + dyn->dt;
	STATUS = gsl_odeiv2_driver_apply (d, &dyn->t, dyn->ti, dyn->y);
	if (STATUS != GSL_SUCCESS) {
	  printf ("ODE solver failed. Error = %d\n", STATUS);
	  return STATUS;
	}
      }
    }

    /** Unpack data */
    dyn->r      = y[EOB_EVOLVE_RAD];
    dyn->phi    = y[EOB_EVOLVE_PHI];
    dyn->prstar = y[EOB_EVOLVE_PRSTAR];
    dyn->pphi   = y[EOB_EVOLVE_PPHI];
        
    /** Checking whether the dynamics produces NaN values
	this can happen if radius r becomes too small */
    if (!(isfinite(dyn->r))) {
	printf ("ODE solver returned NaN radius.\n");
	return ERROR;      
    }

    /** Waveform computation 
	Need a r.h.s. evaluation */
    dyn->store = 1;
    if (usespins) {
      s_rhs(dyn->t, y, dy, dyn);
    } else {
      rhs(dyn->t, y, dy, dyn);
    }
    dyn->store = 0;
    hlm(dyn, hlm); 
   
    /** Update size and push arrays (if needed) */
    iter++;
    if (iter>size) {
      size = iter;
      par_set_i("size", size);
      Waveform_push (&hpp, size);
      Waveform_lm_push (&hlm, size);
      Dynamics_push (&dyn, size);

    }
    
    /** Append dynamics and waveform to vectors */
    dyn->time[iter] = t;
    dyn->data[EOB_RAD][iter]    = dyn->r;
    dyn->data[EOB_MOMG][iter]   = dyn->Omg;
    dyn->data[EOB_PPHI][iter]   = dyn->pphi;
    dyn->data[EOB_PRSTAR][iter] = dyn->prstar;
    dyn->data[EOB_OMGORB][iter] = dyn->Omg_orb;
    dyn->data[EOB_DDOTR][iter]  = dyn->ddotr;

    //hlm_rad_vec.push_back(h_form[lm].dat[0]);
    //hlm_phase_vec.push_back(h_form[lm].dat[1]);
    for (k = 0; k < KMAX; k++) {
      hlm->real[k][iter] =0,; // amplitude
      hlm->imag[k][iter] = 0.; //phase
    }
    
    /** Check when to break the computation
	find peak of omega curve and continue for 4 * dt afterwards */
     if (usespins) {
      dyn->MOmg = dyn->Omg_orb;
    } else {
      dyn->MOmg = dyn->Omg;
    }
    
    if (dyn->ode_stop_MOmgpeak == false) {
      if (dyn->MOmg < dyn->MOmg_prev) {	  
	dyn->ode_stop_MOmgpeak = true;
	dyn->t_stop            = dyn->t + 4.*dt;
      } else {
	dyn->MOmg_prev = dyn->MOmg;
      }
    } else {
      if (t >= t_stop) {
	dyn->ode_stop = true;
      }
    }
  
  } /* end time iteration */
  
  /** Free ODE system solver */
  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
  gsl_odeiv2_driver_free (d);

  // SB stops here waiting for improved NQC & Ringdown ...............

  /** Interpolate on uniform grid (if needed) */
  // ...

  /** Compute NQC corrections */
  // ...

  /** Ringdown attachment */
  // ...  

  // HERE 
  // - assume size is updated
  
  /** Computation of (h+,hx) */
  
  /* Init to zero */
  memset(hpp->real, 0, size*sizeof(double));
  memset(hpp->imag, 0, size*sizeof(double));
  
  /* Scale to physical units (if necessary) */
  double M = par_get_d("M");
  double nu = par_get_d("nu");
  double distance = par_get_d("distance");
  double amplitude_prefactor = 1.;    
  if (!(par_get_i("use_geometric_units"))) {
    M *= MSUN_M;
    amplitude_prefactor = nu*M/(distance*MPC_M);
    } 
  
  /* Spherical harmonics projection */
  double Y_real, Y_imag;
  double psi = par_get_d("polarization"); 
  double iota = par_get_d("inclination");
  double Aki, cosPhi, sinPhi;
  for (k = 0; k < KMAX; k++ ) {
    spinsphericalharm(&Y_real, &Y_imag, -2, L[k], M[k], psi,iota);
    for (iter = 0; iter < size; i++) {
      Aki    = 0.;//hlm_ampl_g[k][i]*amplitude_prefactor;
      cosPhi = 0;//cos(hlm_phase_g[k][i]);
      sinPhi = 0;//-sin(hlm_phase_g[k][i]);
      *hpp->real[iter] += Aki*(cosPhi*Y_real - sinPhi*Y_imag);
      *hpp->imag[iter] -= Aki*(cosPhi*Y_imag + sinPhi*Y_real);
    }
  }

  /** Output dynamics */
  // ...

  /** Free memory for dynamical vars */
  Dynamics_free(dyn);
  
  return OK;
}
