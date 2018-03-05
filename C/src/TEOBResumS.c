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

/**
 * GSL routines for ODE integration
 * https://www.gnu.org/software/gsl/doc/html/ode-initval.html
 * http://www.csse.uwa.edu.au/programming/gsl-1.0/gsl-ref_24.html
 */

/** Following macros can be set during compilation for special feats*/
#ifndef VERBOSE
#define VERBOSE 0 /* verbose mode is off by default */
#endif

#ifndef POSTPOSTCIRCULAR
#define POSTPOSTCIRCULAR 1 /* use post-post-circular initial conditions by default */
#endif

#ifndef EXCLUDESPINSPINTIDES
#define EXCLUDESPINSPINTIDES 0 /* use tidally deformed centr. radius with self-spin and tides by default */
#endif

/** Global vars, defined as external in header */
const int LINDEX[KMAX] = {
    2,2,
    3,3,3,
    4,4,4,4,
    5,5,5,5,5,
    6,6,6,6,6,6,
    7,7,7,7,7,7,7,
    8,8,8,8,8,8,8,8};
const int MINDEX[KMAX] = {
    1,2,
    1,2,3,
    1,2,3,4,
    1,2,3,4,5,
    1,2,3,4,5,6,
    1,2,3,4,5,6,7,
    1,2,3,4,5,6,7,8};

/** TEOBResumS main */
int main (int argc, char* argv) 
{
  char soutdir[STRLEN];

  /** Input parameters */
  if (argc == 1) {
    printf(TEOBResumS_Info);
    print_date_time();
    eob_set_params(argv, argc, VERBOSE);
  } else {
    TEOBResumS_Usage(&argv[0]);
    exit(OK);
  }

  /** Make output dir */
  strcpy(soutdir,par_get_s("output_dir"));
  system_mkdir(soutdir);
  par_db_write_file(strcat(soutdir,"/params.txt"));

  /** Switch to mass-rescaled geometric units (if needed)*/
  double M = par_get_d("M"); /* Msun */ 
  const double f0 = par_get_d("initial_frequency");
  double time_unit_fact = 1;
  double r0;
  if (!(par_get_i("use_geometric_units"))) {
    /* Input given in physical units, 
       rescale to geometric units and mass rescaled quantities
       compute r0 from the initial GW frequency in Hz */
    time_unit_fact = time_units_factor(M);
    r0 = radius0(M, f0);
  } else {
    /* Input given in geometric units, 
       rescale to geometric units and mass rescaled quantities
       compute r0 from the initial GW frequency in geometric units and mass rescaled */
    par_set_d("M", 1.);
    r0 = pow(f0*Pi, -2./3.);
  }

  //TODO: CHECK PAR RANGES AND FIX PARAMETERS

  int size = par_get_i("size"); /* note: size can vary */

  /** Alloc memory for dynamics and multipolar waveform */
  Dynamics *dyn;
  Waveform_lm *hlm; /* h_lm */ 
  Waveform_lm_t *hlm_t;
  Waveform_lm *hlm_nqc; /* NQC */

  Dynamics_alloc (&dyn, size, strcat(soutdir,"/dyn.txt"));
  Waveform_lm_alloc (&hlm, size, strcat(soutdir,"/hlm")); 
  Waveform_lm_t_alloc (&hlm_t);

  /** Set useful pars/vars */
  const double q    = par_get_d("q");
  const double nu   = par_get_d("nu");
  const double chi1 = par_get_d("chi1");
  const double chi2 = par_get_d("chi2");
  const int interp_uniform_grid = par_get_i("interp_uniform_grid");
  int check_status;
  int store_dynamics = par_get_i("output_dynamics");
  const int use_tidal = par_get_i("use_tidal");
  const int use_spins = par_get_i("use_spins");
  if (!(use_tidal)) store_dynamics = 1; 

  Dynamics_set_params(dyn);
  
  /** Compute light-ring and LSO (if needed) */
  if (use_tidal) {
    /* Compute rLR_tidal for NNLO potential and without spin part */
    dyn->use_tidal = TIDES_NNLO;
    dyn->use_spins = 0;
    ROOTFINDER(check_status, eob_dyn_adiabLR(dyn, &(dyn->rLR)));
    par_set_d("rLR_tidal", dyn->rLR_tidal);
    /* Reset options */
    dyn->use_tidal = par_get_i("use_tidal");
    dyn->use_spins = par_get_i("use_spins");
  }
  if (par_get_i("compute_LR")) {
    ROOTFINDER(check_status, eob_dyn_adiabLR(dyn, &(dyn->rLR)));
    par_set_d("rLR", dyn->rLR);
  }
  if (par_get_i("compute_LSO")) {
    //TODO: LSO COMPUTATION IS CORRECT ONLY FOR NOSPIN. IMPLEMENT SPIN VERSION IN eob_dyn_adiabLSO()
    ROOTFINDER(check_status, eob_dyn_adiabLSO(dyn, &(dyn->rLSO)));
    par_set_d("rLSO", dyn->rLSO);
  }
  
  /** Computing the initial conditions */
  int (*p_eob_dyn_rhs)();
  if (use_spins) {
    p_eob_dyn_rhs = &eob_dyn_rhs_s;
    eob_dyn_ic_s(r0, dyn, dyn->y0);
  } else {
    p_eob_dyn_rhs = &eob_dyn_rhs;
    eob_dyn_ic(r0, dyn, dyn->y0);
  }
  gsl_odeiv2_system sys = {p_eob_dyn_rhs, NULL , EOB_EVOLVE_NVARS, dyn};
  
  /** Initial conditions: t, r, phi, prstar, pphi */
  dyn->t                    = 0.;
  dyn->y[EOB_EVOLVE_RAD]    = dyn->y0[EOB_ID_RAD];
  dyn->y[EOB_EVOLVE_PHI]    = 0.;
  dyn->y[EOB_EVOLVE_PRSTAR] = dyn->y0[EOB_ID_PRSTAR];
  dyn->y[EOB_EVOLVE_PPHI]   = dyn->y0[EOB_ID_PPHI];

  /** Final BH */
  if (!(dyn->use_tidal)) {
    HealyBBHFitRemnant(chi1, chi2, q, &(dyn->Mbhf), &(dyn->abhf));
    if (PR) printf("BH_final_mass[Healy] = %e\nBH_final_spin[Healy] = %e",dyn->Mbhf,dyn->abhf);
    dyn->abhf = JimenezFortezaRemnantSpin(dyn->nu, dyn->X1, dyn->X2, chi1, chi2);
    if (PR) printf("BH_final_spin[JimenezForteza] = %e",dyn->abhf);
    par_set_d("BH_final_mass", dyn->Mbhf);
    par_set_d("BH_final_spin", dyn->abhf);
  }
    
  /** Initialize ODE system solver */
  dyn->dt            = par_get_d("dt")       * time_unit_fact;
  dyn->t1            = par_get_d("ode_t1")   * time_unit_fact;
  dyn->t_stop        = par_get_d("ode_tmax") * time_unit_fact;
  par_set_d("dt",       dyn->dt);
  par_set_d("ode_t1",   dyn->t1);
  par_set_d("ode_tmax", dyn->t_stop);
  dyn->ode_stop          = false;
  dyn->ode_stop_MOmgpeak = false;
  dyn->ode_stop_radius   = false;
  const double rstop = par_get_d("ode_stop_at_radius");
  if (rstop>=0.) {
    dyn->ode_stop_radius   = true;
  }

  int j;
  for (j=0; j<ODE_TSTEP_NOPT; j++) {
    if (STREQUAL(par_get_s("ode_timestep"),ode_tstep_opt[j])) {
      if (PR) printf("ode_timestep = %s\n",ode_tstep_opt[j]);
      break;
    }
  }
  if (j==ODE_TSTEP_NOPT) {
    if (PR) printf("ode_timestep '%s' undefined, set to default\n",par_get_s("ode_timestep"));
    j = ODE_TSTEP_ADAPTIVE;
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
    
  /** Solve ODE */
  int STATUS = OK;
  int iter = 0;
  int k;
  while (dyn->ode_stop) {
    
    if (ode_tstep == ODE_TSTEP_UNIFORM) {
      /* Uniform timestepping */
      dyn->ti = dyn->t + dyn->dt;
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
      if (dyn->y[EOB_EVOLVE_RAD]>dyn->rLSO) {
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
    dyn->r      = dyn->y[EOB_EVOLVE_RAD];
    dyn->phi    = dyn->y[EOB_EVOLVE_PHI];
    dyn->prstar = dyn->y[EOB_EVOLVE_PRSTAR];
    dyn->pphi   = dyn->y[EOB_EVOLVE_PPHI];
        
    /** Checking whether the dynamics produces NaN values
	this can happen if radius r becomes too small */
    if (!(isfinite(dyn->r))) {
	printf ("ODE solver returned NaN radius.\n");
	return ERROR;      
    }

    /** Waveform computation 
	Needs a r.h.s. evaluation */
    dyn->store = 1;
    p_eob_dyn_rhs(dyn->t, dyn->y, dyn->dy, dyn); 
    dyn->store = 0;
    eob_wav_hlm(dyn, hlm_t); 
   
    /** Update size and push arrays (if needed) */
    iter++;
    if (iter>size) {
      size = iter;
      par_set_i("size", size);
      Waveform_lm_push (&hlm, size);
      Dynamics_push (&dyn, size);
    }
    
    /** Append dynamics and waveform to vectors */
    for (k = 0; k < KMAX; k++) {
      hlm->ampli[k][iter] = hlm_t->ampli[k]; 
      hlm->phase[k][iter] = hlm_t->phase[k]; 
    }
    if (store_dynamics) {
      dyn->time[iter]             = dyn->t; //FIXME: Is this current time?
      dyn->data[EOB_RAD][iter]    = dyn->r;
      dyn->data[EOB_MOMG][iter]   = dyn->Omg;
      dyn->data[EOB_PPHI][iter]   = dyn->pphi;
      dyn->data[EOB_PRSTAR][iter] = dyn->prstar;
      dyn->data[EOB_OMGORB][iter] = dyn->Omg_orb;
      dyn->data[EOB_DDOTR][iter]  = dyn->ddotr;
    }

    /** Stop integration at given radius (if rstop >= 0) */    
    if ((dyn->ode_stop_radius) && (dyn->r < rstop) ) {
      dyn->ode_stop = true;
    }

    /** Check when to break the computation
	find peak of omega curve and continue for 4 * dt afterwards */
    if (use_spins) {
      dyn->MOmg = dyn->Omg_orb;
    } else {
      dyn->MOmg = dyn->Omg;
    }
    if (dyn->ode_stop_MOmgpeak == false) {
      if (dyn->MOmg < dyn->MOmg_prev) {	  
	dyn->ode_stop_MOmgpeak = true;
	dyn->t_stop            = dyn->t + 4.*dyn->dt;
      } else {
	dyn->MOmg_prev = dyn->MOmg;
      }
    } else {
      if (dyn->t >= dyn->t_stop) {
	dyn->ode_stop = true;
      }
    }
  
  } /* end time iteration */
  
  /** Free ODE system solver */
  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
  gsl_odeiv2_driver_free (d);

  /** Update waveform size */
  par_set_i("size", size); 

  /** Uniform grid */
  if ((!use_tidal) || (interp_uniform_grid)) {

    Waveform_lm *hlm_vecg; 
    Dynamics *dyn_vecg; 
  
    /* Build uniform grid of width dt and alloc tmp memory */
    
    if (DEBUG) printf("iter=%d size=%d (%d)\n",iter,size,(iter==size));    
    //WARNING: is this rounding under control ?!
    const int size_vecg = (int)((dyn->time[size] - dyn->time[0])/dyn->dt + 1);

    Waveform_lm_alloc (&hlm_vecg, size_vecg, "hlm_tmp");
    
    int i;    
    for (i = 0; i < size_vecg; i++) {
      hlm_vecg->time[i] = i*dyn->dt;
    }

    /* Interpolate on uniform grid */
    for (k = 0; k < KMAX; k++) {
      interp_grid(hlm->time, hlm->ampli[k], hlm->size, hlm_vecg->time, size_vecg, hlm_vecg->ampli[k]);
    }
    for (k = 0; k < KMAX; k++) {
      interp_grid(hlm->time, hlm->phase[k], hlm->size, hlm_vecg->time, size_vecg, hlm_vecg->phase[k]);
    }

    /* Swap pointers and free old memory */
    SWAPTRS(hlm_vecg, hlm);

    Waveform_lm_free (hlm_vecg);
    strcpy(hlm->name, strcat(soutdir,"/hlm.txt"));

    if (store_dynamics) {

      /* Same for dynamics */
      Dynamics_alloc(&dyn_vecg, size_vecg, "dyn_vecg");
      for (i = 0; i < size_vecg; i++) {
	dyn_vecg->time[i] = i*dyn->dt;
      }      
      for (k = 0; k < EOB_DYNAMICS_NVARS; k++) {
	interp_grid(dyn->time, dyn->data[k], size, dyn_vecg->time, size_vecg, dyn_vecg->data[k]);
      }

      /* Swap array data pointers and structure pointers
	 maintains old scalar data and parameters */
      for (k = 0; k < EOB_DYNAMICS_NVARS; k++) {
	SWAPTRS(dyn_vecg->data[k], dyn->data[k]); //CHECKME!
      }
      SWAPTRS(dyn_vecg, dyn);

      Dynamics_free (dyn_vecg);
      strcpy(dyn->name, strcat(soutdir,"/dyn.txt"));

    }
    
  }
  
  /** NQC and ringdown for BBH */
  if (!(use_tidal)) {

    Waveform_lm_alloc (&hlm_nqc, size, strcat(soutdir,"/hlm_nqc.txt"));

    /** Compute NQC corrections */
    eob_wav_hlmNQC_find_a1a2a3(size, dyn, hlm, hlm_nqc);

    /** Extend arrays */
    size += par_get_i("ringdown_extend_array");
    Waveform_lm_push (&hlm, size);
    Dynamics_push (&dyn, size );

    //TODO: store/output also ringdown wf
    //Waveform_lm_alloc (&hlm_ringdown, par_get_i("ringdown_extend_array"), "hlm_ringdown");
    
    /** Ringdown attachment */
    eob_wav_ringdown(dyn, hlm);
    //eob_wav_ringdown(dyn, hlm, hlm_ringdown);

    if (par_get_i("output_nqc")) 
      Waveform_lm_output (hlm_nqc);
    //if (par_get_i("output_ringdwon")) 
    //Waveform_lm_output (hlm_ringdown);

    Waveform_lm_free (hlm_nqc);
    //Waveform_lm_free (hlm_ringdown);

  }
    
  /** Update waveform size */
  size = hlm->size;
  par_set_i("size", size); 

  /** Alloc memory for (h+,hx) */
  Waveform *hpc; 
  Waveform_alloc (&hpc, size, "hpc");
  
  /** Scale to physical units (if necessary) */
  const double distance = par_get_d("distance");
  double amplitude_prefactor = 1.;    
  if (!(par_get_i("use_geometric_units"))) {
    M *= MSUN_M;
    amplitude_prefactor = nu*M/(distance*MPC_M);
  } 
  //const double psi = par_get_d("polarization");
  const double psi = par_get_d("coalescence angle"); 
  const double iota = par_get_d("inclination");
  
  /** Computation of (h+,hx) */
  compute_hpc(&hlm, nu, M, distance, amplitude_prefactor, psi, iota, hpc);
    
  /** Output */
  Waveform_output (hpc);
  if (par_get_i("output_multipoles"))
    Waveform_lm_output (hlm);
  if (par_get_i("output_dynamics")) 
    Dynamics_output(dyn);
  
  /** Free memory */
  Dynamics_free (dyn);
  Waveform_lm_free (hlm);
  Waveform_lm_t_free (hlm_t);
  Waveform_free (hpc);
  eob_free_params();

  return OK;
}
