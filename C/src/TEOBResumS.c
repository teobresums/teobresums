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

/**
 * GSL routines for ODE integration
 * https://www.gnu.org/software/gsl/doc/html/ode-initval.html
 * http://www.csse.uwa.edu.au/programming/gsl-1.0/gsl-ref_24.html
 */

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

/** Global var for NQC coefficient */
NQCdata *NQC;

/** TEOBResumS main */
int main (int argc, char* argv[]) 
{

  /* *****************************************
   * Init 
   * *****************************************
   */
  
  /** Input parameters */
  if (argc == 2) {
    eob_set_params(argv[1], argc);    
    if (VERBOSE) {
      printf(TEOBResumS_Info);
      print_date_time();
      PRSECTN("Parameters");
      par_db_screen (VERBOSE);
    }
  } else {
    TEOBResumS_Usage(argv[0]);
    exit(OK);
  }

  /** Make output dir */
  system_mkdir(par_get_s("output_dir"));
  par_db_write_file("params.txt");
  
#ifdef _OPENMP
  openmp_init();
#endif
  
  /** Switch to mass-rescaled geometric units (if needed)*/
  double M = par_get_d("M"); /* Msun */ 
  double time_unit_fact = 1;
  if (!(par_get_i("use_geometric_units"))) {
    /* Input given in physical units, 
       rescale to geometric units and mass rescaled quantities
       compute r0 from the initial GW frequency in Hz */
    time_unit_fact = time_units_factor(M);
  } else {
    /* Input given in geometric units, 
       rescale to geometric units and mass rescaled quantities
       compute r0 from the initial GW frequency in geometric units and mass rescaled */
    par_set_d("M", 1.);
  }

  /** Set useful pars/vars */
  const double q    = par_get_d("q");
  const double nu   = par_get_d("nu");
  const double chi1 = par_get_d("chi1");
  const double chi2 = par_get_d("chi2");
  const int use_spins = par_get_i("use_spins");
  const int use_tidal = par_get_i("use_tidal");
  int store_dynamics = par_get_i("output_dynamics"); 
  if (!(use_tidal)) store_dynamics = 1; /* NQC determination need dynamical variables */
  int use_postadiab_dyn = STREQUAL(par_get_s("postadiabatic_dynamics"),"yes");
  if (use_postadiab_dyn) store_dynamics = 1;
  const double dt = par_get_d("dt");

  /* *****************************************
   * Set Memory & do preliminary computations
   * *****************************************
   */
  
  /** Alloc memory for dynamics and multipolar waveform */
  Dynamics *dyn;
  Waveform_lm *hlm; /* h_lm */ 
  Waveform_lm_t *hlm_t;
  Waveform_lm *hlm_nqc; /* NQC */
  Waveform_lm *hlm_mrg; /* merger chunk */
  Dynamics *dyn_mrg;
  double ytmp[EOB_EVOLVE_NVARS], dytmp[EOB_EVOLVE_NVARS], ttmp; /* Additional buffer for post-Omegapeak ev */
  
  /* Set quick-access parameters dyn (be careful here) */
  Dynamics_alloc (&dyn, 0, "dyn"); 
  Dynamics_set_params(dyn);
  dyn->store = dyn->noflx = 0; /* Default: do not store vars, flux on */
  
  /** Compute initial radius */
  const double f0 = par_get_d("initial_frequency")/time_unit_fact;
  double r0 = eob_dyn_r0_Kepler(f0);
  //double r0 = eob_dyn_r0_eob(f0, dyn); /* TODO: Radius from EOB equations. This is what should be used. */

  /* If f_min is too high fall back to a minimum acceptable initial radius */
  if (r0 < TEOB_R0_THRESHOLD) r0 = TEOB_R0_THRESHOLD;

  /* Saving initial radius */
  par_set_d("r0", r0);

  const int chunk = par_get_i("size");
  int size = chunk; /* note: size can vary */
  
  /* Change size when using PA dynamics */
  if (use_postadiab_dyn) {
    size = par_get_i("postadiabatic_dynamics_size"); 
    double rmin = par_get_d("postadiabatic_dynamics_rmin");
    if(use_tidal) rmin = par_get_d("postadiabatic_dynamics_rmin_BNS");
    size = floor(fabs(r0 - rmin)/POSTADIABATIC_DR) + 1;
    
    /* If initial radius is too close to PA limit then skip PA and go directly to ODE */
    if (size - 1 < POSTADIABATIC_NSTEP_MIN) {
        size = chunk;
        use_postadiab_dyn = 0;
    }
    par_set_i("size",size);
  }

  /* Allocating memory for dynamics and waveform */
  Dynamics_push (&dyn, size);
  Waveform_lm_alloc (&hlm, size, "hlm"); 
  Waveform_lm_t_alloc (&hlm_t); 

  /** Set r.h.s. fun pointer */
  int (*p_eob_dyn_rhs)();
  if (use_spins) p_eob_dyn_rhs = &eob_dyn_rhs_s;
  else           p_eob_dyn_rhs = &eob_dyn_rhs;
  
  /** NQC data */  
  NQCdata_alloc (&NQC);
  eob_nqc_setcoefs(NQC);
  
  /** Compute light-ring and LSO (if needed) */
  int check_status;
  if (use_tidal) {

    /* Compute rLR_tidal for NNLO potential and without spin part */
    dyn->use_tidal = TIDES_NNLO; 
    dyn->use_spins = 0;
    ROOTFINDER(check_status, eob_dyn_adiabLR(dyn, &(dyn->rLR_tidal)));
    par_set_d("rLR_tidal", dyn->rLR_tidal);
    double LambdaAl2  = par_get_d("LambdaAl2");
    if( fabs(LambdaAl2) < TEOB_LAMBDA_TOL ) LambdaAl2 = 0.0;
    double LambdaBl2 = par_get_d("LambdaBl2");
    if( fabs(LambdaBl2) < TEOB_LAMBDA_TOL ) LambdaBl2 = 0.0;
    double q = par_get_d("q");

    /* Reset options */
    dyn->use_tidal = par_get_i("use_tidal");
    dyn->use_spins = par_get_i("use_spins");
    if (VERBOSE) PRFORMd("rLR_tidal",dyn->rLR_tidal); 
  }
  if (par_get_i("compute_LR")) {
    //TODO: LR COMPUTATION IS CORRECT ONLY FOR NOSPIN. IMPLEMENT SPIN VERSION IN eob_dyn_adiabLSO()
    ROOTFINDER(check_status, eob_dyn_adiabLR(dyn, &(dyn->rLR)));
    par_set_d("rLR", dyn->rLR);
    if (VERBOSE) PRFORMd("rLR",dyn->rLR);
  }
  if (par_get_i("compute_LSO")) {
    //TODO: LSO COMPUTATION IS CORRECT ONLY FOR NOSPIN. IMPLEMENT SPIN VERSION IN eob_dyn_adiabLSO()
    ROOTFINDER(check_status, eob_dyn_adiabLSO(dyn, &(dyn->rLSO)));
    par_set_d("rLSO", dyn->rLSO);
    if (VERBOSE) PRFORMd("rLSO",dyn->rLSO);
  }   

  /** Final BH */
  if (!(dyn->use_tidal)) {
    HealyBBHFitRemnant(chi1, chi2, q, &(dyn->Mbhf), &(dyn->abhf));
    dyn->abhf = JimenezFortezaRemnantSpin(dyn->nu, dyn->X1, dyn->X2, chi1, chi2);
    if (VERBOSE) {
      PRSECTN("Final black hole");
      PRFORMd("BH_final_mass[Healy]",dyn->Mbhf); 
      PRFORMd("BH_final_spin[Healy]",dyn->abhf);
      PRFORMd("BH_final_spin[JimenezForteza]",dyn->abhf);
    }
    par_set_d("BH_final_mass", dyn->Mbhf);
    par_set_d("BH_final_spin", dyn->abhf);
  }

  /* Iteration index */
  int iter = 0;  
  
  if (use_postadiab_dyn) {

    /* *****************************************
     * Post-adiabatic dynamics
     * *****************************************
     */
    
    if (VERBOSE) PRSECTN("Post-adiabatic dynamics");

    /** Calculate dynamics */
    eob_dyn_Npostadiabatic(dyn, r0); 
 
    /** Calculate waveform */
    for (int i = 0; i < size; i++) 
      hlm->time[i] = dyn->time[i];
    
    dyn->store = dyn->noflx = 1;

    for (int i = 0; i < size; i++) {
      dyn->y[EOB_EVOLVE_RAD]    = dyn->data[EOB_RAD][i];
      dyn->y[EOB_EVOLVE_PHI]    = dyn->data[EOB_PHI][i];
      dyn->y[EOB_EVOLVE_PRSTAR] = dyn->data[EOB_PRSTAR][i]; 
      dyn->y[EOB_EVOLVE_PPHI]   = dyn->data[EOB_PPHI][i];
      p_eob_dyn_rhs(dyn->t, dyn->y, dyn->dy, dyn); 
      eob_wav_hlm(dyn, hlm_t); 
      for (int k = 0; k < KMAX; k++) {
	hlm->ampli[k][i] = hlm_t->ampli[k];
	hlm->phase[k][i] = hlm_t->phase[k]; 
      }
    }

    dyn->store = dyn->noflx = 0;
    
    if (STREQUAL(par_get_s("postadiabatic_dynamics_stop"),"yes")) {
      if (VERBOSE) printf("Post-adiabatic dynamics: skip evolution.\n");
      /* SKIP ODE EVOLUTION */
      goto END_ODE_EVOLUTION; 
    }
    
    /** Prepare for evolution */
    /* start counting from here */
    iter = size-1; 
    dyn->dt = 0.5*(dyn->time[iter]-dyn->time[iter-1]);
    
    /* Set arrays with initial conditions 
       Note current time is already set in dyn->t */
    dyn->y0[EOB_ID_RAD]  = dyn->r;
    dyn->y0[EOB_ID_PHI]  = dyn->phi;
    dyn->y0[EOB_ID_PPHI] = dyn->pphi;
    dyn->y0[EOB_ID_OMGJ] = dyn->Omg;
    dyn->y0[EOB_ID_PRSTAR] = dyn->prstar;
    //dyn->Omg_orb = 0.;//FIXME ?
    dyn->y[EOB_EVOLVE_RAD]    = dyn->r;
    dyn->y[EOB_EVOLVE_PHI]    = dyn->phi;
    dyn->y[EOB_EVOLVE_PRSTAR] = dyn->prstar; 
    dyn->y[EOB_EVOLVE_PPHI]   = dyn->pphi;
    
  } else {

    /* *****************************************
     * Initial conditions for the evolution
     * *****************************************
     */

    /** Compute the initial conditions */
    if (use_spins) eob_dyn_ic_s(r0, dyn, dyn->y0);
    else           eob_dyn_ic(r0, dyn, dyn->y0);
    
    /** Se arrays with initial conditions */
    dyn->t       = 0.;
    dyn->r       = dyn->y0[EOB_ID_RAD];
    dyn->phi     = 0.;
    dyn->pphi    = dyn->y0[EOB_ID_PPHI];
    dyn->Omg     = dyn->y0[EOB_ID_OMGJ];
    dyn->ddotr   = 0.; 
    dyn->prstar  = dyn->y0[EOB_ID_PRSTAR];
    dyn->Omg_orb = 0.;//FIXME 
    dyn->y[EOB_EVOLVE_RAD]    = dyn->r;
    dyn->y[EOB_EVOLVE_PHI]    = dyn->phi;
    dyn->y[EOB_EVOLVE_PRSTAR] = dyn->prstar; 
    dyn->y[EOB_EVOLVE_PPHI]   = dyn->pphi;
    if (store_dynamics) {
      dyn->time[0]             = dyn->t; 
      dyn->data[EOB_RAD][0]    = dyn->r;
      dyn->data[EOB_PHI][0]    = dyn->phi;
      dyn->data[EOB_PPHI][0]   = dyn->pphi;
      dyn->data[EOB_MOMG][0]   = dyn->Omg;
      dyn->data[EOB_DDOTR][0]  = dyn->ddotr;
      dyn->data[EOB_PRSTAR][0] = dyn->prstar;
      dyn->data[EOB_OMGORB][0] = dyn->Omg_orb;
      dyn->data[EOB_E0][0]     = dyn->E;
    }
    
    /** Waveform computation at t = 0 
	Needs a r.h.s. evaluation for some vars (no flux) */
    dyn->store = dyn->noflx = 1;
    p_eob_dyn_rhs(dyn->t, dyn->y, dyn->dy, dyn); 
    dyn->store = dyn->noflx = 0;
    eob_wav_hlm(dyn, hlm_t); 
    
    /** Append waveform to arrays */
    hlm->time[0] = 0.;
    for (int k = 0; k < KMAX; k++) {
      hlm->ampli[k][0] = hlm_t->ampli[k];
    }
    for (int k = 0; k < KMAX; k++) {
      hlm->phase[k][0] = hlm_t->phase[k]; 
    }

    /** Prepare for evolution */
    dyn->dt = dt;
    
  }
  
  if (VERBOSE) {
    /* Print initial conditions */
    PRSECTN("Initial conditions");
    for (int i = 0; i < EOB_ID_NVARS; i++)
      PRFORMd(eob_id_var[i], dyn->y0[i]);
  }

  /* *****************************************
   * ODE Evolution
   * *****************************************
   */
     
  /** Initialize ODE system solver */
  dyn->t_stop = par_get_d("ode_tmax") * time_unit_fact;
  par_set_d("ode_tmax", dyn->t_stop);
  dyn->ode_stop          = false;
  dyn->ode_stop_MOmgpeak = false;
  dyn->ode_stop_radius   = false;
  const double rstop   = par_get_d("ode_stop_at_radius"); 
  if (rstop>0.) {
    dyn->ode_stop_radius   = true;
  }

  for (dyn->ode_timestep=0; dyn->ode_timestep<ODE_TSTEP_NOPT; dyn->ode_timestep++) {
    if (STREQUAL(par_get_s("ode_timestep"),ode_tstep_opt[ dyn->ode_timestep ])) {
      if (VERBOSE) printf("%-40s = %s\n","ode_timestep",ode_tstep_opt[ dyn->ode_timestep ]);
      break;
    }
  }
  if (dyn->ode_timestep==ODE_TSTEP_NOPT) {
    if (VERBOSE) printf("ode_timestep '%s' undefined, set to default\n",par_get_s("ode_timestep"));
    dyn->ode_timestep = ODE_TSTEP_ADAPTIVE;
  }
  const int ode_tstep = dyn->ode_timestep;   
  const double ode_abstol = par_get_d("ode_abstol");
  const double ode_reltol = par_get_d("ode_reltol");
  if (ode_tstep == ODE_TSTEP_UNIFORM) dyn->dt = dt;
  par_set_d("dt", dyn->dt);  

  /* GSL integrator memory */
  gsl_odeiv2_system sys          = {p_eob_dyn_rhs, NULL , EOB_EVOLVE_NVARS, dyn};
#if (USERK45)
  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_driver * d          = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, dyn->dt, ode_abstol, ode_reltol);    
#else
  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_driver * d          = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, dyn->dt, ode_abstol, ode_reltol);    
#endif
  gsl_odeiv2_step * s            = gsl_odeiv2_step_alloc (T, EOB_EVOLVE_NVARS);
  gsl_odeiv2_control * c         = gsl_odeiv2_control_y_new (ode_abstol, ode_reltol);
  gsl_odeiv2_evolve * e          = gsl_odeiv2_evolve_alloc (EOB_EVOLVE_NVARS);

  /* Set optimized dt around merger */
  const double dt_tuned_mrg = get_mrg_timestep(q, chi1, chi2);
  
  /** Solve ODE */
  if (VERBOSE) PRSECTN("ODE Evolution");
  int STATUS = OK;
  while (!(dyn->ode_stop)) {
    if (VERBOSE) printf("iter %09d | t = %.9e h = %.9e | r = %.9e\n", iter, dyn->t, dyn->dt, dyn->r); 
    iter++;

    if (ode_tstep == ODE_TSTEP_UNIFORM) {
      /* Uniform timestepping  */
      dyn->ti = dyn->t + dyn->dt;
      STATUS = gsl_odeiv2_driver_apply (d, &dyn->t, dyn->ti, dyn->y);
    } 
    
    if (ode_tstep == ODE_TSTEP_ADAPTIVE) {
      /* Adaptive timestepping */
      if ( dyn->ode_stop_MOmgpeak == true ) 
	/* if we are after the peak, slow down and fix the last steps ! */
	STATUS = gsl_odeiv2_evolve_apply_fixed_step (e, c, s, &sys, &dyn->t, dyn->dt, dyn->y);
      else
	STATUS = gsl_odeiv2_evolve_apply (e, c, s, &sys, &dyn->t, dyn->t_stop, &dyn->dt, dyn->y);
    }
    
    if (ode_tstep == ODE_TSTEP_ADAPTIVE_UNIFORM_AFTER_LSO) {
      /* Adaptive timestepping until LSO ... */
      if (dyn->r > dyn->rLSO) { 
	STATUS = gsl_odeiv2_evolve_apply (e, c, s, &sys, &dyn->t, dyn->t_stop, &dyn->dt, dyn->y);
      } else {
	/* ... uniform afterwards */
	dyn->dt = dt_tuned_mrg;
	dyn->ti = dyn->t + dyn->dt;
	STATUS = gsl_odeiv2_evolve_apply_fixed_step (e, c, s, &sys, &dyn->t, dyn->dt, dyn->y);
      }
    }

    /** Check for failures ... */
    if (dyn->ode_stop_MOmgpeak == true) {
      /* ... if after the Omega_orb peak, stop integration */
      if ( (STATUS != GSL_SUCCESS) || (!isfinite(dyn->y[EOB_EVOLVE_RAD])) ) {
	if (VERBOSE) printf("Stop: Peak of Omega reached; 2M not reached.\n");
	iter--; /* do count this iter! */
	dyn->ode_stop = true;
	break; /* (while) stop */
      }
    }

    /* ... if before the Omega_orb peak, this is an actual error */
    if (STATUS != GSL_SUCCESS) {
      printf("Error = %d", STATUS);
      errorexit("ODE solver returned error.\n");
    }
    
    /** Checking whether the dynamics produces NaN values
       this can happen if radius r becomes too small */
    if (!(isfinite(dyn->r))) {
      printf("%.1f\t%.3f\t%.3f\n", q, chi1, chi2);	
      errorexit("ODE solver returned NaN radius.\n");
    }

    /** Unpack data */
    dyn->r      = dyn->y[EOB_EVOLVE_RAD];
    dyn->phi    = dyn->y[EOB_EVOLVE_PHI];
    dyn->prstar = dyn->y[EOB_EVOLVE_PRSTAR];
    dyn->pphi   = dyn->y[EOB_EVOLVE_PPHI];
    
    /** Waveform computation 
	Needs a r.h.s. evaluation for some vars (but no flux) */
    dyn->store = dyn->noflx = 1;
    p_eob_dyn_rhs(dyn->t, dyn->y, dyn->dy, dyn); 
    dyn->store = dyn->noflx = 0;
    eob_wav_hlm(dyn, hlm_t); 
    
    /** Update size and push arrays (if needed) */
    if (iter==size) {
      /* if (DEBUG)  printf("Push memory\n"); */ 
      size += chunk;
      par_set_i("size", size);
      Waveform_lm_push (&hlm, size);
      Dynamics_push (&dyn, size);
    }
    
    /** Append waveform and dynamics to arrays */
    hlm->time[iter] = hlm_t->time;
    for (int k = 0; k < KMAX; k++) {
      hlm->ampli[k][iter] = hlm_t->ampli[k];
      hlm->phase[k][iter] = hlm_t->phase[k]; 
    }
      
    if (store_dynamics) {
      dyn->time[iter]             = dyn->t; 
      dyn->data[EOB_RAD][iter]    = dyn->r;
      dyn->data[EOB_PHI][iter]    = dyn->phi;
      dyn->data[EOB_PPHI][iter]   = dyn->pphi;
      dyn->data[EOB_MOMG][iter]   = dyn->Omg;
      dyn->data[EOB_DDOTR][iter]  = dyn->ddotr;
      dyn->data[EOB_PRSTAR][iter] = dyn->prstar;
      dyn->data[EOB_OMGORB][iter] = dyn->Omg_orb;
      dyn->data[EOB_E0][iter] 	  = dyn->E;
    }

    /** Stop integration if reached max time */    
    if (dyn->t > dyn->t_stop) {
      if (VERBOSE) printf("Stop: Max integration time reached.\n");
      dyn->ode_stop = true;
    }

    /** Stop integration at given radius (if rstop >= 0) */    
    if ((dyn->ode_stop_radius) && (dyn->r < rstop) ) {
      if (VERBOSE) printf("Stop: radius reached.\n");
      dyn->ode_stop = true;
    }

    /** Check when to break the computation
	find peak of omega curve and continue for 2M */
    if (use_spins) {
      dyn->MOmg = dyn->Omg_orb;
    } else {
      dyn->MOmg = dyn->Omg;
    }

    if (dyn->ode_stop_MOmgpeak == false) {
      /* Before the Omega_orb peak */      
      if (dyn->MOmg < dyn->MOmg_prev) {
	/* This is the first step after the peak
	   Set things for uniform tstep evolution */
	dyn->tMOmgpeak = dyn->t; // = dyn->t-0.5*dyn->dt;
	dyn->ode_stop_MOmgpeak = true;
	dyn->dt = MIN(dyn->dt, dt_tuned_mrg); 
	dyn->t_stop = dyn->t + 2.;
	if (VERBOSE) printf("Peak of Omega reached, doing extra steps with h = %e\n",dyn->dt);
      } else {
	/* Peak not reached, update the max */
	dyn->MOmg_prev = dyn->MOmg;
      }      
    } else {
      if (dyn->t >= dyn->t_stop) {
	if (VERBOSE) printf("Stop: Peak of Omega reached.\n");
	dyn->ode_stop = true;
      }
    }
  
  } /* end time iteration */
  
  /** Free ODE system solver */
  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
  gsl_odeiv2_driver_free (d);

  /** Update waveform and dynamics size 
      resize to actual size */
  size = iter+1;
  par_set_i("size", size); 
  Waveform_lm_push (&hlm, size);
  Dynamics_push (&dyn, size);
  
 END_ODE_EVOLUTION:;
  
#if (DEBUG) 
  /* Output wave and dynamics */
  if(par_get_i("output_multipoles")) {
    strcat(hlm->name,"_insplunge");
    Waveform_lm_output (hlm);
  }
  if (par_get_i("output_dynamics"))
    Dynamics_output(dyn);
#endif

  if (!(use_tidal)) {
    
    /* *****************************************
     * Following is for BBH : NQC & Ringdown
     * *****************************************
     */
    
    /* This is a BBH run.
       NQC and ringdown attachment currently assume uniform grids.
       Do we need to interpolate ? */
    int merger_interp = 1; /* In general, yes ... */
    if (ode_tstep != ODE_TSTEP_ADAPTIVE) merger_interp = 0; /* ... except if merger is covered by uniform tstep */

    /* Make sure merger ptrs point to something */
    hlm_mrg = hlm; 
    dyn_mrg = dyn;

    if (merger_interp) {
      
      /** NQC and ringdown attachment is done around merger 
	  using auxiliary variables defined around [tmin,tmax] 
	  Recall that parameters are NOT stored into these auxiliary vars */

      const double tmin = hlm->time[size-1] - 20; /* Use last 20M points */
      const double tmax = hlm->time[size-1] +  2*dt; /* Make sure to use or get last point */
      hlm_mrg = NULL;
      dyn_mrg = NULL;
      
      Waveform_lm_extract (hlm, tmin, tmax, &hlm_mrg, "hlm_mrg");
      Dynamics_extract (dyn, tmin, tmax, &dyn_mrg, "dyn_mrg");

      /**  Interpolate mrg on uniform grid */
      
      /* Build uniform grid of width dt and alloc tmp memory */
      double dt_merger_interp = par_get_d("dt_merger_interp");
      dt_merger_interp = MIN(dt_merger_interp, (dyn->time[size-1] - dyn->tMOmgpeak)/4 ); /* Make sure to have always 3 points */
      //dt_merger_interp = MIN(dt_merger_interp, dyn->dt);

      const int size_mrg = get_uniform_size(hlm_mrg->time[hlm_mrg->size-1], hlm_mrg->time[0], dt_merger_interp);
      if (VERBOSE) {
	PRSECTN("Interpolation of merger to uniform grid");
	PRFORMi("interpolation_grid_size",size_mrg);
	PRFORMd("interpolation_grid_dt",dt_merger_interp);
	PRFORMd("interpolation_grid_t0",hlm_mrg->time[0]);
	PRFORMd("interpolation_grid_tN",hlm_mrg->time[hlm_mrg->size-1]);
      }
      
      /* Interp Waveform */ 
      Waveform_lm_interp (hlm_mrg, size_mrg, hlm_mrg->time[0], dt_merger_interp, "hlm_mrg_interp");
      
      /* Interp Dynamics */
      Dynamics_interp (dyn_mrg, size_mrg, dyn_mrg->time[0], dt_merger_interp, "dyn_mrg_interp");

#if (DEBUG) 
      /* Output post-interpolation wave and dynamics */
      if(par_get_i("output_multipoles")) {
	Waveform_lm_output (hlm_mrg);
	Waveform_lm_output_reim (hlm_mrg);
      }
      if (par_get_i("output_dynamics")) 
	Dynamics_output(dyn_mrg);
#endif
       
    } /* End of merger interp */
    
    
    if (STREQUAL(par_get_s("nqc_coefs_hlm"),"compute")) {

      /** BBH : compute and add NQC */

      if (VERBOSE) PRSECTN("NQC Calculation");
      
      if (merger_interp) {

	/* Compute NQC only around merger, 
	   add to both merger and full waveform */
	Waveform_lm_alloc (&hlm_nqc, hlm_mrg->size, "hlm_nqc"); 
	eob_wav_hlmNQC_find_a1a2a3_mrg(dyn_mrg, hlm_mrg, hlm_nqc, dyn, hlm);
	strcat(hlm_mrg->name,"_nqc");

	/* Join merger to full waveform */
	Waveform_lm_join (hlm, hlm_mrg, hlm_mrg->time[0]);
	Dynamics_join (dyn, dyn_mrg, dyn_mrg->time[0]);
	size = hlm->size;
	par_set_i("size", size);
	
      } else {
	
	/* Compute NQC and add them to full waveform */
	Waveform_lm_alloc (&hlm_nqc, size, "hlm_nqc"); 
	eob_wav_hlmNQC_find_a1a2a3(dyn, hlm, hlm_nqc);

      }
      
      strcat(hlm->name,"_nqc");      

#if (DEBUG) 
      if (par_get_i("output_nqc"))  {
	Waveform_lm_output (hlm_nqc);
	Waveform_lm_output (hlm_mrg);
      }
      if (par_get_i("output_multipoles")) 
	Waveform_lm_output (hlm);
#endif
      
      Waveform_lm_free (hlm_nqc);
      
    }
    
    /** BBH : add Ringdown */
    
    if (VERBOSE) PRSECTN("Ringdown");
    
    /* Extend arrays */    
    const int size_ringdown = par_get_i("ringdown_extend_array");    
    double dt_rngdn = dt;
    if (merger_interp)
      dt_rngdn = par_get_d("dt_merger_interp"); 
    
#if (DEBUG) 
    printf("Push memory for ringdown (%d + %d):",size,par_get_i("ringdown_extend_array"));
    printf(" tend = %e + %d * %e (%e) = %e\n",hlm->time[size-1],size_ringdown,dt_rngdn,dt_rngdn*size_ringdown,hlm->time[size-1]+dt_rngdn*size_ringdown);
#endif

    Waveform_lm_push (&hlm, (size+size_ringdown));
    for (int i = size; i < (size+size_ringdown); i++) 
      hlm->time[i] = hlm->time[i-1] + dt_rngdn;
    size += size_ringdown;
    par_set_i("size", size);
    
    /* Ringdown attachment */
    eob_wav_ringdown(dyn, hlm);
    strcat(hlm->name,"_ringdown");
    
  } /* End of BBH section */

  /* *****************************************
   * Compute h+, hx 
   * *****************************************
   */
  
  /** Scale to physical units (if necessary) */
  const double distance = par_get_d("distance");
  double amplitude_prefactor = 1.;   
  if (!(par_get_i("use_geometric_units"))) {
    amplitude_prefactor = nu*M*MSUN_M/(distance*MPC_M);    
    M *= MSUN_S;   
  } 
  /** Consistent with LAL master conventions for input parameters */
  /** Azimuthal angle phi follows LAL convention of LIGO-T1800226 for master, 
    * where the polarization basis is defined with \f$\Omega=\pi/2\f$ and
    * \f$Z = \sin{\iota}\sin{\Phi}x + \sin{\iota}\cos{\Phi}y + \cos{\iota}z\f$
    */
  const double phi = Pi/2.0 - par_get_d("coalescence_angle"); 
  const double iota = par_get_d("inclination");

  /** Computation of (h+,hx) */

  const int interp_uniform_grid = par_get_i("interp_uniform_grid");

  if (interp_uniform_grid == INTERP_UNIFORM_GRID_HLM) {
    /* Interp to uniform grid the multipoles before hpc computation */
    const double dt_interp = par_get_d("dt_interp");
    const int size_interp = get_uniform_size(hlm->time[size-1], hlm->time[0], dt_interp); 
    Waveform_lm_interp (hlm, size_interp, hlm->time[0], dt_interp, "hlm_interp");  
    size = size_interp;  
  }

  /** Alloc memory for (h+,hx) */
  Waveform *hpc; 
  Waveform_alloc (&hpc, size, "waveform"); 

  /* h+, hx */  
  compute_hpc(hlm, nu, M, distance, amplitude_prefactor, phi, iota, hpc);
  
  if (interp_uniform_grid == INTERP_UNIFORM_GRID_HPC) {
    /* Interp to uniform grid phase and amplitude of h+, hx  */
    const double dt_interp_hpc = par_get_d("dt_interp") * M;
    const int size_interp_hpc = get_uniform_size(hpc->time[size-1], hpc->time[0], dt_interp_hpc); 
    Waveform_rmap (hpc, 1, 0); /* do not unwrap here ... */
    unwrap_proxy(hpc->phase, hlm->phase[1], hpc->size, 1); /* ... but use phi22 as unwrap proxy */
    Waveform_interp_ap (hpc, size_interp_hpc, hpc->time[0], dt_interp_hpc, "waveform_interp");
    /* Waveform_interp (hpc, size_interp_hpc, hpc->time[0], dt_interp_hpc, "waveform_interp"); */ /* this interp real/imag */
    if (par_get_i("output_multipoles")) {
      const double dt_interp_hlm = par_get_d("dt_interp");
      const int size_interp_hlm = get_uniform_size(hlm->time[size-1], hlm->time[0], dt_interp_hlm); 
      Waveform_lm_interp (hlm, size_interp_hlm, hlm->time[0], dt_interp_hlm, "hlm_interp");
    }
    
  }
  
  if ( (interp_uniform_grid) && (par_get_i("output_dynamics")) ) {
    /* Interp to uniform grid the dynamics, rem the dyn size can be different from wf size */
    const double dt_interp_dyn = par_get_d("dt_interp");
    const int size_interp_dyn = get_uniform_size(dyn->time[dyn->size-1], dyn->time[0], dt_interp_dyn);
    Dynamics_interp (dyn, size_interp_dyn, dyn->time[0], dt_interp_dyn, "dyn_interp");  
  }
        
  /** Output */
  if (par_get_i("output_hpc"))
    Waveform_output (hpc);
  if (par_get_i("output_multipoles")) {
    Waveform_lm_output (hlm); 
    Waveform_lm_output_reim (hlm);
  }
  if (par_get_i("output_dynamics"))
    Dynamics_output(dyn);

  /* *****************************************
   * Finalize 
   * *****************************************
   */
  
#ifdef _OPENMP
  openmp_free(); 
#endif
  
  /** Free memory */
  Dynamics_free (dyn);
  Waveform_lm_free (hlm);
  Waveform_lm_t_free (hlm_t);
  Waveform_free (hpc);
  NQCdata_free (NQC);
  eob_free_params();

  return OK;
}
