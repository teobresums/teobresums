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
  const int use_tidal_gravitomagnetic = par_get_i("use_tidal_gravitomagnetic");
  int interp_uniform_grid = par_get_i("interp_uniform_grid");  
  int store_dynamics = par_get_i("output_dynamics");
  if (!(use_tidal)) interp_uniform_grid = 1; /* NQC and ringdown attachment assume uniform grids */
  if (!(use_tidal) && (use_spins)) store_dynamics = 1; /* NQC need dynamical variables */
  const int use_postadiab_dyn = STREQUAL(par_get_s("postadiabatic_dynamics"),"yes");
  if (use_postadiab_dyn) store_dynamics = 1;
  const double dt = par_get_d("dt"); 
  /** Alloc memory for dynamics and multipolar waveform */
  Dynamics *dyn;
  Waveform_lm *hlm; /* h_lm */ 
  Waveform_lm_t *hlm_t;
  Waveform_lm *hlm_nqc; /* NQC */

  const int chunk = par_get_i("size");
  int size = chunk; /* note: size can vary */
    
  if (use_postadiab_dyn) {
    size = par_get_i("postadiabatic_dynamics_size"); 
    par_set_i("size",size);
    Dynamics_alloc (&dyn, size, "dyn");
    Waveform_lm_alloc (&hlm, size, "hlm"); 
  } else {
    Dynamics_alloc (&dyn, size, "dyn"); 
    Waveform_lm_alloc (&hlm, size, "hlm"); 
  }
  Waveform_lm_t_alloc (&hlm_t);

  /* Set quick-access parameters dyn (be careful here) */
  Dynamics_set_params(dyn);
  dyn->store = dyn->noflx = 0; /* Default: do not store vars, flux on */
  
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

  /** Compute initial radius */
  const double f0 = par_get_d("initial_frequency")/time_unit_fact;
  const double r0 = eob_dyn_r0_Kepler(f0);
  //const double r0 = eob_dyn_r0_eob(f0, dyn);

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

    /*
     * Post-adiabatic dynamics
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
      goto END_ODE_EVOLUTION; 
    }
    
    /** Prepare for evolution */
    /* start counting from here */
    iter = size-1; 
    
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

    /*
     * Initial conditions for the evolution
     */
    
    /** Compute the initial conditions */
    if (use_spins) eob_dyn_ic_s(r0, dyn, dyn->y0);
    else           eob_dyn_ic(r0, dyn, dyn->y0);
    
    /** Se arrays with initial conditions */
    dyn->t       = 0.;
    dyn->r       = dyn->y0[EOB_ID_RAD];
    dyn->phi     = 0.;
    dyn->pphi    = dyn->y0[EOB_ID_PPHI];
    dyn->Omg     = dyn->y0[EOB_ID_OMGJ];	//CHECKME
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
      hlm->phase[k][0] = hlm_t->phase[k]; 
    }
    
  }
  
  if (VERBOSE) {
    /* Print initial conditions */
    PRSECTN("Initial conditions");
    for (int i = 0; i < EOB_ID_NVARS; i++)
      PRFORMd(eob_id_var[i], dyn->y0[i]);
  }
   
  /** Initialize ODE system solver */
  dyn->dt     = dt;
  dyn->t_stop = par_get_d("ode_tmax") * time_unit_fact;
  par_set_d("dt",       dyn->dt);
  par_set_d("ode_tmax", dyn->t_stop);
  dyn->ode_stop          = false;
  dyn->ode_stop_MOmgpeak = false;
  dyn->ode_stop_radius   = false;
  const double rstop = par_get_d("ode_stop_at_radius");
  const int nstep_stop = par_get_i("ode_stop_afterNdt");
  if (rstop>0.) {
    dyn->ode_stop_radius   = true;
  }

  int j;
  for (j=0; j<ODE_TSTEP_NOPT; j++) {
    if (STREQUAL(par_get_s("ode_timestep"),ode_tstep_opt[j])) {
      if (VERBOSE) printf("%-40s = %s\n","ode_timestep",ode_tstep_opt[j]);
      break;
    }
  }
  if (j==ODE_TSTEP_NOPT) {
    if (VERBOSE) printf("ode_timestep '%s' undefined, set to default\n",par_get_s("ode_timestep"));
    j = ODE_TSTEP_ADAPTIVE;
  }
  dyn->ode_timestep  = j;
  const int ode_tstep = dyn->ode_timestep;

  /* Adjust interpolation option */
  if (ode_tstep == ODE_TSTEP_UNIFORM)        interp_uniform_grid = 0; /* Not needed with uniform tstep */
  if ((use_postadiab_dyn) && (!(use_tidal))) interp_uniform_grid = 1; /* Always needed with post-adiab and BBH */
  
  const double ode_abstol = par_get_d("ode_abstol");
  const double ode_reltol = par_get_d("ode_reltol");

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
    
  /** Solve ODE */
  if (VERBOSE) PRSECTN("ODE Evolution");
  int STATUS = OK;
  while (!(dyn->ode_stop)) {
    if (VERBOSE) printf("iter %09d | t = %.9e h = %.9e | r = %.9e\n", iter, dyn->t, dyn->dt, dyn->r); 
    iter++;

    if (ode_tstep == ODE_TSTEP_UNIFORM) {
      /*  Uniform timestepping  */
      dyn->ti = dyn->t + dyn->dt;
      STATUS = gsl_odeiv2_driver_apply (d, &dyn->t, dyn->ti, dyn->y);
      if (STATUS != GSL_SUCCESS) {
	printf ("ODE solver failed. Error = %d\n", STATUS);
	return STATUS;
      }
    } 
    
    if (ode_tstep == ODE_TSTEP_ADAPTIVE) {
      /* Adaptive timestepping */
      STATUS = gsl_odeiv2_evolve_apply (e, c, s, &sys, &dyn->t, dyn->t_stop, &dyn->dt, dyn->y);
      if (STATUS != GSL_SUCCESS) {
	printf ("ODE solver failed. Error = %d\n", STATUS);
	return STATUS;	
      }
    }
    
    if (ode_tstep == ODE_TSTEP_ADAPTIVE_UNIFORM_AFTER_LSO) {
      /* Adaptive timestepping until LSO ... */
      if (dyn->r >dyn->rLSO ) { 
	STATUS = gsl_odeiv2_evolve_apply (e, c, s, &sys, &dyn->t, dyn->t_stop, &dyn->dt, dyn->y); 
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
      errorexit("ODE solver returned NaN radius.\n");
    }

    /** Waveform computation 
	Needs a r.h.s. evaluation for some vars (but no flux) */
    dyn->store = dyn->noflx = 1;
    p_eob_dyn_rhs(dyn->t, dyn->y, dyn->dy, dyn); 
    dyn->store = dyn->noflx = 0;
    eob_wav_hlm(dyn, hlm_t); 
    
    /** Update size and push arrays (if needed) */
    if (iter>size) {
    //  if (DEBUG)  printf("Push memory\n"); 
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
      printf("Stop: Max integration time reached.\n");
      dyn->ode_stop = true;
    }

    /** Stop integration at given radius (if rstop >= 0) */    
    if ((dyn->ode_stop_radius) && (dyn->r < rstop) ) {
      if (VERBOSE) printf("Stop: radius reached.\n");
      dyn->ode_stop = true;
    }

    /** Check when to break the computation
	find peak of omega curve and continue for nstep_stop iters */
    if (use_spins) {
      dyn->MOmg = dyn->Omg_orb;
    } else {
      dyn->MOmg = dyn->Omg;
    }
    if (dyn->ode_stop_MOmgpeak == false) {
      if (dyn->MOmg < dyn->MOmg_prev) {	  
	dyn->ode_stop_MOmgpeak = true;
        dyn->t_stop            = dyn->t + nstep_stop*dyn->dt; 
      } else {
	dyn->MOmg_prev = dyn->MOmg;
      }
    } else {
      if (dyn->t >= dyn->t_stop) {
	printf("Stop: Peak of Omega reached.\n");
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
  size = iter;
  par_set_i("size", size); 
  Waveform_lm_push (&hlm, size);
  Dynamics_push (&dyn, size);
  
 END_ODE_EVOLUTION:;
  
  if (DEBUG) {
    /* Output pre-interpolation wave and dynamics */
    if(par_get_i("output_multipoles")) {
      strcat(hlm->name,"_insplunge");
      Waveform_lm_output (hlm);
    }
    if (par_get_i("output_dynamics"))
      Dynamics_output(dyn);
  }
  
  if (interp_uniform_grid) { 

    /* 
     * Interpolate on uniform grid
     */

    /* Auxiliary vars */
    Waveform_lm *hlm_aux; 
    Dynamics *dyn_aux; 
  
    /* Build uniform grid of width dt and alloc tmp memory */
    //CHECKME: is this rounding under control ?!
    const int size_new = (int)((hlm->time[size-1] - hlm->time[0])/dt + 1); /* use dt from parfile */
    if (DEBUG) printf("iter=%d size=%d (%d)\n",iter,size,(iter==size));      
    if (VERBOSE) {
      PRSECTN("Interpolation to uniform grid");
      PRFORMi("interpolation_grid_size",size_new);
      PRFORMd("interpolation_grid_dt",dt);
      PRFORMd("interpolation_grid_t0",hlm->time[0]);
      PRFORMd("interpolation_grid_tN",hlm->time[size-1]);
    }

    /** Waveform */ 
    Waveform_lm_alloc (&hlm_aux, size, "");
    memcpy(hlm_aux, hlm, sizeof(Waveform_lm));
    strcpy(hlm->name, "hlm_insplunge_interp");
    hlm->size = size_new;
    hlm->time = malloc ( size_new * sizeof(double) );
    for (int k = 0; k < KMAX; k++) {
      hlm->ampli[k] = malloc ( size_new * sizeof(double) );
      hlm->phase[k] = malloc ( size_new * sizeof(double) );
    } 
   
    for (int i = 0; i < size_new; i++) {
      hlm->time[i] = i*dt;
    }    
    for (int k = 0; k < KMAX; k++) {
      interp_spline(hlm_aux->time, hlm_aux->ampli[k], size, hlm->time, size_new, hlm->ampli[k]);
    }
    for (int k = 0; k < KMAX; k++) {
      interp_spline(hlm_aux->time, hlm_aux->phase[k], size, hlm->time, size_new, hlm->phase[k]);
    }
    
    Waveform_lm_free (hlm_aux);
    
    if (store_dynamics) {

      /** Dynamics, 
	  Similar, but need to keep the non-array fields */

      Dynamics_alloc(&dyn_aux, size, "");
      memcpy(dyn_aux, dyn, sizeof(Dynamics));
      strcpy(dyn->name, "dyn_interp");
      dyn->dt   = dt;
      dyn->size = size_new; 
      dyn->time = malloc ( size_new * sizeof(double) );
      for (int v = 0; v < EOB_DYNAMICS_NVARS; v++) {
	dyn->data[v] = malloc ( size_new * sizeof(double) );
	memset(dyn->data[v], 0., size_new*sizeof(double));
      }      

      for (int i = 0; i < size_new; i++) {
	dyn->time[i] = hlm->time[i];
      }      
      for (int k = 0; k < EOB_DYNAMICS_NVARS; k++) {
	interp_spline(dyn_aux->time, dyn_aux->data[k], size, dyn->time, size_new, dyn->data[k]);
      }
      
      Dynamics_free (dyn_aux);
      
    }

    /** Update size */
    size = size_new;
    par_set_i("size", size); 
    
  }

  if (DEBUG) {
    /* Output post-interpolation wave and dynamics */
    if(par_get_i("output_multipoles")) {
      Waveform_lm_output (hlm);
      Waveform_lm_output_reim (hlm);
    }
    if (par_get_i("output_dynamics")) {
      Dynamics_output(dyn);
    }
  }
  
  if (!(use_tidal)) {

    /* 
     * BBH : add NQC 
     */
    
    if (STREQUAL(par_get_s("nqc_coefs_hlm"),"compute")) {
      
      if (VERBOSE) PRSECTN("NQC Calculation");
      
      Waveform_lm_alloc (&hlm_nqc, size, "hlm_nqc"); 
      
      /* Compute NQC corrections */
      eob_wav_hlmNQC_find_a1a2a3(dyn, hlm, hlm_nqc);
      strcat(hlm->name,"_nqc");      
      
      if (DEBUG) {
	if (par_get_i("output_nqc")) { 
	  Waveform_lm_output (hlm_nqc);
	}
	if (par_get_i("output_multipoles")) {
	  Waveform_lm_output (hlm);
	  /* Waveform_lm_output_reim (hlm); */
	}	
      }
    
      Waveform_lm_free (hlm_nqc);
    
    }
    
    /* 
     * BBH : add Ringdown 
     */
    
    if (VERBOSE) PRSECTN("Ringdown");

    /* Extend arrays */    
    const int size_ringdown = par_get_i("ringdown_extend_array");
    Waveform_lm_push (&hlm, (size+size_ringdown));
    if (DEBUG) {
      printf("Push memory for ringdown (%d + %d):",size,par_get_i("ringdown_extend_array"));
      printf(" tend = %e + %d * %e (%e) = %e\n",hlm->time[size-1],size_ringdown,dt,dt*size_ringdown,hlm->time[size-1]+dt*size_ringdown);
    }
    for (int i = size; i < (size+size_ringdown); i++) {
      hlm->time[i] = hlm->time[i-1] + dt;
    }
    size += size_ringdown;
    par_set_i("size", size);
    
    /** Ringdown attachment */
    eob_wav_ringdown(dyn, hlm);
    strcat(hlm->name,"_ringdown");
    
    if (DEBUG) {
      if (par_get_i("output_multipoles")) {
	Waveform_lm_output (hlm);
	Waveform_lm_output_reim (hlm);
      }
    }

  }
    
  /** Alloc memory for (h+,hx) */
  Waveform *hpc; 
  Waveform_alloc (&hpc, size, "waveform"); 
  
  /** Scale to physical units (if necessary) */
  const double distance = par_get_d("distance");
  double amplitude_prefactor = 1.;    
  if (!(par_get_i("use_geometric_units"))) {
    M *= MSUN_M;
    amplitude_prefactor = nu*M/(distance*MPC_M);
  } 
  const double psi = par_get_d("coalescence_angle"); 
  const double iota = par_get_d("inclination");

  /** Computation of (h+,hx) */
  compute_hpc(hlm, nu, M, distance, amplitude_prefactor, psi, iota, hpc);
  
  /** Output */
  Waveform_output (hpc);
    
  /** Free memory */
  Dynamics_free (dyn);
  Waveform_lm_free (hlm);
  Waveform_lm_t_free (hlm_t);
  Waveform_free (hpc);
  NQCdata_free (NQC);
  eob_free_params();

  return OK;
}
