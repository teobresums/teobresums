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

/** TEOBResumS v2.* main */
int main (int argc, char* argv[])
{   
  if (VERBOSE) PRSECTN(TEOBResumS_Info);
  
  Waveform *hpc = NULL; /* TD wvf */
  Waveform_lm *hmodes = NULL; /* modes */

  WaveformFD *hfpc = NULL; /* FD wvf */
  WaveformFD_lm *hfmodes = NULL; /* modes */

  Dynamics *dynf = NULL;

  int fc = 1; /* firstcall, set to 1 for now */
  int dc = DEFAULT_PARS_BBH; /* default_choice, set to BBH */
  
  /* Init parameters & set defaults */
  EOBParameters_alloc( &EOBPars );
  EOBParameters_defaults (dc, EOBPars);
  
  if (argv[1]!=NULL) {
    /* Deal with input parfile if necessary 
       (Kepr for backward compatibility) */
    par_db_init ();
    par_db_from_EOBPar (EOBPars);
    par_file_parse_merge (argv[1]);
    par_db_screen (VERBOSE);
    EOBParameters_set_from_db (EOBPars);
    const int output = EOBPars->output_dynamics + EOBPars->output_multipoles + EOBPars->output_hpc + EOBPars->output_nqc;    
    if (output) {
      if (system_mkdir(EOBPars->output_dir)) {
	printf("ERROR(TEOBResumS): %s\n",eob_error_msg[ERROR_MKDIR]);
	return ERROR_MKDIR;
      }
      par_db_write_file("params.txt");
    }
    par_db_free();    
    /* RG: if input parfile specifies BNS runs, change default_choice */ 
    if (EOBPars->LambdaAl2 > 1. && EOBPars->LambdaBl2 >1) dc = DEFAULT_PARS_BNS;
  }
  
  /* Set all firstcalls = 1 */
  for (int k=0; k < NFIRSTCALL; k++){ 
    EOBPars->firstcall[k] = 1;
  }

  /* set domain */
  eob_set_params(dc, fc);
  
  /* TD hpc, FD hpc, TD modes, FD modes, default_choice, firstcall */
  int status = EOBRun(&hpc, &hfpc, 
		      &hmodes, &hfmodes,&dynf,
		      dc, fc);
  if (status) printf("ERROR(TEOBResumS): %s\n",eob_error_msg[status]);

  Waveform_free (hpc);
  WaveformFD_free (hfpc);
  Waveform_lm_free (hmodes);
  WaveformFD_lm_free (hfmodes);
  Dynamics_free (dynf);
  
  EOBParameters_free (EOBPars);
  return status;
}

/**
 * Main routine for TD & FD EOB
 * assumes the EOBPars are allocated and filled.
 */    

int EOBRun(Waveform **hpc, WaveformFD **hfpc, 
	   Waveform_lm **hmodes, WaveformFD_lm **hfmodes, Dynamics **dynf,
	   int default_choice, int firstcall)
{
  int status = OK;

#ifdef _OPENMP
  openmp_init(VERBOSE);
#endif
  
  /* *****************************************
   * Init 
   * *****************************************
   */

  /** Do any output ? */
  const int output = EOBPars->output_dynamics + EOBPars->output_multipoles + EOBPars->output_hpc + EOBPars->output_nqc;    
  if (output) {
    /* Make output dir */
    if (system_mkdir(EOBPars->output_dir)) {
      status = ERROR_MKDIR;
      goto EXIT_POINT;
    }
  }

  /** Switch to mass-rescaled geometric units (if needed)*/
  double M = EOBPars->M; /* Msun */ 
  double time_unit_fact = 1;
  if (!(EOBPars->use_geometric_units)) {
    /* Input given in physical units, 
       rescale to geometric units and mass rescaled quantities
       compute r0 from the initial GW frequency in Hz */
    time_unit_fact = time_units_factor(M);
  } else {
    /* Input given in geometric units, 
       rescale to geometric units and mass rescaled quantities
       compute r0 from the initial GW frequency in geometric units and mass rescaled */
    EOBPars->M=1.;
  }
  
  /** Set useful pars/vars */
  const double q    = EOBPars->q;
  const double nu   = EOBPars->nu;
  const double chi1 = EOBPars->chi1;
  const double chi2 = EOBPars->chi2;
  const double ecc  = EOBPars->ecc;
  const double r_hyp = EOBPars->r_hyp;
  const double H_hyp = EOBPars->H_hyp;
  const double j_hyp = EOBPars->j_hyp;
  const int use_spins = EOBPars->use_spins;
  const int use_tidal = EOBPars->use_tidal;
  int store_dynamics = EOBPars->output_dynamics; 
  if (!(use_tidal)) store_dynamics = 1; /* NQC determination need dynamical variables */
  if ( ecc != 0. )  store_dynamics = 1; /* Eccentric waveform computation needs dynamical variables (sigmoid) */
  if ((ecc != 0. ) || (r_hyp != 0.)) EOBPars->postadiabatic_dynamics = 0;
  int use_postadiab_dyn = EOBPars->postadiabatic_dynamics;
  if (use_postadiab_dyn) store_dynamics = 1;
  const double dt = EOBPars->dt;
  
  /* *****************************************
   * Set Memory & do preliminary computations
   * *****************************************
   */
  
  /** Alloc memory for dynamics and multipolar waveform */
  Dynamics *dyn = NULL;
  Waveform_lm *hlm = NULL; /* h_lm */ 
  WaveformFD_lm *hflm = NULL; /* hf_lm (FD) */ 
  
  Waveform_lm_t *hlm_t = NULL;
  Waveform_lm *hlm_nqc = NULL; /* NQC */
  Waveform_lm *hlm_mrg = NULL; /* merger chunk */
  Dynamics *dyn_mrg = NULL;
  double ytmp[EOB_EVOLVE_NVARS], dytmp[EOB_EVOLVE_NVARS], ttmp; /* Additional buffer for post-Omegapeak ev */
  
  /* Set quick-access parameters dyn (be careful here) */
  Dynamics_alloc (&dyn, 0, "dyn"); 
  Dynamics_set_params(dyn); 
  dyn->store = dyn->noflx = 0; /* Default: do not store vars, flux on */
  
  /** Compute initial radius */
  const double f0 = EOBPars->initial_frequency/time_unit_fact;
  double r0;
  if (r_hyp != 0.) {
    // hyp case
    r0 = r_hyp;
  } else if (ecc != 0.) {
    // eccentric case
    if(EOBPars->ecc_ics == ECCICS_1PA)
      r0 = eob_dyn_r0_ecc(f0, dyn);
    else if (EOBPars->ecc_ics == ECCICS_0PA){
      if(ecc > 1e-4)
        r0 = eob_dyn_r0_ecc(f0, dyn);
      else
        r0 = eob_dyn_r0_Kepler(f0);  // fall back to QC radius
    }
  } else {
    // QC case
    r0 = eob_dyn_r0_Kepler(f0);
  }
  //r0 = eob_dyn_r0_eob(f0, dyn); /* TODO: Radius from EOB equations. This is what should be used. */

  /* If f_min is too high and ecc=0 fall back to a minimum acceptable initial radius */
  if(ecc == 0)
    if (r0 < TEOB_R0_THRESHOLD) r0 = TEOB_R0_THRESHOLD;

  /* Saving initial radius */
  EOBPars->r0 = r0;

  const int chunk = EOBPars->size;
  int size = chunk; /* note: size can vary */

  /* Change size when using PA dynamics */
  if (use_postadiab_dyn) {
    size = EOBPars->postadiabatic_dynamics_size; 
    double rmin = EOBPars->postadiabatic_dynamics_rmin;
    if(use_tidal) rmin = 14.; //EOBPars->postadiabatic_dynamics_rmin;//SB: use EOBPar do not hardcode
    size = floor(fabs(r0 - rmin)/POSTADIABATIC_DR) + 1;
    
    /* If initial radius is too close to PA limit then skip PA and go directly to ODE */
    if (size - 1 < POSTADIABATIC_NSTEP_MIN) {
        size = chunk;
        use_postadiab_dyn = 0;
    }
    EOBPars->size = size;
  }
  
  /* Allocating memory for dynamics and waveform */
  Dynamics_push (&dyn, size); 
  Waveform_lm_alloc (&hlm, size, "hlm"); 
  Waveform_lm_t_alloc (&hlm_t); 

  /** Set r.h.s. fun pointer */
  int (*p_eob_dyn_rhs)();
  if ((ecc != 0.) || (r_hyp != 0.)) {
    p_eob_dyn_rhs = &eob_dyn_rhs_ecc;
  } else if (use_spins) {
    p_eob_dyn_rhs = &eob_dyn_rhs_s;
  } else {
    p_eob_dyn_rhs = &eob_dyn_rhs;
  }

  
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
    if (check_status) {
      status = ERROR_ROOTFINDER;
      goto EXIT_POINT;
    }
    EOBPars->rLR_tidal = dyn->rLR_tidal;
    double LambdaAl2  = EOBPars->LambdaAl2;
    if( fabs(LambdaAl2) < TEOB_LAMBDA_TOL ) LambdaAl2 = 0.0;
    double LambdaBl2 = EOBPars->LambdaBl2;
    if( fabs(LambdaBl2) < TEOB_LAMBDA_TOL ) LambdaBl2 = 0.0;
    double q = EOBPars->q;

    /* Reset options */
    dyn->use_tidal = EOBPars->use_tidal;
    dyn->use_spins = EOBPars->use_spins;
    if (VERBOSE) PRFORMd("rLR_tidal",dyn->rLR_tidal); 
  }

  if (EOBPars->compute_LR && !(use_tidal)) {
    //TODO: LR COMPUTATION IS CORRECT ONLY FOR NOSPIN. IMPLEMENT SPIN VERSION IN eob_dyn_adiabLSO()
    ROOTFINDER(check_status, eob_dyn_adiabLR(dyn, &(dyn->rLR)));
    if (check_status) {
      status = ERROR_ROOTFINDER;
      goto EXIT_POINT;
    }
    EOBPars->rLR = dyn->rLR;
    if (VERBOSE) PRFORMd("rLR",dyn->rLR);
  }

  if (EOBPars->compute_LSO) {
    //TODO: LSO COMPUTATION IS CORRECT ONLY FOR NOSPIN. IMPLEMENT SPIN VERSION IN eob_dyn_adiabLSO()
    ROOTFINDER(check_status, eob_dyn_adiabLSO(dyn, &(dyn->rLSO)));
    if (check_status) {
      status = ERROR_ROOTFINDER;
      goto EXIT_POINT;
    }
    EOBPars->rLSO = dyn->rLSO;
    if (VERBOSE) PRFORMd("rLSO",dyn->rLSO);
  }
  
  /** Final BH */
  /* If (Mbhf, abhf) are user-specified, do not use the fits */
  if (!(dyn->use_tidal)) {

    // Mbhf
    if (EOBPars->Mbhf >  0.) dyn->Mbhf = EOBPars->Mbhf;
    else HealyBBHFitRemnant(chi1, chi2, q, &(dyn->Mbhf), &(dyn->abhf)); // TODO: mode to JF final mass?
    // abhf
    if (EOBPars->abhf > -2.) dyn->abhf = EOBPars->abhf;
    else dyn->abhf = JimenezFortezaRemnantSpin(dyn->nu, dyn->X1, dyn->X2, chi1, chi2);

    if (VERBOSE) {
      PRSECTN("Final black hole");
      PRFORMd("BH_final_mass[Healy/Input]",dyn->Mbhf); 
      PRFORMd("BH_final_spin[JimenezForteza/Input]",dyn->abhf);
    }
    EOBPars->Mbhf = dyn->Mbhf;
    EOBPars->abhf = dyn->abhf;
  }

  /* Iteration index */
  int iter = 0;  
  int pasize = 0;

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
        if((hlm->kmask[k])){
	  hlm->ampli[k][i] = hlm_t->ampli[k];
	  hlm->phase[k][i] = hlm_t->phase[k]; 
        }
      }

      if (dyn->time[size-1] > EOBPars->ode_tmax)
	EOBPars->postadiabatic_dynamics_stop = 1;
    }

    dyn->store = dyn->noflx = 0;
    
    if (EOBPars->postadiabatic_dynamics_stop) {
      if (VERBOSE) printf("Post-adiabatic dynamics: skip evolution.\n");
      /* SKIP ODE EVOLUTION */
      goto END_ODE_EVOLUTION; 
    }
    
    /** Prepare for evolution */
    /* start counting from here */
    pasize = size;
    iter = size-1; 
    dyn->dt = 0.5*(dyn->time[iter]-dyn->time[iter-1]);

    if (dyn->dt > 10.) dyn->dt = 10.;

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
    if (r_hyp != 0.) {
      // hyp case
      // SB 08/2020 workaround 
      // FIXME!!!!
      if (eob_dyn_ic_hyp(r_hyp, H_hyp, j_hyp, dyn, dyn->y0)) {
	      status = 1;
	      goto EXIT_POINT;
      }

    } else if (ecc !=0) {
      // eccentric case
      if(EOBPars->ecc_ics == ECCICS_1PA)
        eob_dyn_ic_ecc_PA(r0, dyn, dyn->y0);   // 1PA ICs
      else if (EOBPars->ecc_ics == ECCICS_0PA){
        if(ecc > 1e-4)
          eob_dyn_ic_ecc(r0, dyn, dyn->y0);    // adiabatic ICs
        else
          eob_dyn_ic_s(r0, dyn, dyn->y0);      // Quasi-circular ICs ("nospin" option is deprecated)
      } else
        errorexit("Unrecognized eccentric_ic flag.\n");

    } else if (use_spins) {
      // quasi-circular ICs with spins
      eob_dyn_ic_s(r0, dyn, dyn->y0);
    } else {
      // quasi-circular ICs without spins (deprecated)
      eob_dyn_ic(r0, dyn, dyn->y0);
    }
    
    /** Se arrays with initial conditions */
    dyn->t       = 0.;
    dyn->r       = dyn->y0[EOB_ID_RAD];
    dyn->phi     = 0.;
    dyn->pphi    = dyn->y0[EOB_ID_PPHI];
    dyn->Omg     = dyn->y0[EOB_ID_OMGJ];
    dyn->ddotr   = 0.; 
    dyn->prstar  = dyn->y0[EOB_ID_PRSTAR];
    dyn->E       = dyn->y0[EOB_ID_E0];
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
      if(hlm->kmask[k]) hlm->ampli[k][0] = hlm_t->ampli[k];
    }
    for (int k = 0; k < KMAX; k++) {
      if(hlm->kmask[k]) hlm->phase[k][0] = hlm_t->phase[k]; 
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
  dyn->t_stop = EOBPars->ode_tmax; // time_unit_fact;
  //EOBPars->ode_tmax = dyn->t_stop;
  dyn->ode_stop          = false;
  dyn->ode_stop_MOmgpeak = false;
  dyn->ode_stop_radius   = false;
  double rstop = EOBPars->ode_stop_radius;
  /* Avoiding useless computations for hyperbolic cases */
  if ((r_hyp != 0.) || (ecc != 0.)) {
    double r_hor = horizon_radius(dyn->nu)+0.2;
    rstop = MAX(rstop,r_hor);
  }
  if (rstop>0.) {
    dyn->ode_stop_radius   = true;
  }

  for (dyn->ode_timestep=0; dyn->ode_timestep<ODE_TSTEP_NOPT; dyn->ode_timestep++) {
    if ( EOBPars->ode_timestep == dyn->ode_timestep ) {
      if (VERBOSE) printf("%-40s = %s\n","ode_timestep",ode_tstep_opt[ dyn->ode_timestep ]);
      break;
    }
  }
  if (dyn->ode_timestep==ODE_TSTEP_NOPT) {
    if (VERBOSE) printf("ode_timestep '%s' undefined, set to default\n",ode_tstep_opt[EOBPars->ode_timestep]);
    dyn->ode_timestep = ODE_TSTEP_ADAPTIVE;
  }
  const int ode_tstep = dyn->ode_timestep;   
  const double ode_abstol = EOBPars->ode_abstol;
  const double ode_reltol = EOBPars->ode_reltol;
  if (ode_tstep == ODE_TSTEP_UNIFORM) dyn->dt = dt;
  EOBPars->dt = dyn->dt;  

  /** GSL routines for ODE integration
   * https://www.gnu.org/software/gsl/doc/html/ode-initval.html
   * http://www.csse.uwa.edu.au/programming/gsl-1.0/gsl-ref_24.html
   */
  
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
  int GSLSTATUS = OK;
  while (!(dyn->ode_stop)) {
    if (VERBOSE) printf("iter %09d | t = %.9e h = %.9e | r = %.9e\n", iter, dyn->t, dyn->dt, dyn->r); 
    iter++;

    if (ode_tstep == ODE_TSTEP_UNIFORM) {
      /* Uniform timestepping  */
      dyn->ti = dyn->t + dyn->dt;
      GSLSTATUS = gsl_odeiv2_driver_apply (d, &dyn->t, dyn->ti, dyn->y);
    } 
    
    if (ode_tstep == ODE_TSTEP_ADAPTIVE) {
      /* Adaptive timestepping */
      if ( dyn->ode_stop_MOmgpeak == true ) {
	/* if we are after the peak, slow down and fix the last steps ! */
	//GSLSTATUS = gsl_odeiv2_evolve_apply_fixed_step (e, c, s, &sys, &dyn->t, dyn->dt, dyn->y);
	dyn->ti = dyn->t + dyn->dt;
	GSLSTATUS = gsl_odeiv2_driver_apply (d, &dyn->t, dyn->ti, dyn->y);
      } else
	GSLSTATUS = gsl_odeiv2_evolve_apply (e, c, s, &sys, &dyn->t, dyn->t_stop, &dyn->dt, dyn->y);
    }
    
    if (ode_tstep == ODE_TSTEP_ADAPTIVE_UNIFORM_AFTER_LSO) {
      /* Adaptive timestepping until LSO ... */
      if (dyn->r > dyn->rLSO) { 
	GSLSTATUS = gsl_odeiv2_evolve_apply (e, c, s, &sys, &dyn->t, dyn->t_stop, &dyn->dt, dyn->y);
      } else {
	/* ... uniform afterwards */
	dyn->dt = dt_tuned_mrg;
	dyn->ti = dyn->t + dyn->dt;
	GSLSTATUS = gsl_odeiv2_driver_apply (d, &dyn->t, dyn->ti, dyn->y);
      }
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

    /** Check for failures ... */
    if (use_spins) {
      dyn->MOmg = dyn->Omg_orb;
    } else {
      dyn->MOmg = dyn->Omg;
    }

    if ((ecc != 0.) || (r_hyp != 0.)) {
      if ( (GSLSTATUS != GSL_SUCCESS) || (!isfinite(dyn->y[EOB_EVOLVE_RAD])) ) {
	if (VERBOSE) printf("Stop: Orbit crossed event horizon.\n");
	iter--; /* do count this iter! */
	dyn->ode_stop = true;
	break; /* (while) stop */
      }
    }
	
    if (dyn->ode_stop_MOmgpeak == true) {
      /* ... if after the Omega_orb peak, stop integration */
      if ( (GSLSTATUS != GSL_SUCCESS) || (!isfinite(dyn->y[EOB_EVOLVE_RAD])) ) {
	if (VERBOSE) printf("Stop: Peak of Omega reached; 2M not reached.\n");
	iter--; /* do count this iter! */
	dyn->ode_stop = true;
	break; /* (while) stop */
      }

      if (dyn->MOmg > dyn->MOmg_prev) {
	if (VERBOSE) printf("Stop: Peak of Omega reached; 2M not reached.\n");
	iter--; /* do count this iter! */
	dyn->ode_stop = true;
	break; /* (while) stop */
      } else {
	/* Mininum not reached, update the max */
	dyn->MOmg_prev = dyn->MOmg;
    } 
    }

    /* ... if before the Omega_orb peak, this is an actual error */
    if (GSLSTATUS != GSL_SUCCESS) {
      printf("GSL Error = %d", GSLSTATUS);
      /* errorexit("ODE solver returned error.\n"); */
      status = ERROR_ODEINT;
      goto EXIT_POINT;
    }
    
    /** Checking whether the dynamics produces NaN values
	this can happen if radius r becomes too small */
    if (!(isfinite(dyn->r))) {
      printf("%.1f\t%.3f\t%.3f\n", q, chi1, chi2);	
      /* errorexit("ODE solver returned NaN radius.\n"); */
      printf("ODE solver returned NaN radius.\n");
      status = ERROR_ODEINT;
      goto EXIT_POINT;
    }
    
    /** Update size and push arrays (if needed) */
    if (iter==size) {
      /* if (DEBUG)  printf("Push memory\n"); */ 
      size += chunk;
      EOBPars->size = size;
      Waveform_lm_push (&hlm, size);
      Dynamics_push (&dyn, size);
    }
    
    /** Append waveform and dynamics to arrays */
    hlm->time[iter] = hlm_t->time;
    for (int k = 0; k < KMAX; k++) {
      if(hlm->kmask[k]){
        hlm->ampli[k][iter] = hlm_t->ampli[k];
        hlm->phase[k][iter] = hlm_t->phase[k]; 
      }
    }
      
    if (store_dynamics) {
      dyn->time[iter]              = dyn->t; 
      dyn->data[EOB_RAD][iter]     = dyn->r;
      dyn->data[EOB_PHI][iter]     = dyn->phi;
      dyn->data[EOB_PPHI][iter]    = dyn->pphi;
      dyn->data[EOB_MOMG][iter]    = dyn->Omg;
      dyn->data[EOB_DDOTR][iter]   = dyn->ddotr;
      dyn->data[EOB_PRSTAR][iter]  = dyn->prstar;
      dyn->data[EOB_OMGORB][iter]  = dyn->Omg_orb;
      dyn->data[EOB_E0][iter] 	   = dyn->E;
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
    if ((dyn->ode_stop_MOmgpeak == false) && (ecc == 0.) && (r_hyp == 0.)) {
      /* Before the Omega_orb peak */
      if (dyn->MOmg < dyn->MOmg_prev) {
	/* This is the first step after the peak
	   Set things for uniform tstep evolution */
	dyn->tMOmgpeak = dyn->t; // = dyn->t-0.5*dyn->dt;
	dyn->ode_stop_MOmgpeak = true;
	dyn->dt = MIN(dyn->dt, dt_tuned_mrg); 
	dyn->t_stop = dyn->t + 2.;
	
	if (EOBPars->use_flm == USEFLM_HM) {
	  dyn->dt     = 0.1;
	  dyn->t_stop = dyn->t + 10.;
	}
	
	if (VERBOSE) printf("Peak of Omega reached, doing extra steps with h = %e\n",dyn->dt);
      } else {
	// Peak not reached, update the max 
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
  EOBPars->size = size; 
  Waveform_lm_push (&hlm, size);
  Dynamics_push (&dyn, size);
      
  /** Find peak of Omega */
  int index_pk = dyn->size-1;
  double Omega_pk = dyn->data[EOB_OMGORB][index_pk];
  for (int j = dyn->size-2; j-- ; ) {
    if (dyn->data[EOB_OMGORB][j] < Omega_pk) 
      break;
    index_pk = j;
    Omega_pk = dyn->data[EOB_OMGORB][j]; 
  }
  double *t_ptr     = &dyn->time[index_pk-2];
  double *Omega_ptr = &dyn->data[EOB_OMGORB][index_pk-2];
  double tOmg_pk    = find_max_grid(t_ptr, Omega_ptr);

  dyn->tOmg_pk = tOmg_pk;

  /** Over-writing waveform in the eccentric case - adding sigmoid */
  if ((ecc != 0) && (r_hyp == 0.)) {
    
    for (int i = 0; i < size; i++) {
      dyn->store = 1;
      dyn->t = dyn->time[i];
      dyn->y[EOB_EVOLVE_PHI] = dyn->data[EOB_PHI][i];
      dyn->y[EOB_EVOLVE_RAD] = dyn->data[EOB_RAD][i];
      dyn->y[EOB_EVOLVE_PPHI] = dyn->data[EOB_PPHI][i];
      dyn->y[EOB_EVOLVE_PRSTAR] = dyn->data[EOB_PRSTAR][i];
      eob_dyn_rhs_ecc(dyn->t, dyn->y, dyn->dy, dyn);
      
      eob_wav_hlm_ecc_sigmoid(dyn, hlm_t);

      for (int k = 0; k < KMAX; k++) {
        if((hlm->kmask[k])){
          hlm->ampli[k][i] = hlm_t->ampli[k];
          hlm->phase[k][i] = hlm_t->phase[k]; 
        }
      }
    }
  }
  
  /** Over-writing waveform in the eccentric case - computing numerical derivatives */
  /*
    if ((ecc != 0) || (r_hyp != 0.)) {
    double *r_omg     = (double*) calloc (size,sizeof(double));
    double *rdot      = (double*) calloc (size,sizeof(double));
    double *r2dot     = (double*) calloc (size,sizeof(double));
    double *r3dot     = (double*) calloc (size,sizeof(double));
    double *r4dot     = (double*) calloc (size,sizeof(double));
    double *r5dot     = (double*) calloc (size,sizeof(double));
    double *Omegadot  = (double*) calloc (size,sizeof(double));
    double *Omega2dot = (double*) calloc (size,sizeof(double));
    double *Omega3dot = (double*) calloc (size,sizeof(double));
    double *Omega4dot = (double*) calloc (size,sizeof(double));
    
    for (int i = 0; i < size; i++) {
    r_omg[i] = eob_dyn_get_romg(dyn->data[EOB_RAD][i],
    dyn->data[EOB_PRSTAR][i],
    dyn->data[EOB_PPHI][i],
    dyn);
    }
    
    D0_x_4(r_omg, dyn->time, size, rdot);
    D0_x_4(rdot,  dyn->time, size, r2dot);
    D0_x_4(r2dot, dyn->time, size, r3dot);
    D0_x_4(r3dot, dyn->time, size, r4dot);
    D0_x_4(r4dot, dyn->time, size, r5dot);
    D0_x_4(dyn->data[EOB_MOMG], dyn->time, size, Omegadot);
    D0_x_4(Omegadot,  dyn->time, size, Omega2dot);
    D0_x_4(Omega2dot, dyn->time, size, Omega3dot);
    D0_x_4(Omega3dot, dyn->time, size, Omega4dot);

    for (int i = 0; i < size; i++) {
    dyn->store = 1;
    dyn->y[EOB_EVOLVE_PHI] = dyn->data[EOB_PHI][i];
    dyn->y[EOB_EVOLVE_RAD] = dyn->data[EOB_RAD][i];
    dyn->y[EOB_EVOLVE_PPHI] = dyn->data[EOB_PPHI][i];
    dyn->y[EOB_EVOLVE_PRSTAR] = dyn->data[EOB_PRSTAR][i];
    eob_dyn_rhs_ecc(0., dyn->y, dyn->dy, dyn);
    dyn->rdot    = rdot[i];
    dyn->r2dot   = r2dot[i];
    dyn->r3dot   = r3dot[i];
    dyn->r4dot   = r4dot[i];
    dyn->r5dot   = r5dot[i];
    dyn->Omegadot  = Omegadot[i];
    dyn->Omega2dot = Omega2dot[i];
    dyn->Omega3dot = Omega3dot[i];
    dyn->Omega4dot = Omega4dot[i];
      
    eob_wav_hlm_ecc(dyn, hlm_t);

    for (int k = 0; k < KMAX; k++) {
    if((hlm->kmask[k])){
    hlm->ampli[k][i] = hlm_t->ampli[k];
    hlm->phase[k][i] = hlm_t->phase[k]; 
    }
    }
    }

    // free some memory 
    free(r_omg);
    free(rdot);
    free(r2dot);
    free(r3dot);
    free(r4dot);
    free(r5dot);
    free(Omegadot);
    free(Omega2dot);
    free(Omega3dot);
    free(Omega4dot);
    }
  */
  
 END_ODE_EVOLUTION:;
  
  /** Unwrap phase for higher modes.
      Skip PA phases: they can jump 2Pi by construction */
  if ((EOBPars->use_flm == USEFLM_HM) || (dyn->ecc != 0.)) {
    for (int k = 0; k < KMAX; k++) {
      if(hlm->kmask[k]){
	unwrap_HM(&hlm->phase[k][pasize],size-pasize);
      }
    }
  }
  
#if (DEBUG)
  /* Output wave and dynamics */
  if(EOBPars->output_multipoles) {
    strcat(hlm->name,"_insplunge");
    Waveform_lm_output (hlm);
  }
  if (EOBPars->output_dynamics)
    Dynamics_output(dyn);
#endif
  
  if (!(use_tidal) && (dyn->r < 3.)) {
    
    /* *****************************************
     * Following is for BBH : NQC & Ringdown
     * *****************************************
     */
    
    /* This is a BBH run.
       NQC and ringdown attachment currently assume uniform grids.
       Do we need to interpolate ? */
    int merger_interp = 1; /* In general, yes ... */
    if ((ode_tstep != ODE_TSTEP_ADAPTIVE) && (EOBPars->use_flm != USEFLM_HM)) merger_interp = 0; /* ... except if merger is covered by uniform tstep */

    /** NQC and ringdown attachment is done around merger 
	using auxiliary variables defined around [tmin,tmax] 
	Recall that parameters are NOT stored into these auxiliary vars */
    
    if (merger_interp) {

      /** Extract the waveform and dynamics around merger */
      double tmin = hlm->time[size-1] - 20.;         /* Use last 20M points */
      const double tmax = hlm->time[size-1] +  2*dt; /* Make sure to use or get last point */
      
      /* For hyperbolic orbits, find time when r = 10 */
      if (r_hyp != 0.) {
	      for (int j = dyn->size-1; j-- ; ) {
	        tmin = dyn->time[j];
	        if (dyn->data[EOB_RAD][j] > 10.)
	          break;
  	    }
      }
      
      /* The following routines alloc memory for the *_mrg ptrs */
      Waveform_lm_extract (hlm, tmin, tmax, &hlm_mrg, "hlm_mrg");
      Dynamics_extract (dyn, tmin, tmax, &dyn_mrg, "dyn_mrg");
      
      /** Build uniform grid of width dt and alloc tmp memory */
      double dt_merger_interp;
      if (EOBPars->use_flm == USEFLM_HM) {
	dt_merger_interp = 0.5;
      } else {
	dt_merger_interp = MIN(EOBPars->dt_merger_interp, dyn->dt);
      }
      double tstart_mrg = tOmg_pk - 8.;
      while (tstart_mrg < hlm_mrg->time[0]) /** Make sure it does not exrapolate */
	tstart_mrg += dt_merger_interp;
      const int size_mrg = get_uniform_size(hlm_mrg->time[hlm_mrg->size-1], tstart_mrg, dt_merger_interp);
      
      if (VERBOSE) {
	PRSECTN("Interpolation of merger to uniform grid");
	PRFORMi("interpolation_grid_size",size_mrg);
	PRFORMd("interpolation_grid_dt",dt_merger_interp);
	PRFORMd("interpolation_grid_t0",hlm_mrg->time[0]);
	PRFORMd("interpolation_grid_tN",hlm_mrg->time[hlm_mrg->size-1]);
      }

#if (DEBUG) 
      // Output post-interpolation wave and dynamics 
      if(EOBPars->output_multipoles) {
	Waveform_lm_output (hlm_mrg);
	Waveform_lm_output_reim (hlm_mrg);
      }
      if (EOBPars->output_dynamics) 
	Dynamics_output(dyn_mrg);
#endif
      
      /* Interp Waveform */ 
      Waveform_lm_interp (hlm_mrg, size_mrg, tstart_mrg, dt_merger_interp, "hlm_mrg_interp");
      
      /* Interp Dynamics */
      Dynamics_interp (dyn_mrg, size_mrg, tstart_mrg, dt_merger_interp, "dyn_mrg_interp");	       
      
    } /* End of merger interp */
    
    
    if (EOBPars->nqc_coefs_hlm == NQC_HLM_COMPUTE) {
      
      /** BBH : compute and add NQC */
      
      if (VERBOSE) PRSECTN("NQC Calculation");
      
      if (merger_interp) { 
	
        /* Compute NQC only around merger, 
	   add to both merger and full waveform */
        Waveform_lm_alloc (&hlm_nqc, hlm_mrg->size, "hlm_nqc"); 
        /* eob_wav_hlmNQC_find_a1a2a3_mrg_22(dyn_mrg, hlm_mrg, hlm_nqc, dyn, hlm); */
        eob_wav_hlmNQC_find_a1a2a3_mrg(dyn_mrg, hlm_mrg, hlm_nqc, dyn, hlm);
        strcat(hlm_mrg->name,"_nqc");
	
        /* Join merger to full waveform */
        Waveform_lm_join (hlm, hlm_mrg, hlm_mrg->time[0]);
        Dynamics_join (dyn, dyn_mrg, dyn_mrg->time[0]);
        size = hlm->size;
        EOBPars->size = size;
	
	/* Free the *_mrg buffers 
	   Note these were allocated in the *_extract() calls if 'merger_interp = 1'  */
      Waveform_lm_free(hlm_mrg);
      Dynamics_free(dyn_mrg);
	
      } else {
	
        /* Compute NQC and add them to full waveform */
        Waveform_lm_alloc (&hlm_nqc, size, "hlm_nqc"); 
        eob_wav_hlmNQC_find_a1a2a3(dyn, hlm, hlm_nqc);
	
      }
      
      strcat(hlm->name,"_nqc");      
      
#if (DEBUG) 
      if (EOBPars->output_nqc)  {
        Waveform_lm_output (hlm_nqc);
        Waveform_lm_output (hlm_mrg);
      }
      if (EOBPars->output_multipoles) 
     	  Waveform_lm_output (hlm);
#endif
      
      Waveform_lm_free (hlm_nqc);
      
    } // NQC_HLM_COMPUTE
    
    
    /** BBH : add Ringdown (can be omitted with flag) */
    if (EOBPars->compute_ringdown) {
      
      if (VERBOSE) PRSECTN("Ringdown");
      
      /* Extend arrays */    
      const int size_ringdown = EOBPars->ringdown_extend_array;    
      double dt_rngdn = dt;
      if (merger_interp)
	      dt_rngdn = EOBPars->dt_merger_interp; 
      
#if (DEBUG) 
      printf("Push memory for ringdown (%d + %d):",size,EOBPars->ringdown_extend_array);
      printf(" tend = %e + %d * %e (%e) = %e\n",hlm->time[size-1],size_ringdown,dt_rngdn,dt_rngdn*size_ringdown,hlm->time[size-1]+dt_rngdn*size_ringdown);
#endif
      
      Waveform_lm_push (&hlm, (size+size_ringdown));
      for (int i = size; i < (size+size_ringdown); i++) 
	hlm->time[i] = hlm->time[i-1] + dt_rngdn;
      size += size_ringdown;
      EOBPars->size = size;
      
      /* Ringdown attachment */
      eob_wav_ringdown(dyn, hlm);
      
#if (DEBUG) 
      // Output wave and dynamics 
      if(EOBPars->output_multipoles) {
	strcat(hlm->name,"_ringdown");
	Waveform_lm_output (hlm);
	Waveform_lm_output_reim (hlm);
      }
#endif
    }
    
  } /* End of BBH section */
  


  /* *****************************************
   * Compute h+, hx 
   * *****************************************
   */
  
  /** Set tc */
  /** FIXME: tc (=merger time) is correct for BNS, not for BBH */
  EOBPars->tc = hlm->time[size-1];


  /** Scale to physical units (if necessary) */
  const double distance = EOBPars->distance;
  double amplitude_prefactor = 1.;   
  if (!(EOBPars->use_geometric_units)) {
    if(EOBPars->domain==DOMAIN_TD) {
      amplitude_prefactor = nu*M*MSUN_M/(distance*MPC_M); 
    } else {
      amplitude_prefactor = nu*M*M*MSUN_M*MSUN_S/(distance*MPC_M);
    }
    M *= MSUN_S;
    EOBPars->tc*= M;   
  } 
  const double phi = Pi/2.0 - EOBPars->coalescence_angle; 
  const double iota = EOBPars->inclination;

  /** Computation of (h+,hx) */
  
  if (EOBPars->domain == DOMAIN_TD) {

    /** TIME DOMAIN */
    
    if (EOBPars->interp_uniform_grid) {
      /* Interp to uniform grid the multipoles before hpc computation */
      const double dt_interp = EOBPars->dt_interp;
      const long int size_interp = get_uniform_size(hlm->time[size-1], hlm->time[0], dt_interp);
      Waveform_lm_interp (hlm, size_interp, hlm->time[0], dt_interp, "hlm_interp");  
      size = size_interp;
    }

    /** Alloc memory for (h+,hx) */
    Waveform_alloc (hpc, size, "waveform");   
  
    /* h+, hx */  
    compute_hpc(hlm, nu, M, distance, amplitude_prefactor, phi, iota, *hpc);
         
  } else {
    
    /** FREQUENCY DOMAIN */
    
    /** Frequency domain multipolar waveform */
    WaveformFD_lm_alloc (&hflm, size, "hflm");

    /** Calculate the SPA for the multipolar waveform */
    SPA(hlm, hflm);
    
    /* The SPA performs an interpolation, 
       need to update the size */
    size = hflm->size;
    
    /** Alloc memory for (h+,hx) */
    WaveformFD_alloc (hfpc, size, "waveform");   
        
    /* h+, hx */  
    compute_hpc_FD(hflm, nu, M, distance, amplitude_prefactor, phi, iota, *hfpc);

    /* time-shift the FD waveforms */
    time_shift_FD(*hfpc, EOBPars->tc);
    
  }

  /* *****************************************
   * Output
   * *****************************************
   */

  if (output) {
    
    if (EOBPars->output_hpc)
      if (EOBPars->domain == DOMAIN_TD) Waveform_output (*hpc); 
      else                              WaveformFD_output (*hfpc);

    //FIXME: To be reinstated
    //if (EOBPars->output_multipoles) Waveform_lm_output (hlm);

    //if (EOBPars->output_multipoles) Waveform_lm_output_reim (hlm);

    //if (EOBPars->output_multipoles_fd) WaveformFD_lm_output (hflm);  //TODO: add this parameter

    if (EOBPars->output_dynamics) {
      if (EOBPars->interp_uniform_grid) {
	/* Interp to uniform grid the dynamics, rem the dyn size can be different from wf size */
	const double dt_interp_dyn = EOBPars->dt_interp;
	const int size_interp_dyn = get_uniform_size(dyn->time[dyn->size-1], dyn->time[0], dt_interp_dyn);
	Dynamics_interp (dyn, size_interp_dyn, dyn->time[0], dt_interp_dyn, "dyn_interp");  
      }
      Dynamics_output(dyn);
    }
  }

  /* *****************************************
   * Finalize 
   * *****************************************
   */
 EXIT_POINT:;

  if (status) {
    /** There was an error, return zeros waveform */
    if (EOBPars->domain == DOMAIN_TD) {
      //SB: the size here needs to be fixed to the required sampling frequency.
      //    if not, the code jumps here and size is still the one from default...
      Waveform_alloc (hpc, size, "waveform");
      Waveform_lm_alloc (hmodes, size, "hlm");
    } else  {                             
      const int interp_fd_size = get_uniform_size(EOBPars->initial_frequency, EOBPars->initial_frequency, EOBPars->df);
      WaveformFD_alloc (hfpc, interp_fd_size, "waveform_fd");
      WaveformFD_lm_alloc (hfmodes, interp_fd_size, "hlm_fd");
    }
  }
  
  *hmodes  = hlm;  /* do not free these! */
  *hfmodes = hflm; /* do not free these! */
  *dynf = dyn;     /* do not free these! */
  
#ifdef _OPENMP
  openmp_free(); 
#endif
  
  /** Free memory */
  /* Dynamics_free (dyn);       */
  /* Waveform_lm_free (hlm);    */
  /* WaveformFD_lm_free (hflm); */
  Waveform_lm_t_free (hlm_t);
  NQCdata_free (NQC);

  return OK;
}

