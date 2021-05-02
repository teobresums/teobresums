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

/**
 * @file TEOBResumPars.h
 * @brief Parameter manager
 *
 * Parameters are managed using the type EOBParameters.
 * Parameters can be input via a parfile as versions < v2.0.
 * 
 * This file contains routines to manage
 * - EOBParameters
 * - the old parameters db (libconfig wrappers, see below)
 * - copying between the two
 * - set EOB parameters for the run
 * 
 * History
 * - SB 11/2019 Added the par structure type EOBParameters and routines to work with it.
 * - v0.0 and v1.0 worked only with the parameter db and the parfile.
 *
 *  TODO
 * - SB 04/2020 we want to simplify the logic for the input file, and also have command-line inputs for the PE parameters (similar to the C++).
 *
 * OLD Documentation about the parameter db:
 *
 * Parameters are managed using a database (db) ' key = value '
 * The db is initialized by the file $TEOBRESUMS/par/default.par
 * Use libconfig API but in a simplified way (http://hyperrealm.com/libconfig/libconfig_manual.html)
 * The db entries can be set/accessed in other parts of the code with simple wrapper functions
 * Note autoconversion int/float is disabled, type must be specified in the *.par
 */

#include "TEOBResumS.h"

#define DEBUG_THIS_FILE (0) /* = 1 to compile and debug various routines in this file */
#if (DEBUG_THIS_FILE)
#undef errorexit
#undef errorexits
#define errorexit(s) {printf("%s\n",s); exit(ERROR);}
#define errorexits(s,t) {printf("%s %s\n",s,t); exit(ERROR);}
#endif


/*
 * routines to work with EOBParameters
 */

/** Global var for EOB parameters 
    (there is an extern in the header) */
EOBParameters *EOBPars;

void EOBParameters_alloc (EOBParameters **eobp)
{
  *eobp = (EOBParameters *) calloc(1, sizeof(EOBParameters));
  if (eobp == NULL)
    errorexit("Out of memory");
  /* the arrays below are allocated in EOBParameters_defaults */
  /* (*eobp)->use_mode_lm_size = 1;  */
  /* (*eobp)->use_mode_lm = malloc ( 1 * sizeof(int) ); */
  /* (*eobp)->use_mode_lm [0] = -1; */
  /* (*eobp)->output_lm_size = 1; */
  /* (*eobp)->output_lm = malloc ( 1 * sizeof(int) ); */
  /* (*eobp)->output_lm [0] = -1; */
} 

void EOBParameters_free (EOBParameters *eobp)
{
  if (!eobp) return;
  if (eobp->use_mode_lm) free (eobp->use_mode_lm);
  if (eobp->output_lm) free (eobp->output_lm);
  free(eobp);
}

/* Following default parameters should match those for production runs */
void EOBParameters_defaults (int choose, EOBParameters *eobp)
{

  eobp->domain = DOMAIN_TD;
  
  eobp->M = 1.;
  eobp->q = 1.;

  eobp->chi1 = 0.;
  eobp->chi2 = 0.;
  
  eobp->distance = 1.;
  eobp->inclination = 0.;
  eobp->coalescence_angle = 0.;
  eobp->polarization = 0.;

  eobp->r0 = 0.;
  eobp->initial_frequency = 0.004;

  eobp->df = 1.;

  eobp->LambdaAl2 = 0.; // Tidal gravitoelectric parameter Lambda for star A ell=2
  eobp->LambdaBl2 = 0.;
  eobp->LambdaAl3 = 0.; 
  eobp->LambdaAl4 = 0.; 
  eobp->LambdaBl3 = 0.; 
  eobp->LambdaBl4 = 0.; 
  eobp->SigmaAl2  = 0.; // Tidal gravitomagnetic parameter Sigma for star A ell=2
  eobp->SigmaBl2  = 0.;
  eobp->use_Yagi_fits = 0;
  eobp->pGSF_tidal =4.0;// p-power in GSF tidal potential model

  eobp->use_spins=1; // use spins ?
  eobp->project_spins=0;

  /* options */

  eobp->use_geometric_units = 1;

  eobp->use_speedytail = 0;

  eobp->dt_merger_interp = 0.5;

  eobp->interp_uniform_grid = 0; //INTERP_UNIFORM_GRID_HPC;
  eobp->dt_interp = 0.5;
  eobp->srate_interp = 4096.;

  int hlm[] = {1};      //indexes of multipoles to use
  eobp->use_mode_lm_size = 1;
  eobp->use_mode_lm = malloc (eobp->use_mode_lm_size * sizeof(int) );
  memcpy(eobp->use_mode_lm, hlm, eobp->use_mode_lm_size * sizeof(int));

  /* EOB Settings */
  
  eobp->postadiabatic_dynamics=1;
  eobp->postadiabatic_dynamics_N=8; // post-adiabatic order
  eobp->postadiabatic_dynamics_size=800; // grid size 
  eobp->postadiabatic_dynamics_rmin=14.; // minimum radius (end of PA dynamics)
  eobp->postadiabatic_dynamics_stop=0; // stop after post-adiabatic dynamics 

  eobp->centrifugal_radius=CENTRAD_NLO; // {LO, NLO, NNLO, NNLOS4, NOSPIN, NOTIDES}
  eobp->use_flm=USEFLM_HM; // "SSLO", "SSNLO", "HM"
  
  eobp->compute_LR=0; // calculate LR ?
  eobp->compute_LSO=0; // calculate LSO ?
  eobp->compute_LR_guess=3.;
  eobp->compute_LSO_guess=6.;

  eobp->nqc=NQC_AUTO; // {"no", "auto", "manual"}
  eobp->nqc_coefs_flx=NQC_FLX_NONE; // {"none", "nrfit_nospin20160209", "nrfit_spin202002", "fromfile"}
  eobp->nqc_coefs_hlm=NQC_HLM_NONE; // {"compute", "none", "nrfit_nospin20160209", "nrfit_spin202002", "fromfile"}
  strcpy(eobp->nqc_coefs_flx_file,"");
  strcpy(eobp->nqc_coefs_hlm_file,"");

  /* Output */
  
  strcpy(eobp->output_dir, "./data");  // output dir
  eobp->output_hpc= 1; // output h+,hx
  eobp->output_multipoles= 0; // output multipoles
  eobp->output_dynamics=0; // output dynamics
  eobp->output_nqc=0; // output NQC waveform
  eobp->output_nqc_coefs=0; // output multipolar NQC coefs (if determined)
  eobp->output_ringdown=0; // output ringdown waveform
  
  int klm[] = {1, 2};      //indexes of multipoles to ouput
  eobp->output_lm_size = 2;
  eobp->output_lm = malloc (eobp->output_lm_size * sizeof(int) );
  memcpy(eobp->output_lm, klm, eobp->output_lm_size * sizeof(int));
  
  /* Evolution settings */
  
  eobp->srate=4096.; // sampling rate, used if input is given in physical unit, reset based on tstep otherwise
  eobp->dt=0.5; // timestep, used if input is given in geometric unit, reset based on srate otherwise
  eobp->size=500; // size of the arrays (chunks, dynamically extended)
  eobp->ringdown_extend_array=500; // grid points to extend arrays for ringdown attachment
  eobp->ode_timestep=ODE_TSTEP_ADAPTIVE; // specify ODE solver timestep "uniform","adaptive","adaptive+uniform_after_LSO","undefined"
  eobp->ode_abstol=1e-13; // ODE solver absolute accuracy
  eobp->ode_reltol=1e-11; //  ODE solver relative accuracy
  eobp->ode_tmax=1e9; // max integration time
  eobp->ode_stop_radius=1.; // stop ODE integration at this radius (if > 0)
  eobp->ode_stop_afterNdt=4; // stop ODE N iters after the Omega peak

  /* Spin dynamics */
  
  eobp->spin_dyn_size = 500; // size for spin dynamics arrays
  eobp->spin_odes_omg_stop = -1; // stop ODE integration at this Momega set by NR BBH mrg freq fit 
  eobp->spin_odes_t_stop = -1;
  eobp->spin_odes_dt = 1; // timestep for spin dynamics
  eobp->spin_interp_domain = 0; // FD or TD interpolation
  eobp->chi1x = eobp->chi1y = eobp->chi1z = 0.;
  eobp->chi2x = eobp->chi2y = eobp->chi2z = 0.;
  eobp->spin_flx=SPIN_FLX_EOB;

  /* OMP settings */

  eobp->openmp_threads=1; // OpenMP threads
  eobp->openmp_timeron=0; // OpenMP timers

   
  /* following pars are set later by the code */

  
  eobp->nu= 0. ; // symmetric mass ratio
  eobp->X1= 0. ; // mass ratio M1/M
  eobp->X2= 0. ; // mass ratio M2/M
  eobp->S1= 0. ;
  eobp->S2= 0. ;
  eobp->S= 0. ;
  eobp->Sstar= 0. ;
  eobp->a1= 0. ;
  eobp->a2= 0. ; // 
  eobp->aK= 0. ; //
  eobp->aK2= 0. ; // 
  eobp->C_Q1= 0. ; //
  eobp->C_Q2= 0. ; //
  eobp->C_Oct1= 0. ; //
  eobp->C_Oct2= 0. ; //
  eobp->C_Hex1= 0. ; //
  eobp->C_Hex2= 0. ; //
  eobp->a6c= 0. ; //
  eobp->cN3LO= 0. ; //
  
  eobp->kapA2= 0. ; // gravitoelectric kappa star A
  eobp->kapA3= 0. ; //
  eobp->kapA4= 0. ; //
  eobp->kapB2= 0. ; //
  eobp->kapB3= 0. ; //
  eobp->kapB4= 0. ; //
  eobp->kapT2= 0. ; //
  eobp->kapT3= 0. ; //
  eobp->kapT4= 0. ; //

  eobp->japA2= 0. ; // gravitomagnetic kappa star A
  //eobp->japA3= 0. ; 
  //eobp->japA4= 0. ; 
  eobp->japB2= 0. ;
  //eobp->japB3= 0. ;
  //eobp->japB4= 0. ;
  eobp->japT2= 0. ;
  //eobp->japT3= 0. ;
  //eobp->japT4= 0. ;  

  eobp->bar_alph2_1= 0. ; //
  eobp->bar_alph2_2= 0. ; //
  eobp->bar_alph3_1= 0. ; //
  eobp->bar_alph3_2= 0. ; //
  eobp->bar_alph2j_1= 0. ; //
  eobp->khatA2= 0. ; //
  eobp->khatB2= 0. ; //
  eobp->rLR_tidal= 0. ; // radius of light-ring for NNLO tidal model
  
  eobp->Mbhf= 0. ; // final BH mass
  eobp->abhf= 0. ; // final BH spin
  
  eobp->rLR= 0. ; // radius of light-ring
  eobp->rLSO= 0. ; // radius of last stable orbit 
  
  eobp->use_tidal= TIDES_OFF ; // index for tidal modus
  eobp->use_tidal_gravitomagnetic= TIDES_GM_OFF ; // index for gravitomagnetic tide

  if (choose == DEFAULT_PARS_BBH) {

    eobp->centrifugal_radius = CENTRAD_NLO;
    eobp->use_flm = USEFLM_HM;
    eobp->use_tidal = TIDES_OFF;
    eobp->use_tidal_gravitomagnetic = TIDES_GM_OFF;

    eobp->nqc = NQC_AUTO; // {"no", "auto", "manual"}
    eobp->nqc_coefs_flx = NQC_FLX_NRFIT_SPIN_202002; // {"none", "nrfit_nospin20160209", "nrfit_spin20202","fromfile"}
    eobp->nqc_coefs_hlm = NQC_HLM_COMPUTE; // {"compute", "none", "nrfit_nospin20160209", "nrfit_spin20202", "fromfile"}
     
  } else if (choose == DEFAULT_PARS_BNS) {

    eobp->use_tidal = TIDES_TEOBRESUM3;
    eobp->use_tidal_gravitomagnetic = TIDES_GM_PN;
    eobp->pGSF_tidal = 4.0;
    eobp->use_Yagi_fits = 1;

    eobp->centrifugal_radius = CENTRAD_NNLO;
    eobp->use_flm = USEFLM_SSNLO;
    eobp->nqc = NQC_NO; // {"no", "auto", "manual"}
    eobp->nqc_coefs_flx = NQC_FLX_NONE; // {"none", "nrfit_nospin20160209", "nrfit_spin20202", "fromfile"}
    eobp->nqc_coefs_hlm = NQC_HLM_NONE;

  }

  /*
  else if (choose == DEFAULT_PARS_BHNS) {
     
    //Still to be decided & set I guess?
    //Let's have a BH with m = 5 and a NS with m = 1.4
    eobp->M = 6.4;
    eobp->q = 3.6;
    eobp->chi1 = 0.;
    eobp->chi2 = 0.;
    eobp->use_spins = 0;
  }
  */

  else errorexit("unknown default parameter choice.");

}

/** overwrites the EOBParameters entries from the db values */
// This is the 'inverse' of par_db_default
void EOBParameters_set_from_db (EOBParameters *eobp)
{
  eobp->M = par_get_d("M"); // total binary mass  
  eobp->q = par_get_d("q"); // mass ratio
  eobp->chi1 = par_get_d("chi1"); // dimensionless spin 1 z-component 
  eobp->chi2 = par_get_d("chi2"); // dimensionless spin 2 z-component 
  eobp->distance = par_get_d("distance"); 
  eobp->inclination = par_get_d("inclination");
  eobp->coalescence_angle = par_get_d("coalescence_angle");
  eobp->polarization = par_get_d("polarisation");
  eobp->r0 = par_get_d("r0"); //initial radius	
  eobp->initial_frequency = par_get_d("initial_frequency"); // initial GW frequency
  
  eobp->LambdaAl2 = par_get_d("LambdaAl2"); // Tidal gravitoelectric parameter Lambda for star A ell=2
  eobp->LambdaBl2 = par_get_d("LambdaBl2");
  eobp->LambdaAl3 = par_get_d("LambdaAl3"); 
  eobp->LambdaBl3 = par_get_d("LambdaBl3"); 
  eobp->LambdaAl4 = par_get_d("LambdaAl4"); 
  eobp->LambdaBl4 = par_get_d("LambdaBl4"); 
  eobp->SigmaAl2 = par_get_d("SigmaAl2"); // Tidal gravitomagnetic parameter Sigma for star A ell=2
  eobp->SigmaBl2 = par_get_d("SigmaBl2");
  eobp->pGSF_tidal = par_get_d("pGSF_tidal");// p-power in GSF tidal potential model

  /* EOB Settings */
  eobp->use_spins = par_get_i("use_spins"); // use spins ? 

  for (eobp->use_tidal=0; eobp->use_tidal<TIDES_NOPT; eobp->use_tidal++) {
    if (STREQUAL(par_get_s("tides"),tides_opt[eobp->use_tidal])) {
      break;
    }
  }
  if (eobp->use_tidal == TIDES_NOPT) {
    eobp->use_tidal = TIDES_OFF;
    if (VERBOSE) printf("tides '%s' undefined, set to '%s'\n",
			par_get_s("tides"),tides_opt[eobp->use_tidal]);
  }

  for (eobp->use_tidal_gravitomagnetic=0; eobp->use_tidal_gravitomagnetic<TIDES_GM_NOPT; eobp->use_tidal_gravitomagnetic++) {
    if (STREQUAL(par_get_s("tides_gravitomagnetic"), tides_gravitomagnetic_opt[eobp->use_tidal_gravitomagnetic])) {
      break;
    }
  }
  if (eobp->use_tidal_gravitomagnetic == TIDES_GM_NOPT) {
    eobp->use_tidal_gravitomagnetic = TIDES_GM_OFF;
    if (VERBOSE) printf("tides GM '%s' undefined, set to '%s'\n",
			par_get_s("tides_gravitomagnetic"), tides_gravitomagnetic_opt[eobp->use_tidal_gravitomagnetic]);
  }
  eobp->use_Yagi_fits = par_get_i("use_Yagi_fits"); // use Yagi fit to obtain tidal parameters Lambda_3,4 from Lambda_2 ?
  eobp->use_geometric_units = par_get_i("use_geometric_units"); // use geometric units for I/O ?
  eobp->use_speedytail = par_get_i("use_speedytail"); // use special routine to speed up tail computation ?
  eobp->dt_merger_interp = par_get_d("dt_merger_interp"); // dt for interpolating merger waveform and NQC/ringdown attachment

  eobp->interp_uniform_grid = par_get_i("interp_uniform_grid"); // interpolate on uniform grid final result ?
  eobp->dt_interp = par_get_d("dt_interp"); // timestep to be used for final interpolation (used if input is given in geometrical unit, unused otherwise)
  eobp->srate_interp = par_get_d("srate_interp"); // sampling rate to be used for final interpolation (used if input is given in physical unit, unused otherwise)

  int *klm = par_get_arrayi("use_mode_lm", &eobp->use_mode_lm_size);
  free(eobp->use_mode_lm);
  eobp->use_mode_lm = malloc ( eobp->use_mode_lm_size * sizeof(int) );
  memcpy(eobp->use_mode_lm, klm, eobp->use_mode_lm_size * sizeof(int));
  free(klm);

  eobp->postadiabatic_dynamics = YESNO2INT(par_get_s("postadiabatic_dynamics")); //CHECK THIS MACRO ITS NEW... no=0, else yes=1
  eobp->postadiabatic_dynamics_N = par_get_i("postadiabatic_dynamics_N"); // post-adiabatic order
  eobp->postadiabatic_dynamics_size = par_get_i("postadiabatic_dynamics_size"); // grid size 
  eobp->postadiabatic_dynamics_rmin = par_get_d("postadiabatic_dynamics_rmin"); // minimum radius (end of PA dynamics)
  //eobp->postadiabatic_dynamics_rmin_BNS = par_get_d("postadiabatic_dynamics_rmin_BNS"); // minimum radius (end of PA dynamics)
  eobp->postadiabatic_dynamics_stop = YESNO2INT(par_get_s("postadiabatic_dynamics_stop")); // stop after post-adiabatic dynamics //FIXME: make bool

  for (eobp->centrifugal_radius=0; eobp->centrifugal_radius<CENTRAD_NOPT; eobp->centrifugal_radius++) {
    if (STREQUAL(par_get_s("centrifugal_radius"), centrifugal_radius_opt[eobp->centrifugal_radius])) {
      break;
    }
  }
  if (eobp->centrifugal_radius == CENTRAD_NOPT) {
    eobp->centrifugal_radius = CENTRAD_NLO;
    if (VERBOSE) printf("centrifugal_radius '%s' undefined, set to '%s'\n",
			par_get_s("centrifugal_radius"), centrifugal_radius_opt[eobp->centrifugal_radius]);
  }

  for (eobp->use_flm=0; eobp->use_flm<USEFLM_NOPT; eobp->use_flm++) {
    if (STREQUAL(par_get_s("use_flm"), use_flm_opt[eobp->use_flm])) {
      break;
    }
  }
  if (eobp->use_flm == USEFLM_NOPT) {
    eobp->use_flm = USEFLM_HM;
    if (VERBOSE) printf("use_flm '%s' undefined, set to '%s'\n",
			par_get_s("use_flm"), use_flm_opt[eobp->use_flm]);
  }

  eobp->compute_LR = par_get_i("compute_LR"); // calculate LR ?
  eobp->compute_LSO = par_get_i("compute_LSO"); // calculate LSO ?
  eobp->compute_LR_guess = par_get_d("compute_LR_guess");
  eobp->compute_LSO_guess = par_get_d("compute_LSO_guess");

  for (eobp->nqc=0; eobp->nqc<NQC_NOPT; eobp->nqc++) {
    if (STREQUAL(par_get_s("nqc"), nqc_opt[eobp->nqc])) {
      break;
    }
  }
  if (eobp->nqc == NQC_NOPT) {
    eobp->nqc = NQC_AUTO;
    if (VERBOSE) printf("nqc '%s' undefined, set to '%s'\n",
			par_get_s("use_flm"), nqc_opt[eobp->nqc]);
  }

  for (eobp->nqc_coefs_flx=0; eobp->nqc_coefs_flx<NQC_FLX_NOPT; eobp->nqc_coefs_flx++) {
    if (STREQUAL(par_get_s("nqc_coefs_flx"), nqc_flx_opt[eobp->nqc_coefs_flx])) {
      break;
    }
  }
  if (eobp->nqc_coefs_flx == NQC_FLX_NOPT) {
    eobp->nqc_coefs_flx = NQC_FLX_NONE;
    if (VERBOSE) printf("nqc '%s' undefined, set to '%s'\n",
			par_get_s("nqc_coefs_flx"), nqc_flx_opt[eobp->nqc_coefs_flx]);
  }

  for (eobp->nqc_coefs_hlm=0; eobp->nqc_coefs_hlm<NQC_HLM_NOPT; eobp->nqc_coefs_hlm++) {
    if (STREQUAL(par_get_s("nqc_coefs_hlm"), nqc_hlm_opt[eobp->nqc_coefs_hlm])) {
      break;
    }
  }
  if (eobp->nqc_coefs_hlm == NQC_HLM_NOPT) {
    eobp->nqc_coefs_hlm = NQC_HLM_NONE;
    if (VERBOSE) printf("nqc '%s' undefined, set to '%s'\n",
			par_get_s("nqc_coefs_hlm"), nqc_hlm_opt[eobp->nqc_coefs_hlm]);
  }

  strcpy(eobp->nqc_coefs_flx_file, par_get_s("nqc_coefs_flx_file"));
  strcpy(eobp->nqc_coefs_hlm_file, par_get_s("nqc_coefs_hlm_file"));

  /* Output */

  strcpy(eobp->output_dir, par_get_s("output_dir")); // output dir
  eobp->output_hpc = par_get_i("output_hpc"); // output h+,hx
  eobp->output_multipoles = par_get_i("output_multipoles"); // output multipoles

  int *olm = par_get_arrayi("output_lm", &eobp->output_lm_size);// indexes of multipoles to ouput
  free(eobp->output_lm);
  eobp->output_lm = malloc ( eobp->output_lm_size * sizeof(int) );
  memcpy(eobp->output_lm, olm, eobp->output_lm_size * sizeof(int));
  free(olm);

  eobp->output_dynamics = par_get_i("output_dynamics"); // output dynamics
  eobp->output_nqc = par_get_i("output_nqc"); // output NQC waveform
  eobp->output_nqc_coefs = par_get_i("output_nqc_coefs"); // output multipolar NQC coefs (if determined)
  eobp->output_ringdown = par_get_i("output_ringdown"); // output ringdown waveform
  
  /* Evolution settings */
  
  eobp->srate = par_get_d("srate"); // sampling rate, used if input is given in physical unit, reset based on tstep otherwise
  eobp->dt = par_get_d("dt"); // timestep, used if input is given in geometric unit, reset based on srate otherwise
  //eobp->size = par_get_i("size"); // size of the arrays (chunks, dynamically extended)
  eobp->ringdown_extend_array = par_get_i("ringdown_extend_array"); // grid points to extend arrays for ringdown attachment


  eobp->ode_timestep = -1; // specify ODE solver timestep "uniform","adaptive","adaptive+uniform_after_LSO","undefined"
  for (eobp->ode_timestep=0; eobp->ode_timestep<ODE_TSTEP_NOPT; eobp->ode_timestep++) {
    if (STREQUAL(par_get_s("ode_timestep"),ode_tstep_opt[eobp->ode_timestep])) {
      break;
    }
  }
  if (eobp->ode_timestep==ODE_TSTEP_NOPT) {
    eobp->ode_timestep = ODE_TSTEP_ADAPTIVE;
    if (VERBOSE) printf("ode_timestep '%s' undefined, set to default %s\n",par_get_s("ode_timestep"),
			ode_tstep_opt[eobp->ode_timestep]);
  }

  eobp->ode_abstol = par_get_d("ode_abstol"); // ODE solver absolute accuracy
  eobp->ode_reltol = par_get_d("ode_reltol"); //  ODE solver relative accuracy
  eobp->ode_tmax = par_get_d("ode_tmax"); // max integration time
  eobp->ode_stop_radius = par_get_d("ode_stop_at_radius"); // stop ODE integration at this radius (if > 0)
  eobp->ode_stop_afterNdt = par_get_i("ode_stop_afterNdt"); // stop ODE N iters after the Omega peak

  /* Spin dynamics */

  eobp->spin_dyn_size = par_get_i("spin_dyn_size"); 
  eobp->spin_odes_omg_stop = par_get_d("spin_odes_omg_stop"); // stop ODE integration at this Momega set by NR BBH mrg freq fit 
  eobp->spin_odes_t_stop = par_get_d("spin_odes_t_stop");
  eobp->spin_odes_dt = par_get_d("spin_odes_dt"); // timestep for spin dynamics
  eobp->spin_interp_domain = par_get_i("spin_interp_domain"); 
  eobp->project_spins = par_get_i("project_spins");           // spin projection in EOB dynamics
  eobp->chi1x = par_get_d("chi1x"); 
  eobp->chi1y = par_get_d("chi1y"); 
  eobp->chi1z = par_get_d("chi1z"); 
  eobp->chi2x = par_get_d("chi2x"); 
  eobp->chi2y = par_get_d("chi2y"); 
  eobp->chi2z = par_get_d("chi2z"); 
  for (eobp->spin_flx=0; eobp->spin_flx<SPIN_FLX_NOPT; eobp->spin_flx++) {
    if (STREQUAL(par_get_s("spin_flx"), spin_flx_opt[eobp->spin_flx])) {
      break;
    }
  }
  /* OMP settings */
  
  eobp->openmp_threads = par_get_i("openmp_threads");
  eobp->openmp_timeron = par_get_i("openmp_timeron"); 
}


/*
 * libconfig wrappers
 */

/** database for parameters */
config_t cfg, *cf;
config_setting_t *csroot, *csetting;

void par_db_init ()
{
  cf = &cfg;
  config_init(cf);
  csroot = config_root_setting(cf); 
}

void par_db_free ()
{
  config_destroy(cf);
}

void par_file_parse (const char *fname)
{    
  if (!(config_read_file(cf, fname))) {
    fprintf(stderr, "%s = %d - %s\n",
	    config_error_file(cf),
	    config_error_line(cf),
	    config_error_text(cf));
    config_destroy(cf);
    errorexit("Problem reading file.");
  }
}

void par_file_parse_merge (const char *fname)
{
  /* parse the file in a new cfg */
  config_t cfg1, *cf1;
  cf1 = &cfg1;
  config_init(cf1);
  if (!(config_read_file(cf1, fname))) {
    fprintf(stderr, "%s = %d - %s\n",
	    config_error_file(cf1),
	    config_error_line(cf1),
	    config_error_text(cf1));
    config_destroy(cf1);
    errorexit("Problem reading file.");
  }
  /* set parameter values in db */
  config_setting_t * csroot1 = config_root_setting(cf1);
  config_setting_t * cs1, * array;
  int np = config_setting_length(csroot1);
  int i, type, n, j;
  char key[STRLEN];
  for (i = 0; i < np; i++) {
    cs1  = config_setting_get_elem(csroot1, i);
    strcpy(key, config_setting_name(cs1));
    type = config_setting_type(cs1);
    if      (type==CONFIG_TYPE_INT)    par_set_i(key, config_setting_get_int(cs1));
    else if (type==CONFIG_TYPE_BOOL)   par_set_b(key, config_setting_get_bool(cs1));
    else if (type==CONFIG_TYPE_FLOAT)  par_set_d(key, config_setting_get_float(cs1));
    else if (type==CONFIG_TYPE_STRING) par_set_s(key, config_setting_get_string(cs1));
    else if (type==CONFIG_TYPE_ARRAY) {
      /* this is a hack, only for int type */
      if (cs1 != NULL) {
	n = config_setting_length(cs1);
	int * array = (int *) malloc (n * sizeof(int));
	if (!array) errorexit("out of memory");
	for (j = 0; j < n; j++) {
	  array[j] = config_setting_get_int_elem(cs1, j);
	}
	par_set_arrayi (key, array, n);
	free(array);
      } 
    } else errorexit("unkown parameter type.");
  }
  /* free extra cfg */
  config_destroy(cf1);
}

/* Set the libconfig db from EOBParamaters */
void par_db_from_EOBPar (EOBParameters *eobp)
{
  
  /* Physical input pars */
  
  par_add_d("M", eobp->M); // total binary mass  
  par_add_d("q", eobp->q); // mass ratio
  par_add_d("chi1", eobp->chi1); // dimensionless spin 1 z-component 
  par_add_d("chi2", eobp->chi2); // dimensionless spin 2 z-component 
  par_add_d("distance", eobp->distance); 
  par_add_d("inclination", eobp->inclination);
  par_add_d("coalescence_angle", eobp->coalescence_angle);
  par_add_d("polarisation", eobp->polarization);
  par_add_d("r0", eobp->r0); //initial radius	
  par_add_d("initial_frequency", eobp->initial_frequency); // initial GW frequency
  
  par_add_d("LambdaAl2", eobp->LambdaAl2); // Tidal gravitoelectric parameter Lambda for star A ell=2
  par_add_d("LambdaBl2",eobp->LambdaBl2);
  par_add_d("LambdaAl3",eobp->LambdaAl3); 
  par_add_d("LambdaAl4", eobp->LambdaAl4); 
  par_add_d("LambdaBl3", eobp->LambdaBl3); 
  par_add_d("LambdaBl4", eobp->LambdaBl4); 
  par_add_d("SigmaAl2", eobp->SigmaAl2); // Tidal gravitomagnetic parameter Sigma for star A ell=2
  par_add_d("SigmaBl2", eobp->SigmaBl2);
  par_add_d("pGSF_tidal", eobp->pGSF_tidal);// p-power in GSF tidal potential model
  
  /* EOB Settings */

  //int indexeslm[1] = {-1};  // indexes of multipoles to use in h+,hx (if [-1], use all) and to output
  char int2YESNO[2][STRLEN] = {"no", "yes"};

  par_add_b("use_spins", eobp->use_spins); // use spins ?
  par_add_s("tides", tides_opt[eobp->use_tidal]);
  par_add_s("tides_gravitomagnetic", tides_gravitomagnetic_opt[eobp->use_tidal_gravitomagnetic]);
  
  par_add_b("use_Yagi_fits", eobp->use_Yagi_fits); // use Yagi fit to obtain tidal parameters Lambda_3,4 from Lambda_2 ?
  par_add_b("use_geometric_units", eobp->use_geometric_units); // use geometric units for I/O ?
  par_add_b("use_speedytail", eobp->use_speedytail); // use special routine to speed up tail computation ?
  
  par_add_d("dt_merger_interp", eobp->dt_merger_interp); // dt for interpolating merger waveform and NQC/ringdown attachment
  par_add_b("interp_uniform_grid", eobp->interp_uniform_grid); // interpolate on uniform grid final result ?
  par_add_d("dt_interp", eobp->dt_interp); // timestep to be used for final interpolation (used if input is given in geometrical unit, unused otherwise)
  par_add_d("srate_interp", eobp->srate_interp); // sampling rate to be used for final interpolation (used if input is given in physical unit, unused otherwise)
  
  par_add_arrayi("use_mode_lm", eobp->use_mode_lm, eobp->use_mode_lm_size); // indexes of multipoles to use in h+,hx (if [-1], use all)

  par_add_s("postadiabatic_dynamics", int2YESNO[eobp->postadiabatic_dynamics]);
  par_add_i("postadiabatic_dynamics_N", eobp->postadiabatic_dynamics_N); // post-adiabatic order
  par_add_i("postadiabatic_dynamics_size", eobp->postadiabatic_dynamics_size); // grid size 
  par_add_d("postadiabatic_dynamics_rmin",eobp->postadiabatic_dynamics_rmin); // minimum radius (end of PA dynamics)
  par_add_d("postadiabatic_dynamics_rmin_BNS",eobp->postadiabatic_dynamics_rmin); // minimum radius (end of PA dynamics)
  par_add_s("postadiabatic_dynamics_stop",int2YESNO[eobp->postadiabatic_dynamics_stop]); // stop after post-adiabatic dynamics //FIXME: make bool

  par_add_s("centrifugal_radius", centrifugal_radius_opt[eobp->centrifugal_radius]); // {LO, NLO, NNLO, NNLOS4, NOSPIN, NOTIDES}
  par_add_s("use_flm", use_flm_opt[eobp->use_flm]); // "SSLO", "SSNLO", "HM"
  par_add_b("compute_LR", eobp->compute_LR); // calculate LR ?
  par_add_b("compute_LSO",eobp->compute_LSO); // calculate LSO ?
  par_add_d("compute_LR_guess", eobp->compute_LR_guess);
  par_add_d("compute_LSO_guess",eobp->compute_LSO_guess);

  par_add_s("nqc", nqc_opt[eobp->nqc]); // {"auto", "manual"}
  par_add_s("nqc_coefs_flx", nqc_flx_opt[eobp->nqc_coefs_flx]); // {"none", "nrfit_nospin20160209", "fromfile"}
  par_add_s("nqc_coefs_hlm", nqc_hlm_opt[eobp->nqc_coefs_hlm]); // {"compute", "none", "nrfit_nospin20160209", "fromfile"}
  par_add_s("nqc_coefs_flx_file", eobp->nqc_coefs_flx_file);
  par_add_s("nqc_coefs_hlm_file", eobp->nqc_coefs_hlm_file);

  /* Output */

  par_add_s("output_dir",eobp->output_dir); // output dir
  par_add_b("output_hpc", eobp->output_hpc); // output h+,hx
  par_add_b("output_multipoles", eobp->output_multipoles); // output multipoles
  par_add_arrayi("output_lm", eobp->output_lm, eobp->output_lm_size); // indexes of multipoles to ouput
  par_add_b("output_dynamics", eobp->output_dynamics); // output dynamics
  par_add_b("output_nqc",eobp->output_nqc); // output NQC waveform
  par_add_b("output_nqc_coefs", eobp->output_nqc_coefs); // output multipolar NQC coefs (if determined)
  par_add_b("output_ringdown", eobp->output_ringdown); // output ringdown waveform
  
  /* Evolution settings */
  
  par_add_d("srate", eobp->srate); // sampling rate, used if input is given in physical unit, reset based on tstep otherwise
  par_add_d("dt", eobp->dt); // timestep, used if input is given in geometric unit, reset based on srate otherwise
  par_add_i("size", eobp->size); // size of the arrays (chunks, dynamically extended)
  par_add_i("ringdown_extend_array", eobp->ringdown_extend_array); // grid points to extend arrays for ringdown attachment
  par_add_s("ode_timestep", ode_tstep_opt[eobp->ode_timestep]); // specify ODE solver timestep "uniform","adaptive","adaptive+uniform_after_LSO","undefined"
  par_add_d("ode_abstol",eobp->ode_abstol); // ODE solver absolute accuracy
  par_add_d("ode_reltol",eobp->ode_reltol); //  ODE solver relative accuracy
  par_add_d("ode_tmax", eobp->ode_tmax); // max integration time
  par_add_d("ode_stop_at_radius",eobp->ode_stop_radius); // stop ODE integration at this radius (if > 0)
  par_add_i("ode_stop_afterNdt",eobp->ode_stop_afterNdt); // stop ODE N iters after the Omega peak

  /* Spin dynamics */

  par_add_i("spin_dyn_size", eobp->spin_dyn_size);
  par_add_d("spin_odes_omg_stop", eobp->spin_odes_omg_stop);
  par_add_d("spin_odes_t_stop", eobp->spin_odes_t_stop);
  par_add_d("spin_odes_dt", eobp->spin_odes_dt);
  par_add_i("spin_interp_domain", eobp->spin_interp_domain);
  par_add_i("project_spins", eobp->project_spins);
  par_add_s("spin_flx", spin_flx_opt[eobp->spin_flx]); // {"PN", "EOB"}
  par_add_d("chi1x", eobp->chi1x);
  par_add_d("chi1y", eobp->chi1y);
  par_add_d("chi1z", eobp->chi1z);
  par_add_d("chi2x", eobp->chi2x);
  par_add_d("chi2y", eobp->chi2y);
  par_add_d("chi2z", eobp->chi2z);
    
  /* OMP settings */
  
  par_add_i("openmp_threads", eobp->openmp_threads); // OpenMP threads
  par_add_b("openmp_timeron", eobp->openmp_timeron); // OpenMP timers
  
  /* following pars are set later by the code */
  //TODO: they will be removed from the db and only kept in EOBParameters

  par_add_d("nu",  eobp->nu); // symmetric mass ratio
  par_add_d("X1",  eobp->X1); // mass ratio M1/M
  par_add_d("X2",  eobp->X2); // mass ratio M2/M
  par_add_d("S1",  eobp->S1);
  par_add_d("S2",  eobp->S2);
  par_add_d("S",  eobp->S);
  par_add_d("Sstar",  eobp->Sstar);
  par_add_d("a1",  eobp->a1);
  par_add_d("a2",  eobp->a2); // 
  par_add_d("aK",  eobp->aK); //
  par_add_d("aK2",  eobp->aK2); // 
  par_add_d("C_Q1",  eobp->C_Q1); //
  par_add_d("C_Q2",  eobp->C_Q2); //
  par_add_d("C_Oct1",  eobp->C_Oct1); //
  par_add_d("C_Oct2",  eobp->C_Oct2); //
  par_add_d("C_Hex1",  eobp->C_Hex1); //
  par_add_d("C_Hex2",  eobp->C_Hex2); //
  par_add_d("a6c",  eobp->a6c); //
  par_add_d("cN3LO",  eobp->cN3LO); //
  
  par_add_d("kappaAl2",  eobp->kapA2); // gravitoelectric kappa star A
  par_add_d("kappaAl3",  eobp->kapA3); //
  par_add_d("kappaAl4",  eobp->kapA4); //
  par_add_d("kappaBl2",  eobp->kapB2); //
  par_add_d("kappaBl3",  eobp->kapB3); //
  par_add_d("kappaBl4",  eobp->kapB4); //
  
  par_add_d("kappajAl2",  eobp->japA2); // gravitomagnetic kappa star A
  par_add_d("kappajBl2",  eobp->japB2);
  par_add_d("kappajTl2",  eobp->japT2);
  
  par_add_d("kappaTl2",  eobp->kapT2); //
  par_add_d("kappaTl3",  eobp->kapT3); //
  par_add_d("kappaTl4",  eobp->kapT4); //
  par_add_d("bar_alph2_1", eobp->bar_alph2_1); //
  par_add_d("bar_alph2_2", eobp->bar_alph2_2); //
  par_add_d("bar_alph3_1", eobp->bar_alph3_1); //
  par_add_d("bar_alph3_2", eobp->bar_alph3_2); //
  par_add_d("bar_alph2j_1", eobp->bar_alph2j_1); //
  par_add_d("khatAl2",  eobp->khatA2); //
  par_add_d("khatBl2",  eobp->khatB2); //
  par_add_d("rLR_tidal",  eobp->rLR_tidal); // radius of light-ring for NNLO tidal model
  
  par_add_d("BH_final_mass",  eobp->Mbhf); // final BH mass
  par_add_d("BH_final_spin",  eobp->abhf); // final BH spin
  
  par_add_d("rLR",  eobp->rLR); // radius of light-ring
  par_add_d("rLSO",  eobp->rLSO); // radius of last stable orbit 
  
  par_add_i("use_tidal", eobp->use_tidal); // index for tidal modus
  par_add_i("use_tidal_gravitomagnetic", eobp->use_tidal_gravitomagnetic); // index for gravitomagnetic tide

  return;    
}

/** default values for parameters */
// This sets parameters without using the external default.par
// See https://github.com/hyperrealm/libconfig/blob/master/examples/c/example3.c
// See par/default.par for what default parameters are needed
// Tested with current code
//DEPRECATED: this routine is not used in v2.0
//            not all the parameters of v2.0 are included here.
#if (0)
void par_db_default ()
{
  
  /* Physical input pars */
  
  par_add_d("M", 1.); // total binary mass  
  par_add_d("q", 1.); // mass ratio
  par_add_d("chi1", 0.); // dimensionless spin 1 z-component 
  par_add_d("chi2", 0.); // dimensionless spin 2 z-component 
  par_add_d("distance", 1.); 
  par_add_d("inclination", 0.);
  par_add_d("coalescence_angle", 0.);
  par_add_d("polarisation", 0.);
  par_add_d("r0", 0.); //initial radius	
  par_add_d("initial_frequency", 0.004); // initial GW frequency
  
  par_add_d("LambdaAl2", 0.); // Tidal gravitoelectric parameter Lambda for star A ell=2
  par_add_d("LambdaBl2",0.);
  par_add_d("LambdaAl3",0.); 
  par_add_d("LambdaAl4", 0.); 
  par_add_d("LambdaBl3", 0.); 
  par_add_d("LambdaBl4", 0.); 
  par_add_d("SigmaAl2", 0.); // Tidal gravitomagnetic parameter Sigma for star A ell=2
  par_add_d("SigmaBl2", 0.);
  par_add_d("pGSF_tidal", 4.0);// p-power in GSF tidal potential model

  /* EOB Settings */

  int indexeslm[1] = {-1};  // indexes of multipoles to use in h+,hx (if [-1], use all) and to output
  
  par_add_b("use_spins", 0); // use spins ?
  par_add_s("tides", "no");
  par_add_s("tides_gravitomagnetic","no");
  par_add_b("use_Yagi_fits", 0); // use Yagi fit to obtain tidal parameters Lambda_3,4 from Lambda_2 ?
  par_add_b("use_geometric_units", 1); // use geometric units for I/O ?
  par_add_b("use_speedytail", 0); // use special routine to speed up tail computation ?
  
  par_add_d("dt_merger_interp", 0.5); // dt for interpolating merger waveform and NQC/ringdown attachment
  par_add_b("interp_uniform_grid", 0); // interpolate on uniform grid final result ?
  par_add_d("dt_interp", 0.5); // timestep to be used for final interpolation (used if input is given in geometrical unit, unused otherwise)
  par_add_d("srate_interp", 4096.); // sampling rate to be used for final interpolation (used if input is given in physical unit, unused otherwise)
  
  par_add_arrayi("use_mode_lm", indexeslm, 1); // indexes of multipoles to use in h+,hx (if [-1], use all)
  
  par_add_s("postadiabatic_dynamics", "no");
  par_add_i("postadiabatic_dynamics_N", 8); // post-adiabatic order
  par_add_i("postadiabatic_dynamics_size", 1000); // grid size 
  par_add_d("postadiabatic_dynamics_rmin", 14.); // minimum radius (end of PA dynamics)
  par_add_d("postadiabatic_dynamics_rmin_BNS", 14.); // minimum radius (end of PA dynamics for BNS) //FIXME:this should be removed, use only above.
  par_add_s("postadiabatic_dynamics_stop","yes"); // stop after post-adiabatic dynamics //FIXME: make bool

  par_add_s("centrifugal_radius", "NLO"); // {LO, NLO, NNLO, NNLOS4, NOSPIN, NOTIDES}
  par_add_s("use_flm", "SSLO"); // "SSLO", "SSNLO", "HM"
  par_add_b("compute_LR", 0); // calculate LR ?
  par_add_b("compute_LSO", 0); // calculate LSO ?
  par_add_d("compute_LR_guess", 3.);
  par_add_d("compute_LSO_guess", 6.);

  par_add_s("nqc", "auto"); // {"auto", "manual"}
  par_add_s("nqc_coefs_flx", "none"); // {"none", "nrfit_nospin20160209", "fromfile"}
  par_add_s("nqc_coefs_hlm", "none"); // {"compute", "none", "nrfit_nospin20160209", "fromfile"}
  par_add_s("nqc_coefs_flx_file", "");
  par_add_s("nqc_coefs_hlm_file", "");

  /* Output */

  par_add_s("output_dir", "./data/"); // output dir
  par_add_b("output_hpc", 1); // output h+,hx
  par_add_b("output_multipoles", 0); // output multipoles
  par_add_arrayi("output_lm", indexeslm, 1); // indexes of multipoles to ouput
  par_add_b("output_dynamics", 0); // output dynamics
  par_add_b("output_nqc", 0); // output NQC waveform
  par_add_b("output_nqc_coefs", 0); // output multipolar NQC coefs (if determined)
  par_add_b("output_ringdown", 0); // output ringdown waveform
  
  /* Evolution settings */
  
  par_add_d("srate", 4096.); // sampling rate, used if input is given in physical unit, reset based on tstep otherwise
  par_add_d("dt", 0.5); // timestep, used if input is given in geometric unit, reset based on srate otherwise
  par_add_i("size", 1); // size of the arrays (chunks, dynamically extended)
  par_add_i("ringdown_extend_array", 500); // grid points to extend arrays for ringdown attachment
  par_add_s("ode_timestep", "uniform"); // specify ODE solver timestep "uniform","adaptive","adaptive+uniform_after_LSO","undefined"
  par_add_d("ode_abstol", 1e-13); // ODE solver absolute accuracy
  par_add_d("ode_reltol", 1e-11); //  ODE solver relative accuracy
  par_add_d("ode_tmax", 1e9); // max integration time
  par_add_d("ode_stop_at_radius", 2.); // stop ODE integration at this radius (if > 0)
  par_add_i("ode_stop_afterNdt", 4); // stop ODE N iters after the Omega peak
  
  /* OMP settings */
  
  par_add_i("openmp_threads", 1); // OpenMP threads
  par_add_b("openmp_timeron", 0); // OpenMP timers
  
  /* following pars are set later by the code */
  // TODO: they will be removed from the db and only kept in EOBParameters

  par_add_d("nu", 0.); // symmetric mass ratio
  par_add_d("X1", 0.); // mass ratio M1/M
  par_add_d("X2", 0.); // mass ratio M2/M
  par_add_d("S1", 0.);
  par_add_d("S2", 0.);
  par_add_d("S", 0.);
  par_add_d("Sstar", 0.);
  par_add_d("a1", 0.);
  par_add_d("a2", 0.); // 
  par_add_d("aK", 0.); //
  par_add_d("aK2", 0.); // 
  par_add_d("C_Q1", 0.); //
  par_add_d("C_Q2", 0.); //
  par_add_d("C_Oct1", 0.); //
  par_add_d("C_Oct2", 0.); //
  par_add_d("C_Hex1", 0.); //
  par_add_d("C_Hex2", 0.); //
  par_add_d("a6c", 0.); //
  par_add_d("cN3LO", 0.); //
  
  par_add_d("kappaAl2", 0.); // gravitoelectric kappa star A
  par_add_d("kappaAl3", 0.); //
  par_add_d("kappaAl4", 0.); //
  par_add_d("kappaBl2", 0.); //
  par_add_d("kappaBl3", 0.); //
  par_add_d("kappaBl4", 0.); //
  
  par_add_d("kappajAl2", 0.); // gravitomagnetic kappa star A
  par_add_d("kappajBl2", 0.);
  par_add_d("kappajTl2", 0.);
  
  par_add_d("kappaTl2", 0.); //
  par_add_d("kappaTl3", 0.); //
  par_add_d("kappaTl4", 0.); //
  par_add_d("bar_alph2_1", 0.); //
  par_add_d("bar_alph2_2", 0.); //
  par_add_d("bar_alph3_1", 0.); //
  par_add_d("bar_alph3_2", 0.); //
  par_add_d("bar_alph2j_1", 0.); //
  par_add_d("khatAl2", 0.); //
  par_add_d("khatBl2", 0.); //
  par_add_d("rLR_tidal", 0.); // radius of light-ring for NNLO tidal model
  
  par_add_d("BH_final_mass", 0.); // final BH mass
  par_add_d("BH_final_spin", 0.); // final BH spin
  
  par_add_d("rLR", 0.); // radius of light-ring
  par_add_d("rLSO", 0.); // radius of last stable orbit 
  
  par_add_i("use_tidal", 0); // index for tidal modus
  par_add_i("use_tidal_gravitomagnetic", 0); // index for gravitomagnetic tide

  return;    
}
#endif 

/** default values for parameters are expected in file
    $TEOBRESUMS/par/default.par */
void par_db_default_fromfile ()
{
  static const char *eobcodeenvv = "TEOBRESUMS";
  char *eobcodepath = NULL;
  eobcodepath = getenv (eobcodeenvv);
  if (!eobcodepath) {
    printf ("Enviroment variable '%s' is not set\n", eobcodeenvv);
    errorexit("Default values for parameters cannot be set.\n");
  }
  else {
    par_file_parse (strcat(eobcodepath,"/par/default.par"));
  }
}

void par_db_write_file (const char *name)
{
  char fname[STRLEN];
  strcpy(fname,par_get_s("output_dir"));
  strcat(fname,"/");
  strcat(fname,name);
  if (!(config_write_file(cf, fname))) {
    fprintf(stderr, "Error writing file %s\n", fname);
    config_destroy(cf);
    errorexit("Problem writing file");
  }
}

void par_db_screen (const int pr)
{
  if (pr) config_write(cf, stderr);
}

/* 'add' routines */

void par_add_i (const char *key, int val)
{
  csetting = config_setting_add(csroot, key, CONFIG_TYPE_INT);
  config_setting_set_int(csetting, val);
  return;
}

void par_add_b (const char *key, int val)
{
  csetting = config_setting_add(csroot, key, CONFIG_TYPE_INT);
  config_setting_set_int(csetting, val);
  /* csetting = config_setting_add(csroot, key, CONFIG_TYPE_BOOL); */
  /* config_setting_set_bool(csetting, val); */
  return;
}

void par_add_d (const char *key, double val)
{
  csetting = config_setting_add(csroot, key, CONFIG_TYPE_FLOAT);
  config_setting_set_float(csetting, val);
  return;
}

void par_add_s (const char *key, const char *val)
{
  csetting = config_setting_add(csroot, key, CONFIG_TYPE_STRING);
  config_setting_set_string(csetting, val);
  return;
}

void par_add_arrayi (const char *key, int *a, int size)
{
  config_setting_t * array;
  array = config_setting_add(csroot, key, CONFIG_TYPE_ARRAY);
  for(int i = 0; i < size; ++i) {
    csetting = config_setting_add(array, NULL, CONFIG_TYPE_INT);
    config_setting_set_int(csetting, a[i]);
  }
  return;
}

void par_add_arrayd (const char *key, double *a, int size)
{
  config_setting_t * array;
  array = config_setting_add(csroot, key, CONFIG_TYPE_ARRAY);
  for(int i = 0; i < size; ++i) {
    csetting = config_setting_add(array, NULL, CONFIG_TYPE_FLOAT);
    config_setting_set_float(csetting, a[i]);
  }
  return;
}

/* 'get' routines */

int par_get_i (const char *key)
{
  int val;
  if (config_lookup_int(cf, key, &val))
    return val;
  else 
    errorexits("unknown parameter", key);
}

int par_get_b (const char *key)
{
  int val;
  if (config_lookup_bool(cf, key, &val))
    return val;
  else 
    errorexits("unknown parameter", key);
} 

double par_get_d (const char *key)
{
  double val;
  if (config_lookup_float(cf, key, &val))
    return val;
  else 
    errorexits("unknown parameter",key);
}

const char * par_get_s (const char *key)
{
  const char *val;
  if (config_lookup_string(cf, key, &val))
    return val;
  else 
    errorexits("unknown parameter",key);
}

int * par_get_arrayi (const char *key, int *n)
{
  const config_setting_t *cs = config_lookup(cf, key);
  if (cs != NULL) {
    int l = config_setting_length(cs);
    int * array = (int *) malloc (l * sizeof(int));
    if (!array) errorexit("out of memory");
    int i;
    for (i = 0; i < l; i++) {
      array[i] = config_setting_get_int_elem(cs, i);
    }
    *n = l;
    return array;
  } else
    errorexits("unknown parameter",key);
}

double * par_get_arrayd (const char *key, int *n)
{
  const config_setting_t *cs = config_lookup(cf, key);
  if (cs != NULL) {
    int l = config_setting_length(cs);
    double * array = (double *) malloc (l * sizeof(int));
    if (!array) errorexit("out of memory");
    int i;
    for (i = 0; i < l; i++) {
      array[i] = config_setting_get_int_elem(cs, i);
    }
    *n = l;
    return array;
  } else
    errorexits("unknown parameter",key);
}

/* 'set' routines, assume all settings are root */

void par_set_i (const char *key, int val)
{
  config_setting_t * cs = config_lookup(cf, key);
  if (!(config_setting_set_int(cs, val)))
    errorexits("unknown parameter/wrong type for",key);
}

void par_set_b (const char *key, int val)
{
  config_setting_t * cs = config_lookup(cf, key);
  if (!(config_setting_set_bool(cs, val)))
    errorexits("unknown parameter/wrong type for",key);
} 

void par_set_d (const char *key, double val)
{
  config_setting_t * cs = config_lookup(cf, key);
  if (!(config_setting_set_float(cs, val)))
    errorexits("unknown parameter/wrong type for",key);
}

void par_set_s (const char *key, const char *val)
{
  config_setting_t * cs = config_lookup(cf, key);
  if (!(config_setting_set_string(cs, val)))
    errorexits("unknown parameter/wrong type for",key);
}

void par_set_arrayi (const char *key, int *array, int n)
{
  config_setting_t *cs = config_lookup(cf, key);
  config_setting_t *parent = config_setting_parent(cs);
  int i;
  /* remove what is there */
  if (!config_setting_remove(parent, key))
    errorexits("error removing array",key);
  /* add the new */
  config_setting_t *sarray, *a;
  sarray = config_setting_add(parent, key, CONFIG_TYPE_ARRAY);
  for(i = 0; i < n; ++i) {
    a = config_setting_add(sarray, NULL, CONFIG_TYPE_INT);
    config_setting_set_int(a, array[i]);
  }
}

void par_set_arrayd (const char *key, double *array, int n)
{
  config_setting_t *cs = config_lookup(cf, key);
  config_setting_t *parent = config_setting_parent(cs);
  int i;
  /* remove what is there */
  if (!config_setting_remove(parent, key))
    errorexits("error removing array",key);
  /* add the new */
  config_setting_t *sarray, *a;
  sarray = config_setting_add(parent, key, CONFIG_TYPE_ARRAY);
  for(i = 0; i < n; ++i) {
    a = config_setting_add(sarray, NULL, CONFIG_TYPE_FLOAT);
    config_setting_set_int(a, array[i]);
  }
}


/*
 * main routine to set parameters for the run
 */

void eob_set_params(EOBParameters *eobp, int default_choice, int firstcall)
{
  
  /* Set intrinsic parameters as given by user */

  const double M =  eobp->M;
  const double fmin = eobp->initial_frequency;
  const double q =  eobp->q;
  eobp->nu = q_to_nu(q);
  eobp->X1 = nu_to_X1(eobp->nu);
  eobp->X2 = 1. -  eobp->X1;
  const double XA = eobp->X1; /* tidal part used different notation, used here for simplicity */
  const double XB = eobp->X2;

  const int usespins = eobp->use_spins;
  const int usetidal = eobp->use_tidal;
  const int usetidalGM =  eobp->use_tidal_gravitomagnetic;

  /* Spin parameters */
  if (usespins==MODE_SPINS_GENERIC) {
    eobp->chi1 = eobp->chi1z;
    eobp->chi2 = eobp->chi2z;
    eobp->ringdown_eulerangles = RD_EULERANGLES_CONSTANT;
  }
  
  const double chi1 = eobp->chi1;
  const double chi2 = eobp->chi2;
  set_spin_vars(XA,XB, eobp->chi1,eobp->chi2, 
		&eobp->S1, &eobp->S2,
		&eobp->a1, &eobp->a2,
		&eobp->aK, &eobp->aK2,
		&eobp->S, &eobp->Sstar);
  /*
    eobp->S1 = SQ(XA) * chi1;
    eobp->S2 = SQ(XB) * chi2;
    eobp->a1 = XA*chi1;
    eobp->a2 = XB*chi2;
    eobp->aK = eobp->a1 +  eobp->a2;
    eobp->aK2 = SQ(eobp->aK);   
    eobp->S = eobp->S1 +  eobp->S2;            // in the EMRL this becomes the spin of the BH 
    eobp->Sstar = XB* eobp->a1 + XA* eobp->a2;  // in the EMRL this becomes the spin of the particle 
  */

  /* Tidal parameters */
  
  if (eobp->use_Yagi_fits) {
    eobp->LambdaAl3 = Yagi13_fit_barlamdel(eobp->LambdaAl2, 3);
    eobp->LambdaBl3 = Yagi13_fit_barlamdel(eobp->LambdaBl2, 3);
    eobp->LambdaAl4 = Yagi13_fit_barlamdel(eobp->LambdaAl2, 4);
    eobp->LambdaBl4 = Yagi13_fit_barlamdel(eobp->LambdaBl2, 4);
  }

#if(USEGRAVITOMAGNETICTERMS)
 eobp->SigmaAl2 = JFAPG_fit_Sigma_Irrotational(eobp->LambdaAl2);
 eobp->SigmaBl2 = JFAPG_fit_Sigma_Irrotational(eobp->LambdaBl2);
#endif

  /* Tidal coupling constants */    
  eobp->kapA2 = 3.   * eobp->LambdaAl2 * XA*XA*XA*XA*XA / q; 
  eobp->kapA3 = 15.  * eobp->LambdaAl3 * XA*XA*XA*XA*XA*XA*XA / q;
  eobp->kapA4 = 105. * eobp->LambdaAl4 * XA*XA*XA*XA*XA*XA*XA*XA*XA / q;
  
  eobp->kapB2 = 3.   * eobp->LambdaBl2 * XB*XB*XB*XB*XB * q;
  eobp->kapB3 = 15.  * eobp->LambdaBl3 * XB*XB*XB*XB*XB*XB*XB * q;
  eobp->kapB4 = 105. * eobp->LambdaBl4 * XB*XB*XB*XB*XB*XB*XB*XB*XB * q;

  /* gravitomagnetic tidal coupling constants el = 2 only */    
  eobp->japA2 = 24.   * eobp->SigmaAl2 * XA*XA*XA*XA*XA / q; 
  eobp->japB2 = 24.   * eobp->SigmaBl2 * XB*XB*XB*XB*XB * q;

  eobp->kapT2 = eobp->kapA2 + eobp->kapB2;
  eobp->kapT3 = eobp->kapA3 + eobp->kapB3;
  eobp->kapT4 = eobp->kapA4 + eobp->kapB4;
  
  eobp->japT2 = eobp->japA2 + eobp->japB2;

  if (usetidal) {
    if (!(eobp->kapT2 > 0.)) errorexit("kappaT2 must be >0");
    if (!(eobp->kapT3 > 0.)) errorexit("kappaT3 must be >0");
    if (!(eobp->kapT4 > 0.)) errorexit("kappaT4 must be >0");
  } 

  /* Tidal coefficients cons dynamics
     \bar{\alpha}_n^{(\ell)}, Eq.(37) of Damour&Nagar, PRD 81, 084016 (2010) */
  if(usetidal){
    eobp->bar_alph2_1 = (5./2.*XA*eobp->kapA2 + 5./2.*XB*eobp->kapB2)/eobp->kapT2;
    eobp->bar_alph2_2 = ((3.+XA/8.+ 337./28.*XA*XA)*eobp->kapA2 + (3.+XB/8.+ 337./28.*XB*XB)*eobp->kapB2)/eobp->kapT2;
    eobp->bar_alph3_1 = ((-2.+15./2.*XA)*eobp->kapA3 + (-2.+15./2.*XB)*eobp->kapB3)/eobp->kapT3;
    eobp->bar_alph3_2 = ((8./3.-311./24.*XA+110./3.*XA*XA)*eobp->kapA3 + (8./3.-311./24.*XB+110./3.*XB*XB)*eobp->kapB3)/eobp->kapT3;
  /* Gravitomagnetic term, see Eq.(6.27) of Bini-Damour-Faye 2012 */
    eobp->bar_alph2j_1 = ( eobp->japA2*(1. + (11./6.)*XA + XA*XA) + eobp->japB2*(1. + (11./6.)*XB + XB*XB) )/eobp->japT2;
  }

  /* Tidal coefficients for the amplitude */
  eobp->khatA2  = 3./2. * eobp->LambdaAl2 * XB/XA * gsl_pow_int(XA,5);
  eobp->khatB2  = 3./2. * eobp->LambdaBl2 * XA/XB * gsl_pow_int(XB,5);
  
  /* Self-spin coefficients */
  eobp->C_Q1   = 1.;
  eobp->C_Q2   = 1.;
  eobp->C_Oct1 = 1.;
  eobp->C_Oct2 = 1.;
  eobp->C_Hex1 = 1.;
  eobp->C_Hex2 = 1.;
  if (eobp->LambdaAl2>0.) {
    double logC_Q1 = logQ(log(eobp->LambdaAl2));
    eobp->C_Q1           = exp(logC_Q1);
    eobp->C_Oct1         = Yagi14_fit_Coct(eobp->C_Q1);
    eobp->C_Hex1         = Yagi14_fit_Chex(eobp->C_Q1);
  }
  if (eobp->LambdaBl2>0.) {
    double logC_Q2 = logQ(log(eobp->LambdaBl2));
    eobp->C_Q2           = exp(logC_Q2);
    eobp->C_Oct2         = Yagi14_fit_Coct(eobp->C_Q2);
    eobp->C_Hex2         = Yagi14_fit_Chex(eobp->C_Q2);
  }

  /* Default settings for NQC */
  // NOTE: The defaults are different from v0.0 and v1.0
  if (eobp->nqc == NQC_AUTO) {
    if (usetidal) {
      eobp->nqc_coefs_flx = NQC_FLX_NONE;
      eobp->nqc_coefs_hlm = NQC_HLM_NONE;
    } else {
      if (usespins) {
	eobp->nqc_coefs_flx = NQC_FLX_NRFIT_SPIN_202002;
	eobp->nqc_coefs_hlm = NQC_HLM_COMPUTE;
      } else {
	eobp->nqc_coefs_flx = NQC_FLX_NRFIT_NOSPIN_201602;
	eobp->nqc_coefs_hlm = NQC_HLM_NRFIT_NOSPIN_201602;
      }
    }
  } 
  
  /** Set more as needed ... */
  eobp->a6c = 0.;
  if (eobp->use_flm == USEFLM_HM) {
    /* Higher modes */
    eobp->a6c = eob_a6c_fit_HM(eobp->nu);
  } else {
    eobp->a6c = eob_a6c_fit(eobp->nu);
  }

  eobp->cN3LO = 0.;
  if (usetidal) eobp->cN3LO = 0.0;
  else if (eobp->use_flm == USEFLM_HM) {
    eobp->cN3LO = eob_c3_fit_HM(eobp->nu,eobp->a1,eobp->a2);
  } else {
    eobp->cN3LO = eob_c3_fit_global(eobp->nu,eobp->a1,eobp->a2);
  }

  double dt = eobp->dt;
  if (eobp->use_geometric_units) {
    /* input given in geometric units, 
       rescale to geometric units and mass rescaled quantities
       compute r0 from the initial GW frequency in geometric units and mass rescaled 
       reset sample rate using dt
    */
    if (VERBOSE) printf("Assume geometric units for pars values\n");
    eobp->r0 = pow(fmin*Pi, -2./3.);
    eobp->srate = 1./dt;
    eobp->distance = 1.;
    eobp->M = 1.;
  } else {
    /* input given in physical units, 
       rescale to geometric units and mass rescaled quantities
       compute r0 from the initial GW frequency in Hz 
    */
    if (VERBOSE) printf("Assume physical units for pars values\n");
    /* Set interpolation dt */
    dt = 1./eobp->srate_interp;
    dt = time_units_conversion(M, dt);
    eobp->dt_interp = dt;
    /* Set dt */
    eobp->r0 = radius0(M, fmin);
    dt = 1./eobp->srate;
    dt = time_units_conversion(M, dt);
    eobp->dt = dt;
    if (VERBOSE) PRFORMd("dt",eobp->dt);
    if (VERBOSE)
      if (eobp->interp_uniform_grid)
	PRFORMd("dt_interp",eobp->dt_interp);
  }

  /* Function pointers */
  
  /** Set f_lm fun pointer */
  if (eobp->use_flm == USEFLM_HM) {
    eob_wav_hlmNewt  = &eob_wav_hlmNewt_HM;
    eob_wav_flm      = &eob_wav_flm_HM;
    eob_wav_flm_s    = &eob_wav_flm_s_HM;
    eob_wav_deltalm  = &eob_wav_deltalm_HM;
    eob_wav_hlmNQC_find_a1a2a3_mrg = &eob_wav_hlmNQC_find_a1a2a3_mrg_HM;
    eob_wav_ringdown = &eob_wav_ringdown_HM;
  } else if (eobp->use_flm == USEFLM_SSLO) {
    /* eob_wav_flm_s = &eob_wav_flm_s_old; */
    eob_wav_hlmNewt  = &eob_wav_hlmNewt_v1;
    eob_wav_flm      = &eob_wav_flm_v1;
    eob_wav_flm_s    = &eob_wav_flm_s_SSLO;
    eob_wav_deltalm  = &eob_wav_deltalm_v1;
    eob_wav_hlmNQC_find_a1a2a3_mrg = &eob_wav_hlmNQC_find_a1a2a3_mrg_22;
    eob_wav_ringdown = &eob_wav_ringdown_v1;
  } else if (eobp->use_flm == USEFLM_SSNLO) {
    eob_wav_hlmNewt  = &eob_wav_hlmNewt_v1;
    eob_wav_flm      = &eob_wav_flm_v1;
    eob_wav_flm_s    = &eob_wav_flm_s_SSNLO;
    eob_wav_deltalm  = &eob_wav_deltalm_v1;
    eob_wav_hlmNQC_find_a1a2a3_mrg = &eob_wav_hlmNQC_find_a1a2a3_mrg_22;
    eob_wav_ringdown = &eob_wav_ringdown_v1;
    /*
      } else if (eobp->use_flm == USEFLM_SSNNLO) {
      eob_wav_hlmNewt = &eob_wav_hlmNewt_v1;
      eob_wav_flm     = &eob_wav_flm_v1;
      eob_wav_flm_s   = &eob_wav_flm_s_SSNNLO;
    */
  } else errorexit("unknown option for use_flm");

  /** Set rc fun pointer */
  if (eobp->centrifugal_radius == CENTRAD_LO) {
    eob_dyn_s_get_rc = &eob_dyn_s_get_rc_LO;
  } else if (eobp->centrifugal_radius == CENTRAD_NLO) {
    eob_dyn_s_get_rc = &eob_dyn_s_get_rc_NLO;
  } else if (eobp->centrifugal_radius == CENTRAD_NNLO) {
    eob_dyn_s_get_rc = &eob_dyn_s_get_rc_NNLO;
  } else if (eobp->centrifugal_radius == CENTRAD_NNLOS4) {
    eob_dyn_s_get_rc = &eob_dyn_s_get_rc_NNLO_S4;
  } else if (eobp->centrifugal_radius == CENTRAD_NOSPIN) {
    eob_dyn_s_get_rc = &eob_dyn_s_get_rc_NOSPIN;
  } else if (eobp->centrifugal_radius == CENTRAD_NOTIDES) {
    eob_dyn_s_get_rc = &eob_dyn_s_get_rc_NOTIDES;
  } else errorexit("unknown option for centrifugal_radius");

}

/** Set parameters */
/* Old routine working with the parfile/libconfig/parameter DB */
//DEPRECATED
#if (0)
void eob_set_params_old(char *s, int n)
{
  /* init db */
  par_db_init ();

  /* Set default values */
  par_db_default_fromfile ();

  /* Parse input parfile */
  if (s!=NULL) {
    /* Deal with input parfile if necessary 
    (this is the current logic, we'll keep it for compatibility) */
    par_file_parse_merge (s);
  }
  const int usespins = par_get_i("use_spins");

  /* Set auxiliary parameters */
  double M = par_get_d("M");
  double fmin = par_get_d("initial_frequency");
  double q =  par_get_d("q");  
  double nu = q_to_nu(q);
  double X1 = nu_to_X1(nu);
  double X2 = 1. - X1;
  double XA = X1; /* tidal part used different notation, used here for simplicity */
  double XB = X2;
  
  double chi1 = par_get_d("chi1");
  double chi2 = par_get_d("chi2");
  double S1 = SQ(X1) * chi1;
  double S2 = SQ(X2) * chi2;
  double a1  = X1*chi1;
  double a2  = X2*chi2;
  double aK  = a1 + a2;
  double aK2 = aK*aK;   
  double S = S1 + S2;            /* in the EMRL this becomes the spin of the BH */
  double Sstar = X2*a1 + X1*a2;  /* in the EMRL this becomes the spin of the particle */

  int j;

  for (j=0; j<=TIDES_NOPT; j++) {
    if (STREQUAL(par_get_s("tides"),tides_opt[j])) {
      /* if (DEBUG) printf("%-40s = %s\n","tides",tides_opt[j]); */
      break;
    }
  }
  if (j==TIDES_NOPT) {
    if (VERBOSE) printf("tides '%s' undefined, set to default\n",par_get_s("tides"));
    j = TIDES_OFF;
  }
  par_set_i("use_tidal",j);
  const int usetidal = j;

  int k;

  for (k=0; k<=TIDES_GM_NOPT; k++) {
    if (STREQUAL(par_get_s("tides_gravitomagnetic"), tides_gravitomagnetic_opt[k])) {
      /* if (DEBUG) printf("%-40s = %s\n","tides_gravitomagnetic", tides_gravitomagnetic_opt[k]); */
      break;
    }
  }
  if (k==TIDES_GM_NOPT) {
    if (VERBOSE) printf("gravitomagnetic tides '%s' undefined, set to default\n",par_get_s("tides_gravitomagnetic"));
    k = TIDES_GM_OFF;
  }
  par_set_i("use_tidal_gravitomagnetic",k);
  const int usetidalGM = k;  
  
  double LambdaAl2 = par_get_d("LambdaAl2");
  double LambdaBl2 = par_get_d("LambdaBl2");
  double LambdaAl3 = par_get_d("LambdaAl3");
  double LambdaBl3 = par_get_d("LambdaBl3");
  double LambdaAl4 = par_get_d("LambdaAl4");
  double LambdaBl4 = par_get_d("LambdaBl4");
  double SigmaAl2 = par_get_d("SigmaAl2");
  double SigmaBl2 = par_get_d("SigmaBl2");

  if (par_get_i("use_Yagi_fits")) {
    LambdaAl3 = Yagi13_fit_barlamdel(LambdaAl2, 3);
    LambdaBl3 = Yagi13_fit_barlamdel(LambdaBl2, 3);
    LambdaAl4 = Yagi13_fit_barlamdel(LambdaAl2, 4);
    LambdaBl4 = Yagi13_fit_barlamdel(LambdaBl2, 4);
  }

#if(USEGRAVITOMAGNETICTERMS)
    SigmaAl2 = JFAPG_fit_Sigma_Irrotational(LambdaAl2);
    SigmaBl2 = JFAPG_fit_Sigma_Irrotational(LambdaBl2);
#endif

  /* Tidal coupling constants */    
  double kapA2 = 3.   * LambdaAl2 * XA*XA*XA*XA*XA / q; 
  double kapA3 = 15.  * LambdaAl3 * XA*XA*XA*XA*XA*XA*XA / q;
  double kapA4 = 105. * LambdaAl4 * XA*XA*XA*XA*XA*XA*XA*XA*XA / q;
  
  double kapB2 = 3.   * LambdaBl2 * XB*XB*XB*XB*XB * q;
  double kapB3 = 15.  * LambdaBl3 * XB*XB*XB*XB*XB*XB*XB * q;
  double kapB4 = 105. * LambdaBl4 * XB*XB*XB*XB*XB*XB*XB*XB*XB * q;

  /* gravitomagnetic tidal coupling constants el = 2 only */    
  double kapA2j = 24.   * SigmaAl2 * XA*XA*XA*XA*XA / q; 
  double kapB2j = 24.   * SigmaBl2 * XB*XB*XB*XB*XB * q;

  double kapT2 = kapA2 + kapB2;
  double kapT3 = kapA3 + kapB3;
  double kapT4 = kapA4 + kapB4;
  double kapT2j = kapA2j + kapB2j;

  if (usetidal) {
    if (!(kapT2 > 0.)) errorexit("kappaT2 must be >0");
    if (!(kapT3 > 0.)) errorexit("kappaT3 must be >0");
    if (!(kapT4 > 0.)) errorexit("kappaT4 must be >0");
  } 
  
  /* Tidal coefficients cons dynamics
     \bar{\alpha}_n^{(\ell)}, Eq.(37) of Damour&Nagar, PRD 81, 084016 (2010) */
  double bar_alph2_1  = 0.;
  double bar_alph2_2  = 0.;
  double bar_alph3_1  = 0.;
  double bar_alph3_2  = 0.;
  double bar_alph2j_1 = 0.;

  if (usetidal){
    bar_alph2_1 = (5./2.*XA*kapA2 + 5./2.*XB*kapB2)/kapT2;
    bar_alph2_2 = ((3.+XA/8.+ 337./28.*XA*XA)*kapA2 + (3.+XB/8.+ 337./28.*XB*XB)*kapB2)/kapT2;
    bar_alph3_1 = ((-2.+15./2.*XA)*kapA3 + (-2.+15./2.*XB)*kapB3)/kapT3;
    bar_alph3_2 = ((8./3.-311./24.*XA+110./3.*XA*XA)*kapA3 + (8./3.-311./24.*XB+110./3.*XB*XB)*kapB3)/kapT3;
    /* Gravitomagnetic term, see Eq.(6.27) if Bini-Damour-Faye 2012 */
    bar_alph2j_1 = ( kapA2j*(1. + (11./6.)*XA + XA*XA) + kapB2j*(1. + (11./6.)*XB + XB*XB) )/kapT2j;
  }

  /* Tidal coefficients for the amplitude */
  double khatA_2  = 3./2. * LambdaAl2 * XB/XA * gsl_pow_int(XA,5);
  double khatB_2  = 3./2. * LambdaBl2 * XA/XB * gsl_pow_int(XB,5);
  
  /* Self-spin coefficients */
  double C_Q1   = 1.;
  double C_Q2   = 1.;
  double C_Oct1 = 1.;
  double C_Oct2 = 1.;
  double C_Hex1 = 1.;
  double C_Hex2 = 1.;
  if (LambdaAl2>0.) {
    double logC_Q1 = logQ(log(LambdaAl2));
    C_Q1           = exp(logC_Q1);
    C_Oct1         = Yagi14_fit_Coct(C_Q1);
    C_Hex1         = Yagi14_fit_Chex(C_Q1);
  }
  if (LambdaBl2>0.) {
    double logC_Q2 = logQ(log(LambdaBl2));
    C_Q2           = exp(logC_Q2);
    C_Oct2         = Yagi14_fit_Coct(C_Q2);
    C_Hex2         = Yagi14_fit_Chex(C_Q2);
  }

  /* Default settings for NQC */
  // NOTE: The defaults are different from v0.0 and v1.0
  if (STREQUAL(par_get_s("nqc"),"auto")) {
    if (usetidal) {
      par_set_s("nqc_coefs_flx","none");
      par_set_s("nqc_coefs_hlm","none");
    } else {
      if (usespins) {
	//par_set_s("nqc_coefs_flx","none");
	par_set_s("nqc_coefs_flx", "nrfit_spin202002");
	par_set_s("nqc_coefs_hlm","compute");
      } else {
	par_set_s("nqc_coefs_flx","nrfit_nospin201602");
	par_set_s("nqc_coefs_hlm","nrfit_nospin201602");
      }
    }
  } 
  
  /* Check the NQC option is known */
  for (j = 0; j < 3; j++)
    if (STREQUAL(par_get_s("nqc_coefs_flx"),nqc_flx_opt[j]))
      break;
  if(j>2) errorexit("unknown option for nqc_coefs_flx");
  for (j = 0; j < 4; j++)
    if (STREQUAL(par_get_s("nqc_coefs_hlm"),nqc_hlm_opt[j]))
      break;
  if(j>3) errorexit("unknown option for nqc_coefs_hlm");
  
  /* Set aux parfiles in database */
  
  par_set_d("nu",nu);
  par_set_d("X1",X1);
  par_set_d("X2",X2);

  par_set_d("S1",S1);
  par_set_d("S2",S2);
  par_set_d("S",S);
  par_set_d("Sstar",Sstar);
  par_set_d("a1",a1);
  par_set_d("a2",a2);
  par_set_d("aK",aK);
  par_set_d("aK2",aK2);

  par_set_d("LambdaAl3",LambdaAl3);
  par_set_d("LambdaAl4",LambdaAl4);
  par_set_d("LambdaBl3",LambdaBl3);
  par_set_d("LambdaBl4",LambdaBl4);

  par_set_d("SigmaAl2",SigmaAl2);
  par_set_d("SigmaBl2",SigmaBl2);

  par_set_d("kappaAl2", kapA2 );
  par_set_d("kappaAl3", kapA3 );
  par_set_d("kappaAl4", kapA4 );
  par_set_d("kappaBl2", kapB2 );
  par_set_d("kappaBl3", kapB3 );
  par_set_d("kappaBl4", kapB4 );
  par_set_d("kappaTl2", kapT2 );
  par_set_d("kappaTl3", kapT3 );
  par_set_d("kappaTl4", kapT4 );

  par_set_d("kappajAl2", kapA2j );
  par_set_d("kappajBl2", kapB2j );
  par_set_d("kappajTl2", kapT2j );

  par_set_d("bar_alph2_1",bar_alph2_1);
  par_set_d("bar_alph2_2",bar_alph2_2);
  par_set_d("bar_alph3_1",bar_alph3_1);
  par_set_d("bar_alph3_2",bar_alph3_2);
  par_set_d("bar_alph2j_1",bar_alph2j_1);

  par_set_d("khatAl2",khatA_2);
  par_set_d("khatBl2",khatB_2);

  par_set_d("C_Q1",C_Q1);
  par_set_d("C_Q2",C_Q2);

  par_set_d("C_Oct1",C_Oct1);
  par_set_d("C_Oct2",C_Oct2);
  par_set_d("C_Hex1",C_Hex1);
  par_set_d("C_Hex2",C_Hex2);
  
  /** Set more as needed ... */
  double a6c = 0.;
  if (EOBPars->use_flm == USEFLM_HM) {
    /* Higher modes */
    a6c = eob_a6c_fit_HM(nu);
  } else {
    a6c = eob_a6c_fit(nu);
  }
  par_set_d("a6c", a6c);

  double c3 = 0.;
  if(usetidal) c3 = 0.0;
  else if (EOBPars->use_flm == USEFLM_HM) {
    /* Higher modes */
    c3 = eob_c3_fit_HM(nu,a1,a2);
  } else {
    c3 = eob_c3_fit_global(nu,a1,a2);
  }
  par_set_d("cN3LO", c3 );

  double dt = par_get_d("dt");
  if (par_get_i("use_geometric_units")) {
    /* input given in geometric units, 
       rescale to geometric units and mass rescaled quantities
       compute r0 from the initial GW frequency in geometric units and mass rescaled 
       reset sample rate using dt
    */
    if (VERBOSE) printf("Assume geometric units for pars values\n");
    par_set_d("r0", pow(fmin*Pi, -2./3.) );
    par_set_d("srate", 1./dt );
    par_set_d("distance", 1. );
    par_set_d("M", 1. );
    /* par_set_d("dt_interp", dt ); */ /* Do not reset interp_dt, just use parfile value */
  } else {
    /* input given in physical units, 
       rescale to geometric units and mass rescaled quantities
       compute r0 from the initial GW frequency in Hz 
    */
    if (VERBOSE) printf("Assume physical units for pars values\n");
    /* Set interpolation dt */
    dt = 1./par_get_d("srate_interp");
    dt = time_units_conversion(M, dt);
    par_set_d("dt_interp", dt );
    /* Set dt */
    par_set_d("r0",  radius0(M, fmin) );
    dt = 1./par_get_d("srate");
    dt = time_units_conversion(M, dt);
    par_set_d("dt", dt);
    if (VERBOSE) PRFORMd("dt",dt);
    if (VERBOSE)
      if (par_get_i("interp_uniform_grid"))
	PRFORMd("dt_interp",par_get_d("dt_interp"));
  }

  /* Function pointers */
  
  /** Set f_lm fun pointer */
  if ((STREQUAL(par_get_s("use_flm"),"HM"))) {
    /* Higher modes */
    eob_wav_hlmNewt  = &eob_wav_hlmNewt_HM;
    eob_wav_flm      = &eob_wav_flm_HM;
    eob_wav_flm_s    = &eob_wav_flm_s_HM;
    eob_wav_deltalm  = &eob_wav_deltalm_HM;
    eob_wav_ringdown = &eob_wav_ringdown_HM;
  }  else if ((STREQUAL(par_get_s("use_flm"),"SSLO"))) {
    /* eob_wav_flm_s = &eob_wav_flm_s_old; */
    eob_wav_hlmNewt  = &eob_wav_hlmNewt_v1;
    eob_wav_flm      = &eob_wav_flm_v1;
    eob_wav_flm_s    = &eob_wav_flm_s_SSLO;
    eob_wav_deltalm  = &eob_wav_deltalm_v1;
    eob_wav_ringdown = &eob_wav_ringdown_v1;
  } else {
    eob_wav_hlmNewt  = &eob_wav_hlmNewt_v1;
    eob_wav_flm      = &eob_wav_flm_v1;
    eob_wav_flm_s    = &eob_wav_flm_s_SSNLO;
    eob_wav_deltalm  = &eob_wav_deltalm_v1;
    eob_wav_ringdown = &eob_wav_ringdown_v1;
  }

  /** Set rc fun pointer */
  if ((STREQUAL(par_get_s("centrifugal_radius"),"LO"))) {
    eob_dyn_s_get_rc = &eob_dyn_s_get_rc_LO;
  } else if ((STREQUAL(par_get_s("centrifugal_radius"),"NLO"))) {
    eob_dyn_s_get_rc = &eob_dyn_s_get_rc_NLO;
  } else if ((STREQUAL(par_get_s("centrifugal_radius"),"NNLO"))) {
    eob_dyn_s_get_rc = &eob_dyn_s_get_rc_NNLO;
  } else if ((STREQUAL(par_get_s("centrifugal_radius"),"NNLOS4"))) {
    eob_dyn_s_get_rc = &eob_dyn_s_get_rc_NNLO_S4;
  } else if ((STREQUAL(par_get_s("centrifugal_radius"),"NOSPIN"))) {
    eob_dyn_s_get_rc = &eob_dyn_s_get_rc_NOSPIN;
  } else if ((STREQUAL(par_get_s("centrifugal_radius"),"NOTIDES"))) {
    eob_dyn_s_get_rc = &eob_dyn_s_get_rc_NOTIDES;
  }
}
#endif

#if (DEBUG_THIS_FILE)

/* test 
   gcc TEOBResumSPars.c -lconfig -o testpars.x */
int main (int argc, char* argv[])
{
  int idx,n;
  /* init db */
  par_db_init ();
  /* parse default params from file */
  par_db_default ();
  /* display default params */
  par_db_screen ();
  /* test 'get' */
  printf("get-test: size = %d\n",par_get_i("size"));
  printf("get-test: M = %e\n",par_get_d("M"));
  int *i = par_get_arrayi ("output_lm", &n);
  printf("get-test: array_length n=%d\n",n);
  for (idx=0; idx<n; idx++) printf("get-test: array[%d]=%d\n",idx,i[idx]);
  /* test 'set' */
  par_set_i("size",10);
  par_set_d("M",100);
  printf("set-test: size = %d (=10)\n",par_get_i("size"));
  printf("set-test: M = %e (=100)\n",par_get_d("M"));
  n=6;
  int j[] = {1,2,3,4,5,6};
  par_set_arrayi ("output_lm", j, n);
  int *jj = par_get_arrayi ("output_lm", &n);
  printf("get-test: array_length n=%d\n",n);
  for (idx=0; idx<n; idx++) printf("get-test: array[%d]=%d\n",idx,jj[idx]);
  /* parse other params */
  par_file_parse_merge ("../par/test.par");
  /* more operations */
  printf("get-test: dt = %e\n",par_get_d("dt"));
  printf("get-test: size = %d\n",par_get_i("size"));
  par_set_i("size",20);
  printf("set-test: size = %d (=20)\n",par_get_i("size"));
  /* display new db */
  par_db_screen ();  
  /* write out db */
  par_db_write_file("../par/test1.par");
  /* free db */
  par_db_free ();
  return OK;
}

#endif
