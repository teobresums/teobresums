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
 * 
 * To add a parameter:
 * - add the relative variable in the EOBParameters structure
 * - add its default value in EOBParameters_defaults() (if taken as input) ...
 * - ... or in eob_set_params() (if computed from other input parameters)
 * - (if taken as input) make sure it is parsed by EOBParameterse_parse_file()
 *   and in the analogue routine of the python wrapper
 *
 * History
 * - SB 09/2021 simplified the logic, added minimal custum routines to handle parfile input. Removed libconfig use.
 * - SB 11/2019 Added the par structure type EOBParameters and routines to work with it.
 * - v0.0 and v1.0 worked only with the parameter db and the parfile, using libconfig wrappers.
 *
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

/** Global vars for function pointers
    (extern in the header) */
void (*eob_wav_hlmNewt)();
void (*eob_wav_flm)();
void (*eob_wav_flm_s)();
void (*eob_wav_deltalm)();
void (*eob_wav_hlmNQC_find_a1a2a3_mrg)();
void (*eob_wav_ringdown)();
void (*eob_dyn_s_get_rc)();
int (*p_eob_spin_dyn_rhs)();

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
  if (eobp->use_mode_lm_inertial) free (eobp->use_mode_lm_inertial);
  if (eobp->use_mode_lm_nqc) free (eobp->use_mode_lm_nqc);
  if (eobp->kpostpeak) free (eobp->kpostpeak);
  if (eobp->knqcpeak22) free (eobp->knqcpeak22);
  if (eobp->output_lm) free (eobp->output_lm);
  if (eobp->freqs) free(eobp->freqs);
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

  eobp->LambdaAl2 = 0.; // Tidal gravitoelectric parameter Lambda for star A ell=2
  eobp->LambdaBl2 = 0.;
  eobp->LambdaAl3 = 0.; 
  eobp->LambdaAl4 = 0.; 
  eobp->LambdaAl5 = 0.; 
  eobp->LambdaAl6 = 0.; 
  eobp->LambdaAl7 = 0.; 
  eobp->LambdaAl8 = 0.; 
  eobp->LambdaBl3 = 0.; 
  eobp->LambdaBl4 = 0.; 
  eobp->LambdaBl5 = 0.; 
  eobp->LambdaBl6 = 0.; 
  eobp->LambdaBl7 = 0.; 
  eobp->LambdaBl8 = 0.; 
  eobp->SigmaAl2  = 0.; // Tidal gravitomagnetic parameter Sigma for star A ell=2
  eobp->SigmaBl2  = 0.;
  eobp->use_lambda234_fits = Lambda234_fits_NO;
  eobp->pGSF_tidal = 4.0;// p-power in GSF tidal potential model

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


  /* Modes to use in the intertial frame to compute hpc */
  int hlm_inertial[] = {1};
  eobp->use_mode_lm_inertial_size = 1;
  eobp->use_mode_lm_inertial = malloc (eobp->use_mode_lm_inertial_size * sizeof(int) );
  memcpy(eobp->use_mode_lm_inertial, hlm_inertial, eobp->use_mode_lm_inertial_size * sizeof(int));

  // TODO: intersect with use_mode_lm so to compute NQCs only for active modes
  int hlm_nqc[] = {0,1,3,4,6,7,8,13};      //indexes of multipoles to use
  eobp->use_mode_lm_nqc_size = 8;
  eobp->use_mode_lm_nqc = malloc (eobp->use_mode_lm_nqc_size * sizeof(int) );
  memcpy(eobp->use_mode_lm_nqc, hlm_nqc, eobp->use_mode_lm_nqc_size * sizeof(int));

  int kpostpeak[] = {1,3,6,7,8};      //indexes of multipoles to use
  eobp->kpostpeak_size = 5;
  eobp->kpostpeak = malloc (eobp->kpostpeak_size * sizeof(int) );
  memcpy(eobp->kpostpeak, kpostpeak, eobp->kpostpeak_size * sizeof(int));
  
  int knqcpeak22[]      = {0};      //indexes of multipoles to use
  eobp->knqcpeak22_size = 1;
  eobp->knqcpeak22      = malloc (eobp->knqcpeak22_size * sizeof(int) );
  memcpy(eobp->knqcpeak22, knqcpeak22, eobp->knqcpeak22_size * sizeof(int));

  /* FD options */
  
  eobp->tc = 0;
  eobp->time_shift_FD=1;
  eobp->df = 1.;
  
  eobp->interp_freqs=0;
  double fr[] = {30.};      //indexes of multipoles to use
  eobp->freqs_size = 1;
  eobp->freqs = malloc (eobp->freqs_size * sizeof(double));
  memcpy(eobp->freqs, fr, eobp->freqs_size * sizeof(double));

  /* EOB Settings */
  
  eobp->postadiabatic_dynamics=1;
  eobp->postadiabatic_dynamics_N=8;      // post-adiabatic order
  eobp->postadiabatic_dynamics_size=800; // grid size 
  eobp->postadiabatic_dynamics_rmin=14.; // minimum radius (end of PA dynamics)
  eobp->postadiabatic_dynamics_stop=0;   // stop after post-adiabatic dynamics 
  eobp->postadiabatic_dynamics_dr=0.1;            // PA step dr = 0.1 


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
  eobp->output_hpc= 0; // output h+,hx
  eobp->output_multipoles= 0; // output multipoles
  eobp->output_dynamics=0; // output dynamics
  eobp->output_nqc=0; // output NQC waveform
  eobp->output_nqc_coefs=0; // output multipolar NQC coefs (if determined)
  eobp->output_ringdown=0; // output ringdown waveform
  
  int klm[] = {1};      //indexes of multipoles to ouput
  eobp->output_lm_size = 1;
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
  eobp->ode_tmax=1e12; // max integration time
  eobp->ode_stop_radius  =1.; // stop ODE integration at this radius (if > 0)
  eobp->ode_stop_afterNdt=4;  // stop ODE N iters after the Omega peak
  eobp->ode_stop_after_peak=0;

  /* Spin dynamics */
  
  eobp->spin_dyn_size      = 500;// size for spin dynamics arrays
  eobp->spin_odes_omg_stop = -1; // stop ODE integration at this Momega set by NR BBH mrg freq fit 
  eobp->spin_odes_t_stop   = -1;
  eobp->spin_odes_dt       = 1;  // timestep for spin dynamics
  eobp->spin_interp_domain = 1;  // FD or TD interpolation
  eobp->chi1x = eobp->chi1y = eobp->chi1z = 0.;
  eobp->chi2x = eobp->chi2y = eobp->chi2z = 0.;
  eobp->spin_flx             = SPIN_FLX_PN;
  eobp->ringdown_eulerangles = RD_EULERANGLES_QNMs;

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
  eobp->kapA5= 0. ; //
  eobp->kapA6= 0. ; //
  eobp->kapA7= 0. ; //
  eobp->kapA8= 0. ; //
  eobp->kapB2= 0. ; //
  eobp->kapB3= 0. ; //
  eobp->kapB4= 0. ; //
  eobp->kapB5= 0. ; //
  eobp->kapB6= 0. ; //
  eobp->kapB7= 0. ; //
  eobp->kapB8= 0. ; //
  eobp->kapT2= 0. ; //
  eobp->kapT3= 0. ; //
  eobp->kapT4= 0. ; //
  eobp->kapT5= 0. ; //
  eobp->kapT6= 0. ; //
  eobp->kapT7= 0. ; //
  eobp->kapT8= 0. ; //

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

  // f-mode frequencies
  for (int l=0; l<6; l++) {
    eobp->bomgfA[l] = eobp->bomgfB[l] = 0.;
  }
	  
  eobp->Mbhf= 0. ; // final BH mass
  eobp->abhf= 0. ; // final BH spin
  
  eobp->rLR= 0. ; // radius of light-ring
  eobp->rLSO= 0. ; // radius of last stable orbit 

  
  /* Choose the default for the binary type */
  eobp->binary = BINARY_BBH;
  eobp->use_tidal = TIDES_OFF ; // index for tidal modus
  eobp->use_tidal_gravitomagnetic = TIDES_GM_OFF ; // index for gravitomagnetic tide
  eobp->use_tidal_fmode_model = 0; // do not 
  
  if (choose == BINARY_BBH) {
    eobp->binary=BINARY_BBH;
    eobp->centrifugal_radius = CENTRAD_NLO;
    eobp->use_flm = USEFLM_HM;
    eobp->use_tidal = TIDES_OFF;
    eobp->use_tidal_gravitomagnetic = TIDES_GM_OFF;
    eobp->use_a6c_fits   = a6c_fits_HM_2023;
    eobp->use_cN3LO_fits = cN3LO_fits_HM_2023_432; 

    eobp->nqc = NQC_AUTO; // {"no", "auto", "manual"}
    eobp->nqc_coefs_flx = NQC_FLX_NRFIT_SPIN_202002; // {"none", "nrfit_nospin20160209", "nrfit_spin20202","fromfile"}
    eobp->nqc_coefs_hlm = NQC_HLM_COMPUTE; // {"compute", "none", "nrfit_nospin20160209", "nrfit_spin20202", "fromfile"}

  } else if (choose == BINARY_BNS) {
    eobp->binary=BINARY_BNS;
    eobp->use_tidal = TIDES_TEOBRESUM3;
    eobp->use_tidal_gravitomagnetic = TIDES_GM_PN;
    eobp->pGSF_tidal = 4.0;
    eobp->use_lambda234_fits = Lambda234_fits_YAGI13;
    eobp->use_a6c_fits   = a6c_fits_V0;
    eobp->use_cN3LO_fits = cN3LO_fits_NO;

    eobp->centrifugal_radius = CENTRAD_NNLO;
    eobp->use_flm = USEFLM_SSNLO;
    eobp->nqc = NQC_NO; // {"no", "auto", "manual"}
    eobp->nqc_coefs_flx = NQC_FLX_NONE; // {"none", "nrfit_nospin20160209", "nrfit_spin20202", "fromfile"}
    eobp->nqc_coefs_hlm = NQC_HLM_NONE;

  } else if (choose == BINARY_BHNS) {
    eobp->binary=BINARY_BHNS;
    eobp->centrifugal_radius = CENTRAD_NLO;
    eobp->use_flm = USEFLM_HM;
    eobp->nqc = NQC_AUTO; // {"no", "auto", "manual"}                                                                                                                                                                                         
    eobp->nqc_coefs_flx = NQC_FLX_NRFIT_SPIN_202002; // {"none", "nrfit_nospin20160209", "nrfit_spin20202","fromfile"}                                                                                                                        
    eobp->nqc_coefs_hlm = NQC_HLM_COMPUTE; // {"compute", "none", "nrfit_nospin20160209", "nrfit_spin20202", "fromfile"}                                                                                                                      
    eobp->use_a6c_fits   = a6c_fits_HM_2023;
    eobp->use_cN3LO_fits = cN3LO_fits_HM_2023_432; 

    eobp->use_tidal = TIDES_TEOBRESUM;
    eobp->use_tidal_gravitomagnetic = TIDES_GM_OFF;//TIDES_GM_PN;
    eobp->use_lambda234_fits = Lambda234_fits_YAGI13;

  }


  else errorexit("unknown default parameter choice.");

}

/*
 * main routine to set parameters for the run
 */

void eob_set_params(int default_choice, int firstcall)
{
  
  /* Set intrinsic parameters as given by user */
  
  const double M    =  EOBPars->M;
  const double fmin = EOBPars->initial_frequency;
  double q          =  EOBPars->q;

  if(q < 1.){
    /* Ensure that the swap is performed correctly 
       According to LAL conventions, if m1<->m2 the x axis is flipped.
       Therefore, the waveform has to remain identical when:
        - label_1<->label_2;
        - coalescence_angle -> coalescence_angle + Pi;
        - In plane spins are rotated by Pi
       Here we enforce this convention.
    */
    q  =  1./q;
    EOBPars->q = q;
    SWAPTRS(EOBPars->chi1z, EOBPars->chi2z);
    SWAPTRS(EOBPars->chi1x, EOBPars->chi2x);
    SWAPTRS(EOBPars->chi1y, EOBPars->chi2y);
    SWAPTRS(EOBPars->chi1,  EOBPars->chi2);
    EOBPars->chi1x *= -1; EOBPars->chi2x *= -1;
    EOBPars->chi1y *= -1; EOBPars->chi2y *= -1;
    SWAPTRS(EOBPars->LambdaAl2, EOBPars->LambdaBl2);
    if (VERBOSE) printf("WARNING: q<1, swapping bodies!\n");
    EOBPars->coalescence_angle -= Pi;
  }

  /* Check: if q is closer to 1 than 1e-8, then q=1 to avoid floating points issues */
  if (DEQUAL(q, 1., 1e-8)){
    if(VERBOSE) printf("WARNING: manually fixing q = 1\n");
    q          = 1.;
    EOBPars->q = 1.;
  }

  
  /* Check: if in-plane components of spins are < 1e-4, then spin aligned, else spin precessing */
  /* Note: the "NOSPIN" option is deprecated, and currently never used */

  double chip_1 = sqrt(SQ(EOBPars->chi1x) + SQ(EOBPars->chi1y));
  double chip_2 = sqrt(SQ(EOBPars->chi2x) + SQ(EOBPars->chi2y));
  if (chip_1 + chip_2 > 1e-4){
    EOBPars->use_spins = MODE_SPINS_GENERIC;
    EOBPars->chi1 = EOBPars->chi1z;
    EOBPars->chi2 = EOBPars->chi2z;
    if (EOBPars->domain == DOMAIN_TD)
      EOBPars->interp_uniform_grid = 1; // for TD twist we require interpolation
  } else {
    if (DUNEQUAL(EOBPars->chi1z, EOBPars->chi1, 1e-8) || DUNEQUAL(EOBPars->chi2z, EOBPars->chi2, 1e-8)){
      /* User specified only chi_{12} or chi_{12}z. The unspecified value is 0 by default. 
         Two possibilities:
         * chi_iz nonzero, chi_i zero   : spin aligned limit of precessing --> overwrite EOBPars->chi{12}
         * chi_iz zero,    chi_i nonzero: spin aligned                     --> do not overwrite EOBPars->chi{12}
      */
      if DISZERO(EOBPars->chi1, 1e-8) EOBPars->chi1 = EOBPars->chi1z;
      if DISZERO(EOBPars->chi2, 1e-8) EOBPars->chi2 = EOBPars->chi2z;
    }
    EOBPars->use_spins = MODE_SPINS_ALIGNED;
  }

  EOBPars->nu = q_to_nu(q);
  EOBPars->X1 = nu_to_X1(EOBPars->nu);
  EOBPars->X2 = 1. -  EOBPars->X1;
  const double XA = EOBPars->X1; /* tidal part used different notation, used here for simplicity */
  const double XB = EOBPars->X2;
  
  const int usespins = EOBPars->use_spins;
  const int usetidal = EOBPars->use_tidal;
  const int usetidalGM =  EOBPars->use_tidal_gravitomagnetic;

  const double chi1 = EOBPars->chi1;
  const double chi2 = EOBPars->chi2;
  set_spin_vars(XA,XB, EOBPars->chi1,EOBPars->chi2, 
		&EOBPars->S1, &EOBPars->S2,
		&EOBPars->a1, &EOBPars->a2,
		&EOBPars->aK, &EOBPars->aK2,
		&EOBPars->S, &EOBPars->Sstar);

  if (usetidal) {
    
    /* Set the tidal parameters */
    
    if (EOBPars->use_lambda234_fits == Lambda234_fits_YAGI13) {
      EOBPars->LambdaAl3 = Yagi13_fit_barlamdel(EOBPars->LambdaAl2, 3);
      EOBPars->LambdaBl3 = Yagi13_fit_barlamdel(EOBPars->LambdaBl2, 3);
      EOBPars->LambdaAl4 = Yagi13_fit_barlamdel(EOBPars->LambdaAl2, 4);
      EOBPars->LambdaBl4 = Yagi13_fit_barlamdel(EOBPars->LambdaBl2, 4);
    } else if (EOBPars->use_lambda234_fits == Lambda234_fits_GODZIEBA20) {
      EOBPars->LambdaAl3 = Godzieba20_fit_barlamdel(EOBPars->LambdaAl2, 3);
      EOBPars->LambdaBl3 = Godzieba20_fit_barlamdel(EOBPars->LambdaBl2, 3);
      EOBPars->LambdaAl4 = Godzieba20_fit_barlamdel(EOBPars->LambdaAl2, 4);
      EOBPars->LambdaBl4 = Godzieba20_fit_barlamdel(EOBPars->LambdaBl2, 4);
    } else if (EOBPars->use_lambda234_fits == Lambda2345678_fits_GODZIEBA20) {
      EOBPars->LambdaAl3 = Godzieba20_fit_barlamdel(EOBPars->LambdaAl2, 3);
      EOBPars->LambdaBl3 = Godzieba20_fit_barlamdel(EOBPars->LambdaBl2, 3);
      EOBPars->LambdaAl4 = Godzieba20_fit_barlamdel(EOBPars->LambdaAl2, 4);
      EOBPars->LambdaBl4 = Godzieba20_fit_barlamdel(EOBPars->LambdaBl2, 4);
      EOBPars->LambdaAl5 = Godzieba20_fit_barlamdel(EOBPars->LambdaAl2, 5);
      EOBPars->LambdaBl5 = Godzieba20_fit_barlamdel(EOBPars->LambdaBl2, 5);
      EOBPars->LambdaAl6 = Godzieba20_fit_barlamdel(EOBPars->LambdaAl2, 6);
      EOBPars->LambdaBl6 = Godzieba20_fit_barlamdel(EOBPars->LambdaBl2, 6);
      EOBPars->LambdaAl7 = Godzieba20_fit_barlamdel(EOBPars->LambdaAl2, 7);
      EOBPars->LambdaBl7 = Godzieba20_fit_barlamdel(EOBPars->LambdaBl2, 7);
      EOBPars->LambdaAl8 = Godzieba20_fit_barlamdel(EOBPars->LambdaAl2, 8);
      EOBPars->LambdaBl8 = Godzieba20_fit_barlamdel(EOBPars->LambdaBl2, 8);
    }
    if(EOBPars->use_tidal_gravitomagnetic){
      EOBPars->SigmaAl2 = JFAPG_fit_Sigma_Irrotational(EOBPars->LambdaAl2);
      EOBPars->SigmaBl2 = JFAPG_fit_Sigma_Irrotational(EOBPars->LambdaBl2);
    }

    /* Tidal coupling constants */
    tidal_kappa_of_Lambda(q, XA, XB, EOBPars->LambdaAl2,EOBPars->LambdaBl2, 2,  &(EOBPars->kapA2), &(EOBPars->kapB2));
    tidal_kappa_of_Lambda(q, XA, XB, EOBPars->LambdaAl3,EOBPars->LambdaBl3, 3,  &(EOBPars->kapA3), &(EOBPars->kapB3));
    tidal_kappa_of_Lambda(q, XA, XB, EOBPars->LambdaAl4,EOBPars->LambdaBl4, 4,  &(EOBPars->kapA4), &(EOBPars->kapB4));
    tidal_kappa_of_Lambda(q, XA, XB, EOBPars->LambdaAl5,EOBPars->LambdaBl5, 5,  &(EOBPars->kapA5), &(EOBPars->kapB5));
    tidal_kappa_of_Lambda(q, XA, XB, EOBPars->LambdaAl6,EOBPars->LambdaBl6, 6,  &(EOBPars->kapA6), &(EOBPars->kapB6));
    tidal_kappa_of_Lambda(q, XA, XB, EOBPars->LambdaAl7,EOBPars->LambdaBl7, 7,  &(EOBPars->kapA7), &(EOBPars->kapB7));
    tidal_kappa_of_Lambda(q, XA, XB, EOBPars->LambdaAl8,EOBPars->LambdaBl8, 8,  &(EOBPars->kapA8), &(EOBPars->kapB8));
    
    /* gravitomagnetic tidal coupling constants el = 2 only */
    EOBPars->japA2 = 24.   * EOBPars->SigmaAl2 * XA*XA*XA*XA*XA / q;
    EOBPars->japB2 = 24.   * EOBPars->SigmaBl2 * XB*XB*XB*XB*XB * q;
    
    EOBPars->kapT2 = EOBPars->kapA2 + EOBPars->kapB2;
    EOBPars->kapT3 = EOBPars->kapA3 + EOBPars->kapB3;
    EOBPars->kapT4 = EOBPars->kapA4 + EOBPars->kapB4;
    EOBPars->kapT5 = EOBPars->kapA5 + EOBPars->kapB5;
    EOBPars->kapT6 = EOBPars->kapA6 + EOBPars->kapB6;
    EOBPars->kapT7 = EOBPars->kapA7 + EOBPars->kapB7;
    EOBPars->kapT8 = EOBPars->kapA8 + EOBPars->kapB8;
    
    EOBPars->japT2 = EOBPars->japA2 + EOBPars->japB2;
    
    if (!(EOBPars->kapT2 > 0.)) errorexit("kappaT2 must be >0");
    if (!(EOBPars->kapT3 > 0.)) errorexit("kappaT3 must be >0");
    if (!(EOBPars->kapT4 > 0.)) errorexit("kappaT4 must be >0");
    
    /* Tidal coefficients cons dynamics
       \bar{\alpha}_n^{(\ell)}, Eq.(37) of Damour&Nagar, PRD 81, 084016 (2010) */
    EOBPars->bar_alph2_1 = (5./2.*XA*EOBPars->kapA2 + 5./2.*XB*EOBPars->kapB2)/EOBPars->kapT2;
    EOBPars->bar_alph2_2 = ((3.+XA/8.+ 337./28.*XA*XA)*EOBPars->kapA2 + (3.+XB/8.+ 337./28.*XB*XB)*EOBPars->kapB2)/EOBPars->kapT2;
    EOBPars->bar_alph3_1 = ((-2.+15./2.*XA)*EOBPars->kapA3 + (-2.+15./2.*XB)*EOBPars->kapB3)/EOBPars->kapT3;
    EOBPars->bar_alph3_2 = ((8./3.-311./24.*XA+110./3.*XA*XA)*EOBPars->kapA3 + (8./3.-311./24.*XB+110./3.*XB*XB)*EOBPars->kapB3)/EOBPars->kapT3;
    /* Gravitomagnetic term, see Eq.(6.27) of Bini-Damour-Faye 2012 */
    EOBPars->bar_alph2j_1 = ( EOBPars->japA2*(1. + (11./6.)*XA + XA*XA) + EOBPars->japB2*(1. + (11./6.)*XB + XB*XB) )/EOBPars->japT2;
    
    
    /* Self-spin coefficients */
    EOBPars->C_Q1   = 1.;
    EOBPars->C_Q2   = 1.;
    EOBPars->C_Oct1 = 1.;
    EOBPars->C_Oct2 = 1.;
    EOBPars->C_Hex1 = 1.;
    EOBPars->C_Hex2 = 1.;
    if (EOBPars->LambdaAl2>0.) {
      double logC_Q1 = YagiYunes13_fit_logQ(log(EOBPars->LambdaAl2));
      EOBPars->C_Q1           = exp(logC_Q1);
      EOBPars->C_Oct1         = Yagi14_fit_Coct(EOBPars->C_Q1);
      EOBPars->C_Hex1         = Yagi14_fit_Chex(EOBPars->C_Q1);
    }
    if (EOBPars->LambdaBl2>0.) {
      double logC_Q2 = YagiYunes13_fit_logQ(log(EOBPars->LambdaBl2));
      EOBPars->C_Q2           = exp(logC_Q2);
      EOBPars->C_Oct2         = Yagi14_fit_Coct(EOBPars->C_Q2);
      EOBPars->C_Hex2         = Yagi14_fit_Chex(EOBPars->C_Q2);
    }
    
    /* Set the f-mode frequencies from fits, if needed */
    if (EOBPars->use_tidal_fmode_model) {
            
      const int lmax = 4;
      double LamAl[] = {0,0,EOBPars->LambdaAl2,EOBPars->LambdaAl3,EOBPars->LambdaAl4};
      double LamBl[] = {0,0,EOBPars->LambdaBl2,EOBPars->LambdaBl3,EOBPars->LambdaBl4};
      
      for (int l=2; l<=lmax; l++) {
	if (LamAl[l] > 0) {
	  EOBPars->bomgfA[l] = Chang14_fit_omegaf(LamAl[l], l);	 
	  if (EOBPars->bomgfA[l]<=0.)
	    errorexit("f-mode frequency of star A cannot be zero or negative");
	  EOBPars->bomgfA[l] /= XA;
	}
	if (LamBl[l] > 0) {
	  EOBPars->bomgfB[l] = Chang14_fit_omegaf(LamBl[l], l);
	  if (EOBPars->bomgfB[l]<=0.)
	    errorexit("f-mode frequency of star B cannot be zero or negative");
	  EOBPars->bomgfB[l] /= XB;
	}
      } 
      
    }  /* EOBPars->use_tidal_fmode_model */

  }  /* use_tidal */
  
  /* Default settings for NQC */
  // NOTE: The defaults are different from v0.0 and v1.0
  if (EOBPars->nqc == NQC_AUTO) {
    if (EOBPars->binary == BINARY_BNS) {
        EOBPars->nqc_coefs_flx = NQC_FLX_NONE;
        EOBPars->nqc_coefs_hlm = NQC_HLM_NONE;
    } else {
      if (usespins) {
        EOBPars->nqc_coefs_flx = NQC_FLX_NRFIT_SPIN_202002;
        EOBPars->nqc_coefs_hlm = NQC_HLM_COMPUTE;
      } else {
        EOBPars->nqc_coefs_flx = NQC_FLX_NRFIT_NOSPIN_201602;
        EOBPars->nqc_coefs_hlm = NQC_HLM_NRFIT_NOSPIN_201602;
      }
    }
  } 
  
  /** Set more as needed ... */
  EOBPars->a6c = 0.;
  switch(EOBPars->use_a6c_fits)
  {
    case(a6c_fits_HM_2023):
      if (EOBPars->use_flm != USEFLM_HM)
        errorexit("a6c_fits_HM should be used with USEFLM_HM\n.");
      EOBPars->a6c = eob_a6c_fit_HM_2023(EOBPars->nu);
      break;
    case(a6c_fits_HM):
      // check compatibility with FLM
      if (EOBPars->use_flm != USEFLM_HM)
        errorexit("a6c_fits_HM should be used with USEFLM_HM\n.");
      EOBPars->a6c = eob_a6c_fit_HM(EOBPars->nu);
      break;
    case(a6c_fits_V0):
      EOBPars->a6c = eob_a6c_fit(EOBPars->nu);
      break;
    case(a6c_fits_NO):
      break;
    default:
      if (VERBOSE) printf("WARINING: No option specified for a6c.\n");
      break; 
  }

  EOBPars->cN3LO = 0.;
  switch(EOBPars->use_cN3LO_fits)
  {
    case(cN3LO_fits_HM_2023_432):
    case(cN3LO_fits_HM_2023_431):
    case(cN3LO_fits_HM_2023_430):
    case(cN3LO_fits_HM_2023_420):
      if (EOBPars->use_flm != USEFLM_HM)
        errorexit("cN3LO_fits_HM_2023 should be used with USEFLM_HM\n.");
      EOBPars->cN3LO = eob_c3_fit_HM_2023(EOBPars->nu,EOBPars->a1,EOBPars->a2);
      break;
    case(cN3LO_fits_HM):
      if (EOBPars->use_flm != USEFLM_HM)
        errorexit("cN3LO_fits_HM should be used with USEFLM_HM\n.");
      EOBPars->cN3LO = eob_c3_fit_HM(EOBPars->nu,EOBPars->a1,EOBPars->a2);
      break;
    case(cN3LO_fits_V0):
      EOBPars->cN3LO = eob_c3_fit_global(EOBPars->nu,EOBPars->a1,EOBPars->a2);
      break;
    case(cN3LO_fits_NO):
      break;
    default:
      if (VERBOSE) printf("WARINING: No option specified for cN3LO.\n");
      break; 
  }

  double dt = EOBPars->dt;
  if (EOBPars->use_geometric_units) {
    /* input given in geometric units, 
       rescale to geometric units and mass rescaled quantities
       compute r0 from the initial GW frequency in geometric units and mass rescaled 
       reset sample rate using dt
    */
    if (VERBOSE) printf("Assume geometric units for pars values\n");
    
    /* set r0 based on fmin if it is not specified, else set fmin according to r0 */
    if (EOBPars->r0 == 0.){
      EOBPars->r0 = pow(fmin*Pi, -2./3.);    
    } else {
      EOBPars->initial_frequency = pow(EOBPars->r0, -1.5)/Pi;
    }
    
    EOBPars->srate = 1./dt;
    EOBPars->distance = 1.;
    EOBPars->M = 1.;
  } else {
    /* input given in physical units, 
       rescale to geometric units and mass rescaled quantities
       compute r0 from the initial GW frequency in Hz 
    */
    if (VERBOSE) printf("Assume physical units for pars values\n");
    
    /* set r0 based on fmin if it is not specified, else set fmin according to r0 */
    if (EOBPars->r0 == 0.){
      EOBPars->r0 = radius0(M, fmin);
    } else {
      EOBPars->initial_frequency = pow(EOBPars->r0, -1.5)/(Pi*M*MSUN_S);
    }
    
    /* Set interpolation dt */
    dt = 1./EOBPars->srate_interp;
    dt = time_units_conversion(M, dt);
    EOBPars->dt_interp = dt;
    /* Set dt */
    dt = 1./EOBPars->srate;
    dt = time_units_conversion(M, dt);
    EOBPars->dt = dt;
    if (VERBOSE) PRFORMd("dt",EOBPars->dt);
    if (VERBOSE)
      if (EOBPars->interp_uniform_grid)
	PRFORMd("dt_interp",EOBPars->dt_interp);
  }

  /* Function pointers */
  
  /** Set f_lm fun pointer */
  if (EOBPars->use_flm == USEFLM_HM) {
    eob_wav_hlmNewt  = &eob_wav_hlmNewt_HM;
    eob_wav_flm      = &eob_wav_flm_HM;
    eob_wav_flm_s    = &eob_wav_flm_s_HM;
    eob_wav_deltalm  = &eob_wav_deltalm_HM;
    eob_wav_hlmNQC_find_a1a2a3_mrg = &eob_wav_hlmNQC_find_a1a2a3_mrg_HM;
    eob_wav_ringdown = &eob_wav_ringdown_HM; 
  } else if (EOBPars->use_flm == USEFLM_SSLO) {
    /* eob_wav_flm_s = &eob_wav_flm_s_old; */
    eob_wav_hlmNewt  = &eob_wav_hlmNewt_v1;
    eob_wav_flm      = &eob_wav_flm_v1;
    eob_wav_flm_s    = &eob_wav_flm_s_SSLO;
    eob_wav_deltalm  = &eob_wav_deltalm_v1;
    eob_wav_hlmNQC_find_a1a2a3_mrg = &eob_wav_hlmNQC_find_a1a2a3_mrg_22;
    eob_wav_ringdown = &eob_wav_ringdown_v1;
  } else if (EOBPars->use_flm == USEFLM_SSNLO) {
    eob_wav_hlmNewt  = &eob_wav_hlmNewt_v1;
    eob_wav_flm      = &eob_wav_flm_v1;
    eob_wav_flm_s    = &eob_wav_flm_s_SSNLO;
    eob_wav_deltalm  = &eob_wav_deltalm_v1;
    eob_wav_hlmNQC_find_a1a2a3_mrg = &eob_wav_hlmNQC_find_a1a2a3_mrg_22;
    eob_wav_ringdown = &eob_wav_ringdown_v1;
    /*
      } else if (EOBPars->use_flm == USEFLM_SSNNLO) {
      eob_wav_hlmNewt = &eob_wav_hlmNewt_v1;
      eob_wav_flm     = &eob_wav_flm_v1;
      eob_wav_flm_s   = &eob_wav_flm_s_SSNNLO;
    */
  } else errorexit("unknown option for use_flm");

  /** Set rc fun pointer */
  if (EOBPars->centrifugal_radius == CENTRAD_LO) {
    eob_dyn_s_get_rc = &eob_dyn_s_get_rc_LO;
  } else if (EOBPars->centrifugal_radius == CENTRAD_NLO) {
    eob_dyn_s_get_rc = &eob_dyn_s_get_rc_NLO;
  } else if (EOBPars->centrifugal_radius == CENTRAD_NNLO) {
    eob_dyn_s_get_rc = &eob_dyn_s_get_rc_NNLO;
  } else if (EOBPars->centrifugal_radius == CENTRAD_NNLOS4) {
    eob_dyn_s_get_rc = &eob_dyn_s_get_rc_NNLO_S4;
  } else if (EOBPars->centrifugal_radius == CENTRAD_NOSPIN) {
    eob_dyn_s_get_rc = &eob_dyn_s_get_rc_NOSPIN;
  } else if (EOBPars->centrifugal_radius == CENTRAD_NOTIDES) {
    eob_dyn_s_get_rc = &eob_dyn_s_get_rc_NOTIDES;
  } else errorexit("unknown option for centrifugal_radius");

}

void update_params(int binary)
{
  /* Updated function pointers */
  
  /** Set f_lm fun pointer */
  if (EOBPars->use_flm == USEFLM_HM) {
    eob_wav_hlmNewt  = &eob_wav_hlmNewt_HM;
    eob_wav_flm      = &eob_wav_flm_HM;
    eob_wav_flm_s    = &eob_wav_flm_s_HM;
    eob_wav_deltalm  = &eob_wav_deltalm_HM;
    if(binary==BINARY_BBH){
      eob_wav_hlmNQC_find_a1a2a3_mrg = &eob_wav_hlmNQC_find_a1a2a3_mrg_HM;
      eob_wav_ringdown = &eob_wav_ringdown_HM; 
    }else { // for BHNS cases I & III
      eob_wav_hlmNQC_find_a1a2a3_mrg = &eob_wav_hlmNQC_find_a1a2a3_mrg_BHNS_HM;
      eob_wav_ringdown = &eob_wav_ringdown_bhns;
    }
  }
}

int EOBParameters_parse_commandline(EOBParameters *eobp, int argc, char **argv)
{
  optind=1; //in order to parse twice, this needs to be 1

  /*always leave the "p:" last to have parfiles overwrite input*/
  char args[] = "hg:o:R:M:q:X:x:Y:y:Z:z:L:l:f:d:i:T:S:F:p:";
  int opt;
  while((opt=getopt(argc, argv, args)) != -1){
    switch(opt){
      case 'h':
        if (!VERBOSE) PRSECTN(TEOBResumS_Info);
        printf(TEOBResumS_Usage);
        exit(0);
        break;
      case 'g':
        eobp->use_geometric_units = YESNO2INT(optarg);
        break;
      case 'o':
        eobp->output_hpc = YESNO2INT(optarg);
        break;
      case 'R':
        eobp->srate_interp = par_get_d(optarg);
        eobp->interp_uniform_grid = 1;
        break;
      case 'M':
        eobp->M = par_get_d(optarg);
        break;
      case 'q':
        eobp->q = par_get_d(optarg);
        break;
      case 'X':
        eobp->chi1x = par_get_d(optarg);
        break;
      case 'x':
        eobp->chi2x = par_get_d(optarg);
        break;
      case 'Y':
        eobp->chi1y = par_get_d(optarg);
        break;
      case 'y':
        eobp->chi2y = par_get_d(optarg);
        break;
      case 'Z':
        eobp->chi1 = par_get_d(optarg);
        eobp->chi1z = par_get_d(optarg);
        break;
      case 'z':
        eobp->chi2 = par_get_d(optarg);
        eobp->chi2z = par_get_d(optarg);
        break;
      case 'L':
        eobp->LambdaAl2 = par_get_d(optarg);
        /*set the defaults for BNS to avoid having to re-parse command line*/
        eobp->use_tidal = TIDES_TEOBRESUM3;
        eobp->use_tidal_gravitomagnetic = TIDES_GM_PN;
        eobp->pGSF_tidal = 4.0;
        eobp->use_lambda234_fits = Lambda234_fits_YAGI13;
        eobp->use_a6c_fits = a6c_fits_V0;
        eobp->use_cN3LO_fits = cN3LO_fits_NO;
        eobp->centrifugal_radius = CENTRAD_NNLO;
        eobp->use_flm = USEFLM_SSNLO;
        eobp->nqc = NQC_NO; 
        eobp->nqc_coefs_flx = NQC_FLX_NONE;
        eobp->nqc_coefs_hlm = NQC_HLM_NONE;
        break;
      case 'l':
        eobp->LambdaBl2 = par_get_d(optarg);
        /*set the defaults for BNS to avoid having to re-parse command line*/
        eobp->use_tidal = TIDES_TEOBRESUM3;
        eobp->use_tidal_gravitomagnetic = TIDES_GM_PN;
        eobp->pGSF_tidal = 4.0;
        eobp->use_lambda234_fits = Lambda234_fits_YAGI13;
        eobp->use_a6c_fits   = a6c_fits_V0;
        eobp->use_cN3LO_fits = cN3LO_fits_NO;
        eobp->centrifugal_radius = CENTRAD_NNLO;
        eobp->use_flm = USEFLM_SSNLO;
        eobp->nqc = NQC_NO;
        eobp->nqc_coefs_flx = NQC_FLX_NONE;
        eobp->nqc_coefs_hlm = NQC_HLM_NONE;
        break;
      case 'f':
        eobp->initial_frequency = par_get_d(optarg);
        break;
      case 'd':
        eobp->distance = par_get_d(optarg);
        break;
      case 'i':
        eobp->inclination = par_get_d(optarg);
        break;
      case 'p':
        EOBParameters_parse_file(optarg, eobp);
        break;
      default:
        errorexit("Invalid option specified, for more information run\n\
        ./TEOBResumS -h\n");
    }
  }
  if (eobp->LambdaAl2 > 1. && eobp->LambdaBl2 > 1.)
    eobp->binary = BINARY_BNS;
  else if (eobp->LambdaAl2 == 0. && eobp->LambdaBl2 > 1.)
    eobp->binary = BINARY_BHNS;
  else
    eobp->binary = BINARY_BBH;
  return eobp->binary;
}

/* Parse an input parfile */
void EOBParameters_parse_file(char *fname, EOBParameters *eobp)
{
  const char DELIMITERS_FOR_COMMENTS[] = "#";
  char line[2*STRLEN];
  char *key, *val;
  
  FILE *fp = fopen(fname,"r");
  if (!fp) errorexit(eob_error_msg[ERROR_FILEOPEN]);
  
  while (fgets(line,2*STRLEN,fp)!=NULL) {    
    if (is_blank(line)) continue;
    if (line[0] == '#') continue;
    remove_comments(line, DELIMITERS_FOR_COMMENTS);
    remove_white_spaces(line);
    if (getkv(line,&key,&val)) continue; 

    EOBParameters_set_key_val(eobp, key, val);

  } // while/fgets
  fclose(fp);
  
}

void EOBParameters_set_key_val(EOBParameters *eobp, char *key, char *val)
{

  if (STREQUAL(key,"use_geometric_units")) {    
    eobp->use_geometric_units = YESNO2INT(string_trim(val));
  } 
  /* Binary parameters */
  
  if (STREQUAL(key,"M")) {
    eobp->M = atof(val); 
  }
  if (STREQUAL(key,"q")) {
    eobp->q = par_get_d(val); 
  }
  if (STREQUAL(key,"chi1")) {
    eobp->chi1 = par_get_d(val);
  }
  if (STREQUAL(key,"chi1x")) {
    eobp->chi1x = par_get_d(val);
  }
  if (STREQUAL(key,"chi1y")) {
    eobp->chi1y = par_get_d(val);
  }
  if (STREQUAL(key,"chi1z")) {
    eobp->chi1z = par_get_d(val);
  }
  if (STREQUAL(key,"chi2")) {
    eobp->chi2 = par_get_d(val);
  }
  if (STREQUAL(key,"chi2x")) {
    eobp->chi2x = par_get_d(val);
  }
  if (STREQUAL(key,"chi2y")) {
    eobp->chi2y = par_get_d(val);
  }
  if (STREQUAL(key,"chi2z")) {
    eobp->chi2z = par_get_d(val);
  }
  if (STREQUAL(key,"distance")) {
    eobp->distance = par_get_d(val);
  }
  if (STREQUAL(key,"inclination")) {
    eobp->inclination = par_get_d(val);
  }
  if (STREQUAL(key,"coalescence_angle")) {
    eobp->coalescence_angle = par_get_d(val);
  }
  if (STREQUAL(key,"polarisation")) {
    eobp->polarization = par_get_d(val);
  }
  if (STREQUAL(key,"r0")) {
    eobp->r0 = par_get_d(val);
  }
  if (STREQUAL(key,"initial_frequency")) {
    eobp->initial_frequency = par_get_d(val);
  }
  if (STREQUAL(key,"LambdaAl2")) {
    eobp->LambdaAl2 = par_get_d(val);
  }
  if (STREQUAL(key,"LambdaBl2")) {
    eobp->LambdaBl2 = par_get_d(val);
  }
  if (STREQUAL(key,"LambdaAl3")) {
    eobp->LambdaAl3 = par_get_d(val);
  }
  if (STREQUAL(key,"LambdaBl3")) {
    eobp->LambdaBl3 = par_get_d(val);
  }
  if (STREQUAL(key,"LambdaAl4")) {
    eobp->LambdaAl4 = par_get_d(val);
  }
  if (STREQUAL(key,"LambdaBl4")) {
    eobp->LambdaBl4 = par_get_d(val);
  }
  if (STREQUAL(key,"SigmaAl2")) {
    eobp->SigmaAl2 = par_get_d(val);
  }
  if (STREQUAL(key,"SigmaBl2")) {
    eobp->SigmaBl2 = par_get_d(val);
  }

  /* EOB Settings */

  if (STREQUAL(key,"use_spins")) {
    eobp->use_spins = par_get_i(val); //FIXME: this was a bool, but could be used as integer to switch between ALIGNED/PRECESSING
  }

  if (STREQUAL(key,"tides")) {     
    val = string_trim(val);
    for (eobp->use_tidal=0; eobp->use_tidal<=TIDES_NOPT; eobp->use_tidal++) {
      if (eobp->use_tidal == TIDES_NOPT) {
        eobp->use_tidal = TIDES_OFF;
          if (VERBOSE) printf("tides '%s' undefined, set to '%s'\n",
          val,tides_opt[eobp->use_tidal]);
          break;
      }
    if (STREQUAL(val,tides_opt[eobp->use_tidal])) break;
    }
  }
  
  if (STREQUAL(key,"tides_gravitomagnetic")) {
    val = string_trim(val);
    for (eobp->use_tidal_gravitomagnetic=0; eobp->use_tidal_gravitomagnetic<=TIDES_GM_NOPT; eobp->use_tidal_gravitomagnetic++) {
      if (eobp->use_tidal_gravitomagnetic == TIDES_GM_NOPT) {
        eobp->use_tidal_gravitomagnetic = TIDES_GM_OFF;
        if (VERBOSE) printf("tides GM '%s' undefined, set to '%s'\n", val, tides_gravitomagnetic_opt[eobp->use_tidal_gravitomagnetic]);
        break;
      }
      if (STREQUAL(val, tides_gravitomagnetic_opt[eobp->use_tidal_gravitomagnetic])) break;
    }
  }
  
  if (STREQUAL(key,"pGSF_tidal")) {
    eobp->pGSF_tidal = par_get_d(val);
  }

  if (STREQUAL(key,"use_lambda234_fits")) {
    val = string_trim(val);
    for (eobp->use_lambda234_fits=0; eobp->use_lambda234_fits<=Lambda234_fits_NOPT; eobp->use_lambda234_fits++) {
      if (eobp->use_lambda234_fits == Lambda234_fits_NOPT) {
        eobp->use_lambda234_fits = Lambda234_fits_YAGI13;
        if (VERBOSE) printf("use_lambda234_fits '%s' undefined, set to '%s'\n",
        val, use_lambda234_fits_opt[eobp->use_lambda234_fits]);
      break;
      }
      if (STREQUAL(val, use_lambda234_fits_opt[eobp->use_lambda234_fits])) break;
    }
  }

  if (STREQUAL(key,"use_a6c_fits")) {
    val = string_trim(val);
    for (eobp->use_a6c_fits=0; eobp->use_a6c_fits<=a6c_fits_NOPT; eobp->use_a6c_fits++) {
      if (eobp->use_a6c_fits == a6c_fits_NOPT) {
        eobp->use_a6c_fits = a6c_fits_HM;
        if (VERBOSE) printf("use_a6c_fits '%s' undefined, set to '%s'\n",
        val, use_a6c_fits_opt[eobp->use_a6c_fits]);
      break;
      }
      if (STREQUAL(val, use_a6c_fits_opt[eobp->use_a6c_fits])) break;
    }
  }

  if (STREQUAL(key,"use_cN3LO_fits")) {
    val = string_trim(val);
    for (eobp->use_cN3LO_fits=0; eobp->use_cN3LO_fits<=cN3LO_fits_NOPT; eobp->use_cN3LO_fits++) {
      if (eobp->use_cN3LO_fits == cN3LO_fits_NOPT) {
        eobp->use_cN3LO_fits = cN3LO_fits_HM_2023_432;
        if (VERBOSE) printf("use_cN3LO_fits '%s' undefined, set to '%s'\n",
        val, use_cN3LO_fits_opt[eobp->use_cN3LO_fits]);
      break;
      }
      if (STREQUAL(val, use_cN3LO_fits_opt[eobp->use_cN3LO_fits])) break;
    }
  }

  if (STREQUAL(key,"use_tidal_fmode_model")) {
    eobp->use_tidal_fmode_model = YESNO2INT(string_trim(val));
  }
  
  if (STREQUAL(key,"use_speedytail")) {    
    eobp->use_speedytail = YESNO2INT(string_trim(val));
  }    
  if (STREQUAL(key,"dt_merger_interp")) {    
    eobp->dt_merger_interp = par_get_d(val);
  }

  if (STREQUAL(key,"interp_uniform_grid")) {    
    eobp->interp_uniform_grid = YESNO2INT(string_trim(val)); //FIXME: take from str comparisons // RG: made yes/no, for now
  }
  if (STREQUAL(key,"dt_interp")) {    
    eobp->dt_interp = par_get_d(val);
  }
  if (STREQUAL(key,"srate_interp")) {    
    eobp->srate_interp = par_get_d(val);
  }
  
  if (STREQUAL(key,"use_mode_lm")) {
    free(eobp->use_mode_lm);
    eobp->use_mode_lm_size = str2iarray(val, &eobp->use_mode_lm);
    /* Over-write use_mode_lm_inertial to use_mode_lm (default behavior if
    use_mode_lm_inertial is not specified) */
    free(eobp->use_mode_lm_inertial);
    eobp->use_mode_lm_inertial_size = str2iarray(val, &eobp->use_mode_lm_inertial);
  }  

  if (STREQUAL(key,"use_mode_lm_inertial")) {
    free(eobp->use_mode_lm_inertial);
    eobp->use_mode_lm_inertial_size = str2iarray(val, &eobp->use_mode_lm_inertial);
  }  

  if (STREQUAL(key,"kpostpeak")) {
    free(eobp->kpostpeak);
    eobp->kpostpeak_size = str2iarray(val, &eobp->kpostpeak);
  }

  if (STREQUAL(key,"knqcpeak22")) {
    free(eobp->knqcpeak22);
    eobp->knqcpeak22_size = str2iarray(val, &eobp->knqcpeak22);
  }

  if (STREQUAL(key,"centrifugal_radius")) {
    val = string_trim(val);
    for (eobp->centrifugal_radius=0; eobp->centrifugal_radius<=CENTRAD_NOPT; eobp->centrifugal_radius++) {
if (eobp->centrifugal_radius == CENTRAD_NOPT) {
  eobp->centrifugal_radius = CENTRAD_NLO;
  if (VERBOSE) printf("centrifugal_radius '%s' undefined, set to '%s'\n",
          val, centrifugal_radius_opt[eobp->centrifugal_radius]);
  break;
}
if (STREQUAL(val, centrifugal_radius_opt[eobp->centrifugal_radius])) break;
    }
  }

  if (STREQUAL(key,"use_flm")) {
    val = string_trim(val);
    for (eobp->use_flm=0; eobp->use_flm<=USEFLM_NOPT; eobp->use_flm++) {
if (eobp->use_flm == USEFLM_NOPT) {
eobp->use_flm = USEFLM_HM;
if (VERBOSE) printf("use_flm '%s' undefined, set to '%s'\n",
        val, use_flm_opt[eobp->use_flm]);
break;
}
if (STREQUAL(val, use_flm_opt[eobp->use_flm])) break;
    }
  }
  
  if (STREQUAL(key,"compute_LR")) {
    eobp->compute_LR = YESNO2INT(string_trim(val));
  }
  if (STREQUAL(key,"compute_LSO")) {    
    eobp->compute_LSO = YESNO2INT(string_trim(val));
  }
  if (STREQUAL(key,"compute_LR_guess")) {    
    eobp->compute_LR_guess = par_get_d(val);
  }
  if (STREQUAL(key,"compute_LSO_guess")) {    
    eobp->compute_LSO_guess = par_get_d(val);
  }
  
  /* Precession settings */

  if (STREQUAL(key,"spin_flx")) {     
    val = string_trim(val);
    for (eobp->spin_flx=0; eobp->spin_flx<=SPIN_FLX_NOPT; eobp->spin_flx++) {
      if (eobp->spin_flx == SPIN_FLX_NOPT) {
        eobp->spin_flx = SPIN_FLX_PN;
        if (VERBOSE) printf("spin flux '%s' undefined, set to '%s'\n",
        val,spin_flx_opt[eobp->spin_flx]);
        break;
      }
    if (STREQUAL(val,spin_flx_opt[eobp->spin_flx])) break;
    }
  }    

  if (STREQUAL(key,"ringdown_eulerangles")) {     
    val = string_trim(val);
    for (eobp->ringdown_eulerangles=0; eobp->ringdown_eulerangles<=RD_EULERANGLES_NOPT; eobp->ringdown_eulerangles++) {
      if (eobp->ringdown_eulerangles == RD_EULERANGLES_NOPT) {
        eobp->ringdown_eulerangles = RD_EULERANGLES_QNMs;
        if (VERBOSE) printf("ringdown euler angles '%s' undefined, set to '%s'\n",
        val,ringdown_eulerangles_opt[eobp->ringdown_eulerangles]);
        break;
      }
    if (STREQUAL(val,ringdown_eulerangles_opt[eobp->ringdown_eulerangles])) break;
    }
  }   

  if (STREQUAL(key,"spin_dyn_size")) {
    eobp->spin_dyn_size = par_get_i(val); //FIXME: this was a bool, but could be used as integer to switch between ALIGNED/PRECESSING
  }
  if (STREQUAL(key,"spin_odes_omg_stop")) {
    eobp->spin_odes_omg_stop = par_get_d(val); //FIXME: this was a bool, but could be used as integer to switch between ALIGNED/PRECESSING
  }
  if (STREQUAL(key,"spin_odes_t_stop")) {
    eobp->spin_odes_t_stop = par_get_d(val); //FIXME: this was a bool, but could be used as integer to switch between ALIGNED/PRECESSING
  }
  if (STREQUAL(key,"spin_odes_dt")) {
    eobp->spin_odes_dt = par_get_d(val); //FIXME: this was a bool, but could be used as integer to switch between ALIGNED/PRECESSING
  }

  /* NQC */
  
  if (STREQUAL(key,"nqc")) {
    val = string_trim(val);
    for (eobp->nqc=0; eobp->nqc<=NQC_NOPT; eobp->nqc++) {
      if (STREQUAL(val,nqc_opt[eobp->nqc])) break;
    }
    if (eobp->nqc == NQC_NOPT) {
      eobp->nqc = NQC_AUTO;
      if (VERBOSE) printf("nqc '%s' undefined, set to '%s'\n",
        val, nqc_opt[eobp->nqc]);
    }
  }

  if (STREQUAL(key,"nqc_coefs_flx")) {
    val = string_trim(val);
    for (eobp->nqc_coefs_flx=0; eobp->nqc_coefs_flx<=NQC_FLX_NOPT; eobp->nqc_coefs_flx++) {
if (eobp->nqc_coefs_flx == NQC_FLX_NOPT) {
  eobp->nqc_coefs_flx = NQC_FLX_NONE;
  if (VERBOSE) printf("nqc '%s' undefined, set to '%s'\n",
          val, nqc_flx_opt[eobp->nqc_coefs_flx]);
  break;
}
if (STREQUAL(val, nqc_flx_opt[eobp->nqc_coefs_flx])) break;
    }
  }

  if (STREQUAL(key,"nqc_coefs_hlm")) {
    val = string_trim(val);
    for (eobp->nqc_coefs_hlm=0; eobp->nqc_coefs_hlm<=NQC_HLM_NOPT; eobp->nqc_coefs_hlm++) {
if (eobp->nqc_coefs_hlm == NQC_HLM_NOPT) {
  eobp->nqc_coefs_hlm = NQC_HLM_NONE;
  if (VERBOSE) printf("nqc '%s' undefined, set to '%s'\n",
          val, nqc_hlm_opt[eobp->nqc_coefs_hlm]);
  break;
}
if (STREQUAL(val, nqc_hlm_opt[eobp->nqc_coefs_hlm])) break;
    }
  }

  if (STREQUAL(key,"nqc_coefs_flx_file")) {      
    if (par_get_s(eobp->nqc_coefs_flx_file, val))
      errorexit("Error converting input string");
  }
  if (STREQUAL(key,"nqc_coefs_hlm_file")) {
    if (par_get_s(eobp->nqc_coefs_hlm_file, val))
      errorexit("Error converting input string");
  }

  /* Post-adiabatic */
  
  if (STREQUAL(key,"postadiabatic_dynamics")) {    
    eobp->postadiabatic_dynamics = YESNO2INT(string_trim(val));
  }
  if (STREQUAL(key,"postadiabatic_dynamics_N")) {      
    eobp->postadiabatic_dynamics_N = par_get_i(val);
  }
  if (STREQUAL(key,"postadiabatic_dynamics_size")) {      
    eobp->postadiabatic_dynamics_size = par_get_i(val);
  }
  if (STREQUAL(key,"postadiabatic_dynamics_rmin")) {
    eobp->postadiabatic_dynamics_rmin = par_get_d(val);
  }
  if (STREQUAL(key,"postadiabatic_dynamics_stop")) {
    eobp->postadiabatic_dynamics_stop = YESNO2INT(string_trim(val));
  }
  
  /* Evolution settings */

  if (STREQUAL(key,"srate")) {
    eobp->srate = par_get_d(val);
  }
  if (STREQUAL(key,"dt")) {
    eobp->dt = par_get_d(val);
  }
  if (STREQUAL(key,"size")) {
    eobp->size = par_get_i(val);
  }
  if (STREQUAL(key,"ringdown_extend_array")) {
    eobp->ringdown_extend_array = par_get_i(val);
  }

  if (STREQUAL(key,"ode_timestep")) {
    val = string_trim(val);
    for (eobp->ode_timestep=0; eobp->ode_timestep<=ODE_TSTEP_NOPT; eobp->ode_timestep++) {
if (eobp->ode_timestep==ODE_TSTEP_NOPT) {
  eobp->ode_timestep = ODE_TSTEP_ADAPTIVE;
  if (VERBOSE) printf("ode_timestep '%s' undefined, set to default %s\n",
          val, ode_tstep_opt[eobp->ode_timestep]);
  break;
}
if (STREQUAL(val,ode_tstep_opt[eobp->ode_timestep])) break;
    }
  }

  if (STREQUAL(key,"ode_abstol")) {
    eobp->ode_abstol = par_get_d(val);
  }
  if (STREQUAL(key,"ode_reltol")) {
    eobp->ode_reltol = par_get_d(val);
  }
  if (STREQUAL(key,"ode_tmax")) {
    eobp->ode_tmax = par_get_d(val);
  }
  if (STREQUAL(key,"ode_stop_at_radius")) {
    eobp->ode_stop_radius = par_get_d(val);
  }
  if (STREQUAL(key,"ode_stop_afterNdt")) {
    eobp->ode_stop_afterNdt = par_get_i(val);
  }
  
  /* Output */
  
  if (STREQUAL(key,"output_dir")) {
    if (par_get_s(eobp->output_dir, val))
      errorexit("Error converting input string");
  }
  if (STREQUAL(key,"output_hpc")) {
    eobp->output_hpc = YESNO2INT(string_trim(val));
  }
  if (STREQUAL(key,"output_multipoles")) {
    eobp->output_multipoles = YESNO2INT(string_trim(val));
  }
  if (STREQUAL(key,"output_lm")) {
    free(eobp->output_lm);
    eobp->output_lm_size = str2iarray(val, &eobp->output_lm);
  }
  if (STREQUAL(key,"output_dynamics")) {
    eobp->output_dynamics = YESNO2INT(string_trim(val));
  }
  if (STREQUAL(key,"output_nqc")) {
    eobp->output_nqc = YESNO2INT(string_trim(val));
  }
  if (STREQUAL(key,"output_nqc_coefs")) {
    eobp->output_nqc_coefs = YESNO2INT(string_trim(val));
  }
  if (STREQUAL(key,"output_ringdown")) {
    eobp->output_ringdown = YESNO2INT(string_trim(val));
  }

  /* FD */
  if (STREQUAL(key,"domain")) {
    eobp->domain = par_get_i(val);
  }
  if (STREQUAL(key,"time_shift_FD")) {
    eobp->time_shift_FD = par_get_i(val);
  }      
  if (STREQUAL(key,"df")) {
    eobp->df = par_get_d(val);
   } 
}

void EOBParameters_tofile (EOBParameters *eobp, char *fname)
{
  //Dump everything to a file for reproducibility
  FILE *f = fopen(fname, "w");
  
  fprintf(f,"%s = \"%s\"\n", "use_geometric_units", INT2YESNO(eobp->use_geometric_units));

  /* intrinsic parameters */
  fprintf(f,"%s = %.16f\n", "M", eobp->M);
  fprintf(f,"%s = %.16f\n", "q", eobp->q);
  fprintf(f,"%s = %.16f\n", "chi1", eobp->chi1);
  fprintf(f,"%s = %.16f\n", "chi2", eobp->chi2);
  fprintf(f,"%s = %.16f\n", "chi1x", eobp->chi1x);
  fprintf(f,"%s = %.16f\n", "chi1y", eobp->chi1y);
  fprintf(f,"%s = %.16f\n", "chi1z", eobp->chi1z);
  fprintf(f,"%s = %.16f\n", "chi2x", eobp->chi2x);
  fprintf(f,"%s = %.16f\n", "chi2y", eobp->chi2y);
  fprintf(f,"%s = %.16f\n", "chi2z", eobp->chi2z);
  fprintf(f,"%s = %.16f\n", "distance", eobp->distance);
  fprintf(f,"%s = %.16f\n", "inclination", eobp->inclination);
  fprintf(f,"%s = %.16f\n", "coalescence_angle", eobp->coalescence_angle);
  fprintf(f,"%s = %.16f\n", "polarization", eobp->polarization);
  fprintf(f,"%s = %.16f\n", "r0", eobp->r0);
  fprintf(f,"%s = %.16f\n", "initial_frequency", eobp->initial_frequency);
  fprintf(f,"%s = %.16f\n", "LambdaAl2", eobp->LambdaAl2);
  fprintf(f,"%s = %.16f\n", "LambdaBl2", eobp->LambdaBl2);
  fprintf(f,"%s = %.16f\n", "LambdaAl3", eobp->LambdaAl3);
  fprintf(f,"%s = %.16f\n", "LambdaBl3", eobp->LambdaBl3);
  fprintf(f,"%s = %.16f\n", "LambdaAl4", eobp->LambdaAl4);
  fprintf(f,"%s = %.16f\n", "LambdaBl4", eobp->LambdaBl4);
  fprintf(f,"%s = %.16f\n", "SigmaAl2", eobp->SigmaAl2);
  fprintf(f,"%s = %.16f\n", "SigmaBl2", eobp->SigmaBl2);

  /* Derived parameters */
  fprintf(f,"%s = %.16f\n", "nu", eobp->nu);
  fprintf(f,"%s = %.16f\n", "X1", eobp->X1);
  fprintf(f,"%s = %.16f\n", "X2", eobp->X2);
  fprintf(f,"%s = %.16f\n", "S2", eobp->S1);
  fprintf(f,"%s = %.16f\n", "S2", eobp->S2);
  fprintf(f,"%s = %.16f\n", "S",  eobp->S);
  fprintf(f,"%s = %.16f\n", "Sstar",  eobp->Sstar);
  fprintf(f,"%s = %.16f\n", "a1",  eobp->a1);
  fprintf(f,"%s = %.16f\n", "a2",  eobp->a2);
  fprintf(f,"%s = %.16f\n", "aK",  eobp->aK);
  fprintf(f,"%s = %.16f\n", "aK2", eobp->aK2);
  fprintf(f,"%s = %.16f\n", "C_Q1", eobp->C_Q1);
  fprintf(f,"%s = %.16f\n", "C_Q2", eobp->C_Q2);
  fprintf(f,"%s = %.16f\n", "C_Oct1", eobp->C_Oct1);
  fprintf(f,"%s = %.16f\n", "C_Oct2", eobp->C_Oct2);
  fprintf(f,"%s = %.16f\n", "C_Hex1", eobp->C_Hex1);
  fprintf(f,"%s = %.16f\n", "C_Hex2", eobp->C_Hex2);
  fprintf(f,"%s = %.16f\n", "a6c", eobp->a6c);
  fprintf(f,"%s = %.16f\n", "cN3LO", eobp->cN3LO);

  fprintf(f,"%s = %.16f\n","kappaAl2",  eobp->kapA2); // gravitoelectric kappa star A
  fprintf(f,"%s = %.16f\n","kappaAl3",  eobp->kapA3); //
  fprintf(f,"%s = %.16f\n","kappaAl4",  eobp->kapA4); //
  fprintf(f,"%s = %.16f\n","kappaAl5",  eobp->kapA5); //
  fprintf(f,"%s = %.16f\n","kappaAl6",  eobp->kapA6); //
  fprintf(f,"%s = %.16f\n","kappaAl7",  eobp->kapA7); //
  fprintf(f,"%s = %.16f\n","kappaAl8",  eobp->kapA8); //

  fprintf(f,"%s = %.16f\n","kappaBl2",  eobp->kapB2); //
  fprintf(f,"%s = %.16f\n","kappaBl3",  eobp->kapB3); //
  fprintf(f,"%s = %.16f\n","kappaBl4",  eobp->kapB4); //
  fprintf(f,"%s = %.16f\n","kappaBl5",  eobp->kapB5); //
  fprintf(f,"%s = %.16f\n","kappaBl6",  eobp->kapB6); //
  fprintf(f,"%s = %.16f\n","kappaBl7",  eobp->kapB7); //
  fprintf(f,"%s = %.16f\n","kappaBl8",  eobp->kapB8); //

  fprintf(f,"%s = %.16f\n","kappajAl2",  eobp->japA2); // gravitomagnetic kappa star A
  fprintf(f,"%s = %.16f\n","kappajBl2",  eobp->japB2);
  fprintf(f,"%s = %.16f\n","kappajTl2",  eobp->japT2);
  
  fprintf(f,"%s = %.16f\n","kappaTl2",  eobp->kapT2); //
  fprintf(f,"%s = %.16f\n","kappaTl3",  eobp->kapT3); //
  fprintf(f,"%s = %.16f\n","kappaTl4",  eobp->kapT4); //
  fprintf(f,"%s = %.16f\n","bar_alph2_1", eobp->bar_alph2_1); //
  fprintf(f,"%s = %.16f\n","bar_alph2_2", eobp->bar_alph2_2); //
  fprintf(f,"%s = %.16f\n","bar_alph3_1", eobp->bar_alph3_1); //
  fprintf(f,"%s = %.16f\n","bar_alph3_2", eobp->bar_alph3_2); //
  fprintf(f,"%s = %.16f\n","bar_alph2j_1", eobp->bar_alph2j_1); //
  fprintf(f,"%s = %.16f\n","khatAl2",  eobp->khatA2); //
  fprintf(f,"%s = %.16f\n","khatBl2",  eobp->khatB2); //
  fprintf(f,"%s = %.16f\n","rLR_tidal",  eobp->rLR_tidal); // radius of light-ring for NNLO tidal model
  
  fprintf(f,"%s = %.16f\n","BH_final_mass",  eobp->Mbhf); // final BH mass
  fprintf(f,"%s = %.16f\n","BH_final_spin",  eobp->abhf); 

  /* EOB Settings */

  fprintf(f,"%s = %d\n"    , "use_spins", eobp->use_spins);
  fprintf(f,"%s = \"%s\"\n", "tides", tides_opt[eobp->use_tidal]);
  fprintf(f,"%s = \"%s\"\n", "tides_gravitomagnetic", tides_gravitomagnetic_opt[eobp->use_tidal_gravitomagnetic]);
  fprintf(f,"%s = %.16f\n" , "pGSF_tidal", eobp->pGSF_tidal);
  fprintf(f,"%s = \"%s\"\n", "use_lambda234_fits", use_lambda234_fits_opt[eobp->use_lambda234_fits]);
  fprintf(f,"%s = \"%s\"\n", "use_a6c_fits_opt", use_a6c_fits_opt[eobp->use_a6c_fits]);
  fprintf(f,"%s = \"%s\"\n", "use_cN3LO_fits_opt", use_cN3LO_fits_opt[eobp->use_cN3LO_fits]);
  fprintf(f,"%s = \"%s\"\n", "use_tidal_fmode_model", INT2YESNO(eobp->use_tidal_fmode_model));
  fprintf(f,"%s = \"%s\"\n", "use_speedytail", INT2YESNO(eobp->use_speedytail));
  fprintf(f,"%s = %.16f\n" , "dt_merger_interp", eobp->dt_merger_interp);
  fprintf(f,"%s = %.16f\n" , "dt_interp", eobp->dt_interp);
  fprintf(f,"%s = %.16f\n" , "srate_interp", eobp->srate_interp);
  fprintf(f,"%s = \"%s\"\n", "interp_uniform_grid", INT2YESNO(eobp->interp_uniform_grid));

  fprintf(f,"%s = [", "use_mode_lm");
  for(int i=0; i<eobp->use_mode_lm_size-1;i++)
    fprintf(f,"%d,", eobp->use_mode_lm[i]);
  fprintf(f,"%d]\n", eobp->use_mode_lm[eobp->use_mode_lm_size-1]);

  fprintf(f,"%s = [", "use_mode_lm_inertial");
  for(int i=0; i<eobp->use_mode_lm_inertial_size-1;i++)
    fprintf(f,"%d,", eobp->use_mode_lm_inertial[i]);
  fprintf(f,"%d]\n", eobp->use_mode_lm_inertial[eobp->use_mode_lm_inertial_size-1]);
  fprintf(f,"%s = [", "kpostpeak");
  for(int i=0; i<eobp->kpostpeak_size-1;i++)
    fprintf(f,"%d,", eobp->kpostpeak[i]);
  fprintf(f,"%d]\n", eobp->kpostpeak[eobp->kpostpeak_size-1]);
  
  fprintf(f,"%s = [", "knqcpeak22");
  for(int i=0; i<eobp->knqcpeak22_size-1;i++)
    fprintf(f,"%d,", eobp->knqcpeak22[i]);
  fprintf(f,"%d]\n", eobp->knqcpeak22[eobp->knqcpeak22_size-1]);

  fprintf(f,"%s = \"%s\"\n", "centrifugal_radius", centrifugal_radius_opt[eobp->centrifugal_radius]);
  fprintf(f,"%s = \"%s\"\n", "use_flm", use_flm_opt[eobp->use_flm]);
  fprintf(f,"%s = \"%s\"\n", "compute_LR", INT2YESNO(eobp->compute_LR));
  fprintf(f,"%s = %d\n"    , "compute_LR_guess", eobp->compute_LR_guess);
  fprintf(f,"%s = \"%s\"\n", "compute_LSO", INT2YESNO(eobp->compute_LSO));
  fprintf(f,"%s = %d\n"    , "compute_LSO_guess", eobp->compute_LSO_guess);

  /* NQC */
  fprintf(f,"%s = \"%s\"\n", "nqc", nqc_opt[eobp->nqc]);
  fprintf(f,"%s = \"%s\"\n", "nqc_coefs_flx", nqc_flx_opt[eobp->nqc_coefs_flx]);
  fprintf(f,"%s = \"%s\"\n", "nqc_coefs_hlm", nqc_hlm_opt[eobp->nqc_coefs_flx]);
  fprintf(f,"%s = \"%s\"\n", "nqc_coefs_flx_file", eobp->nqc_coefs_flx_file);
  fprintf(f,"%s = \"%s\"\n", "nqc_coefs_hlm_file", eobp->nqc_coefs_hlm_file);

  /* Post-adiabatic */
  fprintf(f,"%s = \"%s\"\n", "postadiabatic_dynamics", INT2YESNO(eobp->postadiabatic_dynamics));
  fprintf(f,"%s = %d\n"    , "postadiabatic_dynamics_N", eobp->postadiabatic_dynamics_N);
  fprintf(f,"%s = %d\n"    , "postadiabatic_dynamics_size", eobp->postadiabatic_dynamics_size);
  fprintf(f,"%s = %.16f\n" , "postadiabatic_dynamics_rmin", eobp->postadiabatic_dynamics_rmin);
  fprintf(f,"%s = \"%s\"\n", "postadiabatic_dynamics_stop", INT2YESNO(eobp->postadiabatic_dynamics_stop));

  /* Evolution settings */
  fprintf(f,"%s = %.16f\n" , "srate", eobp->srate);
  fprintf(f,"%s = %.16f\n" , "dt", eobp->dt);
  fprintf(f,"%s = %d\n"    , "size", eobp->size);
  fprintf(f,"%s = %d\n"    , "ringdown_extend_array", eobp->ringdown_extend_array);
  fprintf(f,"%s = \"%s\"\n", "ode_timestep", ode_tstep_opt[eobp->ode_timestep]);
  fprintf(f,"%s = %E\n"    , "ode_abstol", eobp->ode_abstol);
  fprintf(f,"%s = %E\n"    , "ode_reltol", eobp->ode_reltol);
  fprintf(f,"%s = %.16f\n" , "ode_tmax", eobp->ode_tmax);
  fprintf(f,"%s = %f\n"    , "ode_stop_at_radius", eobp->ode_stop_radius);
  fprintf(f,"%s = %d\n"    , "ode_stop_afterNdt", eobp->ode_stop_afterNdt);

  /* Output */
  fprintf(f,"%s = \"%s\"\n", "output_dir", eobp->output_dir);
  fprintf(f,"%s = \"%s\"\n", "output_hpc", INT2YESNO(eobp->output_hpc));
  fprintf(f,"%s = \"%s\"\n", "output_multipoles", INT2YESNO(eobp->output_multipoles));

  fprintf(f,"%s = [", "output_lm");
  for(int i=0; i<eobp->output_lm_size-1;i++)
    fprintf(f,"%d,", eobp->output_lm[i]);
  fprintf(f,"%d]\n", eobp->output_lm[eobp->output_lm_size-1]);

  fprintf(f,"%s = \"%s\"\n", "output_dynamics", INT2YESNO(eobp->output_dynamics));
  fprintf(f,"%s = \"%s\"\n", "output_nqc", INT2YESNO(eobp->output_nqc));
  fprintf(f,"%s = \"%s\"\n", "output_nqc_coefs", INT2YESNO(eobp->output_nqc_coefs));
  fprintf(f,"%s = \"%s\"\n", "output_ringdown", INT2YESNO(eobp->output_ringdown));
  
  fclose(f);
}

#if (DEBUG_THIS_FILE)

/* test 
   gcc TEOBResumSPars.c -lconfig -o testpars.x */
int main (int argc, char* argv[])
{
  EOBParameters_parse_file (argv[1], EOBPars);
}

#endif
