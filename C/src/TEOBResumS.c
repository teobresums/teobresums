/** \file TEOBResumS.c
 *  \brief TEOBResumS main routines
 *
 * EOBRun and main routines for TEOBResumS.
 * 
 * Copyright (C) 2017-2023 See AUTHORS file
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
 * LINDEX
 * ------
 * Global array, given k \in [0, 35] returns associated ell, 
 * defined as external in header.
 */
const int LINDEX[KMAX] = {
  2,2,
  3,3,3,
  4,4,4,4,
  5,5,5,5,5,
  6,6,6,6,6,6,
  7,7,7,7,7,7,7,
  8,8,8,8,8,8,8,8};

/** 
 * MINDEX
 * ------
 * Global array, given k \in [0, 35] returns associated emm, 
 * defined as external in header.
 */
const int MINDEX[KMAX] = {
  1,2,
  1,2,3,
  1,2,3,4,
  1,2,3,4,5,
  1,2,3,4,5,6,
  1,2,3,4,5,6,7,
  1,2,3,4,5,6,7,8};

/** 
 * KINDEX
 * ------
 * Global array, given \ell \in [2, 8] and emm \in [1, 8] returns 
 * associated k, defined as external in header. 
 */
const int KINDEX[9][9] = {  // l (m = 1 ...l)
  {-1,-1,-1,-1,-1,-1,-1,-1},// 0 
  {-1,-1,-1,-1,-1,-1,-1,-1},// 1 
  { 0, 1,-1,-1,-1,-1,-1,-1},// 2 (1 2)
  { 2, 3, 4,-1,-1,-1,-1,-1},// 3 (1 2 3)
  { 5, 6, 7, 8,-1,-1,-1,-1},// 4 (1 2 3 4)
  { 9,10,11,12,13,-1,-1,-1},// 5 (1 2 3 4 5)
  {14,15,16,17,18,19,-1,-1},// 6 (1 ... 6)
  {20,21,22,23,24,25,26,-1},// 7 (1 ... 7)
  {27,28,29,30,31,32,33,34},// 8 (1 ... 8)
};

/** 
 * NQCdata
 * ---------------
 * Global var for NQC coefficient 
 */
NQCdata *NQC;

/** 
 *  Function: main
 *  ---------------
 *   TEOBResumS main called by executasble TEOBResumS.x.
 *   Reads user input, set parameters accordingly and calls EOBRun
 *   
 *   @param[in] argc: argument count, to be passed to getopt via the command-line parser
 *   @param[in] argv: argument array, to be passed to getopt via the command-line parser
 *
 *   @return  zero if the waveform was evaluated correctly
 *              one if any error was raised
 */
int main (int argc, char* argv[])
{   
  if (VERBOSE) PRSECTN(TEOBResumS_Info);
  
  Waveform *hpc = NULL;         /* TD wvf */
  Waveform_lm *hmodes = NULL;   /* modes */
  Waveform_lm *hTmodes = NULL;  /* twisted modes, m>0*/
  Waveform_lm *hTmmodes = NULL; /* twisted modes, m<0 */
  Waveform_lm *hT0modes = NULL; /* twisted modes, m<0 */

  WaveformFD *hfpc = NULL;       /* FD wvf */
  WaveformFD_lm *hfmodes = NULL; /* modes */

  Dynamics *dynf = NULL;

  WaveformFD_lm *hfTmodes = NULL; /* twisted modes */
  
  int fc    = 1;          /* firstcall, set to 1 for now */
  int dc    = BINARY_BBH; /* default_choice, set to BBH */
  int model = MODEL_DALI; /* default_choice, set to Dalì */
  /* Init parameters & set defaults */
  EOBParameters_alloc( &EOBPars );
  EOBParameters_defaults (dc, model, EOBPars);

  if (argv[1]!=NULL) {
    /* Deal with input parfile or command line arguments */
    dc = EOBParameters_parse_commandline(EOBPars,argc,argv);
  }

  /* Properly reset the default & re-read the parameters */
  if ((EOBPars->ecc != 0.) || (EOBPars->r_hyp != 0.))
    EOBPars->model = MODEL_DALI;
    
  if (EOBPars->use_mode_lm)          free (EOBPars->use_mode_lm);
  if (EOBPars->use_mode_lm_inertial) free (EOBPars->use_mode_lm_inertial);
  if (EOBPars->use_mode_lm_nqc)      free (EOBPars->use_mode_lm_nqc);
  if (EOBPars->kpostpeak)            free (EOBPars->kpostpeak);
  if (EOBPars->knqcpeak22)           free (EOBPars->knqcpeak22);
  if (EOBPars->output_lm)            free (EOBPars->output_lm);
  if (EOBPars->freqs)                free (EOBPars->freqs);

  EOBParameters_defaults (dc, EOBPars->model, EOBPars);
  if (argv[1]!=NULL) {
    EOBParameters_parse_commandline(EOBPars,argc,argv);
  }

  const int output = EOBPars->output_dynamics
    + EOBPars->output_multipoles
    + EOBPars->output_hpc
    + EOBPars->output_nqc;  
  
  if (output) {
    if (system_mkdir(EOBPars->output_dir)) {
      printf("ERROR(TEOBResumS): %s\n",eob_error_msg[ERROR_MKDIR]);
      return ERROR_MKDIR;
    }
  }
    
  /* Set all firstcalls = 1 */
  for (int k=0; k < NFIRSTCALL; k++){
    EOBPars->firstcall[k] = 1;
  }
  
  int status = OK;
  /* set domain */
  if (eob_set_params(dc, fc)) {
    printf("ERROR(TEOBResumS): %s\n",eob_error_msg[ERROR_SET_PARAMS]);
    status = ERROR_SET_PARAMS;
    goto EXIT_POINT_MAIN;
  }
  if (output){
    char outpar[STRLEN];
    strcpy(outpar,EOBPars->output_dir);
    EOBParameters_tofile(EOBPars,strcat(outpar,"/params.txt"));
  }
  /* TD hpc, FD hpc, TD modes, FD modes, default_choice, firstcall */
  status = EOBRun(&hpc, &hfpc, 
		      &hmodes, &hfmodes,&dynf,
		      &hTmodes, &hTmmodes, &hT0modes,
		      &hfTmodes,
		      dc, fc);
  if (status) printf("ERROR(TEOBResumS): %s\n",eob_error_msg[status]);

EXIT_POINT_MAIN:;

  Waveform_free (hpc);
  WaveformFD_free (hfpc);
  Waveform_lm_free (hmodes);
  Waveform_lm_free (hTmodes);
  Waveform_lm_free (hTmmodes);
  Waveform_lm_free (hT0modes);
  WaveformFD_lm_free (hfmodes);
  Dynamics_free (dynf);
  WaveformFD_lm_free (hfTmodes);

  EOBParameters_free (EOBPars);
  return status;
}

/**
 *  Function : EOBRun
 *  ---------------
 *    Main routine for TD & FD EOB, assumes the EOBPars are allocated and 
 *    filled.
 *    It is divided into 4 main blocks:
 *      1. Initialization and preliminary computations
 *      2. Postadiabatic dynamics and/or solution of the EOB ODE system
 *      3. NQCs and ringdown attachment
 *      4. h+, hx computation
 *
 *    @param[out] hpc            : struct for waveform polarizations (time domain)  
 *    @param[out] hfpc           : struct for waveform polarizations (frequency domain)  
 *    @param[out] hmodes         : struct for co-precessing waveform multipoles    (time domain)  
 *    @param[out] hfmodes        : struct for co-precessing waveform multipoles    (frequency domain)  
 *    @param[out] dynf           : struct for the EOB dynamics  
 *    @param[out] hTmodes        : struct for inertial waveform multipoles    (time domain, m>0)  
 *    @param[out] hTmmodes       : struct for inertial waveform multipoles    (time domain, m<0)  
 *    @param[out] hT0odes        : struct for inertial waveform multipoles    (time domain, m=0)  
 *    @param[in]  default_choice : Type of binary (UNUSED)  
 *    @param[in]  firstcall      : firstcall parameters (UNUSED) 
 *    @return                      zero if the waveform was evaluated correctly, one if any error was raised  
 */ 
int EOBRun(Waveform **hpc, WaveformFD **hfpc, 
	   Waveform_lm **hmodes, WaveformFD_lm **hfmodes, Dynamics **dynf,
	   Waveform_lm **hTmodes,  Waveform_lm **hTmmodes, Waveform_lm **hT0modes,
	   WaveformFD_lm **hfTmodes,
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

  /* Do any output ? */
  const int output = EOBPars->output_dynamics + EOBPars->output_multipoles + EOBPars->output_hpc + EOBPars->output_nqc;    
  if (output) {
    /* Make output dir */
    if (system_mkdir(EOBPars->output_dir)) {
      status = ERROR_MKDIR;
      goto EXIT_POINT;
    }
  }

  /* Switch to mass-rescaled geometric units (if needed)*/
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
  
  /* Set useful pars/vars */
  const double q      = EOBPars->q;
  const double nu     = EOBPars->nu;
  const double chi1   = EOBPars->chi1;
  const double chi2   = EOBPars->chi2;
  const double ecc    = EOBPars->ecc;
  const double r_hyp  = EOBPars->r_hyp;
  const double H_hyp  = EOBPars->H_hyp;
  const double j_hyp  = EOBPars->j_hyp;
  const int use_spins = EOBPars->use_spins;
  int use_tidal = EOBPars->use_tidal;
  int store_dynamics  = EOBPars->output_dynamics; 
  if (!(EOBPars->binary == BINARY_BNS)) store_dynamics = 1; /* NQC determination need dynamical variables */
  if (use_spins == MODE_SPINS_GENERIC) store_dynamics = 1; /* Precession needs dynamical variables */
  if (ecc != 0.)  store_dynamics = 1; /* Eccentric waveform computation needs dynamical variables (sigmoid) */
  if (EOBPars->model == MODEL_DALI) EOBPars->postadiabatic_dynamics = 0;
  int use_postadiab_dyn = EOBPars->postadiabatic_dynamics;
  if (use_postadiab_dyn) store_dynamics = 1;
  const double dt = EOBPars->dt;
  
  /* *****************************************
   * Set Memory & do preliminary computations
   * *****************************************
   */

  /* Alloc memory for dynamics and multipolar waveform */
  Dynamics *dyn         = NULL;
  DynamicsSpin *spindyn = NULL;

  Waveform_lm *hlm    = NULL; /* h_lm */ 
  WaveformFD_lm *hflm = NULL; /* hf_lm (FD) */
  
  Waveform_lm *hTlm     = NULL; /* h_lm twisted */
  Waveform_lm *hTlm_neg = NULL; /* h_lm twisted */
  Waveform_lm *hTl0     = NULL;

  WaveformFD_lm *hfTlm = NULL; /* h_lm twisted (FD) */
  
  Waveform_lm_t *hlm_t = NULL;
  Waveform_lm *hlm_nqc = NULL; /* NQC */
  Waveform_lm *hlm_mrg = NULL; /* merger chunk */
  Dynamics *dyn_mrg = NULL;
  double ytmp[EOB_EVOLVE_NVARS], dytmp[EOB_EVOLVE_NVARS], ttmp; /* Additional buffer for post-Omegapeak ev */

  /* Alloc dynamics, 
     Set quick-access parameters dyn (be careful here) */
  Dynamics_alloc (&dyn, 0, "dyn"); 
  //Dynamics_set_params(dyn); //TODO: REMOVE always use EOBPars
  dyn->store = dyn->noflx = 0; /* Default: do not store vars, flux on */
  dyn->dt = EOBPars->dt;
  dyn->t_stop = EOBPars->ode_tmax;

  /* Alloc spin dynamics, if needed */
  if (use_spins == MODE_SPINS_GENERIC) {   
    
    DynamicsSpin_alloc(&spindyn, EOBPars->spin_dyn_size); 
    
    /* Set up a reference to spin dynamics in the dynamics structure */
    dyn->spins = spindyn;
    
  }

  /* Compute initial radius */
  double f0 = EOBPars->initial_frequency/time_unit_fact;
  
  EOBPars->f0 = f0;
  double r0;
  if (r_hyp != 0.) {
    // hyp case
    r0 = r_hyp;
  } else {
    r0 = eob_dyn_r0_eob(f0,dyn);
  }
  //double r0 = eob_dyn_r0_Kepler(f0); /* Kepler radius: no longer used */
    
  /* If f_min is too high fall back to a minimum acceptable initial radius */
  if ((EOBPars->model == MODEL_GIOTTO) && (r0 < TEOB_R0_THRESHOLD)) r0 = TEOB_R0_THRESHOLD;

  /* Saving initial radius */
  EOBPars->r0 = r0;
  
  const int chunk = EOBPars->size;
  int size = chunk; /* note: size can vary */

  /* Change size when using PA dynamics */
  if (use_postadiab_dyn) {
    size = EOBPars->postadiabatic_dynamics_size; 
    double rmin = EOBPars->postadiabatic_dynamics_rmin;
    size = floor(fabs(r0 - rmin)/EOBPars->postadiabatic_dynamics_dr) + 1;
    
    /* If initial radius is too close to PA limit then skip PA and go directly to ODE */
    if (size - 1 < POSTADIABATIC_NSTEP_MIN) {
      size = chunk;
      use_postadiab_dyn = 0;
    }
    EOBPars->size = size;
  }

  /* Allocating memory for dynamics and waveform */
  int use_nqc_tmp[KMAX];
  EOBPars->use_mode_lm_nqc_size = intersect_int(EOBPars->use_mode_lm,
                                                EOBPars->use_mode_lm_size,	
                                                EOBPars->use_mode_lm_nqc,
                                                EOBPars->use_mode_lm_nqc_size, 
                                                use_nqc_tmp);

  for(int i = 0; i < EOBPars->use_mode_lm_nqc_size; i++) 
    EOBPars->use_mode_lm_nqc[i] = use_nqc_tmp[i];

  Dynamics_push (&dyn, size); 
  Waveform_lm_alloc (&hlm, size, "hlm", EOBPars->use_mode_lm,EOBPars->use_mode_lm_size); 
  Waveform_lm_t_alloc (&hlm_t); 

  /* Integrate spin dynamics before EOB dyn if projecting */
  if (use_spins == MODE_SPINS_GENERIC && EOBPars->project_spins) {
    int tmp_nqc = EOBPars->nqc_coefs_flx;
    EOBPars->nqc_coefs_flx = NQC_FLX_NONE;
    if (eob_spin_dyn(spindyn, NULL, NULL, Pi * EOBPars->initial_frequency/time_unit_fact))
      errorexit("problem during spin dynamics");
    EOBPars->nqc_coefs_flx = tmp_nqc;

    spindyn->data[EOB_EVOLVE_SPIN_alp][0] = spindyn->data[EOB_EVOLVE_SPIN_alp][1];
    spindyn->data[EOB_EVOLVE_SPIN_gam][0] = spindyn->data[EOB_EVOLVE_SPIN_gam][1];

    if (EOBPars->spin_interp_domain == 0){
      if (VERBOSE) printf("Compute splines in t\n");
      for(int v=0; v < EOB_EVOLVE_SPIN_NVARS; v++)
        gsl_spline_init (spindyn->spline[v], spindyn->time, spindyn->data[v], spindyn->size); 
    } else {
      if (VERBOSE) printf("Compute splines in Momg\n");
      for(int v=0; v < EOB_EVOLVE_SPIN_NVARS; v++)
        gsl_spline_init (spindyn->spline[v], spindyn->data[EOB_EVOLVE_SPIN_Momg], spindyn->data[v], spindyn->size); 
    }  
  }

  /* NQC data */  
  NQCdata_alloc (&NQC); 
  eob_nqc_setcoefs(NQC);
  

  /* Final BH */
  if (!(EOBPars->binary == BINARY_BNS)) {
    // The aligned-spins Mbhf and abhf are computed in eob_set_params() inside
    // TEOBResumSPars.c. Here, we just update the values if the binary is precessing.

    if (use_spins == MODE_SPINS_GENERIC && EOBPars->project_spins) {   
      // (4.17) of https://arxiv.org/abs/2004.06503 
      EOBPars->abhf = PrecessingRemnantSpin(dyn);
    }
    if (VERBOSE) {
      PRSECTN("Final black hole");
      PRFORMd("BH_final_mass",EOBPars->Mbhf); 
      PRFORMd("BH_final_spin",EOBPars->abhf);
    }
  }
  /* Compute light-ring and LSO (if needed) */
  int check_status;
  //if (EOBPars->binary == BINARY_BNS) {
  if (use_tidal) { 
    /* Compute rLR_tidal for NNLO potential and without spin part */
    int tidal_tmp      = EOBPars->use_tidal;
    int spins_tmp      = EOBPars->use_spins;
    int use_fmode_tmp  = EOBPars->use_tidal_fmode_model;
    EOBPars->use_tidal = TIDES_NNLO; 
    EOBPars->use_spins = 0;
    EOBPars->use_tidal_fmode_model = 0;
    ROOTFINDER(check_status, eob_dyn_adiabLR(dyn, &(EOBPars->rLR_tidal)));
    if (check_status) {
      status = ERROR_ROOTFINDER;
      goto EXIT_POINT;
    }
    double LambdaAl2  = EOBPars->LambdaAl2;
    if( fabs(LambdaAl2) < TEOB_LAMBDA_TOL ) LambdaAl2 = 0.0;
    double LambdaBl2 = EOBPars->LambdaBl2;
    if( fabs(LambdaBl2) < TEOB_LAMBDA_TOL ) LambdaBl2 = 0.0;
    double q = EOBPars->q;
    /* Reset options */
    EOBPars->use_tidal = tidal_tmp;
    EOBPars->use_spins = spins_tmp;
    EOBPars->use_tidal_fmode_model = use_fmode_tmp;
    if (VERBOSE) PRFORMd("rLR_tidal",EOBPars->rLR_tidal);
    /* Set ODE stop to LR */
    EOBPars->ode_stop_radius = 1.01*EOBPars->rLR_tidal;    
  }
  if (EOBPars->compute_LR && !(use_tidal)) {
    //TODO: LR COMPUTATION IS CORRECT ONLY FOR NOSPIN. IMPLEMENT SPIN VERSION IN eob_dyn_adiabLSO()
    ROOTFINDER(check_status, eob_dyn_adiabLR(dyn, &(EOBPars->rLR)));
    if (check_status) {
      status = ERROR_ROOTFINDER;
      goto EXIT_POINT;
    }
    if (VERBOSE) PRFORMd("rLR",EOBPars->rLR);
  }
  if (EOBPars->compute_LSO) {
    ROOTFINDER(check_status, eob_dyn_adiabLSO(dyn, &(EOBPars->rLSO)));
    //Spin version
    //TODO: Decide what to do when LSO does not exist (large alignes spins) 
    //double pphiLSO = 0.;
    //ROOTFINDER(check_status, eob_dyn_LSO_s(dyn, &(dyn->rLSO), &pphiLSO));
    if (check_status) {
      status = ERROR_ROOTFINDER;
      goto EXIT_POINT;
    }
    if (VERBOSE) PRFORMd("rLSO",EOBPars->rLSO);
  }   
  /* Compute the dressing factors for the f-mode resonances at r0 */
  if ((EOBPars->use_tidal)&&(EOBPars->use_tidal_fmode_model))
    fmode_resonance_dressing_factors(r0, dyn);  

  /* Iteration index */
  int iter   = 0;  
  int pasize = 0;
 
  if (use_postadiab_dyn) {

    /* *****************************************
     * Post-adiabatic dynamics
     * *****************************************
     */

    if (VERBOSE) PRSECTN("Post-adiabatic dynamics");

    if (use_spins == MODE_SPINS_GENERIC && EOBPars->project_spins){
      if(pow(r0+0.1, -1.5) < spindyn->data[EOB_EVOLVE_SPIN_Momg][0]){
        if (VERBOSE) printf("integrate backwards\n");
        eob_spin_dyn_integrate_backwards(spindyn, NULL, NULL, pow(r0+0.1, -1.5));
      }
    }

    /* Calculate dynamics */
    eob_dyn_Npostadiabatic(dyn, r0, spindyn); 

    /* Calculate waveform */
    for (int i = 0; i < size; i++) 
      hlm->time[i] = dyn->time[i];
    
    dyn->store = dyn->noflx = 1;

    for (int i = 0; i < size; i++) {
      dyn->y[EOB_EVOLVE_RAD]    = dyn->data[EOB_RAD][i];
      dyn->y[EOB_EVOLVE_PHI]    = dyn->data[EOB_PHI][i];
      dyn->y[EOB_EVOLVE_PRSTAR] = dyn->data[EOB_PRSTAR][i]; 
      dyn->y[EOB_EVOLVE_PPHI]   = dyn->data[EOB_PPHI][i];

      if (use_spins == MODE_SPINS_GENERIC && EOBPars->project_spins){
        double SA, SB;
        eob_spin_dyn_Sproj_interp(spindyn, dyn->data[EOB_MOMG][i], &SA, &SB, NULL, NULL, NULL, NULL, 1);
        EOBPars->chi1 = SA/SQ(EOBPars->X1);
        EOBPars->chi2 = SB/SQ(EOBPars->X2);
        set_spin_vars(EOBPars->X1,EOBPars->X2,EOBPars->chi1,EOBPars->chi2,&EOBPars->S1, &EOBPars->S2,&EOBPars->a1, &EOBPars->a2,&EOBPars->aK, &EOBPars->aK2,&EOBPars->S, &EOBPars->Sstar);
      }
      
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
    
    /* Prepare for evolution */
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

    /* Compute the initial conditions */
    if (eob_dyn_ic(r0, dyn, dyn->y0)) {
      status = ERROR_INITIAL_CONDITIONS;
      goto EXIT_POINT;
    }
      
    /* check that hyperbolic orbits ic work */
    // TODO: check if this workaround is needed  
    if ((r_hyp != 0.) && (dyn->y0[EOB_ID_PRSTAR] == 0.)) {
      status = ERROR_INITIAL_CONDITIONS;
      goto EXIT_POINT;
    }
    
    /* Se arrays with initial conditions */
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

    /* Check: print warning if initial separation < 10 */
    if (dyn->r < 10.) {
      if (VERBOSE) printf("WARNING: Initial separation < 10\n");
    }

    if (use_spins == MODE_SPINS_GENERIC && EOBPars->project_spins){
      if(dyn->data[EOB_MOMG][0] < spindyn->data[EOB_EVOLVE_SPIN_Momg][0]){
        eob_spin_dyn_integrate_backwards(spindyn, NULL, NULL, dyn->data[EOB_MOMG][0]);
      }
      /* if we are projecting the spins, we also need to update the initial chi1z and chi2z to be consistent*/
      double SA, SB; // projections of the spin parallel to hatL
      if (EOBPars->spin_interp_domain == 0)
        eob_spin_dyn_Sproj_interp(spindyn, 0., &SA, &SB, NULL, NULL, NULL, NULL, EOBPars->spin_interp_domain);
      else 
        eob_spin_dyn_Sproj_interp(spindyn, dyn->data[EOB_MOMG][0], &SA, &SB, NULL, NULL, NULL, NULL, EOBPars->spin_interp_domain);
      double XA = EOBPars->X1;
      double XB = EOBPars->X2;
      EOBPars->chi1 = SA /SQ(XA);
      EOBPars->chi2 = SB /SQ(XB);
      set_spin_vars(XA,XB,EOBPars->chi1,EOBPars->chi2, 
		    &EOBPars->S1, &EOBPars->S2,
		    &EOBPars->a1, &EOBPars->a2,
		    &EOBPars->aK, &EOBPars->aK2,
		    &EOBPars->S, &EOBPars->Sstar);
    }    
    
    /* Waveform computation at t = 0 
	Needs a r.h.s. evaluation for some vars (no flux) */
    dyn->store = dyn->noflx = 1;
    p_eob_dyn_rhs(dyn->t, dyn->y, dyn->dy, dyn); 
    dyn->store = dyn->noflx = 0;
    eob_wav_hlm(dyn, hlm_t); 
    
    /* Append waveform to arrays */
    hlm->time[0] = 0.;
    for (int k = 0; k < KMAX; k++) {
      if(hlm->kmask[k]) hlm->ampli[k][0] = hlm_t->ampli[k];
    }
    for (int k = 0; k < KMAX; k++) {
      if(hlm->kmask[k]) hlm->phase[k][0] = hlm_t->phase[k]; 
    }

    /* Prepare for evolution */
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

  /* Initialize ODE system solver */
  dyn->t_stop = EOBPars->ode_tmax * time_unit_fact;
  EOBPars->ode_tmax      = dyn->t_stop;
  dyn->ode_stop          = false;
  dyn->ode_stop_MOmgpeak = false;
  dyn->ode_stop_radius   = false;
  double rstop     = EOBPars->ode_stop_radius;
  /* Avoiding useless computations for hyperbolic cases */
  if (EOBPars->model == MODEL_DALI) {
    double r_hor = horizon_radius(EOBPars->nu)+0.2;
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

  /* GSL routines for ODE integration
     https://www.gnu.org/software/gsl/doc/html/ode-initval.html
     http://www.csse.uwa.edu.au/programming/gsl-1.0/gsl-ref_24.html
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
  const double deltat_tuned_mrg = get_mrg_timestop(q, chi1, chi2);
  
  /* Solve ODE */
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
      if (dyn->r > EOBPars->rLSO) { 
	GSLSTATUS = gsl_odeiv2_evolve_apply (e, c, s, &sys, &dyn->t, dyn->t_stop, &dyn->dt, dyn->y);
      } else {
        /* ... uniform afterwards */
        dyn->dt = dt_tuned_mrg;
        dyn->ti = dyn->t + dyn->dt;
        GSLSTATUS = gsl_odeiv2_driver_apply (d, &dyn->t, dyn->ti, dyn->y);
      }
    }
    
    /* Unpack data */
    dyn->r      = dyn->y[EOB_EVOLVE_RAD];
    dyn->phi    = dyn->y[EOB_EVOLVE_PHI];
    dyn->prstar = dyn->y[EOB_EVOLVE_PRSTAR];
    dyn->pphi   = dyn->y[EOB_EVOLVE_PPHI];
    
    /* Waveform computation 
	Needs a r.h.s. evaluation for some vars (but no flux) */
    dyn->store = dyn->noflx = 1;
    p_eob_dyn_rhs(dyn->t, dyn->y, dyn->dy, dyn); 
    dyn->store = dyn->noflx = 0;
    eob_wav_hlm(dyn, hlm_t); 

    if (use_spins) {
      dyn->MOmg = dyn->Omg_orb;
    } else {
      dyn->MOmg = dyn->Omg;
    }

    /* Update spins parallel to L, if required */
    if (use_spins == MODE_SPINS_GENERIC && EOBPars->project_spins) {
      double SA, SB; // projections of the spin parallel to hatL
      if(EOBPars->spin_interp_domain==0){
        eob_spin_dyn_Sproj_interp(dyn->spins, dyn->t, &SA, &SB, NULL, 
				  NULL, NULL, NULL, 
				  EOBPars->spin_interp_domain);
      } else {
        eob_spin_dyn_Sproj_interp(dyn->spins, dyn->Omg, &SA, &SB, NULL, 
				  NULL, NULL, NULL, 
				  EOBPars->spin_interp_domain);
      }
      const double XA = EOBPars->X1;
      const double XB = EOBPars->X2;
      
      EOBPars->chi1 = SA /SQ(XA);
      EOBPars->chi2 = SB /SQ(XB);
      
      set_spin_vars(XA,XB,EOBPars->chi1,EOBPars->chi2, 
		    &EOBPars->S1, &EOBPars->S2,
		    &EOBPars->a1, &EOBPars->a2,
		    &EOBPars->aK, &EOBPars->aK2,
		    &EOBPars->S, &EOBPars->Sstar);
    }
	
    /* Check for failures ... */
    if (EOBPars->model == MODEL_DALI) {
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
    
    /* Checking whether the dynamics produces NaN values
	  this can happen if radius r becomes too small */
    if (!(isfinite(dyn->r))) {
      printf("%.1f\t%.3f\t%.3f\n", q, chi1, chi2);	
      /* errorexit("ODE solver returned NaN radius.\n"); */
      printf("ODE solver returned NaN radius.\n");
      status = ERROR_ODEINT;
      goto EXIT_POINT;
    }
    
    /* Update size and push arrays (if needed) */
    if (iter==size) {
      /* if (DEBUG)  printf("Push memory\n"); */ 
      size += chunk;
      EOBPars->size = size;
      Waveform_lm_push (&hlm, size);
      Dynamics_push (&dyn, size);
    }
    
    /* Append waveform and dynamics to arrays */
    hlm->time[iter] = hlm_t->time;
    for (int k = 0; k < KMAX; k++) {
      if(hlm->kmask[k]){
        hlm->ampli[k][iter] = hlm_t->ampli[k];
        hlm->phase[k][iter] = hlm_t->phase[k]; 
      }
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

    /* Stop integration if reached max time */    
    if (dyn->t > dyn->t_stop) {
      if (VERBOSE) printf("Stop: Max integration time reached.\n");
      dyn->ode_stop = true;
    }
    
    /* Stop integration at given radius (if rstop >= 0) */    
    if ((dyn->ode_stop_radius) && (dyn->r < rstop) ) {
      if (VERBOSE) printf("Stop: radius reached.\n");
      dyn->ode_stop = true;
    }

    /* If it could be a scattering (E > 1.), stop integration at large radii (if rstop >= 0) */
    if ((dyn->ode_stop_radius) && (dyn->r > 5.*r0) && r_hyp != 0. && dyn->E > 1.) {
      if (VERBOSE) printf("Stop: maximum radius reached.\n");
      dyn->ode_stop = true;
    }

    /* Stop integration if r is unchanged */
    double r_prev = dyn->data[EOB_RAD][iter-1];
    if (dyn->r == r_prev && !dyn->noflx) {
      if (VERBOSE) printf("Stop: radius unchanged.\n");
      dyn->ode_stop = true;
    }

    /* Check when to break the computation
	  find peak of omega curve and continue for 2M */
    if ((dyn->ode_stop_MOmgpeak == false) && (EOBPars->model != MODEL_DALI)) {
      /* Before the Omega_orb peak */      
      if (dyn->MOmg < dyn->MOmg_prev) {
	/* This is the first step after the peak
	   Set things for uniform tstep evolution */
        dyn->tMOmgpeak = dyn->t;
        dyn->ode_stop_MOmgpeak = true;

        dyn->dt = MIN(dyn->dt, dt_tuned_mrg); 
        dyn->t_stop = dyn->t + deltat_tuned_mrg;
        
        if(EOBPars->ode_stop_after_peak){
          dyn->ode_stop = true;
          iter--;
        }
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
  
  /* Free ODE system solver */
  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
  gsl_odeiv2_driver_free (d);

  /* Update waveform and dynamics size 
      resize to actual size */
  size = iter+1;
  EOBPars->size = size;

  /* Check: is the dynamics long enough? */
  if (size < 10){
    // CHECKME: 10 points is somewhat arbitrary
    printf("ERROR(TEOBResumS): ODE dynamics size < 10\n");
    status = ERROR_ODEINT;
    goto EXIT_POINT;
  }

  Waveform_lm_push (&hlm, size);
  Dynamics_push (&dyn, size);

  /* Find peak of Omega */
  int index_pk = dyn->size-1;
  double Omega_pk = dyn->data[EOB_OMGORB][index_pk];
  for (int j = dyn->size-2; j>0; j--) {
    if (dyn->data[EOB_OMGORB][j] < Omega_pk) 
      break;
    index_pk = j;
    Omega_pk = dyn->data[EOB_OMGORB][j]; 
  }

  double tOmg_pk;
  if (dyn->size-1 - index_pk < 2) {
    /* If peak is too close to the end, 
       use the point itself */
    tOmg_pk = dyn->time[index_pk];
  } else {
    double *t_ptr     = &dyn->time[index_pk-2];
    double *Omega_ptr = &dyn->data[EOB_OMGORB][index_pk-2];
    tOmg_pk    = find_max_grid(t_ptr, Omega_ptr);
  }
  dyn->tOmg_pk = tOmg_pk;
  
 END_ODE_EVOLUTION:;
 
  /* Unwrap phase for higher modes */
  if ((EOBPars->use_flm == USEFLM_HM) || (EOBPars->use_flm == USEFLM_HM_4PN22) || (ecc != 0.)) {
    for (int k = 0; k < KMAX; k++) {
      if(hlm->kmask[k]){
	      unwrap_HM(hlm->phase[k],size);
      }
    }
  }

#if (DEBUG) 
  // Output wave and dynamics 
  if(EOBPars->output_multipoles) {
    strcat(hlm->name,"_insplunge");
    Waveform_lm_output (hlm);
    Waveform_lm_output_reim (hlm);
  }
  if (EOBPars->output_dynamics) {
    Dynamics_output(dyn);
  }
#endif
  
  /* Precessing BNS (or BBH with no RD) + EOB flux: */
  if( (EOBPars->binary == BINARY_BNS || dyn->data[EOB_RAD][size-1] > 3.) && use_spins == MODE_SPINS_GENERIC && !(EOBPars->project_spins)){
    if (eob_spin_dyn(spindyn, dyn, hlm, Pi*EOBPars->initial_frequency/time_unit_fact))
      errorexit("problem during spin dynamics");
    spindyn->data[EOB_EVOLVE_SPIN_alp][0] = spindyn->data[EOB_EVOLVE_SPIN_alp][1];
    spindyn->data[EOB_EVOLVE_SPIN_gam][0] = spindyn->data[EOB_EVOLVE_SPIN_gam][1];
    
    if(dyn->data[EOB_MOMG][0] < spindyn->data[EOB_EVOLVE_SPIN_Momg][0] && EOBPars->spin_flx != SPIN_FLX_EOB)
      eob_spin_dyn_integrate_backwards(spindyn, dyn, hlm, dyn->data[EOB_MOMG][0]);  
  }

  if (!(EOBPars->binary == BINARY_BNS) && (dyn->data[EOB_RAD][size-1] < 3.)) {
    
    /* ******************************************************
     * Following is for merging BBH and BHNS: NQC & Ringdown
     * ******************************************************
     */
    
    /* Over-writing waveform in the eccentric case - adding sigmoid */
    if ((EOBPars->nqc_coefs_hlm != NQC_HLM_NONE) && (EOBPars->model == MODEL_DALI)) {
   
      for (int i = 0; i < size; i++) {
        dyn->store = 1;
        dyn->t = dyn->time[i];
        dyn->y[EOB_EVOLVE_PHI]    = dyn->data[EOB_PHI][i];
        dyn->y[EOB_EVOLVE_RAD]    = dyn->data[EOB_RAD][i];
        dyn->y[EOB_EVOLVE_PPHI]   = dyn->data[EOB_PPHI][i];
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

    
    /* This is a BBH run.
       NQC and ringdown attachment currently assume uniform grids.
       Do we need to interpolate ? */
    int merger_interp = 1; /* In general, yes ... */
    if ((ode_tstep != ODE_TSTEP_ADAPTIVE) && (EOBPars->use_flm != USEFLM_HM) && (EOBPars->use_flm != USEFLM_HM_4PN22)) merger_interp = 0; /* ... except if merger is covered by uniform tstep */

    /* NQC and ringdown attachment is done around merger 
	using auxiliary variables defined around [tmin,tmax] 
	Recall that parameters are NOT stored into these auxiliary vars */
    
    if (merger_interp) {
      
      /* Extract the waveform and dynamics around merger */
      // Note : 20M might be too little
      double tmin = hlm->time[size-1] - 20; /* Use last 20M points */
      const double tmax = hlm->time[size-1] +  2*dt; /* Make sure to use or get last point */

      /* For hyperbolic or eccentric orbits, find time when r = 10 */
      if (EOBPars->model == MODEL_DALI) {
	      for (int j = dyn->size-1; j>0; j--) {
          tmin = dyn->time[j];
          if (dyn->data[EOB_RAD][j] > 10.)
	          break;
	      }
      }
      
      // Check that the time interval chosen contains the peak of Omega
      if (dyn->tOmg_pk < tmin || dyn->tOmg_pk > tmax)
        errorexit("The peak orbital frequency is not contained in the NQC/RD attachment region.");
      
      /* The following routines alloc memory for the *_mrg ptrs */
      Waveform_lm_extract (hlm, tmin, tmax, &hlm_mrg, "hlm_mrg");
      Dynamics_extract    (dyn, tmin, tmax, &dyn_mrg, "dyn_mrg");
      
      /* Build uniform grid of width dt and alloc tmp memory */
      double dt_merger_interp;
      if ( (EOBPars->use_flm == USEFLM_HM) || (EOBPars->use_flm == USEFLM_HM_4PN22)) {
	      dt_merger_interp = 0.5;
      } else {
	      dt_merger_interp = MIN(EOBPars->dt_merger_interp, dyn->dt);
      }
      double tstart_mrg = dyn->tOmg_pk - 8.;
      while (tstart_mrg < hlm_mrg->time[0]) /** Make sure it does not extrapolate */
	      tstart_mrg += dt_merger_interp;
      const long int size_mrg = get_uniform_size(hlm_mrg->time[hlm_mrg->size-1], tstart_mrg, dt_merger_interp);
      
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
      
      /**  Interpolate mrg on uniform grid */
      
      /* Interp Waveform */ 
      Waveform_lm_interp (hlm_mrg, size_mrg, tstart_mrg, dt_merger_interp, "hlm_mrg_interp");
      
      /* Interp Dynamics */
      Dynamics_interp (dyn_mrg, size_mrg, tstart_mrg, dt_merger_interp, "dyn_mrg_interp");	       
      
      
    } /* End of merger interp */
    
    
    if ((EOBPars->nqc_coefs_hlm == NQC_HLM_COMPUTE)) {
      
      /** BBH : compute and add NQC */
    
      if (VERBOSE) PRSECTN("NQC Calculation");
      
      if (merger_interp) { 
	
        /* Compute NQC only around merger, 
	   add to both merger and full waveform */
        Waveform_lm_alloc (&hlm_nqc, hlm_mrg->size, "hlm_nqc", EOBPars->use_mode_lm, EOBPars->use_mode_lm_size); 
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
        Waveform_lm_alloc (&hlm_nqc, size, "hlm_nqc", EOBPars->use_mode_lm,EOBPars->use_mode_lm_size); 
        eob_wav_hlmNQC_find_a1a2a3(dyn, hlm, hlm_nqc);
      }
      
#if (DEBUG) 
      if (EOBPars->output_nqc)  {
	Waveform_lm_output (hlm_nqc);
	Waveform_lm_output (hlm_mrg);
      }
      if (EOBPars->output_multipoles) {
	strcat(hlm->name,"_nqc");      
	Waveform_lm_output_reim (hlm);
      }
#endif
      
      Waveform_lm_free (hlm_nqc);
      
    }  /* NQC_HLM_COMPUTE */

    /* Precessing BBH + EOB flux: 
    */
    if(use_spins == MODE_SPINS_GENERIC && !(EOBPars->project_spins)){
      if (eob_spin_dyn(spindyn, dyn, hlm, Pi*EOBPars->initial_frequency/time_unit_fact))
        errorexit("problem during spin dynamics");
      spindyn->data[EOB_EVOLVE_SPIN_alp][0] = spindyn->data[EOB_EVOLVE_SPIN_alp][1];
      spindyn->data[EOB_EVOLVE_SPIN_gam][0] = spindyn->data[EOB_EVOLVE_SPIN_gam][1];
      
      if(dyn->data[EOB_MOMG][0] < spindyn->data[EOB_EVOLVE_SPIN_Momg][0] && EOBPars->spin_flx != SPIN_FLX_EOB)
        eob_spin_dyn_integrate_backwards(spindyn, dyn, hlm, dyn->data[EOB_MOMG][0]);  

      /* final state */
      EOBPars->abhf = PrecessingRemnantSpin(dyn); 
    }

    /* BBH : add Ringdown */
    
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
    if (eob_wav_ringdown(dyn, hlm)){
      printf("ERROR(TEOBResumS): %s\n",eob_error_msg[ERROR_RINGDOWN]);
      status = ERROR_RINGDOWN;
      goto EXIT_POINT;
    }
    
  } /* End of BBH section */

#if (DEBUG) 
  if (EOBPars->output_multipoles) {
    strcat(hlm->name,"_ringdown");
    Waveform_lm_output (hlm);
    Waveform_lm_output_reim (hlm);
  }
#endif
      
  /* *****************************************
   * Compute h+, hx 
   * *****************************************
   */
  
  /* Set tc */
  /* Note: tc (=merger time) is correct for BNS, not for BBH */
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
  const double phi  = Pi/2.0 - EOBPars->coalescence_angle; 
  const double iota = EOBPars->inclination;

  /* Computation of (h+,hx) */
  
  if (EOBPars->domain == DOMAIN_TD) {

    /* TIME DOMAIN */
    
    /* Find the time shift to have merger at 0*/
    if(EOBPars->time_shift_TD)
      EOBPars->tc = time_shift_mrg_to_0(hlm);

    if (EOBPars->interp_uniform_grid) {
      /* Interp to uniform grid the multipoles before hpc computation */
      const double dt_interp = EOBPars->dt_interp;
      const long int size_interp = get_uniform_size(hlm->time[size-1], hlm->time[0], dt_interp); 
      printf("Interpolating to uniform grid with dt = %e\n",dt_interp);
      Waveform_lm_interp (hlm, size_interp, hlm->time[0], dt_interp, "hlm_interp");
      printf("Interpolateds\n");
      size = size_interp;
    }

    /* Alloc memory for (h+,hx) */
    Waveform_alloc (hpc, size, "waveform");   
  
    /* h+, hx */  
    if (use_spins == MODE_SPINS_GENERIC){
      if (VERBOSE) PRSECTN("Twisting");
      Waveform_lm_alloc (&hTlm, size, "hTlm",EOBPars->use_mode_lm_inertial,EOBPars->use_mode_lm_inertial_size); 
      Waveform_lm_alloc (&hTlm_neg, size, "hTlm_neg",EOBPars->use_mode_lm_inertial,EOBPars->use_mode_lm_inertial_size); 
      Waveform_lm_alloc (&hTl0, size, "hTl0",EOBPars->use_mode_lm_inertial,EOBPars->use_mode_lm_inertial_size); 
      twist_hlm_TD(dyn, hlm, dyn->spins, 1, hTlm, hTlm_neg, hTl0);
      compute_hpc(hTlm, hTlm_neg, hTl0, nu, M, distance, amplitude_prefactor, phi, iota, *hpc);
    } else
      compute_hpc(hlm, NULL, NULL, nu, M, distance, amplitude_prefactor, phi, iota, *hpc);
         
    /* time-shift the TD waveforms */
    if(EOBPars->time_shift_TD) 
      time_shift_TD((*hpc)->time, EOBPars->tc/time_unit_fact, (*hpc)->size);

    /* Enforce the LAL tetrad convention: hlm --> -1*hlm, 
       implying that h_{+,x} -->-1*h_{+,x} 
    */
    if (EOBPars->lal_tetrad_conventions){
      lal_tetrad_conventions(*hpc);
    }
  } else {
    
    /* FREQUENCY DOMAIN */
    
    /* Frequency domain multipolar waveform */
    WaveformFD_lm_alloc (&hflm, size, "hflm");
    
    /* Calculate the SPA for the multipolar waveform */
    if (VERBOSE) PRSECTN("SPA");
    SPA(hlm, hflm);
    /* The SPA performs an interpolation, 
       need to update the size */
    size = hflm->size;
    
    /* Alloc memory for (h+,hx) */
    WaveformFD_alloc (hfpc, size, "waveform");   

    /* h+, hx */  
    if (use_spins == MODE_SPINS_GENERIC){
      
      if (VERBOSE) PRSECTN("Twist");
      twist_hlm_FD(hflm, dyn->spins, M, amplitude_prefactor, phi, iota, *hfpc); //this also computes hpc
    } else {
      compute_hpc_FD(hflm, nu, M, distance, amplitude_prefactor, phi, iota, *hfpc);
    }

    /* time-shift the FD waveforms */
    if(EOBPars->time_shift_FD) time_shift_FD(*hfpc, EOBPars->tc);
  }

  /* *****************************************
   * Output
   * *****************************************
   */

  if (output) {
    // hpc
    if (EOBPars->output_hpc)
      if (EOBPars->domain == DOMAIN_TD) Waveform_output (*hpc); 
      else                              WaveformFD_output (*hfpc); 
    // modes
    if (EOBPars->output_multipoles) {
      if (EOBPars->domain == DOMAIN_TD) {
	Waveform_lm_output (hlm); 
	Waveform_lm_output_reim (hlm);
	if (use_spins == MODE_SPINS_GENERIC){
	  Waveform_lm_output (hTlm);
	  Waveform_lm_output (hTlm_neg);
	  Waveform_lm_output (hTl0);
	}
      } else {
	WaveformFD_lm_output (hflm);
	if (use_spins == MODE_SPINS_GENERIC)
	  WaveformFD_lm_output (hfTlm);
      }
    }
    // dynamics
    if (EOBPars->output_dynamics) {
      if (EOBPars->interp_uniform_grid) {
	/* Interp to uniform grid the dynamics, 
	   rem the dyn size can be different from wf size */
	const double dt_interp_dyn = EOBPars->dt_interp;
	const int size_interp_dyn = get_uniform_size(dyn->time[dyn->size-1], dyn->time[0], dt_interp_dyn);
	Dynamics_interp (dyn, size_interp_dyn, dyn->time[0], dt_interp_dyn, "dyn_interp");  
      }
      Dynamics_output(dyn);
      if (use_spins == MODE_SPINS_GENERIC)
	DynamicsSpin_output (spindyn);
    }
  }


  /* *****************************************
   * Finalize 
   * *****************************************
   */
 EXIT_POINT:;

  if (status) {
    /* There was an error, return zeros waveform */
    if (EOBPars->domain == DOMAIN_TD) {
      //SB: the size here needs to be fixed to the required sampling frequency.
      //    if not, the code jumps here and size is still the one from default...
      Waveform_alloc (hpc, size, "waveform");
      Waveform_lm_free (hlm);
      Waveform_lm_alloc (&hlm, size, "hlm", EOBPars->use_mode_lm,EOBPars->use_mode_lm_size);
    } else  {                             
      const int interp_fd_size = get_uniform_size(EOBPars->initial_frequency, EOBPars->initial_frequency, EOBPars->df);
      WaveformFD_alloc (hfpc, interp_fd_size, "waveform_fd");
      WaveformFD_lm_alloc (hfmodes, interp_fd_size, "hlm_fd");
    }
  }
  
  *hmodes = hlm;        /* do not free these! */
  *hfmodes = hflm;      /* do not free these! */
  *hTmodes = hTlm;      /* do not free these! */
  *hTmmodes = hTlm_neg; /* do not free these! */
  *hT0modes = hTl0;
  *hfTmodes = hfTlm;    /* do not free these! */
  *dynf = dyn;          /* do not free these! */

#ifdef _OPENMP
  openmp_free(); 
#endif
  
  /* Free memory */
  /* Dynamics_free (dyn);       */
  if (use_spins == MODE_SPINS_GENERIC)
    DynamicsSpin_free (spindyn);  

  Waveform_lm_t_free (hlm_t);
  NQCdata_free (NQC);

  return status;
}

