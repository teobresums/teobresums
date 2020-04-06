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
 * @file TEOBResumS..h
 * @brief Header file of the TEOBResumS C code
 *
 * This file contains all the macros, typdef, and routine prototype.
 * Doxygen documentation should go here.
 *
 */

#define _GNU_SOURCE /* This is just for getline() */
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <libconfig.h> /* library to manage parameters */
#include <complex.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#ifdef _OPENMP
#include <omp.h>
#define omp ignore
#endif

/** Following macros can be set during compilation for special feats */
/* NOTE: below are the defaults; macros should be changed exclusively from makefile */

#ifndef USETIMERS
#define USETIMERS (0) /* This requires OMP otherwise do nothing! */
#endif

#ifndef VERBOSE
#define VERBOSE (0) /* verbose mode is off by default */
#endif

#ifndef POSTPOSTCIRCULAR
#define POSTPOSTCIRCULAR (1) /* use post-post-circular initial conditions by default */
#endif

#ifndef EXCLUDESPINSPINTIDES
#define EXCLUDESPINSPINTIDES (0) /* use tidally deformed centr. radius with self-spin and tides by default */
#endif

#ifndef USEGRAVITOMAGNETICTERMS
#define USEGRAVITOMAGNETICTERMS (1) /* use gravitomagnetic terms in tidal potential and waveform */
#endif

#ifndef USEBTIDALPOTENTIAL
#define USEBTIDALPOTENTIAL (1) /* add B LO tidal potential */
#endif

#ifndef USERK45
#define USERK45 (1) /* use GSL rkf45 instead of rk8pd */
#endif

#ifndef DEBUG 
#define DEBUG (0) /* global debug option */ 
#endif

/** Macros */
#define ERROR 1 /** generic error int */
#define STRLEN 1024 /** Standard string length */
#define TEOBResumS_Info "TEOBResumS code (C) 2017\n"
#define TEOBResumS_Usage(x) {printf("%sUSAGE:\t%s <parfile>\n", TEOBResumS_Info, x);} 
#define SIGN(x,y) ((y) >= 0.0 ? fabs(x) : -fabs(x)) 
#define typeof __typeof__
#define MAX(a,b)				\
  ({ typeof (a) _a = (a);			\
    typeof (b) _b = (b);			\
    _a > _b ? _a : _b; })       
#define MIN(a,b)				\
  ({ typeof (a) _a = (a);			\
    typeof (b) _b = (b);			\
    _a < _b ? _a : _b; })       
#define MAX3(a,b,c) (((a) > (b)) ? MAX(a,c) : MAX(b,c))
#define MIN3(a,b,c) (((a) < (b)) ? MIN(a,c) : MIN(b,c))
#define SQ(a) ((a)*(a))
#define DISZERO(a,eps) (fabs((a))<(eps)) /** Is double = 0 to precision eps ?*/
#define DEQUAL(a,b,eps) (fabs((a)-(b))<(eps)) /** double compare */
#define DUNEQUAL(a,b,eps) (fabs((a)-(b))>(eps))
#define STREQUAL(s,t) ((strcmp((s),(t))==0)) /** string compare */  
#define SWAPTRS(a,b)				\
  ({						\
    typeof(a) temp;				\
    temp = a;					\
    a = b;					\
    b = temp;					\
  })
/* helpers for formatted printf */
#define PRSECTN(s) {printf("#\n# %s\n#\n",s);} /* Print section */
#define PRFORMd(s,x) {printf("%-40s = %.16e\n", s,x);} /* Print double */
#define PRFORMi(s,x) {printf("%-40s = %d\n", s,x);} /* Print int */
#define PRWARN(s) {printf("# WARNING: %s\n",s);} 
#define INT2YESNO(i)((i)?"yes":"no")
#define YESNO2INT(s)((strcmp(s,"no")==0)?0:1)
/* helpers for debug */
#define DBGPR(s) printf("DEBUG: %s\n",s);
#define DBGSTOP errorexit("DEBUG: STOP");

/* Useful constants */
#define Pi (3.1415926535897932384626433832795028)
#define TwoPi (6.283185307179586476925286766559)
#define Sqrt2 (1.41421356237309504880168872420969808)
#define Sqrt3 (1.73205080756887729352744634150587237)
#define ooSqrt2 (0.707106781186547524400844362104849039284836)
#define Log1  (0.)
#define Log2  (0.693147180559945309417232)
#define Log3  (1.09861228866810969139525)
#define Log4  (1.38629436111989061883446)
#define Log5  (1.60943791243410037460076)
#define Log6  (1.79175946922805500081248)
#define Log7  (1.94591014905531330510535)
#define MSUN_M (1.476625061404649406193430731479084713e3)  /* G/c^3 */
#define MSUN_S (4.925491025543575903411922162094833998e-6) /* G/c^3 */
#define MPC_M  (3.085677581491367278913937957796471611e22) 
#define EulerGamma (0.5772156649015328606065121)
#define EulerGamma_Log2 (1.27036284546147817002374) /** EulerGamma + Log2 */
#define TEOB_LAMBDA_TOL (1.0)   /* Minimum tidal Lambda_2 value */
#define POSTADIABATIC_DR (0.1) /* PA step dr = 0.1 */
#define TEOB_R0_THRESHOLD (14) /* PA minimum tolerated radius */
#define POSTADIABATIC_NSTEP_MIN (10) /* Minimum requires PA steps, any less than this, the code switches off PA */

/** Simple/generic error handler */
enum{OK,
     ERROR_OUTOFMEM,
     ERROR_FILEOPEN,
     ERROR_MKDIR,
     ERROR_ROOTFINDER,
     ERROR_ODEINT,
     NERROR
};
static const char* eob_error_msg[] = {
  "ok",
  "out of memory", "error opening file", "error while making directory",
  "root finder failed.", "ODE solver failed."
};

/** Index list of EOB evolved variables */
enum{
  EOB_EVOLVE_RAD, 
  EOB_EVOLVE_PHI,
  EOB_EVOLVE_PRSTAR,
  EOB_EVOLVE_PPHI,
  EOB_EVOLVE_NVARS
};
static const char* eob_evolve_var[] = {"r","phi","Prstar","Pphi"};

/** Index list of EOB variables for initial data */
enum{
  EOB_ID_RAD,
  EOB_ID_PHI,
  EOB_ID_PPHI,
  EOB_ID_PRSTAR,
  EOB_ID_PR,
  EOB_ID_J,
  EOB_ID_E0,
  EOB_ID_OMGJ,
  EOB_ID_NVARS
};
static const char* eob_id_var[] = {"r","phi","Pphi","Prstar","Pr","j","E0","Omega"};

/** Index list of EOB dynamical variables (to be stored in arrays) */ 
enum{
  EOB_RAD, 
  EOB_PHI,
  EOB_PPHI,
  EOB_MOMG,
  EOB_DDOTR,
  EOB_PRSTAR,
  EOB_OMGORB,
  EOB_E0,
  EOB_DYNAMICS_NVARS
};
static const char* eob_var[] = {"r","phi","Pphi","MOmega","ddor","Prstar","MOmega_orb","E"};

#define KMAX (35) /** Multipolar linear index, max value */
#define PMTERMS_eps (1) /** Switch on Fujita-Iyer point-mass terms. This is hard-coded here */

/** List for default parameter options */
enum{
  DEFAULT_PARS_BBH,
  DEFAULT_PARS_BNS,
  DEFAULT_PARS_BHNS,
  DEFAULT_PARS_NOPT
};

/** List of options for tidal potential */
enum{
  TIDES_OFF,  /* = 0 , keep first to allow syntax: if(use_tidal) { ... */
  TIDES_NNLO, 
  TIDES_TEOBRESUM,
  TIDES_TEOBRESUM3,
  TIDES_TEOBRESUM_BHNS,
  TIDES_NOPT
};
static const char* const tides_opt[] = {"no","NNLO","TEOBRESUM", "TEOBRESUM3", "TEOBRESUM_BHNS","undefined"};

/** List of options for the gravitomagnetic tidal potential */
enum{
  TIDES_GM_OFF, /* = 0 , keep first to allow syntax: if(use_tidal_gravitomagnetic) { ... */
  TIDES_GM_PN, 
  TIDES_GM_GSF,
  TIDES_GM_NOPT
};
static const char* const tides_gravitomagnetic_opt[] = {"no","PN","GSF","undefined"};

/** List of options for centrifugal radius */
enum{
  CENTRAD_LO,
  CENTRAD_NLO,
  CENTRAD_NNLO,
  CENTRAD_NNLOS4,
  CENTRAD_NOSPIN,
  CENTRAD_NOTIDES,
  CENTRAD_NOPT
};
static const char* const centrifugal_radius_opt[] = {"LO", "NLO", "NNLO", "NNLOS4", "NOSPIN", "NOTIDES"};

/** List of options for flm amplitudes */
enum{
  USEFLM_SSLO,
  USEFLM_SSNLO,
  USEFLM_SSNNLO,
  USEFLM_HM,
  USEFLM_NOPT
};
static const char* const use_flm_opt[] = {"SSLO", "SSNLO", "SSNNLO", "HM"};

/** List of options for ODE timestepping */
enum{
  ODE_TSTEP_UNIFORM, 
  ODE_TSTEP_ADAPTIVE,
  ODE_TSTEP_ADAPTIVE_UNIFORM_AFTER_LSO,
  ODE_TSTEP_NOPT
};
static const char* const ode_tstep_opt[] = {"uniform","adaptive","adaptive+uniform_after_LSO","undefined"};

/** List of options for interp_uniform_grid */
enum{
  INTERP_UNIFORM_GRID_OFF,
  INTERP_UNIFORM_GRID_HPC,
  INTERP_UNIFORM_GRID_HLM,
  INTERP_UNIFORM_GRID_NOPT
};

/** Error handler for root finders */
enum{ 
  ROOT_ERRORS_NO,
  ROOT_ERRORS_BRACKET,
  ROOT_ERRORS_MAXITS,
  ROOT_ERRORS_NOSUCC,
  ROOT_ERRORS
};
static const char* const root_errors[] = {"none","root is not bracketed.","root finder did not converged.","root finder failed."};
//#define ROOTFINDER(i, x) {if ( ((i) = (x)) && ((i)>ROOT_ERRORS_NO) )  { errorexit(root_errors[(i)]); }} 
#define ROOTFINDER(i, x) {if ( ((i) = (x)) && ((i)>ROOT_ERRORS_NO) )  { printf("%s\n",root_errors[(i)]); }}

/** Maps between linear index and the corresponding (l, m) multipole indices */
extern const int LINDEX[KMAX]; /* defined in TEOBResumS.c */
extern const int MINDEX[KMAX]; /* defined in TEOBResumS.c */

/** Multipolar coefficients for NQC waveform */
typedef struct tagNQCcoefs
{
  double a1[KMAX];
  double a2[KMAX];
  double a3[KMAX];
  double b1[KMAX];
  double b2[KMAX];
  double b3[KMAX];
  double n[KMAX][6];
  int activemode[KMAX]; /* mask for modes with nonzero NQC */
  int maxk; /* index of maximum active mode */
  int add; /* flag */
} NQCcoefs;

/** NQC data for flux and waveform */
typedef struct tagNQCdata
{
  NQCcoefs *flx;
  NQCcoefs *hlm;
} NQCdata;

extern NQCdata *NQC; /* defined in TEOBResumS.c */

enum{
  NQC_NO,
  NQC_AUTO,
  NQC_MANUAL,
  NQC_NOPT
};
static const char* const nqc_opt[] = {"no", "auto", "manual"};

enum{
  NQC_FLX_NONE,
  NQC_FLX_NRFIT_NOSPIN201602,
  NQC_FLX_FIT_SPIN_202002,
  NQC_FLX_FROMFILE,
  NQC_FLX_NOPT
};
static const char* const nqc_flx_opt[] = {"none", "nrfit_nospin201602", "fit_spin_202002", "fromfile"};

enum{
  NQC_HLM_NONE,
  NQC_HLM_NRFIT_NOSPIN201602,
  NQC_HLM_FIT_SPIN_202002,
  NQC_HLM_FROMFILE,
  NQC_HLM_COMPUTE,
  NQC_HLM_NOPT
};
static const char* const nqc_hlm_opt[] = {"none", "nrfit_nospin201602", "fit_spin_202002", "fromfile", "compute"};

enum { 
  FIRSTCALL_EOBWAVFLMV1, 
  FIRSTCALL_EOBWAVFLMHM, 
  FIRSTCALL_EOBDYNSGS, 
  NFIRSTCALL
};

enum {
  DOMAIN_TD,
  DOMAIN_FD_22,
  DOMAIN_FD_HM,
};

/** Waveform data type */
typedef struct tagWaveform
{
  int size;
  double *time;
  double *real; 
  double *imag; 
  double *ampli;
  double *phase;  
  char name[STRLEN];
}  Waveform;

/** FD Waveform data type */
typedef struct tagWaveformFD
{
  int size;
  double *freq;
  double *preal; 
  double *pimag; 
  double *creal; 
  double *cimag; 
  double *ampli;
  double *phase;  
  char name[STRLEN];
}  WaveformFD;

/** Multipolar waveform data type */
typedef struct tagWaveform_lm
{
  int size;
  double *time;
  double *ampli[KMAX]; /* amplitude */
  double *phase[KMAX]; /* phase */
  char name[STRLEN];
  int kmask[KMAX]; /* mask for multipoles */
}  Waveform_lm;

/** Multipolar FD waveform data type */
typedef struct tagWaveformFD_lm
{
  int size;
  double *freq; /* uniform frequency array */
  double *ampli[KMAX]; /* amplitude */
  double *phase[KMAX]; /* phase */
  double *F[KMAX], *Fdot[KMAX]; /* Freq, and drvts for SPA */
  char name[STRLEN];
  int kmask[KMAX]; /* mask for multipoles */
}  WaveformFD_lm;


/** Multipolar waveform at given time or frequency, comes at handy */
typedef struct tagWaveform_lm_t
{
  double time;
  double freq;
  double ampli[KMAX]; /* amplitude */
  double phase[KMAX]; /* phase */
  int kmask[KMAX]; /* mask for multipoles */
}  Waveform_lm_t;

/** Dynamics data type */
typedef struct tagDynamics
{
  char name[STRLEN];
  /* various pointwise variables */
  int store; /* store following values? */
  int noflx; /* compute rhs without flux */
  double t, r, phi, pphi, prstar, ddotr, Omg, Omg_orb;
  double H, Heff, Heff_orb, E, jhat, r_omega, psi, v_phi;
  double A,dA,d2A, B,dB;
  double MOmg, MOmg_prev, tMOmgpeak;
  /* stuff for ODE solver */
  double y[EOB_EVOLVE_NVARS]; /* rhs storage */
  double dy[EOB_EVOLVE_NVARS];
  double y0[EOB_ID_NVARS]; /* ID storage */
  double dt, t_stop, ti;
  int ode_timestep;
  bool ode_stop, ode_stop_MOmgpeak, ode_stop_radius;
  /* arrays */
  int size;
  double *time;
  double *data[EOB_DYNAMICS_NVARS];
  
  /* key parameters for quick access */
  // TODO: REMOVE THEM FROM HERE, put them in EOBParameters
  double M, nu, q, X1, X2;
  double chi1, chi2, S1,S2, S,Sstar, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, a6c, cN3LO;
  double rLR, rLSO;
  double kapA2,kapA3,kapA4, kapB2,kapB3,kapB4, kapT2,kapT3,kapT4;
  double khatA2,khatB2; //FIXME: redundant, =0.5*kapB2,  should be removed and defined locally
  double bar_alph2_1, bar_alph2_2, bar_alph3_1, bar_alph3_2, bar_alph2j_1; //FIXME: these coefficients should be set at first call of metric routine (consistently with other PN coefs), and not used here
  double kapA2j, kapB2j, kapT2j;
  double rLR_tidal, pGSF_tidal;
  double Mbhf, abhf; /* final BH */
  int use_tidal, use_spins, use_tidal_gravitomagnetic;
} Dynamics;

/** Parameter data type */
typedef struct tagEOBParameters
{
  double M, nu, q, X1, X2;
  double chi1, chi2, S1,S2, S,Sstar, a1, a2, aK, aK2;
  double C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, a6c, cN3LO;
  double rLR, rLSO;
  double LambdaAl2,LambdaAl3,LambdaAl4, LambdaBl2,LambdaBl3,LambdaBl4, SigmaAl2,SigmaBl2;
  double kapA2,kapA3,kapA4, kapB2,kapB3,kapB4, kapT2,kapT3,kapT4;
  double japA2,japA3,japA4, japB2,japB3,japB4, japT2,japT3,japT4;//new names!
  
  double khatA2,khatB2; //FIXME: redundant, =0.5*kapB2,  should be removed and defined locally
  double bar_alph2_1, bar_alph2_2, bar_alph3_1, bar_alph3_2, bar_alph2j_1; //FIXME: these coefficients should be set at first call of metric routine (consistently with other PN coefs), and not used here
  double kapA2j, kapB2j, kapT2j;
  double rLR_tidal, pGSF_tidal;
  double Mbhf, abhf; // final BH 

  double r0, initial_frequency;

  double distance, inclination, polarization, coalescence_angle;
  
  int use_tidal, use_spins, use_tidal_gravitomagnetic;
  int use_Yagi_fits;
  int use_geometric_units;
  int use_speedytail;

  double dt_merger_interp, dt_interp, srate_interp;
  int interp_uniform_grid;

  int *use_mode_lm, use_mode_lm_size;

  int postadiabatic_dynamics, postadiabatic_dynamics_stop;
  int postadiabatic_dynamics_N;
  int postadiabatic_dynamics_size;
  double postadiabatic_dynamics_rmin;

  int centrifugal_radius; // NEW, INDEX FOR # {LO, NLO, NNLO, NNLOS4, NOSPIN, NOTIDES}
  int use_flm; //NEW, INDEX FOR  # "SSLO", "SSNLO", "SSNNLO", "HM"

  int compute_LR, compute_LSO, compute_LR_guess, compute_LSO_guess;

  int nqc, nqc_coefs_flx, nqc_coefs_hlm; // NEW, INDEXES
  char nqc_coefs_flx_file[STRLEN], nqc_coefs_hlm_file[STRLEN];

  char output_dir[STRLEN];
  int output_hpc, output_multipoles, output_dynamics, output_nqc, output_nqc_coefs, output_ringdown;
  int *output_lm, output_lm_size; 

  double srate, dt;
  int size;
  int ringdown_extend_array;
  int ode_timestep;
  double ode_abstol, ode_reltol;
  double ode_tmax;
  int ode_stop_afterNdt;
  int ode_stop, ode_stop_MOmgpeak, ode_stop_radius;

  int openmp_threads, openmp_timeron;

  int firstcall[NFIRSTCALL];

  int domain; //Time or frequency domain
  int interp_FD_waveform;
  double df;  //frequency interp df, can be set from srate (?)
  
} EOBParameters;

extern EOBParameters *EOBPars; /* defined in TEOBResumSPars.c */ 

/* Function protoypes grouped based on file */

/* TEOBResumS.c */
int EOBRun(Waveform **hpc, WaveformFD **hfpc, 
	   int default_choice, int firstcall);

/* TEOBResumSPars.c */
void par_db_init ();
void par_db_free ();
void par_db_default_fromfile ();
void par_db_default ();
void par_db_from_EOBPar (EOBParameters *EOBPars);
void par_file_parse (const char *fname);
void par_file_parse_merge (const char *fname);
void par_db_write_file (const char *fname);
void par_db_screen (const int pr);
void par_add_i (const char *key, int val);
void par_add_b (const char *key, int val);
void par_add_d (const char *key, double val);
void par_add_s (const char *key, const char *val);
void par_add_arrayi (const char *key, int *a, int size);
void par_add_arrayd (const char *key, double *a, int size);
void par_set_i(const char *key, int val);
void par_set_b(const char *key, int val);
void par_set_d(const char *key, double val);
void par_set_s(const char *key, const char *val);
void par_set_arrayi (const char *key, int *array, int n);
void par_set_arrayd (const char *key, double *array, int n);
int par_get_i(const char *key);
int par_get_b(const char *key);
double par_get_d(const char *key);
const char * par_get_s(const char *key);
int * par_get_arrayi(const char *key, int *n);
double * par_get_arrayd(const char *key, int *n);
void eob_set_params(char *s, int n);
void eob_set_params_new(char *s, int n, int default_choice);
void eob_set_params_EOBRun(int default_choice, int firstcall);
void eob_free_params();
void EOBParameters_alloc (EOBParameters **eobp);
void EOBParameters_free (EOBParameters *eobp);
void EOBParameters_defaults (int choose, EOBParameters *eobp);
void EOBParameters_set_from_db (EOBParameters *eobp);

/* TEOBResumSUtil.c */
double q_to_nu(const double q);
double nu_to_X1(const double nu);
double Eulerlog(const double x,const int m);
void interp_spline(double *t, double *y, int n, double *ti, int ni, double *yi);
void interp_spline_checklim(double *t, double *y, int n, double *ti, int ni, double *yi);
void interp_spline_omp(double *t, double *y, int n, double *ti, int ni, double *yi);
int find_point_bisection(double x, int n, double *xp, int o);
double baryc_f(double xx, int n, double *f, double *x);
void baryc_weights(int n, double *x, double *omega);
double baryc_f_weights(double xx, int n, double *f, double *x, double *omega);
double interp1d (const int order, double xx, int nx, double *f, double *x);
double find_max (const int n, double dx, double x0, double *f, double *fmax);
double fact(int n);
double wigner_d_function(int l, int m, int s, double i);
int spinsphericalharm(double *rY, double *iY, int s, int l, int m, double phi, double i);
int D0(double *f, double dx, int n, double *df);
int D2(double *f, double dx, int n, double *d2f);
int D0_x(double *f, double *x, int n, double *df);
int D0_x_4(double *f, double *x, int n, double *df);
double d4(double *x, double *y, int i);
double l_deriv(double *x,int i,int j);
double cumtrapz(double *f, double *x, const int n, double *sum);
double cumint3(double *f, double *x, const int n, double *sum);
void unwrap(double *p, const int size);
void unwrap_proxy(double *p, double *r, const int size, const int shift0);
void set_multipolar_idx_mask_old(int *kmask, int n);
void set_multipolar_idx_mask(int *kmask, int n, const int *idx, int m, int on);
int get_uniform_size(const double tf, const double t0, const double dt);
void Waveform_alloc (Waveform **wav, const int size, const char *name);
void Waveform_push (Waveform **wav, int size);
void Waveform_rmap (Waveform *h, const int mode, const int unw);
void Waveform_interp (Waveform *h, const int size, const double t0, const double dt, const char *name);                                                        
void Waveform_interp_ap (Waveform *h, const int size, const double t0, const double dt, const char *name);
void Waveform_output (Waveform *wav);
void Waveform_free (Waveform *wav);
void Waveform_lm_alloc (Waveform_lm **wav, int size, const char *name);
void Waveform_lm_push (Waveform_lm **wav, int size);
void Waveform_lm_output (Waveform_lm *wav);
void Waveform_lm_output_reim (Waveform_lm *wav);
void Waveform_lm_free (Waveform_lm *wav);
void Waveform_lm_interp (Waveform_lm *hlm, const int size, const double t0, const double dt, const char *name);
/* void Waveform_lm_alloc_interp (Waveform_lm *hlm, Waveform_lm **hlm_new, const int size, const double t0, const double dt, const char *name); */
void Waveform_lm_extract (Waveform_lm *hlma, const double to, const double tn, Waveform_lm **hlmb, const char *name);
void Waveform_lm_join (Waveform_lm *hlma, Waveform_lm *hlmb, double to); 
void WaveformFD_alloc (WaveformFD **wav, const int size, const char *name);
void WaveformFD_push (WaveformFD **wav, int size);
void WaveformFD_interp_ap (WaveformFD *h, const int size, const double f0, const double df, const char *name);
void WaveformFD_output (WaveformFD *wav);
void WaveformFD_free (WaveformFD *wav);
void WaveformFD_lm_alloc (WaveformFD_lm **wav, int size, const char *name);
void WaveformFD_lm_push (WaveformFD_lm **wav, int size);
void WaveformFD_lm_output (WaveformFD_lm *wav);
void WaveformFD_lm_output_reim (WaveformFD_lm *wav);
void WaveformFD_lm_free (WaveformFD_lm *wav);
void WaveformFD_lm_interp_ap (WaveformFD_lm *hlm, const int size, const double f0, const double df, const char *name);
void Waveform_lm_t_alloc (Waveform_lm_t **wav);
void Waveform_lm_t_free (Waveform_lm_t *wav);
void Dynamics_alloc (Dynamics **dyn, int size, const char *name);
void Dynamics_push (Dynamics **dyn, int size);
void Dynamics_output (Dynamics *dyn);
void Dynamics_free (Dynamics *dyn);
void Dynamics_interp (Dynamics *dyn, const int size, const double t0, const double dt, const char *name);
void Dynamics_extract (Dynamics *dyna, const double to, const double tn, Dynamics **dynb, const char *name);
void Dynamics_join (Dynamics *dyna, Dynamics *dynb, double to);
void Dynamics_set_params (Dynamics *dyn);
void NQCdata_alloc (NQCdata **nqc);
void NQCdata_free (NQCdata *nqc);
void SPA(Waveform_lm *TDlm, WaveformFD_lm *FDlm);
void compute_hpc(Waveform_lm *hlm, double nu, double M, double distance, double amplitude_prefactor, double psi, double iota, Waveform *hpc);
void compute_hpc_FD(WaveformFD_lm *hlm, double nu, double M, double distance, double amplitude_prefactor, double phi, double iota, WaveformFD *hpc);
double time_units_factor(double M);
double time_units_conversion(double M, double t);
double radius0(double M, double fHz);
int system_mkdir(const char *name);
void print_date_time();
void errorexit(char *file, int line, const char *s);
#define errorexit(s) errorexit(__FILE__, __LINE__, (s))
void errorexits(char *file, int line, const char *s, const char *t);
#define errorexits(s,t) errorexits(__FILE__, __LINE__, (s), (t))

/* TEOBResumSFits.c */
double eob_a6c_fit(double nu);
double eob_a6c_fit_HM(double nu);
double eob_c3_fit_global(double nu, double a1, double a2);
double eob_c3_fit_HM(double nu, double a1, double a2);
void eob_nqc_point(Dynamics *dyn, double *A_tmp, double *dA_tmp, double *omg_tmp, double *domg_tmp);
void eob_nqc_point_HM(Dynamics *dyn, double *A_tmp, double *dA_tmp, double *omg_tmp, double *domg_tmp);
void eob_nqc_point_postpeak(double Mbh, double c1A, double c2A, double c3A, double c4A, 
			    double c1phi, double c2phi, double c3phi, double c4phi,
			    double alpha1, double omega1,
			    double *A_tmp, double *dA_tmp, double *omg_tmp, double *domg_tmp);
double eob_nqc_dtfit(const double chi, const double chi0);
double eob_nqc_timeshift(double nu, double chi1);
void eob_nqc_deltat_lm(Dynamics *dyn, double *Dt_lm);
void eob_nqc_setcoefs(NQCdata *nqc);
void eob_nqc_setcoefs_nospin201602(NQCcoefs *nqc);
void eob_nqc_setcoefs_fromfile(NQCcoefs *nqc, const char *fname);
void eob_nqc_setcoefs_spin202002(NQCcoefs *nqc);
double get_a1_fit_22(double nu, double chi1, double chi2);
double get_a2_fit_22(double nu, double chi1, double chi2);
double logQ(double x);
double Yagi13_fit_barlamdel(double barlam2, int ell);
double Yagi13_fit_barsigmalambda(double barlam2);
double Yagi14_fit_Coct(double C_Q);
double Yagi14_fit_Chex(double C_Q);
double JFAPG_fit_Sigma_Irrotational(double barlam2);
double JFAPG_fit_Sigma_Static(double barlam2);
void HealyBBHFitRemnant(double chi1,double chi2, double q, double *mass, double *spin);
double JimenezFortezaRemnantMass(double nu, double X1, double X2, double chi1, double chi2);
double JimenezFortezaRemnantSpin(double nu, double X1, double X2, double chi1, double chi2);
void QNMHybridFitCab(double nu, double X1, double X2, double chi1, double chi2, double aK,
		     double Mbh, double abh,  
		     double *a1, double *a2, double *a3, double *a4, double *b1, double *b2, double *b3, double *b4, 
		     double *sigmar, double *sigmai);
void QNMHybridFitCab_HM(double nu, double X1, double X2, double chi1, double chi2, double aK,
			double Mbh, double abh,  
			double *a1, double *a2, double *a3, double *a4, double *b1, double *b2, double *b3, double *b4, 
			double *sigmar, double *sigmai);
void QNM_coefs(double af, double *alpha21, double *alpha1, double *omega1);
double eob_approxLR(const double nu);
double get_mrg_timestep(double q, double chi1, double chi2);

/* TEOBResumSDynamics.c */
int eob_dyn_rhs(double t, const double y[], double dy[], void *params);
void eob_ham(double nu, double r, double pph, double prstar, double A, double dA,
	     double *H, double *Heff, double *dHeff_dr, double *dHeff_dprstar, double *dHeff_dpphi);
int eob_dyn_rhs_s(double t, const double y[], double dy[], void *params);
void eob_ham_s(double nu, double r, double rc, double drc_dr, double pphi, double prstar, double S, double Sstar, double chi1, double chi2, double X1, double X2, double aK2, double c3, double A, double dA, double *H, double *Heff, double *Heff_orb, double *dHeff_dr, double *dHeff_dprstar, double *dHeff_dpphi, double *d2Heff_dprstar20);
void eob_dyn_s_GS(double r, double rc, double drc_dr, double aK2, double prstar, double pph, double nu, double chi1, double chi2, double X1, double X2, double cN3LO, double *ggm);
void (*eob_dyn_s_get_rc)();
void eob_dyn_s_get_rc_LO(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, double C_Oct1, double C_Oct2, double C_Hex1, double C_Hex2, int usetidal, double *rc, double *drc_dr, double *d2rc_dr2);
void eob_dyn_s_get_rc_NLO(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, double C_Oct1, double C_Oct2, double C_Hex1, double C_Hex2, int usetidal, double *rc, double *drc_dr, double *d2rc_dr2);
void eob_dyn_s_get_rc_NNLO(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, double C_Oct1, double C_Oct2, double C_Hex1, double C_Hex2, int usetidal, double *rc, double *drc_dr, double *d2rc_dr2);
void eob_dyn_s_get_rc_NNLO_S4(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, double C_Oct1, double C_Oct2, double C_Hex1, double C_Hex2, int usetidal, double *rc, double *drc_dr, double *d2rc_dr2);
void eob_dyn_s_get_rc_NOSPIN(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, double C_Oct1, double C_Oct2, double C_Hex1, double C_Hex2, int usetidal, double *rc, double *drc_dr, double *d2rc_dr2);
void eob_dyn_s_get_rc_NOTIDES(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, double C_Oct1, double C_Oct2, double C_Hex1, double C_Hex2, int usetidal, double *rc, double *drc_dr, double *d2rc_dr2);
double eob_dyn_fLR(double r, void * params);
int eob_dyn_adiabLR(Dynamics *dyn, double *rLR);
double eob_dyn_fLSO(double r, void * params);
int eob_dyn_adiabLSO(Dynamics *dyn, double *rLSO);

/* TEOBResumSPostAdiabatic.c */
int eob_dyn_Npostadiabatic(Dynamics *dyn, double r0);

/* TEOBResumSInitialCondition.c */
void eob_dyn_ic(double r0, Dynamics *dyn, double y_init[]);
void eob_dyn_ic_s(double r0, Dynamics *dyn, double y_init[]);
double eob_dyn_bisecHeff0_s(double nu, double chi1, double chi2, double X1, double X2, double c3, double pph, double rorb, double A, double dA, double rc, double drc_dr, double ak2, double S, double Ss);
double eob_dyn_DHeff0(double x, void *params);
double eob_dyn_r0_Kepler (double f0);
double eob_dyn_r0_eob (double f0, Dynamics *dyn);
double eob_dyn_Omegaorb0(double r, void *params);
double eob_dyn_bisecOmegaorb0(Dynamics *dyn, double omg_orb0,double r0_kepl);

/* TEOBResumSMetric.c */
void eob_metric_A5PNlog(double r, double nu, double *A, double *dA, double *d2A);
void eob_metric_Atidal(double r, Dynamics *dyn, double *AT, double *dAT, double *d2AT);
void eob_metric(double r, Dynamics *dyn, double *A, double *B, double *dA, double *d2A, double *dB);
void eob_metric_s(double r, Dynamics *dyn, double *A, double *B, double *dA, double *d2A, double *dB);

/* TEOBResumSFlux.c */
double eob_flx_Flux(double x, double Omega, double r_omega, double E, double Heff, double jhat, double r, double pr_star, double ddotr, Dynamics *dyn);
double eob_flx_Flux_s(double x, double Omega, double r_omega, double E, double Heff, double jhat, double r, double pr_star, double ddotr, Dynamics *dyn);
void eob_flx_Tlm(double w, double *MTlm);
void eob_flx_FlmNewt(double x, double nu, double *Nlm);
double eob_flx_HorizonFlux(double x, double Heff, double jhat, double nu);
double eob_flx_HorizonFlux_s(double x, double Heff, double jhat, double nu, double X1, double X2, double chi1, double chi2);

/* TEOBResumSWaveform.c */
void eob_wav_hlm(Dynamics *dyn, Waveform_lm_t *hlm);
void (*eob_wav_deltalm)();
void eob_wav_deltalm_v1(double Hreal,double Omega,double nu, double *dlm);
void eob_wav_deltalm_HM(double Hreal,double Omega,double nu, double *dlm);
void eob_wav_hhatlmTail(double Omega,double Hreal,double bphys, Waveform_lm_t *tlm);
void eob_wav_speedyTail(double Omega, double Hreal, double bphys, Waveform_lm_t *tlm);
void (*eob_wav_hlmNewt)();
void eob_wav_hlmNewt_v1(double r, double Omega, double phi, double nu, Waveform_lm_t *hNewt);
void eob_wav_hlmNewt_HM(double r, double Omega, double phi, double nu, Waveform_lm_t *hNewt);
void eob_wav_hlmTidal(double x, Dynamics *dyn, double *hTidallm);
void (*eob_wav_flm)();
void eob_wav_flm_v1(double x,double nu, double *rholm, double *flm);
void eob_wav_flm_old(double x,double nu, double *rholm, double *flm);
void eob_wav_flm_HM(double x,double nu, double *rholm, double *flm);
void (*eob_wav_flm_s)();
void eob_wav_flm_s_SSNLO(double x, double nu, double X1, double X2, double chi1, double chi2, double a1, double a2, double C_Q1, double C_Q2, int usetidal, double *rholm, double *flm);
void eob_wav_flm_s_SSLO(double x, double nu, double X1, double X2, double chi1, double chi2, double a1, double a2, double C_Q1, double C_Q2, int usetidal, double *rholm, double *flm);
void eob_wav_flm_s_old(double x, double nu, double X1, double X2, double chi1, double chi2, double a1, double a2, double C_Q1, double C_Q2, int usetidal, double *rholm, double *flm);
void eob_wav_flm_s_HM(double x, double nu, double X1, double X2, double chi1, double chi2, double a1, double a2, double C_Q1, double C_Q2, int usetidal, double *rholm, double *flm);
void eob_wav_hlmNQC_find_a1a2a3(Dynamics *dyn, Waveform_lm *h, Waveform_lm *hnqc);
void eob_wav_hlmNQC_find_a1a2a3_mrg(Dynamics *dyn_mrg, Waveform_lm *hlm_mrg, Waveform_lm *hnqc,
				    Dynamics *dyn, Waveform_lm *hlm);
void eob_wav_hlmNQC(double  nu, double  r, double  prstar, double  Omega, double  ddotr, NQCcoefs *nqc, Waveform_lm_t *hlmnqc);
void eob_wav_hlmNQC_nospin201602(double  nu, double  r, double  prstar, double  Omega, double  ddotr, Waveform_lm_t *hlmnqc);
void eob_wav_ringdown_template(double x, double a1, double a2, double a3, double a4, double b1, double b2, double b3, double b4, double sigmar, double sigmai, double *psi);
void (*eob_wav_ringdown)();
void eob_wav_ringdown_v1(Dynamics *dyn, Waveform_lm *hlm);
void eob_wav_ringdown_HM(Dynamics *dyn, Waveform_lm *hlm);

#ifdef _OPENMP
/* TEOBResumSOMP.c */
void openmp_init(const int verbose);
void openmp_timer_start(char *name);
void openmp_timer_stop(char *name);
void openmp_free();
#endif
