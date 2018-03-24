/**
 *  Copyright (C) 2017 Sebastiano Bernuzzi, Gregorio Carullo, Walter Del Pozzo, Alessandro Nagar, Ka Wa Tsang
 *  This file is part of TEOBResumS
 *
 *  TEOBResumS is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  TEOBResumS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

/**
 * @file TEOBResumS..h
 * @brief Header file of the TEOBResumS C code
 *
 * This file contains all the macros, typdef, and routine prototype.
 * Doxygen documentation should go here.
 *
 */

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

/** Following macros can be set during compilation for special feats */
#ifndef VERBOSE
#define VERBOSE 0 /* verbose mode is off by default */
#endif

#ifndef POSTPOSTCIRCULAR
#define POSTPOSTCIRCULAR 1 /* use post-post-circular initial conditions by default */
#endif

#ifndef EXCLUDESPINSPINTIDES
#define EXCLUDESPINSPINTIDES 0 /* use tidally deformed centr. radius with self-spin and tides by default */
#endif

#ifndef DEBUG 
#define DEBUG 0 /* global debug option */ 
#endif

/** Macros */
#define ERROR 1 /** generic error int */
#define OK 0 /** generic go int */
#define STRLEN 128 /** Standard string length */
#define TEOBResumS_Info "TEOBResumS code (C) 2017\n"
#define TEOBResumS_Usage(x) {printf("%sUSAGE:\t%s <parfile>\n", TEOBResumS_Info, x);} 
#define SIGN(x,y) ((y) >= 0.0 ? fabs(x) : -fabs(x)) 
#define typeof __typeof__
#define MAX(a,b)                \
  ({ typeof (a) _a = (a);       \
    typeof (b) _b = (b);        \
    _a > _b ? _a : _b; })       
#define MIN(a,b)                \
  ({ typeof (a) _a = (a);       \
    typeof (b) _b = (b);        \
    _a < _b ? _a : _b; })       
#define MAX3(a,b,c) (((a) > (b)) ? MAX(a,c) : MAX(b,c))
#define MIN3(a,b,c) (((a) < (b)) ? MIN(a,c) : MIN(b,c))
#define SQ(a) ((a)*(a))
#define DEQUAL(a,b,eps) (fabs((a)-(b))<(eps)) /** double compare */
#define DUNEQUAL(a,b,eps) (fabs((a)-(b))>(eps))
#define STREQUAL(s,t) ((strcmp((s),(t))==0)) /** string compare */  
#define SWAPTRS(a,b)   \
  ({		    \
    typeof(a) temp; \
    temp = a;	    \
    a = b;	    \
    b = temp;	    \
  })
#define PRFORM(s,x) {printf("%-40s = %.12e\n", s,x);}//TODO: handle also int
/* helpers for debug */
#define DBGPR(s) printf("DEBUG: %s\n",s);
#define DBGSTOP errorexit("DEBUG: STOP");

/* Useful constants */
#define Pi 3.1415926535897932384626433832795028
#define Sqrt2 1.41421356237309504880168872420969808
#define Sqrt3 1.73205080756887729352744634150587237
#define ooSqrt2 0.707106781186547524400844362104849039284836
#define Log2  0.6931471805599453094172321
#define MSUN_M 1.47662504e3; /** G/c^3 */
#define MSUN_S 4.925491025543575903411922162094833998e-6 /** G/c^2 */
#define MPC_M  3.086e22
#define EulerGamma 0.5772156649015328606065121
  
/** Index list of EOB evolved variables */
enum{
    EOB_EVOLVE_RAD, 
    EOB_EVOLVE_PHI,
    EOB_EVOLVE_PRSTAR,
    EOB_EVOLVE_PPHI,
    EOB_EVOLVE_NVARS
};
static const char* eob_evolve_var[] = {"r","phi","prstar","pphi"};

/** Index list of EOB variables for initial data */
enum{
  EOB_ID_RAD, 
  EOB_ID_PPHI,
  EOB_ID_PRSTAR,
  EOB_ID_PR,
  EOB_ID_J,
  EOB_ID_E0,
  EOB_ID_OMGJ,
  EOB_ID_NVARS
};
static const char* eob_id_var[] = {"r","pphi","prstar","pr","j","E0","Omega"};

/** Index list of EOB dynamical variables (to be stored in arrays) */ 
enum{
  EOB_RAD, 
  EOB_PHI,
  EOB_PPHI,
  EOB_MOMG,
  EOB_DDOTR,
  EOB_PRSTAR,
  EOB_OMGORB,
  EOB_DYNAMICS_NVARS
};
static const char* eob_var[] = {"r","phi","pphi","momega","ddor","prstar","Omgorb"};

#define KMAX (35) /** Multipolar linear index, max value */
#define PMTERMS_eps (1) /** Switch on Fujita-Iyer point-mass terms. This is hard-coded here */

/** List of options for tidal potential */
enum{
  TIDES_OFF,  /* = 0 , keep first to allow syntax: if(use_tidal) { ... */
  TIDES_NNLO, 
  TIDES_TEOBRESUM,
  TIDES_NOPT
};
static const char* const tides_opt[] = {"no","NNLO","TEOBRESUM","undefined"};

/** List of options for ODE timestepping */
enum{
  ODE_TSTEP_UNIFORM, 
  ODE_TSTEP_ADAPTIVE,
  ODE_TSTEP_ADAPTIVE_UNIFORM_AFTER_LSO,
  ODE_TSTEP_NOPT
};
static const char* const ode_tstep_opt[] = {"uniform","adaptive","adaptive+uniform_after_LSO","undefined"};

/** Error handler for root finders */
enum{ 
  ROOT_ERRORS_NO,
  ROOT_ERRORS_BRACKET,
  ROOT_ERRORS_MAXITS,
  ROOT_ERRORS_NOSUCC,
  ROOT_ERRORS
};
static const char* const root_errors[] = {"none","root is not bracketed.","root finder did not converged.","root finder failed."};
#define ROOTFINDER(i, x) {if ( ((i) = (x)) && ((i)>ROOT_ERRORS_NO) )  { errorexit(root_errors[(i)]); }} //TODO: CHECK THIS MACRO (LOGIC INVOLVED)

/** Maps between linear index and the corresponding (l, m) multipole indices */
extern const int LINDEX[KMAX]; /* defined in TEOBResumS.c */
extern const int MINDEX[KMAX]; /* defined in TEOBResumS.c */

/** Type for complex double */
typedef double complex cdouble;

/** Waveform data type */
typedef struct tagWaveform
{
  int size;
  double *time;
  double *real;
  double *imag;
  char name[STRLEN];
}  Waveform;

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

/** Multipolar waveform at given time, comes at handy */
typedef struct tagWaveform_lm_t
{
  double time;
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
  double t, r, phi, pphi, prstar, ddotr, Omg, Omg_orb;
  double H, Heff, Heff_orb, E, jhat, r_omega, psi, v_phi;
  double A,dA,d2A, B,dB;
  double MOmg, MOmg_prev;
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
  double M, nu, q, X1, X2;
  double chi1, chi2, S1,S2, S,Sstar, a1, a2, aK2, C_Q1,C_Q2, cN3LO;
  double rLR, rLSO;
  double kapA2,kapA3,kapA4, kapB2,kapB3,kapB4, kapT2,kapT3,kapT4, khatA2,khatB2;
  double bar_alph2_1, bar_alph2_2, bar_alph3_1, bar_alph3_2;
  double rLR_tidal, pGSF_tidal;
  double Mbhf, abhf; /* final BH */
  int use_tidal, use_spins;
} Dynamics;

/* Function protoypes grouped based on file */

/* TEOBResumSPars.c */
void par_db_init ();
void par_db_free ();
void par_db_default ();
void par_file_parse (const char *fname);
void par_file_parse_merge (const char *fname);
void par_db_write_file (const char *fname);
void par_db_screen ();
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
void eob_free_params();

/* TEOBResumSUtil.c */
double q_to_nu(const double q);
double nu_to_X1(const double nu);
double Eulerlog(const double x,const double m);
void interp_spline(double *t, double *y, int n, double *ti, int ni, double *yi);
int find_point_bisection(double x, int n, double *xp, int o);
double baryc_f(double xx, int n, double *f, double *x);
void baryc_weights(int n, double *x, double *omega);
double baryc_f_weights(double xx, int n, double *f, double *x, double *omega);
double interp1d (const int order, double xx, int nx, double *f, double *x);
double fact(int n);
double wigner_d_function(int l, int m, int s, double i);
int spinsphericalharm(double *rY, double *iY, int s, int l, int m, double phi, double i);
void compute_hpc(Waveform_lm **hlm, double nu, double M, double distance, double amplitude_prefactor, double psi, double iota, Waveform *hpc);
int D0(double *f, double dx, int n, double *df);
int D2(double *f, double dx, int n, double *d2f);
int D0_x(double *f, double *x, int n, double *df);
void set_multipolar_idx_mask(int *kmask, int n);
void Waveform_alloc (Waveform **wav, int size, const char *name);
void Waveform_push (Waveform **wav, int size);
void Waveform_output (Waveform *wav);
void Waveform_free (Waveform *wav);
void Waveform_lm_alloc (Waveform_lm **wav, int size, const char *name);
void Waveform_lm_push (Waveform_lm **wav, int size);
void Waveform_lm_output (Waveform_lm *wav);
void Waveform_lm_free (Waveform_lm *wav);
void Waveform_lm_t_alloc (Waveform_lm_t **wav);
void Waveform_lm_t_free (Waveform_lm_t *wav);
void Dynamics_alloc (Dynamics **dyn, int size, const char *name);
void Dynamics_push (Dynamics **dyn, int size);
void Dynamics_output (Dynamics *dyn);
void Dynamics_free (Dynamics *dyn);
void Dynamics_set_params (Dynamics *dyn);
double time_units_factor(double M);
double time_units_conversion(double M, double t);
double radius0(double M, double fHz);
void system_mkdir(const char *name);
void print_date_time();
void errorexit(char *file, int line, const char *s);
#define errorexit(s) errorexit(__FILE__, __LINE__, (s))
void errorexits(char *file, int line, const char *s, const char *t);
#define errorexits(s,t) errorexits(__FILE__, __LINE__, (s), (t))

/* TEOBResumSFits.c */
double c3_fit_global(double nu, double chi1, double chi2, double X1, double X2, double a1, double a2);
double logQ(double x);
double Yagi13_fit_barlamdel(double barlam2, int ell);
void HealyBBHFitRemnant(double chi1,double chi2, double q, double *mass, double *spin);
double JimenezFortezaRemnantSpin(double nu, double X1, double X2, double chi1, double chi2);
void QNMHybridFitCab(double nu, double X1, double X2, double chi1, double chi2, double aK,
		     double Mbh, double abh,  
		     double *a1, double *a2, double *a3, double *a4, double *b1, double *b2, double *b3, double *b4, 
		     double **sigma);
double eob_approxLR(const double nu);

/* TEOBResumSDynamics.c */
int eob_dyn_rhs(double t, const double y[], double dy[], void *params);
void eob_ham(double nu, double r, double pph, double prstar, double A, double dA,
	     double *H, double *Heff, double *dHeff_dr, double *dHeff_dprstar, double *dHeff_dpphi);
int eob_dyn_rhs_s(double t, const double y[], double dy[], void *params);
void eob_ham_s(double nu, double r, double rc, double drc_dr, double pphi, double prstar, double S, double Sstar, double chi1, double chi2, double X1, double X2, double aK2, double c3, double A, double dA, double *H, double *Heff, double *Heff_orb, double *dHeff_dr, double *dHeff_dprstar, double *dHeff_dpphi, double *d2Heff_dprstar20);
void eob_dyn_s_GS(double r, double rc, double drc_dr, double aK2, double prstar, double pph, double nu, double chi1, double chi2, double X1, double X2, double cN3LO, double *ggm);
void eob_dyn_s_get_rc(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, int usetidal, double *rc, double *drc_dr, double *d2rc_dr2); 
double eob_dyn_fLR(double r, void * params);
int eob_dyn_adiabLR(Dynamics *dyn, double *rLR);
double eob_dyn_fLSO(double r, void * params);
int eob_dyn_adiabLSO(Dynamics *dyn, double *rLSO);

/* TEOBResumSInitialCondition.c */
void eob_dyn_ic(double r0, Dynamics *dyn, double y_init[]);
void eob_dyn_ic_s(double r0, Dynamics *dyn, double y_init[]);
double eob_dyn_bisecHeff0_s(double nu, double chi1, double chi2, double X1, double X2, double c3, double pph, double rorb, double A, double dA, double rc, double drc_dr, double ak2, double S, double Ss);
double eob_dyn_DHeff0(double x, void *params);

/* TEOBResumSMetric.c */
void eob_metric_A5PNlog(double r, double nu, double *A, double *dA, double *d2A);
void eob_metric_Atidal(double r, Dynamics *dyn, double *AT, double *dAT, double *d2AT);
void eob_metric(double r, Dynamics *dyn, double *A, double *B, double *dA, double *d2A, double *dB);
void eob_metric_s(double r, Dynamics *dyn, double *A, double *B, double *dA, double *d2A, double *dB);

/* TEOBResumSFlux.c */
double eob_flx_Flux(double x, double Omega, double r_omega, double E, double Heff, double jhat, double r, double pr_star, double ddotr, Dynamics *dyn);
double eob_flx_Flux_s(double x, double Omega, double r_omega, double E, double Heff, double jhat, double r, double pr_star, double ddotr, Dynamics *dyn);
void eob_flx_Tlm(double w, double *MTlm);
void eob_flx_FlmNewt(double x, double nu, int usetidal, int usespins, double *Nlm);
double eob_flx_HorizonFlux(double x, double Heff, double jhat, double nu);
double eob_flx_HorizonFlux_s(double x, double Heff, double jhat, double nu, double X1, double X2, double chi1, double chi2);

/* TEOBResumSWaveform.c */
void eob_wav_hlm(Dynamics *dyn, Waveform_lm_t *hlm);
void eob_wav_deltalm(double Hreal,double Omega,double nu, double *dlm);
void eob_wav_hhatlmTail(double Omega,double Hreal,double bphys, Waveform_lm_t *tlm);
void eob_wav_speedyTail(double Omega, double Hreal, double bphys, Waveform_lm_t *tlm);
void eob_wav_hlmNewt(double r, double Omega, double phi, double nu, int usetidal, Waveform_lm_t *hNewt);
void eob_wav_hlmTidal(double x, Dynamics *dyn, double *hTidallm);
void eob_wav_flm(double x,double nu, double *rholm, double *flm);
void eob_wav_flm_s(double x, double nu, double X1, double X2, double chi1, double chi2, double a1, double a2, double C_Q1, double C_Q2, int usetidal, double *rholm, double *flm);
void eob_wav_hlmNQC_find_a1a2a3(const int size, Dynamics *dyn, Waveform_lm *h, Waveform_lm *hnqc);
void eob_wav_hlmNQC(double  nu, double  r, double  prstar, double  Omega, double  ddotr, Waveform_lm_t *psilmnqc);
void eob_wav_ringdown_template(double x, double a1, double a2, double a3, double a4, double b1, double b2, double b3, double b4, double *sigma, double *psi);
void eob_wav_ringdown(Dynamics *dyn, Waveform_lm *hlm);

