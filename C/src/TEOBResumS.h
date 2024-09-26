/** \file TEOBResumS.h
 *  \brief Header file of the TEOBResumS C code
 *
 * This file contains all the macros, typdef, and routine prototype.
 *
 */

#define _GNU_SOURCE /* This is just for getline() */
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <unistd.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>
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
#define TEOBResumS_Info "TEOBResumS code (c) 2017-2022"
#define TEOBResumS_Usage "COMMAND LINE USAGE:\n\
./TEOBResumS.x [options] or ./TEOBResumS.x -p [parfile_name]\n\n\
Options:\n\
  -M [total mass]\n\
  -q [mass ratio]\n\
  -Z [chi1z]\n\
  -z [chi2z]\n\
  -X [chi1x]\n\
  -x [chi2x]\n\
  -Y [chi1y]\n\
  -y [chi2y]\n\
  -L [lambda1]\n\
  -l [lambda2]\n\
  -d [distance]\n\
  -i [inclination]\n\
  -f [initial frequency]\n\
  -R [sampling rate]\n\
  -o [output:yes/no]\n\
  -g [geometric_units:yes/no]\n\n\
Note: if a parfile is specified together with other options, the latter are overwritten by the parfile.\n"

#define STRLEN 1024 /** Standard string length */
 
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
#define YESNO2INT(s)((strcmp(s,"yes")==0)?1:0)
/* helpers for debug */
#define DBGPR(s) printf("DEBUG: %s\n",s);
#define DBGSTOP errorexit("DEBUG: STOP");

/* indexes of 3-vectors */
enum{Ix, Iy, Iz, IN3};

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
#define TEOB_R0_THRESHOLD (14)  /* PA minimum tolerated radius */
#define POSTADIABATIC_NSTEP_MIN (10) /* Minimum requires PA steps, any less than this, the code switches off PA */

/** Simple/generic error handler */
#define ERROR (1) /** generic error int */
enum{OK,                       /**< No error */
     ERROR_OUTOFMEM,           /**< Out of memory */
     ERROR_FILEOPEN,           /**< Error opening file */
     ERROR_MKDIR,              /**< Error while making directory */
     ERROR_ROOTFINDER,         /**< Root finder failed */
     ERROR_ODEINT,             /**< ODE solver failed */
     ERROR_INITIAL_CONDITIONS, /**< Initial conditions failed */
     ERROR_SET_PARAMS,         /**< Error setting parameters */
     ERROR_RINGDOWN,           /**< Ringdown failed */
     NERROR                    /**< Number of errors */
};
static const char* eob_error_msg[] = {
  "ok",
  "out of memory", "error opening file", "error while making directory",
  "root finder failed.", "ODE solver failed.", "initial conditions failed.",
  "error setting parameters.", "ringdown failed."
};

/** Index list of EOB evolved variables */
enum{
  EOB_EVOLVE_RAD,          /**< index of r */
  EOB_EVOLVE_PHI,          /**< index of phi */
  EOB_EVOLVE_PRSTAR,       /**< index of Prstar */
  EOB_EVOLVE_PPHI,         /**< index of Pphi */
  EOB_EVOLVE_NVARS         /**< number of EOB evolved variables */
};
static const char* eob_evolve_var[] = {"r","phi","Prstar","Pphi"};

/** Index list of EOB variables for initial data */
enum{
  EOB_ID_RAD,          /**< index of r */
  EOB_ID_PHI,          /**< index of phi */
  EOB_ID_PPHI,         /**< index of Pphi */
  EOB_ID_PRSTAR,       /**< index of Prstar */
  EOB_ID_PR,           /**< index of Pr */
  EOB_ID_J,            /**< index of j */
  EOB_ID_E0,           /**< index of E0 */
  EOB_ID_OMGJ,         /**< index of Omega */
  EOB_ID_NVARS         /**< number of EOB variables for initial data */
};
static const char* eob_id_var[] = {"r","phi","Pphi","Prstar","Pr","j","E0","Omega"};

/** Index list of EOB dynamical variables (to be stored in arrays) */ 
enum{
  EOB_RAD,               /**< index of r */
  EOB_PHI,               /**< index of phi */
  EOB_PPHI,              /**< index of Pphi */
  EOB_MOMG,              /**< index of M Omega */
  EOB_DDOTR,             /**< index of ddotr */
  EOB_PRSTAR,            /**< index of Prstar */
  EOB_OMGORB,            /**< index of pure orbital frequency Omega orb */
  EOB_E0,                /**< index of E0 */
  EOB_DYNAMICS_NVARS     /**< number of EOB dynamical variables */
};
static const char* eob_var[] = {"r","phi","Pphi","MOmega","ddor","Prstar","MOmega_orb","E"};

#define KMAX (35) /** Multipolar linear index, max value */
#define PMTERMS_eps (1) /** Switch on Fujita-Iyer point-mass terms. This is hard-coded here */

/** List for binary type */
enum{
  BINARY_BBH,            /**< binary black hole */
  BINARY_BNS,            /**< binary neutron star */
  BINARY_BHNS,           /**< black hole-neutron star */
  BINARY_BHNS_TD,        /**< black hole-neutron star, tidal disruption */
  BINARY_NOPT            /**< number of binary types */
};

/** List for fits Lambda_ell(Lambda2) */
enum{
  Lambda234_fits_NO,                /**< no fits */
  Lambda234_fits_YAGI13,            /**< Yagi et al. 2013, ell < 5 */
  Lambda234_fits_GODZIEBA20,        /**< Godzieba et al. 2020, ell < 5 */
  Lambda2345678_fits_GODZIEBA20,    /**< Godzieba et al. 2020, ell < 9 */
  Lambda234_fits_NOPT               /**< number of fits */
};
static const char* const use_lambda234_fits_opt[] = {"no","YAGI13","GODZIEBA20_ell4", "GODZIEBA20_ell8", "undefined"};

/** List of options for orbital A potential */
enum{
  A_5PNlog,        /**< 5PNlog resummed with P15 */
  A_GSF,           /**< GSF A, Pade' resummed */
  A_5PNlogP33,     /**< 5PNlog resummed with P33 */
  A_NOPT           /**< number of options */
};
static const char* const A_opt[] = {"PN", "GSF", "5PNlogP33", "undefined"};

/** List of options for orbital D potential */
enum{
  D_3PN,          /**< 3PN, Pade' resummed */
  D_GSF,          /**< GSF D, Pade' resummed */
  D_5PNP32,       /**< 5PN, Pade' resummed with P32*/ 
  D_NOPT          /**< number of options */
};
static const char* const D_opt[] = {"PN", "GSF", "5PNP32", "undefined"};

/** List of options for orbital Q potential */
enum{
  Q_3PN,          /**< 3PN, Taylor exp */
  Q_GSF,          /**< GSF Q, Pade' resummed */
  Q_5PNloc,       /**< 5PN local, Taylor exp */
  Q_NOPT          /**< number of options */
};
static const char* const Q_opt[] = {"PN", "GSF", "5PNloc", "undefined"};
enum{
  a6c_fits_NO,                /**< no fits */
  a6c_fits_V0,                /**< V0 fits, Nagar et al 2018 */
  a6c_fits_HM,                /**< HM fits, Nagar et al 2020 */
  a6c_fits_HM_2023,           /**< HM fits, Nagar et al 2023 */
  a6c_fits_ecc,               /**< ecc fits, Nagar et al TODO add ref */
  a6c_fits_P33_HM4PN22,       /**< ecc fits, Nagar et al in prep */
  a6c_fits_NOPT               /**< number of fits */
};
static const char* const use_a6c_fits_opt[] = {"no", "v0", "HM", "HM_2023", "ecc", "HM4PN22", "undefined"};

enum{
  cN3LO_fits_NO,              /**< no fits */
  cN3LO_fits_V0,              /**< V0 fits, Nagar et al 2018 */
  cN3LO_fits_HM,              /**< HM fits, Nagar et al 2020 */
  cN3LO_fits_HM_2023_420,     /**< HM fits, Nagar et al 2023, v420 */
  cN3LO_fits_HM_2023_430,     /**< HM fits, Nagar et al 2023, v430 */
  cN3LO_fits_HM_2023_431,     /**< HM fits, Nagar et al 2023, v431 */
  cN3LO_fits_HM_2023_432,     /**< HM fits, Nagar et al 2023, v432 */
  cN3LO_fits_ecc,             /**< ecc fits, Nagar et al TODO add ref */
  cN3LO_fits_P33_HM4PN22,     /**< ecc fits, Nagar et al in prep */
  cN3LO_fits_NOPT             /**< number of fits */
};
static const char* const use_cN3LO_fits_opt[] = {"no", "v0", "HM", "HM_420", "HM_430", "HM_431", "HM_432", "ecc", "HM4PN22", "undefined"};

/** List of options for tidal potential */
enum{
  TIDES_OFF,               /**< no tides , keep first to allow syntax: if(use_tidal) { ... */
  TIDES_NNLO,              /**< NNLO PN tides */
  TIDES_TEOBRESUM,         /**< TEOBResum tides */
  TIDES_TEOBRESUM3,        /**< TEOBResum3 tides, Akcay et al 2018 */
  TIDES_TEOBRESUM_BHNS,    /**< TEOBResum BHNS tides */
  TIDES_NOPT               /**< number of tides options */
};
static const char* const tides_opt[] = {"no","NNLO","TEOBRESUM", "TEOBRESUM3", "TEOBRESUM_BHNS","undefined"};

/** List of options for the gravitomagnetic tidal potential */
enum{
  TIDES_GM_OFF,            /**< no GM tides, keep first to allow syntax: if(use_tidal_gravitomagnetic) { ... */
  TIDES_GM_PN,             /**< PN GM tides */
  TIDES_GM_GSF,            /**< GSF GM tides */
  TIDES_GM_NOPT            /**< number of GM tides options */
};
static const char* const tides_gravitomagnetic_opt[] = {"no","PN","GSF","undefined"};

/** List of options for centrifugal radius */
enum{
  CENTRAD_LO,             /**< LO centrifugal radius */
  CENTRAD_NLO,            /**< NLO centrifugal radius */
  CENTRAD_NNLO,           /**< NNLO centrifugal radius */
  CENTRAD_NNLOS4,         /**< NNLO centrifugal radius with spin^4 */
  CENTRAD_NOSPIN,         /**< No spin in centrifugal radius */
  CENTRAD_NOTIDES,        /**< No tides in centrifugal radius */
  CENTRAD_NOPT            /**< number of centrifugal radius options */
};
static const char* const centrifugal_radius_opt[] = {"LO", "NLO", "NNLO", "NNLOS4", "NOSPIN", "NOTIDES"};

/** List of options for eccentric initial frequency */
enum{
  ECCFREQ_PERIASTRON,     /**< Initial frequency specified at periastron */
  ECCFREQ_AVERAGE,        /**< Initial frequency specified at average between r+ and r- */
  ECCFREQ_APASTRON,       /**< Initial frequency specified at apastron */
  ECCFREQ_ORBAVGD,        /**< Initial orbit-averaged frequency  */
  ECCFREQ_NOPT            /**< number of eccentric initial frequency options */
};
static const char* const ecc_freq_opt[] = {"periastron", "average", "apastron", "orbitaveraged"};

/** List of options for eccentric initial conditions */
enum{
  ECCICS_0PA,      /**< Adiabatic ICs */
  ECCICS_1PA,      /**< 1PA ICs */
  ECCICS_MA,       /**< Anomaly ICs */
  ECCICS_NOPT      /**< number of eccentric initial conditions options */
};
static const char* const ecc_ics_opt[] = {"0PA", "1PA", "MA", "undefined"};

/** List of options for flm amplitudes */
enum{
  USEFLM_SSLO,            /**< SSLO amplitudes */
  USEFLM_SSNLO,           /**< SSNLO amplitudes */
  USEFLM_SSNNLO,          /**< SSNNLO amplitudes */
  USEFLM_HM,              /**< HM amplitudes */
  USEFLM_HM_4PN22,        /**< HM amplitudes with 4PN in 22 */
  USEFLM_KERR,            /**< Kerr amplitudes */
  USEFLM_NOPT             /**< number of flm amplitudes options */
};
static const char* const use_flm_opt[] = {"SSLO", "SSNLO", "SSNNLO", "HM", "HM4PN22", "Kerr"};

/** List of options for ODE timestepping */
enum{
  ODE_TSTEP_UNIFORM,                       /**< uniform timestep */
  ODE_TSTEP_ADAPTIVE,                      /**< adaptive timestep */
  ODE_TSTEP_ADAPTIVE_UNIFORM_AFTER_LSO,    /**< adaptive timestep, uniform after LSO */
  ODE_TSTEP_NOPT                           /**< number of ODE timestepping options */
};
static const char* const ode_tstep_opt[] = {"uniform","adaptive","adaptive+uniform_after_LSO","undefined"};

/** List of options for 'usespin' parameter */
enum{
MODE_SPINS_NOSPIN,                       /**< no spin, deprecated */     
  MODE_SPINS_ALIGNED,                    /**< aligned spin */
  MODE_SPINS_GENERIC,                    /**< Precessing spin */
  MODE_SPINS_NOPT,                       /**< number of 'usespin' options */
};
static const char* const mode_spin_opt[] = {"nospin","aligned","generic","undefined"};

/** List for precesing vars indexes */
enum{
  EOB_EVOLVE_SPIN_SxA,                     /**< index of Sx for body 1 */
  EOB_EVOLVE_SPIN_SyA,                     /**< index of Sy for body 1 */
  EOB_EVOLVE_SPIN_SzA,                     /**< index of Sz for body 1 */
  EOB_EVOLVE_SPIN_SxB,                     /**< index of Sx for body 2 */
  EOB_EVOLVE_SPIN_SyB,                     /**< index of Sy for body 2 */
  EOB_EVOLVE_SPIN_SzB,                     /**< index of Sz for body 2 */
  EOB_EVOLVE_SPIN_Lx,                      /**< index of Lx */
  EOB_EVOLVE_SPIN_Ly,                      /**< index of Ly */
  EOB_EVOLVE_SPIN_Lz,                      /**< index of Lz */
  EOB_EVOLVE_SPIN_alp,                     /**< index of alpha euler angle */
  EOB_EVOLVE_SPIN_bet,                     /**< index of beta euler angle */
  EOB_EVOLVE_SPIN_gam,                     /**< index of gamma euler angle */
  EOB_EVOLVE_SPIN_Momg,                    /**< index of M Omega */
  EOB_EVOLVE_SPIN_NVARS,                   /**< number of precessing variables */
};
static const char* eob_prec_var[] = {
  "SxA","SyA","SzA",
  "SxB","SyB","SzB",
  "Lx","Ly","Lz",
  "alpha","beta","gamma",
  "Momega"};

/** List of options for beyond-merger euler angles extension */
enum{
  RD_EULERANGLES_CONSTANT,            /**< constant euler angles */
  RD_EULERANGLES_QNMs,                /**< QNMs euler angles */
  RD_EULERANGLES_NOPT,                /**< number of beyond-merger euler angles extension options */
};
static const char* ringdown_eulerangles_opt[] = {"constant", "QNMs", "undefined"};

/** Omegadot for spin dynamics */
enum{
  SPIN_FLX_PN,              /**< PN Omegadot */
  SPIN_FLX_EOB,             /**< EOB Omegadot */
  SPIN_FLX_EOB_HYBRIDv1,    /**< EOB+PN Omegadot, hybrid, v1 */
  SPIN_FLX_EOB_HYBRIDv2,    /**< EOB+PN Omegadot, hybrid, v2 */
  SPIN_FLX_NOPT,            /**< number of Omegadot options */
};
static const char* spin_flx_opt[] = {"PN", "EOB", "HYBv1","HYBv2","undefined"};

/** Error handler for root finders */
enum{ 
  ROOT_ERRORS_NO,            /**< no error */
  ROOT_ERRORS_BRACKET,       /**< root is not bracketed */
  ROOT_ERRORS_MAXITS,        /**< root finder did not converge */
  ROOT_ERRORS_NOSUCC,        /**< root finder failed */
  ROOT_ERRORS                /**< number of root finder errors */
};
static const char* const root_errors[] = {"none","root is not bracketed.","root finder did not converged.","root finder failed."};
/* #define ROOTFINDER(i, x) {if ( ((i) = (x)) && ((i)>ROOT_ERRORS_NO) )  { errorexit(root_errors[(i)]); }} */
#define ROOTFINDER(i, x) {if ( ((i) = (x)) && ((i)>ROOT_ERRORS_NO) )  { printf("%s\n",root_errors[(i)]); }}

/** Maps between linear index and the corresponding (l, m) multipole indices */
extern const int LINDEX[KMAX]; /* defined in TEOBResumS.c */
extern const int MINDEX[KMAX]; /* defined in TEOBResumS.c */
extern const int KINDEX[9][9]; /* defined in TEOBResumS.c */ //FIXME: hardcoded for KMAX=35

/** Multipolar coefficients for NQC waveform */
typedef struct tagNQCcoefs
{
  double a1[KMAX];       /**< array of a1 amplitude coefficient */
  double a2[KMAX];       /**< array of a2 amplitude coefficient */
  double a3[KMAX];       /**< array of a3 amplitude coefficient */
  double b1[KMAX];       /**< array of b1 phase coefficient */
  double b2[KMAX];       /**< array of b2 phase coefficient */
  double b3[KMAX];       /**< array of b3 phase coefficient */
  double n[KMAX][6];
  int activemode[KMAX];  /**< mask for modes with nonzero NQC */
  int maxk;              /**< index of maximum active mode */
  int add;               /**< flag */
} NQCcoefs;

/** NQC data for flux and waveform */
typedef struct tagNQCdata
{
  NQCcoefs *flx;     /**< NQC coefficients for flux */
  NQCcoefs *hlm;     /**< NQC coefficients for waveform */
} NQCdata;

extern NQCdata *NQC; /* defined in TEOBResumS.c */

enum{
  NQC_NO,               /**< NQCs off */
  NQC_AUTO,             /**< NQCs automatically set */
  NQC_MANUAL,           /**< NQCs manually set */
  NQC_NOPT              /**< number of NQC options */
};
static const char* const nqc_opt[] = {"no", "auto", "manual"};

enum{
  NQC_FLX_NONE,                    /**< no NQCs */
  NQC_FLX_NRFIT_NOSPIN_201602,     /**< NQCs from NRfit, no spin  */
  NQC_FLX_NRFIT_SPIN_202002,       /**< NQCs from NRfit, spin  */
  NQC_FLX_FROMFILE,                /**< NQCs from file (iterative method)*/
  NQC_FLX_NOPT                     /**< number of NQC options */
};
static const char* const nqc_flx_opt[] = {"none", "nrfit_nospin201602", "nrfit_spin202002", "fromfile"};

enum{
  NQC_HLM_NONE,                    /**< no NQCs */
  NQC_HLM_NRFIT_NOSPIN_201602,     /**< NQCs from NRfit, no spin  */
  NQC_HLM_NRFIT_SPIN_202002,       /**< NQCs from NRfit, spin  */
  NQC_HLM_FROMFILE,                /**< NQCs from file (iterative method) */
  NQC_HLM_COMPUTE,                 /**< NQCs computed during merger attachment */
  NQC_HLM_NOPT                     /**< number of NQC options */
};
static const char* const nqc_hlm_opt[] = {"none", "nrfit_nospin201602", "nrfit_spin202002", "fromfile", "compute"};

enum { 
  FIRSTCALL_EOBWAVFLMV1,          /**< first call to eob_wav_flm_v1 */
  FIRSTCALL_EOBWAVFLMHM,          /**< first call to eob_wav_flm_HM */
  FIRSTCALL_EOBWAVFLMHM4PN22,     /**< first call to eob_wav_flm_HM_4PN22 */
  FIRSTCALL_EOBWAVFLMKERR,        /**< first call to eob_wav_flm_Kerr */
  FIRSTCALL_EOBWAVFLMKERRS,       /**< first call to eob_wav_flm_s_Kerr */ 
  FIRSTCALL_EOBDYNSGS,            /**< first call to eob_dyn_s_gs */
  FIRSTCALL_SPINDYN,              /**< first call to spin_dyn */
  FIRSTCALL_PNMOMG,               /**< first call to PN omega in spin dyn rhs */
  NFIRSTCALL                      /**< number of firstcalls */
};

enum {
  DOMAIN_TD,                   /**< time domain */
  DOMAIN_FD,                   /**< frequency domain */
};

enum {
  MODEL_DALI,                   /**< Dalì (eccentric) model */
  MODEL_GIOTTO,                 /**< Giotto (quasi-circ) model */
  MODEL_NOPT,                   /**< number of options for model */
};
static const char* const model_opt[] = {"Dali", "Giotto", "undefined"};

/** Waveform data type */
typedef struct tagWaveform
{
  int size;          /**< size of arrays */
  double *time;      /**< time array */
  double *real;      /**< real part */
  double *imag;      /**< imaginary part */
  double *ampli;     /**< amplitude */
  double *phase;     /**< phase */
  char name[STRLEN]; /**< name of waveform */
}  Waveform;

/** FD Waveform data type */
typedef struct tagWaveformFD
{
  int size;           /**< size of arrays */
  double *freq;       /**< frequency array */
  double *preal;      /**< real part of plus polarization */
  double *pimag;      /**< imaginary part of plus polarization */
  double *creal;      /**< real part of cross polarization */
  double *cimag;      /**< imaginary part of cross polarization */
  double *ampli;      /**< amplitude */
  double *phase;      /**< phase */
  char name[STRLEN];  /**< name of waveform */
}  WaveformFD;

/** Multipolar waveform data type */
typedef struct tagWaveform_lm
{
  int size;              /**< size of arrays */
  double *time;          /**< time array */
  double *ampli[KMAX];   /**< amplitude */
  double *phase[KMAX];   /**< phase */
  char name[STRLEN];     /**< name of waveform */
  int kmask[KMAX];       /**< mask for multipoles */
  int kmask_nqc[KMAX];   /**< mask for NQCs */
}  Waveform_lm;

/** Multipolar FD waveform data type */
typedef struct tagWaveformFD_lm
{
  int size;
  double *freq;                  /**< frequency array */
  double *ampli[KMAX];           /**< amplitude */
  double *phase[KMAX];           /**< phase */
  double *F[KMAX], *Fdot[KMAX];  /**< Freq, and drvts for SPA */
  char name[STRLEN];             /**< name of waveform */
  int kmask[KMAX];               /**< mask for multipoles */
  int kmask_nqc[KMAX];           /**< mask for NQCs */
}  WaveformFD_lm;


/** Multipolar waveform at given time or frequency, comes at handy */
typedef struct tagWaveform_lm_t
{
  double time;         /**< time */
  double freq;         /**< frequency */
  double ampli[KMAX];  /**< amplitude */
  double phase[KMAX];  /**< phase */
  int kmask[KMAX];     /**< mask for multipoles */
  int kmask_nqc[KMAX]; /**< mask for NQCs */
}  Waveform_lm_t;

/** Data type for spin dynamics */
typedef struct tagDynamicsSpin
{
  int size;                                                    /**< size of arrays */
  double *time;                                                /**< time array */
  double *data[EOB_EVOLVE_SPIN_NVARS];                         /**< spin dynamics data array */
  double y[EOB_EVOLVE_SPIN_NVARS], dy[EOB_EVOLVE_SPIN_NVARS];  /**< spin eqns dynamics variables and rhs */
  double t, dt;                                                /**< current time and timestep */
  double t_stop;                                               /**< stopping time, if >0 */
  double omg_stop;                                             /**< stopping frequency Momega */
  double omg_backward;                                         /**< omega which divides backward from forward spin evolution */
  double time_backward;                                        /**< time  which divides backward from forward spin evolution */
  gsl_spline *spline[EOB_EVOLVE_SPIN_NVARS];                   /**< splines for spin dynamics */
  gsl_interp_accel *accel[EOB_EVOLVE_SPIN_NVARS];              /**< accelerators for spin dynamics */
} DynamicsSpin;

/** Dynamics data type */
typedef struct tagDynamics
{
  char name[STRLEN];
  
  /* various pointwise variables */
  int store; /* store following values? */
  int noflx; /* compute rhs without flux */
  double t, r, phi, pphi, prstar, ddotr, Omg, Omg_orb;
  double rdot, r2dot, r3dot, r4dot, r5dot, Omegadot, Omega2dot, Omega3dot, Omega4dot;
  double tOmg_pk;
  double H, Heff, Heff_orb, E, jhat, r_omega, psi, v_phi;     /**< current Hamiltonian, angular momentum and derived variables */
  double A,dA,d2A, B,dB;                                      /**< current A,B and derived variables */
  double MOmg, MOmg_prev, tMOmgpeak;                          /**< current M Omega and previous value, time of the Omega peak */
  double dress_tides_fmode_A[6], dress_tides_fmode_B[6];      /**< dressing factors for adiabatic tidal pars  */
  double dress_tides_fmode_A_u[6], dress_tides_fmode_B_u[6];  /**< dressing factors for adiabatic tidal pars, drvts wrt u   */
  double dressed_C_Q1, dressed_C_Q1_u, dressed_C_Q1_uu;       /**< dressing factors for spin-induced quadrupole and drvts   */
  double dressed_C_Q2, dressed_C_Q2_u, dressed_C_Q2_uu;       /**< dressing factors for spin-induced quadrupole and drvts   */
  double dressed_C_Oct1, dressed_C_Oct1_u, dressed_C_Oct1_uu; /**< dressing factors for spin-induced octupole and drvts   */
  double dressed_C_Oct2, dressed_C_Oct2_u, dressed_C_Oct2_uu; /**< dressing factors for spin-induced octupole and drvts   */
  double dressed_C_Hex1, dressed_C_Hex1_u, dressed_C_Hex1_uu; /**< dressing factors for spin-induced hexadecapole and drvts   */
  double dressed_C_Hex2, dressed_C_Hex2_u, dressed_C_Hex2_uu; /**< dressing factors for spin-induced hexadecapole and drvts   */
  
  /* stuff for ODE solver */
  double y[EOB_EVOLVE_NVARS];   /**< rhs storage */
  double dy[EOB_EVOLVE_NVARS];  /**< rhs storage */
  double y0[EOB_ID_NVARS];      /**< ID storage */
  double dt, t_stop, ti;        /**< timestep, stopping time, current time */
  int ode_timestep;             /**< ODE timestep type */
  bool ode_stop, ode_stop_MOmgpeak, ode_stop_radius;  /**< flag to stop ode, stop after Omgpeak, stop after radius */

  /* arrays */
  int size;
  double *time;                        /**< time array */
  double *data[EOB_DYNAMICS_NVARS];    /**< dynamics data array */

  DynamicsSpin *spins;  /**< Ptr to reference spin dynamics, do not allocate mem to this one */

  
} Dynamics;

/** Parameter data type */
typedef struct tagEOBParameters
{

/* binary parameters */
  
  /**@{*/
  /** Mass parameters */
  double M, nu, q, X1, X2;                            
  /**@}*/ 

  /**@{*/
  /** Spin parameters */
  double chi1, chi2, S1,S2, S,Sstar, a1, a2, aK, aK2;
  double chi1x,chi1y,chi1z;     /* Spin components (body 1) */                   
  double chi2x,chi2y,chi2z;     /* Spin components (body 2) */       
  /**@}*/ 

  /**@{*/
  /** Tidal parameters */
  double C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2;                            /* Spin-induced multipoles*/
  double LambdaAl2,LambdaAl3,LambdaAl4,LambdaAl5,LambdaAl6,LambdaAl7,LambdaAl8; /* Multipolar electric adiabatic tidal parameters (body 1)*/
  double LambdaBl2,LambdaBl3,LambdaBl4,LambdaBl5,LambdaBl6,LambdaBl7,LambdaBl8; /* Multipolar electric adiabatic tidal parameters (body 2)*/
  double SigmaAl2,SigmaBl2;                                                     /* ell = 2 magnetic adiabatic tidal parameters*/
  double kapA2,kapA3,kapA4,kapA5,kapA6,kapA7,kapA8;                             /* EOB Multipolar electric adiabatic tidal parameters (body 1)*/
  double kapB2,kapB3,kapB4,kapB5,kapB6,kapB7,kapB8;                             /* EOB Multipolar electric adiabatic tidal parameters (body 2)*/
  double kapT2,kapT3,kapT4,kapT5,kapT6,kapT7,kapT8;                             /* EOB Multipolar electric effective adiabatic tidal parameters */
  double japA2,japA3,japA4, japB2,japB3,japB4, japT2,japT3,japT4;               /* EOB Multipolar magnetic adiabatic tidal parameters */
  double khatA2,khatB2;                                                         //FIXME: redundant, =0.5*kapB2,  should be removed and defined locally
  double bar_alph2_1, bar_alph2_2, bar_alph3_1, bar_alph3_2, bar_alph2j_1;      /* Tidal parameters combinations entering the EOB A potential*/
  double rLR_tidal, pGSF_tidal;                                                 /* Tidal light ring and GSF exponent */
  double bomgfA[6], bomgfB[6];                                                  /* f-mode frequencies star A,B (ell=2,3,4; indexes 0,1 not used) */
  /**@}*/

  /**@{*/
  /** NR-informed conservative variables */
  double a6c, cN3LO;         
  /**@}*/ 

  double r0;                    /**< Initial radial separation */
  double initial_frequency;     /**< initial frequency */
  double f0;                    /**< initial frequency in geom units */
  double distance;              /**< distance */
  double inclination;           /**< inclination */
  double polarization;          /**< polarization */
  double coalescence_angle;     /**< coalescence angle */
  double rLR;                   /**< Location of the light ring */ 
  double rLSO;                  /**< Location of the LSO */
  int compute_LR;               /**< Flag to compute the LR */
  int compute_LSO;              /**< Flag to compute the LSO */
  double compute_LR_guess;      /**< Guess for the LR */
  double compute_LSO_guess;     /**< Guess for the LSO */
  double ecc, anomaly, r_hyp, H_hyp, j_hyp;
  int ecc_freq, ecc_ics;
  double alpha_sigmoid_NQC, delta_t0_sigmoid_NQC;
  double alpha_sigmoid_Newt, delta_t0_sigmoid_Newt;
  int A_pot, D_pot, Q_pot;
  int compute_ringdown;
  
  /* options/settings */
  int binary;                                           /**< binary type (BBH, BNS, BHNS) */
  int model;                                            /**< model (Dalì, Giotto) */
  int centrifugal_radius;                               /**< NEW, INDEX FOR # {LO, NLO, NNLO, NNLOS4, NOSPIN, NOTIDES} */
  int use_flm;                                          /**< NEW, INDEX FOR  # "SSLO", "SSNLO", "SSNNLO", "HM" */
  int use_tidal, use_spins, use_tidal_gravitomagnetic;  /**< Flag for tides, spins and gravito-magnetic tides */
  int use_geometric_units;                              /**< Flag for geometric vs SI units */
  int use_speedytail;                                   /**< Flag for fast computation of tail (speedytail) */
  int use_lambda234_fits;                               /**< Flag for fits Lambda_ell(Lambda2) */
  int use_tidal_fmode_model;                            /**< Flag for f-mode tidal model */
  int use_a6c_fits;                                     /**< Flag for fits a6c */
  int use_cN3LO_fits;                                   /**< Flag for fits cN3LO */
  
  double dt_merger_interp;                              /**< timestep for merger interpolation */ 
  double dt_interp;                                     /**< timestep for interpolation */
  double srate_interp;                                  /**< sampling rate for interpolation */
  int interp_uniform_grid;                              /**< Flag for uniform grid interpolation (yes/no) */

  int *use_mode_lm, use_mode_lm_size;                   /**< Co-precessing modes array and size of co-precessing modes array*/
  int *use_mode_lm_inertial, use_mode_lm_inertial_size; /**< Inertial modes assay and size of Inertial modes assay*/
  int *use_mode_lm_nqc, use_mode_lm_nqc_size;           /**< multipoles to attach NQCs */
  int *kpostpeak, kpostpeak_size;                       /**< multipoles where NQCs are extracted at t_peak_lm + 2 from the RD template */
  int *knqcpeak22, knqcpeak22_size;                     /**< multipoles where NQCs are extracted at t_peak_22  from global fits */

  int size;                                             /**< size of arrays */
  int ode_timestep;                                     /**< ODE timestep type */
  double srate;                                         /**< sampling rate */
  double dt;                                            /**< timestep */
  double ode_abstol;                                    /**< ODE solver absolute tolerance */
  double ode_reltol;                                    /**< ODE solver relative tolerance */
  double ode_tmax;                                      /**< ODE solver max time */
  int ode_stop_afterNdt;                                /**< Number of steps to take after peal of omega */
  int ode_stop_after_peak;                              /**< Flag to check if solver is after the peak (yes/no) */
  int ode_stop;                                         /**< Flag to stop ODE solver (yes/no) */
  int ode_stop_MOmgpeak;                                /**< Flag to stop ODE solver after MOmg peak (yes/no) */
  double ode_stop_radius;                               /**< Radius to stop ODE solver,note: this is *different* in meaning from dyn->ode_stop_radius! */

  /* post-adiabatic */
  int postadiabatic_dynamics;                           /**< Flag to switch on/off post-adiabatic dynamics */
  int postadiabatic_dynamics_stop;                      /**< Flag to stop post-adiabatic dynamics */
  int postadiabatic_dynamics_N;                         /**< Number of PA iterations to take (default=8)*/
  int postadiabatic_dynamics_size;                      /**< Size of PA arrays */
  double postadiabatic_dynamics_dr;                     /**< Radial step for PA */
  double postadiabatic_dynamics_rmin;                   /**< Minimum radius to reach for PA (default=14)*/

  /* spin dynamics */
  int spin_dyn_size;                                    /**< Size of spin dynamics arrays */
  double spin_odes_omg_stop, spin_odes_t_stop;          /**< Stopping frequency and time for spin dynamics */
  double spin_odes_dt;                                  /**< Timestep for spin dynamics */
  int spin_interp_domain;                               /**< Domain for spin dynamics interpolation (time/frequency) */
  int project_spins;                                    /**< Flag to project spins during EOB dynamics (yes/no) */
  int ringdown_eulerangles;                             /**< Flag to choose which beyond-merger euler angles extension to employ */
  int spin_flx;                                         /**< Flag to choose which Omegadot for spin dynamics */

  /* final state */
  double Mbhf; /**< Final BH mass */
  double abhf; /**< Final BH spin */

  /* waveform */
  int ringdown_extend_array; /**< Extend ringdown array */

  /**@{*/
  /** NQC variables and options*/
  int nqc, nqc_coefs_flx, nqc_coefs_hlm; // NEW, INDEXES
  char nqc_coefs_flx_file[STRLEN], nqc_coefs_hlm_file[STRLEN];
  /**@}*/

  /**@{*/
  /** output options*/
  char output_dir[STRLEN];
  int output_hpc, output_multipoles, output_dynamics, output_nqc, output_nqc_coefs, output_ringdown;
  int *output_lm, output_lm_size; 
  /**@}*/

  /* msc */
  int lal_tetrad_conventions;/**< Use the LALSimulation Tetrad convention */
  int openmp_threads;        /**< Number of OpenMP threads */
  int openmp_timeron;        /**< Flag to switch on/off OpenMP timers */
  int firstcall[NFIRSTCALL]; /**< Flag to check if first call to functions */
  int errors_to_warnings;    /**< Flag to change some errors to warnings */

  int domain;                /**< Time or frequency domain */
  double tc;                 /**< Coalescence time */
  int time_shift_FD;         /**< Time shift FD waveform to have merger at t=0 */
  int time_shift_TD;         /**< Time shift TD waveform to have merger at t=0 */
  double df;                 /**< Frequency step */
  int interp_freqs;          /**< Flag to interpolate on a user-given frequency array */
  double *freqs;             /**< Frequency array */
  int freqs_size;            /**< Size of frequency array */
} EOBParameters;

extern EOBParameters *EOBPars; /* defined in TEOBResumSPars.c */ 

/* Function protoypes grouped based on file */

/* TEOBResumS.c */
int EOBRun(Waveform **hpc, WaveformFD **hfpc, 
	   Waveform_lm **hmodes, WaveformFD_lm **hfmodes, Dynamics **dynf,
	   Waveform_lm **hT_modes, Waveform_lm **hTneg_modes, Waveform_lm **hT0_modes,
     WaveformFD_lm **hfTmodes,
	   int default_choice, int firstcall);

/* TEOBResumSPars.c */
int eob_set_params(int default_choice, int firstcall);
void EOBParameters_alloc (EOBParameters **eobp);
void EOBParameters_free (EOBParameters *eobp);
void EOBParameters_defaults (int binary, int orbit, EOBParameters *eobp);
void EOBParameters_parse_file(char *fname, EOBParameters *eobp);
int EOBParameters_parse_commandline(EOBParameters *eobp, int argc, char **argv);
void EOBParameters_set_key_val(EOBParameters *eobp, char *key, char *val);
void EOBParameters_tofile(EOBParameters *eobp,char *fname);

/* TEOBResumSUtils.c */
double q_to_nu(const double q);
double nu_to_X1(const double nu);
double tidal_kappa_of_Lambda(double q, double XA, double XB, double LamA, double LamB, int ell, double *kapA, double *kapB);
void set_spin_vars(double X1, double X2, double chi1, double chi2, double *S1, double *S2, double *a1, double *a2, double *aK, double *aK2, double *S, double *Sstar);
double Eulerlog(const double x,const int m);
double Pade02(double x, double *a);
double Pade32(double x, double *a);
double Pade22(double x, double *a);
double Pade23(double x, double *a);
double Pade51(double x, double *a);
double Pade42(double x, double *a);
double Pade33(double x, double *a);
double Pade15(double x, double *a);
double Pade62(double x, double *a);
void Pade33_forGSF(double coeffs[6], double Dcoeffs[6], double D2coeffs[6], double u, double *P, double *DP, double *D2P);
void Pade76v1_forGSF(double coeffs[4], double Dcoeffs[2], double D2coeffs[2], double u, double *P, double *DP, double *D2P);
void Pade76v2_forGSF(double coeffs[5], double Dcoeffs[3], double D2coeffs[3], double u, double *P, double *DP, double *D2P);
double Taylorseries(double x, double *a, int N);
void vect_dot(double ax, double ay, double az, double bx, double by, double bz, double *s);
void vect_dot3(double *a, double *b, double *s);
void vect_cross(double ax, double ay, double az, double bx, double by, double bz, double *cx, double *cy, double *cz);
void vect_cross3(double *a, double *b, double *c);
void vect_rotate(int axis, double angle, double *vx_p, double *vy_p, double *vz_p);
void vect_rotate3(int axis, double angle, double *v);
double interp_spline_pt(double *t, double *y, int n, double ti);
void interp_spline(double *t, double *y, int n, double *ti, int ni, double *yi);
void interp_spline_checklim(double *t, double *y, int n, double *ti, int ni, double *yi);
void interp_spline_omp(double *t, double *y, int n, double *ti, int ni, double *yi);
int find_point_bisection(double x, int n, double *xp, int o);
double baryc_f(double xx, int n, double *f, double *x);
void baryc_weights(int n, double *x, double *omega);
double baryc_f_weights(double xx, int n, double *f, double *x, double *omega);
double interp1d (const int order, double xx, int nx, double *f, double *x);
double find_max (const int n, double dx, double x0, double *f, double *fmax);
double find_max_grid (double *x, double *f);
int D0(double *f, double dx, int n, double *df);
int D2(double *f, double dx, int n, double *d2f);
int D0_x(double *f, double *x, int n, double *df);
int D0_x_4(double *f, double *x, int n, double *df);
double d4(double *x, double *y, int i);
double l_deriv(double *x,int i,int j);
double cumtrapz(double *f, double *x, const int n, double *sum);
double cumint3(double *f, double *x, const int n, double *sum);
void unwrap(double *p, const int size);
void unwrap_HM(double *p, const int size);
void unwrap_euler(double *p, const int size);
void unwrap_proxy(double *p, double *r, const int size, const int shift0);
void rmap (double *re, double *im, double *p, double *a, const int mode);
void rmap_twist (double *re, double *im, double *p, double *a, const int mode);
void set_multipolar_idx_mask_old(int *kmask, int n);
void set_multipolar_idx_mask(int *kmask, int n, const int *idx, int m, int on);
long int get_uniform_size(const double tf, const double t0, const double dt);
int intersect_int (int *a, int size_a, int *b, int size_b, int *result);
void lal_tetrad_conventions(Waveform *hpc);
void Waveform_alloc (Waveform **wav, const int size, const char *name);
void Waveform_push (Waveform **wav, int size);
void Waveform_rmap (Waveform *h, const int mode, const int unw);
void Waveform_interp (Waveform *h, const int size, const double t0, const double dt, const char *name);                                                        
void Waveform_interp_ap (Waveform *h, const int size, const double t0, const double dt, const char *name);
void Waveform_output (Waveform *wav);
void Waveform_free (Waveform *wav);
void Waveform_lm_alloc (Waveform_lm **wav, int size, const char *name, int *use_mode_lm, int use_mode_lm_size);
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
void WaveformFD_lm_interp_ap_freqs (WaveformFD_lm *hlm, const char *name);
void Waveform_lm_t_alloc (Waveform_lm_t **wav);
void Waveform_lm_t_free (Waveform_lm_t *wav);
void Dynamics_alloc (Dynamics **dyn, int size, const char *name);
void Dynamics_push (Dynamics **dyn, int size);
void Dynamics_output (Dynamics *dyn);
void Dynamics_free (Dynamics *dyn);
void Dynamics_interp (Dynamics *dyn, const int size, const double t0, const double dt, const char *name);
void Dynamics_extract (Dynamics *dyna, const double to, const double tn, Dynamics **dynb, const char *name);
void Dynamics_join (Dynamics *dyna, Dynamics *dynb, double to);
void DynamicsSpin_alloc (DynamicsSpin **dyn, int size);
void DynamicsSpin_push (DynamicsSpin **dyn, int size);
void DynamicsSpin_free (DynamicsSpin *dyn);
void DynamicsSpin_output (DynamicsSpin *dyn);
void DynamicsSpin_join (DynamicsSpin *dyn, DynamicsSpin *dynb, double to);
void Dynamics_set_params (Dynamics *dyn);
void NQCdata_alloc (NQCdata **nqc);
void NQCdata_free (NQCdata *nqc);
double time_units_factor(double M);
double time_units_conversion(double M, double t);
double radius0(double M, double fHz);
int system_mkdir(const char *name);
int is_blank(const char *line);
void remove_white_spaces(char *str);
void remove_comments(char *line, const char *delimiters);
char *trim(char *str);
int getkv(char *line, char **key, char **val);
int is_string(const char *str);
char *string_trim(char *str);
int par_get_i (char *val);
int par_get_b (char *val);
double par_get_d (char *val);
int par_get_s (char * dest, char *src);
int noentries(const char *string);
int str2iarray(const char *string, int **a);
int str2darray(const char *string, double **a);
void print_date_time();
void errorexit(char *file, int line, const char *s);
#define errorexit(s) errorexit(__FILE__, __LINE__, (s))
void errorexits(char *file, int line, const char *s, const char *t);
#define errorexits(s,t) errorexits(__FILE__, __LINE__, (s), (t))

/* TEOBResumSFits.c */
double eob_a6c_fit(double nu);
double eob_a6c_fit_HM(double nu);
double eob_a6c_fit_HM_2023(double nu);
double eob_a6c_fit_ecc(double nu);
double eob_a6c_fit_next(double nu);
double eob_a6c_fit_ecc_P33_4PNh22(double nu);
double eob_c3_fit_global(double nu, double a1, double a2);
double eob_c3_fit_HM(double nu, double a1, double a2);
double eob_c3_fit_HM_2023(double nu, double a1, double a2);
double eob_c3_fit_ecc(double nu, double a1, double a2);
double eob_c3_fit_ecc_P33_4PNh22(double nu, double a1, double a2);
double eob_mrg_momg(double nu, double X1, double X2, double chi1, double chi2);
void eob_nqc_point(Dynamics *dyn, double *A_tmp, double *dA_tmp, double *omg_tmp, double *domg_tmp);
void eob_nqc_point_HM(Dynamics *dyn, double *A_tmp, double *dA_tmp, double *omg_tmp, double *domg_tmp);
void eob_nqc_point_HM_peak22(Dynamics *dyn, double *A_tmp, double *dA_tmp, double *omg_tmp, double *domg_tmp);
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
void YagiYunes13_fit_logQ_coefs(double *c);
double YagiYunes13_fit_logQ(double x);
void YagiYunes13_fit_logQ_drvts(double Lam, double Lam_u, double Lam_uu, double *Q, double *Q_u, double *Q_uu);
double Yagi13_fit_barlamdel(double barlam2, int ell);
double Yagi13_fit_barsigmalambda(double barlam2);
double Yagi14_fit_Coct(double C_Q);
double Yagi14_fit_Chex(double C_Q);
double JFAPG_fit_Sigma_Irrotational(double barlam2);
double JFAPG_fit_Sigma_Static(double barlam2);
double Godzieba20_fit_barlamdel(double barlam2, int ell);
double Chang14_fit_omegaf(double lam, int ell);
void HealyBBHFitRemnant(double chi1,double chi2, double q, double *mass, double *spin);
double JimenezFortezaRemnantMass(double nu, double X1, double X2, double chi1, double chi2);
double JimenezFortezaRemnantSpin(double nu, double X1, double X2, double chi1, double chi2);
double PrecessingRemnantSpin(Dynamics *dyn);
void QNMHybridFitCab(double nu, double X1, double X2, double chi1, double chi2, double aK,
		     double Mbh, double abh,  
		     double *ca1, double *ca2, double *ca3, double *ca4, double *cb1, double *cb2, double *cb3, double *cb4, 
		     double *sigmar, double *sigmai);
void QNMHybridFitCab_HM(double nu, double X1, double X2, double chi1, double chi2, double aK,
			double Mbh, double abh,  
			double *ca1, double *ca2, double *ca3, double *ca4, double *cb1, double *cb2, double *cb3, double *cb4, 
			double *sigmar, double *sigmai);
void QNMHybridFitCab_HM_Pompili23(double nu, double X1, double X2, double chi1, double chi2, double aK, double Mbh, double abh, 
			double *ca1, double *ca2, double *ca3, double *ca4, double *cb1, double *cb2, double *cb3, double *cb4,
      double *sigmar, double *sigmai);
void QNM_coefs(double af, double *alpha21, double *alpha1, double *omega1);
double eob_approxLR(const double nu);
double get_mrg_timestep(double q, double chi1, double chi2);
double get_mrg_timestop(double q, double chi1, double chi2);
double fmode_resonance_dress_Love(double nu, double r, double bomgf, int ell, double *dtides, double *dtides_u);
void fmode_resonance_dressing_factors(double r, Dynamics *dyn);
void fmode_resonance_dress_QOH(Dynamics *dyn);
void eob_bhns_fit(double a, double nu, double *mass, double *spin, double lambda, double m_bh, double a_bh);
void bhns_cases(double nu, double af, double chi1, double lambda, int *binary, int *usetidal);
void QNMHybridFitCab_BHNS_HM(double nu, double X1, double X2, double chi1, double chi2, double aK, 
			double Mbh, double abh,  
			double *ca1, double *ca2, double *ca3, double *ca4, double *cb1, double *cb2, double *cb3, double *cb4, 
			double *sigmar, double *sigmai);
void postpeak_coef(double *ca1, double *ca2, double *ca3, double *ca4, double *cb1, double *cb2, double *cb3, double *cb4, 
                    double *sigmar, double *sigmai, double nu, double chi1, double chi2, double X1, double X2, double aK, double Mbh, double abh,
                    double *Apeak, double *alpha2);
void peak_bhns(double nu, double kt2, double chi1, double X1, double X2, double abh, double *Apeak, double *Opeak);
void QNM_bhns_td(double af, double *alpha1, double *alpha2, double *omega1, double *omega2, double *alpha21, double lambda, double nu, double chi1);
void eob_nqc_point_BHNS_HM(Dynamics *dyn, double *A_tmp, double *dA_tmp, double *omg_tmp, double *domg_tmp, double abh, double kt2);
void eob_nqc_point_test(double Mbh, double c1A, double c2A, double c3A, double c4A, 
			    double c1phi, double c2phi, double c3phi, double c4phi,
			    double alpha1, double omega1,
			    double *A_tmp, double *dA_tmp, double *d2A_tmp, double *omg_tmp, double *domg_tmp, double *d2omg_tmp);
double horizon_radius(const double nu);
double pph_lso_orbital(const double nu);
double pph_lso_spin(const double nu, const double a0);

/* TEOBResumSDynamics.c */
extern int (*p_eob_dyn_rhs)(); /* defined in TEOBResumSPars.c */
int eob_dyn_rhs(double t, const double y[], double dy[], void *params);
void eob_ham(double nu, double r, double pphi, double prstar, double A, double dA, double Q, double dQ, double dQ_dprstar,
	           double *H, double *Heff, double *dHeff_dr, double *dHeff_dprstar, double *dHeff_dpphi);
int eob_dyn_rhs_s(double t, const double y[], double dy[], void *params);
int eob_dyn_rhs_ecc(double t, const double y[], double dy[], void *params);
void eob_ham_s(double nu, double r, double rc, double drc_dr, double d2rc_dr2, double pphi, double prstar, double S, double Sstar, double chi1, double chi2, double X1, double X2, double aK2, double c3, double A, double dA, double d2A, double Q, double dQ, double dQ_dprstar, double d2Q, double d2Q_dprstar2,
               double *H, double *Heff, double *Heff_orb, double *dHeff_dr, double *dHeff_dprstar, double *dHeff_dpphi, double *d2Heff_dprstar20, double *d2Heff_dr2);
void eob_dyn_s_GS(double r, double rc, double drc_dr, double d2rc_rd2, double aK2, double prstar, double pph, double nu, double chi1, double chi2, double X1, double X2, double cN3LO, double *ggm);

extern void (*eob_dyn_s_get_rc)(); /* defined in TEOBResumSPars.c*/
void eob_dyn_s_get_rc_LO(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, double C_Oct1, double C_Oct2, double C_Hex1, double C_Hex2, int usetidal, double *rc, double *drc_dr, double *d2rc_dr2);
void eob_dyn_s_get_rc_NLO(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, double C_Oct1, double C_Oct2, double C_Hex1, double C_Hex2, int usetidal, double *rc, double *drc_dr, double *d2rc_dr2);
void eob_dyn_s_get_rc_NNLO(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, double C_Oct1, double C_Oct2, double C_Hex1, double C_Hex2, int usetidal, double *rc, double *drc_dr, double *d2rc_dr2);
void eob_dyn_s_get_rc_NNLO_S4(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, double C_Oct1, double C_Oct2, double C_Hex1, double C_Hex2, int usetidal, double *rc, double *drc_dr, double *d2rc_dr2);
void eob_dyn_s_get_rc_NOSPIN(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, double C_Oct1, double C_Oct2, double C_Hex1, double C_Hex2, int usetidal, double *rc, double *drc_dr, double *d2rc_dr2);
void eob_dyn_s_get_rc_NOTIDES(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, double C_Oct1, double C_Oct2, double C_Hex1, double C_Hex2, int usetidal, double *rc, double *drc_dr, double *d2rc_dr2);
void eob_dyn_s_rc_add_QOH_drvts(Dynamics *dyn, double rc, double u, double at1, double at2, double *drc_dr, double *d2rc_dr2);
double eob_dyn_get_romg(double r, double prstar, double pphi, Dynamics *dyn);
double eob_dyn_fLR(double r, void * params);
int eob_dyn_adiabLR(Dynamics *dyn, double *rLR);
double eob_dyn_fLSO(double r, void * params);
int eob_dyn_fLSO_s(const gsl_vector *x, void * params, gsl_vector *f);
int eob_dyn_adiabLSO(Dynamics *dyn, double *rLSO);
int eob_dyn_LSO_s(Dynamics *dyn, double *rLSO, double *pphiLSO);
double eob_spin_dyn_alpha(double Lhx, double Lhy, double Lhz);
double alpha_initial_condition(EOBParameters *eobp); 
double eob_spin_dyn_beta(double Lhx, double Lhy, double Lhz);
extern int (*p_eob_spin_dyn_rhs)(); /* defined in TEOBResumSPars.c */
int eob_spin_dyn_rhs_PN(double t, const double y[], double dy[], void *d);
int eob_spin_dyn_rhs_PN_abc(double t, const double y[], double dy[], void *d);
int eob_spin_dyn_rhs_EOB(double t, const double y[], double dy[], void *d);
int eob_spin_dyn_integrate(DynamicsSpin *dyn, Dynamics *eobdyn, Waveform_lm *hlm);
int eob_spin_dyn_integrate_backwards(DynamicsSpin *dyn, Dynamics *eobdyn, Waveform_lm *hlm, double omega_stop);
int eob_spin_dyn(DynamicsSpin *dyn, Dynamics *eobdyn, Waveform_lm *hlm, double omg0);
void eob_spin_dyn_abc_interp(DynamicsSpin *dyn, Dynamics *eobdyn, Waveform_lm *hlm, double time, double *alpha_p, double *beta_p, double *gamma_p, int continue_integration);
void eob_spin_dyn_Sproj_interp(DynamicsSpin *dyn, double time,
			       double *SApara, double *SBpara, double *Spara,
			       double *SAperp, double *SBperp, double *Sperp,
			       int continue_integration);

/* TEOBResumSPostAdiabatic.c */
int eob_dyn_Npostadiabatic(Dynamics *dyn, double r0, DynamicsSpin *spin);

/* TEOBResumSInitialCondition.c */
extern int (*eob_dyn_ic)(); /* defined in TEOBResumSPars.c*/
int eob_dyn_ic_circ(double r0, Dynamics *dyn, double y_init[]);
int eob_dyn_ic_circ_s(double r0, Dynamics *dyn, double y_init[]);
int eob_dyn_ic_ecc(double r0, Dynamics *dyn, double y_init[]);
int eob_dyn_ic_ecc_ma(double r0, Dynamics *dyn, double y_init[]);
int eob_dyn_ic_ecc_ma_split(double r0_kepl, Dynamics *dyn, double y_init[]);
int eob_dyn_ic_ecc_PA(double r0, Dynamics *dyn, double y_init[]);
int eob_dyn_ic_hyp(double r0, Dynamics *dyn, double y_init[]);
double eob_dyn_j0(double r0, Dynamics *dyn);
double eob_dyn_ecc_j0(double r0, Dynamics *dyn);
double eob_dyn_circ_j0(double r0, Dynamics *dyn);
double eob_dyn_bisecHeff0_s(double nu, double chi1, double chi2, double X1, double X2, double c3, double pph, double rorb, double A, double dA, double rc, double drc_dr, double ak2, double S, double Ss);
double eob_dyn_DHeff0(double x, void *params);
double eob_dyn_omg_from_omgbar(double omgbar, double zeta, double ecc);
double eob_dyn_r0_Kepler (double f0);
extern double (*eob_dyn_r0_eob)(); /* defined in TEOBResumSPars.c*/
double eob_dyn_r0_circ (double f0, Dynamics *dyn);
double eob_dyn_r0_ecc (double f0, Dynamics *dyn);
double eob_dyn_Omegaorb0(double r, void *params);
double eob_dyn_Omegaecc0(double r, void *params);
double eob_dyn_Ham0(double pr, void *params);
double eob_dyn_bisecOmegaorb0(Dynamics *dyn, double omg_orb0,double r0_kepl);
double eob_dyn_bisecOmegaecc0(Dynamics *dyn, double omg_orb0,double r0_kepl);
double eob_dyn_bisecHam0(Dynamics *dyn, double pr0PN, double j0, double Hap, double rma);
int eob_dyn_rpr(const gsl_vector *x, void * params, gsl_vector *f);
int eob_dyn_rootfind_rpr(Dynamics *dyn, double *r0, double *pr0, double omg_orb0, double rKepl, double pr0PN);


/* TEOBResumSMetric.c */
extern void (*eob_metric_Apotential)(); /* defined in TEOBResumSPars.c*/
extern void (*eob_metric_Dpotential)(); /* defined in TEOBResumSPars.c*/
extern void (*eob_metric_Qpotential)(); /* defined in TEOBResumSPars.c*/
void eob_metric_A5PNlog(double r, double nu, double *A, double *dA, double *d2A);
void eob_metric_AGSF(double r, double nu, double *A, double *dA, double *d2A);
void eob_metric_A5PNlogP33(double r, double nu, double *A, double *dA, double *d2A);
void eob_metric_D3PN(double r, double nu, double *D, double *dD, double *d2D);
void eob_metric_DGSF(double r, double nu, double *D, double *dD, double *d2D);
void eob_metric_D5PNP32(double r, double nu, double *D, double *dD, double *d2D);
void eob_metric_Q3PN(double r, double prstar, double nu, double *Q, double *dQ_du, double *dQ_dprstar, double *d2Q_du2, double *ddQ_drdprstar, double *d2Q_dprstar2, double *d3Q_du2dprstar, double *d3Q_dudprstar2, double *d3Q_dprstar3);
void eob_metric_QGSF(double r, double prstar, double nu, double *Q, double *dQ_du, double *dQ_dprstar, double *d2Q_du2, double *ddQ_drdprstar, double *d2Q_dprstar2, double *d3Q_du2dprstar, double *d3Q_dudprstar2, double *d3Q_dprstar3);
void eob_metric_Q5PNloc(double r, double prstar, double nu, double *Q, double *dQ_du, double *dQ_dprstar, double *d2Q_du2, double *ddQ_drdprstar, double *d2Q_dprstar2, double *d3Q_du2dprstar, double *d3Q_dudprstar2, double *d3Q_dprstar3);
void eob_metric_Atidal(double r, Dynamics *dyn, double *AT, double *dAT, double *d2AT);
void eob_metric_Btidal(double r, Dynamics *dyn, double *BT, double *dBT, double *d2BT);
void eob_metric(double r, double prstar, Dynamics *dyn, double *A, double *B, double *dA, double *d2A, double *dB, double *d2B,
                double *Q, double *dQ, double *dQ_dprstar, double *d2Q, double *ddQ_drdprstar, double *d2Q_dprstar2,
                double *d3Q_dr2dprstar, double *d3Q_drdprstar2, double *d3Q_dprstar3);
void eob_metric_s(double r, double prstar, Dynamics *dyn, double *A, double *B, double *dA, double *d2A, double *dB, double *d2B,
                  double *Q, double *dQ, double *dQ_dprstar, double *d2Q, double *ddQ_drdprstar, double *d2Q_dprstar2,
                  double *d3Q_dr2dprstar, double *d3Q_drdprstar2, double *d3Q_dprstar3);

/* TEOBResumSFlux.c */
extern double (*eob_flx_Fr)(); /* defined in TEOBResumSPars.c*/
double eob_flx_Flux(double x, double Omega, double r_omega, double E, double Heff, double jhat, double r, double pr_star, double ddotr, Dynamics *dyn);
double eob_flx_Flux_s(double x, double Omega, double r_omega, double E, double Heff, double jhat, double r, double pr_star, double ddotr, Dynamics *dyn);
void eob_flx_Flux_ecc(double x, double Omega, double r_omega, double E, double Heff, double jhat, double r, double pr_star, double pphi, double rdot, double ddotr, double *Fphi, double *Fr, Dynamics *dyn);
void eob_flx_Fphi_ecc(double r, double prstar, double pphi, double Omg, double rdot, double Fphi, double Fphi_lo, double *FlmNewt, double *Flm, double Fphi_H, double Fr, Dynamics *dyn, double *hatflm_NC);
double eob_flx_Fr_ecc(double r, double prstar, double pphi, Dynamics *dyn, double Fphi_qc);
double eob_flx_Fr_ecc_BD(double r, double prstar, double pphi, Dynamics *dyn, double Fphi_qc);
double eob_flx_Fr_ecc_next(double r, double prstar, double pphi, Dynamics *dyn, double Fphi_qc);
double Fphi_NewtPref(double r, double Omg, double rdot, double r2dot, double r3dot, double Omgdot, double Omg2dot);
void eob_flx_Tlm(double w, double *MTlm);
void eob_flx_FlmNewt(double x, double nu, double *Nlm);
double eob_flx_HorizonFlux(double x, double Heff, double jhat, double nu);
double eob_flx_HorizonFlux_s(double x, double Heff, double jhat, double nu, double X1, double X2, double chi1, double chi2);

/* TEOBResumSWaveform.c */
extern void (*eob_wav_hlm)(); /* defined in TEOBResumSPars.c */
void eob_wav_hlm_circ(Dynamics *dyn, Waveform_lm_t *hlm);
void eob_wav_hlm_ecc(Dynamics *dyn, Waveform_lm_t *hlm);
void eob_wav_hlm_ecc_sigmoid(Dynamics *dyn, Waveform_lm_t *hlm);
extern void (*eob_wav_deltalm)(); /* defined in TEOBResumSPars.c */
void eob_wav_deltalm_v1(double Hreal,double Omega,double nu, double *dlm);
void eob_wav_deltalm_HM(double Hreal,double Omega,double nu, double *dlm);
void eob_wav_hlmNewt_ecc(Dynamics *dyn,Waveform_lm_t *hlmNewt);
void eob_wav_hlmNewt_ecc_sigmoid(Dynamics *dyn,Waveform_lm_t *hlmNewt);
void eob_wav_hhatlmTail(double Omega,double Hreal,double bphys, Waveform_lm_t *tlm);
void eob_wav_speedyTail(double Omega, double Hreal, double bphys, Waveform_lm_t *tlm);
extern void (*eob_wav_hlmNewt)(); /* defined in TEOBResumSPars.c */
void eob_wav_hlmNewt_v1(double r, double Omega, double phi, double nu, Waveform_lm_t *hNewt);
void eob_wav_hlmNewt_HM(double r, double Omega, double phi, double nu, Waveform_lm_t *hNewt);
void eob_wav_hlmTidal(double x, Dynamics *dyn, double *hTidallm);
extern void (*eob_wav_flm)(); /* defined in TEOBResumSPars.c */
void eob_wav_flm_v1(double x,double nu, double *rholm, double *flm);
void eob_wav_flm_old(double x,double nu, double *rholm, double *flm);
void eob_wav_flm_HM(double x,double nu, double *rholm, double *flm);
void eob_wav_flm_HM_4PN22(double x,double nu, double *rholm, double *flm);
void eob_wav_flm_Kerr(double x,double nu, double *rholm, double *flm);
extern void (*eob_wav_flm_s)(); /* defined in TEOBResumSPars.c */
void eob_wav_flm_s_SSNLO(double x, double nu, double X1, double X2, double chi1, double chi2, double a1, double a2, double C_Q1, double C_Q2, int usetidal, double *rholm, double *flm);
void eob_wav_flm_s_SSLO(double x, double nu, double X1, double X2, double chi1, double chi2, double a1, double a2, double C_Q1, double C_Q2, int usetidal, double *rholm, double *flm);
void eob_wav_flm_s_old(double x, double nu, double X1, double X2, double chi1, double chi2, double a1, double a2, double C_Q1, double C_Q2, int usetidal, double *rholm, double *flm);
void eob_wav_flm_s_HM(double x, double nu, double X1, double X2, double chi1, double chi2, double a1, double a2, double C_Q1, double C_Q2, int usetidal, double *rholm, double *flm);
void eob_wav_flm_s_HM_4PN22(double x, double nu, double X1, double X2, double chi1, double chi2, double a1, double a2, double C_Q1, double C_Q2, int usetidal, double *rholm, double *flm);
void eob_wav_flm_s_Kerr(double x, double nu, double X1, double X2, double chi1, double chi2, double a1, double a2, double C_Q1, double C_Q2, int usetidal, double *rholm, double *flm);
extern void (*eob_wav_hlmNQC_find_a1a2a3)(); /* defined in TEOBResumSPars.c */
void eob_wav_hlmNQC_find_a1a2a3_circ(Dynamics *dyn, Waveform_lm *h, Waveform_lm *hnqc);
void eob_wav_hlmNQC_find_a1a2a3_ecc(Dynamics *dyn, Waveform_lm *h, Waveform_lm *hnqc);
extern void (*eob_wav_hlmNQC_find_a1a2a3_mrg)(); /* defined in TEOBResumSPars.c */
void eob_wav_hlmNQC_find_a1a2a3_mrg_HM(Dynamics *dyn_mrg, Waveform_lm *hlm_mrg, Waveform_lm *hnqc,
				       Dynamics *dyn, Waveform_lm *hlm);
void eob_wav_hlmNQC_find_a1a2a3_mrg_22(Dynamics *dyn_mrg, Waveform_lm *hlm_mrg, Waveform_lm *hnqc,
				       Dynamics *dyn, Waveform_lm *hlm);
void eob_wav_hlmNQC_find_a1a2a3_mrg_ecc(Dynamics *dyn_mrg, Waveform_lm *hlm_mrg, Waveform_lm *hnqc,
				       Dynamics *dyn, Waveform_lm *hlm);
void eob_wav_hlmNQC(double  nu, double  r, double  prstar, double  Omega, double  ddotr, NQCcoefs *nqc, Waveform_lm_t *hlmnqc);
void eob_wav_hlmNQC_ecc_sigmoid(double  nu, double  r, double  prstar, double  Omega, double  ddotr, double t, double tOmg_pk, NQCcoefs *nqc, Waveform_lm_t *hlmnqc);
void eob_wav_hlmNQC_nospin201602(double  nu, double  r, double  prstar, double  Omega, double  ddotr, Waveform_lm_t *hlmnqc);
void eob_wav_ringdown_template(double x, double a1, double a2, double a3, double a4, double b1, double b2, double b3, double b4, double sigmar, double sigmai, double *psi);
extern int (*eob_wav_ringdown)(); /* defined in TEOBResumSPars.c */
int eob_wav_ringdown_v1(Dynamics *dyn, Waveform_lm *hlm);
int eob_wav_ringdown_HM(Dynamics *dyn, Waveform_lm *hlm);
double eob_wav_hlmTidal_fmode_fact22A(double x, double alpha, double bomgf, double XB);
void SPA(Waveform_lm *TDlm, WaveformFD_lm *FDlm);
void twist_hlm_TD(Dynamics *dyn, Waveform_lm *hlm, DynamicsSpin *spin, int interp_spin_abc, Waveform_lm *hTlm, Waveform_lm *hTlm_neg, Waveform_lm *hTl0);
void twist_hlm_FD(WaveformFD_lm *hlm, DynamicsSpin *spin, double M, double amplitude_prefactor, double phi, double iota, WaveformFD *hpc);
void compute_hpc_old(Waveform_lm *hlm, double nu, double M, double distance, double amplitude_prefactor, double psi, double iota, Waveform *hpc);
void compute_hpc(Waveform_lm *hlm,Waveform_lm *hlm_neg, Waveform_lm *hl0, double nu, double M, double distance, double amplitude_prefactor, double psi, double iota, Waveform *hpc);
void SPA(Waveform_lm *TDlm, WaveformFD_lm *FDlm);
void prolong_euler_angles_TD(double *alpha, double *beta, double *gamma, Dynamics *dyn, DynamicsSpin *spin, Waveform_lm *hlm);
void prolong_euler_angles(double *alpha, double *beta, double *gamma, Dynamics *dyn, DynamicsSpin *spin, Waveform_lm *hlm);
void compute_hpc_FD(WaveformFD_lm *hlm, double nu, double M, double distance, double amplitude_prefactor, double phi, double iota, WaveformFD *hpc);
double time_shift_mrg_to_0(Waveform_lm *hlm);
void time_shift_TD(double *t, double tc, int size);
void time_shift_FD(WaveformFD *hpc, double tc);
void eob_wav_hlmNQC_test_bhns(Dynamics *dyn_mrg, Waveform_lm *hlm_mrg, Waveform_lm *hnqc,
				       Dynamics *dyn, Waveform_lm *hlm, bool *bhns);
void eob_wav_hlmNQC_find_a1a2a3_mrg_BHNS_HM(Dynamics *dyn_mrg, Waveform_lm *hlm_mrg, Waveform_lm *hnqc,
				       Dynamics *dyn, Waveform_lm *hlm);
void eob_wav_ringdown_template_td(double x, double a1, double a2, double a3, double a4, double b1, double b2, double b3, double b4, double sigmai, double *psi, double alpha2, double Amrg);
int eob_wav_ringdown_bhns(Dynamics *dyn, Waveform_lm *hlm);

/* SpecialFuns.c */
double fact(int n);
double doublefact(int n);
double wigner_d_function(int l, int m, int s, double i);
int spinsphericalharm(double *rY, double *iY, int s, int l, int m, double phi, double i);
double Fresnel_Sine_Integral( double x );
double Fresnel_Cosine_Integral( double x );

#ifdef _OPENMP
/* TEOBResumSOMP.c */
void openmp_init(const int verbose);
void openmp_timer_start(char *name);
void openmp_timer_stop(char *name);
void openmp_free();
#endif
