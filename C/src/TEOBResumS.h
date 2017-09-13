/**
 * Copyright (C) 2017 Alessandro Nagar, Gregorio Carullo, Ka Wa Tsang, Philipp Fleig, Sebastiano Bernuzzi, Walter Del Pozzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#include <cmath>
#include <math.h>
#include <vector>
#include <limits>
#include <stdio.h>
#include <cstring> 
#include <libconfig.h> /* library to manage parameters */

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_spline.h>

/** Macros */
#define ERROR 1 /** generic error int */
#define OK 0 /** generic go int */
#define STRLEN 128 /** Standard string length */
#define TEOBResumS_Info "TEOBResumS code (C) 2017\n"
#define TEOBResumS_Usage {printf("USAGE:\n\t%s parfile\n", argv[0]);exit(OK);} 
#define DEBUG 1 /** Flag for debug mode */
#ifndef PR /** Flag for print option (control at compiling time) */
#define PR 0 
#endif

#define SIGN(x,y) ((y) >= 0.0 ? fabs(x) : -fabs(x)) 
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
  
/** List of EOB evolved variables */
enum{
    EOB_EVOLVE_RAD, 
    EOB_EVOLVE_PHI,
    EOB_EVOLVE_PRSTAR,
    EOB_EVOLVE_PPHI,
    EOB_EVOLVE_VARS
  };

/** List of EOB variables for initial data */
enum{
    EOB_ID_RAD, 
    EOB_ID_PPH,
    EOB_ID_PRSTAR,
    EOB_ID_PR,
    EOB_ID_J,
    EOB_ID_E0,
    EOB_ID_OMGJ,
    EOB_ID_VARS
};

/** List of EOB dynamical variables (to be stored in arrays) */ 
enum{
  EOB_RAD, 
  EOB_PHI,
  EOB_PPH,
  EOB_MOMG,
  EOB_DDOTR,
  EOB_PRSTAR,
  EOB_OMGORB,
  EOB_DYNAMICS_VARS
};

#define KMAX 35 /** Multipolar linear index, max value */

/** Maps between linear index and the corresponding (l, m) multipole indices */
const int L[KMAX] = {
    2,2,
    3,3,3,
    4,4,4,4,
    5,5,5,5,5,
    6,6,6,6,6,6,
    7,7,7,7,7,7,7,
    8,8,8,8,8,8,8,8};

const int M[KMAX] = {
    1,2,
    1,2,3,
    1,2,3,4,
    1,2,3,4,5,
    1,2,3,4,5,6,
    1,2,3,4,5,6,7,
    1,2,3,4,5,6,7,8};

/** Waveform data type */
typedef struct tagWaveform
{
  int size;
  double *time;
  double *real;
  double *imag;
  //complex *data; // what's best to use?
  char name[STRLEN];
}  Waveform;

/** Multipolar waveform data type */
typedef struct tagWaveform_lm
{
  int size;
  double *time;
  double *real[KMAX];
  double *imag[KMAX];
  //complex *data[KMAX]; // what's best to use?
  char name[KMAX][STRLEN];
  int *kmask[KMAX]; /* mask for multipoles */
}  Waveform_lm;

/** Dynamics data type */
typedef struct tagDynamics
{
  double t, r, prstar, phi, pphi, ddotr, Omg, Omg_orb, A;
  double y[EOB_EVOLVE_VARS]; /* rhs storage */
  double y0[EOB_ID_VARS]; /* ID storage */
  double t1, dt, t_stop, ti;
  double rLR, rLSO, MOmg, MOmg_prev;
  bool stop_flag, MOmgpeak_flag;
  int size;
  double *time;
  double *data[EOB_DYNAMICS_VARS]; 
} Dynamics;


/* Function protoypes grouped based on file */


/* TEOBResumSPars.c */
void par_db_init ();
void par_db_free ();
void par_db_default ();
void par_file_parse (char *fname);
void par_db_write_file (char *fname);
void par_db_screen ();
void par_set_i(const char *key, int *val);
void par_set_b(const char *key, int *val);
void par_set_d(const char *key, double *val);
void par_set_s(const char *key, const char **val);
int par_get_i(const char *key);
int par_get_b(const char *key);
double par_get_d(const char *key);
const char * par_get_s(const char *key);
int * par_get_arrayi(const char *key, int *n);
double * par_get_arrayd(const char *key, int *n);
//void par_commandline_parse(char *s, int n);
//void TEOBResumSSetParameters(char *s, int n, int mode, int pr);
void TEOBResumSSetParameters(char *s, int pr);

/* TEOBResumSUtil.c */
double q_to_nu(const double q);
double nu_to_X1(const double nu);
double Eulerlog(const double x,const double m);
int find_point_bisection(double x, int n, double *xp, int o);
double baryc_f(double xx, int n, double *f, double *x);
void baryc_weights(int n, double *x, double *omega);
double baryc_f_weights(double xx, int n, double *f, double *x, double *omega);
double interp1d (const int order, double xx, int nx, double *f, double *x);
double fact(int n);
double wigner_d_function(int l, int m, int s, double i);
int spinsphericalharm(double *rY, double *iY, int s, int l, int m, double phi, double i);
int D0(double *f, double dx, int n, double *df);
int D2(double *f, double dx, int n, double *d2f);
int D0_nux(double *f, double *x, int n, double *df);
void set_multipolar_idx_mask(int *kmask, int n);
void Waveform_alloc (Waveform **wav, int size, char *name);
void Waveform_push (Waveform **wav, int size);
void Waveform_free (Waveform *wav);
void Waveform_lm_alloc (Waveform_lm **wav, int size, char **name);
void Waveform_lm_push (Waveform **wav, int size);
void Waveform_lm_free (Waveform_lm *wav);
void Dynamics_alloc (Dynamics **dyn, int size);
void Dynamics_push (Dynamics **dyn, int size);
void Dynamics_output (Dynamics *dyn);
void Dynamics_free (Dynamics *dyn);
void errorexit(char *file, int line, char *s);
#define errorexit(s) errorexit(__FILE__, __LINE__, (s))
void errorexits(char *file, int line, char *s, char *t);
#define errorexits(s,t) errorexits(__FILE__, __LINE__, (s), (t))

/* TEOBResumSFits.c */
double c3_fit_global(double nu, double chi1, double chi2, double X1, double X2, double a1, double a2);
double logQ(double x);
double Yagi13_fit_barlamdel(double barlam2, int ell);

/* TEOBResumSDynamics.c */
int rhs(double t, const double y[], double dy], void *params);
int s_RHS(double t, const double y[], double dy[], void *params);
void s_GS(double r, double rc, double drc_dr, double aK2, double prstar, double pph, double nu, double chi1, double chi2, double X1, double X2, double cN3LO,
	  double *ggm);
void s_get_rc(double r, double nu, double at1,double at2, double aK2, double C_Q1, double C_Q2, int usetidal, double *rcout);
void get_Omg_orb(double *r, double *pph, double *pr_star, double *A, double *B, int size, void *params, double *Omg_orb);

/* TEOBResumSMetric.c */

/* TEOBResumSFlux.c */

/* TEOBResumSWaveform.c */





// ********************************************
// ********************************************
// ********************************************
// ********************************************
// ********************************************
// ********************************************
// ********************************************

// ARRIVED HERE 

















/** Structure of flags to control various physics ingredients */
typedef struct tagTEOBResumFlags
{
    int solver_scheme;  /** Solver scheme to utilise */
    int spin;           /** Spinning dynamics */
    int tidal;          /** Tidal deformability dynamics */
    int RWZ;            /** Regge-Wheeler-Zerilli potential */
    int speedy;         /** Faster tails calculation */
    int dynamics;       /** Output dynamics to file */
    int Yagi_fits;      /** Use 'universal relations' for higher l lambdas */
    int multipoles;     /** Output single multipole waveforms */
    int geometric_units; /** use geometric units */
    int set;            /** flags set bby default or not? */
}   TEOBResumFlags;

/** Structure of parameters */
typedef struct tagTEOBResumParams
{
    int    lm;                  /**                                                          */
    double mtot;                /** Total mass of the binary                                 */
    double q;                   /** Mass ratio of the binary                                 */
    double iota;                /** inclination angle                                        */
    double psi;                 /** polarisation angle                                       */
    double distance;            /** sources distance [Mpc]                                   */
    double nu;                  /** Symmetric mass ratio of the binary, nu = m1*m2/(m1+m2)^2 */
    double r0;                  /** Initial radial separation of the objects                 */
    double f_min;               /** Initial frequency                                        */
    double dt;                  /** Time step of the differential equation evolution         */
    double rLR;                 /** Light ring radius                                        */
    double chi1;                /** Dimensionless spin of the first object                   */
    double chi2;                /** Dimensionless spin of the second object                  */
    double S1;                  /**                                                          */
    double S2;                  /**                                                          */
    double S;                   /**                                                          */
    double Sstar;               /**                                                          */
    double X1;                  /**                                                          */
    double X2;                  /**                                                          */
    double a1;                  /**                                                          */
    double a2;                  /**                                                          */
    double aK;                  /**                                                          */
    double aK2;                 /**                                                          */
    double cN3LO;               /**                                                          */
    double LambdaAl2;           /** Tidal coupling constant of the first object  for l=2     */
    double LambdaAl3;           /** Tidal coupling constant of the first object  for l=3     */
    double LambdaAl4;           /** Tidal coupling constant of the first object  for l=4     */
    double LambdaBl2;           /** Tidal coupling constant of the second object for l=2     */
    double LambdaBl3;           /** Tidal coupling constant of the second object for l=3     */
    double LambdaBl4;           /** Tidal coupling constant of the second object for l=4     */
    double kappaAl2;            /**                                                          */
    double kappaAl3;            /**                                                          */
    double kappaAl4;            /**                                                          */
    double kappaBl2;            /**                                                          */
    double kappaBl3;            /**                                                          */
    double kappaBl4;            /**                                                          */
    double kappaTl2;            /**                                                          */
    double kappaTl3;            /**                                                          */
    double kappaTl4;            /**                                                          */
    double bar_alph2_1;         /**                                                          */
    double bar_alph2_2;         /**                                                          */
    double bar_alph3_1;         /**                                                          */
    double bar_alph3_2;         /**                                                          */
    double C_Q1;                /**                                                          */
    double C_Q2;                /**                                                          */
    double Mbh;                 /** Mass of the final BH                                     */
    double Abh;                 /** Spin of the final BH                                     */
    TEOBResumFlags flags;       /** Flag structure                                           */
}   TEOBResumParams;





/** Sets params flags to their default value */
void SetDefaultFlagsValues(TEOBResumFlags *flags);

/** Copy params flags */
void CopyTEOBResumSFlags(TEOBResumFlags *out, TEOBResumFlags *in);

/** _A_NumDenom_h */

int A_NumDenom(const double r, const vector<double> a, const double nu);

/* _AdiabLR_h */

double fLR(double r, void *params);
double AdiabLR(void *params);

/** Eulerlog */

double Eulerlog(const double x, const double m);


/* _FlmNewt_h */

vector<double> FlmNewt(const double x, void *params);

/* _HealyBBHFitRemnant_h */

double HealyBBHFitRemnant(double chi1, double chi2, double q);

double JimenezFortezaRemnantSpin(TEOBResumParams params);

/* _HorizonFlux_h */

double HorizonFlux(const double x, const double Heff, const double jhat, const double nu);

/* _LALSimIMRTEOBIHES_h */

void TEOBResumS(Waveform **hplus,               /** h+ return array **/
                        Waveform **hcross,      /** hx return array **/
                        double m1,              /** m1(Msun) **/
                        double m2,              /** m2(Msun) **/
                        double spin1x,          /** dimensionless s1x **/
                        double spin1y,          /** dimensionless s1y **/
                        double spin1z,          /** dimensionless s1z **/
                        double spin2x,          /** dimensionless s2x **/
                        double spin2y,          /** dimensionless s2y **/
                        double spin2z,          /** dimensionless s2z **/
                        double inclination,     /** inclination angle (rad) **/
                        double polarisation,    /** polarisation angle (rad) **/
                        double f_min,           /** starting frequency(Hz) **/
                        double dt,              /** sampling interval (s) **/
                        double LambdaAl2,       /** l=2 (tidal deformation of body 1)/(mass of body 1)^5 **/
                        double LambdaBl2,       /** l=2 (tidal deformation of body 2)/(mass of body 2)^5 **/
                        double LambdaAl3,       /** l=3 (tidal deformation of body 1)/(mass of body 1)^5 **/
                        double LambdaBl3,       /** l=3 (tidal deformation of body 2)/(mass of body 2)^5 **/
                        double LambdaAl4,       /** l=4 (tidal deformation of body 1)/(mass of body 1)^5 **/
                        double LambdaBl4,       /** l=4 (tidal deformation of body 2)/(mass of body 2)^5 **/
                        double distance,        /** distance(Mpc) **/
                        int    lm,              /** multipole index for output **/
                        TEOBResumFlags *flags   /** flags **/
);

/* _Metric_h */
    
vector<double> acoeffs(const double r, const double nu);
vector<double> Metric(const double r, void *params, bool nnlo_flag);
vector<double> A5pnP15_dd(const double r, void *params);

    
/* _QNMHybridFitCab_h */
    
void QNMHybridFitCab(TEOBResumParams params, vector<double> &a1, vector<double> &a2, vector<double> &a3, vector<double> &a4, vector<double> &b1, vector<double> &b2, vector<double> &b3, vector<double> &b4, vector<gsl_complex> &sigma);
    
/* _RHS_h */

int rhs(double t, const double y[], double f[], void *params);

/* _SpeedyTail_h */

vector<gsl_complex> speedyTail(const double Omega, const double Hreal, const double bphys, const int L[], const int M[]);

/* _Tlm_h */

vector<double> Tlm(const double w);

/* _deltalm_h */

vector<double> deltalm(const double Hreal, const double Omega, const double nu);

/* _dtnqc_fit_h */

double dtnqc_fit(double chi, double chi0);

/* _flm_h */

vector<double> f_lm(const double x,const double nu);

/* _find_a1a2a3_h */

vector<vector<gsl_complex> > find_a1a2a3(vector<double> T, vector<double> r, vector<double> w, vector<double> pph, vector<double> pr_star, vector<vector<double> > hlm_phase,vector<double> Omg_orb, vector<vector<double> > A, vector<double> ddotr, void *params);

/* _flux_h */

double flux(const double x,const double Omega,const double r_omega,const double E, const double Heff,const double jhat,const double r,const double prstar, const double ddotr, double source[],void *params);

/* _get_Omg_orb_h */

vector<double> get_Omg_orb(vector<double> r, vector<double> pph, vector<double> pr_star, vector<double> A, vector<double> B, void *params);

/* _hhatlmtail_h */

vector<gsl_complex> hhatlmtail(const double Omega, const double Hreal, const double bphys, const int L[], const int M[]);

/* _hlm_h */

vector<gsl_complex> hlm(double t, const double phi, const double r, const double pph, const double prstar, double Omega, const double ddotr, const double H, const double Heff,const double jhat, const double rw,void *params);

/* _hlm_NQC_h */

vector<gsl_complex> hlmNQC(double nu, double r, double prstar, double  Omega, double ddotr);

/* _hlmNewt_h */

vector<gsl_complex> hlmNewt(const double r, const double Omega, const double phi, const double nu, bool tidal_flag);

/* _hlm_Tidal_h */

vector<double> hlm_Tidal(double x,void *params);

/* _initial_h */

vector<double> initial(TEOBResumParams *params);

/* _input_struc_h */



/* _interp_h */

vector<double> interp_grid(vector<double> t_vec, vector<double> data, double dt);

/* _interpolator_h */

double interpolate(double dt, vector<gsl_complex> grid);

double interp1d (const int order, double xx, int nx, double *f, double *x);
double baryc_f_weights(double xx, int n, double *f, double *x, double *omega);
void baryc_weights(int n, double *x, double *omega);
double baryc_f(double xx, int n, double *f, double *x);
int find_point_bisection(double x, int n, double *xp, int o);

/* _read_config_h */

double logQ(double x);

TEOBResumParams process_input_parameters(double m1,
                                double m2,
                                double chi1,
                                double chi2,
                                double f_min,
                                double sampling_rate,
                                double LambdaAl2,
                                double LambdaBl2,
                                double LambdaAl3,
                                double LambdaBl3,
                                double LambdaAl4,
                                double LambdaBl4,
                                TEOBResumFlags *flags
                                );
TEOBResumParams read_config(char *fname);

double Yagi13_fit_barlamdel(double barlam2, int ell);
double time_units_conversion(double M, double dt);
double radius0(double M, double f_start);

/* _ringdown_h */

int ringdown(TEOBResumParams params, vector<vector<double> > &t_vec, vector<double> Omega_vec, vector<vector<double> > &hlm_rad, vector<vector<double> > &hlm_phase);

/* _ringdown_match_h */

gsl_complex ringdown_match(double x, int k, vector<double> a1, vector<double> a2, vector<double> a3, vector<double> a4, vector<double> b1, vector<double> b2, vector<double> b3, vector<double> b4, vector<gsl_complex> sigma);

/* _s_A5PNlog_h */

vector<double> s_A5PNlog(double r, void *params, bool nnlo_flag);

/* _s_D1_h */

vector<double> s_D1(vector<double> f, vector<double> x, int Nmax);
vector<double> u_D1(vector<double> f, vector<double> x, int Nmax);

/* _s_Flux_h */

double s_Flux(double x, double Omega, double r_omega, double E, double Heff, double jhat, double r, double pr_star, double ddotr, void *params);

/* _s_GS_h */

double c3_fit_global(double nu, double chi1, double chi2, double X1, double X2, double a1, double a2, bool tidal_flag);
vector<double> s_GS(double r, double rc, double drc_dr, double aK2, double prstar, double pph, double nu, double chi1, double chi2, double X1, double X2, double cN3LO);

/* _s_HorizonFlux_h */

double s_HorizonFlux(double x, double Heff, double jhat, double nu, double X1, double X2, double chi1, double chi2);

/* _s_Metric_h */

vector<double> s_Metric(double r, void *params, bool nnlo_flag);

/* _s_RHS_h */

int s_RHS(double t, const double y[], double f[], void *params);

/* _s_bisec_h */

struct energy_params
{
    double rorb, A, dA, rc, drc_dr, ak2, S, Ss, nu, chi1, chi2, X1, X2, c3;
};


double energy_params (double x, void *params);
double DHeff0(double x, void *DHeff_params);
double s_bisec(double pph, double rorb, double A, double dA, double rc, double drc_dr, double ak2, double S, double Ss, void *params);

/* _s_flm_h */

vector<double> s_flm(double x, void *params);

/* _s_get_rc_h */

vector <double> s_get_rc(double r, void *params);

/* _s_initial_h */

vector<double> s_initial(TEOBResumParams *params);

/* _s_waveform_h */

vector<gsl_complex> s_waveform(double t, const double y[], void *params, double &Omg, double &Omg_orb, double &A, double &ddotr);

/* _spinsphericalharm_h */

double fact(int n);
double wigner_d_function(int l, int m, int s, double x);
void spinsphericalharm(double *rY, double *iY, int s, int l, int m, double phi, double x);






