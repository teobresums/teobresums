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

#ifndef _TEOBRESUM_H
#define _TEOBRESUM_H
//#ifdef __cplusplus
//extern "C" {
//#endif

#define DEBUG 0

#include <vector>
#include "cmath"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

using namespace::std;

/** Useful constants */
const double EulerGamma = M_EULER;
const double Log2       = M_LN2;
const double pi         = M_PI;
const double MSUN_M     = 1.47662504e3; /** G/c^3 */
const double MSUN_S     = 4.925491025543575903411922162094833998e-6; /** G/c^2 */
const double MPC_M      = 3.086e22;

/** Maps between 1D index and the corresponding l and m multipole indices */
const int L[35] = {
    2,2,
    3,3,3,
    4,4,4,4,
    5,5,5,5,5,
    6,6,6,6,6,6,
    7,7,7,7,7,7,7,
    8,8,8,8,8,8,8,8};

const int M[35] = {
    1,2,
    1,2,3,
    1,2,3,4,
    1,2,3,4,5,
    1,2,3,4,5,6,
    1,2,3,4,5,6,7,
    1,2,3,4,5,6,7,8};

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

/** Algorithm control structure */
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

/** Waveform structure, container for an array of doubles containing h and its length */
typedef struct tagWaveform
{
    double       *data;
    unsigned int length;
}   Waveform;

/** Sets the dynamics controlling flags to their default value */
void SetDefaultFlagsValues(TEOBResumFlags *flags);

void CopyTEOBResumSFlags(TEOBResumFlags *out, TEOBResumFlags *in);

/* _A_NumDenom_h */

vector<double> A_NumDenom(const double r, const vector<double> a, const double nu);

/* _AdiabLR_h */

double fLR(double r, void *params);
double AdiabLR(void *params);

/* _Eulerlog_h */

double Eulerlog(const double x, const double m);

/* _FDdrvt_h */

vector<double> FDdrvt(const vector<double> f, const vector<double> t, int varargin, const int grid_length);

/* _FDdrvt_omega_h */

vector<double> FDdrvt_omega(vector<double> f, double dt);

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

//#ifdef __cplusplus
//}
//#endif
#endif /* _TEOBRESUM_H */

