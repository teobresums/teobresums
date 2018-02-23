/**
 * Copyright (C) 2017 Sebastiano Bernuzzi, Gregorio Carullo, Walter Del Pozzo, Alessandro Nagar, Ka Wa Tsang
 * This file is part of TEOBResumS
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
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 */

#include "TEOBResumS.h"

/** Return symm mass ratio from q */
double q_to_nu(const double q)
{
  double nu = 0;
  if (q>0.)
    nu = q/((q+1.)*(q+1.));
  return nu;
}

/** Return mass ratio M1/M from nu */
double nu_to_X1(const double nu)
{
  if ( (nu<0.) || (nu>0.25) )
    errorexit("sym mass ratio is defined 0<= nu <= 1/4.");
  return 0.5*(1.+sqrt(1.-4.*nu));
}

/** Eulerlog function */
double Eulerlog(const double x,const double m)
{
  /* const double EulerGamma = 0.5772156649015328606065121; */
  /* const double Log2       = 0.6931471805599453094172321; */
  /* above constants are defined in header */
  return EulerGamma + Log2 + log(m) + 0.5*log(x);
}

/** Spline interpolation with GSL routines */
void interp_grid(double *t, double *y, int n, double *ti, int ni, double *yi)
{
  int k;
  double yi_pt;
 
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n);
  gsl_spline_init (spline, t, y, n);
    
  for (k = 0; k < ni; k++) {
    yi[k] = gsl_spline_eval (spline, ti[k], acc);
  }
  
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  
}

/** Find nearest point index in 1d array */
int find_point_bisection(double x, int n, double *xp, int o)
{
  int i0 = o-1, i1 = n-o;
  int i;
  if (n < 2*o) {
    errorexit("not enough point to interpolate");
  }  
  if (x <= xp[i0]) return 0;
  if (x >  xp[i1]) return n-2*o;
  while (i0 != i1-1) {
    i = (i0+i1)/2;
    if (x < xp[i]) i1 = i; else i0 = i;
  }
  return i0-o+1;
}

/** Barycentric Lagrange interpolation at xx with n points of f(x), 
   equivalent to standard Lagrangian interpolation */   
#define tiny 1e-12
double baryc_f(double xx, int n, double *f, double *x)
{
  double omega[n];
  double o, num, den, div, ci;
  int i, j;
  for (i = 0; i < n; i++) {
    if (fabs(xx - x[i]) <= tiny) return f[i];
    o = 1.;
    for (j = 0; j < n; j++) {
      if (j != i) {
	o /= (x[i] - x[j]);
      }
    }
    omega[i] = o;
  }
  num = den = 0.;
  for (i = 0; i < n; i++) {
    div  = xx - x[i];
    ci   = omega[i]/div;
    den += ci;
    num += ci * f[i];
  }
  return( num/den );
}

/** Barycentric Lagrange interpolation at xx with n points of f(x), 
    compute weights */
void baryc_weights(int n, double *x, double *omega)
{  
  double o;
  int i, j;
  for (i = 0; i < n; i++) {
    o = 1.;
    for (j = 0; j < n; j++) {
      if (j != i) { 
	o /= (x[i] - x[j]);
      }
    }
    omega[i] = o;
  }
}

/** Barycentric Lagrange interpolation at xx with n points of f(x), 
    use precomputed weights */
double baryc_f_weights(double xx, int n, double *f, double *x, double *omega)
{
  int i;
  double num, den, div, ci;
  num = den = 0.;
  for (i = 0; i < n; i++) {
    div  = xx - x[i];
    if (fabs(div) <= tiny) return f[i];
    ci   = omega[i]/div;
    den += ci;
    num += ci * f[i];
  }  
  return( num/den );
}

/** 1d Lagrangian barycentric interpolation */
double interp1d (const int order, double xx, int nx, double *f, double *x)
{
  double ff;
  int ix;
  int ox = order > nx ? nx : order;
  ix = find_point_bisection(xx, nx, x, ox/2);
  ff = baryc_f(xx, ox, &f[ix], &x[ix]);  
  return( ff );
}

/** Factorial */
double fact(int n)
{
    double f[] = {1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880.,
        3628800., 39916800., 479001600., 6227020800., 87178291200.};
    if (n < 0){
      errorexit(" computing a negative factorial.\n");
    } else if (n <= 14){
      return f[n];
    } else {
      return n*fact(n-1);
    }
}

/** Wigner d-function */
double wigner_d_function(int l, int m, int s, double i)
{
  double dWig = 0.;
  double costheta = cos(i*0.5);
  double sintheta = sin(i*0.5);
  int ki = MAX( 0  , m-s );
  int kf = MIN( l+m, l-s );
  int k;
  for( k = ki; k <= kf; k++ ){
    dWig +=
      ( pow(-1.,k) * pow(costheta,2*l+m-s-2*k) * pow(sintheta,2*k+s-m) )/
      ( fact(k) * fact(l+m-k) * fact(l-s-k) * fact(s-m+k) );
  }
  return (sqrt(fact(l+m) * fact(l-m) * fact(l+s) * fact(l-s)) * dWig);
}

/** Spin-weighted spherical harmonic 
    Ref: https://arxiv.org/pdf/0709.0093.pdf */
int spinsphericalharm(double *rY, double *iY, int s, int l, int m, double phi, double i)
{
    if ((l<0) || (m<-l) || (m>l)) {
        errorexit(" wrong (l,m) inside spinspharmY\n");
    }
    double c = pow(-1.,-s) * sqrt( (2.*l+1.)/(4.*M_PI) );
    double dWigner = c * wigner_d_function(l,m,-s,i);
    *rY = cos((double)(m)*phi) * dWigner;
    *iY = sin((double)(m)*phi) * dWigner;
    return OK;
}

/** (h+, hx) polarizations from the multipolar waveform */
void compute_hpc(Waveform_lm **hlm, double nu, double M, double distance, double amplitude_prefactor, double psi, double iota, Waveform *hpc)
{
  static const int mneg = 1; /* m>0 modes only, add m<0 modes afterwards */
  double Y_real, Y_imag;
  double Aki, cosPhi, sinPhi;
  int k,i;
  for (k = 0; k < KMAX; k++ ) {
    spinsphericalharm(&Y_real, &Y_imag, -2, LINDEX[k], MINDEX[k], psi,iota);
    for (i = 0; i < (*hlm)->size; i++) {
      Aki    = amplitude_prefactor * (*hlm)->ampli[k][i];
      cosPhi = cos( (*hlm)->phase[k][i] );
      sinPhi = sin( (*hlm)->phase[k][i] );
      hpc->real[i] += Aki*(cosPhi*Y_real + sinPhi*Y_imag);
      hpc->imag[i] -= Aki*(sinPhi*Y_real + cosPhi*Y_imag); // overall check sign
      hpc->time[i] *= M; 
    }
    if ( (mneg) && (MINDEX[k]!=0) ) { 
      /* add m<0 modes */
      spinsphericalharm(&Y_real, &Y_imag, -2, LINDEX[k], -MINDEX[k], psi,iota);
      for (i = 0; i < (*hlm)->size; i++) {
	Aki    = amplitude_prefactor * (*hlm)->ampli[k][i];
	cosPhi = cos( (*hlm)->phase[k][i] ); 
	sinPhi = sin( (*hlm)->phase[k][i] ); 
	hpc->real[i] += Aki*(cosPhi*Y_real - sinPhi*Y_imag);
	hpc->imag[i] += Aki*(sinPhi*Y_real - cosPhi*Y_imag); // overall check sign
      }     
    } /* mneg */
  }
}

/** 4th order centered stencil first derivative, uniform grids */
int D0(double *f, double dx, int n, double *df)
{
  const double oo12dx  = 1./(12*dx);
  int i;
  for (i=2; i<n-2; i++) {
    df[i] = (8.*(f[i+1]-f[i-1]) - f[i+2] + f[i-2])*oo12dx;
  }
  i = 0;
  df[i] = (-25.*f[i] + 48.*f[i+1] - 36.*f[i+2] + 16.*f[i+3] - 3.*f[i+4])*oo12dx;
  i = 1;
  df[i] = (-3.*f[i-1] - 10.*f[i] + 18.*f[i+1] - 6.*f[i+2] + f[i+3])*oo12dx;
  i = n-2;
  df[i] = - (-3.*f[i+1] - 10.*f[i] + 18.*f[i-1] - 6.*f[i-2] + f[i-3])*oo12dx;
  i = n-1;
  df[i] = - (-25.*f[i] + 48.*f[i-1] - 36.*f[i-2] + 16.*f[i-3] - 3.*f[i-4])*oo12dx;
  return OK;
}

/** 4th order centered stencil second derivative, uniform grids */
int D2(double *f, double dx, int n, double *d2f)
{
  const double oo12dx2  = 1./(dx*dx*12);
  int i;
  for (i=2; i<n-2; i++) {
    d2f[i] = (-30*f[i]+16*(f[i+1]+f[i-1])-(f[i+2]+f[i-2]))*oo12dx2;
  }
  i= 0;
  d2f[i] = (45*f[i]-154*f[i+1]+214*f[i+2]-156*f[i+3]+61*f[i+4]-10*f[i+5])*oo12dx2;
  i= 1;
  d2f[i] = (10*f[i-1]-15*f[i]-4*f[i+1]+14*f[i+2]-6*f[i+3]+f[i+4])*oo12dx2;
  i = n-2;
  d2f[i] = (10*f[i+1]-15*f[i]-4*f[i-1]+14*f[i-2]-6*f[i-3]+f[i-4])*oo12dx2;
  i = n-1;
  d2f[i] = (45*f[i]-154*f[i-1]+214*f[i-2]-156*f[i-3]+61*f[i-4]-10*f[i-5])*oo12dx2;
  return OK;
}

/** 4th order centered stencil first derivative, nonuniform grids */
int D0_nux(double *f, double *x, int n, double *df)
{
  int i, n0 = 0;
  for(i=2; i<=n-2; i++)
    {
      df[i] = 1./3.*(8.*f[1+i] - f[2+i] - 8.*f[i-1] + f[i-2])/(x[2+i]-x[i-2]);
    }
  /* df[n0]   = (-24./17.*f[n0] + 59./34.*f[n0+1] - 4./17.*f[n0+2] - 3./34.*f[n0+3])/(x[n0+1]-x[n0]); */
  /* df[n0+1] = (-1./2.*f[n0] + 1./2.*f[n0+2] )/(x[n0+1]-x[n0]); */
  /* df[n]    = -(-24./17.*f[n] + 59./34.*f[n-1] - 4./17.*f[n-2] - 3./34.*f[n-3])/(x[n]-x[n-1]); */
  /* df[n-1]  = -(-1./2.*f[n] + 1./2.*f[n-2])/(x[n]-x[n-1]); */
  df[n0]   = (-25*f[n0] + 48*f[n0+1] - 36*f[n0+2] +16*f[n0+3]-3*f[n0+4])/(3*(x[n0+4]-x[n0]));
  df[n0+1] = (-3*f[n0]-10*f[n0+1]+18*f[n0+2]-6*f[n0+3]+f[n0+4])/(3*(x[n0+4]-x[n0]));
  df[n]    = (25*f[n] - 48*f[n-1] + 36*f[n-2] - 16*f[n-3]+3*f[n-4])/(3*(x[n]-x[n-4]));
  df[n-1]  = (-f[n-4]+6*f[n-3]-18*f[n-2]+10*f[n-1]+3*f[n])/(3*(x[n]-x[n-4]));  
  return OK;
}

/** This routine sets a 0/1 mask for the multipolar linear index */
void set_multipolar_idx_mask(int *kmask, int n)
{
  int m, k,j;
  for (k = 0; k<n; k++) kmask[k] = 0; /* all off */
  int *idx = par_get_arrayi("output_lm", &m);
  if (m==0) return;
  if (m==1 && idx[0]==-1) return;
  for (k = 0; k<n; k++)
    for (j = 0; j<m; j++)
      if (idx[j] == k) kmask[k] = 1; 
}

/* Alloc/Free data type routines */

/** Waveform (complex) */
void Waveform_alloc (Waveform **wav, int size, char *name)
{
  *wav = (Waveform *) calloc(1, sizeof(Waveform)); 
  if (wav == NULL)
    errorexit("Out of memory");
  (*wav)->real = malloc ( size * sizeof(double) );
  (*wav)->imag = malloc ( size * sizeof(double) );
  (*wav)->time = malloc ( size * sizeof(double) );
  memset( (*wav)->real, 0, size * sizeof(double) );
  memset( (*wav)->imag, 0, size * sizeof(double) );
  memset( (*wav)->time, 0, size * sizeof(double) );
  (*wav)->size = size; 
  strcpy((*wav)->name,name);
}

void Waveform_push (Waveform **wav, int size)
{
  if ((*wav)->real) (*wav)->real = realloc ( (*wav)->real, size * sizeof(double) );
  if ((*wav)->imag) (*wav)->imag = realloc ( (*wav)->imag, size * sizeof(double) );
  if ((*wav)->time) (*wav)->time = realloc ( (*wav)->time, size * sizeof(double) );
  const int n  = (*wav)->size;
  const int dn = size - (*wav)->size;
  if (dn>0) {
    memset( (*wav)->real + n, 0, dn * sizeof(double) );
    memset( (*wav)->imag + n, 0, dn * sizeof(double) );
    memset( (*wav)->time + n, 0, dn * sizeof(double) );
  }
  (*wav)->size = size; 
}

void Waveform_output (Waveform *wav)
{
  int i;
  FILE* fp = fopen(wav->name, "w"); 
  for (i = 0; i < wav->size; i++) {
    fprintf(fp, "%.9e %.12e %.12e\n", wav->time[i], wav->real[i], wav->imag[i]);
  }
  fclose(fp);
}

void Waveform_free (Waveform *wav)
{
  if (wav->real) free(wav->real);
  if (wav->imag) free(wav->imag);
  free(wav);
}

/** Multipolar waveform (complex) */
void Waveform_lm_alloc (Waveform_lm **wav, int size, char *name)
{
  *wav = (Waveform_lm *) calloc(1, sizeof(Waveform_lm)); 
  if (wav == NULL)
    errorexit("Out of memory");
  (*wav)->size = size; 
  set_multipolar_idx_mask((*wav)->kmask, KMAX); 
  (*wav)->time = malloc ( size * sizeof(double) );
  memset((*wav)->time, 0, size*sizeof(double));
  int k;
  for (k=0; k<KMAX; k++) {
    (*wav)->ampli[k] = malloc ( size * sizeof(double) );
    (*wav)->phase[k] = malloc ( size * sizeof(double) );
    memset((*wav)->ampli[k], 0, size*sizeof(double));
    memset((*wav)->phase[k], 0, size*sizeof(double));
    strcpy((*wav)->name[k],&name[k]);
  }
}

void Waveform_lm_push (Waveform_lm **wav, int size)
{
  int k, i;
  const int n  = (*wav)->size;
  const int dn = size - (*wav)->size;
  (*wav)->time = realloc ( (*wav)->time, size * sizeof(double) );
  if ((*wav)->time == NULL) errorexit("Out of memory.");
  if (dn>0) memset( (*wav)->time + n, 0, dn * sizeof(double) );
  for (k=0; k<KMAX; k++) {
    (*wav)->ampli[k] = realloc ( (*wav)->ampli[k], size * sizeof(double) );
    if ((*wav)->ampli[k] == NULL) errorexit("Out of memory.");
    (*wav)->phase[k] = realloc ( (*wav)->phase[k], size * sizeof(double) );
    if ((*wav)->phase[k] == NULL) errorexit("Out of memory.");
    if (dn>0) {
      memset( (*wav)->ampli[k] + n, 0, dn * sizeof(double) );
      memset( (*wav)->phase[k] + n, 0, dn * sizeof(double) );
    } 
  }
  (*wav)->size = size;
}

void Waveform_lm_output (Waveform_lm *wav)
{
  int k,i;
  const int n = wav->size;
  for (k=0; k<KMAX; k++) {
    if (wav->kmask[k]) {
      FILE* fp = fopen(wav->name[k], "w"); 
      for (i = 0; i < n; i++) {
	fprintf(fp, "%.9e %.12e %.12e\n", wav->time[i], wav->ampli[k][i], wav->phase[k][i]);
      }
      fclose(fp);
    }
  }
}

void Waveform_lm_free (Waveform_lm *wav)
{
  int k;
  for (k=0; k<KMAX; k++) {
    if (wav->kmask[k]) {
      if (wav->ampli[k]) free(wav->ampli[k]);
      if (wav->phase[k]) free(wav->phase[k]);
    }
  }
  free(wav);
}

/** Dynamics */
void Dynamics_alloc (Dynamics **dyn, int size, char *name)
{
  (*dyn) = (Dynamics *) calloc(1, sizeof(Dynamics)); 
  if (dyn == NULL)
    errorexit("Out of memory");
  strcpy((*dyn)->name,name);
  (*dyn)->size = size; 
  (*dyn)->time = malloc ( size * sizeof(double) );
  memset((*dyn)->time, 0, size*sizeof(double));
  int v;
  for (v = 0; v < EOB_DYNAMICS_VARS; v++) {
    (*dyn)->data[v] = malloc ( size * sizeof(double) );
    memset((*dyn)->data[v], 0, size*sizeof(double));
  }
}

void Dynamics_push (Dynamics **dyn, int size)
{
  int v,i;
  const int n  = (*dyn)->size;
  const int dn = size - (*dyn)->size;
  (*dyn)->time = realloc ( (*dyn)->time, size * sizeof(double) );
  for (v = 0; v < EOB_DYNAMICS_VARS; v++) {
    (*dyn)->data[v] = realloc ( (*dyn)->data[v], size * sizeof(double) );
    if ((*dyn)->data[v] == NULL) errorexit("Out of memory.");
    if (dn>0) memset( (*dyn)->data[v] + n, 0, dn * sizeof(double) );
  }
  (*dyn)->size = size; 
}

void Dynamics_output (Dynamics *dyn)
{
  int v, i;
  FILE* fp = fopen(dyn->name, "w"); 
  for (i = 0; i < dyn->size; i++) {
    fprintf(fp, "%.9e", dyn->time[i]);
    for (v = 0; v < EOB_DYNAMICS_VARS; v++)
      fprintf(fp, " %.12e", dyn->data[v][i]);
    fprintf(fp, "\n"); 
  }
  fclose(fp);
}

void Dynamics_free (Dynamics *dyn)
{
  if (dyn->time) free(dyn->time);
  int v;
  for (v = 0; v < EOB_DYNAMICS_VARS; v++)
    if (dyn->data[v]) free(dyn->data[v]);
  free(dyn);
}

/** Sync some quick access parameters in dyn with parameter database 
    to be used carefully */
void Dynamics_set_params (Dynamics *dyn)
{
  dyn->store = 0;
  dyn->size  = par_get_i("size");
  dyn->M     = par_get_d("M");
  dyn->nu    = par_get_d("nu");
  dyn->q     = par_get_d("q");
  dyn->X1    = par_get_d("X1");  
  dyn->X2    = par_get_d("X2");
  dyn->chi1  = par_get_d("chi1");
  dyn->chi2  = par_get_d("chi2");
  dyn->S1    = par_get_d("S1");
  dyn->S2    = par_get_d("S2");
  dyn->S     = par_get_d("S");
  dyn->Sstar = par_get_d("Sstar");
  dyn->a1    = par_get_d("a1");
  dyn->a2    = par_get_d("a2"); 
  dyn->aK2   = par_get_d("aK2"); 
  dyn->C_Q1  = par_get_d("C_Q1");
  dyn->C_Q2  = par_get_d("C_Q2"); 
  dyn->cN3LO = par_get_d("cN3LO");
  dyn->rLR   = par_get_d("rLR");
  dyn->rLR_tidal = par_get_d("rLR_tidal");
  dyn->rLSO  = par_get_d("rLSO");
  dyn->kapA2 = par_get_d("kappaAl2");
  dyn->kapA3 = par_get_d("kappaAl3");
  dyn->kapA4 = par_get_d("kappaAl4");
  dyn->kapB2 = par_get_d("kappaBl2");
  dyn->kapB3 = par_get_d("kappaBl3");
  dyn->kapB4 = par_get_d("kappaBl4");
  dyn->kapT2 = par_get_d("kappaTl2");
  dyn->kapT3 = par_get_d("kappaTl3");
  dyn->kapT4 = par_get_d("kappaTl4");
  dyn->khatA2 = par_get_d("khatAl2");
  dyn->khatB2 = par_get_d("khatBl2");
  dyn->rLR_tidal  = par_get_d("rLR_tides");
  dyn->pGSF_tidal = par_get_d("pGSF_tides");
  dyn->Mbhf = par_get_d("BH_final_mass");
  dyn->abhf = par_get_d("BH_final_spin");
  dyn->use_tidal = par_get_i("use_tidal");
  dyn->use_spins = par_get_i("use_spins");
  dyn->dt     = par_get_d("dt");
  dyn->t1     = par_get_d("ode_t1");
  dyn->t_stop = par_get_d("ode_tmax");
}

/** Convert time in sec to dimensionless and mass-rescaled units */
double time_units_factor(double M)
{
  return 1./(M*MSUN_S);
}

double time_units_conversion(double M, double t)
{
  return t/(M*MSUN_S);
}

/** Convert frequency in Hz to dimensionless radius */
double radius0(double M, double fHz)
{
  double x = (M*fHz*MSUN_S*2.*Pi)/2.;
  return cbrt(1/(x*x));
}

/** Make dir */
void system_mkdir(const char *name)
{
  char s[STRLEN];
  sprintf(s,"mkdir -p %s",name);
  system(s); 
}

/** Date and time */
void print_date_time() 
{
  time_t rawtime;
  struct tm * timeinfo;  
  time( &rawtime );
  timeinfo = localtime( &rawtime );
  printf("%s", asctime(timeinfo));
}

/** Errorexit routines */
#undef errorexit
#undef errorexits

void errorexit(char *file, int line, const char *s) 
{
  fprintf(stderr, "Error: %s  ", s);
  fprintf(stderr, "(%s, line %d)\n", file, line);
  fflush(stdout);
  fflush(stderr);
  exit(ERROR);
}

void errorexits(char *file, int line, const char *s, const char *t)
{
  fprintf(stderr, "Error: ");
  fprintf(stderr, s, t);
  fprintf(stderr, "  (%s, line %d)\n", file, line);
  fflush(stdout);
  fflush(stderr);
  exit(ERROR);
}

/* do not write functions beyond this line: errorexit undef/define */

























