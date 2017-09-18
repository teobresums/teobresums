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

#include "TEOBResumS.h"

/** Return symm mass ratio from q */
double q_to_nu(const double q)
{
  double nu = 0;
  if (q>0.)
    nu = q/((q+1.)*(q+1.));
  return nu
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
  /* above constants are deinfed in header */
  return EulerGamma + Log2 + log(m) + 0.5*log(x);
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
  for( int k = ki; k <= kf; k++ ){
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
        errorexits(" wrong l (%d) or m (%d) inside spinspharmY\n", l, m);
    }
    double c = pow(-1.,-s) * sqrt( (2.*l+1.)/(4.*M_PI) );
    double dWigner = c * wigner_d_function(l,m,-s,i);
    *rY = cos((double)(m)*phi) * dWigner;
    *iY = sin((double)(m)*phi) * dWigner;
    return OK;
}

/** 4th order centered stencil first derivative, uniform grids */
int D0(double *f, double dx, int n, double *df)
{
  const double oodx  = 1./dx;
  const double c = 1./12.;
  int i;
  for (i=2; i<n-2; i++)
    {
      d1f[i] = c*(8.*(f[i+1]-f[i-1]) - f[i+2] + f[i-2])*oodx;
    }
  i = 0;
  d1f[i] = c*(-25.*f[i] + 48.*f[i+1] - 36.*f[i+2] + 16.*f[i+3] - 3.*f[i+4])*oodx;
  i = 1;
  d1f[i] = c*(-3.*f[i-1] - 10.*f[i] + 18.*f[i+1] - 6.*f[i+2] + f[i+3])*oodx;
  i = n-2;
  d1f[i] = - c*(-3.*f[i+1] - 10.*f[i] + 18.*f[i-1] - 6.*f[i-2] + f[i-3])*oodx;
  i = n-1;
  d1f[i] = - c*(-25.*f[i] + 48.*f[i-1] - 36.*f[i-2] + 16.*f[i-3] - 3.*f[i-4])*oodx;
  return OK;
}

/** 4th order centered stencil second derivative, uniform grids */
int D2(double *f, double dx, int n, double *d2f)
{
  const double oodx  = 1./dx;
  const double oodx2  = oodx*oodx;
  const double c = 1./12.;
  int i;
  for (i=2; i<n-2; i++)
    {
      d2f[i] = c*(-30*f[i]+16*(f[i+1]+f[i-1])-(f[i+2]+f[i-2]))*oodx2;
    }
  i= 0;
  d2f[i] = c*(45*f[i]-154*f[i+1]+214*f[i+2]-156*f[i+3]+61*f[i+4]-10*f[i+5])*oodx2;
  i= 1;
  d2f[i] = c*(10*f[i-1]-15*f[i]-4*f[i+1]+14*f[i+2]-6*f[i+3]+f[i+4])*oodx2;
  i = n-2;
  d2f[i] = c*(10*f[i+1]-15*f[i]-4*f[i-1]+14*f[i-2]-6*f[i-3]+f[i-4])*oodx2;
  i = n-1;
  d2f[i] = c*(45*f[i]-154*f[i-1]+214*f[i-2]-156*f[i-3]+61*f[i-4]-10*f[i-5])*oodx2;
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
  df[n0] = (-24./17.*f[n0] + 59./34.*f[n0+1] - 4./17.*f[n0+2] - 3./34.*f[n0+3])/(x[n0+1]-x[n0]);
  df[n0+1] = (-1./2.*f[n0] + 1./2.*f[n0+2] )/(x[n0+1]-x[n0]);
  df[n]   = -(-24./17.*f[n] + 59./34.*f[n-1] - 4./17.*f[n-2] - 3./34.*f[n-3])/(x[n]-x[n-1]);
  df[n-1] = -(-1./2.*f[n] + 1./2.*f[n-2])/(x[n]-x[n-1]);
  return OK;
}

/** This routine sets a 0/1 mask for the multipolar linear index */
void set_multipolar_idx_mask(int *kmask, int n)
{
  int m, k,j;
  for (k = 0; k<n; k++)
    kmask = 0; /* all off */
  int *idx = par_get_arrayi("lm", &m);
  for (j = 0; j<m; j++)
    for (k = 0; k<n; k++)
      if (idx[j] == kmask[k]) kmask[k] = 1; 
}

/* Alloc/Free data type routines */

/** Waveform (complex) */
void Waveform_alloc (Waveform **wav, int size, char *name)
{
  *wav = (Waveform *) calloc(1, sizeof(Waveform)); 
  if (*wav == NULL)
    errorexit("Out of memory");
  *wav->real = (double*) malloc ( size * sizeof(double) );
  *wav->imag = (double*) malloc ( size * sizeof(double) );
  //*wav->data = // for complex data
  *wav->size = size; 
  strcpy(name,*wav->name);
}

void Waveform_push (Waveform **wav, int size)
{
  if (*wav->real) *wav->real = (double*) realloc ( size * sizeof(double) );
  if (*wav->imag) *wav->imag = (double*) realloc ( size * sizeof(double) );
  //if (*wav->data) *wav->data = // for complex data
  *wav->size = size; 
}

void Waveform_output (Waveform *wav)
{
  int i;
  FILE* fp = fopen(wav->name, "w"); 
  for (i = 0; i < wav->size; i++) {
    fprintf(fp, "%.9e %.12e %.12e\n", wav->time, wav->real[i], wave->imag[i]);
    //fprintf(fp, "%.9e ??? ???\n", wav->time, wav->data[i]);
  }
  fclose(fp);
}

void Waveform_free (Waveform *wav)
{
  if (wav->real) free(wav->real);
  if (wav->imag) free(wav->imag);
  //if (wav->data) free(wav->data);
  free(wav);
}

/** Multipolar waveform (complex) */
void Waveform_lm_alloc (Waveform_lm **wav, int size, char **name)
{
  *wav = (Waveform_lm *) calloc(1, sizeof(Waveform_lm)); 
  if (*wav == NULL)
    errorexit("Out of memory");
  *wav->size = size; 
  set_multipolar_idx_mask(*wav->kmask, KMAX); 
  int k;
  for (k=0; k<KMAX; k++) {
    if (kmask[k]) {
      //*wav->real[k] = (double*) malloc ( size * sizeof(double) );
      //*wav->imag[k] = (double*) malloc ( size * sizeof(double) );
      *wav->ampil[k] = (double*) malloc ( size * sizeof(double) );
      *wav->phase[k] = (double*) malloc ( size * sizeof(double) );
      strcpy(name[k],*wav->name[k]);
    } else {
      //*wav->real[k] = NULL;
      //*wav->imag[k] = NULL;
      *wav->ampli[k] = NULL;
      *wav->phase[k] = NULL;
    }
  }
}

void Waveform_lm_push (Waveform **wav, int size, int *kmask)
{
  int k;
  for (k=0; k<KMAX; k++) {
    if (wav->kmask[k]) {
      //*wav->real[k] = (double*) realloc ( size * sizeof(double) );
      //*wav->imag[k] = (double*) realloc ( size * sizeof(double) );
      *wav->ampli[k] = (double*) realloc ( size * sizeof(double) );
      *wav->phase[k] = (double*) realloc ( size * sizeof(double) );
    }
  }
}

void Waveform_lm_output (Waveform *wav, int *kmask)
{
  int k,i;
  const int n = wav->size;
  for (k=0; k<KMAX; k++) {
    if (wav->kmask[k]) {
      FILE* fp = fopen(wav->name[k], "w"); 
      for (i = 0; i < n; i++) {
	fprintf(fp, "%.9e %.12e %.12e\n", wav->time, wav->real[k][i], wave->imag[k][i]);
	//fprintf(fp, "%.9e ??? ???\n", wav->time, wav->data[k][i]);
      }
      fclose(fp);
    }
  }
}

void Waveform_lm_free (Waveform_lm *wav)
{
  for (k=0; k<KMAX; k++) {
    if (wav->kmask[k]) {
      if (wav->real[k]) free(wav->real[k]);
      if (wav->imag[k]) free(wav->imag[k]);
      if (wav->ampli[k]) free(wav->ampli[k]);
      if (wav->phase[k]) free(wav->phase[k]);
      strcpy(name[k],wav->name[k]);
    }
  }
  free(wav);
}

/** Dynamics */
void Dynamics_alloc (Dynamics **dyn, int size)
{
  *dyn = (Dynamics *) calloc(1, sizeof(Dynamics)); 
  if (*dyn == NULL)
    errorexit("Out of memory");
  *dyn->size = size; 
  *dyn->time = (double*) malloc ( size * sizeof(double) );
  int v;
  for (v = 0; v < EOB_DYNAMICS_VARS; v++)
    *dyn->data[v] = (double*) malloc ( size * sizeof(double) );
}

void Dynamics_push (Dynamics **dyn, int size)
{
  *dyn->time = (double*) realloc ( size * sizeof(double) );
  int v;
  for (v = 0; v < EOB_DYNAMICS_VARS; v++)
    *dyn->data[v] = (double*) realloc ( size * sizeof(double) );
  *dyn->size = size; 
}

void Dynamics_output (Dynamics *dyn)
{
  int v, i;
  FILE* fp = fopen("dynamics.txt", "w"); 
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

/** Errorexit routines */
#undef errorexit
#undef errorexits

void errorexit(char *file, int line, char *s) 
{
  fprintf(stderr, "Error: %s  ", s);
  fprintf(stderr, "(%s, line %d)\n", file, line);
  fflush(stdout);
  fflush(stderr);
  exit(ERROR);
}

void errorexits(char *file, int line, char *s, char *t)
{
  fprintf(stderr, "Error: ");
  fprintf(stderr, s, t);
  fprintf(stderr, "  (%s, line %d)\n", file, line);
  fflush(stdout);
  fflush(stderr);
  exit(ERROR);
}

/* do not write functions beyond this line: errorexit undef/define */






// ********************************************
// ********************************************
// ********************************************
// ********************************************
// ********************************************
// ********************************************
// ********************************************

/* FOLLOWING ROUTINES SHOULD BE CHECKED/IMPROVED/DELETED USED? */

// ********************************************
// ********************************************
// ********************************************
// ********************************************
// ********************************************
// ********************************************
// ********************************************


static void swap_variables(double *v1, double *v2)
{
    double tmp;
    tmp = *v1;
    *v1 = *v2;
    *v2 = tmp;
}



/** */
double interpolate(double dt,vector<gsl_complex> grid)
{
    double xi, yi;
    double x[]      = {0.,0.,0.,0.,0.,0.,0.};
    double y[]      = {0.,0.,0.,0.,0.,0.,0.};
    double step     = 0.01; //dt/5.;
    double omeg_max = 0.;
    double t_max    = 0.;
    bool peak_flag  = false;
    
    for (int i=0; i<=6; i++)
    {
        x[i] = grid[i].dat[0]; //time
        y[i] = grid[i].dat[1]; //omega
    }
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline    = gsl_spline_alloc (gsl_interp_cspline, 7);
    gsl_spline_init (spline, x, y, 7);
    
    for (xi = x[0]; xi < x[6]; xi += step)
    {
        yi = gsl_spline_eval (spline, xi, acc);
        if (peak_flag==false)
        {
            if (yi<omeg_max)
            {
                peak_flag=true;
            }
            else
            {
                omeg_max=yi;
                t_max=xi;
            }
        }
    }
    
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
    
    return t_max;
}

/** */
vector<double> interp_grid(vector<double> t_vec, vector<double> data, double dt)
{
    int i = 0;
    int t_length = t_vec.size();
    int grid_length = (int)((t_vec.back() - t_vec[0])/dt + 1);
    
    double xi, yi;
    vector<double>     data_g(grid_length);
    vector<double> omg_interp(grid_length);
    
    /** Convert all vectors to an array */
    double* t        = &t_vec[0];
    double* data_arr = &data[0];
    double step      = dt;
    
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline    = gsl_spline_alloc (gsl_interp_cspline, t_length);
    gsl_spline_init(spline, t, data_arr, t_length);
    
    for (xi = t_vec[0]; xi < t_vec.back(); xi += step)
    {
        yi = gsl_spline_eval (spline, xi, acc);
        data_g[i] = yi;
        i++;
    }
    
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
    
    return data_g;
}








// ********************************************
// ********************************************
// ********************************************
// ********************************************
// ********************************************
// ********************************************
// ********************************************

/* FOLLOWING ROUTINES SHOULD GO IN DYNAMICS */

// ********************************************
// ********************************************
// ********************************************
// ********************************************
// ********************************************
// ********************************************
// ********************************************


double time_units_conversion(double M, double dt)
{
  return t/(M*MSUN_S);
}



/** */
double radius0(double M, double f_start)
{
  double x = (M*f_start*MSUN_S*2.*M_PI)/2.;
  return cbrt(1/(x*x));
}

/** The root of this function defines the light ring */
double fLR(double r, void *params)
{
    
    vector<double> metric=Metric(r, params,true);
    
    double A  = metric[0];
    double dA = metric[2]; /**derivative w.r.t. u*/
    double f  = A +(0.5*dA)/r;
    
    return f;
}

/** Find adiabatic light-ring, 
    takes nu as input */
double AdiabLR(void *params)
{
    
    int status;
    int iter = 0, max_iter = 200;
    const gsl_root_fsolver_type *T;
    double rLR;
    double x_lo = 0.1, x_hi = 15.;
    
    gsl_root_fsolver *s;
    gsl_function F;
    F.function = &fLR;
    F.params = params;
    T = gsl_root_fsolver_bisection;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
    
    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        rLR    = gsl_root_fsolver_root (s);
        x_lo   = gsl_root_fsolver_x_lower (s);
        x_hi   = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi, 0, 0.000000000000001);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (s);
    
    return rLR;
}



// ********************************************
// ********************************************
// ********************************************
// ********************************************
// ********************************************
// ********************************************
// ********************************************

/* FOLLOWING ROUTINES SHOULD GO IN WAVEFORM */

// ********************************************
// ********************************************
// ********************************************
// ********************************************
// ********************************************
// ********************************************
// ********************************************



/** Alternative implementation of the phase of the tail factor */
vector<gsl_complex> speedyTail(const double Omega, const double Hreal, const double bphys, const int L[], const int M[])
{
    int kmax = 35;
    
    double ratio_rad;
    double ratio_ang;
    double tlm_rad;
    double tlm_ang;
    
    double x;
    double x2;
    double x3;
    double x4;
    double x5;
    
    vector<double> num_ang(kmax);
    vector<gsl_complex> tlm(kmax);
    
    /** Fit coefficients*/
    /*
     double a1[] =
     {
     0.3060234167900050, 0.3060234167900050, 0.0424759238428813, 0.0424489015884926, 0.0424717446800903, 0.0215953972500844, 0.0215873812155663, 0.0215776183122621, 0.0216017621863542, 0.0128123696874894, 0.0128097056242375, 0.0128038943888768, 0.0128025242617949, 0.0128202485907368, 0.0083762045692408, 0.0083751913886140, 0.0083724067460769, 0.0083694435961860, 0.0083710364141552, 0.0083834483913443, 0.0058540393221396, 0.0058536069384738, 0.0058522594457692, 0.0058502436535615, 0.0058491157293566, 0.0058514875071582, 0.0058602498033381, 0.0042956812356573, 0.0042954784390887, 0.0042947951664056, 0.0042935886137697, 0.0042923691461384, 0.0042922256848799, 0.0042945927126022, 0.0043009106861259};
     
     double a2[] =
     {
     -0.0248640052699995, -0.0248640052699995, 0.0000597134489198, 0.0002551406918111, 0.0001741036904709, 0.0000124649041611, 0.0000685496215625, 0.0001131160409390, 0.0000419907542591, 0.0000035218982282, 0.0000219211271097, 0.0000473186962874, 0.0000524142634057, 0.0000106823372552, 0.0000012237574387, 0.0000081742188269, 0.0000201940563214, 0.0000295722761753, 0.0000260539631956, 0.0000018994753518, 0.0000004932942990, 0.0000034477210351, 0.0000092294406360, 0.0000155143073237, 0.0000183386499818, 0.0000137922469695, -0.0000007075155453, 0.0000002223410995, 0.0000016045317657, 0.0000045260028113, 0.0000082655700107, 0.0000112393599417, 0.0000115758243113, 0.0000076838709956, -0.0000014020591745
     };
     
     double a3[] = {
     0.1947531537583291, 0.1947531537583291, -0.0042911639711832, -0.0047431560217121, -0.0046577314472149, -0.0013089557502947, -0.0014343968205390, -0.0014978542575474, -0.0014329302934532, -0.0005167994164556, -0.0005573939123058, -0.0005921030407223, -0.0005978284714483, -0.0005673965369076, -0.0002409269302708, -0.0002561516055118, -0.0002723768586352, -0.0002815958312453, -0.0002792078156272, -0.0002646630240693, -0.0001261183503407, -0.0001325622938779, -0.0001403198638518, -0.0001464084186977, -0.0001485971591029, -0.0001459023931717, -0.0001384829633836, -0.0000719062974278, -0.0000749128468013, -0.0000788187384314, -0.0000824202283094, -0.0000846673495936, -0.0000849054394951, -0.0000829269749240, -0.0000788883333858
     };
     
     double a4[] = {
     -0.0254800838696074, -0.0254800838696074, 0.0006914295465364, 0.0010322294603561, 0.0010057563135650, 0.0001394203795507, 0.0002309706405978, 0.0002596611624417, 0.0002409588083156, 0.0000386949167221, 0.0000679154947896, 0.0000830199015202, 0.0000850120755064, 0.0000780125513602, 0.0000133034384660, 0.0000241813441339, 0.0000311573885555, 0.0000340233089866, 0.0000335167900637, 0.0000307571022927, 0.0000053305073331, 0.0000099143129290, 0.0000132296989826, 0.0000150959309402, 0.0000156304390748, 0.0000151274875147, 0.0000139320508803, 0.0000023959090314, 0.0000045285807761, 0.0000061918979830, 0.0000072894226381, 0.0000078251853305, 0.0000078772667984, 0.0000075606242809, 0.0000069956215270
     };
     
     double a5[] = {
     -0.0000140714300659, -0.0000140714300659, 0.0006914295465364, 0.0010322294603561, 0.0010057563135650, 0.0001394203795507, 0.0002309706405978, 0.0002596611624417, 0.0002409588083156, 0.0000386949167221, 0.0000679154947896, 0.0000830199015202, 0.0000850120755064, 0.0000780125513602, 0.0000133034384660, 0.0000241813441339, 0.0000311573885555, 0.0000340233089866, 0.0000335167900637, 0.0000307571022927, 0.0000053305073331, 0.0000099143129290, 0.0000132296989826, 0.0000150959309402, 0.0000156304390748, 0.0000151274875147, 0.0000139320508803, 0.0000023959090314, 0.0000045285807761, 0.0000061918979830, 0.0000072894226381, 0.0000078251853305, 0.0000078772667984, 0.0000075606242809, 0.0000069956215270
     };
     */
    
    double b1[] =
    {
        0.1113090643348557, 0.1112593821157397, 0.0424759238428813, 0.0424489015884926, 0.0424717446800903, 0.0215953972500844, 0.0215873812155663, 0.0215776183122621, 0.0216017621863542, 0.0128123696874894, 0.0128097056242375, 0.0128038943888768, 0.0128025242617949, 0.0128202485907368, 0.0083762045692408, 0.0083751913886140, 0.0083724067460769, 0.0083694435961860, 0.0083710364141552, 0.0083834483913443, 0.0058540393221396, 0.0058536069384738, 0.0058522594457692, 0.0058502436535615, 0.0058491157293566, 0.0058514875071582, 0.0058602498033381, 0.0042956812356573, 0.0042954784390887, 0.0042947951664056, 0.0042935886137697, 0.0042923691461384, 0.0042922256848799, 0.0042945927126022, 0.0043009106861259};
    
    double b2[] =
    {
        0.0004643273300862, 0.0009375605440004, 0.0000597134489198, 0.0002551406918111, 0.0001741036904709, 0.0000124649041611, 0.0000685496215625, 0.0001131160409390, 0.0000419907542591, 0.0000035218982282, 0.0000219211271097, 0.0000473186962874, 0.0000524142634057, 0.0000106823372552, 0.0000012237574387, 0.0000081742188269, 0.0000201940563214, 0.0000295722761753, 0.0000260539631956, 0.0000018994753518, 0.0000004932942990, 0.0000034477210351, 0.0000092294406360, 0.0000155143073237, 0.0000183386499818, 0.0000137922469695, -0.0000007075155453, 0.0000002223410995, 0.0000016045317657, 0.0000045260028113, 0.0000082655700107, 0.0000112393599417, 0.0000115758243113, 0.0000076838709956, -0.0000014020591745
    };
    
    double b3[] =
    {
        -0.0221835462237291, -0.0235386333304348, -0.0042911639711832, -0.0047431560217121, -0.0046577314472149, -0.0013089557502947, -0.0014343968205390, -0.0014978542575474, -0.0014329302934532, -0.0005167994164556, -0.0005573939123058, -0.0005921030407223, -0.0005978284714483, -0.0005673965369076, -0.0002409269302708, -0.0002561516055118, -0.0002723768586352, -0.0002815958312453, -0.0002792078156272, -0.0002646630240693, -0.0001261183503407, -0.0001325622938779, -0.0001403198638518, -0.0001464084186977, -0.0001485971591029, -0.0001459023931717, -0.0001384829633836, -0.0000719062974278, -0.0000749128468013, -0.0000788187384314, -0.0000824202283094, -0.0000846673495936, -0.0000849054394951, -0.0000829269749240, -0.0000788883333858
    };
    
    double b4[] =
    {
        0.0058366730167965, 0.0070452306758401, 0.0006914295465364, 0.0010322294603561, 0.0010057563135650, 0.0001394203795507, 0.0002309706405978, 0.0002596611624417, 0.0002409588083156, 0.0000386949167221, 0.0000679154947896, 0.0000830199015202, 0.0000850120755064, 0.0000780125513602, 0.0000133034384660, 0.0000241813441339, 0.0000311573885555, 0.0000340233089866, 0.0000335167900637, 0.0000307571022927, 0.0000053305073331, 0.0000099143129290, 0.0000132296989826, 0.0000150959309402, 0.0000156304390748, 0.0000151274875147, 0.0000139320508803, 0.0000023959090314, 0.0000045285807761, 0.0000061918979830, 0.0000072894226381, 0.0000078251853305, 0.0000078772667984, 0.0000075606242809, 0.0000069956215270
    };
        
    const vector<double> Tlm_real = Tlm(Omega*Hreal);
    
    /** Pre-computed psi */
    double psi[] = {0.9227843350984671394, 0.9227843350984671394,
        1.256117668431800473, 1.256117668431800473, 1.256117668431800473,
        1.506117668431800473, 1.506117668431800473, 1.506117668431800473, 1.506117668431800473,
        1.706117668431800473, 1.706117668431800473, 1.706117668431800473, 1.706117668431800473, 1.706117668431800473,
        1.872784335098467139, 1.872784335098467139, 1.872784335098467139, 1.872784335098467139, 1.872784335098467139, 1.872784335098467139,
        2.015641477955609997, 2.015641477955609997, 2.015641477955609997, 2.015641477955609997, 2.015641477955609997, 2.015641477955609997, 2.015641477955609997,
        2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997};
    
    double k;
    for (int i=kmax; i--;)
    {
        k  = M[i] * Omega;
        x  = k * Hreal; //hathatk
        x2 = x * x;
        x3 = x2 * x;
        x4 = x3 * x;
        x5 = x4 * x;
        // Why no linear term in the current implementation?
        
        //num_ang[i] = (1. + a1[i]*x2 + a2[i]*x3)/(1. + a3[i]*x2 + a4[i]*x3 + a5[i]*x4); //A.N. fit
        num_ang[i]    = 1. + b1[i]*x2 + b2[i]*x3 + b3[i]*x4 + b4[i]*x5; //P.F. fit
        
        ratio_rad     = Tlm_real[i];
        ratio_ang     = - 2. * psi[i] * x * num_ang[i];
        tlm_rad       = ratio_rad;
        tlm_ang       = ratio_ang + 2. * x * log(2. * k * bphys);
        
        tlm[i].dat[0] = tlm_rad;
        tlm[i].dat[1] = tlm_ang;
    }
    
    return tlm;
}

















