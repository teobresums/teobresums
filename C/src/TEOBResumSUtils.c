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

/** Eulerlog function (constants are defined in header) */
static const double Logm[] = {0.,Log1,Log2,Log3,Log4,Log5,Log6,Log7};
double Eulerlog(const double x,const int m)
{
  double logm = 0.;
  if ((m>0) & (m<8)) logm = Logm[m];
  else logm = log((double)m);
  return EulerGamma_Log2 + logm + 0.5*log(x);
}

/** Spline interpolation with GSL routines */
void interp_spline(double *t, double *y, int n, double *ti, int ni, double *yi)
{
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n);
  gsl_spline_init (spline, t, y, n);    
  for (int k = 0; k < ni; k++) {
    yi[k] = gsl_spline_eval (spline, ti[k], acc);
  }
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
}

/** Spline interpolation with GSL routines 
    Check limits and set to zero points outside interp range */
void interp_spline_checklim(double *t, double *y, int n, double *ti, int ni, double *yi)
{
  const double ta = t[0];
  const double tb = t[n-1];
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n);
  gsl_spline_init (spline, t, y, n);    
  for (int k = 0; k < ni; k++) {
    if ((ti[k]<=ta)||(ti[k]>=tb)) yi[k] = 0.;
    else yi[k] = gsl_spline_eval (spline, ti[k], acc);
  }
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
}


/* An OpenMP version. We keep two versions because we might want to introduce the 
   thread-parallelism at different levels */
/*
  Notes for SIMD
  * SIMD Functions need to be declared, see e.g.
    P.25 https://info.ornl.gov/sites/publications/files/Pub69214.pdf
  * How does it work with GSL? Is wrapper enough?
  * 'uniform' means 'arg does not change'
  * 'linear' means 'arg will be increased between each successive call fo the fun'
  GSL pointers:
  https://github.com/ampl/gsl/blob/master/interpolation/spline.c#L118
  https://github.com/ampl/gsl/blob/master/interpolation/interp.c#L140
*/
#pragma omp declare simd uniform(x) linear(i: 1)
double gsl_spline_eval_omp(const gsl_spline * spline, double *x, gsl_interp_accel * acc, int i)
{
  return gsl_spline_eval(spline, x[i], acc);
}

void interp_spline_omp(double *t, double *y, int n, double *ti, int ni, double *yi)
{
#ifdef _OPENMP
  if (USETIMERS) openmp_timer_start("interp_spline");
#endif
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n);
  gsl_spline_init (spline, t, y, n);    
#pragma omp simd 
  for (int k = 0; k < ni; k++) 
    yi[k] = gsl_spline_eval_omp (spline, ti, acc, k);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
#ifdef _OPENMP  
  if (USETIMERS) openmp_timer_stop("interp_spline");
#endif
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

/** Find max location by poynomial interpolation around x0 (uniform grid) */
double find_max (const int n, double dx, double x0, double *f, double *fmax)
{
  const int i = (n-1)/2; /* To centre the grid for n = N, e.g., n=7, i-3 = 0. */
  double xmax = x0;
  double d1f = 0., d2f = 0.;
  if (n==3) {
    d1f = 0.5*(f[i+1]-f[i-1]);
    d2f = (f[i-1]-2*f[i]+f[i+1]); 
  } else if (n==5) {
    d1f = (8.*(f[i+1]-f[i-1]) - f[i+2] + f[i-2]); 
    d2f = (-30*f[i]+16*(f[i+1]+f[i-1])-(f[i+2]+f[i-2]));
  } else if (n==7) { 
    d1f = ( 45.0*(f[i+1]-f[i-1]) - 9.0*(f[i+2] - f[i-2]) + f[i+3] - f[i-3] )/(60.0);
    d2f = ( -490.0*f[i]+270.0*(f[i+1]+f[i-1])-27.0*(f[i+2]+f[i-2])+2.0*(f[i+3]+f[i-3]) )/(180.0); 
  }
  else errorexit("Implemented only n = 3,5,7");    
  
  if (d2f != 0.) {    
    xmax -= dx*d1f/d2f;
  }

  /* Eval function */
  if (fmax!=NULL) {
    if (n==3) {
      *fmax = ( -((dx + x0 - xmax)*(-x0 + xmax)*f[-1 + i])
		+ (dx - x0 + xmax)*(2*(dx + x0 - xmax)*f[i] + (-x0 + xmax)*f[1 + i]) );
      *fmax /= (2.*SQ(dx));
    } else if (n==5) {
      *fmax = ((dx + x0 - xmax)*(2*dx + x0 - xmax)*(-x0 + xmax)*(dx - x0 + xmax)*f[-2 + i]
	       + (2*dx - x0 + xmax)*(-4*(dx + x0 - xmax)*(2*dx + x0 - xmax)*(-x0 + xmax)*f[-1 + i]+(dx - x0 + xmax)*(6*(dx + x0 - xmax)*(2*dx + x0 - xmax)*f[i] + (-x0 + xmax)*(4*(2*dx + x0 - xmax)*f[1 + i]- (dx + x0 - xmax)*f[2 + i]))));
      *fmax /= (24.*pow(dx,4));
      
    } else errorexit("Implemented only n = 3,5");    
  }
  
  return xmax;
}

/** Find max location around x0 using 5 points (non-uniform grid) */
double find_max_grid (double *x, double *f)
{
  const int i = 2; /* centre the grid around third point */
  
  double dx[7];
  for (int j = 0; j < 4; j++){
    dx[j] = x[j+1] - x[j];
  }
  
  double d1f = 0., d2f = 0.;
  d1f = 0.5*(f[i+1]-f[i])/dx[i] + 0.5*(f[i]-f[i-1])/dx[i-1];
  d2f = 0.25*((f[i+2]-f[i+1])/dx[i+1] - (f[i]-f[i-1])/dx[i-1])/dx[i]
    + 0.25*((f[i+1]-f[i])/dx[i] - (f[i-1]-f[i-2])/dx[i-2])/dx[i-1];  
      
  double xmax = x[i] - d1f/d2f;

  return xmax;
}

/** Factorial */
static const double f35[] = {1.,
			     1.,
			     2.,
			     6.,
			     24.,
			     120.,
			     720.,
			     5040.,
			     40320.,
			     362880.,
			     3628800., 
			     39916800.,
			     479001600.,
			     6227020800.,
			     87178291200.,
			     1307674368000.,
			     20922789888000.,
			     355687428096000.,
			     6402373705728000.,
			     121645100408832000.,
			     2432902008176640000.,
			     51090942171709440000.,
			     1124000727777607680000.,
			     25852016738884976640000.,
			     620448401733239439360000.,
			     15511210043330985984000000.,
			     403291461126605635584000000.,
			     10888869450418352160768000000.,
			     304888344611713860501504000000.,
			     8841761993739701954543616000000.,
			     265252859812191058636308480000000.,
			     8222838654177922817725562880000000.,
			     263130836933693530167218012160000000.,
			     8683317618811886495518194401280000000.,
			     295232799039604140847618609643520000000.,
			     10333147966386144929666651337523200000000.};
double fact(int n)
{
  if (n < 0){
    errorexit(" computing a negative factorial.\n");
  } else if (n <= 35){
    return f35[n];
  } else {
    return n*fact(n-1);
  }
}

/** Wigner d-function */
double wigner_d_function(int l, int m, int s, double i)
{
  const double costheta = cos(i*0.5);
  const double sintheta = sin(i*0.5);
  const double norm = sqrt( (fact(l+m) * fact(l-m) * fact(l+s) * fact(l-s)) );
  const int ki = MAX( 0  , m-s );
  const int kf = MIN( l+m, l-s );
  double dWig = 0.;  
  double div;
  for (int k = ki; k <= kf; k++ ) {
    div = 1.0/( fact(k) * fact(l+m-k) * fact(l-s-k) * fact(s-m+k) );
    dWig += div*( pow(-1.,k) * pow(costheta,2*l+m-s-2*k) * pow(sintheta,2*k+s-m) );
  }
  return (norm * dWig);
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
void compute_hpc(Waveform_lm *hlm, double nu, double M, double distance, double amplitude_prefactor, double phi, double iota, Waveform *hpc)
{  
#ifdef _OPENMP
  if (USETIMERS) openmp_timer_start("compute_hpc");
#endif
#pragma omp parallel 
  {
    double Y_real[KMAX], Y_imag[KMAX];
    static const int mneg = 1; /* m>0 modes only, add m<0 modes afterwards */
    double Y_real_mneg[KMAX], Y_imag_mneg[KMAX];
    double Aki, cosPhi, sinPhi;
    double sumr, sumi;
    int activemode[KMAX];
    double Msun = M;
    set_multipolar_idx_mask (activemode, KMAX, EOBPars->use_mode_lm, EOBPars->use_mode_lm_size, 1);
    if (!(EOBPars->use_geometric_units)) Msun = M/MSUN_S;
#if (DEBUG)
    printf("h+,x: nu = %e M = %e D = %e Mpc phi = %e iota = %e prefactor = %e\n",
	   nu,Msun,distance,phi,iota,amplitude_prefactor);
#endif
    
    /* Precompute Ylm */
    for (int k = 0; k < KMAX; k++ ) {
      if (!activemode[k]) continue;
      spinsphericalharm(&Y_real[k], &Y_imag[k], -2, LINDEX[k], MINDEX[k], phi, iota);
      /* add m<0 modes */
      if ( (mneg) && (MINDEX[k]!=0) )
	spinsphericalharm(&Y_real_mneg[k], &Y_imag_mneg[k], -2, LINDEX[k], -MINDEX[k], phi, iota); 
      }
      
    /* Sum up  hlm * Ylm 
     * Note because EOB code defines phase>0, 
     * but the convention is h_lm = A_lm Exp[-I phi_lm] we have
     * h_{l,m>=0} = A_lm ( cos(phi_lm) - I*sin(phi_lm) ) for m>0 and
     * h_{l,m<0}  = (-)^l A_l|m| ( cos(phi_l|m|) + I*sin(phi_l|m|) ) for m<0 below
     * We now agree with, e.g., LALSimSphHarmMode.c: 64-74
     */
#pragma omp for
    for (int i = 0; i < hlm->size; i++) {
        hpc->time[i] = hlm->time[i]*M; 
        sumr = sumi = 0.;
	
        /* Loop over modes */
        for (int k = 0; k < KMAX; k++ ) {
	  if (!activemode[k]) continue;
	  Aki  = amplitude_prefactor * hlm->ampli[k][i];
	  cosPhi = cos( hlm->phase[k][i] );
	  sinPhi = sin( hlm->phase[k][i] );
	  sumr += Aki*(cosPhi*Y_real[k] + sinPhi*Y_imag[k]);
	  sumi += Aki*(cosPhi*Y_imag[k] - sinPhi*Y_real[k]); 
          
	  /* add m<0 modes */
	  if ( (mneg) && (MINDEX[k]!=0) ) { 
	    /* H_{l-m} = (-)^l H^{*}_{lm} */
	    if (LINDEX[k] % 2) {
	      sumr -= Aki*(cosPhi*Y_real_mneg[k] - sinPhi*Y_imag_mneg[k]); 
	      sumi -= Aki*(cosPhi*Y_imag_mneg[k] + sinPhi*Y_real_mneg[k]); 
	    }
	    else { 
	      sumr += Aki*(cosPhi*Y_real_mneg[k] - sinPhi*Y_imag_mneg[k]);
	      sumi += Aki*(cosPhi*Y_imag_mneg[k] + sinPhi*Y_real_mneg[k]); 
	    }
	  }    
        }
        /* h = h+ - i hx */
        hpc->real[i] = sumr;
        hpc->imag[i] = -sumi;
    }
  }
#ifdef _OPENMP
  if (USETIMERS) openmp_timer_stop("compute_hpc");
#endif
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

/** 2nd order centered stencil first derivative, nonuniform grids */
int D0_x_2(double *f, double *x, int n, double *df)
{
  int i;
  for(i=1; i<n-1; i++) {
    df[i] = (f[i+1]-f[i-1])/(x[i+1]-x[i-1]);
  }
  i = 0;
  df[i] = (f[i]-f[i+1])/(x[i]-x[i+1]);
  i = n-1;
  df[i] = (f[i-1]-f[i])/(x[i-1]-x[i]);
  return OK;
}

/** 4th order first derivative, nonuniform grid */
/** FIXME: can be further optimized! */
int D0_x_4(double *f, double *x, int n, double *df)
{
  double *ix = x;
  double *iy = f;

  /* left boundary */
  df[0] = d4(ix, iy, 0);
  df[1] = d4(ix, iy, 1);

  /* central stencil */
  for(int k = 2; k < n-2; k++){
    df[k] = d4(ix, iy, 2);
    ix++; iy++;
  }

  /*right boundary*/
  df[n-2] = d4(ix, iy, 3);
  df[n-1] = d4(ix, iy, 4);
  
  return OK;
}

double d4(double *x, double *y, int i)
{
  double dy_i =0.; 
  double ld_j;

  for (int j = 0; j < 5; j++){
    ld_j = l_deriv(x, i, j);
    dy_i = *(y+j) *ld_j + dy_i;
  }
  return dy_i;
}

double l_deriv(double *x,int i,int j)
{
  double xi = *(x +i);
  double xj = *(x +j);
  double ld_ji;
  /*compute denominator*/
  double D = 1.;
  for (int a = 0; a<5; a++){
    if (a != j) D = D*(xj - *(x +a));
  }
  /*compute N and N/D */
  double N = 1.;
  if (i != j){                            //case 1: i != j
    for (int b = 0; b < 5 ; b++){
      if (b != j && b!= i) 
        N = N*(xi - *(x +b));
    }
    ld_ji = N/D;
  } else {                                //case 2: 1 ==j
    double tmp = 0.;
    for (int c = 0; c < 5; c++){
      double N = 1;
      if (c != j){
        for (int b = 0; b < 5; b++){
          if (b !=j && b != c) N = N*(xj - *(x + b));
        }
      tmp = tmp + N;  
      }
    }
    ld_ji = tmp/D;
  }

  return ld_ji;
}

/** Trapezoidal rule */
double cumtrapz(double *f, double *x, const int n, double *sum)
{
  sum[0] = 0.;
  if (x != NULL) {
    for (int i=1; i < n; i++) 
      sum[i] = sum[i-1] + 0.5*(x[i] - x[i-1])*(f[i] + f[i-1]);
  } else {
    for (int i=1; i < n; i++) 
      sum[i] = sum[i-1] + 0.5*(f[i] + f[i-1]);
  }
  return sum[n-1];
}

/** Third-order polynomial integration */
double cumint3(double *f, double *x, const int n, double *sum)
{
  double xe[n+2], fe[n+2];
  double *x0,*x1,*x2,*x3;
  double *f0,*f1,*f2,*f3;
  double a,b,c,d,e,h,g,z;
  const double oo12 = 0.08333333333333333;
  
  for (int i=1; i < n+1; i++) {
    xe[i] = x[i-1];
    fe[i] = f[i-1];
  }
  xe[0] = x[3];
  xe[n+1] = x[n-4];
  fe[0] = f[3];
  fe[n+1] = f[n-4];

  x0 = &xe[0];
  x1 = &xe[1];
  x2 = &xe[2];  
  x3 = &xe[3];

  f0 = &fe[0];
  f1 = &fe[1];
  f2 = &fe[2];  
  f3 = &fe[3];
  
  sum[0] = 0.;
  for (int i=0; i < n-1; i++) {
    a = x1[i]-x0[i];
    b = x2[i]-x1[i];
    c = x3[i]-x2[i];
    d = f1[i]-f0[i];
    e = f2[i]-f1[i];
    h = f3[i]-f2[i];
    g = 0.5*(f1[i]+f2[i]);
    z = b*g + oo12*b*b*(c*b*(2*c+b)*(c+b)*d-a*c*(c-a)*(2*c+2*a+3*b)*e-a*b*(2*a+b)*(a+b)*h)/(a*c*(a+b)*(c+b)*(c+a+b));
    sum[i+1] = sum[i] + z;
  }

  return sum[n-1];
}

/* Simple unwrap for phase angles */ // Do NOT mess with it
void unwrap(double *p, const int size)
{
 if (size < 1) return;
  int j;
  int fact = 0;  // For making the initial negative phase angles positive
  double curr, prev;
  double corr = 0.0;
  double dphi = 0.0;

  prev = p[0];
  if( p[0] < 0 ) fact = 1;
  if( p[1] < p[0] ) 
    dphi = TwoPi;
  //printf("phi[0]= %.3f\tfact=%d\n", p[0], fact);
  for (j = 1; j < size; j++){
    //if(j<10)printf("j=%d\tp[j]= %.3f\t", j, p[j]);
    p[j] += fact*TwoPi;
    //if(j<10)printf("p[j]+/-fact= %.3f\t", p[j]);
    curr = p[j];
    if( curr < prev ) 
      dphi = TwoPi;
    corr += dphi;
    p[j] += corr - fact*TwoPi;      
    //if(j<10)printf("p[j]corr= %.3f\n", p[j]);   
    prev = curr;
    dphi = 0.0;
  }
}

/* Slightly modified unwrap function for HM phases */
void unwrap_HM(double *p, const int size)
{
  if (size < 1) return;
  int j;
  double dphi, corr, curr, prev;
  dphi = 0.;
  corr = 0.;
  
  prev = *p;  
  for (j = 1; j < size; j++){
    p++;
    curr = *p;
    
    if(curr < prev - Pi) 
      dphi = TwoPi;
    if(curr > prev + TwoPi)
      dphi = -TwoPi;

    corr += dphi;
    *p += corr;    

    prev = curr;
    dphi = 0.0;    
  } 
  
}

#define dbg_unwrap_proxy (0)  /* stops after routine, use: ./TEOBResumS.x test.par > out */ 
/* Unwrap unsign number of cycles from reference phase as proxy */
void unwrap_proxy(double *p, double *r, const int size, const int shift0)
{
  if (size < 1) return;

  int *n = (int *) calloc(size, sizeof(int));
  if (n == NULL)
    errorexit("Out of memory");

  const double ooTwoPi = 1.0/(TwoPi);
  const double r0 = r[0];

  /* no cycles from reference phase, r>0 */
  for (int i = 0; i < size; i++)
    n[i] = floor ( (r[i] - r0) * ooTwoPi );

  if (shift0) {
    /* shift phase : p(0) = r(0) */
    const double dp = (r0 - p[0]);
    for (int i = 0; i < size; i++) 
      p[i] += dp; 
  }
  
  /* unwrap based on no cycles */
  const double p0 = p[0];  
  int np = 0;
  if (dbg_unwrap_proxy) printf("#i:0 r:1 n:2 p:3 (r-p_unwrap):4 np:5 (n-np):6 p_unwrap_correct:7 (r-p_unwrap_corrected):8\n");
  for (int i = 0; i < size; i++) {
    if (dbg_unwrap_proxy) printf("%06d %+.16e %06d %+.16e ",i,r[i],n[i],p[i]);
    p[i] += n[i]*TwoPi;
    if (dbg_unwrap_proxy) printf("%+.16e", r[i]-p[i]);
    /* correct cases p = - Pi + epsilon */
    np = floor ( ( p[i] - p0 )*ooTwoPi );
    p[i] += (n[i]-np)*TwoPi;
    if (dbg_unwrap_proxy) printf(" %06d %06d %+.16e %+.16e\n",np, n[i]-np, p[i], r[i]-p[i]);
  }

  free(n);
  if (dbg_unwrap_proxy) DBGSTOP;
}

/** This routine sets a 0/1 mask for the multipolar linear index */
void set_multipolar_idx_mask_old(int *kmask, int n)
{
  int m, k,j;
  for (k = 0; k<n; k++) kmask[k] = 0; /* all off */
  int *idx = par_get_arrayi("use_mode_lm", &m);
  if (m==0) return;
  if (m==1 && idx[0]==-1) return;
  for (k = 0; k<n; k++)
    for (j = 0; j<m; j++)
      if (idx[j] == k) kmask[k] = 1;
  free(idx);
}

/** This routine sets a 0/1 mask for the multipolar linear index 
    work for any parameter and can specify default all on/off */
void set_multipolar_idx_mask (int *kmask, int n, const int *idx, int m, int on)
{
  int k,j;
  for (k = 0; k<n; k++) kmask[k] = 0; /* all off */
  if ( (m<=0) || (m==1 && idx[0]==-1) ) {
    if(on) for (k = 0; k<n; k++) kmask[k] = 1; /* all on */
    return;
  } 
  for (k = 0; k<n; k++)
    for (j = 0; j<m; j++)
      if (idx[j] == k) kmask[k] = 1; 
}

/** Compute size of a uniform grid t0:dt:tf */
int get_uniform_size(const double tN, const double t0, const double dt)
{
  return ((int)((tN - t0)/dt + 1)); 
}

/* Alloc/Free data type routines */

/** Waveform (complex) */
void Waveform_alloc (Waveform **wav, const int size, const char *name)
{
  *wav = (Waveform *) calloc(1, sizeof(Waveform)); 
  if (wav == NULL)
    errorexit("Out of memory");
  (*wav)->real =  malloc ( size * sizeof(double) );
  (*wav)->imag =  malloc ( size * sizeof(double) );
  (*wav)->ampli = malloc ( size * sizeof(double) );
  (*wav)->phase = malloc ( size * sizeof(double) );
  (*wav)->time =  malloc ( size * sizeof(double) );
  memset( (*wav)->real, 0, size * sizeof(double) );
  memset( (*wav)->imag, 0, size * sizeof(double) );
  memset((*wav)->ampli, 0, size * sizeof(double) );
  memset((*wav)->phase, 0, size * sizeof(double) );
  (*wav)->size = size; 
  strcpy((*wav)->name,name);
}

void Waveform_push (Waveform **wav, int size)
{
  if ((*wav)->real)  (*wav)->real  = realloc ( (*wav)->real,  size * sizeof(double) );
  if ((*wav)->imag)  (*wav)->imag  = realloc ( (*wav)->imag,  size * sizeof(double) );
  if ((*wav)->ampli) (*wav)->ampli = realloc ( (*wav)->ampli, size * sizeof(double) );
  if ((*wav)->phase) (*wav)->phase = realloc ( (*wav)->phase, size * sizeof(double) );
  if ((*wav)->time)  (*wav)->time  = realloc ( (*wav)->time,  size * sizeof(double) );
  const int n  = (*wav)->size;
  const int dn = size - (*wav)->size;
  (*wav)->size = size; 
}

/* Compute real/imag <-> amplitude/phase */
void Waveform_rmap (Waveform *h, const int mode, const int unw)
{
  const int size = h->size;
  if (mode) {
    /** (Re, Im) -> (Amplitude, phase) */
#pragma omp simd
    for (int i = 0; i < size; i++)
      h->ampli[i] = sqrt( SQ(h->real[i]) + SQ(h->imag[i]) );
#pragma omp simd
    for (int i = 0; i < size; i++)
      h->phase[i] = Pi - atan2(h->imag[i], h->real[i]); /* exp(- i phi) => Pi  */
    if (unw) unwrap(h->phase, h->size); 
  } else {
    /** (Amplitude, phase) -> (Re, Im) */
    /* h =  A exp( -i phi) */
#pragma omp simd
    for (int i = 0; i < size; i++) 
      h->real[i] = + h->ampli[i] * cos(h->phase[i]);  
#pragma omp simd
    for (int i = 0; i < size; i++)
      h->imag[i] = - h->ampli[i] * sin(h->phase[i]); 
  }
}

void Waveform_interp (Waveform *h, const int size, const double t0, const double dt, const char *name)
{
  /* Alloc and init aux memory */  
  Waveform *h_aux;
  const int oldsize = h->size;
  Waveform_alloc(&h_aux, oldsize, "");
  memcpy(h_aux->time, h->time, oldsize * sizeof(double));
  memcpy(h_aux->real, h->real, oldsize * sizeof(double));
  memcpy(h_aux->imag, h->imag, oldsize * sizeof(double));
 
  /* Realloc arrays */
  h->size = size;
  if (strcmp(name, "")) strcpy(h->name, name);
  if (h->time) free(h->time);
  if (h->real) free(h->real);
  if (h->imag) free(h->imag);
  h->time = malloc ( size * sizeof(double) );
  h->real = malloc ( size * sizeof(double) );
  h->imag = malloc ( size * sizeof(double) );

  /* Fill new time array */
#pragma omp simd
  for (int i = 0; i < size; i++)
    h->time[i] = i*dt + t0;
  
  /* Interp real/imag*/
  interp_spline_omp(h_aux->time, h_aux->real, h_aux->size, h->time, size, h->real);
  interp_spline_omp(h_aux->time, h_aux->imag, h_aux->size, h->time, size, h->imag);

  /* Compute phase and amplitude */
  /* Waveform_rmap(h, 1, 1); */ /* unwrap */
  Waveform_rmap(h, 1, 0);  /* do not unwrap */

  /* Free aux memory */
  Waveform_free (h_aux);
}

void Waveform_interp_ap (Waveform *h, const int size, const double t0, const double dt, const char *name)
{
  /* Alloc and init aux memory */  
  Waveform *h_aux;
  const int oldsize = h->size;
  Waveform_alloc(&h_aux, oldsize, "");
  memcpy(h_aux->time, h->time, oldsize * sizeof(double));
  memcpy(h_aux->real, h->real, oldsize * sizeof(double));
  memcpy(h_aux->imag, h->imag, oldsize * sizeof(double));
  memcpy(h_aux->ampli, h->ampli, oldsize * sizeof(double)); 
  memcpy(h_aux->phase, h->phase, oldsize * sizeof(double)); 

  /* Realloc arrays */
  h->size = size;
  if (strcmp(name, "")) strcpy(h->name, name);
  if (h->time) free(h->time);
  if (h->real) free(h->real);
  if (h->imag) free(h->imag);
  if (h->ampli) free(h->ampli);
  if (h->phase) free(h->phase);
  h->time  = malloc ( size * sizeof(double) );
  h->real  = malloc ( size * sizeof(double) );
  h->imag  = malloc ( size * sizeof(double) );
  h->ampli = malloc ( size * sizeof(double) );
  h->phase = malloc ( size * sizeof(double) );

  /* Fill new time array */
#pragma omp simd
  for (int i = 0; i < size; i++)
    h->time[i] = i*dt + t0;

  /* Interp phase and amplitude */
  interp_spline_omp(h_aux->time, h_aux->ampli, h_aux->size, h->time, size, h->ampli);
  interp_spline_omp(h_aux->time, h_aux->phase, h_aux->size, h->time, size, h->phase);
  
  /* Compute Real/Imag */
  Waveform_rmap (h, 0, 0); /* do not unwrap */
  
  /* Free aux memory */
  Waveform_free (h_aux);
}

void Waveform_output (Waveform *wav)
{
  FILE* fp;
  char fname[STRLEN];
  strcpy(fname,EOBPars->output_dir);
  strcat(fname,"/");
  strcat(fname,wav->name);
  strcat(fname,".txt");
  if ((fp = fopen(fname, "w+")) == NULL)
    errorexits("error opening file",wav->name);
  fprintf(fp, "# q=%e chizA=%e chizB=%e f0=%e\n",EOBPars->q,EOBPars->chi1,EOBPars->chi2,EOBPars->initial_frequency);
  fprintf(fp, "# M=%e LambdaA=[%e,%e,%e] LambdaBl2=[%e,%e,%e]\n",EOBPars->M,
	  EOBPars->LambdaAl2,EOBPars->LambdaAl3,EOBPars->LambdaAl4,
	  EOBPars->LambdaBl2,EOBPars->LambdaBl3,EOBPars->LambdaBl4);
  fprintf(fp, "# D=%e phi=%e iota=%e\n",EOBPars->distance,EOBPars->coalescence_angle,EOBPars->inclination);
  fprintf(fp, "# t:0 real:1 imag:2 Ampli:3 Phase:4\n");
  for (int i = 0; i < wav->size; i++) {
    fprintf(fp, "%.9e %.12e %.12e %.12e %.12e\n", wav->time[i], wav->real[i], wav->imag[i], wav->ampli[i], wav->phase[i]);
  }
  fclose(fp);
}

void Waveform_free (Waveform *wav)
{
  if (!wav) return;
  if (wav->real) free(wav->real);
  if (wav->imag) free(wav->imag);
  if (wav->ampli) free(wav->ampli);
  if (wav->phase) free(wav->phase);
  if (wav->time)   free(wav->time);
  free(wav);
}

/** WaveformFD (complex) */
void WaveformFD_alloc (WaveformFD **wav, const int size, const char *name)
{
  *wav = (WaveformFD *) calloc(1, sizeof(WaveformFD)); 
  if (wav == NULL)
    errorexit("Out of memory");
  (*wav)->preal = malloc ( size * sizeof(double) );
  (*wav)->pimag = malloc ( size * sizeof(double) );
  (*wav)->creal = malloc ( size * sizeof(double) );
  (*wav)->cimag = malloc ( size * sizeof(double) );
  (*wav)->ampli = malloc ( size * sizeof(double) );
  (*wav)->phase = malloc ( size * sizeof(double) );
  (*wav)->freq = malloc ( size * sizeof(double) );
  memset( (*wav)->preal, 0, size * sizeof(double) );
  memset( (*wav)->pimag, 0, size * sizeof(double) );
  memset( (*wav)->creal, 0, size * sizeof(double) );
  memset( (*wav)->cimag, 0, size * sizeof(double) );
  memset((*wav)->ampli, 0, size * sizeof(double) );
  memset((*wav)->phase, 0, size * sizeof(double) );
  (*wav)->size = size; 
  strcpy((*wav)->name,name);
}

void WaveformFD_push (WaveformFD **wav, int size)
{
  if ((*wav)->preal)  (*wav)->preal  = realloc ( (*wav)->preal,  size * sizeof(double) );
  if ((*wav)->pimag)  (*wav)->pimag  = realloc ( (*wav)->pimag,  size * sizeof(double) );
  if ((*wav)->creal)  (*wav)->creal  = realloc ( (*wav)->creal,  size * sizeof(double) );
  if ((*wav)->cimag)  (*wav)->cimag  = realloc ( (*wav)->cimag,  size * sizeof(double) );
  if ((*wav)->ampli) (*wav)->ampli = realloc ( (*wav)->ampli, size * sizeof(double) );
  if ((*wav)->phase) (*wav)->phase = realloc ( (*wav)->phase, size * sizeof(double) );
  if ((*wav)->freq)  (*wav)->freq  = realloc ( (*wav)->freq,  size * sizeof(double) );
  const int n  = (*wav)->size;
  const int dn = size - (*wav)->size;
  (*wav)->size = size; 
}

/* Compute real/imag <-> amplitude/phase */
void WaveformFD_rmap (WaveformFD *h, const int mode, const int unw)
{
  const int size = h->size;
  if (mode) {
    /** (Re, Im) -> (Amplitude, phase) */
#pragma omp simd
    for (int i = 0; i < size; i++) {
      double real = h->preal[i]+h->creal[i];
      double imag = h->pimag[i]+h->cimag[i];
      h->ampli[i] = sqrt( SQ(real) + SQ(imag) );
    }
#pragma omp simd
    for (int i = 0; i < size; i++) {
      double real = h->preal[i]+h->creal[i];
      double imag = h->cimag[i]+h->cimag[i];
      h->phase[i] = Pi - atan2(imag, real); /* exp(- i phi) => Pi  */
      if (unw) unwrap(h->phase, h->size); 
    }   
  } else {
    /** (Amplitude, phase) -> (Re, Im) */
    errorexit("cannot fill real/imag of cross/plus polarization.");
  }
}

void WaveformFD_interp_ap (WaveformFD *h, const int size, const double f0, const double df, const char *name)
{
  /* Alloc and init aux memory */  
  WaveformFD *h_aux;
  const int oldsize = h->size;
  WaveformFD_alloc(&h_aux, oldsize, "");
  memcpy(h_aux->freq, h->freq, oldsize * sizeof(double));
  memcpy(h_aux->preal, h->preal, oldsize * sizeof(double));
  memcpy(h_aux->pimag, h->pimag, oldsize * sizeof(double));
  memcpy(h_aux->creal, h->creal, oldsize * sizeof(double));
  memcpy(h_aux->cimag, h->cimag, oldsize * sizeof(double));
  memcpy(h_aux->ampli, h->ampli, oldsize * sizeof(double)); 
  memcpy(h_aux->phase, h->phase, oldsize * sizeof(double)); 

  /* Realloc arrays */
  h->size = size;
  if (strcmp(name, "")) strcpy(h->name, name);
  if (h->freq) free(h->freq);
  if (h->preal) free(h->preal);
  if (h->pimag) free(h->pimag);
  if (h->creal) free(h->creal);
  if (h->cimag) free(h->cimag);
  if (h->ampli) free(h->ampli);
  if (h->phase) free(h->phase);
  h->freq  = malloc ( size * sizeof(double) );
  h->preal  = malloc ( size * sizeof(double) );
  h->pimag  = malloc ( size * sizeof(double) );
  h->creal  = malloc ( size * sizeof(double) );
  h->cimag  = malloc ( size * sizeof(double) );
  h->ampli = malloc ( size * sizeof(double) );
  h->phase = malloc ( size * sizeof(double) );

  /* Fill new time array */
#pragma omp simd
  for (int i = 0; i < size; i++)
    h->freq[i] = i*df + f0;

  /* Interp phase and amplitude */
  interp_spline_omp(h_aux->freq, h_aux->ampli, h_aux->size, h->freq, size, h->ampli);
  interp_spline_omp(h_aux->freq, h_aux->phase, h_aux->size, h->freq, size, h->phase);
  
  /* Free aux memory */
  WaveformFD_free (h_aux);
}

void WaveformFD_output (WaveformFD *wav)
{
  FILE* fp;
  char fname[STRLEN];
  strcpy(fname,EOBPars->output_dir);
  strcat(fname,"/");
  strcat(fname,wav->name);
  strcat(fname,".txt");
  if ((fp = fopen(fname, "w+")) == NULL)
    errorexits("error opening file",wav->name);
  fprintf(fp, "# q=%e chizA=%e chizB=%e f0=%e\n",EOBPars->q,EOBPars->chi1,EOBPars->chi2,EOBPars->initial_frequency);
  fprintf(fp, "# M=%e LambdaA=[%e,%e,%e] LambdaBl2=[%e,%e,%e]\n",EOBPars->M,
	  EOBPars->LambdaAl2,EOBPars->LambdaAl3,EOBPars->LambdaAl4,
	  EOBPars->LambdaBl2,EOBPars->LambdaBl3,EOBPars->LambdaBl4);
  fprintf(fp, "# D=%e phi=%e iota=%e\n",EOBPars->distance,EOBPars->coalescence_angle,EOBPars->inclination);
  fprintf(fp, "# f:0 real:1 imag:2 Ampli:3 Phase:4\n");
  for (int i = 0; i < wav->size; i++) {
    fprintf(fp, "%.9e %.12e %.12e %.12e %.12e\n", wav->freq[i], wav->preal[i], wav->pimag[i], wav->creal[i], wav->cimag[i]);
  }
  fclose(fp);
}

void WaveformFD_free (WaveformFD *wav)
{
  if (!wav) return;
  if (wav->preal) free(wav->preal);
  if (wav->pimag) free(wav->pimag);
  if (wav->creal) free(wav->creal);
  if (wav->cimag) free(wav->cimag);
  if (wav->ampli) free(wav->ampli);
  if (wav->phase) free(wav->phase);
  if (wav->freq) free(wav->freq);
  free(wav);
}

/** Multipolar waveform (complex) */
void Waveform_lm_alloc (Waveform_lm **wav, int size, const char *name)
{
  *wav = (Waveform_lm *) calloc(1, sizeof(Waveform_lm)); 
  if (wav == NULL)
    errorexit("Out of memory");
  (*wav)->size = size; 
  set_multipolar_idx_mask((*wav)->kmask, KMAX, EOBPars->use_mode_lm, EOBPars->use_mode_lm_size, 0);
  (*wav)->time = malloc ( size * sizeof(double) );
  memset((*wav)->time, 0, size*sizeof(double));
  for (int k=0; k<KMAX; k++) {
    if ((*wav)->kmask[k]){
      (*wav)->ampli[k] = malloc ( size * sizeof(double) );
      (*wav)->phase[k] = malloc ( size * sizeof(double) );
      memset((*wav)->ampli[k], 0, size*sizeof(double));
      memset((*wav)->phase[k], 0, size*sizeof(double));
      strcpy((*wav)->name,name);
    }
  }
}

void Waveform_lm_push (Waveform_lm **wav, int size)
{
  const int n  = (*wav)->size;
  const int dn = size - (*wav)->size;
  (*wav)->time = realloc ( (*wav)->time, size * sizeof(double) );
  if ((*wav)->time == NULL) errorexit("Out of memory.");
  for (int k=0; k<KMAX; k++) {
    if ((*wav)->kmask[k]){
      (*wav)->ampli[k] = realloc ( (*wav)->ampli[k], size * sizeof(double) );
      if ((*wav)->ampli[k] == NULL) errorexit("Out of memory.");
      (*wav)->phase[k] = realloc ( (*wav)->phase[k], size * sizeof(double) );
      if ((*wav)->phase[k] == NULL) errorexit("Out of memory.");
    }
  }
  (*wav)->size = size;
}

/* Interp on uniform time array and overwrite a multipolar waveform */
void Waveform_lm_interp (Waveform_lm *hlm, const int size, const double t0, const double dt, const char *name)
{
  /* Alloc and init aux memory */  
  Waveform_lm *hlm_aux;
  const int oldsize = hlm->size;
  Waveform_lm_alloc(&hlm_aux, oldsize, "");
  memcpy(hlm_aux->time, hlm->time, oldsize * sizeof(double));
  for (int k = 0; k < KMAX; k++) {
    if (hlm->kmask[k]){
      memcpy(hlm_aux->ampli[k], hlm->ampli[k], oldsize * sizeof(double));
      memcpy(hlm_aux->phase[k], hlm->phase[k], oldsize * sizeof(double));
    }
  }

  /* Realloc arrays */
  hlm->size = size;
  if (strcmp(name, "")) strcpy(hlm->name, name);
  if (hlm->time) free(hlm->time); 
  hlm->time = malloc ( size * sizeof(double) );
  for (int k = 0; k < KMAX; k++) {
    if (hlm->kmask[k]){
      if (hlm->ampli[k]) free(hlm->ampli[k]);
      if (hlm->phase[k]) free(hlm->phase[k]);
      hlm->ampli[k] = malloc ( size * sizeof(double) );
      hlm->phase[k] = malloc ( size * sizeof(double) );
    }
  } 
  
  /* Fill new time array */
#pragma omp simd
  for (int i = 0; i < size; i++) 
    hlm->time[i] = i*dt + t0;
  
  /* Interp */
  for (int k = 0; k < KMAX; k++) 
    if (hlm->kmask[k])
      interp_spline_omp(hlm_aux->time, hlm_aux->ampli[k], hlm_aux->size, hlm->time, size, hlm->ampli[k]);
  for (int k = 0; k < KMAX; k++)
    if (hlm->kmask[k]) 
      interp_spline_omp(hlm_aux->time, hlm_aux->phase[k], hlm_aux->size, hlm->time, size, hlm->phase[k]);
  
  /* Free aux memory */
  Waveform_lm_free (hlm_aux);
}

void Waveform_lm_output (Waveform_lm *wav)
{
  char fname[STRLEN*2];
  const int n = wav->size;
  for (int k=0; k<KMAX; k++) {
    if (wav->kmask[k]) {      
      sprintf(fname,"%s/%s_l%01d_m%01d.txt",EOBPars->output_dir,wav->name,LINDEX[k],MINDEX[k]);
      FILE* fp;
      if ((fp = fopen(fname, "w+")) == NULL)
	      errorexits("error opening file",fname);
      for (int i = 0; i < n; i++) {
	      fprintf(fp, "%.9e %.16e %.16e\n", wav->time[i], wav->ampli[k][i], wav->phase[k][i]); 
      }
      fclose(fp);
    }
  }
}

void Waveform_lm_output_reim (Waveform_lm *wav)
{
  char fname[STRLEN*2];
  double re,im;
  const int n = wav->size;
  for (int k=0; k<KMAX; k++) {
    if (wav->kmask[k]) {
      sprintf(fname,"%s/%s_l%01d_m%01d_reim.txt",EOBPars->output_dir,wav->name,LINDEX[k],MINDEX[k]);
      FILE* fp;
      if ((fp = fopen(fname, "w+")) == NULL)
	errorexits("error opening file",fname);
      for (int i = 0; i < n; i++) {
	re = + wav->ampli[k][i] * cos(wav->phase[k][i]);
	im = - wav->ampli[k][i] * sin(wav->phase[k][i]);
	fprintf(fp, "%.9e %.12e %.12e\n", wav->time[i], re,im);
      }
      fclose(fp);
    }
  }
}

/* Extract a multipolar waveforms at times t >= to and t < tn, Alloc a new Waveform_lm var */
void Waveform_lm_extract (Waveform_lm *hlma, const double to, const double tn, Waveform_lm **hlmb, const char *name)
{
  /* Check limits */
  if (tn<to)
    errorexit("Bad choice of times: tn < to");
  if (to > hlma->time[hlma->size-1]) 
    errorexit("Nothing to extract, to > time[size-1]");
  if (tn < hlma->time[0]) 
    errorexit("Nothing to extract, tn < time[0]");

  /* Find indexes of closest elements to  (to, tn) */
  int io = 0;
  int in = hlma->size-1;
  if (to > hlma->time[0])
    io = find_point_bisection(to, hlma->size, hlma->time, 1);
  if (tn < hlma->time[hlma->size-1])
    in = find_point_bisection(tn, hlma->size, hlma->time, 0);
  
  /* Calculate the new size */ 
  const int N  = in-io;

#if (0)
  printf("size = %d t[0] = %.6e t[size-1]= %.6e to = %.6e tn = %.6e\n",
	 hlma->size,hlma->time[0],hlma->time[hlma->size-1],to,tn);
  printf("io = %d in = %d (N=%d) t[io]= %.6e t[in] = %.6e\n",
	 io, in, N, hlma->time[io],hlma->time[in]);
#endif
  
  /* Alloc output waveform b */
  Waveform_lm_alloc (hlmb, N, name);
  /* TODO: Parameters are not copied in the new wf !*/
  
  /* Copy the relevant part of a into b */
  for (int i = 0; i < N; i++) 
    (*hlmb)->time[i] = hlma->time[io + i]; 
  
  for (int k=0; k<KMAX; k++) {
    if ((*hlmb)->kmask[k]){ 
#pragma omp simd
      for (int i = 0; i < N; i++) {
	(*hlmb)->ampli[k][i] = hlma->ampli[k][io + i];
	(*hlmb)->phase[k][i] = hlma->phase[k][io + i];
      }
    }
  }
  
}

/* Join two multipolar waveforms at t = to */
void Waveform_lm_join (Waveform_lm *hlma, Waveform_lm *hlmb, double to)
{
  /* Time arrays are suppose to be ordered as
     hlma->time:  x x x x x x x x x 
     hlmb->time:       o o o o o o o o o 
     to        :                |
     But they do not need to overlap or be uniformly spaced.
     Note to can be 
     to > hlma->time[hlma->size-1] => extend the a waveform
     to < hlmb->time[0]            => join the whole b waveform
     Following checks enforce the above structure, if possible.
  */
  if (hlma->time[0] > hlmb->time[0]) {
    SWAPTRS( hlma, hlmb );
    if ((DEBUG) || (VERBOSE)) PRWARN("Swapped waveforms while joining.");
  }
  if (to > hlmb->time[hlmb->size-1]) {
    /* Nothing to join */
    if ((DEBUG) || (VERBOSE)) PRWARN("Joining time outside range. Waveforms not joined.");
    return;
  }
  if (to <= hlma->time[0]) {
    /* Nothing to join */
    if ((DEBUG) || (VERBOSE)) PRWARN("Joining time outside range. Waveforms not joined.");
    return;
  }

  /* Find indexes of closest elements to "to" 
     If the two arrays exactly overlap at "to" this should give: 
     time_a[ioa] = time_b[iob] */
  const int ioa = find_point_bisection(to, hlma->size, hlma->time, 1);
  int iob = find_point_bisection(to, hlmb->size, hlmb->time, 1);
  if ( DEQUAL(hlma->time[ioa], hlmb->time[iob], 1e-10) ) iob++;
  
  /* Calculate the new size N */
  const int Nb = hlmb->size - iob;
  const int N  = ioa + Nb;

  /* Resize a */
  Waveform_lm_push (&hlma, N);
  
  /* Copy the relevant part of b into a */
#pragma omp simd
  for (int i = 0; i < Nb; i++) {
    hlma->time[ioa + i] = hlmb->time[iob + i]; 
  }

  for (int k=0; k<KMAX; k++) {
    if (hlma->kmask[k]){ 
#pragma omp simd
      for (int i = 0; i < Nb; i++) {
	hlma->ampli[k][ioa + i] = hlmb->ampli[k][iob + i];
	hlma->phase[k][ioa + i] = hlmb->phase[k][iob + i];
      }
    }
  }
    
#if (0)
  printf("Waveform (a) i = %d time[i] = %.6e \n",ioa, hlma->time[ioa]);
  printf("Waveform (b) i = %d time[i] = %.6e \n",iob, hlmb->time[iob]);
  printf("Total size (a)+(b) = %d + %d = %d (%d)\n",Nb, iob, N, hlma->size);
  printf("%.6e - %.6e = %.6e\n",hlma->time[hlma->size-1], hlmb->time[hlmb->size-1],hlma->time[hlma->size-1]-hlmb->time[hlmb->size-1]);
#endif  
  
}

void Waveform_lm_free (Waveform_lm *wav)
{
  if (!wav) return;
  for (int k=0; k<KMAX; k++) {
    if (wav->kmask[k]) { 
      if (wav->ampli[k]) free(wav->ampli[k]);
      if (wav->phase[k]) free(wav->phase[k]);
    }
  }
  if (wav->time) free(wav->time);
  free(wav);
}

/** Multipolar waveform (complex) */
void WaveformFD_lm_alloc (WaveformFD_lm **wav, int size, const char *name)
{
  *wav = (WaveformFD_lm *) calloc(1, sizeof(WaveformFD_lm)); 
  if (wav == NULL)
    errorexit("Out of memory");
  (*wav)->size = size; 
  set_multipolar_idx_mask((*wav)->kmask, KMAX, EOBPars->use_mode_lm, EOBPars->use_mode_lm_size, 0);
  (*wav)->freq = malloc ( size * sizeof(double) );
  memset((*wav)->freq, 0, size*sizeof(double));
  for (int k=0; k<KMAX; k++) {
    if ((*wav)->kmask[k]) { 
      (*wav)->ampli[k] = malloc ( size * sizeof(double) );
      (*wav)->phase[k] = malloc ( size * sizeof(double) );
      (*wav)->F[k] = malloc ( size * sizeof(double) );
      (*wav)->Fdot[k] = malloc ( size * sizeof(double) );
      memset((*wav)->ampli[k], 0, size*sizeof(double));
      memset((*wav)->phase[k], 0, size*sizeof(double));
      memset((*wav)->F[k], 0, size*sizeof(double));
      memset((*wav)->Fdot[k], 0, size*sizeof(double));
      strcpy((*wav)->name,name);
    }
  }
}

void WaveformFD_lm_push (WaveformFD_lm **wav, int size)
{
  const int n  = (*wav)->size;
  const int dn = size - (*wav)->size;
  (*wav)->freq = realloc ( (*wav)->freq, size * sizeof(double) );
  if ((*wav)->freq == NULL) errorexit("Out of memory.");
  for (int k=0; k<KMAX; k++) {
    if ((*wav)->kmask[k]) { 
      (*wav)->ampli[k] = realloc ( (*wav)->ampli[k], size * sizeof(double) );
      if ((*wav)->ampli[k] == NULL) errorexit("Out of memory.");
      (*wav)->phase[k] = realloc ( (*wav)->phase[k], size * sizeof(double) );
      if ((*wav)->phase[k] == NULL) errorexit("Out of memory.");
      (*wav)->F[k] = realloc ( (*wav)->F[k], size * sizeof(double) );
      if ((*wav)->F[k] == NULL) errorexit("Out of memory.");
      (*wav)->Fdot[k] = realloc ( (*wav)->Fdot[k], size * sizeof(double) );
      if ((*wav)->Fdot[k] == NULL) errorexit("Out of memory.");
    }
  }
  (*wav)->size = size;
}

/* Interp on uniform freq array and overwrite a multipolar waveform 
   Notes 
   - each mode has its own frequency array hlm->F[k]
   - the uniform array is stored in hlm->freq 
   - use interp_spline_checklim() and set to zero points outside interp range */
void WaveformFD_lm_interp_ap (WaveformFD_lm *hlm, const int size, const double f0, const double df, const char *name)
{
  /* Alloc and init aux memory */  
  WaveformFD_lm *hlm_aux;
  const int oldsize = hlm->size;
  WaveformFD_lm_alloc(&hlm_aux, oldsize, "");
  for (int k = 0; k < KMAX; k++) {
    if (hlm->kmask[k]) {  
      memcpy(hlm_aux->ampli[k], hlm->ampli[k], oldsize * sizeof(double));
      memcpy(hlm_aux->phase[k], hlm->phase[k], oldsize * sizeof(double));
      memcpy(hlm_aux->F[k], hlm->F[k], oldsize * sizeof(double));
      //memcpy(hlm_aux->Fdot[k], hlm->Fdot[k], oldsize * sizeof(double));
    }
  }

  /* Realloc arrays */
  hlm->size = size;
  if (strcmp(name, "")) strcpy(hlm->name, name);
  if (hlm->freq) free(hlm->freq); 
  hlm->freq = malloc ( size * sizeof(double) );
  for (int k = 0; k < KMAX; k++) {
    if (hlm->kmask[k]) {  
      if (hlm->ampli[k]) free(hlm->ampli[k]);
      if (hlm->phase[k]) free(hlm->phase[k]);
      if (hlm->F[k]) free(hlm->F[k]);
      if (hlm->Fdot[k]) free(hlm->Fdot[k]);
      
      hlm->ampli[k] = malloc ( size * sizeof(double) );
      hlm->phase[k] = malloc ( size * sizeof(double) );
      hlm->F[k] = NULL; // this info is lost in the interp.
      hlm->Fdot[k] = NULL;// this info is lost in the interp.
    } 
  }

  /* Fill new freq array */
#pragma omp simd
  for (int i = 0; i < size; i++) 
    hlm->freq[i] = i*df + f0;
  
  /* Interp */
  for (int k = 0; k < KMAX; k++){
    if (hlm->kmask[k]){
      //interp_spline_omp(hlm_aux->F[k], hlm_aux->ampli[k], hlm_aux->size, hlm->freq, size, hlm->ampli[k]);
      interp_spline_checklim(hlm_aux->F[k], hlm_aux->ampli[k], hlm_aux->size, hlm->freq, size, hlm->ampli[k]);
    }
  }
  for (int k = 0; k < KMAX; k++){
    if (hlm->kmask[k]) {
      //interp_spline_omp(hlm_aux->F[k], hlm_aux->phase[k], hlm_aux->size, hlm->freq, size, hlm->phase[k]);
      interp_spline_checklim(hlm_aux->F[k], hlm_aux->phase[k], hlm_aux->size, hlm->freq, size, hlm->phase[k]);
    }
  }
  /* Free aux memory */
  WaveformFD_lm_free (hlm_aux);
}

void WaveformFD_lm_output (WaveformFD_lm *wav)
{
  int kmask[KMAX];
  set_multipolar_idx_mask(kmask, KMAX, EOBPars->output_lm, EOBPars->output_lm_size, 0);  
  char fname[STRLEN*2];
  const int n = wav->size;
  for (int k=0; k<KMAX; k++) {
    if (kmask[k]) {      
      sprintf(fname,"%s/FD_%s_l%01d_m%01d.txt",EOBPars->output_dir,wav->name,LINDEX[k],MINDEX[k]);
      FILE* fp;
      if ((fp = fopen(fname, "w+")) == NULL)
	errorexits("error opening file",fname);
      for (int i = 0; i < n; i++) {
	fprintf(fp, "%.9e %.16e %.16e\n", wav->freq[i], wav->ampli[k][i], wav->phase[k][i]); 
      }
      fclose(fp);
    }
  }
}

void WaveformFD_lm_output_reim (WaveformFD_lm *wav)
{
  int kmask[KMAX];
  set_multipolar_idx_mask(kmask, KMAX, EOBPars->output_lm, EOBPars->output_lm_size, 0);  
  char fname[STRLEN*2];
  double re,im;
  const int n = wav->size;
  for (int k=0; k<KMAX; k++) {
    if (wav->kmask[k]) {
      sprintf(fname,"%s/FD_%s_l%01d_m%01d_reim.txt",EOBPars->output_dir,wav->name,LINDEX[k],MINDEX[k]);
      FILE* fp;
      if ((fp = fopen(fname, "w+")) == NULL)
	errorexits("error opening file",fname);
      for (int i = 0; i < n; i++) {
	re = + wav->ampli[k][i] * cos(wav->phase[k][i]);
	im = - wav->ampli[k][i] * sin(wav->phase[k][i]);
	fprintf(fp, "%.9e %.12e %.12e\n", wav->freq[i], re,im);
      }
      fclose(fp);
    }
  }
}

void WaveformFD_lm_free (WaveformFD_lm *wav)
{
  if (!wav) return;
  for (int k=0; k<KMAX; k++) {
    if (wav->kmask[k]) {
      if (wav->ampli[k]) free(wav->ampli[k]);
      if (wav->phase[k]) free(wav->phase[k]);
      if (wav->F[k]) free(wav->F[k]);
      if (wav->Fdot[k]) free(wav->Fdot[k]);
    }
  }
  if (wav->freq) free(wav->freq);
  free(wav);
}

/** Multipolar waveform at time point (complex) */
void Waveform_lm_t_alloc (Waveform_lm_t **wav)
{
  *wav = (Waveform_lm_t *) calloc(1, sizeof(Waveform_lm_t)); 
  if (wav == NULL) errorexit("Out of memory");
  (*wav)->time = 0.;
  (*wav)->freq = 0.;
  set_multipolar_idx_mask ((*wav)->kmask, KMAX, EOBPars->use_mode_lm, EOBPars->use_mode_lm_size, 0); 
}

void Waveform_lm_t_free (Waveform_lm_t *wav)
{
  if (!wav) return;
  free(wav);
}

/** Dynamics */
void Dynamics_alloc (Dynamics **dyn, int size, const char *name)
{
  (*dyn) = (Dynamics *) calloc(1, sizeof(Dynamics)); 
  if (dyn == NULL)
    errorexit("Out of memory");
  strcpy((*dyn)->name,name);
  (*dyn)->size = size; 
  (*dyn)->time = malloc ( size * sizeof(double) );
  memset((*dyn)->time, 0, size*sizeof(double));
  for (int v = 0; v < EOB_DYNAMICS_NVARS; v++) {
    (*dyn)->data[v] = malloc ( size * sizeof(double) );
    memset((*dyn)->data[v], 0, size*sizeof(double));
  }
}

void Dynamics_push (Dynamics **dyn, int size)
{
  const int n  = (*dyn)->size;
  const int dn = size - (*dyn)->size;
  (*dyn)->time = realloc ( (*dyn)->time, size * sizeof(double) );
  for (int v = 0; v < EOB_DYNAMICS_NVARS; v++) {
    (*dyn)->data[v] = realloc ( (*dyn)->data[v], size * sizeof(double) );
    if ((*dyn)->data[v] == NULL) errorexit("Out of memory.");
    /* if (dn>0) memset( (*dyn)->data[v] + n, 0, dn * sizeof(double) ); */
  }
  (*dyn)->size = size; 
}

/* Interp and overwrite a multipolar waveform */
void Dynamics_interp (Dynamics *dyn, const int size, const double t0, const double dt, const char *name)
{
  /* Alloc and init aux memory */  
  Dynamics *dyn_aux;
  const int oldsize = dyn->size;
  Dynamics_alloc(&dyn_aux, oldsize, "");
  memcpy(dyn_aux->time, dyn->time, oldsize * sizeof(double));
  for (int v = 0; v < EOB_DYNAMICS_NVARS; v++)
    memcpy(dyn_aux->data[v], dyn->data[v], oldsize * sizeof(double));
  
  /* Overwrite and realloc arrays */
  dyn->dt   = dt;
  dyn->size = size; 
  if (strcmp(name, "")) strcpy(dyn->name, name);
  if (dyn->time) free(dyn->time);
  dyn->time = malloc ( size * sizeof(double) );
  for (int v = 0; v < EOB_DYNAMICS_NVARS; v++) {
    if (dyn->data[v]) free(dyn->data[v]);
    dyn->data[v] = malloc ( size * sizeof(double) );
    /* memset(dyn->data[v], 0., size*sizeof(double)); */
  }        

  /* Fill new time array */
  for (int i = 0; i < size; i++)
    dyn->time[i] = i*dt + t0;
  
  /* Interp */
  for (int k = 0; k < EOB_DYNAMICS_NVARS; k++) 
    interp_spline_omp(dyn_aux->time, dyn_aux->data[k], dyn_aux->size, dyn->time, size, dyn->data[k]);

  /* Free aux memory */
  Dynamics_free (dyn_aux);
}

void Dynamics_output (Dynamics *dyn)
{
  FILE* fp; 
  char fname[STRLEN];
  strcpy(fname,EOBPars->output_dir);
  strcat(fname,"/");
  strcat(fname,dyn->name);
  strcat(fname,".txt");
  if ((fp = fopen(fname, "w+")) == NULL)
    errorexits("error opening file",dyn->name);
  fprintf(fp, "# t:0");
  for (int v = 0; v < EOB_DYNAMICS_NVARS; v++)
    fprintf(fp, " %s:%d",eob_var[v],v+1);
  fprintf(fp, "\n");
  for (int i = 0; i < dyn->size; i++) {
    fprintf(fp, "%.9e", dyn->time[i]);
    for (int v = 0; v < EOB_DYNAMICS_NVARS; v++)
      fprintf(fp, " %.12e", dyn->data[v][i]);
    fprintf(fp, "\n");
  }
  fclose(fp);
}

/* Extract Dynamics at times t >= to and t < tn, Alloc a new Dynamics var */
void Dynamics_extract (Dynamics *dyna, const double to, const double tn, Dynamics **dynb, const char *name)
{
  /* Check limits */
  if (tn<to)
    errorexit("Bad choice of times: tn < to");
  if (to > dyna->time[dyna->size-1]) 
    errorexit("Nothing to extract, to > time[size-1]");
  if (tn < dyna->time[0]) 
    errorexit("Nothing to extract, tn < time[0]");

  /* Find indexes of closest elements to  (to, tn) */
  int io = 0;
  int in = dyna->size-1;
  if (to > dyna->time[0])
    io = find_point_bisection(to, dyna->size, dyna->time, 1);
  if (tn < dyna->time[dyna->size-1])
    in = find_point_bisection(tn, dyna->size, dyna->time, 0);
  
  /* Calculate the new size */ 
  const int N  = in-io;

#if (0)
  printf("size = %d t[0] = %.6e t[size-1]= %.6e to = %.6e tn = %.6e\n",
	 dyna->size,dyna->time[0],dyna->time[dyna->size-1],to,tn);
  printf("io = %d in = %d (N=%d) t[io]= %.6e t[in] = %.6e\n",
	 io, in, N, dyna->time[io],dyna->time[in]);
#endif
  
  /* Alloc output waveform b */
  Dynamics_alloc (dynb, N, name); 

  /* TODO: Parameters are not copied in the new wf !*/
  /*
    (*dynb) = (Dynamics *) calloc(1, sizeof(Dynamics)); 
    if (dynb == NULL)
    errorexit("Out of memory");
    strcpy((*dynb)->name,name);
    (*dynb)->size = N;
    memcpy(*dynb, dyna, sizeof(Dynamics)); // copy parameters 
    (*dynb)->time = malloc ( N * sizeof(double) );
    for (int v = 0; v < EOB_DYNAMICS_NVARS; v++) 
    (*dynb)->data[v] = malloc ( N * sizeof(double) );
  */
    
  /* Copy the relevant part of a into b */
  for (int i = 0; i < N; i++) 
    (*dynb)->time[i] = dyna->time[io + i]; 
  
  for (int v = 0; v < EOB_DYNAMICS_NVARS; v++) {
    for (int i = 0; i < N; i++) {
      (*dynb)->data[v][i] = dyna->data[v][io + i];
    }
  }
  
}

/* Join two dynamics time series at t = to */
void Dynamics_join (Dynamics *dyna, Dynamics *dynb, double to)
{
  /* Time arrays are suppose to be ordered as
     dyna->time:  x x x x x x x x x 
     dynb->time:       o o o o o o o o o 
     to        :                |
     But they do not need to overlap or be uniformly spaced.
     Note to can be 
     to > dyna->time[hlma->size-1] => extend the dynamics data
     to < dynb->time[0]            => join the whole b dynamics
     Following checks enforce the above structure, if possible.
  */
  if (dyna->time[0] > dynb->time[0]) {
    SWAPTRS( dyna, dynb );
    if ((DEBUG) || (VERBOSE)) PRWARN("Swapped dynamics while joining.");
  }
  if (to > dynb->time[dynb->size-1]) {
    /* Nothing to join */
    if ((DEBUG) || (VERBOSE)) PRWARN("Joining time outside range. Dynamics not joined.");
    return;
  }
  if (to <= dyna->time[0]) {
    /* Nothing to join */
    if ((DEBUG) || (VERBOSE)) PRWARN("Joining time outside range. Dynamics not joined.");
    return;
  }

  /* Find indexes of closest elements to to */
  const int ioa = find_point_bisection(to, dyna->size, dyna->time, 1);
  int iob = find_point_bisection(to, dynb->size, dynb->time, 1);
  if ( DEQUAL(dyna->time[ioa], dynb->time[iob], 1e-10) ) iob++;
  
  /* Calculate the new size */
  const int Nb = dynb->size - iob;
  const int N  = ioa + Nb;

  /* Resize a */
  Dynamics_push (&dyna, N);

  /* Copy the relevant part of b into a */
#pragma omp simd
  for (int i = 0; i < Nb; i++) 
    dyna->time[ioa + i] = dynb->time[iob + i]; 
  
  for (int v = 0; v < EOB_DYNAMICS_NVARS; v++) {
#pragma omp simd
    for (int i = 0; i < Nb; i++) {
      dyna->data[v][ioa + i] = dynb->data[v][iob + i];
    }
 }

#if (0)
  printf("Dynamics (a) i = %d time[i] = %.6e \n",ioa, dyna->time[ioa]);
  printf("Dynamics (b) i = %d time[i] = %.6e \n",iob, dynb->time[iob]);
  printf("Total size (a)+(b) = %d + %d = %d (%d)\n",Nb, iob, N, dyna->size);
  printf("%.6e - %.6e = %.6e\n",dyna->time[dyna->size-1], dynb->time[dynb->size-1],dyna->time[dyna->size-1]-dynb->time[dynb->size-1]);
#endif 
    
}

void Dynamics_free (Dynamics *dyn)
{
  if (!dyn) return;
  if (dyn->time) free(dyn->time);
  for (int v = 0; v < EOB_DYNAMICS_NVARS; v++)
    if (dyn->data[v]) free(dyn->data[v]);
  free(dyn);
}

/** Sync some quick access parameters in dyn with parameter database 
    to be used carefully */
void Dynamics_set_params (Dynamics *dyn)
{
  dyn->store = 0;
  dyn->size  = EOBPars->size;
  dyn->M     = EOBPars->M;
  dyn->nu    = EOBPars->nu;
  dyn->q     = EOBPars->q;
  dyn->X1    = EOBPars->X1;  
  dyn->X2    = EOBPars->X2;
  dyn->chi1  = EOBPars->chi1;
  dyn->chi2  = EOBPars->chi2;
  dyn->S1    = EOBPars->S1;
  dyn->S2    = EOBPars->S2;
  dyn->S     = EOBPars->S;
  dyn->Sstar = EOBPars->Sstar;
  dyn->a1    = EOBPars->a1;
  dyn->a2    = EOBPars->a2; 
  dyn->aK2   = EOBPars->aK2; 
  dyn->C_Q1  = EOBPars->C_Q1;
  dyn->C_Q2  = EOBPars->C_Q2;
  dyn->C_Oct1 = EOBPars->C_Oct1;
  dyn->C_Oct2 = EOBPars->C_Oct2;
  dyn->C_Hex1 = EOBPars->C_Hex1;
  dyn->C_Hex2 = EOBPars->C_Hex2;
  dyn->a6c   = EOBPars->a6c;
  dyn->cN3LO = EOBPars->cN3LO;
  dyn->rLR   = EOBPars->rLR;
  dyn->rLR_tidal = EOBPars->rLR_tidal;
  dyn->rLSO  = EOBPars->rLSO;
  dyn->kapA2 = EOBPars->kapA2;
  dyn->kapA3 = EOBPars->kapA3;
  dyn->kapA4 = EOBPars->kapA4;
  dyn->kapB2 = EOBPars->kapB2;
  dyn->kapB3 = EOBPars->kapB3;
  dyn->kapB4 = EOBPars->kapB4;
  dyn->kapT2 = EOBPars->kapT2;
  dyn->kapT3 = EOBPars->kapT3;
  dyn->kapT4 = EOBPars->kapT4;
  dyn->khatA2 = EOBPars->khatA2;
  dyn->khatB2 = EOBPars->khatB2;
  dyn->kapA2j = EOBPars->japA2;
  dyn->kapB2j = EOBPars->japB2;
  dyn->kapT2j = EOBPars->japT2;
  dyn->bar_alph2_1 = EOBPars->bar_alph2_1;
  dyn->bar_alph2_2 = EOBPars->bar_alph2_2;
  dyn->bar_alph3_1 = EOBPars->bar_alph3_1;
  dyn->bar_alph3_2 = EOBPars->bar_alph3_2;
  dyn->bar_alph2j_1 = EOBPars->bar_alph2j_1;
  dyn->pGSF_tidal = EOBPars->pGSF_tidal;
  dyn->Mbhf = EOBPars->Mbhf;
  dyn->abhf = EOBPars->abhf;
  dyn->use_tidal = EOBPars->use_tidal;
  dyn->use_tidal_gravitomagnetic = EOBPars->use_tidal_gravitomagnetic;
  dyn->use_spins = EOBPars->use_spins;
  dyn->dt        = EOBPars->dt;
  dyn->t_stop    = EOBPars->ode_tmax;
}

/** NQC data */
void NQCdata_alloc (NQCdata **nqc)
{
  *nqc = (NQCdata *) calloc(1, sizeof(NQCdata));
  if (nqc == NULL)
    errorexit("Out of memory");
  (*nqc)->flx = (NQCcoefs *) calloc(1, sizeof(NQCcoefs));
  if ((*nqc)->flx == NULL)
    errorexit("Out of memory");
  (*nqc)->hlm = (NQCcoefs *) calloc(1, sizeof(NQCcoefs)); 
  if ((*nqc)->hlm == NULL)
    errorexit("Out of memory");
}

void NQCdata_free (NQCdata *nqc)
{
  if (!nqc) return;
  if (nqc->flx) free (nqc->flx);
  if (nqc->hlm) free (nqc->hlm);
  if (nqc)      free (nqc);
}

/** Stationary Phase Approximation of multipolar wvf */
void SPA(Waveform_lm *TDlm, WaveformFD_lm *FDlm)
{
  const int size = TDlm->size;
  double tmpf0 = EOBPars->initial_frequency;  
  double tmpdf = EOBPars->df;
  double tmpsrate = EOBPars->srate_interp/2.;

  /* if necessary, transform f0, df and srate_interp to geom units */
  if (!(EOBPars->use_geometric_units)){
    double Msun = EOBPars->M;
    double conv = time_units_factor(Msun);
    tmpsrate = tmpsrate/conv; //FIXME: this must be transformed in geom units 
    tmpf0    = tmpf0/conv;
    tmpdf    = tmpdf/conv;
  }
  const double half_srate_interp = tmpsrate;
  const double f0 = tmpf0;
  const double df = tmpdf;
  //double Fmin=0., Fmax=0.;
  //int nmax[KMAX]; // Fmax index
  int nmin[KMAX]; // Fmin > f0 index

  const double Pio4 = Pi/4.;
  int *activemode = TDlm->kmask;
  //int activemode[KMAX];
  //set_multipolar_idx_mask (activemode, KMAX, 
  //EOBPars->use_mode_lm, EOBPars->use_mode_lm_size, 1);
  
  /* For each active mode... */
  for (int k = 0; k < KMAX; k++ ) {
    if (!activemode[k]) continue;
    
    /* Compute frequencies */
    /* D0_x_2(phaset, time, size, *F);  
       D0_x_2(*F, time, size, Fdot); */
    D0_x_4(TDlm->phase[k], TDlm->time, size, FDlm->F[k]);
    D0_x_4(FDlm->F[k], TDlm->time, size, FDlm->Fdot[k]);
    
    for (int i=0; i < size; i++){
      FDlm->F[k][i]     = FDlm->F[k][i]/(TwoPi);
      FDlm->Fdot[k][i]  = FDlm->Fdot[k][i]/(TwoPi);
      FDlm->phase[k][i] = (TwoPi * FDlm->F[k][i] * TDlm->time[i] - TDlm->phase[k][i]- Pio4);
      FDlm->ampli[k][i] = TDlm->ampli[k][i]/sqrt(fabs(FDlm->Fdot[k][i])); 
    }

    /* Track multipoles w\ Fmin > f0 */
    nmin[k] = (FDlm->F[k][0]>f0)?(1):(0);

    /* Get last index until Fdot is monotonically increasing (for attachment) */
    //int i_aux = 0;
    //while(FDlm->Fdot[i_aux+1] > FDlm->Fdot[i_aux]) i_aux++;
    //*newsize = i_aux +1; // i_axu = n
    //SB: I changed the logic above to avoid overflow of the array Fdot.
    //    please check n is correctly set (could be offset of 1...)
    int n = 1;
    while(FDlm->Fdot[k][n] > FDlm->Fdot[k][n-1] && n<size) n++;
    //nmax[k] = n; 

    /* Prolong the waveform, if necessary  
       - Fill the points between Fmax and half_srate_interp such that 
         FDlm->F[k][n] > half_srate_interp (strictly larger)
	 for later interp;
       - The amplitude is expected to behave as 1./f asymptotically 
         (see eg 3.31a of arXiv:gr-qc/0001023);
       - For the phase, we express it as phi(f) = (a + b*f);
    */
    if (FDlm->F[k][n] < half_srate_interp) { 
      const int n1 = n-1;
      double Fn1 = FDlm->F[k][n-1];
      double An1 = FDlm->ampli[k][n-1];
      double pn1 =  FDlm->phase[k][n-1];
      //double c = 0.;
      double b = TwoPi * TDlm->time[n-1];
      double dfk = (half_srate_interp - FDlm->F[k][n1])/(size-n-1); // one point more
      for (int i=n; i < size; i++){
        FDlm->F[k][i] = Fn1 + (i-n1)*dfk;
        FDlm->ampli[k][i] = An1/(FDlm->F[k][i])*Fn1;  // \approx 1/f
        //FDlm->phase[k][i] = (pn1 + b*(FDlm->F[k][i] - Fn1))/(1 + c*(FDlm->F[k][i] - Fn1));
        FDlm->phase[k][i] = (pn1 + b*(FDlm->F[k][i] - Fn1));
      }
    } else{
      /*update values with index >= n to assure F is formally increasing (for interpolation).
        Note that we do not care about frequencies higher than srate_interp/2, so
        we can fill F however we want. We also fill A, to avoid inf which may appear
        and would mess up the spline 
      */
      double Fn1 = FDlm->F[k][n-1];
      const int n1 = n-1;
      for (int i=n; i < size; i++){
        FDlm->F[k][i] = Fn1 + (i-n1);
        FDlm->ampli[k][i] = FDlm->ampli[k][n1];
      }
    }

    /* Absolute min/max */
    //Fmin = MIN(Fmin,FDlm->F[k][0]);
    //Fmax = MIN(Fmax,FDlm->F[k][n]);
  }

  /* Determine the output frequency array */
  const int interp_size = get_uniform_size(half_srate_interp,f0, df);
  
  /* Interpolate each mode */
  WaveformFD_lm_interp_ap (FDlm, interp_size, f0, df, "");

  /* Correct Fmin > f0 */
  for (int k = 0; k < KMAX; k++ ) {
    if (!activemode[k]) continue;
    if (nmin[k]) {
      /* f0 < Fmin , the interpolation returned 0s at those points
	 -> fill points with linear extrapolation starting from innermost */
      int i0 = 0;
      while(FDlm->ampli[k][i0]==0. && i0<size) i0++;
      for (int i = i0-1; i>=0; i--){
        FDlm->ampli[k][i] = 2*FDlm->ampli[k][i+1] - FDlm->ampli[k][i+2]; //TODO: check sign!
        FDlm->phase[k][i] = 2*FDlm->phase[k][i+1] - FDlm->phase[k][i+2]; 
      }
    }
  }
  
}

/** (h+, hx) polarizations from the multipolar waveform, FD, all active modes (interpolate after SPA, needed to correctly add modes together) */
void compute_hpc_FD(WaveformFD_lm *hflm, double nu, double M, double distance, double amplitude_prefactor, double phi, double iota, WaveformFD *hpc)
{  
#ifdef _OPENMP
  if (USETIMERS) openmp_timer_start("compute_hpc_FD");
#endif
  
  double Y_real[KMAX] = {0};
  double Y_imag[KMAX] = {0};
  double Y_real_mneg[KMAX] = {0};
  double Y_imag_mneg[KMAX] = {0};
  static const int mneg = 1; /* m>0 modes only, add m<0 modes afterwards */
  double conv=1.;
  double sumr=0., sumi=0.;
  double sumpr=0., sumpi=0., sumcr=0., sumci=0.;
  
  int *activemode = hflm->kmask;
  
  double Msun = M;
  if (!(EOBPars->use_geometric_units)) {
    Msun = M/MSUN_S;
    conv = time_units_factor(Msun);
    for (int i = 0; i < hflm->size; i++)
	    hflm->freq[i] = hflm->freq[i]*conv; //CHECKME!
  }
  
    /* Precompute Ylm */
  for (int k = 0; k < KMAX; k++ ) {
    if (!activemode[k]) continue;
    spinsphericalharm(&Y_real[k], &Y_imag[k], -2, LINDEX[k], MINDEX[k], phi, iota);
    /* add m<0 modes */
    if ( (mneg) && (MINDEX[k]!=0) ) 
      spinsphericalharm(&Y_real_mneg[k], &Y_imag_mneg[k], -2, LINDEX[k], -MINDEX[k], phi, iota); 
  }
  
  /* Sum up  hlm * Ylm 
   * Note because EOB code defines phase>0, 
   * but the convention is h_lm = A_lm Exp[-I phi_lm] we have
   * h_{l,m>=0} = A_lm ( cos(phi_lm) - I*sin(phi_lm) ) for m>0 and
   * h_{l,m<0}  = (-)^l A_l|m| ( cos(phi_l|m|) + I*sin(phi_l|m|) ) for m<0 below
   * We now agree with, e.g., LALSimSphHarmMode.c: 64-74
   */
#if (0)
  
  //SB: Taken from TD (?) to be checked
  
  for (int i = 0; i < hflm->size; i++) {
    hpc->freq[i] = hflm->freq[i]*M; 
    sumr = sumi = 0.;
    
    /* Loop over modes */
    for (int k = 0; k < KMAX; k++ ) {
      if (!activemode[k]) continue;
      double Aki  = amplitude_prefactor * hflm->ampli[k][i];
      double cosPhi = cos( hflm->phase[k][i] );
      double sinPhi = sin( hflm->phase[k][i] );
      sumr += Aki*(cosPhi*Y_real[k] + sinPhi*Y_imag[k]);
      sumi += Aki*(cosPhi*Y_imag[k] - sinPhi*Y_real[k]); 
      
      /* add m<0 modes */
      if ( (mneg) && (MINDEX[k]!=0) ) { 
	/* H_{l-m} = (-)^l H^{*}_{lm} */
	if (LINDEX[k] % 2) {
	  sumr -= Aki*(cosPhi*Y_real_mneg[k] - sinPhi*Y_imag_mneg[k]); 
	  sumi -= Aki*(cosPhi*Y_imag_mneg[k] + sinPhi*Y_real_mneg[k]); 
	}
	else { 
	  sumr += Aki*(cosPhi*Y_real_mneg[k] - sinPhi*Y_imag_mneg[k]);
	  sumi += Aki*(cosPhi*Y_imag_mneg[k] + sinPhi*Y_real_mneg[k]); 
	}
      }    
    }

    /* h = h+ - i hx */
    hpc->real[i] = sumr;
    hpc->imag[i] = -sumi;
  }
  
  /* Compute amplitude & Phase */
  //Waveform_rmap(hpc,1,1); // with unwrap
  Waveform_rmap(hpc,1,0); // without
  
#else
  
  const double pm3 = 3.*Pi/2.;
  for (int i = 0; i < hflm->size; i++) {
    hpc->freq[i] = hflm->freq[i]; 
    sumpr = sumpi = sumcr = sumci = 0.;
    
    for (int k=0; k< KMAX; k++){
      if (!activemode[k]) continue;
      double Aki  = 0.5*amplitude_prefactor * hflm->ampli[k][i];
      double cosPhi = cos( hflm->phase[k][i] );
      double sinPhi = sin( hflm->phase[k][i] );
      double cosPhipm3 = cos( hflm->phase[k][i] + pm3 );
      double sinPhipm3 = sin( hflm->phase[k][i] + pm3 );
      
      /* H_{l-m} = (-)^l H^{*}_{lm} */
      if (LINDEX[k] % 2) {
        sumpr +=  Aki * (cosPhi*(Y_real[k] - Y_real_mneg[k]) + sinPhi* (Y_imag[k] + Y_imag_mneg[k]));
        sumpi +=  Aki * (cosPhi*(-Y_imag[k] - Y_imag_mneg[k]) + sinPhi* (Y_real[k] - Y_real_mneg[k]));
        sumcr += -Aki * (cosPhipm3*(Y_real[k] + Y_real_mneg[k]) + sinPhipm3* (Y_imag[k] - Y_imag_mneg[k]));
        sumci += -Aki * (-cosPhipm3*(Y_real[k] - Y_real_mneg[k]) + sinPhipm3* (Y_imag[k] + Y_imag_mneg[k]));
      } else {
        sumpr +=  Aki * (cosPhi* (Y_real[k] + Y_real_mneg[k]) + sinPhi* (Y_imag[k] - Y_imag_mneg[k]));
        sumpi +=  Aki * (cosPhi* (-Y_imag[k] + Y_imag_mneg[k]) + sinPhi* (Y_real[k] + Y_real_mneg[k]));
        sumcr += -Aki * (cosPhipm3* (Y_real[k] - Y_real_mneg[k]) + sinPhipm3* (Y_imag[k] + Y_imag_mneg[k]));
        sumci += -Aki * (-cosPhipm3*(Y_imag[k] + Y_imag_mneg[k]) + sinPhipm3* (Y_real[k] - Y_real_mneg[k]));
      }
      
      hpc->preal[i] = sumpr; 
      hpc->pimag[i] = sumpi;
      hpc->creal[i] = sumcr;
      hpc->cimag[i] = sumci;	  
    }
  }
    
#endif
    
#ifdef _OPENMP
    if (USETIMERS) openmp_timer_stop("compute_hpc_FD");
#endif
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
int system_mkdir(const char *name)
{
  char s[STRLEN];
  sprintf(s,"mkdir -p %s",name);
  return system(s);
  /* if (system(s)) errorexit("Error during system call to make directory.");  */
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
