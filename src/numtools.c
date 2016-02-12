/**
 * \file numtools.c
 *
 * \brief Numerical routines.
 *
 * \author S.Bernuzzi
 */

#include "eob.h"

/* alloc Pade coefficients */
tPade *Pade_alloc (int N, int M)
{
  tPade *P = calloc(1, sizeof(tPade));
  P->N = N;
  P->M = M;
  P->n = calloc(N+1, sizeof(double));
  P->d = calloc(M+1, sizeof(double));
  if (N<1) P->d1n = calloc(1, sizeof(double));
  else     P->d1n = calloc(N, sizeof(double));
  if (M<1) P->d1d = calloc(1, sizeof(double));
  else     P->d1d = calloc(M, sizeof(double));
  if (N<2) P->d2n = calloc(1, sizeof(double));
  else     P->d2n = calloc(N-1, sizeof(double)); 
  if (M<2) P->d2d = calloc(1, sizeof(double)); 
  else     P->d2d = calloc(M-1, sizeof(double)); 
  return P;
}

/* set Pade coefficients based on PN Taylor expansion */
void Pade_PNCoefs (tPade *P, double *c)
{
  double c1,c3,c4,c5,c6,c7; /* c0 = 1, c2 = 0 */

  if ( (P->N==1) && (P->M==3) ) {

    c1 = c[1];
    c3 = c[2];
    c4 = c[3];

    P->n[0] = 1;
    P->n[1] = (pow(c1,4) + 2*c1*c3 - c4)/(pow(c1,3) + c3);

    P->d1n[0] = (pow(c1,4) + 2*c1*c3 - c4)/(pow(c1,3) + c3);

    P->d2n[0] = 0.;
    
    P->d[0] = 1;
    P->d[1] = (c1*c3 - c4)/(pow(c1,3) + c3);
    P->d[2] = -((c1*(c1*c3 - c4))/(pow(c1,3) + c3));
    P->d[3] = (-pow(c3,2) - pow(c1,2)*c4)/(pow(c1,3) + c3);
    
    P->d1d[0] = (c1*c3 - c4)/(pow(c1,3) + c3);
    P->d1d[1] = (-2*c1*(c1*c3 - c4))/(pow(c1,3) + c3);
    P->d1d[2] = (3*(-pow(c3,2) - pow(c1,2)*c4))/(pow(c1,3) + c3);

    P->d2d[0] = (-2*c1*(c1*c3 - c4))/(pow(c1,3) + c3);
    P->d2d[1] = (6*(-pow(c3,2) - pow(c1,2)*c4))/(pow(c1,3) + c3);

    return;
  }

  if ( (P->N==2) && (P->M==2) ) {

    c1 = c[1];
    c3 = c[2];
    c4 = c[3];
    
    P->n[0] = 1;
    P->n[1] = (c1*c3 - c4)/c3;
    P->n[2] = (-pow(c3,2) - pow(c1,2)*c4)/(c1*c3);
    
    P->d1n[0] = (c1*c3 - c4)/c3;
    P->d1n[1] = (2*(-pow(c3,2) - pow(c1,2)*c4))/(c1*c3);
    
    P->d2n[0] = (2*(-pow(c3,2) - pow(c1,2)*c4))/(c1*c3);
  
    P->d[0] = 1;
    P->d[1] = -(c4/c3);
    P->d[2] = -(c3/c1);
    
    P->d1d[0] = -(c4/c3);
    P->d1d[1] = (-2*c3)/c1;
    
    P->d2d[0] = (-2*c3)/c1;
    
    return;
  }

  if ( (P->N==1) && (P->M==4) ) {

    c1 = c[1];
    c3 = c[2];
    c4 = c[3];
    c5 = c[4];

    P->n[0] = 1;
    P->n[1] = (pow(c1,5) + 3*pow(c1,2)*c3 - 2*c1*c4 + c5)/(pow(c1,4) + 2*c1*c3 - c4);

    P->d1n[0] = (pow(c1,5) + 3*pow(c1,2)*c3 - 2*c1*c4 + c5)/(pow(c1,4) + 2*c1*c3 - c4);
    
    P->d2n[0] = 0.;
    
    P->d[0] = 1;
    P->d[1] = (pow(c1,2)*c3 - c1*c4 + c5)/(pow(c1,4) + 2*c1*c3 - c4);
    P->d[2] = -((c1*(pow(c1,2)*c3 - c1*c4 + c5))/(pow(c1,4) + 2*c1*c3 - c4));
    P->d[3] = (-2*c1*pow(c3,2) - pow(c1,3)*c4 + c3*c4 + pow(c1,2)*c5)/(pow(c1,4) + 2*c1*c3 - c4);
    P->d[4] = (pow(c1,2)*pow(c3,2) - 2*c1*c3*c4 + pow(c4,2) - pow(c1,3)*c5 - c3*c5)/(pow(c1,4) + 2*c1*c3 - c4);
    
    P->d1d[0] = (pow(c1,2)*c3 - c1*c4 + c5)/(pow(c1,4) + 2*c1*c3 - c4);
    P->d1d[1] = (-2*c1*(pow(c1,2)*c3 - c1*c4 + c5))/(pow(c1,4) + 2*c1*c3 - c4);
    P->d1d[2] = (3*(-2*c1*pow(c3,2) - pow(c1,3)*c4 + c3*c4 + pow(c1,2)*c5))/(pow(c1,4) + 2*c1*c3 - c4);
    P->d1d[3] = (4*(pow(c1,2)*pow(c3,2) - 2*c1*c3*c4 + pow(c4,2) - pow(c1,3)*c5 - c3*c5))/
      (pow(c1,4) + 2*c1*c3 - c4);
    
    P->d2d[0] = (-2*c1*(pow(c1,2)*c3 - c1*c4 + c5))/(pow(c1,4) + 2*c1*c3 - c4);
    P->d2d[1] =  (6*(-2*c1*pow(c3,2) - pow(c1,3)*c4 + c3*c4 + pow(c1,2)*c5))/(pow(c1,4) + 2*c1*c3 - c4);
    P->d2d[2] = (12*(pow(c1,2)*pow(c3,2) - 2*c1*c3*c4 + pow(c4,2) - pow(c1,3)*c5 - c3*c5))/
      (pow(c1,4) + 2*c1*c3 - c4);
    
    return;
  }

  if ( (P->N==1) && (P->M==5) ) {

    c1 = c[1];
    c3 = c[2];
    c4 = c[3];
    c5 = c[4];
    c6 = c[5];
    
    P->n[0] = 1;
    P->n[1] = (pow(c1,6) + 4*pow(c1,3)*c3 + pow(c3,2) - 3*pow(c1,2)*c4 + 2*c1*c5 - c6)/
      (pow(c1,5) + 3*pow(c1,2)*c3 - 2*c1*c4 + c5);

    P->d1n[0] = (pow(c1,6) + 4*pow(c1,3)*c3 + pow(c3,2) - 3*pow(c1,2)*c4 + 2*c1*c5 - c6)/
      (pow(c1,5) + 3*pow(c1,2)*c3 - 2*c1*c4 + c5);
    
    P->d2n[0] = 0.;
    
    P->d[0] = 1;
    P->d[1] = (pow(c1,3)*c3 + pow(c3,2) - pow(c1,2)*c4 + c1*c5 - c6)/
      (pow(c1,5) + 3*pow(c1,2)*c3 - 2*c1*c4 + c5);
    P->d[2] = -((c1*(pow(c1,3)*c3 + pow(c3,2) - pow(c1,2)*c4 + c1*c5 - c6))/
		(pow(c1,5) + 3*pow(c1,2)*c3 - 2*c1*c4 + c5));
    P->d[3] = (-2*pow(c1,2)*pow(c3,2) - pow(c1,4)*c4 + 2*c1*c3*c4 + pow(c1,3)*c5 - c3*c5 - pow(c1,2)*c6)/
      (pow(c1,5) + 3*pow(c1,2)*c3 - 2*c1*c4 + c5);
    P->d[4] = (pow(c1,3)*pow(c3,2) - pow(c3,3) - 4*pow(c1,2)*c3*c4 + 2*c1*pow(c4,2) - pow(c1,4)*c5 - c4*c5 + 
	       pow(c1,3)*c6 + c3*c6)/(pow(c1,5) + 3*pow(c1,2)*c3 - 2*c1*c4 + c5);
    P->d[5] = (2*c1*pow(c3,3) + 2*pow(c1,3)*c3*c4 - pow(c3,2)*c4 - pow(c1,2)*pow(c4,2) - 2*pow(c1,2)*c3*c5 + 
	       2*c1*c4*c5 - pow(c5,2) - pow(c1,4)*c6 - 2*c1*c3*c6 + c4*c6)/
      (pow(c1,5) + 3*pow(c1,2)*c3 - 2*c1*c4 + c5);

    P->d1d[0] = (pow(c1,3)*c3 + pow(c3,2) - pow(c1,2)*c4 + c1*c5 - c6)/
      (pow(c1,5) + 3*pow(c1,2)*c3 - 2*c1*c4 + c5);
    P->d1d[1] = (-2*c1*(pow(c1,3)*c3 + pow(c3,2) - pow(c1,2)*c4 + c1*c5 - c6))/
      (pow(c1,5) + 3*pow(c1,2)*c3 - 2*c1*c4 + c5);
    P->d1d[2] = (3*(-2*pow(c1,2)*pow(c3,2) - pow(c1,4)*c4 + 2*c1*c3*c4 + pow(c1,3)*c5 - c3*c5 - pow(c1,2)*c6))/
      (pow(c1,5) + 3*pow(c1,2)*c3 - 2*c1*c4 + c5);
    P->d1d[3] = (4*(pow(c1,3)*pow(c3,2) - pow(c3,3) - 4*pow(c1,2)*c3*c4 + 2*c1*pow(c4,2) - pow(c1,4)*c5 - c4*c5 + 
		    pow(c1,3)*c6 + c3*c6))/(pow(c1,5) + 3*pow(c1,2)*c3 - 2*c1*c4 + c5);
    P->d1d[4] = (5*(2*c1*pow(c3,3) + 2*pow(c1,3)*c3*c4 - pow(c3,2)*c4 - pow(c1,2)*pow(c4,2) - 2*pow(c1,2)*c3*c5 + 
		    2*c1*c4*c5 - pow(c5,2) - pow(c1,4)*c6 - 2*c1*c3*c6 + c4*c6))/
      (pow(c1,5) + 3*pow(c1,2)*c3 - 2*c1*c4 + c5);
    
    P->d2d[0] = (-2*c1*(pow(c1,3)*c3 + pow(c3,2) - pow(c1,2)*c4 + c1*c5 - c6))/
      (pow(c1,5) + 3*pow(c1,2)*c3 - 2*c1*c4 + c5); 
    P->d2d[1] = (6*(-2*pow(c1,2)*pow(c3,2) - pow(c1,4)*c4 + 2*c1*c3*c4 + pow(c1,3)*c5 - c3*c5 - pow(c1,2)*c6))/
      (pow(c1,5) + 3*pow(c1,2)*c3 - 2*c1*c4 + c5);
    P->d2d[2] = (12*(pow(c1,3)*pow(c3,2) - pow(c3,3) - 4*pow(c1,2)*c3*c4 + 2*c1*pow(c4,2) - pow(c1,4)*c5 - 
		     c4*c5 + pow(c1,3)*c6 + c3*c6))/(pow(c1,5) + 3*pow(c1,2)*c3 - 2*c1*c4 + c5);
    P->d2d[3] = (20*(2*c1*pow(c3,3) + 2*pow(c1,3)*c3*c4 - pow(c3,2)*c4 - pow(c1,2)*pow(c4,2) - 
		     2*pow(c1,2)*c3*c5 + 2*c1*c4*c5 - pow(c5,2) - pow(c1,4)*c6 - 2*c1*c3*c6 + c4*c6))/
      (pow(c1,5) + 3*pow(c1,2)*c3 - 2*c1*c4 + c5);
    
    return;
  }

}

/* free Pade Coefficients */
void Pade_free (tPade *P)
{
  if (P->n) free(P->n);
  if (P->d) free(P->d);
  if (P->d1n) free(P->d1n);
  if (P->d1d) free(P->d1d);
  if (P->d2n) free(P->d2n);
  if (P->d2d) free(P->d2d);
  free(P);
}

/* Euler log */
double Eulerlog (double x, int m)
{
  static double EG = EulerGamma;
  return ( EG + log(2.*m*sqrt(x)) );
}

/* Gamma function for real argument */
double GammaReal (double z)
{

#define nc 15 /* elements of c */
#define SMALL 1e-13
#define INFINITY 666e+300

  /* Lanczos approximation for the complex plane calculated using vpa
     digits(256) the best set of coeffs was selected from a solution space of
     g=0 to 32 with 1 to 32 terms these coeffs really give superb performance
     of 15 sig. digits for 0<=real(z)<=171 coeffs should sum to about
     g*g/2+23/24 */
  static const double g = 607./128.; /* best results when 4<=g<=5 */
  static double c[] = { 0.99999999999999709182,
			57.156235665862923517,
			-59.597960355475491248,
			14.136097974741747174,
			-0.49191381609762019978,
			.33994649984811888699e-4,
			.46523628927048575665e-4,
			-.98374475304879564677e-4,
			.15808870322491248884e-3,
			-.21026444172410488319e-3,
			.21743961811521264320e-3,
			-.16431810653676389022e-3,
			.84418223983852743293e-4,
			-.26190838401581408670e-4,
			.36899182659531622704e-5 };
  
  double zz = z;
  double zh;
  double zgh;
  double zp;
  double ss;
  double Gam;
  double zr;
  
  int pp;

  /* Some checks */
  int iszneg  = 0;    
  int iszzero = 0;    
  int iszone  = 0;
  
  if (z<0.) iszneg = 1; 
  if (DEQUAL(z,0.,SMALL)) iszzero = 1; 
  if (DEQUAL(z,1.,SMALL)) iszone = 1;     

  /* Special case: z = 0,1 */ 
  if ( (iszzero) || (iszone) ) 
    return 1.;
  
  /* Flip sign ? */
  if (iszneg) z = - z; 
    
  /* Special case: negative poles */ 
  if ( ((iszneg) || (iszzero)) && 
       (DEQUAL(round(z),z,SMALL)) ) { 
    return INFINITY;
  }
  
  /* aux quantities */ 
  z   = z - 1;
  zh  = z + 0.5;
  zgh = zh + g;
  
  /* trick for avoiding FP overflow above z=141 */ 
  zp = pow( zgh, zh*0.5 );
  
  /* Sum  pp = 14 ... 1 */ 
  ss = 0.;
  for( pp = nc-1; pp>=1; pp-- ) 
    ss += c[pp]/(z+pp); 
  
  /* Result */ 
  Gam = (sq2Pi*(c[0]+ss)) * ((zp * exp(-zgh)) * zp); 
  
  /* Adjust for negative real parts */
  if (iszneg) 
    Gam = -Pi/(zz * Gam * sin(Pi*zz)); 
  
  return Gam; 
}

/* Gamma function for complex argument */
domplex GammaComplex (domplex z)
{

#define nc 15 /* elements of c */
#define SMALL 1e-13
#define INFINITY 666e+300
  
  /* Lanczos approximation for the complex plane calculated using vpa
     digits(256) the best set of coeffs was selected from a solution space of
     g=0 to 32 with 1 to 32 terms these coeffs really give superb performance
     of 15 sig. digits for 0<=real(z)<=171 coeffs should sum to about
     g*g/2+23/24 */
  static const double g = 607./128.; /* best results when 4<=g<=5 */
  static double c[] = { 0.99999999999999709182,
			57.156235665862923517,
			-59.597960355475491248,
			14.136097974741747174,
			-0.49191381609762019978,
			.33994649984811888699e-4,
			.46523628927048575665e-4,
			-.98374475304879564677e-4,
			.15808870322491248884e-3,
			-.21026444172410488319e-3,
			.21743961811521264320e-3,
			-.16431810653676389022e-3,
			.84418223983852743293e-4,
			-.26190838401581408670e-4,
			.36899182659531622704e-5 };
  
  double zr = creal(z);
  double zi = cimag(z);
  
  if (DEQUAL(zi,0.,SMALL))  {
    /* entry is real */
    return GammaReal (zr);
  }
  
  domplex zz = z;
  domplex zh;
  domplex zgh;
  domplex zp;
  domplex ss;
  domplex Gam;

  int iszrneg = 0;
  int pp;
  
  /* Is z<0 ? => flip sign */ 
  if (zr<0.) {
    iszrneg = 1; 
    z = - z;
  } 

  /* aux quantities */ 
  z   = z - 1;
  zh  = z + 0.5;
  zgh = zh + g;

  /* trick for avoiding FP overflow above z=141 */ 
  zp = cpow( zgh, zh*0.5 );

  /* Sum  pp = 14 ... 1 */ 
  ss = 0.;
  for( pp = nc-1; pp>=1; pp-- ) 
    ss += c[pp]/(z+pp); 

  /* Result */ 
  Gam = (Sq2Pi*(c[0]+ss)) * ((zp * cexp(-zgh)) * zp);   

  /* Adjust for negative real parts */
  if (iszrneg) 
    Gam = -Pi/(zz * Gam * csin(Pi*zz)); 
 
  return Gam; 
}

/* factorial */
double ffact (double n)
{
  static double f[] = {1,         1,          2,     
		       6,         24,         120,
		       720,       5040,       40320, 
		       362880,    3628800,    39916800,   
		       479001600, 6227020800, 87178291200};
  if(n<=14.) return(f[(int)n]);
  else       return(fact(n));
}

double fact (double n)
{
  if(n < 0.) {
    errorexit("not possible computing a negative factorial");
  }
  if (n<=14.) {
    return(ffact(n));
  } else {
    n = n*fact(n-1);
    return(n);
  }
}

/* Find nearest point index in 1d array */
int find_point_bisection(double x, int n, double *xp, int o)
{
  int i0 = o-1, i1 = n-o;
  int i;

  if (n < 2*o) {
    errorexit(" not enough point to interpolate");
  }
  
  if (x <= xp[i0]) return 0;
  if (x >  xp[i1]) return n-2*o;

  while (i0 != i1-1) {
    i = (i0+i1)/2;
    if (x < xp[i]) i1 = i; else i0 = i;
  }

  return i0-o+1;
}

/* Barycentric Lagrange interpolation at xx with n points of f(x), 
   equivalent to standard Lagrangian interpolation */   
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

/* Barycentric Lagrange interpolation at xx with n points of f(x), 
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

/* Barycentric Lagrange interpolation at xx with n points of f(x), 
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

/* 1d Lagrangian barycentric interpolation */
double interp1d (const int order, double xx, int nx, double *f, double *x)
{
  double ff;
  int ix;
  int ox = order > nx ? nx : order;
  ix = find_point_bisection(xx, nx, x, ox/2);
  ff = baryc_f(xx, ox, &f[ix], &x[ix]);  
  return( ff );
}


