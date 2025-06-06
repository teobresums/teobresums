/** \file TEOBResumSUtils.c
 *  \brief Useful functions for TEOBResumS
 * 
 *  This file contains useful functions for TEOBResumS, as well as 
 *  the functions to alloc/free/push/interpolate many of the structures
 *  employed in the code.
 *
 */

#include "TEOBResumS.h"

/**
 * Function: q_to_nu
 * -----------------
 * Return symm mass ratio from q
 * 
 * @param[in] q: mass ratio
 * 
 * @return nu: symm mass ratio 
 */
double q_to_nu(const double q)
{
  double nu = 0;
  if (q>0.)
    nu = q/((q+1.)*(q+1.));
  if (nu > 0.25)
    nu = 0.25;
  
  return nu;
}

/**
 * Function: nu_to_X1
 * -----------------
 * Return mass fraction M1/M from nu
 * 
 * @param[in] nu: symm mass ratio
 * 
 * @return X1: mass ratio 
 */
double nu_to_X1(const double nu)
{
  if ( (nu<0.) || (nu>0.25) )
    errorexit("sym mass ratio is defined 0<= nu <= 1/4.");
  return 0.5*(1.+sqrt(1.-4.*nu));
}

/**
 * Function: tidal_kappa_of_Lambda
 * -------------------------------
 * Return tidal coupling constants from tidal polarizability parameters
 * 
 * @param[in] q    : mass ratio
 * @param[in] XA   : mass fraction M1/M
 * @param[in] XB   : mass fraction M2/M
 * @param[in] LamA : tidal polarizability parameter of star A
 * @param[in] LamB : tidal polarizability parameter of star B
 * @param[in] ell  : multipole index
 * @param[out] kapA: tidal coupling constant of star A
 * @param[out] kapB: tidal coupling constant of star B
*/
double tidal_kappa_of_Lambda(double q, double XA, double XB, double LamA, double LamB, int ell,
			     double *kapA, double *kapB)
{
  if (ell<2) errorexit("ell must be >2");
  double f2lm1 = doublefact(2*ell-1);
  int p = 2*ell + 1;
  *kapA = f2lm1 * LamA * pow(XA,p) / q;
  *kapB = f2lm1 * LamB * pow(XB,p) * q;
}

/**
  * Function: set_spin_vars
  * -----------------------
  * Set spin variables
  * 
  * @param[in] X1   : mass fraction M1/M
  * @param[in] X2   : mass fraction M2/M
  * @param[in] chi1 : dimensionless spin of body 1
  * @param[in] chi2 : dimensionless spin of body 2
  * @param[out] S1  : spin of body 1
  * @param[out] S2  : spin of body 2
  * @param[out] a1  : X1*chi1
  * @param[out] a2  : X2*chi2
  * @param[out] aK  : X1*chi1 + X2*chi2
  * @param[out] aK2 : (X1*chi1 + X2*chi2)^2
  * @param[out] S   : S1 + S2
  * @param[out] Sstar : X2*a1 + X1*a2
*/
void set_spin_vars(double X1, double X2, double chi1, double chi2,
		   double *S1, double *S2, double *a1, double *a2,
		   double *aK, double *aK2,
		   double *S, double *Sstar)
{
  *S1 = SQ(X1) * chi1;
  *S2 = SQ(X2) * chi2;
  *a1 = X1 * chi1;
  *a2 = X2 * chi2;
  *aK = (*a1) +  (*a2);
  *aK2 = SQ((*aK));   
  *S = (*S1) +  (*S2);           
  *Sstar = X2 * (*a1) + X1 * (*a2);
}

/** Eulerlog function (constants are defined in header) */
static const double Logm[] = {0.,Log1,Log2,Log3,Log4,Log5,Log6,Log7};
/**
 * Function: Eulerlog
 * ------------------
 * Return Eulerlog function
 * 
 * @param[in] x: argument
 * @param[in] m: integer
 * 
 * @return Eulerlog(x,m): Eulerlog function
*/
double Eulerlog(const double x,const int m)
{
  double logm = 0.;
  if ((m>0) & (m<8)) logm = Logm[m];
  else logm = log((double)m);
  return EulerGamma_Log2 + logm + 0.5*log(x);
}

/**
 * Function: Pade02
 * ----------------
 * Return Pade02 function
 * 
 * @param[in] x: argument
 * @param[in] a: array of coefficients
 * 
 * @return Pade02(x,a): Pade02 function
*/
double Pade02(double x, double *a){
  double x2  = x*x;
  double d1  = -a[1];
  double d2  = SQ(a[1]) - a[2];
  double den = (1. + d1*x + d2*x2);
  return 1./den;
}

/**
 * Function: Pade22
 * ----------------
 * Return Pade22 function
 * 
 * @param[in] x: argument
 * @param[in] a: array of coefficients
 * 
 * @return Pade22(x,a): Pade22 function
*/
double Pade22(double x, double *a){
  double x2 = x*x;

  double cden =  a[1]*a[3]-SQ(a[2]);
  double n1   = -a[1]*SQ(a[2]) + SQ(a[1])*a[3]+a[2]*a[3] - a[1]*a[4];
  double n2   = -SQ(a[2])*a[2]+2*a[1]*a[2]*a[3]-SQ(a[3])-SQ(a[1])*a[4]+a[2]*a[4];
  double d1   =  a[2]*a[3] - a[1]*a[4];
  double d2   = -SQ(a[3]) + a[2]*a[4];

  double pade = (cden + n1*x + n2*x2)/(cden + d1*x + d2*x2);
  
  return pade;
}

/**
 * Function: Pade32
 * ----------------
 * Return Pade32 function
 * 
 * @param[in] x: argument
 * @param[in] a: array of coefficients
 * 
 * @return Pade32(x,a): Pade32 function
*/
double Pade32(double x, double *a){
  double x2 = x*x;
  double x3 = x2*x;

  double cden = SQ(a[3]) - a[2]*a[4];

  double n1 = a[1]*SQ(a[3]) - a[1]*a[2]*a[4] - a[3]*a[4] + a[2]*a[5];
  double n2 = a[2]*SQ(a[3]) - SQ(a[2])*a[4] - a[1]*a[3]*a[4] + SQ(a[4]) + a[1]*a[2]*a[5] - a[3]*a[5];
  double n3 = SQ(a[3])*a[3] - 2.*a[2]*a[3]*a[4] + a[1]*SQ(a[4]) + SQ(a[2])*a[5] - a[1]*a[3]*a[5];

  double d1 = - a[3]*a[4] + a[2]*a[5];
  double d2 = SQ(a[4]) - a[3]*a[5];

  double pade = (cden + n1*x + n2*x2 + n3*x3)/(cden + d1*x + d2*x2);

  return pade;
}

/**
 * Function: Pade23
 * ----------------
 * Return Pade23 function
 * 
 * @param[in] x: argument
 * @param[in] a: array of coefficients
 * 
 * @return Pade23(x,a): Pade23 function
*/
double Pade23(double x, double *a){
  double x2 = x*x;
  double x3 = x2*x;

  double cden = SQ(a[2])*a[2] - 2.*a[1]*a[2]*a[3] + SQ(a[3]) + SQ(a[1])*a[4] - a[2]*a[4];

  double n1 = a[1]*SQ(a[2])*a[2] - 2.*SQ(a[1])*a[2]*a[3] - SQ(a[2])*a[3] + 2.*a[1]*SQ(a[3]) + SQ(a[1])*a[1]*a[4] - a[3]*a[4] - SQ(a[1])*a[5] + a[2]*a[5];
  double n2 = SQ(SQ(a[2])) - 3.*a[1]*SQ(a[2])*a[3] + SQ(a[1])*SQ(a[3]) + 2.*a[2]*SQ(a[3]) + 2.*SQ(a[1])*a[2]*a[4] - 2.*SQ(a[2])*a[4] - 2.*a[1]*a[3]*a[4] + SQ(a[4]) - SQ(a[1])*a[1]*a[5] + 2.*a[1]*a[2]*a[5] - a[3]*a[5];

  double d1 = -SQ(a[2])*a[3] + a[1]*SQ(a[3]) + a[1]*a[2]*a[4] - a[3]*a[4] - SQ(a[1])*a[5] + a[2]*a[5];
  double d2 = a[2]*SQ(a[3]) - SQ(a[2])*a[4] - a[1]*a[3]*a[4] + SQ(a[4]) + a[1]*a[2]*a[5] - a[3]*a[5];
  double d3 = -SQ(a[3])*a[3] + 2.*a[2]*a[3]*a[4] - a[1]*SQ(a[4]) - SQ(a[2])*a[5] + a[1]*a[3]*a[5];
    
  double pade = (cden + n1*x + n2*x2)/(cden + d1*x + d2*x2 +d3*x3);

  return pade;
}

/**
 * Function: Pade51
 * ----------------
 * Return Pade51 function
 * 
 * @param[in] x: argument
 * @param[in] a: array of coefficients
 * 
 * @return Pade51(x,a): Pade51 function
 */
double Pade51(double x, double *a){
  double x2 = x*x;
  double x3 = x2*x;
  double x4 = x3*x;
  double x5 = x4*x;

  double cden = a[5];

  double n1 = a[1]*a[5] -      a[6];
  double n2 = a[2]*a[5] - a[1]*a[6];
  double n3 = a[3]*a[5] - a[2]*a[6];
  double n4 = a[4]*a[5] - a[3]*a[6];
  double n5 = a[5]*a[5] - a[4]*a[6];

  double d1 = - a[6];

  double pade = (cden + n1*x + n2*x2 + n3*x3 + n4*x4 + n5*x5)/(cden + d1*x);

  return pade;
}

/**
 * Function: Pade42
 * ----------------
 * Return Pade42 function
 * 
 * @param[in] x: argument
 * @param[in] a: array of coefficients
 * 
 * @return Pade42(x,a): Pade42 function
 */
double Pade42(double x, double *a){
  double x2 = x*x;
  double x3 = x2*x;
  double x4 = x3*x;

  double cden = SQ(a[4]) - a[3]*a[5];

  double n1 = a[1]*SQ(a[4]) - a[1]*a[3]*a[5] - a[4]*a[5] + a[3]*a[6];
  double n2 = a[2]*SQ(a[4]) - a[2]*a[3]*a[5] - a[1]*a[4]*a[5] + SQ(a[5]) + a[1]*a[3]*a[6] - a[4]*a[6];
  double n3 = a[3]*SQ(a[4]) - SQ(a[3])*a[5] - a[2]*a[4]*a[5] + a[1]*SQ(a[5]) + a[2]*a[3]*a[6] - a[1]*a[4]*a[6];
  double n4 = SQ(a[4])*a[4] - 2*a[3]*a[4]*a[5] + a[2]*SQ(a[5]) + SQ(a[3])*a[6] - a[2]*a[4]*a[6];
  
  double d1 = - a[4]*a[5] + a[3]*a[6];
  double d2 = SQ(a[5]) - a[4]*a[6];

  double pade = (cden + n1*x + n2*x2 + n3*x3 + n4*x4)/(cden + d1*x + d2*x2);

  return pade;
}

/**
 * Function: Pade33
 * ----------------
 * Return Pade33 function
 * 
 * @param[in] x: argument
 * @param[in] a: array of coefficients
 * 
 * @return Pade33(x,a): Pade33 function
 */
double Pade33(double x, double *a){
  double x2 = x*x;
  double x3 = x2*x;

  double cden = SQ(a[3])*a[3] - 2.*a[2]*a[3]*a[4] + a[1]*SQ(a[4]) + SQ(a[2])*a[5] - a[1]*a[3]*a[5];

  double n1 = a[1]*a[3]*SQ(a[3]) - 2.*a[1]*a[2]*a[3]*a[4] - SQ(a[3])*a[4] + SQ(a[1])*SQ(a[4]) + a[2]*SQ(a[4]) + a[1]*SQ(a[2])*a[5] - SQ(a[1])*a[3]*a[5] + a[2]*a[3]*a[5] - a[1]*a[4]*a[5] - SQ(a[2])*a[6] + a[1]*a[3]*a[6];
  double n2 = a[2]*a[3]*SQ(a[3]) - 2.*SQ(a[2])*a[3]*a[4] - a[1]*SQ(a[3])*a[4] + 2.*a[1]*a[2]*SQ(a[4]) + a[3]*SQ(a[4]) + a[2]*SQ(a[2])*a[5] - SQ(a[3])*a[5] - SQ(a[1])*a[4]*a[5] - a[2]*a[4]*a[5] + a[1]*SQ(a[5]) - a[1]*SQ(a[2])*a[6] + SQ(a[1])*a[3]*a[6] + a[2]*a[3]*a[6] - a[1]*a[4]*a[6];
  double n3 = SQ(SQ(a[3])) -3.*a[2]*SQ(a[3])*a[4] + SQ(a[2])*SQ(a[4]) + 2.*a[1]*a[3]*SQ(a[4]) - a[4]*SQ(a[4]) + 2.*SQ(a[2])*a[3]*a[5] - 2.*a[1]*SQ(a[3])*a[5] - 2.*a[1]*a[2]*a[4]*a[5] + 2.*a[3]*a[4]*a[5] + SQ(a[1])*SQ(a[5]) - a[2]*SQ(a[5]) - a[2]*SQ(a[2])*a[6] + 2.*a[1]*a[2]*a[3]*a[6] - SQ(a[3])*a[6] - SQ(a[1])*a[4]*a[6] + a[2]*a[4]*a[6];
  
  double d1 = -SQ(a[3])*a[4] + a[2]*SQ(a[4]) + a[2]*a[3]*a[5] - a[1]*a[4]*a[5] - SQ(a[2])*a[6] + a[1]*a[3]*a[6];
  double d2 = a[3]*SQ(a[4]) - SQ(a[3])*a[5] - a[2]*a[4]*a[5] + a[1]*SQ(a[5]) + a[2]*a[3]*a[6] - a[1]*a[4]*a[6];
  double d3 = -SQ(a[4])*a[4] + 2.*a[3]*a[4]*a[5] - a[2]*SQ(a[5]) - SQ(a[3])*a[6] + a[2]*a[4]*a[6];

  double pade = (cden + n1*x + n2*x2 + n3*x3)/(cden + d1*x + d2*x2 + d3*x3);

  return pade;
}

/**
 * Function: Pade15
 * ----------------
 * Return Pade15 function
 * 
 * @param[in] x: argument
 * @param[in] a: array of coefficients
 * 
 * @return Pade15(x,a): Pade15 function
 */
double Pade15(double x, double *a){
  double x2 = x*x;
  double x3 = x2*x;
  double x4 = x3*x;
  double x5 = x4*x;

  double cden = SQ(SQ(a[1]))*a[1] - 4.*SQ(a[1])*a[1]*a[2] + 3.*a[1]*SQ(a[2]) + 3.*SQ(a[1])*a[3] - 2.*a[2]*a[3] - 2.*a[1]*a[4] + a[5];

  double n1 = SQ(SQ(a[1]))*SQ(a[1]) - 5.*SQ(SQ(a[1]))*a[2] + 6.*SQ(a[1])*SQ(a[2]) - SQ(a[2])*a[2] + 4.*SQ(a[1])*a[1]*a[3] - 6.*a[1]*a[2]*a[3] + SQ(a[3]) - 3.*SQ(a[1])*a[4] + 2.*a[2]*a[4] + 2.*a[1]*a[5] - a[6];

  double d1 = -SQ(SQ(a[1]))*a[2] + 3.*SQ(a[1])*SQ(a[2]) - SQ(a[2])*a[2] + SQ(a[1])*a[1]*a[3] - 4.*a[1]*a[2]*a[3] + SQ(a[3]) - SQ(a[1])*a[4] + 2.*a[2]*a[4] + a[1]*a[5] - a[6];  
  double d2 = SQ(a[1])*a[1]*SQ(a[2]) - 2.*a[1]*SQ(a[2])*a[2] - SQ(SQ(a[1]))*a[3] + SQ(a[1])*a[2]*a[3] + 2.*SQ(a[2])*a[3] - a[1]*SQ(a[3]) + SQ(a[1])*a[1]*a[4] - SQ(a[1])*a[5] - a[2]*a[5] + a[1]*a[6];  
  double d3 = -SQ(a[1])*SQ(a[2])*a[2] + SQ(SQ(a[2])) + 2.*SQ(a[1])*a[1]*a[2]*a[3] - a[1]*SQ(a[2])*a[3] - 2.*SQ(a[1])*SQ(a[3]) +  a[2]*SQ(a[3]) - SQ(SQ(a[1]))*a[4] + SQ(a[1])*a[2]*a[4] - 2.*SQ(a[2])*a[4] + 2.*a[1]*a[3]*a[4] + SQ(a[1])*a[1]*a[5] - a[3]*a[5] - SQ(a[1])*a[6] + a[2]*a[6];  
  double d4 = a[1]*SQ(SQ(a[2])) - 3.*SQ(a[1])*SQ(a[2])*a[3] - SQ(a[2])*a[2]*a[3] + SQ(a[1])*a[1]*SQ(a[3]) + 4.*a[1]*a[2]*SQ(a[3]) - SQ(a[3])*a[3] + 2.*SQ(a[1])*a[1]*a[2]*a[4] - a[1]*SQ(a[2])*a[4] - 4.*SQ(a[1])*a[3]*a[4] + 2.*a[1]*SQ(a[4]) - SQ(SQ(a[1]))*a[5] + SQ(a[1])*a[2]*a[5] + SQ(a[2])*a[5] - a[4]*a[5] + SQ(a[1])*a[1]*a[6] - 2.*a[1]*a[2]*a[6] + a[3]*a[6];  
  double d5 = -SQ(SQ(a[2]))*a[2] + 4.*a[1]*SQ(a[2])*a[2]*a[3] - 3.*SQ(a[1])*a[2]*SQ(a[3]) - 3.*SQ(a[2])*SQ(a[3]) + 2.*a[1]*SQ(a[3])*a[3] - 3.*SQ(a[1])*SQ(a[2])*a[4] + 3.*SQ(a[2])*a[2]*a[4] + 2.*SQ(a[1])*a[1]*a[3]*a[4] + 2.*a[1]*a[2]*a[3]*a[4] - SQ(a[3])*a[4] - SQ(a[1])*SQ(a[4]) - 2.*a[2]*SQ(a[4]) + 2.*SQ(a[1])*a[1]*a[2]*a[5] - 4.*a[1]*SQ(a[2])*a[5] - 2.*SQ(a[1])*a[3]*a[5] + 4.*a[2]*a[3]*a[5] + 2.*a[1]*a[4]*a[5] - SQ(a[5]) - SQ(SQ(a[1]))*a[6] + 3.*SQ(a[1])*a[2]*a[6] - SQ(a[2])*a[6] - 2.*a[1]*a[3]*a[6] + a[4]*a[6];

  double pade = (cden + n1*x)/(cden + d1*x + d2*x2 + d3*x3 + d4*x4 + d5*x5);

  return pade;
}

/**
 * Function: Pade62
 * ----------------
 * Return Pade62 function
 * 
 * @param[in] x: argument
 * @param[in] a: array of coefficients
 * 
 * @return Pade62(x,a): Pade62 function
 */
double Pade62(double x, double *a){
  double x2 = x*x;
  double x3 = x2*x;
  double x4 = x3*x;	
  double x5 = x4*x;
  double x6 = x5*x;

  double cden = SQ(a[6]) - a[5]*a[7];

  double n1 = a[1]*SQ(a[6]) - a[1]*a[5]*a[7] - a[6]*a[7] + a[5]*a[8];
  double n2 = a[2]*SQ(a[6]) - a[2]*a[5]*a[7] - a[1]*a[6]*a[7] + SQ(a[7]) + a[1]*a[5]*a[8] - a[6]*a[8];
  double n3 = a[3]*SQ(a[6]) - a[3]*a[5]*a[7] - a[2]*a[6]*a[7] + a[1]*SQ(a[7]) + a[2]*a[5]*a[8] - a[1]*a[6]*a[8];
  double n4 = a[4]*SQ(a[6]) - a[4]*a[5]*a[7] - a[3]*a[6]*a[7] + a[2]*SQ(a[7]) + a[3]*a[5]*a[8] - a[2]*a[6]*a[8];
  double n5 = a[5]*SQ(a[6]) - SQ(a[5])*a[7] - a[4]*a[6]*a[7] + a[3]*SQ(a[7]) + a[4]*a[5]*a[8] - a[3]*a[6]*a[8];
  double n6 = a[6]*SQ(a[6]) - 2.*a[5]*a[6]*a[7] + a[4]*SQ(a[7]) + SQ(a[5])*a[8] - a[4]*a[6]*a[8];

  double d1 = -a[6]*a[7] + a[5]*a[8];
  double d2 = SQ(a[7]) - a[6]*a[8];

  double pade = (cden + n1*x + n2*x2 + n3*x3 + n4*x4 + n5*x5 + n6*x6)/(cden + d1*x + d2*x2);

  return pade;
}

/**
 * Function: Pade33_forGSF
 * ----------------
 * Return Pade33_forGSF function
 * 
 * @param[in] coeffs: coefficients
 * @param[in] Dcoeffs: derivatives of coefficients
 * @param[in] D2coeffs: second derivatives of coefficients
 * @param[in] u: argument
 * @param[out] P: Pade33_forGSF function
 * @param[out] DP: derivative of Pade33_forGSF function
 * @param[out] D2P: second derivative of Pade33_forGSF function
 * 
 */
void Pade33_forGSF(double coeffs[6], double Dcoeffs[5], double D2coeffs[5], double u, double *P, double *DP, double *D2P)
{
  // Padé 3,3 of f[u] = 1 + c1 u + c2[u] u^2 + c3[u] u^3 + c4[u] u^4 + c5[u] u^5 + c6[u] u^6 
  // (functional form of the integer part of a1SF, d1SF, q1SF)
  // + notice that c1 does *not* depend on u so it has no derivatives.

  double u2 = u*u;
  double u3 = u2*u;
  double u4 = u3*u;

  // Coefficients of the Taylor-expanded function + derivatives (except from c1, they are functions of log(u))
  double c1 = coeffs[0];
  double c2 = coeffs[1];
  double c3 = coeffs[2];
  double c4 = coeffs[3];
  double c5 = coeffs[4];
  double c6 = coeffs[5];

  double dc2 = Dcoeffs[0];
  double dc3 = Dcoeffs[1];
  double dc4 = Dcoeffs[2];
  double dc5 = Dcoeffs[3];
  double dc6 = Dcoeffs[4];

  double d2c2 = D2coeffs[0];
  double d2c3 = D2coeffs[1];
  double d2c4 = D2coeffs[2];
  double d2c5 = D2coeffs[3];
  double d2c6 = D2coeffs[4];

  // Powers of the coefficients
  double c1to2 = c1*c1;
  double c2to2 = c2*c2;
  double c2to3 = c2to2*c2;
  double c3to2 = c3*c3;
  double c3to3 = c3to2*c3;
  double c3to4 = c3to3*c3;
  double c4to2 = c4*c4;
  double c4to3 = c4to2*c4;
  double c5to2 = c5*c5;
  // Powers of the derivatives
  double dc2to2 = dc2*dc2;
  double dc3to2 = dc3*dc3;
  double dc4to2 = dc4*dc4;
  double dc5to2 = dc5*dc5;

  // Coefficients of the numerator + derivatives
  double n0   = c3to3 - 2.*c2*c3*c4 + c1*c4to2 + c2to2*c5 - c1*c3*c5;
  double Dn0  = - 2.*c3*c4*dc2 + 2.*c2*c5*dc2 + 3.*c3to2*dc3 - 2.*c2*c4*dc3 - c1*c5*dc3 - 2.*c2*c3*dc4 + 2.*c1*c4*dc4 
                + c2to2*dc5 - c1*c3*dc5;
  double D2n0 = - 2.*c3*c4*d2c2 + 2.*c2*c5*d2c2 + 3.*c3to2*d2c3 - 2.*c2*c4*d2c3 - c1*c5*d2c3 - 2.*c2*c3*d2c4 + 2.*c1*c4*d2c4 
                + c2to2*d2c5 - c1*c3*d2c5 + 2.*c5*dc2to2 - 4.*c4*dc2*dc3 + 6.*c3*dc3to2 - 4.*c3*dc2*dc4 - 4.*c2*dc3*dc4 
                + 2.*c1*dc4to2 + 4.*c2*dc2*dc5 - 2.*c1*dc3*dc5;
  double n1   = c1*c3to3 - 2.*c1*c2*c3*c4 - c3to2*c4 + c1to2*c4to2 + c2*c4to2 + c1*c2to2*c5 - c1to2*c3*c5 + c2*c3*c5 
                - c1*c4*c5 - c2to2*c6 + c1*c3*c6;
  double Dn1  = - 2.*c1*c3*c4*dc2 + c4to2*dc2 + 2.*c1*c2*c5*dc2 + c3*c5*dc2 - 2.*c2*c6*dc2 + 3.*c1*c3to2*dc3 - 2.*c1*c2*c4*dc3 
                - 2.*c3*c4*dc3 - c1to2*c5*dc3 + c2*c5*dc3 + c1*c6*dc3 - 2.*c1*c2*c3*dc4 - c3to2*dc4 + 2.*c1to2*c4*dc4 
                + 2.*c2*c4*dc4 - c1*c5*dc4 + c1*c2to2*dc5 - c1to2*c3*dc5 + c2*c3*dc5 - c1*c4*dc5 - c2to2*dc6 + c1*c3*dc6;
  double D2n1 = - 2.*c1*c3*c4*d2c2 + c4to2*d2c2 + 2.*c1*c2*c5*d2c2 + c3*c5*d2c2 - 2.*c2*c6*d2c2 + 3.*c1*c3to2*d2c3 
                - 2.*c1*c2*c4*d2c3 - 2.*c3*c4*d2c3 - c1to2*c5*d2c3 + c2*c5*d2c3 + c1*c6*d2c3 - 2.*c1*c2*c3*d2c4 
                - c3to2*d2c4 + 2.*c1to2*c4*d2c4 + 2.*c2*c4*d2c4 - c1*c5*d2c4 + c1*c2to2*d2c5 - c1to2*c3*d2c5 
                + c2*c3*d2c5 - c1*c4*d2c5 - c2to2*d2c6 + c1*c3*d2c6 + 2.*c1*c5*dc2to2 - 2.*c6*dc2to2 - 4.*c1*c4*dc2*dc3 
                + 2.*c5*dc2*dc3 + 6.*c1*c3*dc3to2 - 2.*c4*dc3to2 - 4.*c1*c3*dc2*dc4 + 4.*c4*dc2*dc4 - 4.*c1*c2*dc3*dc4 
                - 4.*c3*dc3*dc4 + 2.*c1to2*dc4to2 + 2.*c2*dc4to2 + 4.*c1*c2*dc2*dc5 + 2.*c3*dc2*dc5 - 2.*c1to2*dc3*dc5 
                + 2.*c2*dc3*dc5 - 2.*c1*dc4*dc5 - 4.*c2*dc2*dc6 + 2.*c1*dc3*dc6;
  double n2   = c2*c3to3 - 2.*c2to2*c3*c4 - c1*c3to2*c4 + 2.*c1*c2*c4to2 + c3*c4to2 + c2to3*c5 - c3to2*c5 - c1to2*c4*c5 
                - c2*c4*c5 + c1*c5to2 - c1*c2to2*c6 + c1to2*c3*c6 + c2*c3*c6 - c1*c4*c6;
  double Dn2  = c3to3*dc2 - 4.*c2*c3*c4*dc2 + 2.*c1*c4to2*dc2 + 3.*c2to2*c5*dc2 - c4*c5*dc2 - 2.*c1*c2*c6*dc2 + c3*c6*dc2 
                + 3.*c2*c3to2*dc3 - 2.*c2to2*c4*dc3 - 2.*c1*c3*c4*dc3 + c4to2*dc3 - 2.*c3*c5*dc3 + c1to2*c6*dc3 + c2*c6*dc3 
                - 2.*c2to2*c3*dc4 - c1*c3to2*dc4 + 4.*c1*c2*c4*dc4 + 2.*c3*c4*dc4 - c1to2*c5*dc4 - c2*c5*dc4 - c1*c6*dc4 
                + c2to3*dc5 - c3to2*dc5 - c1to2*c4*dc5 - c2*c4*dc5 + 2.*c1*c5*dc5 - c1*c2to2*dc6 + c1to2*c3*dc6 + c2*c3*dc6 
                - c1*c4*dc6;
  double D2n2 = c3to3*d2c2 - 4.*c2*c3*c4*d2c2 + 2.*c1*c4to2*d2c2 + 3.*c2to2*c5*d2c2 - c4*c5*d2c2 - 2.*c1*c2*c6*d2c2 + c3*c6*d2c2 
                + 3.*c2*c3to2*d2c3 - 2.*c2to2*c4*d2c3 - 2.*c1*c3*c4*d2c3 + c4to2*d2c3 - 2.*c3*c5*d2c3 + c1to2*c6*d2c3 + c2*c6*d2c3 
                - 2.*c2to2*c3*d2c4 - c1*c3to2*d2c4 + 4.*c1*c2*c4*d2c4 + 2.*c3*c4*d2c4 - c1to2*c5*d2c4 - c2*c5*d2c4 - c1*c6*d2c4 
                + c2to3*d2c5 - c3to2*d2c5 - c1to2*c4*d2c5 - c2*c4*d2c5 + 2.*c1*c5*d2c5 - c1*c2to2*d2c6 + c1to2*c3*d2c6 
                + c2*c3*d2c6 - c1*c4*d2c6 - 4.*c3*c4*dc2to2 + 6.*c2*c5*dc2to2 - 2.*c1*c6*dc2to2 + 6.*c3to2*dc2*dc3 - 8.*c2*c4*dc2*dc3 
                + 2.*c6*dc2*dc3 + 6.*c2*c3*dc3to2 - 2.*c1*c4*dc3to2 - 2.*c5*dc3to2 - 8.*c2*c3*dc2*dc4 + 8.*c1*c4*dc2*dc4 - 2.*c5*dc2*dc4 
                - 4.*c2to2*dc3*dc4 - 4.*c1*c3*dc3*dc4 + 4.*c4*dc3*dc4 + 4.*c1*c2*dc4to2 + 2.*c3*dc4to2 + 6.*c2to2*dc2*dc5 - 2.*c4*dc2*dc5 
                - 4.*c3*dc3*dc5 - 2.*c1to2*dc4*dc5 - 2.*c2*dc4*dc5 + 2.*c1*dc5to2 - 4.*c1*c2*dc2*dc6 + 2.*c3*dc2*dc6 + 2.*c1to2*dc3*dc6 
                + 2.*c2*dc3*dc6 - 2.*c1*dc4*dc6;
  double n3   = c3to4 - 3.*c2*c3to2*c4 + c2to2*c4to2 + 2.*c1*c3*c4to2 - c4to3 + 2.*c2to2*c3*c5 - 2.*c1*c3to2*c5 - 2.*c1*c2*c4*c5 
                + 2.*c3*c4*c5 + c1to2*c5to2 - c2*c5to2 - c2to3*c6 + 2.*c1*c2*c3*c6 - c3to2*c6 - c1to2*c4*c6 + c2*c4*c6;
  double Dn3  = - 3.*c3to2*c4*dc2 + 2.*c2*c4to2*dc2 + 4.*c2*c3*c5*dc2 - 2.*c1*c4*c5*dc2 - c5to2*dc2 - 3.*c2to2*c6*dc2 + 2.*c1*c3*c6*dc2 
                + c4*c6*dc2 + 4.*c3to3*dc3 - 6.*c2*c3*c4*dc3 + 2.*c1*c4to2*dc3 + 2.*c2to2*c5*dc3 - 4.*c1*c3*c5*dc3 + 2.*c4*c5*dc3 
                + 2.*c1*c2*c6*dc3 - 2.*c3*c6*dc3 - 3.*c2*c3to2*dc4 + 2.*c2to2*c4*dc4 + 4.*c1*c3*c4*dc4 - 3.*c4to2*dc4 - 2.*c1*c2*c5*dc4 
                + 2.*c3*c5*dc4 - c1to2*c6*dc4 + c2*c6*dc4 + 2.*c2to2*c3*dc5 - 2.*c1*c3to2*dc5 - 2.*c1*c2*c4*dc5 + 2.*c3*c4*dc5 
                + 2.*c1to2*c5*dc5 - 2.*c2*c5*dc5 - c2to3*dc6 + 2.*c1*c2*c3*dc6 - c3to2*dc6 - c1to2*c4*dc6 + c2*c4*dc6;
  double D2n3 = - 3.*c3to2*c4*d2c2 + 2.*c2*c4to2*d2c2 + 4.*c2*c3*c5*d2c2 - 2.*c1*c4*c5*d2c2 - c5to2*d2c2 - 3.*c2to2*c6*d2c2 
                + 2.*c1*c3*c6*d2c2 + c4*c6*d2c2 + 4.*c3to3*d2c3 - 6.*c2*c3*c4*d2c3 + 2.*c1*c4to2*d2c3 + 2.*c2to2*c5*d2c3 
                - 4.*c1*c3*c5*d2c3 + 2.*c4*c5*d2c3 + 2.*c1*c2*c6*d2c3 - 2.*c3*c6*d2c3 - 3.*c2*c3to2*d2c4 + 2.*c2to2*c4*d2c4 
                + 4.*c1*c3*c4*d2c4 - 3.*c4to2*d2c4 - 2.*c1*c2*c5*d2c4 + 2.*c3*c5*d2c4 - c1to2*c6*d2c4 + c2*c6*d2c4 + 2.*c2to2*c3*d2c5 
                - 2.*c1*c3to2*d2c5 - 2.*c1*c2*c4*d2c5 + 2.*c3*c4*d2c5 + 2.*c1to2*c5*d2c5 - 2.*c2*c5*d2c5 - c2to3*d2c6 + 2.*c1*c2*c3*d2c6 
                - c3to2*d2c6 - c1to2*c4*d2c6 + c2*c4*d2c6 + 2.*c4to2*dc2to2 + 4.*c3*c5*dc2to2 - 6.*c2*c6*dc2to2 - 12*c3*c4*dc2*dc3 
                + 8.*c2*c5*dc2*dc3 + 4.*c1*c6*dc2*dc3 + 12*c3to2*dc3to2 - 6.*c2*c4*dc3to2 - 4.*c1*c5*dc3to2 - 2.*c6*dc3to2 
                - 6.*c3to2*dc2*dc4 + 8.*c2*c4*dc2*dc4 - 4.*c1*c5*dc2*dc4 + 2.*c6*dc2*dc4 - 12*c2*c3*dc3*dc4 + 8.*c1*c4*dc3*dc4 
                + 4.*c5*dc3*dc4 + 2.*c2to2*dc4to2 + 4.*c1*c3*dc4to2 - 6.*c4*dc4to2 + 8.*c2*c3*dc2*dc5 - 4.*c1*c4*dc2*dc5 - 4.*c5*dc2*dc5 
                + 4.*c2to2*dc3*dc5 - 8.*c1*c3*dc3*dc5 + 4.*c4*dc3*dc5 - 4.*c1*c2*dc4*dc5 + 4.*c3*dc4*dc5 + 2.*c1to2*dc5to2 - 2.*c2*dc5to2 
                - 6.*c2to2*dc2*dc6 + 4.*c1*c3*dc2*dc6 + 2.*c4*dc2*dc6 + 4.*c1*c2*dc3*dc6 - 4.*c3*dc3*dc6 - 2.*c1to2*dc4*dc6 + 2.*c2*dc4*dc6;

  // Coefficients of the denominator + derivatives
  double d0   = c3to3 - 2.*c2*c3*c4 + c1*c4to2 + c2to2*c5 - c1*c3*c5;
  double Dd0  = - 2.*c3*c4*dc2 + 2.*c2*c5*dc2 + 3.*c3to2*dc3 - 2.*c2*c4*dc3 - c1*c5*dc3 - 2.*c2*c3*dc4 + 2.*c1*c4*dc4 + c2to2*dc5 
                - c1*c3*dc5;
  double D2d0 = - 2.*c3*c4*d2c2 + 2.*c2*c5*d2c2 + 3.*c3to2*d2c3 - 2.*c2*c4*d2c3 - c1*c5*d2c3 - 2.*c2*c3*d2c4 + 2.*c1*c4*d2c4 
                + c2to2*d2c5 - c1*c3*d2c5 + 2.*c5*dc2to2 - 4.*c4*dc2*dc3 + 6.*c3*dc3to2 - 4.*c3*dc2*dc4 - 4.*c2*dc3*dc4 
                + 2.*c1*dc4to2 + 4.*c2*dc2*dc5 - 2.*c1*dc3*dc5;
  double d1   = - c3to2*c4 + c2*c4to2 + c2*c3*c5 - c1*c4*c5 - c2to2*c6 + c1*c3*c6;
  double Dd1  = c4to2*dc2 + c3*c5*dc2 - 2.*c2*c6*dc2 - 2.*c3*c4*dc3 + c2*c5*dc3 + c1*c6*dc3 - c3to2*dc4 + 2.*c2*c4*dc4 
                - c1*c5*dc4 + c2*c3*dc5 - c1*c4*dc5 - c2to2*dc6 + c1*c3*dc6;
  double D2d1 = c4to2*d2c2 + c3*c5*d2c2 - 2.*c2*c6*d2c2 - 2.*c3*c4*d2c3 + c2*c5*d2c3 + c1*c6*d2c3 - c3to2*d2c4 + 2.*c2*c4*d2c4 
                - c1*c5*d2c4 + c2*c3*d2c5 - c1*c4*d2c5 - c2to2*d2c6 + c1*c3*d2c6 - 2.*c6*dc2to2 + 2.*c5*dc2*dc3 - 2.*c4*dc3to2 
                + 4.*c4*dc2*dc4 - 4.*c3*dc3*dc4 + 2.*c2*dc4to2 + 2.*c3*dc2*dc5 + 2.*c2*dc3*dc5 - 2.*c1*dc4*dc5 - 4.*c2*dc2*dc6 + 2.*c1*dc3*dc6;
  double d2   = c3*c4to2 - c3to2*c5 - c2*c4*c5 + c1*c5to2 + c2*c3*c6 - c1*c4*c6;
  double Dd2  = - c4*c5*dc2 + c3*c6*dc2 + c4to2*dc3 - 2.*c3*c5*dc3 + c2*c6*dc3 + 2.*c3*c4*dc4 - c2*c5*dc4 - c1*c6*dc4 - c3to2*dc5 
                - c2*c4*dc5 + 2.*c1*c5*dc5 + c2*c3*dc6 - c1*c4*dc6;
  double D2d2 = - c4*c5*d2c2 + c3*c6*d2c2 + c4to2*d2c3 - 2.*c3*c5*d2c3 + c2*c6*d2c3 + 2.*c3*c4*d2c4 - c2*c5*d2c4 - c1*c6*d2c4 
                - c3to2*d2c5 - c2*c4*d2c5 + 2.*c1*c5*d2c5 + c2*c3*d2c6 - c1*c4*d2c6 + 2.*c6*dc2*dc3 - 2.*c5*dc3to2 - 2.*c5*dc2*dc4 
                + 4.*c4*dc3*dc4 + 2.*c3*dc4to2 - 2.*c4*dc2*dc5 - 4.*c3*dc3*dc5 - 2.*c2*dc4*dc5 + 2.*c1*dc5to2 + 2.*c3*dc2*dc6 
                + 2.*c2*dc3*dc6 - 2.*c1*dc4*dc6;
  double d3   = - c4to3 + 2.*c3*c4*c5 - c2*c5to2 - c3to2*c6 + c2*c4*c6;
  double Dd3  = - c5to2*dc2 + c4*c6*dc2 + 2.*c4*c5*dc3 - 2.*c3*c6*dc3 - 3.*c4to2*dc4 + 2.*c3*c5*dc4 + c2*c6*dc4 + 2.*c3*c4*dc5 
                - 2.*c2*c5*dc5 - c3to2*dc6 + c2*c4*dc6;
  double D2d3 = - c5to2*d2c2 + c4*c6*d2c2 + 2.*c4*c5*d2c3 - 2.*c3*c6*d2c3 - 3.*c4to2*d2c4 + 2.*c3*c5*d2c4 + c2*c6*d2c4 + 2.*c3*c4*d2c5 
                - 2.*c2*c5*d2c5 - c3to2*d2c6 + c2*c4*d2c6 - 2.*c6*dc3to2 + 2.*c6*dc2*dc4 + 4.*c5*dc3*dc4 - 6.*c4*dc4to2 - 4.*c5*dc2*dc5 
                + 4.*c4*dc3*dc5 + 4.*c3*dc4*dc5 - 2.*c2*dc5to2 + 2.*c4*dc2*dc6 - 4.*c3*dc3*dc6 + 2.*c2*dc4*dc6;

  double Num   = n0 + n1*u + n2*u2 + n3*u3;
  double dNum  = Dn0 + n1 + Dn1*u + 2.*n2*u + Dn2*u2 + 3.*n3*u2 + Dn3*u3;
  double d2Num = D2n0 + 2.*Dn1 + 2.*n2 + D2n1*u + 4.*Dn2*u + 6.*n3*u + D2n2*u2 + 6.*Dn3*u2 + D2n3*u3;
  double Den   = d0 + d1*u + d2*u2 + d3*u3;
  double dDen  = d1 + Dd0 + 2.*d2*u + Dd1*u + 3.*d3*u2 + Dd2*u2 + Dd3*u3;
  double d2Den = 2.*d2 + D2d0 + 2.*Dd1 + D2d1*u + 6.*d3*u + 4.*Dd2*u + D2d2*u2 + 6.*Dd3*u2 + D2d3*u3;

  double dNumDen = dNum*Den - Num*dDen; 
  double Den2    = Den*Den;
  double Den3    = Den2*Den;

  *P   = Num/Den;
  *DP  = dNumDen/Den2;
  *D2P = (Den*(d2Num*Den - Num*d2Den) - 2.*dDen*dNumDen)/Den3;

}

/**
  * Function: Pade33_forGSF function
  * --------------------------------
  *   Padé 7,6 of f[v] = 1 + c1 v^7 + c2 v^9 + c3[v] v^11 + c4[v] v^13  (where v = sqrt(u))
  *   (functional form of the tail part of a1SF & d1SF)
  *   + notice that c1 and c2 do *not* depend on u so they have no derivatives.
  * 
  *    @param[in] coeffs: coefficients of the Taylor-expanded function + derivatives (except from c1, they are functions of log(u))
  *    @param[in] Dcoeffs: derivatives of the Taylor-expanded function + derivatives (except from c1, they are functions of log(u))
  *    @param[in] D2coeffs: second derivatives of the Taylor-expanded function + derivatives (except from c1, they are functions of log(u))
  *    @param[in] u: argument
  *    @param[out] P: Pade33_forGSF function
  *    @param[out] DP: derivative of Pade33_forGSF function
  *    @param[out] D2P: second derivative of Pade33_forGSF function
  */
void Pade76v1_forGSF(double coeffs[4], double Dcoeffs[2], double D2coeffs[2], double u, double *P, double *DP, double *D2P)
{ 
  double u2    = u*u;
  double u3    = u2*u;
  double u4    = u3*u;
  double u5    = u4*u;
  double u6    = u5*u;
  double u3by2 = sqrt(u3);
  double u5by2 = sqrt(u5);
  double u7by2 = sqrt(u6*u);

  double c1 = coeffs[0];
  double c2 = coeffs[1];
  double c3 = coeffs[2];
  double c4 = coeffs[3];

  double dc3 = Dcoeffs[0];
  double dc4 = Dcoeffs[1];

  double d2c3 = D2coeffs[0];
  double d2c4 = D2coeffs[1];

  // Powers
  double c1to2 = c1*c1;
  double c1to3 = c1to2*c1;
  double c1to4 = c1to3*c1;
  double c2to2 = c2*c2;
  double c2to3 = c2to2*c2;

  // Coefficients of the numerator + derivatives
  double n0    =  c1to3;
  double n1    = -c1to2*c2;
  double n2    =  c1*c2to2 - c1to2*c3;
  double Dn2   = -c1to2*dc3;
  double D2n2  = -c1to2*d2c3;
  double n3    = -c2to3 + 2*c1*c2*c3 - c1to2*c4;
  double Dn3   =  2*c1*c2*dc3 - c1to2*dc4;
  double D2n3  =  2*c1*c2*d2c3 - c1to2*d2c4;
  double n7by2 =  c1to4;

  // Coefficients of the denominator + derivatives
  double d0   = n0;
  double d1   = n1;
  double d2   = n2;
  double Dd2  = Dn2;
  double D2d2 = D2n2;
  double d3   = n3;
  double Dd3  = Dn3;
  double D2d3 = D2n3;

  double Num   = n0 + n1*u + n2*u2 + n3*u3 + n7by2*u7by2;
  double dNum  = n1 + 2.*n2*u + Dn2*u2 + 3.*n3*u2 + Dn3*u3 + 0.5*(7.*n7by2*u5by2);
  double d2Num = 2.*n2 + 4.*Dn2*u + 6.*n3*u + D2n2*u2 + 6.*Dn3*u2 + D2n3*u3 + 0.25*(35.*n7by2*u3by2);
  double Den   = d0 + d1*u + d2*u2 + d3*u3;
  double dDen  = d1 + 2.*d2*u + 3.*d3*u2 + Dd2*u2 + Dd3*u3;
  double d2Den = 2.*d2 + 6.*d3*u + 4*Dd2*u + D2d2*u2 + 6.*Dd3*u2 + D2d3*u3;

  double dNumDen = dNum*Den - Num*dDen; // (Num*Den)'
  double Den2    = Den*Den;
  double Den3    = Den2*Den;

  *P   = Num/Den;
  *DP  = dNumDen/Den2;
  *D2P = (Den*(d2Num*Den - Num*d2Den) - 2.*dDen*dNumDen)/Den3;

}

/**
 * Function: Pade76v2_forGSF function
 * ----------------------------------
 *   Padé 7,6 of f[v] = 1 + c1 v^5 + c2 v^7 + c3[v] v^9 + c4[v] v^11 + c5[v] v^13  (where v = sqrt(u))
 *   (functional form of the tail part of q1SF)
 *   + notice that c1 and c2 do *not* depend on u so they have no derivatives.
 *    
 *   @param[in] coeffs: coefficients
 *   @param[in] Dcoeffs: derivatives of the coefficients
 *   @param[in] D2coeffs: second derivatives of the coefficients
 *   @param[in] u: argument
 *   @param[out] P: Pade76v2_forGSF function
 *   @param[out] DP: derivative of Pade76v2_forGSF function
 *   @param[out] D2P: second derivative of Pade76v2_forGSF function
 *   
*/
void Pade76v2_forGSF(double coeffs[5], double Dcoeffs[3], double D2coeffs[3], double u, double *P, double *DP, double *D2P)
{
  double sqrtu = sqrt(u);
  double u2    = u*u;
  double u3    = u2*u;
  double u4    = u3*u;
  double u5    = u4*u;
  double u6    = u5*u;
  double u3by2 = sqrt(u3);
  double u5by2 = sqrt(u5);
  double u7by2 = sqrt(u6*u);

  // Coefficients of the Taylor-expanded function + derivatives of c3, c4, c5
  double c1 = coeffs[0];
  double c2 = coeffs[1];
  double c3 = coeffs[2];
  double c4 = coeffs[3];
  double c5 = coeffs[4];

  double dc3 = Dcoeffs[0];
  double dc4 = Dcoeffs[1];
  double dc5 = Dcoeffs[2];

  double d2c3 = D2coeffs[0];
  double d2c4 = D2coeffs[1];
  double d2c5 = D2coeffs[2];

  // Powers of the coefficients
  double c1to2 = c1*c1;
  double c2to2 = c2*c2;
  double c2to3 = c2to2*c2;
  double c2to4 = c2to3*c2;
  double c3to2 = c3*c3;
  double c3to3 = c3to2*c3;
  double c4to2 = c4*c4;

  // Powers of the derivatives
  double dc3to2 = dc3*dc3;
  double dc4to2 = dc4*dc4;

  // Coefficients of the numerator + derivatives
  double n0      =  c2to3 - 2.*c1*c2*c3 + c1to2*c4;
  double Dn0     = -2.*c1*c2*dc3 + c1to2*dc4;
  double D2n0    = -2.*c1*c2*d2c3 + c1to2*d2c4;
  double n1      = -c2to2*c3 + c1*c2*c4 + c1*(c3to2 - c1*c5);
  double Dn1     = -c2to2*dc3 + c1*c2*dc4 + c1*(2.*c3*dc3 - c1*dc5);
  double D2n1    = -c2to2*d2c3 + c1*c2*d2c4 + c1*(2.*c3*d2c3 - c1*d2c5 + 2.*dc3to2);
  double n2      = -c2to2*c4 - c1*c3*c4 + c2*(c3to2 + c1*c5);
  double Dn2     = -c1*c4*dc3 - c2to2*dc4 - c1*c3*dc4 + c2*(2.*c3*dc3 + c1*dc5);
  double D2n2    = -c1*c4*d2c3 - c2to2*d2c4 - c1*c3*d2c4 + c2*(2.*c3*d2c3 + c1*d2c5 + 2.*dc3to2) - 2.*c1*dc3*dc4;
  double n5by2   =  c1*(c2to3 - 2.*c1*c2*c3 + c1to2*c4);
  double Dn5by2  =  c1*(-2.*c1*c2*dc3 + c1to2*dc4);
  double D2n5by2 =  c1*(-2.*c1*c2*d2c3 + c1to2*d2c4);
  double n3      = -c3to3 - c1*c4to2 - c2to2*c5 + c3*(2*c2*c4 + c1*c5);
  double Dn3     = -3.*c3to2*dc3 + (2.*c2*c4 + c1*c5)*dc3 - 2.*c1*c4*dc4 - c2to2*dc5 + c3*(2.*c2*dc4 + c1*dc5);
  double D2n3    = -3.*c3to2*d2c3 + (2.*c2*c4 + c1*c5)*d2c3 - 2.*c1*c4*d2c4 - c2to2*d2c5 + c3*(2.*c2*d2c4 + c1*d2c5) 
                   -6.*c3*dc3to2 - 2.*c1*dc4to2 + 2.*dc3*(2.*c2*dc4 + c1*dc5);
  double n7by2   = c2to4 - 3.*c1*c2to2*c3 + 2.*c1to2*c2*c4 + c1to2*(c3to2 - c1*c5);
  double Dn7by2  = -3.*c1*c2to2*dc3 + 2.*c1to2*c2*dc4 + c1to2*(2.*c3*dc3 - c1*dc5);
  double D2n7by2 = -3.*c1*c2to2*d2c3 + 2.*c1to2*c2*d2c4 + c1to2*(2.*c3*d2c3 - c1*d2c5 + 2.*dc3to2);

  // Coefficients of the denominator + derivatives
  double d0   = n0;
  double Dd0  = Dn0;
  double D2d0 = D2n0;
  double d1   = n1;
  double Dd1  = Dn1;
  double D2d1 = D2n1;
  double d2   = n2;
  double Dd2  = Dn2;
  double D2d2 = D2n2;
  double d3   = n3;
  double Dd3  = Dn3;
  double D2d3 = D2n3;

  double Num   = n0 + n1*u + n2*u2 + n3*u3 + n5by2*u5by2 + n7by2*u7by2;
  double dNum  = Dn0 + n1 + Dn1*u + 2.*n2*u + Dn2*u2 + 3.*n3*u2 + Dn3*u3 + 0.5*(5.*n5by2*u3by2) 
               + Dn5by2*u5by2 + 0.5*(7.*n7by2*u5by2) + Dn7by2*u7by2;
  double d2Num = D2n0 + 2.*Dn1 + 2.*n2 + 0.25*(15.*n5by2*sqrtu) + D2n1*u + 4.*Dn2*u + 6.*n3*u + D2n2*u2 
               + 6.*Dn3*u2 + D2n3*u3 + 5.*Dn5by2*u3by2 + 0.25*(35.*n7by2*u3by2) + D2n5by2*u5by2 + 7.*Dn7by2*u5by2 + D2n7by2*u7by2;
  double Den   = d0 + d1*u + d2*u2 + d3*u3;
  double dDen  = d1 + Dd0 + 2.*d2*u + Dd1*u + 3.*d3*u2 + Dd2*u2 + Dd3*u3;
  double d2Den = 2.*d2 + D2d0 + 2.*Dd1 + D2d1*u + 6.*d3*u + 4.*Dd2*u + D2d2*u2 + 6.*Dd3*u2 + D2d3*u3;

  double dNumDen = dNum*Den - Num*dDen; // (Num*Den)'
  double Den2    = Den*Den;
  double Den3    = Den2*Den;

  *P   = Num/Den;
  *DP  = dNumDen/Den2;
  *D2P = (Den*(d2Num*Den - Num*d2Den) - 2.*dDen*dNumDen)/Den3;

}

/**
 * Function: Taylorseries
 * ----------------------
 *   Return Taylor series
 * 
 *   @param[in] x: argument
 *   @param[in] a: array of coefficients
 *   @param[in] N: order of the series
 * 
 *   @return Taylorseries(x,a,N): Taylor series
*/
double Taylorseries(double x, double *a, int N){
  double xn[N+1];
  xn[0] = 1.;
  for (int n=1; n<N+1; n++) {	
    xn[n] = xn[n-1]*x;
  }
  
  double sum = 0.;
#if (1)
  sum = a[0];
  for (int n=1; n<N+1; n++) {	
    sum += a[n]*xn[n];
  }
#else  
  sum = a[N]*xn[N];
  for (int n=N-1; n-- > 1;) {
    sum += a[n]*xn[n];
  } 
  sum rholm[k] += a[0];
#endif

  return sum;
}

/**
 * Function: vect_dot
 * ------------------
 * Return vector scalar product
 * 
 * @param[in] ax: x component of vector a
 * @param[in] ay: y component of vector a
 * @param[in] az: z component of vector a
 * @param[in] bx: x component of vector b
 * @param[in] by: y component of vector b
 * @param[in] bz: z component of vector b
 * 
 * @return s: vector scalar product
 */
void vect_dot(double ax, double ay, double az,
	      double bx, double by, double bz,
	      double *s)
{
  *s = ax*bx + ay*by + az*bz;
}

/**
 * Function: vect_dot3
 * -------------------
 * Return vector scalar product for 3d vectors
 * 
 * @param[in] a: vector a
 * @param[in] b: vector b
 * 
 * @return s: vector scalar product
*/
void vect_dot3(double *a, double *b, 
	       double *s)
{
  *s = a[Ix]*b[Ix] + a[Iy]*b[Iy] + a[Iz]*b[Iz];
}

/**
 * Function: vect_cross
 * --------------------
 * Return vector cross product
 * 
 * @param[in] ax: x component of vector a
 * @param[in] ay: y component of vector a
 * @param[in] az: z component of vector a
 * @param[in] bx: x component of vector b
 * @param[in] by: y component of vector b
 * @param[in] bz: z component of vector b
 * @param[out] cx: x component of vector c
 * @param[out] cy: y component of vector c
 * @param[out] cz: z component of vector c
*/
void vect_cross(double ax, double ay, double az,
		double bx, double by, double bz,
		double *cx, double *cy, double *cz)
{
  *cx = ay*bz - az*by;
  *cy = az*bx - ax*bz;
  *cz = ax*by - ay*bx;
}

/**
 * Function: vect_cross3
 * ---------------------
 * Return vector cross product for 3d vectors
 * 
 * @param[in] a: vector a
 * @param[in] b: vector b
 * @param[out] c: vector c
*/
void vect_cross3(double *a, double *b, 
		 double *c)
{
  c[Ix] = a[Iy]*b[Iz] - a[Iz]*b[Iy];
  c[Iy] = a[Iz]*b[Ix] - a[Ix]*b[Iz];
  c[Iz] = a[Ix]*b[Iy] - a[Iy]*b[Ix];
}

/**
 * Function: vect_rotate
 * ---------------------
 * Rotate vector about an axis
 * 
 * @param[in] axis: axis of rotation
 * @param[in] angle: angle of rotation
 * @param[in,out] vx_p: x component of vector
 * @param[in,out] vy_p: y component of vector
 * @param[in,out] vz_p: z component of vector
*/
void vect_rotate(int axis, double angle, double *vx_p, double *vy_p, double *vz_p)
{
  double tmp1, tmp2;
  const double ca = cos(angle);
  const double sa = sin(angle);

  double vx = *vx_p;
  double vy = *vy_p;
  double vz = *vz_p;  

  if (axis==3) {
    // z-axis
    tmp1 = vx * ca - vy * sa;	    
    tmp2 = vx * sa + vy * ca;
    vx = tmp1;
    vy = tmp2;
  } else if (axis==2) {
    // y-axis
    tmp1 = vx * ca + vz * sa;
    tmp2 = - vx * sa + vz * ca;
    vx = tmp1;
    vz = tmp2;
  } else if (axis==1) {
    // x-axis
    tmp1 = vy * ca - vz * sa;
    tmp2 = vy * sa + vz * ca;
    vy = tmp1;
    vz = tmp2;
  } else
    errorexit("unkown rotation axis");
  
  *vx_p = vx;
  *vy_p = vy;
  *vz_p = vz;  
}

/**
 * Function: vect_rotate3
 * ----------------------
 * Rotate 3d vector about an axis
 * 
 * @param[in] axis: axis of rotation
 * @param[in] angle: angle of rotation
 * @param[in,out] v: vector
*/
void vect_rotate3(int axis, double angle, double *v)
{
  double tmp1, tmp2;
  const double ca = cos(angle);
  const double sa = sin(angle);

  double vx = v[Ix];
  double vy = v[Iy];
  double vz = v[Iz];  

  if (axis==3) {
    // z-axis
    tmp1 = vx * ca - vy * sa;	    
    tmp2 = vx * sa + vy * ca;
    vx = tmp1;
    vy = tmp2;
  } else if (axis==2) {
    // y-axis
    tmp1 = vx * ca + vz * sa;
    tmp2 = - vx * sa + vz * ca;
    vx = tmp1;
    vz = tmp2;
  } else if (axis==1) {
    // x-axis
    tmp1 = vy * ca - vz * sa;
    tmp2 = vy * sa + vz * ca;
    vy = tmp1;
    vz = tmp2;
  } else
    errorexit("unkown rotation axis");
  
  v[Ix] = vx;
  v[Iy] = vy;
  v[Iz] = vz;  
}

/**
 * Function: interp_spline_pt
 * --------------------------
 *   Spline interpolation with GSL routines at single point
 *   
 *   @param[in] t: array of x values
 *   @param[in] y: array of y values
 *   @param[in] n: number of points
 *   @param[in] ti: x value to interpolate
 * 
 *   @return yi: interpolated value
*/
double interp_spline_pt(double *t, double *y, int n, double ti)
{
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n);
  gsl_spline_init (spline, t, y, n);    
  double yi = gsl_spline_eval (spline, ti, acc);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return yi;
}

/**
  * Function: interp_spline
  * -----------------------
  *   Spline interpolation with GSL routines
  *   
  *   @param[in] t: array of x values
  *   @param[in] y: array of y values
  *   @param[in] n: number of points
  *   @param[in] ti: array of x values to interpolate
  *   @param[in] ni: number of points to interpolate
  *   @param[out] yi: array of interpolated values
*/
void interp_spline(double *t, double *y, int n, double *ti, int ni, double *yi)
{
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n);
  gsl_spline_init (spline, t, y, n);    
  for (int k = 0; k < ni; k++) {
    yi[k] = gsl_spline_eval (spline, *(ti + k), acc);
    /* yi[k] = gsl_spline_eval (spline, ti[k], acc); */
  }
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
}

/**
  * Function: interp_spline_checklim
  * --------------------------------
  *   Spline interpolation with GSL routines
  *   Check limits and set to zero points outside interp range
  *   
  *   @param[in] t: array of x values
  *   @param[in] y: array of y values
  *   @param[in] n: number of points
  *   @param[in] ti: array of x values to interpolate
  *   @param[in] ni: number of points to interpolate
  *   @param[out] yi: array of interpolated values
*/
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

/* An OpenMP SIMD version - experimental
 * SIMD Functions need to be enabled, see e.g.
 * P.25 https://info.ornl.gov/sites/publications/files/Pub69214.pdf
 * GSL pointers:
 * https://github.com/ampl/gsl/blob/master/interpolation/spline.c#L118
 * https://github.com/ampl/gsl/blob/master/interpolation/interp.c#L140
*/
#pragma omp declare simd uniform(_spline) linear(_x) // linear(_acc) ?
double gsl_spline_eval_simd_enabled(const gsl_spline * _spline, double * _x, gsl_interp_accel * _acc)
{
  return gsl_spline_eval(_spline, *_x, _acc);
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
  for (int k = 0; k < ni; k++) {
    /* yi[k] = gsl_spline_eval_simd_enabled(spline, ti + k, NULL); */ // bsearch, 2x slower
    yi[k] = gsl_spline_eval_simd_enabled(spline, ti + k, acc);
  }
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
#ifdef _OPENMP  
  if (USETIMERS) openmp_timer_stop("interp_spline");
#endif
}
 
/**
  * Function: find_point_bisection
  * ------------------------------
  *   Find nearest point index in 1d array
  *   
  *   @param[in] x: value to find
  *   @param[in] n: number of points
  *   @param[in] xp: array of x values
  *   @param[in] o: offset
  * 
  *   @return i: index of nearest point
*/
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

#define tiny 1e-12
/**
  * Function: baryc_f
  * -----------------
  *   Barycentric Lagrange interpolation at xx with n points of f(x), 
  *   equivalent to standard Lagrangian interpolation
  *   
  *   @param[in] xx: x value to interpolate
  *   @param[in] n: number of points
  *   @param[in] f: array of y values
  *   @param[in] x: array of x values
  * 
  *   @return f(xx): interpolated value
*/
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

/**
 * Function: baryc_weights
 * -----------------------
 *   Barycentric Lagrange interpolation at xx with n points of f(x), 
 *   compute weights
 *   
 *   @param[in] n: number of points
 *   @param[in] x: array of x values
 *   @param[out] omega: array of weights
 */
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

/**
 * Function: baryc_f_weights
 * -------------------------
 *   Barycentric Lagrange interpolation at xx with n points of f(x), 
 *   use precomputed weights
 *   
 *   @param[in] xx: x value to interpolate
 *   @param[in] n: number of points
 *   @param[in] f: array of y values
 *   @param[in] x: array of x values
 *   @param[in] omega: array of weights
 * 
 *   @return f(xx): interpolated value
 */
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

/**
 * Function: interp1d
 * ------------------
 *   1d Lagrangian barycentric interpolation
 *   
 *   @param[in] order: order of interpolation
 *   @param[in] xx: x value to interpolate
 *   @param[in] nx: number of points
 *   @param[in] f: array of y values
 *   @param[in] x: array of x values
 * 
 *   @return f(xx): interpolated value
 */
double interp1d (const int order, double xx, int nx, double *f, double *x)
{
  double ff;
  int ix;
  int ox = order > nx ? nx : order;
  ix = find_point_bisection(xx, nx, x, ox/2);
  ff = baryc_f(xx, ox, &f[ix], &x[ix]);  
  return( ff );
}

/**
 * Function: find_max
 * ------------------
 *   Find max location by poynomial interpolation around x0 (uniform grid)
 *   
 *   @param[in] n: number of points
 *   @param[in] dx: grid spacing
 *   @param[in] x0: centre of grid
 *   @param[in] f: array of y values
 *   @param[out] fmax: interpolated value at max
 * 
 *   @return xmax: max location
*/
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

/**
 * Function: find_max_grid
 * -------------------
 *   Find max location around x0 using 5 points (non-uniform grid)
 *   
 *   @param[in] x: array of x values
 *   @param[in] f: array of y values
 *   @param[in] x0: centre of grid
 *   @param[out] fmax: interpolated value at max
 * 
 *   @return xmax: max location
 */
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

/**
 * Function: D0
 * ------------
 *   4th order centered stencil first derivative, uniform grids
 *   
 *   @param[in]  f: array of y values
 *   @param[in]  dx: grid spacing
 *   @param[in]  n: number of points
 *   @param[out] df: array of derivatives
 * 
 *   @return OK: exit code
 */
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

/**
 * Function: D2
 * ------------
 *   4th order centered stencil second derivative, uniform grids
 *   
 *   @param[in]  f: array of y values
 *   @param[in]  dx: grid spacing
 *   @param[in]  n: number of points
 *   @param[out] d2f: array of derivatives
 * 
 *   @return OK: exit code
 */
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

/**
 * Function: D0_x_2
 * ----------------
 *   2nd order centered stencil first derivative, nonuniform grids
 *   
 *   @param[in]  f: array of y values
 *   @param[in]  x: array of x values
 *   @param[in]  n: number of points
 *   @param[out] df: array of derivatives
 * 
 *   @return OK: exit code
 */
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

/**
 * Function: D0_x_4
 * ----------------
 *   4th order centered stencil first derivative, nonuniform grids
 *   Based on Lagrange polynomial
 *   @note: can be further optimized
 * 
 *   @param[in]  f: array of y values
 *   @param[in]  x: array of x values
 *   @param[in]  n: number of points
 *   @param[out] df: array of derivatives
 * 
 *   @return OK: exit code
 */
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
  ix--;iy--;
  /*right boundary*/
  df[n-2] = d4(ix, iy, 3);
  df[n-1] = d4(ix, iy, 4);
  
  return OK;
}

/**
 * Function: d4
 * ------------
 *   Helper function for D0_x_4 
 *   By definition, we consider the polynomials l_j(x) such that
 *   p(x)    = \sum_j y_j \, l_j(x), 
 *   hence the derivative at x_i will be
 *   p'(x_i) = \sum_j y_j \, l'_j(x_i) 
 *   
 *   @param[in] x: array of x values
 *   @param[in] y: array of y values
 *   @param[in] i: index in [0,4]
 *   
 *   @return dy_i: derivative of y at x[i]
*/
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

/**
 * Function: l_deriv
 * -----------------
 *   Compute the derivative of the interpolating
 *   polynomial j at x[i]
 * 
 *  @param[in] x: array of x values
 *  @param[in] i: position of x at which l_j is evaluated
 *  @param[in] j: index of l_j
 * 
 *  @return ld_ji: derivative of l_j at x[i]
*/
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

/**
  * Function: cumtrapz
  * ------------------
  *   Trapezoidal rule
  *   
  *   @param[in] f: array of y values
  *   @param[in] x: array of x values
  *   @param[in] n: number of points
  *   @param[out] sum: array of cumulative sums
  * 
  *   @return sum[n-1]: total integral
*/
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

/**
  * Function: cumint3
  * -----------------
  *   Third-order polynomial integration
  *   
  *   @param[in] f: array of y values
  *   @param[in] x: array of x values
  *   @param[in] n: number of points
  *   @param[out] sum: array of cumulative sums
  * 
  *   @return sum[n-1]: total integral
*/
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

/**
 * Function: unwrap
 * ----------------
 *   Unwrap phase angles
 *   
 *   @param[in,out] p: array of phase angles
 *   @param[in] size: number of points
*/
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

/**
 * Function: unwrap_HM
 * -------------------
 *    Modified unwrap function for HM phases
 *   
 *   @param[in,out] p: array of phase angles
 *   @param[in] size: number of points
*/
void unwrap_HM(double *p, const int size)
{
  if (size < 1) return;
  
  int j;
  double delta, dphi, corr, curr, prev;

  dphi = 0.;
  corr = 0.;
  
  prev  = p[0];
  delta = p[1]-p[0];
  
  for (j = 1; j < size; j++){
    // Setting current data point
    p[j] += corr;
    curr = p[j];

    // Check if decreasing too much - adding 2Pi
    if((curr < prev - Pi) && (curr - prev < delta - Pi))
      dphi = TwoPi;

    // Check if increasing too much - removing 2Pi
    if((curr > prev + TwoPi) && (curr - prev > delta + Pi))
      dphi = -TwoPi;

    // Adding corrections
    corr += dphi;
    p[j] += dphi;

    // Resetting for next iteration
    prev = p[j];
    delta = p[j]-p[j-1];
    dphi = 0.0;
  }
  
}

/**
 * Function: unwrap_euler
 * ----------------------
 *   Unwrap function for euler angles
 *   
 *   @param[in,out] p: array of euler angles
 *   @param[in] size: number of points
*/
void unwrap_euler(double *p, const int size)
{
  if (size < 1) return;
  
  int j;
  double delta, dphi, corr, curr, prev;

  dphi = 0.;
  corr = 0.;
  
  prev  = p[0];
  delta = p[1]-p[0];
  
  for (j = 1; j < size; j++){
    // Setting current data point
    p[j] += corr;
    curr = p[j];

    // Check if decreasing too much - adding 2Pi
    if((curr < prev - Pi) && (curr - prev < delta - Pi))
      dphi = TwoPi;

    // Check if increasing too much - removing 2Pi
    if((curr > prev + Pi) && (curr - prev > delta + Pi))
      dphi = -TwoPi;

    // Adding corrections
    corr += dphi;
    p[j] += dphi;

    // Resetting for next iteration
    prev = p[j];
    delta = p[j]-p[j-1];
    dphi = 0.0;
  }
}

#define dbg_unwrap_proxy (0)  /* stops after routine, use: ./TEOBResumS.x test.par > out */ 
/**
 * Function: unwrap_proxy
 * ----------------------
 *   Unwrap using number of cycles from reference phase as proxy
 *   
 *   @param[in,out] p: array of phase angles
 *   @param[in] r: array of reference phase angles
 *   @param[in] size: number of points
 *   @param[in] shift0: shift phase to r[0]
*/
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

/**
  * Function: rmap
  * --------------
  *   Compute real/imag <-> amplitude/phase
  *   
  *   @param[in,out] re: array of real values
  *   @param[in,out] im: array of imag values
  *   @param[in,out] p: array of phase angles
  *   @param[in,out] a: array of amplitudes
  *   @param[in] mode: 0: (Re, Im) -> (Amplitude, phase)
  *                    1: (Amplitude, phase) -> (Re, Im)
*/
void rmap (double *re, double *im, double *p, double *a, const int mode)
{
  /* h =  A exp( -i phi) */
  if (mode) {
    /** (Re, Im) -> (Amplitude, phase) */
    *a = sqrt( SQ((*re)) + SQ((*im)) );
    *p = Pi - atan2((*im), (*re)); /* exp(- i phi) => Pi  */
    //TODO: unwrap?
  } else {
    /** (Amplitude, phase) -> (Re, Im) */
    *re = + (*a) * cos((*p));  
    *im = - (*a) * sin((*p));  
  }
}

/**
  * Function: rmap_twist
  * --------------------
  *   Compute real/imag <-> amplitude/phase
  *   for twisted waveforms. 
  *   @note: `rmap` added a "Pi" shift, hence we add a new function
  *   
  *   @param[in,out] re: array of real values
  *   @param[in,out] im: array of imag values
  *   @param[in,out] p: array of phase angles
  *   @param[in,out] a: array of amplitudes
  *   @param[in] mode: 0: (Re, Im) -> (Amplitude, phase)
  *                    1: (Amplitude, phase) -> (Re, Im)
*/
void rmap_twist (double *re, double *im, double *p, double *a, const int mode)
{
  /* h =  A exp( -i phi) */
  if (mode) {
    /** (Re, Im) -> (Amplitude, phase) */
    *a = sqrt( SQ((*re)) + SQ((*im)) );
    *p = - atan2((*im), (*re)); /* exp(- i phi) => Pi  */
  } else {
    /** (Amplitude, phase) -> (Re, Im) */
    *re = + (*a) * cos((*p));  
    *im = - (*a) * sin((*p));  
  }
}

/**
  * Function: set_multipolar_idx_mask
  * ---------------------------------
  *   This routine sets a 0/1 mask for the multipolar linear index 
  *   work for any parameter and can specify default all on/off
  *   
  *   @param[in,out] kmask: array of 0/1 mask
  *   @param[in] n: number of points
  *   @param[in] idx: array of indices
  *   @param[in] m: number of indices
  *   @param[in] on: 0: all off
  *                  1: all on
*/
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

/**
 * Function: get_uniform_size
 * --------------------------
 *   Compute size of a uniform grid t0:dt:tf
 *   
 *   @param[in] tN: final time
 *   @param[in] t0: initial time
 *   @param[in] dt: time step
 * 
 *   @return size: number of points
*/
long int get_uniform_size(const double tN, const double t0, const double dt)
{
  return ((long int)((tN - t0)/dt + 1)); 
}

/**
 * Function: intersect_int
 * -----------------------
 *   Intersection of two int arrays
 *   
 *   @param[in] a: array of int
 *   @param[in] size_a: size of a
 *   @param[in] b: array of int
 *   @param[in] size_b: size of b
 *   @param[out] result: array of int
 * 
 *   @return k: size of result
*/
int intersect_int (int *a, int size_a, int *b, int size_b, int *result)
{
  int k = 0;
  for (int i = 0; i < size_a; i++){
    for (int j = 0; j < size_b; j++){
      if (a[i]==b[j]){
	      result[k] = a[i];
	      k++;
      }
    }
  }
  return k;
}

/**
 * Function: Waveform_alloc
 * ------------------------
 *  Alloc memory for Waveform
 * 
 * @param[in,out] wav: pointer to Waveform
 * @param[in] size: number of points
 * @param[in] name: name of waveform
*/void Waveform_alloc (Waveform **wav, const int size, const char *name)
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

/**
 * Function: Waveform_push
 * -----------------------
 * Realloc memory for Waveform changing size
 * 
 * @param[in,out] wav: pointer to Waveform
 * @param[in] size: new number of points
*/
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

/**
 * Function: Waveform_rmap
 * -----------------------
 *   Compute real/imag <-> amplitude/phase
 *   
 *   @param[in,out] h: pointer to Waveform
 *   @param[in] mode: 1: (Re, Im) -> (Amplitude, phase)
 *                    0: (Amplitude, phase) -> (Re, Im)
 *   @param[in] unw: 0: do not unwrap
 *                   1: unwrap
*/
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
    /* (Amplitude, phase) -> (Re, Im) */
    /* h =  A exp( -i phi) */
#pragma omp simd
    for (int i = 0; i < size; i++) 
      h->real[i] = + h->ampli[i] * cos(h->phase[i]);  
#pragma omp simd
    for (int i = 0; i < size; i++)
      h->imag[i] = - h->ampli[i] * sin(h->phase[i]); 
  }
}

/**
 * Function: Waveform_interp
 * -------------------------
 *   Interpolate Real/imag waveform to new time array
 * 
 *   @param[in,out] h: pointer to Waveform
 *   @param[in] size: new number of points
 *   @param[in] t0: initial time
 *   @param[in] dt: time step
 *   @param[in] name: name of waveform
*/
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
  interp_spline(h_aux->time, h_aux->real, h_aux->size, h->time, size, h->real);
  interp_spline(h_aux->time, h_aux->imag, h_aux->size, h->time, size, h->imag);

  /* Compute phase and amplitude */
  /* Waveform_rmap(h, 1, 1); */ /* unwrap */
  Waveform_rmap(h, 1, 0);  /* do not unwrap */

  /* Free aux memory */
  Waveform_free (h_aux);
}

/**
  * Function: Waveform_interp_ap
  * ----------------------------
  *   Interpolate amplitude/phase waveform to new time array
  * 
  *   @param[in,out] h: pointer to Waveform
  *   @param[in] size: new number of points
  *   @param[in] t0: initial time
  *   @param[in] dt: time step
  *   @param[in] name: name of waveform
  */
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
  interp_spline(h_aux->time, h_aux->ampli, h_aux->size, h->time, size, h->ampli);
  interp_spline(h_aux->time, h_aux->phase, h_aux->size, h->time, size, h->phase);
  
  /* Compute Real/Imag */
  Waveform_rmap (h, 0, 0); /* do not unwrap */
  
  /* Free aux memory */
  Waveform_free (h_aux);
}

/**
 * Function: Waveform_output
 * -------------------------
 *   Output Waveform to file
 *   The output file is named as: output_dir/wav->name.txt
 * 
 *   @param[in] wav: pointer to Waveform
*/
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

/**
 * Function: Waveform_free
 * ------------------------
 *   Free memory for Waveform
 * 
 *   @param[in,out] wav: pointer to Waveform
*/
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

/**
  * Function: WaveformFD_alloc
  * ---------------------------
  *  Alloc memory for WaveformFD
  * 
  * @param[in,out] wav: pointer to WaveformFD
  * @param[in] size: number of points
  * @param[in] name: name of waveform
*/
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

/**
 * Function: WaveformFD_push:
 * --------------------------
 *   Realloc memory for WaveformFD changing size
 * 
 *   @param[in,out] wav: pointer to WaveformFD
 *   @param[in] size: new number of points
*/
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

/**
 * Function: WaveformFD_rmap
 * -------------------------
 *   Compute real/imag <-> amplitude/phase
 *   
 *   @param[in,out] h: pointer to WaveformFD
 *   @param[in] mode: 1: (Re, Im) -> (Amplitude, phase)
 *                    0: (Amplitude, phase) -> (Re, Im)
 *   @param[in] unw: 0: do not unwrap
 *                   1: unwrap
*/
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

/**
 * Function: WaveformFD_interp_ap
 * ------------------------------
 *   Interpolate amp/phase waveform to new time array
 * 
 *   @param[in,out] h: pointer to WaveformFD
 *   @param[in] size: new number of points
 *   @param[in] f0: initial frequency
 *   @param[in] df: frequency step
 *   @param[in] name: name of waveform
*/
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
  interp_spline(h_aux->freq, h_aux->ampli, h_aux->size, h->freq, size, h->ampli);
  interp_spline(h_aux->freq, h_aux->phase, h_aux->size, h->freq, size, h->phase);
  
  /* Free aux memory */
  WaveformFD_free (h_aux);
}

/**
 * Function: WaveformFD_output
 * ---------------------------
 *   Output WaveformFD to file
 *   The output file is named as: output_dir/wav->name.txt
 * 
 *   @param[in] wav: pointer to WaveformFD
*/
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

/**
 * Function: WaveformFD_free
 * --------------------------
 *   Free memory for WaveformFD
 * 
 *   @param[in,out] wav: pointer to WaveformFD
*/
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

/**
  * Function: Waveform_lm_alloc
  * ---------------------------
  *  Alloc memory for Waveform_lm
  * 
  * @param[in,out] wav: pointer to Waveform_lm
  * @param[in] size: number of points
  * @param[in] name: name of waveform
  * @param[in] use_mode_lm: array of multipolar indices
  * @param[in] use_mode_lm_size: size of use_mode_lm
*/
void Waveform_lm_alloc (Waveform_lm **wav, int size, const char *name, int *use_mode_lm, int use_mode_lm_size)
{
  *wav = (Waveform_lm *) calloc(1, sizeof(Waveform_lm)); 
  if (wav == NULL)
    errorexit("Out of memory");
  (*wav)->size = size; 
  set_multipolar_idx_mask((*wav)->kmask,     KMAX, use_mode_lm, use_mode_lm_size, 0);
  set_multipolar_idx_mask((*wav)->kmask_nqc, KMAX, EOBPars->use_mode_lm_nqc, EOBPars->use_mode_lm_nqc_size, 0);
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

/**
 * Function: Waveform_lm_push
 * --------------------------
 *   Realloc memory for Waveform_lm changing size
 * 
 *   @param[in,out] wav: pointer to Waveform_lm
 *   @param[in] size: new number of points
*/
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

/**
 * Function: Waveform_lm_interp
 * ----------------------------
 *   Interpolate multipolar waveform to new time array
 * 
 *   @param[in,out] hlm: pointer to Waveform_lm
 *   @param[in] size: new number of points
 *   @param[in] t0: initial time
 *   @param[in] dt: time step
 *   @param[in] name: name of waveform
*/
void Waveform_lm_interp (Waveform_lm *hlm, const int size, const double t0, const double dt, const char *name)
{
  /* Alloc and init aux memory */  
  Waveform_lm *hlm_aux;
  const int oldsize = hlm->size;
  Waveform_lm_alloc(&hlm_aux, oldsize, "", EOBPars->use_mode_lm, EOBPars->use_mode_lm_size);
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
      interp_spline(hlm_aux->time, hlm_aux->ampli[k], hlm_aux->size, hlm->time, size, hlm->ampli[k]);
  for (int k = 0; k < KMAX; k++)
    if (hlm->kmask[k]) 
      interp_spline(hlm_aux->time, hlm_aux->phase[k], hlm_aux->size, hlm->time, size, hlm->phase[k]);
  
  /* Free aux memory */
  Waveform_lm_free (hlm_aux);
}

/**
 * Function: Waveform_lm_output
 * ----------------------------
 *   Output multipolar Waveform to file
 *   The output file is named as: output_dir/wav->name_l%d_m%d.txt
 *   
 *   @param[in] wav: pointer to Waveform_lm
*/
void Waveform_lm_output (Waveform_lm *wav)
{
  int kmask[KMAX];
  set_multipolar_idx_mask(kmask, KMAX, EOBPars->output_lm, EOBPars->output_lm_size, 0);  
  char fname[STRLEN*2];
  const int n = wav->size;
  
  for (int k=0; k<KMAX; k++) {
    if (kmask[k]) {      
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
/**
 * Function: Waveform_lm_output_reim
 * ---------------------------------
 *   Output multipolar Waveform to file as real/imag
 *   The output file is named as: output_dir/wav->name_l%d_m%d_reim.txt
 *   
 *   @param[in] wav: pointer to Waveform_lm
*/
void Waveform_lm_output_reim (Waveform_lm *wav)
{
  int kmask[KMAX];
  set_multipolar_idx_mask(kmask, KMAX, EOBPars->output_lm, EOBPars->output_lm_size, 0);  
  char fname[STRLEN*2];
  double re,im;
  const int n = wav->size;
  for (int k=0; k<KMAX; k++) {
    if (kmask[k]) {
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

/**
  * Function: Waveform_lm_extract
  * -----------------------------
  *   Extract a multipolar waveforms at times t >= to and t < tn, Alloc a new Waveform_lm var
  * 
  *   @param[in] hlma: pointer to Waveform_lm
  *   @param[in] to: initial time
  *   @param[in] tn: final time
  *   @param[in,out] hlmb: pointer to Waveform_lm
  *   @param[in] name: name of waveform
*/
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
  Waveform_lm_alloc (hlmb, N, name, EOBPars->use_mode_lm, EOBPars->use_mode_lm_size);
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

/**
  * Function: Waveform_lm_join
  * --------------------------
  *   Join two multipolar waveforms at t = to
  * 
  *   @note Time arrays are supposed to be ordered as
  *   hlma->time:  x x x x x x x x x 
  *   hlmb->time:       o o o o o o o o o 
  *   to        :                |
  *   But they do not need to overlap or be uniformly spaced.
  *   to can be 
  *   to > hlma->time[hlma->size-1] => extend the a waveform
  *   to < hlmb->time[0]            => join the whole b waveform
  *   A number of checks enforce the above structure, if possible.
  * 
  *   @param[in,out] hlma: pointer to Waveform_lm
  *   @param[in,out] hlmb: pointer to Waveform_lm
  *   @param[in] to: time to join
*/
void Waveform_lm_join (Waveform_lm *hlma, Waveform_lm *hlmb, double to)
{
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

/**
 * Function: Waveform_lm_free
 * ---------------------------
 *   Free memory for Waveform_lm
 * 
 *   @param[in,out] wav: pointer to Waveform_lm
*/
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

/**
  * Function: WaveformFD_lm_alloc
  * ------------------------------
  *  Alloc memory for WaveformFD_lm
  * 
  * @param[in,out] wav: pointer to WaveformFD_lm
  * @param[in] size: number of points
  * @param[in] name: name of waveform
*/
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

/**
 * Function: WaveformFD_lm_push
 * -----------------------------
 *   Realloc memory for WaveformFD_lm changing size
 * 
 *   @param[in,out] wav: pointer to WaveformFD_lm
 *   @param[in] size: new number of points
*/
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

/**
 * Function: WaveformFD_lm_interp_ap
 * ---------------------------------
 *   Interpolate amp/phase waveform to new freq array
 *   @note each mode has its own frequency array hlm->F[k]
 *   @note the uniform array is stored in hlm->freq
 *   @note use interp_spline_checklim() and set to zero points outside interp range
 * 
 *   @param[in,out] hlm: pointer to WaveformFD_lm
 *   @param[in] size: new number of points
 *   @param[in] f0: initial frequency
 *   @param[in] df: frequency step
 *   @param[in] name: name of waveform
*/
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
      /* interp_spline_omp(hlm_aux->F[k], hlm_aux->ampli[k], hlm_aux->size, hlm->freq, size, hlm->ampli[k]); */
      interp_spline_checklim(hlm_aux->F[k], hlm_aux->ampli[k], hlm_aux->size, hlm->freq, size, hlm->ampli[k]);
    }
  }
  for (int k = 0; k < KMAX; k++){
    if (hlm->kmask[k]) {
      /* interp_spline_omp(hlm_aux->F[k], hlm_aux->phase[k], hlm_aux->size, hlm->freq, size, hlm->phase[k]); */
      interp_spline_checklim(hlm_aux->F[k], hlm_aux->phase[k], hlm_aux->size, hlm->freq, size, hlm->phase[k]);
    }
  }
  /* Free aux memory */
  WaveformFD_lm_free (hlm_aux);
}

/**
 * Function: WaveformFD_lm_interp_ap_freqs
 * ---------------------------------------
 *   Interpolate amp/phase waveform to new freq array
 *   @note each mode has its own frequency array hlm->F[k]
 *   @note the uniform array is stored in hlm->freq
 *   @note use interp_spline_checklim() and set to zero points outside interp range
 * 
 *   @param[in,out] hlm: pointer to WaveformFD_lm
 *   @param[in] name: name of waveform
*/
void WaveformFD_lm_interp_ap_freqs (WaveformFD_lm *hlm, const char *name)
{
  /* Alloc and init aux memory */  
  WaveformFD_lm *hlm_aux;
  int size = EOBPars->freqs_size;
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
    hlm->freq[i] = EOBPars->freqs[i];
  
  /* Interp */
  for (int k = 0; k < KMAX; k++){
    if (hlm->kmask[k]){
      /* interp_spline(hlm_aux->F[k], hlm_aux->ampli[k], hlm_aux->size, hlm->freq, size, hlm->ampli[k]); */
      interp_spline_checklim(hlm_aux->F[k], hlm_aux->ampli[k], hlm_aux->size, hlm->freq, size, hlm->ampli[k]);
    }
  }
  for (int k = 0; k < KMAX; k++){
    if (hlm->kmask[k]) {
      /* interp_spline(hlm_aux->F[k], hlm_aux->phase[k], hlm_aux->size, hlm->freq, size, hlm->phase[k]); */
      interp_spline_checklim(hlm_aux->F[k], hlm_aux->phase[k], hlm_aux->size, hlm->freq, size, hlm->phase[k]);
    }
  }
  /* Free aux memory */
  WaveformFD_lm_free (hlm_aux);

}

/**
 * Function: WaveformFD_lm_output
 * ------------------------------
 *  Output multipolar Waveform to file
 *  The output file is named as: output_dir/FD_wav->name_l%d_m%d.txt
 * 
 *  @param[in] wav: pointer to WaveformFD_lm
*/
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

/**
 * Function: WaveformFD_lm_output_reim
 * -----------------------------------
 *   Output multipolar Waveform to file as real/imag
 *   The output file is named as: output_dir/FD_wav->name_l%d_m%d_reim.txt
 *   
 *   @param[in] wav: pointer to WaveformFD_lm
*/
void WaveformFD_lm_output_reim (WaveformFD_lm *wav)
{
  int kmask[KMAX];
  set_multipolar_idx_mask(kmask, KMAX, EOBPars->output_lm, EOBPars->output_lm_size, 0);  
  char fname[STRLEN*2];
  double re,im;
  const int n = wav->size;
  for (int k=0; k<KMAX; k++) {
    if (kmask[k]) {
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

/**
 * Function: WaveformFD_lm_free
 * ----------------------------
 *  Free memory for WaveformFD_lm
 * 
 *  @param[in,out] wav: pointer to WaveformFD_lm
*/
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

/**
 * Function: Waveform_lm_t_alloc
 * -----------------------------
 *   Alloc memory for Waveform_lm_t (waveform at time point)
 * 
 *   @param[in,out] wav: pointer to Waveform_lm_t
*/
void Waveform_lm_t_alloc (Waveform_lm_t **wav)
{
  *wav = (Waveform_lm_t *) calloc(1, sizeof(Waveform_lm_t)); 
  if (wav == NULL) errorexit("Out of memory");
  (*wav)->time = 0.;
  (*wav)->freq = 0.;
  set_multipolar_idx_mask ((*wav)->kmask, KMAX, EOBPars->use_mode_lm, EOBPars->use_mode_lm_size, 0); 
  set_multipolar_idx_mask ((*wav)->kmask_nqc, KMAX, EOBPars->use_mode_lm_nqc, EOBPars->use_mode_lm_nqc_size, 0);
}

/**
 * Function: Waveform_lm_t_free
 * ----------------------------
 *   Free memory for Waveform_lm_t
 * 
 *   @param[in,out] wav: pointer to Waveform_lm_t
*/
void Waveform_lm_t_free (Waveform_lm_t *wav)
{
  if (!wav) return;
  free(wav);
}

/**
 * Function: Dynamics_alloc
 * ------------------------
 *   Alloc memory for Dynamics
 * 
 *   @param[in,out] dyn: pointer to Dynamics
 *   @param[in] size: number of points
 *   @param[in] name: name of waveform
*/
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
  (*dyn)->spins = NULL;
}

/**
 * Function: Dynamics_push
 * -----------------------
 *   Realloc memory for Dynamics changing size
 * 
 *   @param[in,out] dyn: pointer to Dynamics
 *   @param[in] size: new number of points
*/
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

/**
 * Function: Dynamics_interp
 * -------------------------
 *   Interp and overwrite the EOB dynamics
 * 
 *   @param[in,out] dyn: pointer to Dynamics
 *   @param[in] size: new number of points
 *   @param[in] t0: initial time
 *   @param[in] dt: time step
 *   @param[in] name: name of dynamics
*/
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
    interp_spline(dyn_aux->time, dyn_aux->data[k], dyn_aux->size, dyn->time, size, dyn->data[k]);

  /* Free aux memory */
  Dynamics_free (dyn_aux);
}

/**
 * Function: Dynamics_output
 * -------------------------
 *   Output Dynamics to file
 *   The output file is named as: output_dir/dyn->name.txt
 * 
 *   @param[in] dyn: pointer to Dynamics
*/
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

/**
 * Function: Dynamics_extract
 * --------------------------
 *   Extract Dynamics at times t >= to and t < tn, Alloc a new Dynamics var
 * 
 *   @param[in] dyna: pointer to Dynamics
 *   @param[in] to: initial time
 *   @param[in] tn: final time
 *   @param[in,out] dynb: pointer to Dynamics
 *   @param[in] name: name of dynamics
*/
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

/**
  * Function: Dynamics_join
  * -----------------------
  *   Join two dynamics time series at t = to
  * 
  *   @note Time arrays are suppose to be ordered as
  *   dyna->time:  x x x x x x x x x 
  *   dynb->time:       o o o o o o o o o 
  *   to        :                |
  *   But they do not need to overlap or be uniformly spaced.
  *   Note to can be 
  *   to > dyna->time[dyna->size-1] => extend the dynamics data
  *   to < dynb->time[0]            => join the whole b dynamics
  *   Following checks enforce the above structure, if possible.
  * 
  *   @param[in,out] dyna: pointer to Dynamics
  *   @param[in,out] dynb: pointer to Dynamics
  *   @param[in] to: time to join
*/
void Dynamics_join (Dynamics *dyna, Dynamics *dynb, double to)
{
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

/**
 * Function: Dynamics_free
 * -----------------------
 *   Free memory for Dynamics
 * 
 *   @param[in,out] dyn: pointer to Dynamics
*/
void Dynamics_free (Dynamics *dyn)
{
  if (!dyn) return;
  if (dyn->time) free(dyn->time);
  for (int v = 0; v < EOB_DYNAMICS_NVARS; v++)
    if (dyn->data[v]) free(dyn->data[v]);
  if (dyn->spins) dyn->spins = NULL;
  free(dyn);
}

/**
 * Function: DynamicsSpin_alloc
 * ----------------------------
 *   Alloc memory for DynamicsSpin
 *  
 *   @param[in,out] dyn: pointer to DynamicsSpin
*/
void DynamicsSpin_alloc (DynamicsSpin **dyn, int size)
{
  *dyn = (DynamicsSpin *) calloc(1, sizeof(DynamicsSpin)); 
  if (dyn == NULL) errorexit("Out of memory");
  (*dyn)->size = size; 
  (*dyn)->time = malloc (size * sizeof(double));
  memset((*dyn)->time, 0, size*sizeof(double));
  for (int v=0; v<EOB_EVOLVE_SPIN_NVARS; v++) {
    (*dyn)->data[v] = malloc (size * sizeof(double));  
    memset((*dyn)->data[v], 0, size*sizeof(double));
    (*dyn)->spline[v] = gsl_spline_alloc (gsl_interp_cspline, size);
    (*dyn)->accel[v]  = gsl_interp_accel_alloc();

  }
  (*dyn)->omg_stop=-1; // set from EOBPars or by NR merger
  (*dyn)->t_stop=-1; // use Momg as stopping criterion, if not otherwise specified.
  (*dyn)->time_backward=0;
  (*dyn)->omg_backward= 0; 
}

/**
 * Function: DynamicsSpin_push
 * ---------------------------
 *  Realloc memory for DynamicsSpin changing size
 * 
 *  @param[in,out] dyn: pointer to DynamicsSpin
 *  @param[in] size: new number of points
*/
void DynamicsSpin_push (DynamicsSpin **dyn, int size)
{
  const int n  = (*dyn)->size;
  const int dn = size - (*dyn)->size;
  (*dyn)->time = realloc ( (*dyn)->time, size * sizeof(double) );
  for (int v = 0; v < EOB_EVOLVE_SPIN_NVARS; v++) {
    (*dyn)->data[v] = realloc ( (*dyn)->data[v], size * sizeof(double) );
    if ((*dyn)->data[v] == NULL) errorexit("Out of memory.");
    gsl_spline_free ((*dyn)->spline[v]);
    (*dyn)->spline[v] = gsl_spline_alloc (gsl_interp_cspline, size);
    gsl_interp_accel_free((*dyn)->accel[v]);
    (*dyn)->accel[v]  = gsl_interp_accel_alloc();
  }
  (*dyn)->size = size; 
}

/**
 * Function: DynamicsSpin_free
 * ---------------------------
 *   Free memory for DynamicsSpin
 * 
 *   @param[in,out] dyn: pointer to DynamicsSpin
*/
void DynamicsSpin_free (DynamicsSpin *dyn)
{
  if (!dyn) return;
  if (dyn->time) free(dyn->time);
  if (dyn->data)
    for (int v=0; v<EOB_EVOLVE_SPIN_NVARS; v++){
      if (dyn->data[v]) free(dyn->data[v]);
      gsl_spline_free (dyn->spline[v]);
      gsl_interp_accel_free(dyn->accel[v]);
    }
  free(dyn);
}

/**
 * Function: DynamicsSpin_output
 * -----------------------------
 *   Output DynamicsSpin to file
 *   The output file is named as: output_dir/dynspin.txt
 * 
 *   @param[in] dyn: pointer to DynamicsSpin
*/
void DynamicsSpin_output (DynamicsSpin *dyn)
{  
  char fname[STRLEN*2];
  const int n = dyn->size;
  sprintf(fname,"%s/dynspin.txt",EOBPars->output_dir);
  FILE* fp;
  if ((fp = fopen(fname, "w+")) == NULL) errorexits("error opening file",fname);
  fprintf(fp, "# t "); 
  for (int v=0; v<EOB_EVOLVE_SPIN_NVARS; v++)
    fprintf(fp, " %s",  eob_prec_var[v]); 
  fprintf(fp, "\n"); 
  for (int i = 0; i < n; i++) {
    fprintf(fp, "%.9e", dyn->time[i]);
    for (int v=0; v<EOB_EVOLVE_SPIN_NVARS; v++)
      fprintf(fp, " %.16e", dyn->data[v][i]); 
    fprintf(fp, "\n");
  }
  fclose(fp);
}


/**
 * Function: DynamicsSpin_join
 * ---------------------------
 *   Join two dynamics time series at t = to
 * 
 *   @note Time arrays are suppose to be ordered as
 *   dyna->time:  x x x x x x x x x 
 *   dynb->time:       o o o o o o o o o 
 *   to        :                |
 *   But they do not need to overlap or be uniformly spaced.
 *   Note to can be 
 *   to > dyna->time[dyna->size-1] => extend the dynamics data
 *   to < dynb->time[0]            => join the whole b dynamics
 *   Following checks enforce the above structure, if possible.
 * 
 *   @param[in,out] dyna: pointer to DynamicsSpin
 *   @param[in,out] dynb: pointer to DynamicsSpin
 *   @param[in] to: time to join
*/
void DynamicsSpin_join (DynamicsSpin *dyna, DynamicsSpin *dynb, double to)
{
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
  DynamicsSpin_push (&dyna, N);

  /* Copy the relevant part of b into a */
#pragma omp simd
  for (int i = 0; i < Nb; i++) 
    dyna->time[ioa + i] = dynb->time[iob + i]; 
  
  for (int v = 0; v < EOB_EVOLVE_SPIN_NVARS; v++) {
#pragma omp simd
    for (int i = 0; i < Nb; i++) {
      dyna->data[v][ioa + i] = dynb->data[v][iob + i];
    }
 }
  for(int v=0; v < EOB_EVOLVE_SPIN_NVARS; v++)
    gsl_spline_init (dyna->spline[v], dyna->time, dyna->data[v], dyna->size);   

  

#if (0)
  printf("Dynamics (a) i = %d time[i] = %.6e \n",ioa, dyna->time[ioa]);
  printf("Dynamics (b) i = %d time[i] = %.6e \n",iob, dynb->time[iob]);
  printf("Total size (a)+(b) = %d + %d = %d (%d)\n",Nb, iob, N, dyna->size);
  printf("%.6e - %.6e = %.6e\n",dyna->time[dyna->size-1], dynb->time[dynb->size-1],dyna->time[dyna->size-1]-dynb->time[dynb->size-1]);
#endif 
    
}

/**
 * Function: DynamicsSpin_set_params
 * ---------------------------------
 *   Sync some quick access parameters in dyn with parameter database
 *   @note to be used carefully
 * 
 *   @param[in,out] dyn: pointer to Dynamics
*/
void Dynamics_set_params (Dynamics *dyn)
{
  dyn->store = 0;
  dyn->size  = EOBPars->size;
  dyn->dt        = EOBPars->dt;
  dyn->t_stop    = EOBPars->ode_tmax;
  for (int l=2; l<6; l++) {
    dyn->dress_tides_fmode_A[l] = dyn->dress_tides_fmode_B[l] = 1;
    dyn->dress_tides_fmode_A_u[l] = dyn->dress_tides_fmode_B_u[l] = 0;
  }

}

/**
 * Function: NQCdata_alloc
 * -----------------------
 *   Alloc memory for NQCdata
 * 
 *   @param[in,out] nqc: pointer to NQCdata
*/
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

/**
 * Function: NQCdata_free
 * ----------------------
 *   Free memory for NQCdata
 * 
 *   @param[in,out] nqc: pointer to NQCdata
*/
void NQCdata_free (NQCdata *nqc)
{
  if (!nqc) return;
  if (nqc->flx) free (nqc->flx);
  if (nqc->hlm) free (nqc->hlm);
  if (nqc)      free (nqc);
}
  

/**
 * Function: lal_tetrad_conventions
 * --------------------------------
 *   Enforce the LALSimulation tetrad conventions.
 *   With respect to SXS, hlm->-1*hlm.
 *   This implies that the +,x polarizations are mapped in
 *   -1*(h+,hx)
 *
 *   @param[in,out] hpc: pointer to the +,x polarizations
*/
void lal_tetrad_conventions(Waveform *hpc)
{
  for(int i=0; i<hpc->size; i++){
    hpc->real[i] *= -1.;
    hpc->imag[i] *= -1.;
  }
}

/**
 * Function: time_shift_mrg_to_0
 * -----------------------------
 *   Find the time shift to move the merger to t=0
 *   Define "merger" as the last peak of
 *   Q = sqrt(A); A = sum_{ell, emm} |A_{ell,emm}|^2
 * 
 *   @param[in,out] hlm: pointer to coprecessing hlm  
*/
double time_shift_mrg_to_0(Waveform_lm *hlm)
{
  int i, imrg;
  double Q, Qmrg;
  Qmrg  = 0;
  imrg  = 0;

#if(0)
  /* Scan the entire waveform */
  for (i=0; i<hlm->size; i++){
    Q = 0.;
    for (int k=0; k<KMAX; k++){
      if (hlm->kmask[k]) 
        Q += hlm->ampli[k][i]*hlm->ampli[k][i];

    if (Q > Qmrg){
      Qmrg = Q;
      imrg = i;
    }
  }
#endif

  /* Start looking for the last peak from the end */
  Qmrg  = 0.;
  imrg  = hlm->size-1;
  for (i=hlm->size-1; i>0; i--){
    Q = 0.;
    for (int k=0; k<KMAX; k++){
      if (hlm->kmask[k]) 
        Q += hlm->ampli[k][i]*hlm->ampli[k][i];
    }
    if (Q < Qmrg){
      imrg = i;
      Qmrg = Q;
      break;
    } else {
      Qmrg = Q;
    }
  }
  return hlm->time[imrg];

}

/**
 * Function: time_shift_TD
 * ------------------------
 *  Time shift the waveform polarizations
 * 
 *   @param[in,out] t: pointer to Waveform->t or Waveformlm->t
 *   @param[in] tc: time shift
*/
void time_shift_TD(double *t, double tc, int size)
{
  for(int i=0; i<size; i++){
    t[i] -= tc;
  }
}


/**
 * Function: time_shift_FD
 * -----------------------
 *   Perform a time shift in FD
 * 
 *   @param[in,out] hpc: pointer to WaveformFD
 *   @param[in] tc: time shift
*/
void time_shift_FD(WaveformFD *hpc, double tc){

  double ts = -tc;
  double twopits = 2*Pi*ts;
  //printf("Time shifting\n");
  //printf("%.2e %.2e %. \n", tc, EOBPars->df, hpc->size);

  for(int i=0; i < hpc->size; i++){
    double f = hpc->freq[i];
    double cn = cos(-twopits*f);
    double sn = sin(-twopits*f);
    /** Shift: h-->h*e^{-2 Pi ts f} */
    double tmpr = hpc->preal[i]*cn + hpc->pimag[i]*sn;
    double tmpi = hpc->pimag[i]*cn - hpc->preal[i]*sn;
    hpc->preal[i] = tmpr;
    hpc->pimag[i] = tmpi;

    tmpr  = hpc->creal[i]*cn + hpc->cimag[i]*sn;
    tmpi  = hpc->cimag[i]*cn - hpc->creal[i]*sn;
    hpc->creal[i] = tmpr;
    hpc->cimag[i] = tmpi;
  }
}

/**
  * Function: time_units_factor
  * ---------------------------
  *   Convert time in sec to dimensionless and mass-rescaled units
  * 
  *   @param[in] M: mass
  *   @return conversion factor
  */
double time_units_factor(double M)
{
  return 1./(M*MSUN_S);
}
/**
 * 
  * Function: time_units_conversion
  * -------------------------------
  *   Convert time in sec to dimensionless and mass-rescaled units
  * 
  *   @param[in] M: mass
  *   @param[in] t: time in sec
  *   @return time in dimensionless and mass-rescaled units
  */
double time_units_conversion(double M, double t)
{
  return t/(M*MSUN_S);
}

/**
  * Function: radius0
  * -----------------
  *   Convert frequency in Hz to dimensionless radius using Kepler law
  * 
  *   @param[in] M: mass
  *   @param[in] fHz: frequency in Hz
  *   @return dimensionless radius
  */
double radius0(double M, double fHz)
{
  double x = (M*fHz*MSUN_S*2.*Pi)/2.;
  return cbrt(1/(x*x));
}

/**
  * Function: system_mkdir
  * ----------------------
  *   Make dir
  * 
  *   @param[in] name: name of dir
  *   @return 0 if success
*/
int system_mkdir(const char *name)
{
  char s[STRLEN];
  sprintf(s,"mkdir -p %s",name);
  return system(s);
}

/**
  * Function: is_blank
  * ------------------
  *   Check if string is blank
  * 
  *   @param[in] line: string
  *   @return 1 if blank
*/
int is_blank(const char *line) 
{
  const char accept[]=" \t\r\n"; 
  return (strspn(line, accept) == strlen(line));
}

/**
  * Function: remove_white_spaces
  * -----------------------------
  *   Remove white spaces from string
  * 
  *   @param[in,out] str: string
*/
void remove_white_spaces(char *str)
{
  int i = 0, j = 0;
  while (str[i]) {
    if (str[i] != ' ')
      str[j++] = str[i];
    i++;
  }
  str[j] = '\0';
}

/**
  * Function: remove_comments
  * --------------------
  *   Cut string to the first delimiter
  * 
  *   @param[in,out] line  : string
  *   @param[in] delimiters: delimiter
*/
void remove_comments(char *line, const char *delimiters)
{
  int sz = strcspn(line,delimiters);
  char *newline = (char *) calloc(sizeof(char), sz+1);
  strncpy(newline,line,sz);
  newline[sz] = '\0';
  strcpy(line, newline);
  free(newline);
}

/**
  * Function: trim
  * --------------
  *   Get rid of trailing and leading whitespace
  * 
  *   @param[in,out] str: string
  *   @return pointer to beginning of string
*/
char *trim(char *str)
{
  char *start = str;
  char *end = str + strlen(str);  
  while(*start && isspace(*start))
    start++;
  while(end > start && isspace(*(end - 1)))
    end--;
  *end = '\0';
  return start;
}

/**
 * Function: startswith
 * --------------------
 *  Check if string starts with beg
 * 
 *  @param[in] str: string
 *  @param[in] beg: beginning of string
 * 
 *  @return 1 if true
*/
int startswith(const char *str, const char *beg)
{
  return (strncmp(beg, str, strlen(beg)) == 0);
}
/**
 * Function: endswith
 * ------------------
 *  Check if string ends with end
 * 
 *  @param[in] str: string
 *  @param[in] end: end of string
 * 
 *  @return 1 if true
*/
int endswith(const char *str, const char *end)
{
  return (strncmp(str+strlen(str)-strlen(end), end, strlen(end)) == 0);
}

/**
 * Function: getkv
 * ---------------
 *  Get key value separated by a '='
 * 
 *  @param[in] line: string
 *  @param[in,out] key: key
 *  @param[in,out] val: value
 * 
 *  @return 0 if success
*/
int getkv(char *line, char **key, char **val) 
{
  char *ptr = strchr(line,'=');
  if (ptr == NULL)
    return 1;
  *ptr++ = '\0';
  *key = trim(line);
  *val = trim(ptr);
  return 0;
}

/**
 * Function: is_string
 * -------------------
 *   Identify string in parfile
 * 
 *   @param[in] str: string
 * 
 *   @return 1 if true
*/
int is_string(const char *str)
{
  if (startswith(str,"\"") && endswith(str,"\""))
    return 1;
  return 0;
}

/**
 * Function: string_trim
 * ---------------------
 *   Trim string
 * 
 *   @param[in] str: string
 * 
 *   @return pointer to beginning of string
*/
char *string_trim(char *str)
{
  char *start = str;
  char *end = str + strlen(str);  
  if (is_string(str)) {
    start++;
    end--;
    *end = '\0';
    return start;
  }
  return NULL;
}

/* Helpers conversion str to type */

/**
 * Function: par_get_i
 * -------------------
 *   Convert string to int
 * 
 *   @param[in] val: string
 * 
 *   @return int
*/
int par_get_i (char *val)
{
  return atoi(val);
}

/**
 * Function: par_get_b
 * -------------------
 *   Convert string to bool
 * 
 *   @param[in] val: string
 * 
 *   @return bool
*/
int par_get_b (char *val)
{
  return atoi(val)?1:0;
}
/**
 * Function: par_get_d
 * -------------------
 *   Convert string to double
 * 
 *   @param[in] val: string
 * 
 *   @return double
*/
double par_get_d (char *val)
{
  return atof(val);
}
/**
  * Function: par_get_s
  * -------------------
  *   Copy string
  * 
  *   @param[in,out] dest: destination
  *   @param[in] src: source
  * 
  *   @return 0 if successful, 1 src is not a string
  */
int par_get_s (char * dest, char *src)
{
  if (is_string(src)) {    
    strncpy(dest, src+1, strlen(src)-2);
    return 0;
  }
  return 1;
}

/** return number of entries in a string according to "delimiters" */ 
const char * ARRAY_DELIMITER[] = { "[", ",", "]" };
/**
 * Function: noentries
 * -------------------
 *   Return number of entries in a string according to "delimiters"
 * 
 *   @param[in] string: string
 * 
 *   @return number of entries
*/
int noentries(const char *string)
{
  int n = 0;
  char *s, *t; 
  int len = strlen(string);
  
  if (! ((strncmp(ARRAY_DELIMITER[0], string, 1) == 0) &&
	 (strncmp(string+len-1, ARRAY_DELIMITER[2], 1) == 0)) ) {
    /* printf(" %c...%c\n",string[0],string[len-1]);  */
    errorexit("Parsing string for array: must start with [ and end with ].");
  }

  s = strndup(string+1,len-2);
  t = strtok(s,ARRAY_DELIMITER[1]);
  while (t != NULL) {
    n++;
    t = strtok(NULL,ARRAY_DELIMITER[1]);
  }
  free(s);
  return n;
}

/**
 * Function: str2iarray
 * --------------------
 *   Convert a string to an array of int (alloc mem)
 * 
 *   @param[in] string: string
 *   @param[in,out] a: pointer to array
 * 
 *   @return number of entries
*/
int str2iarray(const char *string, int **a)
{
  char *s, *t;
  int n = noentries(string);
  (*a) = (int *) malloc (n * sizeof(int));
  if (!(*a)) errorexit(eob_error_msg[ERROR_OUTOFMEM]);
  const int len = strlen(string);
  s = strndup(string+1,len-2);
  t = strtok(s,ARRAY_DELIMITER[1]);
  int j = 0;
  while ((t != NULL) && (j<n)) {
    (*a)[j] = atoi(t);     
    t = strtok(NULL,ARRAY_DELIMITER[1]);
    j++;
  }
  free(s);  
  return n;
}

/**
 * Function: str2darray
 * --------------------
 *   Convert a string to an array of double (alloc mem)
 * 
 *   @param[in] string: string
 *   @param[in,out] a: pointer to array
 * 
 *   @return number of entries
*/
int str2darray(const char *string, double **a)
{
  char *s, *t;
  int n = noentries(string);
  (*a) = (double*) malloc (n * sizeof(double));
  if (!(*a)) errorexit(eob_error_msg[ERROR_OUTOFMEM]);
  const int len = strlen(string);
  s = strndup(string+1,len-2);
  t = strtok(s,ARRAY_DELIMITER[1]);
  int j = 0;
  while ((t != NULL) && (j<n)) {
    (*a)[j] = atof(t);
    t = strtok(NULL,ARRAY_DELIMITER[1]);
    j++;
  }
  free(s);  
  return n;
}

/**
 * Function: print_date_time
 * -------------------------
 *   Print date and time
*/
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
/**
 * Function: errorexit
 * -------------------
 *   Print error message and exit
 * 
 *   @param[in] file: file name
 *   @param[in] line: line number
 *   @param[in] s: error message
*/
void errorexit(char *file, int line, const char *s) 
{
  fprintf(stderr, "Error: %s  ", s);
  fprintf(stderr, "(%s, line %d)\n", file, line);
  fflush(stdout);
  fflush(stderr);
  exit(ERROR);
}

/**
 * Function: errorexits
 * --------------------
 *   Print error message and exit
 * 
 *   @param[in] file: file name
 *   @param[in] line: line number
 *   @param[in] s: error message
 *   @param[in] t: error message
*/
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
