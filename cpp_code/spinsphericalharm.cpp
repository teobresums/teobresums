#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include "TEOBResum.h"

/* factorial */
double fact(int n){
  double f[] = {1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880.,
                3628800., 39916800., 479001600., 6227020800., 87178291200.};
  if (n < 0){
   fprintf(stderr, " computing a negative factorial. \n");
   exit(1);
  } else if (n <= 14){
   return f[n];
  } else {
   return n*fact(n-1);
  }
}

/* Wigner d-function */
double wigner_d_function(int l, int m, int s, double i){
  double dWig = 0.;
  
  double costheta = cos(i*0.5);
  double sintheta = sin(i*0.5);
 
  int ki = std::max( 0  , m-s );
  int kf = std::min( l+m, l-s );
  
  for( int k = ki; k <= kf; k++ ){
    dWig += 
      ( pow(-1.,k) * pow(costheta,2*l+m-s-2*k) * pow(sintheta,2*k+s-m) )/
      ( fact(k) * fact(l+m-k) * fact(l-s-k) * fact(s-m+k) );
  }

  return (sqrt(fact(l+m) * fact(l-m) * fact(l+s) * fact(l-s)) * dWig);
}

/* spin-weighted spherical harmonic */
void spinsphericalharm(double *rY, double *iY, int s, int l, int m, double phi, double i){
  /* Following the Ref.: https://arxiv.org/pdf/0709.0093.pdf */
  if ((l<0) || (m<-l) || (m>l)) {
    fprintf(stderr, " wrong l (%d) or m (%d) inside spinspharmY\n", l, m);
    exit(1);
  }

  double c = pow(-1.,-s) * sqrt( (2.*l+1.)/(4.*M_PI) );
  double dWigner = c * wigner_d_function(l,m,-s,i); 

  *rY = cos((double)(m)*phi) * dWigner;
  *iY = sin((double)(m)*phi) * dWigner;
}
