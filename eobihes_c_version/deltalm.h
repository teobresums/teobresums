//
//  deltalm.h
//
//  Created by Philipp Fleig on 17/03/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef deltalm_h
#define deltalm_h

#include <gsl/gsl_math.h>

vector<double> EOBdeltalm(const double Hreal,const double Omega,const double nu){


/*
EOBdeltalm Residual phase corrections delta_{lm} up to l=m=5.

   EOBdeltalm(Hreal,Omega,nu, EOBopt);


   nu=0, the delta_lmm are written in Taylor-expanded form and all terms
   up to 4.5PN accuracy are included.
   nu~=0,
         (i) the 4.5PN, test-mass terms are excluded because of the too
         large PN-gap with the nu-dependent terms
         (ii) the deltalmm are replaced by suitable Pade' approximants for
         multipoles (2,2), (2,1), (3,3), (3,1)
   The l=m=2 residual phase includes the 3.5PN (nu-dependent) correction
   obtained in Faye et al.

   Reference(s)
    Damour, Iyer & Nagar, PRD 79, 064004 (2008)
    Fujita & Iyer, PRD 82 044051 (2010)
    Faye et al., Class. Q. Grav. 29 175004 (2012)
    Damour, Nagar & Bernuzzi, PRD 87, 084035 (2013)


 TODO: this routine requires optimization
 - precompute coefficients c(nu)
 - evaluate efficiently polynomials
*/

// Shorthands
double pi2    = pi*pi;
double nu2    = nu*nu;

double y      = gsl_pow_int(Hreal*Omega,2./3.);
double sqrt_y = sqrt(y);
double y3     = y*y*y;
double y32    = Hreal*Omega;

// Residual phases in Pade-resummed form when possible

// Leading order contributions
double delta22LO = 7./3.   * y32;
double delta21LO = 2./3.   * y32;
double delta33LO = 13./10. * y32;
double delta31LO = 13./30. * y32;
    
int kmax=35;
    
vector<double> deltalm(35);
    
for (int i=kmax; i--; ) {
    deltalm[i]=0.;
}

// l=2 ------------------------------------------------------------------

// Pade(2,2) approximant
double num        = (808920.*nu*pi*sqrt(y) + 137388.*pi2*y + 35.*nu2*(136080. + (154975. - 1359276.*nu)*y));
double den        = (808920.*nu*pi*sqrt(y) + 137388.*pi2*y + 35.*nu2*(136080. + (154975. + 40404.*nu)*y));
deltalm[1] = delta22LO*num/den;

// Pade(1,2) approximant
num        = 69020.*nu + 5992.*pi*sqrt_y;
den        = 5992.*pi*sqrt_y + 2456.*nu*(28.+493.*nu*y);
deltalm[0] = delta21LO*num/den;

// l=3 ------------------------------------------------------------------

// Pade(1,2) approximant
num        = 1.   + 94770.*pi/(566279.*nu)*sqrt_y;
den        = num + 80897.*nu/3159.*y;
deltalm[4] = delta33LO*num/den;

// Taylor-expanded form
deltalm[3] = (10.+33.*nu)/(15.*(1.-3.*nu))*y32    +  52./21.*pi*y3;

// Pade(1,2) approximant
num        = 4641.*nu + 1690.*pi*sqrt_y;
den        = num + 18207.*nu2*y;
deltalm[2] = delta31LO*num/den;

// l=4 ------------------------------------------------------------------

deltalm[8] =  (112.+219.*nu)/(120.*(1.-3.*nu)) * y32 + 25136./3465.*pi*y3;
deltalm[7] = (486.+4961.*nu)/(810.*(1.-2.*nu)) * y32 +   1571./385.*pi*y3;
deltalm[6] =  7.*(1.+6.*nu)/(15.*(1.-3.*nu))  * y32 +  6284./3465.*pi*y3;
deltalm[5] =   (2.+507.*nu)/(10.*(1.-2.*nu))   * y32 +  1571./3465.*pi*y3;

// l=5 ------------------------------------------------------------------

deltalm[9] = (96875. + 857528.*nu)/(131250.*(1.-2.*nu)) * y32;
    
return deltalm;

}

#endif /* deltalm_h */
