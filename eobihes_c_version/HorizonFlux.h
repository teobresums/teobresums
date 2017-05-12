//
//  HorizonFlux.h
//
//  Created by Philipp Fleig on 13/03/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef HorizonFlux_h
#define HorizonFlux_h

#include <gsl/gsl_math.h>

double HorizonFlux(const double x,const double Heff,const double jhat,const double nu){

//EOBFluxHorizon Compute horizon-absorbed fluxes.
//
//   [hatF FlmH] = EOBFluxHorizon(x,Heff,jhat,nu, EOBopt)
//
//   Nagar & Akcay, PRD 85, 044025 (2012)
//   Bernuzzi, Nagar & Zenginoglu, PRD 86, 104038 (2012)

const int kmax   = 2;
vector<double> rhoHlm(kmax);
vector<double> FlmHLO(kmax);
vector<double> FlmH(kmax);

// Shorthands
double nu2 = nu*nu;
double nu3 = nu*nu2;

double x2  = x*x;
double x3  = x*x2;
double x4  = x*x3;
double x5  = x*x4;
double x9  = x4*x5;
double x10 = x*x9;

const int k22 = 1;//LM2K(2,2);
const int k21 = 0;//LM2K(2,1);

// The Newtonian asymptotic contribution
const double FNewt22 = 32./5.*x5;


// Compute leading-order part (nu-dependent)
FlmHLO[1] = 32./5.*(1-4*nu+2*nu2)*x9;
FlmHLO[0] = 32./5.*(1-4*nu+2*nu2)*x10;


// Compute rho_lm
vector<double> c1(kmax);
vector<double> c2(kmax);
vector<double> c3(kmax);
vector<double> c4(kmax);
    
c1[1] = (4.-21.*nu + 27.*nu2 - 8.*nu3)/(4.*(1.-4.*nu+2.*nu2));
c2[1] =  4.78752;
c3[1] = 26.760136;
c4[1] = 43.861478;
    
c1[0] = 0.58121;
c2[0] = 1.01059;
c3[0] = 7.955729;
c4[0] = 1.650228;
    
// NOTE: the following is a polynomial evaluation, can be optimized
rhoHlm[1] = 1. + c1[1]*x + c2[1]*x2 + c3[1]*x3 + c4[1]*x4;
rhoHlm[0] = 1. + c1[0]*x + c2[0]*x2 + c3[0]*x3 + c4[0]*x4;

// Compute horizon multipolar flux (only l=2)
const double Heff2 = Heff*Heff;
const double jhat2 = jhat*jhat;

FlmH[k22] = FlmHLO[k22] * Heff2 * gsl_pow_int(rhoHlm[k22],4);
FlmH[k21] = FlmHLO[k21] * jhat2 * gsl_pow_int(rhoHlm[k21],4);

//Sum over multipoles and normalize to the 22 Newtonian multipole
double const hatFH = (FlmH[0]+FlmH[1])/FNewt22;

return hatFH;
}
#endif /* HorizonFlux_h */
