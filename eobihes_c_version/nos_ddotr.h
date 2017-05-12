//
//  ddotr.h
//  EOB_IHES_cpp_bitbucket
//
//  Created by Philipp Fleig on 10/04/2017.
//  Copyright © 2017 Philipp Fleig. All rights reserved.
//

#ifndef ddotr_h
#define ddotr_h

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include "Metric.h"
#include <vector>
#include "cmath"
#include "Flux.h"

#include "hlm.h"

#include <ios>
#include <fstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "input_struc.h"
#include "multipole_index.h"


#include <limits>

typedef std::numeric_limits< double > dbl;

using namespace::std;

double nos_ddotr(double t,double r,double pphi,double prstar,void *params) {

    
    vector<double> f(4);
    double nu = (*(input *)params).nu;
    
    const double z3=2.0*nu*(4.0-3.0*nu);

    const double prstar2=prstar*prstar;
    const double prstar3=prstar2*prstar;
    const double prstar4=prstar3*prstar;

    const double pphi2=pphi*pphi;

    const double r2 = r*r;
    const double u=1./r;
    const double u2=u*u;
    const double u3=u2*u;
    
    const vector<double> metric=Metric(r, params,false);
    const double A=metric[0];
    const double dA=metric[1];
    const double B=metric[3];
    const double dB = metric[4];
    const double one_A    = 1./A;
    const double one_B    = 1./B;
    
    //const double j = pphi;
    
    double Heff=sqrt(prstar2+A*(1. + pphi2*u2 +  z3*prstar4*u2));
    const double H = sqrt( 1. + 2.*nu*(Heff - 1.) )/nu;
    const double E = H*nu;
    
    const double sqrAB = sqrt(A/B);
    
    //r evol eqn rhs
    f[0] = (prstar+2.0*z3*A*prstar3*u2)/Heff;
    f[0] *= (sqrAB / E);
    
    //phi evol eqn rhs
    double Omega   = A*pphi*u2/Heff;
    Omega *= 1./E;
    f[1] = Omega;
    
    //prstar evol eqn rhs
    f[2] = (dA + ( pphi2 + z3*prstar4 )*(dA*u2-2.0*A*u3))/Heff;
    f[2] *= - 0.5 * sqrAB / E;
    
    //pphi evol eqn rhs
    const double sqrW = sqrt(A*(1. + pphi2*u2));
    double psi = 0.;
    bool psi_flag = false; // flag for an alternative way of computing psi
    if (psi_flag==false) {
        psi = 2.*(1.0 + 2.0*nu*(sqrW - 1.0))/(r2*dA);
    } if (psi_flag==true) {
        psi = 2.*(1.0 + 2.0*nu*(Heff - 1.0))/(r2*dA);
    }
    
    const double E2       = E*E;
    const double tmpE     = 1./Heff+nu/E2;
    const double denE     = E*Heff;
    const double one_denE = 1./denE;
    
    const double dHeff_dr      = 0.5*(dA + (pphi2 + z3*prstar4)*(dA*u2 - 2*A*u3))/Heff;
    const double dHeff_dprstar = (prstar + z3*2.0*A*u2*prstar3)/Heff;
    const double dprstar_dt = f[2];
    const double dr_dt      = f[0];
    const double ddotr_dr = sqrAB*( (prstar + z3*2.*A*u2*prstar3)*(0.5*(dA*one_A-dB*one_B)-dHeff_dr*tmpE)+ 2.0*z3*(dA*u2 - 2*A*u3)*prstar3)/E;
    const double ddotr_dprstar = sqrAB*( 1+z3*6.*A*u2*prstar2-(prstar + z3*2*A*u2*prstar3)*dHeff_dprstar*tmpE)*one_denE;
    
    double ddotr = dprstar_dt*ddotr_dprstar + dr_dt*ddotr_dr;

    return ddotr;
}

#endif /* ddotr_h */
