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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <vector>
#include <cmath>
#include <ios>
#include <fstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <limits>

#include "TEOBResum.h"

typedef std::numeric_limits< double > dbl;

using namespace::std;

int rhs(double t, const double y[], double f[], void *params){
    (void)(t); /* avoid unused parameter warning */

    double nu = (*(TEOBResumParams *)params).nu;
    
    const double z3 = 2.0*nu*(4.0-3.0*nu);
    const double prstar  = y[2];
    const double prstar2 = prstar*prstar;
    const double prstar3 = prstar2*prstar;
    const double prstar4 = prstar3*prstar;
    
    const double pphi  = y[3];
    const double pphi2 = pphi*pphi;
    
    const double r  = y[0];
    const double r2 = r*r;
    const double u  = 1./r;
    const double u2 = u*u;
    const double u3 = u2*u;
    
    const vector<double> metric = Metric(r, params, false);
    const double A     = metric[0];
    const double dA    = metric[1];
    const double B     = metric[3];
    const double dB    = metric[4];
    const double one_A = 1./A;
    const double one_B = 1./B;

    //const double j = pphi;
    
    double Heff = sqrt(prstar2+A*(1. + pphi2*u2 +  z3*prstar4*u2));
    const double H = sqrt( 1. + 2.*nu*(Heff - 1.) )/nu;
    const double E = H*nu;
    
    const double sqrAB = sqrt(A/B);
    
    //r evol eqn rhs
    f[0]  = (prstar+2.0*z3*A*prstar3*u2)/Heff;
    f[0] *= (sqrAB / E);

    //phi evol eqn rhs
    double Omega   = A*pphi*u2/Heff;
    Omega *= 1./E;
    f[1] = Omega;
    
    //prstar evol eqn rhs
    f[2]  = (dA + ( pphi2 + z3*prstar4 )*(dA*u2-2.0*A*u3))/Heff;
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
    
    
    const double r_omega = r*pow(psi,1.0/3.0);
    const double v_phi   = r_omega*Omega;
    const double x       = v_phi * v_phi;
    const double jhat    = pphi/(r_omega*v_phi);
    
    const double E2       = E*E;
    const double tmpE     = 1./Heff+nu/E2;
    const double denE     = E*Heff;
    const double one_denE = 1./denE;
    
    const double dHeff_dr      = 0.5*(dA + (pphi2 + z3*prstar4)*(dA*u2 - 2.*A*u3))/Heff;
    const double dHeff_dprstar = (prstar + z3*2.0*A*u2*prstar3)/Heff;
    const double dprstar_dt    = f[2];
    const double dr_dt         = f[0];
    const double ddotr_dr      = sqrAB*( (prstar + z3*2.*A*u2*prstar3)*(0.5*(dA*one_A-dB*one_B)-dHeff_dr*tmpE)+ 2.0*z3*(dA*u2 - 2.*A*u3)*prstar3)/E;
    const double ddotr_dprstar = sqrAB*( 1.+z3*6.*A*u2*prstar2-(prstar + z3*2.*A*u2*prstar3)*dHeff_dprstar*tmpE)*one_denE;
    
    
    // Approximate ddot(r) without Flux <= ???
    const double ddotr = dprstar_dt*ddotr_dprstar + dr_dt*ddotr_dr;
    
    double source[] = {
        jhat,Heff,
        Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,
        Heff,jhat,Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,jhat,Heff,
        Heff,jhat,Heff,jhat,Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,jhat,Heff,jhat,Heff};
    
    f[3] = flux(x,Omega,r_omega,E,Heff,jhat,r, prstar,ddotr,source,params);

    
    return GSL_SUCCESS;
}

