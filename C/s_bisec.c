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
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "s_GS.h"
#include "input_struc.h"
#include "s_bisec.h"

struct energy_params
{
    double rorb, A, dA, rc, drc_dr, ak2, S, Ss, nu, chi1, chi2, X1, X2, c3;
};


double DHeff0(double x, void *DHeff_params)
{
    
    struct energy_params *p
    = (struct energy_params *) DHeff_params;

    double rorb   = p->rorb;
    double A      = p->A;
    double dA     = p->dA;
    double rc     = p->rc;
    double drc_dr = p->drc_dr;
    double ak2    = p->ak2;
    double S      = p->S;
    double Ss     = p->Ss;
    double nu     = p->nu;
    double chi1   = p->chi1;
    double chi2   = p->chi2;
    double X1     = p->X1;
    double X2     = p->X2;
    double c3     = p->c3;
    
    double ggm0[14];
    s_GS(ggm0, rorb, rc, drc_dr, ak2, 0., x, nu, chi1, chi2, X1, X2, c3);
    double dGS_dr  = ggm0[6];
    double dGSs_dr = ggm0[7];

    double x2 = x*x;
    
    double uc  = 1./rc;
    double uc2 = uc*uc;
    double uc3 = uc2*uc;
    
    // Orbital circular effective Hamiltonian
    double Horbeff0 = sqrt(A*(1. + x2*uc2));

    double dHeff_dr = x*(dGS_dr*S + dGSs_dr*Ss) + 1./(2.*Horbeff0)*( dA*(1. + x2*uc2) - 2.*A*uc3*drc_dr*x2);
    
    return dHeff_dr;
}

double s_bisec(
               double pph,
               double rorb,
               double A,
               double dA,
               double rc,
               double drc_dr,
               double ak2,
               double S,
               double Ss,
               void *params
){
    double nu   = (*(input *)params).nu;
    double chi1 = (*(input *)params).chi1;
    double chi2 = (*(input *)params).chi2;
    double X1   = (*(input *)params).X1;
    double X2   = (*(input *)params).X2;
    double c3   = (*(input *)params).cN3LO;

    int status;
    int iter = 0, max_iter = 200;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    
    double r;
    double x_lo = 0.5*pph, x_hi = 1.5*pph;
    gsl_function F;
    struct energy_params DHeff_params = {rorb,A,dA,rc,drc_dr,ak2,S,Ss,nu,chi1,chi2,X1,X2,c3};
    F.function = &DHeff0;
    
    F.params = &DHeff_params;
    T = gsl_root_fsolver_bisection;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        r      = gsl_root_fsolver_root (s);
        x_lo   = gsl_root_fsolver_x_lower (s);
        x_hi   = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi,
                                         0, 0.000000000000001);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (s);
    
    return r;
}

