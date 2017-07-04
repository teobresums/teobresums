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

#include <math.h>
#include <vector>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "Metric.h"
#include "input_struc.h"
#include "AdiabLR.h"

using namespace::std;

double fLR(double r, void *params)
{
    
    vector<double> metric=Metric(r, params,true);
    
    double A  = metric[0];
    double dA = metric[2]; /**derivative w.r.t. u*/
    double f  = A +(0.5*dA)/r;
    
    return f;
}
/** Takes nu as input */
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