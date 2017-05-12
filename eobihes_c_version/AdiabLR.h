//
//  AdiabLR.h
//
//  Created by Philipp Fleig on 28/10/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef EOBAdiabLR_h
#define EOBAdiabLR_h

#include <math.h>
#include <vector>
#include "Metric.h"
#include "input_struc.h"

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>


using namespace::std;

double fLR(double r, void *params){
    
    vector<double> metric=Metric(r, params,true);
    double A = metric[0];
    double dA = metric[2]; //derivative w.r.t. u

    double f = A + 0.5 * dA / r;
    
    return f;
}

double AdiabLR(void *params){//double nu){
    
    int status;
    int iter = 0, max_iter = 200;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    
    double rLR;
    double x_lo = 0.1, x_hi = 15.;
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
        rLR = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi,
                                         0, 0.000000000000001);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (s);
    
    return rLR;
}


#endif /* AdiabLR_h */
