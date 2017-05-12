//
//  interpolator.h
//
//  Created by Philipp Fleig on 11/09/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef interpolator_h
#define interpolator_h

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include "Metric.h"
#include <vector>
#include "cmath"
#include "flux.h"

#include "hlm.h"

#include <ios>
#include <fstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "hlm.h"

#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


#include <limits>

using namespace::std;

double interpolate(double dt,vector<gsl_complex> grid)
{
    double xi, yi;
    double x[]={0.,0.,0.,0.,0.,0.,0.};
    double y[]={0.,0.,0.,0.,0.,0.,0.};
    
    for (int i=0; i<=6; i++) {
        x[i] = grid[i].dat[0]; //time
        y[i] = grid[i].dat[1]; //omega
    }
    
    double step = 0.01;//dt/5.;
    double omeg_max = 0.;
    double t_max = 0.;
    bool peak_flag=false;
    //printf("%.12e %.12e %.12e %.12e %.12e %.12e %.12e \n",x[0],x[1],x[2],x[3],x[4],x[5],x[6]);
    
        {
            gsl_interp_accel *acc = gsl_interp_accel_alloc ();
            gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, 7);
            gsl_spline_init (spline, x, y, 7);
            for (xi = x[0]; xi < x[6]; xi += step)
            {
                yi = gsl_spline_eval (spline, xi, acc);
                if (peak_flag==false) {
                    if (yi<omeg_max) {
                        peak_flag=true;
                    }else{
                        omeg_max=yi;
                        t_max=xi;
                    }
                }
            }
            gsl_spline_free (spline);
            gsl_interp_accel_free (acc);
        }
    
    return t_max;
}

#endif /* interpolator_h */
