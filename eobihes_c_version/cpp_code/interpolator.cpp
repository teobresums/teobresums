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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#include <ios>
#include <cmath>
#include <math.h>
#include <vector>
#include <limits>
#include <stdio.h>
#include <fstream>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_complex_math.h>

#include "hlm.h"
#include "flux.h"
#include "Metric.h"
#include "interpolator.h"


using namespace::std;

double interpolate(double dt,vector<gsl_complex> grid)
{
    double xi, yi;
    double x[]      = {0.,0.,0.,0.,0.,0.,0.};
    double y[]      = {0.,0.,0.,0.,0.,0.,0.};
    double step     = 0.01; //dt/5.;
    double omeg_max = 0.;
    double t_max    = 0.;
    bool peak_flag  = false;
    
    for (int i=0; i<=6; i++)
    {
        x[i] = grid[i].dat[0]; //time
        y[i] = grid[i].dat[1]; //omega
    }
    

    
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline    = gsl_spline_alloc (gsl_interp_cspline, 7);
    gsl_spline_init (spline, x, y, 7);
    
    for (xi = x[0]; xi < x[6]; xi += step)
    {
        yi = gsl_spline_eval (spline, xi, acc);
        if (peak_flag==false)
        {
            if (yi<omeg_max)
            {
                peak_flag=true;
            }
            else
            {
                omeg_max=yi;
                t_max=xi;
            }
        }
    }
    
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
    
    return t_max;
}
