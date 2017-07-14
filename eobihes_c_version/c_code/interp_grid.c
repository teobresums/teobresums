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

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "interp_grid.h"

void interp_grid(
    double data_g[],             /** OUTPUT Dimension: 35*/
    double t_vec[],
    const int size_t_vec,
    double data[],
    double dt){
    
    int i = 0;
    const int t_length = size_t_vec;
    const int grid_length = (int)(t_vec[size_t_vec-1]-t_vec[0])/dt + 2;

    double xi, yi;
    double omg_interp[grid_length];
    double t_interp[grid_length];

    double step = dt;

    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, t_length);
    gsl_spline_init (spline, t_vec, data, t_length);

    for (xi = t_vec[0]; xi < t_vec[size_t_vec-1]; xi += step)
    {
        yi          = gsl_spline_eval (spline, xi, acc);
        data_g[i]   = yi;
        t_interp[i] = xi;
        i++;
    }

    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

}
