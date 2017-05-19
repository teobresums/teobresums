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
#include <fstream>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <math.h>
#include <cmath>
#include <limits>
#include <vector>

#include "Metric.h"
#include "flux.h"
#include "Q_omega.h"
#include "hlm.h"

using namespace::std;

int interpolate_wf(double dt, vector<double> t_vec, vector<double> hlm_rad, vector<double> hlm_phase, vector<double> Momg_vec, bool waveform_flag, std::ofstream& wave, std::ofstream& Qomg, double Mbh)
{

    /** NOTE: BEFORE, THE Momg_vec HAD A Mbh multiplied onto it! */

    double xi, yi;
    int i=0;
    int t_length     = t_vec.size();
    int grid_length  = (int)(t_vec.back()-t_vec[0])/dt + 2;
    double* t        = &t_vec[0];
    double* radial   = &hlm_rad[0];
    double* phase    = &hlm_phase[0];
    double Momg      = Momg_vec[0];
    double step      = dt;
    bool Qomega_flag = true;
    
    vector<gsl_complex> hlm_interp(grid_length);
    vector<double> omg_interp(grid_length);
    vector<double> t_interp(grid_length);
    
    /** Convert all vectors to an array */
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline    = gsl_spline_alloc (gsl_interp_cspline, t_length);
    gsl_spline_init (spline, t, radial, t_length);
    
    for (xi = t_vec[0]; xi < t_vec.back(); xi += step)
    {
        yi                   = gsl_spline_eval (spline, xi, acc);
        hlm_interp[i].dat[0] = yi;
        t_interp[i]          = xi;
        i++;
    }
    
    gsl_spline_init (spline, t, phase, t_length);
    i=0;
    
    for (xi = t_vec[0]; xi < t_vec.back(); xi += step)
    {
        yi                   = gsl_spline_eval (spline, xi, acc);
        hlm_interp[i].dat[1] = yi;
        i++;
    }

    gsl_spline_init (spline, t, &Momg, t_length);
    i=0;
    
    for (xi = t_vec[0]; xi < t_vec.back(); xi += step)
    {
        yi            = gsl_spline_eval (spline, xi, acc);
        omg_interp[i] = yi;
        i++;
    }
    
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
    
    wave.precision(dbl::max_digits10);
    gsl_complex z;
    gsl_complex hlm_rect;
    
    for (int i=0; i < grid_length-1; i++)
    {
        /** Convert to rectangular complex form*/
        z               = gsl_complex_polar(hlm_interp[i].dat[0],hlm_interp[i].dat[1]);
        hlm_rect.dat[0] = GSL_REAL(z);
        hlm_rect.dat[1] = GSL_IMAG(z);
        
        if (waveform_flag==true)
        {
            wave << t_interp[i] << "\t" << hlm_interp[i].dat[0] << "\t" << hlm_interp[i].dat[1] << "\t" << hlm_rect.dat[0] << "\t" << hlm_rect.dat[1] << endl;
        }
    }
    
    if (Qomega_flag==true)
    {
        Q_omega(t_interp,omg_interp,Qomg);
    }
    
    return 0;
}