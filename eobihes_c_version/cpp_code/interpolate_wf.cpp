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
#include "cmath"
#include <math.h>
#include <limits>
#include <vector>
#include <fstream>
#include <stdio.h>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "hlm.h"
#include "flux.h"
#include "Metric.h"
#include "Q_omega.h"
#include "interpolate_wf.h"

typedef std::numeric_limits< double > dbl;

using namespace::std;

//int interpolate_wf(double dt,vector<double> t_vec,vector<double> hlm_rad, vector<double> hlm_phase,bool waveform_flag,std::ofstream& wave,double Mbh)
int interpolate_wf(double dt,vector<vector<double> > t_vec,vector<vector<double> > hlm_rad, vector<vector<double> > hlm_phase,bool waveform_flag, double distance)
{

    /** Note: before the Momg_vec had a Mbh multiplied onto it! */

    for (int k=35; k--; )
    {
        int t_length = t_vec[k].size();
        int grid_length = (int)(t_vec[k].back()-t_vec[k][0])/dt + 2;
        vector<double> t_interp(grid_length);
        double xi, yi;
        int i=0;

        /** Convert all vectors to an array */
        if (k==1)
        {
            double* t      = &t_vec[k][0];
            double* radial = &hlm_rad[k][0];
            double* phase  = &hlm_phase[k][0];
            double step    = dt;

            gsl_interp_accel *acc = gsl_interp_accel_alloc ();
            gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, t_length);
            gsl_spline_init (spline, t, radial, t_length);

            for (xi = t_vec[k][0]; xi < t_vec[k].back(); xi += step)
            {
                yi = gsl_spline_eval (spline, xi, acc);
                //hlm_interp[i].dat[0] = yi;
                t_interp[i] = xi;
                i++;
            }
            gsl_spline_init (spline, t, phase, t_length);
            i=0;
            for (xi = t_vec[k][0]; xi < t_vec[k].back(); xi += step)
            {
                yi = gsl_spline_eval (spline, xi, acc);
                //hlm_interp[i].dat[1] = yi;
                i++;
            }
    
            gsl_spline_free (spline);
            gsl_interp_accel_free (acc);
//            for (int j=0; j < grid_length-2; j++)
//            {
//                waveform[j] += hlm_interp[j].dat[0]*cos(hlm_interp[j].dat[1])/distance;
//            }
        }
        else continue;

    }

    return 0;
}
