//
//  interpolator_waveform.h
//
//  Created by Philipp Fleig on 17/10/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef interpolator_waveform_h
#define interpolator_waveform_h

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
#include "Q_omega.h"

#include <limits>

using namespace::std;

int interpolate_wf(double dt,vector<double> t_vec,vector<double> hlm_rad, vector<double> hlm_phase,vector<double> Momg_vec,bool waveform_flag,std::ofstream& wave,std::ofstream& Qomg,double Mbh)
{

    //NOTE: BEFORE, THE Momg_vec HAD A Mbh multiplied onto it!!!!!!!!!

    int t_length = t_vec.size();
    double xi, yi;
    int grid_length = (int)(t_vec.back()-t_vec[0])/dt + 2;
    vector<gsl_complex> hlm_interp(grid_length);
    vector<double> omg_interp(grid_length);
    vector<double> t_interp(grid_length);
    
    // Convert all vectors to an array
    double* t = &t_vec[0];
    double* radial = &hlm_rad[0];
    double* phase = &hlm_phase[0];
    double Momg = Momg_vec[0];
        
    double step = dt;
    {
        gsl_interp_accel *acc = gsl_interp_accel_alloc ();
        gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, t_length);
        
        gsl_spline_init (spline, t, radial, t_length);
        int i=0;
        for (xi = t_vec[0]; xi < t_vec.back(); xi += step)
        {
            yi = gsl_spline_eval (spline, xi, acc);
            hlm_interp[i].dat[0] = yi;
            t_interp[i] = xi;
            i++;
        }
        
        gsl_spline_init (spline, t, phase, t_length);
        i=0;
        for (xi = t_vec[0]; xi < t_vec.back(); xi += step)
        {
            yi = gsl_spline_eval (spline, xi, acc);
            hlm_interp[i].dat[1] = yi;
            i++;
        }

        //(gsl spline * spline, const double xa[], const double ya[], size t size)
        gsl_spline_init (spline, t, &Momg, t_length);
        i=0;
        for (xi = t_vec[0]; xi < t_vec.back(); xi += step)
        {
            yi = gsl_spline_eval (spline, xi, acc);
            omg_interp[i] = yi;
            i++;
        }
        
        gsl_spline_free (spline);
        gsl_interp_accel_free (acc);
        
    }
    
        wave.precision(dbl::max_digits10);
        gsl_complex z;
        gsl_complex hlm_rect;
        for (int i=0; i < grid_length-1; i++) {
        //for (int i=0; i < hlm_rad.size(); i++) {
            //convert to rectangular complex form
            z = gsl_complex_polar(hlm_interp[i].dat[0],hlm_interp[i].dat[1]);
            //z = gsl_complex_polar(hlm_rad[i],hlm_phase[i]);
            hlm_rect.dat[0]=GSL_REAL(z);
            hlm_rect.dat[1]=GSL_IMAG(z);
            
            if (waveform_flag==true) {
                wave << t_interp[i] << "\t" << hlm_interp[i].dat[0] << "\t" << hlm_interp[i].dat[1] << "\t" << hlm_rect.dat[0] << "\t" << hlm_rect.dat[1] << endl;
            }
        }
    
    bool Qomega_flag = true;
    if (Qomega_flag==true) {
        Q_omega(t_interp,omg_interp,Qomg);
    }

    
    return 0;
}

#endif /* interpolator_waveform_h */
