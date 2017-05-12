//
//  interp.h
//  EOB_IHES_cpp_bitbucket
//
//  Created by Philipp Fleig on 27/04/2017.
//  Copyright © 2017 Philipp Fleig. All rights reserved.
//

#ifndef interp_h
#define interp_h

#include <gsl/gsl_errno.h>
#include <vector>

using namespace::std;

vector<double> interp_grid(vector<double> t_vec,vector<double> data,double dt)
{
    
    //NOTE: BEFORE, THE Momg_vec HAD A Mbh multiplied onto it!!!!!!!!!
    
    int t_length = t_vec.size();
    double xi, yi;
    int grid_length = (int)(t_vec.back()-t_vec[0])/dt + 2;
    vector<double> data_g(grid_length);
    vector<double> omg_interp(grid_length);
    vector<double> t_interp(grid_length);
    
    // Convert all vectors to an array
    double* t = &t_vec[0];
    double* data_arr = &data[0];
    
    double step = dt;
    {
        gsl_interp_accel *acc = gsl_interp_accel_alloc ();
        gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, t_length);
        
        gsl_spline_init (spline, t, data_arr, t_length);
        int i=0;
        for (xi = t_vec[0]; xi < t_vec.back(); xi += step)
        {
            yi = gsl_spline_eval (spline, xi, acc);
            data_g[i] = yi;
            t_interp[i] = xi;
            i++;
        }
        
        gsl_spline_free (spline);
        gsl_interp_accel_free (acc);
        
    }
    
    return data_g;
}

#endif /* interp_h */
