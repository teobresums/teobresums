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

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

#include "Array.h"
#include "AdiabLR.h"
#include "HealyBBHFitRemnant.h"
#include "file_names.h"
#include "find_a1a2a3.h"
#include "initial.h"
#include "input_struc.h"
#include "interp_grid.h"
//#include "interpolate_wf.h"
#include "multipole_index.h"
#include "QNMHybridFitCab.h"
#include "read_config.h"
//#include "ringdown.h"
#include "RHS.h"
#include "s_initial.h"
#include "s_RHS.h"
#include "s_waveform.h"


int main (int argc, char* argv[])
{
    
    int lm, solver_scheme, i;
    double q, r0, dt, chi1, chi2, r_min, rLR, nu, r, prstar, phi, pphi, MOmg, t, y[4], t1, h, r_LSO, MOmg_prev, t_stop, Omg, Omg_orb, A, ddotr, ti;
    bool stop_flag, MOmgpeak_flag;
    
    
    q    = atof(argv[1]);
    chi1 = atof(argv[2]);
    chi2 = atof(argv[3]);
    r0   = atof(argv[4]);
    
    input params  = read_config(q,chi1,chi2,r0);

    // lm = params.lm;
    q             = params.q;
    nu            = params.nu;
    r0            = params.r0;
    dt            = params.dt;
    solver_scheme = params.solver_scheme;
    
    if (params.tidal==true) {
        rLR        = AdiabLR(&params);
        params.rLR = rLR;
        r_min      = rLR;
        printf("%s %.16e \n","rLR",params.rLR);
    }
    
    /* Creating folder with permission to read, write and execute*/
    mkdir("data",0777);
    
    /* Output file definitions */
    //vector<string> fnames = file_names(&params);
    
    char fnames_data[512];
    char fnames_init[512];
    char wavenames[35][512];
    file_names(fnames_data, fnames_init, wavenames, &params);

    FILE *data = fopen(fnames_data, "w");
    FILE *init = fopen(fnames_init, "w");

    /* Defining data vectors and variables */
    
    //std::vector<gsl_complex> hlm_vec={};
    
    Double_Array t_vec;
    Double_Array r_vec;

    Double_Array MOmg_vec;
    Double_Array Omg_orb_vec;
    Double_Array ddotr_vec;
    Double_Array prstar_vec;
    Double_Array pph_vec;

    Double_Array hlm_rad_vec;
    Double_Array hlm_phase_vec;
    Double_Array hlm_ampl[35];
    Double_Array hlm_phase[35];

    // Memory has to be free at the end.
    initArray(&t_vec, 5);
    initArray(&r_vec, 5);

    initArray(&MOmg_vec, 5);
    initArray(&Omg_orb_vec, 5);
    initArray(&ddotr_vec, 5);
    initArray(&prstar_vec, 5);
    initArray(&pph_vec, 5);

    initArray(&hlm_rad_vec, 5);
    initArray(&hlm_phase_vec, 5);
    for (int k=35; k--; ){
      initArray(&hlm_ampl[k], 5);
      initArray(&hlm_phase[k], 5);
    }

    /** Computing the initial conditions */
    double initial_data[7];
    // FIXME sys cannot be modified once it is initialized. Maybe we need to put rhs and s_RHS in one function.
//    gsl_odeiv2_system sys = {rhs, NULL, 4, &params};
    gsl_odeiv2_system sys = {s_RHS, NULL, 4, &params};
    
    if (params.spin==true) {
//        sys = {s_RHS, NULL, 4, &params};
        s_initial(initial_data, &params);
    } else if (params.spin==false) {
//        sys = {rhs, NULL, 4, &params};
        initial(initial_data, &params);
    }

    fprintf( init, "r p_phi p_r* p_r j E0 Omega_j\n" );
    fprintf( init, "%.17f %.17f %.17f %.17f %.17f %.17f %.17f\n", initial_data[0], initial_data[1], initial_data[2], initial_data[3], initial_data[4], initial_data[5], initial_data[6]);
    printf ("%.20e %.20e %.20e %.20e %.20e %.20e %.20e \n", initial_data[0], initial_data[1], initial_data[2], initial_data[3], initial_data[4], initial_data[5], initial_data[6]);
    
    /** Initial conditions: t, r, phi, prstar, pphi */
    t = 0.0;
    y[0] = initial_data[0];
    y[1] = 0.;
    y[2] = initial_data[2];
    y[3] = initial_data[1];

    double final_mass = HealyBBHFitRemnant(chi1, chi2, q);
    printf("%s %.8e \n","final BBH mass", final_mass);
    
    clock_t start, end;
    start = clock();
    
    /** Initialize ODE system solver */
    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_step * s            = gsl_odeiv2_step_alloc (T, 4);
    gsl_odeiv2_control * c         = gsl_odeiv2_control_y_new (1.e-13, 1.e-11);
    gsl_odeiv2_evolve * e          = gsl_odeiv2_evolve_alloc (4);
    gsl_odeiv2_driver * d          = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,1e-2, 1000., 1000.);
    
    t     = 0.0;
    r_LSO = 6.0;
    t1    = 1.e15;
    h     = 0.0001;
    
    MOmg_prev     = 0.;
    t_stop        = 0.;
    Omg           = 0.;
    Omg_orb       = 0.;
    A             = 0.;
    ddotr         = 0.;
    stop_flag     = false;
    MOmgpeak_flag = false;
    
    gsl_complex h_form[35];
    while (stop_flag == false)
    {
        //while (y[0]>r_min) {
        
        switch (solver_scheme)
        {
            case 0:
            {
                int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, y);
                //time_step << t << "\t" << h << "\t" << y[3] << endl;
                if (status != GSL_SUCCESS)
                {
                    break;
                }
                break;
            }
            case 1:
            {
                if (y[0]>r_LSO)
                {
                    int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, y);
                    if (status != GSL_SUCCESS)
                    {
                        break;
                    }
                }
                else
                {
                    double ti = t+dt;
                    int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
                    if (status != GSL_SUCCESS)
                    {
                        printf ("error, return value=%d\n", status);
                        break;
                    }
                }
                break;
            }
            case 2:
            {
                double ti = t+dt;
                int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
                if (status != GSL_SUCCESS)
                {
                    printf ("error, return value=%d\n", status);
                    break;
                }
                break;
            }
            default:
            {
                int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, y);
                if (status != GSL_SUCCESS)
                {
                    break;
                }
                break;
            }
        }
        
        int status  = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, y);
        if (status != GSL_SUCCESS) {
            break;
        }
        
        /** Read out computation */
        r      = y[0];
        phi    = y[1];
        prstar = y[2];
        pphi   = y[3];
        
        /** Checking whether the dynamics produces NaN values; this can happen if radius r becomes too small */
        if (r!=r) {
            printf("%s \n","dynamics is producing NaN values");
        } else {
            
            /** Write and print dynamics*/
            if (params.dynamics==true) {
              fprintf( data, "%.17f %.17f %.17f %.17f %.17f\n", t, r, phi, prstar, pphi );
            }
            
            /** Waveform computation*/
//            printf( "%.17f %.17f %.17f %.17f %.17f\n", t, r, phi, prstar, pphi );
            s_waveform(h_form, &Omg, &Omg_orb, &A, &ddotr, t, y, &params);
            
            /** Append dynamics and waveform to vectors */
            // FIXME: lm doesnt have any value, get it from input params? (KaWa)
            insertArray( &hlm_rad_vec, h_form[lm].dat[0]);
            insertArray( &hlm_phase_vec, h_form[lm].dat[1]);
            
            for (int k=35; k--; ) {
              insertArray( &hlm_ampl[k], h_form[k].dat[0] );
              insertArray( &hlm_phase[k], h_form[k].dat[1] );
            }
        
            insertArray( &t_vec, t);
            insertArray( &MOmg_vec, Omg);
            insertArray( &r_vec, r);
            insertArray( &pph_vec, pphi);
            insertArray( &prstar_vec, prstar);
            insertArray( &Omg_orb_vec, Omg_orb);
            insertArray( &ddotr_vec, ddotr);
        }
        
        /** Check when to break the computation;find peak of omega curve and continue for delta_t=10. afterwards */
        //MOmg = Omg; //NOTE: was MOmg = Omg_orb; before!!! (only for the spinning case)
        if (params.spin==true) {
            MOmg = Omg_orb;
        } else {
            MOmg = Omg;
        }

        if (MOmgpeak_flag==false) {
          if (MOmg < MOmg_prev) {
              MOmgpeak_flag = true;
              t_stop        = t + 10.;
          } else {
              MOmg_prev = MOmg;
          }
        } else {
          if (t >= t_stop) {
            stop_flag = true;
          }
        }
    }
    gsl_odeiv2_evolve_free (e); gsl_odeiv2_control_free (c); gsl_odeiv2_step_free (s);gsl_odeiv2_driver_free (d);
    
    /** Interpolate quantities on a grid of width dt */
    const int grid_length = (int)(t_vec.array[t_vec.used-1]-t_vec.array[0])/dt + 1;
    
    double t_vecg[grid_length];
    i  = 0;
    ti = 0.;
                        
    for (ti = t_vec.array[0]; ti < t_vec.array[t_vec.used-1]; ti += dt) {
      t_vecg[i] = ti;
      i++;
    }

    /* Allocate needed arrays*/
    // FIXME why t_length != grid_length? In find_a1a2a3, They are assummed to be the same.
    const int t_length = (int)(t_vec.array[t_vec.used-1]-t_vec.array[0])/dt + 2;
    double          r_vecg[t_length];
    double       MOmg_vecg[t_length];
    double        pph_vecg[t_length];
    double     prstar_vecg[t_length];
    double  hlm_phase_vecg[t_length];
    double    hlm_rad_vecg[t_length];
    double      ddotr_vecg[t_length];
    double     OmgOrb_vecg[t_length];
    double  hlm_ampl_g[35][t_length];
    double hlm_phase_g[35][t_length];
    
    /* Fill in the arrays passing them as input */
    interp_grid(r_vecg,         t_vec, r_vec.array,         dt);
    interp_grid(MOmg_vecg,      t_vec, MOmg_vec.array,      dt);
    interp_grid(pph_vecg,       t_vec, pph_vec.array,       dt);
    interp_grid(prstar_vecg,    t_vec, prstar_vec.array,    dt);
    interp_grid(hlm_phase_vecg, t_vec, hlm_phase_vec.array, dt);
    interp_grid(hlm_rad_vecg,   t_vec, hlm_rad_vec.array,   dt);
    interp_grid(ddotr_vecg,     t_vec, ddotr_vec.array,     dt);
    interp_grid(OmgOrb_vecg,    t_vec, Omg_orb_vec.array,   dt);
    
 
    for (int k=35; k--; )
    {
        interp_grid( hlm_ampl_g[k], t_vec, hlm_ampl[k].array,  dt);
        interp_grid(hlm_phase_g[k], t_vec, hlm_phase[k].array, dt);
    }

    /** NQCs */
    if (params.tidal==false && params.spin==true)
    { 
        gsl_complex nqc[35][grid_length];
        find_a1a2a3(grid_length, nqc, t_vecg, r_vecg, MOmg_vecg, pph_vecg, prstar_vecg, hlm_phase_g, OmgOrb_vecg, hlm_ampl_g, ddotr_vecg, &params);
      
        for (int k=35; k--; )
        {
            for (int i=grid_length; i--; )
            {
                hlm_ampl_g[k][i]  = hlm_ampl_g[k][i]  * nqc[k][i].dat[0];
                hlm_phase_g[k][i] = hlm_phase_g[k][i] + nqc[k][i].dat[1];
            }
        }
    }
    
    /** Define a time vector for each multipole */
    double t_g[35][grid_length];
    for (int k=35; k--; )
    {
        for (int j=0; j<grid_length;j++)
        {
            t_g[k][j] = t_vecg[j];
        }
    }
    
    /** Ringdown attachment */
    if (params.tidal==false && params.spin==false)
    {
//        ringdown(nu,q,dt,final_mass,t_g,MOmg_vecg,hlm_ampl_g,hlm_phase_g);
    }
    
    /** Compute interpolation of waveform on grid and write to output file */
//    interpolate_wf(dt, t_g, hlm_ampl_g,hlm_phase_g,params.waveform,wavenames,final_mass);
    
    end = clock();
 
    printf("Time required for execution: %e seconds.\n\n", (double)(end-start)/CLOCKS_PER_SEC);
    
    return 0;
}
