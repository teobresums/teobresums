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

#include <iostream>
#include <list>
#include <fstream>
#include <stdio.h>
#include <time.h>
#include <cmath>
#include <string>
#include <cstring>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <tuple>
#include <vector>
#include <complex.h>
#include <sys/stat.h>

#include "TEOBResumS.h"

using namespace::std;

static void swap_variables(double *v1, double *v2)
{
    double tmp;
    tmp = *v1;
    *v1 = *v2;
    *v2 = tmp;
}

void TEOBResumS(Waveform **hplus,       /** h+ return array                                      **/
                Waveform **hcross,      /** hx return array                                      **/
                double m1,              /** m1 (Msun)                                            **/
                double m2,              /** m2 (Msun)                                            **/
                double spin1x,          /** dimensionless s1x                                    **/
                double spin1y,          /** dimensionless s1y                                    **/
                double spin1z,          /** dimensionless s1z                                    **/
                double spin2x,          /** dimensionless s2x                                    **/
                double spin2y,          /** dimensionless s2y                                    **/
                double spin2z,          /** dimensionless s2z                                    **/
                double inclination,     /** inclination angle (rad)                              **/
                double polarisation,    /** polarisation angle (rad)                             **/
                double f_min,           /** starting frequency (Hz)                              **/
                double dt,              /** sampling interval (s)                                **/
                double LambdaAl2,       /** l=2 (tidal deformation of body 1)/(mass of body 1)^5 **/
                double LambdaBl2,       /** l=2 (tidal deformation of body 2)/(mass of body 2)^5 **/
                double LambdaAl3,       /** l=3 (tidal deformation of body 1)/(mass of body 1)^5 **/
                double LambdaBl3,       /** l=3 (tidal deformation of body 2)/(mass of body 2)^5 **/
                double LambdaAl4,       /** l=4 (tidal deformation of body 1)/(mass of body 1)^5 **/
                double LambdaBl4,       /** l=4 (tidal deformation of body 2)/(mass of body 2)^5 **/
                double distance,        /** distance(Mpc)                                        **/
                int    lm,              /** multipole index for output                           **/
                TEOBResumFlags *flags   /** flags **/
){
    
    int i = 0;
    int grid_length = 0;
    double rLR, r, prstar, phi, pphi, MOmg, t, y[4], t1, h, r_LSO, MOmg_prev, t_stop, Omg, Omg_orb, A, ddotr, ti;
    bool stop_flag, MOmgpeak_flag;
    
    if ((spin1x!=0)||(spin1y!=0)||(spin2x!=0)||(spin2y!=0))
    {
        printf("ERROR! Non-aligned spins not supported! Aborting.\n");
        exit(-1);
    }

    if (m2 > m1)
    {
        printf("Warning! m1 > m2, swapping component masses, spins and tidal coefficients\n");
        swap_variables(&m1, &m2);
        swap_variables(&spin1x, &spin2x);
        swap_variables(&spin1y, &spin2y);
        swap_variables(&spin1z, &spin2z);
        swap_variables(&LambdaAl2, &LambdaBl2);
        swap_variables(&LambdaAl3, &LambdaBl3);
        swap_variables(&LambdaAl4, &LambdaBl4);

    }
    
    TEOBResumParams params = process_input_parameters(
                                                    m1,
                                                    m2,
                                                    spin1z,
                                                    spin2z,
                                                    f_min,
                                                    dt,
                                                    LambdaAl2,
                                                    LambdaBl2,
                                                    LambdaAl3,
                                                    LambdaBl3,
                                                    LambdaAl4,
                                                    LambdaBl4,
                                                    flags);

    double q      = params.q;
    dt            = params.dt;

    if (params.flags.tidal==1)
    {
        rLR        = AdiabLR(&params);
        params.rLR = rLR;
    }
    
    /** Defining data vectors and variables */

    std::vector<double> t_vec={};
    std::vector<double> r_vec={};
    std::vector<double> pph_vec={};
    std::vector<double> MOmg_vec={};
    std::vector<double> ddotr_vec={};
    std::vector<double> prstar_vec={};
    std::vector<double> hlm_rad_vec={};
    std::vector<double> Omg_orb_vec={};
    std::vector<double> hlm_phase_vec={};
    std::vector<vector<double> > hlm_ampl(35);
    std::vector<vector<double> > hlm_phase(35);
    
    /** Computing the initial conditions */
    vector<double> initial_data(7);
    gsl_odeiv2_system sys = {rhs, NULL , 4, &params};
    
    if (params.flags.spin==1)
    {
        sys          = {s_RHS, NULL , 4, &params};
        initial_data = s_initial(&params);
    }
    else
    {
        sys          = {rhs, NULL , 4, &params};
        initial_data = initial(&params);
    }
    
    /** Initial conditions: t, r, phi, prstar, pphi */
    t    = 0.0;
    y[0] = initial_data[0];
    y[1] = 0.;
    y[2] = initial_data[2];
    y[3] = initial_data[1];


    double final_mass = HealyBBHFitRemnant(spin1z, spin2z, q);
    params.Mbh = final_mass;
    
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
    
    while (stop_flag == false)
    {
        switch (params.flags.solver_scheme)
        {
            case 0:
            {
                int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, y);

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
        
        /** Read out computation */
        r      = y[0];
        phi    = y[1];
        prstar = y[2];
        pphi   = y[3];
        
        /** Checking whether the dynamics produces NaN values
            this can happen if radius r becomes too small */
        if (std::isfinite(r))
        {
            /** Waveform computation*/
            vector<gsl_complex> h_form = s_waveform(t, y, &params, Omg, Omg_orb, A, ddotr);

            /** Append dynamics and waveform to vectors */
            hlm_rad_vec.push_back(h_form[lm].dat[0]);
            hlm_phase_vec.push_back(h_form[lm].dat[1]);
            
            for (int k=35; k--; )
            {
                hlm_ampl[k].push_back(h_form[k].dat[0]);
                hlm_phase[k].push_back(h_form[k].dat[1]);
            }
        
            t_vec.push_back(t);
            MOmg_vec.push_back(Omg);
            r_vec.push_back(r);
            pph_vec.push_back(pphi);
            prstar_vec.push_back(prstar);
            Omg_orb_vec.push_back(Omg_orb);
            ddotr_vec.push_back(ddotr);
        }
        
        /** Check when to break the computation;find peak of omega curve and continue for delta_t=10. afterwards */
        //MOmg = Omg; //NOTE: was MOmg = Omg_orb; before!!! (only for the spinning case)
        if (params.flags.spin==1)
        {
            MOmg = Omg_orb;
        }
        else
        {
            MOmg = Omg;
            
        }
        if (MOmgpeak_flag==false)
        {
            if (MOmg < MOmg_prev)
            {
                MOmgpeak_flag = true;
                t_stop        = t + 4.*dt;
            }
            else
            {
                MOmg_prev = MOmg;
            }
        }
        else
        {
            if (t >= t_stop)
            {
                stop_flag = true;
            }
        }
    }
    
    gsl_odeiv2_evolve_free (e);
    gsl_odeiv2_control_free (c);
    gsl_odeiv2_step_free (s);
    gsl_odeiv2_driver_free (d);
    

    /** To compute the NQC corrections, a precise determination of the time at which Omega_orb peaks is required,
        in order to construct a grid which passes from that t_peak.
        The first step is to find the t_peak on the grid. */

    double t_max_grid, omg_max_grid;
    int index_max;

    for (int index=1; index<Omg_orb_vec.size(); index++)
    {
        if(Omg_orb_vec[index] < Omg_orb_vec[index-1])
        {
            index_max    = index-1;
            omg_max_grid = Omg_orb_vec[index_max];
            t_max_grid   = t_vec[index_max];
            index        = Omg_orb_vec.size();
        }

    }

    cout << "t_max_grid:\t" <<  t_max_grid << "\nomg_max_grid:\t" <<  omg_max_grid << endl;

    /* Then take a few points around the peak and analytically interpolate these points
        to find a better approximation for t_peak. */
    double x1 = t_vec[index_max - 1];
    double x2 = t_vec[index_max    ];
    double x3 = t_vec[index_max + 1];
    double y1 = Omg_orb_vec[index_max - 1];
    double y2 = Omg_orb_vec[index_max    ];
    double y3 = Omg_orb_vec[index_max + 1];

    double c1 = (pow(x3,2)*(y1 - y2) + pow(x1,2)*(y2 - y3) +
                 pow(x2,2)*(-y1 + y3))/((x1 - x2)*(x1 - x3)*(x2 - x3));
    double c2 = (x3*(-y1 + y2) + x2*(y1 - y3) + x1*(-y2 + y3))/
    ((x1 - x2)*(x1 - x3)*(x2 - x3));

    double t_max = (-c1)/(2.*c2);

    /** This is just a temporary check to compare the omega computed on the grid
        and the one coming from the interpolation.*/
    const int omg_size = Omg_orb_vec.size();
    double Omg_array[omg_size];
    const int t_size = t_vec.size();
    double t_array[t_size];

    for (int k=0; k < omg_size; k++)
    {
        Omg_array[k] = Omg_orb_vec[k];
    }

    for (int k=0; k <t_size; k++)
    {
        t_array[k] = t_vec[k];
    }
    
    double omega_max = interp1d (3, t_max, omg_size, Omg_array, t_array);
    cout << "t_max_interp:\t" << t_max << "\nomg_max_interp:\t" <<  omega_max << endl;

    int N_before = int((t_max - t_vec[0])/dt);

    /** Interpolate quantities on a grid of width dt */
    grid_length = (int)((t_vec.back()-t_vec[0])/dt + 1);

    vector<double> t_vecg(grid_length);
    i  = 0;
    ti = 0;
    
    for (int k = N_before; k >= 0 ; k--)
    {
        t_vecg[k] = t_max - (N_before-k)*dt;
    }

    for (int k = N_before+1; k < grid_length; k++)
    {
        t_vecg[k] = t_max + (k-N_before)*dt;
    }

    for(int k = 0; k < grid_length; k++)
    {
        cout << "i:\t" << k << "\tt_vec[i]:\t" << t_vecg[k] << endl;
    }

    
    vector<double> hlm_phase_vecg = interp_grid(t_vec,hlm_phase_vec,dt);
    vector<double> hlm_rad_vecg   = interp_grid(t_vec,hlm_rad_vec,dt);
    vector<double> r_vecg         = r_vec;
    vector<double> MOmg_vecg      = MOmg_vec;
    vector<double> pph_vecg       = pph_vec;
    vector<double> prstar_vecg    = prstar_vec;
    vector<double> ddotr_vecg     = ddotr_vec;
    vector<double> OmgOrb_vecg    = Omg_orb_vec;
    if (params.flags.tidal == 0)
    {
        r_vecg         = interp_grid(t_vec,r_vec,dt);
        MOmg_vecg      = interp_grid(t_vec,MOmg_vec,dt);
        pph_vecg       = interp_grid(t_vec,pph_vec,dt);
        prstar_vecg    = interp_grid(t_vec,prstar_vec,dt);
        ddotr_vecg     = interp_grid(t_vec,ddotr_vec,dt);
        OmgOrb_vecg    = interp_grid(t_vec,Omg_orb_vec,dt);
    }
    std::vector<vector<double> > hlm_ampl_g(35);
    std::vector<vector<double> > hlm_phase_g(35);
    
    if (DEBUG)
    {
        char   outputr[256]   = "waveform_preIntepolation.dat";
        std::FILE* waveform_preint   = std::fopen(outputr, "w");
        int j                 = 0;
        int N                 = hlm_ampl[1].size();
        
        for (j=0;j<N;j++)
        {
            std::fprintf(waveform_preint, "%f\t%e\t%e\n", t_vec[j], hlm_ampl[1][j], hlm_phase[1][j]);
        }
        std::fclose(waveform_preint);
    }
    
    
    for (int k=35; k--; )
    {
        vector<double> amplitude = hlm_ampl[k];
        vector<double> phase     = hlm_phase[k];
        hlm_ampl_g[k]            = interp_grid(t_vec,amplitude,dt);
        hlm_phase_g[k]           = interp_grid(t_vec,phase,dt);
    }
    
    /** NQCs corrections */
    /** NOTE THAT IF YOU REMOVE PARAMS.SPIN==TRUE EVERYTHING IS FUCKED UP FOR SOME REASON */
    if (DEBUG)
    {
        char   outputr[256]   = "waveform_preNQC.dat";
        std::FILE* waveform_preNQC   = std::fopen(outputr, "w");
        int j                 = 0;
        int N                 = hlm_ampl_g[1].size();
        
        for (j=0;j<N;j++)
        {
            std::fprintf(waveform_preNQC, "%f\t%e\t%e\n", (double)j*dt, hlm_ampl_g[1][j], hlm_phase_g[1][j]);
        }
        std::fclose(waveform_preNQC);
    }
    
//    if (DEBUG)
//    {
        char   outputr[256]   = "waveform_nqc.dat";
        std::FILE* nqcs       = std::fopen(outputr, "w");
//    }
    
    if (params.flags.tidal==0 && params.flags.spin==1)
    {
        
        vector<vector<gsl_complex> > nqc = find_a1a2a3(t_vecg,r_vecg,MOmg_vecg,pph_vecg,prstar_vecg,hlm_phase_g,OmgOrb_vecg,hlm_ampl_g,ddotr_vecg,&params);
        
        for (int k=35; k--; )
        {
            if (k==1)
            {
                for (int i=0; i<grid_length; i++ )
                {
                    hlm_ampl_g[k][i]  = hlm_ampl_g[k][i]  * nqc[k][i].dat[0];
                    hlm_phase_g[k][i] = hlm_phase_g[k][i] + nqc[k][i].dat[1];
                    std::fprintf(nqcs, "%f\t%e\t%e\n", t_vecg[i], nqc[k][i].dat[0], nqc[k][i].dat[1]);
                }
            }
        }
    }
//    if (DEBUG)
//    {
        std::fclose(nqcs);
//    }
    if (DEBUG)
    {
        char   outputr[256]   = "waveform_postNQC.dat";
        std::FILE* waveform_postNQC   = std::fopen(outputr, "w");
        int j                 = 0;
        int N                 = hlm_ampl_g[1].size();
        
        for (j=0;j<N;j++)
        {
            std::fprintf(waveform_postNQC, "%f\t%e\t%e\n", (double)j*dt, hlm_ampl_g[1][j], hlm_phase_g[1][j]);
        }
        std::fclose(waveform_postNQC);
    }
    /** Define a time vector for each multipole
        These will be cut by the ringdown, where
        each multipole has its own starting time */
    
    vector<vector<double> > t_g(35);
    for (int k=35; k--; )
    {
        t_g[k] = t_vecg;
    }
    
    /** Ringdown attachment */
    if (params.flags.tidal==0)
    {
        ringdown(params,t_g,OmgOrb_vecg,hlm_ampl_g,hlm_phase_g);
    }
    /** All multipoles will now have size N+Nringdown */
    /** Multipole for which no ringdown model is available will be filled with 0s */
    /** We pick the index 1 since it is the 22 mode and it is always computed */
    
    int N = hlm_ampl_g[1].size();

    /** Allocate hplus and hcross */
    Waveform *hplus_out = (Waveform *)malloc(sizeof(Waveform));
    if (hplus_out == NULL)
    {
        printf("ERROR allocating hplus.\n");
        exit(-1);
    }
    hplus_out->data = (double *)malloc(N*sizeof(double));
    if (hplus_out->data == NULL)
    {
        printf("ERROR allocating hplus->data.\n");
        exit(-1);
    }
    hplus_out->length = N;
    Waveform *hcross_out = (Waveform *)malloc(sizeof(Waveform));
    if (hcross_out == NULL)
    {
        printf("ERROR allocating hcross.\n");
        exit(-1);
    }
    hcross_out->data = (double *)malloc(N*sizeof(double));
    if (hcross_out->data == NULL)
    {
        printf("ERROR allocating hcross->data.\n");
        exit(-1);
    }
    hcross_out->length = N;
    
    /** Set them to zero initially */
    
    memset(hplus_out->data, 0, N*sizeof(double));
    memset(hcross_out->data, 0, N*sizeof(double));
    
    /** Spherical harmonics projection **/
    /** construct hplus and hcross **/
    /** h22 = 1/R * (nu*M)*G/c^2 h_code_output */

    double mtot_m = 1.;
    double amplitude_prefactor = 1.;    
    if (!(params.flags.geometric_units))
    {
      mtot_m = (m1+m2)*MSUN_M;
      amplitude_prefactor = params.nu*mtot_m/(distance*MPC_M);
    } 

    if (params.flags.multipoles == 1)
    {
        for (i=0; i<N; i++)
        {
            hplus_out->data[i]   = hlm_ampl_g[lm][i];
            hcross_out->data[i]  = hlm_phase_g[lm][i];
        }

    }
    else
    {
        for (int k=35; k--; )
        {
            if (k==1)
            {
                double Y_real, Y_imag;
                spinsphericalharm(&Y_real, &Y_imag, -2, L[k], M[k], polarisation, inclination);
                
                /** there is a MINUS SIGN in the phase h = A exp(-i phase) **/
                for (i=0; i<N; i++)
                {
                    double Aki = hlm_ampl_g[k][i]*amplitude_prefactor;
                    double cosPhi = cos(hlm_phase_g[k][i]);
                    double sinPhi = -sin(hlm_phase_g[k][i]);
                    hplus_out->data[i] += Aki*(cosPhi*Y_real - sinPhi*Y_imag);
                    hcross_out->data[i] -= Aki*(cosPhi*Y_imag + sinPhi*Y_real);
                }
            }
        }
    }
    *hplus = hplus_out;
    *hcross= hcross_out;
}
