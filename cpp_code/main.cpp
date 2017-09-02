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
#include <cstdio>
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
#include <sys/stat.h>

#include "TEOBResum.h"

using namespace::std;

int main (int argc, char* argv[])
{
    double m1            = 40.0;
    double m2            = 40.0;
    double q             = m1/m2;
    double chi1          = 0.95;
    double chi2          = 0.2;
    double f_min         = 20.;
    double sampling_rate = 4096.;
    double LambdaAl2     = 0.0;
    double LambdaBl2     = 0.0;
    double distance      = 40;
    double inclination   = 0.0;
    double polarisation  = 0.0;
    int   NQC            = 0;
    int   tidal          = 0;
    int   speedy         = 1;
    int   RWZ            = 0;
    int    dynamics      = 0;
    int    waveform      = 0;
    bool   multipoles    = false;
    int    mult_index    = -1;
    int    lm            = 1;
    int    solver_scheme = 0;
    char   output[256]   = "waveform.dat";
    char   parfile[256]  = "";
    
    
    if (argc < 2) printf("using default values\n");
   
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i],"-p")==0)
        {
            sprintf(parfile,"%s",argv[i+1]);
            printf("found parfile: %s\n",parfile);
            printf("Warning! Will use default values:\n");
            printf("m1 = %f\n",m1);
            printf("f_min = %f\n",f_min);
            printf("iota = %f\n",inclination);
            printf("psi = %f\n",polarisation);
            TEOBResumParams params = read_config(parfile);
            q = params.q;
            m2 = m1/q;
            chi1 = params.chi1;
            chi2 = params.chi2;
            NQC  = params.flags.NQC;
            RWZ  = params.flags.RWZ;
            solver_scheme = params.solver_scheme;
            tidal = params.flags.tidal;
            speedy = params.flags.speedy;
            lm = params.lm;
            LambdaAl2 = params.LambdaAl2;
            LambdaBl2 = params.LambdaBl2;
            break;
        }
        else if (strcmp(argv[i],"-m1")==0)
        {
            m1 = atof(argv[i+1]);
            printf("m1: %f\n",m1);
        }
        else if (strcmp(argv[i],"-m2")==0)
        {
            m2 = atof(argv[i+1]);
            printf("m2: %f\n",m2);
        }
        else if (strcmp(argv[i],"-chi1")==0)
        {
            chi1 = atof(argv[i+1]);
            printf("chi1: %f\n",chi1);
        }
        else if (strcmp(argv[i],"-chi2")==0)
        {
            chi2 = atof(argv[i+1]);
            printf("chi2: %f\n",chi2);
        }
        else if (strcmp(argv[i],"-f_min")==0)
        {
            f_min = atof(argv[i+1]);
            printf("f_min: %f\n",f_min);
        }
        else if (strcmp(argv[i],"-srate")==0)
        {
            sampling_rate = atof(argv[i+1]);
            printf("srate: %f\n",sampling_rate);
        }
        else if (strcmp(argv[i],"-tidal")==0)
        {
            tidal = true;
            NQC = false;
            printf("tidal = true\n");
        }
        else if (strcmp(argv[i],"-NQC")==0)
        {
            NQC = true;
            printf("NQC = true\n");
        }
        else if (strcmp(argv[i],"-speedy")==0)
        {
            speedy = true;
            printf("speedy = true\n");
        }
        else if (strcmp(argv[i],"-RWZ")==0)
        {
            RWZ = true;
            printf("RWZ = true\n");
        }
        else if (strcmp(argv[i],"-lambda1")==0)
        {
            LambdaAl2 = atof(argv[i+1]);
            printf("lambda1: %f\n",LambdaAl2);
        }
        else if (strcmp(argv[i],"-lambda2")==0)
        {
            LambdaBl2 = atof(argv[i+1]);
            printf("lambda2: %f\n",LambdaBl2);
        }
        else if (strcmp(argv[i],"-distance")==0)
        {
            distance = atof(argv[i+1]);
            printf("distance: %f\n",distance);
        }
        else if (strcmp(argv[i],"-inclination")==0)
        {
            inclination = atof(argv[i+1]);
            printf("inclination: %f\n",inclination);
        }
        else if (strcmp(argv[i],"-polarisation")==0)
        {
            polarisation = atof(argv[i+1]);
            printf("polarisation: %f\n",polarisation);
        }
        else if (strcmp(argv[i],"-o")==0)
        {
            sprintf(output,"%s",argv[i+1]);
        }
        else if (strcmp(argv[i],"-multipoles")==0)
        {
            multipoles = true;
            printf("multipoles = true\n");
        }
        else if (strcmp(argv[i],"-mult_index")==0)
        {
            mult_index = atoi(argv[i+1]);
            printf("mult_index: %d\n",mult_index);
        }
    }
    

    
    if (multipoles == true && mult_index!=-1)
    {
        Waveform *ampl;
        Waveform *phase;
        
        XLALSimIMRTEOBIHES_single_mode(&ampl,
                                       &phase,
                                       m1,
                                       m2,
                                       0.0,
                                       0.0,
                                       chi1,
                                       0.0,
                                       0.0,
                                       chi2,
                                       inclination,
                                       polarisation,
                                       f_min,
                                       sampling_rate,
                                       LambdaAl2,
                                       LambdaBl2,
                                       distance,
                                       NQC,
                                       tidal,
                                       speedy,
                                       RWZ,
                                       dynamics,
                                       waveform,
                                       lm,
                                       solver_scheme,
                                       mult_index);
        
        std::FILE* f = std::fopen(output, "w");
        int i        = 0;
        int N        = ampl->length;
        double dt    = 1./sampling_rate;
        for (i=0;i<N;i++)
        {
            std::fprintf(f, "%f\t%e\t%e\n", i*dt, ampl->data[i], phase->data[i]);
        }
        std::fclose(f);
        free(ampl->data);
        free(phase->data);
        free(ampl);
        free(phase);
    }
    else if (multipoles == true && mult_index==-1)
    {
        cout << "Need to input also the index of the multipole" << endl;
        exit(-1);
    }
    else
    { 
        
        Waveform *hplus;
        Waveform *hcross;
        
        XLALSimIMRTEOBIHES(&hplus,
                           &hcross,
                           m1,
                           m2,
                           0.0,
                           0.0,
                           chi1,
                           0.0,
                           0.0,
                           chi2,
                           inclination,
                           polarisation,
                           f_min,
                           sampling_rate,
                           LambdaAl2,
                           LambdaBl2,
                           distance,
                           NQC,
                           tidal,
                           speedy,
                           RWZ,
                           dynamics,
                           waveform,
                           lm,
                           solver_scheme);
        
        std::FILE* f = std::fopen(output, "w");
        int i        = 0;
        int N        = hplus->length;
        double dt    = 1./sampling_rate;
        for (i=0;i<N;i++)
        {
            std::fprintf(f, "%f\t%e\t%e\n", i*dt, hplus->data[i], hcross->data[i]);
        }
        std::fclose(f);
        free(hplus->data);
        free(hcross->data);
        free(hplus);
        free(hcross);
        
    }
    


    return 0;
}
