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

#include "TEOBResumS.h"

using namespace::std;

#define USAGE "\n\
USAGE:\n\
\t./TEOBResumS.x -p <parfile>\n\
\t./TEOBResumS.x [OPTIONS]\n\
"
const char *optstr[] =
{
  "-p"           , "<parfile>",  "reads input parameters from parfile. Overrides all other arguments.", "",
  "-m1"          , "<double>",   "mass of the primary [Msun].", "40",
  "-m2"          , "<double>",   "mass of the secondary [Msun]. ", "40",
  "-chi1"        , "<double>",   "dimensionless spin component along the orbital angular momentum of the primary.", "0",
  "-chi2"        , "<double>",   "dimensionless spin component along the orbital angular momentum of the secondary.", "0",
  "-distance"    , "<double>",   "source distance [Mpc].", "100",
  "-inclination" , "<double>",   "(IOTA) inclination angle [rad].", "0",
  "-polarisation", "<double>",   "(PSI) polarisation angle [rad].", "0",
  "-f_min"       , "<double>",   "starting frequency [Hz].", "20",
  "-srate"       , "<double>",   "sampling rate [Hz].", "4096",
  "-lambda1_l2"  , "<double>",   "l=2 tidal deformability for body 1 (Lambda/M^5). Only if tidal corrections are enabled.", "0",
  "-lambda2_l2"  , "<double>",   "l=2 tidal deformability for body 2 (Lambda/M^5). Only if tidal corrections are enabled.", "0",
  "-lambda1_l3"  , "<double>",   "l=3 tidal deformability for body 1 (Lambda/M^5). Only if tidal corrections are enabled.", "0",
  "-lambda2_l3"  , "<double>",   "l=3 tidal deformability for body 2 (Lambda/M^5). Only if tidal corrections are enabled.", "0",
  "-lambda1_l4"  , "<double>",   "l=4 tidal deformability for body 1 (Lambda/M^5). Only if tidal corrections are enabled.", "0",
  "-lambda2_l4"  , "<double>",   "l=4 tidal deformability for body 2 (Lambda/M^5). Only if tidal corrections are enabled.", "0",
  "-tidal"       , "<int>",      "enable tidal corrections.", "0 (false)",
  "-speedy"      , "<int>",      "faster tails calculations.", "1 (true)",
  "-dynamics"    , "<int>",      "output dynamics evolution.", "0 (false)",
  "-RW"          , "<int>",      "Regge-Wheeler-Zerilli potential.", "0 (false)",
  "-multipoles"  , "<int>",      "enable single multipole output, in geometrical units.", "0 (false)",
  "-lm"          , "<int>",      "index for the output multipole. Requires multipoles output format.", "-1",
  "-Yagi_fits"   , "<int>",      "enable Yagi fits for Lambda_l=3,4.", "0 (false)",
  "-output"      , "<filename>", "output file. If multipoles is enable will contain t/M amplitude phase. Otherwise t(s) h+ hx.", "'waveform.dat'"
};

int main (int argc, char* argv[])
{
    double sampling_rate = 4096.;
    double m1            = 40.0;
    double m2            = 40.0;
    double q             = m1/m2;
    double chi1          = 0.0;
    double chi2          = 0.0;
    double f_min         = 20.;
    double dt            = 1./sampling_rate;
    double LambdaAl2     = 0.0;
    double LambdaBl2     = 0.0;
    double LambdaAl3     = 0.0;
    double LambdaBl3     = 0.0;
    double LambdaAl4     = 0.0;
    double LambdaBl4     = 0.0;
    double distance      = 100;
    double inclination   = 0.0;
    double polarisation  = 0.0;
    int    lm            = -1;
    char   output[256]   = "waveform.dat";
    char   parfile[256]  = "";
    TEOBResumFlags flags;
    SetDefaultFlagsValues(&flags);
    flags.set = 1;
    
    if (argc < 2)
    {
        fprintf(stderr,USAGE);
        fprintf(stderr,"\nOPTIONS:\n");
        for (int i = 0; i < (20*4); i=i+4)
            fprintf(stderr,"\t%-20s %-10s %s [%s]\n",optstr[i],optstr[i+1],optstr[i+2],optstr[i+3]);
        exit(0);
    }
   
    for (int i = 1; i < argc; i++)
    {
        if ((strcmp(argv[i],"-h")==0)||(strcmp(argv[i],"-help")==0))
        {
            fprintf(stderr,USAGE);
	    fprintf(stderr,"\nOPTIONS:\n");
	    for (int i = 0; i < (20*4); i=i+4)
	      fprintf(stderr,"\t%-20s %-10s %s [%s]\n",optstr[i],optstr[i+1],optstr[i+2],optstr[i+3]);
            exit(0);
        }
        else if (strcmp(argv[i],"-p")==0)
        {
            sprintf(parfile,"%s",argv[i+1]);
            printf("Reading parfile: %s\n",parfile);

            TEOBResumParams params = read_config(parfile);
            q = params.q;
            m1 = params.mtot * q/(1.+q);
            m2 = params.mtot - m1;
            chi1 = params.chi1;
            chi2 = params.chi2;
            dt = params.dt;
            flags.RWZ  = params.flags.RWZ;
            flags.solver_scheme = params.flags.solver_scheme;
            flags.tidal = params.flags.tidal;
            flags.speedy = params.flags.speedy;
            flags.Yagi_fits = params.flags.Yagi_fits;
            flags.multipoles = params.flags.multipoles;
            flags.geometric_units = params.flags.geometric_units;

//            flags->spin = params.flags.spin; FIX ME, NOW IT WORKS A CDC
            lm = params.lm;
            LambdaAl2 = params.LambdaAl2;
            LambdaBl2 = params.LambdaBl2;
            if (flags.Yagi_fits==0)
            {
                LambdaAl3 = params.LambdaAl3;
                LambdaBl3 = params.LambdaBl3;
                LambdaAl4 = params.LambdaAl4;
                LambdaBl4 = params.LambdaBl4;
            }

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
            dt = 1./sampling_rate;
        }
        else if (strcmp(argv[i],"-tidal")==0)
        {
            flags.tidal = 1;
            printf("tidal = true\n");
        }
        else if (strcmp(argv[i],"-speedy")==0)
        {
            flags.speedy = 1;
            printf("speedy = true\n");
        }
        else if (strcmp(argv[i],"-RWZ")==0)
        {
            flags.RWZ = 1;
            printf("RWZ = true\n");
        }
        else if (strcmp(argv[i],"-Yagi_fits")==0)
        {
            flags.Yagi_fits = 1;
            printf("Yagi_fits = true\n");
        }
        else if (strcmp(argv[i],"-lambda1_l2")==0)
        {
            LambdaAl2 = atof(argv[i+1]);
            printf("lambda1_l2: %f\n",LambdaAl2);
        }
        else if (strcmp(argv[i],"-lambda2_l2")==0)
        {
            LambdaBl2 = atof(argv[i+1]);
            printf("lambda2_l2: %f\n",LambdaBl2);
        }
        else if (strcmp(argv[i],"-lambda1_l3")==0)
        {
            LambdaAl3 = atof(argv[i+1]);
            printf("lambda1_l3: %f\n",LambdaAl3);
        }
        else if (strcmp(argv[i],"-lambda2_l3")==0)
        {
            LambdaBl3 = atof(argv[i+1]);
            printf("lambda2_l3: %f\n",LambdaBl3);
        }
        else if (strcmp(argv[i],"-lambda1_l4")==0)
        {
            LambdaAl4 = atof(argv[i+1]);
            printf("lambda1_l4: %f\n",LambdaAl4);
        }
        else if (strcmp(argv[i],"-lambda2_l4")==0)
        {
            LambdaBl4 = atof(argv[i+1]);
            printf("lambda2_l4: %f\n",LambdaBl4);
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
        else if (strcmp(argv[i],"-output")==0)
        {
            sprintf(output,"%s",argv[i+1]);
        }
        else if (strcmp(argv[i],"-multipoles")==0)
        {
            flags.multipoles = 1;
            printf("multipoles = true\n");
        }
        else if (strcmp(argv[i],"-lm")==0)
        {
            lm = atoi(argv[i+1]);
            printf("lm: %d\n",lm);
        }
        else if (strcmp(argv[i],"-geometric")==0)
        {
            flags.geometric_units = 1;
            printf("geometric = true\n");
        }
    }
    
    if (flags.multipoles == 1 && lm==-1)
    {
        cout << "Need to input also the index of the multipole via the option -lm" << endl;
        exit(-1);
    }
    else
    { 
        
        Waveform *hplus;
        Waveform *hcross;
        
        if (flags.geometric_units==1)
        {
            printf("Will use geometric units and mass rescaled quantities\n");
        }
        else
        {
            double dt_phys = time_units_conversion(m1+m2, dt);
            if (dt_phys > 10.0)
            {
                printf("ERROR! dt = %f is too big and will cause the interpolator to crash when attaching the ringdown.\n",dt_phys);
                printf("Decrease dt in input and retry.\n");
                exit(-1);
            }
            
        }
        
        if (flags.multipoles==1)
        {
            printf("Will output l = %d m = %d waveform\n",L[lm],M[lm]);
        }
        TEOBResumS(&hplus,
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
                    dt,
                    LambdaAl2,
                    LambdaBl2,
                    LambdaAl3,
                    LambdaBl3,
                    LambdaAl4,
                    LambdaBl4,
                    distance,
                    lm,
                    &flags);
        
        std::FILE* f = std::fopen(output, "w");
        int i        = 0;
        int N        = hplus->length;
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
