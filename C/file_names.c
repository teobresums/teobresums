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

#include <stdio.h>
#include <string.h>

#include "multipole_index.h"
#include "input_struc.h"

#include "file_names.h"

void file_names (char data[512],        /** OUTPUT */
                 char init[512],        /** OUTPUT */
                 char waveNames[35][512],  /** OUTPUT */
                 input *params
){

    int solver_scheme = (*params).solver_scheme;
    int lm            = (*params).lm;
    double q          = (*params).q;
    double r0         = (*params).r0;
    double dt         = (*params).dt;
    double chi1       = (*params).chi1;
    double chi2       = (*params).chi2;

    char wave[35][512];
    const char *conf;
    const char *scheme;
    char QOmg[512];
    char Afunc[512];
    char time_step[512];
    char runReport[512];

     if ((*params).tidal==true){conf = "BNS";}
     else {conf = "BBH";}
    
     switch (solver_scheme)
     {
         case 0:
             scheme="adap";
             break;
         case 1:
             scheme="hybrid";
             break;
         case 2:
             scheme="fixed";
             break;
         default:
             scheme="adap";
             break;
     }

     sprintf(data, "data/dynamics_%s_%i%i_q_%g_r0_%g_dt_%g_chi1_%g_chi2_%g_%s.dat", conf, L[lm], M[lm], q, r0, dt, chi1, chi2, scheme);
     sprintf(init, "data/initialData_%s_%i%i_q_%g_r0_%g_dt_%g_chi1_%g_chi2_%g_%s.dat", conf, L[lm], M[lm], q, r0, dt, chi1, chi2, scheme);

    for (int k=35; k--; )
    {
     sprintf(wave[k], "data/h_%s_%i%i_q_%g_r0_%g_dt_%g_chi1_%g_chi2_%g_%s.dat", conf, L[k], M[k], q, r0, dt, chi1, chi2, scheme);
     strcpy(waveNames[k], wave[k]);
    }

     sprintf(QOmg, "data/Qomg_%s_%i%i_q_%g_r0_%g_dt_%g_chi1_%g_chi2_%g_%s.dat", conf, L[lm], M[lm], q, r0, dt, chi1, chi2, scheme);
     sprintf(Afunc, "data/Afunc_%s_%i%i_q_%g_r0_%g_dt_%g_chi1_%g_chi2_%g_%s.dat", conf, L[lm], M[lm], q, r0, dt, chi1, chi2, scheme);
     sprintf(time_step, "data/Tstep_%s_%i%i_q_%g_r0_%g_dt_%g_chi1_%g_chi2_%g_%s.dat", conf, L[lm], M[lm], q, r0, dt, chi1, chi2, scheme);
     sprintf(runReport, "data/runReport_%s_%i%i_q_%g_r0_%g_dt_%g_chi1_%g_chi2_%g_%s.dat", conf, L[lm], M[lm], q, r0, dt, chi1, chi2, scheme);

}


