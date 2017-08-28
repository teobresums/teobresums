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
#include <sys/stat.h>

#include "initial.h"
#include "s_initial.h"
#include "RHS.h"
#include "s_RHS.h"
#include "s_waveform.h"
#include "QNMHybridFitCab.h"
#include "HealyBBHFitRemnant.h"
#include "ringdown.h"
#include "input_struc.h"
#include "multipole_index.h"
#include "interpolate_wf.h"
#include "AdiabLR.h"
#include "read_config.h"
#include "file_names.h"
#include "find_a1a2a3.h"
#include "interp_grid.h"
#include "LALEOB.h"

using namespace::std;

int main (int argc,char* argv[])
{
    double m1 = 1.5;
    double m2 = 1.5;
    double q = m1/m2;
    double chi1 = 0.6;
    double chi2 = 0.6;
    double f_min = 100.;
    double sampling_rate = 2048.;
    double LambdaAl2 = 100.0;
    double LambdaBl2 = 1000.0;
    double distance = 1.0;
    bool NQC = false;
    bool tidal = true;
    bool speedy = true;
    bool RWZ = false;
    int lm = 1;
    int solver_scheme = 0;
    vector<double> waveform = LALEOB( m1,
                           m2,
                           q,
                           chi1,
                           chi2,
                           f_min,
                           sampling_rate,
                           LambdaAl2,
                           LambdaBl2,
                           distance,
                             NQC,
                             tidal,
                             speedy,
                             RWZ,
                              lm,
                             solver_scheme);
    int i=0;
    for (i=0;i<waveform.size();i++) fprintf(stderr,"%d %e\n",i,waveform[i]);
    
    return 0;
}
