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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#ifndef _read_config_h
#define _read_config_h

#include "input_struc.h"

double logQ(double x);

input process_input_parameters(double m1,
                               double m2,
                               double chi1,
                               double chi2,
                               double f_min,
                               double sampling_rate,
                               double LambdaAl2,
                               double LambdaBl2,
                               bool   NQC,
                               bool   tidal,
                               bool   speedy,
                               bool   RWZ,
                               int    lm,
                               int    solver_scheme
                               );
input read_config(char *fname);

double Yagi13_fit_barlamdel(double barlam2, int ell);
double time_units_conversion(double M, double Srate);
double radius0(double M, double f_start);
#endif /* _read_config_h */
