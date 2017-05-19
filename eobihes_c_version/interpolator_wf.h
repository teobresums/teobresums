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

#ifndef _interpolator_waveform_h
#define _interpolator_waveform_h

#include <stdio.h>
#include <vector>
#include <ios>
#include <fstream>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_errno.h>

using namespace::std;

int interpolate_wf(double dt, vector<double> t_vec, vector<double> hlm_rad, vector<double> hlm_phase, vector<double> Momg_vec, bool waveform_flag, std::ofstream& wave, std::ofstream& Qomg, double Mbh);

#endif /* interpolator_waveform_h */
