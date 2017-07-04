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

#ifndef _find_a1a2a3_h
#define _find_a1a2a3_h

#include <ios>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <limits>

#include "input_struc.h"
#include "multipole_index.h"
#include "s_D1.h"
#include "dtnqc_fit.h"

typedef std::numeric_limits< double > dbl;

using namespace::std;

vector<vector<gsl_complex> > find_a1a2a3_hlm(vector<double> T, vector<double> r, vector<double> w, vector<double> pph, vector<double> pr_star, vector<vector<double> > hlm_phase,vector<double> Omg_orb, vector<vector<double> > A, vector<double> ddotr, void *params);

#endif /* _find_a1a2a3_h */
