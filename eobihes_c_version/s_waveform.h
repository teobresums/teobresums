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

#ifndef _s_waveform_h
#define _s_waveform_h

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include "Metric.h"
#include <vector>
#include "cmath"
#include "flux.h"

#include "hlm.h"

#include <ios>
#include <fstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "hlm.h"
#include "interpolator.h"
#include "multipole_index.h"
#include "input_struc.h"
#include "s_GS.h"


#include <limits>

typedef std::numeric_limits< double > dbl;

using namespace::std;

vector<gsl_complex> s_waveform(double t, const double y[], void *params, double &Omg, double &Omg_orb, double &A, double &ddotr);

#endif /* _s_waveform_h */
