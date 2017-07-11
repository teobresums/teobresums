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

#include <gsl/gsl_complex.h>
#include <limits>

typedef std::numeric_limits< double > dbl;

gsl_complex** find_a1a2a3(double* T, const int size_T, double* r, double* w, double* pph, double* pr_star, double** hlm_phase, double* Omg_orb, double** A, double* ddotr, void *params);

#endif /* _find_a1a2a3_h */
