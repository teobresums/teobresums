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

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <cmath>
#include <vector>
#include "ringdown_match.h"

using namespace::std;

gsl_complex ringdown_match(double x, int k, vector<double> a1, vector<double> a2, vector<double> a3, vector<double> a4, vector<double> b1, vector<double> b2, vector<double> b3, vector<double> b4, vector<gsl_complex> sigma){

    gsl_complex psi;
                    
    double amp   =  ( a1[k] * tanh(a2[k]*x +a3[k]) + a4[k] ) ;
    double phase = -b1[k]*log((1. + b3[k]*exp(-b2[k]*x) + b4[k]*exp(-2.*b2[k]*x))/(1.+b3[k]+b4[k]));
                    
    psi.dat[0] = amp * exp(-sigma[k].dat[0]*x);
    psi.dat[1] = - (phase - sigma[k].dat[1]*x); //NOTE: minus sign in front by convention
    
    return psi;
    
}

