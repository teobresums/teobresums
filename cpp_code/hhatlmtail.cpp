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

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <vector>
#include <cmath>

#include "TEOBResum.h"

using namespace::std;

vector<gsl_complex> hhatlmtail(const double Omega, const double Hreal, const double bphys, const int L[], const int M[])
{

/** EOBhhatlmTail Computes the tail contribution to the resummed wave.
  *
  *   tlm = EOBTail(L,M, Omega,E, bphys)
  *
  *   Reference(s)
  *   Damour, Iyer & Nagar, PRD 79, 064004 (2009)
  *
  */
 
    int kmax  = 35;
    const double pi = M_PI;
    double k;
    double hhatk;
    
    gsl_sf_result num_rad;
    gsl_sf_result num_phase;
    gsl_sf_result denom_rad;
    gsl_sf_result denom_phase;
    
    double ratio_rad;
    double ratio_ang;
    double tlm_rad;
    double tlm_phase;
    vector<gsl_complex> tlm(kmax);
    
    for (int i=kmax; i--;)
    {
            k     = M[i] * Omega;
            hhatk = k * Hreal;
        
            gsl_sf_lngamma_complex_e(L[i] + 1., -2.*hhatk, &num_rad, &num_phase);
            gsl_sf_lngamma_complex_e(L[i] + 1., 0., &denom_rad, &denom_phase);
        
            ratio_rad = num_rad.val-denom_rad.val;
            ratio_ang = num_phase.val-0.;
            
            tlm_rad   = ratio_rad + pi * hhatk;
            tlm_phase = ratio_ang + 2.*hhatk*log(2.*k*bphys);
            
            tlm[i].dat[0] = exp(tlm_rad);
            tlm[i].dat[1] = tlm_phase;
    }
        return tlm;
}
