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

#include "math.h"
#include <gsl/gsl_math.h>
#include <stdlib.h>

#include "constants.h"
#include "multipole_index.h"
#include "Tlm.h"

void Tlm(
                double MTlm[],    /** OUTPUT Dimension: 35 */
                const double w
){

    /* Factorials evaluated */
    double f14[] = {1.,         1.,          2.,
                    6.,         24.,         120.,
                    720.,       5040.,       40320.,
                    362880.,    3628800.,    39916800.,
                    479001600., 6227020800., 87178291200.};
    double hhatk, x2, y, prod;
    const int kmax = 35;
    
    for (int i=kmax; i--; )
    {
        hhatk = M[i] * w;   
        x2    = 4.*hhatk*hhatk;
        prod  = 1.;

        for (int j=1; j <= L[i]; j++ )
        {
            prod *= ( j*j + x2 );
        }
        
        y  = 4.*pi*hhatk;
        y /= ( 1. - exp(-y) );
        
        /** Put the three pieces together */
        MTlm[i] = sqrt( 1./(f14[L[i]]*f14[L[i]]) * y * prod );
    }
}
