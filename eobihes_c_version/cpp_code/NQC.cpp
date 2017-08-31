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

#include <cmath>
#include "NQC.h"

vector<double> NQC(const double r, const double prstar, const double Omega, const double ddotr, int i)
{


    vector<double> n(6);
    
    switch (i)
    {
            // l=2 -------------------------------------------------------------------
            // (2,1)
        case 0:
            // NQC corrections to the modulus
            n[0] = (prstar/(r*Omega))*(prstar/(r*Omega));
            n[1] = ddotr/(r*Omega*Omega);
            n[2] = n[0]*prstar*prstar;
            
            
            //% NQC corrections to the phase
            n[3] = prstar/(r*Omega);
            n[4] = n[3]*cbrt(Omega*Omega);
            n[5] = n[4]*prstar*prstar;
            break;
            // (2,2)
        case 1:
            // NQC corrections to the modulus
            n[0] = (prstar/(r*Omega))*(prstar/(r*Omega));
            n[1] = ddotr/(r*Omega*Omega);
            n[2] = n[0]*prstar*prstar;
            
            
            //% NQC corrections to the phase
            n[3] = prstar/(r*Omega);
            n[4] = n[3]*(r*Omega)*(r*Omega);
            n[5] = n[4]*prstar*prstar;
            break;
            // l=3 -------------------------------------------------------------------
            // (3,3)
        case 4:
            // NQC corrections to the modulus
            n[0] = (prstar/(r*Omega))*(prstar/(r*Omega));
            n[1] = ddotr/(r*Omega*Omega);
            n[2] = n[0]*prstar*prstar;
            
            
            //% NQC corrections to the phase
            n[3] = prstar/(r*Omega);
            n[4] = n[3]*cbrt(Omega*Omega);
            n[5] = n[4]*prstar*prstar;
            break;
        default:
            n[0] = 0.;
            n[1] = 0.;
            n[2] = 0.;
            
            n[3] = 0.;
            n[4] = 0.;
            n[5] = 0.;
            break;
    }
    
    return n;
}

