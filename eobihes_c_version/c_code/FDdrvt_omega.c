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

#include "FDdrvt_omega.h"

void FDdrvt_omega(
        double d1f[],          /** OUTPUT Dimension: size_f-2 */
        double f[],
        const int size_f,
        double dt
    ){
    
    const double oodt = 1./dt;
    const double c    = 1./12.;
    
    for (long int i=0; i<size_f-2; i++)
    {
        switch (i)
        {
            case 0:
                d1f[i] = c*(-25.*f[i]+48.*f[i+1]-36.*f[i+2]+16.*f[i+3]-3.*f[i+4])*oodt;
                //d2f[i] = c*(45*f[i]-154*f[i+1]+214*f[i+2]-156*f[i+3]+61*f[i+4]-10*f[i+5])*oodt2;
                break;
            case 1:
                d1f[i] = c*(-3.*f[i-1]-10.*f[i]+18.*f[i+1]-6.*f[i+2]+f[i+3])*oodt;
                //d2f[i] = c*(10*f[i-1]-15*f[i]-4*f[i+1]+14*f[i+2]-6*f[i+3]+f[i+4])*oodt2;
                break;
                /*
                 case 198:
                 d1f[i] = - c*(-3.*f[i+1]-10.*f[i]+18.*f[i-1]-6.*f[i-2]+f[i-3])*oodt;
                 //d2f[i] = c*(10*f[i+1]-15*f[i]-4*f[i-1]+14*f[i-2]-6*f[i-3]+f[i-4])*oodt2;
                 break;
                 case 199:
                 d1f[i] = - c*(-25.*f[i]+48.*f[i-1]-36.*f[i-2]+16.*f[i-3]-3.*f[i-4])*oodt;
                 //d2f[i] = c*(45*f[i]-154*f[i-1]+214*f[i-2]-156*f[i-3]+61*f[i-4]-10*f[i-5])*oodt2;
                 break;
                 */
            default:
                d1f[i] = c*(8.*(f[i+1]-f[i-1]) - f[i+2] + f[i-2])*oodt;
                //d2f[i] = c*(-30*f[i]+16*(f[i+1]+f[i-1])-(f[i+2]+f[i-2]))*oodt2;
                break;
        }
    }
}
