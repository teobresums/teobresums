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

#include <math.h>

#include "s_HorizonFlux.h"

double s_HorizonFlux(
                     double x,
                     double Heff,
                     double jhat,
                     double nu,
                     double X1,
                     double X2,
                     double chi1,
                     double chi2
){

    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x3*x;
    double x5 = x4*x;
    double v5 = sqrt(x5);

    double cv5[2];
    double cv8[2];

    /** Coefficients of the v^5 term (Alvi leading order) */
    cv5[0] = -1./4.*chi1*(1.+3.*chi1*chi1)*X1*X1*X1;
    cv5[1] = -1./4.*chi2*(1.+3.*chi2*chi2)*X2*X2*X2;

    /** Coefficients of the v^8=x^4 term */
    cv8[0] = 0.5*(1.+sqrt(1.-chi1*chi1))*(1.+3.*chi1*chi1)*X1*X1*X1*X1;
    cv8[1] = 0.5*(1.+sqrt(1.-chi2*chi2))*(1.+3.*chi2*chi2)*X2*X2*X2*X2;

    double FH22_S = (cv5[0]+cv5[1])*v5;
    double FH22   = (cv8[0]+cv8[1])*x4;
    double FH21   =  0.0;

    /** Newton-normalized horizon flux: use only l=2 fluxes */
    double hatFH  = FH22_S + FH22 + FH21;
        
    return hatFH;
}

