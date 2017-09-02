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
#include "HorizonFlux.h"

double HorizonFlux(const double x, const double Heff, const double jhat, const double nu)
{

    /** Compute horizon-absorbed fluxes.
     *
     *
     *
     * Nagar & Akcay, PRD 85, 044025 (2012)
     * Bernuzzi, Nagar & Zenginoglu, PRD 86, 104038 (2012)
     */

    const int kmax = 2;
    double rhoHlm[kmax];
    double FlmHLO[kmax];
    double FlmH[kmax];

    /** Shorthands */
    double nu2 = nu*nu;
    double nu3 = nu*nu2;
    double x2  = x*x;
    double x3  = x*x2;
    double x4  = x*x3;
    double x5  = x*x4;
    double x9  = x4*x5;
    double x10 = x*x9;

    const int k22 = 1;//LM2K(2,2);
    const int k21 = 0;//LM2K(2,1);

    /** The Newtonian asymptotic contribution */
    const double FNewt22 = 32./5.*x5;

    /** Compute leading-order part (nu-dependent) */
    FlmHLO[1] = 32./5.*(1-4*nu+2*nu2)*x9;
    FlmHLO[0] = 32./5.*(1-4*nu+2*nu2)*x10;

    /** Compute rho_lm */
    double c1[kmax];
    double c2[kmax];
    double c3[kmax];
    double c4[kmax];

    c1[0] = 0.58121;
    c2[0] = 1.01059;
    c3[0] = 7.955729;
    c4[0] = 1.650228;
    
    c1[1] = (4.-21.*nu + 27.*nu2 - 8.*nu3)/(4.*(1.-4.*nu+2.*nu2));
    c2[1] =  4.78752;
    c3[1] = 26.760136;
    c4[1] = 43.861478;
    
    rhoHlm[1] = 1. + c1[1]*x + c2[1]*x2 + c3[1]*x3 + c4[1]*x4;
    rhoHlm[0] = 1. + c1[0]*x + c2[0]*x2 + c3[0]*x3 + c4[0]*x4;

    /** Compute horizon multipolar flux (only l=2) */
    const double Heff2 = Heff*Heff;
    const double jhat2 = jhat*jhat;

    FlmH[k22] = FlmHLO[k22] * Heff2 *rhoHlm[k22]*rhoHlm[k22]*rhoHlm[k22]*rhoHlm[k22];
    FlmH[k21] = FlmHLO[k21] * jhat2 *rhoHlm[k21]*rhoHlm[k21]*rhoHlm[k21]*rhoHlm[k21];

    /** Sum over multipoles and normalize to the 22 Newtonian multipole */
    const double hatFH = (FlmH[0]+FlmH[1])/FNewt22;

    return hatFH;
}

