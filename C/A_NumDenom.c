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

#include <math.h>

#include "A_NumDenom.h"

void A_NumDenom(
    double        frac[],            /** OUTPUT Dimension: 4;  Fraction ... */
    const double  r,                 /** Radius ... */
    const double  a[],               /** ... */
    const double  nu                 /** Mass ratio */
    ){
    
    /** Shorthands */
    const double u  = 1./r;
    const double u2 = u*u;
    const double u3 = u2*u;
    const double u4 = u3*u;
    const double u5 = u4*u;

    /**  Point-mass PN coefs */
    const double a3  = a[0];
    const double a4  = a[1];
    const double a5  = a[4];
    const double a6  = a[7];
    const double a5l = a[3];
    const double a6l = a[6];

    /**  Coefficients of the denominator and numerator of the Pade' A function */
    const double C1  = (8.*a3 + 3.*a4 + a5);
    const double C2  = (4.*a4 + 2.*a5 + a6 + a3*(8. - a3*nu));
    const double C3  = (12.*a3 + 4.*a4 + a5);
    const double C4  = (-32. + C3*nu);
    const double C6  = (4.*a3 + a4);
    const double C5  = (16.*a6 + (2.*a4 + a5)*(8.*a3 + 2.*a4 + a5)*nu - C6*a6*nu + a3*a3*C6*nu*nu);
    const double C7  = (64. - (4.*C1 + a6)*nu + a3*a3*nu*nu);
    const double C8  = a6*(8. - a3*nu) + nu*(8.*a3*a3 + 16.*a3*a4 + 4.*a4*a4 + a3*a3*a3*nu) + a5*(16. + a4*nu);
    const double C9  = 4.*(4.*a4 + 2.*a5 + a6) + a3*(8.*a3 + 4.*a4 + a5)*nu;
    const double C10 = (-8 + a3*nu);
    const double C11 = (-32. + C3*nu);
    const double C12 = (-4.*a4 - 2.*a5 - a6 + a3*C10);

    frac[0] = 1. + (C7*u)/C11;
    frac[1] = (-32. + nu*(12.*a3 + 4.*a4 + a5 - C2*u - 2.*C2*u2 - (C9)*u3 - (C8)*u4 - C5*u5))/C11;
    frac[2] = (-((4.*a5l + a6l)*nu*C4) + C4*C7 + a5l*nu*(-64 + nu*(12.*a4 + 4.*a5 + a6 + a3*(32 - a3*nu))))/(C4*C4);
    frac[3] = (nu*(-((2.*a5l + a6l)*C4) + a5l*nu*C2 + C4*C12 - 2.*(2.*a5l + a6l)*C4*u + 2.*a5l*nu*C2*u + 4.*C4*C12*u + a5l*nu*(C9)*u2 - 3.*(C9)*C4*u2 - C4*(4.*a6l + a5l*(8 + a3*nu))*u2 + a5l*nu*(C8)*u3 - 4.*C4*(C8)*u3 + C4*(a6l*C10 - a5l*(16 + a4*nu))*u3 + a5l*nu*C5*u4 - 5.*C4*C5*u4 + C4*(-2.*(4.*a3 + 2.*a4 + a5)*a5l*nu + a6l*(-16 + 4.*a3*nu + a4*nu))*u4))/(C11*C11);
}
