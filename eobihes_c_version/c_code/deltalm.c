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
#include <math.h>
#include "deltalm.h"
#include "constants.h"

double* deltalm(const double Hreal, const double Omega, const double nu)
{

    /*
    EOBdeltalm Residual phase corrections delta_{lm} up to l=m=5.
    If nu=0, the delta_lmm are written in Taylor-expanded form and all terms up to 4.5PN accuracy are included.
    nu~=0,
    (i) the 4.5PN, test-mass terms are excluded because of the too
    large PN-gap with the nu-dependent terms
    (ii) the deltalm are replaced by suitable Pade' approximants for
    multipoles (2,2), (2,1), (3,3), (3,1)
    The l=m=2 residual phase includes the 3.5PN (nu-dependent) correction
    obtained in Faye et al.

    Reference(s)
    Damour, Iyer & Nagar, PRD 79, 064004 (2008)
    Fujita & Iyer, PRD 82 044051 (2010)
    Faye et al., Class. Q. Grav. 29 175004 (2012)
    Damour, Nagar & Bernuzzi, PRD 87, 084035 (2013)

     TODO: this routine requires optimization
     - precompute coefficients c(nu)
     - evaluate efficiently polynomials
    */
    
    /** Useful shorthands*/
    double pi2    = pi*pi;
    double nu2    = nu*nu;
    double y      = gsl_pow_int(Hreal*Omega,2./3.);
    double sqrt_y = sqrt(y);
    double y3     = y*y*y;
    double y32    = Hreal*Omega;
    
    /** Leading order contributions*/
    double delta22LO = 7./3.   * y32;
    double delta21LO = 2./3.   * y32;
    double delta33LO = 13./10. * y32;
    double delta31LO = 13./30. * y32;
    
    int kmax = 35;
    static double deltalmvec[35];
    for (int i=kmax; i--; ) {deltalmvec[i]=0.;}
    
    double num;
    double den;
    
    // Residual phases in Pade-resummed form when possible


    /** l=2 ------------------------------------------------------------------*/
    
    /** Pade(1,2) approximant */
    num           = 69020.*nu + 5992.*pi*sqrt_y;
    den           = 5992.*pi*sqrt_y + 2456.*nu*(28.+493.*nu* y);
    deltalmvec[0] = delta21LO*num/den;
    
    /** Pade(2,2) approximant */
    num           = (808920.*nu*pi*sqrt(y) + 137388.*pi2*y + 35.*nu2*(136080. + (154975. - 1359276.*nu)*y));
    den           = (808920.*nu*pi*sqrt(y) + 137388.*pi2*y + 35.*nu2*(136080. + (154975. + 40404.*nu)*y));
    deltalmvec[1] = delta22LO*num/den;

    
    /** l=3 ------------------------------------------------------------------*/

    /** Pade(1,2) approximant */
    num           = 4641.*nu + 1690.*pi*sqrt_y;
    den           = num + 18207.*nu2*y;
    deltalmvec[2] = delta31LO*num/den;
    
    /** Taylor-expanded form */
    num           = 1.  + 94770.*pi/(566279.*nu)*sqrt_y;
    den           = num + 80897.* nu/3159.*y;
    deltalmvec[3] = (10.+33.*nu)/(15.*(1.-3.*nu)) * y32 + 52./21.*pi*y3;
    
    /** Pade(1,2) approximant */
    deltalmvec[4] = delta33LO*num/den;


    /** l=4 ------------------------------------------------------------------*/
    deltalmvec[5] =   (2.+507.*nu)/(10.*(1.-2.*nu))*y32   + 1571./3465.*pi*y3;
    deltalmvec[6] =  7.*(1.+6.*nu)/(15.*(1.-3.*nu))*y32   + 6284./3465.*pi*y3;
    deltalmvec[7] = (486.+4961.*nu)/(810.*(1.-2.*nu))*y32 + 1571./385.*pi*y3;
    deltalmvec[8] =  (112.+219.*nu)/(120.*(1.-3.*nu))*y32 + 25136./3465.*pi*y3;

    /** l=5 ------------------------------------------------------------------*/

    deltalmvec[9] = (96875. + 857528.*nu)/(131250.*(1.-2.*nu))*y32;
        
    return deltalmvec;

}

