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

#include <cmath>
#include "TEOBResum.h"

double HealyBBHFitRemnant(double chi1,double chi2, double q)
{

    /**
    
    * Computes the mass and angular momentum of the final black hole
    * implementing the fits of Healey, Lousto and Zochlower (HLZ),
    * arXiv: 1406.7295, published as PRD 90, 104004 (2014)
    *
    * Usage: [a, Mbh] = EOB_HealyBBHFitRemnant(q, chi1, chi2)
    *
    * WARNING: the formula uses the convention that M2 > M1, so that
    *          chi2 should refer to the black hole with the largest
    *          mass. In the EOB code, this is given by chi1, since
    *          in EOB code we use the convention that M1 > M2
    *
    *          Here it is q=M2/M1, with M2>M1
    *
    * Improved with (Eisco, Jisco) + iterative procedure 23/02/2016
    * parameters (TABLE VI)
    */

    /** Final mass:                    Angular momentum: */
     
    double M0  =  0.951507;            double L0  =  0.686710;
    double K1  = -0.051379;            double L1  =  0.613247;
    double K2a = -0.004804;            double L2a = -0.145427;
    double K2b = -0.054522;            double L2b = -0.115689;
    double K2c = -0.000022;            double L2c = -0.005254;
    double K2d =  1.995246;            double L2d =  0.801838;
    double K3a =  0.007064;            double L3a = -0.073839;
    double K3b = -0.017599;            double L3b =  0.004759;
    double K3c = -0.119175;            double L3c = -0.078377;
    double K3d =  0.025000;            double L3d =  1.585809;
    double K4a = -0.068981;            double L4a = -0.003050;
    double K4b = -0.011383;            double L4b = -0.002968;
    double K4c = -0.002284;            double L4c =  0.004364;
    double K4d = -0.165658;            double L4d = -0.047204;
    double K4e =  0.019403;            double L4e = -0.053099;
    double K4f =  2.980990;            double L4f =  0.953458;
    double K4g =  0.020250;            double L4g = -0.067998;
    double K4h = -0.004091;            double L4h =  0.001629;
    double K4i =  0.078441;            double L4i = -0.066693;



    /** Parameters */
    double nu      = q/((1.+q)*(1.+q));
    /** Masses: convention here is that m2>m1 */
    double X2      = 0.5*(1.+sqrt(1.-4*nu));
    double X1      = 1.-X2;
    /** Spin variables */
    double s1      = X1*X1*chi1;
    double s2      = X2*X2*chi2;
    double S       = s1 + s2;
    double S2      = S*S;
    double S3      = S*S2;
    double S4      = S2*S2;
    double Delta   = X1/X2*s2 - X2/X1*s1 + s2 - s1;
    double Delta2  = Delta*Delta;
    double Delta3  = Delta*Delta2;
    double Delta4  = Delta2*Delta2;
    /** Mass ratio variables */
    double deltam  = -sqrt(1-4*nu); // X1 - X2
    double deltam2 =  deltam*deltam;
    double deltam3 =  deltam*deltam2;
    double deltam4 =  deltam*deltam3;
    double deltam6 =  deltam2*deltam4;

    /** Initialize the angular momentum */
    double a0 = s1 + s2;
    int a0_sign = 0.;
  
    if (a0==0)
    {
        a0_sign=0;
    }
    else if (a0>0)
    {
        a0_sign=1;
    }
    else if (a0<0)
    {
        a0_sign=-1;
    }

    /** Set-up an interative procedure to compute properly the "isco" quantities */

    double a2;
    double Z1;
    double Z2;
    double risco;
    double uisco;
    double Eisco;
    double Jisco;
    double abh;
    double Mbh=0.;
    
    for(int i=20; i--;)
    {
    
        a2     = a0*a0;
        Z1     = 1 + cbrt(1-a2)*(cbrt(1+a0) + cbrt(1-a0));
        Z2     = sqrt(3*a2 + Z1*Z1);
        risco  = 3 + Z2 - a0_sign*sqrt((3-Z1)*(3+Z1+2.*Z2));
        uisco  = 1./risco;
        Eisco  = (1 - 2.*uisco + a0*sqrt(uisco*uisco*uisco))/sqrt(1-3*uisco + 2*a0*sqrt(uisco*uisco*uisco));
        Jisco  = 2./(sqrt(3.*risco))*(3.*sqrt(risco)-2.*a0);
        
        /** Dimensionless spin: J/Mbh^2 */
        abh = (4*nu)*(4*nu)*(L0 + L1*S + L2a*Delta*deltam + L2b*S2 + L2c*Delta2 + L2d*deltam2 + L3a*Delta*S*deltam + L3b*S*Delta2 + L3c*S3 + L3d*S*deltam2 + L4a*Delta*S2*deltam + L4b*Delta3*deltam + L4c*Delta4 + L4d*S4 + L4e*Delta2*S2 + L4f*deltam4 + L4g*Delta*deltam3 + L4h*Delta2*deltam2 + L4i*S2*deltam2) + S*(1+8*nu)*deltam4 + nu*Jisco*deltam6;

        Mbh = (4*nu)*(4*nu)*(M0 + K1*S + K2a*Delta*deltam + K2b*S2 + K2c*Delta2 + K2d*deltam2 + K3a*Delta*S*deltam + K3b*S*Delta2 + K3c*S3 + K3d*S*deltam2 + K4a*Delta*S2*deltam + K4b*Delta3*deltam + K4c*Delta4 + K4d*S4 + K4e*Delta2*S2 + K4f*deltam4 + K4g*Delta*deltam3 + K4h*Delta2*deltam2 + K4i*S2*deltam2) + (1 + nu*(Eisco + 11))*deltam6;

        a0 = abh;
    
    }
    return Mbh;
}
