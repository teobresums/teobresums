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

#include <fstream>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <ios>
#include <limits>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#include "input_struc.h"
#include "Metric.h"
#include "multipole_index.h"
#include "s_A5PNlog.h"
#include "s_ddotr.h"
#include "s_Flux.h"
#include "s_get_rc.h"
#include "s_GS.h"
#include "s_Metric.h"


typedef std::numeric_limits< double > dbl;

double s_ddotr(double t, double r, double pph, double prstar, void *params)
{
    
    double nu         = (*(input *)params).nu;
    double S1         = (*(input *)params).S1;
    double S2         = (*(input *)params).S2;
    double chi1       = (*(input *)params).chi1;
    double chi2       = (*(input *)params).chi2;
    double X1         = (*(input *)params).X1;
    double X2         = (*(input *)params).X2;
    double c3         = (*(input *)params).cN3LO;
    bool   tidal_flag = (*(input *)params).tidal;
    
    // Kerr parameter
    double a1  = (*(input *)params).a1; //a1 = X1*chi1
    double a2  = (*(input *)params).a2; //a2 = X2*chi2
    double aK2 = (*(input *)params).aK2;

    // shorthands
    double prstar2 = prstar*prstar;
    double prstar3 = prstar2*prstar;
    double prstar4 = prstar3*prstar;
    double pphi2   = pph*pph;
    
    double z3 = 2.*nu*(4.-3.*nu);
    
    double A, B, dA;
    double metric[5];
    if (tidal_flag==true)
    {
        metric = Metric(r, params, false);
        A      = metric[0];
        B      = metric[3];
        dA     = metric[1];
    }
    else
    {
        metric = s_Metric(r, params, false); // false is added to correct the number of arguments.
        A  = metric[0];
        B  = metric[1];
        dA = metric[2];
    }
    
    double rc_vec[3];//[rc, drc, d2rc]
    //rc_vec = s_get_rc(r,aK2,params);//[rc, drc, d2rc]
    s_get_rc(rc_vec, r, params);
    double rc     = rc_vec[0];
    double drc_dr = rc_vec[1];
    double uc     = 1./rc;
    double uc2    = uc*uc;
    double uc3    = uc2*uc;
    
    double Heff_orb = sqrt(prstar2+A*(1. + pphi2*uc2 +  z3*prstar4*uc2));
    
    
    // spin variable (with dimensions)
    double S     = S1 + S2;        // => in the EMRL this becomes the spin of the BH
    double Sstar = X2*a1 + X1*a2;  // => in the EMRL this becomes the spin of the particle
    
    double ggm[14];
    s_GS(ggm, r, rc, drc_dr, aK2, prstar, pph, nu, chi1, chi2, X1, X2, c3);//nu,chi1,chi2,X1,X2);
    
    double GS              = ggm[2];
    double GSs             = ggm[3];
    double dGS_dr          = ggm[6];
    double dGSs_dr         = ggm[7];
    double d2GS_dprstar20  = ggm[12];
    double d2GSs_dprstar20 = ggm[13];
    
    
    double Heff  = Heff_orb + (GS*S + GSs*Sstar)*pph;
    double H     = sqrt( 1. + 2.*nu*(Heff - 1.) );// /nu;
    double one_H = 1./H;
    
    //==========================================================
    // Compute same quantities with prstar=0. This to obtain psi.
    // Procedure consistent with the nonspinning case
    //==========================================================
    double ggm0[14];
    s_GS(ggm0, r, rc, drc_dr, aK2, 0., pph, nu, chi1, chi2, X1, X2, c3);//nu,chi1,chi2,X1,X2);
    
    double sqrtAbyB  = sqrt(A/B);
    
    double dHeff_dr =  pph*(dGS_dr*S + dGSs_dr*Sstar) + 1./(2.*Heff_orb)*( dA*(1. + pphi2*uc2 + z3*prstar4*uc2) - 2.*A*uc3*drc_dr*(pphi2 + z3*prstar4) );
    
    double dp_rstar_dt_0 = - sqrtAbyB*one_H*dHeff_dr;
    
    // second derivative of Heff wrt to pr_star neglecting all pr_star^2 terms
    double d2Heff_dprstar20 = pph*(d2GS_dprstar20*S + d2GSs_dprstar20*Sstar) +  (1./Heff_orb)*(1. + 2.*A*uc2*z3*prstar2);
    
    double ddotr_dp_rstar = sqrtAbyB*one_H*d2Heff_dprstar20;
    
    //-------------------------------------------
    // 0.th -- approximate ddot(r)_0 without Fphi
    //-------------------------------------------
    double ddotr  =  dp_rstar_dt_0*ddotr_dp_rstar; //order pr_star^2 neglected
    
    return ddotr;

}

