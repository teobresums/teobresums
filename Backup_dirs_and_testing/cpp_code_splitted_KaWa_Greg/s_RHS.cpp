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


#include <ios>
#include <cmath>
#include <limits>
#include <vector>
#include <fstream>
#include <stdio.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>

#include "input_struc.h"
#include "Metric.h"
#include "multipole_index.h"
#include "s_A5PNlog.h"
#include "s_Flux.h"
#include "s_GS.h"
#include "s_Metric.h"
#include "s_RHS.h"
#include "s_get_rc.h"
typedef std::numeric_limits< double > dbl;

using namespace::std;

int s_RHS(double t, const double y[], double f[], void *params)
{
/*
% This function provides the rhs of the full EOB equations
% written using p*, i.e. the conjugate momentum to  the r*
% generalized tortoise coordinate.
*/
    
/* y content
%---------------------------
%      y(0) = r
%      y(1) = phi
%      y(2) = pr*
%      y(3) = pphi
%---------------------------
*/
    
    (void)(t); /* avoid unused parameter warning */

    double nu       = (*(input *)params).nu;
    double S        = (*(input *)params).S;
    double Sstar    = (*(input *)params).Sstar;
    double chi1     = (*(input *)params).chi1;
    double chi2     = (*(input *)params).chi2;
    double X1       = (*(input *)params).X1;
    double X2       = (*(input *)params).X2;
    double c3       = (*(input *)params).cN3LO;
    bool tidal_flag = (*(input *)params).tidal;
    double aK2      = (*(input *)params).aK2;
 
    
    double r      = y[0];
    double prstar = y[2];
    double pph    = y[3];
    double pphi2  = pph*pph;

    // shorthands
    double prstar2 = prstar*prstar;
    double prstar3 = prstar2*prstar;
    double prstar4 = prstar3*prstar;

    double z3 = 2.*nu*(4.-3.*nu);
    
    double A, B, dA;
    vector<double> metric;
    if (tidal_flag==true)
    {
        metric = Metric(r, params, false);
        A      = metric[0];
        B      = metric[3];
        dA     = metric[1];
    }
    else
    {
      metric = s_Metric(r, params, false); //{A,B,dA,d2A} data[0]=A; data[1]=A_dr; data[2]=A_du; data[3]=B; data[4]=B_dr;
        A    = metric[0];
        B    = metric[1];
        dA   = metric[2];
    }
    
    vector<double> rc_vec;
    rc_vec = s_get_rc(r,params); //[rc, drc, d2rc]
    double rc     = rc_vec[0];
    double drc_dr = rc_vec[1];
    double uc     = 1./rc;
    double uc2    = uc*uc;
    double uc3    = uc2*uc;
    
    double Heff_orb = sqrt( prstar2+A*(1. + pphi2*uc2 +  z3*prstar4*uc2) );
    
    vector<double> ggm = s_GS(r, rc, drc_dr, aK2, prstar, pph, nu, chi1, chi2, X1, X2, c3);

    double GS              = ggm[2];
    double GSs             = ggm[3];
    double dGS_dprstar     = ggm[4];
    double dGSs_dprstar    = ggm[5];
    double dGS_dr          = ggm[6];
    double dGSs_dr         = ggm[7];
    double dGSs_dpph       = ggm[9];
    double d2GS_dprstar20  = ggm[12];
    double d2GSs_dprstar20 = ggm[13];
    
    
    double Heff  = Heff_orb + (GS*S + GSs*Sstar)*pph;
    double H     = sqrt( 1. + 2.*nu*(Heff - 1.) );
    double one_H = 1./H;

    double sqrtAbyB = sqrt(A/B);
    
    double dHeff_dr = pph*(dGS_dr*S + dGSs_dr*Sstar) + 1./(2.*Heff_orb)*( dA*(1. + pphi2*uc2 + z3*prstar4*uc2) - 2.*A*uc3*drc_dr*(pphi2 + z3*prstar4) );
    
    double dp_rstar_dt_0 = - sqrtAbyB*one_H*dHeff_dr;

    double dHeff_dprstar = pph*(dGS_dprstar*S + dGSs_dprstar*Sstar) + (prstar/Heff_orb)*(1. + 2.*A*uc2*z3*prstar2);
    
    // second derivative of Heff wrt to pr_star neglecting all pr_star^2 terms
    double d2Heff_dprstar20 = pph*(d2GS_dprstar20*S + d2GSs_dprstar20*Sstar) +  (1./Heff_orb)*(1. + 2.*A*uc2*z3*prstar2);

    double ddotr_dp_rstar = sqrtAbyB*one_H*d2Heff_dprstar20;
    
    //-------------------------------------------
    // 0.th -- approximate ddot(r)_0 without Fphi
    //-------------------------------------------
    double ddotr = dp_rstar_dt_0*ddotr_dp_rstar; //order pr_star^2 neglected

    
    //------------------ dr/dt ------------------
    f[0] = sqrtAbyB*one_H*dHeff_dprstar;

    //----------------- d\phi/dt ----------------
    double dHeff_dpph = GS*S + (GSs + pph*dGSs_dpph)*Sstar + pph*A*uc2/Heff_orb;
    f[1] = one_H*dHeff_dpph;

    //----------------- dp_{r*}/dt --------------
    f[2] = -sqrtAbyB*one_H*dHeff_dr;

    //------------------ dp_{\phi}/dt -----------
    double Omg = f[1];


    //----------------------------------
    // Compute here the new r_omg radius
    //----------------------------------
    //==========================================================
    // Compute same quantities with prstar=0. This to obtain psi.
    // Procedure consistent with the nonspinning case
    //==========================================================
    vector<double> ggm0 = s_GS(r, rc, drc_dr, aK2, 0., pph, nu, chi1, chi2, X1, X2, c3);//nu,chi1,chi2,X1,X2);
    
    double GS_0      = ggm0[2];
    double GSs_0     = ggm0[3];
    double dGS_dr_0  = ggm0[6];
    double dGSs_dr_0 = ggm0[7];
    double Heff_orb_0 = sqrt(A*(1.0 + pphi2*uc2));   // effective Hamiltonian H_0^eff
    double Heff_0     = Heff_orb_0 + (GS_0*S + GSs_0*Sstar)*pph;
    double H0         = sqrt(1.0 + 2.0*nu*(Heff_0 - 1.0) );
    double one_H0     = 1./H0;
    double Gtilde     = GS_0*S     + GSs_0*Sstar;
    double dGtilde_dr = dGS_dr_0*S + dGSs_dr_0*Sstar;
    double duc_dr     = -uc2*drc_dr;
    double psic       = (duc_dr + dGtilde_dr*rc*sqrt(A/pphi2 + A*uc2)/A)/(-0.5*dA);
    double r_omg      = pow( (pow( gsl_pow_int(rc,3)*psic,-1./2.)+Gtilde )*one_H0 ,-2./3.);
    double v_phi = r_omg*Omg;
    double x     = v_phi*v_phi;
    double jhat  = pph/(r_omg*v_phi);
    double Fphi  = s_Flux(x,Omg,r_omg,H,Heff,jhat,r,prstar,ddotr,params);
    f[3] = Fphi;

    return GSL_SUCCESS;
}

