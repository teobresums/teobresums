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

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include "Metric.h"
#include <vector>
#include "cmath"
#include "flux.h"

#include "hlm.h"

#include <ios>
#include <fstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "hlm.h"
#include "interpolator.h"
#include "multipole_index.h"
#include "input_struc.h"
#include "s_GS.h"
#include "s_Metric.h"
#include "s_get_rc.h"


#include <limits>

#include "s_waveform.h"

typedef std::numeric_limits< double > dbl;

using namespace::std;

vector<gsl_complex> s_waveform(double t, const double y[], void *params, double &Omg, double &Omg_orb, double &A, double &ddotr){
    
    double nu         = (*(input *)params).nu;
    bool   tidal_flag = (*(input *)params).tidal;
    bool   spin_flag  = (*(input *)params).spin;
    double chi1       = (*(input *)params).chi1;
    double chi2       = (*(input *)params).chi2;
    double X1         = (*(input *)params).X1;
    double X2         = (*(input *)params).X2;
    double S          = (*(input *)params).S;
    double Sstar      = (*(input *)params).Sstar;
    double c3         = (*(input *)params).cN3LO;
    double aK2        = (*(input *)params).aK2;
    
    double f[] = {0.,0.,0.};
    
    const double z3      = 2.0*nu*(4.0-3.0*nu);
    const double phi     = y[1];
    const double prstar  = y[2];
    const double prstar2 = prstar*prstar;
    const double prstar3 = prstar2*prstar;
    const double prstar4 = prstar3*prstar;
    
    const double pphi  = y[3];
    const double pphi2 = pphi*pphi;
    
    const double r  = y[0];
    const double r2 = r*r;
    const double u  = 1./r;
    const double u2 = u*u;
    const double u3 = u2*u;

    
    vector<double> metric(5);
    double /*A,*/ dA, B, dB, one_A, one_B;
    double jhat, Omega, /*ddotr,*/ H, Heff, r_omega;
    
    /*
     //{H,Heff,dHeff_dr,dHeff_dprstar,dHeff_dpph,Omega,ddotr,jhat,r_omega,Omg_orb}
     vector<double> Hamiltonian_vars = Hamiltonian(y,params);
     H = Hamiltonian_vars[0];
     Heff = Hamiltonian_vars[1];
     dHeff_dr = Hamiltonian_vars[2];
     dHeff_dprstar = Hamiltonian_vars[3];
     dHeff_dpph = Hamiltonian_vars[4];
     //Omega = energy_variables[5];
     //ddotr = energy_variables[6];
     //jhat = energy_variables[7];
     //r_omega = energy_variables[8];
     */
    
    //const double j = pphi;
    if (spin_flag==false) {
        metric = Metric(r, params,false);
        A      = metric[0];
        dA     = metric[1];
        B      = metric[3];
        dB     = metric[4];
        one_A  = 1./A;
        one_B  = 1./B;
        
        Heff     = sqrt( prstar2 + A*(1. + pphi2*u2 + z3*prstar4*u2) );
        H        = sqrt( 1. + 2.*nu*(Heff - 1.) )/nu; //note the 1/nu here
        double E = H*nu;
        
        double sqrAB = sqrt(A/B);
        
        //r evol eqn rhs
        f[0]  = (prstar + 2.0*z3*A*prstar3*u2)/Heff;
        f[0] *= (sqrAB / E);
        
        //phi evol eqn rhs
        Omega  = A*pphi*u2/Heff;
        Omega *= 1./E;
        f[1]   = Omega;
        
        Omg    = Omega;
        
        //MOmg = Mbh*Omega; //passed back to main and then used to determine peak of MOmg curve
        
        //prstar evol eqn rhs
        f[2]  = (dA + ( pphi2 + z3*prstar4 )*( dA*u2 - 2.0*A*u3 ))/Heff;
        f[2] *= -0.5 * sqrAB / E;
        
        //pphi evol eqn rhs
        double sqrW = sqrt( A*(1. + pphi2*u2) );
        double psi  = 2.*(1.0 + 2.0*nu*(sqrW - 1.0))/(r2*dA);
        
        r_omega      = r*pow(psi,1.0/3.0);
        double v_phi = r_omega*Omega;
        jhat         = pphi/(r_omega*v_phi);
        
        double E2       = E*E;
        double tmpE     = 1./Heff+nu/E2;
        double denE     = E*Heff;
        double one_denE = 1./denE;
        
        double dHeff_dr      = 0.5*(dA + (pphi2 + z3*prstar4)*(dA*u2 - 2*A*u3))/Heff;
        double dHeff_dprstar = (prstar + z3*2.0*A*u2*prstar3)/Heff;
        double dprstar_dt    = f[2];
        double dr_dt         = f[0];
        double ddotr_dr      = sqrAB*( (prstar + z3*2.*A*u2*prstar3)*(0.5*(dA*one_A-dB*one_B)-dHeff_dr*tmpE)+ 2.0*z3*(dA*u2 - 2*A*u3)*prstar3)*one_denE;
        double ddotr_dprstar = sqrAB*( 1+z3*6.*A*u2*prstar2-(prstar + z3*2*A*u2*prstar3)*dHeff_dprstar*tmpE)*one_denE;

        ddotr = dprstar_dt*ddotr_dprstar + dr_dt*ddotr_dr;
        
    }else if (spin_flag==true){
        
        double z3 = 2.*nu*(4.-3.*nu);

        double /*A,*/ B, dA;
        vector<double> metric;
        if (tidal_flag==true) {
            metric = Metric(r, params,false);
            A  = metric[0];
            B  = metric[3];
            dA = metric[1];
        } else {
	          metric = s_Metric(r, params,false); //{A,B,dA,d2A} data[0]=A; data[1]=A_dr; data[2]=A_du; data[3]=B; data[4]=B_dr;
            A      = metric[0];
            B      = metric[1];
            dA     = metric[2];
        }

        vector<double> rc_vec;
        rc_vec        = s_get_rc(r,params);//nu,X1,X2,chi1,chi2); //[rc, drc, d2rc]
        double rc     = rc_vec[0];
        double drc_dr = rc_vec[1];
        double uc     = 1./rc;
        double uc2    = uc*uc;
        double uc3    = uc2*uc;
        
        double Heff_orb = sqrt( prstar2+A*(1. + pphi2*uc2 +  z3*prstar4*uc2) );
        
        vector<double> ggm = s_GS(r,rc,drc_dr,aK2,prstar,pphi,nu,chi1,chi2,X1,X2,c3);
        
        double GS              = ggm[2];
        double GSs             = ggm[3];
        double dGS_dprstar     = ggm[4];
        double dGSs_dprstar    = ggm[5];
        double dGS_dr          = ggm[6];
        double dGSs_dr         = ggm[7];
        double dGSs_dpph       = ggm[9];
        double d2GS_dprstar20  = ggm[12];
        double d2GSs_dprstar20 = ggm[13];
        
        Heff = Heff_orb + (GS*S + GSs*Sstar)*pphi;
        H    = sqrt( 1. + 2.*nu*(Heff - 1.) );// /nu;
        double one_H = 1./H;
        
        double sqrtAbyB = sqrt(A/B);
        
        double dHeff_dr = pphi*(dGS_dr*S + dGSs_dr*Sstar) + 1./(2*Heff_orb)*( dA*(1 + pphi2*uc2 + z3*prstar4*uc2) - 2*A*uc3*drc_dr*(pphi2 + z3*prstar4) );
        
        double dp_rstar_dt_0 = - sqrtAbyB*one_H*dHeff_dr;
        
        double dHeff_dprstar = pphi*(dGS_dprstar*S + dGSs_dprstar*Sstar) + (prstar/Heff_orb)*(1 + 2*A*uc2*z3*prstar2);
        
        // second derivative of Heff wrt to pr_star neglecting all pr_star^2 terms
        double d2Heff_dprstar20 = pphi*(d2GS_dprstar20*S + d2GSs_dprstar20*Sstar) +  (1./Heff_orb)*(1 + 2*A*uc2*z3*prstar2);
        
        
        double ddotr_dp_rstar = sqrtAbyB*one_H*d2Heff_dprstar20;
        
        //-------------------------------------------
        // 0.th -- approximate ddot(r)_0 without Fphi
        //-------------------------------------------
        ddotr = dp_rstar_dt_0*ddotr_dp_rstar;  // + dr_dt.*ddotr_dr; //order pr_star^2 neglected
        
        
        //------------------ dr/dt ------------------
        f[0] = sqrtAbyB*one_H*dHeff_dprstar;
        
        //----------------- d\phi/dt ----------------
        Omg_orb = one_H*pphi*A*uc2/Heff_orb;
        double dHeff_dpph = GS*S + (GSs + pphi*dGSs_dpph)*Sstar + pphi*A*uc2/Heff_orb;
        f[1] = one_H*dHeff_dpph;
        
        //----------------- dp_{r*}/dt --------------
        f[2] = -sqrtAbyB*one_H*dHeff_dr;
        
        //------------------ dp_{\phi}/dt -----------
        Omega = f[1];
        
        Omg   = Omega;
        
        //MOmg = Omega;//Mbh*Omega; //passed back to main and then used to determine peak of MOmg curve
        
        
        //----------------------------------
        // Compute here the new r_omg radius
        //----------------------------------
        //==========================================================
        // Compute same quantities with prstar=0. This to obtain psi.
        // Procedure consistent with the nonspinning case
        //==========================================================
        vector<double> ggm0 = s_GS(r,rc,drc_dr,aK2,0.,pphi,nu,chi1,chi2,X1,X2,c3);
        
        double GS_0      = ggm0[2];
        double GSs_0     = ggm0[3];
        double dGS_dr_0  = ggm0[6];
        double dGSs_dr_0 = ggm0[7];
        
        double Heff_orb_0 = sqrt(A*(1.0 + pphi2*uc2));                     // effective Hamiltonian H_0^eff
        double Heff_0     = Heff_orb_0 + (GS_0*S + GSs_0*Sstar)*pphi;
        double H0         = sqrt(1.0 + 2.0*nu*(Heff_0 - 1.0) );
        double one_H0     = 1./H0;
        
        double Gtilde     = GS_0*S     + GSs_0*Sstar;
        double dGtilde_dr = dGS_dr_0*S + dGSs_dr_0*Sstar;
        double duc_dr     = -uc2*drc_dr;
        double psic       = (duc_dr + dGtilde_dr*rc*sqrt(A/pphi2 + A*uc2)/A)/(-0.5*dA);
        r_omega           = pow( (pow( gsl_pow_int(rc,3)*psic,-1./2.)+Gtilde )*one_H0 ,-2./3.);
        
        double v_phi = r_omega*Omega;
        
        jhat         = pphi/(r_omega*v_phi);
        
        H *= 1./nu; //note the 1/nu
    }
    
    vector<gsl_complex> waveform = hlm(t, phi, r, pphi, prstar, Omega, ddotr, H, Heff, jhat, r_omega, params);
    
    return waveform;
}

