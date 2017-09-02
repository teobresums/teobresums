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


#include <vector>
#include <cmath>

#include <ios>
#include <stdio.h>
#include <fstream>

#include <limits>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "TEOBResumS.h"

typedef std::numeric_limits< double > dbl;

using namespace::std;

vector<gsl_complex> s_waveform(double t, const double y[], void *params, double &Omg, double &Omg_orb, double &A, double &ddotr){
    
    double nu         = (*(TEOBResumParams *)params).nu;
    bool   tidal_flag = (*(TEOBResumParams *)params).flags.tidal;
    bool   spin_flag  = (*(TEOBResumParams *)params).flags.spin;
    double chi1       = (*(TEOBResumParams *)params).chi1;
    double chi2       = (*(TEOBResumParams *)params).chi2;
    double X1         = (*(TEOBResumParams *)params).X1;
    double X2         = (*(TEOBResumParams *)params).X2;
    double S          = (*(TEOBResumParams *)params).S;
    double Sstar      = (*(TEOBResumParams *)params).Sstar;
    double c3         = (*(TEOBResumParams *)params).cN3LO;
    double aK2        = (*(TEOBResumParams *)params).aK2;
    
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
    double dA, B, dB, one_A, one_B;
    double jhat = 0.0, Omega = 0.0, H = 0.0, Heff = 0.0, r_omega = 0.0;
    
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
    
    if (spin_flag==false)
    {
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
        
        
        //prstar evol eqn rhs
        f[2]  = (dA + ( pphi2 + z3*prstar4 )*( dA*u2 - 2.0*A*u3 ))/Heff;
        f[2] *= -0.5 * sqrAB / E;
        
        //pphi evol eqn rhs
        double sqrW = sqrt( A*(1. + pphi2*u2) );
        double psi  = 2.*(1.0 + 2.0*nu*(sqrW - 1.0))/(r2*dA);
        
        r_omega      = r*cbrt(psi);
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
        
    }
    else
    {
        
        double z3 = 2.*nu*(4.-3.*nu);
        
        double /*A,*/ B, dA;
        vector<double> metric;
        if (tidal_flag==true) {
            metric = Metric(r, params,false);
            A  = metric[0];
            B  = metric[3];
            dA = metric[1];
        }
        else
        {
            metric = s_Metric(r, params,false); //{A,B,dA,d2A} data[0]=A; data[1]=A_dr; data[2]=A_du; data[3]=B; data[4]=B_dr;
            A      = metric[0];
            B      = metric[1];
            dA     = metric[2];
        }
        
        vector<double> rc_vec;
        rc_vec        = s_get_rc(r,params); //[rc, drc, d2rc]
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
        
        /** Second derivative of Heff wrt to pr_star neglecting all pr_star^2 terms */
        double d2Heff_dprstar20 = pphi*(d2GS_dprstar20*S + d2GSs_dprstar20*Sstar) +  (1./Heff_orb)*(1 + 2*A*uc2*z3*prstar2);
        
        
        double ddotr_dp_rstar = sqrtAbyB*one_H*d2Heff_dprstar20;
        
        /*
         *-------------------------------------------
         * 0.th -- approximate ddot(r)_0 without Fphi
         *-------------------------------------------
         */
        ddotr = dp_rstar_dt_0*ddotr_dp_rstar;  // + dr_dt.*ddotr_dr; //order pr_star^2 neglected
        
        
        /*------------------ dr/dt ------------------*/
        f[0] = sqrtAbyB*one_H*dHeff_dprstar;
        
        /*----------------- d\phi/dt ----------------*/
        Omg_orb = one_H*pphi*A*uc2/Heff_orb;
        double dHeff_dpph = GS*S + (GSs + pphi*dGSs_dpph)*Sstar + pphi*A*uc2/Heff_orb;
        f[1] = one_H*dHeff_dpph;
        
        /*----------------- dp_{r*}/dt --------------*/
        f[2] = -sqrtAbyB*one_H*dHeff_dr;
        
        /*------------------ dp_{\phi}/dt -----------*/
        Omega = f[1];
        
        Omg   = Omega;
        
        
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
        
        double Heff_orb_0 = sqrt(A*(1.0 + pphi2*uc2));    /** Effective Hamiltonian H_0^eff*/
        double Heff_0     = Heff_orb_0 + (GS_0*S + GSs_0*Sstar)*pphi;
        double H0         = sqrt(1.0 + 2.0*nu*(Heff_0 - 1.0) );
        double one_H0     = 1./H0;
        
        double Gtilde     = GS_0*S     + GSs_0*Sstar;
        double dGtilde_dr = dGS_dr_0*S + dGSs_dr_0*Sstar;
        double duc_dr     = -uc2*drc_dr;
        double psic       = (duc_dr + dGtilde_dr*rc*sqrt(A/pphi2 + A*uc2)/A)/(-0.5*dA);
        r_omega           = pow( (1.0/sqrt( rc*rc*rc*psic)+Gtilde )*one_H0 ,-2./3.);
        
        double v_phi = r_omega*Omega;
        
        jhat         = pphi/(r_omega*v_phi);
        
        H *= 1./nu; /** Note the 1/nu */
    }
    
    vector<gsl_complex> waveform = hlm(t, phi, r, pphi, prstar, Omega, ddotr, H, Heff, jhat, r_omega, params);
    
    return waveform;
}

vector<double> deltalm(const double Hreal,const double Omega,const double nu)
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
    
    int kmax=35;
    vector<double> deltalmvec(35);
    for (int i=kmax; i--; ) {deltalmvec[i]=0.;}
    
    double num;
    double den;
    
    // Residual phases in Pade-resummed form when possible
    
    
    /** l=2 ------------------------------------------------------------------*/
    
    /** Pade(1,2) approximant */
    num        = 69020.*nu + 5992.*pi*sqrt_y;
    den        = 5992.*pi*sqrt_y + 2456.*nu*(28.+493.*nu* y);
    deltalmvec[0] = delta21LO*num/den;
    
    /** Pade(2,2) approximant */
    num        = (808920.*nu*pi*sqrt(y) + 137388.*pi2*y + 35.*nu2*(136080. + (154975. - 1359276.*nu)*y));
    den        = (808920.*nu*pi*sqrt(y) + 137388.*pi2*y + 35.*nu2*(136080. + (154975. + 40404.*nu)*y));
    deltalmvec[1] = delta22LO*num/den;
    
    
    /** l=3 ------------------------------------------------------------------*/
    
    /** Pade(1,2) approximant */
    num        = 4641.*nu + 1690.*pi*sqrt_y;
    den        = num + 18207.*nu2*y;
    deltalmvec[2] = delta31LO*num/den;
    
    /** Taylor-expanded form */
    num        = 1.  + 94770.*pi/(566279.*nu)*sqrt_y;
    den        = num + 80897.* nu/3159.*y;
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

