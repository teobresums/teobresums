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
#include "QNMHybridFitCab.h"

using namespace::std;

/** from https://arxiv.org/abs/1611.00332 */
static double JimenezFortezaRemnantSpin(input params)
{
    double nu      = params.nu;
    double xnu     = sqrt(1.0-4.0*nu);
    double X1      = params.X1;
    double X2      = params.X2;
    double chi1    = params.chi1;
    double chi2    = params.chi2;
    double Dchi    = chi1-chi2;
    double S       = (X1*X1*chi1+X2*X2*chi2)/(X1*X1+X2*X2);
    double a2      = 3.833;
    double a3      = -9.49;
    double a5      = 2.513;
    
    /** The functional form is taken from eq. (7), page 5. */
    double Lorb_spin_zero  = (1.3*a3*nu*nu*nu + 5.24*a2*nu*nu + 2.*sqrt(3)*nu)/(2.88*a5*nu + 1);
    
    /** Coeffcients taken from Table II, page 6: */
    double b1      = 1.00096;
    double b2      = 0.788;
    double b3      = 0.654;
    double b5      = 0.840;
    
    /** These values are taken from Table III, page 7: */
    double f21     = 8.774;
    double f31     = 22.83;
    double f50     = 1.8805;
    double f11     = 0.345225*f21 + 0.0321306*f31 - 3.66556*f50 + 7.5397;
    
    /** These values are taken from Table IV, page 10 */
    double f12     = 0.512;
    double f22     = -32.1;
    double f32     = -154;
    double f51     = -4.77;
    
    /** The following quantities were taken from the relation given in eq. (11), */
    /** page 7: fi3 = 64 - 64.*fi0 - 16.*fi1 - 4.*fi2; */
    double f13     = 64 - 16.*f11 - 4.*f12;
    double f23     = 64 - 16.*f21 - 4.*f22;
    double f33     = 64 - 16.*f31 - 4.*f32;
    double f53     = 64 - 64.*f50 - 16.*f51;
    
    /** this transformation is given in eq. (9), page (7) */
    double b1t     = b1*(f11*nu + f12*nu*nu + f13*nu*nu*nu);
    double b2t     = b2*(f21*nu + f22*nu*nu + f23*nu*nu*nu);
    double b3t     = b3*(f31*nu + f32*nu*nu + f33*nu*nu*nu);
    double b5t     = b5*(f50 + f51*nu + f53*nu*nu*nu);
    
    /** The functional form is taken from eq. (8), page 6. */
    double Lorb_eq_spin  = (0.00954*b3t*S*S*S + 0.0851*b2t*S*S - 0.194*b1t*S)/(1 - 0.579*b5t*S);
    
    /** These values are taken from Table IV, page 10: */
    double d10     = 0.322;
    double d11     = 9.33;
    double d20     = -0.0598;
    double d30     = 2.32;
    double d31     = -3.26;
    /** The functional form is taken from eq. (19a-c), page 10.*/
    double A1      = d10*xnu*nu*nu*(d11*nu+1);
    double A2      = d20*nu*nu*nu;
    double A3      = d30*xnu*nu*nu*nu*(d31*nu+1);
    
    /** The functional form is taken from eq. (15), page 9. */
    double Lorb_uneq_mass  = A1*Dchi + A2*Dchi*Dchi + A3*S*Dchi;
    
    return X1*X1*chi1+X2*X2*chi2 + Lorb_spin_zero + Lorb_eq_spin + Lorb_uneq_mass;
}

void QNMHybridFitCab(input params, vector<double> &a1, vector<double> &a2, vector<double> &a3, vector<double> &a4, vector<double> &b1, vector<double> &b2, vector<double> &b3, vector<double> &b4)
{

    // Shorthands
    int k22 = 1;
    int k21 = 0;
    int k33 = 4;

    double nu  = params.nu;
    double nu2 = nu*nu;

    for (int i=35; i--; )
    {
        a1[i]=0.;
    }
    a2 = a1;
    a3 = a1;
    a4 = a1;
    b1 = a1;
    b2 = a1;
    b3 = a1;
    b4 = a1;

    vector<double> alpha21 = a1;
    vector<double> alpha1  = a1;
    vector<double> c3A     = a1;
    vector<double> c3phi   = a1;
    vector<double> c4phi   = a1;
    vector<double> Domg    = a1;
    vector<double> Amrg    = a1;
    vector<double> c2A     = a1;

    bool   spin_flag    = params.spin;
    double af           = JimenezFortezaRemnantSpin(params);
    double a12          = params.X1*params.chi1 - params.X2*params.chi2;
    double X12          = params.X1 - params.X2;
    double aeff         = params.aK + 1./3.*a12*X12;
    double aeff2        = aeff*aeff;
    double aeff3        = aeff2*aeff;
    double aeff_omg     = params.aK + a12*X12;
    double af2          = af*af;
    double af3          = af2*af;
    double aeff_omg2    = aeff_omg  * aeff_omg;
    double aeff_omg3    = aeff_omg2 * aeff_omg;
    double aeff_omg4    = aeff_omg2 * aeff_omg2;
    double X12_2        = X12*X12;
    double Mbh          = params.Mbh;
    
    if (spin_flag == false)
    {
        // l=2 -------------------------------------------------------------------
        
        alpha21[k22] = -0.339603474900442   * nu2 + 0.016503807342446  * nu + 0.181732224707156 ;
        alpha1[k22]  = -0.180901799858457   * nu2 + 0.021618562848427  * nu + 0.087248294522483 ;
        c3A[k22]     =  0.921101179605961   * nu  - 0.584774362431865                           ;
        c3phi[k22]   = -2.225295573520468   * nu  + 4.612106599603013                           ;
        c4phi[k22]   = -7.295803229089186   * nu  + 3.970238231730844                           ;
        Domg[k22]    =  1.133871014007213   * nu2 - 0.006827404737053  * nu + 0.114279202150220 ;
        Amrg[k22]    =  1.232248048103415   * nu2 + 0.311304775361839  * nu + 1.420821306854634 ;
        
        alpha21[k21] = -0.2820300392081096  * nu2 + 0.0102479949004937 * nu + 0.1834031745850074;
        alpha1[k21]  = -0.1315990000203947  * nu2 + 0.0104928325168123 * nu + 0.0882125138092582;
        c3A[k21]     = -0.5287168076930696  * nu  + 0.1849855679575291                          ;
        c3phi[k21]   = -4.5223371288433771  * nu  + 2.6118602299368048                          ;
        c4phi[k21]   = -3.4863108298513281  * nu  + 1.0786026182494697                          ;
        Domg[k21]    =  0.3841167347154108  * nu2 + 0.1795259779761766 * nu + 0.0916127596380461;
        Amrg[k21]    = -5.7730289498820806  * nu2 + 0.0313143958615908 * nu + 0.4315195950116911;
        
        
        // l=3 ------------------------------------------------------------------
        alpha21[k33] = -0.3729868089698167  * nu2 + 0.0204125338148693 * nu + 0.1863293206607959;
        alpha1[k33]  = -0.1877130001846801  * nu2 + 0.0150701243580132 * nu + 0.0915671940372841;
        c3A[k33]     =  1.0651201096037879  * nu  - 0.3603063407108281                          ;
        c3phi[k33]   = -2.5337475643757652  * nu  + 2.6347286453107643                          ;
        c4phi[k33]   = -7.9494744361321930  * nu  + 2.8106916948840492                          ;
        Domg[k33]    =  2.5797331166178403  * nu2 - 0.5337830158170729 * nu + 0.1930766531716946;
        Amrg[k33]    = -10.4024985230145610 * nu2 + 1.3517710770250695 * nu + 0.4307642913724235;

    }
    else
    {
        double omega1_c    = -0.0598837831 * af3 + 0.8082136788 * af2 - 1.7408467418 * af + 1;
        double omega1_d    = -0.2358960279 * af3 + 1.3152369374 * af2 - 2.0764065380 * af + 1;
        double omega1      =  0.3736716844 * (omega1_c/omega1_d);
    
        /** alpha1 is alpha1[k22] */
        double alpha1_c    =  0.1211263886 * af3 + 0.7015835813 * af2 - 1.8226060896 * af + 1;
        double alpha1_d    =  0.0811633377 * af3 + 0.7201166020 * af2 - 1.8002031358 * af + 1;
        alpha1[k22]      =  0.0889623157 * (alpha1_c/alpha1_d);
        
        double alpha21_c   =  0.4764196512 * af3 - 0.0593165805 * af2 - 1.4168096833 * af + 1;
        double alpha21_d   =  0.4385578151 * af3 - 0.0763529088 * af2 - 1.3595491146 * af + 1;
        alpha21[k22]     =  0.1849525596 * (alpha21_c/alpha21_d);
        
        double a_c3A 	    =  0.0169543;
        double b_c3A 	    = -0.0799343;
        double c_c3A 	    = -0.115928;
        double c3A_nu      =  0.8298678603 * nu - 0.5615838975;
        double c3A_eq      =  (c_c3A * X12 + 0.0907476903) * aeff3 + (b_c3A * X12 + 0.0227344099) * aeff2 + (a_c3A * X12 - 0.1994944332)*aeff;
        c3A[k22]         =  c3A_nu + c3A_eq;
        
        double a_c3phi  	= -0.462321;
        double b_c3phi  	= -0.904512;
        double c_c3phi  	=  0.437747;
        double d_c3phi  	=  1.8275;
        double c3phi_nu    =  0.4558467286 * nu + 3.8883812141;
        double c3phi_equal =  (d_c3phi*X12-2.0575868122) * aeff_omg4 +(c_c3phi*X12-0.5051534498)*aeff_omg3 +(b_c3phi*X12+2.5742292762)*aeff_omg2 +(a_c3phi*X12+2.5599640181)*aeff_omg;
        c3phi[k22]       = c3phi_nu + c3phi_equal;
        
        double a_c4phi 	= -0.449976;
        double b_c4phi 	= -0.980913;
        double c4phi_nu    =  2.0822327682 * nu + 1.4996868401;
        double c4phi_equal =  (b_c4phi*X12+3.5695199109) * aeff_omg2 + (a_c4phi * X12 + 4.1312404030) * aeff_omg;
        c4phi[k22]       =  c4phi_nu + c4phi_equal;
        
        double a2_omgmx 	= -0.122735;
        double a1_omgmx 	=  0.0857478;
        double b2_omgmx 	= -0.0760023;
        double b1_omgmx 	=  0.0826514;
        double omgmx_eq_c  =  (a2_omgmx*X12_2 +a1_omgmx*X12 -0.1416002395) * aeff_omg + 1;
        double omgmx_eq_d  =  (b2_omgmx*X12_2 +b1_omgmx*X12 -0.3484804901) * aeff_omg + 1;
        double omgmx_eq    =  omgmx_eq_c/omgmx_eq_d;
        double omgmx       =  (0.481958619443355 * nu2 + 0.223976694441952 * nu + 0.273813064427363) * omgmx_eq;
        
        double a2_A_scaled = -0.0820894;
        double a1_A_scaled = 0.176126;
        double b2_A_scaled = -0.150239;
        double b1_A_scaled = 0.20491;
        double A_scaled_eq = ((a2_A_scaled*X12*X12 + a1_A_scaled*X12 -0.2935238329)*aeff + 1)/((b2_A_scaled*X12*X12 + b1_A_scaled*X12 -0.4728707630)*aeff + 1);
        double A_scaled    = (+1.826573640739664*nu2 +0.100709438291872*nu +1.438424467327531)*A_scaled_eq;
        
        Amrg[k22]      = A_scaled*(1-0.5*omgmx*aeff);
        Domg[k22]      = omega1 - Mbh*omgmx;
    }
    for (int i=35; i--; )
    {
        c2A[i] = 0.5*alpha21[i];
    }
    
    for (int i=35; i--;)
    {
        double cosh_c3A = cosh(c3A[i]);
        a1[i] = Amrg[i] * alpha1[i] * cosh_c3A * cosh_c3A / c2A[i];
        a2[i] = c2A[i];
        a3[i] = c3A[i];
        a4[i] = Amrg[i] - a1[i] * tanh(c3A[i]);
        
        b1[i] = Domg[i] * (1+c3phi[i]+c4phi[i]) / (b2[i]*(c3phi[i] + 2.*c4phi[i]));
        b2[i] = alpha21[i];
        b3[i] = c3phi[i];
        b4[i] = c4phi[i];
    }
}

