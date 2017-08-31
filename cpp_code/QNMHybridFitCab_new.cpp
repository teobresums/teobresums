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

void QNMHybridFitCab(double nu, vector<double> &a1, vector<double> &a2, vector<double> &a3, vector<double> &a4, vector<double> &b1, vector<double> &b2, vector<double> &b3, vector<double> &b4)
{

    // Shorthands
    int k22 = 1;
    int k21 = 0;
    int k33 = 4;

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

    
    double af, aeff_omg;
    double af2       = af*af;
    double af3       = af2*af;
    double aeff_omg2 = aeff_omg  * aeff_omg;
    double aeff_omg3 = aeff_omg2 * aeff_omg;
    double X12_2     = X12*X12;
    
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
    else
    {
        omega1_c    = -0.0598837831 * af3 + 0.8082136788 * af2 - 1.7408467418 * af + 1;
        omega1_d    = -0.2358960279 * af3 + 1.3152369374 * af2 - 2.0764065380 * af + 1;
        omega1      =  0.3736716844 * (omega1_c/omega1_d);
        
        alpha1_c    =  0.1211263886 * af3 + 0.7015835813 * af2 - 1.8226060896 * af + 1;
        alpha1_d    =  0.0811633377 * af3 + 0.7201166020 * af2 - 1.8002031358 * af + 1;
        alpha1      =  0.0889623157 * (alpha1_c/alpha1_d);
        
        alpha21_c   =  0.4764196512 * af3 - 0.0593165805 * af2 - 1.4168096833 * af + 1;
        alpha21_d   =  0.4385578151 * af3 - 0.0763529088 * af2 - 1.3595491146 * af + 1;
        alpha21     =  0.1849525596 * (alpha21_c/alpha21_d);
        
        a_c3A 	    =  0.0169543;
        b_c3A 	    = -0.0799343;
        c_c3A 	    = -0.115928;
        c3A_nu      =  0.8298678603 * nu - 0.5615838975;
        c3A_eq      =  (c_c3A * X12 + 0.0907476903) * aeff3 + (b_c3A * X12 + 0.0227344099) * aeff2 + (a_c3A * X12 - 0.1994944332)*aeff;
        c3A         =  c3A_nu + c3A_eq;
        
        a_c3phi  	= -0.462321;
        b_c3phi  	= -0.904512;
        c_c3phi  	=  0.437747;
        d_c3phi  	=  1.8275;
        c3phi_nu    =  0.4558467286 * nu + 3.8883812141;
        c3phi_equal =  (d_c3phi.*X12-2.0575868122) * aeff_omg4 +(c_c3phi.*X12-0.5051534498).*aeff_omg.^3 +(b_c3phi.*X12+2.5742292762).*aeff_omg.^2 +(a_c3phi.*X12+2.5599640181).*aeff_omg;
        c3phi       = c3phi_nu + c3phi_equal;
        
        a_c4phi 	= -0.449976;
        b_c4phi 	= -0.980913;
        c4phi_nu    =  2.0822327682 * nu + 1.4996868401;
        c4phi_equal =  (b_c4phi.*X12+3.5695199109) * aeff_omg2 + (a_c4phi * X12 + 4.1312404030) * aeff_omg;
        c4phi       =  c4phi_nu + c4phi_equal;
        
        a2_omgmx 	= -0.122735;
        a1_omgmx 	=  0.0857478;
        b2_omgmx 	= -0.0760023;
        b1_omgmx 	=  0.0826514;
        omgmx_eq_c  =  (a2_omgmx*X12_2 +a1_omgmx*X12 -0.1416002395) * aeff_omg + 1;
        omgmx_eq_d  =  (b2_omgmx*X12_2 +b1_omgmx*X12 -0.3484804901) * aeff_omg + 1;
        omgmx_eq    =  omgmx_eq_c./omgmx_eq_d;
        omgmx       =  (0.481958619443355 * nu2 + 0.223976694441952 * nu + 0.273813064427363) * omgmx_eq;
        
        a2_A_scaled = -0.0820894;
        a1_A_scaled = 0.176126;
        b2_A_scaled = -0.150239;
        b1_A_scaled = 0.20491;
        A_scaled_eq = ((a2_A_scaled.*X12.^2 + a1_A_scaled.*X12 -0.2935238329).*aeff + 1)./((b2_A_scaled.*X12.^2 + b1_A_scaled.*X12 -0.4728707630).*aeff + 1);
        A_scaled    = (+1.826573640739664.*nu.^2 +0.100709438291872.*nu +1.438424467327531).*A_scaled_eq;
        
        A22mrg      = A_scaled.*(1-0.5.*omgmx.*aeff);
        Domg        = omega1 - Mbh.*omgmx;
    }
}

