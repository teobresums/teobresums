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
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <iostream>


#include "TEOBResumS.h"

using namespace::std;

/** Final spin fitting formula */
/** from https://arxiv.org/abs/1611.00332 */

double JimenezFortezaRemnantSpin(TEOBResumParams params)
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

/** Fits for the 22 mode for spinning systems */
void QNMHybridFitCab(TEOBResumParams params, vector<double> &a1, vector<double> &a2, vector<double> &a3, vector<double> &a4, vector<double> &b1, vector<double> &b2, vector<double> &b3, vector<double> &b4, vector<gsl_complex> &sigma)
{
    
    // Shorthands
    int k22 = 1;
    int k21 = 0;
    int k33 = 4;
    int k44 = 8;
    
    double nu  = params.nu;
    double nu2 = nu*nu;
    double nu3 = nu2*nu;
    
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
    vector<double> omega1  = a1;
    vector<double> c3A     = a1;
    vector<double> c3phi   = a1;
    vector<double> c4phi   = a1;
    vector<double> Domg    = a1;
    vector<double> Amrg    = a1;
    vector<double> c2A     = a1;
    
    int   spin_flag     = params.flags.spin;
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
    
    if (spin_flag == 0)
    {
        /*=========================================================*/
        /* Updated fits: 05/09/2017 from CoM extrapolated SXS data */
        /*=========================================================*/
      
        // l=2 -------------------------------------------------------------------
        
        /* (l=2, m=2)*/
	alpha21[k22] = -0.3025985041156393 *nu2 +  0.0032794155172817 *nu +  0.1828276903682022;
	alpha1[k22]  = -0.1615300454109702 *nu2 +  0.0147030662812516 *nu +  0.0878204175700328;
	c3A[k22]     =  0.8118901739129283 *nu  -  0.5584875090785957;
	c3phi[k22]   =  0.7156419884962878 *nu  +  3.8436474282409803;
	c4phi[k22]   =  2.2336960710670901 *nu  +  1.4736119175780844;
	Domg[k22]    =  0.8846304360111242 *nu2 +  0.0872792137250448 *nu +  0.1058414813686749;
	Amrg[k22]     = 1.4935750287318139 *nu2 +  0.2157497669089671 *nu +  1.4292027468283439;

	/* (l=2, m=1)*/
	alpha21[k21] = -0.2741607253846813 *nu2 +  0.0079342900879431 *nu +  0.1835522430667348;
	alpha1[k21]  = -0.1277546304610336 *nu2 +  0.0093615534859368 *nu +  0.0882855170502398;
	c3A[k21]     = -0.9431151070942140 *nu  +  0.2569989171628133;
	c3phi[k21]   = -3.4479482376671666 *nu  +  2.4755856452648359;
	c4phi[k21]   = -3.4024504071619841 *nu  +  1.0650118588151427;
	Domg[k21]    =  0.2660644668923829 *nu2 +  0.2276854484140649 *nu +  0.0884880283627388;
	Amrg[k21]    = -5.7236432632743952 *nu2 +  0.0390010969627653 *nu +  0.4291847351869338;
	        
        // l=3 ------------------------------------------------------------------
	/* (l=3,m=3)*/
        alpha21[k33] = -0.3620553934265325 *nu2 +  0.0171973908686402 *nu +  0.1865364041200878;
	alpha1[k33]  = -0.1821867653548689 *nu2 +  0.0134440240947561 *nu +  0.0916720214797975;
	c3A[k33]     =  2.7565431398030675 *nu  -  0.5506682334306747;
	c3phi[k33]   = -0.2497526471104979 *nu  +  2.3737675006958683;
	c4phi[k33]   = -2.9538823110315420 *nu  +  1.4483501341373066;
	Domg[k33]    =  1.3341439550896721 *nu2 -  0.1717105341058959 *nu +  0.1694617455660599;
	Amrg[k33]    = -9.3034388918614841 *nu2 +  1.0189351143222705 *nu +  0.4533252110436300;

	// l=4 ------------------------------------------------------------------
	/* (l=4,m=4)*/
	alpha21[k44] = -0.3991680748908423 *nu2 +   0.0287698202159666 *nu +  0.1880112530796091;
	alpha1[k44]  = -0.2003781755488581 *nu2 +   0.0171888841352427 *nu +  0.0930836242032652;
	c3A[k44]     =  3.1899853343683140 *nu  +  -0.4131730594856833;
	c3phi[k44]   = 31.5753575286023747 *nu  +  -1.0375600524681363;
	c4phi[k44]   = 25.4170586178559716 *nu  +  -0.4151371540505313;
	Domg[k44]    = -1.5342842283421341 *nu2 +   1.5224173843877831 *nu +  0.0897013049238634;
	Amrg[k44]    =  0.9438333992719329 *nu2 +  -1.0464153920266663 *nu +  0.2897769169572948;

	        
        for (int i=35; i--; )
        {
            switch (i)
            {
                case 0:
                    sigma[0].dat[0] = -0.208936*nu3 - 0.028103*nu2 - 0.005383*nu + 0.08896;
                    sigma[0].dat[1] =  0.733477*nu3 + 0.188359*nu2 + 0.220659*nu + 0.37367;
                    break;
                case 1:
                    sigma[1].dat[0] = -0.364177*nu3 + 0.010951*nu2 - 0.010591*nu + 0.08896;
                    sigma[1].dat[1] =  2.392808*nu3 + 0.051309*nu2 + 0.449425*nu + 0.37365;
                    break;
                case 4:
                    sigma[4].dat[0] = -0.319703*nu3 - 0.030076*nu2-0.009034*nu + 0.09270;
                    sigma[4].dat[1] =  2.957425*nu3 + 0.178146*nu2 + 0.709560*nu + 0.59944;
                    break;
                default:
                    sigma[i].dat[0] = 0.;
                    sigma[i].dat[1] = 0.;
                    break;
            }
        }

    }
    else
    {

      /* Setting up coefficients from the phenomenological description of the ringdown.
         For notation: Damour&Nagar, PRD 90 (2015), 024054 and Del Pozzo & Nagar, PRD 95 (2017), 124034
         Current global fits are new. See Nagar+ 2017 (in preparation) for a global performance
         and Riemenschneider& Nagar (2017) in preparation for the description of the fits */
       
        /* omg1 - imaginary part of the fundamental mode */
        double omega1_c    = -0.0598837831 * af3 + 0.8082136788 * af2 - 1.7408467418 * af + 1;
        double omega1_d    = -0.2358960279 * af3 + 1.3152369374 * af2 - 2.0764065380 * af + 1;
        omega1[k22]        =  0.3736716844 * (omega1_c/omega1_d);
        
        /* alpha1 - real part (damping time) of the fundamental mode */
        double alpha1_c    =  0.1211263886 * af3 + 0.7015835813 * af2 - 1.8226060896 * af + 1;
        double alpha1_d    =  0.0811633377 * af3 + 0.7201166020 * af2 - 1.8002031358 * af + 1;
        alpha1[k22]        =  0.0889623157 * (alpha1_c/alpha1_d);

	/* alpha2 - alpha1 */
        double alpha21_c   =  0.4764196512 * af3 - 0.0593165805 * af2 - 1.4168096833 * af + 1;
        double alpha21_d   =  0.4385578151 * af3 - 0.0763529088 * af2 - 1.3595491146 * af + 1;
        alpha21[k22]       =  0.1849525596 * (alpha21_c/alpha21_d);

	/* c3A */
        double a_c3A 	    =  0.0169543;
        double b_c3A 	    = -0.0799343;
        double c_c3A 	    = -0.115928;
        double c3A_nu       =  0.8298678603 * nu - 0.5615838975;
        double c3A_eq       =  (c_c3A * X12 + 0.0907476903) * aeff3 + (b_c3A * X12 + 0.0227344099) * aeff2 + (a_c3A * X12 - 0.1994944332)*aeff;
        c3A[k22]            =  c3A_nu + c3A_eq;

	/* c3_phi */
        double a_c3phi      = -0.462321;
        double b_c3phi      = -0.904512;
        double c_c3phi      =  0.437747;
        double d_c3phi      =  1.8275;
        double c3phi_nu     =  0.4558467286 * nu + 3.8883812141;
        double c3phi_equal  =  (d_c3phi*X12-2.0575868122) * aeff_omg4 +(c_c3phi*X12-0.5051534498)*aeff_omg3 +(b_c3phi*X12+2.5742292762)*aeff_omg2 +(a_c3phi*X12+2.5599640181)*aeff_omg;
        c3phi[k22]          = c3phi_nu + c3phi_equal;

	/* c4_phi */
        double a_c4phi      = -0.449976;
        double b_c4phi      = -0.980913;
        double c4phi_nu     =  2.0822327682 * nu + 1.4996868401;
        double c4phi_equal  =  (b_c4phi*X12+3.5695199109) * aeff_omg2 + (a_c4phi * X12 + 4.1312404030) * aeff_omg;
        c4phi[k22]          =  c4phi_nu + c4phi_equal;

	/* omg_mrg: the "merger frequency", i.e. the frequency at the peak of |h22| */
	/* Special scaling and independent variables used for the fit. AN&GR 2017 */
        double a2_omgmx     = -0.122735;
        double a1_omgmx     =  0.0857478;
        double b2_omgmx     = -0.0760023;
        double b1_omgmx     =  0.0826514;
        double omgmx_eq_c   =  (a2_omgmx*X12_2 +a1_omgmx*X12 -0.1416002395) * aeff_omg + 1;
        double omgmx_eq_d   =  (b2_omgmx*X12_2 +b1_omgmx*X12 -0.3484804901) * aeff_omg + 1;
        double omgmx_eq     =  omgmx_eq_c/omgmx_eq_d;
        double omgmx        =  (0.481958619443355 * nu2 + 0.223976694441952 * nu + 0.273813064427363) * omgmx_eq;

	/* the peak of the h22 metric (strain) waveform.*/
	/* Special scaling and independent variables used for the fit. AN& GR 2017*/	
        double a2_A_scaled = -0.0820894;
        double a1_A_scaled = 0.176126;
        double b2_A_scaled = -0.150239;
        double b1_A_scaled = 0.20491;
        double A_scaled_eq = ((a2_A_scaled*X12*X12 + a1_A_scaled*X12 -0.2935238329)*aeff + 1)/((b2_A_scaled*X12*X12 + b1_A_scaled*X12 -0.4728707630)*aeff + 1);
        double A_scaled    = (+1.826573640739664*nu2 +0.100709438291872*nu +1.438424467327531)*A_scaled_eq;
        
        Amrg[k22]      = A_scaled*(1-0.5*omgmx*aeff);
        Domg[k22]      = omega1[k22] - Mbh*omgmx;

	/* renaming real & imaginary part of the QNM complex frequency sigma */
        sigma[k22].dat[0] = alpha1[k22];
        sigma[k22].dat[1] = omega1[k22];
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
        
        b2[i] = alpha21[i];
        b3[i] = c3phi[i];
        b4[i] = c4phi[i];
        b1[i] = Domg[i] * (1+c3phi[i]+c4phi[i]) / (b2[i]*(c3phi[i] + 2.*c4phi[i]));
    }
}

/** Ringdown waveform template */
gsl_complex ringdown_match(double x, int k, vector<double> a1, vector<double> a2, vector<double> a3, vector<double> a4, vector<double> b1, vector<double> b2, vector<double> b3, vector<double> b4, vector<gsl_complex> sigma){
    
    gsl_complex psi;
    
    double amp   =  ( a1[k] * tanh(a2[k]*x +a3[k]) + a4[k] ) ;
    double phase = -b1[k]*log((1. + b3[k]*exp(-b2[k]*x) + b4[k]*exp(-2.*b2[k]*x))/(1.+b3[k]+b4[k]));
    
    psi.dat[0] = amp * exp(-sigma[k].dat[0]*x);
    psi.dat[1] = - (phase - sigma[k].dat[1]*x); //NOTE: minus sign in front by convention
    
    return psi;
    
}

/* ringdown calculation and match to the dynamics */ 
int ringdown(TEOBResumParams params, vector<vector<double> > &t_vec, vector<double> Omega_vec, vector<vector<double> > &hlm_rad, vector<vector<double> > &hlm_phase){
    
    //EOBhlm Compute the multipolar resummed waveform.
    //
    //   [hlm,phi,omglm,domglm,d2omglm, psilm,Alm,dpsilm] = ...
    //   EOBhlm(nu,t,phi,r,pph,prstar, Omega,E,Heff, EOBopt, EOBmet)
    //   return (complex) multipolar wave, phase, frequency, and derivatives,
    //   and the RWZ normalized (complex) wave, its amplitude and (complex)
    //   derivative.
    //
    //   WAVE = EOBhlm( .... ) return a structure with the wave
    //
    
    //NOTE: before the Omega_vec had a Mbh multiplied onto it!!!!!!!!!!!!!!!!!!!!!!
    
    
    double dt   = params.dt;
    double nu   = params.nu;
    double Mbh  = params.Mbh;
    double chi1 = params.chi1;
    
    
    long int pk_index = Omega_vec.size()-1;
    double Omega_pk   = Omega_vec[pk_index];
    long int i        = pk_index-1;
    while (Omega_vec[i] > Omega_pk)
    {
        pk_index = i;
        Omega_pk = Omega_vec[i];
        i--;
    }
    vector<gsl_complex> Omega_pk_grid(7);
    Omega_pk_grid[0].dat[0] = t_vec[1][pk_index-3];
    Omega_pk_grid[0].dat[1] = Omega_vec[pk_index-3];
    Omega_pk_grid[1].dat[0] = t_vec[1][pk_index-2];
    Omega_pk_grid[1].dat[1] = Omega_vec[pk_index-2];
    Omega_pk_grid[2].dat[0] = t_vec[1][pk_index-1];
    Omega_pk_grid[2].dat[1] = Omega_vec[pk_index-1];
    Omega_pk_grid[3].dat[0] = t_vec[1][pk_index];
    Omega_pk_grid[3].dat[1] = Omega_vec[pk_index];
    Omega_pk_grid[4].dat[0] = t_vec[1][pk_index+1];
    Omega_pk_grid[4].dat[1] = Omega_vec[pk_index+1];
    Omega_pk_grid[5].dat[0] = t_vec[1][pk_index+2];
    Omega_pk_grid[5].dat[1] = Omega_vec[pk_index+2];
    Omega_pk_grid[6].dat[0] = t_vec[1][pk_index+3];
    Omega_pk_grid[6].dat[1] = Omega_vec[pk_index+3];
    
    double tOmg_pk = 0.;
    double DeltaT_nqc = 0.;
    vector<double> tmrg(35);
    vector<double> tmatch(35);

    //compute true peak by interpolation from the grid
    tOmg_pk  = interpolate(dt, Omega_pk_grid);
    tOmg_pk *= 1./Mbh;

    //calculate tmatch
    double xnu =(1.-4.*nu);
    
    vector<double> dtmrg(2);
    int k21 = 0;
    int k22 = 1;
    int k33 = 4;
    
    
    /** time-shift needed when the largest object is highly spinning*/
    if (chi1 >= 0.8498)
    {
        /* Interpolating fit for Deltat_NQC. See Eq.(21) of arXiv:1506.08457 */
        DeltaT_nqc = dtnqc_fit(chi1,0.8498);
    }
    else
    {
        DeltaT_nqc = 1.;
    }
    
    tmrg[k22] = tOmg_pk-(DeltaT_nqc + 2)/Mbh;     //t_max(A22) => MERGER
    
    /** nonspinning case - old */
    /*tmrg[k22]  = tOmg_pk-3./Mbh; */             // t_max(A22) => MERGER
    
    
    /** only nonspinning case */
    dtmrg[0]   = 5.70364338 + 1.85804796*xnu  + 4.0332262*xnu*xnu; //k21
    dtmrg[1]   = 4.29550934 - 0.85938*xnu;                         //k33
    tmrg[k21]  = tmrg[k22] + dtmrg[0]/Mbh;     // t_max(A21) => peak of 21 mode
    tmrg[k33]  = tmrg[k22] + dtmrg[1]/Mbh;     // t_max(A33) => peak of 33 mode


    /* postmerger-ringdown matching time */
    tmatch = tmrg;
    
    int kmax = 35;
    vector<gsl_complex> sigma(kmax); //move this to another place
    
    for (int i=35; i--; ) {
        sigma[i].dat[0]=0.;
        sigma[i].dat[1]=0.;
        tmatch[i] += 2./Mbh;
    }
    
    vector<double> a1(35);
    vector<double> a2(35);
    vector<double> a3(35);
    vector<double> a4(35);
    vector<double> b1(35);
    vector<double> b2(35);
    vector<double> b3(35);
    vector<double> b4(35);
    
    QNMHybridFitCab(params,a1,a2,a3,a4,b1,b2,b3,b4,sigma);
    
    /*deleting data points up to tmatch (starting from the back)*/
    vector<long> I(35);
    vector<long> Size(35);
    for (int k = 35; k--; )
    {
        i = t_vec[k].size()-1;
        Size[k] = t_vec[k].size();
        switch (k)
        {
            case 0:
                while (t_vec[k][i]/Mbh>tmatch[k])
                {
                    t_vec[k].pop_back();
                    hlm_rad[k].pop_back();
                    hlm_phase[k].pop_back();
                    i--;
                }
                i++;
                break;
            case 1:
                while (t_vec[k][i]/Mbh>tmatch[k])
                {
                    t_vec[k].pop_back();
                    hlm_rad[k].pop_back();
                    hlm_phase[k].pop_back();
                    i--;
                }
                i++;
                break;
            case 4:
                while (t_vec[k][i]/Mbh>tmatch[k]) {
                    t_vec[k].pop_back();
                    hlm_rad[k].pop_back();
                    hlm_phase[k].pop_back();
                    i--;
                }
                i++;
                break;
            default:
                break;
        }
        I[k] = i;
    }
//    if (DEBUG)
//    {
//        char   outputr[256]   = "tmatch.dat";
//        std::FILE* match_file   = std::fopen(outputr, "w");
//        std::fprintf(match_file,"%e\t%e\n",tmatch[1]*Mbh,tmrg[1]*Mbh);
//        std::fclose(match_file);
//    }
    //Calculate deltaphi
    vector<gsl_complex> psi(35);
    vector<double> Deltaphi(35);
    for (int k=35; k--; )
    {
        double x    = t_vec[k][I[k]]/Mbh-tmrg[k];
        psi[k]      = ringdown_match(x, k, a1, a2, a3, a4, b1, b2, b3, b4, sigma);
        Deltaphi[k] = psi[k].dat[1] - hlm_phase[k][I[k]];
    }
    
    /** add 500 points of ringdown attachment */
    /** if we select the number of points in each mode to be such that the number of output points
     is the same, we do not need any interpolation, since the vectors are defined on the same time grid. */
    int Nringdown = 500;
    for (int k=35; k--; )
    {
        double t = t_vec[k][I[k]];
        int n_removed = Size[k]-I[k];

        switch (k)
        {
            case 0:
                for (int j=0; j < Nringdown+n_removed; j++)
                {
                    double x = t/Mbh-tmrg[k];
                    psi[k] = ringdown_match(x, k, a1, a2, a3, a4, b1, b2, b3, b4, sigma);
                    
                    psi[k].dat[1]  = psi[k].dat[1] - Deltaphi[k];
                    hlm_rad[k].push_back(psi[k].dat[0]);
                    hlm_phase[k].push_back(psi[k].dat[1]);
                    t_vec[k].push_back(t);
                    t += dt;
                }
                break;
            case 1:
                for (int j=0; j < Nringdown+n_removed; j++)
                {
                    double x = t/Mbh-tmrg[k];
                    psi[k] = ringdown_match(x, k, a1, a2, a3, a4, b1, b2, b3, b4, sigma);
                    
                    psi[k].dat[1]  = psi[k].dat[1] - Deltaphi[k];
                    hlm_rad[k].push_back(psi[k].dat[0]);
                    hlm_phase[k].push_back(psi[k].dat[1]);
                    t_vec[k].push_back(t);
                    t += dt;

                }
                break;
            case 4:
                for (int j=0; j < Nringdown+n_removed; j++)
                {
                    double x = t/Mbh-tmrg[k];
                    psi[k] = ringdown_match(x, k, a1, a2, a3, a4, b1, b2, b3, b4, sigma);
                    
                    psi[k].dat[1]  = psi[k].dat[1] - Deltaphi[k];
                    hlm_rad[k].push_back(psi[k].dat[0]);
                    hlm_phase[k].push_back(psi[k].dat[1]);
                    t_vec[k].push_back(t);
                    t += dt;
                }
                break;
            default:
                for (int j=0; j < Nringdown+n_removed-1; j++)
                {
                    hlm_rad[k].push_back(0.0);
                    hlm_phase[k].push_back(0.0);
                    t_vec[k].push_back(t);
                    t += dt;
                }
                break;
        }
    }
//    if (DEBUG)
//    {
//        char   outputr[256]   = "ringdown.dat";
//        std::FILE* ringfile   = std::fopen(outputr, "w");
//        int j                 = 0;
//        int N                 = hlm_rad[1].size();
//
//        for (j= N - Nringdown+(Size[1]-I[1]);j<N;j++)
//        {
//            std::fprintf(ringfile, "%f\t%e\t%e\n", (double)j*dt, hlm_rad[1][j], hlm_phase[1][j]);
//        }
//        std::fclose(ringfile);
//    }
    return 0;
}
