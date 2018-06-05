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

#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <cmath>
#include <iostream>


#include "TEOBResumS.h"

using namespace::std;

/** This routine computes the waveform, according to eq. (76) and subsequent references of https://arxiv.org/abs/1406.6913v1 */


vector<gsl_complex> hlmNewt(const double r,
                            const double Omega,
                            const double phi,
                            const double nu,
                            bool         tidal_flag,
                            bool         spin_flag
                            ){
    
    /******************************************************************************************
     * Computes the leading-order (Newtonian) prefactor  of the multipolar resummed waveform. *
     * Reference: Damour, Iyer & Nagar, PRD 79, 064004 (2009)                                 *
     ******************************************************************************************/
    
    
    /** Shorthands */
    double nu2   = nu*nu;
    double nu3   = nu*nu2;
    
    double vphi  = r*Omega;
    double vphi2 = vphi*vphi;
    double vphi3 = vphi*vphi2;
    double vphi4 = vphi*vphi3;
    double vphi5 = vphi*vphi4;
    double vphi6 = vphi*vphi5;
    double vphi7 = vphi*vphi6;
    double vphi8 = vphi*vphi7;
    double vphi9 = vphi*vphi8;
    
    
    /** Polynomials in nu */
    const double X12 = sqrt(1.-4.*nu);
    const double p1  = 1.;
    double p2 = 0.0;      /*(2,1), (3,3) & (3,1)*/
    double p4 = 0.0;      /*(4,3), (4,1), (5,5), (5,3), (5,1) */
    const double p3 = (3.*nu-1.);
    //double p2 = sqrt(1.-4.*nu);                   /*(2,1), (3,3) & (3,1)*/
    //double p4 = (2.*nu-1.)*sqrt(1.-4.*nu);        /*(4,3), (4,1), (5,5), (5,3), (5,1) */
    const double p5 = 1.-5.*nu+5.*nu2;
    const double p6 = (1.-4.*nu+3.*nu2)*X12;
    const double p7 = 7.*nu3 - 14.*nu2 + 7.*nu -1.;
    const double p8 = (4*nu3-10*nu2+6*nu-1)*X12;

    /* Special treatment when spin is on because of the singularity in the sqrt(1-4*nu) factor whn m=odd when nu=1/4.
       The nu-dependence is not factored out in the spin part, while in the nonspinning part it is re-introduced
       when the full waveform is computed. See discussion in Damour & Nagar, PRD 90, 044018, Sec. 4, Eq.(89).
       This is not done for multipole with l>4 because no spinning information is included in the waveform. */

    if (spin_flag==true)      
      {
	p2 = 1.0;
	p4 = 2.*nu-1;
      }
    else
      { /* Nonspinning case*/
	if (tidal_flag==false)
	  {
	    //p2 = sqrt(1.-4.*nu);
	    //p4 = (2.*nu-1)*sqrt(1.-4.*nu);

	    p2 = X12;
	    p4 = (2.*nu-1)*X12;
	    
	      }
	else
	  {
	    p2 = 1.;
	    //p4 = (2.*nu-1)*sqrt(1.-4.*nu);
	    p4 = (2.*nu-1)*X12;
	      }
      }

    /* old conditions: incorrect and inconsistent with the flux */
    /*if (spin_flag==true)
    {
	p2 = 1;
	p4 = 2.*nu-1;
    }
    	
    if (tidal_flag==true)
    {
	p2 = 1.;
	p4 = 2.*nu-1;
	}*/
    
    int kmax = 35;
    const double phix2 = 2. * phi;
    const double phix3 = 3. * phi;
    const double phix4 = 4. * phi;
    const double phix5 = 5. * phi;
    const double phix6 = 6. * phi;
    const double phix7 = 7. * phi;
    const double phix8 = 8. * phi;
    
    double M[] = {
        phi,phix2,
        phi,phix2,phix3,
        phi,phix2,phix3,phix4,
        phi,phix2,phix3,phix4,phix5,
        phi,phix2,phix3,phix4,phix5,phix6,
        phi,phix2,phix3,phix4,phix5,phix6,phix7,
        phi,phix2,phix3,phix4,phix5,phix6,phix7,phix8};
    
    vector<gsl_complex> hlmNewt(kmax);
    
    const double pv23 = p2 * vphi3;
    const double pv34 = p3 * vphi4;
    const double pv45 = p4 * vphi5;
    const double pv56 = p5 * vphi6;
    const double pv67 = p6 * vphi7;
    const double pv78 = p7 * vphi8;
    const double pv89 = p8 * vphi9;
    
    
    /** Compute hlmNewt (without phase factor) in complex Polar coords
     * l=2 ------------------------------------------------------------------ */
    
    hlmNewt[0].dat[1] = 3.*pi/2. - M[0];
    hlmNewt[0].dat[0] = 8./3.*sqrt(pi/5.)                * pv23;
    
    hlmNewt[1].dat[1] = pi - M[1];
    hlmNewt[1].dat[0] = 8.*sqrt(pi/5.)                   * p1 * vphi2;
    
    /** l=3 ------------------------------------------------------------------ */
    
    hlmNewt[2].dat[1] = 3.*pi/2. - M[2];
    hlmNewt[2].dat[0] = 1./3.*sqrt(2.*pi/35.)            * pv23;
    
    hlmNewt[3].dat[1] =  - M[3];
    hlmNewt[3].dat[0] =  8./3.*sqrt(pi/7.)                * pv34;
    
    hlmNewt[4].dat[1] = pi/2. - M[4];
    hlmNewt[4].dat[0] = 3.*sqrt(6.*pi/7.)                * pv23;
    
    /** l=4 ------------------------------------------------------------------ */
    
    hlmNewt[5].dat[1] = pi/2. - M[5];
    hlmNewt[5].dat[0] = 1./105.*sqrt(2.*pi)              * pv45;
    
    hlmNewt[6].dat[1] = - M[6];
    hlmNewt[6].dat[0] = 8./63.*sqrt(pi)                  * pv34;
    
    hlmNewt[7].dat[1] = 3.*pi/2. - M[7];
    hlmNewt[7].dat[0] = 9./5*sqrt(2*pi/7.)               * pv45;
    
    hlmNewt[8].dat[1] = pi - M[8];
    hlmNewt[8].dat[0] = 64./9.*sqrt(pi/7.)               * pv34;
    
    /** l=5 ------------------------------------------------------------------ */
    
    hlmNewt[9].dat[1] = pi/2. - M[9];
    hlmNewt[9].dat[0] = 1./180.*sqrt(pi/77.)             * pv45;
    
    hlmNewt[10].dat[1] = pi - M[10];
    hlmNewt[10].dat[0] = 16./135.*sqrt(pi/11.)           * pv56;
    
    hlmNewt[11].dat[1] = 3.*pi/2. - M[11];
    hlmNewt[11].dat[0] = 9./20.*sqrt(3*pi/22.)           * pv45;
    
    hlmNewt[12].dat[1] = - M[12];
    hlmNewt[12].dat[0] = 256./45.*sqrt(pi/33.)           * pv56;
    
    hlmNewt[13].dat[1] = pi/2. - M[13];
    hlmNewt[13].dat[0] = 125./12.*sqrt(5.*pi/66.)        * pv45;
    if (spin_flag==true)
      /* Need to correct the l=5, m=odd modes when spin is present
         because p4 is defined without the factor sqrt(1-4*nu) that 
	 is reintroduced in the calculation of the flm Fixed, 17-04-2018 */
      {
	hlmNewt[9].dat[0] = 1./180.*sqrt(pi/77.)         * sqrt(1.-4.*nu)  * pv45;
	hlmNewt[11].dat[0] = 9./20.*sqrt(3*pi/22.)       * sqrt(1.-4.*nu)  * pv45;
	hlmNewt[13].dat[0] = 125./12.*sqrt(5.*pi/66.)    * sqrt(1.-4.*nu)  * pv45;
      }

    
    /** l=6 ------------------------------------------------------------------ */
    
    hlmNewt[14].dat[1] = 3.*pi/2. - M[14];
    hlmNewt[14].dat[0] = 1./2079.*sqrt(2.*pi/65.)        * pv67;
    
    hlmNewt[15].dat[1] = pi - M[15];
    hlmNewt[15].dat[0] = 16./1485.*sqrt(pi/13.)          * pv56;
    
    hlmNewt[16].dat[1] = pi/2. - M[16];
    hlmNewt[16].dat[0] = 81./385.*sqrt(pi/13.)           * pv67;
    
    hlmNewt[17].dat[1] = - M[17];
    hlmNewt[17].dat[0] = 1024./495.*sqrt(2.*pi/195.)     * pv56;
    
    hlmNewt[18].dat[1] = 3.*pi/2. - M[18];
    hlmNewt[18].dat[0] = 625./63*sqrt(5.*pi/429.)        * pv67;
    
    hlmNewt[19].dat[1] = pi - M[19];
    hlmNewt[19].dat[0] = 432./5*sqrt(pi/715.)            * pv56;
    
    /** l=7 ------------------------------------------------------------------ */
    
    hlmNewt[20].dat[1] = 3.*pi/2. - M[20];
    hlmNewt[20].dat[0] = 1./108108.*sqrt(pi/10.)         * pv67;
    
    hlmNewt[21].dat[1] = - M[21];
    hlmNewt[21].dat[0] = 8./3003.*sqrt(pi/15.)           * pv78;
    
    hlmNewt[22].dat[1] = pi/2. - M[22];
    hlmNewt[22].dat[0] = 243./20020.*sqrt(3.*pi/10.)     * pv67;
    
    hlmNewt[23].dat[1] = pi - M[23];
    hlmNewt[23].dat[0] = 1024./1365.*sqrt(2.*pi/165.)    * pv78;
    
    hlmNewt[24].dat[1] = 3.*pi/2. - M[24];
    hlmNewt[24].dat[0] = 3125./3276.*sqrt(5.*pi/66.)     * pv67;
    
    hlmNewt[25].dat[1] = - M[25];
    hlmNewt[25].dat[0] = 648./35.*sqrt(3.*pi/715.)       * pv78;
    
    hlmNewt[26].dat[1] = pi/2. - M[26];
    hlmNewt[26].dat[0] = 16807./180.*sqrt(7.*pi/4290.)   * pv67;
    
    /** l=8 ------------------------------------------------------------------ */
    hlmNewt[27].dat[1] = 3.*pi/2. - M[27];
    hlmNewt[27].dat[0] = ( 5.54485779151375621e-7 ) 	* pv89;
    
    hlmNewt[28].dat[1] = 0.0 - M[28];
    hlmNewt[28].dat[0] = ( 0.0000763473331250837455 ) 	* pv78;
    
    hlmNewt[29].dat[1] = pi/2. - M[29];
    hlmNewt[29].dat[0] = ( 0.00353250998285463003 )  	* pv89;
    
    hlmNewt[30].dat[1] = pi - M[30];
    hlmNewt[30].dat[0] = ( 0.0204988821800401766 )	* pv78;
    
    hlmNewt[31].dat[1] = 3.*pi/2. - M[31];
    hlmNewt[31].dat[0] = ( 0.19579402814926015 ) 	* pv89;
    
    hlmNewt[32].dat[1] = 0.0 - M[32];
    hlmNewt[32].dat[0] = ( 0.584571015778149663 ) 	* pv78;
    
    hlmNewt[33].dat[1] = pi/2. - M[33];
    hlmNewt[33].dat[0] = ( 2.44207899966355693 ) * pv89;

    hlmNewt[34].dat[1] = pi - M[34];
    hlmNewt[34].dat[0] = 131072./315.*sqrt(2.*pi/17017.) * pv78;
    
    
    return hlmNewt;
}

vector<gsl_complex> hhatlmtail(const double Omega,
                               const double Hreal,
                               const double bphys,
                               const int    L[],
                               const int    M[])
{
    
  /******************************************************************
     * Computes the tail contribution to the resummed wave.         *
     * Reference: Damour, Iyer & Nagar, PRD 79, 064004 (2009) [DNV] *
     ****************************************************************/
    
    int kmax        = 35;
    const double pi = M_PI;
    double k;
    double hhatk;
    
    gsl_sf_result num_rad;
    gsl_sf_result num_phase;
    gsl_sf_result denom_rad;
    gsl_sf_result denom_phase;
    
    double ratio_rad;
    double ratio_ang;
    double tlm_rad;
    double tlm_phase;
    vector<gsl_complex> tlm(kmax);
    
    for (int i=kmax; i--;)
    {
        k     = M[i] * Omega;
        hhatk = k * Hreal;
        
        gsl_sf_lngamma_complex_e(L[i] + 1., -2.*hhatk, &num_rad, &num_phase);
        gsl_sf_lngamma_complex_e(L[i] + 1., 0., &denom_rad, &denom_phase);
        
        ratio_rad     = num_rad.val-denom_rad.val;
        ratio_ang     = num_phase.val-0.;
        
        tlm_rad       = ratio_rad + pi * hhatk;
        tlm_phase     = ratio_ang + 2.*hhatk*log(2.*k*bphys);
        
        tlm[i].dat[0] = exp(tlm_rad);
        tlm[i].dat[1] = tlm_phase;
    }
    return tlm;
}

vector<double> hlm_Tidal(double x,
                         void *params
                         ){
    
    /********************************************************
     * Calculate tidal correction to multipolar waveform.   *
     * Damour, Nagar & Villain, Phys.Rev. D85 (2012) 123007 *
     ********************************************************/
    
    int kmax  = 35;
    double x5 = gsl_pow_int(x,5);
    
    vector<double> hA(kmax);
    
    double lambdaA2 = (*(TEOBResumParams *)params).LambdaAl2;
    double lambdaB2 = (*(TEOBResumParams *)params).LambdaBl2;
    double XA       = (*(TEOBResumParams *)params).X1;
    double XB       = (*(TEOBResumParams *)params).X2;
    double khatA_2  = 3./2. * lambdaA2 * XB/XA * gsl_pow_int(XA,5);
    double khatB_2  = 3./2. * lambdaB2 * XA/XB * gsl_pow_int(XB,5);
    
    for (int i=kmax; i--; )
    {
        hA[i]=0.;
    }
    
    vector<double> hB     = hA;
    vector<double> betaA1 = hA;
    vector<double> betaB1 = hA;
    vector<double> hTidallm(kmax);
    
    /* Individual tidal multipoles. The l=2 also brings the NLO (relative 1PN)
       tidal correction, \beta_1^{22}, as given in Eq.(A15) of DNV.
       The other modes only have the LO contribution and the beta's don't exist. */

    /*(2,2) mode */
    hA[1]     = 2. * khatA_2 *(XA/XB+3.);
    hB[1]     = 2. * khatB_2 *(XB/XA+3.);
    
    betaA1[1] = (-202. + 560*XA - 340*XA*XA + 45*XA*XA*XA)/(42*(3-2*XA));
    betaB1[1] = (-202. + 560*XB - 340*XB*XB + 45*XB*XB*XB)/(42*(3-2*XB));

    /*(2,1) mode */
    hA[0]     = 3. * khatA_2 * XB * (3.-4.*XA)/XA;
    hB[0]     = 3. * khatB_2 * XA * (3.-4.*XB)/XB;
    
    /* (3,1) mode */
    hA[2] = 12. * khatA_2 * XB*XB/XA;
    hB[2] = 12. * khatB_2 * XA*XA/XB;

    /* (3,3) mode */
    hA[4] = hA[2];
    hB[4] = hB[2];
    
    /**********************************************************************
     Combining the pieces together to construct the tidal waveform. Note 
     that the signs are correct and consistent with Eqs.(A14) and (A16) and
     (A17) of DNV, PRD 85, 123007 (2012). The - sign in front of hA for 
     subdominant multipoles matches Eqs.(A16)-(A17) of DNV.
    **********************************************************************/
    
    /* (2,2) */
    hTidallm[1] = ( hA[1]*(1. + betaA1[1]*x) + hB[1]*(1. + betaB1[1]*x) )*x5;
    
    /* (2,1) */
    hTidallm[0] = ( -hA[0] + hB[0] )*x5;
        
    /* (3,3) */
    hTidallm[4] = ( -hA[4] + hB[4] )*x5;

    /* (3,1) */
    hTidallm[2] = ( -hA[2] + hB[2] )*x5;
    
        
    return hTidallm;
    
}

vector<double> f_lm(const double x,const double nu)
{
    
    /**EOBflm Compute the resummed amplitudes in the general nu-dependent case.
     *
     *   x   :: PN parameter
     *   nu  :: symmetric mass ratio. EMRL case is nu=0
     *   eps :: eps=1 switches on Iyer-Fujita terms; eps=0 switches them off.
     *
     *   Reference(s)
     *   Damour, Iyer & Nagar, PRD 79, 064004 (2009)     [theory]
     *   Fujita & Iyer, PRD 82, 044051 (2010)            [test-mass 5.5PN]
     *   Damour, Nagar & Bernuzzi, PRD 87, 084035 (2013) [complete information]
     *
     *   TODO: this routine requires optimization
     *   - precompute coefficients c(nu)
     *   - evaluate efficiently polynomials
     */
    
    const int    kmax =35; /** Length of vector needed to store the 8 quadropoles, 35=8+7+6+5+4+3+2*/
    const double eps  = 1.0;
    const double Pi2  = pi*pi;
    
    /** Shorthands */
    const double x2  = x*x;
    const double x3  = x*x2;
    const double x4  = x*x3;
    const double x5  = x*x4;
    const double nu2 = nu*nu;
    const double nu3 = nu*nu2;
    const double nu4 = nu*nu3;
    
    /** Compute EulerLogs */
    const double EulerGamma = 0.5772156649015328606065121;
    const double Log1       = 0.0;
    const double Log2       = 0.6931471805599453094172321;
    const double Log3       = 1.0986122886681098;
    const double Log4       = 1.3862943611198906;
    const double Log5       = 1.6094379124341003;
    const double Log6       = 1.791759469228055;
    const double Log7       = 1.9459101490553132;

    const double half_logx = 0.5*log(x)+EulerGamma+Log2;
    const double el1 = half_logx+Log1;//Eulerlog(x,1);
    const double el2 = half_logx+Log2;//Eulerlog(x,2);
    const double el3 = half_logx+Log3;//Eulerlog(x,3);
    const double el4 = half_logx+Log4;//Eulerlog(x,4);
    const double el5 = half_logx+Log5;//Eulerlog(x,5);
    const double el6 = half_logx+Log6;//Eulerlog(x,6);
    const double el7 = half_logx+Log7;//Eulerlog(x,7);
    
    vector<double> rholm(kmax);
    vector<double> flm(kmax);
    
    
    /** l=2 ------------------------------------------------------------------
     *  (2,1) */
    rholm[0] = 1.                                                                           +
    (-1.0535714285714286 + 0.27380952380952384 *nu                          )*x  +
    (-0.8327841553287982 - 0.7789824263038548  *nu + 0.13116496598639457*nu2)*x2 +
    (2.9192806270460925  - 1.019047619047619   *el1                         )*x3 +
    (-1.28235780892213   + 1.073639455782313   *el1                         )*x4 +
    (-3.8466571723355227 + 0.8486467106683944  *el1                     )*eps*x5 ;
    
    /** (2,2) */
    rholm[1] = 1.                                                                                                                         +
    (-1.0238095238095237 + 0.6547619047619048*nu                                                                          )*x  +
    (-1.94208238851096   - 1.5601379440665155*nu + 0.4625614134542706*nu2                                                 )*x2 +
    (12.736034731834051  - 2.902228713904598 *nu - 1.9301558466099282*nu2 + 0.2715020968103451*nu3 - 4.076190476190476*el2)*x3 +
    (-2.4172313935587004 + 4.173242630385488 *el2                                                                         )*x4 +
    (-30.14143102836864  + 7.916297736025627 *el2                                                                         )*x5 ;
    
    flm[0]   = gsl_pow_int(rholm[0], 2);
    flm[1]   = gsl_pow_int(rholm[1], 2);
    
    
    /** l=3 ------------------------------------------------------------------
     *  (3,1) */
    rholm[2] = 1. + (-0.7222222222222222 - 0.2222222222222222*nu)*x + (0.014169472502805836 - 0.9455667789001122*nu - 0.46520763187429853*nu2)*x2 + x3*(1.9098284139598072 - 0.4126984126984127*el1+ (-4.646868015386534 + (0.21354166666666666)*Pi2)*nu + 2.3020866307903347*nu2 - 0.5813492634480288*nu3) + x4*(0.5368150316615179 + 0.2980599647266314*el1) + eps*x5*(1.4497991763035063 - 0.0058477188106817735*el1);
    /** (3,2) */
    rholm[3] = 1. + (0.003703703703703704*(328. - 1115.*nu + 320.*nu2)*x)/(-1. + 3.*nu) + (6.235191420376606e-7*(-1.444528e6 + 8.050045e6*nu - 4.725605e6*nu2 - 2.033896e7*nu3 + 3.08564e6*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + x3*(6.220997955214429 - 1.6507936507936507*el2) + eps*x4*(-3.4527288879001268 + 2.005408583186361*el2);
    /** (3,3) */
    rholm[4] = 1. + (-1.1666666666666667 + 0.6666666666666666*nu)*x + (-1.6967171717171716 - 1.8797979797979798*nu + 0.45151515151515154*nu2)*x2 + x3*(14.10891386831863 - 3.7142857142857144*el3 + (-5.031429681429682 + (0.21354166666666666)*Pi2)*nu - 1.7781727531727531*nu2 + 0.25923767590434255*nu3) + x4*(-6.723375314944128 + 4.333333333333333*el3) + eps*x5*(-29.568699895427518 + 6.302092352092352*el3);
    
    flm[2] = gsl_pow_int(rholm[2], 3);
    flm[4] = gsl_pow_int(rholm[4], 3);
    flm[3] = gsl_pow_int(rholm[3], 3);
    
    /** l=4 ------------------------------------------------------------------
     *  (4,1) */
    rholm[5] = 1. + (0.001893939393939394*(602. - 1385.*nu + 288.*nu2)*x)/(-1. + 2.*nu) - 0.36778992787515513*x2 + x3*(0.6981550175535535 - 0.2266955266955267*el1) + eps*x4*(-0.7931524512893319 + 0.2584672482399755*el1);
    /** (4,2) */
    rholm[6] = 1. + (0.0007575757575757576*(1146. - 3530.*nu + 285.*nu2)*x)/(-1. + 3.*nu) - (3.1534122443213353e-9*(1.14859044e8 - 2.95834536e8*nu - 1.204388696e9*nu2 + 3.04798116e9*nu3 + 3.79526805e8*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + 4.550378418934105e-12*x3*(8.48238724511e11 - 1.9927619712e11*el2) + eps*x4*(-0.6621921297263365 + 0.787251738160829*el2);
    /** (4,3) */
    rholm[7] = 1. + (0.005681818181818182*(222. - 547.*nu + 160.*nu2)*x)/(-1. + 2.*nu) - 0.9783218202252293*x2 + eps*(x3*(8.519456157072423 - 2.0402597402597404*el3) +      x4*(-5.353216984886716 + 2.5735094451003544*el3));
    /** (4,4) */
    rholm[8] = 1. + (0.0007575757575757576*(1614. - 5870.*nu + 2625.*nu2)*x)/(-1. + 3.*nu) + (3.1534122443213353e-9*(-5.11573572e8 + 2.338945704e9*nu - 3.13857376e8*nu2 - 6.733146e9*nu3 + 1.252563795e9*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + x3*(15.108111214795123 - 3.627128427128427*el4) + eps*x4*(-8.857121657199649 + 4.434988849534304*el4);
    
    flm[5] = gsl_pow_int(rholm[5], 4);
    flm[6] = gsl_pow_int(rholm[6], 4);
    flm[7] = gsl_pow_int(rholm[7], 4);
    flm[8] = gsl_pow_int(rholm[8], 4);
    
    /** l=5 ------------------------------------------------------------------
     *  (5,1) */
    rholm[9] = 1. + (0.002564102564102564*(319. - 626.*nu + 8.*nu2)*x)/(-1. + 2.*nu) - 0.1047896120973044*x2 + eps*(x3*(0.642701885362399 - 0.14414918414918415*el1) +      x4*(-0.07651588046467575 + 0.11790664036817883*el1));
    /** (5,2) */
    rholm[10] = 1. + (0.00007326007326007326*(-15828. + 84679.*nu - 104930.*nu2 + 21980.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + eps*(- 0.4629337197600934*x2 +      x3*(2.354458371550237 - 0.5765967365967366*el2));
    /** (5,3) */
    rholm[11] = 1. + (0.002564102564102564*(375. - 850.*nu + 176.*nu2)*x)/(-1. + 2.*nu) - 0.5788010707241477*x2 + eps*(x3*(5.733973288504755 - 1.2973426573426574*el3) +      x4*(-1.9573287625526001 + 1.2474448628294783*el3));
    /** (5,4) */
    rholm[12] = 1. + (0.00007326007326007326*(-17448. + 96019.*nu - 127610.*nu2 + 33320.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + eps*(- 1.0442142414362194*x2 +   x3*(10.252052781721588 - 2.3063869463869464*el4));
    /**  (5,5) */
    rholm[13] = 1. + (0.002564102564102564*(487. - 1298.*nu + 512.*nu2)*x)/(-1. + 2.*nu) - 1.5749727622804546*x2 + eps*(x3*(15.939827047208668 - 3.6037296037296036*el5) +      x4*(-10.272578060123237 + 4.500041838503377*el5));
    
    flm[9]  = gsl_pow_int(rholm[9], 5);
    flm[10] = gsl_pow_int(rholm[10], 5);
    flm[11] = gsl_pow_int(rholm[11], 5);
    flm[12] = gsl_pow_int(rholm[12], 5);
    flm[13] = gsl_pow_int(rholm[13], 5);
    
    
    /** l=6 ------------------------------------------------------------------
     *  (6,1) */
    rholm[14] = 1. + (0.006944444444444444*(-161. + 694.*nu - 670.*nu2 + 124.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + eps*(- 0.29175486850885135*x2 + x3*(0.21653486654395454 - 0.10001110001110002*el1));
    /** (6,2) */
    rholm[15] = 1. + (0.011904761904761904*(-74. + 378.*nu - 413.*nu2 + 49.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + eps*( - 0.24797525070634313*x2 + x3*(1.7942694138754138 - 0.40004440004440006*el2));
    /** (6,3) */
    rholm[16] = 1. + (0.006944444444444444*(-169. + 742.*nu - 750.*nu2 + 156.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + eps*(- 0.5605554442947213*x2 + x3*(4.002558222882566 - 0.9000999000999002*el3));
    /** 6,4)  */
    rholm[17] = 1. + (0.011904761904761904*(-86. + 462.*nu - 581.*nu2 + 133.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + eps*(- 0.7228451986855349*x2 + x3*(7.359388663371044 - 1.6001776001776002*el4));
    /** 6,5)  */
    rholm[18] = 1. + (0.006944444444444444*(-185. + 838.*nu - 910.*nu2 + 220.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + eps*(- 1.0973940686333457*x2 + x3*(11.623366217471297 - 2.5002775002775004*el5));
    /** (6,6) */
    rholm[19] = 1. + (0.011904761904761904*(-106. + 602.*nu - 861.*nu2 + 273.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + eps*(- 1.5543111183867486*x2 + x3*(16.645950799433503 - 3.6003996003996006*el6));
    
    flm[14] = gsl_pow_int(rholm[14], 6);
    flm[15] = gsl_pow_int(rholm[15], 6);
    flm[16] = gsl_pow_int(rholm[16], 6);
    flm[17] = gsl_pow_int(rholm[17], 6);
    flm[18] = gsl_pow_int(rholm[18], 6);
    flm[19] = gsl_pow_int(rholm[19], 6);
    
    /** l=7 ------------------------------------------------------------------
     *  (7,1) */
    rholm[20] = 1. + (0.0014005602240896359*(-618. + 2518.*nu - 2083.*nu2 + 228.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + eps*( - 0.1508235111143767*x2 + x3*(0.2581280702019663 - 0.07355557607658449*el1));
    /** (7,2) */
    rholm[21] = 1. + (0.00006669334400426837*(16832. - 123489.*nu + 273924.*nu2 - 190239.*nu3 + 32760.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 0.351319484450667*eps*x2;
    /** (7,3) */
    rholm[22] = 1. + (0.0014005602240896359*(-666. + 2806.*nu - 2563.*nu2 + 420.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + eps*(- 0.37187416047628863*x2 + x3*(3.0835293524055283 - 0.6620001846892604*el3));
    /** (7,4) */
    rholm[23] = 1. + (0.00006669334400426837*(17756. - 131805.*nu + 298872.*nu2 - 217959.*nu3 + 41076.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 0.6473746896670599*eps*x2;
    /** (7,5) */
    rholm[24] = 1. + (0.0014005602240896359*(-762. + 3382.*nu - 3523.*nu2 + 804.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + eps*(- 0.8269193364414116*x2 + x3*(8.750589067052443 - 1.838889401914612*el5));
    /** (7,6) */
    rholm[25] = 1. + (0.0006002400960384153*(2144. - 16185.*nu + 37828.*nu2 - 29351.*nu3 + 6104.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - 1.1403265020692532*eps*x2;
    /** (7,7) */
    rholm[26] = 1. + (0.0014005602240896359*(-906. + 4246.*nu - 4963.*nu2 + 1380.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + eps*(- 1.5418467934923434*x2 + x3*(17.255875091408523 - 3.6042232277526396*el7));
    
    flm[20] = gsl_pow_int(rholm[20], 7);
    flm[21] = gsl_pow_int(rholm[21], 7);
    flm[22] = gsl_pow_int(rholm[22], 7);
    flm[23] = gsl_pow_int(rholm[23], 7);
    flm[24] = gsl_pow_int(rholm[24], 7);
    flm[25] = gsl_pow_int(rholm[25], 7);
    flm[26] = gsl_pow_int(rholm[26], 7);
    
    /** l=8 ------------------------------------------------------------------
     *  (8,1) */
    rholm[27] = 1. + (0.00005482456140350877*(20022. - 126451.*nu + 236922.*nu2 - 138430.*nu3 + 21640.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - 0.26842133517043704*eps*x2;
    /** (8,2) */
    rholm[28] = 1. + (0.0003654970760233918*(2462. - 17598.*nu + 37119.*nu2 - 22845.*nu3 + 3063.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - 0.2261796441029474*eps*x2;
    /** (8,3) */
    rholm[29] = 1. + (0.00005482456140350877*(20598. - 131059.*nu + 249018.*nu2 - 149950.*nu3 + 24520.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - 0.4196774909106648*eps*x2;
    /** (8,4) */
    rholm[30] = 1. + (0.0003654970760233918*(2666. - 19434.*nu + 42627.*nu2 - 28965.*nu3 + 4899.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - 0.47652059150068155*eps*x2;
    /** (8,5) */
    rholm[31] = 1. + (0.00027412280701754384*(4350. - 28055.*nu + 54642.*nu2 - 34598.*nu3 + 6056.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3)- 0.7220789990670207*eps*x2;
    /** (8,6) */
    rholm[32] = 1. + (0.0010964912280701754*(1002. - 7498.*nu + 17269.*nu2 - 13055.*nu3 + 2653.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 0.9061610303170207*eps*x2;
    /** (8,7) */
    rholm[33] = 1. + (0.00005482456140350877*(23478. - 154099.*nu + 309498.*nu2 - 207550.*nu3 + 38920.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3)- 1.175404252991305*eps*x2;
    /** (8,8) */
    rholm[34] = 1. + (0.0003654970760233918*(3482. - 26778.*nu + 64659.*nu2 - 53445.*nu3 + 12243.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 1.5337092502821381*eps*x2;
    
    flm[27] = gsl_pow_int(rholm[27], 8);
    flm[28] = gsl_pow_int(rholm[28], 8);
    flm[29] = gsl_pow_int(rholm[29], 8);
    flm[30] = gsl_pow_int(rholm[30], 8);
    flm[31] = gsl_pow_int(rholm[31], 8);
    flm[32] = gsl_pow_int(rholm[32], 8);
    flm[33] = gsl_pow_int(rholm[33], 8);
    flm[34] = gsl_pow_int(rholm[34], 8);
    
    return flm;
}


vector<double> s_flm(double x, void *params){
    
    /* This function computes the residual amplitude corrections flm's as 
       introduced in Damour, Iyer & Nagar, PRD 79, 064004 (2008).
       The orbital part is taken at the usual 3^{+2} PN order, i.e. 3PN terms
       are integrated by the 4PN and 5PN test-particle terms, with the higher
       modes obtained by Fujita & Iyer.

     */

    /* declaring the spin-dependent terms */
    double rho22S;
    double rho32S;
    double rho44S;
    double rho42S;
    double f21S;
    double f33S;
    double f31S;
    double f43S;
    double f41S;

    /* defining parameters */
    double nu         = (*(TEOBResumParams *)params).nu;
    double X1         = (*(TEOBResumParams *)params).X1;
    double X2         = (*(TEOBResumParams *)params).X2;
    double chi1       = (*(TEOBResumParams *)params).chi1;
    double chi2       = (*(TEOBResumParams *)params).chi2;
    double a1         = (*(TEOBResumParams *)params).a1;
    double a2         = (*(TEOBResumParams *)params).a2;
    double C_Q1       = (*(TEOBResumParams *)params).C_Q1;
    double C_Q2       = (*(TEOBResumParams *)params).C_Q2;
    bool   tidal_flag = (*(TEOBResumParams *)params).flags.tidal;

    /* Switching on Fujita-Iyer terms. This is hard-coded here */
    const double eps = 1.0;
    const double Pi2  = pi*pi;

    /* Defining useful shorthands. Note that the variables called here (a1,a2)
       are what we usually cal tilde{a}_1 and tilde{a}_2 and are defined as
       a1 = X1*chi1, a2=X2*chi2 and are passed here as parameters. Special 
       combinations of these quantities are used here to write the spin-dependent
       part of the waveform in particularly compact form, so that the (spinning)
       test-particle limit is recovered just by visual inspection of the equations */
      
    double a0      = a1+a2;
    double a12     = a1-a2;
    double X12     = X1-X2;
    double a0X12   = a0*X12;
    double a12X12  = a12*X12;
          
    /* velocity variables */
    double x2 = x*x;
    double x3 = x*x2;
    double x4 = x*x3;
    double x5 = x*x4;
    double v  = sqrt(x);
    double v2 = x;
    double v3 = v*v2;
    double v4 = v3*v;
    double v5 = v4*v;
        
    // Shorthands
    const double nu2 = nu*nu;
    const double nu3 = nu*nu2;
    const double nu4 = nu*nu3;
    
    // Compute EulerLogs
    const double el1 = Eulerlog(x,1);
    const double el2 = Eulerlog(x,2);
    const double el3 = Eulerlog(x,3);
    const double el4 = Eulerlog(x,4);
    const double el5 = Eulerlog(x,5);
    const double el6 = Eulerlog(x,6);
    const double el7 = Eulerlog(x,7);
    
    /*********************
     l=m=2 multipole
    *********************/

    /***********************************************************************************************************
     spin-orbit coefficients: the Damour & Nagar pieces
     are rewritten in compact form
     double deltam  = X12;
     double chiS    = 0.5*(chi1 + chi2);
     double chiA    = 0.5*(chi1 - chi2);
     double cSO_lo    = -2./3.*(chiS*(1-nu)+chiA*deltam);
     double cSO_nlo   = (  (-34./21. + 49./18.*nu + 209./126.*nu2)*chiS + (-34./21. - 19./42.*nu)*deltam*chiA );
    ************************************************************************************************************/
    double cSO_lo    = (-0.5*a0 - a12X12/6.);
    double cSO_nlo   = (-52./63.-19./504.*nu)*a0 - (50./63.+209./504.*nu)*a12X12;

    /*---------------------------------------
       SPIN-SPIN contribution (both BBH & BNS
    -----------------------------------------*/
    
    double cSS_lo;
    if (tidal_flag==true)
    {
        cSS_lo = 0.5*(C_Q1*a1*a1 + 2.*a1*a2 + C_Q2*a2*a2);
    }
    else
    {
        cSS_lo = 0.5*a0*a0; 
    }

    /***************************************************************/
    /* rho_22^S: Eq. (80) of Damour & Nagar, PRD 90, 044018 (2014) */
    /***************************************************************/
    rho22S = cSO_lo*v3 + cSS_lo*v4 + cSO_nlo*v5 ;
    
    /************************************************************************
     l>=3, m=even: multipoles rewritten in compact and self-explanatory form
    *************************************************************************/
    rho32S = (a0-a12X12)/(3.*(1.-3.*nu))*v;
    rho44S = (-19./30.*a0 -  (1.-21.*nu)/(30.-90.*nu)*a12X12)*v3;
    rho42S = ( -1./30.*a0 - (19.-39.*nu)/(30.-90.*nu)*a12X12)*v3;
    
    /************************************************************************
     l>=2, m=odd: multipoles rewritten in compact and self-explanatory form
    *************************************************************************/
    f21S = -1.5*a12*v + ((110./21. + 79./84.*nu)*a12 - 13./84.*a0X12)*v3;
    f33S = ((-0.25 + 2.5*nu)*a12 - 1.75*a0X12)*v3;
    f31S = ((-2.25 + 6.5*nu)*a12 + 0.25*a0X12)*v3;
    f43S = (( 5. -10.*nu)*a12 - 5.*a0X12)/(-4.+8.*nu)*v;
    f41S = f43S;
    
    /**********************************************************************************************
     Old writing using precisely the expressions in the DN paper using the variables chiA and chiS
     These expressions are fully equivalent to the ones above, but are less compact and less self
     explanatory, since the spinning test-mass limit is not recognizable at sight
  
     rho32S  = -4.*nu/(3.*(3.*nu-1.))*chiS*v;
     rho44S  = -1./(15.*(1.-3.*nu))*((42.*nu2-41.*nu+10.)*chiS + (10.-39.*nu)*deltam*chiA)*v3;
     rho42S  = -1./(15.*(1.-3.*nu))*((78.*nu2-59.*nu+10.)*chiS + (10.-21.*nu)*deltam*chiA)*v3;
     f21S = -1.5*(deltam*chiS + chiA)*v+ v3*(  (61./12. + 79./84.*nu)*deltam*chiS + (61./12. + 131./84.*nu)*chiA );
     f33S = -( deltam*chiS*(2. -  5./2.*nu) + chiA*(2. - 19./2.*nu) )*v3;
     f31S = -( deltam*chiS*(2. - 13./2.*nu) + chiA*(2. - 11./2.*nu) )*v3;
     f43S = -5.*nu/(2.*(2.*nu-1.))*(deltam*chiS - chiA)*v;
     f41S =  f43S;
    *********************************************************************************************/
        
    const int kmax = 35; //length of vector needed to store all the multipoles, 35=8+7+...+2
    vector<double> rholm(kmax);
    vector<double> flm(kmax);
    
    
    // l=2 ------------------------------------------------------------------
    //(2,2)
    rholm[1] = 1. + (-1.0238095238095237 + 0.6547619047619048*nu)*x + (-1.94208238851096 - 1.5601379440665155*nu + 0.4625614134542706*nu2)*x2 + x3*(12.736034731834051 - 2.902228713904598*nu - 1.9301558466099282*nu2 + 0.2715020968103451*nu3 - 4.076190476190476*el2) + x4*(-2.4172313935587004 + 4.173242630385488*el2) + x5*(-30.14143102836864 + 7.916297736025627*el2);
    
    flm[1] = gsl_pow_int(rholm[1]+ rho22S, 2);
    
    //(2,1)
    rholm[0] = 1. + (-1.0535714285714286 + 0.27380952380952384*nu)*x + (-0.8327841553287982 - 0.7789824263038548*nu + 0.13116496598639457*nu2)*x2 + x3*(2.9192806270460925 - 1.019047619047619*el1) + x4*(-1.28235780892213 + 1.073639455782313*el1) + eps*x5*(-3.8466571723355227 + 0.8486467106683944*el1);
    
    flm[0] = gsl_pow_int(rholm[0], 2);
    flm[0] = (X12*flm[0] + f21S);
    
    
    // l=3 ------------------------------------------------------------------
    //(3,3)
    rholm[4] = 1. + (-1.1666666666666667 + 0.6666666666666666*nu)*x + (-1.6967171717171716 - 1.8797979797979798*nu + 0.45151515151515154*nu2)*x2 + x3*(14.10891386831863 - 3.7142857142857144*el3 + (-5.031429681429682 + (0.21354166666666666)*Pi2)*nu - 1.7781727531727531*nu2 + 0.25923767590434255*nu3) + x4*(-6.723375314944128 + 4.333333333333333*el3) + eps*x5*(-29.568699895427518 + 6.302092352092352*el3);
    
    flm[4] = gsl_pow_int(rholm[4], 3);
    flm[4] = (X12*flm[4] + f33S);
    
    //(3,2)
    rholm[3] = 1. + (0.003703703703703704*(328. - 1115.*nu + 320.*nu2)*x)/(-1. + 3.*nu) + (6.235191420376606e-7*(-1.444528e6 + 8.050045e6*nu - 4.725605e6*nu2 - 2.033896e7*nu3 + 3.08564e6*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + x3*(6.220997955214429 - 1.6507936507936507*el2) + eps*x4*(-3.4527288879001268 + 2.005408583186361*el2);
    
    flm[3]=gsl_pow_int(rholm[3]+ rho32S, 3);
    
    //(3,1)
    rholm[2] = 1. + (-0.7222222222222222 - 0.2222222222222222*nu)*x + (0.014169472502805836 - 0.9455667789001122*nu - 0.46520763187429853*nu2)*x2 + x3*(1.9098284139598072 - 0.4126984126984127*el1+ (-4.646868015386534 + (0.21354166666666666)*Pi2)*nu + 2.3020866307903347*nu2 - 0.5813492634480288*nu3) + x4*(0.5368150316615179 + 0.2980599647266314*el1) + eps*x5*(1.4497991763035063 - 0.0058477188106817735*el1);
    
    flm[2] = gsl_pow_int(rholm[2], 3);
    flm[2] = (X12*flm[2] + f31S);
    
    // l=4 ------------------------------------------------------------------
    //(4,4)
    rholm[8] = 1. + (0.0007575757575757576*(1614. - 5870.*nu + 2625.*nu2)*x)/(-1. + 3.*nu) + (3.1534122443213353e-9*(-5.11573572e8 + 2.338945704e9*nu - 3.13857376e8*nu2 - 6.733146e9*nu3 + 1.252563795e9*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + x3*(15.108111214795123 - 3.627128427128427*el4) + eps*x4*(-8.857121657199649 + 4.434988849534304*el4);
    
    flm[8] = gsl_pow_int(rholm[8] + rho44S, 4);
    
    //(4,3)
    rholm[7] = 1. + (0.005681818181818182*(222. - 547.*nu + 160.*nu2)*x)/(-1. + 2.*nu) - 0.9783218202252293*x2 + eps*(x3*(8.519456157072423 - 2.0402597402597404*el3) +      x4*(-5.353216984886716 + 2.5735094451003544*el3));
    
    flm[7] = gsl_pow_int(rholm[7], 4);
    flm[7] = (X12*flm[7] + f43S);
    
    //(4,2)
    rholm[6] = 1. + (0.0007575757575757576*(1146. - 3530.*nu + 285.*nu2)*x)/(-1. + 3.*nu) - (3.1534122443213353e-9*(1.14859044e8 - 2.95834536e8*nu - 1.204388696e9*nu2 + 3.04798116e9*nu3 + 3.79526805e8*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + 4.550378418934105e-12*x3*(8.48238724511e11 - 1.9927619712e11*el2) + eps*x4*(-0.6621921297263365 + 0.787251738160829*el2);
    
    flm[6]=gsl_pow_int(rholm[6] + rho42S, 4);
    
    //(4,1)
    rholm[5] = 1. + (0.001893939393939394*(602. - 1385.*nu + 288.*nu2)*x)/(-1. + 2.*nu) - 0.36778992787515513*x2 + x3*(0.6981550175535535 - 0.2266955266955267*el1) + eps*x4*(-0.7931524512893319 + 0.2584672482399755*el1);
    
    flm[5] = gsl_pow_int(rholm[5], 4);
    flm[5] = (X12*flm[5] + f41S);
    
    // l=5 ------------------------------------------------------------------
    //(5,5)
    rholm[13] = 1. + (0.002564102564102564*(487. - 1298.*nu + 512.*nu2)*x)/(-1. + 2.*nu) - 1.5749727622804546*x2 + eps*(x3*(15.939827047208668 - 3.6037296037296036*el5) +      x4*(-10.272578060123237 + 4.500041838503377*el5));
    //(5,4)
    rholm[12] = 1. + (0.00007326007326007326*(-17448. + 96019.*nu - 127610.*nu2 + 33320.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + eps*(- 1.0442142414362194*x2 +   x3*(10.252052781721588 - 2.3063869463869464*el4));
    //(5,3)
    rholm[11] = 1. + (0.002564102564102564*(375. - 850.*nu + 176.*nu2)*x)/(-1. + 2.*nu) - 0.5788010707241477*x2 + eps*(x3*(5.733973288504755 - 1.2973426573426574*el3) +      x4*(-1.9573287625526001 + 1.2474448628294783*el3));
    //(5,2)
    rholm[10] = 1. + (0.00007326007326007326*(-15828. + 84679.*nu - 104930.*nu2 + 21980.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + eps*(- 0.4629337197600934*x2 +      x3*(2.354458371550237 - 0.5765967365967366*el2));
    //(5,1)
    rholm[9] = 1. + (0.002564102564102564*(319. - 626.*nu + 8.*nu2)*x)/(-1. + 2.*nu) - 0.1047896120973044*x2 + eps*(x3*(0.642701885362399 - 0.14414918414918415*el1) +      x4*(-0.07651588046467575 + 0.11790664036817883*el1));
    
    flm[13] = gsl_pow_int(rholm[13], 5);
    flm[12] = gsl_pow_int(rholm[12], 5);
    flm[11] = gsl_pow_int(rholm[11], 5);
    flm[10] = gsl_pow_int(rholm[10], 5);
    flm[9]  = gsl_pow_int(rholm[9], 5);
    
    // l=6 ------------------------------------------------------------------
    //(6,6)
    rholm[19] = 1. + (0.011904761904761904*(-106. + 602.*nu - 861.*nu2 + 273.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + eps*(- 1.5543111183867486*x2 + x3*(16.645950799433503 - 3.6003996003996006*el6));
    //(6,5)
    rholm[18] = 1. + (0.006944444444444444*(-185. + 838.*nu - 910.*nu2 + 220.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + eps*(- 1.0973940686333457*x2 + x3*(11.623366217471297 - 2.5002775002775004*el5));
    //(6,4)
    rholm[17] = 1. + (0.011904761904761904*(-86. + 462.*nu - 581.*nu2 + 133.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + eps*(- 0.7228451986855349*x2 + x3*(7.359388663371044 - 1.6001776001776002*el4));
    //(6,3)
    rholm[16] = 1. + (0.006944444444444444*(-169. + 742.*nu - 750.*nu2 + 156.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + eps*(- 0.5605554442947213*x2 + x3*(4.002558222882566 - 0.9000999000999002*el3));
    //(6,2)
    rholm[15] = 1. + (0.011904761904761904*(-74. + 378.*nu - 413.*nu2 + 49.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + eps*( - 0.24797525070634313*x2 + x3*(1.7942694138754138 - 0.40004440004440006*el2));
    //(6,1)
    rholm[14] = 1. + (0.006944444444444444*(-161. + 694.*nu - 670.*nu2 + 124.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + eps*(- 0.29175486850885135*x2 + x3*(0.21653486654395454 - 0.10001110001110002*el1));
    
    flm[19] = gsl_pow_int(rholm[19], 6);
    flm[18] = gsl_pow_int(rholm[18], 6);
    flm[17] = gsl_pow_int(rholm[17], 6);
    flm[16] = gsl_pow_int(rholm[16], 6);
    flm[15] = gsl_pow_int(rholm[15], 6);
    flm[14] = gsl_pow_int(rholm[14], 6);
    
    // l=7 ------------------------------------------------------------------
    //(7,7)
    rholm[26] = 1. + (0.0014005602240896359*(-906. + 4246.*nu - 4963.*nu2 + 1380.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + eps*(- 1.5418467934923434*x2 + x3*(17.255875091408523 - 3.6042232277526396*el7));
    //(7,6)
    rholm[25] = 1. + (0.0006002400960384153*(2144. - 16185.*nu + 37828.*nu2 - 29351.*nu3 + 6104.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - 1.1403265020692532*eps*x2;
    //(7,5)
    rholm[24] = 1. + (0.0014005602240896359*(-762. + 3382.*nu - 3523.*nu2 + 804.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + eps*(- 0.8269193364414116*x2 + x3*(8.750589067052443 - 1.838889401914612*el5));
    //(7,4)
    rholm[23] = 1. + (0.00006669334400426837*(17756. - 131805.*nu + 298872.*nu2 - 217959.*nu3 + 41076.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 0.6473746896670599*eps*x2;
    //(7,3)
    rholm[22] = 1. + (0.0014005602240896359*(-666. + 2806.*nu - 2563.*nu2 + 420.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + eps*(- 0.37187416047628863*x2 + x3*(3.0835293524055283 - 0.6620001846892604*el3));
    //(7,2)
    rholm[21] = 1. + (0.00006669334400426837*(16832. - 123489.*nu + 273924.*nu2 - 190239.*nu3 + 32760.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 0.351319484450667*eps*x2;
    //(7,1)
    rholm[20] = 1. + (0.0014005602240896359*(-618. + 2518.*nu - 2083.*nu2 + 228.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + eps*( - 0.1508235111143767*x2 + x3*(0.2581280702019663 - 0.07355557607658449*el1));
    
    flm[26] = gsl_pow_int(rholm[26], 7);
    flm[25] = gsl_pow_int(rholm[25], 7);
    flm[24] = gsl_pow_int(rholm[24], 7);
    flm[23] = gsl_pow_int(rholm[23], 7);
    flm[22] = gsl_pow_int(rholm[22], 7);
    flm[21] = gsl_pow_int(rholm[21], 7);
    flm[20] = gsl_pow_int(rholm[20], 7);
    
    
    // l=8 ------------------------------------------------------------------
    //(8,8)
    rholm[34] = 1. + (0.0003654970760233918*(3482. - 26778.*nu + 64659.*nu2 - 53445.*nu3 + 12243.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 1.5337092502821381*eps*x2;
    //(8,7)
    rholm[33] = 1. + (0.00005482456140350877*(23478. - 154099.*nu + 309498.*nu2 - 207550.*nu3 + 38920.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3)- 1.175404252991305*eps*x2;
    //(8,6)
    rholm[32] = 1. + (0.0010964912280701754*(1002. - 7498.*nu + 17269.*nu2 - 13055.*nu3 + 2653.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 0.9061610303170207*eps*x2;
    //(8,5)
    rholm[31] = 1. + (0.00027412280701754384*(4350. - 28055.*nu + 54642.*nu2 - 34598.*nu3 + 6056.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3)- 0.7220789990670207*eps*x2;
    //(8,4)
    rholm[30] = 1. + (0.0003654970760233918*(2666. - 19434.*nu + 42627.*nu2 - 28965.*nu3 + 4899.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - 0.47652059150068155*eps*x2;
    //(8,3)
    rholm[29] = 1. + (0.00005482456140350877*(20598. - 131059.*nu + 249018.*nu2 - 149950.*nu3 + 24520.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - 0.4196774909106648*eps*x2;
    //(8,2)
    rholm[28] = 1. + (0.0003654970760233918*(2462. - 17598.*nu + 37119.*nu2 - 22845.*nu3 + 3063.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - 0.2261796441029474*eps*x2;
    //(8,1)
    rholm[27] = 1. + (0.00005482456140350877*(20022. - 126451.*nu + 236922.*nu2 - 138430.*nu3 + 21640.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - 0.26842133517043704*eps*x2;
    
    flm[34] = gsl_pow_int(rholm[34], 8);
    flm[33] = gsl_pow_int(rholm[33], 8);
    flm[32] = gsl_pow_int(rholm[32], 8);
    flm[31] = gsl_pow_int(rholm[31], 8);
    flm[30] = gsl_pow_int(rholm[30], 8);
    flm[29] = gsl_pow_int(rholm[29], 8);
    flm[28] = gsl_pow_int(rholm[28], 8);
    flm[27] = gsl_pow_int(rholm[27], 8);
    
    return flm;
    
}

                            /*******************************************************
                             * Checkpoint Charlie: You are entering the NQC sector *
                             *******************************************************/



vector<vector<gsl_complex> > find_a1a2a3(
                                         vector<double>          T,
                                         vector<double>          r,
                                         vector<double>          w,
                                         vector<double>          pph,
                                         vector<double>          pr_star,
                                         vector<vector<double>>  hlm_phase,
                                         vector<double>          Omg_orb,
                                         vector<vector<double>>  A,
                                         vector<double>          ddotr,
                                         void                    *params
                                         ){
    
    
    /**************************************************************
     * Computes the factors and the coefficients) that build the  *
     * NQC corrections to the waveform in the spinning case.      *
     **************************************************************/
    
    double A_tmp, dA_tmp, omg_tmp, domg_tmp;
    
    double nu         = (*(TEOBResumParams *)params).nu;
    double nu2        = nu*nu;
    double nu3        = nu*nu*nu;
    double X1         = (*(TEOBResumParams *)params).X1;
    double X2         = (*(TEOBResumParams *)params).X2;
    double X12        = X1 - X2;
    double chi1       = (*(TEOBResumParams *)params).chi1;
    double chi2       = (*(TEOBResumParams *)params).chi2;
    
    double aK         = (*(TEOBResumParams *)params).aK;
    double aK2        = aK*aK;
    double aK3        = aK*aK*aK;
    double aK4        = aK*aK*aK*aK;
    double a12        = X1*chi1 - X2*chi2;
    
    long int t_length = T.size();
    double c_p1,     c_p2,     c_p3,   c_p4;
    double c_pdA1,   c_pdA2,   c_pdA3, c_pdA4;
    double c_pdomg1, c_pdomg2;
    double n0, d1;
    double a0_omg_tmp, a1_omg_tmp, a2_omg_tmp, b0_omg_tmp, b1_omg_tmp, b2_omg_tmp;
    double a0_domg_tmp, a1_domg_tmp, a2_domg_tmp, b0_domg_tmp,b1_domg_tmp, b2_domg_tmp;
    double a0_A_tmp, a1_A_tmp , a2_A_tmp, b0_A_tmp,b1_A_tmp, b2_A_tmp;
    double a0_dA_tmp, a1_dA_tmp, a2_dA_tmp, b0_dA_tmp,b1_dA_tmp, b2_dA_tmp;
    double omg_tmp_nu, omg_tmp_equal, domg_tmp_nu, domg_tmp_equal, A_tmp_scale_nu, A_tmp_scale_equal, dA_tmp_scale_nu, dA_tmp_scale_equal ;
    double aeff        = aK + 1/3*a12*X12;
    double aeff_omg    = aK + a12*X12;
    vector<double> P(2);
    vector<double> M(4);
    
    vector<double> p1(2);
    vector<double> p2(2);
    vector<double> p3(2);
    vector<double> p4(2);
    vector<double> pA(5);
    vector<double> pdA(5);
    
    vector<double> pomg(5);
    vector<double> pdomg(5);
    vector<double> pn0(2),       pd1(2);
    vector<double> ppdomg1(2),   ppdomg2(2);
    vector<double> pdA1(2),      pdA2(2),      pdA3(2),      pdA4(2);
    vector<double> n1(t_length), n2(t_length), n3(t_length), n4(t_length), n5(t_length), n6(t_length);
    vector<double> max_A(35),    max_dA(35),   d2max(35),    d3max(35),    max_omg(35),  max_domg(35), maxd2omg(35), DeltaT(35);
    
    vector<vector<double> > ai(   35, vector<double>(3));
    vector<vector<double> > bi(   35, vector<double>(3));
    vector<vector<double> > omg(  35, vector<double>(t_length));
    vector<vector<double> > domg( 35, vector<double>(t_length));
    vector<vector<double> > m11(  35, vector<double>(t_length));
    vector<vector<double> > m12(  35, vector<double>(t_length));
    vector<vector<double> > m13(  35, vector<double>(t_length));
    vector<vector<double> > m21(  35, vector<double>(t_length));
    vector<vector<double> > m22(  35, vector<double>(t_length));
    vector<vector<double> > p1tmp(35, vector<double>(t_length));
    vector<vector<double> > p2tmp(35, vector<double>(t_length));
    
    vector<vector<gsl_complex>> o(35, vector<gsl_complex>(t_length));
    
    for (int k=35; k--;)
    {
        for (int j=t_length;j--;)
        {
	  /* Regge-Wheeler-Zerilli normalized amplitude. The ringdown
             coefficient refer to this normalization.
             Check Nagar & Rezzolla, CQG 22 (2005) R167 */           
	  A[k][j] = A[k][j]/sqrt( (L[k]+2)*(L[k]+1)*L[k]*(L[k]-1) );
        }
    }
    
    
    for (int k=35; k--;)
    {
        omg[k] = s_D1(hlm_phase[k], T, t_length-1);
        domg[k] = s_D1(omg[k],       T, t_length-1);
    }
    
    /**  Case 'NQC_fit_hybrid' */
    if (nu == 0.25)
    {
        pA[0]    =  0.00178195;
        pA[1]    =  0.00435589;
        pA[2]    =  0.00344489;
        pA[3]    = -0.00076165;
        pA[4]    =  0.31973334;
        A_tmp    =  pA[0]*aK4    + pA[1]*aK3   + pA[2]*aK2    + pA[3]*aK     + pA[4];
        
        pdA[0]   =  0.00000927;
        pdA[1]   = -0.00024550;
        pdA[2]   =  0.00012469;
        pdA[3]   =  0.00123845;
        pdA[4]   = -0.00195014;
        dA_tmp   =  pdA[0]*aK4   + pdA[1]*aK3   + pdA[2]*aK2   + pdA[3]*aK   + pdA[4];
        
        pomg[0]  =  0.00603482;
        pomg[1]  =  0.01604555;
        pomg[2]  =  0.02290799;
        pomg[3]  =  0.07084587;
        pomg[4]  =  0.38321834;
        omg_tmp  =  pomg[0]*aK4  + pomg[1]*aK3  + pomg[2]*aK2  + pomg[3]*aK  + pomg[4];
        
        pdomg[0] =  0.00024066;
        pdomg[1] =  0.00038123;
        pdomg[2] = -0.00049714;
        pdomg[3] =  0.00041219;
        pdomg[4] =  0.01190548;
        domg_tmp =  pdomg[0]*aK4 + pdomg[1]*aK3 + pdomg[2]*aK2 + pdomg[3]*aK + pdomg[4];
    }
    
    else if( nu > 0.16)
    {
        p1[0]      =  0.04680896;
        p1[1]      = -0.00632114;
        p2[0]      =  0.06586192;
        p2[1]      = -0.01180039;
        p3[0]      = -0.11617413;
        p3[1]      =  0.02704959;
        p4[0]      =  0.15597465;
        p4[1]      =  0.28034978;
        c_p1       =  p1[0]*nu + p1[1];
        c_p2       =  p2[0]*nu + p2[1];
        c_p3       =  p3[0]*nu + p3[1];
        c_p4       =  p4[0]*nu + p4[1];
        A_tmp      =  c_p1*aK3 + c_p2*aK2 + c_p3*aK + c_p4;
        
        pdA1[0]    = -0.00130824;
        pdA1[1]    =  0.00006202;
        pdA2[0]    =  0.00199855;
        pdA2[1]    = -0.00027474;
        pdA3[0]    =  0.00218838;
        pdA3[1]    =  0.00071540;
        pdA4[0]    = -0.00362779;
        pdA4[1]    = -0.00105397;
        c_pdA1     =  pdA1[0]*nu + pdA1[1];
        c_pdA2     =  pdA2[0]*nu + pdA2[1];
        c_pdA3     =  pdA3[0]*nu + pdA3[1];
        c_pdA4     =  pdA4[0]*nu + pdA4[1];
        dA_tmp     =  c_pdA1*aK3   + c_pdA2*aK2 + c_pdA3*aK+ c_pdA4;
        
        pn0[0]     =  0.46908067;
        pn0[1]     =  0.27022141;
        pd1[0]     =  0.64131115;
        pd1[1]     = -0.37878384;
        n0         =  pn0[0]*nu + pn0[1];
        d1         =  pd1[0]*nu + pd1[1];
        omg_tmp    =  n0/(1 + d1*aK);
        
        ppdomg1[0] =  0.00061175;
        ppdomg1[1] =  0.00074001;
        ppdomg2[0] =  0.02504442;
        ppdomg2[1] =  0.00548217;
        c_pdomg1   =  ppdomg1[0]*nu + ppdomg1[1];
        c_pdomg2   =  ppdomg2[0]*nu + ppdomg2[1];
        domg_tmp   =  c_pdomg1*aK   + c_pdomg2;
    }
    
    else
    {
      /* This fit is the most updated and recent that Gunnar did by 
         incorporating in the fits also the test-particle NQC point
	 obtained from the most-recent Teukoslky waveforms done by
	 M. Colleoni using the 6PN-accurare iResum-radiation reaction.
	 These points assure a smooth connection between merger and
	 ringdown also outside the "calibration" domain, notably for
	 large-mass ratios (though q<=20) and large (negative) spins
         Updated, 28/09/2017 */
      
	a0_omg_tmp    = -0.1460961247;
	a1_omg_tmp    =  0.0998056;
	a2_omg_tmp    = -0.118098;
	b0_omg_tmp    = -0.3430184009;
	b1_omg_tmp    =  0.0921551;
	b2_omg_tmp    = -0.0740285;
	omg_tmp_nu    = +0.5427169903*nu2 +0.2512395608*nu +0.2863992248;
	omg_tmp_equal =((a2_omg_tmp*X12*X12 + a1_omg_tmp*X12 + a0_omg_tmp)*aeff_omg+1)/((b2_omg_tmp*X12*X12 +b1_omg_tmp*X12 + b0_omg_tmp)*aeff_omg+1);
	omg_tmp       = omg_tmp_nu*omg_tmp_equal;


	a0_domg_tmp    = +0.0604556289;
	b0_domg_tmp    = -0.0299583285;
	a1_domg_tmp    = 0.0711715;
	a2_domg_tmp    = -0.0500886;
	b1_domg_tmp    = 0.0461239;
	b2_domg_tmp    = -0.0153068;
	
	domg_tmp_nu    = ( +0.0045213831*nu +0.0064934920)/( -1.4466409969*nu+1);
	domg_tmp_equal = (a2_domg_tmp*X12*X12 +a1_domg_tmp*X12 +b0_domg_tmp)*aeff_omg*aeff_omg +(b2_domg_tmp*X12*X12 +b1_domg_tmp*X12+a0_domg_tmp)*aeff_omg+1;
	domg_tmp       = domg_tmp_nu*domg_tmp_equal;

	a0_A_tmp 	= -0.2750516062;
	b0_A_tmp 	= -0.4693776065;
	a1_A_tmp 	=  0.143066;
	a2_A_tmp 	= -0.0425947;
	b1_A_tmp 	=  0.176955;
	b2_A_tmp 	= -0.111902;
	
	A_tmp_scale_nu    = -0.9862040409*nu3 +0.8167558040*nu2 -0.0427442282*nu+0.2948879452;
	A_tmp_scale_equal = ((a2_A_tmp*X12*X12 + a1_A_tmp*X12 +a0_A_tmp)*aeff+1)/((b2_A_tmp*X12*X12 + b1_A_tmp*X12 +b0_A_tmp)*aeff+1);
	A_tmp             = A_tmp_scale_nu*A_tmp_scale_equal*(1-0.5*omg_tmp*aeff);

	a0_dA_tmp 	= +0.0037461628;
	b0_dA_tmp 	= +0.0636082543;
	a1_dA_tmp 	=  0.00129393;
	a2_dA_tmp 	= -0.00239069;
	b1_dA_tmp 	= -0.0534209;
	b2_dA_tmp 	= -0.186101;
	
	dA_tmp_scale_nu    = ( -0.0847947167*nu -0.0042142765)/( +16.1559461812*nu+1);
	dA_tmp_scale_equal = ((a2_dA_tmp*X12*X12 + a1_dA_tmp*X12+ a0_dA_tmp)*aeff)/((b2_dA_tmp*X12*X12 + b1_dA_tmp*X12 + b0_dA_tmp)*aeff+1);
	dA_tmp             = (dA_tmp_scale_nu +dA_tmp_scale_equal)*omg_tmp;
        
    }
    
    for (int k=35; k--;)
    {
        max_A[k]    = 0.;
        max_dA[k]   = 0.;
        max_omg[k]  = 0.;
        max_domg[k] = 0.;
    }
    /** Switch on the 22  values (only) */
    max_A[1]    = A_tmp;
    max_dA[1]   = dA_tmp;
    max_omg[1]  = omg_tmp;
    max_domg[1] = domg_tmp;

    if (DEBUG)
    {
        printf("--------------------------------\n");
        printf("NR values for NQC determination:\n");
        printf("--------------------------------\n");
        printf("Amrg    =%10.6f\n",max_A[1]);
        printf("dAmrg   =%10.6f\n",max_dA[1]);
        printf("omg_mrg =%10.6f\n",max_omg[1]);
        printf("domg_mrg=%10.6f\n",max_domg[1]);
    }

    /** NQC corrections to AMPLITUDE (n1,n2) and PHASE (n4,n5)
     * NQC basis for (2,2) waveform : AMPLITUDE
     * note: n3 and n6 are not used
     */
    
    for (int j=t_length-1;j--;) //THIS WAS BAD  (int j=t_length;j--;)
    {
        double pr_star2 = pr_star[j] * pr_star[j];
        double r2       = r[j] * r[j];
        double w2       = w[j] * w[j];
        
        n1[j]  = pr_star2/(r2*w2);        // [pr*/(r Omg)]^2
        n2[j]  = ddotr[j]/(r[j]*w2);      // [ddot{r}/(r Omg^2)]
        
        // NQC basis for (2,2) waveform: PHASE
        
        n4[j]  = pr_star[j]/(r[j]*w[j]);           //  pr*/(r Omg)
        n5[j]  = n4[j]*r2*w2;                      // (pr*)*(r Omg)
    }

    if (DEBUG)
    {
        for (int j=t_length-1;j--;)
        {
            char   outputNQC[256]   = "NQC_func.dat";
            std::FILE* NQCfile   = std::fopen(outputNQC, "w");
            std::fprintf(NQCfile, "%20.12f\t%20.12f\t%20.12f\t%20.12f\t%20.12f\n", T[j], n1[j], n2[j], n4[j], n5[j]);
            std::fclose(NQCfile);
        }
    }
    
    /** Take the needed derivatives for the phase */
    
    vector<double>  d_n4 = s_D1(n4,T,t_length-1);
    vector<double>  d_n5 = s_D1(n5,T,t_length-1);
    vector<double> d2_n4 = s_D1(d_n4,T,t_length-1);
    vector<double> d2_n5 = s_D1(d_n5,T,t_length-1);

    
    /* A check: output derivatives*/
    if (DEBUG)
    {
        char   outputdNQC[256]   = "dNQC_func.dat";
        std::FILE* dNQCfile   = std::fopen(outputdNQC, "w");
        for (int j=t_length-1;j--;)
        {                
            std::fprintf(dNQCfile, "%20.12f\t%20.12f\t%20.12f\t%20.12f\t%20.12f\n", T[j], d_n4[j], d_n5[j], d2_n4[j], d2_n5[j]);
        }
        
        std::fclose(dNQCfile);
    }
     
    int Omgmax_index = 0;
    double Omg_max   = Omg_orb[0];
    int i            = 1;
    
    while (Omg_orb[i] > Omg_max)
    {
        Omg_max      = Omg_orb[i];
        Omgmax_index = i;
        i++;
    }
    
    double tOmgOrb_pk = T[Omgmax_index];
    double DeltaT_nqc = 0.;

    /*Additional usual time-shift for matching/NQC calculation with respect to the
      peak of the pure orbital frequency. It is 1 always, except for very large mass
      ratios and high, negative, spins */

    if ((chi1 <-0.85) && (nu <= 14./225.))
      {
	DeltaT_nqc = 4.0;
      }
    else
    {
        DeltaT_nqc = 1.0;
    }

    /****************************************************************************
     All this commented block belongs to pre May 10, 2018. It was the hack needed
     from the very beginning, 2013, to account for high-positive spin. This was
     actually related to a factor X1-X2 missing in the l=5, m=odd modes. Fixed
     this, all this is now useless, though it needs a new representation of c3
     - THIS COMMENTED BLOCK TO BE DELETED
     ****************************************************************************
    
    /* Additional time-shift needed when the largest object  is  highly spinning
       and/or for high, anti-aligned spins. Quick hack that will be removed once
       the iResum waveform will be available */
    //if (chi1 >= 0.8498)
    //  {
       
	/* Interpolating fit for Deltat_NQC. See Eq.(21) of arXiv:1506.08457 
	   This is a formula that was obtained in the equal-mass, equal-spin
	   case and promoted also to any other case where the spin on the larger
	   BH is larger than 0.8498. This is a guess to extrapolate the model
	   outside the domain of calibration */
      
    //  DeltaT_nqc = dtnqc_fit(chi1,0.8498);
    //  }
    //else if ((chi1 <=-0.80) && (nu <= 8./81.))
      /* This condition was a simple hack to avoid unphysical features in the
         modulus amplitude when one (or both) the spins are large and negative
         and the mass ratio is large. This little modification in the location
         of the NQC point guarantees that the determination of the NQC parameters
         guarantees just a small perturbation of the non-NQC EOB waveform. The
         iResum waveform will be robust enough that this hack will not be needed*/  
    //  {
    //	DeltaT_nqc = 3.0;
    //  }    
    //else
      /* standard choice inspired by test-particle results */  
    //  {
    //    DeltaT_nqc = 1.;
    //  }
    //=============================================================================

    
    double tNQC = tOmgOrb_pk - DeltaT_nqc;
    if (DEBUG)
    {
        printf("-------------------------------\n");
        printf("Check: NQC related information:\n");
        printf("-------------------------------\n");
        printf("DeltaT_tNQC = %f\n",DeltaT_nqc);
        printf("tNQC [bare] = %f\n",tNQC);
    }
    
    i        = 0;
    int jmax = 0;
    while (T[i] < tNQC)
    {
        jmax = i;
        i++;
    }
    
    /** Determination of NQC correction: solving a linear systems */
    
    for (int k=35; k--; )
    {
        for (int j=t_length; j--;)
        {
            /** Matrix elements: waveform amplitude at all points */
            m11[k][j] = n1[j]*A[k][j];
            m12[k][j] = n2[j]*A[k][j];
            
        }
	/* Attention here: one is taking FD derivatives */
        m21[k]   = s_D1(m11[k],T,t_length-1);
        m22[k]   = s_D1(m12[k],T,t_length-1);
        
        p1tmp[k] = A[k];                          // amplitude
        p2tmp[k] = s_D1(p1tmp[k],T,t_length-1);   // derivative of amplitude
    }

    /* A check: output derivatives*/
    if (DEBUG)
    {
        char   outputA[256]   = "Amp_func.dat";
        std::FILE* Afile      = std::fopen(outputA, "w");
        for (int j=t_length;j--;)
        {                
            std::fprintf(Afile, "%f\t%f\t%f\n", T[j], p1tmp[1][j], p2tmp[1][j]);
        }
        
        std::fclose(Afile);
    }
  
    // similar test for the frequencies

    
    double detM = 1.;
    for (int k=35;k--;)
    {
        /** Computation of ai coefficients */
        P[0]     = max_A[k]  - p1tmp[k][jmax];
        P[1]     = max_dA[k] - p2tmp[k][jmax];
        
        M[0]     = m11[k][jmax];
        M[1]     = m12[k][jmax];
        M[2]     = m21[k][jmax];
        M[3]     = m22[k][jmax];
        
        detM     = M[0]*M[3]-M[1]*M[2];
        ai[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
        ai[k][1] = (M[0]*P[1] - M[2]*P[0])/detM;
        
        /** Computation of bi coefficients */
        P[0]     = omg[k][jmax]   - max_omg[k];
        P[1]     = domg[k][jmax]  - max_domg[k];
        
        M[0]     = d_n4[jmax];
        M[1]     = d_n5[jmax];
        M[2]     = d2_n4[jmax];
        M[3]     = d2_n5[jmax];
        
        detM     =  M[0]*M[3] - M[1]*M[2];
        bi[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
        bi[k][1] = (M[0]*P[1] - M[2]*P[0])/detM;
    }

    if (DEBUG)
    {
        /*
        printf("m11 = %f\n",m11[1][jmax]);
        printf("m12 = %f\n",m12[1][jmax]);
        printf("m21 = %f\n",m21[1][jmax]);
        printf("m22 = %f\n",m22[1][jmax]);
        printf("P[0] = %f\n", max_A[1]  - p1tmp[1][jmax]);
        printf("P[1] = %f\n", max_dA[1] - p2tmp[1][jmax]);
        */
        printf("-------------------\n");
        printf("NQC coefficients:  \n");
        printf("-------------------\n");
        printf("a1 = %f\n",ai[1][0]);
        printf("a2 = %f\n",ai[1][1]);
        printf("b1 = %f\n",bi[1][0]);
        printf("b2 = %f\n",bi[1][1]);
    }
    
    for (int k=35;k--;)
    {
        for (int j=0; j<t_length-1;j++)
        {
	  o[k][j].dat[0] = 1. + ai[k][0]*n1[j] + ai[k][1]*n2[j];
	  o[k][j].dat[1] = 0. + bi[k][0]*n4[j] + bi[k][1]*n5[j];
        }
    }

    return o;
    
}

vector<gsl_complex> hlmNQC(double  nu,
                           double  r,
                           double  prstar,
                           double  Omega,
                           double  ddotr
                           ){
    
    /********************************************************************************
     * This function computes the NQC corrections to the RWZ multipolar waveform.   *
     * References: Nagar, Damour, Reisswig, Pollney http://arxiv.org/abs/1506.08457 *
     ********************************************************************************/
    
    vector<double> n(6);
    const int kmax = 35;
    double a1;
    double a2;
    double a3;
    double b1;
    double b2;
    double b3;
    
    /** FITS: possibly to be improved further. Current fits: 9/02/2016 */
    const double xnu  = 1-4*nu;
    const double xnu2 = (1-4*nu)*(1-4*nu);
    
    /** NQC multipolar correction factor */
    vector<gsl_complex> psilmnqc(kmax);
    
    for (int i=kmax;i--;)
    {
        switch (i)
        {
                /** l=2 -------------------------------------------------------------------
                 *(2,1) */
            case 0:
                
                a1 =  0.0162387198*( 7.32653082*xnu2 + 1.19616248*xnu + 0.73496656);
                a2 =                -1.80492460*xnu2 + 1.78172686*xnu + 0.30865284 ;
                a3 =                                                            0.0;
                
                b1 =  -0.0647955017*(3.59934444*xnu2 - 4.08628784*xnu + 1.37890907);
                b2 =   1.3410693180*(0.38491989*xnu2 + 0.10969453*xnu + 0.97513971);
                b3 =                                                            0.0;
                
                // NQC corrections to the modulus
                n[0] = (prstar/(r*Omega))*(prstar/(r*Omega));
                n[1] = ddotr/(r*Omega*Omega);
                n[2] = n[0]*prstar*prstar;
                
                
                //% NQC corrections to the phase
                n[3] = prstar/(r*Omega);
                n[4] = n[3]*cbrt(Omega*Omega);
                n[5] = n[4]*prstar*prstar;
                    
                break;
                
                /* (2,2) */
            case 1:
                
                a1   = -0.0805236959*( 1 - 2.00332326*xnu2)/( 1 + 3.08595088*xnu2);
                a2   =  1.5299534255*( 1 + 1.16438929*xnu2)/( 1 + 1.92033923*xnu2);
                a3   =  0.0;
                
                b1   = 0.146768094955*( 0.07417121*xnu + 1.01691256);
                b2   = 0.896911234248*(-0.61072011*xnu + 0.94295129);
                b3   = 0.0;
                
                // NQC corrections to the modulus
                n[0] = (prstar/(r*Omega))*(prstar/(r*Omega));
                n[1] = ddotr/(r*Omega*Omega);
                n[2] = n[0]*prstar*prstar;
                
                
                //% NQC corrections to the phase
                n[3] = prstar/(r*Omega);
                n[4] = n[3]*(r*Omega)*(r*Omega);
                n[5] = n[4]*prstar*prstar;
                
                break;
                
                /** l=3 -------------------------------------------------------------------
                 * (3,3) */
            case 4:
                
                a1   = -0.0377680000*(1 - 14.61548907*xnu2)/( 1 + 2.44559263*xnu2);
                a2   =  1.9898000000*(1 + 2.09750346 *xnu2)/( 1 + 2.57489466*xnu2);
                a3   =  0.0;
                
                b1   = 0.1418400000*(1.07430512 - 1.23906804*xnu + 4.44910652*xnu2);
                b2   = 0.6191300000*(0.80672432 + 4.07432829*xnu - 7.47270977*xnu2);
                b3   = 0.0;
                
                // NQC corrections to the modulus
                n[0] = (prstar/(r*Omega))*(prstar/(r*Omega));
                n[1] = ddotr/(r*Omega*Omega);
                n[2] = n[0]*prstar*prstar;
                
                
                //% NQC corrections to the phase
                n[3] = prstar/(r*Omega);
                n[4] = n[3]*cbrt(Omega*Omega);
                n[5] = n[4]*prstar*prstar;
                
                break;
                
            default:
                
                a1   = 0.;
                a2   = 0.;
                a3   = 0.;
                
                b1   = 0.;
                b2   = 0.;
                b3   = 0.;
                
                n[0] = 0.;
                n[1] = 0.;
                n[2] = 0.;
                
                n[3] = 0.;
                n[4] = 0.;
                n[5] = 0.;
                
                break;
        }
        
        psilmnqc[i].dat[0] = 1. + a1*n[0] + a2*n[1] + a3*n[2];
        psilmnqc[i].dat[1] = 0. + b1*n[3] + b2*n[4] + b3*n[5];
    }
    
    return psilmnqc;
}

double dtnqc_fit(const double chi,
                 const double chi0
                 ){
    
    /** Function providing a fit of Deltat_NQC vs chi, via a simple rational function. */
    
    const double n1    = -16.06288206;
    const double d1    = -4.04266459;
    double x     = chi-chi0;
    double dtnqc = (1.+n1*x)/(1.+d1*x);
    
    return dtnqc;
    
}

                            /*******************************************************
                             * Checkpoint Charlie: You are leaving the NQC sector  *
                             *******************************************************/

vector<gsl_complex> hlm(double       t,
                        const double phi,
                        const double r,
                        const double pph,
                        const double prstar,
                        double       Omega,
                        const double ddotr,
                        const double H,
                        const double Heff,
                        const double jhat,
                        const double rw,
                        void         *params
                        ){
    int kmax = 35;
    
    vector<gsl_complex> hlm(kmax);
    
    double nu            = (*(TEOBResumParams *)params).nu;
    bool tidal_flag      = (*(TEOBResumParams *)params).flags.tidal;
    bool spin_flag       = (*(TEOBResumParams *)params).flags.spin;
    bool speedytail_flag = (*(TEOBResumParams *)params).flags.speedy;
    
    double source[]      = {
        jhat,Heff,
        Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,
        Heff,jhat,Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,jhat,Heff,
        Heff,jhat,Heff,jhat,Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,jhat,Heff,jhat,Heff};
    
    /** Newtonian waveform */
    vector<gsl_complex> hNewt = hlmNewt(rw,Omega,phi, nu, tidal_flag, spin_flag);
    
    /** Compute corrections */
    double x = rw*Omega*rw*Omega;
    vector<double> flm(35);
    if (spin_flag==1)
    {
        flm = s_flm(x, params);
    }
    else
    {
        flm = f_lm(x, nu);
    }
    
    /** Computing the tail */
    const double r0    = 1.213061319425267e+00;
    const double Hreal = H * nu;
    vector<gsl_complex> tlm(kmax);
    if (speedytail_flag==0)
    {
        tlm = hhatlmtail(Omega,Hreal, r0, L, M);
    }
    else
    {
        tlm = speedyTail(Omega,Hreal, r0, L, M);
    }
    
    /** Residual phase corrections delta_{lm} */
    const vector<double> EOBdeltalm = deltalm(Hreal, Omega, nu);
    
    vector<gsl_complex> h_NQC(kmax);
    if (tidal_flag==0)
    {
        h_NQC = hlmNQC(nu,r,prstar,Omega,ddotr);
    }
    
    for (int k=35; k--;)
    {
        tlm[k].dat[1] += EOBdeltalm[k];
        
        /** Compute \hat{h}_lm */
        hlm[k].dat[0] =   hNewt[k].dat[0] * flm[k] * source[k] * tlm[k].dat[0] ;
        hlm[k].dat[1] =   hNewt[k].dat[1] + tlm[k].dat[1] ;
        hlm[k].dat[1] = - hlm[k].dat[1]; /** Minus sign by convention */
        
        /** NQC correction */
        if (tidal_flag==0 && spin_flag==0)
        {
	    hlm[k].dat[0] *= h_NQC[k].dat[0];
            hlm[k].dat[1] -= h_NQC[k].dat[1];
        }
    }
    
    if (tidal_flag==1)
    {
        
        /** Compute tidal contribution */
        vector<double> hlmtidal   = hlm_Tidal(x, params);
        
        /** Update waveform */
        double p2 = sqrt(1-4*nu);
        for (int k=35; k--;)
        {
            switch (k)
            {
                case 0:
                    hlm[0].dat[0] *= p2;
                    break;
                case 2:
                    hlm[2].dat[0] *= p2;
                    break;
                case 4:
                    hlm[4].dat[0] *= p2;
                    break;
                default:
                    break;
            }
            
            double rad_tid = hNewt[k].dat[0] * tlm[k].dat[0] * hlmtidal[k] ;
            hlm[k].dat[0] += rad_tid;
        }
    }
    
    return hlm;
}



