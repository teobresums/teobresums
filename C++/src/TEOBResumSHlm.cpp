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

#include "TEOBResumS.h"

using namespace::std;

/** This routine computes the waveform, according to eq. (76) of https://arxiv.org/abs/1406.6913v1 */


vector<gsl_complex> hlmNewt(const double r,
                            const double Omega,
                            const double phi,
                            const double nu,
                            bool         tidal_flag)
{
    
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
    
    
    /** Polynomials in nu */
    const double p1 = 1.;
    double       p2 = sqrt(1.-4.*nu);
    const double p3 = (3.*nu-1.);
    const double p4 = (2.*nu-1.)*sqrt(1.-4.*nu);
    const double p5 = 1.-5.*nu+5.*nu2;
    const double p6 = (1.-4.*nu+3.*nu2)*sqrt(1.-4.*nu);
    const double p7 = 7.*nu3 - 14.*nu2 + 7.*nu -1.;
    
    if (tidal_flag==true)
    {
        p2 = 1.;
    }
    
    int kmax = 35;
    const double phix2 = 2. * phi;
    const double phix3 = 3. * phi;
    const double phix4 = 4. * phi;
    const double phix5 = 5. * phi;
    const double phix6 = 6. * phi;
    const double phix7 = 7. * phi;
    
    double M[] = {
        phi,phix2,
        phi,phix2,phix3,
        phi,phix2,phix3,phix4,
        phi,phix2,phix3,phix4,phix5,
        phi,phix2,phix3,phix4,phix5,phix6,
        phi,phix2,phix3,phix4,phix5,phix6,phix7,
        phi,phix2,phix3,phix4,phix5,phix6,phix7,8.*phi};
    
    vector<gsl_complex> hlmNewt(kmax);
    
    const double pv23 = p2 * vphi3;
    const double pv34 = p3 * vphi4;
    const double pv45 = p4 * vphi5;
    const double pv56 = p5 * vphi6;
    const double pv67 = p6 * vphi7;
    const double pv78 = p7 * vphi8;
    
    
    /** Compute hlmNewt (without phase factor) in complex Polar coords
     * l=2 ------------------------------------------------------------------ */
    
    hlmNewt[0].dat[1] = 3.*pi/2. - M[0];
    hlmNewt[0].dat[0] = 8./3.*sqrt(pi/5.) * pv23;
    
    hlmNewt[1].dat[1] = pi - M[1];
    hlmNewt[1].dat[0] = 8.*sqrt(pi/5.) * p1 * vphi2;
    
    /** l=3 ------------------------------------------------------------------ */
    
    hlmNewt[2].dat[1] = 3.*pi/2. - M[2];
    hlmNewt[2].dat[0] = 1./3.*sqrt(2.*pi/35.) * pv23;
    
    hlmNewt[3].dat[1] = - M[3];
    hlmNewt[3].dat[0] = 8./3.*sqrt(pi/7.) * pv34;
    
    hlmNewt[4].dat[1] = pi/2. - M[4];
    hlmNewt[4].dat[0] = 3.*sqrt(6.*pi/7.) * pv23;
    
    /** l=4 ------------------------------------------------------------------ */
    
    hlmNewt[5].dat[1] = pi/2. - M[5];
    hlmNewt[5].dat[0] = 1./105.*sqrt(2.*pi) * pv45;
    
    hlmNewt[6].dat[1] = - M[6];
    hlmNewt[6].dat[0] = 8./63.*sqrt(pi) * pv34;
    
    hlmNewt[7].dat[1] = 3.*pi/2. - M[7];
    hlmNewt[7].dat[0] = 9./5*sqrt(2*pi/7.) * pv45;
    
    hlmNewt[8].dat[1] = pi - M[8];
    hlmNewt[8].dat[0] = 64./9.*sqrt(pi/7.) * pv34;
    
    /** l=5 ------------------------------------------------------------------ */
    
    hlmNewt[9].dat[1] = pi/2. - M[9];
    hlmNewt[9].dat[0] = 1./180.*sqrt(pi/77.) *pv45;
    
    hlmNewt[10].dat[1] = pi - M[10];
    hlmNewt[10].dat[0] = 16./135.*sqrt(pi/11.) *pv56;
    
    hlmNewt[11].dat[1] = 3.*pi/2. - M[11];
    hlmNewt[11].dat[0] = 9./20.*sqrt(3*pi/22.) *pv45;
    
    hlmNewt[12].dat[1] = - M[12];
    hlmNewt[12].dat[0] = 256./45.*sqrt(pi/33.) *pv56;
    
    hlmNewt[13].dat[1] = pi/2. - M[13];
    hlmNewt[13].dat[0] = 125./12.*sqrt(5.*pi/66.) * pv45;
    
    /** l=6 ------------------------------------------------------------------ */
    
    hlmNewt[14].dat[1] = 3.*pi/2. - M[14];
    hlmNewt[14].dat[0] = 1./2079.*sqrt(2.*pi/65.) *pv67;
    
    hlmNewt[15].dat[1] = pi - M[15];
    hlmNewt[15].dat[0] = 16./1485.*sqrt(pi/13.) *pv56;
    
    hlmNewt[16].dat[1] = pi/2. - M[16];
    hlmNewt[16].dat[0] = 81./385.*sqrt(pi/13.)  *pv67;
    
    hlmNewt[17].dat[1] = - M[17];
    hlmNewt[17].dat[0] = 1024./495.*sqrt(2.*pi/195.) *pv56;
    
    hlmNewt[18].dat[1] = 3.*pi/2. - M[18];
    hlmNewt[18].dat[0] = 625./63*sqrt(5.*pi/429.) * pv67;
    
    hlmNewt[19].dat[1] = pi - M[19];
    hlmNewt[19].dat[0] = 432./5*sqrt(pi/715.)  *pv56;
    
    /** l=7 ------------------------------------------------------------------ */
    
    hlmNewt[20].dat[1] = 3.*pi/2. - M[20];
    hlmNewt[20].dat[0] = 1./108108.*sqrt(pi/10.)  *pv67;
    
    hlmNewt[21].dat[1] = - M[21];
    hlmNewt[21].dat[0] = 8./3003.*sqrt(pi/15.)    *pv78;
    
    hlmNewt[22].dat[1] = pi/2. - M[22];
    hlmNewt[22].dat[0] = 243./20020.*sqrt(3.*pi/10.) *pv67;
    
    hlmNewt[23].dat[1] = pi - M[23];
    hlmNewt[23].dat[0] = 1024./1365.*sqrt(2.*pi/165.)  *pv78;
    
    hlmNewt[24].dat[1] = 3.*pi/2. - M[24];
    hlmNewt[24].dat[0] = 3125./3276.*sqrt(5.*pi/66.) *pv67;
    
    hlmNewt[25].dat[1] = - M[25];
    hlmNewt[25].dat[0] = 648./35.*sqrt(3.*pi/715.)  *pv78;
    
    hlmNewt[26].dat[1] = pi/2. - M[26];
    hlmNewt[26].dat[0] = 16807./180.*sqrt(7.*pi/4290.) *pv67;
    
    /** l=8 ------------------------------------------------------------------ */
    
    hlmNewt[27].dat[1] = pi - M[27];
    hlmNewt[27].dat[0] = 131072./315.*sqrt(2.*pi/17017.) *pv78;
    
    hlmNewt[28].dat[1] = - M[28];
    hlmNewt[28].dat[0] = 0.;
    
    hlmNewt[29].dat[1] = - M[29];
    hlmNewt[29].dat[0] = 0.;
    
    hlmNewt[30].dat[1] = - M[30];
    hlmNewt[30].dat[0] = 0.;
    
    hlmNewt[31].dat[1] = - M[31];
    hlmNewt[31].dat[0] = 0.;
    
    hlmNewt[32].dat[1] = - M[32];
    hlmNewt[32].dat[0] = 0.;
    
    hlmNewt[33].dat[1] = - M[33];
    hlmNewt[33].dat[0] = 0.;
    
    hlmNewt[34].dat[1] = - M[34];
    hlmNewt[34].dat[0] = 0.;
    
    return hlmNewt;
}

vector<gsl_complex> hhatlmtail(const double Omega,
                               const double Hreal,
                               const double bphys,
                               const int    L[],
                               const int    M[])
{
    
    /***********************************************************
     *   Computes the tail contribution to the resummed wave.  *
     *   Reference: Damour, Iyer & Nagar, PRD 79, 064004 (2009)*
     ***********************************************************/
    
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
        
        ratio_rad = num_rad.val-denom_rad.val;
        ratio_ang = num_phase.val-0.;
        
        tlm_rad   = ratio_rad + pi * hhatk;
        tlm_phase = ratio_ang + 2.*hhatk*log(2.*k*bphys);
        
        tlm[i].dat[0] = exp(tlm_rad);
        tlm[i].dat[1] = tlm_phase;
    }
    return tlm;
}



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
    
    double nu           = (*(TEOBResumParams *)params).nu;
    int tidal_flag      = (*(TEOBResumParams *)params).flags.tidal;
    int spin_flag       = (*(TEOBResumParams *)params).flags.spin;
    int NQC_flag        = (*(TEOBResumParams *)params).flags.NQC;
    int speedytail_flag = (*(TEOBResumParams *)params).flags.speedy;
    
    double source[]     = {
        jhat,Heff,
        Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,
        Heff,jhat,Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,jhat,Heff,
        Heff,jhat,Heff,jhat,Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,jhat,Heff,jhat,Heff};
    
    /** Newtonian waveform */
    vector<gsl_complex> hNewt = hlmNewt( rw,Omega,phi, nu,tidal_flag);
    
    /** Compute corrections */
    double x = gsl_pow_int(rw*Omega,2);
    vector<double> flm(35);
    if (spin_flag==1)
    {
        flm = s_flm(x,params);
    }
    else
    {
        flm = f_lm(x,nu);
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
    if (NQC_flag==1)
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
        if (NQC_flag==1)
        {
            hlm[k].dat[0] *= h_NQC[k].dat[0];
            hlm[k].dat[1] -= h_NQC[k].dat[1];
        }
    }
    
    if (tidal_flag==1)
    {
        
        /** Compute tidal contribution */
        vector<double> hlmtidal = hlm_Tidal(x, params);
        
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
    
    
    // l=2 -------------------------------------------------------------------*/
    
    hA[1]     = 2 * khatA_2 *(XA/XB+3);
    hB[1]     = 2 * khatB_2 *(XB/XA+3);
    
    betaA1[1] = (-202. + 560*XA - 340*XA*XA + 45*XA*XA*XA)/(42*(3-2*XA));
    betaB1[1] = (-202. + 560*XB - 340*XB*XB + 45*XB*XB*XB)/(42*(3-2*XB));
    
    hA[0]     = 3 * khatA_2 * XB * (3-4*XA)/XA;
    hB[0]     = 3 * khatB_2 * XA * (3-4*XB)/XB;
    
    
    /** l=3 ------------------------------------------------------------------*/
    
    hA[2] = hA[4];
    hB[2] = hB[4];
    
    hA[4] = 12 * khatA_2 * XB*XB/XA;
    hB[4] = 12 * khatB_2 * XA*XA/XB;
    
    
    
    /** l=2 ------------------------------------------------------------------
     * (2,1) */
    hTidallm[0] = ( -hA[0] + hB[0] )*x5;
    
    /* (2,2) */
    hTidallm[1] = ( hA[1]*(1. + betaA1[1]*x) + hB[1]*(1. + betaB1[1]*x) )*x5;
    
    /** l=3 ------------------------------------------------------------------
     * (3,1) */
    hTidallm[2] = ( -hA[2] + hB[2] )*x5;
    
    /* (3,3) */
    hTidallm[4] = ( -hA[4] + hB[4] )*x5;
    
    return hTidallm;
    
}

