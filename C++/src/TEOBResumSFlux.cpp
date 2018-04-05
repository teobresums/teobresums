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

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "TEOBResumS.h"

double flux(const double x,
            const double Omega,
            const double r_omega,
            const double E,
            const double Heff,
            const double jhat,
            const double r,
            const double prstar,
            const double ddotr,
            double source[],
            void *params)
 /*    This function computes the Newton-normalized energy flux according
       to the usual Damour-Iyer-Nagar resummation procedure, see also
       Damour, Nagar & Bernuzzi, PRD 87 (2013). This function is used
       only for NONSPINNING systems.
       Arguments of the function are:
       x       => PN variable
       Omega   => Orbital frequency
       r_omega => Kepler-modified radial separation
       E       => Energy
       Heff    => Effective energy
       jhat    => Newton-Normalized angular momentum       
       r       => radial separation
       pr_star => radial momentum (conjugate to r*)
       ddotr   => second time-derivative of the radial separation
       source[]=> effective energy (l+m=even) or angular momentum (l+m=odd)
     */
  
  {
    
    bool   tidal_flag   = (*(TEOBResumParams *)params).flags.tidal;
    double nu           = (*(TEOBResumParams *)params).nu;
    double X1           = (*(TEOBResumParams *)params).X1;
    double X2           = (*(TEOBResumParams *)params).X2;
    double X12          = X1-X2;
    double               Flm;
    double               Modhhatlm;
    const vector<double> flm          = f_lm(x,nu);
    const vector<double> FNewtlm      = FlmNewt(x,params);
    const double         FNewt22      = FNewtlm[1];
    vector<double>       MTlm         = Tlm(E*Omega);
    vector<double>       hlmTidal     = hlm_Tidal(x,params);
    
    /** Compute NQC correction to the modulus of the (l,m) waveform */
    vector<gsl_complex> hlm_NQC = hlmNQC(nu,r,prstar,Omega,ddotr);
    double SFlm=0.; 

    for (int k=35; k--;)
    {
        /** Compute modulus of hhat_lm */
        Modhhatlm = source[k] * MTlm[k] * flm[k];
        
        /** NQCs for the non tidal case*/
        if (tidal_flag==false)
        {
            switch (k)
            {
//                case 0:
//                    Modhhatlm *= hlm_NQC[k].dat[0];
//                    break;
                case 1:
                    Modhhatlm *= hlm_NQC[k].dat[0];
                    break;
//                case 4:
//                    Modhhatlm *= hlm_NQC[k].dat[0];
//                    break;
                default:
                    Modhhatlm *= 1.;
                    break;
            }
        }
        
        /* Re-introduce the nu-dependent prefactor in the tidal case */
        if (tidal_flag==true)
        {
            switch (k)
            {
                case 0: // (2,1)
		    Modhhatlm *=X12;
                    break;
                case 2: // (3,1)
		    Modhhatlm *=X12;
                    break;
                case 4: // (3,3)
		    Modhhatlm *=X12;
                    break;
                default: Modhhatlm *= 1.;
                    break;
            }
            
            Modhhatlm += MTlm[k]*hlmTidal[k];
            
        }

        /** Total flux multipoles */
        Flm = (Modhhatlm * Modhhatlm) * FNewtlm[k];
        SFlm += Flm;
    }
    
    /** Sum over multipoles and normalize to the 22 Newtonian multipole */
    double hatf = SFlm/(FNewt22);

    if (tidal_flag==false)
    {
      /* Add horizon flux in the BBH case */
        double hatFH = HorizonFlux(x,Heff,jhat,nu);
        hatf += hatFH;
    }
    
    double Fphi = -32.0/5.0* nu * gsl_pow_int(r_omega,4) * gsl_pow_int(Omega,5) * hatf;
    
    return Fphi;
}

/** Flux calculation for spinning systems */
double s_Flux(double x,
              double Omega,
              double r_omega,
              double E,
              double Heff,
              double jhat,
              double r,
              double pr_star,
              double ddotr,
              void *params){
    /* This function computes the Newton-normalized energy flux according
       to the usual Damour-Iyer-Nagar resummation procedure.
       Arguments of the function are:
       x       => PN variable
       r_omega => Kepler-modified radial separation
       Omega   => Orbital frequency
       E       => Energy
       Heff    => Effective energy
       jhat    => Newton-Normalized angular momentum       
       r       => radial separation
       pr_star => radial momentum (conjugate to r*)
       ddotr   => second time-derivative of the radial separation
     */
    
    double nu       = (*(TEOBResumParams *)params).nu;
    double chi1     = (*(TEOBResumParams *)params).chi1;
    double chi2     = (*(TEOBResumParams *)params).chi2;
    double X1       = (*(TEOBResumParams *)params).X1;
    double X2       = (*(TEOBResumParams *)params).X2;
    bool tidal_flag = (*(TEOBResumParams *)params).flags.tidal;

    /* effective source: jhat for l+m=odd and Heff when l+m=even */
    double prefact[] = {
        jhat, Heff,
        Heff, jhat, Heff,
        jhat, Heff, jhat, Heff,
        Heff, jhat, Heff, jhat, Heff,
        jhat, Heff, jhat, Heff, jhat, Heff,
        Heff, jhat, Heff, jhat, Heff, jhat, Heff,
        jhat, Heff, jhat, Heff, jhat, Heff, jhat, Heff};
    
    /* Temporary fix: NQC function put to 1. This will change with iterations */
    double hnqclm = 1.;
    double Flm;
    double Modhhatlm;
    double X12 = X1-X2;
    
    vector<double> flm     = s_flm(x,params);
    vector<double> FNewtlm = FlmNewt(x,params);
    double FNewt22         = FNewtlm[1];
    
    vector<double> MTlm     = Tlm(E*Omega);
    vector<double> hlmTidal = hlm_Tidal(x,params);
    
    double SFlm = 0.;
    for (int k=35; k--;) {
      /** Compute modulus of hhat_lm */
        Modhhatlm = prefact[k] * MTlm[k] * flm[k];
        
        /** NQCs for the non tidal case */
        if (tidal_flag==false)
        {
            Modhhatlm *= hnqclm;
        }
        
        /* Adding the tidal waveform amplitude to the point-mass baseline */
        if (tidal_flag==true)
        {            
            Modhhatlm += MTlm[k]*hlmTidal[k];
        }
        
        
        /** Total flux multipoles */
        Flm = (Modhhatlm * Modhhatlm) * FNewtlm[k];
        
        SFlm += Flm;
    }
    
    // Sum over multipoles and normalize to the 22 Newtonian multipole
    double hatf = SFlm/(FNewt22);
    
    if (tidal_flag==false)
    {
        /* Add horizon flux in the BBH case */
        double hatFH = s_HorizonFlux(x, Heff, jhat, nu, X1, X2, chi1, chi2);
        hatf += hatFH;
    }
    
    double Fphi = -32.0/5.0* nu * gsl_pow_int(r_omega,4) * gsl_pow_int(Omega,5) * hatf;
    
    return Fphi;
}

vector<double> FlmNewt(const double x, void *params)
{
  /** This function computes the Newtonian prefactors in the flux.
      The multipolar Newtonian prefactors can be obtained from 
      Eq.(4) of Damour-Iyer-Nagar, PRD 79, 064004 (2009) [DIN]*/
    double nu       = (*(TEOBResumParams *)params).nu;
    bool tidal_flag = (*(TEOBResumParams *)params).flags.tidal;
    bool spin_flag  = (*(TEOBResumParams *)params).flags.spin;
    
    /** Shorthands*/
    const double nu2 = nu*nu;
    const double nu3 = nu2*nu;
    const double x5  = x*x*x*x*x;
    const double x6  = x*x5;
    const double x7  = x*x6;
    const double x8  = x*x7;
    const double x9  = x*x8;
    const double x10 = x*x9;
    const double x11 = x*x10;
    const double x12 = x*x11;

    /* setting up numerical constants, the c_ell+epsilon, Eq.(7) of DIN */
    double sp2 = 0.0;
    double sp4 = 0.0;
    const double sp3 = (1.-3.*nu)*(1.-3.*nu);
    const double sp5 = (1.-5.*nu+5.*nu2)*(1.-5.*nu+5.*nu2);
    const double sp6 = (1-4*nu)*(3*nu2-4*nu +1)*(3*nu2-4*nu +1);    
    const double sp7 = (1 - 7*nu + 14*nu2 - 7*nu3)*(1 - 7*nu + 14*nu2 - 7*nu3);
    const double sp8 = (1 - 4*nu)*(1 - 6*nu + 10*nu2 - 4*nu3)*(1 - 6*nu + 10*nu2 - 4*nu3);
    
    vector<double> Nlm(35);
    
    if (spin_flag==true)
    {
        sp2 = 1.;
        sp4 = (2.*nu-1.)*(2.*nu-1.);	
    }
    else
    {   /** Nonspinning case*/
        if (tidal_flag==false)
        {
            sp2 = 1.-4.*nu;
            sp4 = (1.-4.*nu)*(1.-2.*nu)*(1.-2.*nu);
        }
        else
        {
            sp2 = 1.;
            sp4 = (1.-4.*nu)*(1.-2.*nu)*(1.-2.*nu);
        }
    }
    
    /** Newtonian partial fluxes*/

    /* l=2 */
    Nlm[0]  = 8./45                                 * sp2 * x6 ; /*(2,1)*/
    Nlm[1]  = 32./5.                                *       x5 ; /*(2,1)*/

    /* l=3 */
    Nlm[2]  = 1./1260                               * sp2 * x6 ; /*(3,1)*/
    Nlm[3]  = 32./63                                * sp3 * x7 ; /*(3,2)*/
    Nlm[4]  = 243./28                               * sp2 * x6 ; /*(3,3)*/

    /* l=4 */
    Nlm[5]  = 1./44100                              * sp4 * x8 ; /*(4,1)*/
    Nlm[6]  = 32./3969                              * sp3 * x7 ; /*(4,2)*/
    Nlm[7]  = 729./700                              * sp4 * x8 ; /*(4,3)*/
    Nlm[8]  = 8192./567                             * sp3 * x7 ; /*(4,4)*/

    /* l=5 */
    Nlm[9]  = 1./19958400                           * sp4 * x8 ; /*(5,1)*/
    Nlm[10] = 256./400950                           * sp5 * x9 ; /*(5,2)*/
    Nlm[11] = 2187./70400                           * sp4 * x8 ; /*(5,3)*/
    Nlm[12] = 131072./66825                         * sp5 * x9 ; /*(5,4)*/
    Nlm[13] = 1953125./76032                        * sp4 * x8 ; /*(5,5)*/

    /* l=6 */
    Nlm[14] = 1./1123782660                         * sp6 * x10; /*(6,1)*/
    Nlm[15] = 128./28667925                         * sp5 * x9 ; /*(6,2)*/
    Nlm[16] = 59049./15415400                       * sp6 * x10; /*(6,3)*/
    Nlm[17] = 4194304./47779875                     * sp5 * x9 ; /*(6,4)*/
    Nlm[18] = 48828125./13621608                    * sp6 * x10; /*(6,5)*/
    Nlm[19] = 839808./17875                         * sp5 * x9 ; /*(6,6)*/

    /* l=7 */
    Nlm[20] = 1./9.3498717312e11                    * sp6 * x10; /*(7,1)*/
    Nlm[21] = 32./1.35270135e8                      * sp7 * x11; /*(7,2)*/
    Nlm[22] = 1594323./3.2064032e10                 * sp6 * x10; /*(7,3)*/
    Nlm[23] = 4194304./3.07432125e8                 * sp7 * x11; /*(7,4)*/
    Nlm[24] = 1220703125./5.666588928e9             * sp6 * x10; /*(7,5)*/
    Nlm[25] = 5668704./875875                       * sp7 * x11; /*(7,6)*/
    Nlm[26] = 96889010407./1111968000               * sp6 * x10; /*(7,7)*/

    /* l=8 */
    Nlm[27] = 1./8.174459284992e13                  * sp8 * x12; /*(8,1)*/
    Nlm[28] = 32./3.4493884425e10                   * sp7 * x11; /*(8,2)*/
    Nlm[29] = 177147./3.96428032e10                 * sp8 * x12; /*(8,3)*/
    Nlm[30] = 4194304./1.5679038375e10              * sp7 * x11; /*(8,4)*/
    Nlm[31] = 30517578125./8.00296713216e11         * sp8 * x12; /*(8,5)*/
    Nlm[32] = 51018336./1.04229125e8                * sp7 * x11; /*(8,6)*/
    Nlm[33] = 4747561509943./4.083146496e11         * sp8 * x12; /*(8,7)*/
    Nlm[34] = 274877906944./1688511825.             * sp7 * x11; /*(8,8)*/
    
    return Nlm;
}


double HorizonFlux(const double x, const double Heff, const double jhat, const double nu)
{
    
    /** Compute horizon-absorbed fluxes.
     *
     *
     *
     * Nagar & Akcay, PRD 85, 044025 (2012)
     * Bernuzzi, Nagar & Zenginoglu, PRD 86, 104038 (2012)
     */
    
    const int kmax   = 2;
    vector<double> rhoHlm(kmax);
    vector<double> FlmHLO(kmax);
    vector<double> FlmH(kmax);
    
    /** Shorthands */
    double nu2 = nu*nu;
    double nu3 = nu*nu2;
    double x2  = x*x;
    double x3  = x*x2;
    double x4  = x*x3;
    double x5  = x*x4;
    double x9  = x4*x5;
    double x10 = x*x9;
    
    const int k22 = 1;
    const int k21 = 0;
    
    /** The Newtonian asymptotic contribution */
    const double FNewt22 = 32./5.*x5;
    
    /** Compute leading-order part (nu-dependent) */
    FlmHLO[1] = 32./5.*(1-4*nu+2*nu2)*x9;
    FlmHLO[0] = 32./5.*(1-4*nu+2*nu2)*x10;
    
    /** Compute rho_lm */
    vector<double> c1(kmax);
    vector<double> c2(kmax);
    vector<double> c3(kmax);
    vector<double> c4(kmax);
    
    c1[0] = 0.58121;
    c2[0] = 1.01059;
    c3[0] = 7.955729;
    c4[0] = 1.650228;
    
    c1[1] = (4.-21.*nu + 27.*nu2 - 8.*nu3)/(4.*(1.-4.*nu+2.*nu2));
    c2[1] =  4.78752;
    c3[1] = 26.760136;
    c4[1] = 43.861478;
    
    rhoHlm[1] = 1. + c1[1]*x + c2[1]*x2 + c3[1]*x3 + c4[1]*x4;
    rhoHlm[0] = 1. + c1[0]*x + c2[0]*x2 + c3[0]*x3 + c4[0]*x4;
    
    /** Compute horizon multipolar flux (only l=2) */
    const double Heff2 = Heff*Heff;
    const double jhat2 = jhat*jhat;
    
    FlmH[k22] = FlmHLO[k22] * Heff2 * gsl_pow_int(rhoHlm[k22],4);
    FlmH[k21] = FlmHLO[k21] * jhat2 * gsl_pow_int(rhoHlm[k21],4);
    
    /** Sum over multipoles and normalize to the 22 Newtonian multipole */
    double const hatFH = (FlmH[0]+FlmH[1])/FNewt22;
    
    return hatFH;
}


double s_HorizonFlux(double x, double Heff, double jhat, double nu, double X1, double X2, double chi1, double chi2)
{
    /** Higher order terms can be found in ref Chatziiaonnou, Poisson, Yunes PRD 94, 084043 (2016) */
    
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x3*x;
    double x5 = x4*x;
    double v5 = sqrt(x5);
    
    vector<double> cv5(2);
    vector<double> cv8(2);
    
    /** Coefficients of the v^5 term (Alvi leading order) */
    cv5[0] = -1./4.*chi1*(1.+3.*chi1*chi1)*X1*X1*X1;
    cv5[1] = -1./4.*chi2*(1.+3.*chi2*chi2)*X2*X2*X2;
    
    /** Coefficients of the v^8=x^4 term
        eq.46 of https://arxiv.org/pdf/1112.2840.pdf */
    cv8[0] = 0.5*(1.+sqrt(1.-chi1*chi1))*(1.+3.*chi1*chi1)*X1*X1*X1*X1;
    cv8[1] = 0.5*(1.+sqrt(1.-chi2*chi2))*(1.+3.*chi2*chi2)*X2*X2*X2*X2;
    
    double FH22_S = (cv5[0]+cv5[1])*v5;
    double FH22   = (cv8[0]+cv8[1])*x4;
    double FH21   =  0.0;
    
    /** Newton-normalized horizon flux: use only l=2 fluxes */
    double hatFH  = FH22_S + FH22 + FH21;
    
    return hatFH;
}

vector<double> Tlm(const double w)
{
    
    /* Factorials evaluated */
    double f14[] = {1.,         1.,          2.,
        6.,         24.,         120.,
        720.,       5040.,       40320.,
        362880.,    3628800.,    39916800.,
        479001600., 6227020800., 87178291200.};
    double hhatk, x2, y, prod;
    vector<double> MTlm(35);
    int kmax = 35;
    
    for (int i=kmax; i--; )
    {
        hhatk = M[i] * w;
        x2    = 4.*hhatk*hhatk;
        prod  = 1.;
        
        for (int j=1; j <= L[i]; j++ )
        {
            prod *= ( j*j + x2 );
        }
        
        y  = 4.*pi*hhatk;
        y /= ( 1. - exp(-y) );
        
        /** Put the three pieces together */
        MTlm[i] = sqrt( 1./(f14[L[i]]*f14[L[i]]) * y * prod );
    }
    
    return MTlm;
}
