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

#include "TEOBResum.h"

/** Flux calculation for non-spinning systems */
double flux(const double x,const double Omega,const double r_omega,const double E, const double Heff,const double jhat,const double r,const double prstar, const double ddotr, double source[],void *params)
{
    
    bool tidal_flag = (*(TEOBResumParams *)params).flags.tidal;
    bool NQC_flag   = (*(TEOBResumParams *)params).flags.NQC;
    double nu       = (*(TEOBResumParams *)params).nu;
    double Flm;
    double Modhhatlm;
    double sqrt_one_4nu = sqrt(1.-4.*nu);
    const vector<double> flm     = f_lm(x,nu);
    const vector<double> FNewtlm = FlmNewt(x,params);
    const double FNewt22         = FNewtlm[1];
    vector<double> MTlm          = Tlm(E*Omega);
    vector<double> hlmTidal      = hlm_Tidal(x,params);
    
    /** Compute NQC correction to the modulus of the (l,m) waveform */
    vector<gsl_complex> hlm_NQC = hlmNQC(nu,r,prstar,Omega,ddotr);
    double SFlm=0.;
    for (int k=35; k--;)
    {
        /** Compute modulus of hhat_lm */
        Modhhatlm = source[k] * MTlm[k] * flm[k];
        
        /** NQCs */
        if (NQC_flag==true)
        {
            switch (k)
            {
                case 0:
                    Modhhatlm *= hlm_NQC[k].dat[0];
                    break;
                case 1:
                    Modhhatlm *= hlm_NQC[k].dat[0];
                    break;
                case 4:
                    Modhhatlm *= hlm_NQC[k].dat[0];
                    break;
                default:
                    Modhhatlm *= 1.;
                    break;
            }
        }
        
        //Make tidal corrections
        if (tidal_flag==true)
        {
            switch (k)
            {
                case 0: // (2,1)
                    Modhhatlm *=sqrt_one_4nu;
                    break;
                case 2: // (3,1)
                    Modhhatlm *=sqrt_one_4nu;
                    break;
                case 4: // (3,3)
                    Modhhatlm *=sqrt_one_4nu;
                    break;
                default: Modhhatlm *= 1.;
                    break;
            }
            
            Modhhatlm += MTlm[k]*hlmTidal[k];
            
        }
        //if (k==1)
        //{
        //    printf(" %s %.8e \n","hlmTidal",hlmTidal[k]);
        //}
        /** Total flux multipoles */
        Flm = (Modhhatlm * Modhhatlm) * FNewtlm[k];
        SFlm += Flm;
    }
    
    /** Sum over multipoles and normalize to the 22 Newtonian multipole */
    double hatf = SFlm/(FNewt22);
    
    if (tidal_flag==false)
    {
        double hatFH = HorizonFlux(x,Heff,jhat,nu);
        hatf += hatFH;
    }
    
    double Fphi = -32.0/5.0* nu * gsl_pow_int(r_omega,4) * gsl_pow_int(Omega,5) * hatf;
    
    return Fphi;
}

/** Flux calculation for spinning systems */
double s_Flux(double x, double Omega, double r_omega, double E, double Heff, double jhat, double r, double pr_star, double ddotr, void *params){
    /*
     % DINFLUX This function computes the Newton.Normalized energy flux according to
     %         the DIN resummation procedure. It is also designed so to add non-QC
     %         and non-K corrections to  (2,2) partial flux.
     %
     %         USAGE:
     %
     %         [Flm F hatF hatF_resum]=DINFlux(x,Omega,E,Heff,jhat,nu,lmax,r,pr_star,ddotr)
     %
     %         where:
     %
     %         x       :: PN argument
     %         Omega   :: Orbital frequency
     %         E       :: Energy
     %         Heff    :: Effective energy
     %         jhat    :: Newton-Normalized angular momentum
     %         nu      :: symmetric mass ratio
     %         lmax    :: maximum l
     %         r       :: EOB radius
     %         pr_star :: radial momentum
     %         ddotr   :: \ddot{r}
     */
    
    double nu       = (*(TEOBResumParams *)params).nu;
    double chi1     = (*(TEOBResumParams *)params).chi1;
    double chi2     = (*(TEOBResumParams *)params).chi2;
    double X1       = (*(TEOBResumParams *)params).X1;
    double X2       = (*(TEOBResumParams *)params).X2;
    bool tidal_flag = (*(TEOBResumParams *)params).flags.tidal;
    bool NQC_flag   = (*(TEOBResumParams *)params).flags.NQC;
    
    double prefact[] = {
        jhat,Heff,
        Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,
        Heff,jhat,Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,jhat,Heff,
        Heff,jhat,Heff,jhat,Heff,jhat,Heff,
        jhat,Heff,jhat,Heff,jhat,Heff,jhat,Heff};
    
    
    double hnqclm = 1.;
    double Flm;
    double Modhhatlm;
    double sqrt_one_4nu = sqrt(1.-4.*nu);
    
    vector<double> flm     = s_flm(x,params);
    vector<double> FNewtlm = FlmNewt(x,params);
    double FNewt22         = FNewtlm[1];
    
    vector<double> MTlm     = Tlm(E*Omega);
    vector<double> hlmTidal = hlm_Tidal(x,params);
    
    double SFlm = 0.;
    for (int k=35; k--;) {
        // Compute modulus of hhat_lm
        Modhhatlm = prefact[k] * MTlm[k] * flm[k];
        
        //Include NQC with flag
        if (NQC_flag==true) {
            Modhhatlm *= hnqclm;
        }
        
        //Make tidal corrections
        if (tidal_flag==true) {
            switch (k) {
                case 0: // (2,1)
                    Modhhatlm *= sqrt_one_4nu;
                    break;
                case 2: // (3,1)
                    Modhhatlm *= sqrt_one_4nu;
                    break;
                case 4: // (3,3)
                    Modhhatlm *= sqrt_one_4nu;
                    break;
                    
                default: Modhhatlm *= 1.;
                    break;
            }
            
            Modhhatlm += MTlm[k]*hlmTidal[k];
            
        }
        
        
        // Total flux multipoles
        Flm = (Modhhatlm * Modhhatlm) * FNewtlm[k];
        
        SFlm += Flm;
    }
    
    // Sum over multipoles and normalize to the 22 Newtonian multipole
    double hatf = SFlm/(FNewt22);
    
    if (tidal_flag==false) {
        double hatFH = s_HorizonFlux(x, Heff, jhat, nu, X1, X2, chi1, chi2);
        hatf += hatFH;
    }
    
    double Fphi = -32.0/5.0* nu * gsl_pow_int(r_omega,4) * gsl_pow_int(Omega,5) * hatf;
    
    return Fphi;
}

vector<double> FlmNewt(const double x, void *params)
{
    
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
    
    double sp2 = 0.0;
    double sp4 = 0.0;
    const double sp3 = (1.-3.*nu)*(1.-3.*nu);
    const double sp5 = (1.-5.*nu+5.*nu2)*(1.-5.*nu+5.*nu2);
    const double sp6 = (1-4*nu)*(3*nu2-4*nu +1)*(3*nu2-4*nu +1);
    const double sp7 = (1 - 7*nu + 14*nu2 - 7*nu3)*(1 - 7*nu + 14*nu2 - 7*nu3);
    const double sp8 = (1 - 4*nu)*(1 - 6*nu + 10*nu2 - 4*nu3)*(1 - 6*nu + 10*nu2 - 4*nu3);
    const double sp9 = (1 - 4*nu)*(1 - 4*nu + 3*nu2)*(1 - 4*nu + 3*nu2);
    
    vector<double> Nlm(35);
    
    if (spin_flag==true)
    {
        sp2 = 1.;
        sp4 = (2*nu-1)*(2*nu-1);
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
    
    Nlm[0]  = 8./45                                 * sp2 * x6 ;
    Nlm[1]  = 32./5.                                *       x5 ;
    
    Nlm[2]  = 1./1260                               * sp2 * x6 ;
    Nlm[3]  = 32./63                                * sp3 * x7 ;
    Nlm[4]  = 243./28                               * sp2 * x6 ;
    
    Nlm[5]  = 1./44100                              * sp4 * x8 ;
    Nlm[6]  = 32./3969                              * sp3 * x7 ;
    Nlm[7]  = 729./700                              * sp4 * x8 ;
    Nlm[8]  = 8192./567                             * sp3 * x7 ;
    
    Nlm[9]  = 1./19958400                           * sp4 * x8 ;
    Nlm[10] = 256./400950                           * sp5 * x9 ;
    Nlm[11] = 2187./70400                           * sp4 * x8 ;
    Nlm[12] = 131072./66825                         * sp5 * x9 ;
    Nlm[13] = 1953125./76032                        * sp4 * x8 ;
    
    Nlm[14] = 1./1123782660                         * sp6 * x10;
    Nlm[15] = 128./28667925                         * sp5 * x9 ;
    Nlm[16] = 59049./15415400                       * sp6 * x10;
    Nlm[17] = 4194304./47779875                     * sp5 * x9 ;
    Nlm[18] = 48828125./13621608                    * sp6 * x10;
    Nlm[19] = 839808./17875                         * sp5 * x9 ;
    
    Nlm[20] = 1./9.3498717312e11                    * sp9 * x10;
    Nlm[21] = 32./1.35270135e8                      * sp7 * x11;
    Nlm[22] = 1594323./3.2064032e10                 * sp9 * x10;
    Nlm[23] = 4194304./3.07432125e8                 * sp7 * x11;
    Nlm[24] = 1220703125./5.666588928e9             * sp6 * x10;
    Nlm[25] = 5668704./875875                       * sp7 * x11;
    Nlm[26] = 96889010407./1111968000               * sp6 * x10;
    
    Nlm[27] = 1./8.174459284992e13                  * sp8 * x12;
    Nlm[28] = 32./3.4493884425e10                   * sp7 * x11;
    Nlm[29] = 177147./3.96428032e10                 * sp8 * x12;
    Nlm[30] = 4194304./1.5679038375e10              * sp7 * x11;
    Nlm[31] = 30517578125./8.00296713216e11         * sp8 * x12;
    Nlm[32] = 51018336./1.04229125e8                * sp7 * x11;
    Nlm[33] = 4747561509943./4.083146496e11         * sp8 * x12;
    Nlm[34] = (7*nu3-14*nu2+7*nu-1)*(7*nu3-14*nu2+7*nu-1)*274877906944./1688511825.*x11;
    
    return Nlm;
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
    
    const double eps = 1.0;
    
    /** Shorthands */
    const double x2  = x*x;
    const double x3  = x*x2;
    const double x4  = x*x3;
    const double x5  = x*x4;
    const double nu2 = nu*nu;
    const double nu3 = nu*nu2;
    const double nu4 = nu*nu3;
    
    /** Compute EulerLogs */
    const double el1 = Eulerlog(x,1);
    const double el2 = Eulerlog(x,2);
    const double el3 = Eulerlog(x,3);
    const double el4 = Eulerlog(x,4);
    const double el5 = Eulerlog(x,5);
    const double el6 = Eulerlog(x,6);
    const double el7 = Eulerlog(x,7);
    
    const int kmax=35; /** Length of vector needed to store the 8 quadropoles, 35=8+7+6+5+4+3+2*/
    vector<double> rholm(kmax);
    vector<double> flm(kmax);
    
    
    /** l=2 ------------------------------------------------------------------
     *  (2,1) */
    rholm[0] = 1. + (-1.0535714285714286 + 0.27380952380952384*nu)*x + (-0.8327841553287982 - 0.7789824263038548*nu + 0.13116496598639457*nu2)*x2 + x3*(2.9192806270460925 - 1.019047619047619*el1) + x4*(-1.28235780892213 + 1.073639455782313*el1) + eps*x5*(-3.8466571723355227 + 0.8486467106683944*el1);
    /** (2,2) */
    rholm[1] = 1. + (-1.0238095238095237 + 0.6547619047619048*nu)*x + (-1.94208238851096 - 1.5601379440665155*nu + 0.4625614134542706*nu2)*x2 + x3*(12.736034731834051 - 2.902228713904598*nu - 1.9301558466099282*nu2 + 0.2715020968103451*nu3 - 4.076190476190476*el2) + x4*(-2.4172313935587004 + 4.173242630385488*el2) + x5*(-30.14143102836864 + 7.916297736025627*el2);
    
    flm[0] = gsl_pow_int(rholm[0], 2);
    flm[1] = gsl_pow_int(rholm[1], 2);
    
    
    /** l=3 ------------------------------------------------------------------
     *  (3,1) */
    rholm[2] = 1. + (-0.7222222222222222 - 0.2222222222222222*nu)*x + (0.014169472502805836 - 0.9455667789001122*nu - 0.46520763187429853*nu2)*x2 + x3*(1.9098284139598072 - 0.4126984126984127*el1) + x4*(0.5368150316615179 + 0.2980599647266314*el1) + eps*x5*(1.4497991763035063 - 0.0058477188106817735*el1);
    /** (3,2) */
    rholm[3] = 1. + (0.003703703703703704*(328. - 1115.*nu + 320.*nu2)*x)/(-1. + 3.*nu) + (6.235191420376606e-7*(-1.444528e6 + 8.050045e6*nu - 4.725605e6*nu2 - 2.033896e7*nu3 + 3.08564e6*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + x3*(6.220997955214429 - 1.6507936507936507*el2) + eps*x4*(-3.4527288879001268 + 2.005408583186361*el2);
    /** (3,3) */
    rholm[4] = 1. + (-1.1666666666666667 + 0.6666666666666666*nu)*x + (-1.6967171717171716 - 1.8797979797979798*nu + 0.45151515151515154*nu2)*x2 + x3*(14.10891386831863 - 3.7142857142857144*el3) + x4*(-6.723375314944128 + 4.333333333333333*el3) + eps*x5*(-29.568699895427518 + 6.302092352092352*el3);
    
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
    
    const int k22 = 1;//LM2K(2,2);
    const int k21 = 0;//LM2K(2,1);
    
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

vector<double> s_flm(double x, void *params){
    
    /*
     % Function EOB_flm(x). This function explicitly computes
     % the resummed flm in the general nu-dependent case
     %
     % Generated with EOB_rholm_complete.nb mathematica notebook
     %
     % USAGE: s=EOB_flm(x,nu,eps)
     %
     %        x   :: PN parameter
     %        nu  :: symmetric mass ratio. EMRL case is nu=0
     %        eps :: eps=1 switches on Iyer-Fujita terms; eps=0 switches them off.
     %
     %        s is the output structure of the form:
     */
    
    //vector<double> flm = EOBflm(x,nu);
    double nu         = (*(TEOBResumParams *)params).nu;
    double X1         = (*(TEOBResumParams *)params).X1;
    double X2         = (*(TEOBResumParams *)params).X2;
    double chi1       = (*(TEOBResumParams *)params).chi1;
    double chi2       = (*(TEOBResumParams *)params).chi2;
    bool   tidal_flag = (*(TEOBResumParams *)params).flags.tidal;
    
    double rho22S;
    
    
    double deltam  = X1 - X2;
    double chiS    = 0.5*(chi1 + chi2);
    double chiA    = 0.5*(chi1 - chi2);
    
    // velocity variables
    double v  = sqrt(x);
    double v2 = x;
    double v3 = v*v2;
    double v4 = v3*v;
    double v5 = v4*v;
    
    double nu2 = nu*nu;
    
    const double eps = 1.0;
    
    // Shorthands
    const double x2  = x*x;
    const double x3  = x*x2;
    const double x4  = x*x3;
    const double x5  = x*x4;
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
    
    //================
    // l=m=2 multipole
    //================
    
    // spin-orbit terms (even-parity)
    double cSO_lo    = -2./3.*(chiS*(1-nu)+chiA*deltam);
    double cSO_nlo   = (  (-34./21. + 49./18.*nu + 209./126.*nu2)*chiS + (-34./21. - 19./42.*nu)*deltam*chiA );
    
    //========================================
    // SPIN-SPIN contribution (even-parity):
    // put it to zero when tides are present
    //========================================
    double a1      = (*(TEOBResumParams *)params).a1;
    double a2      = (*(TEOBResumParams *)params).a2;
    double a0      = a1+a2;
    double C_Q1    = (*(TEOBResumParams *)params).C_Q1;
    double C_Q2    = (*(TEOBResumParams *)params).C_Q2;
    
    double cSS_lo;
    if (tidal_flag==true) {
        cSS_lo = 0.5*(C_Q1*a1*a1 + 2.*a1*a2 + C_Q2*a2*a2);
    } else {
        cSS_lo = 0.5*a0*a0; // spin-spin contribution to zero for BNS
    }
    
    //the spin-dependent part in taylor-expanded form (bad towards merger)
    rho22S = cSO_lo*v3 + cSS_lo*v4 + cSO_nlo*v5 ;//+ cSS_nlo*v6 + cSO_nnlo*v7;
    
    
    //===============
    // l>2 multipoles
    //===============
    double rho32S  = -4.*nu/(3.*(3.*nu-1.))*chiS*v;
    double rho44S  = -1./(15.*(1.-3.*nu))*((42.*nu2-41.*nu+10.)*chiS + (10.-39.*nu)*deltam*chiA)*v3;
    double rho42S  = -1./(15.*(1.-3.*nu))*((78.*nu2-59.*nu+10.)*chiS + (10.-21.*nu)*deltam*chiA)*v3;
    
    double sqrt_one_4nu = sqrt(1.-4.*nu);
    
    // Note that these are deltam*flm of Eq. (A15a)-(A15d) of Taracchini et al.
    // Typo in f31 corrected
    double f21S = -1.5*(deltam*chiS + chiA)*v+ v3*(  (61./12. + 79./84.*nu)*deltam*chiS + (61./12. + 131./84.*nu)*chiA );
    double f33S = -( deltam*chiS*(2. -  5./2.*nu) + chiA*(2. - 19./2.*nu) )*v3;
    double f31S = -( deltam*chiS*(2. - 13./2.*nu) + chiA*(2. - 11./2.*nu) )*v3;
    double f43S = -5.*nu/(2.*(2.*nu-1.))*(deltam*chiS - chiA)*v;
    double f41S =  f43S;
    
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
    flm[0] = (sqrt_one_4nu*flm[0] + f21S);
    
    
    // l=3 ------------------------------------------------------------------
    //(3,3)
    rholm[4] = 1. + (-1.1666666666666667 + 0.6666666666666666*nu)*x + (-1.6967171717171716 - 1.8797979797979798*nu + 0.45151515151515154*nu2)*x2 + x3*(14.10891386831863 - 3.7142857142857144*el3) + x4*(-6.723375314944128 + 4.333333333333333*el3) + eps*x5*(-29.568699895427518 + 6.302092352092352*el3);
    
    flm[4] = gsl_pow_int(rholm[4], 3);
    flm[4] = (sqrt_one_4nu*flm[4] + f33S);
    
    //(3,2)
    rholm[3] = 1. + (0.003703703703703704*(328. - 1115.*nu + 320.*nu2)*x)/(-1. + 3.*nu) + (6.235191420376606e-7*(-1.444528e6 + 8.050045e6*nu - 4.725605e6*nu2 - 2.033896e7*nu3 + 3.08564e6*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + x3*(6.220997955214429 - 1.6507936507936507*el2) + eps*x4*(-3.4527288879001268 + 2.005408583186361*el2);
    
    flm[3]=gsl_pow_int(rholm[3]+ rho32S, 3);
    
    //(3,1)
    rholm[2] = 1. + (-0.7222222222222222 - 0.2222222222222222*nu)*x + (0.014169472502805836 - 0.9455667789001122*nu - 0.46520763187429853*nu2)*x2 + x3*(1.9098284139598072 - 0.4126984126984127*el1) + x4*(0.5368150316615179 + 0.2980599647266314*el1) + eps*x5*(1.4497991763035063 - 0.0058477188106817735*el1);
    
    flm[2] = gsl_pow_int(rholm[2], 3);
    flm[2] = (sqrt_one_4nu*flm[2] + f31S);
    
    // l=4 ------------------------------------------------------------------
    //(4,4)
    rholm[8] = 1. + (0.0007575757575757576*(1614. - 5870.*nu + 2625.*nu2)*x)/(-1. + 3.*nu) + (3.1534122443213353e-9*(-5.11573572e8 + 2.338945704e9*nu - 3.13857376e8*nu2 - 6.733146e9*nu3 + 1.252563795e9*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + x3*(15.108111214795123 - 3.627128427128427*el4) + eps*x4*(-8.857121657199649 + 4.434988849534304*el4);
    
    flm[8] = gsl_pow_int(rholm[8] + rho44S, 4);
    
    //(4,3)
    rholm[7] = 1. + (0.005681818181818182*(222. - 547.*nu + 160.*nu2)*x)/(-1. + 2.*nu) - 0.9783218202252293*x2 + eps*(x3*(8.519456157072423 - 2.0402597402597404*el3) +      x4*(-5.353216984886716 + 2.5735094451003544*el3));
    
    flm[7] = gsl_pow_int(rholm[7], 4);
    flm[7] = (sqrt_one_4nu*flm[7] + f43S);
    
    //(4,2)
    rholm[6] = 1. + (0.0007575757575757576*(1146. - 3530.*nu + 285.*nu2)*x)/(-1. + 3.*nu) - (3.1534122443213353e-9*(1.14859044e8 - 2.95834536e8*nu - 1.204388696e9*nu2 + 3.04798116e9*nu3 + 3.79526805e8*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + 4.550378418934105e-12*x3*(8.48238724511e11 - 1.9927619712e11*el2) + eps*x4*(-0.6621921297263365 + 0.787251738160829*el2);
    
    flm[6]=gsl_pow_int(rholm[6] + rho42S, 4);
    
    //(4,1)
    rholm[5] = 1. + (0.001893939393939394*(602. - 1385.*nu + 288.*nu2)*x)/(-1. + 2.*nu) - 0.36778992787515513*x2 + x3*(0.6981550175535535 - 0.2266955266955267*el1) + eps*x4*(-0.7931524512893319 + 0.2584672482399755*el1);
    
    flm[5] = gsl_pow_int(rholm[5], 4);
    flm[5] = (sqrt_one_4nu*flm[5] + f41S);
    
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

double s_HorizonFlux(double x, double Heff, double jhat, double nu, double X1, double X2, double chi1, double chi2)
{
    
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
    
    /** Coefficients of the v^8=x^4 term */
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
