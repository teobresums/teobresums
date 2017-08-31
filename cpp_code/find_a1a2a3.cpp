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
#include <vector>
#include <stdio.h>
#include <fstream>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <limits>

#include "find_a1a2a3.h"
#include "input_struc.h"
#include "multipole_index.h"
#include "s_D1.h"
#include "dtnqc_fit.h"

typedef std::numeric_limits< double > dbl;

using namespace::std;

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

    double A_tmp, dA_tmp, omg_tmp, domg_tmp;
    
    double nu         = (*(input *)params).nu;
    double nu2 = nu*nu;
    double nu3 = nu*nu*nu;
    double X1       = (*(input *)params).X1;
    double X2       = (*(input *)params).X2;
    double X12      = X1 - X2;
    double chi1       = (*(input *)params).chi1;
    double chi2       = (*(input *)params).chi2;

    double aK         = (*(input *)params).aK;
    double aK2        = aK*aK;
    double aK3        = aK*aK*aK;
    double aK4        = aK*aK*aK*aK;
    double a12        = X1*chi1 - X2*chi2;

    long int t_length = T.size();
    double c_p1,     c_p2,     c_p3,   c_p4;
    double c_pdA1,   c_pdA2,   c_pdA3, c_pdA4;
    double c_pdomg1, c_pdomg2;
    double n0, d1;
    double a1_omg_tmp, a2_omg_tmp, b1_omg_tmp, b2_omg_tmp, a1_domg_tmp, a2_domg_tmp, b1_domg_tmp, b2_domg_tmp, a1_A_tmp , a2_A_tmp, b1_A_tmp, b2_A_tmp, a1_dA_tmp, a2_dA_tmp, b1_dA_tmp, b2_dA_tmp, omg_tmp_nu, omg_tmp_equal, domg_tmp_nu, domg_tmp_equal,  A_tmp_scale_nu, A_tmp_scale_equal, dA_tmp_scale_nu, dA_tmp_scale_equal ;
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
            /** In general divide by sqrt( (l+2) (l+1) l (l-1) ). Use the multipole structure to get the correct L. */
            A[k][j] = A[k][j]/sqrt(24.);
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
        A_tmp    =  pA[0]*aK4    + pA[1]*aK3   + pA[2]*aK2    + pA[4]*aK     + pA[4];

        pdA[0]   =  0.00000927;
        pdA[1]   = -0.00024550;
        pdA[2]   =  0.00012469;
        pdA[3]   =  0.00123845;
        pdA[4]   = -0.00195014;
        dA_tmp   =  pdA[0]*aK4   + pdA[1]*aK3   + pdA[2]*aK2   + pdA[4]*aK   + pdA[4];
        
        pomg[0]  =  0.00603482;
        pomg[1]  =  0.01604555;
        pomg[2]  =  0.02290799;
        pomg[3]  =  0.07084587;
        pomg[4]  =  0.38321834;
        omg_tmp  =  pomg[0]*aK4  + pomg[1]*aK3  + pomg[2]*aK2  + pomg[4]*aK  + pomg[4];

        pdomg[0] =  0.00024066;
        pdomg[1] =  0.00038123;
        pdomg[2] = -0.00049714;
        pdomg[3] =  0.00041219;
        pdomg[4] =  0.01190548;
        domg_tmp =  pdomg[0]*aK4 + pdomg[1]*aK3 + pdomg[2]*aK2 + pdomg[4]*aK + pdomg[4];
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
        a1_omg_tmp         =  0.205958;
        a2_omg_tmp         = -0.282734;
        b1_omg_tmp         =  0.186073;
        b2_omg_tmp         = -0.217723;
        omg_tmp_nu         =  0.6383186929*nu*nu + 0.2198527359*nu+ 0.2886403943;
        omg_tmp_equal      = ((a2_omg_tmp*X12*X12 + a1_omg_tmp*X12 - 0.1401748476)*aeff_omg + 1)/((b2_omg_tmp*X12*X12 + b1_omg_tmp*X12 - 0.3375083723)*aeff_omg + 1);
        omg_tmp            = omg_tmp_nu*omg_tmp_equal;
        
        a1_domg_tmp        =  0.0709177;
        a2_domg_tmp        = -0.0505505;
        b1_domg_tmp        =  0.033916;
        b2_domg_tmp        = -0.00755181;
        domg_tmp_nu        =  0.0449367831*nu*nu + 0.0097045815*nu + 0.0066911252;
        domg_tmp_equal     = (a2_domg_tmp*X12*X12 + a1_domg_tmp*X12 - 0.0277484292)*aeff_omg*aeff_omg + (b2_domg_tmp*X12*X12 + b1_domg_tmp*X12 + 0.0603634961)*aeff_omg + 1;
        domg_tmp           = domg_tmp_nu*domg_tmp_equal;
        
        a1_A_tmp           =  0.0905463;
        a2_A_tmp           =  0.0381341;
        b1_A_tmp           =  0.111952;
        b2_A_tmp           = -0.00790612;
        A_tmp_scale_nu     = -1.4938817908*nu3 +1.0576568105*nu2 - 0.0779048897*nu+0.2964517117;
        A_tmp_scale_equal  = ((a2_A_tmp*X12*X12 + a1_A_tmp*X12 - 0.2764889288)*aeff+1)/((b2_A_tmp*X12*X12 + b1_A_tmp*X12 -0.4706843028)*aeff+1);

        A_tmp              = A_tmp_scale_nu*A_tmp_scale_equal*(1-0.5*omg_tmp*aeff);
        
        a1_dA_tmp          =  0.00143545;
        a2_dA_tmp          = -0.00162301;
        b1_dA_tmp          =  0.00271927;
        b2_dA_tmp          = -0.00490688;
        dA_tmp_scale_nu    = -0.0017246790*nu-0.0046671920;
        dA_tmp_scale_equal = (a2_dA_tmp*X12*X12 + a1_dA_tmp*X12-0.0001583384)*aeff*aeff + (b2_dA_tmp*X12*X12 + b1_dA_tmp*X12+0.0037503520)*aeff;
        dA_tmp             = (dA_tmp_scale_nu  + dA_tmp_scale_equal)*omg_tmp;
        
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
    
    /** NQC corrections to AMPLITUDE (n1,n2,n3) and PHASE (n4,n5,n6)
     * NQC basis for (2,2) waveform : AMPLITUDE
     * note: n3 and n6 are not used
     */
    
    for (int j=t_length;j--;)
    {
        n1[j]  = pow(pr_star[j]/(r[j]*w[j]),2);  // [pr*/(r Omg)]^2
        n2[j]  = ddotr[j]/(r[j]*w[j]*w[j]);      // [ddot{r}/(r Omg^2)]
        n3[j]  = n1[j]*pr_star[j]*pr_star[j];    // [pr*/(r Omg)]^2 *(pr*)^2

        // NQC basis for (2,2) waveform: PHASE

        n4[j]  = pr_star[j]/(r[j]*w[j]);           //  pr*/(r Omg)
        n5[j]  = n4[j]*pow(r[j]*w[j],2);           // (pr*)*(r Omg)
        n6[j]  = n5[j]*pow(pr_star[j],2);          // (pr*^3)*(r Omg)
    }
                                  
    /** Take the needed derivatives for the phase */
    vector<double>  d_n4 = s_D1(n4,T,t_length-1);
    vector<double>  d_n5 = s_D1(n5,T,t_length-1);
    vector<double> d2_n4 = s_D1(d_n4,T,t_length-1);
    vector<double> d2_n5 = s_D1(d_n5,T,t_length-1);

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
    
    if (chi1 >= 0.8498)
    {
        /* Interpolating fit for Deltat_NQC. See Eq.(21) of arXiv:1506.08457 */
        DeltaT_nqc = dtnqc_fit(chi1,0.8498);
    }
    else
    {
        DeltaT_nqc = 1.;
    }

    double tNQC = tOmgOrb_pk - DeltaT_nqc;

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
        m21[k]   = s_D1(m11[k],T,t_length-1);
        m22[k]   = s_D1(m12[k],T,t_length-1);
        
        p1tmp[k] = A[k];
        p2tmp[k] = s_D1(p1tmp[k],T,t_length-1);
    }


        
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
        ai[k][2] = 0.;

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
        bi[k][2] =  0.;
    }

    for (int k=35;k--;)
    {
        for (int j=0; j<t_length-1;j++)
        {
            o[k][j].dat[0] = 1. + ai[k][0]*n1[j] + ai[k][1]*n2[j] + ai[k][2]*n3[j];
            o[k][j].dat[1] = 0. + bi[k][0]*n4[j] + bi[k][1]*n5[j] + bi[k][2]*n6[j];
        }
    }

    return o;

}
