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

gsl_complex** find_a1a2a3(double* T, const int size_T, double* r, double* w, double* pph, double* pr_star, double** hlm_phase, double* Omg_orb, double** A, double* ddotr, void *params)
{

    /* determine NQC parameters */
    double nu    = (*(input *)params).nu;
    double chi1  = (*(input *)params).chi1;
    double aKerr = (*(input *)params).aK;
    double A_tmp, dA_tmp, omg_tmp, domg_tmp;

    const long int t_length = size_T;
    
    double pA[5];
    double pdA[5];
    double pomg[5];
    double pdomg[5];
    double p1v[2], p2v[2], p3v[2];
    double pdA1v[2], pdA2v[2], pdA3v[2];
    double pn0[2], pd1[2];
    double ppdomg1[2], ppdomg2[2];
    double * omg[35];
    double *domg[35];
    double n1[t_length], n2[t_length], n3[t_length], n4[t_length], n5[t_length], n6[t_length];
    double max_A[35], max_dA[35], d2max[35], d3max[35], max_omg[35], max_domg[35], maxd2omg[35], DeltaT[35];
    //double ai[3];
    //double bi[3];
    double ai[35][3];
    double bi[35][3];
    double P[2];
    double M[4];
    // double m11[t_length], m12[t_length], m13[t_length];
    double *  m11[35];
    double *  m12[35];
    double *  m13[35];
    double *  m21[35];
    double *  m22[35];
    double *p1tmp[35];
    double *p2tmp[35];

    static gsl_complex o[35][t_length];

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
        domg[k] = s_D1(omg[k], T, t_length-1);
    }

    /**  Case 'NQC_fit_hybrid' */
    if (nu==0.25)
    {
        
        /** Amplitude */
        
        pA[0] =  0.00178195;
        pA[1] =  0.00435589;
        pA[2] =  0.00344489;
        pA[3] = -0.00076165;
        pA[4] =  0.31973334;
        A_tmp = pA[0]*pow(aKerr,4) + pA[1]*pow(aKerr,3) + pA[2]*pow(aKerr,2) + pA[3]*aKerr + pA[4];
        
        /** Derivative of amplitude */
        
        pdA[0] =  0.00000927;
        pdA[1] = -0.00024550;
        pdA[2] =  0.00012469;
        pdA[3] =  0.00123845;
        pdA[4] = -0.00195014;
        dA_tmp = pdA[0]*pow(aKerr,4) + pdA[1]*pow(aKerr,3) + pdA[2]*pow(aKerr,2) + pdA[3]*aKerr + pdA[4];
        
        /** Frequency */
        
        pomg[0] = 0.00603482;
        pomg[1] = 0.01604555;
        pomg[2] = 0.02290799;
        pomg[3] = 0.07084587;
        pomg[4] = 0.38321834;
        omg_tmp = pomg[0]*pow(aKerr,4) + pomg[1]*pow(aKerr,3) + pomg[2]*pow(aKerr,2) + pomg[3]*aKerr + pomg[4];
        
        /** Derivative of the frequency */
        
        pdomg[0] = 0.00024066;
        pdomg[1] = 0.00038123;
        pdomg[2] = -0.00049714;
        pdomg[3] = 0.00041219;
        pdomg[4] = 0.01190548;
        domg_tmp = pdomg[0]*pow(aKerr,4) + pdomg[1]*pow(aKerr,3) + pdomg[2]*pow(aKerr,2) + pdomg[3]*aKerr + pdomg[4];
        
    }
    else
    {
        /** Amplitude */
        
        p1v[0]    =  0.05385059;
        p1v[1]    = -0.00890942;
        p2v[0]    = -0.07942102;
        p2v[1]    =  0.02152423;
        p3v[0]    =  0.14805262;
        p3v[1]    =  0.28210487;
        double p1 =  p1v[0]*nu + p1v[1];
        double p2 =  p2v[0]*nu + p2v[1];
        double p3 =  p3v[0]*nu + p3v[1];
        A_tmp     =  p1*pow(aKerr,2) + p2*aKerr + p3;
        
        /** Derivative of amplitude */
        
        pdA1v[0]    =  0.00248472;
        pdA1v[1]    = -0.00033422;
        pdA2v[0]    =  0.00105298;
        pdA2v[1]    =  0.00085160;
        pdA3v[0]    = -0.00339257;
        pdA3v[1]    = -0.00110932;
        double pdA1 =  pdA1v[0]*nu + pdA1v[1];
        double pdA2 =  pdA2v[0]*nu + pdA2v[1];
        double pdA3 =  pdA3v[0]*nu + pdA3v[1];
        dA_tmp      =  pdA1*pow(aKerr,2) + pdA2*aKerr + pdA3;
        
        /** Frequency */
        pn0[0]    =  0.45584139;
        pn0[1]    =  0.27315247;
        pd1[0]    =  0.75276414;
        pd1[1]    = -0.40081625;
        double n0 =  pn0[0]*nu + pn0[1];
        double d1 =  pd1[0]*nu + pd1[1];
        omg_tmp   =  n0/(1. + d1*aKerr);

        /** Derivative of the frequency */
        
        ppdomg1[0]    = -0.00177362;
        ppdomg1[1]    =  0.00123900;
        ppdomg2[0]    =  0.02424739;
        ppdomg2[1]    =  0.00566504;
        double pdomg1 =  ppdomg1[0]*nu + ppdomg1[1];
        double pdomg2 =  ppdomg2[0]*nu + ppdomg2[1];
        domg_tmp      =  pdomg1*aKerr + pdomg2;
        
    }

    for (int k=35;k--;)
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
    double * d_n4 = s_D1(n4,T,t_length-1);
    double * d_n5 = s_D1(n5,T,t_length-1);
    double *d2_n4 = s_D1(d_n4,T,t_length-1);
    double *d2_n5 = s_D1(d_n5,T,t_length-1);

    int Omgmax_index = 0;
    double Omg_max = Omg_orb[0];
    int i = 1;
    
    while (Omg_orb[i] > Omg_max)
    {
        Omg_max = Omg_orb[i];
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

    i = 0;
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
        m21[k] = s_D1(m11[k],T,t_length-1);
        m22[k] = s_D1(m12[k],T,t_length-1);
        
        p1tmp[k] = A[k];
        p2tmp[k] = s_D1(p1tmp[k],T,t_length-1);
    }


        
    double detM = 1.;
    for (int k=35;k--;)
    {
        //int k=1;
        
        /** Computation of ai coefficients */
        P[0] = max_A[k]  - p1tmp[k][jmax];
        P[1] = max_dA[k] - p2tmp[k][jmax];
        
        M[0] = m11[k][jmax];
        M[1] = m12[k][jmax];
        M[2] = m21[k][jmax];
        M[3] = m22[k][jmax];
        
        detM     = M[0]*M[3]-M[1]*M[2];
        ai[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
        ai[k][1] = (M[0]*P[1] - M[2]*P[0])/detM;
        ai[k][2] = 0.;

        /** Computation of bi coefficients */
        P[0] = omg[k][jmax]   - max_omg[k];
        P[1] = domg[k][jmax]  - max_domg[k];
        
        M[0] = d_n4[jmax];
        M[1] = d_n5[jmax];
        M[2] = d2_n4[jmax];
        M[3] = d2_n5[jmax];

        detM = M[0]*M[3]-M[1]*M[2];
        bi[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
        bi[k][1] = (M[0]*P[1] - M[2]*P[0])/detM;
        bi[k][2] = 0.;
    }

    for (int k=35;k--;)
    {
        for (int j=0; j<t_length-1;j++)
        {
            o[k][j].dat[0] = 1. + ai[k][0]*n1[j] + ai[k][1]*n2[j] + ai[k][2]*n3[j];
            o[k][j].dat[1] = bi[k][0]*n4[j] + bi[k][1]*n5[j] + bi[k][2]*n6[j];
        }
    }

    //printf("%s %.16e %.16e %.16e %.16e \n","ai and bi for k=1",ai[1][0],ai[1][1],bi[1][0],bi[1][1]);
    return o;

}
