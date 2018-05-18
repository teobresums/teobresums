/**
 * This file is part of TEOBResumS
 *
 * Copyright (C) 2017-2018 See AUTHORS file
 *
 * TEOBResumS is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * TEOBResumS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/.       
 *
 */

#include "TEOBResumS.h"

/** post-adiabatic dynamics */
int eob_dyn_Npostadiabatic(Dynamics *dyn, double r0)
{
    /* Unpack values */
    const double nu    = dyn->nu;
    const double S     = dyn->S;
    const double Sstar = dyn->Sstar;
    const double chi1  = dyn->chi1;
    const double chi2  = dyn->chi2;
    const double X1    = dyn->X1;
    const double X2    = dyn->X2;
    const double c3    = dyn->cN3LO;
    const double aK2   = dyn->aK2;
    const double a1    = dyn->a1;
    const double a2    = dyn->a2;
    const double C_Q1  = dyn->C_Q1;
    const double C_Q2  = dyn->C_Q2;
    const int usetidal = dyn->use_tidal;
    const int usespins = dyn->use_spins;

    const int size = dyn->size;
    const int Npa = par_get_i("postadiabatic_dynamics_N");

    /* Mem for quantities to be used. FIXME: Can we avoid allocating the memory?
     Can we use already-defined C-structures here? */
    double *A_vec                  = (double*)malloc(size * sizeof (double));
    double *dA_vec                 = (double*)malloc(size * sizeof (double));
    double *d2A_vec                = (double*)malloc(size * sizeof (double));
    double *B_vec                  = (double*)malloc(size * sizeof (double));
    double *dB_vec                 = (double*)malloc(size * sizeof (double));
    double *sqrtAbyB_vec           = (double*)malloc(size * sizeof (double));
    double *rc_vec                 = (double*)malloc(size * sizeof (double));
    double *drc_dr_vec             = (double*)malloc(size * sizeof (double));
    double *d2rc_dr_vec            = (double*)malloc(size * sizeof (double));
    double *uc2_vec                = (double*)malloc(size * sizeof (double));
    double *duc_dr_vec             = (double*)malloc(size * sizeof (double));
    double *dAuc2_dr_vec           = (double*)malloc(size * sizeof (double));
    double *G_vec                  = (double*)malloc(size * sizeof (double));
    double *dG_dr_vec              = (double*)malloc(size * sizeof (double));
    double *dG_dprstar_vec         = (double*)malloc(size * sizeof (double));
    double *dG_dprstarbyprstar_vec = (double*)malloc(size * sizeof (double));
    double *dpphi_dr_vec           = (double*)malloc(size * sizeof (double));
    double *dprstar_dr_vec         = (double*)malloc(size * sizeof (double));
    double *dHdPrs                 = (double*)malloc(size * sizeof (double));
    double *dHdPphi_o_dHdPrs       = (double*)malloc(size * sizeof (double));
    double ggm[14];

    double a_coeff, b_coeff, c_coeff, Delta, sol_p, sol_m, j02, uc;

    /** Compute radius of inflection point of Pr* */
    //TODO implementing a robust stopping condition (which)
    double ri = 2.;

    /** Build a uniform grid and compute circular dynamics */
    const double dr = (r0 - ri)/(size-1);
    for (int i = 0; i < size; i++)
    {
        //FIXME dyn->data[EOB_RAD] forse non e' quello da usare, in altri punti usa dyn->r. Io non voglio che dyn->data cambi, quindi non mi piace passargli il puntatore.
        dyn->r = r0 - i*dr;
        /** Pre-compute metric functions and centrifugal radius */
        if(usespins)
        {
            eob_metric_s(dyn->r,dyn, &A_vec[i], &B_vec[i], &dA_vec[i], &d2A_vec[i], &dB_vec[i]);
            eob_dyn_s_get_rc(dyn->r, nu, a1, a2, aK2, C_Q1, C_Q2, usetidal, &rc_vec[i], &drc_dr_vec[i], &d2rc_dr_vec[i]);
            eob_dyn_s_GS(dyn->r, rc_vec[i], drc_dr_vec[i], aK2, 0.0, 0.0, nu, chi1, chi2, X1, X2, c3, ggm);
            
            G_vec[i]                  = ggm[2] *S+ggm[3] *Sstar;    // Tilde G
            dG_dr_vec[i]              = ggm[6] *S+ggm[7] *Sstar;
            dG_dprstar_vec[i]         = ggm[4] *S+ggm[5] *Sstar;
            dG_dprstarbyprstar_vec[i] = ggm[10]*S+ggm[11]*Sstar;
        }
        else
        {
            eob_metric(dyn->r ,dyn, &A_vec[i], &B_vec[i], &dA_vec[i], &d2A_vec[i], &dB_vec[i]);
            rc_vec[i]                 = dyn->r ;
            drc_dr_vec[i]             = 1;
            d2rc_dr_vec[i]            = 0;
            G_vec[i]                  = 0.0;    // Tilde G
            dG_dr_vec[i]              = 0.0;
            dG_dprstar_vec[i]         = 0.0;
            dG_dprstarbyprstar_vec[i] = 0.0;
        }
        sqrtAbyB_vec[i] = sqrt(A_vec[i]/B_vec[i]);
        
        /** Auxiliary functions */
        uc              = 1./rc_vec[i];
        uc2_vec[i]      = uc*uc;
        duc_dr_vec[i]   = -uc2_vec[i]*drc_dr_vec[i];
        dAuc2_dr_vec[i] = uc2_vec[i]*(dA_vec[i]-2*A_vec[i]*uc*drc_dr_vec[i]);
        
        /*===================#
         Circular dynamics #
         ===================*/

        // Computing the circular angular momentum by solving eq. (A15) of TEOBResum paper (which is equivalent to solve eq.(4)=0 of arXiv:1805.03891). The procedure to choose the physical solution of the quadratic equation is effective but not understood.  FIXME

        a_coeff    = dAuc2_dr_vec[i]*dAuc2_dr_vec[i]   - 4*A_vec[i]*uc2_vec[i]*dG_dr_vec[i]*dG_dr_vec[i];  // First coefficient of the quadratic equation
        b_coeff    = 2*dA_vec[i]*dAuc2_dr_vec[i] - 4*A_vec[i]*dG_dr_vec[i]*dG_dr_vec[i];    // Second coefficient of the quadratic equation
        c_coeff    = dA_vec[i]*dA_vec[i];                             // Third coefficient of the quadratic equation
        Delta      = b_coeff*b_coeff - 4*a_coeff*c_coeff ;            // Delta of the quadratic equation
        sol_p      = (-b_coeff + sqrt(Delta))/(2*a_coeff);            // Plus  solution of the quadratic equation
        sol_m      = (-b_coeff - sqrt(Delta))/(2*a_coeff);            // Minus solution of the quadratic equation
        
        // Effective prescription: If the Tilde G function is negative, take the positive solution and vice versa.
        if (G_vec[i] < 0)
        {
            j02 = sol_p;
        }
        else
        {
            j02 = sol_m;
        }
        // Define momenta in the circular orbit approximation
        dyn->pphi                = sqrt(j02);
        dyn->data[EOB_PRSTAR][i] = 0.0;
        dprstar_dr_vec[i]        = 0.0;
        
        // Circular Hamiltonians, ref: arXiv: 1406.6913
        if(usespins)
        {

            eob_ham_s(nu,dyn->r,rc_vec[i],drc_dr_vec[i],dyn->pphi,0.0,S,Sstar,chi1,chi2,X1,X2,aK2,c3,A_vec[i],dA_vec[i],
                      &H,             real EOB Hamiltonian divided by mu=m1m2/(m1+m2)
                      &Heff,          effective EOB Hamiltonian (divided by mu)
               &Heff_orb,
               &dHeff_dr,      drvt Heff,r
               &dHeff_dprstar, drvt Heff,prstar
               &dHeff_dpphi,    drvt Heff,pphi
               &d2Heff_dprstar20
               )

            Horbeff[0,:] = sqrt(A*(1 + pphi[0,:]**2*uc2))
            Heff[0,:]    = G[0,:]*pphi[0,:] + Horbeff[0,:]
            H[0,:]       = sqrt(1 + 2*nu*(Heff[0,:] - 1))  //Convention: H = H_EOB * nu
            one_H[0,:]   = 1./H[0,:]

            // Circular orbital frequency
            dHeff_dpph   = G[0,:] + pphi[0,:]*A*uc2./Horbeff[0,:]
            Omg[0,:]     = one_H[0,:]*dHeff_dpph
        }
        else
        {
            //NON spinning hamiltonian
        }
        dyn->data[EOB_RAD][i]    = dyn->r;
        dyn->data[EOB_PPHI][i]    = dyn->pphi;
    }
    // END r-GRID FOR
    
    // Initialize derivatives. FIXME check which derivative to use
//    D0(dyn->data[EOB_PPHI], dr, 12, dphidr_vec);
    
    /** Compute Pphi and Pr* by iteration */
    for (int n = 0; n < Npa; n++)
    {
        for (int i = 0; i < size; i++)
        {
          
          // TODO PA iteration
          // maximize use of routines already coded e.g. flux routines, hamiltonian, etc
          // check code in *Dynamics.c and *InitialConditions.c
          // store final result into dyn->data[<variable_index>][i]
          // note Dynamics structure has 'storage' arrays dyn->data and point-wise vars,
          /*
          dyn->data[EOB_PHI][i]    = dyn->phi;
          dyn->data[EOB_PPHI][i]   = dyn->pphi;
          dyn->data[EOB_MOMG][i]   = dyn->Omg;
          dyn->data[EOB_DDOTR][i]  = dyn->ddotr;
          dyn->data[EOB_PRSTAR][i] = dyn->prstar;
          dyn->data[EOB_OMGORB][i] = dyn->Omg_orb;
          */
          
        }
    }
    /** Compute time */
    cumint3(dHdPrs, dyn->data[EOB_RAD], size, dyn->time);
    /** FIXME dHdPrs^-1*/

    /* Set last value for evolution */
    dyn->t = dyn->time[size-1];

    /** Compute orbital phase */
    cumint3(dHdPphi_o_dHdPrs, dyn->data[EOB_RAD], size, dyn->data[EOB_PHI]);

    /** Print su file */
    /* Free memory */
    free(dHdPrs);
    free(dHdPphi_o_dHdPrs);
    printf("\n\nDAJECHEGIRO\n\n");
    return OK;
}
