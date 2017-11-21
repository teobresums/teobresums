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

#include <gsl/gsl_math.h>

#include "TEOBResumS.h"

using namespace std;

/** Initial conditions calculation for non-spinning systems */
vector<double> initial(TEOBResumParams *params)
{
  /* This routine computes the initial configuration for EOB
     evolution in the nonspinning case. It determines 
     post-post-circular initial data as described in Sec.IIE
     of DNB, PRD 87, 084035 */
    double nu = (*params).nu;
    double r0 = (*params).r0;
    /**************************************/
    /* Tidal parameters */
    /**************************************/   
    double lambdaAl2 = (*(TEOBResumParams *)params).LambdaAl2;
    double lambdaAl3 = (*(TEOBResumParams *)params).LambdaAl3;
    double lambdaAl4 = (*(TEOBResumParams *)params).LambdaAl4;
    double lambdaBl2 = (*(TEOBResumParams *)params).LambdaBl2;
    double lambdaBl3 = (*(TEOBResumParams *)params).LambdaBl3;
    double lambdaBl4 = (*(TEOBResumParams *)params).LambdaBl4;
    double kappa2T   = (*(TEOBResumParams *)params).kappaTl2;
    double tidal_flag= (*(TEOBResumParams *)params).flags.tidal;
    /* Keplerian orbital frequency */
    double OmgKepler = 1./(sqrt(r0)*sqrt(r0)*sqrt(r0));
    
    vector<double> y_init(7);
    int N  = 6;
    const double dr = 1.e-8;
    
    vector<double> r(2*N), dA(2*N), j(2*N), j2(2*N), djdr(2*N); /** j:angular momentum */
    vector<double> E0(2*N), Omega_j(2*N);
    vector<double> Fphi(2*N), Ctmp(2*N), prstar(2*N), pr(2*N), pph(2*N), dprstardt(2*N);
    vector<double> metric(5);
    
    double r2, r3, A, B, d2A, j3;
    double z3 = 2.0*nu*(4.0-3.0*nu);
    double H0eff, H0, psi, r_omega, v_phi, jhat, x;
    
    for (int i=2*N; i--;)
    {
        
        r[i]  = r0+(i-N+1)*dr;
        r2 = r[i]*r[i];
        r3 = r[i]*r[i]*r[i];
        
        /** Compute metric  */
        metric = Metric(r[i],params,false);
        A    =  metric[0];
        dA[i]=  metric[1];
        B    =  metric[3];
        d2A  =  A5pnP15_dd(r[i],params)[0];
        
        /** Angular momentum for circular orbit: circular ID  */
        j2[i]   =  r3*dA[i]/(2.*A-r[i]*dA[i]);
        j[i]    =  sqrt(j2[i]);
        j3      =  j2[i]*j[i];
        djdr[i] = -j3/r3*( 2.0 - 3.0*A/(r[i]*dA[i]) - A*d2A/(dA[i]*dA[i]) );
        
        /** For circular orbit at r0=r(N)  */
        H0eff      = sqrt(A*(1.0 + j2[i]/r2));                     /** effective circular energy per unit mass E_0^eff  */
        E0[i]      = sqrt(1.0 + 2.0*nu*(H0eff - 1.0) );            /** real energy for circular orbit E_0               */
        H0         = E0[i]/nu;                                     /** H_0/nu: this is the usual hatH_EOB               */
        Omega_j[i] = A*j[i]/(nu*r2*H0*H0eff);                      /** Orbital frequency (from Hamilton's equation)     */
        psi        = 2.*(1.0 + 2.0*nu*(H0eff - 1.0))/(r2*dA[i]);   /** correction factor to the radius  */
        r_omega    = r[i]*cbrt(psi);                               /** EOB-corrected radius  */
        v_phi      = Omega_j[i]*r_omega;                           /** "corrected" azimuthal velocity such that Kepler's law is satisfied, r_omg^3 Omg_i^2 = 1  */
        
        x    = v_phi * v_phi;
        jhat = j[i]/(r_omega*v_phi); /** Newton-normalized angular momentum  */
        
        double prefact[] = {
            jhat,H0eff,
            H0eff,jhat,H0eff,
            jhat,H0eff,jhat,H0eff,
            H0eff,jhat,H0eff,jhat,H0eff,
            jhat,H0eff,jhat,H0eff,jhat,H0eff,
            H0eff,jhat,H0eff,jhat,H0eff,jhat,H0eff,
            jhat,H0eff,jhat,H0eff,jhat,H0eff,jhat,H0eff};
        
        Fphi[i] = flux(x, Omega_j[i], r_omega, E0[i], H0eff, jhat, r[i], 0., 0., prefact, params);
        
        /** Radial momentum conjugate to r*: post-circular ID  */
        Ctmp[i]   = sqrt(B/A)*nu*H0*H0eff;
        prstar[i] = Ctmp[i]*Fphi[i]/djdr[i];
        
        /** Radial momentum conjugate to r  */
        pr[i] = prstar[i]*sqrt(B/A);
        
    }
    
    dprstardt = s_D1(prstar, r, 2*N);
    
    for(int i=2*N; i--;)
    {
        dprstardt[i] *= Fphi[i]/djdr[i];
        pph[i] = j[i]*sqrt(1. + 2.*Ctmp[i]/dA[i]*dprstardt[i] - z3*gsl_pow_int(prstar[i],4)/j2[i]);
    }
    
    y_init[0] = r[N-1];
    y_init[1] = pph[N-1];
    y_init[2] = prstar[N-1];
    y_init[3] = pr[N-1];
    y_init[4] = j[N-1];
    y_init[5] = E0[N-1];
    y_init[6] = Omega_j[N-1];
    if (DEBUG)
    {
        printf("-----------------------------------\n");
        printf("Initial configuration:\n");
        printf("-----------------------------------\n");
        printf("r[0]         = %18.16f\n",y_init[0]);
        printf("pphi[0]      = %18.16f\n",y_init[1]);
        printf("j[0]         = %18.16f\n",y_init[4]);
        printf("pr[0]        = %18.16f\n",y_init[3]);
        printf("prstar[0]    = %18.16f\n",y_init[2]);
        printf("E[0]         = %18.16f\n",y_init[5]);
        printf("Omega[0]     = %18.16f\n",y_init[6]);
        printf("Omega-Kepler = %18.16f\n",OmgKepler);
    }
    if (tidal_flag==1)
    {
        if (DEBUG)
        {
            printf("lambda_l2[A] = %18.16f\n",lambdaAl2);
            printf("lambda_l3[A] = %18.16f\n",lambdaAl3);
            printf("lambda_l4[A] = %18.16f\n",lambdaAl4);
            printf("lambda_l2[B] = %18.16f\n",lambdaBl2);
            printf("lambda_l3[B] = %18.16f\n",lambdaBl3);
            printf("lambda_l4[B] = %18.16f\n",lambdaBl4);
            printf("kappa2T      = %18.16f\n",kappa2T);
        }
    }
    
    return y_init;
}

/** Initial conditions calculation for spinning systems */
vector<double> s_initial(TEOBResumParams *params){
    
    /*
     s_initial.
     This function computes the post-circular initial data
     in the presence of spin for initial relative separation
     r0. NOTE: this is just post-circular and not post-post
     circular as in the nonspinning case.
     Reference: TD, AN & SB, PRD87 (2013) 084035, 
     Eq. (51) and (52)
     Notation is as follows:
     r0       => relative separation
     pph0     => circular angular momentum
     pr       => post-circular radial momentum
     pprstar0 => post-circular r*-conjugate radial momentum
     
     Here the post-circular conditions are extended to the 
     spinning case and one follows TD & AN, PRD 90, 044018 (2014)

     It consists of two step:
     
     0. Compute j0            =>   circular ID, j!=0, pr =0
     1. From j0, compute pr*  =>   post circular ID, j!=0, pr !=0
     
     */
    
    double nu   = (*(TEOBResumParams *)params).nu;
    double r0   = (*(TEOBResumParams *)params).r0;
    double chi1 = (*(TEOBResumParams *)params).chi1;
    double chi2 = (*(TEOBResumParams *)params).chi2;
    double S1   = (*(TEOBResumParams *)params).S1;
    double S2   = (*(TEOBResumParams *)params).S2;
    double c3   = (*(TEOBResumParams *)params).cN3LO;
    double X1   = (*(TEOBResumParams *)params).X1;
    double X2   = (*(TEOBResumParams *)params).X2;
    double a1   = (*(TEOBResumParams *)params).a1;
    double a2   = (*(TEOBResumParams *)params).a2;
    double aK2  = (*(TEOBResumParams *)params).aK2;
    /**************************************/
    /* Tidal parameters for output */
    /**************************************/   
    double CQ1       = (*(TEOBResumParams *)params).C_Q1;
    double CQ2       = (*(TEOBResumParams *)params).C_Q2;
    double lambdaAl2 = (*(TEOBResumParams *)params).LambdaAl2;
    double lambdaAl3 = (*(TEOBResumParams *)params).LambdaAl3;
    double lambdaAl4 = (*(TEOBResumParams *)params).LambdaAl4;
    double lambdaBl2 = (*(TEOBResumParams *)params).LambdaBl2;
    double lambdaBl3 = (*(TEOBResumParams *)params).LambdaBl3;
    double lambdaBl4 = (*(TEOBResumParams *)params).LambdaBl4;
    double kappa2T   = (*(TEOBResumParams *)params).kappaTl2;
    double tidal_flag= (*(TEOBResumParams *)params).flags.tidal;
    
    //-----------------------------------------------------------------
    // Build  a small  grid (2*N points) around the initial position r0
    // This is used later to compute a spatial derivative numerically
    //-----------------------------------------------------------------
    
    vector<double> y_init(7);
    
    int N = 10;
    const double dr = 1.e-8;
    vector<double> r(2*N);
    
    // Angular momentum for circular orbit: circular ID
    vector<double> j2(2*N);
    vector<double> j(2*N);
    vector<double> djdr(2*N);
    
    // For circular orbit at r0=r(N)
    vector<double> E0(2*N);         // real Hamiltonian      H_0
    vector<double> Omega_j(2*N);    // Orbital frequency (from Hamilton's equation)        
    vector<double> Fphi(2*N);       // radial momentum
    vector<double> prstar(2*N);
    vector<double> pr(2*N);
    vector<double> pph(2*N);    
    vector<double> dprstardt(2*N);  // NOTE: Fr* here
    vector<double> rc(2*N);
    vector<double> drc(2*N);
    vector<double> d2rc(2*N);
    double rorb;
    double pphorb;
    vector<double> metric(5);
    vector<double> A(2*N);
    vector<double> dA(2*N);
    vector<double> B(2*N);
    vector<double> d2A(2*N);

    
    /* Dimensionful spin variables entering the EOB spin-orbit sector.
       Weighted averages of the dimensionful spins */
    double S  = S1 + S2;        // => in the extreme-mass-ratio-limit this becomes the spin of the BH
    double Ss = X2*a1 + X1*a2;  // => in the extreme-mass-ratio-limit this becomes the spin of the particle
    
    for (int i=2*N; i--;) {
        
        r[i] = r0+(i-N+1)*dr;
        
        // metric is here. The tidal parameters are within this routine
        metric = s_Metric(r[i], params,false);
        A[i]   = metric[0];
        B[i]   = metric[1];
        dA[i]  = metric[2];
        d2A[i] = metric[3];
        
        vector<double> rc_rad = s_get_rc(r[i],params);
	/** centrifugal radius (i.e., what is at the denominator of pphi^2 **/
        rc[i]  = rc_rad[0];
        drc[i] = rc_rad[1];
        
        //Compute minimum of Heff0 using bisection method
        rorb   = r[i];
        pphorb = rorb/sqrt(rorb-3.);
        pph[i] = s_bisec(pphorb,rorb,A[i],dA[i],rc[i],drc[i],aK2,S,Ss,params);
    }
    
    vector<double> dpph_dr = s_D1(pph,r,2*N); 
    
    for (int i=2*N; i--;) {
        
        double sqrtAbyB     = sqrt(A[i]/B[i]);
        
        double uc  = 1./rc[i];
        double uc2 = uc*uc;
        
        // circular angular momentum
        double pph2 = pph[i]*pph[i];
        
        // Orbital effective Hamiltonian
        double Horbeff0 = sqrt(A[i]*(1. + pph2*uc2));
        
        // Compute gyro-gravitomagnetic coupling functions
        vector<double> ggm0 = s_GS(r[i], rc[i], drc[i], aK2, 0., pph[i], nu, chi1, chi2, X1, X2, c3);
        double GS_0                   = ggm0[2];
        double GSs_0                  = ggm0[3];
        double dGS_dr_0               = ggm0[6];
        double dGSs_dr_0              = ggm0[7];
        double dGSs_dpph_0            = ggm0[9];
        double dGS_dprstarbyprstar_0  = ggm0[10];
        double dGSs_dprstarbyprstar_0 = ggm0[11];
        
        // Final effective Hamiltonian
        double Heff0 = (GS_0*S + GSs_0*Ss)*pph[i] + Horbeff0;
        
        // Real Hamiltonian: beware that this is NOT divided by nu
        double H0     = sqrt( 1. + 2.*nu*(Heff0 - 1.));
        double one_H0 = 1./H0;
        
        // get gyro-gravitomagnetic (derivative) functions
        double dHeff_dprstarbyprstar = pph[i]*(dGS_dprstarbyprstar_0*S + dGSs_dprstarbyprstar_0*Ss) + 1./Horbeff0;
        
        double C0 = sqrtAbyB*one_H0*dHeff_dprstarbyprstar;
        
        // orbital frequency for circular orbit
        double dHeff_dpph = GS_0*S + (GSs_0 + pph[i]*dGSs_dpph_0)*Ss + pph[i]*A[i]*uc2/Horbeff0;
        double Omg        = one_H0*dHeff_dpph;
        
        //================
        // computing F_phi
        //================
        
        // get the psi correctly with the spin contribution
        double Gtilde     =  GS_0*S     + GSs_0*Ss;
        double dGtilde_dr =  dGS_dr_0*S + dGSs_dr_0*Ss;
        double duc_dr     = -uc2*drc[i];
        double psic       = (duc_dr + dGtilde_dr*rc[i]*sqrt(A[i]/pph2 + A[i]*uc2)/A[i])/(-0.5*dA[i]);
        double r_omg      =  pow((pow(rc[i]*rc[i]*rc[i]*psic,-1./2)+Gtilde)*one_H0,-2./3.);
        double v_phi      =  r_omg*Omg;
        double x          =  v_phi*v_phi;
        double jhat       =  pph[i]/(r_omg*v_phi); // Newton-normalized angular momentum
        
        Fphi[i] = s_Flux(x, Omg, r_omg, H0, Heff0, jhat, r[i], 0., 0., params);
        
        //=================================================================
        // 1. The radial momentum conjugate to r*: post-circular correction
        //=================================================================
        
        prstar[i] = Fphi[i]/(dpph_dr[i]*C0);
        pr[i]     = prstar[i]* sqrt(B[i]/A[i]);
        
        j[i]       = pph[i];
        E0[i]      = H0;
        Omega_j[i] = Omg;
        
    }
    
    y_init[0] = r[N-1];
    y_init[1] = pph[N-1];
    y_init[2] = prstar[N-1];
    y_init[3] = pr[N-1];
    y_init[4] = j[N-1];
    y_init[5] = E0[N-1];
    y_init[6] = Omega_j[N-1];

    if (DEBUG)
    {
        printf("-----------------------------------\n");
        printf("Initial configuration:\n");
        printf("-----------------------------------\n");
        printf("r[0]         = %18.16f\n",y_init[0]);
        printf("pphi[0]      = %18.16f\n",y_init[1]);
        printf("j[0]         = %18.16f\n",y_init[4]);
        printf("pr[0]        = %18.16f\n",y_init[3]);
        printf("prstar[0]    = %18.16f\n",y_init[2]);
        printf("E[0]         = %18.16f\n",y_init[5]);
        printf("Omega[0]     = %18.16f\n",y_init[6]);
    }

    if (tidal_flag==1)
      {
          if (DEBUG)
          {
            printf("CQ1          = %18.16f\n",CQ1);
            printf("CQ2          = %18.16f\n",CQ2);
            printf("lambda_l2[A] = %18.16f\n",lambdaAl2);
            printf("lambda_l3[A] = %18.16f\n",lambdaAl3);
            printf("lambda_l4[A] = %18.16f\n",lambdaAl4);
            printf("lambda_l2[B] = %18.16f\n",lambdaBl2);
            printf("lambda_l3[B] = %18.16f\n",lambdaBl3);
            printf("lambda_l4[B] = %18.16f\n",lambdaBl4);
            printf("kappa2T      = %18.16f\n",kappa2T);
          }
      }
    
    return y_init;
}

