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
#include "input_struc.h"
#include "initial.h"
#include "FDdrvt.h"
#include "Metric.h"
#include "flux.h"

double* initial(input *params)
{

    double nu = (*params).nu;
    double r0 = (*params).r0;

    static double y_init[7];
    const int N = 6;
    const double dr = 1.e-8;
    
    double r[2*N], dA[2*N], j[2*N], j2[2*N], djdr[2*N]; /** j:angular momentum */
    double E0[2*N], Omega_j[2*N];
    double Fphi[2*N], Ctmp[2*N], prstar[2*N], pr[2*N], pph[2*N], dprstardt[2*N];
    double metric[5];

    double r2, r3, u, A, B, d2A, j3;
    double z3 = 2.0*nu*(4.0-3.0*nu);
    double H0eff, H0, psi, r_omega, v_phi, jhat, x;
    double c1 = (-227./140.*nu + 1957./1680.);
    double c2 = (753./560.*nu*nu + 165703./70560.*nu - 25672541./5080320.);
    double Frhat, Frstar;

    for (int i=2*N; i--;)
    {
        
        r[i]  = r0+(i-N+1)*dr;
        r2 = r[i]*r[i];
        r3 = r[i]*r[i]*r[i];
        u  = 1./r[i];
        
        /** Compute metric  */
        metric = Metric(r[i], params, false);
        A      = metric[0];
        dA[i]  = metric[1];
        B      = metric[3];
        d2A    = A5pnP15_dd(r[i],params)[0];
        
        /** Angular momentum for circular orbit: circular ID  */
        j2[i]   =  r3*dA[i]/(2.*A-r[i]*dA[i]);
        j[i]    =  sqrt(j2[i]);
        j3      =  j2[i]*j[i];
        djdr[i] = -j3/r3*( 2.0 - 3.0*A/(r[i]*dA[i]) - A*d2A/(dA[i]*dA[i]) );
        
        /** For circular orbit at r0=r(N)  */
        H0eff      = sqrt(A*(1.0 + j2[i]/r2));                     /** ffective Hamiltonian H_0^eff  */
        E0[i]      = sqrt(1.0 + 2.0*nu*(H0eff - 1.0) );            /** real Hamiltonian      H_0  */
        H0         = E0[i]/nu;                                     /** H_0/nu  */
        Omega_j[i] = A*j[i]/(nu*r2*H0*H0eff);                      /** Orbital frequency (from Hamilton's equation)  */
        psi        = 2.*(1.0 + 2.0*nu*(H0eff - 1.0))/(r2*dA[i]);   /** correction factor to the radius  */
        r_omega    = r[i]*pow(psi,1.0/3.0);                        /** EOB-corrected radius  */
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
        
        /** Angular momentum again: post-post-circular ID  */
        Frhat  = 1. + c1*u + c2*u*u;
        Frstar = -5./3.*Fphi[i]*prstar[i]/j[i] * Frhat;
        
        /** Radial momentum conjugate to r  */
        pr[i] = prstar[i]*sqrt(B/A);
        
    }
    
    dprstardt = FDdrvt(prstar,r,4,2*N);

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
    
    return y_init;
}

