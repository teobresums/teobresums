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
#include "get_Omg_orb.h"

using namespace::std;

vector<double> get_Omg_orb(vector<double> r, vector<double> pph, vector<double> pr_star, vector<double> A, vector<double> B, void *params)
{

    double nu   = (*(input *)params).nu;
    double aK2  = (*(input *)params).aK2;
    double S1   = (*(input *)params).S1;
    double S2   = (*(input *)params).S2;
    double a1   = (*(input *)params).a1;
    double a2   = (*(input *)params).a2;
    double X1   = (*(input *)params).X1;
    double X2   = (*(input *)params).X2;
    double chi1 = (*(input *)params).chi1;
    double chi2 = (*(input *)params).chi2;
    double c3   = (*(input *)params).cN3LO;
    
    double S     = S1 + S2;
    double Sstar = X2*a1 + X1*a2;
    double z3    = 2.*nu*(4.-3.*nu);
    double rc;
    double drc_dr;
    double uc    = 1./rc;
    double uc2   = uc*uc;
    
    double prstar2;
    double prstar4;
    double pphi2;
    double Horbeff;
    double Heff;
    double H;
    double one_H
    double GS;
    double GSs;
    
    long int r_length = r.size();
    vector<double> Omg_orb(r_length);
    vector<double> rc_vec;
    vector<double> ggm;
    
    for (int i=r_length; i--;)
    {

        rc_vec = s_get_rc(r[i],aK2,params);
        rc     = rc_vec[0];
        drc_dr = rc_vec[1];

        ggm = s_GS(r[i],rc,drc_dr,aK2,pr_star[i],pph[i],nu,chi1,chi2,X1,X2,c3);
        GS  = ggm[2];
        GSs = ggm[3];
        
        prstar2 = pr_star[i]*pr_star[i];
        prstar4 = prstar2*prstar2;
        pphi2   = pph[i]*pph[i];
        Horbeff = sqrt(prstar2 + A[i]*(1. + pphi2*uc2 +  z3*prstar4*uc2));
        Heff    = Horbeff + (GS*S + GSs*Sstar)*pph[i];
        H       = sqrt( 1. + 2.*nu*(Heff - 1.) );
        one_H   = 1./H;

        Omg_orb[i] = one_H * pph[i] * A[i] * uc2/Horbeff;
    }

    return Omg_orb;
    
}
