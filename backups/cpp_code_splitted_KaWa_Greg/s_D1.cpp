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
#include "s_D1.h"

using namespace::std;

vector<double> s_D1(vector<double> f, vector<double> x, int Nmax){
/* Computes the first derivative of the function. Centered but at the edges. USAGE: df = EOB_D1(f,x) */
int Nmin = 0;
    
    vector<double> df(Nmax+1);
    for(int i=2;i<=Nmax-2;i++){
        df[i] = 1./3.*(8.*f[1+i] - f[2+i] - 8.*f[i-1] + f[i-2])/(x[2+i]-x[i-2]);
    }

    // 4th order boundaries
    df[0] = (-24./17.*f[Nmin] + 59./34.*f[Nmin+1] - 4./17.*f[Nmin+2] - 3./34.*f[Nmin+3])/(x[Nmin+1]-x[Nmin]);
    df[1] = (-1./2.*f[Nmin] + 1./2.*f[Nmin+2] )/(x[Nmin+1]-x[Nmin]);

    df[Nmax]   = -(-24./17.*f[Nmax] + 59./34.*f[Nmax-1] - 4./17.*f[Nmax-2] - 3./34.*f[Nmax-3])/(x[Nmax]-x[Nmax-1]);
    df[Nmax-1] = -(-1./2.*f[Nmax] + 1./2.*f[Nmax-2])/(x[Nmax]-x[Nmax-1]);

    return df;
}

