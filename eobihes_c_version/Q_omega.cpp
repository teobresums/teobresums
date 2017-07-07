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
#include <iostream>
#include <fstream>
#include "Q_omega.h"

using namespace::std;

//NOTE: routine not in use
int Q_omega(vector<double> t, vector<double> f, std::ofstream& Qomg){
    
    double dt   = t[1]-t[0];
    double oodt = 1./dt;
    
    double d1f;
    
    const double c = 1./12.;
    
    for (long int i=0; i<=t.size()-3; i++) {
        switch (i) {
            case 0:
                d1f = c*(-25.*f[i]+48.*f[i+1]-36.*f[i+2]+16.*f[i+3]-3.*f[i+4])*oodt;
                break;
            case 1:
                d1f = c*(-3.*f[i-1]-10.*f[i]+18.*f[i+1]-6.*f[i+2]+f[i+3])*oodt;
                break;
            default:
                d1f = c*(8.*(f[i+1]-f[i-1]) - f[i+2] + f[i-2])*oodt;
                break;
        }
        double Qomega = (f[i]*f[i])/d1f;
        Qomg << t[i] << "\t"<< f[i] << "\t" << d1f << "\t" << Qomega << endl;
    }
    
    
    return 0;
    
}

