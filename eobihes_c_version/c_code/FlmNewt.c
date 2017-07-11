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

#include <stdbool.h>
#include "FlmNewt.h"
#include "input_struc.h"

//To do: pre-calculate the sp * x products

double* FlmNewt(const double x, void *params)
{

    double nu       = (*(input *)params).nu;
    bool tidal_flag = (*(input *)params).tidal;
    bool spin_flag  = (*(input *)params).spin;
    
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

    double sp2;
    double sp4;
    const double sp3 = (1.-3.*nu)*(1.-3.*nu);
    const double sp5 = (1.-5.*nu+5.*nu2)*(1.-5.*nu+5.*nu2);
    const double sp6 = (1-4*nu)*(3*nu2-4*nu +1)*(3*nu2-4*nu +1);
    const double sp7 = (1 - 7*nu + 14*nu2 - 7*nu3)*(1 - 7*nu + 14*nu2 - 7*nu3);
    const double sp8 = (1 - 4*nu)*(1 - 6*nu + 10*nu2 - 4*nu3)*(1 - 6*nu + 10*nu2 - 4*nu3);
    const double sp9 = (1 - 4*nu)*(1 - 4*nu + 3*nu2)*(1 - 4*nu + 3*nu2);
    
    static double Nlm[35];

    if (spin_flag==true)
    {
        sp2 = 1.;
        sp4 = (2*nu-1)*(2*nu-1);
    }
    else if (spin_flag==false)
    {   /** Nonspinning case*/
        if (tidal_flag==false)
        {
            sp2 = 1.-4.*nu;
            sp4 = (1.-4.*nu)*(1.-2.*nu)*(1.-2.*nu);
        }
        else if (tidal_flag==true)
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

