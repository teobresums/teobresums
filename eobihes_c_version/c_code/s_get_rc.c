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

#include <math.h>
#include <stdbool.h>

#include "input_struc.h"
#include "s_get_rc.h"


double* s_get_rc(double r, void *params)
{
    
    
    double nu         = (*(input *)params).nu;
    double at1        = (*(input *)params).a1;
    double at2        = (*(input *)params).a2;
    bool   tidal_flag = (*(input *)params).tidal;
    double aK2        = (*(input *)params).aK2;
    
    double C_Q1       = (*(input *)params).C_Q1;
    double C_Q2       = (*(input *)params).C_Q2;
    
    double rc, drc_dr, d2rc_dr2;
    
    double u   = 1./r;
    double u2  = u*u;
    double u3  = u*u2;
    double r2  = r*r;
    
    
    if (tidal_flag==true)
    {
        
         /* inclusion of LO spin-square coupling. The S1*S1 term coincides with the BBH one, no effect of structure.
         The self-spin couplings, S1*S1 and S2*S2 get a EOS-dependent coefficient, CQ, that describe the quadrupole
         deformation due to spin. Notation of Levi-Steinhoff, JCAP 1412 (2014), no.12, 003. Notation analogous to
         the parameter a of Poisson, PRD 57, (1998) 5287-5290 or C_ES^2 in Porto & Rothstein, PRD 78 (2008), 044013
        
         NS quadrupoles due to rotation
         These are parameters that should be specified in the parameter file

         Inclusion of LO spin-square coupling. The S1*S1 term coincides with the BBH one, no effect of structure.
         The self-spin couplings, S1*S1 and S2*S2 get a EOS-dependent coefficient, CQ, that describe the quadrupole
         deformation due to spin. Notation of Levi-Steinhoff, JCAP 1412 (2014), no.12, 003. Notation analogous to
         the parameter a of Poisson, PRD 57, (1998) 5287-5290 or C_ES^2 in Porto & Rothstein, PRD 78 (2008), 044013
        
         The implementation uses the I-Love-Q fits of Table I of Yunes-Yagi
         paper, PRD 88, 023009, the bar{Q}(bar{\lambda)^{tid}) relation, line 3 of the table. The dimensionless bar{\lambda} love number is related to our apsidal constant as lambda = 2/3 k2/(C^5) so that both quantities have to appear here.*/
        
        //BNS effective spin parameter
        double a02      = C_Q1*at1*at1 + 2.*at1*at2 + C_Q2*at2*at2;
	                      
        //tidally-modified centrifugal radius
        double rc2 = r2 + a02*(1.+2.*u);
        rc         = sqrt(rc2);
        drc_dr     = r/rc*(1.-a02*u3);
        d2rc_dr2   = 1./rc*( 1.-drc_dr*r/rc*(1.-a02*u3)+2.*a02*u3);

	//NO spin-spin-tidal couplings
	/*double rc2 = r2;
	  rc = r;
	  drc_dr = 1;
	  d2rc_dr2 = 0;*/
	
    }
    else
    {

        double X12      = sqrt(1.-4.*nu);        
        double alphanu2 = 1. + 0.5/aK2*(- at2*at2*(5./4. + 5./4.*X12 + nu/2.) - at1*at1*(5./4. - 5./4.*X12 +nu/2.) + at1*at2*(-2.+nu));
        
        double rc2 = r2 + aK2*(1. + 2.*alphanu2/r);
        rc         = sqrt(rc2);
        drc_dr     = r/rc*(1.+aK2*(-alphanu2*u3 ));
        d2rc_dr2   = 1./rc*(1.-drc_dr*r/rc*(1.-alphanu2*aK2*u3)+ 2.*alphanu2*aK2*u3);
    }
    
    static double result[]= {rc, drc_dr, d2rc_dr2};
    
    return result;
}

