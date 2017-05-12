//
//  EOB_Metric.h
//
//  Created by Philipp Fleig on 20/04/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef EOB_Metric_h
#define EOB_Metric_h

#include <gsl/gsl_math.h>
#include <math.h>
#include <vector>

#include "s_get_rc.h"
#include "s_A5PNlog.h"

vector<double> s_Metric(double r, void *params){
/*
%                       This function computes the EOB metric potentials
%                       A(r), B(r), and their derivatives, as functions
%                       of the EOB (Boyer-Lindquist) radius r.
%                       USAGE:
%
%                       [A dA d2A B dB] = EOB_Metric(r,d2A_flag)
%
%                       where r is the EOB radius; dA,dB and d2A are the
%                       first and second derivatives of A and B. The
%                       default value for d2A_flag is 'no'. One must set
%                       it up to 'yes' if the computation of the second
%                       radial derivative of A is needed.
*/
    
double nu = (*(input *)params).nu;
double aK2 = (*(input *)params).aK2;

vector<double> rc_vec = s_get_rc(r,aK2,params); //[rc, drc, d2rc]
double rc = rc_vec[0];
double drc = rc_vec[1];
double d2rc = rc_vec[2];

double r2  = r*r;
double u   = 1./r;
double u2  = u*u;
double u3  = u2*u;
double u4  = u2*u2;
double uc  = 1./rc;
double uc2 = uc*uc;
double uc3 = uc2*uc;
    
//vector<double> metric = A5pnP15(rc, nu);
vector<double> metric = s_A5PNlog(rc,nu);
double Aorb = metric[0];
double dAorb = metric[1];
double d2Aorb = metric[2];
double Dorb = metric[3];

double AKerr_Multipole   = (1+2.*uc)/(1+2.*u);

double fss   = 1.;

// Correction factor
double A = Aorb*AKerr_Multipole*fss;

double dA = dAorb*drc*(1+2*uc)/(1+2*u) - 2.*Aorb*drc*uc2/(1+2*u) + 2.*Aorb*(1.+2*uc)*u2/((1+2*u)*(1+2*u));

double d2A = d2Aorb*(1+2*uc)/(1+2*u) + 4*dAorb*( u2*(1+2*uc)/((1+2*u)*(1+2*u)) - uc2/(1+2*u)*drc) + Aorb*(-4*u3*(1+2*uc)/((1+2*u)*(1+2*u)) + 8*u4*(1+2*uc)/((1+2*u)*(1+2*u)*(1+2*u))+4*uc3*(1+2*u)*drc*drc - 2*uc2/(1+2*u)*d2rc);


// The B function
double B  = r2*uc2*Dorb/A;

    return {A,B,dA,d2A};
}



#endif /* EOB_Metric_h */
