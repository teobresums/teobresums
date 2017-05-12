//
//  EOB_HorizonFlux.h
//
//  Created by Philipp Fleig on 25/04/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef EOB_HorizonFlux_h
#define EOB_HorizonFlux_h

#include "math.h"

double s_HorizonFlux(double x,double Heff,double jhat,double nu,double X1,double X2,double chi1,double chi2){

double x2   = x*x;
double x3   = x2*x;
double x4   = x3*x;
double x5   = x4*x;
double v5   = sqrt(x5);

vector<double> cv5(2);
vector<double> cv8(2);

// coefficients of the v^5 term (Alvi leading order)
cv5[0] = -1./4.*chi1*(1+3*chi1*chi1)*X1*X1*X1;
cv5[1] = -1./4.*chi2*(1+3*chi2*chi2)*X2*X2*X2;

// coefficients of the v^8=x^4 term
cv8[0] = 0.5*(1+sqrt(1-chi1*chi1))*(1+3*chi1*chi1)*X1*X1*X1*X1;
cv8[1] = 0.5*(1+sqrt(1-chi2*chi2))*(1+3*chi2*chi2)*X2*X2*X2*X2;

double FH22_S = (cv5[0]+cv5[1])*v5;
double FH22  = (cv8[0]+cv8[1])*x4;
double FH21  =  0.0;

/* Newton-normalized horizon flux: use only l=2 fluxes */
double hatFH   = FH22_S + FH22 + FH21;
    
    return hatFH;
}

#endif /* EOB_HorizonFlux_h */
