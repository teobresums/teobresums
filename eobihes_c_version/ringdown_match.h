//
//  ringdown_match.h
//
//  Created by Philipp Fleig on 10/09/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef ringdown_match_h
#define ringdown_match_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "cmath"

gsl_complex ringdown_match(double x, int k, vector<double> a1,vector<double> a2,vector<double> a3,vector<double> a4,vector<double> b1,vector<double> b2,vector<double> b3,vector<double> b4,vector<gsl_complex> sigma){

    gsl_complex psi;
                    
    double amp  =  ( a1[k] * tanh(a2[k]*x +a3[k]) + a4[k] ) ;
    double phase  = -b1[k]*log((1. + b3[k]*exp(-b2[k]*x) + b4[k]*exp(-2.*b2[k]*x))/(1.+b3[k]+b4[k]));
                    
    psi.dat[0] = amp * exp(-sigma[k].dat[0]*x);
    psi.dat[1] = - (phase - sigma[k].dat[1]*x); //NOTE: minus sign in front by convention
    
    return psi;
    
}

#endif /* ringdown_match_h */
