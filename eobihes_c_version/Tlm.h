//
//  Tlm.h
//
//  Created by Philipp Fleig on 03/03/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef Tlm_h
#define Tlm_h

#include <stdlib.h>
#include "math.h"
#include <gsl/gsl_math.h>
#include "multipole_index.h"

/* factorials evaluated */
double f14[] = {1.,         1.,          2.,
    6.,         24.,         120.,
    720.,       5040.,       40320.,
    362880.,    3628800.,    39916800.,
    479001600., 6227020800., 87178291200.};

vector<double> Tlm(const double w){
    
    double hhatk,x2, y,prod;
    vector<double> MTlm(35);
    int kmax=35;
    
    for (int i=kmax; i--;) {
        
        hhatk = M[i]*w;
        
        x2   = 4.*hhatk*hhatk;
        prod = 1.;
        for(int j = 1; j <= L[i]; j++ ) {
            prod *= (j*j+x2);
        }
        
        y = 4.*pi*hhatk;
        y = y/( 1. - exp(-y) );
        
        // put the three pieces together
        MTlm[i]=sqrt( 1./(f14[L[i]]*f14[L[i]]) * y * prod );
    }

    return MTlm;
}

#endif /* Tlm_h */
