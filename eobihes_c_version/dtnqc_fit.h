//
//  dtnqc_fit.h
//  EOB_IHES_cpp_bitbucket
//
//  Created by Philipp Fleig on 10/04/2017.
//  Copyright © 2017 Philipp Fleig. All rights reserved.
//

#ifndef dtnqc_fit_h
#define dtnqc_fit_h

double dtnqc_fit(double chi,double chi0){
    // function providing a fit of Deltat_NQC vs chi
    // via a simple rational function
    double n1    = -16.06288206;
    double d1    =  -4.04266459;
    double x     = chi-chi0;
    double dtnqc = (1.+n1*x)/(1.+d1*x);
    
    return dtnqc;
    
}

#endif /* dtnqc_fit_h */
