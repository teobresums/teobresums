//
//  read_config.h
//  EOB_ihes
//
//  Created by Philipp Fleig on 23/01/2017.
//  Copyright © 2017 Philipp Fleig. All rights reserved.
//

#ifndef read_config_h
#define read_config_h

#include <iostream>
#include <list>
#include <fstream>
#include <stdio.h>
#include <string>

#include "s_GS.h"
#include "input_struc.h"

using namespace::std;

input read_config (double q,double chi1,double chi2,double r0) {
    
    input params;            /* declaration of variable of struture type */
    
    string param_name;
    double param_value;
    ifstream fin ("file_run.par");
    if (!fin) //checks to see if file opens properly
    {
        cerr << "error: Could not find the parameters file.";
    }

    int i = 0;
    
    while ( fin >> param_name >> param_value )
    {
        //param_values.push_back(param_value);
        cout << param_name <<"\t"<< param_value << endl;
        switch (i) {
            case 0:
                if (param_value==0) {
                    params.NQC = false;
                } else {
                    params.NQC = true;
                }
                break;
            case 1:
                if (param_value==0) {
                    params.tidal = false;
                } else {
                    params.tidal = true;
                }
                break;
            case 2:
                if (param_value==0) {
                    params.RWZ = false;
                } else {
                    params.RWZ = true;
                }
                break;
            case 3:
                if (param_value==0) {
                    params.speedy = false;
                } else {
                    params.speedy = true;
                }
                break;
            case 4:
                if (param_value==0) {
                    params.dynamics = false;
                } else {
                    params.dynamics = true;
                }
                break;
            case 5:
                if (param_value==0) {
                    params.waveform = false;
                } else {
                    params.waveform = true;
                }
                break;
            case 6:
                params.lm = param_value;
                break;
            case 7:
                params.dt = param_value;
                break;
            case 8:
                params.solver_scheme = param_value;
                break;
            default:
                break;
        }
        
        i++;
    }
    fin.close();

    /*set spin flag*/
    if (chi1 != .0 || chi2 != .0) {
        params.spin = true;
        params.NQC=false;
    } else {params.spin =false;}
    
    double nu = q/((q+1.)*(q+1.));
    
    params.nu = nu;
    params.q = q;
    params.r0 = r0;
    
    double X1 = 0.5*(1+sqrt(1.-4.*nu));
    double X2 = 1. - X1;
    
    params.chi1 = chi1;
    params.chi2 = chi2;
    params.X1 = X1;
    params.X2 = X2;
    
    double S1 = pow(params.X1, 2.) * params.chi1;
    double S2 = pow(params.X2, 2.) * params.chi2;
    
    params.S1 = S1;
    params.S2 = S2;
    
    double a1  = X1*chi1;
    double a2  = X2*chi2;
    double aK  = a1 + a2;
    double aK2 = aK*aK;
    
    params.S     = S1 + S2;        // => in the EMRL this becomes the spin of the BH
    params.Sstar = X2*a1 + X1*a2;  // => in the EMRL this becomes the spin of the particle
    
    params.a1  = X1*chi1;
    params.a2  = X2*chi2;
    params.aK  = a1 + a2;
    params.aK2 = aK2;
    
    params.rLR = 0.;
    
    if (params.tidal==true) {
        params.NQC = false;
    }
    
    double cN3LO = c3_fit_global(nu,chi1,chi2,X1,X2,a1,a2,params.tidal);
    params.cN3LO = cN3LO;
    
    return params;
}



#endif /* read_config_h */
