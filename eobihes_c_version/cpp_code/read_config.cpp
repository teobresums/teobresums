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

#include <iostream>
#include <list>
#include <fstream>
#include <stdio.h>
#include <string>

#include "s_GS.h"
#include "input_struc.h"

#include <math.h>
#include "read_config.h"

using namespace::std;

double logQ(double x){
    
    // implements the logQ-vs-log(lambda) fit of Table I of Yunes-Yagi
    // here x = log(lambda) and the output is the log of the coefficient
    // that describes the quadrupole deformation due to spin.
    double ai = 0.194;
    double bi = 0.0936;
    double ci = 0.0474;
    double di = -4.21e-3;
    double ei = 1.23e-4;
    double x2 = x*x;
    double x3 = x*x2;
    double x4 = x*x3;
    
    return ai + bi*x + ci*x2 + di*x3 + ei*x4;
}

//input read_config(double q, double chi1, double chi2, double r0){
input read_config(char *fname){

    input params;            /* declaration of variable of struture type */
    
    string param_name;
    double param_value;
    //ifstream fin ("file_run.par");
    ifstream fin (fname);
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
	  params.q = param_value;
	  break;
        case 1:
	  params.chi1 = param_value;
	  break;
        case 2:
	  params.chi2 = param_value;
	  break;
	case 3:
	  params.r0 = param_value;
	  break;
	case 4:
                if (param_value==0) {
                    params.NQC = false;
                } else {
                    params.NQC = true;
                }
                break;
	case 5:
                if (param_value==0) {
                    params.tidal = false;
                } else {
                    params.tidal = true;
                }
                break;
            case 6:
                if (param_value==0) {
                    params.RWZ = false;
                } else {
                    params.RWZ = true;
                }
                break;
            case 7:
                if (param_value==0) {
                    params.speedy = false;
                } else {
                    params.speedy = true;
                }
                break;
            case 8:
                if (param_value==0) {
                    params.dynamics = false;
                } else {
                    params.dynamics = true;
                }
                break;
            case 9:
                if (param_value==0) {
                    params.waveform = false;
                } else {
                    params.waveform = true;
                }
                break;
            case 10:
                params.lm = param_value;
                break;
            case 11:
                params.dt = param_value;
                break;
            case 12:
                params.solver_scheme = param_value;
                break;
            case 13:
                params.LambdaAl2 = param_value;
                break;
            case 14:
                params.LambdaAl3 = param_value;
                break;
            case 15:
	      params.LambdaAl4 = param_value;
                break;
            case 16:
                params.LambdaBl2 = param_value;
                break;
            case 17:
                params.LambdaBl3 = param_value;
                break;
            case 18:
                params.LambdaBl4 = param_value;
                break;
            default:
                break;
        }
        
        i++;
    }
    fin.close();

    // tmp
    double chi1 = params.chi1;
    double chi2 = params.chi2;
    double q = params.q;

    /*set spin flag*/
    if (chi1 != .0 || chi2 != .0) {
        params.spin = true;
        params.NQC  = false;
    } else {params.spin = false;}
    
    double nu = q/((q+1.)*(q+1.));
    params.nu = nu;
    //params.q  = q;
    //params.r0 = r0;
    
    double X1 = 0.5*(1.+sqrt(1.-4.*nu));
    double X2 = 1. - X1;
    
    params.X1   = X1;
    params.X2   = X2;
    
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

    double cN3LO = c3_fit_global(nu,chi1,chi2,X1,X2,a1,a2,params.tidal);
    if (params.tidal==true) {
        params.NQC = false;
    }        
    double lambda1 = params.LambdaAl2;
    double lambda2 = params.LambdaBl2;
    double logC_Q1 = logQ(log(lambda1));
    double logC_Q2 = logQ(log(lambda2));
    double C_Q1    = exp(logC_Q1);
    double C_Q2    = exp(logC_Q2);
    
    params.C_Q1 = C_Q1;
    params.C_Q2 = C_Q2;
    
    return params;
}

