//
//  NQC.h
//
//  Created by Philipp Fleig on 23/02/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef NQC_h
#define NQC_h

vector<double> NQC(const double r,const double prstar, const double Omega, const double ddotr, int i){


vector<double> n(6);
    
    switch (i) {
            // l=2 -------------------------------------------------------------------
            // (2,1)
        case 0:
            // NQC corrections to the modulus
            n[0] = (prstar/(r*Omega))*(prstar/(r*Omega));
            n[1] = ddotr/(r*Omega*Omega);
            n[2] = n[0]*prstar*prstar;
            
            
            //% NQC corrections to the phase
            n[3] = prstar/(r*Omega);
            n[4] = n[3]*pow(Omega, 2./3.);
            n[5] = n[4]*prstar*prstar;
            break;
            // (2,2)
        case 1:
            // NQC corrections to the modulus
            n[0] = (prstar/(r*Omega))*(prstar/(r*Omega));
            n[1] = ddotr/(r*Omega*Omega);
            n[2] = n[0]*prstar*prstar;
            
            
            //% NQC corrections to the phase
            n[3] = prstar/(r*Omega);
            n[4] = n[3]*(r*Omega)*(r*Omega);
            n[5] = n[4]*prstar*prstar;
            break;
            // l=3 -------------------------------------------------------------------
            // (3,3)
        case 4:
            // NQC corrections to the modulus
            n[0] = (prstar/(r*Omega))*(prstar/(r*Omega));
            n[1] = ddotr/(r*Omega*Omega);
            n[2] = n[0]*prstar*prstar;
            
            
            //% NQC corrections to the phase
            n[3] = prstar/(r*Omega);
            n[4] = n[3]*pow(Omega, 2./3.);
            n[5] = n[4]*prstar*prstar;
            break;
        default:
            n[0] = 0.;
            n[1] = 0.;
            n[2] = 0.;
            
            n[3] = 0.;
            n[4] = 0.;
            n[5] = 0.;
            break;
    }
    
return n;
}


#endif /* NQC_h */
