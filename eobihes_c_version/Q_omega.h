//
//  Q_omega.h
//
//  Created by Philipp Fleig on 24/10/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef Q_omega_h
#define Q_omega_h

//NOTE: routine not in use
int Q_omega(vector<double> t,vector<double> f,std::ofstream& Qomg) {
    
    double dt    = t[1]-t[0];
    double oodt  = 1./dt;
    
    double d1f;
    
    const double c = 1./12.;
    
    for (long int i=0; i<=t.size()-3; i++) {
        switch (i) {
            case 0:
                d1f = c*(-25.*f[i]+48.*f[i+1]-36.*f[i+2]+16.*f[i+3]-3.*f[i+4])*oodt;
                break;
            case 1:
                d1f = c*(-3.*f[i-1]-10.*f[i]+18.*f[i+1]-6.*f[i+2]+f[i+3])*oodt;
                break;
            default:
                d1f = c*(8.*(f[i+1]-f[i-1]) - f[i+2] + f[i-2])*oodt;
                break;
        }
        double Qomega = (f[i]*f[i])/d1f;
        Qomg << t[i] << "\t"<< f[i] << "\t" << d1f << "\t" << Qomega << endl;
    }
    
    
    return 0;
    
}

#endif /* Q_omega_h */
