//
//  input_struc.h
//
//  Created by Philipp Fleig on 16/10/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef input_struc_h
#define input_struc_h

struct input {         /* declaration of structure for input parameters */
    int   lm;
    double q;
    double nu;
    double r0;
    double dt;
    double chi1;
    double chi2;
    double S1;
    double S2;
    double S;
    double Sstar;
    double X1;
    double X2;
    double a1;
    double a2;
    double aK;
    double aK2;
    bool NQC;
    bool tidal;
    bool RWZ;
    bool speedy;
    bool dynamics;
    bool waveform;
    bool spin;
    double rLR;
    int solver_scheme;
    double cN3LO;
    
};

struct input_flags {bool NQC;bool spin;bool tidal;bool RWZ;bool speedy;bool dynamics;bool waveform;};

#endif /* input_struc_h */
