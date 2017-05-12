//
//  A_NumDenom.h
//  grav_waves
//
//  Created by Philipp Fleig on 16/02/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef A_NumDenom_h
#define A_NumDenom_h

#include <math.h>
#include <vector>

using namespace::std;

vector<double> A_NumDenom(const double r, const vector<double> a, const double nu){
    
    vector<double> frac(4);
    
    // Shorthands
    const double u = 1./r;
    const double u2=u*u;
    const double u3=u2*u;
    const double u4=u3*u;
    const double u5=u4*u;
    
    const double sm = nu;
    
    // Point-mass PN coefs
    const double a3=a[0];
    const double a4=a[1];
    const double a5=a[4];
    const double a6=a[7];
    const double a5l=a[3];
    const double a6l=a[6];
    
    const double C1 = (8.*a3 + 3.*a4 + a5);
    const double C2 = (4.*a4 + 2.*a5 + a6 + a3*(8 - a3*sm));
    const double C3 = (12.*a3 + 4.*a4 + a5);
    const double C4 = (-32 + C3*sm);
    const double C6 = (4.*a3 + a4);
    const double C5 = (16.*a6 + (2.*a4 + a5)*(8.*a3 + 2.*a4 + a5)*sm - C6*a6*sm + a3*a3*C6*sm*sm);
    const double C7 = (64. - (4.*C1 + a6)*sm + a3*a3*sm*sm);
    const double C8 = a6*(8 - a3*sm) + sm*(8.*a3*a3 + 16.*a3*a4 + 4.*a4*a4 + a3*a3*a3*sm) + a5*(16 + a4*sm);
    const double C9 = 4.*(4.*a4 + 2.*a5 + a6) + a3*(8.*a3 + 4.*a4 + a5)*sm;
    const double C10 = (-8 + a3*sm);
    const double C11 = (-32. + C3*sm);
    const double C12 = (-4.*a4 - 2.*a5 - a6 + a3*C10);
    
    frac[0]=1. + (C7*u)/C11;
    frac[1]=(-32. + sm*(12.*a3 + 4.*a4 + a5 - C2*u - 2.*C2*u2 - (C9)*u3 - (C8)*u4 - C5*u5))/C11;
    
    frac[2]=(-((4.*a5l + a6l)*sm*C4) + C4*C7 + a5l*sm*(-64 + sm*(12.*a4 + 4.*a5 + a6 + a3*(32 - a3*sm))))/(C4*C4);
    frac[3]=(sm*(-((2.*a5l + a6l)*C4) + a5l*sm*C2 + C4*C12 - 2.*(2.*a5l + a6l)*C4*u + 2.*a5l*sm*C2*u + 4.*C4*C12*u + a5l*sm*(C9)*u2 - 3.*(C9)*C4*u2 - C4*(4.*a6l + a5l*(8 + a3*sm))*u2 + a5l*sm*(C8)*u3 - 4.*C4*(C8)*u3 + C4*(a6l*C10 - a5l*(16 + a4*sm))*u3 + a5l*sm*C5*u4 - 5.*C4*C5*u4 + C4*(-2.*(4.*a3 + 2.*a4 + a5)*a5l*sm + a6l*(-16 + 4.*a3*sm + a4*sm))*u4))/(C11*C11);
    
    
    
    return frac;
}

#endif /* A_NumDenom_h */


/*
frac[0]=1. + ((64. - (4.*(8.*a3 + 3.*a4 + a5) + a6)*sm + a3*a3*sm*sm)*u)/(-32. + (12.*a3 + 4.*a4 + a5)*sm);
frac[1]=(-32. + sm*(12.*a3 + 4.*a4 + a5 - (4.*a4 + 2.*a5 + a6 + a3*(8 - a3*sm))*u - 2.*(4.*a4 + 2.*a5 + a6 + a3*(8 - a3*sm))*u2 - (4.*(4.*a4 + 2.*a5 + a6) + a3*(8.*a3 + 4.*a4 + a5)*sm)*u3 - (a6*(8 - a3*sm) + sm*(8.*a3*a3 + 16.*a3*a4 + 4.*a4*a4 + a3*a3*a3*sm) + a5*(16 + a4*sm))*u4 - (16.*a6 + (2.*a4 + a5)*(8.*a3 + 2.*a4 + a5)*sm - (4.*a3 + a4)*a6*sm + a3*a3*(4.*a3 + a4)*sm*sm)*u5))/(-32. + (12.*a3 + 4.*a4 + a5)*sm);

frac[2]=(-((4.*a5l + a6l)*sm*(-32 + (12.*a3 + 4.*a4 + a5)*sm)) + (-32 + (12.*a3 + 4.*a4 + a5)*sm)*(64 - (4.*(8.*a3 + 3.*a4 + a5) + a6)*sm + a3*a3*sm*sm) + a5l*sm*(-64 + sm*(12.*a4 + 4.*a5 + a6 + a3*(32 - a3*sm))))/((-32 + (12.*a3 + 4.*a4 + a5)*sm)*(-32 + (12.*a3 + 4.*a4 + a5)*sm));
frac[3]=(sm*(-((2.*a5l + a6l)*(-32 + (12.*a3 + 4.*a4 + a5)*sm)) + a5l*sm*(4.*a4 + 2.*a5 + a6 + a3*(8 - a3*sm)) + (-32 + (12.*a3 + 4.*a4 + a5)*sm)*(-4.*a4 - 2.*a5 - a6 + a3*(-8 + a3*sm)) - 2.*(2.*a5l + a6l)*(-32 + (12.*a3 + 4.*a4 + a5)*sm)*u + 2.*a5l*sm*(4.*a4 + 2.*a5 + a6 + a3*(8 - a3*sm))*u + 4.*(-32 + (12.*a3 + 4.*a4 + a5)*sm)*(-4.*a4 - 2.*a5 - a6 + a3*(-8 + a3*sm))*u + a5l*sm*(4.*(4.*a4 + 2.*a5 + a6) + a3*(8.*a3 + 4.*a4 + a5)*sm)*u2 - 3.*(4.*(4.*a4 + 2.*a5 + a6) + a3*(8.*a3 + 4.*a4 + a5)*sm)*(-32 + (12.*a3 + 4.*a4 + a5)*sm)*u2 - (-32 + (12.*a3 + 4.*a4 + a5)*sm)*(4.*a6l + a5l*(8 + a3*sm))*u2 + a5l*sm*(a6*(8 - a3*sm) + sm*(8.*a3*a3 + 16.*a3*a4 + 4.*a4*a4 + a3*a3*a3*sm) + a5*(16 + a4*sm))*u3 - 4.*(-32 + (12.*a3 + 4.*a4 + a5)*sm)*(a6*(8 - a3*sm) + sm*(8.*a3*a3 + 16.*a3*a4 + 4.*a4*a4 + a3*a3*a3*sm) + a5*(16 + a4*sm))*u3 + (-32 + (12.*a3 + 4.*a4 + a5)*sm)*(a6l*(-8 + a3*sm) - a5l*(16 + a4*sm))*u3 + a5l*sm*(16.*a6 + (2.*a4 + a5)*(8.*a3 + 2.*a4 + a5)*sm - (4.*a3 + a4)*a6*sm + a3*a3*(4.*a3 + a4)*sm*sm)*u4 - 5.*(-32 + (12.*a3 + 4.*a4 + a5)*sm)*(16.*a6 + (2.*a4 + a5)*(8.*a3 + 2.*a4 + a5)*sm - (4.*a3 + a4)*a6*sm + a3*a3*(4.*a3 + a4)*sm*sm)*u4 + (-32 + (12.*a3 + 4.*a4 + a5)*sm)*(-2.*(4.*a3 + 2.*a4 + a5)*a5l*sm + a6l*(-16 + 4.*a3*sm + a4*sm))*u4))/((-32. + (12.*a3 + 4.*a4 + a5)*sm)*(-32. + (12.*a3 + 4.*a4 + a5)*sm));
*/
