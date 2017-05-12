//
//  EOB_get_rc.h
//
//  Created by Philipp Fleig on 20/04/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef EOB_get_rc_h
#define EOB_get_rc_h

vector <double> s_get_rc(double r, double aK2, void *params){

    
    double nu = (*(input *)params).nu;
    double chi1 = (*(input *)params).chi1;
    double chi2 = (*(input *)params).chi2;
    double X1 = (*(input *)params).X1;
    double X2 = (*(input *)params).X2;
    bool tidal_flag = (*(input *)params).tidal;
    
    double rc,drc_dr, d2rc_dr2;
    
    if (tidal_flag==true) {
        rc = r;
        drc_dr = 1.;
        d2rc_dr2 = 0.;
    } else {
        // shorthands
        double nu2 = nu*nu;
        double u   = 1./r;
        double u2  = u*u;
        double u3  = u*u2;
        double r2  = r*r;

        //ratios
        double X1_X2 = X1/X2;
        double X2_X1 = X2/X1;

        //'equatorial'
        double a11 = nu/16.*(-32*nu - 22*nu2 + X1_X2*(21-44*nu-44*nu2) + (X1_X2*X1_X2)*(-21*nu-22*nu2));
        double a12 = nu/8.*(24 - 53*nu-44*nu2 + (X1_X2 + X2_X1)*(-32*nu-22*nu2));
        double a22 = nu/16.*(-32*nu - 22*nu2 + X2_X1*(21-44*nu-44*nu2) + (X2_X1*X2_X1)*(-21*nu-22*nu2));
        double c11 = nu/16.*(88*nu + 14*nu2 + X1_X2*(-117+196*nu+28*nu2) + (X1_X2*X1_X2)*(117*nu+14*nu2));
        double c12 = nu/8.*(-120+229*nu+28*nu2 + (X1_X2 + X2_X1)*(112*nu + 14*nu2));
        double c22 = nu/16.*(88*nu + 14*nu2 + X2_X1*(-117+196*nu+28*nu2) + (X2_X1*X2_X1)*(117*nu+14*nu2));

        // alphanu2 from Balmelli-Jetzer, general case
        double alphanu2 = 1. + 0.5/aK2*( (a11+c11)*chi1*chi1 + (a22+c22)*chi2*chi2 + (a12+c12)*chi1*chi2);

        double rc2      = r2 + aK2*(1. + 2*alphanu2/r);
        rc       = sqrt(rc2);
        drc_dr   = r/rc*(1+aK2*(-alphanu2*u3 ));
        d2rc_dr2 = 1./rc*(1-drc_dr*r/rc*(1.-alphanu2*aK2*u3)+ 2*alphanu2*aK2*u3);
    }

    return {rc, drc_dr, d2rc_dr2};
}

#endif /* EOB_get_rc_h */
