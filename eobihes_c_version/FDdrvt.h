//
//  FDdrvt.h
//
//  Created by Philipp Fleig on 07/03/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef FDdrvt_h
#define FDdrvt_h

vector<double> FDdrvt(const vector<double> f,const vector<double> t, int varargin, const int grid_length) {

const double dt    = t[1]-t[0];
const double oodt  = 1./dt;
    
vector<double> d1f(grid_length);
    
    
    const double c = 1./12.;
    
    for (int i=grid_length; i--; ) {
        switch (i) {
            case 0://0
                d1f[i] = c*(-25.*f[i]+48.*f[i+1]-36.*f[i+2]+16.*f[i+3]-3.*f[i+4])*oodt;
                //d2f[i] = c*(45*f[i]-154*f[i+1]+214*f[i+2]-156*f[i+3]+61*f[i+4]-10*f[i+5])*oodt2;
                break;
            case 1://1
                d1f[i] = c*(-3.*f[i-1]-10.*f[i]+18.*f[i+1]-6.*f[i+2]+f[i+3])*oodt;
                //d2f[i] = c*(10*f[i-1]-15*f[i]-4*f[i+1]+14*f[i+2]-6*f[i+3]+f[i+4])*oodt2;
                break;
            case 10://10
                d1f[i] = - c*(-3.*f[i+1]-10.*f[i]+18.*f[i-1]-6.*f[i-2]+f[i-3])*oodt;
                //d2f[i] = c*(10*f[i+1]-15*f[i]-4*f[i-1]+14*f[i-2]-6*f[i-3]+f[i-4])*oodt2;
                break;
            case 11://11
                d1f[i] = - c*(-25.*f[i]+48.*f[i-1]-36.*f[i-2]+16.*f[i-3]-3.*f[i-4])*oodt;
                //d2f[i] = c*(45*f[i]-154*f[i-1]+214*f[i-2]-156*f[i-3]+61*f[i-4]-10*f[i-5])*oodt2;
                break;
            default: d1f[i] = c*(8.*(f[i+1]-f[i-1]) - f[i+2] + f[i-2])*oodt;
                //d2f[i] = c*(-30*f[i]+16*(f[i+1]+f[i-1])-(f[i+2]+f[i-2]))*oodt2;
                break;
        }
    }

    
return d1f;
    
}


#endif /* FDdrvt_h */
