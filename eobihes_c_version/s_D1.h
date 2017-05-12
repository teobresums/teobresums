//
//  EOB_D1.h
//
//  Created by Philipp Fleig on 27/04/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef EOB_D1_h
#define EOB_D1_h

vector<double> s_D1(vector<double> f,vector<double> x,int Nmax){
/* Computes the first derivative of the function. Centered but at the edges. USAGE: df = EOB_D1(f,x) */
int Nmin = 0;
//int Nmax = 11;
    
    vector<double> df(Nmax+1);
    for(int i=2;i<=Nmax-2;i++){
        df[i] = 1./3.*(8.*f[1+i] - f[2+i] - 8.*f[i-1] + f[i-2])/(x[2+i]-x[i-2]);
    }

    // 4th order boundaries
    df[0]   = (-24./17.*f[Nmin] + 59./34.*f[Nmin+1] - 4./17.*f[Nmin+2] - 3./34.*f[Nmin+3])/(x[Nmin+1]-x[Nmin]);
    df[1] = (-1./2.*f[Nmin] + 1./2.*f[Nmin+2] )/(x[Nmin+1]-x[Nmin]);

    df[Nmax]   = -(-24./17.*f[Nmax] + 59./34.*f[Nmax-1] - 4./17.*f[Nmax-2] - 3./34.*f[Nmax-3])/(x[Nmax]-x[Nmax-1]);
    df[Nmax-1] = -(-1./2.*f[Nmax] + 1./2.*f[Nmax-2])/(x[Nmax]-x[Nmax-1]);

    return df;

}
#endif /* EOB_D1_h */
