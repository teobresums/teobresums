//
//  Eulerlog.h
//
//  Created by Philipp Fleig on 26/02/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef Eulerlog_h
#define Eulerlog_h

#include <math.h>
#include <gsl/gsl_math.h>

double Eulerlog(const double x,const double m){
    
    const double EulerGamma=0.5772156649015328606065121;
    const double Log2=0.6931471805599453094172321;

    return EulerGamma + Log2 + log(m) + 0.5*log(x);
}


#endif /* Eulerlog_h */
