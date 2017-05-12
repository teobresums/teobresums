//
//  hlmNewt.h
//
//  Created by Philipp Fleig on 17/03/2016.
//  Copyright © 2016 Philipp Fleig. All rights reserved.
//

#ifndef hlmNewt_h
#define hlmNewt_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "cmath"

vector<gsl_complex> hlmNewt(const double r, const double Omega, const double phi, const double nu,bool tidal_flag){

/*
EOB_hlmNewt Computes the leading-order (Newtonian) prefactor  of the
            multipolar resummed waveform.

   hlmNewt = EOB_hlmNewt(r,Omega,phi, nu, EOBopt)

   Reference: Damour, Iyer & Nagar, PRD 79, 064004 (2009)

 TODO: this routine requires optimization
 - precompute coefficients c(nu)
 - evaluate efficiently polynomials
*/


// Shorthands
double nu2   = nu*nu;
double nu3   = nu*nu2;

double vphi  = r*Omega;
double vphi2 = vphi*vphi;
double vphi3 = vphi*vphi2;
double vphi4 = vphi*vphi3;
double vphi5 = vphi*vphi4;
double vphi6 = vphi*vphi5;
double vphi7 = vphi*vphi6;
double vphi8 = vphi*vphi7;


// Poly in nu
const double p1 = 1.;
double p2 = sqrt(1.-4.*nu);
const double p3 = (3.*nu-1.);
const double p4 = (2.*nu-1.)*sqrt(1.-4.*nu);
const double p5 = 1.-5.*nu+5.*nu2;
const double p6 = (1.-4.*nu+3.*nu2)*sqrt(1.-4.*nu);
const double p7 = 7.*nu3 - 14.*nu2 + 7.*nu -1.;
    
    if (tidal_flag==true) {
        p2 = 1.;
    }
    
    int kmax = 35;
    
    const double phix2 = 2. * phi;
    const double phix3 = 3. * phi;
    const double phix4 = 4. * phi;
    const double phix5 = 5. * phi;
    const double phix6 = 6. * phi;
    const double phix7 = 7. * phi;
    
    
    double M[] = {
        phi,phix2,
        phi,phix2,phix3,
        phi,phix2,phix3,phix4,
        phi,phix2,phix3,phix4,phix5,
        phi,phix2,phix3,phix4,phix5,phix6,
        phi,phix2,phix3,phix4,phix5,phix6,phix7,
        phi,phix2,phix3,phix4,phix5,phix6,phix7,8.*phi};


vector<gsl_complex> hlmNewt(kmax);
    
    const double pv23 = p2 * vphi3;
    const double pv34 = p3 * vphi4;
    const double pv45 = p4 * vphi5;
    const double pv56 = p5 * vphi6;
    const double pv67 = p6 * vphi7;
    const double pv78 = p7 * vphi8;
    

// Compute hlmNewt (without phase factor) in complex Polar coords

// l=2 ------------------------------------------------------------------

hlmNewt[1].dat[1] = pi - M[1];
hlmNewt[1].dat[0] = 8.*sqrt(pi/5.)         * p1 * vphi2; //... .*exp(-2*1i.*phi);
    
hlmNewt[0].dat[1] = 3.*pi/2. - M[0];
hlmNewt[0].dat[0] = 8./3.*sqrt(pi/5.)  * pv23;  // * p2 * vphi3; //... .*exp(-1i.*phi);

// l=3 ------------------------------------------------------------------

hlmNewt[4].dat[1] = pi/2. - M[4];
hlmNewt[4].dat[0] = 3.*sqrt(6.*pi/7.)  * pv23; //  * p2 * vphi3; //... .*exp(-3*1i.*phi);

hlmNewt[3].dat[1] = - M[3];
hlmNewt[3].dat[0] = 8./3.*sqrt(pi/7.)     *pv34;//   * p3 * vphi4; //... .*exp(-2*1i.*phi);

hlmNewt[2].dat[1] = 3.*pi/2. - M[2];
hlmNewt[2].dat[0] = 1./3.*sqrt(2.*pi/35.) *pv23;// * p2 * vphi3; //... .*exp(-1i.*phi);

// l=4 ------------------------------------------------------------------

hlmNewt[8].dat[1] = pi - M[8];
hlmNewt[8].dat[0] = 64./9.*sqrt(pi/7.)   *pv34;//    * p3 * vphi4; //... .*exp(-4*1i.*phi);

hlmNewt[7].dat[1] = 3.*pi/2. - M[7];
hlmNewt[7].dat[0] = 9./5*sqrt(2*pi/7.)  *pv45;// * p4 * vphi5; //... .*exp(-3*1i.*phi);

hlmNewt[6].dat[1] = - M[6];
hlmNewt[6].dat[0] = 8./63.*sqrt(pi) * pv34;//        * p3 * vphi4; //... .*exp(-2*1i.*phi);

hlmNewt[5].dat[1] = pi/2. - M[5];
hlmNewt[5].dat[0] = 1./105.*sqrt(2.*pi) * pv45;//   * p4 * vphi5; //... .*exp(-1i.*phi);

// l=5 ------------------------------------------------------------------

hlmNewt[13].dat[1] = pi/2. - M[13];
hlmNewt[13].dat[0] = 125./12.*sqrt(5.*pi/66.) * pv45;//* p4 * vphi5; //... .*exp(-5*1i.*phi);

hlmNewt[12].dat[1] = - M[12];
hlmNewt[12].dat[0] = 256./45.*sqrt(pi/33.) *pv56;//    * p5 * vphi6; //... .*exp(-4*1i.*phi);

hlmNewt[11].dat[1] = 3.*pi/2. - M[11];
hlmNewt[11].dat[0] = 9./20.*sqrt(3*pi/22.) *pv45;// * p4 * vphi5; //... .*exp(-3*1i.*phi);

hlmNewt[10].dat[1] = pi - M[10];
hlmNewt[10].dat[0] = 16./135.*sqrt(pi/11.) *pv56;//   * p5 * vphi6; //... .*exp(-2*1i.*phi);

hlmNewt[9].dat[1] = pi/2. - M[9];
hlmNewt[9].dat[0] = 1./180.*sqrt(pi/77.) *pv45;//    * p4 * vphi5; //... .*exp(-1i.*phi);

// l=6 ------------------------------------------------------------------

hlmNewt[19].dat[1] = pi - M[19];
hlmNewt[19].dat[0] = 432./5*sqrt(pi/715.)  *pv56;//    * p5 * vphi6; //... .*exp(-6.*1i.*phi);

hlmNewt[18].dat[1] = 3.*pi/2. - M[18];
hlmNewt[18].dat[0] = 625./63*sqrt(5.*pi/429.) * pv67;//* p6 * vphi7; //... .*exp(-5.*1i.*phi);

hlmNewt[17].dat[1] = - M[17];
hlmNewt[17].dat[0] = 1024./495.*sqrt(2.*pi/195.) *pv56;// * p5 * vphi6; //... .*exp(-4.*1i.*phi);

hlmNewt[16].dat[1] = pi/2. - M[16];
hlmNewt[16].dat[0] = 81./385.*sqrt(pi/13.)  *pv67;//  * p6 * vphi7; //... .*exp(-3.*1i.*phi);

hlmNewt[15].dat[1] = pi - M[15];
hlmNewt[15].dat[0] = 16./1485.*sqrt(pi/13.) *pv56;//    * p5 * vphi6; //... .*exp(-2.*1i.*phi);

hlmNewt[14].dat[1] = 3.*pi/2. - M[14];
hlmNewt[14].dat[0] = 1./2079.*sqrt(2.*pi/65.) *pv67;//  * p6 * vphi7; //... .*exp(-1i.*phi);

// l=7 ------------------------------------------------------------------

hlmNewt[26].dat[1] = pi/2. - M[26];
hlmNewt[26].dat[0] = 16807./180.*sqrt(7.*pi/4290.) *pv67;// * p6 * vphi7; //... .*exp(-7.*1i.*phi);

hlmNewt[25].dat[1] = - M[25];
hlmNewt[25].dat[0] = 648./35.*sqrt(3.*pi/715.)  *pv78;//       * p7 * vphi8; //... .*exp(-6.*1i.*phi);

hlmNewt[24].dat[1] = 3.*pi/2. - M[24];
hlmNewt[24].dat[0] = 3125./3276.*sqrt(5.*pi/66.) *pv67;//  * p6 * vphi7; //... .*exp(-5.*1i.*phi);

hlmNewt[23].dat[1] = pi - M[23];
hlmNewt[23].dat[0] = 1024./1365.*sqrt(2.*pi/165.)  *pv78;//   * p7 * vphi8; //... .*exp(-4.*1i.*phi);

hlmNewt[22].dat[1] = pi/2. - M[22];
hlmNewt[22].dat[0] = 243./20020.*sqrt(3.*pi/10.) *pv67;//   * p6 * vphi7; //... .*exp(-3.*1i.*phi);

hlmNewt[21].dat[1] = - M[21];
hlmNewt[21].dat[0] = 8./3003.*sqrt(pi/15.)    *pv78;//        * p7 * vphi8; //... .*exp(-2.*1i.*phi);

hlmNewt[20].dat[1] = 3.*pi/2. - M[20];
hlmNewt[20].dat[0] = 1./108108.*sqrt(pi/10.)  *pv67;//      * p6 * vphi7; //... .*exp(-1i.*phi);

// l=8 ------------------------------------------------------------------

hlmNewt[34].dat[1] = - M[34];
hlmNewt[34].dat[0] = 0.; //... .*exp(-7.*1i.*phi);

hlmNewt[33].dat[1] = - M[33];
hlmNewt[33].dat[0] = 0.; //... .*exp(-6.*1i.*phi);

hlmNewt[32].dat[1] = - M[32];
hlmNewt[32].dat[0] = 0.; //... .*exp(-5.*1i.*phi);

hlmNewt[31].dat[1] = - M[31];
hlmNewt[31].dat[0] = 0.; //... .*exp(-4.*1i.*phi);

hlmNewt[30].dat[1] = - M[30];
hlmNewt[30].dat[0] = 0.; //... .*exp(-3.*1i.*phi);

hlmNewt[29].dat[1] = - M[29];
hlmNewt[29].dat[0] = 0.; //... .*exp(-2.*1i.*phi);

hlmNewt[28].dat[1] = - M[28];
hlmNewt[28].dat[0] = 0.; //... .*exp(-1i.*phi);
    
hlmNewt[27].dat[1] = pi - M[27];
hlmNewt[27].dat[0] = 131072./315.*sqrt(2.*pi/17017.) *pv78;// * p7 * vphi8; //... .*exp(-8.*1i.*phi);

return hlmNewt;
}


#endif /* hlmNewt_h */
