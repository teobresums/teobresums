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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#include <ios>
#include <cmath>
#include <math.h>
#include <vector>
#include <limits>
#include <stdio.h>
#include <fstream>
#include <iostream>
//#include <string.h> 
#include <cstring> 

#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_complex_math.h>

#include "TEOBResumS.h"


using namespace::std;

double Eulerlog(const double x,const double m)
{
    
    const double EulerGamma = 0.5772156649015328606065121;
    const double Log2       = 0.6931471805599453094172321;
    
    return EulerGamma + Log2 + log(m) + 0.5*log(x);
}

double interpolate(double dt,vector<gsl_complex> grid)
{
    double xi, yi;
    double x[]      = {0.,0.,0.,0.,0.,0.,0.};
    double y[]      = {0.,0.,0.,0.,0.,0.,0.};
    double step     = 0.01; //dt/5.;
    double omeg_max = 0.;
    double t_max    = 0.;
    bool peak_flag  = false;
    
    for (int i=0; i<=6; i++)
    {
        x[i] = grid[i].dat[0]; //time
        y[i] = grid[i].dat[1]; //omega
    }
    
    
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline    = gsl_spline_alloc (gsl_interp_cspline, 7);
    gsl_spline_init (spline, x, y, 7);
    
    for (xi = x[0]; xi < x[6]; xi += step)
    {
        yi = gsl_spline_eval (spline, xi, acc);
        if (peak_flag==false)
        {
            if (yi<omeg_max)
            {
                peak_flag=true;
            }
            else
            {
                omeg_max=yi;
                t_max=xi;
            }
        }
    }
    
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
    
    return t_max;
}

vector<double> interp_grid(vector<double> t_vec, vector<double> data, double dt)
{
    int i = 0;
    int t_length = t_vec.size();
    int grid_length = (int)((t_vec.back()-t_vec[0])/dt + 1);
    
    double xi, yi;
    vector<double> data_g(grid_length);
    vector<double> omg_interp(grid_length);
    
    /** Convert all vectors to an array */
    double* t = &t_vec[0];
    double* data_arr = &data[0];
    double step = dt;
    
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, t_length);
    gsl_spline_init (spline, t, data_arr, t_length);
    
    for (xi = t_vec[0]; xi < t_vec.back(); xi += step)
    {
        yi = gsl_spline_eval (spline, xi, acc);
        data_g[i] = yi;
        i++;
    }
    
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
    
    return data_g;
}

vector<double> FDdrvt(const vector<double> f,const vector<double> t, int varargin, const int grid_length)
{
    
    const double dt    = t[1]-t[0];
    const double oodt  = 1./dt;
    const double c = 1./12.;
    vector<double> d1f(grid_length);
    
    for (int i=grid_length; i--;)
    {
        switch (i)
        {
            case 0:
                d1f[i] = c*(-25.*f[i] + 48.*f[i+1] - 36.*f[i+2] + 16.*f[i+3] - 3.*f[i+4])*oodt;
                //d2f[i] = c*(45*f[i]-154*f[i+1]+214*f[i+2]-156*f[i+3]+61*f[i+4]-10*f[i+5])*oodt2;
                break;
            case 1:
                d1f[i] = c*(-3.*f[i-1] - 10.*f[i] + 18.*f[i+1] - 6.*f[i+2] + f[i+3])*oodt;
                //d2f[i] = c*(10*f[i-1]-15*f[i]-4*f[i+1]+14*f[i+2]-6*f[i+3]+f[i+4])*oodt2;
                break;
            case 10:
                d1f[i] = - c*(-3.*f[i+1] - 10.*f[i] + 18.*f[i-1] - 6.*f[i-2] + f[i-3])*oodt;
                //d2f[i] = c*(10*f[i+1]-15*f[i]-4*f[i-1]+14*f[i-2]-6*f[i-3]+f[i-4])*oodt2;
                break;
            case 11:
                d1f[i] = - c*(-25.*f[i] + 48.*f[i-1] - 36.*f[i-2] + 16.*f[i-3] - 3.*f[i-4])*oodt;
                //d2f[i] = c*(45*f[i]-154*f[i-1]+214*f[i-2]-156*f[i-3]+61*f[i-4]-10*f[i-5])*oodt2;
                break;
            default: d1f[i] = c*(8.*(f[i+1]-f[i-1]) - f[i+2] + f[i-2])*oodt;
                //d2f[i] = c*(-30*f[i]+16*(f[i+1]+f[i-1])-(f[i+2]+f[i-2]))*oodt2;
                break;
        }
    }
    
    return d1f;
}

vector<double> FDdrvt_omega(vector<double> f, double dt)
{
    vector<double> d1f(f.size()-2);
    const double oodt  = 1./dt;
    const double c     = 1./12.;
    
    for (unsigned long int i=0; i<f.size()-2; i++)
    {
        switch (i)
        {
            case 0:
                d1f[i] = c*(-25.*f[i]+48.*f[i+1]-36.*f[i+2]+16.*f[i+3]-3.*f[i+4])*oodt;
                //d2f[i] = c*(45*f[i]-154*f[i+1]+214*f[i+2]-156*f[i+3]+61*f[i+4]-10*f[i+5])*oodt2;
                break;
            case 1:
                d1f[i] = c*(-3.*f[i-1]-10.*f[i]+18.*f[i+1]-6.*f[i+2]+f[i+3])*oodt;
                //d2f[i] = c*(10*f[i-1]-15*f[i]-4*f[i+1]+14*f[i+2]-6*f[i+3]+f[i+4])*oodt2;
                break;
                /*
                 case 198:
                 d1f[i] = - c*(-3.*f[i+1]-10.*f[i]+18.*f[i-1]-6.*f[i-2]+f[i-3])*oodt;
                 //d2f[i] = c*(10*f[i+1]-15*f[i]-4*f[i-1]+14*f[i-2]-6*f[i-3]+f[i-4])*oodt2;
                 break;
                 case 199:
                 d1f[i] = - c*(-25.*f[i]+48.*f[i-1]-36.*f[i-2]+16.*f[i-3]-3.*f[i-4])*oodt;
                 //d2f[i] = c*(45*f[i]-154*f[i-1]+214*f[i-2]-156*f[i-3]+61*f[i-4]-10*f[i-5])*oodt2;
                 break;
                 */
            default:
                d1f[i] = c*(8.*(f[i+1]-f[i-1]) - f[i+2] + f[i-2])*oodt;
                //d2f[i] = c*(-30*f[i]+16*(f[i+1]+f[i-1])-(f[i+2]+f[i-2]))*oodt2;
                break;
        }
    }
    
    return d1f;
    
}

vector<gsl_complex> speedyTail(const double Omega, const double Hreal, const double bphys, const int L[], const int M[])
{
    int kmax = 35;
    
    double ratio_rad;
    double ratio_ang;
    double tlm_rad;
    double tlm_ang;
    
    double x;
    double x2;
    double x3;
    double x4;
    double x5;
    
    vector<double> num_ang(kmax);
    vector<gsl_complex> tlm(kmax);
    
    /** Fit coefficients*/
    /*
     double a1[] =
     {
     0.3060234167900050, 0.3060234167900050, 0.0424759238428813, 0.0424489015884926, 0.0424717446800903, 0.0215953972500844, 0.0215873812155663, 0.0215776183122621, 0.0216017621863542, 0.0128123696874894, 0.0128097056242375, 0.0128038943888768, 0.0128025242617949, 0.0128202485907368, 0.0083762045692408, 0.0083751913886140, 0.0083724067460769, 0.0083694435961860, 0.0083710364141552, 0.0083834483913443, 0.0058540393221396, 0.0058536069384738, 0.0058522594457692, 0.0058502436535615, 0.0058491157293566, 0.0058514875071582, 0.0058602498033381, 0.0042956812356573, 0.0042954784390887, 0.0042947951664056, 0.0042935886137697, 0.0042923691461384, 0.0042922256848799, 0.0042945927126022, 0.0043009106861259};
     
     double a2[] =
     {
     -0.0248640052699995, -0.0248640052699995, 0.0000597134489198, 0.0002551406918111, 0.0001741036904709, 0.0000124649041611, 0.0000685496215625, 0.0001131160409390, 0.0000419907542591, 0.0000035218982282, 0.0000219211271097, 0.0000473186962874, 0.0000524142634057, 0.0000106823372552, 0.0000012237574387, 0.0000081742188269, 0.0000201940563214, 0.0000295722761753, 0.0000260539631956, 0.0000018994753518, 0.0000004932942990, 0.0000034477210351, 0.0000092294406360, 0.0000155143073237, 0.0000183386499818, 0.0000137922469695, -0.0000007075155453, 0.0000002223410995, 0.0000016045317657, 0.0000045260028113, 0.0000082655700107, 0.0000112393599417, 0.0000115758243113, 0.0000076838709956, -0.0000014020591745
     };
     
     double a3[] = {
     0.1947531537583291, 0.1947531537583291, -0.0042911639711832, -0.0047431560217121, -0.0046577314472149, -0.0013089557502947, -0.0014343968205390, -0.0014978542575474, -0.0014329302934532, -0.0005167994164556, -0.0005573939123058, -0.0005921030407223, -0.0005978284714483, -0.0005673965369076, -0.0002409269302708, -0.0002561516055118, -0.0002723768586352, -0.0002815958312453, -0.0002792078156272, -0.0002646630240693, -0.0001261183503407, -0.0001325622938779, -0.0001403198638518, -0.0001464084186977, -0.0001485971591029, -0.0001459023931717, -0.0001384829633836, -0.0000719062974278, -0.0000749128468013, -0.0000788187384314, -0.0000824202283094, -0.0000846673495936, -0.0000849054394951, -0.0000829269749240, -0.0000788883333858
     };
     
     double a4[] = {
     -0.0254800838696074, -0.0254800838696074, 0.0006914295465364, 0.0010322294603561, 0.0010057563135650, 0.0001394203795507, 0.0002309706405978, 0.0002596611624417, 0.0002409588083156, 0.0000386949167221, 0.0000679154947896, 0.0000830199015202, 0.0000850120755064, 0.0000780125513602, 0.0000133034384660, 0.0000241813441339, 0.0000311573885555, 0.0000340233089866, 0.0000335167900637, 0.0000307571022927, 0.0000053305073331, 0.0000099143129290, 0.0000132296989826, 0.0000150959309402, 0.0000156304390748, 0.0000151274875147, 0.0000139320508803, 0.0000023959090314, 0.0000045285807761, 0.0000061918979830, 0.0000072894226381, 0.0000078251853305, 0.0000078772667984, 0.0000075606242809, 0.0000069956215270
     };
     
     double a5[] = {
     -0.0000140714300659, -0.0000140714300659, 0.0006914295465364, 0.0010322294603561, 0.0010057563135650, 0.0001394203795507, 0.0002309706405978, 0.0002596611624417, 0.0002409588083156, 0.0000386949167221, 0.0000679154947896, 0.0000830199015202, 0.0000850120755064, 0.0000780125513602, 0.0000133034384660, 0.0000241813441339, 0.0000311573885555, 0.0000340233089866, 0.0000335167900637, 0.0000307571022927, 0.0000053305073331, 0.0000099143129290, 0.0000132296989826, 0.0000150959309402, 0.0000156304390748, 0.0000151274875147, 0.0000139320508803, 0.0000023959090314, 0.0000045285807761, 0.0000061918979830, 0.0000072894226381, 0.0000078251853305, 0.0000078772667984, 0.0000075606242809, 0.0000069956215270
     };
     */
    
    double b1[] =
    {
        0.1113090643348557, 0.1112593821157397, 0.0424759238428813, 0.0424489015884926, 0.0424717446800903, 0.0215953972500844, 0.0215873812155663, 0.0215776183122621, 0.0216017621863542, 0.0128123696874894, 0.0128097056242375, 0.0128038943888768, 0.0128025242617949, 0.0128202485907368, 0.0083762045692408, 0.0083751913886140, 0.0083724067460769, 0.0083694435961860, 0.0083710364141552, 0.0083834483913443, 0.0058540393221396, 0.0058536069384738, 0.0058522594457692, 0.0058502436535615, 0.0058491157293566, 0.0058514875071582, 0.0058602498033381, 0.0042956812356573, 0.0042954784390887, 0.0042947951664056, 0.0042935886137697, 0.0042923691461384, 0.0042922256848799, 0.0042945927126022, 0.0043009106861259};
    
    double b2[] =
    {
        0.0004643273300862, 0.0009375605440004, 0.0000597134489198, 0.0002551406918111, 0.0001741036904709, 0.0000124649041611, 0.0000685496215625, 0.0001131160409390, 0.0000419907542591, 0.0000035218982282, 0.0000219211271097, 0.0000473186962874, 0.0000524142634057, 0.0000106823372552, 0.0000012237574387, 0.0000081742188269, 0.0000201940563214, 0.0000295722761753, 0.0000260539631956, 0.0000018994753518, 0.0000004932942990, 0.0000034477210351, 0.0000092294406360, 0.0000155143073237, 0.0000183386499818, 0.0000137922469695, -0.0000007075155453, 0.0000002223410995, 0.0000016045317657, 0.0000045260028113, 0.0000082655700107, 0.0000112393599417, 0.0000115758243113, 0.0000076838709956, -0.0000014020591745
    };
    
    double b3[] =
    {
        -0.0221835462237291, -0.0235386333304348, -0.0042911639711832, -0.0047431560217121, -0.0046577314472149, -0.0013089557502947, -0.0014343968205390, -0.0014978542575474, -0.0014329302934532, -0.0005167994164556, -0.0005573939123058, -0.0005921030407223, -0.0005978284714483, -0.0005673965369076, -0.0002409269302708, -0.0002561516055118, -0.0002723768586352, -0.0002815958312453, -0.0002792078156272, -0.0002646630240693, -0.0001261183503407, -0.0001325622938779, -0.0001403198638518, -0.0001464084186977, -0.0001485971591029, -0.0001459023931717, -0.0001384829633836, -0.0000719062974278, -0.0000749128468013, -0.0000788187384314, -0.0000824202283094, -0.0000846673495936, -0.0000849054394951, -0.0000829269749240, -0.0000788883333858
    };
    
    double b4[] =
    {
        0.0058366730167965, 0.0070452306758401, 0.0006914295465364, 0.0010322294603561, 0.0010057563135650, 0.0001394203795507, 0.0002309706405978, 0.0002596611624417, 0.0002409588083156, 0.0000386949167221, 0.0000679154947896, 0.0000830199015202, 0.0000850120755064, 0.0000780125513602, 0.0000133034384660, 0.0000241813441339, 0.0000311573885555, 0.0000340233089866, 0.0000335167900637, 0.0000307571022927, 0.0000053305073331, 0.0000099143129290, 0.0000132296989826, 0.0000150959309402, 0.0000156304390748, 0.0000151274875147, 0.0000139320508803, 0.0000023959090314, 0.0000045285807761, 0.0000061918979830, 0.0000072894226381, 0.0000078251853305, 0.0000078772667984, 0.0000075606242809, 0.0000069956215270
    };
    
    
    
    const vector<double> Tlm_real = Tlm(Omega*Hreal);
    
    /** Pre-computed psi */
    double psi[] = {0.9227843350984671394, 0.9227843350984671394,
        1.256117668431800473, 1.256117668431800473, 1.256117668431800473,
        1.506117668431800473, 1.506117668431800473, 1.506117668431800473, 1.506117668431800473,
        1.706117668431800473, 1.706117668431800473, 1.706117668431800473, 1.706117668431800473, 1.706117668431800473,
        1.872784335098467139, 1.872784335098467139, 1.872784335098467139, 1.872784335098467139, 1.872784335098467139, 1.872784335098467139,
        2.015641477955609997, 2.015641477955609997, 2.015641477955609997, 2.015641477955609997, 2.015641477955609997, 2.015641477955609997, 2.015641477955609997,
        2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997};
    
    double k;
    for (int i=kmax; i--;)
    {
        k  = M[i] * Omega;
        x  = k * Hreal; //hathatk
        x2 = x * x;
        x3 = x2 * x;
        x4 = x3 * x;
        x5 = x4 * x;
        // Why no linear term in the current implementation?
        
        //num_ang[i] = (1. + a1[i]*x2 + a2[i]*x3)/(1. + a3[i]*x2 + a4[i]*x3 + a5[i]*x4); //A.N. fit
        num_ang[i]    = 1. + b1[i]*x2 + b2[i]*x3 + b3[i]*x4 + b4[i]*x5; //P.F. fit
        
        ratio_rad     = Tlm_real[i];
        ratio_ang     = - 2. * psi[i] * x * num_ang[i];
        tlm_rad       = ratio_rad;
        tlm_ang       = ratio_ang + 2. * x * log(2. * k * bphys);
        
        tlm[i].dat[0] = tlm_rad;
        tlm[i].dat[1] = tlm_ang;
    }
    
    return tlm;
}

/* factorial */
double fact(int n)
{
    double f[] = {1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880.,
        3628800., 39916800., 479001600., 6227020800., 87178291200.};
    if (n < 0){
        fprintf(stderr, " computing a negative factorial. \n");
        exit(1);
    } else if (n <= 14){
        return f[n];
    } else {
        return n*fact(n-1);
    }
}

/* Wigner d-function */
double wigner_d_function(int l, int m, int s, double i)
{
    double dWig = 0.;
    
    double costheta = cos(i*0.5);
    double sintheta = sin(i*0.5);
    
    int ki = std::max( 0  , m-s );
    int kf = std::min( l+m, l-s );
    
    for( int k = ki; k <= kf; k++ ){
        dWig +=
        ( pow(-1.,k) * pow(costheta,2*l+m-s-2*k) * pow(sintheta,2*k+s-m) )/
        ( fact(k) * fact(l+m-k) * fact(l-s-k) * fact(s-m+k) );
    }
    
    return (sqrt(fact(l+m) * fact(l-m) * fact(l+s) * fact(l-s)) * dWig);
}

/* spin-weighted spherical harmonic */
void spinsphericalharm(double *rY, double *iY, int s, int l, int m, double phi, double i)
{
    /* Following the Ref.: https://arxiv.org/pdf/0709.0093.pdf */
    if ((l<0) || (m<-l) || (m>l)) {
        fprintf(stderr, " wrong l (%d) or m (%d) inside spinspharmY\n", l, m);
        exit(1);
    }
    
    double c = pow(-1.,-s) * sqrt( (2.*l+1.)/(4.*M_PI) );
    double dWigner = c * wigner_d_function(l,m,-s,i);
    
    *rY = cos((double)(m)*phi) * dWigner;
    *iY = sin((double)(m)*phi) * dWigner;
}

double fLR(double r, void *params)
{
    
    vector<double> metric=Metric(r, params,true);
    
    double A  = metric[0];
    double dA = metric[2]; /**derivative w.r.t. u*/
    double f  = A +(0.5*dA)/r;
    
    return f;
}
/** Takes nu as input */
double AdiabLR(void *params)
{
    
    int status;
    int iter = 0, max_iter = 200;
    const gsl_root_fsolver_type *T;
    double rLR;
    double x_lo = 0.1, x_hi = 15.;
    
    gsl_root_fsolver *s;
    gsl_function F;
    F.function = &fLR;
    F.params = params;
    T = gsl_root_fsolver_bisection;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
    
    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        rLR    = gsl_root_fsolver_root (s);
        x_lo   = gsl_root_fsolver_x_lower (s);
        x_hi   = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi, 0, 0.000000000000001);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (s);
    
    return rLR;
}

vector<double> s_D1(vector<double> f, vector<double> x, int Nmax)
{
    /* Computes the first derivative of the function. Centered but at the edges. USAGE: df = EOB_D1(f,x) */
    int Nmin = 0;
    
    vector<double> df(Nmax+1);
    for(int i=2;i<=Nmax-2;i++)
    {
        df[i] = 1./3.*(8.*f[1+i] - f[2+i] - 8.*f[i-1] + f[i-2])/(x[2+i]-x[i-2]);
    }
    
    // 4th order boundaries
    df[0] = (-24./17.*f[Nmin] + 59./34.*f[Nmin+1] - 4./17.*f[Nmin+2] - 3./34.*f[Nmin+3])/(x[Nmin+1]-x[Nmin]);
    df[1] = (-1./2.*f[Nmin] + 1./2.*f[Nmin+2] )/(x[Nmin+1]-x[Nmin]);
    
    df[Nmax]   = -(-24./17.*f[Nmax] + 59./34.*f[Nmax-1] - 4./17.*f[Nmax-2] - 3./34.*f[Nmax-3])/(x[Nmax]-x[Nmax-1]);
    df[Nmax-1] = -(-1./2.*f[Nmax] + 1./2.*f[Nmax-2])/(x[Nmax]-x[Nmax-1]);
    
    return df;
}

/** Sets the dynamics controlling flags to their default value */
void SetDefaultFlagsValues(TEOBResumFlags *flags)
{
    flags->spin             = 0;
    flags->tidal            = 0;
    flags->RWZ              = 0;
    flags->speedy           = 1;
    flags->dynamics         = 0;
    flags->multipoles       = 0;
    flags->Yagi_fits        = 0;
    flags->geometric_units  = 0;
    flags->solver_scheme    = 0;
    flags->set              = 0;
}

double logQ(double x)
{
    /** logQ-vs-log(lambda) fit of Table I of Yunes-Yagi
     here x = log(lambda) and the output is the log of the coefficient
     that describes the quadrupole deformation due to spin. */
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

double Yagi13_fit_barlamdel(double barlam2, int ell)
{
    /*
     Yagi 2013 fits for multipolar
     $\bar{\lambda}_\ell$ = 2 k_\ell/(C^{2\ell+1} (2\ell-1)!!)$
     Eq.(10),(61); Tab.I; Fig.8 http://arxiv.org/abs/1311.0872
     */
    double lnx = log(barlam2);
    double coeffs[5] = {0.0};
    if (ell == 3)
    {
        coeffs[0] = 2.52e-5;
        coeffs[1] = -1.31e-3;
        coeffs[2] = 2.51e-2;
        coeffs[3] = 1.18;
        coeffs[4] = -1.15;
    }
    else if (ell == 4)
    {
        coeffs[0] = 2.8e-5;
        coeffs[1] =-1.81e-3;
        coeffs[2] =3.95e-2;
        coeffs[3] =1.43;
        coeffs[4] =-2.45;
    }
    else return 0.0;
    
    double lny = coeffs[0]*lnx*lnx*lnx*lnx+coeffs[1]*lnx*lnx*lnx+coeffs[2]*lnx*lnx+coeffs[3]*lnx+coeffs[4];
    return exp(lny);
}

double radius0(double M, double f_start)
{
    double x = (M*f_start*MSUN_S*2.*M_PI)/2.;
    return cbrt(1/(x*x));
}

double time_units_conversion(double M, double dt)
{
    return dt/(M*MSUN_S);
}

TEOBResumParams read_config(char *fname)
{

    TEOBResumParams params;
    
    SetDefaultFlagsValues(&params.flags);
    
    string param_name;
    double param_value;
    
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

	if(param_name=="Mtot") {
	  params.mtot = param_value;
	}
	if(param_name=="iota") {
	  params.iota = param_value;
	}
	if(param_name=="psi") {
	  params.psi = param_value;
	}
	if(param_name=="distance") {
	  params.psi = param_value;
	}
	if(param_name=="q") {
	  params.q = param_value;
	}
	if(param_name=="chi1") {
	  params.chi1 = param_value;
	}
	if(param_name=="chi2") {
	  params.chi2 = param_value;
	}
	// if(param_name=="r0") {
	//   params.r0 = param_value;
	// }
	if(param_name=="f_min") {
	  params.f_min = param_value;
	}
	if(param_name=="tidal") {
	  params.flags.tidal = param_value;
	}
	if(param_name=="RWZ") {
	  params.flags.RWZ = param_value;
	}
	if(param_name=="speedy") {
	  params.flags.speedy = param_value;
	}
	if(param_name=="dynamics") {
	  params.flags.dynamics = param_value;
	}
	if(param_name=="Yagi_fit") {
	  params.flags.Yagi_fits = param_value;
	}
	if(param_name=="multipoles") {
	  params.flags.multipoles = param_value;
	}
	if(param_name=="lm") {
	  params.lm = param_value;
	}
	if(param_name=="dt") {
	  params.dt = param_value;
	}
	if(param_name=="solver_scheme") {
	  params.flags.solver_scheme = param_value;
	}
	if(param_name=="LambdaAl2") {
	  params.LambdaAl2 = param_value;
	}
	if(param_name=="LambdaBl2") {
	  params.LambdaBl2 = param_value;
	}
	if(param_name=="LambdaAl3") {
	  params.LambdaAl3 = param_value;
	}
	if(param_name=="LambdaBl3") {
	  params.LambdaBl3 = param_value;
	}
	if(param_name=="LambdaAl4") {
	  params.LambdaAl4 = param_value;
	}
	if(param_name=="LambdaBl4") {
	  params.LambdaBl4 = param_value;
	}
	if(param_name=="geometric_units") {
	  params.flags.geometric_units = param_value;
	}

        i++;
    }
    fin.close();
    
    // calculate parameters
    
    double chi1 = params.chi1;
    double chi2 = params.chi2;
    double q = params.q;
    if (params.flags.Yagi_fits == 1)
    {
        params.LambdaAl3 = Yagi13_fit_barlamdel(params.LambdaAl2, 3);
        params.LambdaBl3 = Yagi13_fit_barlamdel(params.LambdaBl2, 3);
        params.LambdaAl4 = Yagi13_fit_barlamdel(params.LambdaAl2, 4);
        params.LambdaBl4 = Yagi13_fit_barlamdel(params.LambdaBl2, 4);
    }
    
    /** Override spin settings if spins are given in input */
    if (chi1 != .0 || chi2 != .0) params.flags.spin = 1;
    
    double nu = q/((q+1.)*(q+1.));
    params.nu = nu;
    
    double X1 = 0.5*(1.+sqrt(1.-4.*nu));
    double X2 = 1. - X1;
    double XA = X1; // a different notation used in tidal part, keep here for simplicity
    double XB = X2;
    params.X1 = X1;
    params.X2 = X2;
    
    double S1 = params.X1 * params.X1 * params.chi1;
    double S2 = params.X2 * params.X2 * params.chi2;
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
    params.cN3LO = c3_fit_global(nu,chi1,chi2,X1,X2,a1,a2,params.flags.tidal);
    
    
    // tidal params
    double LambdaAl2 = params.LambdaAl2;
    double LambdaAl3 = params.LambdaAl3;
    double LambdaAl4 = params.LambdaAl4;
    
    double LambdaBl2 = params.LambdaBl2;
    double LambdaBl3 = params.LambdaBl3;
    double LambdaBl4 = params.LambdaBl4;
    
    /** Computing the tidal coupling constants */
    
    double kapA2 = 3.   * LambdaAl2 * XA*XA*XA*XA*XA / q; //Note: kap stands for kappa; see eqn(1) of REF
    double kapA3 = 15.  * LambdaAl3 * XA*XA*XA*XA*XA*XA*XA / q;
    double kapA4 = 105. * LambdaAl4 * XA*XA*XA*XA*XA*XA*XA*XA*XA / q;
    
    double kapB2 = 3.   * LambdaBl2 * XB*XB*XB*XB*XB * q;
    double kapB3 = 15.  * LambdaBl3 * XB*XB*XB*XB*XB*XB*XB * q;
    double kapB4 = 105. * LambdaBl4 * XB*XB*XB*XB*XB*XB*XB*XB*XB * q;
    
    double kapT2 = kapA2 + kapB2;
    double kapT3 = kapA3 + kapB3;
    double kapT4 = kapA4 + kapB4;
    
    params.kappaAl2 = kapA2;
    params.kappaAl3 = kapA3;
    params.kappaAl4 = kapA4;
    
    params.kappaBl2 = kapB2;
    params.kappaBl3 = kapB3;
    params.kappaBl4 = kapB4;
    
    params.kappaTl2 = kapT2;
    params.kappaTl3 = kapT3;
    params.kappaTl4 = kapT4;
    
    //-----------------------------------------------------------------------------------
    // Definition of the conservative tidal coefficients \bar{\alpha}_n^{(\ell)}, Eq.(37)
    // of Damour&Nagar, PRD 81, 084016 (2010)
    //-----------------------------------------------------------------------------------
    params.bar_alph2_1 = (5./2.*XA*kapA2 + 5./2.*XB*kapB2)/kapT2;
    params.bar_alph2_2 = ((3.+XA/8.+ 337./28.*XA*XA)*kapA2 + (3.+XB/8.+ 337./28.*XB*XB)*kapB2)/kapT2;
    params.bar_alph3_1 = ((-2.+15./2.*XA)*kapA3 + (-2.+15./2.*XB)*kapB3)/kapT3;
    params.bar_alph3_2 = ((8./3.-311./24.*XA+110./3.*XA*XA)*kapA3 + (8./3.-311./24.*XB+110./3.*XB*XB)*kapB3)/kapT3;
    
    
    //double lambda1 = params.LambdaAl2;
    //double lambda2 = params.LambdaBl2;
    double logC_Q1 = logQ(log(LambdaAl2));
    double logC_Q2 = logQ(log(LambdaBl2));
    double C_Q1    = exp(logC_Q1);
    double C_Q2    = exp(logC_Q2);
    params.C_Q1 = C_Q1;
    params.C_Q2 = C_Q2;
    
    return params;
}

void CopyTEOBResumSFlags(TEOBResumFlags *out, TEOBResumFlags *in)
{
    out->spin  = in->spin;
    out->tidal = in->tidal;
    out->RWZ   = in->RWZ;
    out->speedy = in->speedy;
    out->dynamics = in->dynamics;
    out->multipoles = in->multipoles;
    out->Yagi_fits  = in->Yagi_fits;
    out->geometric_units = in->geometric_units;
    out->solver_scheme = in->solver_scheme;
    out->set = in->set;
}

TEOBResumParams process_input_parameters(
                                         double m1,
                                         double m2,
                                         double chi1,
                                         double chi2,
                                         double f_min,
                                         double dt,
                                         double LambdaAl2,
                                         double LambdaBl2,
                                         double LambdaAl3,
                                         double LambdaBl3,
                                         double LambdaAl4,
                                         double LambdaBl4,
                                         TEOBResumFlags *flags   /** flags **/
                                         )
{
    TEOBResumParams params;
    
    CopyTEOBResumSFlags(&params.flags,flags);

    double mtot = m1+m2;
    double q = m1/m2;

    params.mtot = mtot;
    params.q = q;
    
    if (params.flags.geometric_units==0)
    {
      // input given in physical units, 
      // rescale to geometric units and mass rescaled quantities
      // compute r0 from the initial GW frequency in Hz
        params.dt = time_units_conversion(mtot, dt);
        params.r0 = radius0(mtot, f_min);
    }
    else
    {
      // input given in geometric units, 
      // rescale to geometric units and mass rescaled quantities
      // compute r0 from the initial GW frequency in geometric units and mass rescaled
        params.dt = dt;
        params.r0 = pow(f_min*M_PI, -2./3.);
    }

    printf(" dt = %e r0 = %e\n",params.dt, params.r0);
    
    params.chi1 = chi1;
    params.chi2 = chi2;
    
    if (params.flags.tidal == 1 && params.flags.Yagi_fits==1)
    {
        LambdaAl3 = Yagi13_fit_barlamdel(LambdaAl2, 3);
        LambdaBl3 = Yagi13_fit_barlamdel(LambdaBl2, 3);
        LambdaAl4 = Yagi13_fit_barlamdel(LambdaAl2, 4);
        LambdaBl4 = Yagi13_fit_barlamdel(LambdaBl2, 4);
    }
    
    /** Override spin settings if spins are given in input */
    if (chi1 != .0 || chi2 != .0) params.flags.spin = 1;
    
    double nu = q/((q+1.)*(q+1.));
    params.nu = nu;
    
    double X1 = 0.5*(1.+sqrt(1.-4.*nu));
    double X2 = 1. - X1;
    double XA = X1; // a different notation used in tidal part, keep here for simplicity
    double XB = X2;
    params.X1   = X1;
    params.X2   = X2;
    
    double S1 = params.X1*params.X1 * params.chi1;
    double S2 = params.X2*params.X2 * params.chi2;
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
    
    params.cN3LO = c3_fit_global(nu,chi1,chi2,X1,X2,a1,a2,params.flags.tidal);
    
    /** Computing the tidal coupling constants */
    
    double kapA2 = 3.   * LambdaAl2 * XA*XA*XA*XA*XA / q; //Note: kap stands for kappa; see eqn(1) of REF
    double kapA3 = 15.  * LambdaAl3 * XA*XA*XA*XA*XA*XA*XA / q;
    double kapA4 = 105. * LambdaAl4 * XA*XA*XA*XA*XA*XA*XA*XA*XA / q;
    
    double kapB2 = 3.   * LambdaBl2 * XB*XB*XB*XB*XB * q;
    double kapB3 = 15.  * LambdaBl3 * XB*XB*XB*XB*XB*XB*XB * q;
    double kapB4 = 105. * LambdaBl4 * XB*XB*XB*XB*XB*XB*XB*XB*XB * q;
    
    double kapT2 = kapA2 + kapB2;
    double kapT3 = kapA3 + kapB3;
    double kapT4 = kapA4 + kapB4;
    
    params.kappaAl2 = kapA2;
    params.kappaAl3 = kapA3;
    params.kappaAl4 = kapA4;
    
    params.kappaBl2 = kapB2;
    params.kappaBl3 = kapB3;
    params.kappaBl4 = kapB4;
    
    params.kappaTl2 = kapT2;
    params.kappaTl3 = kapT3;
    params.kappaTl4 = kapT4;
    
    //-----------------------------------------------------------------------------------
    // Definition of the conservative tidal coefficients \bar{\alpha}_n^{(\ell)}, Eq.(37)
    // of Damour&Nagar, PRD 81, 084016 (2010)
    //-----------------------------------------------------------------------------------
    params.bar_alph2_1 = (5./2.*XA*kapA2 + 5./2.*XB*kapB2)/kapT2;
    params.bar_alph2_2 = ((3.+XA/8.+ 337./28.*XA*XA)*kapA2 + (3.+XB/8.+ 337./28.*XB*XB)*kapB2)/kapT2;
    params.bar_alph3_1 = ((-2.+15./2.*XA)*kapA3 + (-2.+15./2.*XB)*kapB3)/kapT3;			     			   
    params.bar_alph3_2 = ((8./3.-311./24.*XA+110./3.*XA*XA)*kapA3 + (8./3.-311./24.*XB+110./3.*XB*XB)*kapB3)/kapT3;
    
    double logC_Q1 = logQ(log(LambdaAl2));
    double logC_Q2 = logQ(log(LambdaBl2));
    double C_Q1    = exp(logC_Q1);
    double C_Q2    = exp(logC_Q2);    
    params.C_Q1 = C_Q1;
    params.C_Q2 = C_Q2;
    
    return params;
}
