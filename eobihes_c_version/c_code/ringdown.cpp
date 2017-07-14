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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <math.h>

#include "FDdrvt_omega.h"
#include "interpolator.h"
#include "QNMHybridFitCab.h"
#include "ringdown.h"
#include "ringdown_match.h"

//int ringdown(double nu,double q,int lm,double dt,double Mbh,vector<double> &t_vec,vector<double> &Omega_vec,vector<double> &hlm_rad,vector<double> &hlm_phase){

int ringdown(double nu, double q, double dt, double Mbh, vector<vector<double> > &t_vec, vector<double> Omega_vec, vector<vector<double> > &hlm_rad, vector<vector<double> > &hlm_phase){
    
    /*  =====================================================================
     	 EOBhlm Compute the multipolar resummed waveform.
      
        [hlm,phi,omglm,domglm,d2omglm, psilm,Alm,dpsilm] = ...
        EOBhlm(nu,t,phi,r,pph,prstar, Omega,E,Heff, EOBopt, EOBmet)
        return (complex) multipolar wave, phase, frequency, and derivatives,
        and the RWZ normalized (complex) wave, its amplitude and (complex)
        derivative.
     
        WAVE = EOBhlm( .... ) return a structure with the wave
     
     
         NOTE: before the Omega_vec had a Mbh multiplied onto it!
        ===================================================================== */
//    for (long j=Omega_vec.size(); j--; ) {
//        Omega_vec[j] = Mbh*Omega_vec[j];
//    }
    double nu2 = nu*nu;
    double nu3 = nu2*nu;
    long int pk_index = Omega_vec.size()-1;
    double Omega_pk   = Omega_vec[pk_index];
    long int i        = pk_index-1;
    while (Omega_vec[i] > Omega_pk)
    {
        pk_index = i;
        Omega_pk = Omega_vec[i];
        i--;
    }
    vector<gsl_complex> Omega_pk_grid(7);
    Omega_pk_grid[0].dat[0] = t_vec[1][pk_index-3];
    Omega_pk_grid[0].dat[1] = Omega_vec[pk_index-3];
    Omega_pk_grid[1].dat[0] = t_vec[1][pk_index-2];
    Omega_pk_grid[1].dat[1] = Omega_vec[pk_index-2];
    Omega_pk_grid[2].dat[0] = t_vec[1][pk_index-1];
    Omega_pk_grid[2].dat[1] = Omega_vec[pk_index-1];
    Omega_pk_grid[3].dat[0] = t_vec[1][pk_index];
    Omega_pk_grid[3].dat[1] = Omega_vec[pk_index];
    Omega_pk_grid[4].dat[0] = t_vec[1][pk_index+1];
    Omega_pk_grid[4].dat[1] = Omega_vec[pk_index+1];
    Omega_pk_grid[5].dat[0] = t_vec[1][pk_index+2];
    Omega_pk_grid[5].dat[1] = Omega_vec[pk_index+2];
    Omega_pk_grid[6].dat[0] = t_vec[1][pk_index+3];
    Omega_pk_grid[6].dat[1] = Omega_vec[pk_index+3];
    
    double tOmg_pk = 0.;
    vector<double> tmrg(35);
    vector<double> tmatch(35);
    
    //compute true peak by interpolation from the grid
    tOmg_pk  = interpolate(dt, Omega_pk_grid);
    tOmg_pk *= 1./Mbh;
    
    //claculate tmatch
    double xnu =(1.-4.*nu);
        
    vector<double> dtmrg(2);
    int k21 = 0;
    int k22 = 1;
    int k33 = 4;
    dtmrg[0] = 5.70364338 + 1.85804796*xnu  + 4.0332262*xnu*xnu; //k21
    dtmrg[1] = 4.29550934 - 0.85938*xnu; //k33
        
    tmrg[k22]  = tOmg_pk-3./Mbh;               // t_max(A22) => MERGER
    tmrg[k21]  = tmrg[k22] + dtmrg[0]/Mbh;     // t_max(A21) => peak of 21 mode
    tmrg[k33]  = tmrg[k22] + dtmrg[1]/Mbh;     // t_max(A33) => peak of 33 mode
        
    tmatch = tmrg;
    
    int kmax = 35;
    vector<gsl_complex> sigma(kmax); //move this to another place
    
    for (int i=35; i--; ) {
        sigma[i].dat[0]=0.;
        sigma[i].dat[1]=0.;
        tmatch[i] += 2./Mbh;
    }
    
    /*
    switch ((int)q) {
        case 1:
            sigma[1].dat[0] = 0.081279587;
            sigma[1].dat[1] = 0.52680957;
            break;
        case 2:
            sigma[1].dat[0] = 0.083200939;
            sigma[1].dat[1] = 0.5020357;
            break;
        case 3:
            sigma[1].dat[0] = 0.084967983;
            sigma[1].dat[1] = 0.47554712;
            break;
        case 4:
            sigma[1].dat[0] = 0.086036548;
            sigma[1].dat[1] = 0.45686417;
            break;
        case 6:
            sigma[1].dat[0] = 0.087149214;
            sigma[1].dat[1] = 0.43387484;
            break;
        case 7://7.1429
            sigma[1].dat[0] = 0.087496124;
            sigma[1].dat[1] = 0.42556176;
        case 8:
            sigma[1].dat[0] = 0.087685918;
            sigma[1].dat[1] = 0.42068061;
            break;
        case 9://9.9891:
            sigma[1].dat[0] = 0.087989127;
            sigma[1].dat[1] = 0.41225268;
            break;
        case 20:
            sigma[1].dat[0] = 0.088536933;
            sigma[1].dat[1] = 0.39399845;
            break;
        default:
            break;
    }
    */

    for (i=35; i--; ) {
        switch (i) {
            case 0:
                sigma[0].dat[0] = -0.208936*nu3-0.028103*nu2-0.005383*nu + 0.08896;
                sigma[0].dat[1] = 0.733477*nu3 + 0.188359*nu2 + 0.220659*nu + 0.37367;
                break;
            case 1:
                sigma[1].dat[0] = -0.364177*nu3 + 0.010951*nu2-0.010591*nu + 0.08896;
                sigma[1].dat[1] = 2.392808*nu3 + 0.051309*nu2 + 0.449425*nu + 0.37365;
                break;
            case 4:
                sigma[4].dat[0] = -0.319703*nu3-0.030076*nu2-0.009034*nu + 0.09270;
                sigma[4].dat[1] = 2.957425*nu3 + 0.178146*nu2 + 0.709560*nu + 0.59944;
                break;
            default:
                sigma[i].dat[0] = 0.;
                sigma[i].dat[1] = 0.;
                break;
        }
    }
    
    /*deleting data points up to tmatch (starting from the back)*/
    vector<long> I(35);
    for (int k = 35; k--; ) {
        i = t_vec[k].size()-1;
        switch (k) {
            case 0:
                while (t_vec[k][i]/Mbh>tmatch[k]) {
                    t_vec[k].pop_back();
                    hlm_rad[k].pop_back();
                    hlm_phase[k].pop_back();
                    i--;
                }
                i++;
                break;
            case 1:
                while (t_vec[k][i]/Mbh>tmatch[k]) {
                    t_vec[k].pop_back();
                    hlm_rad[k].pop_back();
                    hlm_phase[k].pop_back();
                    i--;
                }
                i++;
                break;
            case 4:
                while (t_vec[k][i]/Mbh>tmatch[k]) {
                    t_vec[k].pop_back();
                    hlm_rad[k].pop_back();
                    hlm_phase[k].pop_back();
                    i--;
                }
                i++;
                break;
            default:
                break;
        }
        I[k] = i;
    }
    
    vector<double> a1(35);
    vector<double> a2(35);
    vector<double> a3(35);
    vector<double> a4(35);
    vector<double> b1(35);
    vector<double> b2(35);
    vector<double> b3(35);
    vector<double> b4(35);
    QNMHybridFitCab(nu,a1,a2,a3,a4,b1,b2,b3,b4);
    
    //Calculate deltaphi
    vector<gsl_complex> psi(35);
    vector<double> Deltaphi(35);
    for (int k=35; k--; ) {
        double x    = t_vec[k][I[k]]/Mbh-tmrg[k];
        psi[k]      = ringdown_match(x, k, a1, a2, a3, a4, b1, b2, b3, b4, sigma);
        Deltaphi[k] = psi[k].dat[1] - hlm_phase[k][I[k]];
    }
    
    /*add 200 points of ringdown attachment*/
    for (int k=35; k--; ) {
        double t = t_vec[k][I[k]];
        switch (k) {
            case 0:
                for (int j=0; j < 200; j++) {
                    double x = t/Mbh-tmrg[k];
                    psi[k] = ringdown_match(x, k, a1, a2, a3, a4, b1, b2, b3, b4, sigma);
                    
                    psi[k].dat[1]  = psi[k].dat[1] - Deltaphi[k];
                    hlm_rad[k].push_back(psi[k].dat[0]);
                    hlm_phase[k].push_back(psi[k].dat[1]);
                    t_vec[k].push_back(t);
                    t += dt;
                }
                break;
            case 1:
                for (int j=0; j < 200; j++) {
                    double x = t/Mbh-tmrg[k];
                    psi[k] = ringdown_match(x, k, a1, a2, a3, a4, b1, b2, b3, b4, sigma);
                    
                    psi[k].dat[1]  = psi[k].dat[1] - Deltaphi[k];
                    hlm_rad[k].push_back(psi[k].dat[0]);
                    hlm_phase[k].push_back(psi[k].dat[1]);
                    t_vec[k].push_back(t);
                    t += dt;
                }
                break;
            case 4:
                for (int j=0; j < 200; j++) {
                    double x = t/Mbh-tmrg[k];
                    psi[k] = ringdown_match(x, k, a1, a2, a3, a4, b1, b2, b3, b4, sigma);
                    
                    psi[k].dat[1]  = psi[k].dat[1] - Deltaphi[k];
                    hlm_rad[k].push_back(psi[k].dat[0]);
                    hlm_phase[k].push_back(psi[k].dat[1]);
                    t_vec[k].push_back(t);
                    t += dt;
                }
                break;
            default:
                break;
        }
    }

    return 0;
}


//    //calculate deltaphi
//    double x=t_vec[i]/Mbh-tmrg[lm];
//    gsl_complex psi = ringdown_match(x, lm, a1, a2, a3, a4, b1, b2, b3, b4, sigma);
//    double Deltaphi  = psi.dat[1] - hlm_phase[i];

//double t = t_vec[i];
//vector<double> phase;
//    for (int j=0; j < 200; j++) {
//        x = t/Mbh-tmrg[lm];
//        psi = ringdown_match(x, lm, a1, a2, a3, a4, b1, b2, b3, b4, sigma);
//
//        psi.dat[1]  = psi.dat[1] - Deltaphi;
//        hlm_rad.push_back(psi.dat[0]);
//        hlm_phase.push_back(psi.dat[1]);
//        //phase.push_back(psi.dat[1]);
//        t_vec.push_back(t);
//        t += dt;
//    }

//calculating the MOmega
//    vector<double> dphase = FDdrvt_omega(hlm_phase, len(hlm_phase), dt);
//    Omega_vec={};
//    for (long int j=0; j < hlm_phase.size()-2; j++) {
//        Omega_vec.push_back( - dphase[j] );
//    }


//    vector<vector<double> > dphase(35);
//    //vector<vector<double> > Omega_vec(35);
//    for (int k=35; k--; ) {
//        dphase[k] = FDdrvt_omega(hlm_phase[k], len(hlm_phase[k]), dt);
//        //Omega_vec={};
//        for (long int j=0; j < hlm_phase[k].size()-2; j++) {
//            Omega_vec.push_back( - dphase[k][j] );
//        }    }

/* alpha is the real and omega the imaginary part
 +\alpha_{22}^{1}  &= -0.364177\nu^{3} + 0.010951\nu^{2}-0.010591\nu + 0.08896\\
 +\nonumber
 +\alpha_{21}^{1} & = -0.208936\nu^{3}-0.028103\nu^{2}-0.005383\nu + 0.08896\\
 +\nonumber
 +\alpha_{33}^{1} & = -0.319703\nu^{3}-0.030076\nu^{2}-0.009034\nu + 0.09270\\
 +\nonumber
 +\omega_{22}^{1} &=  2.392808\nu^{3} + 0.051309\nu^{2} + 0.449425\nu + 0.37365\\
 +\nonumber
 +\omega_{21}^{1}& =  0.733477\nu^{3} + 0.188359\nu^{2} + 0.220659\nu + 0.37367\\
 +\nonumber
 +\omega_{33}^{1}& = 2.957425\nu^{3} + 0.178146\nu^{2} + 0.709560\nu + 0.59944
 */
