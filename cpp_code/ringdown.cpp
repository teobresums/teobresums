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
#include "cmath"
#include <vector>

#include "TEOBResum.h"

using namespace::std;


int ringdown(TEOBResumParams params, vector<vector<double> > &t_vec, vector<double> Omega_vec, vector<vector<double> > &hlm_rad, vector<vector<double> > &hlm_phase){
    
    //EOBhlm Compute the multipolar resummed waveform.
    //
    //   [hlm,phi,omglm,domglm,d2omglm, psilm,Alm,dpsilm] = ...
    //   EOBhlm(nu,t,phi,r,pph,prstar, Omega,E,Heff, EOBopt, EOBmet)
    //   return (complex) multipolar wave, phase, frequency, and derivatives,
    //   and the RWZ normalized (complex) wave, its amplitude and (complex)
    //   derivative.
    //
    //   WAVE = EOBhlm( .... ) return a structure with the wave
    //
    
    //NOTE: before the Omega_vec had a Mbh multiplied onto it!!!!!!!!!!!!!!!!!!!!!!
    
    
    double dt   = params.dt;
    double nu   = params.nu;
    double Mbh  = params.Mbh;
    double chi1 = params.chi1;

    
    long int pk_index = Omega_vec.size()-1;
    double Omega_pk   = Omega_vec[pk_index];
    long int i        = pk_index-1;
    while (Omega_vec[i] > Omega_pk) {
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
    double DeltaT_nqc = 0.;
    vector<double> tmrg(35);
    vector<double> tmatch(35);
    
    //compute true peak by interpolation from the grid
    tOmg_pk  = interpolate(dt, Omega_pk_grid);
    tOmg_pk *= 1./Mbh;
    
    //calculate tmatch
    double xnu =(1.-4.*nu);
        
    vector<double> dtmrg(2);
    int k21 = 0;
    int k22 = 1;
    int k33 = 4;
    

    /** time-shift needed when the largest object is highly spinning*/
    if (chi1 >= 0.8498)
    {
        /* Interpolating fit for Deltat_NQC. See Eq.(21) of arXiv:1506.08457 */
        DeltaT_nqc = dtnqc_fit(chi1,0.8498);
    }
    else
    {
        DeltaT_nqc = 1.;
    }

    tmrg[k22] = tOmg_pk-(DeltaT_nqc + 2)/Mbh;     //t_max(A22) => MERGER

    /** nonspinning case - old */
    /*tmrg[k22]  = tOmg_pk-3./Mbh; */             // t_max(A22) => MERGER


    /** only nonspinning case */
    dtmrg[0]   = 5.70364338 + 1.85804796*xnu  + 4.0332262*xnu*xnu; //k21
    dtmrg[1]   = 4.29550934 - 0.85938*xnu;                         //k33
    tmrg[k21]  = tmrg[k22] + dtmrg[0]/Mbh;     // t_max(A21) => peak of 21 mode
    tmrg[k33]  = tmrg[k22] + dtmrg[1]/Mbh;     // t_max(A33) => peak of 33 mode


    /* postmerger-ringdown matching time */    
    tmatch = tmrg;
    
    int kmax = 35;
    vector<gsl_complex> sigma(kmax); //move this to another place
    
    for (int i=35; i--; ) {
        sigma[i].dat[0]=0.;
        sigma[i].dat[1]=0.;
        tmatch[i] += 2./Mbh;
    }
    
    vector<double> a1(35);
    vector<double> a2(35);
    vector<double> a3(35);
    vector<double> a4(35);
    vector<double> b1(35);
    vector<double> b2(35);
    vector<double> b3(35);
    vector<double> b4(35);
    
    QNMHybridFitCab(params,a1,a2,a3,a4,b1,b2,b3,b4,sigma);
    
    /*deleting data points up to tmatch (starting from the back)*/
    vector<long> I(35);
    vector<long> Size(35);
    for (int k = 35; k--; ) {
        i = t_vec[k].size()-1;
        Size[k] = t_vec[k].size();
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
    
    //Calculate deltaphi
    vector<gsl_complex> psi(35);
    vector<double> Deltaphi(35);
    for (int k=35; k--; ) {
        double x    = t_vec[k][I[k]]/Mbh-tmrg[k];
        psi[k]      = ringdown_match(x, k, a1, a2, a3, a4, b1, b2, b3, b4, sigma);
        Deltaphi[k] = psi[k].dat[1] - hlm_phase[k][I[k]];
    }
    
    /** add 200 points of ringdown attachment */
    /** if we select the number of points in each mode to be such that the number of output points
        is the same, we do not need any interpolation, since the vectors are defined on the same time grid. */
    int Nringdown = 500;
    for (int k=35; k--; ) {
        double t = t_vec[k][I[k]];
        int n_removed = Size[k]-I[k];
        switch (k) {
            case 0:
                for (int j=0; j < Nringdown+n_removed; j++) {
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
                for (int j=0; j < Nringdown+n_removed; j++) {
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
                for (int j=0; j < Nringdown+n_removed; j++) {
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
                for (int j=0; j < Nringdown+n_removed-1; j++) {
                    hlm_rad[k].push_back(0.0);
                    hlm_phase[k].push_back(0.0);
                    t_vec[k].push_back(t);
                    t += dt;
                }
                break;
        }
    }

    return 0;
}

