/**
 * This file is part of TEOBResumS
 *
 * Copyright (C) 2017-2018 See AUTHORS file
 *
 * TEOBResumS is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * TEOBResumS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/.       
 *
 */

#include "TEOBResumS.h"

/** post-adiabatic dynamics */
int eob_dyn_Npostadiabatic(Dynamics *dyn, double r0)
{
      
  /* Unpack values */
  const double nu    = dyn->nu;
  const double S     = dyn->S;
  const double Sstar = dyn->Sstar;
  const double chi1  = dyn->chi1;
  const double chi2  = dyn->chi2;
  const double X1    = dyn->X1;
  const double X2    = dyn->X2;
  const double c3    = dyn->cN3LO;
  const double aK2   = dyn->aK2;
  const double a1    = dyn->a1;
  const double a2    = dyn->a2;
  const double C_Q1  = dyn->C_Q1;
  const double C_Q2  = dyn->C_Q2;
  const int usetidal = dyn->use_tidal;
  const int usespins = dyn->use_spins;

  const int size = dyn->size;
  const int Npa = par_get_i("postadiabatic_dynamics_N");

  /* Mem for integrant of t and phi */
  double *dHdPrs = (double*)malloc(size * sizeof (double));
  double *dHdPphi_o_dHdPrs = (double*)malloc(size * sizeof (double));
  
  /** Compute radius of inflection point of Pr* */
  //TODO
  double ri = 0.;
    
  /** Build a uniform grid */
  const double dr = (r0 - ri)/(size-1);
  for (int i = 0; i < size; i++) 
    dyn->data[EOB_RAD][i] = r0 - i*dr;

  /** Compute Pphi and Pr* by iteration */
  for (int n = 0; n < Npa; n++) {
    for (int i = 0; i < size; i++) {
      
      // TODO PA iteration 
      // maximize use of routines already coded e.g. flux routines, hamiltonian, etc
      // check code in *Dynamics.c and *InitialConditions.c
      // store final result into dyn->data[<variable_index>][i]
      // note Dynamics structure has 'storage' arrays dyn->data and point-wise vars, 
      /* 
      dyn->data[EOB_RAD][i]    = dyn->r;
      dyn->data[EOB_PHI][i]    = dyn->phi;
      dyn->data[EOB_PPHI][i]   = dyn->pphi;
      dyn->data[EOB_MOMG][i]   = dyn->Omg;
      dyn->data[EOB_DDOTR][i]  = dyn->ddotr;
      dyn->data[EOB_PRSTAR][i] = dyn->prstar;
      dyn->data[EOB_OMGORB][i] = dyn->Omg_orb;
      */
      
    }
  }
  
  /** Compute time */
  cumint3(dHdPrs, dyn->data[EOB_RAD], size, dyn->time);

  /* Set last value for evolution */
  dyn->t = dyn->time[size-1];
  
  /** Compute orbital phase */
  cumint3(dHdPphi_o_dHdPrs, dyn->data[EOB_RAD], size, dyn->data[EOB_PHI]);

  /* Free memory */
  free(dHdPrs);
  free(dHdPphi_o_dHdPrs);

  return OK;
}
