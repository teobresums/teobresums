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

#include "TEOBResumS.h"

/** logQ-vs-log(lambda) fit of Table I of Yunes-Yagi
    here x = log(lambda) and the output is the log of the coefficient
    that describes the quadrupole deformation due to spin. */
double logQ(double x)
{
  const double ai = 0.194;
  const double bi = 0.0936;
  const double ci = 0.0474;
  const double di = -4.21e-3;
  const double ei = 1.23e-4;
  const double x2 = x*x;
  const double x3 = x*x2;
  const double x4 = x*x3;
  return ai + bi*x + ci*x2 + di*x3 + ei*x4;
}

/** Yagi 2013 fits for multipolar
    $\bar{\lambda}_\ell$ = 2 k_\ell/(C^{2\ell+1} (2\ell-1)!!)$
    Eq.(10),(61); Tab.I; Fig.8 http://arxiv.org/abs/1311.0872 */
double Yagi13_fit_barlamdel(double barlam2, int ell)
{
  double lnx = log(barlam2);
  double coeffs[5];
  if (ell == 3)
    {
      coeffs[0] = 2.52e-5;
      coeffs[1] = -1.31e-3;
      coeffs[2] = 2.51e-2;
      coeffs[3] = 1.18;
      coeffs[4] = -1.15;
    }
  if (ell == 4)
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

