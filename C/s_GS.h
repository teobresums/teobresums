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

#ifndef _s_GS_h
#define _s_GS_h

#include <stdbool.h>

//double c3_fit_global(void *params){
double c3_fit_global(double nu, double chi1, double chi2, double X1, double X2, double a1, double a2, bool tidal_flag);
void s_GS(double result[], double r, double rc, double drc_dr, double aK2, double prstar, double pph, double nu, double chi1, double chi2, double X1, double X2, double cN3LO);

#endif /* _s_GS_h */
