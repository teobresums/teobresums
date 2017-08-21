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

#ifndef _input_struc_h
#define _input_struc_h

#include <string> //  std::string outputdir; requires this, not sure why

struct input
{         /* declaration of structure for input parameters */
    int   lm;
    int solver_scheme;
    double q;
    double nu;
    double r0;
    double dt;
    double rLR;
    double chi1;
    double chi2;
    double S1;
    double S2;
    double S;
    double Sstar;
    double X1;
    double X2;
    double a1;
    double a2;
    double aK;
    double aK2;
    double cN3LO;
    double LambdaAl2;
    double LambdaAl3;
    double LambdaAl4;
    double LambdaBl2; 
    double LambdaBl3; 
    double LambdaBl4;
    double C_Q1;
    double C_Q2;
    bool NQC;
    bool tidal;
    bool RWZ;
    bool speedy;
    bool dynamics;
  bool waveform;
  bool spin;
  //char outputdir[256];
  std::string outputdir; // = "data";
};

struct input_flags
{
    bool NQC;
    bool spin;
    bool tidal;
    bool RWZ;
    bool speedy;
    bool dynamics;
    bool waveform;
};

#endif /* _input_struc_h */
