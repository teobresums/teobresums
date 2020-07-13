/* program to test spin dynamics */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "TEOBResumS.h"

int main (int argc, char* argv[])
{
  // Defaults
  double f0 = 10; // initial frequency
  double Sx0 = 0; // x-comp spin
  double Sy0 = 0; // y-comp spin
  double Sz0 = 0.5; // z-comp spin
  double ti = 10; // time to interp

  // Read-in pars
  int opt;
  while ((opt = getopt(argc, argv, "fxyzt")) != -1) {
    switch (opt) {
    case 'f': f0 = (double)optarg; break;
    case 'x': Sx0 = (double)optarg; break;
    case 'y': Sy0 = (double)optarg; break;
    case 'x': Sz0 = (double)optarg; break;
    case 't': ti = (double)optarg; break;
    default:      
      fprintf(stderr, "Usage: %s [-fxyzt] <value>\n", argv[0]);
      abort();
      //exit(1);
    }
  }
  
  DynamicsSpin **dp;
  int size = 100; // can vary
  
  // Alloc mem
  DynamicsSpin_alloc(&dp, size); 

  // Init
  dp->f0 = f0;
  dp->Sx0 = Sx0;
  dp->Sy0 = Sy0;
  dp->Sz0 = Sz0;
  
  // Solve 
  if (eob_spin_dyn(dp)) {
    fprintf("ERROR\n");
    return ERROR_ODEINT;
  }

  // Test interp
  double alpha, beta, gamma;
  eob_spin_dyn_abc_interp(dp, ti, &alpha,&beta,&gamma);

  fprintf("angles at %e\n alpha = %e\n beta= %e\n gamma=%e\n",ti,alpha,beta,gamma);
  
  // Output
  DynamicsSpin_output(dp);

  // Free mem
  DynamicsSpin_free(dp);
  return OK;
}


  
