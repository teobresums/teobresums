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
  double m = 1; // mass
  double q = 1; // mass ratio
  double S1x = 0; // x-comp spin 1
  double S1y = 0; // y-comp spin
  double S1z = 0.5; // z-comp spin
  double S2x = 0; // x-comp spin 2
  double S2y = 0; // y-comp spin
  double S2z = 0.5; // z-comp spin
  double ti = 10; // time to interp
  double momg_stop = 0.1; // ODE stop max freq
  double dt = 1; // ODE dt
    
  // Read-in pars
  int opt;
  char opts[256] = "fmqxyzXYZtod";
  while ((opt = getopt(argc, argv, opts)) != -1) {
    switch (opt) {
    case 'f': f0 = (double)optarg; break;
    case 'm': m = (double)optarg; break;
    case 'q': q = (double)optarg; break;
    case 'x': S1x = (double)optarg; break;
    case 'y': S1y = (double)optarg; break;
    case 'x': S1z = (double)optarg; break;
    case 'X': S2x = (double)optarg; break;
    case 'Y': S2y = (double)optarg; break;
    case 'Y': S2z = (double)optarg; break;
    case 't': ti = (double)optarg; break;
    case 'o': momg_stop = (double)optarg; break;
    case 'd': dt = (double)optarg; break;
    default:      
      fprintf(stderr, "Usage: %s [-%s] <value>\n", argv[0],opts);
      abort();
      //exit(1);
    }
  }
  
  // Init parameters & set defaults 
  EOBParameters_alloc( &EOBPars );
  EOBParameters_defaults (dc, EOBPars);

  EOBPars->spin_odes_omg_stop = momg_stop; 
  EOBPars->spin_odes_dt = ; 

  EOBPars->initial_frequency = f0;

  EOBPars->q = q;
  EOBPars->nu = q_to_nu(q);

  double m1 = q*m/2.;
  double m2 = m-m1;
  
  EOBPars->chi1x = S1x/SQ(m1);
  EOBPars->chi1y = S1y/SQ(m1);
  EOBPars->chi1z = S1z/SQ(m1);
  EOBPars->chi2x = S2x/SQ(m2);
  EOBPars->chi2y = S2y/SQ(m2);
  EOBPars->chi2z = S2z/SQ(m2);
  
  // Alloc mem
  DynamicsSpin **dp;
  int size = 100; // can vary
  DynamicsSpin_alloc(&dp, size); 

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
  EOBParameters_free (EOBPars);
  return OK;
}


  
