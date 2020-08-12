/* program to test spin dynamics 
   uses geometric units */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

#include "../../src/TEOBResumS.h"

/** Need to redefine some globals 
    since TEOBResumS.c is not compiled */
const int LINDEX[KMAX] = {
    2,2,
    3,3,3,
    4,4,4,4,
    5,5,5,5,5,
    6,6,6,6,6,6,
    7,7,7,7,7,7,7,
    8,8,8,8,8,8,8,8};
const int MINDEX[KMAX] = {
    1,2,
    1,2,3,
    1,2,3,4,
    1,2,3,4,5,
    1,2,3,4,5,6,
    1,2,3,4,5,6,7,
    1,2,3,4,5,6,7,8};
const int KINDEX[9][9] = {  // l (m = 1 ...l)
  {-1,-1,-1,-1,-1,-1,-1,-1},// 0 
  {-1,-1,-1,-1,-1,-1,-1,-1},// 1 
  { 0, 1,-1,-1,-1,-1,-1,-1},// 2 (1 2)
  { 2, 3, 4,-1,-1,-1,-1,-1},// 3 (1 2 3)
  { 5, 6, 7, 8,-1,-1,-1,-1},// 4 (1 2 3 4)
  { 9,10,11,12,13,-1,-1,-1},// 5 (1 2 3 4 5)
  {14,15,16,17,18,19,-1,-1},// 6 (1 ... 6)
  {20,21,22,23,24,25,26,-1},// 7 (1 ... 7)
  {27,28,29,30,31,32,33,34},// 8 (1 ... 8)
};

NQCdata *NQC;


int main (int argc, char* argv[])
{
  // Defaults
  double f0 = 0.01; // initial frequency
  double m = 1; // mass
  double q = 1; // mass ratio
  double S1x = 0; // x-comp spin 1
  double S1y = 0; // y-comp spin
  double S1z = 0.125; // z-comp spin
  double S2x = 0; // x-comp spin 2
  double S2y = 0; // y-comp spin
  double S2z = 0.125; // z-comp spin
  double ti = 10; // time to interp
  double momg_stop = -1; // ODE stop max freq if >0 , otherwise calculate from NR mrg fit
  double dt = 1; // ODE dt
  int interpint = 0; // test integration in interp
  
  // Read-in pars
  int opt;
  char opts[256] = "f:m:q:x:y:z:X:Y:Z:t:o:d:i:h";
  while ((opt = getopt(argc, argv, opts)) != -1) {
    switch (opt) {
    case 'f': f0 = atof(optarg); break;
    case 'm': m = atof(optarg); break;
    case 'q': q = atof(optarg); break;
    case 'x': S1x = atof(optarg); break;
    case 'y': S1y = atof(optarg); break;
    case 'z': S1z = atof(optarg); break;
    case 'X': S2x = atof(optarg); break;
    case 'Y': S2y = atof(optarg); break;
    case 'Z': S2z = atof(optarg); break;
    case 't': ti = atof(optarg); break;
    case 'o': momg_stop = atof(optarg); break;
    case 'd': dt = atof(optarg); break;
    case 'i': interpint = atoi(optarg); break;
    case '?':
      fprintf(stderr, "Usage: %s -{%s}\n", argv[0],opts);
      return 1;
    default:      
      abort();
    }
  }
  
  // Init parameters & set defaults 
  EOBParameters_alloc( &EOBPars );
  EOBParameters_defaults (DEFAULT_PARS_BBH, EOBPars);

  EOBPars->M = m;
  EOBPars->spin_odes_omg_stop = momg_stop; 
  EOBPars->spin_odes_dt = dt; 
  EOBPars->initial_frequency = f0;
  EOBPars->q = q;
  EOBPars->nu = q_to_nu(q);

  double m1 = m/(1.+q);
  double m2 = m-m1;
  
  EOBPars->chi1x = S1x/SQ(m1);
  EOBPars->chi1y = S1y/SQ(m1);
  EOBPars->chi1z = S1z/SQ(m1);
  EOBPars->chi2x = S2x/SQ(m2);
  EOBPars->chi2y = S2y/SQ(m2);
  EOBPars->chi2z = S2z/SQ(m2);

  printf("Mass  : %e\n",m);
  printf("Spin 1: (Sx,Sy,Sz)=(%e,%e,%e)\n",S1x, S1y, S1z);
  printf("Spin 2: (Sx,Sy,Sz)=(%e,%e,%e)\n",S2x, S2y, S2z);
  
  EOBPars->spin_interp_integrate = interpint;
  
  // Alloc mem
  DynamicsSpin *dp;
  int size = 100; // can vary
  DynamicsSpin_alloc(&dp, size); 

  // Solve 
  if (eob_spin_dyn(dp)) {
    return ERROR_ODEINT;
  }

  printf("final time at %e (tstop = %e):\n\t momg = %e\n\t momg_stop = %e (%e)\n",
	 dp->t, dp->t_stop, dp->y[EOB_EVOLVE_SPIN_Momg],dp->omg_stop, momg_stop);
  
  // Test interp
  double alpha, beta, gamma;
  eob_spin_dyn_abc_interp(dp, ti, &alpha,&beta,&gamma, EOBPars->spin_interp_integrate);

  printf("interp. angles at %e:\n\t alpha = %e\n\t beta = %e\n\t gamma = %e\n",
	 ti,alpha,beta,gamma);
  
  // Output
  strcpy(EOBPars->output_dir, "./output/");
  system_mkdir(EOBPars->output_dir);
  DynamicsSpin_output(dp);

  // Free mem
  DynamicsSpin_free(dp);
  EOBParameters_free (EOBPars);
  return OK;
}


  
