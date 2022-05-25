/*
Small C script to test the implementation of 
the analytical derivatives of the metric potentials
w.r.t. the real numerical values. 

TODO: add tides

*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include "../C/src/TEOBResumS.h"

/** Necessary definitions, or compiler will complain... */ 
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

int main (int argc, char* argv[]){

  //Default
  double q      = 2.;
  double chi1   = 0.9;
  double chi2   = 0.3;

  //radial grid
  double rmin = 2.;
  double rmax = 20.;
  double dr   = 0.1;
  int N = floor((rmax-rmin)/dr) +1;
  
  double r[N];
  for (int i=0; i<N; i++)
    r[i] = rmin + i*dr; 


  // Parser
  int opt;
  char opts[256] = "q:z:Z:h";
  while ((opt = getopt(argc, argv, opts)) != -1) {
    switch (opt) {
    case 'q': q = atof(optarg); break;
    case 'z': chi1 = atof(optarg); break;
    case 'Z': chi2 = atof(optarg); break;
    case '?':
      fprintf(stderr, "Usage: %s -{%s}\n", argv[0],opts);
      return 1;
    default:      
      abort();
    }
  }

  // Init parameters & set defaults 
  EOBParameters_alloc( &EOBPars );
  EOBParameters_defaults (BINARY_BBH, EOBPars);

  EOBPars->q = q;
  EOBPars->chi1=chi1;
  EOBPars->chi2=chi2;

  eob_set_params(BINARY_BBH, 1);
  /* set firstcall */
  for (int k=0; k < NFIRSTCALL; k++){ 
    EOBPars->firstcall[k] = 1;
  }

  Dynamics *dyn;
  Dynamics_alloc (&dyn, 0, "dyn"); 
  Dynamics_set_params(dyn); 

  // Compute the metric
  double A[N], B[N], dA[N], d2A[N], dB[N];
  for(int i=0; i<N; i++){
    eob_metric_s(r[i], dyn, &A[i], &B[i], &dA[i], &d2A[i], &dB[i]);
  }

  // Compute the numerical derivative of the metric over the radial grid
  double dA_num[N], dB_num[N], d2A_num[N];
  D0(A, dr, N, dA_num);
  D2(A, dr, N, d2A_num);
  D0(B, dr, N, dB_num);

  // Compare point by point dA analytic with dA numeric
  double maxdA=0;
  double maxd2A=0;
  double maxdB=0;
  int err = 0;

  for (int i=0; i<N; i++){
    double check_dA  = fabs(dA[i]  - dA_num[i]);
    double check_d2A = fabs(d2A[i] - d2A_num[i]);
    double check_dB  = fabs(dB[i]  - dB_num[i]);

    // dA
    if (check_dA > maxdA && i > 0)
      maxdA = check_dA;
    if (check_dA  > 1e-5 && i > 0){
      printf("ERROR: dA may be wrong! Time to debug!\n");
      printf("|dA_analytic - dA_num| = %.10f\n", check_dA);
      err = 1;
      break;
    }
    //d2A
    if (check_d2A > maxd2A && i > 0)
      maxd2A = check_d2A;
    if (check_d2A  > 1e-5 && i > 0){
      printf("ERROR: d2A may be wrong! Time to debug!\n");
      printf("|d2A_analytic - d2A_num| = %.10f\n", check_d2A);
      err = 1;

      break;
    }
    // dB
    if (check_dB > maxdB && i > 0)
      maxdB = check_dB;
    if (check_dB  > 1e-5 && i > 0){
      printf("ERROR: dB may be wrong! Time to debug!\n");
      printf("|dB_analytic - dB_num| = %.10f\n", check_dB);
      err = 1;
      break;
    }
  }

  printf("Maximum errors for (dA, d2A, dB) = (%.8f, %.8f, %.8f)\n", maxdA, maxd2A, maxdB);

  /* Free */ 
  EOBParameters_free (EOBPars);
  Dynamics_free(dyn);

  if(err)
    errorexit("Derivatives test not successful");
}
