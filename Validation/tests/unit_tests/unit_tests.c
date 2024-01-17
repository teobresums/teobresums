/** \file tests.c
  *  \brief Unit tests for TEOBResumS
  * Testing is performed using the mintunit library
  * see https://jera.com/techinfo/jtns/jtn002
  *     https://github.com/danielgustafsson/mint
*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "mintunit.h"
#include "../../../C/src/TEOBResumS.h"

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
Dynamics *dyn;

/**
* List tests below
*/
int test_metric_derivatives(double q, double chi1, double chi2);
int test_ecc_ICs(double q, double chi1, double chi2, double ecc, double f0);

static char *
test_metric(void)
{ 
  double q[3]    = {1., 1.,  2.};
  double chi1[3] = {0., 0.1, 0.9};
  double chi2[3] = {0., -0.3, 0.};
  
  for(int i=0; i<3;i++)
	  mt_assert("Test metric function", test_metric_derivatives(q[i],chi1[i],chi2[i]) == 0
             );
	return NULL;
}

static char *
test_ICs(void)
{ 
  double q[4]    = {1.     ,   1.,  2., 10.};
  double chi1[4] = {0.01   ,  0.1, 0.9, 0.};
  double chi2[4] = {0.001  , -0.3,  0., -0.5};
  double ecc[4]  = {0.00001,  0.1, 0.9, 0.001};
  double f0      = 0.001;

  for(int i=0; i<4;i++)
	  mt_assert("Test eccentric ICs", test_ecc_ICs(q[i],chi1[i],chi2[i], ecc[i], f0) == 0
             );
	return NULL;
}

static char *
run_tests()
{
	mt_init(2);
	mt_run_test("Testing metric", test_metric);
	mt_run_test("Testing ICs", test_ICs);

	return NULL;
}

/**
 * Function: test_ecc_ICs
 * ----------------------
 * 
 *   Routine to test the initial conditions for the
 *   eccentric code. Compare the "generic" mean-anomaly function 
 *   with anomaly = pi to the 0PA function.
 * 
 *   @param[in] q : mass ratio
 *   @param[in] chi1: dimensionless spin of body 1
 *   @param[in] chi2: dimensionless spin of body 2
 *   @param[in] ecc: eccentricity
 *   @param[in] f0: initial frequency
 * 
 *   @return 0 or 1 depending on whether test passed
*/
int test_ecc_ICs(double q, double chi1, double chi2, double ecc, double f0)
{ 
  //reset relevant parameters
  Dynamics *dyn_0PA;
  Dynamics *dyn_ma;
  Dynamics_alloc (&dyn_0PA, 0, "dyn"); 
  Dynamics_alloc (&dyn_ma, 0, "dyn");
  EOBParameters_defaults(BINARY_BBH, 1, EOBPars);
  EOBPars->q    = q;
  EOBPars->chi1 = chi1;
  EOBPars->chi2 = chi2;
  EOBPars->ecc  = ecc;
  EOBPars->f0   = f0;
  eob_set_params(BINARY_BBH, 1);
  Dynamics_set_params(dyn_0PA);
  Dynamics_set_params(dyn_ma);

  double r0_ecc = eob_dyn_r0_ecc(f0, dyn_0PA);

  eob_dyn_ic_ecc(r0_ecc, dyn_0PA, dyn_0PA->y0);
  eob_dyn_ic_ecc_ma(r0_ecc, dyn_ma, dyn_ma->y0);

  for(int v=0; v< EOB_DYNAMICS_NVARS; v++){
    if (fabs(dyn_0PA->y0[v] - dyn_ma->y0[v]) > 1e-4){
        printf("Error in ICs for %s %lf %lf\n", eob_id_var[v], dyn_0PA->y0[v], dyn_ma->y0[v]);
      	return 1;
    }
  }
  Dynamics_free(dyn_ma);
  Dynamics_free(dyn_0PA);
  return 0;
}

/**
  * Function: test_mtric_derivatives
  * --------------------------------
  *   Routine to test that the A(nu; r) function 
  *   and derivatives are implemented correctly via
  *   comparisons to numerical derivatives
  *   
  *   @param[in] q : mass ratio
  *   @param[in] chi1: dimensionless spin of body 1
  *   @param[in] chi2: dimensionless spin of body 2
  *
  *   @return 0 or 1 depending on whether test passed
  */
int test_metric_derivatives(double q, double chi1, double chi2)
{ 
  //reset relevant parameters
  // TODO: add tidal part

  EOBPars->q = q;
  EOBPars->chi1 = chi1;
  EOBPars->chi2 = chi2;
  eob_set_params(BINARY_BBH, 1);
  //radial grid
  const double rmin = 2.;
  const double rmax = 20.;
  const double dr   = 0.1;
  const int N       = floor((rmax-rmin)/dr) +1;
  double pl_hold; // placeholder for unused variables
  double r[N];
  for (int i=0; i<N; i++)
    r[i] = rmin + i*dr;

  // Compute the metric
  double A[N], B[N], dA[N], d2A[N], dB[N];
  for(int i=0; i<N; i++){
    eob_metric_s(r[i], 0., dyn, &A[i], &B[i], &dA[i], &d2A[i], &dB[i], 
                           &pl_hold, &pl_hold, &pl_hold, &pl_hold, 
                           &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
  }

  // Compute the numerical derivative of the metric over the radial grid
  double dA_num[N], dB_num[N], d2A_num[N];
  D0(A, dr, N, dA_num);
  D2(A, dr, N, d2A_num);
  D0(B, dr, N, dB_num);

  // Compare point by point dA analytic with dA numeric
  double maxdA  = 0;
  double maxd2A = 0;
  double maxdB  = 0;
  int err       = 0;

  for (int i=0; i<N; i++){
    double check_dA  = fabs(dA[i]  - dA_num[i]);
    double check_d2A = fabs(d2A[i] - d2A_num[i]);
    double check_dB  = fabs(dB[i]  - dB_num[i]);
    // dA
    if (check_dA > maxdA && i > 0)
      maxdA = check_dA;
    if (check_dA  > 1e-5 && i > 0){
      err = 1;
      break;
    }
    //d2A
    if (check_d2A > maxd2A && i > 0)
      maxd2A = check_d2A;
    if (check_d2A  > 1e-5 && i > 0){
      err = 1;
      break;
    }
    // dB
    if (check_dB > maxdB && i > 0)
      maxdB = check_dB;
    if (check_dB  > 1e-5 && i > 0){
      err = 1;
      break;
    }
  }
  return err;
}

int main(int argc, char **argv)
{

  // Preliminaries
  EOBParameters_alloc( &EOBPars );
  EOBParameters_defaults (BINARY_BBH, 0, EOBPars);
  eob_set_params(BINARY_BBH, 1);
  Dynamics_alloc (&dyn, 0, "dyn"); 
  Dynamics_set_params(dyn);
  /* set firstcall */
  for (int k=0; k < NFIRSTCALL; k++){ 
    EOBPars->firstcall[k] = 1;
  }

  /* Run tests */
  char *test_result = run_tests();
	if (!test_result)
		printf("\n# All tests passed\n");
	printf("# Test Suites run: %d\n# Unit Tests run : %d\n", tests_run, units_run);

  /* Free */ 
  EOBParameters_free (EOBPars);
  Dynamics_free(dyn);

	return test_result == NULL ? 0 : 1;
}
