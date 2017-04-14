/**
 * \file eob.h
 *
 * \brief Main header file for EOB code.
 *
 * \author S.Bernuzzi
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
// use : #include <tgmath.h> ?

#include <libconfig.h>

#include <gsl/gsl_const.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>

#define eobcode_info "EOB IHES Code"
#define eobcode_usage {printf("usage: \t %s parfile\n", argv[0]);exit(1);} 

#define eobcode_dyn_vars 4
#define eobcode_lmax 8
#define eobcode_kmax 42 /* l=2...8, m>0 */ 
#define eobcode_NQC 6 /* no NQS parameters, n_i  i = 1...6 */

#define eobcode_met_A_PNcoefs 8 /* PN coefs metric fun A, up to 5PN */
#define eobcode_met_D_PNcoefs 3 /* PN coefs metric fun D, up to 3PN */

#define eobcode_wav_rholm_PNcoefs 9 /* PN coefs rholm, up to 5PN */
#define eobcode_wav_deltalm_PNcoefs 5 /* PN coefs deltalm, up to 5PN */

void (*eob_met_A_pt)()=0;
void (*eob_met_B_pt)()=0;
void (*eob_met_D_pt)()=0;

/* include below here various headers files */
#include "utils.h"
#include "numtools.h"
#include "wavtools.h"
#include "parameters.h"
