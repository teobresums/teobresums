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

/**
 * @file TEOBResumPars.h
 * @brief Parameter manager
 *
 * Parameters are managed using a database (db) ' key = value '
 * The db is initialized by the file $TEOBRESUMS/par/default.par
 * Use libconfig API but in a simplified way (http://hyperrealm.com/libconfig/libconfig_manual.html)
 * The db entries can be set/accessed in other parts of the code with simple wrapper functions
 */

#include "TEOBResumS.h"

#define DEBUG_THIS_FILE 1 /* to compile and debug this files */
#if (DEBUG_THIS_FILE)
#undef errorexit
#undef errorexits
#define errorexit(s) {printf("%s\n",s); exit(ERROR);}
#define errorexits(s,t) {printf("%s %s\n",s,t); exit(ERROR);}
#endif

/** database for parameters */
config_t cfg, *cf;
config_setting_t * csroot;

void par_db_init ()
{
  cf = &cfg;
  config_init(cf);
  /* csroot = config_root_setting(cf); */
}

void par_db_free ()
{
  config_destroy(cf);
}

void par_file_parse (char *fname)
{
  if (!(config_read_file(cf, fname))) {
    fprintf(stderr, "%s = %d ; %s\n",
            config_error_file(cf),
            config_error_line(cf),
            config_error_text(cf));
    config_destroy(cf);
    errorexit("Problem reading file.");
  }
}

/** default values for parameters are expected in file
    $TEOBRESUMS/par/default.par */
void par_db_default ()
{
  static const char *eobcodeenvv = "TEOBRESUMS";
  char *eobcodepath = NULL;
  eobcodepath = getenv (eobcodeenvv);
  if (! eobcodepath) {
    printf ("Enviroment variable '%s' is not set\n", eobcodeenvv);
    printf ("Default values for parameters are not set.\n");
    return;
  }
  else {
    par_file_parse (strcat(eobcodepath,"/par/default.par"));
  }
}

void par_db_write_file (char *fname)
{
  if (!(config_write_file(cf, fname))) {
    fprintf(stderr, "Error writing file %s\n", fname);
    config_destroy(cf);
    errorexit("Problem writing file");
  }
}

void par_db_screen ()
{
  config_write(cf, stderr);
}

/* 'get' routines, assume all settings are root */

int par_get_i (const char *key)
{
  int val;
  if (config_lookup_int(cf, key, &val))
    return val;
  else 
    errorexits("unknown parameter", key);
}

int par_get_b (const char *key)
{
  int val;
  if (config_lookup_bool(cf, key, &val))
    return val;
  else 
    errorexits("unknown parameter", key);
} 

double par_get_d (const char *key)
{
  double val;
  if (config_lookup_float(cf, key, &val))
    return val;
  else 
    errorexits("unknown parameter",key);
}

const char * par_get_s (const char *key)
{
  const char *val;
  if (config_lookup_string(cf, key, &val))
    return val;
  else 
    errorexits("unknown parameter",key);
}

int * par_get_arrayi (const char *key, int *n)
{
  const config_setting_t *cs = config_lookup(cf, key);
  if (cs != NULL) {
    int l = config_setting_length(cs);
    int * array = (int *) malloc (l * sizeof(int));
    if (!array) errorexit("out of memory");
    int i;
    for (i = 0; i < l; i++) {
      array[i] = config_setting_get_int_elem(cs, i);
    }
    *n = l;
    return array;
  } else
    errorexits("unknown parameter",key);
}

double * par_get_arrayd (const char *key, int *n)
{
  const config_setting_t *cs = config_lookup(cf, key);
  if (cs != NULL) {
    int l = config_setting_length(cs);
    double * array = (double *) malloc (l * sizeof(int));
    if (!array) errorexit("out of memory");
    int i;
    for (i = 0; i < l; i++) {
      array[i] = config_setting_get_int_elem(cs, i);
    }
    *n = l;
    return array;
  } else
    errorexits("unknown parameter",key);
}

/* 'set' routines, assume all settings are root */

void par_set_i (const char *key, int val)
{
  config_setting_t * cs = config_lookup(cf, key);
  if (!(config_setting_set_int(cs, val)))
    errorexits("unknown parameter/wrong type for",key);
}

void par_set_b (const char *key, int val)
{
  config_setting_t * cs = config_lookup(cf, key);
  if (!(config_setting_set_bool(cs, val)))
    errorexits("unknown parameter/wrong type for",key);
} 

void par_set_d (const char *key, double val)
{
  config_setting_t * cs = config_lookup(cf, key);
  if (!(config_setting_set_float(cs, val)))
    errorexits("unknown parameter/wrong type for",key);
}

void par_set_s (const char *key, const char *val)
{
  config_setting_t * cs = config_lookup(cf, key);
  if (!(config_setting_set_string(cs, val)))
    errorexits("unknown parameter/wrong type for",key);
}

#if (!DEBUG_THIS_FILE)

/** Set parameters */
//void TEOBResumSSetParameters(char *s, int n, int mode, int pr)
void TEOBResumSSetParameters(char *s, int pr)
{

  /* Set default values */
  par_db_default ();

  //if (mode==INPUT_FILE) {
  /* Parse input parfile */
  par_file_parse (s);
  //} 
  //if (mode==COMMAND_LINE) {
  /* Parse command line */
  //par_file_parse (s, n);    
  //}
  //if (mode==NONE) {
  //errorexit("unknown mode for input parameters mode ");
  //}
  
  /* Set auxiliary parameters */
  double dt = par_get_d("dt");
  double M = par_get_d("M");
  double fmin = par_get_d("initial_frequency");
  
  double q =  par_get_d("q");  
  double nu = q_to_nu(q);
  double X1 = nu_to_X1(nu);
  double X2 = 1. - X1;
  double XA = X1; /* tidal part used different notation, used here for simplicity */
  double XB = X2;
  
  double chi1 = par_get_d("chi1");
  double chi2 = par_get_d("chi2");
  double S1 = SQ(X1) * chi1;
  double S2 = SQ(X2) * chi2;
  double a1  = X1*chi1;
  double a2  = X2*chi2;
  double aK  = a1 + a2;
  double aK2 = aK*aK;   
  double S = S1 + S2;            /* in the EMRL this becomes the spin of the BH */
  double Sstar = X2*a1 + X1*a2;  /* in the EMRL this becomes the spin of the particle */

  // lets not do the following. flag_spin should rule!
  // Override spin settings if spins are given in input 
  //if (chi1 != .0 || chi2 != .0) params.flags.spin = 1;

  int flag_tidal = par_get_i("use_tidal");
  
  double LambdaAl2 = par_get_d("LambdaAl2");
  double LambdaBl2 = par_get_d("LambdaBl2");
  double LambdaAl3 = par_get_d("LambdaAl3");
  double LambdaBl3 = par_get_d("LambdaBl3");
  double LambdaAl4 = par_get_d("LambdaAl4");
  double LambdaBl4 = par_get_d("LambdaBl4");
  
  if (par_get_i("use_Yagi_fits")) {
    LambdaAl3 = Yagi13_fit_barlamdel(LambdaAl2, 3);
    LambdaBl3 = Yagi13_fit_barlamdel(LambdaBl2, 3);
    LambdaAl4 = Yagi13_fit_barlamdel(LambdaAl2, 4);
    LambdaBl4 = Yagi13_fit_barlamdel(LambdaBl2, 4);
  }
  
  /* Tidal coupling constants */    
  double kapA2 = 3.   * LambdaAl2 * XA*XA*XA*XA*XA / q; 
  double kapA3 = 15.  * LambdaAl3 * XA*XA*XA*XA*XA*XA*XA / q;
  double kapA4 = 105. * LambdaAl4 * XA*XA*XA*XA*XA*XA*XA*XA*XA / q;
  
  double kapB2 = 3.   * LambdaBl2 * XB*XB*XB*XB*XB * q;
  double kapB3 = 15.  * LambdaBl3 * XB*XB*XB*XB*XB*XB*XB * q;
  double kapB4 = 105. * LambdaBl4 * XB*XB*XB*XB*XB*XB*XB*XB*XB * q;
  
  double kapT2 = kapA2 + kapB2;
  double kapT3 = kapA3 + kapB3;
  double kapT4 = kapA4 + kapB4;

  /* Tidal coefficients cons dynamics
     \bar{\alpha}_n^{(\ell)}, Eq.(37) of Damour&Nagar, PRD 81, 084016 (2010) */
  double bar_alph2_1 = (5./2.*XA*kapA2 + 5./2.*XB*kapB2)/kapT2;
  double bar_alph2_2 = ((3.+XA/8.+ 337./28.*XA*XA)*kapA2 + (3.+XB/8.+ 337./28.*XB*XB)*kapB2)/kapT2;
  double bar_alph3_1 = ((-2.+15./2.*XA)*kapA3 + (-2.+15./2.*XB)*kapB3)/kapT3;
  double bar_alph3_2 = ((8./3.-311./24.*XA+110./3.*XA*XA)*kapA3 + (8./3.-311./24.*XB+110./3.*XB*XB)*kapB3)/kapT3;

  /* Tidal coefficients for the amplitude */
  double khatA_2  = 3./2. * LambdaAl2 * XB/XA * gsl_pow_int(XA,5);
  double khatB_2  = 3./2. * LambdaBl2 * XA/XB * gsl_pow_int(XB,5);

  /* self-spin coefficients */
  double logC_Q1 = logQ(log(LambdaAl2));
  double logC_Q2 = logQ(log(LambdaBl2));
  double C_Q1    = exp(logC_Q1);
  double C_Q2    = exp(logC_Q2);
  
  /* Set aux parfiles in database */
  
  par_set_d("nu",nu);
  par_set_d("X1",X1);
  par_set_d("X2",X2);

  par_set_d("S1",S1);
  par_set_d("S1",S2);
  par_set_d("S",S);
  par_set_d("Sstar",Sstar);
  par_set_d("a1",a1);
  par_set_d("a2",a2);
  par_set_d("aK",aK);
  par_set_d("aK2",aK2);

  par_set_d("LambdaAl3",LambdaAl3);
  par_set_d("LambdaAl4",LambdaAl4);
  par_set_d("LambdaBl3",LambdaBl3);
  par_set_d("LambdaBl4",LambdaBl4);

  par_set_d("kappaAl2", kapA2 );
  par_set_d("kappaAl3", kapA3 );
  par_set_d("kappaAl4", kapA4 );
  par_set_d("kappaBl2", kapB2 );
  par_set_d("kappaBl3", kapB3 );
  par_set_d("kappaBl4", kapB4 );
  par_set_d("kappaTl2", kapT2 );
  par_set_d("kappaTl3", kapT3 );
  par_set_d("kappaTl4", kapT4 );

  par_set_d("bar_alph2_1",bar_alph2_1);
  par_set_d("bar_alph2_2",bar_alph2_2);
  par_set_d("bar_alph3_1",bar_alph3_1);
  par_set_d("bar_alph3_2",bar_alph3_2);

  par_set_d("khatAl2",khatA_2);
  par_set_d("khatBl2",khatB_2);

  par_set_d("C_Q1",C_Q1);
  par_set_d("C_Q2",C_Q2);

  // Set more as needed ...
  double c3  = c3_fit_global(nu,chi1,chi2,X1,X2,a1,a2);
  par_set_d("cN3LO", c3 );
  
  if (par_get_i("use_geometric_units")) {
    /* input given in geometric units, 
       rescale to geometric units and mass rescaled quantities
       compute r0 from the initial GW frequency in geometric units and mass rescaled 
       reset sample rate using dt
    */
    par_set_d("r0", pow(fmin*Pi, -2./3.) );
    par_set_d("srate", 1./dt );
    par_set_d("distance", 1. );
    par_set_d("M", 1. );
  } else {
    /* input given in physical units, 
       rescale to geometric units and mass rescaled quantities
       compute r0 from the initial GW frequency in Hz 
    */
    dt = time_units_conversion(M, dt);
    if (DUNEQUAL(dt,1./par_get_d("srate"),1e-12)) {
      dt  = 1./par_get_d("srate");
      par_set_d("dt", dt);
      printf("warning: input values for dt inconsistent with sample rate, use sample rate.\n");
    } 
    par_set_d("r0",  radius0(M, fmin) );
  }

  /* Print the parameters used for the run */
  if (pr) par_db_screen ();

}

#else

/* test */
int main (int argc, char* argv[])
{
  /* init db */
  par_db_init ();
  /* parse default params from file */
  par_db_default ();
  /* display default params */
  par_db_screen ();
  /* test 'get' */
  printf("get-test: size = %d\n",par_get_i("size"));
  printf("get-test: M = %e\n",par_get_d("M"));
  int idx,n;
  int *i = par_get_arrayi ("output_lm", &n);
  printf("get-test: array_length n=%d\n",n);
  for (idx=0; idx<n; idx++) printf("get-test: array[%d]=%d\n",idx,i[idx]);
  /* test 'set' */
  par_set_i("size",10);
  par_set_d("M",100);
  printf("set-test: size = %d (=10)\n",par_get_i("size"));
  printf("set-test: M = %e (=100)\n",par_get_d("M"));
  /* parse other params */
  par_file_parse ("../par/test.par");
  /* more operations */
  printf("get-test: size = %d\n",par_get_i("size"));
  par_set_i("size",20);
  printf("set-test: size = %d (=20)\n",par_get_i("size"));
  /* display new db */
  //par_db_screen ();  
  /* write out db */
  //par_db_write_file("../par/test1.par");
  /* free db */
  par_db_free ();
  return OK;
}

#endif
