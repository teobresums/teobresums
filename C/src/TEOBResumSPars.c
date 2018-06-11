/**
 * This file is part of TEOBResumS
 *
 * Copyright (C) 2017-2018 See AUTHORS file
 *
 * TEOBResumS is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * TEOBResumS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/.       
 *
 */

/**
 * @file TEOBResumPars.h
 * @brief Parameter manager
 *
 * Parameters are managed using a database (db) ' key = value '
 * The db is initialized by the file $TEOBRESUMS/par/default.par
 * Use libconfig API but in a simplified way (http://hyperrealm.com/libconfig/libconfig_manual.html)
 * The db entries can be set/accessed in other parts of the code with simple wrapper functions
 * Note autoconversion int/float is disabled, type must be specified in the *.par
 */

#include "TEOBResumS.h"

#define DEBUG_THIS_FILE 0 /* = 1 to compile and debug this files */
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

void par_file_parse (const char *fname)
{    
  if (!(config_read_file(cf, fname))) {
    fprintf(stderr, "%s = %d - %s\n",
	    config_error_file(cf),
	    config_error_line(cf),
	    config_error_text(cf));
    config_destroy(cf);
    errorexit("Problem reading file.");
  }
}

void par_file_parse_merge (const char *fname)
{
  /* parse the file in a new cfg */
  config_t cfg1, *cf1;
  cf1 = &cfg1;
  config_init(cf1);
  if (!(config_read_file(cf1, fname))) {
    fprintf(stderr, "%s = %d - %s\n",
	    config_error_file(cf1),
	    config_error_line(cf1),
	    config_error_text(cf1));
    config_destroy(cf1);
    errorexit("Problem reading file.");
  }
  /* set parameter values in db */
  config_setting_t * csroot1 = config_root_setting(cf1);
  config_setting_t * cs1, * array;
  int np = config_setting_length(csroot1);
  int i, type, n, j;
  char key[STRLEN];
  for (i = 0; i < np; i++) {
    cs1  = config_setting_get_elem(csroot1, i);
    strcpy(key, config_setting_name(cs1));
    type = config_setting_type(cs1);
    if      (type==CONFIG_TYPE_INT)    par_set_i(key, config_setting_get_int(cs1));
    else if (type==CONFIG_TYPE_BOOL)   par_set_b(key, config_setting_get_bool(cs1));
    else if (type==CONFIG_TYPE_FLOAT)  par_set_d(key, config_setting_get_float(cs1));
    else if (type==CONFIG_TYPE_STRING) par_set_s(key, config_setting_get_string(cs1));
    else if (type==CONFIG_TYPE_ARRAY) {
      /* this is a hack, only for int type */
      if (cs1 != NULL) {
	n = config_setting_length(cs1);
	int * array = (int *) malloc (n * sizeof(int));
	if (!array) errorexit("out of memory");
	for (j = 0; j < n; j++) {
	  array[j] = config_setting_get_int_elem(cs1, j);
	}
	par_set_arrayi (key, array, n);
	free(array);
      } 
    } else errorexit("unkown parameter type.");
  }
  /* free extra cfg */
  config_destroy(cf1);
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

void par_db_write_file (const char *name)
{
  char fname[STRLEN];
  strcpy(fname,par_get_s("output_dir"));
  strcat(fname,"/");
  strcat(fname,name);
  if (!(config_write_file(cf, fname))) {
    fprintf(stderr, "Error writing file %s\n", fname);
    config_destroy(cf);
    errorexit("Problem writing file");
  }
}

void par_db_screen (const int pr)
{
  if (pr) config_write(cf, stderr);
}

/* 'get' routines */

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

void par_set_arrayi (const char *key, int *array, int n)
{
  config_setting_t *cs = config_lookup(cf, key);
  config_setting_t *parent = config_setting_parent(cs);
  int i;
  /* remove what is there */
  if (!config_setting_remove(parent, key))
    errorexits("error removing array",key);
  /* add the new */
  config_setting_t *sarray, *a;
  sarray = config_setting_add(parent, key, CONFIG_TYPE_ARRAY);
  for(i = 0; i < n; ++i) {
    a = config_setting_add(sarray, NULL, CONFIG_TYPE_INT);
    config_setting_set_int(a, array[i]);
  }
}

void par_set_arrayd (const char *key, double *array, int n)
{
  config_setting_t *cs = config_lookup(cf, key);
  config_setting_t *parent = config_setting_parent(cs);
  int i;
  /* remove what is there */
  if (!config_setting_remove(parent, key))
    errorexits("error removing array",key);
  /* add the new */
  config_setting_t *sarray, *a;
  sarray = config_setting_add(parent, key, CONFIG_TYPE_ARRAY);
  for(i = 0; i < n; ++i) {
    a = config_setting_add(sarray, NULL, CONFIG_TYPE_FLOAT);
    config_setting_set_int(a, array[i]);
  }
}

#if (!DEBUG_THIS_FILE)

/** Set parameters */
void eob_set_params(char *s, int n)
{
  /* init db */
  par_db_init ();

  /* Set default values */
  par_db_default ();

  //if (mode==INPUT_FILE) {
  /* Parse input parfile */
  par_file_parse_merge (s);
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

  int j;

  for (j=0; j<=TIDES_NOPT; j++) {
    if (STREQUAL(par_get_s("tides"),tides_opt[j])) {
      if (DEBUG) printf("%-40s = %s\n","tides",tides_opt[j]);
      break;
    }
  }
  if (j==TIDES_NOPT) {
    if (VERBOSE) printf("tides '%s' undefined, set to default\n",par_get_s("tides"));
    j = TIDES_OFF;
  }
  par_set_i("use_tidal",j);
  const int usetidal = j;
  
  double LambdaAl2 = par_get_d("LambdaAl2");
  double LambdaBl2 = par_get_d("LambdaBl2");
  double LambdaAl3 = par_get_d("LambdaAl3");
  double LambdaBl3 = par_get_d("LambdaBl3");
  double LambdaAl4 = par_get_d("LambdaAl4");
  double LambdaBl4 = par_get_d("LambdaBl4");
  double SigmaAl2 = par_get_d("SigmaAl2");
  double SigmaBl2 = par_get_d("SigmaBl2");
  
  if (par_get_i("use_Yagi_fits")) {
    LambdaAl3 = Yagi13_fit_barlamdel(LambdaAl2, 3);
    LambdaBl3 = Yagi13_fit_barlamdel(LambdaBl2, 3);
    LambdaAl4 = Yagi13_fit_barlamdel(LambdaAl2, 4);
    LambdaBl4 = Yagi13_fit_barlamdel(LambdaBl2, 4);
    SigmaAl2 = Yagi13_fit_barsigmalambda(LambdaAl2);
    SigmaBl2 = Yagi13_fit_barsigmalambda(LambdaBl2);
  }
  
  /* Tidal coupling constants */    
  double kapA2 = 3.   * LambdaAl2 * XA*XA*XA*XA*XA / q; 
  double kapA3 = 15.  * LambdaAl3 * XA*XA*XA*XA*XA*XA*XA / q;
  double kapA4 = 105. * LambdaAl4 * XA*XA*XA*XA*XA*XA*XA*XA*XA / q;
  
  double kapB2 = 3.   * LambdaBl2 * XB*XB*XB*XB*XB * q;
  double kapB3 = 15.  * LambdaBl3 * XB*XB*XB*XB*XB*XB*XB * q;
  double kapB4 = 105. * LambdaBl4 * XB*XB*XB*XB*XB*XB*XB*XB*XB * q;
  
  /* gravitomagnetic tidal coupling constants el = 2 only */    
  double kapA2j = 24.   * SigmaAl2 * XA*XA*XA*XA*XA / q; 
  double kapB2j = 24.   * SigmaBl2 * XB*XB*XB*XB*XB * q;


  double kapT2 = kapA2 + kapB2;
  double kapT3 = kapA3 + kapB3;
  double kapT4 = kapA4 + kapB4;
  double kapT2j = kapA2j + kapB2j;
  
  if (usetidal) {
    if (!(kapT2 > 0.)) errorexit("kappaT2 must be >0");
    if (!(kapT3 > 0.)) errorexit("kappaT3 must be >0");
    if (!(kapT4 > 0.)) errorexit("kappaT4 must be >0");
    if (!(kapT2j > 0.)) errorexit("kappaT2j must be >0");
  } 
  
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
  double C_Q1 = 1.;
  double C_Q2 = 1.;
  if (LambdaAl2>0.) {
    double logC_Q1 = logQ(log(LambdaAl2));
    C_Q1           = exp(logC_Q1);
  }
  if (LambdaBl2>0.) {
    double logC_Q2 = logQ(log(LambdaBl2));
    C_Q2           = exp(logC_Q2);
  }
  
  /* Set aux parfiles in database */
  
  par_set_d("nu",nu);
  par_set_d("X1",X1);
  par_set_d("X2",X2);

  par_set_d("S1",S1);
  par_set_d("S2",S2);
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

  par_set_d("SigmaAl2",SigmaAl2);
  par_set_d("SigmaBl2",SigmaBl2);

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

  /** Set more as needed ... */
  double c3 = 0.;
  if(usetidal) c3 = 0.0;
  else         c3 = eob_c3_fit_global(nu,chi1,chi2,X1,X2,a1,a2); 
  par_set_d("cN3LO", c3 );
  
  if (par_get_i("use_geometric_units")) {
    /* input given in geometric units, 
       rescale to geometric units and mass rescaled quantities
       compute r0 from the initial GW frequency in geometric units and mass rescaled 
       reset sample rate using dt
    */
    if (DEBUG) printf("Assume geometric units for pars values\n");
    par_set_d("r0", pow(fmin*Pi, -2./3.) );
    par_set_d("srate", 1./dt );
    par_set_d("distance", 1. );
    par_set_d("M", 1. );
  } else {
    /* input given in physical units, 
       rescale to geometric units and mass rescaled quantities
       compute r0 from the initial GW frequency in Hz 
    */
    if (DEBUG) printf("Assume physical units for pars values\n");
    dt = time_units_conversion(M, dt);
    if (DUNEQUAL(dt,1./par_get_d("srate"),1e-12)) {
      dt  = 1./par_get_d("srate");
      par_set_d("dt", dt);
      printf("warning: input values for dt inconsistent with sample rate, use sample rate.\n");
    } 
    par_set_d("r0",  radius0(M, fmin) );
  }

}

void eob_free_params()
{
  par_db_free ();
}

#else

/* test 
   gcc TEOBResumSPars.c -lconfig -o testpars.x */
int main (int argc, char* argv[])
{
  int idx,n;
  /* init db */
  par_db_init ();
  /* parse default params from file */
  par_db_default ();
  /* display default params */
  par_db_screen ();
  /* test 'get' */
  printf("get-test: size = %d\n",par_get_i("size"));
  printf("get-test: M = %e\n",par_get_d("M"));
  int *i = par_get_arrayi ("output_lm", &n);
  printf("get-test: array_length n=%d\n",n);
  for (idx=0; idx<n; idx++) printf("get-test: array[%d]=%d\n",idx,i[idx]);
  /* test 'set' */
  par_set_i("size",10);
  par_set_d("M",100);
  printf("set-test: size = %d (=10)\n",par_get_i("size"));
  printf("set-test: M = %e (=100)\n",par_get_d("M"));
  n=6;
  int j[] = {1,2,3,4,5,6};
  par_set_arrayi ("output_lm", j, n);
  int *jj = par_get_arrayi ("output_lm", &n);
  printf("get-test: array_length n=%d\n",n);
  for (idx=0; idx<n; idx++) printf("get-test: array[%d]=%d\n",idx,jj[idx]);
  /* parse other params */
  par_file_parse_merge ("../par/test.par");
  /* more operations */
  printf("get-test: dt = %e\n",par_get_d("dt"));
  printf("get-test: size = %d\n",par_get_i("size"));
  par_set_i("size",20);
  printf("set-test: size = %d (=20)\n",par_get_i("size"));
  /* display new db */
  par_db_screen ();  
  /* write out db */
  par_db_write_file("../par/test1.par");
  /* free db */
  par_db_free ();
  return OK;
}

#endif
