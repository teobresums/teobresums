/** \brief parameters.c
 *         Parameters Manager
 * 
 * Wrappers for libconfig routines
 * http://www.hyperrealm.com/libconfig/
 * http://www.hyperrealm.com/libconfig/libconfig_manual.html 
 */

#include "eob.h"

config_t cfg, *cf;

void par_db_init ()
{
  cf = &cfg;
  config_init(cf);
}

void par_db_free ()
{
  config_destroy(cf);
}

void par_db_default ()
{
  static const char *eobcodeenvv = "EOBCODE";
  char *eobcodepath = NULL;
  eobcodepath = getenv (eobcodeenvv);
  if (! eobcodepath) {
    printf ("Enviroment variable '%s' is not set\n", eobcodeenvv);
    printf ("Default values for parameters are not set.\n");
    return;
  }
  else {
    par_cfgfile_parse (strcat(eobcodepath,"cfg/defaultpars.cfg"));
  }
}

void par_file_parse (char *fname)
{
  if (!(config_read_file(cf, fname))) {
    fprintf(stderr, "%s:%d - %s\n",
            config_error_file(cf),
            config_error_line(cf),
            config_error_text(cf));
    config_destroy(cf);
    errorexit("Problem reading file.");
  }
}

void par_db_write_file (char *fname)
{
  if (!(config_write_file(&cfg, fname))) {
    fprintf(stderr, "Error writing file %s\n", fname);
    config_destroy(&cfg);
    errorexit("Problem writing file");
  }
}

void par_db_screen ()
{
  config_write(&cfg, stderr);
}

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
  char *val;
  if (config_lookup_string(cf, key, &val))
    return val;
  else 
    errorexits("unknown parameter",key);
}

int * par_get_arrayi (const char *key, int *n)
{
  const config_setting_t *a;
  a = config_lookup(cf, key);
  // unsafe. todo: check
  int l = config_setting_length(a);
  int *array = NULL;
  array = (int *) malloc (l * sizeof(int));
  if (!array) errorexit("out of memory");
  int i;
  for (i = 0; i < l; l++) {
    array[i] = config_setting_get_int_elem(a, i);
  }
  *n = l;
  return array;
}

double * par_get_arrayd (const char *key, int *n)
{
  const config_setting_t *a;
  a = config_lookup(cf, key);
  // unsafe. todo: check
  int l = config_setting_length(a);
  double *array = NULL;
  array = (double *) malloc (l * sizeof(double));
  if (!array) errorexit("out of memory");
  int i;
  for (i = 0; i < l; l++) {
    array[i] = config_setting_get_float_elem(a, i);
  }
  *n = l;
  return array;
}

void par_set_i (const char *key, int *val)
{
  if (!(config_setting_lookup_int(cf, key, val)))
    errorexits("unknown parameter/wrong type for",key);
}

void par_set_b (const char *key, int *val)
{
  if (!(config_setting_lookup_bool(cf, key, val)))
    errorexits("unknown parameter/wrong type for",key);
} 

void par_set_d (const char *key, double *val)
{
  if (!(config_setting_lookup_float(cf, key, val)))
    errorexits("unknown parameter/wrong type for",key);
}

void par_set_s (const char *key, const char **val)
{
  if (!(config_setting_lookup_string(cf, key, val)))
    errorexits("unknown parameter/wrong type for",key);
}
