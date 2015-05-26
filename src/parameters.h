
/* parameters.c */
void par_db_init ();
void par_db_free ();
void par_db_default ();
void par_file_parse (char *fname);
void par_db_write_file (char *fname);
void par_db_screen ();
void par_set_i(const char *key, int *val);
void par_set_b(const char *key, int *val);
void par_set_d(const char *key, double *val);
void par_set_s(const char *key, const char **val);
int par_get_i(const char *key);
int par_get_b(const char *key);
double par_get_d(const char *key);
const char * par_get_s(const char *key);
int * par_get_arrayi(const char *key, int *n);
double * par_get_arrayd(const char *key, int *n);
