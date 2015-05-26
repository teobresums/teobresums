/**
 * \file utils.c
 *
 * \brief Functions and macros for basic support 
 *
 * \author S.Bernuzzi
 */

#include "eob.h"

/* basic array operations */
void array_sum_complex (domplex *x, domplex a, domplex *y, int size)
{
  int i;
  for (i = 0; i < size; i++) {
    x[i] += a * y[i];
  }
}

void array_sum_real (double *x, double a, double *y, int size)
{
  
  int i;
  for (i = 0; i < size; i++) {
    x[i] += a * y[i];
  }
}

void array_multiply_complex (domplex *x, domplex a, domplex *y, int size)
{
  int i;
  if (a==1.) {
    for (i = 0; i < size; i++) {
      x[i] *= y[i];
    }
  } else {
    for (i = 0; i < size; i++) {
      x[i] *= a * y[i];
    }
  }
}

void array_multiply_real (double *x, double a, double *y, int size)
{
  int i;
  if (a==1.) {
    for (i = 0; i < size; i++) {
      x[i] *= y[i];
    }
  } else {
   for (i = 0; i < size; i++) {
     x[i] *= a * y[i];
   }
  }
}

/* 1d Lagrangian barycentric interpolation on array */
double *interp1d_array (const int order, double *xx, int nxx, double *f, double *x, int nx)
{
  double **ff;
  ff = calloc( nxx, sizeof(double) );
  int i;
  for (i=0; i<nx; i++)
    ff[i] = interp1d (order, xx[i], nx, f,x);
  return ff;
}

/* errorexit routines */
#undef errorexit
#undef errorexits

void errorexit(char *file, int line, char *s) 
{
  fprintf(stderr, "Error: %s  ", s);
  fprintf(stderr, "(%s, line %d)\n", file, line);
  fflush(stdout);
  fflush(stderr);
  exit(1);
}

void errorexits(char *file, int line, char *s, char *t)
{
  fprintf(stderr, "Error: ");
  fprintf(stderr, s, t);
  fprintf(stderr, "  (%s, line %d)\n", file, line);
  fflush(stdout);
  fflush(stderr);
  exit(1);
}

/* do not write functions beyond this line: errorexit undef/define */
