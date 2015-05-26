/**
 * \file utils.h
 * \brief Function prototypes and macros for basic support.
 *
 * \author S.Bernuzzi
 */

/* service macros */
#define iseq(s,t) ((strcmp((s),(t))==0))

/**
 * \fn void array_sum_complex (domplex *x, domplex a, domplex *y, int size)
 * \brief Sum complex arrays, x += a * y.
 *
 * \param[in/out]  x 
 * \param[in]      a
 * \param[in]      y 
 * \return         void
 */
void array_sum_complex (domplex *x, domplex a, domplex *y, int size);

/**
 * \fn void array_sum_real (double *x, double a, double *y, int size)
 * \brief Sum real arrays, x += a * y.
 *
 * \param[in/out]  x
 * \param[in]      a
 * \param[in]      y 
 * \return         void
 */
void array_sum_real (double *x, double a, double *y, int size);

/**
 * \fn void array_mutliply_complex (domplex *x, domplex a, domplex *y, int size)
 * \brief Multiply complex array, x *= a * y
 *
 * \param[in/out]  x
 * \param[in]      a
 * \param[in]      y 
 * \return         void
 */
void array_multiply_complex (domplex *x, domplex a, domplex *y, int size);

/**
 * \fn void array_mutliply_real (double *x, double a, double *y, int size)
 * \brief Multiply real arrays, x *= a * y
 *
 * \param[in/out]  x
 * \param[in]      a
 * \param[in]      y 
 * \return         void
 */
void array_multiply_real (double *x, double a, double *y, int size);

/**
 * \fn double *interp1d_array (const int order, double *xx, int nxx, double *f, double *x, int nx)
 * \brief Interpolate f(x) on xx
 *
 * Allocate new array ff
 *
 * \param[in] order interpolation order
 * \param[in] xx    points to  interpolate
 * \param[in] nxx   size of xx
 * \param[in] f     array for f(x)
 * \param[in] x     array x
 * \param[in] nx    size of x and f
 * \return    ff    allocated and interpolated array size = nxx
 */
double *interp1d_array (const int order, double *xx, int nxx, double *f, double *x, int nx);

/**
 * \fn void errorexit(char *file, int line, char *s)
 * \brief Basic error function
 */
void errorexit(char *file, int line, char *s);
#define errorexit(s) errorexit(__FILE__, __LINE__, (s))

/**
 * \fn void errorexits(char *file, int line, char *s, char *t) 
 * \brief Basic error function with optional string
 */
void errorexits(char *file, int line, char *s, char *t);
#define errorexits(s,t) errorexits(__FILE__, __LINE__, (s), (t))

