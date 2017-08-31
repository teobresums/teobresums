/**
 * \file numtools.h
 * \brief Function prototypes and macros for numerical routines.
 *
 * \author S.Bernuzzi
 */

/* constants */
#define Pi 3.1415926535897932384626433832795028
#define Sq2Pi 2.5066282746310005024157652848110 /* sqrt(2Pi) */
#define Sqrt2 1.41421356237309504880168872420969808 /* sqrt(2) */
#define Sqrt3 1.73205080756887729352744634150587237 /* sqrt(3) */
#define ooSqrt2 0.707106781186547524400844362104849039284836 /* 1/sqrt(2) */
#define oo2 0.500000000000000000000000000000000000 /* 1/2 */
#define oo3 0.333333333333333333333333333333333333 /* 1/3 */
#define oo12 0.0833333333333333333333333333333333333 /* 1/12 */
#define to3 0.6666666666666666666666666666666 /* 2/3 */ 
#define EulerGamma 0.57721566490153286061
#define twooSqrte 1.213061319425267e+00 /* 2/sqrt(e); */
#define log2 0.693147180559945309417232

/* math macros */
#define SIGN(x,y) ((y) >= 0.0 ? fabs(x) : -fabs(x))
#define MAX(a,b)                \
  ({ typeof (a) _a = (a);       \
    typeof (b) _b = (b);        \
    _a > _b ? _a : _b; })	
#define MIN(a,b)		\
  ({ typeof (a) _a = (a);       \
    typeof (b) _b = (b);        \
    _a < _b ? _a : _b; })	
#define MAX3(a,b,c) (((a) > (b)) ? MAX(a,c) : MAX(b,c))
#define MIN3(a,b,c) (((a) < (b)) ? MIN(a,c) : MIN(b,c))
#define SQ(a) ((a)*(a))
#define DEQUAL(a,b,eps) (fabs((a)-(b))<(eps))
#define DUNEQUAL(a,b,eps) (fabs((a)-(b))>(eps))

/* debug macros */
#define dbg_pr(n) printf("dbg pt %d\n",n);
#define dbg_stop errorexit("dbg stop");

/**
 * \fn double Eulerlog (double x, int m)
 * \brief Euler log function.
 *
 * \param[in] x 
 * \param[in] m 
 * \return    \f$ 0.57721... + \log(2 m \sqrt(x))\f$
 */
double Eulerlog (double x, int m);
 
/**
 * \fn domplex GammaComplex (domplex z)
 * \brief Complex Gamma function valid on the entire complex plane.
 *
 * Complex Gamma function valid in the entire complex plane. 
 * Accuracy is 15 significant digits along the real axis and 
 * 13 significant digits elsewhere. This routine uses a superb 
 * Lanczos series approximation for the complex Gamma function. 
 *
 * \param[in] z complex double
 * \return    complex double \f$\Gamma(z)\f$
 *
 * \author Paul Godfrey
 *         pgodfrey@intersil.com
 *         http://winnie.fit.edu/~gabdo/gamma.txt
 *         Sept 11, 2001
 */
domplex GammaComplex (domplex z);

/**
 * \fn double fact (double n)
 * \brief Factorial.
 *
 * \param[in] n double 
 * \return    \f$n!\f$ factorial
 */
double fact (double n);

/**
 * \fn int find_point_bisection(double x, int n, double *xp, int o)
 * \brief Find nearest point index in 1d array. 
 *
 * Find index such that xp[i] <= x < xp[i+1] uses bisection, which
 * relies on x being ordered o is "offset", number of points smaller
 * than x that are required. Returns j = i-(o-1), i.e. 
 * if o = 2, then xp[j] < xp[j+1] <= x < xp[j+2] < xp[j+3]
 * which is useful for interpolation.
 *
 * \param[in] x  point to find 
 * \param[in] n  number of points in array
 * \param[in] xp array
 * \param[in] o  offset
 * \return    i  index 
 */
int find_point_bisection (double x, int n, double *xp, int o);

/**
 * \fn double interp1d (const int order, double xx, int nx, double *f, double *x)
 * \brief 1d Lagrangian barycentric interpolation.
 *
 * \param[in] order polynomial order
 * \param[in] xx    point to interpolate 
 * \param[in] nx    arary size
 * \param[in] f     array
 * \param[in] x     array
 * \return    ff    \f$f(x)\f$ interpolated at point xx 
 */
double interp1d (const int order, double xx, int nx, double *f, double *x);
