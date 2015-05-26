/**
 * \file types.h
 *
 * \brief Special types for EOB code.
 *
 * Defines the special types used in the EOB code. 
 * Most of them are simple structures collecting the relevant
 * quantities.    
 *
 * \author S.Bernuzzi
 */

/** 
 * \var typedef double complex domplex;
 * \brief Type definition for double complex.
 */
typedef double complex domplex;

/** 
 * \var typedef struct tPade
 * \brief Type definition for Pade coefficients.
 * \f$P^m_n(x) = \frac{n_0 + n_1 x + ... + n_{N-1} x^{N-1}}{d_0 + d_1 x + ... + d_{M-1} x^{M-1}}\f$
 */
typedef struct tPade {
  int N,M;
  double *n,*d;
  double *d1n,*d1d;
  double *d2n,*d2d;
} tPade;

/** 
 * \var typedef teobcode_PNorder
 * \brief Type definition for the formal post-Newtonian order of the model.
 */
typedef teobcode_PNorder enum{
  1PN,
  2PN,
  3PN,
  4PN,
  5PN,
  3PNP13,
  3PNP22,
  4PNP14,
  5PNP15,
  5PNP15NR
};

/** 
 * \var typedef struct tmet
 * \brief Type definition for EOB metric.
 *
 * Contains the storages for EOB metric variables. 
 */
typedef struct tmet {

  int Nt;

  double *A, *dA, *d2A;
  double *dA_u, *d2A_u;
  double *D, *dD;
  double *B, *dB;

  double *Acoefs;
  double *Dcoefs;
  tPade *PadeCoefs;

} tmet;

/** 
 * \var typedef struct tdyn
 * \brief Type definition for EOB dynamics.
 *
 * Contains the storages for EOB dynamical variables. 
 */
typedef struct tdyn {

  int Nt;
  double *t;

  double *r;
  double *phi;
  double *Prs;
  double *Pphi;

  double *Pr;  
  double *Omega;
  double *ddotr;

  double *E;
  double *Heff;
  double *jhat;

  double *NQCn;
  
} tdyn;

/** 
 * \var typedef struct tlmk
 * \brief Type definition for multipolar indexes.
 *
 * Collection of multipolar indexes.
 *
 * The index k runs over all the mulipoles, the physical multipole  
 * index can be found as l[k] or m[k]. Typically, one has
 * k    = 0 1 2 3 4 5 6 7 8 ...
 * l[k] = 2 2 2 3 3 3 3 4 4 ...
 * m[k] = 0 1 2 0 1 2 3 0 1 ...
 */
typedef struct tlmk {
  int *l, *m, *k;
  int Nk;
  int lmin,lmax,Nl;
  int k20,k21,k22,
    k30,k31,k32,k33,
    k40,k41,k42,k43,k44,
    k50,k51,k52,k53,k54,k55,
    k60,k61,k62,k63,k64,k65,k66,
    k70,k71,k72,k73,k74,k75,k76,k77,
    k80,k81,k82,k83,k84,k85,k86,k87,k88;
} tlmk;

/** 
 * \var typedef struct twav
 * \brief Type definition for multipolar waveform.
 *
 * The multipolar waveform is a complex or a real time series of
 * dimension size = Nt x Nk = (times) x (multipoles)   
 */
typedef struct twav {

  int Nt;
  int Nk;
  int size;
  tlmk *lmk;

  domplex *hlm;
  double *amp; 
  double *phi;
  double *omg;

  domplex *hlmNewt;
  double *Seff; 
  double *Tlm; 
  double *flm; 
  double *rholm; 
  double *expdlm; 
  double *deltalm; 

  double *NQCab;
 
} twav;

/** 
 * \var typedef struct teob
 * \brief Type definition for EOB variables.
 *
 * A container for storages and parameters.
 */
typedef struct teob {
  
  int Nt;
  int Nk;
  int size;

  tlmk lmk;

  tmet *met;
  tdyn *dyn;
  twav *wav;

  teobcode_PNorder PNorder;

} teob;
