/**
 * \file wavtools.h
 * \brief Function prototypes for multipolar waveform manipulations.
 *
 * Contains the basic function prototypes to manipulate the multipolar
 * waveforms types, including routines for multipolar indexes.     
 *
 * \author S.Bernuzzi
 */

/**
 * \fn void wav_alloc (twav **wav, int Nt, int lmin, int lmax)
 * \brief Allocate memory for a twav type
 *
 * \param[in] wav multipolar waveform to allocate 
 * \param[in] Nt
 * \param[in] lmin
 * \param[in] lmax
 * \return    void
 *
 * \note Does not allocate data storages !
 */
void wav_alloc (twav **wav, int Nt, int lmin, int lmax);

/**
 * \fn void wav_free (twav *wav)
 * \brief Free memory of a twav type
 *
 * \param[in] wav multipolar waveform to free
 * \return    void
 */
void wav_free (twav *wav);

/**
 * \fn void wav_amplitude (domplex *h, const int size, double *a);
 * \brief Compute amplitude of complex multipolar waveform.
 *
 * \param[in]  h     complex waveform 
 * \param[in]  size  size of h 
 * \param[out] a     real amplitude |h| 
 * \return     void
 */
void wav_amplitude (domplex *h, const int size, double *a);

/**
 * \fn void wav_phase (domplex *h, const int size, double *p);
 * \brief Compute phase of complex multipolar waveform.
 *
 * \param[in]  h     complex waveform 
 * \param[in]  size  size of h 
 * \param[out] p     real phase arg(h) 
 * \return     void 
 */
void wav_phase (domplex *h, const int size, double *p);

/**
 * \fn void wav_sumlm_complex (domplex *h, int Nt, int Nk, int *l, domplex *sum)
 * \brief Sum of complex multipolar waveform
 *        \f$\sum_{\ell=\ell_{min}}^{\ell=\ell_{max}}\sum_{m=0}^\ell
 *        h_{\ell m}\f$.  
 *
 * \param[in]  h   multipolar complex waveform
 * \param[in]  Nt  time size
 * \param[in]  Nk  multipolar size
 * \param[in]  l   multipolar index \f$\ell\f$
 * \param[out] sum complex sums at each time
 * \return     void 
 */
void wav_sumlm_complex (domplex *h, int Nt, int Nk, int *l, domplex *sum);

/**
 * \fn void wav_sumlm_real (double *hlm, int Nt, int Nk, int *l, double *sum)
 * \brief Sum of real multipolar waveform
 *        \f$\sum_{\ell=\ell_{min}}^{\ell=\ell_{max}}\sum_{m=0}^\ell
 *        h_{\ell m}\f$ 
 *
 * \param[in]  h   multipolar real waveform
 * \param[in]  Nt  time size
 * \param[in]  Nk  multipolar size
 * \param[in]  l   multipolar index \f$\ell\f$
 * \param[out] sum real sums at each time
 * \return     void 
 */
void wav_sumlm_real (double *hlm, int Nt, int Nk, int *l, double *sum);

/**
 * \fn
 * \brief 
 *        
 *
 * \return    
 *
 */
tlmk *lmk_alloc ();

/**
 * \fn
 * \brief 
 *        
 *
 * \return    
 *
 */
void lmk_free (tlmk *lmk);
 
/**
 * \fn int lmindex_lmk (int lmin, int lmax, int *l, int *m, int *k) 
 * \brief Given \f$(\ell_{min},\ell_{max})\f$, return index
 *        lists \f$\ell,m,k\f$.  
 *
 * \param[in]  lmin,lmax \f$(\ell_{min},\ell_{max})\f$
 * \param[out] l,m,k     \f$\ell,m,k\f$ int arrays
 * \return     int       size of \f$\ell,m,k\f$ arrays
 *
 * \warning Specific for \f$m>0\f$
 */
int lmk_set_lmk (int lmin, int lmax, int *l, int *m, int *k); 

/**
 * \fn int lmindex_lm2k (int l, int m) 
 * \brief Given \f$(\ell,m)\f$, return single multipolar index
 *        \f$k\f$.  
 *
 * \param[in]  l,m  \f$(\ell,m)\f$
 * \return     k    \f$k\f$ 
 *
 * \warning Hardcoded for \f$\ell=2,...,8\f$
 * \warning Specific for \f$m>0\f$
 */
int lmk_lm2k (int l, int m); 

/**
 * \fn void lmindex_k2lm (int k, int *l, int *m)
 * \brief Given \f$k\f$, return multipolar indexes \f$(\ell,m)\f$.  
 *
 * \param[in]   k   \f$k\f$
 * \param[out]  l,m \f$(\ell,m)\f$ values
 * \return      void
 *
 * \warning Hardcoded for \f$\ell=2,...,8\f$
 * \warning Specific for \f$m>0\f$ // TODO !!!!
 */
void lmk_k2lm (int k, int *l, int *m); 

/* macros for multipolar indexes */
#define eobcode_kmultipolar_index_giveall	\
  int k21 = lmk_lm2k(2,1);			\
  int k22 = lmk_lm2k(2,2);			\
  int k31 = lmk_lm2k(3,1);			\
  int k32 = lmk_lm2k(3,2);			\
  int k33 = lmk_lm2k(3,3);			\
  int k41 = lmk_lm2k(4,1);			\
  int k43 = lmk_lm2k(4,2);			\
  int k43 = lmk_lm2k(4,3);			\
  int k44 = lmk_lm2k(4,4);			\
  int k51 = lmk_lm2k(5,1);			\
  int k52 = lmk_lm2k(5,2);			\
  int k53 = lmk_lm2k(5,3);			\
  int k54 = lmk_lm2k(5,4);			\
  int k55 = lmk_lm2k(5,5);			\
  int k61 = lmk_lm2k(6,1);			\
  int k62 = lmk_lm2k(6,2);			\
  int k63 = lmk_lm2k(6,3);			\
  int k64 = lmk_lm2k(6,4);			\
  int k65 = lmk_lm2k(6,5);			\
  int k66 = lmk_lm2k(6,6);			\
  int k71 = lmk_lm2k(7,1);			\
  int k72 = lmk_lm2k(7,2);			\
  int k73 = lmk_lm2k(7,3);			\
  int k74 = lmk_lm2k(7,4);			\
  int k75 = lmk_lm2k(7,5);			\
  int k76 = lmk_lm2k(7,6);			\
  int k77 = lmk_lm2k(7,7);			\
  int k81 = lmk_lm2k(8,1);			\
  int k82 = lmk_lm2k(8,2);			\
  int k83 = lmk_lm2k(8,3);			\
  int k84 = lmk_lm2k(8,4);			\
  int k85 = lmk_lm2k(8,5);			\
  int k86 = lmk_lm2k(8,6);			\
  int k87 = lmk_lm2k(8,7);			\
  int k88 = lmk_lm2k(8,8);

#define eobcode_kmultipolar_index_lmkgive	\
  int k21 = lmk->k21;				\
  int k22 = lmk->k22;				\
  int k31 = lmk->k31;				\
  int k32 = lmk->k32;				\
  int k33 = lmk->k33;				\
  int k41 = lmk->k41;				\
  int k43 = lmk->k42;				\
  int k43 = lmk->k43;				\
  int k44 = lmk->k44;				\
  int k51 = lmk->k51;				\
  int k52 = lmk->k52;				\
  int k53 = lmk->k53;				\
  int k54 = lmk->k54;				\
  int k55 = lmk->k55;				\
  int k61 = lmk->k61;				\
  int k62 = lmk->k62;				\
  int k63 = lmk->k63;				\
  int k64 = lmk->k64;				\
  int k65 = lmk->k65;				\
  int k66 = lmk->k66;				\
  int k71 = lmk->k71;				\
  int k72 = lmk->k72;				\
  int k73 = lmk->k73;				\
  int k74 = lmk->k74;				\
  int k75 = lmk->k75;				\
  int k76 = lmk->k76;				\
  int k77 = lmk->k77;				\
  int k81 = lmk->k81;				\
  int k82 = lmk->k82;				\
  int k83 = lmk->k83;				\
  int k84 = lmk->k84;				\
  int k85 = lmk->k85;				\
  int k86 = lmk->k86;				\
  int k87 = lmk->k87;				\
  int k88 = lmk->k88;

#define eobcode_kmultipolar_index_lmksetk	\
  lmk->k21 = lmk_lm2k(2,1);			\
  lmk->k22 = lmk_lm2k(2,2);			\
  lmk->k31 = lmk_lm2k(3,1);			\
  lmk->k32 = lmk_lm2k(3,2);			\
  lmk->k33 = lmk_lm2k(3,3);			\
  lmk->k41 = lmk_lm2k(4,1);			\
  lmk->k43 = lmk_lm2k(4,2);			\
  lmk->k43 = lmk_lm2k(4,3);			\
  lmk->k44 = lmk_lm2k(4,4);			\
  lmk->k51 = lmk_lm2k(5,1);			\
  lmk->k52 = lmk_lm2k(5,2);			\
  lmk->k53 = lmk_lm2k(5,3);			\
  lmk->k54 = lmk_lm2k(5,4);			\
  lmk->k55 = lmk_lm2k(5,5);			\
  lmk->k61 = lmk_lm2k(6,1);			\
  lmk->k62 = lmk_lm2k(6,2);			\
  lmk->k63 = lmk_lm2k(6,3);			\
  lmk->k64 = lmk_lm2k(6,4);			\
  lmk->k65 = lmk_lm2k(6,5);			\
  lmk->k66 = lmk_lm2k(6,6);			\
  lmk->k71 = lmk_lm2k(7,1);			\
  lmk->k72 = lmk_lm2k(7,2);			\
  lmk->k73 = lmk_lm2k(7,3);			\
  lmk->k74 = lmk_lm2k(7,4);			\
  lmk->k75 = lmk_lm2k(7,5);			\
  lmk->k76 = lmk_lm2k(7,6);			\
  lmk->k77 = lmk_lm2k(7,7);			\
  lmk->k81 = lmk_lm2k(8,1);			\
  lmk->k82 = lmk_lm2k(8,2);			\
  lmk->k83 = lmk_lm2k(8,3);			\
  lmk->k84 = lmk_lm2k(8,4);			\
  lmk->k85 = lmk_lm2k(8,5);			\
  lmk->k86 = lmk_lm2k(8,6);			\
  lmk->k87 = lmk_lm2k(8,7);			\
  lmk->k88 = lmk_lm2k(8,8);

