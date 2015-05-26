/**
 * \file waveform.h
 *
 * \brief Function prototypes for EOB waveform.
 *
 * \author S.Bernuzzi
 */

/**
 * \fn domplex **eob_wav_Newt_coefs (double nu)
 *
 * \brief Coefficients of the \f$h^{Newt}_{\ell m}\f$ waveform.
 *
 * \param[in] nu symmetric mass-ratio
 * \return    Allocated complex array hlmNewtCoefs of size = eobcode_kmax
 */
domplex **eob_wav_Newt_coefs (double nu);

/**
 * \fn void eob_wav_hlmNewt_n (double r, double phi, double Omega, 
 *           		       tlmk *lmk, 
 *			       domplex *hlmNewtCoefs,   
 *			       domplex *hlmNewt_n)
 *
 * \brief Newtonian waveform \f$h^\text{Newt}_{\ell m}(v_\phi)\f$. 
 * 
 * The Newtonian waveform is evaluated on the dynamics using the azimuthal velocity \f$v_\phi=r\Omega\f$ 
 * See discussions in
 *    <a href="http://arxiv.org/abs/1212.4357">http://arxiv.org/abs/1212.4357</a> 
 *    <a href="http://arxiv.org/abs/1406.5983">http://arxiv.org/abs/1406.5983</a>
 *
 * \param[in]  r            radial coordinate 
 * \param[in]  phi          orbital phase
 * \param[in]  Omega        orbital frequency \f$\Omega\f$
 * \param[in]  lmk          multipolar indexes 
 * \param[in]  hlmNewtCoefs precomputed coefficients
 * \param[out] hlmNewt      \f$h^\text{Newt}_{\ell m}(v_\phi)\f$ waveform
 * \return     void
 */
void eob_wav_hlmNewt_n (double r, double phi, double Omega, 
			tlmk *lmk, 
			domplex *hlmNewtCoefs,   
			domplex *hlmNewt_n);

/**
 * \fn void eob_wav_hlmNewt (double *r, double *phi, double *Omega, 
 *		             tlmk *lmk, int Nt, 
 *		             domplex *hlmNewtCofes,   
 *		             domplex *hlmNewt)
 *
 * \brief Evaluate Newtonian waveform on time series.
 * 
 * \param[in]  r            radial coordinate 
 * \param[in]  phi          orbital phase
 * \param[in]  Omega        orbital frequency \f$\Omega\f$
 * \param[in]  lmk          multipolar indexes 
 * \param[in]  Nt           timesteps (length of r)
 * \param[in]  hlmNewtCoefs precomputed coefficients
 * \param[out] hlmNewt      \f$h^\text{Newt}_{\ell m}(v_\phi)\f$ waveform
 * \return     void
 */
void eob_wav_hlmNewt (double *r, double *phi, double *Omega, 
		      tlmk *lmk, int Nt, 
		      domplex *hlmNewtCofes,   
		      domplex *hlmNewt)

/**
 * \fn double **eob_wav_rholm_coefs (double nu)
 *
 * \brief Coefficients of \f$\rho_{\ell m}\f$ amplitudes.
 *
 * Damour, Nagar & Bernuzzi, PRD 87, 084035 (April 2013) Appendix D
 * <a href="http://arxiv.org/abs/1212.4357">http://arxiv.org/abs/1212.4357</a>
 *
 * \param[in]  nu symmetric mass ratio
 * \return     Allocated double array rholmCoefs of size = eobcode_kmax x eobcode_wav_rholm_PNcoefs
 */
double **eob_wav_rholm_coefs (double nu);

/**
 * \fn void eob_wav_rholm_n (double x, 
 *		             tlmk *lmk,
 *		             double **rholmCoefs,  
 *		             double *rholm_n)
 *
 * \brief Evaluate \f$\rho_{\ell m}(x)\f$ amplitudes at x.
 * 
 * \param[in]  x          PN parameter
 * \param[in]  lmk        multipolar indexes 
 * \param[in]  rholmCoefs precomputed coefficients 
 * \param[out] rholm      \f$\rho_{\ell m}(x)\f$
 * \return     void
 */
void eob_wav_rholm_n (double x, 
		      tlmk *lmk,
		      double **rholmCoefs,  
		      double *rholm_n);

/**
 * \fn double eob_wav_rho22P23(double nu, double x)
 *
 * \brief Pade (2,3) resummed expression of \f$\rho_{22}\f$.
 *
 * First write the 5PN Taylor expansion as
 * \f$\rho_{22}(x) = 1 + f_1 x + f_2 x^2 + f_3 x^3 + f_4 x^4 + f_5 x^5\f$
 * Then do the Pade. Note the use of the eulerlog22 function.
 * Use Nagar's old implementation.
 * 
 * Damour, Nagar & Bernuzzi, PRD 87, 084035 (April 2013) Appendix D
 * <a href="http://arxiv.org/abs/1212.4357">http://arxiv.org/abs/1212.4357</a>
 *
 * \param[in]  nu     symmetric mass-ratio
 * \param[in]  x      PN expansion parameter
 * \return     Resummed version of \f$\rho_{22}(x)\f$
 */
double eob_wav_rho22P23(double nu, double x);

/**
 * \fn void eob_wav_rholm (double *x, 
 *		           tlmk *lmk, int Nt
 *		           double **rholmCoefs,  
 *		           double *rholm)
 *
 * \brief Evaluate \f$\rho_{\ell m}(x)\f$ on time series.
 * 
 * \param[in]  x          PN parameter
 * \param[in]  lmk        multipolar indexes 
 * \param[in]  Nt         timesteps (length of x)
 * \param[in]  rholmCoefs precomputed coefficients 
 * \param[out] rholm      \f$\rho_{\ell m}(x)\f$
 * \return     void
 */
void eob_wav_rholm (double *x, 
		    tlmk *lmk, int Nt, 
		    double **rholmCoefs,  
		    double *rholm);

/**
 * \fn double **eob_wav_deltalm_coefs (double nu)
 *
 * \brief Coefficients of \f$\delta_{\ell m}\f$ residual phase in Taylor form.
 *
 * Damour, Nagar & Bernuzzi, PRD 87, 084035 (April 2013) Appendix D
 * <a href="http://arxiv.org/abs/1212.4357">http://arxiv.org/abs/1212.4357</a>
 *
 * \param[in]  nu symmetric mass ratio
 * \return     Allocated double array deltalmCoefs of size = eobcode_kmax x eobcode_wav_deltalm_PNcoefs
 */
double **eob_wav_deltalm_coefs (double nu);

/**
 * \fn void eob_wav_deltalm_n (double Heob, double Omega, 
 *                             tlmk *lmk, 
 *			       double *deltalmCoefs,
 *			       double *deltalm_n, domplex *expdlm_n)
 *
 * \brief Residual phase \f$\delta_{\ell m}\f$ of the resummed waveform tail term, Taylor form.
 * 
 * Damour, Iyer & Nagar, PRD 79, 064004 (2009)
 * <a href="http://arxiv.org/abs/0811.2069>http://arxiv.org/abs/0811.2069</a>
 *
 * \param[in]  Heob         Hamiltonian
 * \param[in]  Omega        orbital frequency
 * \param[in]  lmk          multipolar indexes 
 * \param[in]  deltalmCoefs precomputed coefficients 
 * \param[out] deltalm      \f$\delta_{\ell m}\f$
 * \param[out] expdlm       \f$\exp(i \delta_{\ell m})\f$
 * \return     void
 */
void eob_wav_deltalm_n (double Heob, double Omega, 
			tlmk *lmk, 
			double *deltalmCoefs,
			double *deltalm_n, domplex *expdlm_n);

/**
 * \fn void eob_wav_deltalm_resum_n (double Heob, double Omega, 
 *			             tlmk *lmk, 
 *			             double *deltalmCoefs,
 *			             double *deltalm_n, domplex *expdlm_n)
 *
 * \brief Residual phase \f$\delta_{\ell m}\f$ of the resummed waveform tail term, resummed expressions.
 * 
 * Damour, Nagar & Bernuzzi, PRD 87, 084035 (April 2013) Appendix D
 * <a href="http://arxiv.org/abs/1212.4357">http://arxiv.org/abs/1212.4357</a>
 *
 * \param[in]  Heob         Hamiltonian
 * \param[in]  Omega        orbital frequency
 * \param[in]  lmk          multipolar indexes 
 * \param[in]  deltalmCoefs precomputed coefficients 
 * \param[out] deltalm      \f$\delta_{\ell m}\f$
 * \param[out] expdlm       \f$\exp(i \delta_{\ell m})\f$
 * \return     void
 *
 * \warninig  The sym mass-ratio is obtained with a workaround.
 */
void eob_wav_deltalm_resum_n (double Heob, double Omega, 
			      tlmk *lmk, 
			      double *deltalmCoefs,
			      double *deltalm_n, domplex *expdlm_n);

/**
 * \fn void eob_wav_deltalm_n (double Heob, double Omega, 
 *                             tlmk *lmk, 
 *			       double *deltalmCoefs,
 *			       double *deltalm_n, domplex *expdlm_n)
 *
 * \brief Evaluate \f$\delta_{\ell m}\f$ on time series.
 * 
 * Damour, Iyer & Nagar, PRD 79, 064004 (2009)
 * <a href="http://arxiv.org/abs/0811.2069>http://arxiv.org/abs/0811.2069</a>
 *
 * \param[in]  Heob         Hamiltonian
 * \param[in]  Omega        orbital frequency
 * \param[in]  lmk          multipolar indexes 
 * \param[in]  Nt           timesteps (length of Omega)
 * \param[in]  deltalmCoefs precomputed coefficients 
 * \param[out] deltalm      \f$\delta_{\ell m}\f$
 * \param[out] expdlm       \f$\exp(i \delta_{\ell m})\f$
 * \return     void
 */
void eob_wav_deltalm (double *Heob, double *Omega, 
		      tlmk *lmk, int Nt, 
		      double *deltalmCoefs,
		      double *deltalm, domplex *expdlm);

/**
 * \fn void eob_wav_deltalm_n (double Heob, double Omega, 
 *                             tlmk *lmk, 
 *			       double *deltalmCoefs,
 *			       double *deltalm_n, domplex *expdlm_n)
 *
 * \brief Evaluate resummed \f$\delta_{\ell m}\f$ on time series.
 * 
 * Damour, Iyer & Nagar, PRD 79, 064004 (2009)
 * <a href="http://arxiv.org/abs/0811.2069>http://arxiv.org/abs/0811.2069</a>
 *
 * \param[in]  Heob         Hamiltonian
 * \param[in]  Omega        orbital frequency
 * \param[in]  lmk          multipolar indexes 
 * \param[in]  Nt           timesteps (length of Omega)
 * \param[in]  deltalmCoefs precomputed coefficients 
 * \param[out] deltalm      \f$\delta_{\ell m}\f$
 * \param[out] expdlm       \f$\exp(i \delta_{\ell m})\f$
 * \return     void
 */
void eob_wav_deltalm_resum (double *Heob, double *Omega, 
			    tlmk *lmk, int Nt, 
			    double *deltalmCoefs,
			    double *deltalm, domplex *expdlm);

/**
 * \fn double ** eob_wav_NQCCoefs (double nu)
 *
 * \brief Coefficients for NQC waveform corrections. 
 *
 * Return the six multipolar coefficients \f$a^{(\ell m)}_i\f$ and \f$b^{(\ell m)_i\f$ \f$i=1...3\f$
 * as given by the NR fits of 
 * Damour, Nagar & Bernuzzi, PRD 87, 084035 (April 2013) 
 * <a href="http://arxiv.org/abs/1212.4357">http://arxiv.org/abs/1212.4357</a>
 *
 * \param[in]  nu symmetric mass-ratio
 * \return     Allocated double array NQCab \f$(a^{(\ell m)_1,a^{(\ell m)_2,a^{(\ell m)_3,b^{(\ell m)_1,b^{(\ell m)_2,b^{(\ell m)_3)\f$ of size = eobcode_kmax x eobcode_NQC 
 */
double ** eob_wav_NQCCoefs (double nu);

/**
  * \fn void eob_wav_hlmNQC_n (double *NQCn_n,
  *                            tlmk *lmk, 
  *		               double **NQCCoefs,
  *		               domplex *hlmNQC_n);
  * 
  * \brief Next-to-quasi-circular correction to the waveform.
  * 
  * Calculate
  * \f$\hat{h}^\text{NQC}_{\ell m} = (1 + a^{(\ell m)}_i n^i)\exp[i (a^{(\ell m)}_i n^{i+3})]\f$
  * 
  * Damour, Nagar & Bernuzzi, PRD 87, 084035 (April 2013) 
  * <a href="http://arxiv.org/abs/1212.4357">http://arxiv.org/abs/1212.4357</a>
  * 
  * \param[in]  NQCn     NQC dynamical quantities
  * \param[in]  lmk      multipolar indexes  
  * \param[in]  NQCCoefs precomputed coefficients 
  * \param[out] hlmNQC   \f$\hat{h}^\text{NQC}_{\ell m}\f$
  * \return     void
  */
void eob_wav_hlmNQC_n (double *NQCn_n,
		       tlmk *lmk, 
		       double **NQCCoefs,
		       domplex *hlmNQC_n);
/**
  * \fn void eob_wav_hlmNQC (double *NQCn,
  *		             tlmk *lmk, int Nt, 
  *		             double **NQCCoefs,
  *		             domplex *hlmNQC)
  * 
  * \brief Evaluate NQC on time series.
  * 
  * \f$\hat{h}^\text{NQC}_{\ell m} = (1 + a^{(\ell m)}_i n^i)\exp[i (a^{(\ell m)}_i n^{i+3})]\f$
  * Damour, Nagar & Bernuzzi, PRD 87, 084035 (April 2013) 
  * <a href="http://arxiv.org/abs/1212.4357">http://arxiv.org/abs/1212.4357</a>
  * 
  * \param[in]  NQCn     NQC dynamical quantities
  * \param[in]  lmk      multipolar indexes  
  * \param[in]  Nt       timesteps 
  * \param[in]  NQCCoefs precomputed coefficients 
  * \param[out] hlmNQC   \f$\hat{h}^\text{NQC}_{\ell m}\f$
  * \return     void
 */
void eob_wav_hlmNQC (double *NQCn,
		     tlmk *lmk, int Nt, 
		     double **NQCCoefs,
		     domplex *hlmNQC);

/**
 * \fn void eob_wav_Tlm_n (double Omega, double E, 
 *    		           tlmk *lmk, 
 *		           double r0, 
 *		           domplex *Tlm_n);
 *
 * \brief Tail contribution \f$T_{\ell m}\f$ of the resummed waveform.
 *
 * Damour, Iyer & Nagar, PRD 79, 064004 (2009)
 * <a href="http://arxiv.org/abs/0811.2069>http://arxiv.org/abs/0811.2069</a>
 *
 * \param[in]  Omega orbital frequency \f$\Omega\f$
 * \param[in]  E     energy
 * \param[in]  lmk   multipolar indexes  
 * \param[in]  r0    parameter
 * \param[out] Tlm   \f$T_{\ell m}\f$ of size Nt x Nk
 * \return     void
 */
void eob_wav_Tlm_n (double Omega, double E, 
		    tlmk *lmk, 
		    double r0, 
		    domplex *Tlm_n);

/**
 * \fn void eob_wav_Tlm (double Omega, double E, 
 *    		         tlmk *lmk, 
 *		         double r0, 
 *		         domplex *Tlm_n);
 *
 * \brief Evaluate \f$T_{\ell m}\f$ on time series.
 *
 * Damour, Iyer & Nagar, PRD 79, 064004 (2009)
 * <a href="http://arxiv.org/abs/0811.2069>http://arxiv.org/abs/0811.2069</a>
 *
 * \param[in]  Omega orbital frequency \f$\Omega\f$
 * \param[in]  E     energy
 * \param[in]  lmk   multipolar indexes
 * \param[in]  Nt    timesteps (length of Omega)
 * \param[in]  r0    parameter
 * \param[out] Tlm   \f$T_{\ell m}\f$ of size Nt x Nk
 * \return     void
 */
void eob_wav_Tlm (double *Omega, double *E, 
		  tlmk *lmk, int Nt, 
		  double r0, 
		  domplex *Tlm);

/**
 * \fn void eob_wav_Seff_n (double Heff, double jhat,
 *	     	            tlmk *lmk,
 *		            double *Seff_n);
 *
 * \brief Source factor \f$S^{\epsilon}_\text{eff}\f$ of the resummed waveform.
 *
 * Damour, Iyer & Nagar, PRD 79, 064004 (2009)
 * <a href="http://arxiv.org/abs/0811.2069>http://arxiv.org/abs/0811.2069</a>
 *
 * \param[in]  Heff  Effective Hamiltonian \f$\hat{H}_\text{eff}\f$
 * \param[in]  jhat  \f$p_\phi/(r_\omega v_\phi)\f$
 * \param[in]  lmk   multipolar indexes
 * \param[out] Seff  \f$S^{\epsilon}_\text{eff}\f$ of size Nt x Nk
 * \return     void
 */
void eob_wav_Seff_n (double Heff, double jhat,
		     tlmk *lmk,
		     double *Seff_n);

/**
 * \fn void eob_wav_Seff (double *Heff, double *jhat,
 *	     	          tlmk *lmk, int Nt,
 *		          double *Seff_n);
 *
 * \brief Evaluate \f$S^{\epsilon}_\text{eff}\f$ on time series.
 *
 * Damour, Iyer & Nagar, PRD 79, 064004 (2009)
 * <a href="http://arxiv.org/abs/0811.2069>http://arxiv.org/abs/0811.2069</a>
 *
 * \param[in]  Heff  Effective Hamiltonian \f$\hat{H}_\text{eff}\f$
 * \param[in]  jhat  \f$p_\phi/(r_\omega v_\phi)\f$
 * \param[in]  lmk   multipolar indexes
 * \param[in]  Nt    timesteps (length of Omega)
 * \param[out] Seff  \f$S^{\epsilon}_\text{eff}\f$ of size Nt x Nk
 * \return     void
 */
void eob_wav_Seff (double *Heff, double *jhat, 
		   tlmk *lmk, int Nt, 
		   double *Seff);

/**
 * \fn void eob_wav_flm_n (double *rholm_n, 
 *     		           tlmk *lmk, 
 *		           double *flm_n);
 *
 * \brief Amplitudes \f$f_{\ell m}(x)\f$ of the resummed waveform.
 * 
 * Damour, Iyer & Nagar, PRD 79, 064004 (2009)
 * <a href="http://arxiv.org/abs/0811.2069>http://arxiv.org/abs/0811.2069</a>
 *
 * \param[in]  rholm \f$\rho_{\ell m}(x)\f$ 
 * \param[in]  lmk   multipolar indexes
 * \param[out] flm   \f$f_{\ell m}(x) = \rho_{\ell m}(x)^\ell\f$
 * \return     void
 */
void eob_wav_flm_n (double *rholm_n, 
		    tlmk *lmk, 
		    double *flm_n);

/**
 * \fn void eob_wav_flm (double *rholm, 
 *     		         tlmk *lmk, 
 *		         double *flm);
 *
 * \brief  Evaluate \f$f_{\ell m}(x)\f$ on time series.
 * 
 * Damour, Iyer & Nagar, PRD 79, 064004 (2009)
 * <a href="http://arxiv.org/abs/0811.2069>http://arxiv.org/abs/0811.2069</a>
 *
 * \param[in]  rholm \f$\rho_{\ell m}(x)\f$ 
 * \param[in]  lmk   multipolar indexes
 * \param[in]  Nt    timesteps 
 * \param[out] flm   \f$f_{\ell m}(x) = \rho_{\ell m}(x)^\ell\f$
 * \return     void
 */
void eob_wav_flm (double *rholm,
		  tlmk *lmk, int Nt, 
		  double *flm);
/**
 * \fn double eob_wav_rholm_xarg (double r, double Omega, double x)
 *     
 * \brief Argument for \f$\rho_{\ell m}(x)\f$.
 *
 * \param[in] r     radial coordinate 
 * \param[in] Omega orbital frequency \f$\Omega\f$
 * \param[in] x     \f$x=r_\omega\Omega=(r\psi^{1/3})\omega\f$
 * \return    xarg  
 */
double eob_wav_rholm_xarg (double r, double Omega, double x);


