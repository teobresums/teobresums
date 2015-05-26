/**
 * \file metric.h
 * \brief Function prototypes for EOB metric.
 *
 * \author S.Bernuzzi
 */

/**
 * \fn void eob_met_A_5PNP15_n (double u, void *coefs, 
 *                              double *A, double *dA_u, double *d2A_u)
 *
 * \brief Metric potential \f$A(u)\f$, 5PN expression resummed Pade(1,5).
 *
 * \param[in]  u      1/r coordinate
 * \param[in]  ACoefs PN coefficients
 * \param[out] A      \f$A(u)\f$
 * \param[out] dA     \f$A'(u)\f$
 * \param[out] d2A    \f$A''(u)\f$
 * \return     void
 */
void eob_met_A_5PNP15_n (double u, void *coefs,
			 double *A, double *dA_u, double *d2A_u);

/**
 * void eob_met_A_4PNP14_n (double u, void *coefs, 
 *                          double *A, double *dA_u, double *d2A_u)
 *
 * \brief Metric potential \f$A(u)\f$, 4PN expression resummed Pade(1,4).
 *
 * \param[in]  u      1/r coordinate
 * \param[in]  ACoefs PN coefficients
 * \param[out] A      \f$A(u)\f$
 * \param[out] dA     \f$A'(u)\f$
 * \param[out] d2A    \f$A''(u)\f$
 * \return     void
 */
void eob_met_A_4PNP14_n (double u, void *coefs,
			 double *A, double *dA_u, double *d2A_u);

/**
 * void eob_met_A_3PNP13_n (double u, void *coefs,
 *   			    double *A, double *dA_u, double *d2A_u)
 *
 * \brief Metric potential \f$A(u)\f$, 3PN expression resummed Pade(1,3).
 *
 * \param[in]  u      1/r coordinate
 * \param[in]  ACoefs PN coefficients
 * \param[out] A      \f$A(u)\f$
 * \param[out] dA     \f$A'(u)\f$
 * \param[out] d2A    \f$A''(u)\f$
 * \return     void
 */
void eob_met_A_3PNP13_n (double u, void *coefs,
			 double *A, double *dA_u, double *d2A_u);

/**
 * \fn void eob_met_A_5PN_n (double u, void *coefs,
 *   		             double *A, double *dA_u, double *d2A_u)
 *
 * \brief Metric potential \f$A(u)\f$, 5PN expression.
 *
 * \param[in]  u      1/r coordinate
 * \param[in]  ACoefs PN coefficients
 * \param[out] A      \f$A(u)\f$
 * \param[out] dA     \f$A'(u)\f$
 * \param[out] d2A    \f$A''(u)\f$
 * \return     void
 */
void eob_met_A_5PN_n (double u, void *coefs,
		      double *A, double *dA_u, double *d2A_u);

/**
 * \fn void eob_met_A_4PN_n (double u, void *coefs,
 *   		             double *A, double *dA_u, double *d2A_u)
 *
 * \brief Metric potential \f$A(u)\f$, 4PN expression.
 *
 * \param[in]  u      1/r coordinate
 * \param[in]  ACoefs PN coefficients
 * \param[out] A      \f$A(u)\f$
 * \param[out] dA     \f$A'(u)\f$
 * \param[out] d2A    \f$A''(u)\f$
 * \return     void
 */
void eob_met_A_4PN_n (double u, void *coefs,
		      double *A, double *dA_u, double *d2A_u);

/**
 * \fn void eob_met_A_3PN_n (double u, void *coefs,
 * 		             double *A, double *dA_u, double *d2A_u)
 *
 * \brief Metric potential \f$A(u)\f$, 3PN expression.
 *
 * \param[in]  u      1/r coordinate
 * \param[in]  ACoefs PN coefficients
 * \param[out] A      \f$A(u)\f$
 * \param[out] dA     \f$A'(u)\f$
 * \param[out] d2A    \f$A''(u)\f$
 * \return     void
 */
void eob_met_A_3PN_n (double u, void *coefs,
		      double *A, double *dA_u, double *d2A_u);

/**
 * \fn void eob_met_A_2PN_n (double u, void *coefs,
 *		      double *A, double *dA_u, double *d2A_u)
 *
 * \brief Metric potential \f$A(u)\f$, 2PN expression.
 *
 * \param[in]  u      1/r coordinate
 * \param[in]  ACoefs PN coefficients
 * \param[out] A      \f$A(u)\f$
 * \param[out] dA     \f$A'(u)\f$
 * \param[out] d2A    \f$A''(u)\f$
 * \return     void
 */
void eob_met_A_2PN_n (double u, void *coefs,
		      double *A, double *dA_u, double *d2A_u);

/**
 * \fn void eob_met_A_1PN_n (double u, void *coefs,
 *		      double *A, double *dA_u, double *d2A_u)
 *
 * \brief Metric potential \f$A(u)\f$, 1PN expression.
 *
 * \param[in]  u      1/r coordinate
 * \param[in]  ACoefs PN coefficients
 * \param[out] A      \f$A(u)\f$
 * \param[out] dA     \f$A'(u)\f$
 * \param[out] d2A    \f$A''(u)\f$
 * \return     void
 */
void eob_met_A_1PN_n (double u, void *coefs,
		      double *A, double *dA_u, double *d2A_u);

/**
 * \fn eob_met_dA_rfromu_n (double u, double dA_u)
 *
 * \brief Return \f$A'(r)\f$ from \f$A'(u)\f$
 *
 * \param[in] u    1/r coordinate
 * \param[in] dA_u \f$A'(u)\f$
 * \return    \f$A'(r)\f$
 */
double eob_met_dA_rfromu_n (double u, double dA_u);

/**
 * \fn double eob_met_d2A_rfromu_n (double u, double dA_u, double d2A_u)
 *
 * \brief Return \f$A''(r)\f$ from \f$A'(u)\f$ and \f$A''(u)\f$
 *
 * \param[in] u     1/r coordinate
 * \param[in] dA_u  \f$A'(u)\f$
 * \param[in] d2A_u \f$A''(u)\f$
 * \return    \f$A''(r)\f$
 */
double eob_met_d2A_rfromu_n (double u, double dA_u, double d2A_u);

/**
 * \fn double *eob_met_ACoefs (double nu)
 *
 * \brief Coefficients of the \f$A(r)\f$ metric potential.
 *
 * \param[in]  nu symmetric mass-ratio
 * \return     Allocated double array ACoefs of size = eobcode_met_A_PNcoefs
 */
double *eob_met_ACoefs (double nu);

/**
 * \fn void eob_met_B_fromAD_n (double A, double dA, double D, double dD, 
 *			        double *B, double *dB)
 *
 * \brief  Metric potential \f$B(u)\f$ from \f$A(u)\f$ and \f$D(u)\f$.
 *
 * \param[in]  A   \f$A(u)\f$
 * \param[in]  dA  \f$A'(u)\f$
 * \param[in]  D   \f$D(u)\f$
 * \param[in]  dD  \f$D'(u)\f$
 * \param[out] B   \f$B(u)\f$
 * \param[out] dB  \f$B'(u)\f$
 */
void eob_met_B_fromAD_n (double A, double dA, double D, double dD, 
			 double *B, double *dB);

/**
 * \fn void eob_met_B_1PN_n (double u, double *BCoefs,
 *           		     double *D, double *dD)
 *
 * \brief Metric potential \f$B(u)\f$, 1PN expression.
 *
 * \param[in]  u      1/r coordinate
 * \param[in]  DCoefs PN coefficients
 * \param[out] B      \f$B(u)\f$
 * \param[out] dB     \f$B'(u)\f$
 * \return     void
 */
void eob_met_B_1PN_n (double u, double *BCoefs,
		      double *D, double *dD);

/**
 * \fn void eob_met_D_3PNP03_n (double u, double *DCoefs, double *D, double *dD)
 *
 * \brief Metric potential \f$D(u)\f$, 3PN expression resummed Pade(0,3).
 *
 * \param[in]  u      1/r coordinate
 * \param[in]  DCoefs PN coefficients
 * \param[out] D      \f$D(u)\f$
 * \param[out] dD     \f$D'(u)\f$
 * \return     void
 */
void eob_met_D_3PNP03_n (double u, double *DCoefs,
			 double *D, double *dD);

/**
 * \fn void eob_met_D_3PN_n (double u, double *DCoefs, double *D, double *dD)
 *
 * \brief Metric potential \f$D(r)\f$, 3PN expression.
 *
 * \param[in]  u      1/r coordinate
 * \param[in]  DCoefs PN coefficients
 * \param[out] D      \f$D(r)\f$
 * \param[out] dD     \f$D'(r)\f$
 * \return     void
 */
void eob_met_D_3PN_n (double u, double *DCoefs,
		      double *D, double *dD);

/**
 * \fn void eob_met_D_2PN (double u, double *DCoefs, double *D, double *dD)
 *
 * \brief Metric potential \f$D(u)\f$, 2PN expression.
 *
 * \param[in]  u      1/r coordinate
 * \param[in]  DCoefs PN coefficients
 * \param[out] D      \f$D(u)\f$
 * \param[out] dD     \f$D'(u)\f$
 * \return     void
 */
void eob_met_D_2PN_n (double u, double *DCoefs,
		      double *D, double *dD);

/**
 * \fn void eob_met_D_1PN (double u, double *DCoefs, double *D, double *dD)
 *
 * \brief Metric potential \f$D(u)\f$, 1PN expression.
 *
 * \param[in]  u      1/r coordinate
 * \param[in]  DCoefs PN coefficients
 * \param[out] D      \f$D(u)\f$
 * \param[out] dD     \f$D'(u)\f$
 * \return     void
 */
void eob_met_D_1PN_n (double u, double *DCoefs,
		      double *D, double *dD);

/**
 * \fn double *eob_met_DCoefs (double nu)
 *
 * \brief Coefficients of the \f$D(r)\f$ metric potential, up ot 3PN.
 *
 * \param[in]  nu symmetric mass ratio
 * \return     Allocated double array DCoefs of size = eobcode_met_D_PNcoefs
 */
double *eob_met_DCoefs (double nu);
