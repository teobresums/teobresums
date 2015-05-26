/**
 * \file metric.c
 *
 * \brief Routines to construct the EOB metric.
 *
 * \author S.Bernuzzi
 */

#include "eob.h"

void eob_met_A_5PNP15_n (double u, void *coefs,
			 double *A, double *dA_u, double *d2A_u)
{
  tPade *P = (tPade *)coefs;
  double u2 = SQ(u);
  double u3 = u2*u; 
  double u4 = u3*u; 
  double u5 = u4*u; 
  double u6 = u5*u; 
  double N,D, dN,dD, d2N,d2D, ooD;

  N = P->n[1]*u + P->n[0];
  D = (((((((P->d[6]*u6)+(P->d[5]*u5))  + (P->d[4]*u4)) + (P->d[3]*u3)) + (P->d[2]*u2)) +  (P->d[1]*u)) + P->d[0]);
  ooD = 1./D;

  dN = P->d1n[0];
  dD = ((((((P->d1d[5]*u5) + (P->d1d[4]*u4)) + (P->d1d[3]*u3)) + (P->d1d[2]*u2)) +  (P->d1d[1]*u)) + P->d1d[0]);
  
  *A = N*ooD;
  *dA_u = (-N*dD + D*dN)*SQ(ooD);
  
  if (d2A_u!=NULL) {	  
    d2N = 0.;
    d2D = (((((P->d2d[4]*u4) + (P->d2d[3]*u3)) + (P->d2d[2]*u2)) +  (P->d2d[1]*u)) + P->d2d[0]);
    *d2A_u = (2*N*SQ(dD) - 2*D*dD*dN - D*N*d2D + SQ(D)*d2N)*SQ(ooD)*ooD;
  }  
}

void eob_met_A_4PNP14_n (double u, void *coefs,
			  double *A, double *dA_u, double *d2A_u)
{
  tPade *P = (tPade *)coefs;
  double u2 = SQ(u);
  double u3 = u2*u; 
  double u4 = u3*u; 
  double u5 = u4*u; 
  double N,D, dN,dD, d2N,d2D, ooD;

  N = P->n[1]*u + P->n[0];
  D = ((((((P->d[5]*u5)+(P->d[4]*u4)) + (P->d[3]*u3)) + (P->d[2]*u2)) +  (P->d[1]*u)) + P->d[0]);
  ooD = 1./D;

  dN = P->d1n[0];
  dD = (((((P->d1d[4]*u4) + (P->d1d[3]*u3)) + (P->d1d[2]*u2)) +  (P->d1d[1]*u)) + P->d1d[0]);
  
  *A = N*ooD;
  *dA_u = (-N*dD + D*dN)*SQ(ooD);
  
  if (d2A_u!=NULL) {	  
    d2N = 0.;
    d2D = ((((P->d2d[3]*u3) + (P->d2d[2]*u2)) +  (P->d2d[1]*u)) + P->d2d[0]);
    *d2A_u = (2*N*SQ(dD) - 2*D*dD*dN - D*N*d2D + SQ(D)*d2N)*SQ(ooD)*ooD;
  }  
}

void eob_met_A_3PNP13_n (double u, void *coefs,
			 double *A, double *dA_u, double *d2A_u)
{
  tPade *P = (tPade *)coefs;
  double u2 = SQ(u);
  double u3 = u2*u; 
  double u4 = u3*u; 
  double N,D, dN,dD, d2N,d2D, ooD;

  N = P->n[1]*u + P->n[0];
  D = (((((P->d[4]*u4) + (P->d[3]*u3)) + (P->d[2]*u2)) +  (P->d[1]*u)) + P->d[0]);
  ooD = 1./D;

  dN = P->d1n[0];
  dD = ((((P->d1d[3]*u3) + (P->d1d[2]*u2)) +  (P->d1d[1]*u)) + P->d1d[0]);
  
  *A = N*ooD;
  *dA_u = (-N*dD + D*dN)*SQ(ooD);
  
  if (d2A_u!=NULL) {	  
    d2N = 0.;
    d2D = (((P->d2d[2]*u2) +  (P->d2d[1]*u)) + P->d2d[0]);
    *d2A_u = (2*N*SQ(dD) - 2*D*dD*dN - D*N*d2D + SQ(D)*d2N)*SQ(ooD)*ooD;
  }
}

void eob_met_A_5PN_n (double u, void *coefs,
		      double *A, double *dA_u, double *d2A_u)
{
  double *ACoefs = (double *)coefs;
  double a0 = ACoefs[0];
  double a1 = ACoefs[1];
  double a3 = ACoefs[2];
  double a4 = ACoefs[3];
  double a5c = ACoefs[5];
  double a5l = ACoefs[6];
  double a6c = ACoefs[7];
  double a6l = ACoefs[8];
  double logu = log(u);
  double a5 = a5c + a5l*logu;
  double a6 = a6c + a6l*logu;
  double u2 = SQ(u);
  double u3 = u2*u; 
  double u4 = u3*u; 
  double u5 = u4*u; 
  double u6 = u5*u; 
  *A = (a0 + 
	((a1*u) + 
	 ((a3*u3) + 
	  ((a4*u4) +
	   ((a5*u5) +
	    (a6*u6))))));
  if (dA_u!=NULL) *dA_u = (a1 + 
			   ((3.*a3*u2) + 
			    ((4*a4*u3) + 
			     (((5.*a5 + a5l)*u4) + 
			      ((6.*a6 + a6l)*u5) ))));
  if (d2A_u!=NULL) *d2A_u = ((6.*a3*u) +
			     ((12.*a4*u2) + 
			      ((((20.*a5) + (9*a5l))*u3) + 
			       (((20.*a6) + (11*a6l))*u4))));
}

void eob_met_A_4PN_n (double u, void *coefs,
		      double *A, double *dA_u, double *d2A_u)
{
  double *ACoefs = (double *)coefs;
  double a0 = ACoefs[0];
  double a1 = ACoefs[1];
  double a3 = ACoefs[2];
  double a4 = ACoefs[3];
  double a5c = ACoefs[5];
  double a5l = ACoefs[6];
  double logu = log(u);
  double a5 = a5c + a5l*logu;
  double u2 = SQ(u);
  double u3 = u2*u; 
  double u4 = u3*u; 
  double u5 = u4*u; 
  *A = (a0 +
	((a1*u) +
	 ((a3*u3) + 
	  ((a4*u4) + 
	   (a5*u5)))));
  if (dA_u!=NULL) *dA_u = (a1 +
			   ((3.*a3*u2) +
			    ((4.*a4*u3) +
			     ((5.*a5 + a5l)*u4))));
  if (d2A_u!=NULL) *d2A_u = ((6.*a3*u) + 
			     ((12.*a4*u2) + 
			      (((20.*a5) + (9*a5l))*u3)));

}

void eob_met_A_3PN_n (double u, void *coefs,
		       double *A, double *dA_u, double *d2A_u)
{
  double *ACoefs = (double *)coefs;
  double a0 = ACoefs[0];
  double a1 = ACoefs[1];
  double a3 = ACoefs[2];
  double a4 = ACoefs[3];
  double u2 = SQ(u);
  double u3 = u2*u; 
  double u4 = u3*u; 
  *A = (a0 +
	((a1*u) +
	 ((a3*u3) +
	  (a4*u4))));
  if (dA_u!=NULL) *dA_u = (a1 + ((4*a4*u3) + (3.*a3*u2)));
  if (d2A_u!=NULL) *d2A_u = ((6.*a3*u) + (12.*a4*u2));
}

void eob_met_A_2PN_n (double u, void *coefs,
		      double *A, double *dA_u, double *d2A_u)
{
  double *ACoefs = (double *)coefs;
  double a0 = ACoefs[0];
  double a1 = ACoefs[1];
  double a3 = ACoefs[2];
  double u2 = SQ(u);
  double u3 = u2*u; 
  *A = (a0 + (a1*u + (a3*u3)));
  if (dA_u!=NULL) *dA_u = a1 + (3.*a3*u2);
  if (d2A_u!=NULL) *d2A_u = 6.*a3*u;
}

void eob_met_A_1PN_n (double u, void *coefs,
		      double *A, double *dA_u, double *d2A_u)
{
  double *ACoefs = (double *)coefs;
  double a0 = ACoefs[0];
  double a1 = ACoefs[1];
  *A = ((a1*u) + a0);
  if (dA_u!=NULL) *dA_u = a1;
  if (d2A_u!=NULL) *d2A_u = 0.;
}

double eob_met_dA_rfromu_n (double u, double dA_u) 
{
  double u2 = SQ(u); 
  return -u2.*dA_u;
}

double eob_met_d2A_rfromu_n (double u, double dA_u, double d2A_u) 
{
  double u3 = pow(u,3); 
  double u4 = u3*u;
  return u4*d2A_u + 2.*u3*dA_u;
}

double *eob_met_ACoefs (double nu)
{ 
  double *ACoefs;
  ACoefs = calloc( eobcode_met_A_PNcoefs, sizeof(double) );
  ACoefs[0] = 1.;
  /* 1PN 
   */
  ACoefs[1] = -2.;
  /* 2PN
     a3 u^3 */
  ACoefs[2] = nu*2;
  /* 3PN
     a4 u^4 */
  ACoefs[3] = nu*(94/3 - 41/32*SQ(Pi));
  /* 4PN
     (a5c + a5l*logu) u^5 */  
  ACoefs[4] = nu*( -4237./60. + 2275./512.*SQ(Pi) + 256./5.*Log2 + 128./5.*EulerGamma
		   +nu*(-221/6+41/32*pi^2) ); /* a5c */
  ACoefs[5] = nu*64./5.;/* a5l */
  /* 5PN
     (a6c + a6l*logu) u^6 */  
  ACoefs[6] = 0.; /* a6c, unknown */
  ACoefs[7] = nu*(-7004/105-144/5*nu);/* a5c */
  return ACoefs;
}

void eob_met_B_fromAD_n (double A, double dA, double D, double dD, 
			 double *B, double *dB)
{
  double ooA = 1./A;
  *B   = D*ooA;
  *dB  = (dD*A - D*dA)*SQ(ooA);
}

void eob_met_B_1PN_n (double u, double *BCoefs,
		       double *B, double *dB)
{
  double mtu = -2*u;
  double A = 1+tu;
  double ooA = 1./A;
  *B = ooA;
  *dB = SQ(ooA)*mtu;
}

void eob_met_D_3PNP03_n (double u, double *DCoefs,
			 double *D, double *dD)
{
  double T,dT;
  eob_met_D_3PN (u, DCoefs, &T,&dT);
  *D  = 1./T;
  *dD = -dT*(SQ(*D)); 
}

void eob_met_D_3PN_n (double u, double *DCoefs,
		      double *D, double *dD)
{
  double u2 = SQ(u);
  double u3 = u2*u;
  double u4 = u3*u;
  double d0 = DCoefs[0];
  double d2 = DCoefs[1];
  double d3 = DCoefs[2];
  *D  = (d0 + ((d2*u2) + (d3*u3)));
  *dD = - ((2.*d2*u3) + (3.*d3*u4));
}

void eob_met_D_2PN_n (double u, double *DCoefs,
		      double *D, double *dD)
{
  double u2 = SQ(u);
  double u3 = u2*u;
  double d0 = DCoefs[0];
  double d2 = DCoefs[1];
  *D  = d0 + d2*u2;
  *dD = -2.*d2*u3;
}

void eob_met_D_1PN_n (double u, double *DCoefs,
		      double *D, double *dD)
{
  *D  = 1.;
  *dD = 0.;
}

double *eob_met_DCoefs (double nu)
{
  double *DCoefs;
  DCoefs = calloc( eobcode_met_D_PNcoefs, sizeof(double) );
  DCoefs[0] = 1.;
  DCoefs[1] = 6.*nu;
  DCoefs[2] = 2*nu*(26-3*nu);
  return DCoefs;
}
