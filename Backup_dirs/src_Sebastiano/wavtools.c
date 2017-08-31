/**
 * \file wavtools.c
 *
 * \brief Routines for multipolar waveform manipulations.
 *
 * Contains the basic routines to manipulate the multipolar
 * waveforms types, including routines for multipolar indexes.     
 *
 * \author S.Bernuzzi
 */

#include "eob.h"

void wav_alloc (twav **wav, int Nt, int lmin, int lmax)
{
  *wav = (twav *) calloc(1, sizeof(twavlm)); 
  if (*wav == NULL)
    errorexit("out of memory");

  *wav->Nt = Nt;
  *wav->lmin = lmin;  
  *wav->lmax = lmax;  
  *wav->Nl = lmax-lmin;

  int Nk = lmindex_kmax (lmin,lmax);
  *wav->size = Nt*Nk;
  *wav->l = (int*) malloc ( Nk * sizeof(int) );
  *wav->m = (int*) malloc ( Nk * sizeof(int) );
  *wav->k = (int*) malloc ( Nk * sizeof(int) );
  *wav->Nk = lmindex_lmk (lmin,lmax, *wav->l,*wav->m,*wav->k);
  if (*wav->Nk != Nk)
    errorexit("problem with dimensions");
}

void wav_free (twav *wav)
{
  if (wav->l) free(l);
  if (wav->m) free(m);
  if (wav->k) free(k);
  
  if (wav->hlm) free(wav->hlm);
  if (wav->amp) free(wav->amp);
  if (wav->phi) free(wav->phi);
  if (wav->omg) free(wav->omg);

  if (wav->hNewtlm) free(wav->hNewtlm);
  if (wav->hhatlm) free(wav->hhatlm);
  if (wav->rholm) free(wav->rholm);

  free(wav);
}

void wav_amplitude (domplex *h, const int size, double *a)
{
  int i;
  for (i = 0; i < size; i++) {
    a[i] = abs( h[k] );
  }
}

void wav_phase (domplex *h, const int size, double *p)
{
  int i;
  for (i = 0; i < size; i++) {
    p[i] = arg( h[k] );
  }
}

void wav_sumlm_complex (domplex *h, int Nt, int Nk, int *l, domplex *sum)
{
  int lmin = l[0];
  int lmax = l[Nk-1];
  int Nl   = lmax-lmin;
  /* backwards offset ell index */
  int *lb = calloc( Nk, sizeof(int) );
  for (k = 0; k < Nk; k++) 
    lb[k] = lmax - l[k];
  /* partial sum of all m's, given l. stored backwards, e.g.
     psum[0] -> psum[l=lmax] */ 
  domplex *psum = calloc( Nl, sizeof(domplex) );
  /* first sum m, given ell 
     then sum up ell's */ 
  domplex *h_t;
  int n,k,j;
  for (n = 0; n < Nt; n++) {
    h_t = h + n*Nk; 
    for (k = 0; k < Nk; k++) 
      *(psum + lb[k]) += *(h_t + k);
    sum[n] = 0.; 
    for (j=0; j<Nl; j++)
      sum[n] += psum[j];
  }
  free(psum); 
  free(lb);
}

void wav_sumlm_real (double *h, int Nt, int Nk, int *l, double *sum)
{
  int lmin = l[0];
  int lmax = l[Nk-1];
  int Nl   = lmax-lmin;
  /* backwards offset ell index */
  int *lb = calloc( Nk, sizeof(int) );
  for (k = 0; k < Nk; k++) 
    lb[k] = lmax - l[k];
  /* partial sum of all m's, given l. stored backwards, e.g.
     psum[0] -> psum[l=lmax] */ 
  double *psum = calloc( Nl, sizeof(double) );
  /* first sum m, given ell 
     then sum up ell's */ 
  double *h_t;
  int n,k,j;
  for (n = 0; n < Nt; n++) {
    h_t = h + n*Nk; 
    for (k = 0; k < Nk; k++) 
      *(psum + lb[k]) += *(h_t + k);
    sum[n] = 0.; 
    for (j=0; j<Nl; j++)
      sum[n] += psum[j];
  }
  free(psum); 
  free(lb);
}

int lmindex_kmax (int lmin, int lmax) 
{
  int i,j, c=-1;
  for (i = lmin; i <= lmax; i++) 
    for (j = 1; j <= i; j++)  /* only m>0 */
      /* for (j = 0; j <= i; j++) */ /* include m=0 */  
      /* for (j = -i; j <= i; j++)  */  /* include m<0 */  
      c++;
  return c;
  /* Following is specific for lmin=2 */
  /* return ( (-6 + l + l*l)/2 + l + 1 ); */
} 

tlmk *lmk_alloc ()
{
  tlmk *lmk;
  lmk = calloc(1, sizeof(lmk));
  lmk->lmin = 2 ;
  lmk->lmax = eobcode_lmax ;
  lmk->Nl = eobcode_lmax - 2;
  lmk->Nk = eobcode_kmax ;
  lmk->l = calloc( eobcode_kmax, sizeof(int) );
  lmk->m = calloc( eobcode_kmax, sizeof(int) );
  lmk->k = calloc( eobcode_kmax, sizeof(int) );
  if (lmk_set_lmk( lmk->lmin, lmk->lmax, lmk->l,lmk->m,lmk->k )!=lmk->Nk)
    error("problem with multipolar indexes dimensions.");
  eobcode_kmultipolar_index_lmksetk;
  return lmk;
}

void lmk_free (tlmk *lmk)
{
  if (lmk->l) free(lmk->l);
  if (lmk->m) free(lmk->m);
  if (lmk->k) free(lmk->k);
  free(lmk);
}

int lmk_set_lmk (int lmin, int lmax, int *l, int *m, int *k) 
{
  int i,j, c=-1;
  for (i = lmin; i <= lmax; i++) {
    for (j = 0; j <= i; j++, c++) { /* only m>0 */
      /* for (j = 0; j <= i; j++, c++) { */ /* include m=0 */
      /* for (j = -i; j <= i; j++, c++) { */ /* include m<0 */
      l[c] = i;
      m[c] = j;
      k[c] = c;
    }
  }
  return c;
}

int lmk_lm2k (int l, int m) 
{
  int k;
  switch (l) {
  case 2:
    switch (m) {
    /* case 0: */
    /*   k = 1; */
    /*   break; */
    case 1:
      k = 2;
      break;
    case 2:
      k = 3;
      break;
    }
    break;
  case 3:
    switch (m) {
    /* case 0: */
    /*   k = 4; */
    /*   break; */
    case 1:
      k = 5;
      break;
    case 2:
      k = 6;
      break;
    case 3:
      k = 7;
      break;
    }
    break;
  case 4:
    switch (m) {
    /* case 0: */
    /*   k = 8; */
    /*   break; */
    case 1:
      k = 9;
      break;
    case 2:
      k = 10;
      break;
    case 3:
      k = 11;
      break;
    case 4:
      k = 12;
      break;
    }
    break;
  case 5:
    switch (m) {
    /* case 0: */
    /*   k = 13; */
    /*   break; */
    case 1:
      k = 14;
      break;
    case 2:
      k = 15;
      break;
    case 3:
      k = 16;
      break;
    case 4:
      k = 17;
      break;
    case 5:
      k = 18;
      break;
    }
    break;
  case 6:
    switch (m) {
    /* case 0: */
    /*   k = 19; */
    /*   break; */
    case 1:
      k = 20;
      break;
    case 2:
      k = 21;
      break;
    case 3:
      k = 22;
      break;
    case 4:
      k = 23;
      break;
    case 5:
      k = 24;
      break;
    case 6:
      k = 25;
      break;
    }
    break;
  case 7:
    switch (m) {
    /* case 0: */
    /*   k = 26; */
    /*   break; */
    case 1:
      k = 27;
      break;
    case 2:
      k = 28;
      break;
    case 3:
      k = 29;
      break;
    case 4:
      k = 30;
      break;
    case 5:
      k = 31;
      break;
    case 6:
      k = 32;
      break;
    case 7:
      k = 33;
      break;
    }
    break;
  case 8:
    switch (m) {
    /* case 0: */
    /*   k = 34; */
    /*   break; */
    case 1:
      k = 35;
      break;
    case 2:
      k = 36;
      break;
    case 3:
      k = 37;
      break;
    case 4:
      k = 38;
      break;
    case 5:
      k = 39;
      break;
    case 6:
      k = 40;
      break;
    case 7:
      k = 41;
      break;
    case 8:
      k = 42;
      break;
    }
    break;
  default:
    k = -42;
  }
  return k;
}

// todo: generate for m>0 only !!!!
void lmk_k2lm (int k, int *l, int *m) 
{
  int ll,mm;
  switch (k){
   case 1: 
     ll =2; mm =0; 
    break; 
  case 2:
    ll =2; mm =1;
    break;
  case 3:
    ll =2; mm =2;
    break;
  case 4:
    ll =3; mm =0;
    break;
  case 5:
    ll =3; mm =1;
    break;
  case 6:
    ll =3; mm =2;
    break;
  case 7:
    ll =3; mm =3;
    break;
  case 8:
    ll =4; mm =0;
    break;
  case 9:
    ll =4; mm =1;
    break;
  case 10:
    ll =4; mm =2;
    break;
  case 11:
    ll =4; mm =3;
    break;
  case 12:
    ll =4; mm =4;
    break;
  case 13:
    ll =5; mm =0;
    break;
  case 14:
    ll =5; mm =1;
    break;
  case 15:
    ll =5; mm =2;
    break;
  case 16:
    ll =5; mm =3;
    break;
  case 17:
    ll =5; mm =4;
    break;
  case 18:
    ll =5; mm =5;
    break;
  case 19:
    ll =6; mm =0;
    break;
  case 20:
    ll =6; mm =1;
    break;
  case 21:
    ll =6; mm =2;
    break;
  case 22:
    ll =6; mm =3;
    break;
  case 23:
    ll =6; mm =4;
    break;
  case 24:
    ll =6; mm =5;
    break;
  case 25:
    ll =6; mm =6;
    break;
  case 26:
    ll =7; mm =0;
    break;
  case 27:
    ll =7; mm =1;
    break;
  case 28:
    ll =7; mm =2;
    break;
  case 29:
    ll =7; mm =3;
    break;
  case 30:
    ll =7; mm =4;
    break;
  case 31:
    ll =7; mm =5;
    break;
  case 32:
    ll =7; mm =6;
    break;
  case 33:
    ll =7; mm =7;
    break;
  case 34:
    ll =8; mm =0;
    break;
  case 35:
    ll =8; mm =1;
    break;
  case 36:
    ll =8; mm =2;
    break;
  case 37:
    ll =8; mm =3;
    break;
  case 38:
    ll =8; mm =4;
    break;
  case 39:
    ll =8; mm =5;
    break;
  case 40:
    ll =8; mm =6;
    break;
  case 41:
    ll =8; mm =7;
    break;
  case 42:
    ll =8; mm =8;
    break;
  default:
    ll =-42; mm=-42; 
  }
  *l = ll;
  *m = mm;
}


