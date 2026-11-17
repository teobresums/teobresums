/** \file TEOBResumSPostAdiabEcc.c
 *  \brief Post-adiabatic dynamics for bound eccentric orbits, in equinoctial elements
 *
 *  Equations cited from eccpa paper as Eq. (number) [label]. Notation, code -> paper:
 *    J                 pphi
 *    U                 r_a^2 = alpha^2 + beta^2 (derivatives in width are taken in U)
 *    e_geo             r_a/r_c
 *    theta             Lambda = X + T and lambda = phi - V, Eq. (44) [eq:TV_maps]
 *    omega[th]         their frequencies n and (1 + k) n
 *    F, F_th           rates at one point, Eqs. (16), (41), App. A
 *    Fbar, Fbar_th     orbit averages, Eq. (B1); Fbar2, Fbar3, Fbarcorr: higher orders
 *    ctilde1, thtilde1 oscillating parts c~(1), theta~(1), Eqs. (26), (29)
 *    kappa_m           m-th Fourier coefficient in chi over r_a^m, Eq. (C6)
 *    nodes             X_j = w + 2 pi (j + 1/2)/N, Eq. (C4)
 *    checkpoints       points of LO evolution where corr2_*, corr3_* are computed (App. C)
 */
#include <string.h>
#include <stdlib.h>
#include "TEOBResumS.h"

/* PostAdDynamics: allocation and release */

/** Free PostAdDynamics and its arrays */
void PostAdDynamics_free(PostAdDynamics **pa)
{
  PostAdDynamics *p = *pa;
  if (p == NULL) return;
  free(p->time);
  for (int v = 0; v < EVO_NVARS; v++) free(p->data[v]);
  free(p->evo_M);
  for (int a = 0; a < OSC_N; a++)
    for (int c = 0; c < 2; c++)
      for (int m = 0; m <= ECCPA_HARMONICS_MAX; m++) free(p->evo_osc1_coef[a][c][m]);
  for (int th = 0; th < TH_N; th++)
    for (int c = 0; c < 2; c++)
      for (int m = 0; m <= ECCPA_HARMONICS_MAX; m++) {
        free(p->evo_kappa[th][c][m]);
        free(p->evo_kappa_dU[th][c][m]);
      }
  free(p->fine_arrays);
  free(p);
  *pa = NULL;
}

/**
 * Function: PostAdDynamics_alloc
 * ------------------------------
 *   Allocate PostAdDynamics (calloc'd)
 *
 *   @param[out] pa   : struct, NULL if out of memory
 *   @param[in]  size : number of evolution samples arrays can hold
 *   @return  0 on success, 1 if out of memory
 */
int PostAdDynamics_alloc(PostAdDynamics **pa, int size)
{
  *pa = (PostAdDynamics *) calloc(1, sizeof(PostAdDynamics));
  if (*pa == NULL) return 1;
  PostAdDynamics *p = *pa;

  // largest number of nodes an orbit may be sampled at.
  // PA may use fewer, initial data keeps it at max
  p->nmax = ECCPA_NODES_MAX;

  p->size = size;
  int ok = 1;
  p->time = (double *) calloc(size, sizeof(double));
  ok = ok && p->time;
  for (int v = 0; v < EVO_NVARS; v++) {
    p->data[v] = (double *) calloc(size, sizeof(double));
    ok = ok && p->data[v];
  }
  p->evo_M = (int *) calloc(size, sizeof(int));
  ok = ok && p->evo_M;
  for (int a = 0; a < OSC_N; a++)
    for (int c = 0; c < 2; c++)
      for (int m = 0; m <= ECCPA_HARMONICS_MAX; m++) {
        p->evo_osc1_coef[a][c][m] = (double *) calloc(size, sizeof(double));
        ok = ok && p->evo_osc1_coef[a][c][m];
      }
  for (int th = 0; th < TH_N; th++)
    for (int c = 0; c < 2; c++)
      for (int m = 0; m <= ECCPA_HARMONICS_MAX; m++) {
        p->evo_kappa[th][c][m]    = (double *) calloc(size, sizeof(double));
        p->evo_kappa_dU[th][c][m] = (double *) calloc(size, sizeof(double));
        ok = ok && p->evo_kappa[th][c][m] && p->evo_kappa_dU[th][c][m];
      }
  if (!ok) {
    PostAdDynamics_free(pa);
    return 1;
  }
  return 0;
}

/** Eight-point Gauss-Legendre rule on [-1, 1] (Abramowitz & Stegun, Table 25.4) */
static const double eccpa_gl8_x[8] = {
  -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
   0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363
};
static const double eccpa_gl8_w[8] = {
   0.1012285362903763,  0.2223810344533745,  0.3137066458778873,  0.3626837833783620,
   0.3626837833783620,  0.3137066458778873,  0.2223810344533745,  0.1012285362903763
};

/* Helpers wrapping Hamiltonian calls */

/** dH/dr at (r, prstar = 0, pphi). */
static double eccpa_cons_Hr(double r, double pphi, Dynamics *dyn)
{
  double dH_dr;
  eob_ecc_ham(r, 0.0, pphi, dyn, NULL, &dH_dr, NULL, NULL, NULL, NULL, NULL,
              NULL, NULL, NULL, NULL, NULL, NULL);
  return dH_dr;
}

/** dH/d(prstar^2) = (dH/dprstar)/(2 prstar). prstar = 0 uses l'Hopital rule */
static double eccpa_cons_dH_dprs2(double r, double prstar, double pphi, Dynamics *dyn)
{
  if (prstar <= 0.0) {   /* callers pass a square root, so this is prstar = 0 exactly */
    double d2H0;
    eob_ecc_ham(r, 0.0, pphi, dyn, NULL, NULL, NULL, NULL, NULL, &d2H0, NULL,
                NULL, NULL, NULL, NULL, NULL, NULL);
    return 0.5*d2H0;
  }
  double dH_dprstar;
  eob_ecc_ham(r, prstar, pphi, dyn, NULL, NULL, &dH_dprstar, NULL, NULL, NULL, NULL,
              NULL, NULL, NULL, NULL, NULL, NULL);
  return dH_dprstar/(2.0*prstar);
}

/**
 * Function: eccpa_cons_deltaH
 * ---------------------------
 *   deltaH(r_c) of Eq. (C1) and its derivative H_r(r_c + r_a) - H_r(r_c - r_a).
 *   Direct difference of two H, or integral of H_r at small e (App. C).
 *
 *   @return  1 if inner turning point (or a node) is not above r = 1, 0 otherwise
 */
static int eccpa_cons_deltaH(double J, double r_c, double r_a, int use_direct, Dynamics *dyn,
                             double *deltaH, double *ddeltaH_drc)
{
  if (use_direct) {
    double Hp, Hm, Hrp, Hrm;
    if (!(r_c - r_a > 1.0)) return 1;
    eob_ecc_ham(r_c + r_a, 0.0, J, dyn, &Hp, &Hrp, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL, NULL);
    eob_ecc_ham(r_c - r_a, 0.0, J, dyn, &Hm, &Hrm, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL, NULL);
    *deltaH      = Hp - Hm;
    *ddeltaH_drc = Hrp - Hrm;
    return 0;
  }
  double sum = 0.0;
  for (int k = 0; k < 8; k++) {
    const double s = r_c + r_a*eccpa_gl8_x[k];
    if (!(s > 1.0)) return 1;
    sum += eccpa_gl8_w[k]*eccpa_cons_Hr(s, J, dyn);
  }
  *deltaH      = sum*r_a;
  *ddeltaH_drc = eccpa_cons_Hr(r_c + r_a, J, dyn) - eccpa_cons_Hr(r_c - r_a, J, dyn);
  return 0;
}

/* Conservative orbit, set by (pphi, r_a) */

/** Direction w of periastron, z = r_a exp(i w); 0 at z = 0, where it is undefined and nothing
 *  depends on it. */
static double eccpa_peri_angle(double alpha, double beta)
{
  return (alpha*alpha + beta*beta > 0.0) ? atan2(beta, alpha) : 0.0;
}

/**
 * Function: eob_eccpa_nodecount
 * -----------------------------
 *   Nodes for a one-period average to reach ECCPA_AVG_TOL, Eq. (C4).
 *   For a Fourier series (series = 1) only N/2 - 1 harmonics are resolved.
 *
 *   @param[in] e_geo  : r_a/r_c
 *   @param[in] series : 0 for an average, 1 for a Fourier series
 *   @return  number of nodes
 */
int eob_eccpa_nodecount(double e_geo, int series)
{
  if (!(e_geo > 0.0)) return series ? ECCPA_NODES_MIN : 1;
  if (e_geo >= 1.0)   return ECCPA_NODES_MAX;
  const double need = log(ECCPA_AVG_TOL/30.0)/log(e_geo);
  const double target = series ? 2.0*need + 2.0 : need;
  if (target > ECCPA_NODES_MAX) return ECCPA_NODES_MAX;   /* before doubling, which overflows as e -> 1 */
  int N = 1;
  while (N < target) N *= 2;
  if (N < ECCPA_NODES_MIN) N = ECCPA_NODES_MIN;
  if (N > ECCPA_NODES_MAX) N = ECCPA_NODES_MAX;
  return N;
}

/** eob_eccpa_nodecount capped at pa->nmax (stage uses a lower cap than initial data). */
static int eccpa_nnodes(const PostAdDynamics *pa, double e_geo, int series)
{
  const int N = eob_eccpa_nodecount(e_geo, series);
  return (N > pa->nmax) ? pa->nmax : N;
}

/**
 * Function: eob_eccpa_nper
 * ------------------------
 *   Samples per orbital period for reconstruction, Eq. (C9).
 *
 *   @param[in] e_geo : r_a/r_c, largest evolution reached
 *   @return  samples per orbital period
 */
int eob_eccpa_nper(double e_geo)
{
  if (!(e_geo > 0.0))  return ECCPA_NPER_MIN;
  if (e_geo >= 1.0)    return 1000000;
  const double need = log(ECCPA_NPER_TOL)/log(e_geo);
  if (need > 1000000) return 1000000;   /* before cast, which overflows as e -> 1 */
  int n = (int) ceil(need);
  if (n < ECCPA_NPER_MIN) n = ECCPA_NPER_MIN;
  if (n > 1000000) n = 1000000;
  return n;
}

/**
 * Function: eccpa_cons_circular
 * -----------------------------
 *   Circular radius r_circ at this pphi (minimum of H(r, 0, pphi)) and H_rr, H_rrr there.
 *   H_rJ only when r_circ is solved for; if rc_known > 0 r_circ is taken as given.
 *
 *   @return  0 on success, 1 if no minimum was found
 */
static int eccpa_cons_circular(PostAdDynamics *pa, double J, Dynamics *dyn, int oslot,
                               double rc_known)
{

  /* Newtonian guess for r_circ */
  double r = SQ(J);
  if (!(r > 2.0)) r = 10.0;
  double Hrr = 0.0;

  /* Newton on H_r = 0 first, bisection if it fails */
  for (int attempt = 0; attempt < 2; attempt++) {
    int fail = 0;
    if (rc_known > 0.0) r = rc_known;
    else if (attempt == 0) {
      fail = 1;   /* unless Newton converges */
      for (int it = 0; it < 12; it++) {
        double f;
        eob_ecc_ham_rr(r, J, dyn, &f, &Hrr);
        if (!isfinite(f) || !isfinite(Hrr) || Hrr == 0.0) break;
        const double step = f/Hrr;
        r -= step;
        if (!(r > 1.0) || !isfinite(r)) break;
        if (fabs(step) < ECCPA_NEWTON_RTOL*r) { fail = 0; break; }
      }
    } else {
      /* H_r > 0 outside minimum, and < 0 between minimum and peak of the
         potential. Step from seed until H_r changes sign, then bisect. */
      double lo = SQ(J), hi;
      if (!(lo > 2.0)) lo = 10.0;
      hi = lo;
      double flo, fhi, dum;
      eob_ecc_ham_rr(lo, J, dyn, &flo, &dum);
      if (!isfinite(flo)) return 1;
      fhi = flo;
      for (int it = 0; ; it++) {
        if (flo < 0.0 && fhi > 0.0) break;
        /* 80 is a cap, should never be reached under normal circumstances... */
        if (it == 80) return 1;
        if (flo > 0.0) { hi = lo; fhi = flo; lo *= 0.9; eob_ecc_ham_rr(lo, J, dyn, &flo, &dum); }
        else           { lo = hi; flo = fhi; hi *= 1.1; eob_ecc_ham_rr(hi, J, dyn, &fhi, &dum); }
        if (!isfinite(flo) || !isfinite(fhi) || !(lo > 1.0)) return 1;
      }
      // Bisection
      for (int it = 0; it < 200; it++) {
        r = 0.5*(lo + hi);
        double f;
        eob_ecc_ham_rr(r, J, dyn, &f, &dum);
        if (!isfinite(f)) return 1;
        if (f < 0.0) lo = r; else hi = r;
        if (hi - lo < ECCPA_NEWTON_RTOL*r) break;
      }
    }
    if (!fail) {
      pa->cons_orb[oslot][PO_R_CIRC] = r;

      /* derivatives of H at r_circ, by finite differences of H_r rather than of H, because H
         contains rest mass, which would cancel in differences, and H_r doesn't. */
      const double h  = ECCPA_DR_REL*r;
      /* dH/dr at r */
      double f = 0.0;
      eob_ecc_ham(r, 0.0, J, dyn, NULL, &f, NULL, NULL, NULL, NULL, NULL,
                  NULL, NULL, NULL, NULL, NULL, NULL);
      const double fp = eccpa_cons_Hr(r + h, J, dyn);
      const double fm = eccpa_cons_Hr(r - h, J, dyn);
      pa->cons_orb[oslot][PO_HRR]  = (fp - fm)/(2.0*h); // Finite difference for H_rr
      pa->cons_orb[oslot][PO_HRRR] = (fp - 2.0*f + fm)/(h*h); // Finite difference for H_rrr
      /* H_rJ is needed only when r_circ is solved for here */
      if (!(rc_known > 0.0)) {
        const double hJ = ECCPA_DR_REL*fabs(J);
        pa->cons_orb[oslot][PO_HRJ] =
          (eccpa_cons_Hr(r, J + hJ, dyn) - eccpa_cons_Hr(r, J - hJ, dyn))/(2.0*hJ);
      }
      if (pa->cons_orb[oslot][PO_HRR] > 0.0) return 0;   /* a minimum of radial potential */
    }
    if (rc_known > 0.0) return 1;   /* nothing to retry with a known radius */
  }
  return 1;
}

/**
 * Function: eccpa_cons_rc_newton
 * ------------------------------
 *   r_c by Newton on deltaH(r_c) = 0, Eq. (C1).
 *   Wide orbits restart from previous root in slot (reconstruction only, cons_orbit clears it).
 *
 *   @param[in]     r_a : half-width of radial oscillation
 *   @param[in,out] r_c : seed in, root out
 *   @return  1 if r_c is trusted, 0 if bracket has to run
 */
static int eccpa_cons_rc_newton(PostAdDynamics *pa, double r_a, Dynamics *dyn, int oslot, double *r_c)
{
  const double J = pa->cons_orb[oslot][PO_PPHI];
  const double r_circ = pa->cons_orb[oslot][PO_R_CIRC];
  const double e_seed = r_a/(*r_c);   /* r_c is in/out, with seed in and root out */
  const double e4     = e_seed*e_seed*e_seed*e_seed;
  int bad_fast = 0;

  /* Newton from seed, unless seed is good enough or orbit is wide */
  if (1.4*e4 > ECCPA_RC_SEED_TOL && !(r_a > ECCPA_RC_WIDE*r_circ)) {
    const int use_direct = (e_seed > ECCPA_E_DIRECT);
    double step_prev = 0.0;
    int converged = 0;
    for (int it = 0; it < 12; it++) {
      double deltaH, ddeltaH_drc;
      if (eccpa_cons_deltaH(J, *r_c, r_a, use_direct, dyn, &deltaH, &ddeltaH_drc)
          || ddeltaH_drc == 0.0 || !isfinite(deltaH) || !isfinite(ddeltaH_drc)) {
        bad_fast = 1;
        break;
      }
      const double step = deltaH/ddeltaH_drc;
      *r_c -= step;
      if (!(*r_c - r_a > 1.0) || !isfinite(*r_c)) {
        bad_fast = 1;
        break;
      }
      /* Stop on next step we expect, step^2/step_prev, rather than on step just
         taken. */
      if ((step_prev > 0.0 && step*step < ECCPA_NEWTON_RTOL*(*r_c)*step_prev)
          || fabs(step) < ECCPA_NEWTON_RTOL*(*r_c)) { converged = 1; break; }
      step_prev = fabs(step);
    }
    if (!converged) bad_fast = 1;   /* out of iterations: bracket decides */
  }

  /* wide orbit: restart from previous root */
  const double seed = pa->cons_orb[oslot][PO_R_C];
  if (!bad_fast && r_a > ECCPA_RC_WIDE*r_circ && isfinite(seed)
      && seed - r_a > 1.0 && seed > r_circ - r_a && seed < r_circ + r_a) {
    double rs = seed, step_prev = 0.0;
    for (int it = 0; it < 8; it++) {
      double deltaH, ddeltaH_drc;
      if (eccpa_cons_deltaH(J, rs, r_a, 1, dyn, &deltaH, &ddeltaH_drc)
          || !(ddeltaH_drc > 0.0) || !isfinite(deltaH)) break;
      const double step = deltaH/ddeltaH_drc;
      rs -= step;
      if (!(rs - r_a > 1.0) || !isfinite(rs)) break;
      if ((step_prev > 0.0 && step*step < ECCPA_NEWTON_RTOL*rs*step_prev)
          || fabs(step) < ECCPA_NEWTON_RTOL*rs) {
        if (fabs(rs - seed) < 0.01*r_a) { *r_c = rs; return 1; }
        break;
      }
      step_prev = fabs(step);
    }
  }

  /* Check that r_c lies within expected range (r_circ - r_a, r_circ + r_a). 
     if not return an error (and proceed with eccpa_cons_rc_bracket) */
  return !(bad_fast || !isfinite(*r_c) || !(*r_c - r_a > 1.0)
           || *r_c < r_circ - r_a || *r_c > r_circ + r_a
           || r_a > ECCPA_RC_WIDE*r_circ);
}

/**
 * Function: eccpa_cons_rc_bracket
 * -------------------------------
 *   r_c by bisection on deltaH in (r_circ - r_a, r_circ + r_a), with Newton steps when they stay inside.
 *   lower end is moved out if deltaH >= 0 there (large r_a).
 *
 *   @param[in]  r_a : half-width of radial oscillation
 *   @param[out] r_c : root
 *   @return  0 on success, 1 if no bound orbit of this width exists
 */
static int eccpa_cons_rc_bracket(PostAdDynamics *pa, double r_a, Dynamics *dyn, int oslot,
                                 double *r_c)
{
  const double J = pa->cons_orb[oslot][PO_PPHI];
  const double r_circ = pa->cons_orb[oslot][PO_R_CIRC];
  double lo = r_circ - r_a, hi = r_circ + r_a;
  if (lo < r_a + 2.0) lo = r_a + 2.0;
  double deltaH_lo = 0.0, deltaH_hi = 0.0, Hp = 0.0, Hm = 0.0, H_at_lo_out = 0.0;
  int lo_ok = 0;
  for (int it = 0; it <= 6; it++) {
    if (!(hi > lo)) return 1;            /* no room for an orbit this wide */
    eob_ecc_ham(lo + r_a, 0.0, J, dyn, &Hp, NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
    H_at_lo_out = Hp;                    /* = H(r_circ) if lo was not moved */
    eob_ecc_ham(lo - r_a, 0.0, J, dyn, &Hm, NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
    deltaH_lo = Hp - Hm;
    if (isfinite(deltaH_lo) && deltaH_lo < 0.0) { lo_ok = 1; break; }
    lo = r_a + 2.0*(lo - r_a);           /* inner turning point twice as far out */
  }
  if (!lo_ok) return 1;
  const int clamped = (r_circ - r_a < lo);   /* lower end was moved */
  eob_ecc_ham(hi + r_a, 0.0, J, dyn, &Hp, NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
  /* H(r_circ) already computed unless lo moved */
  if (clamped)
    eob_ecc_ham(hi - r_a, 0.0, J, dyn, &Hm, NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
  else
    Hm = H_at_lo_out;         /* same radius, carried over rather than recomputed */
  deltaH_hi = Hp - Hm;
  if (!isfinite(deltaH_hi) || !(deltaH_hi > 0.0)) return 1;
  *r_c = 0.5*(lo + hi);
  for (int it = 0; it < 80; it++) {
    /* r_c - r_a doesn't go below 2, so only a non-finite deltaH fails here */
    double deltaH, ddeltaH_drc;
    if (eccpa_cons_deltaH(J, *r_c, r_a, 1, dyn, &deltaH, &ddeltaH_drc) || !isfinite(deltaH))
      return 1;
    if (deltaH > 0.0) hi = *r_c; else lo = *r_c;
    double r_new = (ddeltaH_drc != 0.0 && isfinite(ddeltaH_drc)) ? *r_c - deltaH/ddeltaH_drc
                                                               : 0.5*(lo + hi);
    /* Newton inside bracket, bisection outside it */
    if (!isfinite(r_new) || !(r_new > lo) || !(r_new < hi)) r_new = 0.5*(lo + hi);
    const double move = fabs(r_new - *r_c);
    *r_c = r_new;
    if (move < ECCPA_NEWTON_RTOL*(*r_c)) break;
  }
  return (!isfinite(*r_c) || !(*r_c - r_a > 1.0)) ? 1 : 0;
}

/**
 * Function: eccpa_cons_rc
 * -----------------------
 *   r_c, e_geo, E and H_r at periastron (PO_HR_M) of orbit (pphi, r_a), Eq. (C1).
 *   Newton from expansion r_c = r_circ - H_rrr r_a^2/(6 H_rr), bracket if Newton fails.
 *
 *   @return  0 on success, 1 if no bound orbit of this width was found
 */
static int eccpa_cons_rc(PostAdDynamics *pa, double r_a, Dynamics *dyn, int oslot)
{
  const double J = pa->cons_orb[oslot][PO_PPHI];

  /* seed. At r_a = 0 it is r_circ, two turning points coincide and there is nothing
     to solve, so we go straight to quantities below. */
  const double *orb = pa->cons_orb[oslot];
  double r_c = orb[PO_R_CIRC] - orb[PO_HRRR]*r_a*r_a/(6.0*orb[PO_HRR]);

  if (r_a > 0.0
      && !eccpa_cons_rc_newton(pa, r_a, dyn, oslot, &r_c)
      && eccpa_cons_rc_bracket(pa, r_a, dyn, oslot, &r_c)) return 1;

  pa->cons_orb[oslot][PO_R_C]   = r_c;
  pa->cons_orb[oslot][PO_E_GEO] = r_a/r_c;
  pa->cons_orb[oslot][PO_GAP_OK] = 0.0;

  /* E = H at periastron, and H_r there, from one call */
  double *slot = pa->cons_orb[oslot];
  eob_ecc_ham(r_c - r_a, 0.0, J, dyn, &slot[PO_E], &slot[PO_HR_M], NULL, NULL, NULL, NULL, NULL,
              NULL, NULL, NULL, NULL, NULL, NULL);
  if (!isfinite(slot[PO_E])) return 1;
  return 0;
}

/**
 * Function: eob_eccpa_cons_orbit
 * ------------------------------
 *   Fill an orbit slot from (pphi, r_a). Returns immediately if slot already holds this orbit.
 *
 *   @param[in] J     : pphi
 *   @param[in] r_a   : |z|
 *   @param[in] dyn   : Dynamics structure
 *   @param[in] oslot : orbit slot to fill
 *   @return  0 on success
 */
int eob_eccpa_cons_orbit(PostAdDynamics *pa, double J, double r_a, Dynamics *dyn, int oslot)
{
  double *orb = pa->cons_orb[oslot];
  if (orb[PO_FULL] == 1.0 && orb[PO_PPHI] == J && orb[PO_R_A] == r_a) return 0; // cache hit
  memset(orb, 0, sizeof pa->cons_orb[oslot]);
  orb[PO_PPHI] = J;
  orb[PO_R_A]  = r_a;
  if (eccpa_cons_circular(pa, J, dyn, oslot, 0.0)) return 1;
  if (eccpa_cons_rc(pa, r_a, dyn, oslot)) return 1;

  /* c_ra from differentiating deltaH(r_c) = 0; below ECCPA_E_TAYLOR from seed, -H_rrr/(3 H_rr) */
  if (orb[PO_E_GEO] > ECCPA_E_TAYLOR) {
    const double hp = eccpa_cons_Hr(orb[PO_R_C] + r_a, J, dyn), hm = orb[PO_HR_M];
    if (hp == hm) return 1;
    orb[PO_C_RA] = -(hp + hm)/((hp - hm)*r_a);
  } else {
    orb[PO_C_RA] = -orb[PO_HRRR]/(3.0*orb[PO_HRR]);
  }
  orb[PO_FULL] = 1.0;
  return 0;
}

/**
 * Function: eob_eccpa_cons_drc_dJ
 * -------------------------------
 *   dr_c/dpphi at fixed r_a, -[H_pphi(r_+) - H_pphi(r_-)]/[H_r(r_+) - H_r(r_-)].
 *   Below ECCPA_E_TAYLOR e -> 0 limit -H_rJ/H_rr.
 *
 *   @return  dr_c/dpphi, or 0 if orbit is degenerate
 */
double eob_eccpa_cons_drc_dJ(PostAdDynamics *pa, int oslot, Dynamics *dyn)
{
  if (!(pa->cons_orb[oslot][PO_E_GEO] > ECCPA_E_TAYLOR))
    return (pa->cons_orb[oslot][PO_HRR] != 0.0)
         ? -pa->cons_orb[oslot][PO_HRJ]/pa->cons_orb[oslot][PO_HRR] : 0.0;

  const double J   = pa->cons_orb[oslot][PO_PPHI];
  const double r_a = pa->cons_orb[oslot][PO_R_A], r_c = pa->cons_orb[oslot][PO_R_C];
  double hp, hm, Jp, Jm;
  eob_ecc_ham(r_c + r_a, 0.0, J, dyn, NULL, &hp, NULL, &Jp, NULL, NULL, NULL,
              NULL, NULL, NULL, NULL, NULL, NULL);
  eob_ecc_ham(r_c - r_a, 0.0, J, dyn, NULL, &hm, NULL, &Jm, NULL, NULL, NULL,
              NULL, NULL, NULL, NULL, NULL, NULL);
  if (hp == hm) return 0.0;
  return -(Jp - Jm)/(hp - hm);
}

/** dr_c/dpphi, dr_c/dr_a and dE/dr_a, dE/dpphi, with E = H at periastron (for eccpa_osculate) */
static void eccpa_cons_dE(PostAdDynamics *pa, int oslot, double r_a,
                          Dynamics *dyn, double *c_ra, double *drc_dJ, double *drc_dra,
                          double *dE_dra, double *dE_dJ)
{
  *c_ra    = pa->cons_orb[oslot][PO_C_RA];
  *drc_dJ  = eob_eccpa_cons_drc_dJ(pa, oslot, dyn);
  *drc_dra = (*c_ra)*r_a;
  const double *orb = pa->cons_orb[oslot];
  double Hr = 0.0, Om = 0.0;
  eob_ecc_ham(orb[PO_R_C] - r_a, 0.0, orb[PO_PPHI], dyn, NULL, &Hr, NULL, &Om, NULL, NULL,
              NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  *dE_dra = Hr*((*drc_dra) - 1.0);
  *dE_dJ  = Hr*(*drc_dJ) + Om;
}

/**
 * Function: eccpa_cons_gap_series
 * -------------------------------
 *   gap(chi) = E - H(r(chi), 0, pphi), Eq. (C2), as a Fourier series, from primitive of
 *   d(gap)/d(chi) = -H_r r_a sin(chi); gap = 0 at both turning points.
 */
static void eccpa_cons_gap_series(PostAdDynamics *pa, int oslot, Dynamics *dyn)
{
  double *orb = pa->cons_orb[oslot];
  EccSeries *gs = &pa->cons_gap[oslot];
  const double r_a = orb[PO_R_A];
  int Ng = 2*eccpa_nnodes(pa, orb[PO_E_GEO], 0);
  if (Ng < ECCPA_NODES_MIN_GAP) Ng = ECCPA_NODES_MIN_GAP;
  if (Ng > ECCPA_NODES_DIM) Ng = ECCPA_NODES_DIM;

  /* d(gap)/d(chi) = -H_r dr/dchi on nodes */
  double f[ECCPA_NODES_DIM];
  for (int j = 0; j < Ng; j++) {
    const double c = 2.0*Pi*(j + 0.5)/Ng;
    f[j] = -eccpa_cons_Hr(orb[PO_R_C] - r_a*cos(c), orb[PO_PPHI], dyn)*r_a*sin(c);
  }
  
  /* Fourier series, and its primitive */
  EccSeries fs;
  eob_ecc_series_fit(f, Ng, &fs);
  eob_ecc_series_prim(&fs, gs);
  
  /* constant of integration makes gap(0) = 0 at periastron */
  double gap0 = -eob_ecc_series_eval(gs, 0.0);

  /* fix apastron residual with (g/2)(1 - cos chi): exact at both ends, no new harmonics */
  const double g = eob_ecc_series_eval(gs, Pi) + gap0;
  if (gs->M >= 1) {
    gs->a[1] += 0.5*g;
    gap0     -= 0.5*g;
  }
  orb[PO_GAP0]   = gap0;
  orb[PO_GAP_OK] = 1.0;
}

/**
 * Function: eccpa_cons_state_at_chi
 * ---------------------------------
 *   State at anomaly chi (from periastron), Eqs. (C2)-(C3); written to pa except pa->X.
 *   Uses psi = gap/(dr/dX)^2, finite at turning points and at r_a = 0.
 */
static void eccpa_cons_state_at_chi(PostAdDynamics *pa, int oslot, double chi, int series,
                                    int with_flux, Dynamics *dyn)
{
  const double J = pa->cons_orb[oslot][PO_PPHI], r_a = pa->cons_orb[oslot][PO_R_A];
  const double t1 = 0.5 - 0.5/sqrt(3.0), t2 = 0.5 + 0.5/sqrt(3.0); /* 2-point Gauss on [0,1] */

  pa->r    = pa->cons_orb[oslot][PO_R_C] - r_a*cos(chi);
  pa->drdX = r_a*sin(chi); /* dr/dchi = dr/dX since X = w + chi at fixed w */
  const double r = pa->r;

  /* psi. At r_a = 0 we use its limit. Elsewhere nodes keep sin(chi) >= sin(pi/N), so we
     don't divide by zero. */
  double psi;
  double H0 = 0.0, d2H0_zero = 0.0;
  eob_ecc_ham(r, 0.0, J, dyn, &H0, NULL, NULL, NULL, NULL, &d2H0_zero, NULL,
              NULL, NULL, NULL, NULL, NULL, NULL);

  if (r_a > 0.0) {
    double gap;
    if (!series && pa->cons_orb[oslot][PO_E_GEO] > ECCPA_E_DIRECT) {
      gap = pa->cons_orb[oslot][PO_E] - H0;
    } else {
      if (pa->cons_orb[oslot][PO_GAP_OK] != 1.0) eccpa_cons_gap_series(pa, oslot, dyn);
      gap = eob_ecc_series_eval(&pa->cons_gap[oslot], chi) + pa->cons_orb[oslot][PO_GAP0];
    }
    /* Near a turning point both gap and (dr/dX)^2 go to zero, 
    so we switch on sin^2(chi), not on denominator being zero. */
    const double s2 = (pa->drdX*pa->drdX)/(r_a*r_a);
    if (s2 > ECCPA_TURNING_POINT_TOL && gap > 0.0) {
      psi = gap/(pa->drdX*pa->drdX);
    } else {
      /* At a turning point r_t, gap ~ -H_r(r_t)(r - r_t) and (dr/dX)^2 ~ 2 r_a |r - r_t|, so
         psi -> |H_r(r_t)|/(2 r_a). Also where E - H rounds to zero or below next to one. */
      psi = fabs(eccpa_cons_Hr(r, J, dyn))/(2.0*r_a);
    }
  } else {
    /* at r_a = 0 (circular orbit) psi = H_rr/2 */
    psi = 0.5*pa->cons_orb[oslot][PO_HRR];
  }
  if (!(psi > 0.0) || !isfinite(psi)) psi = 0.5*pa->cons_orb[oslot][PO_HRR];

  /* Fixed-point iteration for Wstar (App. C), starting from Phi at prstar = 0, which is
     dH/d(prstar^2) there, half d2H/dprstar2. */
  double Phi = 0.5*d2H0_zero;
  /* dH/d(prstar^2) at prstar = 0 is positive for EOB Hamiltonian.*/
  if (!(Phi > 0.0)) Phi = 1.0;
  double Wstar = sqrt(psi/Phi);
  double move_prev = 0.0;
  for (int it = 0; it < 6; it++) {
    const double prs_abs = fabs(pa->drdX)*Wstar;
    if (!(prs_abs > 0.0)) break;
    const double prs2 = prs_abs*prs_abs;
    /* Phi of Eq. (C2): finite diff above ECCPA_PRSTAR_DIRECT, 2-pt Gauss below */
    double Phi_new;
    if (prs_abs > ECCPA_PRSTAR_DIRECT) {
      double H_prs;
      eob_ecc_ham(r, prs_abs, J, dyn, &H_prs, NULL, NULL, NULL, NULL, NULL, NULL,
                  NULL, NULL, NULL, NULL, NULL, NULL);
      Phi_new = (H_prs - H0)/prs2;
    } else {
      Phi_new = 0.5*(eccpa_cons_dH_dprs2(r, sqrt(t1*prs2), J, dyn)
                   + eccpa_cons_dH_dprs2(r, sqrt(t2*prs2), J, dyn));
    }
    if (!(Phi_new > 0.0)) break;
    const double W_new = sqrt(psi/Phi_new);
    const double move  = fabs(W_new - Wstar);
    Wstar = W_new;
    /* Stop on next move we expect, move^2/move_prev (iteration shrinks move by
       about same factor each pass). ECCPA_WSTAR_TOL is accuracy Wstar needs. */
    if (move_prev > 0.0 && move*move < ECCPA_WSTAR_TOL*Wstar*move_prev) break;
    if (move < ECCPA_NEWTON_RTOL*Wstar) break;
    move_prev = move;
  }
  pa->Wstar  = Wstar;
  pa->prstar = pa->drdX*Wstar;

  /* Compute various quantities at current phase-space point */
  double H, dH_dr, dH_dprstar, dH_dpphi, sqrtAbyB, d2H0, Heff, rc, drc, d2rc, A, dA, Heff_orb;
  eob_ecc_ham(r, pa->prstar, J, dyn, &H, &dH_dr, &dH_dprstar, &dH_dpphi,
              &sqrtAbyB, &d2H0, &Heff, &rc, &drc, &d2rc, &A, &dA, &Heff_orb);
  pa->H = H;  pa->Heff_orb = Heff_orb;  pa->rc = rc;  pa->A = A;
  /* dH/d(prstar^2) for dt/dX; at prstar = 0 it is half d2H/dprstar2 */
  const double dH_dprs2 = (fabs(pa->prstar) > 0.0) ? dH_dprstar/(2.0*pa->prstar) : 0.5*d2H0;

  pa->dH_dprs2     = dH_dprs2;
  pa->sqrtAbyB = sqrtAbyB;
  pa->rdot  = sqrtAbyB*dH_dprstar;
  pa->Omg   = dH_dpphi;
  pa->ddotr = -sqrtAbyB*dH_dr*sqrtAbyB*d2H0;
  pa->dt_dX = 1.0/(2.0*sqrtAbyB*pa->Wstar*dH_dprs2);
  if (!isfinite(pa->dt_dX) || pa->dt_dX <= 0.0) pa->dt_dX = 0.0;

  if (!with_flux) {
    pa->Fr = pa->FE = pa->F[C_PPHI] = 0.0;
    return;
  }
  eob_ecc_flux_at(r, pa->prstar, J, pa->rdot, pa->Omg, pa->ddotr,
                  H, Heff, sqrtAbyB, dH_dr, dyn, rc, drc, d2rc, A, dA, &pa->F[C_PPHI], &pa->Fr);
  pa->FE = pa->rdot*pa->Fr + pa->Omg*pa->F[C_PPHI];
}

/**
 * Function: eccpa_osculate
 * ------------------------
 *   Rates (Fa, Fb) of (alpha, beta) from Eq. (41) at node j, or at point in pa if j < 0.
 *   Needs dWstar/dr_a and dWstar/dpphi at fixed r, from H = E.
 *
 *   @return  0 on success, 1 if 2x2 system is singular
 */
static int eccpa_osculate(PostAdDynamics *pa, int oslot,
                          double alpha, double beta,
                          int j, double c_ra, double drc_dJ, double drc_dra,
                          double dE_dra, double dE_dJ, double *Fa, double *Fb)
{
  /* node j, or point left in pa when j < 0 */
  const double nd_Wstar = (j >= 0) ? pa->node[PA_WSTAR][j] : pa->Wstar;
  const double nd_drdX = (j >= 0) ? pa->node[PA_DRDX][j] : pa->drdX;
  const double nd_r = (j >= 0) ? pa->node[PA_R][j] : pa->r;
  const double nd_dH_dprs2 = (j >= 0) ? pa->node[PA_DH_DPRS2][j] : pa->dH_dprs2;
  const double nd_Omg = (j >= 0) ? pa->node[PA_OMG][j] : pa->Omg;
  const double nd_X = (j >= 0) ? pa->node[PA_X][j] : pa->X;
  const double nd_FJ = (j >= 0) ? pa->node[PA_F + C_PPHI][j] : pa->F[C_PPHI];
  const double nd_sqrtAbyB = (j >= 0) ? pa->node[PA_SQRTAB][j] : pa->sqrtAbyB;
  const double nd_Fr = (j >= 0) ? pa->node[PA_FR][j] : pa->Fr;
  const double r_a  = pa->cons_orb[oslot][PO_R_A];
  const double Ws   = nd_Wstar;
  const double u2   = nd_drdX*nd_drdX;
  const double rcmr = pa->cons_orb[oslot][PO_R_C] - nd_r;
  double dWs_dJ = 0.0, dWs_dra = 0.0;
  if (Ws > 0.0 && u2 > 0.0 && nd_dH_dprs2 > 0.0) {
    const double du2_dra = 2.0*r_a - 2.0*rcmr*drc_dra;
    const double du2_dJ  =         - 2.0*rcmr*drc_dJ;
    const double den     = 2.0*Ws*u2;
    dWs_dra = (dE_dra/nd_dH_dprs2            - Ws*Ws*du2_dra)/den;
    dWs_dJ  = ((dE_dJ - nd_Omg)/nd_dH_dprs2 - Ws*Ws*du2_dJ)/den;
  }
  /* q is eta of Eq. (42), (dr/dX)(1/r_a)(dWstar/dr_a)/Wstar = sin(chi) (dWstar/dr_a)/Wstar.
     At r_a = 0 it multiplies alpha and beta, which are zero, so we set it to 0. */
  const double sinchi = (r_a > 0.0) ? nd_drdX/r_a : 0.0;
  const double q = (Ws != 0.0 && r_a > 0.0) ? sinchi*dWs_dra/Ws : 0.0;

  const double c = cos(nd_X), s = sin(nd_X);
  const double a11 = c - c_ra*alpha,  a12 = s - c_ra*beta;
  const double a21 = s + q*alpha,     a22 = q*beta - c;
  const double R1  = drc_dJ*nd_FJ;
  const double R2  = (Ws != 0.0)
                   ? (nd_sqrtAbyB*nd_Fr - nd_drdX*dWs_dJ*nd_FJ)/Ws : 0.0;
  const double D   = a11*a22 - a12*a21;
  if (!(fabs(D) > 0.0)) return 1;
  *Fa = (R1*a22 - a12*R2)/D;
  *Fb = (a11*R2 - R1*a21)/D;
  return 0;
}

/**
 * Function: eob_eccpa_cons_sample_nodes
 * -------------------------------------
 *   Sample orbit on N nodes of Eq. (C4) into pa->node, with flux and rates if with_flux.
 *
 *   @param[in,out] pa          : struct
 *   @param[in]     oslot       : slot holding orbit
 *   @param[in]     alpha, beta : eccentricity vector
 *   @param[in]     N           : number of nodes
 *   @param[in]     with_flux   : evaluate flux and rates
 *   @param[in]     dyn         : Dynamics structure
 *   @param[out]    out_dt, out_Omg : if not NULL, only dt/dX and Omega, there
 *   @return  0 on success, -1 on failure
 */
int eob_eccpa_cons_sample_nodes(PostAdDynamics *pa, int oslot, double alpha, double beta,
                                int N, int with_flux, Dynamics *dyn,
                                double *out_dt, double *out_Omg)
{
  if (N < 1 || N > ECCPA_NODES_DIM) return -1;
  const double r_a = pa->cons_orb[oslot][PO_R_A];
  const double w   = eccpa_peri_angle(alpha, beta);
  const int rates = with_flux && !out_dt;
  double c_ra = 0.0, drc_dJ = 0.0, drc_dra = 0.0, dE_dra = 0.0, dE_dJ = 0.0;
  if (rates) eccpa_cons_dE(pa, oslot, r_a, dyn, &c_ra, &drc_dJ, &drc_dra, &dE_dra, &dE_dJ);

  for (int j = 0; j < N; j++) {
    const double chi = 2.0*Pi*(j + 0.5)/N;
    eccpa_cons_state_at_chi(pa, oslot, chi, 1, with_flux, dyn);
    pa->X = chi + w;
    if (out_dt) {
      /* only dt/dX and Omega are wanted, and pa->node[][] must not be overwritten */
      out_dt[j]  = pa->dt_dX;
      out_Omg[j] = pa->Omg;
      continue;
    }
    pa->node[PA_X][j]       = pa->X;
    pa->node[PA_R][j]       = pa->r;
    pa->node[PA_DRDX][j]    = pa->drdX;
    pa->node[PA_PRSTAR][j]  = pa->prstar;
    pa->node[PA_WSTAR][j]   = pa->Wstar;
    pa->node[PA_DTDX][j]    = pa->dt_dX;
    pa->node[PA_DH_DPRS2][j]    = pa->dH_dprs2;
    pa->node[PA_SQRTAB][j]  = pa->sqrtAbyB;
    pa->node[PA_OMG][j]     = pa->Omg;
    pa->node[PA_FR][j]      = pa->Fr;
    pa->node[PA_FE][j]    = pa->FE;
    pa->node[PA_F + C_PPHI][j]    = pa->F[C_PPHI];
    if (rates && eccpa_osculate(pa, oslot, alpha, beta, j, c_ra, drc_dJ, drc_dra, dE_dra, dE_dJ,
                                &pa->node[PA_F + C_ALPHA][j], &pa->node[PA_F + C_BETA][j]))
      return -1;
  }
  if (!out_dt) pa->nnodes = N;
  return 0;
}

/**
 * Function: eob_eccpa_cons_sample_point
 * -------------------------------------
 *   Like eob_eccpa_cons_sample_nodes, at one anomaly X, results left in pa.
 *
 *   @param[in,out] pa          : struct
 *   @param[in]     oslot       : slot holding orbit
 *   @param[in]     alpha, beta : eccentricity vector
 *   @param[in]     X           : anomaly
 *   @param[in]     with_flux   : evaluate flux and rates
 *   @param[in]     dyn         : Dynamics structure
 *   @return  0 on success, 1 if eccpa_osculate fails
 */
int eob_eccpa_cons_sample_point(PostAdDynamics *pa, int oslot, double alpha, double beta,
                                double X, int with_flux, Dynamics *dyn)
{
  const double w = eccpa_peri_angle(alpha, beta);
  eccpa_cons_state_at_chi(pa, oslot, X - w, 0, with_flux, dyn);
  pa->X = X;
  if (!with_flux) return 0;
  double c_ra, drc_dJ, drc_dra, dE_dra, dE_dJ;
  eccpa_cons_dE(pa, oslot, pa->cons_orb[oslot][PO_R_A], dyn, &c_ra, &drc_dJ, &drc_dra, &dE_dra, &dE_dJ);
  return eccpa_osculate(pa, oslot, alpha, beta, -1, c_ra, drc_dJ, drc_dra, dE_dra, dE_dJ,
                        &pa->F[C_ALPHA], &pa->F[C_BETA]);
}

static double eccpa_tp_f(double r, void *params)
{
  eccpa_tp_par *p = (eccpa_tp_par *) params;
  double H;
  eob_ecc_ham(r, 0.0, p->pphi, p->dyn, &H, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  return H - p->E;
}

/**
 * Function: eob_eccpa_cons_turning_points
 * ---------------------------------------
 *   Turning points of orbit (E, pphi), roots of H(r, 0, pphi) = E.
 *
 *   @param[in,out] pa      : orbit written to pa->cons_orb[slot]
 *   @param[in]     E, pphi : constants
 *   @param[in]     rseed   : a radius between turning points
 *   @param[in]     dyn     : Dynamics structure
 *   @param[in]     slot    : orbit slot
 *   @return  0 on success, ECCPA_UNBOUND or ECCPA_NOTURN otherwise
 */
int eob_eccpa_cons_turning_points(PostAdDynamics *pa, double E, double pphi, double rseed,
                                  Dynamics *dyn, int slot)
{
  PROF_START(PROF_ECCPA_ORBIT);
  int status = 0;   /* one exit, so timer is always stopped */
  pa->cons_orb[slot][PO_FULL] = pa->cons_orb[slot][PO_GAP_OK] = 0.0;

  eccpa_tp_par par = { pphi, E, dyn };

  /* seed must be strictly between turning points (H < E there). If it sits on a
     turning point, brackets below don't give right roots. */
  double f_seed = eccpa_tp_f(rseed, &par);
  if (!isfinite(f_seed)) { status = ECCPA_NOTURN; goto done; }
  if (f_seed >= -1e-12 * fabs(E)) {
    /* seed is not inside, so use r_circ instead: minimum of H(r, 0, pphi), which is
       between turning points if a bound orbit with this E exists. */
    if (eccpa_cons_circular(pa, pphi, dyn, slot, 0.0)) {
      status = ECCPA_NOTURN;
      goto done;
    }
    const double r_circ = pa->cons_orb[slot][PO_R_CIRC];

    /* H at circular radius is minimum: above E means no bound orbit here. */
    if (eccpa_tp_f(r_circ, &par) >= 0.0) {
      status = ECCPA_NOTURN;
      goto done;
    }
    rseed = r_circ;
  }

  /* Bracket two turning points by stepping out from seed for apastron and in for
     periastron, stopping at first radius where H(r, 0, pphi) > E. */
  double r_lo = rseed, r_hi = rseed;
  int ok_hi = 0, ok_lo = 0;
  for (int k = 0; k < 60; k++) {
    r_hi *= 1.6;
    const double f = eccpa_tp_f(r_hi, &par);
    if (DEBUG) printf("Iter HI %d: r_hi = %g, f = %g\n", k, r_hi, f);
    if (!isfinite(f)) break;
    if (f > 0.0) {
      ok_hi = 1;
      break;
    }
  }
  for (int k = 0; k < 60; k++) {
    r_lo /= 1.6;
    if (r_lo < 1.0) break;
    const double f = eccpa_tp_f(r_lo, &par);
    if (DEBUG) printf("Iter LO %d: r_lo = %g, f = %g\n", k, r_lo, f);
    if (!isfinite(f)) break;
    if (f > 0.0) {
      ok_lo = 1;
      break;
    }
  }
  if (!ok_hi) { status = ECCPA_UNBOUND; goto done; }   /* no outer turning point: unbound */
  if (!ok_lo) { status = ECCPA_NOTURN; goto done; }    /* no inner turning point: plunging */

  gsl_function F;
  F.function = &eccpa_tp_f;
  F.params   = &par;
  gsl_root_fsolver *s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

  double rmax = 0.0, rmin = 0.0;
  for (int side = 0; side < 2; side++) {
    double a = (side == 0) ? rseed : r_lo;
    double b = (side == 0) ? r_hi  : rseed;
    if (gsl_root_fsolver_set(s, &F, a, b) != GSL_SUCCESS) {
      gsl_root_fsolver_free(s);
      status = ECCPA_NOTURN;
      goto done;
    }
    int gsl_status = GSL_CONTINUE;
    for (int it = 0; it < 100 && gsl_status == GSL_CONTINUE; it++) {
      if (gsl_root_fsolver_iterate(s) != GSL_SUCCESS) break;
      a = gsl_root_fsolver_x_lower(s);
      b = gsl_root_fsolver_x_upper(s);
      gsl_status = gsl_root_test_interval(a, b, 0.0, 1e-13);
    }
    const double root = gsl_root_fsolver_root(s);
    if (side == 0) rmax = root;
    else           rmin = root;
  }
  gsl_root_fsolver_free(s);

  if (!(rmax > rmin)) { status = ECCPA_NOTURN; goto done; }

  /* 2 Newton steps on H - E: Brent's 1e-13 -> roundoff */
  for (int side = 0; side < 2; side++) {
    double *root = (side == 0) ? &rmax : &rmin;
    for (int it = 0; it < 3; it++) {
      double H, dH_dr;
      eob_ecc_ham(*root, 0.0, pphi, dyn, &H, &dH_dr, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
      if (dH_dr == 0.0 || !isfinite(dH_dr)) break;
      const double step = (H - E) / dH_dr;
      if (!isfinite(step)) break;
      *root -= step;
      if (fabs(step) < 1e-15 * fabs(*root)) break;
    }
  }
  if (!(rmax > rmin)) { status = ECCPA_NOTURN; goto done; }

  pa->cons_orb[slot][PO_E]     = E;
  pa->cons_orb[slot][PO_PPHI]  = pphi;
  pa->cons_orb[slot][PO_R_MAX] = rmax;
  pa->cons_orb[slot][PO_R_MIN] = rmin;
  pa->cons_orb[slot][PO_R_C]   = 0.5 * (rmax + rmin);
  pa->cons_orb[slot][PO_R_A]   = 0.5 * (rmax - rmin);
  pa->cons_orb[slot][PO_E_GEO] = pa->cons_orb[slot][PO_R_A] / pa->cons_orb[slot][PO_R_C];

done:
  PROF_STOP(PROF_ECCPA_ORBIT);
  return status;
}

/* Averaged rates and maps M_theta (T and -V); coefficients stored divided by r_a^m */

/**
 * Function: eccpa_kappa
 * ---------------------
 *   kappa_m of dt/dX, Eq. (C6). Below ECCPA_PROBE_RA r_c taken on probe orbit.
 *
 *   @return  0 on success, 1 if probe orbit could not be built
 */
static int eccpa_kappa(PostAdDynamics *pa, int oslot, const double *base,
                       int N, Dynamics *dyn,
                       double *kappa_re, double *kappa_im)
{
  const double probe = ECCPA_PROBE_RA*pa->cons_orb[oslot][PO_R_C];
  const double ra_use = (pa->cons_orb[oslot][PO_R_A] > probe) ? pa->cons_orb[oslot][PO_R_A] : probe;

  const double *src;
  double probe_dtdX[ECCPA_NODES_DIM], probe_Omg[ECCPA_NODES_DIM];
  if (ra_use == pa->cons_orb[oslot][PO_R_A]) {
    /* dt/dX depends only on chi, so caller's nodes will do, whatever their w */
    src = base;
  } else {
    /* probe orbit, in aux slot; no flux */
    if (eob_eccpa_cons_orbit(pa, pa->cons_orb[oslot][PO_PPHI], ra_use, dyn, PO_AUX)) return 1;
    if (eob_eccpa_cons_sample_nodes(pa, PO_AUX, ra_use, 0.0, N, 0, dyn,
                                    probe_dtdX, probe_Omg) < 0) return 1;
    src = probe_dtdX;
  }

  EccSeries s;
  eob_ecc_series_fit(src, N, &s);

  /* eob_ecc_series_fit returns real series f = a_0 + sum_m (a_m cos + b_m sin), whose
     complex coefficient of exp(i m chi) is (a_m - i b_m)/2 for m >= 1. */
  double p = 1.0;
  kappa_re[0] = s.a[0];
  kappa_im[0] = 0.0;
  for (int m = 1; m <= ECCPA_HARMONICS_MAX; m++) {
    p *= ra_use;
    if (m > s.M) {
      kappa_re[m] = kappa_im[m] = 0.0;
      continue;
    }
    kappa_re[m] =  0.5*s.a[m]/p;
    kappa_im[m] = -0.5*s.b[m]/p;
  }
  return 0;
}

/** Powers r_a^0 to r_a^M, to multiply harmonic m by r_a^m. */
static void eccpa_rapow(double *rap, int M, double r_a)
{
  rap[0] = 1.0;
  for (int m = 1; m <= M; m++) rap[m] = rap[m-1]*r_a;
}

/** c[m] = cos(m chi), s[m] = sin(m chi) for m = 0..M, by angle-addition recurrence: two
 *  library calls instead of 2M. drift over M steps is a few epsilon. */
static void eccpa_trig(double chi, int M, double *c, double *s)
{
  const double c1 = cos(chi), s1 = sin(chi);
  c[0] = 1.0; s[0] = 0.0;
  for (int m = 1; m <= M; m++) {
    c[m] = c[m-1]*c1 - s[m-1]*s1;
    s[m] = s[m-1]*c1 + c[m-1]*s1;
  }
}

/** sum_m r_a^m Re[kappa_m exp(i m chi)], and its chi derivative if dsum != NULL */
static double eccpa_ksum(const double *kappa_re, const double *kappa_im, int M,
                         double r_a, double chi, double *dsum)
{
  double t = 0.0, d = 0.0, p = 1.0, c[ECCPA_HARMONICS_MAX+1], s[ECCPA_HARMONICS_MAX+1];
  eccpa_trig(chi, M, c, s);
  for (int m = 1; m <= M; m++) {
    p *= r_a;
    t += p*(kappa_re[m]*c[m] - kappa_im[m]*s[m]);
    if (dsum) d += p*m*(-kappa_re[m]*s[m] - kappa_im[m]*c[m]);
  }
  if (dsum) *dsum = d;
  return t;
}

/** <zdot>, Eq. (C6), from node rates Fa, Fb and kappa_m of dt/dX */
static void eccpa_convolve(const double *Fa, const double *Fb, const double *X, int N,
                           const double *kappa_re, const double *kappa_im,
                           double alpha, double beta, double n,
                           double *Fbar_alpha, double *Fbar_beta)
{
  const double r_a = sqrt(alpha*alpha + beta*beta);
  double Cre[ECCPA_NODES_DIM], Cim[ECCPA_NODES_DIM];
  for (int j = 0; j < N; j++) {
    const double c = cos(X[j]), s = sin(X[j]);
    Cre[j] =  Fa[j]*c + Fb[j]*s;
    Cim[j] = -Fa[j]*s + Fb[j]*c;
  }
  EccSeries Cs_re, Cs_im;
  eob_ecc_series_fit(Cre, N, &Cs_re);
  eob_ecc_series_fit(Cim, N, &Cs_im);

  /* r_a^m for hat tau_m = kappa_|m| r_a^|m| */
  double rap[ECCPA_HARMONICS_MAX+2];
  eccpa_rapow(rap, ECCPA_HARMONICS_MAX+1, r_a);

  /* hat C_m, m in [-M, M], from two real fits */
  double totr = 0.0, toti = 0.0;
  for (int m = -ECCPA_HARMONICS_MAX; m <= ECCPA_HARMONICS_MAX; m++) {
    const int am = (m < 0) ? -m : m;
    if (am > Cs_re.M) continue;
    double cr, ci;
    if (m == 0) {
      cr = Cs_re.a[0];
      ci = Cs_im.a[0];
    }
    else if (m > 0)  { cr = 0.5*(Cs_re.a[am] + Cs_im.b[am]);
                       ci = 0.5*(Cs_im.a[am] - Cs_re.b[am]); }
    else             { cr = 0.5*(Cs_re.a[am] - Cs_im.b[am]);
                       ci = 0.5*(Cs_im.a[am] + Cs_re.b[am]); }

    /* kappa at harmonic -(1 + m); dt/dX is real, so a negative harmonic is conjugate */
    const int k = -(1 + m);
    const int ak = (k < 0) ? -k : k;
    if (ak > ECCPA_HARMONICS_MAX) continue;
    const double wr = kappa_re[ak];
    const double wi = (k >= 0) ? kappa_im[ak] : -kappa_im[ak];

    /* r_a^(|1 + m| - 1), because one r_a of hat tau goes into z below. At m = -1 this is 1/r_a,
       which is fine because hat C_{-1} = O(r_a). */
    double fac;
    if (m == -1) {
      if (!(r_a > 0.0)) continue;
      fac = 1.0/r_a;
    } else {
      fac = rap[ak - 1];
    }
    totr += (cr*wr - ci*wi)*fac;
    toti += (cr*wi + ci*wr)*fac;
  }
  /* n r_a exp(i w) = n z times sum */
  *Fbar_alpha = n*(alpha*totr - beta*toti);
  *Fbar_beta  = n*(alpha*toti + beta*totr);
}

/**
 * Function: eob_eccpa_rates
 * -------------------------
 *   First-order averaged rates Fbar^(1) of (pphi, alpha, beta) and E, Eq. (B1), with n, k, omega
 *   and adiabaticity, Eq. (C8). z average is convolution of Eq. (C6).
 *
 *   @param[in]     J, alpha, beta : state
 *   @param[in]     Nmin : at least this many nodes; 0 = eob_eccpa_nodecount
 *   @param[in]     dyn  : Dynamics structure
 *   @param[in,out] pa   : struct
 *   @return  0 on success
 */
int eob_eccpa_rates(PostAdDynamics *pa, double J, double alpha, double beta, int Nmin,
                    Dynamics *dyn)
{
  pa->FbarE = pa->Fbar[C_PPHI] = pa->Fbar[C_ALPHA] = pa->Fbar[C_BETA] = 0.0;
  pa->dc_dpphi[C_ALPHA] = pa->dc_dpphi[C_BETA] = pa->adiab = 0.0;
  for (int th = 0; th < TH_N; th++) pa->omega[th] = 0.0;
  pa->P_r = 0.0;
  const double r_a = sqrt(alpha*alpha + beta*beta);
  if (eob_eccpa_cons_orbit(pa, J, r_a, dyn, PO_CUR)) return 1;
  int N = eccpa_nnodes(pa, pa->cons_orb[PO_CUR][PO_E_GEO], 0);
  if (Nmin > N) N = (Nmin > ECCPA_NODES_DIM) ? ECCPA_NODES_DIM : Nmin;
  pa->nnodes = N;

  /* Sample orbit on N nodes, with flux and three F_mu. This is only flux
     evaluation here, and everything below reuses node values. */
  if (eob_eccpa_cons_sample_nodes(pa, PO_CUR, alpha, beta, N, 1, dyn, NULL, NULL) < 0) {
    return 1;
  }

  /* sums of Eq. (C4) for P_r, azimuth advance, F_E and F_pphi */
  const double node_spacing = 2.0*Pi/N;
  double P_r = 0.0, dphi = 0.0, sum_FE = 0.0, sum_FJ = 0.0;

  for (int j = 0; j < N; j++) {
    const double dt = pa->node[PA_DTDX][j]*node_spacing;
    P_r  += dt;
    // sum to avg at end
    dphi += pa->node[PA_OMG][j]*dt;
    sum_FE += pa->node[PA_FE][j]*dt;
    sum_FJ += pa->node[PA_F + C_PPHI][j]*dt;
  }

  if (!(P_r > 0.0)) {
    return 1;
  }
  pa->P_r  = P_r;
  pa->omega[TH_LAMBDA_MEAN]    = 2.0*Pi/P_r;
  const double k = dphi/(2.0*Pi) - 1.0;   /* periastron advance */
  pa->omega[TH_AZIMUTH]  = (1.0 + k)*pa->omega[TH_LAMBDA_MEAN];
  pa->FbarE = sum_FE/P_r;
  pa->Fbar[C_PPHI] = sum_FJ/P_r;

  /* convolution. C and dt/dX are on same nodes, so both fits give harmonics in chi, and
     exp(i w) is put back once at end. */
  double kappa_re[ECCPA_HARMONICS_MAX+1], kappa_im[ECCPA_HARMONICS_MAX+1];
  if (eccpa_kappa(pa, PO_CUR, pa->node[PA_DTDX], N, dyn, kappa_re, kappa_im)) {
    return 1;
  }
  eccpa_convolve(pa->node[PA_F + C_ALPHA], pa->node[PA_F + C_BETA], pa->node[PA_X], N, kappa_re, kappa_im,
                 alpha, beta, pa->omega[TH_LAMBDA_MEAN], &pa->Fbar[C_ALPHA], &pa->Fbar[C_BETA]);
  pa->dc_dpphi[C_ALPHA] = (pa->Fbar[C_PPHI] != 0.0) ? pa->Fbar[C_ALPHA]/pa->Fbar[C_PPHI] : 0.0;
  pa->dc_dpphi[C_BETA] = (pa->Fbar[C_PPHI] != 0.0) ? pa->Fbar[C_BETA]/pa->Fbar[C_PPHI] : 0.0;

  const double Ebind = pa->cons_orb[PO_CUR][PO_E] - 1.0/EOBPars->nu;
  pa->adiab = fabs(pa->FbarE)*P_r/fabs(Ebind);
  return 0;
}

/** n = 2 pi/P_r and omega = <Omega> from dt/dX and Omega on N nodes. 2 pi of n cancels
 *  against node spacing 2 pi/N. Returns non-zero if P_r is not positive. */
static int eccpa_mean_freqs(const double *dt, const double *Omg, int N,
                            double *n, double *omg)
{
  double sum_f = 0.0, sum_Of = 0.0;
  for (int j = 0; j < N; j++) {
    sum_f  += dt[j];
    sum_Of += Omg[j]*dt[j];
  }
  if (!(sum_f > 0.0)) return 1;
  *n   = (double) N/sum_f;
  *omg = sum_Of/sum_f;
  return 0;
}

/**
 * Function: eccpa_cons_map_coefs
 * -------------------------------
 *   omega[th] and kappa_m of M_theta = T, -V, Eq. (44), from dT/dX = n dt/dX - 1 and
 *   dV/dX = (Omega - omega) dt/dX. Probe orbit below ECCPA_PROBE_RA.
 *
 *   @param[in]  base_dt, base_Omg : dt/dX and Omega on N nodes
 *   @param[out] map_coefs         : omega[th] and kappa_m
 *   @return  0 on success
 */
static int eccpa_cons_map_coefs(PostAdDynamics *pa, int oslot,
                                const double *base_dt, const double *base_Omg, int N, Dynamics *dyn,
                                EccMapCoefs *map_coefs)
{
  memset(map_coefs, 0, sizeof *map_coefs);

  /* n and omega always from this orbit, never from probe, which would be off by
     O(ECCPA_PROBE_RA^2) */
  if (eccpa_mean_freqs(base_dt, base_Omg, N, &map_coefs->omega[TH_LAMBDA_MEAN],
                       &map_coefs->omega[TH_AZIMUTH])) return 1;

  const double probe  = ECCPA_PROBE_RA*pa->cons_orb[oslot][PO_R_C];
  const double ra_use = (pa->cons_orb[oslot][PO_R_A] > probe) ? pa->cons_orb[oslot][PO_R_A] : probe;

  double f[ECCPA_NODES_DIM], g[ECCPA_NODES_DIM];
  double probe_dtdX[ECCPA_NODES_DIM], probe_Omg[ECCPA_NODES_DIM];
  const double *nd_dt = base_dt, *nd_Omg = base_Omg;
  /* integrands use n and omega of orbit they are sampled on (probe, if used) */
  double n_h = map_coefs->omega[TH_LAMBDA_MEAN], omg_h = map_coefs->omega[TH_AZIMUTH];
  if (ra_use != pa->cons_orb[oslot][PO_R_A]) {
    if (eob_eccpa_cons_orbit(pa, pa->cons_orb[oslot][PO_PPHI], ra_use, dyn, PO_AUX)) return 1;
    if (eob_eccpa_cons_sample_nodes(pa, PO_AUX, ra_use, 0.0, N, 0, dyn,
                                    probe_dtdX, probe_Omg) < 0) {
      return 1;
    }
    nd_dt = probe_dtdX; nd_Omg = probe_Omg;
    if (eccpa_mean_freqs(nd_dt, nd_Omg, N, &n_h, &omg_h)) return 1;
  }

  for (int j = 0; j < N; j++) {
    f[j] = n_h*nd_dt[j] - 1.0;
    g[j] = (nd_Omg[j] - omg_h)*nd_dt[j];
  }

  EccSeries fs, gs;
  eob_ecc_series_fit(f, N, &fs);
  eob_ecc_series_fit(g, N, &gs);
  map_coefs->M = (fs.M > ECCPA_HARMONICS_MAX) ? ECCPA_HARMONICS_MAX : fs.M;

  /* harmonic m divided by i m: (-b_m, -a_m)/m, times 1/r_a^m */
  double p = 1.0;
  for (int m = 1; m <= map_coefs->M; m++) {
    p *= ra_use;
    const double d = m*p;
    map_coefs->kappa[TH_LAMBDA_MEAN][0][m] = -fs.b[m]/d;
    map_coefs->kappa[TH_LAMBDA_MEAN][1][m] = -fs.a[m]/d;
    /* M_lambda = -V, so opposite sign */
    map_coefs->kappa[TH_AZIMUTH][0][m] = gs.b[m]/d;
    map_coefs->kappa[TH_AZIMUTH][1][m] = gs.a[m]/d;
  }
  return 0;
}


/** eccpa_cons_map_coefs on orbit (pphi, r_a) with periastron along (alpha, beta); 1 on failure */
static int eccpa_cons_map_coefs_of_orbit(PostAdDynamics *pa, double J, double r_a, double alpha, double beta,
                                         int oslot, int N, Dynamics *dyn, EccMapCoefs *map_coefs)
{
  double dtdX[ECCPA_NODES_DIM], Omg[ECCPA_NODES_DIM];
  if (eob_eccpa_cons_orbit(pa, J, r_a, dyn, oslot)) return 1;
  if (eob_eccpa_cons_sample_nodes(pa, oslot, alpha, beta, N, 0, dyn, dtdX, Omg) < 0) return 1;
  return eccpa_cons_map_coefs(pa, oslot, dtdX, Omg, N, dyn, map_coefs);
}

/** change of U = alpha^2 + beta^2 for a change (dalpha, dbeta) of z: 2 (alpha dalpha + beta dbeta). */
static double eccpa_dU(double alpha, double beta, double dalpha, double dbeta)
{
  return 2.0*(alpha*dalpha + beta*dbeta);
}

/** change of a quantity f(pphi, U) for a change (dpphi, dalpha, dbeta) of elements. 
 *  Used for frequencies, omega_,mu c~_mu in Eq. (B5). */
static double eccpa_dfreq(double f_J, double f_U, double alpha, double beta,
                          double dpphi, double dalpha, double dbeta)
{
  return f_J*dpphi + f_U*eccpa_dU(alpha, beta, dalpha, dbeta);
}

/**
 * Function: eccpa_phase_forcing
 * -----------------------------
 *   F_theta = M_theta,mu F_mu at one point (App. A). kappa_m derivatives are one-sided
 *   differences against map_coefs_dJ, map_coefs_dU; conj(z)^m is differentiated exactly.
 *   If ctJ != NULL, also M_Lambda,mu c~(1)_mu into M_ct_out.
 */
static void eccpa_phase_forcing(const EccMapCoefs *map_coefs, const EccMapCoefs *map_coefs_dJ,
                                const EccMapCoefs *map_coefs_dU, double hJ, double dU, double r_a,
                                double chi, double X, double alpha, double beta,
                                double FJ, double Fa, double Fb,
                                const double *ctJ, double cta, double ctb,
                                double *F_th_out, double *M_ct_out)
{
  const int M = map_coefs->M;
  double rap[ECCPA_HARMONICS_MAX+2];
  eccpa_rapow(rap, ECCPA_HARMONICS_MAX+1, r_a);
  /* map_coefs_dU is a difference in U, so it is multiplied by rate of U,
     Udot = 2 (alpha Fa + beta Fb), or by change of U along c~(1) */
  const double Udot = eccpa_dU(alpha, beta, Fa, Fb);
  const double dU_ct = ctJ ? eccpa_dU(alpha, beta, cta, ctb) : 0.0;
  double F_th[TH_N] = { 0.0, 0.0 }, M_ct = 0.0;
  const double cX = cos(X), sX = sin(X);
  double c[ECCPA_HARMONICS_MAX+1], s[ECCPA_HARMONICS_MAX+1];
  eccpa_trig(chi, M, c, s);
  for (int m = 1; m <= M; m++) {
    /* conj(z)^m exp(i m X) = r_a^m exp(i m chi) = Pr + i Pi_, and
       conj(z)^(m-1) exp(i m X) = r_a^(m-1) exp(i (m-1) chi) exp(i X) = qr + i qi */
    const double Pr = rap[m]*c[m],  Pi_ = rap[m]*s[m];
    const double qr = rap[m-1]*(c[m-1]*cX - s[m-1]*sX);
    const double qi = rap[m-1]*(c[m-1]*sX + s[m-1]*cX);
    /* conj(zdot) = Fa - i Fb, and Re[ (gr + i gi)(qr + i qi)(Fa - i Fb) ] = gr ur - gi ui */
    const double ur = qr*Fa + qi*Fb, ui = qi*Fa - qr*Fb;

    for (int th = 0; th < TH_N; th++) {
      const double *k0 = map_coefs->kappa[th][0], *k1 = map_coefs->kappa[th][1];
      /* first part: change of kappa_m, by finite differences */
      const double kr = (map_coefs_dJ->kappa[th][0][m] - k0[m])/hJ*FJ
                      + (map_coefs_dU->kappa[th][0][m] - k0[m])/dU*Udot;
      const double ki = (map_coefs_dJ->kappa[th][1][m] - k1[m])/hJ*FJ
                      + (map_coefs_dU->kappa[th][1][m] - k1[m])/dU*Udot;
      F_th[th] += kr*Pr - ki*Pi_;
      /* second part: change of conj(z)^m, m conj(z)^(m-1) conj(zdot) */
      const double gr = k0[m]*m, gi = k1[m]*m;
      F_th[th] += gr*ur - gi*ui;
    }

    if (ctJ) {
      const double *k0 = map_coefs->kappa[TH_LAMBDA_MEAN][0], *k1 = map_coefs->kappa[TH_LAMBDA_MEAN][1];
      const double kr = (map_coefs_dJ->kappa[TH_LAMBDA_MEAN][0][m] - k0[m])/hJ*(*ctJ)
                      + (map_coefs_dU->kappa[TH_LAMBDA_MEAN][0][m] - k0[m])/dU*dU_ct;
      const double ki = (map_coefs_dJ->kappa[TH_LAMBDA_MEAN][1][m] - k1[m])/hJ*(*ctJ)
                      + (map_coefs_dU->kappa[TH_LAMBDA_MEAN][1][m] - k1[m])/dU*dU_ct;
      const double gr = k0[m]*m, gi = k1[m]*m;
      const double uYr = qr*cta + qi*ctb;
      const double uYi = qi*cta - qr*ctb;
      M_ct += kr*Pr - ki*Pi_ + gr*uYr - gi*uYi;
    }
  }
  for (int th = 0; th < TH_N; th++) F_th_out[th] = F_th[th];
  if (M_ct_out) *M_ct_out = M_ct;
}

/** T and dT/dX at anomaly chi, for an orbit of width r_a. */
static void eccpa_Tmap(const EccMapCoefs *map_coefs, double r_a, double chi, double *T, double *dT)
{
  *T = eccpa_ksum(map_coefs->kappa[TH_LAMBDA_MEAN][0], map_coefs->kappa[TH_LAMBDA_MEAN][1], map_coefs->M,
                  r_a, chi, dT);
}

/** X with X + T(X) = Lambda: Newton from X0, bisection on [Lambda - pi, Lambda + pi] if it doesn't settle */
static double eccpa_X_of_Lam(const EccMapCoefs *map_coefs, double r_a, double w,
                             double Lam, double X0)
{
  double X = X0, step = 1.0;
  for (int it = 0; it < 8; it++) {
    double T, dT;
    eccpa_Tmap(map_coefs, r_a, X - w, &T, &dT);
    const double fp = 1.0 + dT;
    if (!(fabs(fp) > 0.0)) break;
    step = (X + T - Lam)/fp;
    X -= step;
    if (fabs(step) < 1e-14) break;
  }
  if (fabs(step) < 1e-10*(1.0 + fabs(X))) return X;
  double lo = Lam - Pi, hi = Lam + Pi;
  for (int it = 0; it < 64; it++) {
    double T;
    X = 0.5*(lo + hi);
    eccpa_Tmap(map_coefs, r_a, X - w, &T, NULL);
    if (X + T > Lam) hi = X; else lo = X;
  }
  return X;
}

/** Primitive of Eq. (C7) on nodes, with zero time average */
static void eccpa_osc_primitive(const double *src, const double *dtdX, int N, double P_r, double *out)
{
  double q[ECCPA_NODES_DIM];
  for (int j = 0; j < N; j++) q[j] = src[j]*dtdX[j];
  EccSeries qs, Gs;
  eob_ecc_series_fit(q, N, &qs);
  eob_ecc_series_prim(&qs, &Gs);
  const double node_spacing = 2.0*Pi/N;
  double tav = 0.0;
  for (int j = 0; j < N; j++) {
    out[j] = eob_ecc_series_eval(&Gs, 2.0*Pi*(j + 0.5)/N);
    tav += out[j]*dtdX[j]*node_spacing;
  }
  tav /= P_r;
  for (int j = 0; j < N; j++) out[j] -= tav;
}

/** Rotate z parts by -w, i.e. measure them from periastron (smooth along evolution) */
static void eccpa_to_peri(const double *za, const double *zb, int N,
                          double alpha, double beta, double *zr, double *zi)
{
  const double w  = eccpa_peri_angle(alpha, beta);
  const double cw = cos(w), sw = sin(w);
  for (int j = 0; j < N; j++) {
    zr[j] =  za[j]*cw + zb[j]*sw;
    zi[j] = -za[j]*sw + zb[j]*cw;
  }
}

/** Inverse of eccpa_to_peri for one value; cw, sw = cos w, sin w */
static void eccpa_from_peri(double zr, double zi, double cw, double sw, double *za, double *zb)
{
  *za = zr*cw - zi*sw;
  *zb = zr*sw + zi*cw;
}

/** Lagrange weights at X on m points x[i0 .. i0+m-1] */
static void eccpa_lagrange(const double *x, int i0, int m, double X, double *w)
{
  for (int a = 0; a < m; a++) {
    double L = 1.0;
    for (int b = 0; b < m; b++)
      if (b != a) L *= (X - x[i0+b])/(x[i0+a] - x[i0+b]);
    w[a] = L;
  }
}

/** a[0] + sum_m a[m] cos(m chi) + b[m] sin(m chi), up to M. For maps, which are stored
 *  divided by r_a^m, use eccpa_ksum. */
static double eccpa_csum(const double *a, const double *b, int M, double chi)
{
  double c[ECCPA_HARMONICS_MAX+1], s[ECCPA_HARMONICS_MAX+1], v = a[0];
  eccpa_trig(chi, M, c, s);
  for (int m = 1; m <= M; m++) v += a[m]*c[m] + b[m]*s[m];
  return v;
}

/** first-order oscillating parts at anomaly chi, from their Fourier coefficients in pa->osc1_coef
 *  (z parts measured from periastron, see eccpa_to_peri). */
static void eccpa_osc1_at(const PostAdDynamics *pa, double chi, double *g)
{
  double c[ECCPA_HARMONICS_MAX+1], s[ECCPA_HARMONICS_MAX+1];
  eccpa_trig(chi, ECCPA_HARMONICS_MAX, c, s);
  for (int a = 0; a < OSC_N; a++) {
    g[a] = pa->osc1_coef[a].a[0];
    for (int m = 1; m <= pa->osc1_coef[a].M; m++) g[a] += pa->osc1_coef[a].a[m]*c[m] + pa->osc1_coef[a].b[m]*s[m];
  }
}

/** Fourier coefficients in chi of one oscillating-part component from its node values; returns harmonics kept */
static int eccpa_osc_coefs(const double *g, int N, double *a, double *b)
{
  EccSeries f;
  eob_ecc_series_fit(g, N, &f);
  for (int m = 0; m <= ECCPA_HARMONICS_MAX; m++) {
    a[m] = (m <= f.M) ? f.a[m] : 0.0;
    b[m] = (m <= f.M) ? f.b[m] : 0.0;
  }
  return (f.M > ECCPA_HARMONICS_MAX) ? ECCPA_HARMONICS_MAX : f.M;
}

/** eccpa_osc_coefs of one quantity on nodes, written to checkpoint k of corr. */
static void eccpa_osc_coefs_to_corr(const double *g, int N,
                                    double corr[2][ECCPA_HARMONICS_MAX+1][ECCPA_CKPT_MAX], int k)
{
  double a[ECCPA_HARMONICS_MAX+1], b[ECCPA_HARMONICS_MAX+1];
  eccpa_osc_coefs(g, N, a, b);
  for (int m = 0; m <= ECCPA_HARMONICS_MAX; m++) {
    corr[0][m][k] = a[m];
    corr[1][m][k] = b[m];
  }
}

/** One component of a higher-order oscillating part, from its values on N nodes, at
 *  anomaly chi. */
static double eccpa_osc2_col(const double *col, int N, double chi)
{
  double a[ECCPA_HARMONICS_MAX+1], b[ECCPA_HARMONICS_MAX+1];
  eccpa_osc_coefs(col, N, a, b);
  return eccpa_csum(a, b, ECCPA_HARMONICS_MAX, chi);
}

/**
 * Function: eob_eccpa_osc1
 * ------------------------
 *   First-order oscillating parts c~(1) and theta~(1), Eqs. (C7), (B2), and averaged angle
 *   rates Fbar^(1)_theta, Eq. (B1), on mean orbit. Fourier coefficients into pa->osc1_coef.
 *   Also leaves in pa map_coefs (and U neighbour), domega_dJ/dU, F_th and M_ctilde1.
 *
 *   @param[in]     J, alpha, beta : state
 *   @param[in]     N_fix : node count to use; 0 = from eccentricity
 *   @param[in]     dyn   : Dynamics structure
 *   @param[in,out] pa    : struct
 *   @return  0 on success
 */
int eob_eccpa_osc1(PostAdDynamics *pa, double J, double alpha, double beta, int N_fix,
                   Dynamics *dyn)
{
  memset(pa->ctilde1, 0, sizeof pa->ctilde1);
  memset(pa->thtilde1,  0, sizeof pa->thtilde1);
  memset(pa->M_ctilde1,   0, sizeof pa->M_ctilde1);
  memset(pa->osc1_coef, 0, sizeof pa->osc1_coef);
  memset(&pa->map_coefs, 0, sizeof pa->map_coefs);
  for (int th = 0; th < TH_N; th++) pa->Fbar_th[th] = 0.0;
  pa->osc1_nnodes = 0;
  const double r_a = sqrt(alpha*alpha + beta*beta);
  if (eob_eccpa_cons_orbit(pa, J, r_a, dyn, PO_CUR)) return 1;

  /* Nodes for oscillating parts (see eob_eccpa_nodecount). This is only a request, since
     eob_eccpa_rates may raise it; N below is what was used. */
  const int N_nodes_osc = (N_fix > 0) ? N_fix : eccpa_nnodes(pa, pa->cons_orb[PO_CUR][PO_E_GEO], 1);
  if (eob_eccpa_rates(pa, J, alpha, beta, N_nodes_osc, dyn)) return 1;
  const int N = pa->nnodes;
  pa->osc1_nnodes = N;

  const double node_spacing   = 2.0*Pi/N;
  const double P_r = pa->P_r;

  /* c~(1), primitive of F_mu - Fbar_mu */
  for (int mu = 0; mu < C_N; mu++) {
    double q[ECCPA_NODES_DIM];
    for (int j = 0; j < N; j++) q[j] = pa->node[PA_F + mu][j] - pa->Fbar[mu];
    eccpa_osc_primitive(q, pa->node[PA_DTDX], N, P_r, pa->ctilde1[mu]);
  }

  /* Fourier coefficients, z parts measured from periastron */
  {
    double zr[ECCPA_NODES_DIM], zi[ECCPA_NODES_DIM];
    eccpa_to_peri(pa->ctilde1[C_ALPHA], pa->ctilde1[C_BETA], N, alpha, beta, zr, zi);
    const double *src[3] = { pa->ctilde1[C_PPHI], zr, zi };
    EccOscCoefs *dst[3] = { &pa->osc1_coef[OSC_PPHI], &pa->osc1_coef[OSC_Z_RE], &pa->osc1_coef[OSC_Z_IM] };
    for (int a = 0; a < 3; a++) dst[a]->M = eccpa_osc_coefs(src[a], N, dst[a]->a, dst[a]->b);
  }

  /* map coefs here and at pphi + hJ, U + dU (one-sided diffs) */
  if (eccpa_cons_map_coefs(pa, PO_CUR, pa->node[PA_DTDX], pa->node[PA_OMG], N, dyn, &pa->map_coefs)) return 1;

  const double hJ    = ECCPA_DPPHI_REL*fabs(J);
  const double U     = r_a*r_a;
  const double U_nbr = U + SQ(ECCPA_DU_RA*pa->cons_orb[PO_CUR][PO_R_C]);
  EccMapCoefs  map_coefs_dJ, map_coefs_dU;
  if (eccpa_cons_map_coefs_of_orbit(pa, J + hJ, r_a, alpha, beta, PO_AUX, N, dyn, &map_coefs_dJ)) return 1;
  /* U neighbour: same pphi and z direction, longer z */
  {
    const double scale = (r_a > 0.0) ? sqrt(U_nbr)/r_a : 0.0;
    const double aA = (r_a > 0.0) ? alpha*scale : sqrt(U_nbr);
    const double bA = (r_a > 0.0) ? beta*scale  : 0.0;
    if (eccpa_cons_map_coefs_of_orbit(pa, J, sqrt(U_nbr), aA, bA, PO_AUX, N, dyn, &map_coefs_dU)) return 1;
  }

  /* width derivatives in U; dU = exact diff of two U */
  const double dU = SQ(sqrt(U_nbr)) - U;
  pa->map_coefs_dU = map_coefs_dU;
  pa->map_coefs_dU_step = dU;
  /* frequencies differentiated in pphi and U; higher-order rates use them too */
  for (int th = 0; th < TH_N; th++) {
    pa->domega_dJ[th] = (map_coefs_dJ.omega[th] - pa->map_coefs.omega[th])/hJ;
    pa->domega_dU[th] = (map_coefs_dU.omega[th] - pa->map_coefs.omega[th])/dU;
  }

  double sum_domega1[TH_N] = { 0.0, 0.0 }, sum_F_th[TH_N] = { 0.0, 0.0 };
  double domega1[TH_N][ECCPA_NODES_DIM];
  for (int j = 0; j < N; j++) {
    const double chi = 2.0*Pi*(j + 0.5)/N;
    const double X   = pa->node[PA_X][j];
    const double dt  = pa->node[PA_DTDX][j]*node_spacing;
    /* instantaneous rates at node */
    const double FJ = pa->node[PA_F + C_PPHI][j], Fa = pa->node[PA_F + C_ALPHA][j], Fb = pa->node[PA_F + C_BETA][j];

    /* domega^(1) at node, Eq. (B5); angle parts below need it at every node */
    for (int th = 0; th < TH_N; th++) {
      domega1[th][j] = eccpa_dfreq(pa->domega_dJ[th], pa->domega_dU[th], alpha, beta,
                                   pa->ctilde1[C_PPHI][j], pa->ctilde1[C_ALPHA][j], pa->ctilde1[C_BETA][j]);
      sum_domega1[th] += domega1[th][j]*dt;
    }

    double F_th[TH_N], M_ct;
    eccpa_phase_forcing(&pa->map_coefs, &map_coefs_dJ, &map_coefs_dU, hJ, dU, r_a, chi, X, alpha, beta,
                        FJ, Fa, Fb, &pa->ctilde1[C_PPHI][j], pa->ctilde1[C_ALPHA][j], pa->ctilde1[C_BETA][j],
                        F_th, &M_ct);
    pa->M_ctilde1[j] = M_ct;
    for (int th = 0; th < TH_N; th++) {
      pa->F_th[th][j] = F_th[th];
      sum_F_th[th] += F_th[th]*dt;
    }
  }

  for (int th = 0; th < TH_N; th++) {
    const double mean_domega1 = sum_domega1[th]/P_r;   /* zero up to roundoff */
    pa->Fbar_th[th] = mean_domega1 + sum_F_th[th]/P_r;
  }

  /* theta~(1), primitive of S^(1) - Fbar^(1), Eq. (B2) */
  for (int th = 0; th < TH_N; th++) {
    double q[ECCPA_NODES_DIM];
    for (int j = 0; j < N; j++) q[j] = domega1[th][j] + pa->F_th[th][j] - pa->Fbar_th[th];
    eccpa_osc_primitive(q, pa->node[PA_DTDX], N, P_r, pa->thtilde1[th]);
    EccOscCoefs *dst = &pa->osc1_coef[OSC_TH + th];
    dst->M = eccpa_osc_coefs(pa->thtilde1[th], N, dst->a, dst->b);
  }
  return 0;
}

/**
 * Function: eob_eccpa_phases
 * --------------------------
 *   t, Lambda and lambda at samples by quadratures of Eq. (45), on a fine grid kept in pa.
 *
 *   @param[in,out] pa : samples
 *   @param[in]     n  : number of samples
 *   @param[in]     n_sub_per_interval : subdivisions per interval; <= 0 = ECCPA_PHASE_NSUB
 *   @param[in]     dyn : Dynamics structure
 *   @return  0 on success
 */
int eob_eccpa_phases(PostAdDynamics *pa, int n, int n_sub_per_interval, Dynamics *dyn)
{
  if (n < 4) return 1;
  if (n_sub_per_interval <= 0) n_sub_per_interval = ECCPA_PHASE_NSUB;
  const int fine_n = (n - 1)*n_sub_per_interval + 1;

  /* seven arrays of length fine_n in one allocation, kept in pa->fine_arrays at end */
  double *fine_arrays = (double *) malloc(7*fine_n*sizeof *fine_arrays);
  if (!fine_arrays) return 1;
  double *pphi_fine = fine_arrays;
  double *dt_dpphi  = fine_arrays + fine_n;
  double *t_fine    = fine_arrays + 4*fine_n;
  double *dth_dpphi[TH_N], *th_fine[TH_N];
  for (int th = 0; th < TH_N; th++) {
    dth_dpphi[th] = fine_arrays + (2 + th)*fine_n;
    th_fine[th]   = fine_arrays + (5 + th)*fine_n;
  }

  const int interp_freqs = EOBPars->postadiabatic_dynamics_ecc_interp_freqs;

  /* stencil width, max number of samples */
  int n_stencil = ECCPA_PHASE_STENCIL;
  if (n_stencil < 2) n_stencil = 2;
  if (n_stencil > n) n_stencil = n;

  int bad = 0;
  for (int i_interval = 0; i_interval < n - 1 && !bad; i_interval++) {
    const double Ja = pa->data[EVO_PPHI][i_interval], Jb = pa->data[EVO_PPHI][i_interval+1];
    const double dpphi = Jb - Ja;                      /* negative: pphi decreases */
    for (int s = 0; s < n_sub_per_interval; s++) {
      const int    k  = i_interval*n_sub_per_interval + s;
      const double frac_in_interval = (double) s/n_sub_per_interval;
      const double J  = Ja + frac_in_interval*dpphi;
      pphi_fine[k] = J;
      /* n, omega interpolated from samples; interp_freqs off (debug): solved on orb */
      double omega[TH_N] = { 0.0, 0.0 };
      if (!interp_freqs) {
        /* cubic Hermite in pphi, exact at both ends */
        const double al = herm_cubic(frac_in_interval, dpphi, pa->data[EVO_ALPHA][i_interval], pa->data[EVO_DADJ][i_interval],
                                     pa->data[EVO_ALPHA][i_interval+1], pa->data[EVO_DADJ][i_interval+1]);
        const double be = herm_cubic(frac_in_interval, dpphi, pa->data[EVO_BETA][i_interval], pa->data[EVO_DBDJ][i_interval],
                                     pa->data[EVO_BETA][i_interval+1], pa->data[EVO_DBDJ][i_interval+1]);
        if (eob_eccpa_cons_orbit(pa, J, sqrt(al*al + be*be), dyn, PO_CUR)) {
          bad = 1;
          break;
        }
        const int N = eccpa_nnodes(pa, pa->cons_orb[PO_CUR][PO_E_GEO], 0);
        if (eob_eccpa_cons_sample_nodes(pa, PO_CUR, al, be, N, 0, dyn, NULL, NULL) < 0) {
          bad = 1;
          break;
        }
        if (eccpa_mean_freqs(pa->node[PA_DTDX], pa->node[PA_OMG], N,
                             &omega[TH_LAMBDA_MEAN], &omega[TH_AZIMUTH])) {
          bad = 1;
          break;
        }
      }

      /* Fbar[C_PPHI] (FbarJ here), Fbar_th (and frequencies) by Lagrange on n_stencil
         samples, stencil moved inwards at two ends */
      int i0 = i_interval - (n_stencil/2 - 1);
      if (i0 < 0) i0 = 0;
      if (i0 > n - n_stencil) i0 = n - n_stencil;
      double FbarJ = 0.0, Fbar_th[TH_N] = { 0.0, 0.0 };
      double lag_w[ECCPA_PHASE_STENCIL];
      eccpa_lagrange(pa->data[EVO_PPHI], i0, n_stencil, J, lag_w);
      for (int a = 0; a < n_stencil; a++) {
        const double L = lag_w[a];
        FbarJ += L*pa->data[EVO_FBAR_PPHI][i0+a];
        for (int th = 0; th < TH_N; th++) {
          Fbar_th[th] += L*pa->data[EVO_FBAR_TH + th][i0+a];
          if (interp_freqs) omega[th] += L*pa->data[EVO_OMEGA + th][i0+a];
        }
      }
      if (!(FbarJ < 0.0)) {
        bad = 1;
        break;
      }
      dt_dpphi[k] = 1.0/FbarJ;
      for (int th = 0; th < TH_N; th++) dth_dpphi[th][k] = (omega[th] + Fbar_th[th])/FbarJ;
    }
  }
  if (!bad) {
    /* last point of fine grid is last sample */
    const int k = fine_n - 1;
    pphi_fine[k] = pa->data[EVO_PPHI][n-1];
    dt_dpphi[k] = 1.0/pa->data[EVO_FBAR_PPHI][n-1];
    for (int th = 0; th < TH_N; th++)
      dth_dpphi[th][k] = (pa->data[EVO_OMEGA + th][n-1] + pa->data[EVO_FBAR_TH + th][n-1])
                         /pa->data[EVO_FBAR_PPHI][n-1];

    cumint3(dt_dpphi, pphi_fine, fine_n, t_fine);
    for (int th = 0; th < TH_N; th++) cumint3(dth_dpphi[th], pphi_fine, fine_n, th_fine[th]);
    for (int i = 0; i < n; i++) {
      const int k_sample = i*n_sub_per_interval;
      pa->time[i] = t_fine[k_sample];
      for (int th = 0; th < TH_N; th++) pa->data[EVO_TH + th][i] = th_fine[th][k_sample];
    }
    /* Kept for reconstruction, which reads angles on this grid; integrands are the
       derivatives its cubic Hermite needs. */
    free(pa->fine_arrays);
    pa->fine_arrays = fine_arrays;  pa->fine_n = fine_n;
    pa->fine_pphi = pphi_fine;  pa->fine_t = t_fine;  pa->fine_dt_dpphi = dt_dpphi;
    for (int th = 0; th < TH_N; th++) {
      pa->fine_th[th] = th_fine[th];
      pa->fine_dth_dpphi[th] = dth_dpphi[th];
    }
    return 0;
  }
  free(fine_arrays);
  return bad;
}

/** Cubic Lagrange of arr[i0 .. i0+3] with weights w from eccpa_lagrange */
static double eccpa_lag(const double *w, int i0, const double *col)
{
  double v = 0.0;
  for (int a = 0; a < 4; a++) v += w[a] * col[i0+a];
  return v;
}

/** Cubic Lagrange weights in pphi on checkpoints, and index they start at. With fewer
 *  than four checkpoints, linear between two nearest. */
static int eccpa_corr_weights(const PostAdDynamics *corr, double J, double *wt)
{
  int i0 = 0;                        /* ckpt_pphi decreases; find interval holding J */
  while (i0 < corr->ckpt_n - 2 && corr->ckpt_pphi[i0+1] > J) i0++;
  if (corr->ckpt_n < 4) {
    const double d = corr->ckpt_pphi[i0+1] - corr->ckpt_pphi[i0];
    const double u = (d != 0.0) ? (J - corr->ckpt_pphi[i0])/d : 0.0;
    wt[0] = 1.0 - u;  wt[1] = u;  wt[2] = 0.0;  wt[3] = 0.0;
    return i0;
  }
  i0--;
  if (i0 < 0) i0 = 0;
  if (i0 > corr->ckpt_n - 4) i0 = corr->ckpt_n - 4;
  eccpa_lagrange(corr->ckpt_pphi, i0, 4, J, wt);
  return i0;
}

/** Second- or third-order oscillating parts at (J, chi), interpolated from checkpoints */
static void eccpa_osc2_at(const PostAdDynamics *corr, int ord, double J, double chi, double *out)
{
  for (int sr = 0; sr < OSC_N; sr++) out[sr] = 0.0;
  /* at NLO checkpoints hold only rate corrections */
  if (corr->ckpt_n < 2 || EOBPars->postadiabatic_dynamics_ecc_map < ord
      || (ord >= 3 && corr->corr3_n < 2)) return;
  /* third order uses same checkpoints as second */
  const double (*ts)[2][ECCPA_HARMONICS_MAX+1][ECCPA_CKPT_MAX] = (ord >= 3) ? corr->corr3_osc_coef : corr->corr2_osc_coef;
  double wt[4];
  const int i0 = eccpa_corr_weights(corr, J, wt);

  /* m = 0 is not zero, because eccpa_osc_primitive makes time average zero, not average
     in chi */
  for (int sr = 0; sr < OSC_N; sr++)
    out[sr] = eccpa_lag(wt, i0, ts[sr][0][0]);

  double c[ECCPA_HARMONICS_MAX+1], s[ECCPA_HARMONICS_MAX+1];
  eccpa_trig(chi, ECCPA_HARMONICS_MAX, c, s);
  for (int m = 1; m <= ECCPA_HARMONICS_MAX; m++)
    for (int sr = 0; sr < OSC_N; sr++) {
      const double am = eccpa_lag(wt, i0, ts[sr][0][m]), bm = eccpa_lag(wt, i0, ts[sr][1][m]);
      out[sr] += am*c[m] + bm*s[m];
    }
}

/** Step (a) of reconstruction: mean pphi and angles at time t, on fine grid; 1 if no grid */
static int eccpa_at_time(const PostAdDynamics *pa, double t, double *J, double *th_bar)
{
  if (pa->fine_n < 2) return 1;
  /* last fine point with t_k <= t (fine_t increases) */
  const int kf = find_point_bisection(t, pa->fine_n, pa->fine_t, 1);
  const double dpphi = pa->fine_pphi[kf+1] - pa->fine_pphi[kf];
  double frac_in_interval = (pa->fine_t[kf+1] > pa->fine_t[kf])
                          ? (t - pa->fine_t[kf])/(pa->fine_t[kf+1] - pa->fine_t[kf]) : 0.0;
  for (int it = 0; it < 4; it++) {
    const double t2 = frac_in_interval*frac_in_interval;
    const double f  = herm_cubic(frac_in_interval, dpphi, pa->fine_t[kf], pa->fine_dt_dpphi[kf],
                                 pa->fine_t[kf+1], pa->fine_dt_dpphi[kf+1]) - t;
    const double fp = ( 6*t2 - 6*frac_in_interval)*pa->fine_t[kf]
                    + (3*t2 - 4*frac_in_interval + 1)*dpphi*pa->fine_dt_dpphi[kf]
                    + (-6*t2 + 6*frac_in_interval)*pa->fine_t[kf+1]
                    + (3*t2 - 2*frac_in_interval)*dpphi*pa->fine_dt_dpphi[kf+1];
    if (!(fabs(fp) > 0.0)) break;
    const double step = f/fp;
    frac_in_interval -= step;
    if (frac_in_interval < -0.5) frac_in_interval = -0.5;
    if (frac_in_interval >  1.5) frac_in_interval =  1.5;
    if (fabs(step) < 1e-14) break;
  }
  *J   = pa->fine_pphi[kf] + frac_in_interval*dpphi;
  for (int th = 0; th < TH_N; th++)
    th_bar[th] = herm_cubic(frac_in_interval, dpphi, pa->fine_th[th][kf], pa->fine_dth_dpphi[th][kf],
                            pa->fine_th[th][kf+1], pa->fine_dth_dpphi[th][kf+1]);
  return 0;
}

/** Map coefficients of osculating orbit, times ra^m, moved from mean U with one Taylor term
 *  (dropped term is O(epsilon^2)). */
static void eccpa_osc_map(const PostAdDynamics *pa, int i0, const double *w_osc, int M,
                          double ra, EccMapCoefs *map_coefs_osc)
{
  double rap[ECCPA_HARMONICS_MAX+1];
  eccpa_rapow(rap, M, ra);
  const double abar = eccpa_lag(w_osc, i0, pa->data[EVO_ALPHA]);
  const double bbar = eccpa_lag(w_osc, i0, pa->data[EVO_BETA]);
  const double dU = ra*ra - (abar*abar + bbar*bbar);
  memset(map_coefs_osc, 0, sizeof *map_coefs_osc);
  map_coefs_osc->M = M;
  for (int k = 1; k <= M; k++)
    for (int th = 0; th < TH_N; th++)
      for (int c = 0; c < 2; c++)
        map_coefs_osc->kappa[th][c][k] = rap[k]*(eccpa_lag(w_osc, i0, pa->evo_kappa[th][c][k])
                                                 + dU*eccpa_lag(w_osc, i0, pa->evo_kappa_dU[th][c][k]));
}

/**
 * Function: eob_eccpa_recon
 * -------------------------
 *   Phase-space point at time t, steps (a)-(d) of reconstruction (App. C).
 *   Meant to be called in increasing t; pa->rec_nseed = 0 resets seeds.
 *
 *   @param[in,out] pa  : evolution; point is written to pa
 *   @param[in]     n   : number of samples
 *   @param[in]     t   : time, in [pa->time[0], pa->time[n-1]]
 *   @param[in]     dyn : Dynamics structure
 *   @return  0 on success
 */
int eob_eccpa_recon(PostAdDynamics *pa, int n, double t, Dynamics *dyn)
{
  pa->t_at = pa->phi = pa->pphi = pa->X = 0.0;
  if (n < 4) return 1;
  pa->t_at = t;

  /* (a) */
  double th_bar[TH_N], J;
  if (eccpa_at_time(pa, t, &J, th_bar)) return 1;
  const double Lam = th_bar[TH_LAMBDA_MEAN];

  /* (b) zbar, cubic Hermite in pphi */
  /* first interval whose lower end is at or below J (pphi decreases along samples) */
  int i = 0, i_hi = n - 2;
  while (i < i_hi) {
    const int mid = (i + i_hi)/2;
    if (pa->data[EVO_PPHI][mid+1] > J) i = mid + 1; else i_hi = mid;
  }
  const double dpphi = pa->data[EVO_PPHI][i+1] - pa->data[EVO_PPHI][i];
  const double frac_in_interval = (dpphi != 0.0) ? (J - pa->data[EVO_PPHI][i])/dpphi : 0.0;
  const double abar = herm_cubic(frac_in_interval, dpphi, pa->data[EVO_ALPHA][i], pa->data[EVO_DADJ][i],
                                 pa->data[EVO_ALPHA][i+1], pa->data[EVO_DADJ][i+1]);
  const double bbar = herm_cubic(frac_in_interval, dpphi, pa->data[EVO_BETA][i], pa->data[EVO_DBDJ][i],
                                 pa->data[EVO_BETA][i+1], pa->data[EVO_DBDJ][i+1]);

  int i0 = i - 1;
  if (i0 < 0) i0 = 0;
  if (i0 > n - 4) i0 = n - 4;
  /* cubic weights at mean pphi */
  double w_mean[4];
  eccpa_lagrange(pa->data[EVO_PPHI], i0, 4, J, w_mean);

  const double rbar = sqrt(abar*abar + bbar*bbar);
  const double wbar = eccpa_peri_angle(abar, bbar);
  /* Use largest M of four samples. Where node count goes up, a neighbour has
     harmonics this sample doesn't, and dropping them would make T and V jump */
  int M = pa->evo_M[i0];
  for (int a = 1; a < 4; a++) if (pa->evo_M[i0+a] > M) M = pa->evo_M[i0+a];

  /* (c) X_m: oscillating parts read on mean orb, Eq. (50) */
  EccMapCoefs map_coefs_mean;
  memset(&map_coefs_mean, 0, sizeof map_coefs_mean);
  map_coefs_mean.M = M;
  for (int k = 1; k <= M; k++) {
    map_coefs_mean.kappa[TH_LAMBDA_MEAN][0][k] = eccpa_lag(w_mean, i0, pa->evo_kappa[TH_LAMBDA_MEAN][0][k]);
    map_coefs_mean.kappa[TH_LAMBDA_MEAN][1][k] = eccpa_lag(w_mean, i0, pa->evo_kappa[TH_LAMBDA_MEAN][1][k]);
  }
  /* seed is a linear extrapolation in Lambda from last two points; after a jump of more
     than a radian, Lambda itself */
  double Xseed = Lam;
  if (pa->rec_nseed >= 2 && fabs(Lam - pa->rec_L1) < 1.0 && pa->rec_L1 != pa->rec_L2) {
    Xseed = pa->rec_X1 + (pa->rec_X1 - pa->rec_X2)*(Lam - pa->rec_L1)/(pa->rec_L1 - pa->rec_L2);
  }
  const double Xq = eccpa_X_of_Lam(&map_coefs_mean, rbar, wbar, Lam, Xseed);
  pa->rec_L2 = pa->rec_L1;
  pa->rec_X2 = pa->rec_X1;
  pa->rec_L1 = Lam;
  pa->rec_X1 = Xq;
  pa->rec_nseed++;
  const double chi = Xq - wbar;

  /* first-order oscillating parts at this pphi and chi, then second and third order
     from checkpoints */
  double sc[OSC_N][2][ECCPA_HARMONICS_MAX+1];
  for (int sr = 0; sr < OSC_N; sr++)
    for (int k = 0; k <= M; k++) {
      sc[sr][0][k] = eccpa_lag(w_mean, i0, pa->evo_osc1_coef[sr][0][k]);
      sc[sr][1][k] = eccpa_lag(w_mean, i0, pa->evo_osc1_coef[sr][1][k]);
    }
  double acc[OSC_N], ck[ECCPA_HARMONICS_MAX+1], sk[ECCPA_HARMONICS_MAX+1];
  eccpa_trig(chi, M, ck, sk);
  for (int sr = 0; sr < OSC_N; sr++) acc[sr] = sc[sr][0][0];
  for (int k = 1; k <= M; k++)
    for (int sr = 0; sr < OSC_N; sr++) acc[sr] += sc[sr][0][k]*ck[k] + sc[sr][1][k]*sk[k];
  double osc2[OSC_N], osc3[OSC_N];
  eccpa_osc2_at(pa, 2, J, chi, osc2);
  eccpa_osc2_at(pa, 3, J, chi, osc3);   /* zero below map order 3 */

  /* (d) osculating elements; z parts are rotated back by wbar */
  const double ct_zr = acc[OSC_Z_RE] + osc2[OSC_Z_RE] + osc3[OSC_Z_RE];
  const double ct_zi = acc[OSC_Z_IM] + osc2[OSC_Z_IM] + osc3[OSC_Z_IM];
  double ct_za, ct_zb;
  eccpa_from_peri(ct_zr, ct_zi, cos(wbar), sin(wbar), &ct_za, &ct_zb);
  const double alpha = abar + ct_za;
  const double beta  = bbar + ct_zb;
  const double pphi  = J + acc[OSC_PPHI] + osc2[OSC_PPHI] + osc3[OSC_PPHI];
  double th_osc[TH_N];
  for (int th = 0; th < TH_N; th++)
    th_osc[th] = th_bar[th] + acc[OSC_TH + th] + osc2[OSC_TH + th] + osc3[OSC_TH + th];
  const double LamP = th_osc[TH_LAMBDA_MEAN];
  const double ra = sqrt(alpha*alpha + beta*beta);
  const double wo = eccpa_peri_angle(alpha, beta);
  /* X from T of osculating orbit, where point is; T of mean orbit would leave an
     O(epsilon) error in X */
  double w_osc[4];
  eccpa_lagrange(pa->data[EVO_PPHI], i0, 4, pphi, w_osc);
  EccMapCoefs map_coefs_osc;
  eccpa_osc_map(pa, i0, w_osc, M, ra, &map_coefs_osc);
  const double X = eccpa_X_of_Lam(&map_coefs_osc, 1.0, wo, LamP, Xq + (LamP - Lam));
  pa->X = X;
  const double M_lam = eccpa_ksum(map_coefs_osc.kappa[TH_AZIMUTH][0], map_coefs_osc.kappa[TH_AZIMUTH][1], M,
                                  1.0, X - wo, NULL);
  /* r_circ interpolated at osculating pphi, H derivatives recomputed */
  {
    const double rc_i = eccpa_lag(w_osc, i0, pa->data[EVO_RCIRC]);
    /* previous point's root, a seed for eccpa_cons_rc */
    const double rc_prev = pa->cons_orb[PO_CUR][PO_R_C];
    memset(pa->cons_orb[PO_CUR], 0, sizeof pa->cons_orb[PO_CUR]);
    pa->cons_orb[PO_CUR][PO_R_C]  = rc_prev;
    pa->cons_orb[PO_CUR][PO_PPHI] = pphi;
    pa->cons_orb[PO_CUR][PO_R_A]  = ra;
    if (eccpa_cons_circular(pa, pphi, dyn, PO_CUR, rc_i)) return 1;
    if (eccpa_cons_rc(pa, pa->cons_orb[PO_CUR][PO_R_A], dyn, PO_CUR)) return 1;
  }
  eob_eccpa_cons_sample_point(pa, PO_CUR, alpha, beta, X, 0, dyn);

  /* r, prstar, Omg, rdot and ddotr are in pa already, from eob_eccpa_cons_sample_point */
  pa->pphi     = pphi;
  pa->phi      = th_osc[TH_AZIMUTH] - M_lam;
  return 0;
}

/** One map coefficient at a displaced state by Taylor expansion on six-orbit stencil */
static double eccpa_map_tay(double f00, double fJp, double fJm, double fUp, double fUm,
                            double fC, double HJ, double HU, double dJ, double dU,
                            double *dfdJ, double *dfdU)
{
  const double gJ  = 0.5*(fJp - fJm)/HJ;
  const double gU  = 0.5*(fUp - fUm)/HU;
  const double gJJ = (fJp - 2.0*f00 + fJm)/(HJ*HJ);
  const double gUU = (fUp - 2.0*f00 + fUm)/(HU*HU);
  const double gJU = (fC - fJp - fUp + f00)/(HJ*HU);   /* from one corner of stencil */
  if (dfdJ) *dfdJ = gJ + gJJ*dJ + gJU*dU;
  if (dfdU) *dfdU = gU + gJU*dJ + gUU*dU;
  return f00 + gJ*dJ + gU*dU + 0.5*(gJJ*dJ*dJ + gUU*dU*dU) + gJU*dJ*dU;
}

/** Map coefficients at a displaced state and at its two neighbours for eccpa_phase_forcing */
static void eccpa_maps_tay(const PostAdDynamics *pa, double J, double Js, double r_a,
                           double r_c, EccMapCoefs *map_coefs, EccMapCoefs *map_coefs_dJ, EccMapCoefs *map_coefs_dU,
                           double *hJ, double *dU)
{
  const double U     = r_a*r_a;
  const double U_nbr = U + SQ(ECCPA_DU_RA*r_c);
  *hJ = ECCPA_DPPHI_REL*fabs(Js);
  *dU = SQ(sqrt(U_nbr)) - U;
  const double dJs = Js - J, dUs = U - pa->map_coefs_fd_U0;
  double dJ_, dU_;
  const EccMapCoefs *st = pa->map_coefs_fd;
  const double hpphi = pa->map_coefs_fd_h_pphi, hU = pa->map_coefs_fd_h_U;
  for (int th = 0; th < TH_N; th++) {
    map_coefs->omega[th] = eccpa_map_tay(st[MAPFD_00].omega[th], st[MAPFD_JP].omega[th], st[MAPFD_JM].omega[th],
                                         st[MAPFD_UP].omega[th], st[MAPFD_UM].omega[th], st[MAPFD_CR].omega[th],
                                         hpphi, hU, dJs, dUs, &dJ_, &dU_);
    map_coefs_dJ->omega[th] = map_coefs->omega[th] + *hJ*dJ_;  map_coefs_dU->omega[th] = map_coefs->omega[th] + *dU*dU_;
  }
  map_coefs->M = map_coefs_dJ->M = map_coefs_dU->M = st[MAPFD_00].M;
  for (int k = 0; k <= map_coefs->M; k++)
    for (int th = 0; th < TH_N; th++)
      for (int c = 0; c < 2; c++) {
        double *f = &map_coefs->kappa[th][c][k];
        *f = eccpa_map_tay(st[MAPFD_00].kappa[th][c][k], st[MAPFD_JP].kappa[th][c][k], st[MAPFD_JM].kappa[th][c][k],
                           st[MAPFD_UP].kappa[th][c][k], st[MAPFD_UM].kappa[th][c][k], st[MAPFD_CR].kappa[th][c][k],
                           hpphi, hU, dJs, dUs, &dJ_, &dU_);
        map_coefs_dJ->kappa[th][c][k] = *f + *hJ*dJ_;  map_coefs_dU->kappa[th][c][k] = *f + *dU*dU_;
      }
}

/**
 * Function: eob_eccpa_rates_higher
 * --------------------------------
 *   Higher-order rates and oscillating parts of App. B on mean orbit, as far as orders ask.
 *   S^(3) terms that need neighbouring checkpoints are added in eob_eccpa_corr_build.
 *   With rates, pa->Fbar = Fbar^(1) + Fbarcorr; c~(2), theta~(2), S^(3) go to pa->mean_node.
 *
 *   @param[in]     J, alpha, beta : state
 *   @param[in]     rates : 1 for rates, 0 for c~(2) and theta~(2) only
 *   @param[in]     dyn   : Dynamics structure
 *   @param[in,out] pa    : struct
 *   @return  0 on success
 */
int eob_eccpa_rates_higher(PostAdDynamics *pa, double J, double alpha, double beta, int rates,
                           Dynamics *dyn)
{
  for (int mu = 0; mu < C_N; mu++) pa->Fbarcorr[mu] = pa->Fbar2[mu] = 0.0;
  for (int th = 0; th < TH_N; th++) pa->Fbar2_th[th] = 0.0;
  pa->s3_ok = 0;
  /* What orders ask for */
  const int map = EOBPars->postadiabatic_dynamics_ecc_map, rate = EOBPars->postadiabatic_dynamics_ecc_rate;
  const int want_osc2  = !rates || map >= 2 || rate >= 3;
  const int want_pass2 = rates && (rate >= 3 || map >= 3);
  int want_S3 = rates && map >= 3;   /* dropped if a third-order piece can't be built */
  if (eob_eccpa_osc1(pa, J, alpha, beta, 0, dyn)) return 1;
  const int N = pa->osc1_nnodes;
  if (N < 4) return 0;                       /* nothing to displace against */

  const double node_spacing = 2.0*Pi/N;
  const double P_r = pa->P_r, n = pa->omega[TH_LAMBDA_MEAN];
  const double omega[TH_N] = { pa->omega[TH_LAMBDA_MEAN], pa->omega[TH_AZIMUTH] };
  const double Fbar1_th[TH_N] = { pa->Fbar_th[TH_LAMBDA_MEAN], pa->Fbar_th[TH_AZIMUTH] };
  const double Fbar1[C_N] = { pa->Fbar[C_PPHI], pa->Fbar[C_ALPHA], pa->Fbar[C_BETA] };
  const double dadJ = pa->dc_dpphi[C_ALPHA], dbdJ = pa->dc_dpphi[C_BETA];
  const double domega_dJ[TH_N] = { pa->domega_dJ[TH_LAMBDA_MEAN], pa->domega_dJ[TH_AZIMUTH] };
  const double domega_dU[TH_N] = { pa->domega_dU[TH_LAMBDA_MEAN], pa->domega_dU[TH_AZIMUTH] };
  const double r_a = sqrt(alpha*alpha + beta*beta);
  const double cbar[C_N] = { J, alpha, beta };

  /* neighbours and displaced states in pa_nbr, pa keeps mean state */
  PostAdDynamics *pa_nbr = NULL;
  if (PostAdDynamics_alloc(&pa_nbr, 1)) return 1;
  pa_nbr->nmax = pa->nmax;
  int status = 1;

  /* Lambdabar at each node; neighbours are compared at same Lambdabar */
  for (int j = 0; j < N; j++) {
    double T;
    eccpa_Tmap(&pa->map_coefs, r_a, 2.0*Pi*(j + 0.5)/N, &T, NULL);
    pa->mean_node[MEAN_NODE_LAM][j] = pa->node[PA_X][j] + T;
  }

  /* mean orbit's r_c and periastron, for third order, and its kappa_m for the
     convolutions of z corrections below */
  const double rc_m = pa->cons_orb[PO_CUR][PO_R_C], w0_m = eccpa_peri_angle(alpha, beta);
  double kappa_mean_re[ECCPA_HARMONICS_MAX+1], kappa_mean_im[ECCPA_HARMONICS_MAX+1];
  if (eccpa_kappa(pa, PO_CUR, pa->node[PA_DTDX], N, dyn, kappa_mean_re, kappa_mean_im)) goto out;

  /* map coefs of displaced states from a 6-orb stencil in (pphi, U) */
  double HJ = 0.0, HU = 0.0;
  pa->map_coefs_fd_U0 = r_a*r_a;
  for (int j = 0; j < N; j++) {
    const double dU_Y = fabs(pa->ctilde1[C_PPHI][j]);
    const double az = fabs(pa->ctilde1[C_ALPHA][j])
                    + fabs(pa->ctilde1[C_BETA][j]);
    if (dU_Y > HJ) HJ = dU_Y;
    if (az > HU) HU = az;
  }
  HU = 2.0*r_a*HU + HU*HU;                       /* largest change of U, not of |z| */
  if (!(HJ > 0.0)) HJ = ECCPA_DPPHI_REL_OSC*fabs(J);
  if (!(HU > 0.0)) HU = ECCPA_DPPHI_REL_OSC*pa->map_coefs_fd_U0;
  /* r_a < displacements: centre stencil at U0 = HU */
  if (!(pa->map_coefs_fd_U0 - HU > 0.0)) {
    pa->map_coefs_fd_U0 = HU;
    HU *= 0.5;
  }
  if (!(HU > 0.0)) goto out;                     /* nothing to build a stencil on */
  const double ra_c = sqrt(pa->map_coefs_fd_U0);          /* r_a at centre of stencil */
  pa->map_coefs_fd_h_pphi = HJ; pa->map_coefs_fd_h_U = HU;
  /* six orbits in order of MAPFD_*; periastron along alpha, since coefficients
     don't depend on w */
  const double J_st[MAPFD_N]  = { J, J + HJ, J - HJ, J, J, J + HJ };
  const double ra_st[MAPFD_N] = { ra_c, ra_c, ra_c, sqrt(pa->map_coefs_fd_U0 + HU),
                                  sqrt(pa->map_coefs_fd_U0 - HU), sqrt(pa->map_coefs_fd_U0 + HU) };
  for (int k = 0; k < MAPFD_N; k++)
    if (eccpa_cons_map_coefs_of_orbit(pa, J_st[k], ra_st[k], ra_st[k], 0.0, PO_AUX, N, dyn, &pa->map_coefs_fd[k]))
      goto out;

  /* first pass: +-c~(1); half-diff = dF^(2), half-sum = quadratic part */
  double dF2[C_N][ECCPA_NODES_DIM];
  double domega2_quad[TH_N][ECCPA_NODES_DIM], dF2_th[TH_N][ECCPA_NODES_DIM];
  double dX1[ECCPA_NODES_DIM];
  double sum_domega2[TH_N] = { 0.0, 0.0 }, sum_dF_th2[TH_N] = { 0.0, 0.0 };
  /* sum_dF2 = P_r Fbar^(2), time average of dF^(2) not yet divided by P_r */
  double sum_dF2[C_N] = { 0.0, 0.0, 0.0 };

  for (int j = 0; j < N; j++) {
    /* dX^(1), Eq. (B6) [eq:dX_eq] */
    const double dX = (pa->node[PA_DTDX][j] > 0.0)
                    ? (pa->thtilde1[TH_LAMBDA_MEAN][j] - pa->M_ctilde1[j])
                      /(n*pa->node[PA_DTDX][j]) : 0.0;
    dX1[j] = dX;
    double F_s[C_N][2], omega_s[TH_N][2], F_th_s[TH_N][2];
    for (int sg = 0; sg < 2; sg++) {
      const double s = sg ? -1.0 : 1.0;
      double c_P[C_N];
      for (int mu = 0; mu < C_N; mu++) c_P[mu] = cbar[mu] + s*pa->ctilde1[mu][j];
      const double X_P = pa->node[PA_X][j] + s*dX;
      if (eob_eccpa_cons_orbit(pa_nbr, c_P[C_PPHI], sqrt(c_P[C_ALPHA]*c_P[C_ALPHA] + c_P[C_BETA]*c_P[C_BETA]), dyn, PO_CUR)
          || eob_eccpa_cons_sample_point(pa_nbr, PO_CUR, c_P[C_ALPHA], c_P[C_BETA], X_P, 1, dyn)) goto out;
      for (int mu = 0; mu < C_N; mu++) F_s[mu][sg] = pa_nbr->F[mu];

      /* n, omega and F_Lambda, F_lambda at this state, from stencil */
      EccMapCoefs map_coefs_P, map_coefs_P_dJ, map_coefs_P_dU;
      const double r_a_P = pa_nbr->cons_orb[PO_CUR][PO_R_A], r_c_P = pa_nbr->cons_orb[PO_CUR][PO_R_C];
      double hJP, dUP;
      const double w_P = eccpa_peri_angle(c_P[C_ALPHA], c_P[C_BETA]);
      eccpa_maps_tay(pa, J, c_P[C_PPHI], r_a_P, r_c_P, &map_coefs_P, &map_coefs_P_dJ, &map_coefs_P_dU, &hJP, &dUP);
      double F_th[TH_N];
      eccpa_phase_forcing(&map_coefs_P, &map_coefs_P_dJ, &map_coefs_P_dU, hJP, dUP, r_a_P, X_P - w_P, X_P,
                          c_P[C_ALPHA], c_P[C_BETA],
                          F_s[C_PPHI][sg], F_s[C_ALPHA][sg], F_s[C_BETA][sg], NULL, 0.0, 0.0, F_th, NULL);
      for (int th = 0; th < TH_N; th++) {
        omega_s[th][sg] = map_coefs_P.omega[th];
        F_th_s[th][sg]  = F_th[th];
      }
    }

    const double dt = pa->node[PA_DTDX][j]*node_spacing;
    for (int mu = 0; mu < C_N; mu++) dF2[mu][j] = 0.5*(F_s[mu][0] - F_s[mu][1]);
    for (int mu = 0; mu < C_N; mu++) sum_dF2[mu] += dF2[mu][j]*pa->node[PA_DTDX][j]*node_spacing;
    for (int th = 0; th < TH_N; th++) {
      domega2_quad[th][j] = 0.5*(omega_s[th][0] + omega_s[th][1]) - omega[th];
      dF2_th[th][j]  = 0.5*(F_th_s[th][0] - F_th_s[th][1]);
      sum_domega2[th] += domega2_quad[th][j]*dt;
      sum_dF_th2[th]  += dF2_th[th][j]*dt;
    }
  }

  /* Fbar^(2)_Lambda and Fbar^(2)_lambda, Eq. (B1) */
  for (int th = 0; th < TH_N; th++) pa->Fbar2_th[th] = (sum_domega2[th] + sum_dF_th2[th])/P_r;

  /* second order of z rate alone, by same convolution as first order */
  double Fbar2_z[2] = { 0.0, 0.0 };
  if (rates)
    eccpa_convolve(dF2[C_ALPHA], dF2[C_BETA], pa->node[PA_X], N, kappa_mean_re, kappa_mean_im, alpha, beta, n,
                   &Fbar2_z[0], &Fbar2_z[1]);

  /* D^(1), D^(2) of first-order parts, Eq. (B3), centred diffs along evolution */
  const double hJ = ECCPA_DPPHI_REL_OSC*fabs(J);
  if (!(hJ > 0.0)) goto out;
  const double Fbar2[C_N] = { sum_dF2[C_PPHI]/P_r, Fbar2_z[0], Fbar2_z[1] };
  /* D^(2)'s step: largest relative change of three elements is ECCPA_DPPHI_REL_OSC */
  const double scale[C_N] = { fabs(J), rc_m, rc_m };
  double rel = fabs(Fbar2[C_PPHI])/scale[C_PPHI];
  for (int mu = 1; mu < C_N; mu++)
    if (fabs(Fbar2[mu])/scale[mu] > rel) rel = fabs(Fbar2[mu])/scale[mu];
  const double h2 = (rel > 0.0) ? ECCPA_DPPHI_REL_OSC/rel : 0.0;
  const double h_d[2] = { hJ, h2 }, fac_d[2] = { Fbar1[C_PPHI], 1.0 };
  const double v_d[2][C_N] = { { 1.0, dadJ, dbdJ }, { Fbar2[C_PPHI], Fbar2[C_ALPHA], Fbar2[C_BETA] } };
  double D1ctilde1[C_N][ECCPA_NODES_DIM] = {{0.0}}, D1thtilde1[TH_N][ECCPA_NODES_DIM] = {{0.0}};
  double D2ctilde1[C_N][ECCPA_NODES_DIM] = {{0.0}}, D2thtilde1[TH_N][ECCPA_NODES_DIM] = {{0.0}};
  double (*Dc[2])[ECCPA_NODES_DIM] = { D1ctilde1, D2ctilde1 };
  double (*Dth[2])[ECCPA_NODES_DIM] = { D1thtilde1, D2thtilde1 };

  for (int d = 0; d < 2; d++) {
    if (d == 0 ? !want_osc2 : !(want_S3 && h2 > 0.0)) break;   /* D^(2) needs D^(1)'s orders */
    for (int sg = 0; sg < 2; sg++) {
      /* 1 */
      const double s = sg ? -1.0 : 1.0;
      double c_n[C_N];
      for (int mu = 0; mu < C_N; mu++) c_n[mu] = cbar[mu] + s*h_d[d]*v_d[d][mu];
      /* 2, checked because a wrong derivative gives a plausible higher-order rate, not an error */
      if (eob_eccpa_osc1(pa_nbr, c_n[C_PPHI], c_n[C_ALPHA], c_n[C_BETA], N, dyn)
          || pa_nbr->osc1_nnodes != N) {
        if (d == 0) goto out;
        want_S3 = 0;
        break;
      }
      const double ran = sqrt(c_n[C_ALPHA]*c_n[C_ALPHA] + c_n[C_BETA]*c_n[C_BETA]);
      const double wn  = eccpa_peri_angle(c_n[C_ALPHA], c_n[C_BETA]);
      const double cw = cos(wn), sw = sin(wn);
      for (int j = 0; j < N; j++) {
        /* 3 */
        const double Xn = eccpa_X_of_Lam(&pa_nbr->map_coefs, ran, wn, pa->mean_node[MEAN_NODE_LAM][j],
                                         pa->node[PA_X][j]);
        double acc[OSC_N];
        eccpa_osc1_at(pa_nbr, Xn - wn, acc);
        double ctilde1_nbr[C_N];
        ctilde1_nbr[C_PPHI] = acc[OSC_PPHI];
        eccpa_from_peri(acc[OSC_Z_RE], acc[OSC_Z_IM], cw, sw, &ctilde1_nbr[C_ALPHA], &ctilde1_nbr[C_BETA]);
        /* 4 */
        for (int mu = 0; mu < C_N; mu++)
          Dc[d][mu][j] += s*fac_d[d]*ctilde1_nbr[mu]/(2.0*h_d[d]);
        for (int th = 0; th < TH_N; th++)
          Dth[d][th][j] += s*fac_d[d]*acc[OSC_TH + th]/(2.0*h_d[d]);
      }
    }
  }

  /* Kept for third order are sources S^(2) (derivative of an oscillating part along
     Lambdabar is its source over n) and n_,mu c~(1), n_,mu c~(2), and same for omega. */
  double S2[C_N][ECCPA_NODES_DIM], S2_th[TH_N][ECCPA_NODES_DIM];
  double domega1[TH_N][ECCPA_NODES_DIM], domega2_lin[TH_N][ECCPA_NODES_DIM];
  /* -- second-order oscillating parts, from S^(2) of Eq. (B2) --------------------------------
     sources average to zero by construction, so primitive only removes roundoff. */
  if (want_osc2) {
    for (int mu = 0; mu < C_N; mu++) {
      const double Fbar2_mu = sum_dF2[mu]/P_r;
      for (int j = 0; j < N; j++)
        S2[mu][j] = dF2[mu][j] - Fbar2_mu - D1ctilde1[mu][j]
                  - (pa->node[PA_F + mu][j] - Fbar1[mu])/n*Fbar1_th[TH_LAMBDA_MEAN];
      eccpa_osc_primitive(S2[mu], pa->node[PA_DTDX], N, P_r, pa->mean_node[MEAN_NODE_CTILDE2 + mu]);
    }
    /* theta~(2), Eq. (B5) */
    for (int th = 0; th < TH_N; th++) {
      for (int j = 0; j < N; j++) {
        domega1[th][j] = eccpa_dfreq(domega_dJ[th], domega_dU[th], alpha, beta,
                                     pa->ctilde1[C_PPHI][j],
                                     pa->ctilde1[C_ALPHA][j],
                                     pa->ctilde1[C_BETA][j]);
        domega2_lin[th][j] = eccpa_dfreq(domega_dJ[th], domega_dU[th], alpha, beta,
                                        pa->mean_node[MEAN_NODE_CTILDE2 + C_PPHI][j],
                                        pa->mean_node[MEAN_NODE_CTILDE2 + C_ALPHA][j],
                                        pa->mean_node[MEAN_NODE_CTILDE2 + C_BETA][j]);
        S2_th[th][j] = domega2_quad[th][j] + domega2_lin[th][j] + dF2_th[th][j] - pa->Fbar2_th[th]
                     - D1thtilde1[th][j]
                     - (domega1[th][j] + pa->F_th[th][j] - Fbar1_th[th])/n*Fbar1_th[TH_LAMBDA_MEAN];
      }
      eccpa_osc_primitive(S2_th[th], pa->node[PA_DTDX], N, P_r, pa->mean_node[MEAN_NODE_THTILDE2 + th]);
    }
  }

  if (rates)
    for (int mu = 0; mu < C_N; mu++) pa->Fbar2[mu] = Fbar2[mu];
  if (want_pass2) {
    /* -- second pass: each node at cbar + c~(1) + c~(2) ----------------------------------------- */
    double F_P[C_N][ECCPA_NODES_DIM], omega_P[TH_N][ECCPA_NODES_DIM];
    double F_th_P[TH_N][ECCPA_NODES_DIM], F_th_stencil[TH_N][ECCPA_NODES_DIM];
    double sum_dF23_pphi = 0.0;
    double dF23_a[ECCPA_NODES_DIM], dF23_b[ECCPA_NODES_DIM];   /* z rates at state minus mean orbit's */
    /* for third order: maps on mean orbit, same at every node */
    EccMapCoefs mc_mean, mc_mean_dJ, mc_mean_dU;
    double hJ_mean = 0.0, dU_mean = 0.0;
    if (want_S3) eccpa_maps_tay(pa, J, J, r_a, rc_m, &mc_mean, &mc_mean_dJ, &mc_mean_dU, &hJ_mean, &dU_mean);
    for (int j = 0; j < N; j++) {
      double c_P[C_N];
      for (int mu = 0; mu < C_N; mu++)
        c_P[mu] = cbar[mu] + pa->ctilde1[mu][j] + pa->mean_node[MEAN_NODE_CTILDE2 + mu][j];
      const double r_a_P = sqrt(c_P[C_ALPHA]*c_P[C_ALPHA] + c_P[C_BETA]*c_P[C_BETA]);
      const double w_P = eccpa_peri_angle(c_P[C_ALPHA], c_P[C_BETA]);
      if (eob_eccpa_cons_orbit(pa_nbr, c_P[C_PPHI], r_a_P, dyn, PO_CUR)) goto out;

      /* X that carries Lambdabar + Lambda~(1) + Lambda~(2), on displaced orbit's own T, so
         nonlinear terms of map are kept */
      double X_P = pa->node[PA_X][j] + dX1[j];
      int map_P_ok = 0;
      EccMapCoefs map_coefs_orbit_P;
      if (!eccpa_cons_map_coefs_of_orbit(pa_nbr, c_P[C_PPHI], r_a_P, c_P[C_ALPHA], c_P[C_BETA], PO_CUR, N, dyn, &map_coefs_orbit_P)) {
        map_P_ok = 1;
        X_P = eccpa_X_of_Lam(&map_coefs_orbit_P, r_a_P, w_P,
                            pa->mean_node[MEAN_NODE_LAM][j] + pa->thtilde1[TH_LAMBDA_MEAN][j]
                            + pa->mean_node[MEAN_NODE_THTILDE2 + TH_LAMBDA_MEAN][j], X_P);
      }
      if (eob_eccpa_cons_sample_point(pa_nbr, PO_CUR, c_P[C_ALPHA], c_P[C_BETA], X_P, 1, dyn)) goto out;
      const double Fa_P = pa_nbr->F[C_ALPHA], Fb_P = pa_nbr->F[C_BETA];

      /* third order: F, n, omega, F_theta at this state */
      if (want_S3) {
        for (int mu = 0; mu < C_N; mu++) F_P[mu][j] = pa_nbr->F[mu];
        if (!map_P_ok) want_S3 = 0;
        else {
          for (int th = 0; th < TH_N; th++) omega_P[th][j] = map_coefs_orbit_P.omega[th];
          double F_th[TH_N], hJs, dUs;
          EccMapCoefs mc, mc_dJ, mc_dU;
          /* displaced state */
          eccpa_maps_tay(pa, J, c_P[C_PPHI], r_a_P, pa_nbr->cons_orb[PO_CUR][PO_R_C], &mc, &mc_dJ, &mc_dU, &hJs, &dUs);
          eccpa_phase_forcing(&mc, &mc_dJ, &mc_dU, hJs, dUs, r_a_P, X_P - w_P, X_P, c_P[C_ALPHA], c_P[C_BETA],
                              pa_nbr->F[C_PPHI], Fa_P, Fb_P, NULL, 0.0, 0.0, F_th, NULL);
          for (int th = 0; th < TH_N; th++) F_th_P[th][j] = F_th[th];
          /* mean orbit, with maps computed once before loop */
          const double X0 = pa->node[PA_X][j];
          eccpa_phase_forcing(&mc_mean, &mc_mean_dJ, &mc_mean_dU, hJ_mean, dU_mean, r_a, X0 - w0_m, X0,
                              alpha, beta, pa->node[PA_F + C_PPHI][j],
                              pa->node[PA_F + C_ALPHA][j], pa->node[PA_F + C_BETA][j],
                              NULL, 0.0, 0.0, F_th, NULL);
          for (int th = 0; th < TH_N; th++) F_th_stencil[th][j] = F_th[th];
        }
      }
      /* F_E and F_pphi are averaged with mean orbit's dt/dX, Fa and Fb go to convolution */
      const double dt = pa->node[PA_DTDX][j]*node_spacing;
      sum_dF23_pphi += (pa_nbr->F[C_PPHI] - pa->node[PA_F + C_PPHI][j])*dt;
      dF23_a[j] = Fa_P - pa->node[PA_F + C_ALPHA][j];
      dF23_b[j] = Fb_P - pa->node[PA_F + C_BETA][j];
    }

    /* S^(3), Eq. (B2), without neighbouring-checkpoint terms */
    if (want_S3) {
      double dF3[C_N][ECCPA_NODES_DIM], domega3[TH_N][ECCPA_NODES_DIM], dF3_th[TH_N][ECCPA_NODES_DIM];
      double Fbar3[C_N] = { 0.0, 0.0, 0.0 }, Fbar3_th[TH_N] = { 0.0, 0.0 };
      for (int j = 0; j < N; j++) {
        const double dt = pa->node[PA_DTDX][j]*node_spacing;
        for (int mu = 0; mu < C_N; mu++) {
          dF3[mu][j] = F_P[mu][j] - pa->node[PA_F + mu][j] - dF2[mu][j];
          Fbar3[mu] += dF3[mu][j]*dt;
        }
        for (int th = 0; th < TH_N; th++) {
          domega3[th][j] = omega_P[th][j] - omega[th] - domega1[th][j] - domega2_lin[th][j] - domega2_quad[th][j];
          dF3_th[th][j]  = (F_th_P[th][j] - F_th_stencil[th][j]) - dF2_th[th][j];
          Fbar3_th[th] += (domega3[th][j] + dF3_th[th][j])*dt;
        }
      }
      for (int mu = 0; mu < C_N; mu++) Fbar3[mu] /= P_r;
      for (int th = 0; th < TH_N; th++) Fbar3_th[th] /= P_r;
      for (int j = 0; j < N; j++) {
        for (int mu = 0; mu < C_N; mu++)
          pa->mean_node[MEAN_NODE_S3_C + mu][j] = dF3[mu][j] - Fbar3[mu] - D2ctilde1[mu][j]
                               - (pa->node[PA_F + mu][j] - Fbar1[mu])/n*pa->Fbar2_th[TH_LAMBDA_MEAN] - S2[mu][j]/n*Fbar1_th[TH_LAMBDA_MEAN];
        for (int th = 0; th < TH_N; th++)
          pa->mean_node[MEAN_NODE_S3_TH + th][j] = domega3[th][j] + dF3_th[th][j] - Fbar3_th[th] - D2thtilde1[th][j]
                            - (domega1[th][j] + pa->F_th[th][j] - Fbar1_th[th])/n*pa->Fbar2_th[TH_LAMBDA_MEAN]
                            - S2_th[th][j]/n*Fbar1_th[TH_LAMBDA_MEAN];
      }
      pa->s3_ok  = 1;
    }

    /* second pass gives Fbar^(2) + Fbar^(3); z part by same convolution as the
       first order */
    pa->Fbarcorr[C_PPHI] = sum_dF23_pphi/P_r;
    eccpa_convolve(dF23_a, dF23_b, pa->node[PA_X], N, kappa_mean_re, kappa_mean_im, alpha, beta, n,
                   &pa->Fbarcorr[C_ALPHA], &pa->Fbarcorr[C_BETA]);
  } else {
    for (int mu = 0; mu < C_N; mu++) pa->Fbarcorr[mu] = pa->Fbar2[mu];
  }
  /* add to first-order rates, kept in Fbar1 since pa->Fbar is overwritten here */
  pa->Fbar[C_PPHI]  = Fbar1[C_PPHI]  + pa->Fbarcorr[C_PPHI];
  pa->Fbar[C_ALPHA] = Fbar1[C_ALPHA] + pa->Fbarcorr[C_ALPHA];
  pa->Fbar[C_BETA]  = Fbar1[C_BETA]  + pa->Fbarcorr[C_BETA];
  pa->dc_dpphi[C_ALPHA] = (pa->Fbar[C_PPHI] != 0.0) ? pa->Fbar[C_ALPHA]/pa->Fbar[C_PPHI] : 0.0;
  pa->dc_dpphi[C_BETA] = (pa->Fbar[C_PPHI] != 0.0) ? pa->Fbar[C_BETA]/pa->Fbar[C_PPHI] : 0.0;
  status = 0;
out:
  PostAdDynamics_free(&pa_nbr);
  return status;
}

/**
 * Function: eob_eccpa_corr_build
 * ------------------------------
 *   Higher-order corrections at checkpoints of LO evolution (App. C); at map order 3,
 *   a second loop builds c~(3).
 *
 *   @param[in,out] pa  : first evolution; corrections into pa->ckpt_*, corr2_*, corr3_*
 *   @param[in]     n   : number of samples
 *   @param[in]     dyn : Dynamics structure
 *   @return  0 if enough checkpoints were built
 */
int eob_eccpa_corr_build(PostAdDynamics *pa, int n, Dynamics *dyn)
{
  pa->ckpt_n = 0;
  pa->corr3_n = 0;
  if (n < 2) return 1;
  int n_ckpt = n/ECCPA_CKPT_EVERY;
  if (n_ckpt < 12) n_ckpt = 12;
  if (n_ckpt > ECCPA_CKPT_MAX) n_ckpt = ECCPA_CKPT_MAX;
  if (n_ckpt > n) n_ckpt = n;

  int n3 = 0;                                /* checkpoints that carry third order */
  for (int k = 0; k < n_ckpt; k++) {
    const int i = (n_ckpt == 1) ? 0 : (int) ((double) k*(n - 1)/(n_ckpt - 1) + 0.5);
    const double J  = pa->data[EVO_PPHI][i];
    const double al = pa->data[EVO_ALPHA][i], be = pa->data[EVO_BETA][i];
    /* first-order rates = Fbar - Fbarcorr; first failed checkpoint ends list */
    const int bad2 = eob_eccpa_rates_higher(pa, J, al, be, 1, dyn);
    if (bad2) break;   /* rest of evolution uses last checkpoint */
    const double Fbar1J = pa->Fbar[C_PPHI] - pa->Fbarcorr[C_PPHI];
    /* Rate order 3 stores whole correction, rate order 2 second order alone */
    const double *corr = (EOBPars->postadiabatic_dynamics_ecc_rate >= 3) ? pa->Fbarcorr : pa->Fbar2;

    pa->ckpt_pphi[pa->ckpt_n]    = J;
    pa->corr2_relJ[pa->ckpt_n] = (Fbar1J != 0.0) ? corr[C_PPHI]/Fbar1J : 0.0;
    for (int th = 0; th < TH_N; th++)
      pa->corr2_relth[th][pa->ckpt_n] = (pa->omega[th] != 0.0) ? pa->Fbar2_th[th]/pa->omega[th] : 0.0;
    /* z correction divided by z (times conj(z)/|z|^2), so it can be interpolated and is
       proportional to z when multiplied back. */
    {
      const double dr = corr[C_ALPHA], di = corr[C_BETA];
      const double r2 = al*al + be*be;
      pa->corr2_zr[pa->ckpt_n] = (r2 > 0.0) ? ( dr*al + di*be)/r2 : 0.0;
      pa->corr2_zi[pa->ckpt_n] = (r2 > 0.0) ? (-dr*be + di*al)/r2 : 0.0;
    }
    /* Fourier coefficients of c~(2), Lambda~(2), lambda~(2) at this checkpoint, z parts measured
       from periastron; only a map of order 2 or more reads them */
    if (EOBPars->postadiabatic_dynamics_ecc_map >= 2) {
      const int Nn = pa->nnodes;
      double zr[ECCPA_NODES_DIM], zi[ECCPA_NODES_DIM];
      eccpa_to_peri(pa->mean_node[MEAN_NODE_CTILDE2 + C_ALPHA], pa->mean_node[MEAN_NODE_CTILDE2 + C_BETA], Nn, al, be,
                    zr, zi);
      const double *src[OSC_N] = { pa->mean_node[MEAN_NODE_CTILDE2 + C_PPHI], zr, zi,
                                   pa->mean_node[MEAN_NODE_THTILDE2 + TH_LAMBDA_MEAN],
                                   pa->mean_node[MEAN_NODE_THTILDE2 + TH_AZIMUTH] };
      for (int sr = 0; sr < OSC_N; sr++) eccpa_osc_coefs_to_corr(src[sr], Nn, pa->corr2_osc_coef[sr], pa->ckpt_n);
    }
    if (EOBPars->postadiabatic_dynamics_ecc_map >= 3 && pa->s3_ok) {
      const int kk = pa->ckpt_n, Nn = pa->nnodes;
      pa->corr3_N[kk] = Nn;  pa->corr3_i[kk] = i;
      pa->corr3_Fbar1J[kk] = Fbar1J;  pa->corr3_Pr[kk] = pa->P_r;
      for (int th = 0; th < TH_N; th++) {
        pa->corr3_domega_dJ[th][kk] = pa->domega_dJ[th];
        pa->corr3_domega_dU[th][kk] = pa->domega_dU[th];
      }
      for (int j = 0; j < Nn; j++) {
        for (int mu = 0; mu < C_N; mu++) pa->corr3_src_c[mu][kk][j] = pa->mean_node[MEAN_NODE_S3_C + mu][j];
        for (int th = 0; th < TH_N; th++) pa->corr3_src_th[th][kk][j] = pa->mean_node[MEAN_NODE_S3_TH + th][j];
        pa->corr3_X[kk][j]    = pa->node[PA_X][j];
        pa->corr3_dtdX[kk][j] = pa->node[PA_DTDX][j];
        pa->corr3_Lam[kk][j]  = pa->mean_node[MEAN_NODE_LAM][j];
      }
      n3++;
    }
    pa->ckpt_n++;
  }

  /* At least three quarters of checkpoints must have been built; with fewer, cubic
     interpolation would span too much of inspiral. */
  const int ok2 = (pa->ckpt_n >= 2 && pa->ckpt_n >= (3*n_ckpt + 3)/4);

  /* third-order map: second loop over checkpoints */
  if (ok2 && EOBPars->postadiabatic_dynamics_ecc_map >= 3 && n3 == pa->ckpt_n && pa->ckpt_n >= 3) {
    const int K = pa->ckpt_n, m = (K >= 4) ? 4 : K;
    for (int k = 0; k < K; k++) {
      const int Nk = pa->corr3_N[k], ik = pa->corr3_i[k];
      int k0 = k - 1;
      if (k0 < 0) k0 = 0;
      if (k0 > K - m) k0 = K - m;
      double dw[4];
      for (int q = 0; q < m; q++) {                     /* d/dJ of q-th Lagrange weight at J_k */
        double sum = 0.0;
        for (int r = 0; r < m; r++) {
          if (r == q) continue;
          double p = 1.0/(pa->ckpt_pphi[k0+q] - pa->ckpt_pphi[k0+r]);
          for (int t = 0; t < m; t++)
            if (t != q && t != r) p *= (pa->ckpt_pphi[k] - pa->ckpt_pphi[k0+t])/(pa->ckpt_pphi[k0+q] - pa->ckpt_pphi[k0+t]);
          sum += p;
        }
        dw[q] = sum;
      }
      double src3_c[C_N][ECCPA_NODES_DIM], ctilde3[C_N][ECCPA_NODES_DIM];
      double src3_th[TH_N][ECCPA_NODES_DIM], thtilde3[TH_N][ECCPA_NODES_DIM];
      for (int j = 0; j < Nk; j++) {
        /* d/dpphi of c~(2) and theta~(2) */
        double dct2[C_N] = { 0.0, 0.0, 0.0 }, dtht2[TH_N] = { 0.0, 0.0 };
        for (int q = 0; q < m; q++) {
          const int kq = k0 + q, iq = pa->corr3_i[kq];
          EccMapCoefs map_coefs_q;
          memset(&map_coefs_q, 0, sizeof map_coefs_q);
          map_coefs_q.M = pa->evo_M[iq];
          for (int mm = 1; mm <= map_coefs_q.M; mm++)
            for (int c = 0; c < 2; c++)
              map_coefs_q.kappa[TH_LAMBDA_MEAN][c][mm] = pa->evo_kappa[TH_LAMBDA_MEAN][c][mm][iq];
          const double aq = pa->data[EVO_ALPHA][iq], bq = pa->data[EVO_BETA][iq];
          const double raq = sqrt(aq*aq + bq*bq), wq = eccpa_peri_angle(aq, bq);
          const double Xq = eccpa_X_of_Lam(&map_coefs_q, raq, wq, pa->corr3_Lam[k][j], pa->corr3_X[k][j]);
          const double chq = Xq - wq;
          double v[OSC_N];
          for (int sr = 0; sr < OSC_N; sr++) {
            v[sr] = pa->corr2_osc_coef[sr][0][0][kq];
            for (int mm = 1; mm <= ECCPA_HARMONICS_MAX; mm++)
              v[sr] += pa->corr2_osc_coef[sr][0][mm][kq]*cos(mm*chq) + pa->corr2_osc_coef[sr][1][mm][kq]*sin(mm*chq);
          }
          const double cq = cos(wq), sq = sin(wq);
          dct2[C_PPHI]  += dw[q]*v[OSC_PPHI];
          double za, zb;
          eccpa_from_peri(v[OSC_Z_RE], v[OSC_Z_IM], cq, sq, &za, &zb);
          dct2[C_ALPHA] += dw[q]*za;
          dct2[C_BETA]  += dw[q]*zb;
          for (int th = 0; th < TH_N; th++) dtht2[th] += dw[q]*v[OSC_TH + th];
        }
        for (int mu = 0; mu < C_N; mu++) src3_c[mu][j] = pa->corr3_src_c[mu][k][j] - pa->corr3_Fbar1J[k]*dct2[mu];
        for (int th = 0; th < TH_N; th++) src3_th[th][j] = pa->corr3_src_th[th][k][j] - pa->corr3_Fbar1J[k]*dtht2[th];
      }
      for (int mu = 0; mu < C_N; mu++)
        eccpa_osc_primitive(src3_c[mu], pa->corr3_dtdX[k], Nk, pa->corr3_Pr[k], ctilde3[mu]);
      const double ak = pa->data[EVO_ALPHA][ik], bk = pa->data[EVO_BETA][ik];
      for (int j = 0; j < Nk; j++)                       /* omega_,mu c~(3) */
        for (int th = 0; th < TH_N; th++)
          src3_th[th][j] += eccpa_dfreq(pa->corr3_domega_dJ[th][k], pa->corr3_domega_dU[th][k], ak, bk,
                                        ctilde3[C_PPHI][j], ctilde3[C_ALPHA][j], ctilde3[C_BETA][j]);
      for (int th = 0; th < TH_N; th++)
        eccpa_osc_primitive(src3_th[th], pa->corr3_dtdX[k], Nk, pa->corr3_Pr[k], thtilde3[th]);
      /* z parts measured from periastron, stored like c~(2) */
      double zr[ECCPA_NODES_DIM], zi[ECCPA_NODES_DIM];
      eccpa_to_peri(ctilde3[C_ALPHA], ctilde3[C_BETA], Nk, ak, bk, zr, zi);
      const double *srcs[OSC_N] = { ctilde3[C_PPHI], zr, zi, thtilde3[TH_LAMBDA_MEAN], thtilde3[TH_AZIMUTH] };
      for (int sr = 0; sr < OSC_N; sr++) eccpa_osc_coefs_to_corr(srcs[sr], Nk, pa->corr3_osc_coef[sr], k);
    }
    pa->corr3_n = K;
  }
  return ok2 ? 0 : 1;
}

/**
 * Function: eob_eccpa_corr_apply
 * ------------------------------
 *   Add checkpoint corrections, interpolated in pphi, to first-order rates.
 *   Also sets pa->Fbar2_th (angle-rate corrections, read by eccpa_store_sample).
 *
 *   @param[in]     corr_pa        : checkpoint corrections
 *   @param[in]     J, alpha, beta : state
 *   @param[in,out] pa             : rates
 */
void eob_eccpa_corr_apply(const PostAdDynamics *corr_pa, PostAdDynamics *pa,
                          double J, double alpha, double beta)
{
  if (corr_pa->ckpt_n < 2) return;

  double wt[4];
  const int i0 = eccpa_corr_weights(corr_pa, J, wt);
  const double relJ = eccpa_lag(wt, i0, corr_pa->corr2_relJ);
  const double zr   = eccpa_lag(wt, i0, corr_pa->corr2_zr);
  const double zi   = eccpa_lag(wt, i0, corr_pa->corr2_zi);

  pa->Fbar[C_PPHI] *= 1.0 + relJ;
  pa->Fbar[C_ALPHA] += zr*alpha - zi*beta;                                       /* (zr + i zi) z */
  pa->Fbar[C_BETA] += zr*beta  + zi*alpha;
  pa->dc_dpphi[C_ALPHA]  = (pa->Fbar[C_PPHI] != 0.0) ? pa->Fbar[C_ALPHA]/pa->Fbar[C_PPHI] : 0.0;
  pa->dc_dpphi[C_BETA]  = (pa->Fbar[C_PPHI] != 0.0) ? pa->Fbar[C_BETA]/pa->Fbar[C_PPHI] : 0.0;
  /* angle corrections are not part of z equations; eccpa_store_sample adds them to the
     stored angle rates */
  for (int th = 0; th < TH_N; th++) pa->Fbar2_th[th] = eccpa_lag(wt, i0, corr_pa->corr2_relth[th])*pa->omega[th];
}

/**
 * Function: eob_eccpa_to_mean
 * ---------------------------
 *   Osculating -> mean elements, Eq. (B7), at map order of reconstruction.
 *
 *   @param[in]     J0, r_a, X0  : osculating state, periastron along +alpha
 *   @param[in]     phi0         : azimuth
 *   @param[in]     order        : map order
 *   @param[in]     dyn          : Dynamics structure
 *   @param[in,out] pa           : struct
 *   @param[out]    Jm, am, bm   : mean elements
 *   @param[out]    th_m         : mean angles Lambda, lambda
 *   @param[in]     corr         : checkpoint corrections, or NULL
 *   @return  0 on success, 1 if oscillating parts could not be built
 */
int eob_eccpa_to_mean(PostAdDynamics *pa, double J0, double r_a, double X0, double phi0,
                      int order, Dynamics *dyn,
                      double *Jm, double *am, double *bm, double *th_m,
                      const PostAdDynamics *corr)
{
  *Jm = J0; *am = r_a; *bm = 0.0; th_m[TH_LAMBDA_MEAN] = X0; th_m[TH_AZIMUTH] = phi0;

  /* Lambda = X + T(X), lambda = phi - V(X) on osculating orbit */
  int osc_ok = !eob_eccpa_osc1(pa, J0, r_a, 0.0, 0, dyn);
  EccMapCoefs map_coefs_osc;
  if (osc_ok) {
    map_coefs_osc = pa->map_coefs;
  } else {
    if (eob_eccpa_cons_orbit(pa, J0, r_a, dyn, PO_CUR)) return 1;
    const int N0 = eccpa_nnodes(pa, pa->cons_orb[PO_CUR][PO_E_GEO], 0);
    if (eccpa_cons_map_coefs_of_orbit(pa, J0, r_a, r_a, 0.0, PO_CUR, N0, dyn, &map_coefs_osc)) return 1;
  }
  double T0;
  eccpa_Tmap(&map_coefs_osc, r_a, X0, &T0, NULL);
  double th_p[TH_N];
  th_p[TH_LAMBDA_MEAN] = X0 + T0;
  th_p[TH_AZIMUTH]     = phi0 + eccpa_ksum(map_coefs_osc.kappa[TH_AZIMUTH][0],
                                           map_coefs_osc.kappa[TH_AZIMUTH][1],
                                           map_coefs_osc.M, r_a, X0, NULL);
  for (int th = 0; th < TH_N; th++) th_m[th] = th_p[th];
  if (!osc_ok) return 1;
  /* one pass at first order, a second with higher orders */
  const int use_corr = (corr && order >= 3 && corr->corr3_n >= 2);
  /* Fourier coefficients of c~(2) (first C_N rows) and theta~(2) (next TH_N), at order 2 */
  double q2a[C_N + TH_N][ECCPA_HARMONICS_MAX+1], q2b[C_N + TH_N][ECCPA_HARMONICS_MAX+1];
  double Xq = X0;
  for (int stage = 0; stage < (order < 2 ? 1 : 2); stage++) {
    if (stage && !use_corr) {
      /* starting an order-2 evolution from first-order mean elements would be inconsistent, so
         fail instead */
      if (eob_eccpa_rates_higher(pa, *Jm, *am, *bm, 0, dyn)) return 1;
      /* second-order parts fitted once, where they were built */
      for (int mu = 0; mu < C_N; mu++)
        eccpa_osc_coefs(pa->mean_node[MEAN_NODE_CTILDE2 + mu], pa->nnodes, q2a[mu], q2b[mu]);
      for (int th = 0; th < TH_N; th++)
        eccpa_osc_coefs(pa->mean_node[MEAN_NODE_THTILDE2 + th], pa->nnodes, q2a[C_N + th], q2b[C_N + th]);
    }
    for (int it = 0; it < 4; it++) {
      /* first-order parts at current mean elements; very first pass reuses build
         above */
      if ((stage || it) && eob_eccpa_osc1(pa, *Jm, *am, *bm, 0, dyn)) break;
      const double rab = sqrt(*am * *am + *bm * *bm);
      const double wb  = eccpa_peri_angle(*am, *bm);
      /* X_m on current mean orbit */
      Xq = eccpa_X_of_Lam(&pa->map_coefs, rab, wb, th_m[TH_LAMBDA_MEAN], Xq);
      const double chi = Xq - wb;
      double g[OSC_N];
      eccpa_osc1_at(pa, chi, g);
      /* checkpoint z parts are from periastron, mean_node ones in (alpha, beta) */
      double dJ = J0 - g[OSC_PPHI], zr = g[OSC_Z_RE], zi = g[OSC_Z_IM], qa = 0.0, qb = 0.0;
      double th_new[TH_N];
      for (int th = 0; th < TH_N; th++) th_new[th] = th_p[th] - g[OSC_TH + th];
      if (stage && use_corr) {
        for (int o = 2; o <= 3; o++) {
          double osc[OSC_N];
          eccpa_osc2_at(corr, o, *Jm, chi, osc);
          dJ -= osc[OSC_PPHI]; zr += osc[OSC_Z_RE]; zi += osc[OSC_Z_IM];
          for (int th = 0; th < TH_N; th++) th_new[th] -= osc[OSC_TH + th];
        }
      } else if (stage) {
        dJ -= eccpa_csum(q2a[C_PPHI], q2b[C_PPHI], ECCPA_HARMONICS_MAX, chi);
        qa  = eccpa_csum(q2a[C_ALPHA], q2b[C_ALPHA], ECCPA_HARMONICS_MAX, chi);
        qb  = eccpa_csum(q2a[C_BETA], q2b[C_BETA], ECCPA_HARMONICS_MAX, chi);
        for (int th = 0; th < TH_N; th++)
          th_new[th] -= eccpa_csum(q2a[C_N + th], q2b[C_N + th], ECCPA_HARMONICS_MAX, chi);
      }
      double za, zb;
      eccpa_from_peri(zr, zi, cos(wb), sin(wb), &za, &zb);
      dJ -= *Jm;
      const double da = r_a - za - qa - *am;
      const double db =     - zb - qb - *bm;
      *Jm += dJ;
      *am += da;
      *bm += db;
      for (int th = 0; th < TH_N; th++) th_m[th] = th_new[th];
      /* each pass shrinks step by a factor O(epsilon) */
      if (fabs(dJ) + fabs(da) + fabs(db) < 1e-14*(fabs(J0) + r_a)) break;
    }
  }
  return 0;
}

/* evolution */

/**
 * Function: eccpa_z_rhs
 * ---------------------
 *   d(alpha, beta)/du with u = pphi(0) - pphi. GSL_EBADFUNC if no bound orbit exists.
 */
static int eccpa_z_rhs(double u, const double y[], double dy[], void *params)
{
  eccpa_evo_par *P = (eccpa_evo_par *) params;
  const double J = P->J0 - u;
  PostAdDynamics *pa = P->pa;
  const int status = eob_eccpa_rates(pa, J, y[0], y[1], 0, P->dyn);
  /* in second pass, add higher-order corrections from checkpoints */
  if (!status && P->corr) eob_eccpa_corr_apply(P->corr, pa, J, y[0], y[1]);
  if (J <= 0.0 || status
      || !(pa->Fbar[C_PPHI] < 0.0) || !isfinite(pa->Fbar[C_ALPHA]) || !isfinite(pa->Fbar[C_BETA])) {
    dy[0] = dy[1] = 0.0;
    return GSL_EBADFUNC;
  }
  P->last_FbarE = pa->FbarE;
  P->last_FbarJ = pa->Fbar[C_PPHI];
  P->last_dc_dpphi[C_ALPHA] = pa->dc_dpphi[C_ALPHA]; P->last_dc_dpphi[C_BETA] = pa->dc_dpphi[C_BETA];
  for (int th = 0; th < TH_N; th++) P->last_omega[th] = pa->omega[th];
  P->last_adiab = pa->adiab;
  P->last_N = pa->nnodes;
  /* d/du = -d/dpphi. pphi must keep decreasing (Fbar[C_PPHI] < 0, checked above, also false for
     NaN), so dz/du has sign of zdot. */
  dy[0] = -pa->dc_dpphi[C_ALPHA];
  dy[1] = -pa->dc_dpphi[C_BETA];
  return GSL_SUCCESS;
}

/**
 * Function: eccpa_store_sample
 * ----------------------------
 *   Store evolution sample q after an accepted step: state, rates, what reconstruction interpolates
 *   (map coefs of T, V, first-order oscillating parts, averaged phase rates), and pa->rmin_last.
 *
 *   @return  0 on success, 1 if no orbit was found
 */
static int eccpa_store_sample(const eccpa_evo_par *P, double u, const double y[], int order,
                        int q)
{
  PostAdDynamics *pa = P->pa;
  const double J   = P->J0 - u;
  const double r_a = sqrt(y[0]*y[0] + y[1]*y[1]);
  if (eob_eccpa_cons_orbit(pa, J, r_a, P->dyn, PO_CUR)) return 1;
  /* periastron of this sample, for rmin stop; PO_CUR is overwritten below */
  pa->rmin_last = pa->cons_orb[PO_CUR][PO_R_C] - pa->cons_orb[PO_CUR][PO_R_A];

  pa->data[EVO_PPHI][q]  = J;
  pa->data[EVO_ALPHA][q] = y[0];
  pa->data[EVO_BETA][q]  = y[1];
  pa->data[EVO_DADJ][q]  = P->last_dc_dpphi[C_ALPHA];
  pa->data[EVO_DBDJ][q]  = P->last_dc_dpphi[C_BETA];
  pa->data[EVO_E][q]     = pa->cons_orb[PO_CUR][PO_E];
  pa->data[EVO_E_GEO][q] = pa->cons_orb[PO_CUR][PO_E_GEO];
  /* r_circ depends on pphi only; reconstruction interpolates it */
  pa->data[EVO_RCIRC][q] = pa->cons_orb[PO_CUR][PO_R_CIRC];
  pa->data[EVO_FBAR_PPHI][q]  = P->last_FbarJ;
  for (int th = 0; th < TH_N; th++) pa->data[EVO_OMEGA + th][q] = P->last_omega[th];
  pa->data[EVO_ADIAB][q] = P->last_adiab;
  pa->time[q] = 0.0;
  for (int th = 0; th < TH_N; th++) pa->data[EVO_TH + th][q] = pa->data[EVO_FBAR_TH + th][q] = 0.0;
  for (int a = 0; a < OSC_N; a++)
    for (int c = 0; c < 2; c++)
      for (int m = 0; m <= ECCPA_HARMONICS_MAX; m++) pa->evo_osc1_coef[a][c][m][q] = 0.0;
  if (order < 1) {
    /* LO has no oscillating parts, but reconstruction still needs T and V, to get X from
       Lambda and phi from lambda. */
    EccMapCoefs map_coefs;
    if (eccpa_cons_map_coefs_of_orbit(pa, J, r_a, y[0], y[1], PO_CUR, P->last_N, P->dyn, &map_coefs))
      return 1;   /* without T and V sample can't be stored, as at higher order */
    for (int th = 0; th < TH_N; th++)
      for (int c = 0; c < 2; c++)
        for (int m = 0; m <= ECCPA_HARMONICS_MAX; m++) {
          pa->evo_kappa[th][c][m][q]    = map_coefs.kappa[th][c][m];
          pa->evo_kappa_dU[th][c][m][q] = 0.0;   /* osculating = mean at LO */
        }
    pa->evo_M[q] = map_coefs.M;
    return 0;
  }

  /* if oscillating parts can't be built, sample can't be stored and evolution stops */
  if (eob_eccpa_osc1(pa, J, y[0], y[1], 0, P->dyn)) return 1;
  /* oscillating parts are built when map or rate needs them. phase order
     only decides whether their averaged phase rates enter phases. */
  if (EOBPars->postadiabatic_dynamics_ecc_phase >= 1)
    for (int th = 0; th < TH_N; th++) pa->data[EVO_FBAR_TH + th][q] = pa->Fbar_th[th];
  if (P->corr) {
    /* second-order phase rates come from checkpoint corrections. call also
       corrects z rates in pa, so it is outside phase-order test. */
    eob_eccpa_corr_apply(P->corr, pa, J, y[0], y[1]);
    if (EOBPars->postadiabatic_dynamics_ecc_phase >= 2)
      for (int th = 0; th < TH_N; th++) pa->data[EVO_FBAR_TH + th][q] += pa->Fbar2_th[th];
  }
  /* Map coefficients, and their derivative in U (one-sided difference with neighbour at
     U + dU). reconstruction uses it at osculating U. */
  const double idU = (pa->map_coefs_dU_step != 0.0) ? 1.0/pa->map_coefs_dU_step : 0.0;
  for (int th = 0; th < TH_N; th++)
    for (int c = 0; c < 2; c++)
      for (int m = 0; m <= ECCPA_HARMONICS_MAX; m++) {
        pa->evo_kappa[th][c][m][q]    = pa->map_coefs.kappa[th][c][m];
        pa->evo_kappa_dU[th][c][m][q] = (pa->map_coefs_dU.kappa[th][c][m]
                                         - pa->map_coefs.kappa[th][c][m])*idU;
      }
  pa->evo_M[q] = pa->map_coefs.M;
  for (int a = 0; a < OSC_N; a++)
    for (int m = 0; m <= ECCPA_HARMONICS_MAX; m++) {
      pa->evo_osc1_coef[a][0][m][q] = pa->osc1_coef[a].a[m];
      pa->evo_osc1_coef[a][1][m][q] = pa->osc1_coef[a].b[m];
    }

  return 0;
}

/**
 * Function: eob_eccpa_evolve
 * --------------------------
 *   Integrate z in pphi, Eq. (45), until a stopping condition of App. C.
 *
 *   @param[in,out] pa : samples in pa->time and pa->data
 *   @param[in]     J0, alpha0, beta0 : initial mean state
 *   @param[in]     order : order of oscillating parts; from 2 an LO pass places checkpoints
 *   @param[in]     dyn   : Dynamics structure
 *   @return  0 if evolution ran
 *
 *   At map order 3 ICs are redone from pa->ic_*, angle offsets into pa->ic_th_off (ic_valid = 2).
 */
int eob_eccpa_evolve(PostAdDynamics *pa, double J0, double alpha0, double beta0,
                     int order, Dynamics *dyn)
{
  const int max_out = pa->size;
  pa->n_samples = 0;
  if (max_out < 2) return 1;

  const double amax  = EOBPars->postadiabatic_dynamics_ecc_adiab;   /* resolved in eob_set_params */
  const double rmin  = EOBPars->postadiabatic_dynamics_ecc_rmin;

  const double hmax_rel = ECCPA_EVO_STEP_MAX;

  eccpa_evo_par P;
  P.dyn = dyn; P.pa = pa; P.J0 = J0; P.corr = NULL;
  double y[2] = { alpha0, beta0 };
  double u = 0.0, dy0[2];

  /* from order 2: an LO pass places checkpoints, second interpolates them */
  if (order >= 2) {
    const int ran = !eob_eccpa_evolve(pa, J0, alpha0, beta0, 0, dyn);
    const int n1  = pa->n_samples;
    /* At least four samples, which phase quadrature needs. With fewer stage is skipped
       anyway (App. C), so corrections are not worth computing. */
    const int ok1 = ran && n1 >= 4 && !eob_eccpa_corr_build(pa, n1, dyn);
    /* second pass starts with no samples; corrections are already stored. */
    pa->n_samples = 0;
    if (!ok1) return 1;
    P.corr = pa;
    /* Redo osculating -> mean map of initial data at third order, now that the
       checkpoint corrections include third-order map. */
    if (EOBPars->postadiabatic_dynamics_ecc_map >= 3 && pa->corr3_n >= 2 && pa->ic_valid == 1) {
      double Jm, am, bm, th_m[TH_N];
      if (!eob_eccpa_to_mean(pa, pa->ic_J0, pa->ic_ra, pa->ic_X0, pa->ic_phi0, 3, dyn,
                             &Jm, &am, &bm, th_m, pa)) {
        J0 = Jm;  P.J0 = Jm;  y[0] = am;  y[1] = bm;
        for (int th = 0; th < TH_N; th++) pa->ic_th_off[th] = th_m[th];
        pa->ic_valid = 2;
      }
    }
  }

  if (eccpa_z_rhs(0.0, y, dy0, &P) != GSL_SUCCESS || eccpa_store_sample(&P, u, y, order, 0)) return 1;
  pa->n_samples = 1;

  /* first step: fraction of pphi change over one RR time */
  double h = ECCPA_EVO_STEP_FIRST*fabs(P.last_FbarJ)
             *fabs((pa->data[EVO_E][0] - 1.0/EOBPars->nu)/P.last_FbarE);
  if (!(h > 0.0) || h > hmax_rel*J0) h = hmax_rel*J0;
  const double u_end = J0;    /* pphi = 0 is never reached; a stopping condition comes first */

  /* abs tolerance relative to initial r_c */
  if (eob_eccpa_cons_orbit(pa, J0, sqrt(y[0]*y[0] + y[1]*y[1]), dyn, PO_CUR)) return 1;
  const double atol = ECCPA_EVO_ATOL*pa->cons_orb[PO_CUR][PO_R_C];

  gsl_odeiv2_system sys = { eccpa_z_rhs, NULL, 2, &P };
  /* rkf45: few RHS calls, no Jacobian */
  const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_step    *stepper = gsl_odeiv2_step_alloc(T, 2);
  gsl_odeiv2_control *ct = gsl_odeiv2_control_y_new(atol, ECCPA_EVO_RTOL);
  gsl_odeiv2_evolve  *ev = gsl_odeiv2_evolve_alloc(2);

  while (u < u_end) {
    double y_prev[2];
    memcpy(y_prev, y, sizeof y_prev);
    const double u_prev = u;
    const double adiab_prev = pa->data[EVO_ADIAB][pa->n_samples - 1];

    const int status = gsl_odeiv2_evolve_apply(ev, ct, stepper, &sys, &u, u_end, &h, y);
    if (status == GSL_EBADFUNC) {
      /* no bound orbit somewhere in step: stop at last sample */
      memcpy(y, y_prev, sizeof y); u = u_prev;
      break;
    }
    if (status != GSL_SUCCESS) break;
    /* Cap step, so phase quadrature has enough samples to interpolate rates
       between. With default tolerances GSL's own steps are smaller and this doesn't bind. */
    {
      const double hmax = hmax_rel*(P.J0 - u);
      if (h > hmax) h = hmax;
    }

    if (pa->n_samples >= max_out) {
      break;
    }
    if (eccpa_store_sample(&P, u, y, order, pa->n_samples)) {
      memcpy(y, y_prev, sizeof y); u = u_prev;
      break;
    }
    pa->n_samples++;

    const int S = pa->n_samples - 1;
    if (rmin > 0.0 && pa->rmin_last < rmin) break;

    if (pa->data[EVO_ADIAB][S] > amax) {
      /* adiabaticity crossed amax in this step: regula falsi on log(adiab/amax) */
      const double span = u - u_prev;
      double f_lo = log(adiab_prev/amax), f_hi = log(pa->data[EVO_ADIAB][S]/amax);
      double lo = 0.0, hi = 1.0;
      for (int it = 0; it < 5; it++) {
        if (f_hi == f_lo) break;
        double tm = lo - f_lo*(hi - lo)/(f_hi - f_lo);
        if (!(tm > lo && tm < hi)) tm = 0.5*(lo + hi);
        if (eob_eccpa_rates(pa, P.J0 - (u_prev + span*tm),
                            y_prev[0] + tm*(y[0] - y_prev[0]),
                            y_prev[1] + tm*(y[1] - y_prev[1]), 0, dyn)) {
              /* if there is no orbit at this t, move upper end down and try again; it can't
                 be crossing */
              hi = tm;
              f_hi = HUGE_VAL;
              continue;
            }
        const double fm = log(pa->adiab/amax);
        if (pa->adiab > amax) {
          hi = tm;
          f_hi = fm;
        } else {
          lo = tm;
          f_lo = fm;
        }
      }
      const double th = (fabs(f_lo) < fabs(f_hi)) ? lo : hi;
      if (th < ECCPA_CROSS_TMIN) {
        /* crossing is at previous sample, or so close that a sample there would
           duplicate it: previous sample is handoff */
        pa->n_samples--;
      } else if (th < 1.0) {
        /* Integrate again from previous sample to crossing, so last state is a
           real solution and not an interpolation. */
        gsl_odeiv2_driver *dd = gsl_odeiv2_driver_alloc_y_new(&sys, T, 0.5*span*th,
                                                              atol, ECCPA_EVO_RTOL);
        memcpy(y, y_prev, sizeof y);
        u = u_prev;
        const int dstatus = gsl_odeiv2_driver_apply(dd, &u, u_prev + span*th, y);
        gsl_odeiv2_driver_free(dd);
        double dyl[2];
        if (dstatus == GSL_SUCCESS && eccpa_z_rhs(u, y, dyl, &P) == GSL_SUCCESS) {
          /* crossing replaces last sample (this call also fills P->last_* for it) */
          pa->n_samples--;
          if (!eccpa_store_sample(&P, u, y, order, pa->n_samples)) pa->n_samples++;
        }
      }
      break;
    }
  }

  gsl_odeiv2_evolve_free(ev);
  gsl_odeiv2_control_free(ct);
  gsl_odeiv2_step_free(stepper);
  return (pa->n_samples >= 2) ? 0 : 1;
}

/**
 * Function: eob_eccpa_advance
 * ---------------------------
 *   Driver hook (Sec. II C): initial point -> mean elements, evolution, reconstruction, handoff.
 *
 *   @param[in,out] dynp  : Dynamics; advanced to handoff on success
 *   @param[out]    t_adv : time advanced (may be NULL)
 *   @param[out]    n_pa  : samples written (may be NULL)
 *   @return  0 if dynamics were advanced, non-zero if stage declined
 */
int eob_eccpa_advance(Dynamics **dynp, double *t_adv, int *n_pa)
{
  if (t_adv) *t_adv = 0.0;
  if (n_pa)  *n_pa  = 0;
  Dynamics *dyn = *dynp;

  const double r0      = dyn->y[EOB_EVOLVE_RAD];
  const double prstar0 = dyn->y[EOB_EVOLVE_PRSTAR];
  const double phi0    = dyn->y[EOB_EVOLVE_PHI];
  const double J0      = dyn->y[EOB_EVOLVE_PPHI];

  /* skip stage if reconstruction needs more than _nper samples (App. C) */
  if (eob_eccpa_nper(EOBPars->ecc) > EOBPars->postadiabatic_dynamics_ecc_nper) {
    if (VERBOSE)
      printf("Eccentric PA (equinoctial): eccentricity %.4f needs %d samples per orbital period "
             "and nper is %d; skipping.\n", EOBPars->ecc, eob_eccpa_nper(EOBPars->ecc),
             EOBPars->postadiabatic_dynamics_ecc_nper);
    return 1;
  }


  PostAdDynamics *pa = NULL;
  if (PostAdDynamics_alloc(&pa, 1024)) {
    if (VERBOSE) printf("Eccentric PA (equinoctial): out of memory; skipping.\n");
    return 1;
  }
  pa->nmax = ECCPA_NODES_MAX_STAGE;   /* for stage; initial data keep ECCPA_NODES_MAX */

  /* osculating orbit through initial point, from (E, pphi); periastron along +alpha, so
     anomaly is chi. */
  double E0;
  eob_ecc_ham(r0, prstar0, J0, dyn, &E0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  if (eob_eccpa_cons_turning_points(pa, E0, J0, r0, dyn, PO_CUR)) {
    if (VERBOSE)
      printf("Eccentric PA (equinoctial): no bound orbit through r = %.4f, prstar = %.3e at "
             "pphi = %.4f; skipping.\n", r0, prstar0, J0);
    goto fail;
  }
  const double r_a = pa->cons_orb[PO_CUR][PO_R_A];
  double cc = (r_a > 0.0) ? (pa->cons_orb[PO_CUR][PO_R_C] - r0)/r_a : 1.0;
  if (cc >  1.0) cc =  1.0;
  if (cc < -1.0) cc = -1.0;
  /* same orbit from (pphi, r_a). two must agree, and r_c - r_a > 1; otherwise z of a
     wrong orbit would grow without bound. */
  if (eob_eccpa_cons_orbit(pa, J0, r_a, dyn, PO_AUX)
      || !(pa->cons_orb[PO_AUX][PO_R_C] - r_a > 1.0)
      || !(fabs(pa->cons_orb[PO_AUX][PO_R_C] - pa->cons_orb[PO_CUR][PO_R_C]) < 1e-6*fabs(r0))) {
    if (VERBOSE)
      printf("Eccentric PA (equinoctial): no bound orbit through r = %.4f at pphi = %.4f "
             "(half-width %.4g, r_c %.6f against %.6f); skipping.\n", r0, J0, r_a,
             pa->cons_orb[PO_AUX][PO_R_C], pa->cons_orb[PO_CUR][PO_R_C]);
    goto fail;
  }

  /* anomaly of initial point by atan2, Eqs. (39)-(40); Wstar depends on X */
  /* orbit from (pphi, r_a), already solved in PO_AUX above */
  memcpy(pa->cons_orb[PO_CUR], pa->cons_orb[PO_AUX], sizeof pa->cons_orb[PO_CUR]);
  pa->cons_gap[PO_CUR] = pa->cons_gap[PO_AUX];
  double X0 = acos(cc);
  if (prstar0 < 0.0) X0 = 2.0*Pi - X0;
  for (int it = 0; it < 3 && r_a > 0.0; it++) {
    eob_eccpa_cons_sample_point(pa, PO_CUR, r_a, 0.0, X0, 0, dyn);
    if (!(pa->Wstar > 0.0)) break;
    const double sX = prstar0/(r_a*pa->Wstar);
    const double Xn = atan2(sX, cc);
    const double dX = Xn - X0;
    X0 = (Xn < 0.0) ? Xn + 2.0*Pi : Xn;
    if (fabs(dX) < 1e-15) break;
  }

  /* osculating -> mean at map order of level; mean angles are phase offsets */
  double J0m, alpha0, beta0, th_off[TH_N];
  if (eob_eccpa_to_mean(pa, J0, r_a, X0, phi0, EOBPars->postadiabatic_dynamics_ecc_map, dyn,
                        &J0m, &alpha0, &beta0, th_off, NULL)) {
    if (VERBOSE) printf("Eccentric PA (equinoctial): the osculating -> mean conversion of the "
                        "initial data failed; skipping.\n");
    goto fail;
  }

  /* osculating initial state, to redo conversion at third order once checkpoint
     corrections exist (eob_eccpa_evolve) */
  pa->ic_J0 = J0;  pa->ic_ra = r_a;  pa->ic_X0 = X0;  pa->ic_phi0 = phi0;  pa->ic_valid = 1;

  /* evolution and phases */
  PROF_START(PROF_ECCPA_SEC);
  int rc = eob_eccpa_evolve(pa, J0m, alpha0, beta0,
                            EOBPars->postadiabatic_dynamics_ecc_osc_order, dyn);
  int n = pa->n_samples;
  if (!rc) rc = eob_eccpa_phases(pa, n, 0, dyn);
  PROF_STOP(PROF_ECCPA_SEC);
  if (rc || n < 4) {
    if (VERBOSE) printf("Eccentric PA (equinoctial): evolution did not run, skipping.\n");
    goto fail;
  }

  if (pa->ic_valid == 2)
    for (int th = 0; th < TH_N; th++) th_off[th] = pa->ic_th_off[th];

  /* add offsets, initial time and mean angles of initial data */
  {
    const double t_off = dyn->t;
    for (int i = 0; i < n; i++) {
      pa->time[i] += t_off;
      for (int th = 0; th < TH_N; th++) pa->data[EVO_TH + th][i] += th_off[th];
    }
    for (int k = 0; k < pa->fine_n; k++) {     /* and fine grid, which reconstruction reads */
      pa->fine_t[k] += t_off;
      for (int th = 0; th < TH_N; th++) pa->fine_th[th][k] += th_off[th];
    }
  }

  /* evolution steps in pphi and gets t only from phases, so it can go past t_stop. Drop
     those samples; ODE starts from last one before t_stop. */
  while (n > 4 && pa->time[n-1] > dyn->t_stop) n--;
  pa->n_samples = n;
  if (pa->time[n-1] > dyn->t_stop) {
    /* skip stage if even first samples are past t_stop */
    if (VERBOSE) printf("Eccentric PA (equinoctial): ode_tmax is inside the first few secular "
                        "steps, skipping.\n");
    goto fail;
  }

  const double span = pa->time[n-1] - pa->time[0];
  const double P_r0 = 2.0*Pi/pa->data[EVO_OMEGA + TH_LAMBDA_MEAN][0];
  if (!(span > 2.0*P_r0)) {
    if (VERBOSE) printf("Eccentric PA (equinoctial): nothing to gain (%.1f radial periods), "
                        "skipping.\n", span/P_r0);
    goto fail;
  }

  /* DEBUG: reconstruct only last export_window SPA widths before handoff */
  double t_rec0 = pa->time[0];
  if (DEBUG && EOBPars->postadiabatic_dynamics_ecc_export && EOBPars->postadiabatic_dynamics_ecc_export_window > 0.0) {
    const int kf = find_point_bisection(pa->time[n-1], pa->fine_n, pa->fine_t, 1);
    const double psdd = 2.0*(pa->fine_dth_dpphi[TH_AZIMUTH][kf+1]/pa->fine_dt_dpphi[kf+1]
                             - pa->fine_dth_dpphi[TH_AZIMUTH][kf]/pa->fine_dt_dpphi[kf])/(pa->fine_t[kf+1] - pa->fine_t[kf]);
    if (psdd > 0.0)
      t_rec0 = pa->time[n-1] - EOBPars->postadiabatic_dynamics_ecc_export_window*sqrt(2.0*Pi/psdd);
    if (t_rec0 < pa->time[0]) t_rec0 = pa->time[0];
  }
  /* reconstruct with step = fixed fraction of orbital period */
  int nrec = 0;
  const int nper = EOBPars->postadiabatic_dynamics_ecc_nper;
  if (nper > 1) {
    /* Samples per orbital period, set once by most eccentric sample (usually first), at
       most nper. */
    int nper_use = ECCPA_NPER_MIN;
    for (int q = 0; q < n; q++) {
      const int nq = eob_eccpa_nper(pa->data[EVO_E_GEO][q]);
      if (nq > nper_use) nper_use = nq;
    }
    if (nper_use > nper) nper_use = nper;
    if (nper_use < 2) nper_use = 2;

    /* Omega grows, so step is smallest at handoff, which bounds number of samples */
    const int want = (int)((pa->time[n-1] - t_rec0)*pa->data[EVO_OMEGA + TH_AZIMUTH][n-1]/(2.0*Pi)*nper_use) + 8;
    if (want > dyn->size) {
      Dynamics_push(dynp, want + 16);
      dyn = *dynp;
    }

    PROF_START(PROF_ECCPA_RECON);
    const int store = dyn->store, noflx = dyn->noflx;
    dyn->store = 1; dyn->noflx = 1;
    double t = t_rec0;
    pa->rec_nseed = 0;   /* no anomaly seed from earlier calls */
    int k = 0;
    for (; t < pa->time[n-1]; k++) {
      /* want is only an estimate: phi oscillates around its mean by O(e) */
      if (k == dyn->size) {
        Dynamics_push(dynp, dyn->size + want/4 + 16);
        dyn = *dynp;
      }
      if (eob_eccpa_recon(pa, n, t, dyn)) break;
      /* into dyn->data only: dyn->y is initial state until handoff, since a stage that
         declines below leaves ODE to start from it */
      /* other quantities from point just computed, without calling right-hand side */
      dyn->time[k]             = t;
      dyn->data[EOB_RAD][k]    = pa->r;
      dyn->data[EOB_PHI][k]    = pa->phi;
      dyn->data[EOB_PPHI][k]   = pa->pphi;
      dyn->data[EOB_PRSTAR][k] = pa->prstar;
      dyn->data[EOB_MOMG][k]   = pa->Omg;
      dyn->data[EOB_DDOTR][k]  = pa->ddotr;
      /* orbital frequency without spin-orbit terms of dH/dpphi, as ODE stores it and
         merger peak search reads it */
      {
        const double ucr = (pa->rc != 0.0) ? 1.0/pa->rc : 0.0;
        dyn->data[EOB_OMGORB][k] = (pa->H != 0.0 && pa->Heff_orb != 0.0)
          ? (pa->pphi * pa->A * ucr * ucr) / (EOBPars->nu * pa->H * pa->Heff_orb) : pa->Omg;
      }
      dyn->data[EOB_E0][k]     = EOBPars->nu * pa->H;
      /* step by a fixed fraction of current orbital period */
      t += (pa->Omg > 0.0) ? 2.0*Pi/(pa->Omg*nper_use) : span/want;
    }
    nrec = k;
    dyn->store = store; dyn->noflx = noflx;
    PROF_STOP(PROF_ECCPA_RECON);
    if (t < pa->time[n-1]) {
      /* a sample failed: a gap before handoff would be bridged silently, so decline */
      if (VERBOSE) printf("Eccentric PA (equinoctial): reconstruction failed at t = %g; skipping.\n", t);
      goto fail;
    }
    if (nrec < 2) nrec = 0;
  }

  /* handoff state. */
  if (eob_eccpa_recon(pa, n, pa->time[n-1], dyn)) {
    if (VERBOSE) printf("Eccentric PA (equinoctial): handoff reconstruction failed.\n");
    goto fail;
  }

  EOBPars->pa_tab3_n    = pa->corr3_n;
  dyn->t      = pa->t_at;
  dyn->r      = pa->r;
  dyn->phi    = pa->phi;
  dyn->prstar = pa->prstar;
  dyn->pphi   = pa->pphi;
  dyn->y[EOB_EVOLVE_RAD]    = pa->r;
  dyn->y[EOB_EVOLVE_PHI]    = pa->phi;
  dyn->y[EOB_EVOLVE_PRSTAR] = pa->prstar;
  dyn->y[EOB_EVOLVE_PPHI]   = pa->pphi;

  /* one call to right-hand side, to fill quantities waveform reads from dyn */
  {
    const int store = dyn->store, noflx = dyn->noflx;
    dyn->store = 1; dyn->noflx = 1;
    p_eob_dyn_rhs(dyn->t, dyn->y, dyn->dy, dyn);
    dyn->store = store; dyn->noflx = noflx;
  }

  dyn->y0[EOB_ID_RAD]    = pa->r;
  dyn->y0[EOB_ID_PHI]    = pa->phi;
  dyn->y0[EOB_ID_PPHI]   = pa->pphi;
  dyn->y0[EOB_ID_PRSTAR] = pa->prstar;
  dyn->y0[EOB_ID_OMGJ]   = dyn->Omg;
  dyn->y0[EOB_ID_E0]     = dyn->E;

  /* DEBUG builds only: evolution, its fine phase grid and, with export_nodes, modes on
     nodes of each mean orbit, written to output_dir for frequency-domain route. */
#if (DEBUG)
  if (EOBPars->postadiabatic_dynamics_ecc_export) {
    char fname[STRLEN];
    FILE *fp;
    sprintf(fname, "%s/eccpa_secular.txt", EOBPars->output_dir);
    if ((fp = fopen(fname, "w"))) {
      fprintf(fp, "#");
      for (int v = 0; v <= EVO_NVARS; v++) fprintf(fp, " %s", pa_sec_var[v]);
      fprintf(fp, "\n");
      for (int q = 0; q < n; q++) {
        fprintf(fp, "%.16e", pa->time[q]);
        for (int v = 0; v < EVO_NVARS; v++) fprintf(fp, " %.16e", pa->data[v][q]);
        fprintf(fp, "\n");
      }
      fclose(fp);
    }
    const double *fcol[PAF_NVARS] = {
      [PAF_PPHI]     = pa->fine_pphi,     [PAF_T]        = pa->fine_t,
      [PAF_LAM_RAD]  = pa->fine_th[TH_LAMBDA_MEAN], [PAF_LAM_AZI] = pa->fine_th[TH_AZIMUTH],
      [PAF_DT]       = pa->fine_dt_dpphi,
      [PAF_DLAM_RAD] = pa->fine_dth_dpphi[TH_LAMBDA_MEAN], [PAF_DLAM_AZI] = pa->fine_dth_dpphi[TH_AZIMUTH],
    };
    sprintf(fname, "%s/eccpa_phase.txt", EOBPars->output_dir);
    if ((fp = fopen(fname, "w"))) {
      fprintf(fp, "#");
      for (int v = 0; v < PAF_NVARS; v++) fprintf(fp, " %s", pa_fine_var[v]);
      fprintf(fp, "\n");
      for (int k = 0; k < pa->fine_n; k++) {
        for (int v = 0; v < PAF_NVARS; v++) fprintf(fp, "%s%.16e", v ? " " : "", fcol[v][k]);
        fprintf(fp, "\n");
      }
      fclose(fp);
    }

    /* modes on nodes of each mean orb and 4 neighbours, for FD route */
    const int Nw = EOBPars->postadiabatic_dynamics_ecc_export_nodes;
    if (Nw > 0) {
      if (Nw > ECCPA_NODES_DIM) {
        fprintf(stderr, "postadiabatic_dynamics_ecc_export_nodes is above ECCPA_NODES_DIM, nodes not written\n");
        goto export_done;
      }
      const int nm = EOBPars->use_mode_lm_size, ncol = PAN_NVARS + 2*nm;
      const long nqi = (long) n*Nw;
      double *out = (double *) malloc(5*ncol*nqi*sizeof(double));
      if (!out) goto export_done;
      Dynamics *dn;
      Dynamics_alloc(&dn, 1, "pa_node");
      dn->store = dn->noflx = 1;
      Waveform_lm_t hn;
      EccMapCoefs map_coefs_n;
      double y[EOB_EVOLVE_NVARS], dy[EOB_EVOLVE_NVARS];
      const int gcol[5][2] = { {PAN_YJ, OSC_PPHI}, {PAN_YZR, OSC_Z_RE}, {PAN_YZI, OSC_Z_IM},
                               {PAN_V1, OSC_TH + TH_LAMBDA_MEAN}, {PAN_W1, OSC_TH + TH_AZIMUTH} };
      for (int q = 0; q < n; q++)
        for (int ob = 0; ob < 5; ob++) {
          const double sa = 1.0 + 1e-4*((ob == 1) - (ob == 2));
          const double J  = pa->data[EVO_PPHI][q]*(1.0 + 1e-4*((ob == 3) - (ob == 4)));
          const double a  = pa->data[EVO_ALPHA][q]*sa, b = pa->data[EVO_BETA][q]*sa;
          const double ra = sqrt(a*a + b*b), w = eccpa_peri_angle(a, b);
          const int bad = eob_eccpa_cons_orbit(pa, J, ra, dn, PO_CUR)
                       || eob_eccpa_cons_sample_nodes(pa, PO_CUR, a, b, Nw, 0, dn, NULL, NULL) < 0
                       || eccpa_cons_map_coefs(pa, PO_CUR, pa->node[PA_DTDX], pa->node[PA_OMG], Nw, dn, &map_coefs_n);
          for (int i = 0; i < Nw; i++) {
            double *o = out + ob*ncol*nqi + (long) q*Nw + i;
            if (bad) {
              for (int v = 0; v < ncol; v++) o[v*nqi] = NAN;
              continue;
            }
            const double X = pa->node[PA_X][i];
            y[EOB_EVOLVE_RAD]    = pa->node[PA_R][i];
            y[EOB_EVOLVE_PHI]    = 0.0;
            y[EOB_EVOLVE_PRSTAR] = pa->node[PA_PRSTAR][i];
            y[EOB_EVOLVE_PPHI]   = J;
            p_eob_dyn_rhs(pa->time[q], y, dy, dn);
            eob_wav_hlm(dn, &hn);
            o[PAN_X*nqi]       = X;
            o[PAN_LAM*nqi]     = X + eccpa_ksum(map_coefs_n.kappa[TH_LAMBDA_MEAN][0],
                                                map_coefs_n.kappa[TH_LAMBDA_MEAN][1],
                                                map_coefs_n.M, ra, X - w, NULL);
            /* V = -M_lambda */
            o[PAN_V*nqi]       = -eccpa_ksum(map_coefs_n.kappa[TH_AZIMUTH][0],
                                             map_coefs_n.kappa[TH_AZIMUTH][1],
                                             map_coefs_n.M, ra, X - w, NULL);
            o[PAN_DTDX*nqi]    = pa->node[PA_DTDX][i];
            o[PAN_R*nqi]       = pa->node[PA_R][i];
            o[PAN_PRSTAR*nqi]  = pa->node[PA_PRSTAR][i];
            o[PAN_OMG*nqi]     = dn->Omg;
            o[PAN_J*nqi]       = J;
            o[PAN_RA*nqi]      = ra;
            o[PAN_RDOT*nqi]    = dn->rdot;
            o[PAN_R2DOT*nqi]   = dn->r2dot;
            o[PAN_R3DOT*nqi]   = dn->r3dot;
            o[PAN_OMGDOT*nqi]  = dn->Omegadot;
            o[PAN_OMG2DOT*nqi] = dn->Omega2dot;
            o[PAN_PRSDOT*nqi]  = dn->prsdot;
            /* sample's first-order oscillating parts, summed as in reconstruction */
            double ck[ECCPA_HARMONICS_MAX+1], sk[ECCPA_HARMONICS_MAX+1];
            eccpa_trig(X - w, pa->evo_M[q], ck, sk);
            for (int g = 0; g < 5; g++) {
              double acc = pa->evo_osc1_coef[gcol[g][1]][0][0][q];
              for (int k = 1; k <= pa->evo_M[q]; k++)
                acc += pa->evo_osc1_coef[gcol[g][1]][0][k][q]*ck[k] + pa->evo_osc1_coef[gcol[g][1]][1][k][q]*sk[k];
              o[gcol[g][0]*nqi] = acc;
            }
            for (int j = 0; j < nm; j++) {
              const int k = EOBPars->use_mode_lm[j];
              o[(PAN_NVARS + 2*j)*nqi]     =  hn.ampli[k]*cos(hn.phase[k]);
              o[(PAN_NVARS + 2*j + 1)*nqi] = -hn.ampli[k]*sin(hn.phase[k]);
            }
          }
        }
      Dynamics_free(dn);
      /* one row per (orbit, sample, node); orbits are mean one, r_a up and down, pphi up
         and down; modes as real and imaginary parts, in order of use_mode_lm */
      sprintf(fname, "%s/eccpa_nodes.txt", EOBPars->output_dir);
      if ((fp = fopen(fname, "w"))) {
        fprintf(fp, "# orbit sample node");
        for (int v = 0; v < PAN_NVARS; v++) fprintf(fp, " %s", pa_node_var[v]);
        for (int j = 0; j < nm; j++) fprintf(fp, " re_hlm_%d im_hlm_%d", EOBPars->use_mode_lm[j], EOBPars->use_mode_lm[j]);
        fprintf(fp, "\n");
        for (int ob = 0; ob < 5; ob++)
          for (long qi = 0; qi < nqi; qi++) {
            fprintf(fp, "%d %ld %ld", ob, qi/Nw, qi%Nw);
            for (int v = 0; v < ncol; v++) fprintf(fp, " %.16e", out[(ob*ncol + v)*nqi + qi]);
            fprintf(fp, "\n");
          }
        fclose(fp);
      }
      free(out);
    }
  }
 export_done:
#endif

  if (VERBOSE) {
    PRSECTN("Eccentric post-adiabatic dynamics (equinoctial)");
    printf("%-40s = %s (map %d, rate %d, phase %d)\n","level",
           eccpa_level_opt[EOBPars->postadiabatic_dynamics_ecc_level],
           EOBPars->postadiabatic_dynamics_ecc_map, EOBPars->postadiabatic_dynamics_ecc_rate,
           EOBPars->postadiabatic_dynamics_ecc_phase);
    printf("%-40s = %d\n","evolution steps",n);
    printf("%-40s = %.6e (%.0f radial orbits)\n","time advanced",span,span/P_r0);
    printf("%-40s = %.6e -> %.6e\n","eccentricity (r_a/r_c)",pa->data[EVO_E_GEO][0],pa->data[EVO_E_GEO][n-1]);
    printf("%-40s = %.6e\n","adiabaticity at handoff",pa->data[EVO_ADIAB][n-1]);
    printf("%-40s = %d\n","reconstructed samples",nrec);
    if (EOBPars->postadiabatic_dynamics_ecc_map >= 3)
      printf("%-40s = %d checkpoints%s\n","third-order map",pa->corr3_n,
             (pa->ic_valid == 2) ? ", initial data at third order" : ", NOT built: second order");
  }

  if (t_adv) *t_adv = span;
  if (n_pa)  *n_pa  = nrec;
  PostAdDynamics_free(&pa);
  return 0;

fail:
  PostAdDynamics_free(&pa);
  return 1;
}

/**
 * Function: eccpa_ic_orbit
 * ------------------------
 *   Mean orbit of initial data, Eq. (48): r_a = e r_c and ecc_freq frequency = omg_tgt.
 *
 *   @param[in,out] pa       : orbit left in PO_CUR
 *   @param[in]     e        : eccentricity r_a/r_c
 *   @param[in]     omg_tgt  : target frequency
 *   @param[in]     J_seed, ra_seed : first guess
 *   @param[in]     dyn      : Dynamics structure
 *   @param[out]    J_out, ra_out   : solution
 *   @return  ICEQ_OK, or reason it failed
 */
static int eccpa_ic_orbit(PostAdDynamics *pa, double e, double omg_tgt,
                          double J_seed, double ra_seed, Dynamics *dyn,
                          double *J_out, double *ra_out)
{
  double J = J_seed, r_a = ra_seed;
  double Jp = 0.0, fp = 0.0;
  double bestf = 1e300, bestJ = J, bestra = r_a;
  int matched = 0;
  for (int it = 0; !matched && it < 48; it++) {
    for (int is = 0; is < 48; is++) {
      if (eob_eccpa_cons_orbit(pa, J, r_a, dyn, PO_CUR)) return ICEQ_ORBIT;
      const double dra = e*pa->cons_orb[PO_CUR][PO_R_C] - r_a;
      r_a += dra;
      if (fabs(dra) <= 1e-14*(1.0 + r_a)) break;
    }
    /* frequency per ecc_freq, as in adiabatic ICs */
    double omg;
    if (EOBPars->ecc_freq == ECCFREQ_ORBAVGD) {
      /* <Omega> on nodes of this orbit, as eccpa_cons_map_coefs computes it */
      const int N = eccpa_nnodes(pa, pa->cons_orb[PO_CUR][PO_E_GEO], 0);
      double n_rad;
      if (eob_eccpa_cons_sample_nodes(pa, PO_CUR, r_a, 0.0, N, 0, dyn, NULL, NULL) < 0
          || eccpa_mean_freqs(pa->node[PA_DTDX], pa->node[PA_OMG], N, &n_rad, &omg))
        return ICEQ_ORBIT;
    } else if (EOBPars->ecc_freq == ECCFREQ_AVERAGE) {
      eob_eccpa_cons_sample_point(pa, PO_CUR, r_a, 0.0, 0.0, 0, dyn);
      const double omg_peri = pa->Omg;
      eob_eccpa_cons_sample_point(pa, PO_CUR, r_a, 0.0, Pi, 0, dyn);
      omg = 0.5*(omg_peri + pa->Omg);
    } else {
      const double Xt = (EOBPars->ecc_freq == ECCFREQ_PERIASTRON) ? 0.0 : Pi;
      eob_eccpa_cons_sample_point(pa, PO_CUR, r_a, 0.0, Xt, 0, dyn);
      omg = pa->Omg;
    }
    const double f = omg - omg_tgt;
    if (fabs(f) < fabs(bestf)) { bestf = f; bestJ = J; bestra = r_a; }
    if (fabs(f) <= ECCPA_IC_FREQ_TOL*omg_tgt) { matched = 1; break; }
    double dJ;
    if (it == 0) {
      dJ = f*J/(3.0*omg_tgt);
    } else {
      const double den = f - fp;
      if (!(fabs(den) > 0.0)) break;
      dJ = -f*(J - Jp)/den;
    }
    Jp = J; fp = f;
    J += dJ;
    if (!(J > 0.0)) return ICEQ_FREQ;
    /* stop if there is no progress, and let test below decide */
    if (it > 2 && fabs(dJ) <= 1e-14*J) break;
  }

  /* accept best iterate within ECCPA_IC_FREQ_TOL_FLOOR (freq noise at high e) */
  if (!matched && fabs(bestf) <= ECCPA_IC_FREQ_TOL_FLOOR*omg_tgt) {
    J = bestJ; r_a = bestra;
    for (int is = 0; is < 48; is++) {
      if (eob_eccpa_cons_orbit(pa, J, r_a, dyn, PO_CUR)) return ICEQ_ORBIT;
      const double dra = e*pa->cons_orb[PO_CUR][PO_R_C] - r_a;
      r_a += dra;
      if (fabs(dra) <= 1e-14*(1.0 + r_a)) break;
    }
    matched = 1;
  }
  *J_out = J; *ra_out = r_a;
  return matched ? ICEQ_OK : ICEQ_FREQ;
}

/**
 * Function: eob_eccpa_ic_state
 * ----------------------------
 *   Initial phase-space point for e and omg_tgt (Sec. III B): mean orbit, Eq. (48), oscillating
 *   parts at X, Eq. (49), osculating elements, Eq. (51).
 *
 *   @param[in]  e       : eccentricity r_a/r_c
 *   @param[in]  omg_tgt : target orbital frequency, per ecc_freq
 *   @param[in]  X       : anomaly from periastron
 *   @param[in]  order   : map order (max 2)
 *   @param[in]  dyn     : Dynamics structure
 *   @param[out] out     : state, and reason on failure
 *   @return  0 on success, 1 on failure
 */
int eob_eccpa_ic_state(double e, double omg_tgt, double X, int order,
                       Dynamics *dyn, EccIcEq *out)
{
  memset(out, 0, sizeof *out);
  const int ok = (e < 1.0) && (e >= 0.0);
  int reason = ok ? ICEQ_OK : ICEQ_INPUT;

  PostAdDynamics *pa = NULL;
  if (PostAdDynamics_alloc(&pa, 8)) {   /* no evolution here, so 8 samples are plenty */
    out->reason = ICEQ_MEMORY;
    return 1;
  }

  out->order = order;
  if (!ok) goto out;

  /* seeds from Kepler's third law */
  double J = 0.0, r_a = 0.0;
  reason = eccpa_ic_orbit(pa, e, omg_tgt, eob_dyn_j0(pow(omg_tgt, -2./3.)*(1. - e*e), dyn),
                          e*pow(omg_tgt, -2./3.), dyn, &J, &r_a);
  if (reason != ICEQ_OK) goto out;

  /* Eq. (51), mean periastron along +alpha */
  double g[OSC_N], osc2[OSC_N];
  for (int a = 0; a < OSC_N; a++) g[a] = osc2[a] = 0.0;
  /* mean orbit's r_c, before PO_CUR is overwritten by osculating orbit */
  const double rc_m = pa->cons_orb[PO_CUR][PO_R_C];
  EccMapCoefs map_coefs_mean;
  memset(&map_coefs_mean, 0, sizeof map_coefs_mean);
  if (order >= 1) {
    if (eob_eccpa_osc1(pa, J, r_a, 0.0, 0, dyn)) { reason = ICEQ_GEN1; goto out; }
    map_coefs_mean = pa->map_coefs;
    eccpa_osc1_at(pa, X, g);
  }
  if (order >= 2) {
    /* fail rather than deliver a lower order than asked */
    if (eob_eccpa_rates_higher(pa, J, r_a, 0.0, 0, dyn)) { reason = ICEQ_GEN2; goto out; }
    /* lambda~(2) is not needed: phi = 0 at start */
    const int Nn = pa->nnodes;
    osc2[OSC_PPHI] = eccpa_osc2_col(pa->mean_node[MEAN_NODE_CTILDE2 + C_PPHI],  Nn, X);
    osc2[OSC_Z_RE] = eccpa_osc2_col(pa->mean_node[MEAN_NODE_CTILDE2 + C_ALPHA], Nn, X);
    osc2[OSC_Z_IM] = eccpa_osc2_col(pa->mean_node[MEAN_NODE_CTILDE2 + C_BETA],  Nn, X);
    osc2[OSC_TH + TH_LAMBDA_MEAN] = eccpa_osc2_col(pa->mean_node[MEAN_NODE_THTILDE2 + TH_LAMBDA_MEAN], Nn, X);
  }

  /* at order 0 there is nothing to invert, Lambda = X */
  double T0 = 0.0;
  double LamP = X;
  if (order >= 1) {
    eccpa_Tmap(&map_coefs_mean, r_a, X, &T0, NULL);
    LamP = X + T0 + g[OSC_TH + TH_LAMBDA_MEAN] + osc2[OSC_TH + TH_LAMBDA_MEAN];
  }
  const double pphi  = J   + g[OSC_PPHI] + osc2[OSC_PPHI];
  const double alpha = r_a + g[OSC_Z_RE] + osc2[OSC_Z_RE];
  const double beta  =       g[OSC_Z_IM] + osc2[OSC_Z_IM];

  /* X' solves X' + T(X') = Lambda on osculating orbit. phi = 0 at start, so lambda~ is
     not needed. */
  const double ra_o = sqrt(alpha*alpha + beta*beta);
  const double w_o  = eccpa_peri_angle(alpha, beta);
  /* at order 0 osculating orbit is mean one, already in PO_CUR, and X is unchanged */
  double Xp = X;
  if (order >= 1) {
    EccMapCoefs map_coefs_osc;
    if (eob_eccpa_cons_orbit(pa, pphi, ra_o, dyn, PO_CUR)) { reason = ICEQ_OSC_ORBIT; goto out; }
    /* on nodes of oscillating parts, where mean orbit's map came from; another node
       count would leave an error in X' as large as second order */
    const int N = (pa->osc1_nnodes > 0) ? pa->osc1_nnodes
                                        : eccpa_nnodes(pa, pa->cons_orb[PO_CUR][PO_E_GEO], 0);
    if (eccpa_cons_map_coefs_of_orbit(pa, pphi, ra_o, alpha, beta, PO_CUR, N, dyn, &map_coefs_osc)) {
      reason = ICEQ_OSC_MAP;
      goto out;
    }
    Xp = eccpa_X_of_Lam(&map_coefs_osc, ra_o, w_o, LamP, LamP - T0);
  }
  eob_eccpa_cons_sample_point(pa, PO_CUR, alpha, beta, Xp, 0, dyn);
  /* r between osculating turning points, Xp close to X (else wrong branch) */
  const double rc_o = pa->cons_orb[PO_CUR][PO_R_C], ra_osc = pa->cons_orb[PO_CUR][PO_R_A];
  if (!(fabs(pa->r - rc_o) <= ra_osc*(1.0 + 1e-6) + 1e-9)) reason = ICEQ_OFF_ORBIT;
  if (!(fabs(Xp - X) < 0.5)) reason = ICEQ_BRANCH;
  if (reason != ICEQ_OK) goto out;

  out->r        = pa->r;
  out->prstar   = pa->prstar;
  out->sqrtAbyB = pa->sqrtAbyB;
  out->pphi     = pphi;
  out->H        = pa->H;
  out->Omg      = pa->Omg;
  out->J        = J;
  out->r_a      = r_a;
  out->e_geo    = (rc_m > 0.0) ? r_a/rc_m : 0.0;
  out->X        = X;
  out->Xp       = Xp;

out:
  out->reason = reason;
  PostAdDynamics_free(&pa);
  return (reason == ICEQ_OK) ? 0 : 1;
}
