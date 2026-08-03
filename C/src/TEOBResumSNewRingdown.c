/** \file TEOBResumSNewRingdown.c
 *  \brief New multi-mode ringdown + NQC model (EOBPars->ringdown_model ==
 *  RINGDOWN_A22, user-facing option string "new_A22" -- see
 *  TEOBResumSPars.c's function-pointer dispatch and TEOBResumSWrap.c's
 *  "ringdown_model" option). Supports (2,1),(2,2),
 *  (3,2),(3,3),(4,3),(4,4),(5,5).
 *
 *  Everything here is anchored at ONE shared attachment point t0, located
 *  EXACTLY as in the old model's tmrg[k22] (true-(2,2)-amplitude-peak
 *  proxy): the refined (7-point local polynomial) peak of the orbital
 *  frequency Omega_orb, minus eob_nqc_timeshift(nu,chi1), minus an
 *  additional NR-calibrated "-2M" (see eob_A22_find_t0_raw) -- NOT a fresh
 *  peak-search on the (2,2) amplitude array itself (tried that first; it
 *  read uninitialized memory in the ringdown-extended tail of the
 *  waveform array and is unnecessary -- the old model's own
 *  Omega-peak+shift proxy is already the right answer here). ALL modes
 *  attach at this SAME t0 -- "usual" for this kind of peak-anchored
 *  template, no per-mode dtmrg offset like the old model's HM path.
 *
 *  (3,2) and (4,3) are not spherical-harmonic modes of a single QNM --
 *  the postpeak template below is intrinsically SPHEROIDAL (Kerr QNMs are
 *  eigenfunctions of the spin-weighted spheroidal, not spherical,
 *  equation), so (3,2)/(4,3) are built by attaching spheroidal S32/S43
 *  templates and MIXING them with the already-attached spherical Y22/Y33
 *  via Berti-Klein (arXiv:1408.1860) coefficients mu_{m,l,l',n=0}, same
 *  formalism as eob_ringdown/py/spheroidal.py's own
 *  spherical_from_spheroidal (l'<=l truncation) -- see A22_attach_mixed_
 *  mode(). The (2,1),(2,2),(3,3),(4,4),(5,5) modes need no mixing (l has
 *  no lower partner at fixed m, or m=l is the leading/only term).
 *
 *  PLACEHOLDER STATE (as of this writing): the (nu,spin) global-fit
 *  amplitude/phase template surfaces (hatA0/hatdA0/hatd2A0/omg0/domg0/
 *  cAmp0/cAmp1/cphf0/cphf1) exist for (2,2) (l2m2 functions below) and
 *  (3,2) (l3m2 functions below -- calibrated on the SPHEROIDAL S32
 *  decomposition, NOT the spherical (3,2) mode, see A22_get_fit_set's own
 *  header comment), both from eob_ringdown/py/fitter.py's fit_2d_to_c;
 *  every other mode ((2,1),(3,3),(4,3),(4,4),(5,5)) uses
 *  A22_fit_placeholder (returns 0 -- see A22_get_fit_set), so ITS OWN
 *  postpeak template is identically zero (smoothly faded out by NQC
 *  approaching t0, then exactly 0 for t>=t0) until real fits are
 *  generated and swapped in. The Berti-Klein mixing coefficients (A22_mu)
 *  are real (18 hardcoded (m,l,l') combinations, both signs of m -- see
 *  A22_mu's own header comment), independent of the amplitude/phase
 *  placeholder state above. QNM complex frequencies (sigmar,sigmai) are
 *  NOT placeholder -- Kerr perturbation theory, not part of the (nu,spin)
 *  fits -- (2,2) keeps its own long-validated hardcoded fit (QNM_coefs'
 *  own (2,2) entry is a DIFFERENT, close but not identical fit --
 *  swapping would be a silent regression for the one mode extensively
 *  tested so far); every other mode uses QNM_coefs (already implemented
 *  in this repo, covers exactly this same mode set).
 *
 *  The nine hatA0_l2m2/hatdA0_l2m2/hatd2A0_l2m2/omg0_l2m2/domg0_l2m2/
 *  cAmp_l2m2_0/cAmp_l2m2_1/cphf_l2m2_0/cphf_l2m2_1 functions below are
 *  generated verbatim by eob_ringdown/py/fitter.py's fit_2d_to_c (see
 *  eob_ringdown/py/generic_q/gfits_A22_v1/l2m2.c) -- x=a0, y=1-X12 is
 *  fit_2d's own boundary convention, NOT raw X12.
 */

#include "TEOBResumS.h"

/* === BEGIN fit_2d_to_c: hatA0_l2m2 === */
double hatA0_l2m2(double x, double y) {
    double x2 = x*x;
    double y2 = y*y;
    double y3 = y2*y;
    return ((1.44959188382681914e+00)*((1.0 + (6.89849337014823316e-01)*((-5.75261729800441213e-01)/(1.0 + 3.61775910818433410e-01*y))*x + (6.89849337014823316e-01)*((-9.37280958978483902e-02)/(1.0 + 9.09090746955315998e-01*y))*x2)/(1.0 + ((-7.71842771502132208e-01)/(1.0 + 6.05283452507252129e-01*y))*x))*(1.0 + (6.89849337014823316e-01)*(-2.48881225900474852e-02 + -5.59213712513535674e-02*x)*y + (6.89849337014823316e-01)*(5.14933452191105534e-01 + -9.42615659347894180e-02*x)*y2 + (6.89849337014823316e-01)*(-3.72149358035482214e-01 + 1.36910594043419293e-01*x)*y3));
}
/* === END fit_2d_to_c: hatA0_l2m2 === */

/* === BEGIN fit_2d_to_c: hatd2A0_l2m2 === */
double hatd2A0_l2m2(double x, double y) {
    double x2 = x*x;
    double y2 = y*y;
    return ((-2.43607606729431304e-03)*((1.0 + (-4.10496212916156708e+02)*(2.00729089641414572e-03 + -2.09867228908063551e-03*y)*x + (-4.10496212916156708e+02)*(7.32092049813945337e-04 + -4.41753499064325968e-04*y)*x2)/(1.0 + (-4.38300163515166097e-01 + 2.42984252628984904e-01*y)*x))*((1.0 + (-4.10496212916156708e+02)*(-7.53541473064639738e-03 + -3.17115963318180245e-04*x)*y + (-4.10496212916156708e+02)*(5.41188839082077538e-04 + 4.22940242772500748e-03*x)*y2)/(1.0 + (1.32873215103040798e+00 + 5.33370232445838588e-01*x)*y)));
}
/* === END fit_2d_to_c: hatd2A0_l2m2 === */

/* === BEGIN fit_2d_to_c: hatdA0_l2m2 === */
double hatdA0_l2m2(double x, double y) {
    return 0.0;
}
/* === END fit_2d_to_c: hatdA0_l2m2 === */

/* === BEGIN fit_2d_to_c: omg0_l2m2 === */
double omg0_l2m2(double x, double y) {
    double x2 = x*x;
    double y2 = y*y;
    return ((2.73371257534512468e-01)*((1.0 + (3.65802904452657174e+00)*(-1.14304064540043701e-01 + 2.17201793541147964e-01*y)*x + (3.65802904452657174e+00)*(-1.36957589116038052e-02 + 4.12461891885490678e-02*y)*x2)/(1.0 + (-7.44128599675536129e-01 + 8.99849754167902249e-01*y)*x))*((1.0 + (3.65802904452657174e+00)*(-3.21684549917094428e-02)*y + (3.65802904452657174e+00)*(-9.40808316024106944e-02)*y2)/(1.0 + (-5.89639038902594170e-01)*y)));
}
/* === END fit_2d_to_c: omg0_l2m2 === */

/* === BEGIN fit_2d_to_c: domg0_l2m2 === */
double domg0_l2m2(double x, double y) {
    double x2 = x*x;
    double y2 = y*y;
    return ((6.03664316737079557e-03)*((1.0 + (1.65654979476870551e+02)*(-1.02840941698208762e-03 + 9.31127691979416819e-03*y)*x + (1.65654979476870551e+02)*(-1.87809567547554072e-03 + 4.68493729526423680e-03*y)*x2)/(1.0 + (-4.38021725085685140e-03 + 7.30617208534026896e-01*y)*x))*((1.0 + (1.65654979476870551e+02)*(5.25259082991305897e-03 + 1.41927898992048169e-03*x)*y + (1.65654979476870551e+02)*(-5.40475178495663025e-03 + -1.68772276752119114e-03*x)*y2)/(1.0 + (-4.76190113524564307e-01 + 3.27342166791816602e-01*x)*y)));
}
/* === END fit_2d_to_c: domg0_l2m2 === */

/* === BEGIN fit_2d_to_c: cAmp_l2m2_0 === */
double cAmp_l2m2_0(double x, double y) {
    double x2 = x*x;
    double y2 = y*y;
    return ((7.92857893446788642e-02)*(1.0 + (1.26126006723941799e+01)*(-1.62516034244934819e-02 + -5.32442242893030934e-03*y)*x + (1.26126006723941799e+01)*(-4.83140284311247073e-03 + -8.43758937150740698e-03*y)*x2)*(1.0 + (1.26126006723941799e+01)*(6.49449223018548020e-03 + 1.02578823374962981e-03*x)*y + (1.26126006723941799e+01)*(-1.04777326975056123e-02 + -2.14079556733814003e-03*x)*y2));
}
/* === END fit_2d_to_c: cAmp_l2m2_0 === */

/* === BEGIN fit_2d_to_c: cAmp_l2m2_1 === */
double cAmp_l2m2_1(double x, double y) {
    double x2 = x*x;
    double y2 = y*y;
    return ((-6.77268917975073714e-01)*((1.0 + (-1.47651837174196743e+00)*(7.30753461689540584e-02 + 2.60491910626122536e-01*y)*x + (-1.47651837174196743e+00)*(-5.49848457011212904e-02 + 5.56214192867383328e-02*y)*x2)/(1.0 + (-6.66666660765624863e-01 + 2.42466943148070657e-01*y)*x))*(1.0 + (-1.47651837174196743e+00)*(5.36950907493955021e-01 + 1.44249988588240036e-01*x)*y + (-1.47651837174196743e+00)*(-3.53423560160593342e-01 + -3.43099440031013725e-01*x)*y2));
}
/* === END fit_2d_to_c: cAmp_l2m2_1 === */

/* === BEGIN fit_2d_to_c: cphf_l2m2_0 === */
double cphf_l2m2_0(double x, double y) {
    double x2 = x*x;
    double y2 = y*y;
    return ((1.54723813969217444e-01)*((1.0 + (6.46312920000118130e+00)*(-1.25950376953620463e-01 + 7.32170371846544848e-02*y)*x + (6.46312920000118130e+00)*(-3.09648950254924442e-02 + 2.34723288226331255e-02*y)*x2)/(1.0 + (-9.09090821424046291e-01 + 8.25880645553469184e-01*y)*x))*((1.0 + (6.46312920000118130e+00)*(8.76874800717417524e-01)*y + (6.46312920000118130e+00)*(-1.73551445103763885e-01)*y2)/(1.0 + (4.65721911571504599e+00)*y)));
}
/* === END fit_2d_to_c: cphf_l2m2_0 === */

/* === BEGIN fit_2d_to_c: cphf_l2m2_1 === */
double cphf_l2m2_1(double x, double y) {
    double x2 = x*x;
    double y2 = y*y;
    return ((1.55047890230702312e+00)*((1.0 + (6.44962016904620739e-01)*(6.56059001787636875e-01 + 9.17559442843698392e-01*y)*x + (6.44962016904620739e-01)*(9.22191005769817718e-01 + -8.15626864654214301e-01*y)*x2)/(1.0 + (-6.66666662114926201e-01 + 1.63373145052067548e+00*y)*x))*(1.0 + (6.44962016904620739e-01)*(1.22911691298600068e+00 + -3.15335366655109473e-01*x)*y + (6.44962016904620739e-01)*(-8.75633452375245924e-01 + 1.08803703408451535e+00*x)*y2));
}
/* === END fit_2d_to_c: cphf_l2m2_1 === */

/* The eight hatA0_l3m2/hatdA0_l3m2/hatd2A0_l3m2/omg0_l3m2/domg0_l3m2/
   cAmp_l3m2_0/cphf_l3m2_0/cphf_l3m2_1 functions below are generated
   verbatim by eob_ringdown/py/fitter.py's fit_2d_to_c (see
   eob_ringdown/py/generic_q/gfits_A22_v1/l3m2.c) -- x=a0, y=1-X12, SAME
   convention as the l2m2 functions above. Naming matches l2m2's own
   (_l3m2, not the earlier _3_2 -- renamed together with this file's
   switch to fitting A0/dA0/d2A0/omg0/domg0 on the SPHERICAL basis, see
   below, to make the (l,m)-in-the-name convention consistent regardless
   of basis).

   hatA0_l3m2/hatdA0_l3m2/hatd2A0_l3m2/omg0_l3m2/domg0_l3m2 target the
   SPHERICAL (3,2) mode (comparable_mass_fitter.py's
   run_primary_fit_multistart, dual-fit change: A0/dA0/d2A0/omg0/domg0 are
   always fit on the spherical wave, even for l>m) -- used by
   A22_attach_mixed_mode_with_target as Y(3,2)'s own "natural" point-value
   target for the spheroidal-spherical mixing (see that function's own
   docstring), NOT as S32's target directly.

   cAmp_l3m2_0/cphf_l3m2_0/cphf_l3m2_1 are different: genuinely SPHEROIDAL
   (cAmp/cphf's own kind_pk='A22' fit was done on the spheroidal
   decomposition, not the spherical (3,2) mode) -- these feed S32's own
   template SHAPE directly (QNMHybridFitCab_A22_lm), only mixed with
   (2,2) at the VALUE/DERIVATIVE-target level, via A22_mu/
   A22_mix_target_at_point/A22_attach_mixed_mode.

   cAmp_l3m2_0 is fit under model_Amp='ext_tanhcosh_A1' (see
   A22_amp_uses_tanhcosh_A1, flipped ON for k32 below): it represents c3,
   the template's SOLE free amplitude shape coefficient, NOT c2 -- c2 is
   fixed analytically at 0.5*alpha21 (QNMHybridFitCab_A22_lm), so there
   is no cAmp_l3m2_1 anymore (a previous ext_tanhcosh_A2-fit version of
   this mode had one; removed along with the old A2 cAmp_l3m2_0 when this
   mode's Python fit was redone under A1). */

/* === BEGIN fit_2d_to_c: hatA0_l3m2 === */
double hatA0_l3m2(double x, double y) {
    double x2 = x*x;
    double y2 = y*y;
    double y3 = y2*y;
    return ((1.60964643045435757e-01)*((1.0 + (6.21254445125398380e+00)*(-6.54242533666669762e-02 + 8.32580754688301561e-02*y)*x + (6.21254445125398380e+00)*(-6.63112809284352976e-03 + 1.49072670789939077e-01*y)*x2)/(1.0 + (-8.24590647113676689e-01 + 5.12643359574348256e-01*y)*x))*(1.0 + (6.21254445125398380e+00)*(-3.15162681727964022e-01 + 4.84628662017182119e-02*x)*y + (6.21254445125398380e+00)*(2.90774191300878726e-01 + 1.89449854137422574e-01*x)*y2 + (6.21254445125398380e+00)*(-7.27890096914212692e-02 + -1.91399744699783769e-01*x)*y3));
}
/* === END fit_2d_to_c: hatA0_l3m2 === */

/* === BEGIN fit_2d_to_c: hatdA0_l3m2 === */
//double hatdA0_l3m2(double x, double y) {
//    double x2 = x*x;
//    double y2 = y*y;
//    double y3 = y2*y;
//    double y4 = y3*y;
//    return ((5.29198763964656461e-03)*((1.0 + (1.88964916038010045e+02)*(-4.97843368909211149e-03 + 1.25072713348087592e-02*y)*x + (1.88964916038010045e+02)*(-1.92826979812946598e-03 + 7.89802554384371434e-03*y)*x2)/(1.0 + (-1.75314087009973529e-01 + 6.19530584127318873e-01*y)*x))*((1.0 + (1.88964916038010045e+02)*(-2.73502182065837214e-02)*y + (1.88964916038010045e+02)*(5.17578533074357638e-02)*y2 + (1.88964916038010045e+02)*(-4.18013741957152449e-02)*y3 + (1.88964916038010045e+02)*(1.31954666057227700e-02)*y4)/(1.0 + (-2.26436033386185809e+00)*y + (1.59992872154669663e+00)*y2)));
//}
//double hatdA0_l3m2(double x, double y) {
//    double x2 = x*x;
//    double y2 = y*y;
//    double y3 = y2*y;
//    return ((5.29198763964656461e-03)*(1.0 + (1.88964916038010045e+02)*((-4.14460558960919219e-03 + 9.15078608862608493e-03*y)/(1.0 + 3.24181201496521798e-01*y))*x + (1.88964916038010045e+02)*((-2.42275705298551809e-03 + 3.58980437540631012e-02*y)/(1.0 + 6.35768349220929707e+00*y))*x2)*((1.0 + (1.88964916038010045e+02)*(-2.36249478943974904e-02)*y + (1.88964916038010045e+02)*(3.45214300576660721e-02)*y2 + (1.88964916038010045e+02)*(-1.46764600924891714e-02)*y3)/(1.0 + (-1.69949216519579882e+00)*y + (1.16661172704800986e+00)*y2)));
//}
double hatdA0_l3m2(double x, double y) {
    double x2 = x*x;
    double y2 = y*y;
    double y3 = y2*y;
    double z_fit = ((1.00529198763964667e+00)*(1.0 + (9.94735870070871719e-01)*(-4.14460558482392279e-03 + -1.35357537331327465e-01*y)*x + (9.94735870070871719e-01)*(-2.42275705175849353e-03 + 2.00286995788515124e-02*y)*x2)*((1.0 + (9.94735870070871719e-01)*(-1.72311644657728746e+00 + 1.69330584831125203e-01*x)*y + (9.94735870070871719e-01)*(1.20112553259372068e+00 + -2.76695531891434754e-01*x)*y2 + (9.94735870070871719e-01)*(-1.46765090011954664e-02 + 1.61082171318875700e-01*x)*y3)/(1.0 + (-1.69949149233591346e+00 + 2.54061921821024453e-02*x)*y + (1.16660407116527964e+00 + -3.65765945228090084e-02*x)*y2)));
    return (z_fit - (1.00000000000000000e+00)) / (1.00000000000000000e+00);
}
/* === END fit_2d_to_c: hatdA0_l3m2 === */

/* === BEGIN fit_2d_to_c: hatd2A0_l3m2 === */
/* NOTE: fit on a shifted/scaled copy of hatd2A0[(3,2)] (z_shift=1,
   z_scale=1) to avoid an ill-posed fit (data too close to 0); the inverse
   transform is already applied in the return statement below, so this
   function returns hatd2A0[(3,2)] itself, not the internal z_fit. */
double hatd2A0_l3m2(double x, double y) {
    double x2 = x*x;
    double y2 = y*y;
    double y3 = y2*y;
    double z_fit = ((1.00015642634349522e+00)*(1.0 + (9.99843598121878729e-01)*(-4.95491249720662365e-04 + -6.30556636403498266e-02*y)*x + (9.99843598121878729e-01)*(-9.47012598033288613e-05 + 2.79618136564744246e-03*y)*x2)*((1.0 + (9.99843598121878729e-01)*(-1.83872526811808346e+00 + 8.11649710659989726e-02*x)*y + (9.99843598121878729e-01)*(1.46954999309500223e+00 + -1.45364512602376761e-01*x)*y2 + (9.99843598121878729e-01)*(-2.07425511266372351e-03 + 9.27475617460521012e-02*x)*y3)/(1.0 + (-1.83653046648004703e+00 + 1.73755670241115198e-02*x)*y + (1.46528972101176747e+00 + -2.85424704511642972e-02*x)*y2)));
    return (z_fit - (1.00000000000000000e+00)) / (1.00000000000000000e+00);
}
/* === END fit_2d_to_c: hatd2A0_l3m2 === */

/* === BEGIN fit_2d_to_c: omg0_l3m2 === */
double omg0_l3m2(double x, double y) {
    double x2 = x*x;
    double y2 = y*y;
    double y3 = y2*y;
    return ((3.21205789913811046e-01)*((1.0 + (3.11326891171024478e+00)*(-1.97411938162358885e-01)*x + (3.11326891171024478e+00)*(-1.87056859446916440e-02)*x2)/(1.0 + (-8.16563602638818553e-01)*x))*((1.0 + (3.11326891171024478e+00)*(-6.02869704012361951e-01 + -3.14435084500446327e-01*x)*y + (3.11326891171024478e+00)*(6.43477261506784304e-01 + 5.42316072914238168e-02*x)*y2 + (3.11326891171024478e+00)*(2.75154083344231692e-01 + -3.48926105710420487e-01*x)*y3)/(1.0 + (-2.37715535046672555e+00 + -3.38357880336099814e-01*x)*y + (3.17823809715805039e+00 + -1.33946686283112171e+00*x)*y2)));
}
/* === END fit_2d_to_c: omg0_l3m2 === */

/* === BEGIN fit_2d_to_c: domg0_l3m2 === */
double domg0_l3m2(double x, double y) {
    double x2 = x*x;
    double y2 = y*y;
    double y3 = y2*y;
    return ((1.13465560786072239e-02)*((1.0 + (8.81324688365483979e+01)*(-1.02113763236537734e-02 + 4.18348992778966571e-03*y)*x + (8.81324688365483979e+01)*(-1.79507043285084779e-03 + 7.63015219370886174e-03*y)*x2)/(1.0 + (-5.63563458762453728e-01 + 3.62730831919867158e-01*y)*x))*((1.0 + (8.81324688365483979e+01)*(-2.55132017697322806e-02)*y + (8.81324688365483979e+01)*(3.37875878419301862e-03)*y2 + (8.81324688365483979e+01)*(2.92914306610257205e-02)*y3)/(1.0 + (-3.29525776151496963e+00)*y + (3.67049811575346085e+00)*y2)));
}
/* === END fit_2d_to_c: domg0_l3m2 === */

// === BEGIN fit_2d_to_c: cAmp_l3m2_0 ===
double cAmp_l3m2_0(double x, double y) {
    double x2 = x*x;
    double y2 = y*y;
    return ((-1.02973256973294536e+00)*((1.0 + (-9.71125930550437610e-01)*(3.90983156866200221e-01)*x + (-9.71125930550437610e-01)*(5.03877517520369578e-01)*x2)/(1.0 + (-9.99999986707884259e-01)*x))*(1.0 + (-9.71125930550437610e-01)*(1.72993370257658063e+00 + 9.76294902366990069e-02*x)*y + (-9.71125930550437610e-01)*(-7.16354935746834398e-01 + -9.11197388071349534e-01*x)*y2));
}
// === END fit_2d_to_c: cAmp_l3m2_0 ===

/* === BEGIN fit_2d_to_c: cphf_l3m2_0 === */
double cphf_l3m2_0(double x, double y) {
    double x2 = x*x;
    double y2 = y*y;
    return ((1.78178382437494953e-01)*((1.0 + (5.61235311669079984e+00)*(-1.51462781177380879e-01 + 3.29078592166769240e-02*y)*x + (5.61235311669079984e+00)*(-2.11022022502822398e-02 + -2.07396649680172074e-01*y)*x2)/(1.0 + (-8.25258729270220592e-01 + -3.45614236769241426e-02*y)*x))*(1.0 + (5.61235311669079984e+00)*(-5.02969109410731294e-02 + -6.45566310047882008e-02*x)*y + (5.61235311669079984e+00)*(1.05199432940434226e-01 + -2.69809169965290474e-01*x)*y2));
}
/* === END fit_2d_to_c: cphf_l3m2_0 === */

/* === BEGIN fit_2d_to_c: cphf_l3m2_1 === */
double cphf_l3m2_1(double x, double y) {
    double x2 = x*x;
    double y2 = y*y;
    return ((3.72480030488304203e+00)*(1.0 + (2.68470768403086202e-01)*((3.15495293408449440e+00)/(1.0 + -2.32097522319127086e-01*y))*x + (2.68470768403086202e-01)*((3.36531552030830872e+00)/(1.0 + 7.81032063711039104e+00*y))*x2)*(1.0 + (2.68470768403086202e-01)*(-5.44494033947331602e+00 + -2.81073791677560392e+00*x)*y + (2.68470768403086202e-01)*(8.31137712495789138e+00 + -4.68375478882972907e+00*x)*y2));
}
/* === END fit_2d_to_c: cphf_l3m2_1 === */


/* ============================================================
 * Placeholder (nu,spin) global-fit surface -- returns 0 for every
 * mode that doesn't have real fits yet (everything except (2,2),(3,2)).
 * ============================================================ */
static double A22_fit_placeholder(double x, double y) { return 0.0; }

typedef double (*A22_fit2d_t)(double, double);

typedef struct {
  A22_fit2d_t hatA0, hatdA0, hatd2A0, omg0, domg0, cAmp0, cAmp1, cphf0, cphf1;
} A22_fit_set_t;

/**
 * Function: A22_get_fit_set
 * --------------------------
 *   Per-mode (nu,spin) global-fit surface lookup: real fits for (2,2),
 *   A22_fit_placeholder (identically 0) for every other supported mode.
 *   Swap in real per-mode functions here (same fit_2d_to_c-generated
 *   shape as the l2m2 ones above) once they exist -- nothing else in
 *   this file needs to change.
 */
/**
 * Function: A22_amp_uses_tanhcosh_A1
 * ------------------------------------
 *   Modes using waveform_fit.py's 'ext_tanhcosh_A1' amplitude template
 *   (ONE free shape coefficient, c3 -- c2 is NOT fit at all, it is FIXED
 *   analytically at 0.5*alpha21 = 0.5*(alpha2-alpha1), the difference of
 *   inverse damping times between the n=1 overtone and the fundamental
 *   mode QNM_coefs already computes -- the SAME c2A the ORIGINAL,
 *   production ringdown model uses, see TEOBResumSFits.c) instead of the
 *   default 'ext_tanhcosh_A2' (two free shape coefficients, c2 AND c3,
 *   both read from this mode's own (nu,spin) global-fit surface).
 *
 *   Returns 0 (A2, the default) for every mode except the ones flipped
 *   on below. Only flip a mode's k index on ONCE its OWN cAmp0 fit
 *   surface has been regenerated (fit_2d_to_c, from a Python primary-fit
 *   batch run under model_Amp='ext_tanhcosh_A1') to represent c3 ALONE
 *   -- NOT before: with a mode's fit surface still under
 *   ext_tanhcosh_A2 (cAmp0=c2, cAmp1=c3), flipping it on here would
 *   silently misinterpret a fitted c2 surface as a fixed-formula c3
 *   input -- wrong physics, not just a placeholder-0. cAmp1 goes unused
 *   for a mode flipped on here (A22_get_fit_set wires it to
 *   A22_fit_placeholder -- see QNMHybridFitCab_A22_lm, which never
 *   reads fs.cAmp1 for a mode this function returns true for).
 */
static int A22_amp_uses_tanhcosh_A1(int k)
{
  const int k32 = 3;
  return k == k32;   /* (3,2): cAmp_l3m2_0 refit under ext_tanhcosh_A1 --
                         see gfits_A22_v1/l3m2.c */
}

static void A22_get_fit_set(int k, A22_fit_set_t *fs)
{
  const int k22 = 1;
  const int k32 = 3;
  if (k == k22) {
    fs->hatA0 = hatA0_l2m2;   fs->hatdA0 = hatdA0_l2m2;   fs->hatd2A0 = hatd2A0_l2m2;
    fs->omg0  = omg0_l2m2;    fs->domg0  = domg0_l2m2;
    fs->cAmp0 = cAmp_l2m2_0;  fs->cAmp1  = cAmp_l2m2_1;
    fs->cphf0 = cphf_l2m2_0;  fs->cphf1  = cphf_l2m2_1;
  } else if (k == k32) {
    fs->hatA0 = hatA0_l3m2;   fs->hatdA0 = hatdA0_l3m2;   fs->hatd2A0 = hatd2A0_l3m2;
    fs->omg0  = omg0_l3m2;    fs->domg0  = domg0_l3m2;
    /* cAmp0 = c3 (ext_tanhcosh_A1 -- A22_amp_uses_tanhcosh_A1(k32)==1
       below); cAmp1 is unused for this mode (QNMHybridFitCab_A22_lm
       never reads fs.cAmp1 when A22_amp_uses_tanhcosh_A1 is true) --
       wired to the placeholder rather than left dangling. */
    fs->cAmp0 = cAmp_l3m2_0;  fs->cAmp1  = A22_fit_placeholder;
    fs->cphf0 = cphf_l3m2_0;  fs->cphf1  = cphf_l3m2_1;
  } else {
    fs->hatA0 = fs->hatdA0 = fs->hatd2A0 = A22_fit_placeholder;
    fs->omg0  = fs->domg0  = A22_fit_placeholder;
    fs->cAmp0 = fs->cAmp1  = A22_fit_placeholder;
    fs->cphf0 = fs->cphf1  = A22_fit_placeholder;
  }
}

/**
 * Function: A22_mu
 * -----------------
 *   Spheroidal-spherical mixing coefficient mu_{m,l,l',n=0}(af), Berti-
 *   Klein fits (arXiv:1408.1860). Hardcoded here for exactly the 18
 *   (m,l,l') combinations this file needs -- the diagonal ones (2,2),
 *   (3,3),(4,4),(5,5) (needed for eob_wav_ringdown_A22's base-mode
 *   normalization) plus (2,1),(3,2),(3,3->2),(4,3),(4,4->3) i.e. every
 *   mu that appears in the (3,2)-from-(2,2) and (4,3)-from-(3,3) mixing
 *   formulas -- TIMES TWO, for BOTH signs of m (mm below can be
 *   NEGATIVE, for af<0): the fit table has GENUINELY SEPARATE rows for
 *   +m and -m (not related by any simple sign flip -- verified directly,
 *   e.g. (m,l,lp)=(1,2,2) vs (-1,2,2) have completely different p1..p4
 *   coefficients), so an af<0 remnant needs its OWN, independently-
 *   hardcoded -m row, not a derived one. Values taken verbatim from
 *   ~/repos/kerrorbits/kerrorbits/modemix_spheroidal_fits.txt (n'=0
 *   rows only), same fit form/sign convention documented in that repo's
 *   modemix_spheroidal_readme.txt and independently confirmed identical
 *   in eob_ringdown/py/spheroidal.py's own
 *   spheroidal_spherical_mixing_coefficient (which performs the SAME
 *   direct (m,l,lp) table lookup with the sign-flipped m, not a
 *   symmetry relation):
 *     mm,spin = (-m,-af) if af<0 else (m,af)      [af here is the ONLY
 *                                                   spin dependence -- a
 *                                                   single remnant spin]
 *     delta   = 1 if l==lp else 0
 *     re      = delta + p1r*spin^p2r + p3r*spin^p4r
 *     im      =         p1i*spin^p2i + p3i*spin^p4i
 *     mu      = (re + i*im) * (-1)^(l+lp)         [applied AGAIN if af<0]
 *     mu      = conj(mu)                          [the value actually used]
 */
static void A22_mu(int m, int l, int lp, double af, double *mu_re, double *mu_im)
{
  int mm = m;
  double spin = af;
  if (af < 0.0) { mm = -m; spin = -af; }

  double p1r,p2r,p3r,p4r,p1i,p2i,p3i,p4i;
  if      (mm== 1 && l==2 && lp==2) { p1r=-6.50602e-03;p2r=2.29605e+00;p3r=-3.08509e-03;p4r=9.09064e+00;p1i=-7.59772e-03;p2i=8.69142e-01;p3i=4.79185e-03;p4i=5.78945e+00; }
  else if (mm== 2 && l==2 && lp==2) { p1r=-7.39789e-03;p2r=2.88933e+00;p3r=-6.61369e-03;p4r=1.71287e+01;p1i= 1.53046e-02;p2i=1.21928e+00;p3i=-9.34293e-03;p4i=2.49915e+01; }
  else if (mm== 2 && l==3 && lp==2) { p1r=-1.03512e-01;p2r=1.22285e+00;p3r=-5.74989e-02;p4r=8.70536e+00;p1i=-1.60040e-02;p2i=9.53385e-01;p3i= 1.00344e-02;p4i=1.47550e+01; }
  else if (mm== 2 && l==3 && lp==3) { p1r=-2.51626e-02;p2r=2.43323e+00;p3r=-1.30029e-02;p4r=1.09648e+01;p1i=-3.23812e-02;p2i=9.21248e-01;p3i= 1.88165e-02;p4i=1.06988e+01; }
  else if (mm== 3 && l==3 && lp==3) { p1r=-1.20858e-02;p2r=2.94266e+00;p3r=-1.12402e-02;p4r=1.72866e+01;p1i=-1.72988e-02;p2i=8.51080e-01;p3i= 1.43925e-02;p4i=5.98677e+00; }
  else if (mm== 3 && l==4 && lp==3) { p1r=-1.32724e-01;p2r=1.24145e+00;p3r=-7.34932e-02;p4r=8.51670e+00;p1i=-1.15656e-02;p2i=9.01101e-01;p3i= 7.32291e-03;p4i=1.17007e+01; }
  else if (mm== 3 && l==4 && lp==4) { p1r=-3.42425e-02;p2r=2.59986e+00;p3r=-2.12509e-02;p4r=1.25368e+01;p1i=-4.54823e-02;p2i=9.25823e-01;p3i= 2.84355e-02;p4i=1.14204e+01; }
  else if (mm== 4 && l==4 && lp==4) { p1r=-1.45512e-02;p2r=2.96007e+00;p3r=-1.40395e-02;p4r=1.75478e+01;p1i=-3.40880e-02;p2i=8.81775e-01;p3i= 2.48362e-02;p4i=8.29363e+00; }
  else if (mm== 5 && l==5 && lp==5) { p1r=-1.56530e-02;p2r=2.96856e+00;p3r=-1.56783e-02;p4r=1.78822e+01;p1i=-4.53719e-02;p2i=9.00055e-01;p3i= 3.20586e-02;p4i=9.31368e+00; }
  else if (mm==-1 && l==2 && lp==2) { p1r=-5.01437e-03;p2r=1.97161e+00;p3r= 3.89867e-04;p4r=2.96189e+00;p1i=-4.98154e-02;p2i=1.00177e+00;p3i= 2.47007e-03;p4i=4.38106e+00; }
  else if (mm==-2 && l==2 && lp==2) { p1r=-1.43708e-02;p2r=2.11765e+00;p3r= 1.03499e-02;p4r=2.22856e+00;p1i=-7.01457e-02;p2i=1.00458e+00;p3i= 6.65951e-04;p4i=3.52711e+00; }
  else if (mm==-2 && l==3 && lp==2) { p1r=-1.34751e-01;p2r=1.08821e+00;p3r= 7.96321e-02;p4r=1.27858e+00;p1i=-1.74388e-02;p2i=1.01083e+00;p3i= 5.15606e-03;p4i=1.82061e+00; }
  else if (mm==-2 && l==3 && lp==3) { p1r=-2.00333e-02;p2r=2.03214e+00;p3r= 7.15335e-03;p4r=2.49949e+00;p1i=-8.94720e-02;p2i=1.00714e+00;p3i= 2.71266e-03;p4i=3.76609e+00; }
  else if (mm==-3 && l==3 && lp==3) { p1r=-1.78203e-01;p2r=2.15584e+00;p3r= 1.73530e-01;p4r=2.16591e+00;p1i=-7.18064e-02;p2i=1.01171e+00;p3i= 2.99296e-04;p4i=5.06507e+00; }
  else if (mm==-3 && l==4 && lp==3) { p1r=-1.55284e-01;p2r=1.08704e+00;p3r= 9.09275e-02;p4r=1.31882e+00;p1i=-1.63145e-02;p2i=1.04521e+00;p3i= 9.11357e-03;p4i=1.65083e+00; }
  else if (mm==-3 && l==4 && lp==4) { p1r=-5.85885e-01;p2r=2.17527e+00;p3r= 5.73950e-01;p4r=2.18532e+00;p1i=-8.41048e-02;p2i=1.00883e+00;p3i= 1.66840e-03;p4i=3.77298e+00; }
  else if (mm==-4 && l==4 && lp==4) { p1r=-1.84411e-01;p2r=2.14629e+00;p3r= 1.79545e-01;p4r=2.15630e+00;p1i=-7.45925e-02;p2i=1.01474e+00;p3i= 1.84831e-04;p4i=3.40093e+01; }
  else if (mm==-5 && l==5 && lp==5) { p1r=-1.77413e-01;p2r=2.13982e+00;p3r= 1.72488e-01;p4r=2.14990e+00;p1i=-7.73863e-02;p2i=1.01807e+00;p3i= 1.67690e-04;p4i=7.54301e+00; }
  else {
    errorexit("A22_mu: no Berti-Klein mixing fit hardcoded for this (m,l,l') -- "
             "only the combinations needed for (2,1),(2,2),(3,2),(3,3),(4,3),"
             "(4,4),(5,5) are available.\n");
    p1r=p2r=p3r=p4r=p1i=p2i=p3i=p4i=0.0; /* unreached, silences -Wmaybe-uninitialized */
  }

  const double delta = (l == lp) ? 1.0 : 0.0;
  const double re_part = delta + p1r*pow(spin,p2r) + p3r*pow(spin,p4r);
  const double im_part =         p1i*pow(spin,p2i) + p3i*pow(spin,p4i);
  const double sign = ((l+lp) % 2 == 0) ? 1.0 : -1.0; /* (-1)^(l+lp) */

  double mu_r = re_part*sign;
  double mu_i = im_part*sign;
  if (af < 0.0) { mu_r *= sign; mu_i *= sign; }

  *mu_re =  mu_r; /* conj() below only flips the imaginary part */
  *mu_im = -mu_i;
}

/**
 *  Function: QNMHybridFitCab_A22_lm
 *  ----------------------------------
 *  Postpeak-template coefficients + NQC target point for mode k
 *  (KMAX index), new ringdown-NQC model. Generalizes the original
 *  (2,2)-only QNMHybridFitCab_A22 to any of the modes this file
 *  supports (see eob_wav_ringdown_A22's own mode-support guard).
 *
 *  QNM complex frequency (sigmar,sigmai = alpha1,omega1): Kerr
 *  perturbation theory, not something the (nu,spin) global-fit surfaces
 *  replace. (2,2) keeps its own long-validated hardcoded fit (UNCHANGED
 *  from QNMHybridFitCab's spinning branch) -- QNM_coefs()'s own (2,2)
 *  entry is a DIFFERENT, close-but-not-identical fit for the same
 *  physical quantity (verified: differs at the 1e-4-1e-3 level across
 *  spin), and this session's whole validation history (dephasing checks,
 *  SXS comparisons, the eob_wav_ringdown_HM bugfix cross-check, etc.)
 *  used THIS exact formula, so switching would be a silent regression
 *  for the one mode that's actually been tested. Every other mode uses
 *  QNM_coefs (already implemented in this repo, TEOBResumSFits.c,
 *  covers exactly this mode set).
 *
 *  A0/dA0/d2A0/omg0/domg0 come from this mode's own global-fit surface
 *  (A22_get_fit_set -- real for (2,2), placeholder-0 otherwise),
 *  evaluated at THIS binary's own (a0, 1-X12), then un-normalized from
 *  hat-prefixed back to plain via the eq. 63-style factor (1-Shat*omg0)
 *  (arXiv:2001.09082 eq. 63, generalized to dA0/d2A0 -- see
 *  comparable_mass_fitter.py's SPIN_AMP_NORM/_spin_amp_norm_factor, this
 *  mode's own kind='S_hat_omega', coeff=1.0).
 *
 *  Amplitude coefficients (ca1..ca5) implement waveform_fit.py's
 *  'ext_tanhcosh_A2' template exactly:
 *    y(x) = ca1*tanh(ca2 x+ca3) + ca4 + ca5/cosh(ca2 x+ca3)
 *  matching A0, dA0, d2A0 at x=0 by construction (ca2,ca3 are the two
 *  free/fitted shape coefficients). Guarded against cA2==0 (the
 *  placeholder value): the algebra below divides by SQ(cA2), which a
 *  genuine fit is never exactly 0 at, but the placeholder always is.
 *
 *  Phase coefficients (cb1..cb3) implement 'ext_ph2' exactly:
 *    phase(x) = -cb1*log((1+cb3*exp(-cb2 x))/(1+cb3))
 *  with cb2,cb3 the two free/fitted coefficients and cb1 fixed by
 *  delta_omg0 = Im(sigma_plus1) - Mbh*omg0 (the Mbh factor is NOT
 *  optional: omg0 here is in TOTAL-mass units, matching the
 *  EOB-literature omg_mrg convention, while sigmai/omega1 are in
 *  REMNANT-mass (Mbh) units, the QNM-table's own natural normalization --
 *  see compare_global_fit.py's build_global_overrides for the identical
 *  relation, verified there to machine precision). Guarded against
 *  d1f==0 or d2f==0 (placeholder) for the same reason as cA2 above.
 *
 *  With every placeholder coefficient 0, both guards trigger and
 *  ca1=ca4=ca5=ca2=ca3=cb1=cb2=cb3=0 -- eob_wav_ringdown_template_A22
 *  then returns amp=0 identically (phase is finite but irrelevant,
 *  multiplied by 0 amplitude everywhere it's used).
 *
 *   @param[in] k : KMAX mode index (see eob_wav_ringdown_A22's own guard
 *     for the supported set)
 *   @param[in] nu,X1,X2,chi1,chi2,Mbh,abh
 *   @param[out] ca1,ca2,ca3,ca4,ca5 : amplitude template coefficients
 *   @param[out] cb1,cb2,cb3         : phase template coefficients
 *   @param[out] sigmar,sigmai       : QNM complex frequency (Mbh units)
 *   @param[out] A0_out,dA0_out,d2A0_out,omg0_out,domg0_out : NQC target
 *     point (plain, un-normalized, at t0) -- see
 *     eob_wav_hlmNQC_find_a1a2a3_mrg_A22. d2A0_out is this mode's OWN
 *     value (needed to combine (3,2)/(4,3)'s own spheroidal target with
 *     (2,2)/(3,3)'s spherical one via A22_mix_target_at_point -- see
 *     there); for modes that don't need mixing it can simply be ignored.
 */
/**
 * Function: A22_build_template_coeffs
 * --------------------------------------
 *   Given a target point (A0,dA0,d2A0,omg0), this mode's own QNM
 *   (alpha1,omega1), and its own template SHAPE (cA2,cA3,d1f,d2f, from
 *   a (nu,spin) fit surface), build the ext_tanhcosh_A2/ext_ph2
 *   coefficients ca1..cb3 -- the exact algebra QNMHybridFitCab_A22_lm
 *   itself uses, factored out so it can ALSO be called with a target
 *   point that is NOT read directly off a fit surface (see
 *   eob_wav_ringdown_A22's (3,2)/(4,3) branches, which derive the target
 *   by properly mixing (2,2)'s/(3,3)'s own target with (3,2)'s/(4,3)'s
 *   own natural spherical inspiral value -- see A22_mix_target_at_point).
 *   See QNMHybridFitCab_A22_lm's own docstring for the ext_tanhcosh_A2/
 *   ext_ph2 formulas and the cA2==0/d1f==0||d2f==0 (placeholder) guards.
 */
static void A22_build_template_coeffs(double A0, double dA0, double d2A0,
                                      double alpha1, double omega1, double Mbh, double omg0,
                                      double cA2, double cA3, double d1f, double d2f,
                                      double *ca1, double *ca2, double *ca3, double *ca4, double *ca5,
                                      double *cb1, double *cb2, double *cb3)
{
  const double delta_omg0 = omega1 - Mbh*omg0;

  if (cA2 == 0.0) {
    *ca1 = 0.0; *ca4 = 0.0; *ca5 = 0.0; *ca2 = 0.0; *ca3 = 0.0;
  } else {
    const double s3   = sinh(cA3);
    const double c3   = cosh(cA3);
    const double t3   = tanh(cA3);
    const double s2c3 = sinh(2.*cA3);
    const double c2c3 = cosh(2.*cA3);
    const double Aa   = A0*alpha1 + dA0;
    const double Ab   = A0*alpha1*alpha1 + d2A0 + 2.*alpha1*dA0;
    *ca1 = -0.5/SQ(cA2) * ( cA2*Aa*(c2c3-3.) + Ab*s2c3 );
    *ca4 = A0 + Ab*SQ(c3)/SQ(cA2) + Aa*c3*s3/cA2;
    *ca5 = -c3/SQ(cA2) * ( Ab + 2.*cA2*Aa*t3 );
    *ca2 = cA2;
    *ca3 = cA3;
  }

  if (d1f == 0.0 || d2f == 0.0) {
    *cb1 = 0.0; *cb2 = 0.0; *cb3 = 0.0;
  } else {
    *cb1 = delta_omg0*(1. + d2f)/(d1f*d2f);
    *cb2 = d1f;
    *cb3 = d2f;
  }
}

void QNMHybridFitCab_A22_lm(int k, double nu, double X1, double X2, double chi1, double chi2,
                            double Mbh, double abh,
                            double *ca1, double *ca2, double *ca3, double *ca4, double *ca5,
                            double *cb1, double *cb2, double *cb3,
                            double *sigmar, double *sigmai,
                            double *A0_out, double *dA0_out, double *d2A0_out,
                            double *omg0_out, double *domg0_out,
                            double *cA2_out, double *cA3_out, double *d1f_out, double *d2f_out)
{
  const int k22 = 1;
  const double af  = abh;

  double alpha1, omega1;
  double alpha21_k = 0.0;  /* only meaningful (and only ever used) for a
                              mode in A22_amp_uses_tanhcosh_A1 -- (2,2)
                              never needs it, its own hardcoded fit below
                              is untouched */
  if (k == k22) {
    const double af2 = SQ(af);
    const double af3 = af2*af;
    const double omega1_c = -0.0598837831*af3 + 0.8082136788*af2 - 1.7408467418*af + 1;
    const double omega1_d = -0.2358960279*af3 + 1.3152369374*af2 - 2.0764065380*af + 1;
    omega1 =  0.3736716844 * (omega1_c/omega1_d);

    const double alpha1_c =  0.1211263886*af3 + 0.7015835813*af2 - 1.8226060896*af + 1;
    const double alpha1_d =  0.0811633377*af3 + 0.7201166020*af2 - 1.8002031358*af + 1;
    alpha1 =  0.0889623157 * (alpha1_c/alpha1_d);
  } else {
    double alpha21[KMAX], alpha1_arr[KMAX], omega1_arr[KMAX];
    QNM_coefs(af, alpha21, alpha1_arr, omega1_arr);
    alpha1 = alpha1_arr[k];
    omega1 = omega1_arr[k];
    alpha21_k = alpha21[k];
  }

  A22_fit_set_t fs;
  A22_get_fit_set(k, &fs);

  /* (nu, spin) global-fit surface coordinates -- x=a0, y=1-X12, and Shat
     for the eq. 63 un-normalization -- same formulas
     QNMHybridFitCab/compute_spin_variables already use */
  const double a1c  = X1*chi1;
  const double a2c  = X2*chi2;
  const double a0   = a1c + a2c;
  const double a12  = a1c - a2c;
  const double X12  = X1 - X2;
  const double y    = 1.0 - X12;
  const double Shat = 0.5*(a0 + X12*a12);

  const double hatA0   = fs.hatA0(a0, y);
  const double hatdA0  = fs.hatdA0(a0, y);
  const double hatd2A0 = fs.hatd2A0(a0, y);
  const double omg0    = fs.omg0(a0, y);
  const double domg0   = fs.domg0(a0, y);
  /* ext_tanhcosh_A1 (A22_amp_uses_tanhcosh_A1): cA2 is FIXED at
     0.5*alpha21_k, not fit -- this mode's fit surface (cAmp0) then holds
     the sole free coefficient c3 directly, and cAmp1 goes unused.
     Every other mode is untouched: both cA2 (cAmp0) and cA3 (cAmp1)
     read from the fit surface, exactly as before this option existed. */
  const int use_tanhcosh_A1 = A22_amp_uses_tanhcosh_A1(k);
  const double cA2 = use_tanhcosh_A1 ? 0.5*alpha21_k : fs.cAmp0(a0, y);
  const double cA3 = use_tanhcosh_A1 ? fs.cAmp0(a0, y) : fs.cAmp1(a0, y);
  const double d1f     = fs.cphf0(a0, y);
  const double d2f     = fs.cphf1(a0, y);

  const double norm_factor = 1.0 - Shat*omg0;
  const double A0   = hatA0   * norm_factor;
  const double dA0  = hatdA0  * norm_factor;
  const double d2A0 = hatd2A0 * norm_factor;

  A22_build_template_coeffs(A0, dA0, d2A0, alpha1, omega1, Mbh, omg0, cA2, cA3, d1f, d2f,
                            ca1, ca2, ca3, ca4, ca5, cb1, cb2, cb3);

  *sigmar = alpha1;
  *sigmai = omega1;

  if (A0_out)    *A0_out    = A0;
  if (dA0_out)   *dA0_out   = dA0;
  if (d2A0_out)  *d2A0_out  = d2A0;
  if (omg0_out)  *omg0_out  = omg0;
  if (domg0_out) *domg0_out = domg0;
  if (cA2_out)   *cA2_out   = cA2;
  if (cA3_out)   *cA3_out   = cA3;
  if (d1f_out)   *d1f_out   = d1f;
  if (d2f_out)   *d2f_out   = d2f;
}

/**
 * Function: eob_wav_ringdown_template_A22
 * ----------------------------------------
 *   New (2,2) postpeak template -- tanh+1/cosh amplitude
 *   (waveform_fit.py's 'ext_tanhcosh_A2'), single-exponential-log phase
 *   ('ext_ph2') -- same psi[0]/psi[1] output convention as
 *   eob_wav_ringdown_template (amplitude, and MINUS(phase-omega*x) for
 *   the phase).
 *
 *   logarg is clamped to a tiny positive floor rather than left free to
 *   go non-positive: verified directly (high-spin, unequal-mass
 *   SXS:BBH:0293-like parameters) that a (nu,spin) fit surface evaluated
 *   near/outside its own calibrated domain can produce a cb2/cb3 pair
 *   for which (1+cb3*exp(-cb2*x))/(1+cb3) crosses 0 almost immediately
 *   past x=0 (not just at extreme x) -- log() of that is NaN, which then
 *   silently corrupts the ENTIRE final waveform once the ringdown-tail
 *   samples feed a downstream cubic-spline resample (interp_uniform_grid),
 *   not just the ringdown itself. This does not change the template for
 *   any (cb1,cb2,cb3) that stay in the physically-intended regime (where
 *   the argument is always positive).
 */
void eob_wav_ringdown_template_A22(double x,
                                   double ca1, double ca2, double ca3, double ca4, double ca5,
                                   double cb1, double cb2, double cb3,
                                   double sigmar, double sigmai, double *psi)
{
  const double z     = ca2*x + ca3;
  const double amp   = ca1*tanh(z) + ca4 + ca5/cosh(z);
  double logarg = (1. + cb3*exp(-cb2*x))/(1. + cb3);
  if (!(logarg > 0.0)) logarg = 1e-300; /* catches NaN (Inf-Inf/...) too */
  const double phase = -cb1*log(logarg);
  psi[0] = amp * exp(-sigmar*x);
  psi[1] = -(phase - sigmai*x);
}

/**
 * Function: eob_A22_find_t0_raw
 * ------------------------------
 *   t0 (RAW, M=1 units) for the new ringdown-NQC model: the EXACT SAME
 *   procedure eob_wav_ringdown_v1 itself uses for tmrg[k22] (refined peak
 *   of the orbital frequency Omega_orb via 7-point local polynomial fit,
 *   including its edge extrapolation for when the peak sits too close to
 *   the end of the dynamics array), minus the same NR-calibrated
 *   eob_nqc_timeshift(nu,chi1) shift AND the same extra "-2M" term
 *   eob_wav_ringdown_v1 itself applies (tmrgA22 = tOmg_pk-(DeltaT_nqc+2)/Mbh)
 *   -- that "-2M" is NOT an attachment-convention artifact, it is part of
 *   the NR-calibrated proxy itself: it's what makes the Omega_orb-peak
 *   proxy line up with the TRUE (2,2) amplitude peak (eob_wav_ringdown_v1
 *   then attaches its own template 2M further on, at tmatch=tmrg+2M, i.e.
 *   at template-argument x=+2M, not x=0). The (nu,spin) global-fit
 *   surfaces above are calibrated at the TRUE |h22| peak directly
 *   (waveform_fit.py's kind_pk='A22': argmax(abs(h22)), see
 *   PrimaryFit/get_peak) -- the same point tmrg proxies for -- so t0 here
 *   must match tmrg, not tmatch. The new model then attaches directly AT
 *   this t0 (template-argument x=0 there), unlike the old model's extra
 *   +2M attachment delay past its own tmrg.
 *
 *   Always called with the FULL dynamics `dyn` (both from
 *   eob_wav_ringdown_A22 and eob_wav_hlmNQC_find_a1a2a3_mrg_A22, which
 *   also receives the merger-segment-only dyn_mrg but deliberately does
 *   NOT use it here) so both consumers agree on IDENTICALLY the same t0.
 */
double eob_A22_find_t0_raw(Dynamics *dyn, double nu, double chi1)
{
  double *Omega = dyn->data[EOB_OMGORB];
  const int dynsize = dyn->size;

  int index_pk = dynsize-1;
  double Omega_pk = Omega[index_pk];
  for (int j = dynsize-2; j--; ) {
    if (Omega[j] < Omega_pk) break;
    index_pk = j;
    Omega_pk = Omega[j];
  }

  /* LOCAL sample spacing around index_pk, NOT dyn->dt: this function is
     called TWICE per EOBRunPy run (once from
     eob_wav_hlmNQC_find_a1a2a3_mrg_A22, on the pre-merger-join `dyn`; once
     from eob_wav_ringdown_A22, on `dyn` AFTER TEOBResumS.c's
     Dynamics_join() has spliced in the uniformly-interpolated
     merger-segment dynamics at dt_merger_interp) -- Dynamics_join() only
     copies in the new tail's time/data arrays, it does NOT update the
     scalar dyn->dt field, which stays stuck at whatever the ORIGINAL
     (pre-join, typically much coarser adaptive-step) value was. Using
     that stale dt in the find_max() 7-point local-polynomial fit below
     (which assumes UNIFORM spacing dt between its input points) silently
     computed a badly wrong t0 in the second (post-join) call -- verified
     directly: dyn->dt stuck at 0.1 while the ACTUAL post-join spacing at
     index_pk was 0.5, producing a >100M spurious difference between the
     NQC-fit call's and the ringdown-attach call's own idea of "t0", which
     showed up as a real amplitude/phase discontinuity at attachment
     (worst for (3,2), whose NQC target/ringdown-template split across
     this same t0 amplifies any such mismatch through the spheroidal-
     spherical mixing). Recomputing dt from the time array directly next
     to index_pk is correct in BOTH calls regardless of any join. */
  const double dt = dyn->time[index_pk] - dyn->time[index_pk-1];

  const int n = 7;
  double tmax = dyn->time[index_pk];
  double *Omega_ptr = &Omega[index_pk-3];
  double tOmg_pk;

  if ( (index_pk + (n-1)/2) > (dynsize-1) ) {
    double Omega_pk_grid[7];
    const int ni = (index_pk + (n-1)/2) - (dynsize-1);
    for (int j = 0; j < (7-ni); j++)
      Omega_pk_grid[j] = Omega_ptr[j];
    if (ni==1) {
      Omega_pk_grid[6] = 2.*Omega_pk_grid[5]-Omega_pk_grid[4];
    } else if (ni==2) {
      Omega_pk_grid[5] = 2.*Omega_pk_grid[4]-Omega_pk_grid[3];
      Omega_pk_grid[6] = 2.*Omega_pk_grid[5]-Omega_pk_grid[4];
    } else if (ni==3) {
      Omega_pk_grid[4] = 2.*Omega_pk_grid[3]-Omega_pk_grid[2];
      Omega_pk_grid[5] = 2.*Omega_pk_grid[4]-Omega_pk_grid[3];
      Omega_pk_grid[6] = 2.*Omega_pk_grid[5]-Omega_pk_grid[4];
    } else {
      errorexit("eob_A22_find_t0_raw: wrong counting (ni)\n");
    }
    tOmg_pk = find_max(n, dt, tmax, Omega_pk_grid, NULL);
  } else {
    tOmg_pk = find_max(n, dt, tmax, Omega_ptr, NULL);
  }

  const double DeltaT_nqc = eob_nqc_timeshift(nu, chi1);
  return tOmg_pk - DeltaT_nqc - 2.0;
}

/* KMAX indices of the modes this file supports (l(l-1)/2+m-2 convention,
   same as everywhere else in TEOBResumS -- see LINDEX/MINDEX). */
#define A22_K21 0
#define A22_K22 1
#define A22_K32 3
#define A22_K33 4
#define A22_K43 7
#define A22_K44 8
#define A22_K55 13

/**
 * Function: A22_attach_mixed_mode
 * ---------------------------------
 *   Attach mode (l,m) in [idx,size) when it needs spheroidal-spherical
 *   MIXING with the already-attached Y(l-1,m) (currently only (3,2)
 *   mixing with (2,2), and (4,3) mixing with (3,3)):
 *
 *     Y(l,m) = [mu(m,l,l-1,0)/mu(m,l-1,l-1,0)] * Y(l-1,m) + mu(m,l,l,0) * S(l,m)
 *
 *   (eob_ringdown/py/spheroidal.py's own spherical_from_spheroidal,
 *   l'<=l truncation, and kerrorbits/ringdown.py's own
 *   get_mixed_ringdown_waveform -- same formula). Array-based, NOT
 *   closed-form: S(l,m) and Y(l-1,m) decay/oscillate at DIFFERENT QNM
 *   rates, so their weighted sum has no simple constant-phase-shift
 *   representation the way a single scaled mode does (see the base-mode
 *   loop in eob_wav_ringdown_A22) -- amp/phase are extracted via atan2
 *   and numerically unwrap()ped, same as spheroidal.py's own
 *   create_wdict.
 *
 *   S(l,m)'s own template (ca1..sigmai) is built by the CALLER
 *   (A22_attach_mixed_mode_with_target) from a target point that already
 *   properly mixes (l-1,m)'s own SPHERICAL (nu,spin)-fit target with
 *   (l,m)'s own SPHERICAL (nu,spin)-fit target -- see A22_mix_target_
 *   at_point there -- so ca1..sigmai already match VALUE *and* DERIVATIVE
 *   of the true target at x=0, not just the value (an earlier version of
 *   this function instead rescaled a raw, unmixed S(l,m) template by a
 *   single constant "mu_eff" solved to match the value alone --
 *   value-continuous but NOT derivative-continuous, producing a visible
 *   frequency kink right at attachment; verified directly, e.g.
 *   SXS:BBH:1155's (3,2): M*omega jumps ~0.31->0.38 in one sample at t0
 *   with that approach; a later version instead used S(l,m)'s own
 *   QNMHybridFitCab_A22_lm(k,...) output DIRECTLY, un-mixed -- wrong once
 *   it became clear the (nu,spin) fit for A0/dA0/d2A0/omg0/domg0 targets
 *   the SPHERICAL mode, not the spheroidal one; only cAmp/cphf are
 *   genuinely spheroidal-basis, see A22_attach_mixed_mode_with_target's
 *   own docstring). What THIS function still must do is apply
 *   phi_own_rot: the template itself is built in
 *   eob_wav_ringdown_template_A22's own "phase(0)=0" gauge (matching only
 *   the target's magnitude/frequency, via ca1..cb3), so its own ABSOLUTE
 *   phase at x=0 still needs a separate rotation to match the target's own
 *   phase there -- exactly kerrorbits' own phases_S[lm]
 *   (rm2pi_from_phase(outSlm['p'])), applied the same way in
 *   get_mixed_ringdown_waveform's own hSlm = ...*exp(-1j*phases_Slm[lm]).
 *
 *   @param[in]  k                 : KMAX index of mode (l,m) being attached
 *   @param[in]  l,m               : this mode's own (l,m)
 *   @param[in]  af                : remnant spin (A22_mu's own argument)
 *   @param[in]  Yprev_re,Yprev_im : already-attached Y(l-1,m) real/imag,
 *                                   arrays of length nring over [idx,size),
 *                                   in Y(l-1,m)'s own REAL (Deltaphi-
 *                                   corrected, physically continuous with
 *                                   the inspiral) phase gauge -- matching
 *                                   kerrorbits' own phases_Ylm[(l-1,m)]
 *   @param[in]  nring             : size-idx
 *   @param[in]  t,ooMbh,t0,idx,size : as in eob_wav_ringdown_A22
 *   @param[in]  ca1..sigmai       : S(l,m)'s own template coefficients,
 *                                   built by the caller from the properly
 *                                   mixed target (see above), NOT
 *                                   QNMHybridFitCab_A22_lm(k,...)'s raw
 *                                   (unmixed) output
 *   @param[in]  phi_own_rot       : S(l,m)'s own target ABSOLUTE phase at
 *                                   x=0 (A22_mix_target_at_point's own
 *                                   phi_out) -- the extra rotation
 *                                   ca1..cb3's "phase(0)=0" gauge needs,
 *                                   see above
 *   @param[in]  mu_diag_lp_re/im  : mu(m,l-1,l-1,0), already computed for
 *                                   Y(l-1,m) by the caller (base-mode loop)
 *   @param[out] hlm               : phase[k]/ampli[k] filled for j in [idx,size)
 */
static void A22_attach_mixed_mode(int k, int l, int m, double af,
                                  const double *Yprev_re, const double *Yprev_im, int nring,
                                  double *t, double ooMbh, double t0, int idx, int size,
                                  double ca1, double ca2, double ca3, double ca4, double ca5,
                                  double cb1, double cb2, double cb3, double sigmar, double sigmai,
                                  double phi_own_rot,
                                  double mu_diag_lp_re, double mu_diag_lp_im,
                                  Waveform_lm *hlm)
{
  double mu_cross_re, mu_cross_im;  /* mu(m,l,l-1,0) */
  double mu_own_re,   mu_own_im;    /* mu(m,l,l,0)   */
  A22_mu(m, l, l-1, af, &mu_cross_re, &mu_cross_im);
  A22_mu(m, l,  l , af, &mu_own_re,   &mu_own_im);

  /* cPrev = mu_cross / mu_diag_lp  (complex division) */
  const double denom_abs2 = mu_diag_lp_re*mu_diag_lp_re + mu_diag_lp_im*mu_diag_lp_im;
  const double cPrev_re = (mu_cross_re*mu_diag_lp_re + mu_cross_im*mu_diag_lp_im)/denom_abs2;
  const double cPrev_im = (mu_cross_im*mu_diag_lp_re - mu_cross_re*mu_diag_lp_im)/denom_abs2;

  const double cos_rot = cos(phi_own_rot), sin_rot = sin(phi_own_rot);

  double *Yre = (double*) malloc(nring*sizeof(double));
  double *Yim = (double*) malloc(nring*sizeof(double));
  double *Yph = (double*) malloc(nring*sizeof(double));

  for (int j = idx; j < size; j++) {
    const int jj = j - idx;
    double psi[2];
    eob_wav_ringdown_template_A22(t[j]*ooMbh - t0, ca1,ca2,ca3,ca4,ca5,cb1,cb2,cb3,sigmar,sigmai, psi);
    /* S_raw = amp*exp(-i*phase), same h=A*exp(-i*phi) convention as
       everywhere else in TEOBResumS (and spheroidal.py's own
       create_wdict) -- in the template's OWN "phase(0)=0" gauge. */
    const double Sraw_re =  psi[0]*cos(psi[1]);
    const double Sraw_im = -psi[0]*sin(psi[1]);
    /* rotate by exp(-i*phi_own_rot) to the target's own absolute phase */
    const double S_re = Sraw_re*cos_rot + Sraw_im*sin_rot;
    const double S_im = -Sraw_re*sin_rot + Sraw_im*cos_rot;
    const double z2_re = mu_own_re*S_re - mu_own_im*S_im;
    const double z2_im = mu_own_re*S_im + mu_own_im*S_re;
    const double zz1_re = cPrev_re*Yprev_re[jj] - cPrev_im*Yprev_im[jj];
    const double zz1_im = cPrev_re*Yprev_im[jj] + cPrev_im*Yprev_re[jj];
    Yre[jj] = zz1_re + z2_re;
    Yim[jj] = zz1_im + z2_im;
    Yph[jj] = -atan2(Yim[jj], Yre[jj]);
  }
  unwrap(Yph, nring);

  /* TEOBResumS's own hlm->ampli[k] is a SIGNED quantity for some (l,m)
     (e.g. (3,2)): the mode's own "canonical" sign, carried unchanged
     throughout inspiral, can be negative (verified directly: -0.15ish
     for (3,2) at this attachment point, for a q=2/chi1=0.3/chi2=-0.1
     config) -- but Yre/Yim above always yields amp=sqrt(Yre^2+Yim^2)>=0
     (a POSITIVE-modulus gauge). Writing that positive modulus directly
     into hlm->ampli[k], as this used to do, silently flips the array's
     own sign convention exactly at attachment -- z=A*exp(-i*phi) itself
     stays perfectly continuous (a POSITIVE amp with phase+pi represents
     the SAME complex value as a NEGATIVE amp with phase), so this was
     invisible in any complex-valued check, but any consumer that treats
     "ampli" as its own continuous real time series (a downstream spline
     resample, e.g. EOBRunPy's own interp_uniform_grid, or a naive
     plot of A(t)) sees a genuine sign flip -- interpolating straight
     through zero and producing a large, spurious dip/spike exactly at
     the join. Fix: pick whichever of the two equivalent (ampli,phase)
     representations continues the sign hlm->ampli[k][idx] ALREADY had
     (read here before being overwritten below), applying the
     compensating +pi phase shift the OTHER representation would need
     BEFORE the Dphi continuity correction (so Dphi still exactly
     matches hlm->phase[k][idx], just in the now sign-consistent gauge). */
  const int neg_amp = (hlm->ampli[k][idx] < 0.0);
  if (neg_amp) {
    for (int jj = 0; jj < nring; jj++) Yph[jj] += M_PI;
  }
  const double Dphi = Yph[0] - hlm->phase[k][idx];
  const double amp_sign = neg_amp ? -1.0 : 1.0;
  for (int j = idx; j < size; j++) {
    const int jj = j - idx;
    hlm->ampli[k][j] = amp_sign * sqrt(Yre[jj]*Yre[jj] + Yim[jj]*Yim[jj]);
    hlm->phase[k][j] = Yph[jj] - Dphi;
  }
  free(Yre); free(Yim); free(Yph);
}

/**
 * Function: eob_wav_ringdown_A22
 * --------------------------------
 *   Ringdown calculation and match to the dynamics for the new multi-
 *   mode model (EOBPars->ringdown_model == RINGDOWN_A22). t0 is located
 *   exactly as in the old model (eob_A22_find_t0_raw -- Omega-orbital-
 *   peak minus eob_nqc_timeshift, no extra +2M) and is SHARED by every
 *   mode -- "usual" for this kind of peak-anchored template, no
 *   per-mode dtmrg offset like the old model's HM path.
 *
 *   (2,1),(2,2),(3,3),(4,4),(5,5) need no spheroidal-spherical mixing --
 *   their (nu,spin) global-fit surfaces (QNMHybridFitCab_A22_lm) are
 *   calibrated directly at the TRUE SPHERICAL |h_lm| peak (kind_pk=
 *   'A22', see that function's own docstring), so the template built
 *   from them IS already the spherical mode: attached AS-IS, with NO
 *   mu(m,l,l,0) rescaling of their own amplitude/phase (mu(l,l,0) is
 *   only meaningful for the SPHEROIDAL S_lm piece that feeds INTO a
 *   mixing sum -- applying it to an already-spherical template would be
 *   a spurious double-transform: verified this actually happens once
 *   real, non-unit Berti-Klein mu(l,l,l,0) values are used, at the
 *   ~0.1-1% amplitude level, silently breaking the (2,2)-exactness this
 *   whole model's validation history depends on). mu_diag_re/im[k] is
 *   still computed here for k22/k33 -- needed downstream by (3,2)'s/
 *   (4,3)'s own mixing (both the template attachment and the NQC
 *   target, see A22_attach_mixed_mode / A22_mix_target_at_point) -- just
 *   not applied to k's OWN output. (3,2) and (4,3) themselves ARE
 *   genuinely spheroidal (no direct spherical fit) and get MIXED with
 *   the already-attached (2,2)/(3,3) via A22_attach_mixed_mode -- see
 *   this file's own header comment and that function's docstring.
 *
 *   errorexit()s if any active mode is outside the supported set, or if
 *   (3,2)/(4,3) is active without (2,2)/(3,3) (needed for its own
 *   mixing).
 *
 *   @param[in]  dyn : dynamics
 *   @param[out] hlm : waveform
 *   @return         : 0 if successful
 */

/**
 * Function: A22_mix_target_at_point
 * -----------------------------------
 *   Combine two independently-known analytic complex pieces
 *   z1 = A0_1*exp(-i*[phi0_1+phi1(x)]), z2 = A0_2*exp(-i*[phi0_2+phi2(x)])
 *   -- each given by its own (A0,dA0,d2A0,omega0,domega0) in the gauge
 *   where its OWN "shape" phase is 0 at x=0 (the natural gauge
 *   QNMHybridFitCab_A22_lm's own A0_out/dA0_out/d2A0_out/omg0_out/
 *   domg0_out are already in, and the SAME gauge a natural, un-fit
 *   inspiral trajectory's own finite-difference derivatives sit in too),
 *   PLUS its own separate absolute phase offset phi0_i (c1,c2 below are
 *   PRE-ROTATED by exp(-i*phi0_i) so this function itself only ever sees
 *   the phi0=0 gauge quantities -- see eob_wav_ringdown_A22's own
 *   callers) -- into the (A,dA,d2A,omega,phi) of Z = c1*z1 + c2*z2 at
 *   x=0, for arbitrary complex constants c1,c2.
 *
 *   Used by A22_attach_mixed_mode_with_target to build S(l,m)'s own
 *   PROPERLY MIXED target: inverting Y(l,m) = cPrev*Y(l-1,m) + mu_own*
 *   S(l,m) for S(l,m), i.e. z1 = Y(l-1,m)'s own SPHERICAL (nu,spin)-fit
 *   target (c1 = -cPrev/mu_own, rotated by Y(l-1,m)'s own REAL absolute
 *   phase at t0) and z2 = Y(l,m) itself, i.e. THIS mode's own SPHERICAL
 *   target too (c2 = 1/mu_own, rotated by ITS OWN real absolute phase at
 *   t0) -- exactly kerrorbits/ringdown.py's own get_ringdown_waveforms_
 *   spheroidal/store(): h1dict/h2dict there are built the identical way
 *   (add_complex_constant_contribution, phi0=phases_Y[...]) from
 *   nqcYlm[lm]/nqcYlm[(l-1,m)] (both SPHERICAL-basis NQC targets), then
 *   sum_of_three_modes combines them into outS['A'/'dA'/'d2A'/'omg'],
 *   which store() writes directly as the S(l,m) template's own A0/dA0/
 *   d2A0/delta_omg0 -- i.e. this function is that same combination, and
 *   its A_out/dA_out/d2A_out/omg_out feed A22_build_template_coeffs the
 *   same way. phi_out is the combination's own absolute phase (matching
 *   kerrorbits' phases_S[lm] = rm2pi_from_phase(outSlm['p'])) -- since
 *   eob_wav_ringdown_template_A22 always builds S(l,m) in ITS OWN
 *   "phase(0)=0" gauge (matching only A_out/dA_out/omg_out), phi_out is
 *   the SEPARATE rotation A22_attach_mixed_mode's own phi_own_rot
 *   parameter needs, exactly like kerrorbits' own hSlm=...*
 *   exp(-1j*phases_Slm[lm]) in get_mixed_ringdown_waveform.
 *
 *   Derivation: z(0)=A0, z'(0)=dA0-i*A0*omega0,
 *   z''(0)=(d2A0-A0*omega0^2)-i*(2*dA0*omega0+A0*domega0), from
 *   d/dx[A(x)e^{-i phi(x)}] with phi(0)=0; derivatives are linear so
 *   Z=c1 z1+c2 z2's own derivatives are the same combination of the
 *   pieces'; then A=|Z|, phi=atan2(-Im(Z),Re(Z)), w=Z'*conj(Z),
 *   dA=Re(w)/A, omega=-Im(w)/A^2, d2A=[Re(Z''*conj(Z))+|Z'|^2-dA^2]/A --
 *   this last relation (and A/dA/omega) verified numerically against
 *   finite-difference differentiation of the directly-summed signal
 *   (matched to 6-7 significant figures).
 */
static void A22_mix_target_at_point(double A0_1, double dA0_1, double d2A0_1, double omg0_1, double domg0_1,
                                    double c1_re, double c1_im,
                                    double A0_2, double dA0_2, double d2A0_2, double omg0_2, double domg0_2,
                                    double c2_re, double c2_im,
                                    double *A_out, double *dA_out, double *d2A_out, double *omg_out, double *phi_out)
{
  const double z1_re = A0_1,  z1_im = 0.0;
  const double z1p_re = dA0_1, z1p_im = -A0_1*omg0_1;
  const double z1pp_re = d2A0_1 - A0_1*omg0_1*omg0_1;
  const double z1pp_im = -(2.*dA0_1*omg0_1 + A0_1*domg0_1);

  const double z2_re = A0_2,  z2_im = 0.0;
  const double z2p_re = dA0_2, z2p_im = -A0_2*omg0_2;
  const double z2pp_re = d2A0_2 - A0_2*omg0_2*omg0_2;
  const double z2pp_im = -(2.*dA0_2*omg0_2 + A0_2*domg0_2);

#define A22_CMUL_RE(ar,ai,br,bi) ((ar)*(br)-(ai)*(bi))
#define A22_CMUL_IM(ar,ai,br,bi) ((ar)*(bi)+(ai)*(br))

  const double Z_re   = A22_CMUL_RE(c1_re,c1_im,z1_re,z1_im)   + A22_CMUL_RE(c2_re,c2_im,z2_re,z2_im);
  const double Z_im   = A22_CMUL_IM(c1_re,c1_im,z1_re,z1_im)   + A22_CMUL_IM(c2_re,c2_im,z2_re,z2_im);
  const double Zp_re  = A22_CMUL_RE(c1_re,c1_im,z1p_re,z1p_im) + A22_CMUL_RE(c2_re,c2_im,z2p_re,z2p_im);
  const double Zp_im  = A22_CMUL_IM(c1_re,c1_im,z1p_re,z1p_im) + A22_CMUL_IM(c2_re,c2_im,z2p_re,z2p_im);
  const double Zpp_re = A22_CMUL_RE(c1_re,c1_im,z1pp_re,z1pp_im) + A22_CMUL_RE(c2_re,c2_im,z2pp_re,z2pp_im);
  const double Zpp_im = A22_CMUL_IM(c1_re,c1_im,z1pp_re,z1pp_im) + A22_CMUL_IM(c2_re,c2_im,z2pp_re,z2pp_im);

#undef A22_CMUL_RE
#undef A22_CMUL_IM

  const double A = sqrt(Z_re*Z_re + Z_im*Z_im);
  *A_out   = A;
  *phi_out = atan2(-Z_im, Z_re);

  if (A == 0.0) {
    /* Both pieces are identically 0 (e.g. (4,3) mixing with (3,3) when
       BOTH are still placeholder fits) -- amplitude derivative/d2A/omega
       are 0/0 here. They multiply/shift a target amplitude that is
       itself exactly 0, so any finite value would do; 0 is the simplest,
       and (unlike NaN) doesn't poison A22_build_template_coeffs below. */
    *dA_out = 0.0; *d2A_out = 0.0; *omg_out = 0.0;
    return;
  }

  const double w_re  = Zp_re*Z_re + Zp_im*Z_im;         /* Re(Z' conj(Z)) */
  const double w_im  = Zp_im*Z_re - Zp_re*Z_im;         /* Im(Z' conj(Z)) */
  const double dA    = w_re/A;
  const double omg   = -w_im/(A*A);
  const double wp_re = (Zpp_re*Z_re + Zpp_im*Z_im) + (Zp_re*Zp_re + Zp_im*Zp_im); /* Re(Z''conj(Z)) + |Z'|^2 */
  const double d2A   = (wp_re - dA*dA)/A;

  *dA_out  = dA;
  *d2A_out = d2A;
  *omg_out = omg;
}

/**
 * Function: A22_attach_mixed_mode_with_target
 * -----------------------------------------------
 *   Build S(l,m)'s own properly-mixed target and template coefficients,
 *   then attach (l,m) via A22_attach_mixed_mode -- the C-code analogue of
 *   kerrorbits/ringdown.py's own get_imr_multipoles(): (1) EVERY mode,
 *   base or mixed, gets NQC-corrected on its own SPHERICAL basis first
 *   (kerrorbits' get_nqc_multipole, using nqc_point_vals_A22pk -- here,
 *   eob_wav_hlmNQC_find_a1a2a3_mrg_A22, using QNMHybridFitCab_A22_lm(k,...)'s
 *   own A0/dA0/d2A0/omg0/domg0 target directly, now that the (nu,spin) fit
 *   surface itself is calibrated on the SPHERICAL decomposition -- see
 *   comparable_mass_fitter.py's run_primary_fit_multistart, dual-fit
 *   change: A0/dA0/... always spherical, cAmp/cphf always spheroidal);
 *   (2) S(l,m)'s own target is recovered by INVERTING Y(l,m)=cPrev*
 *   Y(l-1,m)+mu_own*S(l,m) at x=0 using Y(l-1,m)'s and Y(l,m)'s own
 *   (NQC-corrected, spherical) targets (kerrorbits' get_ringdown_
 *   waveforms_spheroidal/store(), combining nqcYlm[lm]/nqcYlm[(l-1,m)] via
 *   sum_of_three_modes) -- A22_mix_target_at_point, VALUE *and*
 *   DERIVATIVE, not just value (an intermediate version of this function
 *   used only the value, producing a visible frequency kink right at
 *   attachment: SXS:BBH:1155's (3,2) M*omega jumping ~0.31->0.38 in one
 *   sample at t0); (3) S(l,m)'s template is built from this target using
 *   its OWN shape (cAmp/cphf, genuinely spheroidal, from THIS SAME
 *   QNMHybridFitCab_A22_lm(k,...) call -- kerrorbits' get_pfit_spheroidal);
 *   (4) the actual attached array is the array-level SUM cPrev*Y(l-1,m)(t)
 *   + mu_own*S(l,m)(t) (A22_attach_mixed_mode, kerrorbits' own
 *   get_mixed_ringdown_waveform) -- unchanged by this function.
 *
 *   (4,3) has no dedicated (nu,spin) fit yet (hatA0_l4m3 doesn't exist --
 *   A22_get_fit_set falls back to A22_fit_placeholder for k43, detected
 *   here via cA2==0.0): falls back to a local finite-difference of its own
 *   untouched (NQC stays OFF for it, see eob_wav_hlmNQC_find_a1a2a3_mrg_A22)
 *   natural inspiral trajectory instead, same as an earlier version of
 *   this function did for (3,2) too before hatA0_l3m2 existed.
 *
 *   @param[in] k,l,m   : mode (l,m) being attached (KMAX index k)
 *   @param[in] kprev   : KMAX index of Y(l-1,m), already attached
 *   @param[in] Yprev_re,Yprev_im,nring : Y(l-1,m)'s own buffer (real
 *     absolute phase gauge), as built by the base-mode loop
 *   @param[in] mu_diag_lp_re/im : mu(m,l-1,l-1,0), already computed by
 *     the caller
 */
static void A22_attach_mixed_mode_with_target(
    int k, int l, int m, int kprev,
    double nu, double X1, double X2, double chi1, double chi2, double Mbh, double abh,
    const double *Yprev_re, const double *Yprev_im, int nring,
    double *t, double ooMbh, double t0, int idx, int size,
    double mu_diag_lp_re, double mu_diag_lp_im,
    Waveform_lm *hlm)
{
  /* S(l,m)'s own SHAPE (cA2,cA3,d1f,d2f) and QNM (sigmar,sigmai), plus
     (l,m)'s own SPHERICAL point-value target (A0_nat..domg0_nat) -- ALL
     from the SAME (nu,spin) fit surface: the shape feeds S(l,m)'s
     template directly (step 3 above); the point-value target is the
     mode's own "natural" (Y(l,m)) piece in the mixing (step 2), replacing
     eob_wav_hlmNQC_find_a1a2a3_mrg_A22's role of applying it AS a
     correction to the array (which A22_attach_mixed_mode's own summing
     loop bypasses for the mixed-mode output anyway -- see its docstring).
     Placeholder (no real fit, e.g. (4,3) currently) is detected via
     cA2==0.0 (A22_build_template_coeffs's own guard) and handled by
     falling back to a local finite-difference of the untouched array. */
  double ca1_unused,ca2_unused,ca3_unused,ca4_unused,ca5_unused,
        cb1_unused,cb2_unused,cb3_unused,sigmar,sigmai;
  double cA2,cA3,d1f,d2f;
  double A0_nat,dA0_nat,d2A0_nat,omg0_nat,domg0_nat;
  QNMHybridFitCab_A22_lm(k, nu, X1, X2, chi1, chi2, Mbh, abh,
                         &ca1_unused,&ca2_unused,&ca3_unused,&ca4_unused,&ca5_unused,
                         &cb1_unused,&cb2_unused,&cb3_unused,&sigmar,&sigmai,
                         &A0_nat,&dA0_nat,&d2A0_nat,&omg0_nat,&domg0_nat,
                         &cA2,&cA3,&d1f,&d2f);

  /* TEOBResumS's own standard (pre-existing, un-mixed) multipolar
     waveform construction puts Y(l,m)'s OWN absolute phase an extra pi
     out of phase with Y(l-1,m)'s, relative to what real NR shows --
     verified directly for (3,2)/(2,2) on SXS:BBH:1155: tracking
     phi32_EOB-phi22_EOB (TEOBResumS's own raw, pre-NQC, pre-attachment
     trajectory -- so this is NOT an NQC/mixing/attachment artifact, it is
     already present ~100M before merger) against phi32_NR-phi22_NR (real
     NR, both same "+angle" sign convention once TEOBResumS's own
     "-angle" hlm->phase is corrected for), their SUM stays close to -pi
     (const to ~0.02-0.08 rad) across the entire late inspiral -- i.e.
     Y(l,m)'s own phase needs a +pi correction (equivalently, an overall
     sign flip: h_lm -> -h_lm) for its OWN relative phase against Y(l-1,m)
     to agree with NR's -- exactly the (-1)^(l+m) parity BOTH tower
     members implicitly need relative to each other (l+m is ALWAYS odd
     for l=m+1, the only mixing this file implements, so this is
     unconditional here, but written as a parity check for correctness/
     documentation, matching the general (-1)^(l+m) rule -- NOT
     (-1)^l or (-1)^m alone, which would each get one of (2,1)/(3,2)/
     (4,3)/(5,4) wrong). Y(l-1,m) itself (phi0_prev) needs NO such
     correction (l-1+m is even for these pairs). */
  double phi0_nat = hlm->phase[k][idx];
  if ((l+m) % 2 != 0) phi0_nat += M_PI;
  if (cA2 == 0.0) {
    /* No dedicated fit yet -- fall back to the untouched (NQC-OFF) array's
       own local finite-difference derivative, in PLAIN (t, M=1) units (NOT
       Mbh-rescaled "x" units -- QNMHybridFitCab_A22_lm's own dA0_out/
       omg0_out are plain t-derivatives, see its own docstring), matching
       kerrorbits' own phases_Y[lm]/nqcYlm[lm] mechanism minus the
       dedicated fit this mode doesn't have. Falls back further to a
       narrower stencil (no domega/d2A) near the array edges, then to a
       flat (zero-derivative) target -- idx is always deep inside a
       physically reasonable dynamics array in practice, this is
       defensive only. */
    A0_nat = hlm->ampli[k][idx];
    dA0_nat = d2A0_nat = omg0_nat = domg0_nat = 0.0;
    if (idx >= 2 && idx + 2 < size) {
      const double dt_loc = (t[idx+1] - t[idx-1]) * 0.5;
      const double Am1 = hlm->ampli[k][idx-1], Ap1 = hlm->ampli[k][idx+1];
      const double Pm2 = hlm->phase[k][idx-2], Pm1 = hlm->phase[k][idx-1];
      const double Pp1 = hlm->phase[k][idx+1], Pp2 = hlm->phase[k][idx+2];
      dA0_nat   = (Ap1 - Am1) / (2.*dt_loc);
      d2A0_nat  = (Ap1 - 2.*A0_nat + Am1) / (dt_loc*dt_loc);
      const double omg_m1 = (hlm->phase[k][idx] - Pm2) / (2.*dt_loc);
      const double omg_p1 = (Pp2 - hlm->phase[k][idx]) / (2.*dt_loc);
      omg0_nat  = (Pp1 - Pm1) / (2.*dt_loc);
      domg0_nat = (omg_p1 - omg_m1) / (2.*dt_loc);
    } else if (idx >= 1 && idx + 1 < size) {
      const double dt_loc = (t[idx+1] - t[idx-1]) * 0.5;
      dA0_nat  = (hlm->ampli[k][idx+1] - hlm->ampli[k][idx-1]) / (2.*dt_loc);
      d2A0_nat = (hlm->ampli[k][idx+1] - 2.*A0_nat + hlm->ampli[k][idx-1]) / (dt_loc*dt_loc);
      omg0_nat = (hlm->phase[k][idx+1] - hlm->phase[k][idx-1]) / (2.*dt_loc);
    }
  }

  /* Y(l-1,m)'s own SPHERICAL target (always a real fit -- (2,2)/(3,3)) and
     its own REAL absolute phase at idx (the base-mode loop's own
     Deltaphi-corrected value -- matching kerrorbits' phases_Y[(l-1,m)]). */
  double ca1p,ca2p,ca3p,ca4p,ca5p,cb1p,cb2p,cb3p,sigmarp,sigmaip;
  double A0_prev, dA0_prev, d2A0_prev, omg0_prev, domg0_prev;
  QNMHybridFitCab_A22_lm(kprev, nu, X1, X2, chi1, chi2, Mbh, abh,
                         &ca1p,&ca2p,&ca3p,&ca4p,&ca5p,&cb1p,&cb2p,&cb3p,&sigmarp,&sigmaip,
                         &A0_prev,&dA0_prev,&d2A0_prev,&omg0_prev,&domg0_prev, NULL,NULL,NULL,NULL);
  const double phi0_prev = hlm->phase[kprev][idx];

  /* mu coefficients */
  double mu_cross_re, mu_cross_im, mu_diag_re, mu_diag_im, mu_own_re, mu_own_im;
  A22_mu(m, l, l-1, abh, &mu_cross_re, &mu_cross_im);
  A22_mu(m, l-1, l-1, abh, &mu_diag_re, &mu_diag_im);
  A22_mu(m, l, l, abh, &mu_own_re, &mu_own_im);
  (void) mu_diag_lp_re; (void) mu_diag_lp_im; /* mu_diag recomputed locally, same as caller's own -- kept as a parameter for symmetry with A22_attach_mixed_mode's own signature */
  const double mu_diag_abs2 = mu_diag_re*mu_diag_re + mu_diag_im*mu_diag_im;
  const double cPrev_re = (mu_cross_re*mu_diag_re + mu_cross_im*mu_diag_im)/mu_diag_abs2;
  const double cPrev_im = (mu_cross_im*mu_diag_re - mu_cross_re*mu_diag_im)/mu_diag_abs2;
  const double mu_own_abs2 = mu_own_re*mu_own_re + mu_own_im*mu_own_im;

  /* Invert Y(l,m) = cPrev*Y(l-1,m) + mu_own*S(l,m) for S(l,m):
     c1 = (-cPrev/mu_own)*exp(-i*phi0_prev) [piece 1 = Y(l-1,m)'s own target]
     c2 = (1/mu_own)*exp(-i*phi0_nat)       [piece 2 = Y(l,m)'s own natural target] */
  const double negcPrev_over_muown_re = -(cPrev_re*mu_own_re + cPrev_im*mu_own_im)/mu_own_abs2;
  const double negcPrev_over_muown_im = -(cPrev_im*mu_own_re - cPrev_re*mu_own_im)/mu_own_abs2;
  const double inv_muown_re =  mu_own_re/mu_own_abs2;
  const double inv_muown_im = -mu_own_im/mu_own_abs2;

  const double cos_p = cos(phi0_prev), sin_p = sin(phi0_prev);
  const double c1_re = negcPrev_over_muown_re*cos_p + negcPrev_over_muown_im*sin_p;
  const double c1_im = -negcPrev_over_muown_re*sin_p + negcPrev_over_muown_im*cos_p;

  const double cos_n = cos(phi0_nat), sin_n = sin(phi0_nat);
  const double c2_re = inv_muown_re*cos_n + inv_muown_im*sin_n;
  const double c2_im = -inv_muown_re*sin_n + inv_muown_im*cos_n;

  double A_S, dA_S, d2A_S, omega_S, phi_S;
  A22_mix_target_at_point(A0_prev, dA0_prev, d2A0_prev, omg0_prev, domg0_prev, c1_re, c1_im,
                          A0_nat,  dA0_nat,  d2A0_nat,  omg0_nat,  domg0_nat,  c2_re, c2_im,
                          &A_S, &dA_S, &d2A_S, &omega_S, &phi_S);

  /* S(l,m)'s own template from this properly mixed target, using ITS OWN
     shape (cA2,cA3,d1f,d2f) and QNM (sigmar,sigmai). */
  double ca1,ca2,ca3,ca4,ca5,cb1,cb2,cb3;
  A22_build_template_coeffs(A_S, dA_S, d2A_S, sigmar, sigmai, Mbh, omega_S, cA2, cA3, d1f, d2f,
                            &ca1,&ca2,&ca3,&ca4,&ca5,&cb1,&cb2,&cb3);

  A22_attach_mixed_mode(k, l, m, abh, Yprev_re, Yprev_im, nring,
                        t, ooMbh, t0, idx, size,
                        ca1,ca2,ca3,ca4,ca5,cb1,cb2,cb3,sigmar,sigmai,
                        phi_S,
                        mu_diag_re, mu_diag_im, hlm);
}

int eob_wav_ringdown_A22(Dynamics *dyn, Waveform_lm *hlm)
{
  const double Mbh   = EOBPars->Mbhf;
  const double abh   = EOBPars->abhf;
  const double nu    = EOBPars->nu;
  const double chi1  = EOBPars->chi1;
  const double chi2  = EOBPars->chi2;
  const double X1    = EOBPars->X1;
  const double X2    = EOBPars->X2;
  const double ooMbh = 1./Mbh;

  const int size = hlm->size;
  double *t = hlm->time;

  int mode_supported[KMAX];
  for (int k = 0; k < KMAX; k++) mode_supported[k] = 0;
  mode_supported[A22_K21] = mode_supported[A22_K22] = mode_supported[A22_K32] =
    mode_supported[A22_K33] = mode_supported[A22_K43] = mode_supported[A22_K44] =
    mode_supported[A22_K55] = 1;

  for (int k = 0; k < KMAX; k++) {
    if (hlm->kmask[k] && !mode_supported[k]) {
      errorexit("eob_wav_ringdown_A22: mode not yet supported by "
               "ringdown_model=\"A22\" -- only (2,1),(2,2),(3,2),(3,3),(4,3),"
               "(4,4),(5,5) are implemented.\n");
    }
  }
  if (hlm->kmask[A22_K32] && !hlm->kmask[A22_K22]) {
    errorexit("eob_wav_ringdown_A22: (3,2) needs (2,2) also active (its "
             "spheroidal-spherical mixing uses the already-attached (2,2) "
             "spherical mode).\n");
  }
  if (hlm->kmask[A22_K43] && !hlm->kmask[A22_K33]) {
    errorexit("eob_wav_ringdown_A22: (4,3) needs (3,3) also active (its "
             "spheroidal-spherical mixing uses the already-attached (3,3) "
             "spherical mode).\n");
  }

  /* t0, in Mbh units (matching sigma/alpha/ca-cb's own natural scale) --
     SAME anchor for every mode. */
  const double t0raw = eob_A22_find_t0_raw(dyn, nu, chi1);
  const double t0 = t0raw * ooMbh;

  /* Attachment index: NEAREST sample to t0 (RAW units) -- SAME for every
     mode, and CRITICALLY the same "nearest sample" convention (not
     "first sample at/after") that
     eob_wav_hlmNQC_find_a1a2a3_mrg_A22's own jmax search uses for this
     SAME t0 -- verified these two searches, run on hlm's vs hlm_mrg's
     own (overlapping but separately-allocated) time arrays, can land on
     samples up to one dt_merger_interp apart if the conventions don't
     match (a "first at/after" vs "nearest" pick straddling the true t0
     differently), which showed up as a real (though modest, O(dt))
     mismatch between the NQC-fit's own matching instant and the
     ringdown-attach's own instant -- amplified into a large, visible
     discontinuity for (3,2) specifically, since its spheroidal-spherical
     mixing makes its amplitude vary rapidly (a real "notch" feature)
     within just a few M of merger, unlike (2,2)'s comparatively smooth,
     slowly-varying peak. */
  int idx = 0;
  {
    double best = fabs(t[0] - t0raw);
    for (int j = 1; j < size; j++) {
      const double d = fabs(t[j] - t0raw);
      if (d < best) { best = d; idx = j; }
    }
  }
  const int nring = size - idx;

  /* Base ("diagonal", no-mixing) modes: (2,1),(2,2),(3,3),(4,4),(5,5).
     Closed-form scale+constant-phase-shift attach -- also builds the
     Y(2,2)/Y(3,3) buffers (3,2)/(4,3) need for their own mixing below. */
  const int base_modes[5] = {A22_K21, A22_K22, A22_K33, A22_K44, A22_K55};
  double mu_diag_re[KMAX], mu_diag_im[KMAX];
  double *Y22_re = NULL, *Y22_im = NULL, *Y33_re = NULL, *Y33_im = NULL;
  if (hlm->kmask[A22_K32]) {
    Y22_re = (double*) malloc(nring*sizeof(double));
    Y22_im = (double*) malloc(nring*sizeof(double));
  }
  if (hlm->kmask[A22_K43]) {
    Y33_re = (double*) malloc(nring*sizeof(double));
    Y33_im = (double*) malloc(nring*sizeof(double));
  }

  for (int i = 0; i < 5; i++) {
    const int k = base_modes[i];
    if (!hlm->kmask[k]) continue;

    double ca1,ca2,ca3,ca4,ca5,cb1,cb2,cb3,sigmar,sigmai;
    QNMHybridFitCab_A22_lm(k, nu, X1, X2, chi1, chi2, Mbh, abh,
                           &ca1,&ca2,&ca3,&ca4,&ca5,&cb1,&cb2,&cb3,&sigmar,&sigmai,
                           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    const int l = LINDEX[k], m = MINDEX[k];
    double mu_re, mu_im;
    A22_mu(m, l, l, abh, &mu_re, &mu_im);
    mu_diag_re[k] = mu_re; mu_diag_im[k] = mu_im; /* NOT applied below -- see docstring */

    double psi[2];
    eob_wav_ringdown_template_A22(t[idx]*ooMbh - t0, ca1,ca2,ca3,ca4,ca5,cb1,cb2,cb3,
                                  sigmar,sigmai, psi);
    const double Deltaphi = psi[1] - hlm->phase[k][idx];

    for (int j = idx; j < size; j++) {
      eob_wav_ringdown_template_A22(t[j]*ooMbh - t0, ca1,ca2,ca3,ca4,ca5,cb1,cb2,cb3,
                                    sigmar,sigmai, psi);
      const double amp_Y   = psi[0];
      const double phase_Y = psi[1] - Deltaphi;
      hlm->ampli[k][j] = amp_Y;
      hlm->phase[k][j] = phase_Y;
      /* Y22/Y33 buffers feeding the (3,2)/(4,3) MIXING below use amp_Y/
         phase_Y -- the REAL, Deltaphi-corrected phase, physically
         continuous with the inspiral (matching kerrorbits/ringdown.py's
         own get_mixed_ringdown_waveform, whose hYlm1m is likewise
         rotated by phases_Ylm[(l-1,m)], the (l-1,m) mode's own already-
         tracked inspiral phase, not an arbitrary "phase=0 at t0" gauge)
         -- NOT the raw psi[1]. A22_attach_mixed_mode's own mu_eff
         inversion is solved AT idx using whatever gauge Y(l-1,m) is in
         here, so it reproduces the correct target regardless -- but the
         SHAPE of the combination away from idx does depend on this
         choice, and the correct one (matching the reference
         architecture) is Y(l-1,m)'s own real absolute phase. */
      if (k == A22_K22 && Y22_re) {
        const int jj = j - idx;
        Y22_re[jj] =  amp_Y*cos(phase_Y);
        Y22_im[jj] = -amp_Y*sin(phase_Y);
      }
      if (k == A22_K33 && Y33_re) {
        const int jj = j - idx;
        Y33_re[jj] =  amp_Y*cos(phase_Y);
        Y33_im[jj] = -amp_Y*sin(phase_Y);
      }
    }
  }

  /* Mixed modes: (3,2) [needs Y22], (4,3) [needs Y33]. */
  if (hlm->kmask[A22_K32]) {
    A22_attach_mixed_mode_with_target(A22_K32, 3, 2, A22_K22, nu, X1, X2, chi1, chi2, Mbh, abh,
                                      Y22_re, Y22_im, nring, t, ooMbh, t0, idx, size,
                                      mu_diag_re[A22_K22], mu_diag_im[A22_K22], hlm);
  }
  if (hlm->kmask[A22_K43]) {
    A22_attach_mixed_mode_with_target(A22_K43, 4, 3, A22_K33, nu, X1, X2, chi1, chi2, Mbh, abh,
                                      Y33_re, Y33_im, nring, t, ooMbh, t0, idx, size,
                                      mu_diag_re[A22_K33], mu_diag_im[A22_K33], hlm);
  }

  if (Y22_re) { free(Y22_re); free(Y22_im); }
  if (Y33_re) { free(Y33_re); free(Y33_im); }

  return 0;
}

/**
 * Function: eob_wav_hlmNQC_find_a1a2a3_mrg_A22
 * -----------------------------------------------
 *   NQC coefficients for the new multi-mode model -- SAME basis-function
 *   / 2x2-linear-solve machinery as eob_wav_hlmNQC_find_a1a2a3_mrg_22,
 *   independently for EVERY active mode EXCEPT (4,3) (see below) (basis
 *   functions n1,n2,n4,n5 and their derivatives depend only on the
 *   dynamics, so are computed ONCE and shared; the target point and the
 *   resulting (ai,bi) are per-mode), but:
 *     - the target point (A0, dA0, omega0, domega0) comes from the SAME
 *       (nu,spin) global-fit surface the ringdown template for that mode
 *       uses (QNMHybridFitCab_A22_lm -- real fit for (2,2)/(3,2),
 *       placeholder 0 for everything else, see this file's own header
 *       comment), NOT eob_nqc_point's separately-calibrated NR fit -- and,
 *       since comparable_mass_fitter.py's run_primary_fit_multistart dual-
 *       fit change, is the SPHERICAL basis for EVERY mode (l==m already
 *       was; (3,2)/(4,3) now are too, unlike cAmp/cphf, which stay
 *       spheroidal-basis -- see A22_attach_mixed_mode_with_target's own
 *       docstring), exactly matching what THIS correction needs (it acts
 *       on hlm->ampli/phase[k], the SPHERICAL multipole);
 *     - the matching time is the SAME t0 every mode's ringdown attaches
 *       at (eob_A22_find_t0_raw, called with the FULL dynamics `dyn` here
 *       too, then located within the merger-segment's own, shorter time
 *       array), NOT the old model's own tOmgOrb_pk-DeltaT_nqc point on
 *       dyn_mrg alone;
 *     - plain-strain (h/nu) units throughout, matching this pipeline's
 *       own A0/dA0/domega0 convention exactly -- NOT eob_nqc_point's RWZ-
 *       normalized convention (nlm factor), which is specific to that
 *       older, separately-calibrated fit.
 *   For the placeholder ((4,3), currently the only one) mode, a zero
 *   target smoothly fades that mode's own (otherwise perfectly normal
 *   PN/factorized-resummed) inspiral-plunge amplitude down as t0 is
 *   approached (the NQC correction is a bounded function of the dynamics,
 *   applied multiplicatively over the WHOLE trajectory, small away from
 *   merger -- same reasoning as the (2,2)-only model's own inspiral-
 *   negligibility argument), rather than leaving an abrupt jump at t0
 *   where the (currently zero) ringdown template takes over.
 *
 *   (3,2) DOES get a genuine NQC correction here, on hlm->ampli/phase[k32]
 *   itself, matching kerrorbits/ringdown.py's own get_nqc_multipole
 *   (nqc_point_vals_A22pk, also spherical-basis) -- this is READ AGAIN by
 *   A22_attach_mixed_mode_with_target, but only its resulting PHASE
 *   (hlm->phase[A22_K32][idx], the physically-meaningful absolute phase
 *   of the NQC-corrected signal -- matching kerrorbits' own phases_Y[lm])
 *   -- the VALUE+DERIVATIVE target that function needs is read DIRECTLY
 *   from QNMHybridFitCab_A22_lm(k32,...)'s own output, not by re-reading
 *   hlm->ampli[A22_K32][idx] (see that function's own docstring).
 *
 *   (4,3) still gets NO NQC correction (ai0=ai1=bi0=bi1=0 identically,
 *   skipped in the loop below): no dedicated fit exists for it yet
 *   (A22_get_fit_set falls back to A22_fit_placeholder for k43) -- its
 *   spherical trajectory is simply left as whatever the standard, mode-
 *   independent factorized-resummed PN inspiral-plunge machinery already
 *   gives it, and A22_attach_mixed_mode_with_target's own cA2==0.0
 *   fallback reads THIS untouched value/phase as Y(4,3)'s own "natural"
 *   target.
 */
void eob_wav_hlmNQC_find_a1a2a3_mrg_A22(Dynamics *dyn_mrg, Waveform_lm *hlm_mrg, Waveform_lm *hnqc,
                                        Dynamics *dyn, Waveform_lm *hlm)
{
  int mode_supported[KMAX];
  for (int k = 0; k < KMAX; k++) mode_supported[k] = 0;
  mode_supported[A22_K21] = mode_supported[A22_K22] = mode_supported[A22_K32] =
    mode_supported[A22_K33] = mode_supported[A22_K43] = mode_supported[A22_K44] =
    mode_supported[A22_K55] = 1;
  for (int k = 0; k < KMAX; k++) {
    if (hlm->kmask[k] && !mode_supported[k]) {
      errorexit("eob_wav_hlmNQC_find_a1a2a3_mrg_A22: mode not yet supported by "
               "ringdown_model=\"A22\".\n");
    }
  }
  if (hlm->kmask[A22_K32] && !hlm->kmask[A22_K22]) {
    errorexit("eob_wav_hlmNQC_find_a1a2a3_mrg_A22: (3,2) needs (2,2) also active.\n");
  }
  if (hlm->kmask[A22_K43] && !hlm->kmask[A22_K33]) {
    errorexit("eob_wav_hlmNQC_find_a1a2a3_mrg_A22: (4,3) needs (3,3) also active.\n");
  }

  const double nu    = EOBPars->nu;
  const double chi1  = EOBPars->chi1;
  const double chi2  = EOBPars->chi2;
  const double X1    = EOBPars->X1;
  const double X2    = EOBPars->X2;
  const double Mbh   = EOBPars->Mbhf;
  const double abh   = EOBPars->abhf;

  double *t       = hlm_mrg->time;
  double *r       = dyn_mrg->data[EOB_RAD];
  double *w       = dyn_mrg->data[EOB_MOMG];
  double *pr_star = dyn_mrg->data[EOB_PRSTAR];
  double *ddotr   = dyn_mrg->data[EOB_DDOTR];

  const int size = hlm_mrg->size;
  const double dt = t[1]-t[0];

  for (int i = 0; i < size; i++) hnqc->time[i] = t[i];

  /* Basis functions and their derivatives -- MODE-INDEPENDENT (functions
     of the dynamics only), computed ONCE and shared by every mode. */
  double *n1    = (double*) calloc(size, sizeof(double));
  double *n2    = (double*) calloc(size, sizeof(double));
  double *n4    = (double*) calloc(size, sizeof(double));
  double *n5    = (double*) calloc(size, sizeof(double));
  double *d_n4  = (double*) calloc(size, sizeof(double));
  double *d_n5  = (double*) calloc(size, sizeof(double));
  double *d2_n4 = (double*) calloc(size, sizeof(double));
  double *d2_n5 = (double*) calloc(size, sizeof(double));

  double pr_star2, r2, w2;
  for (int j=0; j<size; j++) {
    pr_star2 = SQ(pr_star[j]);
    r2       = SQ(r[j]);
    w2       = SQ(w[j]);
    n1[j] = pr_star2/(r2*w2);
    n2[j] = ddotr[j]/(r[j]*w2);
    n4[j] = pr_star[j]/(r[j]*w[j]);
    n5[j] = n4[j]*r2*w2;
  }
  D0(n4, dt, size, d_n4);
  D0(n5, dt, size, d_n5);
  D0(d_n4, dt, size, d2_n4);
  D0(d_n5, dt, size, d2_n5);

  /* Matching time: the SAME t0 the ringdown attaches at (full dynamics,
     RAW units), located within the merger segment's own time array --
     SAME for every mode. */
  const double t0 = eob_A22_find_t0_raw(dyn, nu, chi1);
  int jmax = 0;
  {
    double best = fabs(t[0]-t0);
    for (int j=1; j<size; j++) {
      double d = fabs(t[j]-t0);
      if (d < best) { best = d; jmax = j; }
    }
  }

  /* Per-mode buffers (mode-dependent: phase derivatives, amplitude and
     its derivative, basis-function*amplitude products). */
  double *omg   = (double*) calloc(size, sizeof(double));
  double *domg  = (double*) calloc(size, sizeof(double));
  double *m11   = (double*) calloc(size, sizeof(double));
  double *m12   = (double*) calloc(size, sizeof(double));
  double *m21   = (double*) calloc(size, sizeof(double));
  double *m22   = (double*) calloc(size, sizeof(double));
  double *p1tmp = (double*) calloc(size, sizeof(double));
  double *p2tmp = (double*) calloc(size, sizeof(double));

  const int all_modes[7] = {A22_K21, A22_K22, A22_K32, A22_K33, A22_K43, A22_K44, A22_K55};
  for (int im = 0; im < 7; im++) {
    const int k = all_modes[im];
    if (!hlm->kmask[k]) continue;

    double ai0=0., ai1=0., bi0=0., bi1=0.;

    /* (3,2) DOES get a genuine NQC correction here, exactly like (2,2):
       its own (nu,spin) fit (hatA0_l3m2/hatdA0_l3m2/hatd2A0_l3m2/
       omg0_l3m2/domg0_l3m2) is real, not placeholder, and (since
       comparable_mass_fitter.py's run_primary_fit_multistart dual-fit
       change) targets the SPHERICAL (3,2) mode -- exactly the basis this
       NQC correction (on
       hlm->ampli/phase[A22_K32], the SPHERICAL multipole) needs, matching
       kerrorbits/ringdown.py's own get_nqc_multipole (using
       nqc_point_vals_A22pk, also spherical-basis). This does NOT double
       up with A22_attach_mixed_mode_with_target's own use of
       QNMHybridFitCab_A22_lm(k32,...)'s A0/dA0/d2A0/omg0/domg0 output as
       Y(3,2)'s "natural" mixing target -- that target is read DIRECTLY
       from the fit, not by re-reading hlm->ampli[A22_K32][idx] (which
       this correction also updates, but only its resulting PHASE is
       still used downstream, via hlm->phase[A22_K32][idx] -- matching
       kerrorbits' own phases_Y[lm] = rm2pi_from_phase(phase of the
       NQC-corrected signal at t_match), which is exactly why this
       correction is needed here, not skippable in favor of the raw
       untouched array).
       (4,3) still gets NO NQC correction (ai0=ai1=bi0=bi1=0 identically,
       skipped below): no dedicated fit exists for it yet (A22_get_fit_set
       falls back to A22_fit_placeholder for k43) -- its spherical
       trajectory is simply left as whatever the standard, mode-
       independent factorized-resummed PN inspiral-plunge machinery
       already gives it, and A22_attach_mixed_mode_with_target's own
       cA2==0.0 fallback reads THIS untouched value/phase as Y(4,3)'s
       own "natural" target. */
    if (k != A22_K43) {
      D0(hlm_mrg->phase[k], dt, size, omg);
      D0(omg, dt, size, domg);

      double ca1,ca2,ca3,ca4,ca5,cb1,cb2,cb3,sigmar,sigmai;
      double max_A, max_dA, max_d2A, max_omg, max_domg;
      QNMHybridFitCab_A22_lm(k, nu, X1, X2, chi1, chi2, Mbh, abh,
                             &ca1,&ca2,&ca3,&ca4,&ca5,&cb1,&cb2,&cb3,&sigmar,&sigmai,
                             &max_A, &max_dA, &max_d2A, &max_omg, &max_domg,
                             NULL, NULL, NULL, NULL);

      double P[2], M[4], oodetM;

      for (int j=0; j<size; j++) p1tmp[j] = fabs(hlm_mrg->ampli[k][j]);
      for (int j=0; j<size; j++) { m11[j] = n1[j]*p1tmp[j]; m12[j] = n2[j]*p1tmp[j]; }
      D0(m11, dt, size, m21);
      D0(m12, dt, size, m22);
      D0(p1tmp, dt, size, p2tmp);

      P[0] = max_A  - p1tmp[jmax];
      P[1] = max_dA - p2tmp[jmax];
      M[0] = m11[jmax]; M[1] = m12[jmax]; M[2] = m21[jmax]; M[3] = m22[jmax];
      oodetM = 1.0/(M[0]*M[3]-M[1]*M[2]);
      if (isfinite(oodetM)) {
        ai0 = (M[3]*P[0] - M[1]*P[1])*oodetM;
        ai1 = (M[0]*P[1] - M[2]*P[0])*oodetM;
      }

      P[0] = omg[jmax]  - max_omg;
      P[1] = domg[jmax] - max_domg;
      M[0] = d_n4[jmax]; M[1] = d_n5[jmax]; M[2] = d2_n4[jmax]; M[3] = d2_n5[jmax];
      oodetM = 1.0/(M[0]*M[3]-M[1]*M[2]);
      if (isfinite(oodetM)) {
        bi0 = (M[3]*P[0] - M[1]*P[1])*oodetM;
        bi1 = (M[0]*P[1] - M[2]*P[0])*oodetM;
      }
    }

    if (VERBOSE) {
      printf("NQC coefficients for mode k=%d (A22 model):\n", k);
      PRFORMd("a1",ai0); PRFORMd("a2",ai1);
      PRFORMd("b1",bi0); PRFORMd("b2",bi1);
    }

    for (int j=0; j<size; j++) {
      hnqc->ampli[k][j] = 1. + ai0*n1[j] + ai1*n2[j];
      hnqc->phase[k][j] =      bi0*n4[j] + bi1*n5[j];
    }
    for (int j=0; j<size; j++) {
      hlm_mrg->ampli[k][j] *= hnqc->ampli[k][j];
      hlm_mrg->phase[k][j] -= hnqc->phase[k][j];
    }

    /* Apply the SAME (ai0,ai1,bi0,bi1) to the full waveform/dynamics --
       skipped entirely when all four are exactly 0 (k32/k43's own no-op,
       see above): NOT just an optimization -- n1f..n5f below are
       evaluated fresh, per-j, straight from dyn's own (FULL-length, not
       dyn_mrg's merger-segment-only) r/prstar/Omega/ddotr, unlike n1/n2/
       n4/n5 above (precomputed ONCE from dyn_mrg, already confirmed
       finite by every other mode's own genuinely-nonzero correction).
       For some configurations dyn's own arrays are 0/ill-defined at
       some j outside the merger segment (e.g. verified: high-spin,
       unequal-mass SXS:BBH:0293-like parameters), which 0.0*(0/0 or
       Inf) turns into a genuine NaN via IEEE-754 (0*NaN=NaN, not 0) --
       invisible for every OTHER mode, whose ai0..bi1 are never exactly
       0.0, but a real, silent corruption of hlm->ampli/phase[k] (from
       j=0, i.e. the WHOLE array, not just the ringdown tail) once k32/
       k43 started deliberately using an EXACT 0.0 no-op. */
    if (ai0 != 0.0 || ai1 != 0.0 || bi0 != 0.0 || bi1 != 0.0) {
    double *rf       = dyn->data[EOB_RAD];
    double *wf       = dyn->data[EOB_MOMG];
    double *pr_starf = dyn->data[EOB_PRSTAR];
    double *ddotrf   = dyn->data[EOB_DDOTR];
    const int fullsize = hlm->size;
    for (int j=0; j<fullsize; j++) {
      const double pr_star2f = SQ(pr_starf[j]);
      const double r2f       = SQ(rf[j]);
      const double w2f       = SQ(wf[j]);
      const double n1f = pr_star2f/(r2f*w2f);
      const double n2f = ddotrf[j]/(rf[j]*w2f);
      const double n4f = pr_starf[j]/(rf[j]*wf[j]);
      const double n5f = n4f*r2f*w2f;
      hlm->ampli[k][j] *= (1. + ai0*n1f + ai1*n2f);
      hlm->phase[k][j] -= (bi0*n4f + bi1*n5f);
    }
    }
  }

  free(n1); free(n2); free(n4); free(n5);
  free(d_n4); free(d_n5); free(d2_n4); free(d2_n5);
  free(omg); free(domg); free(m11); free(m12); free(m21); free(m22);
  free(p1tmp); free(p2tmp);
}
