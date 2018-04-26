/**
 * Copyright (C) 2017 Sebastiano Bernuzzi, Gregorio Carullo, Walter Del Pozzo, Alessandro Nagar, Ka Wa Tsang
 * This file is part of TEOBResumS
 *
 * TEOBResumS is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * TEOBResumS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 */

#include "TEOBResumS.h"

/** hlmNewt coefficients for amplitude */
static const double ChlmNewt_ampli[35] = {2.1137745587232057, 6.341323676169617, 0.1412325034218127, 1.7864655618418102, 4.9229202032627635, 0.023872650234580958, 0.2250735048768909, 1.7053495827316825, 4.763908164911493, 0.001122165903318321, 0.06333806741197714, 0.2945348827200268, 1.755276012972272, 5.0817902739730565, 0.00014954736544380072, 0.005296595280255638, 0.10342548105284892, 0.3713362832603404, 1.8983258440274462, 5.727111757630886, 5.184622059790144e-06, 0.0012191691413436815, 0.011783593824950922, 0.14639129995388936, 0.4653654097044924, 2.125638973894669, 6.685178460621457, 7.9955401278763745, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

/* 
   {8./3.*sqrt(Pi/5.), 8.*sqrt(Pi/5.),
   1./3.*sqrt(2.*Pi/35.), 8./3.*sqrt(Pi/7.), 3.*sqrt(6.*Pi/7.), 
   1./105.*sqrt(2.*Pi), 8./63.*sqrt(Pi), 9./5*sqrt(2*Pi/7.), 64./9.*sqrt(Pi/7.),
   1./180.*sqrt(Pi/77.), 16./135.*sqrt(Pi/11.), 9./20.*sqrt(3*Pi/22.), 256./45.*sqrt(Pi/33.), 125./12.*sqrt(5.*Pi/66.), 
   1./2079.*sqrt(2.*Pi/65.), 16./1485.*sqrt(Pi/13.), 81./385.*sqrt(Pi/13.), 1024./495.*sqrt(2.*Pi/195.), 625./63*sqrt(5.*Pi/429.), 432./5*sqrt(Pi/715.), 
   1./108108.*sqrt(Pi/10.), 8./3003.*sqrt(Pi/15.), 243./20020.*sqrt(3.*Pi/10.), 1024./1365.*sqrt(2.*Pi/165.), 3125./3276.*sqrt(5.*Pi/66.), 648./35.*sqrt(3.*Pi/715.), 16807./180.*sqrt(7.*Pi/4290.), 
   131072./315.*sqrt(2.*Pi/17017.), 0., 0., 0., 0., 0., 0., 0.
   };
*/

/** hlmNewt additive coefficients for phase */
static const double ChlmNewt_phase[35] = {4.71238898038469, 3.141592653589793, 4.71238898038469, 0.0, 1.5707963267948966, 1.5707963267948966, 0.0, 4.71238898038469, 3.141592653589793, 1.5707963267948966, 3.141592653589793, 4.71238898038469, 0.0, 1.5707963267948966, 4.71238898038469, 3.141592653589793, 1.5707963267948966, 0.0, 4.71238898038469, 3.141592653589793, 4.71238898038469, 0.0, 1.5707963267948966, 3.141592653589793, 4.71238898038469, 0.0, 1.5707963267948966, 3.141592653589793, 0, 0, 0, 0, 0, 0, 0};

/*{
  3.*Pi/2., Pi, 
  3.*Pi/2., 0., Pi/2., 
  Pi/2., 0., 3.*Pi/2., Pi, 
  Pi/2., Pi, 3.*Pi/2., 0., Pi/2., 
  3.*Pi/2.,  Pi, Pi/2., 0., 3.*Pi/2., Pi, 
  3.*Pi/2., 0., Pi/2., Pi, 3.*Pi/2., 0., Pi/2., 
  Pi, 0, 0, 0, 0, 0, 0, 0
  };
*/

/** Leading-order (Newtonian) prefactor  of the multipolar resummed waveform. 
    Reference: Damour, Iyer & Nagar, PRD 79, 064004 (2009) */
void eob_wav_hlmNewt(double r,
		     double Omega,
		     double phi,
		     double nu,
		     Waveform_lm_t *hlmNewt)
{
  /** Shorthands */
  double nu2   = nu*nu;
  double nu3   = nu*nu2;
  
  double vphi  = r*Omega;
  double vphi2 = vphi*vphi;
  double vphi3 = vphi*vphi2;
  double vphi4 = vphi*vphi3;
  double vphi5 = vphi*vphi4;
  double vphi6 = vphi*vphi5;
  double vphi7 = vphi*vphi6;
  double vphi8 = vphi*vphi7;
  
  /** Polynomials in nu */
  const double p1 = 1.;
  const double p2 = sqrt(1.-4.*nu); 
  const double p3 = (3.*nu-1.);
  const double p4 = (2.*nu-1.)*sqrt(1.-4.*nu);
  const double p5 = 1.-5.*nu+5.*nu2;
  const double p6 = (1.-4.*nu+3.*nu2)*sqrt(1.-4.*nu);
  const double p7 = 7.*nu3 - 14.*nu2 + 7.*nu -1.;
    
  const double phix2 = 2. * phi;
  const double phix3 = 3. * phi;
  const double phix4 = 4. * phi;
  const double phix5 = 5. * phi;
  const double phix6 = 6. * phi;
  const double phix7 = 7. * phi;
  
  const double pv23 = p2 * vphi3;
  const double pv34 = p3 * vphi4;
  const double pv45 = p4 * vphi5;
  const double pv56 = p5 * vphi6;
  const double pv67 = p6 * vphi7;
  const double pv78 = p7 * vphi8;
  
  double phim[35] = {
    phi,phix2,
    phi,phix2,phix3,
    phi,phix2,phix3,phix4,
    phi,phix2,phix3,phix4,phix5,
    phi,phix2,phix3,phix4,phix5,phix6,
    phi,phix2,phix3,phix4,phix5,phix6,phix7,
    phi,phix2,phix3,phix4,phix5,phix6,phix7,8.*phi
  };

  double Alm[35] = {
    pv23, p1 * vphi2, 
    pv23, pv34, pv23, 
    pv45, pv34, pv45, pv34, 
    pv45, pv56, pv45, pv56, pv45, 
    pv67, pv56, pv67, pv56, pv67, pv56, 
    pv67, pv78, pv67, pv78, pv67, pv78, pv67, 
    pv78, 1., 1., 1., 1., 1., 1., 1. 
  };
    
  /** Compute hlmNewt (without phase factor) in complex Polar coords */
  for (int k = 0; k < KMAX; k++) {
    hlmNewt->phase[k] = - phim[k] + ChlmNewt_phase[k];
    hlmNewt->ampli[k] = ChlmNewt_ampli[k] * Alm[k];
  }
  
}

/** Tail contribution to the resummed wave.   
    Ref. Damour, Iyer & Nagar, PRD 79, 064004 (2009) */
void eob_wav_hhatlmTail(double Omega, double Hreal, double bphys, Waveform_lm_t *tlm)
{
  double k;
  double hhatk;
  
  gsl_sf_result num_rad;
  gsl_sf_result num_phase;
  gsl_sf_result denom_rad;
  gsl_sf_result denom_phase;
  
  double ratio_rad;
  double ratio_ang;
  double tlm_rad;
  double tlm_phase;
  
  int i;
  for (i = 0; i < KMAX; i++) {
    k     = MINDEX[i] * Omega;
    hhatk = k * Hreal;
    
    gsl_sf_lngamma_complex_e(LINDEX[i] + 1., -2.*hhatk, &num_rad, &num_phase);
    gsl_sf_lngamma_complex_e(LINDEX[i] + 1., 0., &denom_rad, &denom_phase);
    
    ratio_rad     = num_rad.val-denom_rad.val;
    ratio_ang     = num_phase.val-0.;
    
    tlm_rad       = ratio_rad + Pi * hhatk;
    tlm_phase     = ratio_ang + 2.*hhatk*log(2.*k*bphys);
    
    tlm->ampli[i] = exp(tlm_rad);
    tlm->phase[i] = tlm_phase;
  }

}

/** Alternative implementation of the phase of the tail factor */
void eob_wav_speedyTail(double Omega, double Hreal, double bphys, Waveform_lm_t *tlm)
{
  double x;
  double x2;
  double x3;
  double x4;
  double x5;
  double tlm_ang;
  double num_ang;
  
  /** Fit coefficients*/
  const double b1[] = {
    0.1113090643348557, 0.1112593821157397, 0.0424759238428813, 0.0424489015884926, 0.0424717446800903, 0.0215953972500844, 0.0215873812155663, 0.0215776183122621, 0.0216017621863542, 0.0128123696874894, 0.0128097056242375, 0.0128038943888768, 0.0128025242617949, 0.0128202485907368, 0.0083762045692408, 0.0083751913886140, 0.0083724067460769, 0.0083694435961860, 0.0083710364141552, 0.0083834483913443, 0.0058540393221396, 0.0058536069384738, 0.0058522594457692, 0.0058502436535615, 0.0058491157293566, 0.0058514875071582, 0.0058602498033381, 0.0042956812356573, 0.0042954784390887, 0.0042947951664056, 0.0042935886137697, 0.0042923691461384, 0.0042922256848799, 0.0042945927126022, 0.0043009106861259};
  
  const double b2[] = {
    0.0004643273300862, 0.0009375605440004, 0.0000597134489198, 0.0002551406918111, 0.0001741036904709, 0.0000124649041611, 0.0000685496215625, 0.0001131160409390, 0.0000419907542591, 0.0000035218982282, 0.0000219211271097, 0.0000473186962874, 0.0000524142634057, 0.0000106823372552, 0.0000012237574387, 0.0000081742188269, 0.0000201940563214, 0.0000295722761753, 0.0000260539631956, 0.0000018994753518, 0.0000004932942990, 0.0000034477210351, 0.0000092294406360, 0.0000155143073237, 0.0000183386499818, 0.0000137922469695, -0.0000007075155453, 0.0000002223410995, 0.0000016045317657, 0.0000045260028113, 0.0000082655700107, 0.0000112393599417, 0.0000115758243113, 0.0000076838709956, -0.0000014020591745};
  
  const double b3[] = {
    -0.0221835462237291, -0.0235386333304348, -0.0042911639711832, -0.0047431560217121, -0.0046577314472149, -0.0013089557502947, -0.0014343968205390, -0.0014978542575474, -0.0014329302934532, -0.0005167994164556, -0.0005573939123058, -0.0005921030407223, -0.0005978284714483, -0.0005673965369076, -0.0002409269302708, -0.0002561516055118, -0.0002723768586352, -0.0002815958312453, -0.0002792078156272, -0.0002646630240693, -0.0001261183503407, -0.0001325622938779, -0.0001403198638518, -0.0001464084186977, -0.0001485971591029, -0.0001459023931717, -0.0001384829633836, -0.0000719062974278, -0.0000749128468013, -0.0000788187384314, -0.0000824202283094, -0.0000846673495936, -0.0000849054394951, -0.0000829269749240, -0.0000788883333858};
  
  const double b4[] = {
    0.0058366730167965, 0.0070452306758401, 0.0006914295465364, 0.0010322294603561, 0.0010057563135650, 0.0001394203795507, 0.0002309706405978, 0.0002596611624417, 0.0002409588083156, 0.0000386949167221, 0.0000679154947896, 0.0000830199015202, 0.0000850120755064, 0.0000780125513602, 0.0000133034384660, 0.0000241813441339, 0.0000311573885555, 0.0000340233089866, 0.0000335167900637, 0.0000307571022927, 0.0000053305073331, 0.0000099143129290, 0.0000132296989826, 0.0000150959309402, 0.0000156304390748, 0.0000151274875147, 0.0000139320508803, 0.0000023959090314, 0.0000045285807761, 0.0000061918979830, 0.0000072894226381, 0.0000078251853305, 0.0000078772667984, 0.0000075606242809, 0.0000069956215270
  };
  
  double Tlm_real[KMAX];
  eob_flx_Tlm(Omega*Hreal, Tlm_real);
  
  /** Pre-computed psi */
  const double psi[] = {0.9227843350984671394, 0.9227843350984671394,
			1.256117668431800473, 1.256117668431800473, 1.256117668431800473,
			1.506117668431800473, 1.506117668431800473, 1.506117668431800473, 1.506117668431800473,
			1.706117668431800473, 1.706117668431800473, 1.706117668431800473, 1.706117668431800473, 1.706117668431800473,
			1.872784335098467139, 1.872784335098467139, 1.872784335098467139, 1.872784335098467139, 1.872784335098467139, 1.872784335098467139,
			2.015641477955609997, 2.015641477955609997, 2.015641477955609997, 2.015641477955609997, 2.015641477955609997, 2.015641477955609997, 2.015641477955609997,
			2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997, 2.140641477955609997};
  
  double k;
  int i;
  for (i=0; i<KMAX; i++) {
    k  = MINDEX[i] * Omega;
    x  = k * Hreal; /* hathatk */
    x2 = x * x;
    x3 = x2 * x;
    x4 = x3 * x;
    x5 = x4 * x;      
    num_ang   = 1. + b1[i]*x2 + b2[i]*x3 + b3[i]*x4 + b4[i]*x5; 
    tlm_ang   = (- 2. * psi[i] * x * num_ang) + 2.*x* log(2. * k * bphys);
    tlm->ampli[i] = Tlm_real[i];
    tlm->phase[i] = tlm_ang;
  }
  
}

/** Residual phase corrections delta_{lm} up to l=m=5.
    Reference(s)
    Damour, Iyer & Nagar, PRD 79, 064004 (2008)
    Fujita & Iyer, PRD 82 044051 (2010)
    Faye et al., Class. Q. Grav. 29 175004 (2012)
    Damour, Nagar & Bernuzzi, PRD 87, 084035 (2013) */   
//TODO: this routine can be optimized: precompute coefficients c(nu)
void eob_wav_deltalm(double Hreal,double Omega,double nu, double *dlm)
{
    
  /** Useful shorthands*/
  double pi     = Pi;
  double pi2    = SQ(Pi);
  double nu2    = SQ(nu);
  double y      = cbrt(Hreal*Omega*Hreal*Omega);
  double sqrt_y = sqrt(y);
  double y3     = y*y*y;
  double y32    = Hreal*Omega;
  
  /** Leading order contributions*/
  double delta22LO = 7./3.   * y32;
  double delta21LO = 2./3.   * y32;
  double delta33LO = 13./10. * y32;
  double delta31LO = 13./30. * y32;
  
  /** Init phase */
  for (int k = 0; k < KMAX; k++) {
    dlm[k] = 0.;
  }
  
  /** Residual phases in Pade-resummed form when possible */
  double num;
  double den;

  /* l=2 */
  /* Pade(1,2) approximant */
  num        = 69020.*nu + 5992.*pi*sqrt_y;
  den        = 5992.*pi*sqrt_y + 2456.*nu*(28.+493.*nu* y);
  dlm[0] = delta21LO*num/den;
  /* Pade(2,2) approximant */
  num        = (808920.*nu*pi*sqrt(y) + 137388.*pi2*y + 35.*nu2*(136080. + (154975. - 1359276.*nu)*y));
  den        = (808920.*nu*pi*sqrt(y) + 137388.*pi2*y + 35.*nu2*(136080. + (154975. + 40404.*nu)*y));
  dlm[1] = delta22LO*num/den;
  
  /* l=3 */
  /* Pade(1,2) approximant */
  num        = 4641.*nu + 1690.*pi*sqrt_y;
  den        = num + 18207.*nu2*y;
  dlm[2] = delta31LO*num/den;
  /* Taylor-expanded form */
  num        = 1.  + 94770.*pi/(566279.*nu)*sqrt_y;
  den        = num + 80897.* nu/3159.*y;
  dlm[3] = (10.+33.*nu)/(15.*(1.-3.*nu)) * y32 + 52./21.*pi*y3;
  /* Pade(1,2) approximant */
  dlm[4] = delta33LO*num/den;
  
  /* l=4 */
  dlm[5] =   (2.+507.*nu)/(10.*(1.-2.*nu))*y32   + 1571./3465.*pi*y3;
  dlm[6] =  7.*(1.+6.*nu)/(15.*(1.-3.*nu))*y32   + 6284./3465.*pi*y3;
  dlm[7] = (486.+4961.*nu)/(810.*(1.-2.*nu))*y32 + 1571./385.*pi*y3;
  dlm[8] =  (112.+219.*nu)/(120.*(1.-3.*nu))*y32 + 25136./3465.*pi*y3;
  
  /* l=5 */
  dlm[9] = (96875. + 857528.*nu)/(131250.*(1.-2.*nu))*y32;
  
}

/** Resummed amplitudes in the general nu-dependent case.
 *  Refs:
 *  . Damour, Iyer & Nagar, PRD 79, 064004 (2009)     [theory]
 *  . Fujita & Iyer, PRD 82, 044051 (2010)            [test-mass 5.5PN]
 *  . Damour, Nagar & Bernuzzi, PRD 87, 084035 (2013) [complete information]
 */
void eob_wav_flm(double x,double nu, double *rholm, double *flm)
{

  /** Coefficients */
  static double clm[KMAX][6];
  
  const double nu2 = nu*nu;
  const double nu3 = nu*nu2;
  const double nu4 = nu*nu3;
  
  static int firstcall = 1;
  if (firstcall) {
    if (0) printf("Precompute some rholm coefs\n");
    firstcall = 0;
    
    for (int k=0; k<KMAX; k++) clm[k][0] = 1.;
    for (int k=0; k<KMAX; k++) for (int n=1; n<6; n++) clm[k][n] = 0.;

    /** (2,1) */
    clm[0][1] = (-1.0535714285714286 + 0.27380952380952384 *nu);
    clm[0][2] = (-0.8327841553287982 - 0.7789824263038548  *nu + 0.13116496598639457*nu2);
    /* clm[0][3] = (2.9192806270460925  - 1.019047619047619   *el1); */
    /* clm[0][4] = (-1.28235780892213   + 1.073639455782313   *el1); */
    /* clm[0][5] = (-3.8466571723355227 + 0.8486467106683944  *el1)*PMTERMS_eps; */

    /** (2,2) */
    clm[1][1] = (-1.0238095238095237 + 0.6547619047619048*nu);
    clm[1][2] = (-1.94208238851096   - 1.5601379440665155*nu + 0.4625614134542706*nu2);
    /* clm[1][3] = (12.736034731834051  - 2.902228713904598 *nu - 1.9301558466099282*nu2 + 0.2715020968103451*nu3 - 4.076190476190476*el2); */
    /* clm[1][4] = (-2.4172313935587004 + 4.173242630385488 *el2); */
    /* clm[1][5] = (-30.14143102836864  + 7.916297736025627 *el2); */

    /** (3,1) */
    clm[2][1] = (-0.7222222222222222 - 0.2222222222222222*nu);
    clm[2][2] = (0.014169472502805836 - 0.9455667789001122*nu - 0.46520763187429853*nu2);
    /* clm[2][3] = (1.9098284139598072 - 0.4126984126984127*el1); */
    /* clm[2][4] = (0.5368150316615179 + 0.2980599647266314*el1); */
    /* clm[2][5] = (1.4497991763035063 - 0.0058477188106817735*el1)*PMTERMS_eps; */

    /** (3,2) */
    clm[3][1] = (0.003703703703703704*(328. - 1115.*nu + 320.*nu2))/(-1. + 3.*nu);
    clm[3][2] = (6.235191420376606e-7*(-1.444528e6 + 8.050045e6*nu - 4.725605e6*nu2 - 2.033896e7*nu3 + 3.08564e6*nu4))/((-1. + 3.*nu)*(-1. + 3.*nu));
    /* clm[3][3] = (6.220997955214429 - 1.6507936507936507*el2); */
    /* clm[3][4] = (-3.4527288879001268 + 2.005408583186361*el2)*PMTERMS_eps; */

    /** (3,3) */
    clm[4][1] =  (-1.1666666666666667 + 0.6666666666666666*nu);
    clm[4][2] = (-1.6967171717171716 - 1.8797979797979798*nu + 0.45151515151515154*nu2);
    /* clm[4][3] = (14.10891386831863 - 3.7142857142857144*el3); */
    /* clm[4][4] = (-6.723375314944128 + 4.333333333333333*el3); */
    /* clm[4][5] = (-29.568699895427518 + 6.302092352092352*el3)*PMTERMS_eps; */

    /** (4,1) */
    clm[5][1] = (0.001893939393939394*(602. - 1385.*nu + 288.*nu2))/(-1. + 2.*nu);
    clm[5][2] = (- 0.36778992787515513);
    /* clm[5][3] = (0.6981550175535535 - 0.2266955266955267*el1); */
    /* clm[5][4] = (-0.7931524512893319 + 0.2584672482399755*el1)*PMTERMS_eps; */

    /** (4,2) */
    clm[6][1] = (0.0007575757575757576*(1146. - 3530.*nu + 285.*nu2))/(-1. + 3.*nu);
    clm[6][2] = - (3.1534122443213353e-9*(1.14859044e8 - 2.95834536e8*nu - 1.204388696e9*nu2 + 3.04798116e9*nu3 + 3.79526805e8*nu4))/((-1. + 3.*nu)*(-1. + 3.*nu));
    /* clm[6][3] = 4.550378418934105e-12*(8.48238724511e11 - 1.9927619712e11*el2); */
    /* clm[6][4] = (-0.6621921297263365 + 0.787251738160829*el2)*PMTERMS_eps; */
    
    /** (4,3) */
    clm[7][1] = (0.005681818181818182*(222. - 547.*nu + 160.*nu2))/(-1. + 2.*nu);
    clm[7][2] = (- 0.9783218202252293);
    /* clm[7][3] = (8.519456157072423 - 2.0402597402597404*el3)*PMTERMS_eps; */
    /* clm[7][4] = (-5.353216984886716 + 2.5735094451003544*el3)*PMTERMS_eps; */

    /** (4,4) */
    clm[8][1] = (0.0007575757575757576*(1614. - 5870.*nu + 2625.*nu2))/(-1. + 3.*nu);
    clm[8][2] = (3.1534122443213353e-9*(-5.11573572e8 + 2.338945704e9*nu - 3.13857376e8*nu2 - 6.733146e9*nu3 + 1.252563795e9*nu4))/((-1. + 3.*nu)*(-1. + 3.*nu));
    /* clm[8][3] = (15.108111214795123 - 3.627128427128427*el4); */
    /* clm[8][4] = (-8.857121657199649 + 4.434988849534304*el4)*PMTERMS_eps; */

    /** (5,1) */
    clm[9][1] = (0.002564102564102564*(319. - 626.*nu + 8.*nu2))/(-1. + 2.*nu);
    clm[9][2] = (- 0.1047896120973044);
    /* clm[9][3] = (0.642701885362399 - 0.14414918414918415*el1)*PMTERMS_eps; */
    /* clm[9][4] = (-0.07651588046467575 + 0.11790664036817883*el1)*PMTERMS_eps; */

    /** (5,2) */
    clm[10][1] = (0.00007326007326007326*(-15828. + 84679.*nu - 104930.*nu2 + 21980.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[10][2] = (- 0.4629337197600934)*PMTERMS_eps; 
    /* clm[10][3] = (2.354458371550237 - 0.5765967365967366*el2)*PMTERMS_eps; */

    /** (5,3) */
    clm[11][1] = (0.002564102564102564*(375. - 850.*nu + 176.*nu2))/(-1. + 2.*nu);
    clm[11][2] = (- 0.5788010707241477);
    /* clm[11][3] = (5.733973288504755 - 1.2973426573426574*el3)*PMTERMS_eps; */
    /* clm[11][4] = (-1.9573287625526001 + 1.2474448628294783*el3)*PMTERMS_eps; */

    /** (5,4) */
    clm[12][1] = (0.00007326007326007326*(-17448. + 96019.*nu - 127610.*nu2 + 33320.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[12][2] = (- 1.0442142414362194)*PMTERMS_eps;
    /* clm[12][3] = (10.252052781721588 - 2.3063869463869464*el4)*PMTERMS_eps; */

    /** (5,5) */
    clm[13][1] = (0.002564102564102564*(487. - 1298.*nu + 512.*nu2))/(-1. + 2.*nu);
    clm[13][2] = (- 1.5749727622804546);
    /* clm[13][3] = (15.939827047208668 - 3.6037296037296036*el5)*PMTERMS_eps; */
    /* clm[13][4] = (-10.272578060123237 + 4.500041838503377*el5)*PMTERMS_eps; */
    
    /** (6,1) */
    clm[14][1] = (0.006944444444444444*(-161. + 694.*nu - 670.*nu2 + 124.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[14][2] = (- 0.29175486850885135)*PMTERMS_eps;
    /* clm[14][3] = (0.21653486654395454 - 0.10001110001110002*el1)*PMTERMS_eps; */

    /** (6,2) */
    clm[15][1] = (0.011904761904761904*(-74. + 378.*nu - 413.*nu2 + 49.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[15][2] = ( - 0.24797525070634313)*PMTERMS_eps;
    /* clm[15][3] = (1.7942694138754138 - 0.40004440004440006*el2)*PMTERMS_eps; */

    /** (6,3) */
    clm[16][1] = (0.006944444444444444*(-169. + 742.*nu - 750.*nu2 + 156.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[16][2] = (- 0.5605554442947213)*PMTERMS_eps;
    /* clm[16][3] = (4.002558222882566 - 0.9000999000999002*el3)*PMTERMS_eps; */

    /** (6,4) */
    clm[17][1] = (0.011904761904761904*(-86. + 462.*nu - 581.*nu2 + 133.*nu3))/(1. - 5.*nu + 5.*nu2);
    clm[17][2] = (- 0.7228451986855349)*PMTERMS_eps;
    /* clm[17][3] = (7.359388663371044 - 1.6001776001776002*el4)*PMTERMS_eps; */

    /** (6,5) */
    clm[18][1] = (0.006944444444444444*(-185. + 838.*nu - 910.*nu2 + 220.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[18][2] = (- 1.0973940686333457)*PMTERMS_eps;
    /* clm[18][3] = (11.623366217471297 - 2.5002775002775004*el5)*PMTERMS_eps; */

    /** (6,6) */
    clm[19][1] = (0.011904761904761904*(-106. + 602.*nu - 861.*nu2 + 273.*nu3))/(1. - 5.*nu + 5.*nu2); 
    clm[19][2] = (- 1.5543111183867486)*PMTERMS_eps;
    /* clm[19][3] = (16.645950799433503 - 3.6003996003996006*el6)*PMTERMS_eps; */

    /** (7,1) */
    clm[20][1] = (0.0014005602240896359*(-618. + 2518.*nu - 2083.*nu2 + 228.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[20][2] = ( - 0.1508235111143767)*PMTERMS_eps;
    /* clm[20][3] = (0.2581280702019663 - 0.07355557607658449*el1)*PMTERMS_eps; */

    /** (7,2) */
    clm[21][1] = (0.00006669334400426837*(16832. - 123489.*nu + 273924.*nu2 - 190239.*nu3 + 32760.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[21][2] = (- 0.351319484450667)*PMTERMS_eps;
    
    /** (7,3) */
    clm[22][1] = (0.0014005602240896359*(-666. + 2806.*nu - 2563.*nu2 + 420.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[22][2] = (- 0.37187416047628863)*PMTERMS_eps;
    /* clm[22][3] = (3.0835293524055283 - 0.6620001846892604*el3)*PMTERMS_eps; */

    /** (7,4) */
    clm[23][1] = (0.00006669334400426837*(17756. - 131805.*nu + 298872.*nu2 - 217959.*nu3 + 41076.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[23][2] = (- 0.6473746896670599)*PMTERMS_eps;
    
    /** (7,5) */
    clm[24][1] = (0.0014005602240896359*(-762. + 3382.*nu - 3523.*nu2 + 804.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[24][2] = (- 0.8269193364414116)*PMTERMS_eps;
    /* clm[24][3] = (8.750589067052443 - 1.838889401914612*el5)*PMTERMS_eps; */

    /** (7,6) */
    clm[25][1] = (0.0006002400960384153*(2144. - 16185.*nu + 37828.*nu2 - 29351.*nu3 + 6104.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[25][2] = (- 1.1403265020692532)*PMTERMS_eps;
    
    /** (7,7) */
    clm[26][1] = (0.0014005602240896359*(-906. + 4246.*nu - 4963.*nu2 + 1380.*nu3))/(1. - 4.*nu + 3.*nu2);
    clm[26][2] = (- 1.5418467934923434)*PMTERMS_eps;
    /* clm[26][3] = (17.255875091408523 - 3.6042232277526396*el7)*PMTERMS_eps; */

    /** (8,1) */
    clm[27][1] = (0.00005482456140350877*(20022. - 126451.*nu + 236922.*nu2 - 138430.*nu3 + 21640.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
    clm[27][2] = (- 0.26842133517043704)*PMTERMS_eps;

    /** (8,2) */
    clm[28][1] = (0.0003654970760233918*(2462. - 17598.*nu + 37119.*nu2 - 22845.*nu3 + 3063.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[28][2] = (- 0.2261796441029474)*PMTERMS_eps;

    /** (8,3) */
    clm[29][1] = (0.00005482456140350877*(20598. - 131059.*nu + 249018.*nu2 - 149950.*nu3 + 24520.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
    clm[29][2] = (- 0.4196774909106648)*PMTERMS_eps;

    /** (8,4) */
    clm[30][1] = (0.0003654970760233918*(2666. - 19434.*nu + 42627.*nu2 - 28965.*nu3 + 4899.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[30][2] = (- 0.47652059150068155)*PMTERMS_eps;

    /** (8,5) */
    clm[31][1] = (0.00027412280701754384*(4350. - 28055.*nu + 54642.*nu2 - 34598.*nu3 + 6056.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
    clm[31][2] = (- 0.7220789990670207)*PMTERMS_eps;

    /** (8,6) */
    clm[32][1] = (0.0010964912280701754*(1002. - 7498.*nu + 17269.*nu2 - 13055.*nu3 + 2653.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[32][2] = (- 0.9061610303170207)*PMTERMS_eps;

    /** (8,7) */
    clm[33][1] = (0.00005482456140350877*(23478. - 154099.*nu + 309498.*nu2 - 207550.*nu3 + 38920.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
    clm[33][2] = (- 1.175404252991305)*PMTERMS_eps;

    /** (8,8) */
    clm[34][1] = (0.0003654970760233918*(3482. - 26778.*nu + 64659.*nu2 - 53445.*nu3 + 12243.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
    clm[34][2] = (- 1.5337092502821381)*PMTERMS_eps;

  }

  /** Compute EulerLogs */
  const double el1 = Eulerlog(x,1);
  const double el2 = Eulerlog(x,2);
  const double el3 = Eulerlog(x,3);
  const double el4 = Eulerlog(x,4);
  const double el5 = Eulerlog(x,5);
  const double el6 = Eulerlog(x,6);
  const double el7 = Eulerlog(x,7);

  /** Coefs with Eulerlogs */
  clm[0][3] = (2.9192806270460925  - 1.019047619047619   *el1);
  clm[0][4] = (-1.28235780892213   + 1.073639455782313   *el1);
  clm[0][5] = (-3.8466571723355227 + 0.8486467106683944  *el1)*PMTERMS_eps;
  
  clm[1][3] = (12.736034731834051  - 2.902228713904598 *nu - 1.9301558466099282*nu2 + 0.2715020968103451*nu3 - 4.076190476190476*el2);
  clm[1][4] = (-2.4172313935587004 + 4.173242630385488 *el2);
  clm[1][5] = (-30.14143102836864  + 7.916297736025627 *el2);
  
  clm[2][3] = (1.9098284139598072 - 0.4126984126984127*el1);
  clm[2][4] = (0.5368150316615179 + 0.2980599647266314*el1);
  clm[2][5] = (1.4497991763035063 - 0.0058477188106817735*el1)*PMTERMS_eps;
    
  clm[3][3] = (6.220997955214429 - 1.6507936507936507*el2);
  clm[3][4] = (-3.4527288879001268 + 2.005408583186361*el2)*PMTERMS_eps;
  
  clm[4][3] = (14.10891386831863 - 3.7142857142857144*el3);
  clm[4][4] = (-6.723375314944128 + 4.333333333333333*el3);
  clm[4][5] = (-29.568699895427518 + 6.302092352092352*el3)*PMTERMS_eps;
  
  clm[5][3] = (0.6981550175535535 - 0.2266955266955267*el1);
  clm[5][4] = (-0.7931524512893319 + 0.2584672482399755*el1)*PMTERMS_eps;
  
  clm[6][3] = 4.550378418934105e-12*(8.48238724511e11 - 1.9927619712e11*el2);
  clm[6][4] = (-0.6621921297263365 + 0.787251738160829*el2)*PMTERMS_eps;
  
  clm[7][3] = (8.519456157072423 - 2.0402597402597404*el3)*PMTERMS_eps;
  clm[7][4] = (-5.353216984886716 + 2.5735094451003544*el3)*PMTERMS_eps;
  
  clm[8][3] = (15.108111214795123 - 3.627128427128427*el4);
  clm[8][4] = (-8.857121657199649 + 4.434988849534304*el4)*PMTERMS_eps;
  
  clm[9][3] = (0.642701885362399 - 0.14414918414918415*el1)*PMTERMS_eps;
  clm[9][4] = (-0.07651588046467575 + 0.11790664036817883*el1)*PMTERMS_eps;
  
  clm[10][3] = (2.354458371550237 - 0.5765967365967366*el2)*PMTERMS_eps;
  
  clm[11][3] = (5.733973288504755 - 1.2973426573426574*el3)*PMTERMS_eps;
  clm[11][4] = (-1.9573287625526001 + 1.2474448628294783*el3)*PMTERMS_eps;
  
  clm[12][3] = (10.252052781721588 - 2.3063869463869464*el4)*PMTERMS_eps;
  
  clm[13][3] = (15.939827047208668 - 3.6037296037296036*el5)*PMTERMS_eps;
  clm[13][4] = (-10.272578060123237 + 4.500041838503377*el5)*PMTERMS_eps;
  
  clm[14][3] = (0.21653486654395454 - 0.10001110001110002*el1)*PMTERMS_eps;
  
  clm[15][3] = (1.7942694138754138 - 0.40004440004440006*el2)*PMTERMS_eps;
  
  clm[16][3] = (4.002558222882566 - 0.9000999000999002*el3)*PMTERMS_eps;
  
  clm[17][3] = (7.359388663371044 - 1.6001776001776002*el4)*PMTERMS_eps;
  
  clm[18][3] = (11.623366217471297 - 2.5002775002775004*el5)*PMTERMS_eps;
  
  clm[19][3] = (16.645950799433503 - 3.6003996003996006*el6)*PMTERMS_eps;
  
  clm[20][3] = (0.2581280702019663 - 0.07355557607658449*el1)*PMTERMS_eps;
  
  clm[22][3] = (3.0835293524055283 - 0.6620001846892604*el3)*PMTERMS_eps;
  
  clm[24][3] = (8.750589067052443 - 1.838889401914612*el5)*PMTERMS_eps;
  
  clm[26][3] = (17.255875091408523 - 3.6042232277526396*el7)*PMTERMS_eps;
    
  /** rho_lm */
  const double x2  = x*x;
  const double x3  = x*x2;
  const double x4  = x*x3;
  const double x5  = x*x4;
  const double xn[] = {1.,x,x2,x3,x4,x5};

  for (int k=0; k<KMAX; k++) {
    /* Note: the two sums give different result */
#if (1)
    rholm[k] = clm[k][0];
    for (int n=1; n<6; n++) {
      rholm[k] += clm[k][n]*xn[n];
    }
#else  
    rholm[k] = x5*clm[k][5];
    for (int n=5; n-- >1; ) { // 4,3,2,1 // 
      rholm[k] += clm[k][n]*xn[n];
    } 
    rholm[k] += clm[k][0]; 
#endif
  }

  /** Amplitudes */
  for (int k = 0; k < KMAX; k++) {
      flm[k] = gsl_pow_int(rholm[k], LINDEX[k]);
  }

}

/** Resummed amplitudes for the spin case. 
    This function computes the residual amplitude corrections flm's as 
    introduced in Damour, Iyer & Nagar, PRD 79, 064004 (2008).
    The orbital part is taken at the usual 3^{+2} PN order, i.e. 3PN terms
    are integrated by the 4PN and 5PN test-particle terms, with the higher
    modes obtained by Fujita & Iyer.
    Note that the variables called here (a1,a2)
    are what we usually cal tilde{a}_1 and tilde{a}_2 and are defined as
    a1 = X1*chi1, a2=X2*chi2 and are passed here as parameters. Special 
    combinations of these quantities are used here to write the spin-dependent
    part of the waveform in particularly compact form, so that the (spinning)
    test-particle limit is recovered just by visual inspection of the equations */
void eob_wav_flm_s(double x, double nu, double X1, double X2, double chi1, double chi2, double a1, double a2, double C_Q1, double C_Q2, int usetidal,
		   double *rholm, double *flm)
{

  /** Orbital part */
  //double rholm_orb[KMAX], flm_orb[KMAX];
  eob_wav_flm(x,nu, rholm, flm);

  /** Spin corrections */
  double rho22S;
  double rho32S;
  double rho44S;
  double rho42S;
  double f21S;
  double f33S;
  double f31S;
  double f43S;
  double f41S;
      
  const double a0      = a1+a2;
  const double a12     = a1-a2;
  const double X12     = X1-X2;
  const double a0X12   = a0*X12;
  const double a12X12  = a12*X12;
  
  const double v  = sqrt(x);
  const double v2 = x;
  const double v3 = v*v2;
  const double v4 = v3*v;
  const double v5 = v4*v;
     
  /** l=m=2 multipole */
  /* spin-orbit */
  const double cSO_lo    = (-0.5*a0 - a12X12/6.);
  const double cSO_nlo   = (-52./63.-19./504.*nu)*a0 - (50./63.+209./504.*nu)*a12X12;
  
  /* SPIN-SPIN contribution */
  double cSS_lo;
  if (usetidal) {
#if (EXCLUDESPINSPINTIDES)
    /* Switch off spin-spin-tidal couplings */
    /* See also: eob_dyn_s_get_rc() */
    cSS_lo = 0.;
    /* Above code switch off everything, 
       Alt. one can set C_Q1=C_Q2=0, but keep the term: */
    /*
      cSS_lo = a1*a2;
    */
#else
    cSS_lo = 0.5*(C_Q1*a1*a1 + 2.*a1*a2 + C_Q2*a2*a2);
#endif
  } else {
    cSS_lo = 0.5*a0*a0; 
  }

  /* rho_22^S: Eq. (80) of Damour & Nagar, PRD 90, 044018 (2014) */
  rho22S = cSO_lo*v3 + cSS_lo*v4 + cSO_nlo*v5 ;
    
  /** l>=3, m=even: multipoles rewritten in compact and self-explanatory form */
  rho32S = (a0-a12X12)/(3.*(1.-3.*nu))*v;
  rho44S = (-19./30.*a0 -  (1.-21.*nu)/(30.-90.*nu)*a12X12)*v3;
  rho42S = ( -1./30.*a0 - (19.-39.*nu)/(30.-90.*nu)*a12X12)*v3;
  
  /** l>=2, m=odd: multipoles rewritten in compact and self-explanatory form */
  f21S = -1.5*a12*v + ((110./21. + 79./84.*nu)*a12 - 13./84.*a0X12)*v3;
  f33S = ((-0.25 + 2.5*nu)*a12 - 1.75*a0X12)*v3;
  f31S = ((-2.25 + 6.5*nu)*a12 + 0.25*a0X12)*v3;
  f43S = (( 5. -10.*nu)*a12 - 5.*a0X12)/(-4.+8.*nu)*v;
  f41S = f43S;
    
  /** Amplitudes (correct with spin terms) */
  flm[0] = gsl_pow_int(rholm[0], 2);
  flm[0] = (X12*flm[0] + f21S);
  
  flm[1] = gsl_pow_int(rholm[1]+ rho22S, 2);
  
  flm[2] = gsl_pow_int(rholm[2], 3);
  flm[2] = (X12*flm[2] + f31S);
  
  flm[3] = gsl_pow_int(rholm[3]+ rho32S, 3);
  
  flm[4] = gsl_pow_int(rholm[4], 3);
  flm[4] = (X12*flm[4] + f33S);
  
  flm[5] = gsl_pow_int(rholm[5], 4);
  flm[5] = (X12*flm[5] + f41S);
  
  flm[6] = gsl_pow_int(rholm[6] + rho42S, 4);
  
  flm[7] = gsl_pow_int(rholm[7], 4);
  flm[7] = (X12*flm[7] + f43S);
  
  flm[8] = gsl_pow_int(rholm[8] + rho44S, 4);
      
}

/** Calculate tidal correction to multipolar waveform amplitude
    Ref. Damour, Nagar & Villain, Phys.Rev. D85 (2012) 123007 */
void eob_wav_hlmTidal(double x, Dynamics *dyn, double *hTidallm)
{
  const double XA       = dyn->X1;
  const double XB       = dyn->X2;
  const double khatA_2  = dyn->khatA2;
  const double khatB_2  = dyn->khatB2;
  
  const double x5 = gsl_pow_int(x,5);
  
  double hA[KMAX], hB[KMAX], betaA1[KMAX],betaB1[KMAX];

  memset(hTidallm, 0., KMAX*sizeof(double));
  memset(hA, 0., KMAX*sizeof(double));
  memset(hB, 0., KMAX*sizeof(double));
  memset(betaA1, 0., KMAX*sizeof(double));
  memset(betaB1, 0., KMAX*sizeof(double));

  /** l=2 */
  hA[1]     = 2 * khatA_2 *(XA/XB+3);
  hB[1]     = 2 * khatB_2 *(XB/XA+3);

  betaA1[1] = (-202. + 560*XA - 340*XA*XA + 45*XA*XA*XA)/(42*(3-2*XA));
  betaB1[1] = (-202. + 560*XB - 340*XB*XB + 45*XB*XB*XB)/(42*(3-2*XB));
    
  hA[0]     = 3 * khatA_2 * XB * (3-4*XA)/XA;
  hB[0]     = 3 * khatB_2 * XA * (3-4*XB)/XB;
    
  /** l=3 */
  hA[2] = hA[4];
  hB[2] = hB[4];
  
  hA[4] = 12 * khatA_2 * XB*XB/XA;
  hB[4] = 12 * khatB_2 * XA*XA/XB;
  
  /** l=2 */
  /* (2,1) */
  hTidallm[0] = ( -hA[0] + hB[0] )*x5;
  /* (2,2) */
  hTidallm[1] = ( hA[1]*(1. + betaA1[1]*x) + hB[1]*(1. + betaB1[1]*x) )*x5;
  
  /** l=3 */
  /* (3,1) */
  hTidallm[2] = ( -hA[2] + hB[2] )*x5;
  /* (3,3) */
  hTidallm[4] = ( -hA[4] + hB[4] )*x5;
  
}

/** Function providing a fit of Deltat_NQC vs chi, via a simple rational function. */
double eob_wav_dtnqc_fit(const double chi, const double chi0)
{
  const double n1 = -16.06288206;
  const double d1 = -4.04266459;
  double x     = chi-chi0;
  double dtnqc = (1.+n1*x)/(1.+d1*x);
  return dtnqc;  
}

/** Time-shift for NQC */
double eob_wav_timeshift_nqc(double nu, double chi1)
{

  double DeltaT_nqc = 1.; /* standard choice inspired by test-particle results */  

  if (chi1 >= 0.8498) {
    
    /* Interpolating fit for Deltat_NQC. See Eq.(21) of arXiv:1506.08457 
       This is a formula that was obtained in the equal-mass, equal-spin
       case and promoted also to any other case where the spin on the larger
       BH is larger than 0.8498. This is a guess to extrapolate the model
       outside the domain of calibration */
    
    DeltaT_nqc = eob_wav_dtnqc_fit(chi1,0.8498);

  } 

  if ((chi1 <=-0.80) && (nu <= 8./81.)) {
    
    /* This condition was a simple hack to avoid unphysical features in the
       modulus amplitude when one (or both) the spins are large and negative
       and the mass ratio is large. This little modification in the location
       of the NQC point guarantees that the determination of the NQC parameters
       guarantees just a small perturbation of the non-NQC EOB waveform. The
       iResum waveform will be robust enough that this hack will not be needed*/  
    
	DeltaT_nqc = 3.0;
  }    
  
  return DeltaT_nqc;  
}

/** Computes the factors and the coefficients that build the  
    NQC corrections to the waveform in the spinning case */
void eob_wav_hlmNQC_find_a1a2a3(const int size, Dynamics *dyn, Waveform_lm *h, Waveform_lm *hnqc)
{
  double A_tmp, dA_tmp, omg_tmp, domg_tmp;
    
  const double nu   = dyn->nu;
  const double X1   = dyn->X1;
  const double X2   = dyn->X2;
  const double chi1 = dyn->chi1;
  const double chi2 = dyn->chi2;
  const double aK   = dyn->a1 + dyn->a2;

  const double nu2  = SQ(nu);
  const double nu3  = nu2*nu;
  const double X12  = X1 - X2;  
  const double aK2  = SQ(aK);
  const double aK3  = aK2*aK;
  const double aK4  = aK2*aK2;
  const double a12  = X1*chi1 - X2*chi2;
  const double aeff     = aK + 1/3*a12*X12;
  const double aeff_omg = aK + a12*X12;
    
  double *T       = dyn->time;
  double *r       = dyn->data[EOB_RAD];
  double *w       = dyn->data[EOB_MOMG];
  double *pph     = dyn->data[EOB_PPHI];
  double *pr_star = dyn->data[EOB_PRSTAR];
  double *Omg_orb = dyn->data[EOB_OMGORB];
  double *ddotr   = dyn->data[EOB_DDOTR];

  double c_p1,     c_p2,     c_p3,   c_p4;
  double c_pdA1,   c_pdA2,   c_pdA3, c_pdA4;
  double c_pdomg1, c_pdomg2;
  double n0, d1;
  double a0_omg_tmp, a1_omg_tmp, a2_omg_tmp, b0_omg_tmp, b1_omg_tmp, b2_omg_tmp, a0_domg_tmp, a1_domg_tmp, a2_domg_tmp, b0_domg_tmp, b1_domg_tmp, b2_domg_tmp, a0_A_tmp, a1_A_tmp , a2_A_tmp, b0_A_tmp, b1_A_tmp, b2_A_tmp, a0_dA_tmp, a1_dA_tmp, a2_dA_tmp, b0_dA_tmp, b1_dA_tmp, b2_dA_tmp, omg_tmp_nu, omg_tmp_equal, domg_tmp_nu, domg_tmp_equal,  A_tmp_scale_nu, A_tmp_scale_equal, dA_tmp_scale_nu, dA_tmp_scale_equal ;
  
  double P[2], M[4], p1[2], p2[2], p3[2], p4[2], pA[5], pdA[5];    
  double pomg[5], pdomg[5], pn0[2], pd1[2], ppdomg1[2], ppdomg2[2], pdA1[2],pdA2[2],pdA3[2],pdA4[2];
  double n1[size], n2[size], n3[size], n4[size], n5[size], n6[size], d_n4[size], d_n5[size], d2_n4[size], d2_n5[size];
  double max_A[35],max_dA[35],d2max[35],d3max[35],max_omg[35],max_domg[35],maxd2omg[35], DeltaT[35];
    
  double A[KMAX][size];
  double phase[KMAX][size];

  double ai[KMAX][size];
  double bi[KMAX][size];
  double omg[KMAX][size];
  double domg[KMAX][size];
  double m11[KMAX][size];
  double m12[KMAX][size];
  double m13[KMAX][size];
  double m21[KMAX][size];
  double m22[KMAX][size];
  double p1tmp[KMAX][size];
  double p2tmp[KMAX][size];

  int k, j;

  /** Regge-Wheeler-Zerilli normalized amplitude. The ringdown
      coefficient refer to this normalization.
      Check Nagar & Rezzolla, CQG 22 (2005) R167 */           
  for (k=0; k<KMAX; k++) {
    for (j=0; j<size; j++) {
      A[k][j]     = h->ampli[k][j]/sqrt( (LINDEX[k]+2)*(LINDEX[k]+1)*LINDEX[k]*(LINDEX[k]-1) );
      phase[k][j] = h->phase[k][j];
    }
  }

  /** derivatives (nonuniform grid) */
  for (k=0; k<KMAX; k++) {
    D0_x(phase[k], T, size, omg[k]);
    D0_x(omg[k]  , T, size, domg[k]);
  }
  
  /**  Case 'NQC_fit_hybrid' */
  if (nu == 0.25) {

    pA[0]    =  0.00178195;
    pA[1]    =  0.00435589;
    pA[2]    =  0.00344489;
    pA[3]    = -0.00076165;
    pA[4]    =  0.31973334;
    A_tmp    =  pA[0]*aK4    + pA[1]*aK3   + pA[2]*aK2    + pA[3]*aK     + pA[4];
    
    pdA[0]   =  0.00000927;
    pdA[1]   = -0.00024550;
    pdA[2]   =  0.00012469;
    pdA[3]   =  0.00123845;
    pdA[4]   = -0.00195014;
    dA_tmp   =  pdA[0]*aK4   + pdA[1]*aK3   + pdA[2]*aK2   + pdA[3]*aK   + pdA[4];
    
    pomg[0]  =  0.00603482;
    pomg[1]  =  0.01604555;
    pomg[2]  =  0.02290799;
    pomg[3]  =  0.07084587;
    pomg[4]  =  0.38321834;
    omg_tmp  =  pomg[0]*aK4  + pomg[1]*aK3  + pomg[2]*aK2  + pomg[3]*aK  + pomg[4];
    
    pdomg[0] =  0.00024066;
    pdomg[1] =  0.00038123;
    pdomg[2] = -0.00049714;
    pdomg[3] =  0.00041219;
    pdomg[4] =  0.01190548;
    domg_tmp =  pdomg[0]*aK4 + pdomg[1]*aK3 + pdomg[2]*aK2 + pdomg[3]*aK + pdomg[4];
  
  }  else if( nu > 0.16) {
    
    p1[0]      =  0.04680896;
    p1[1]      = -0.00632114;
    p2[0]      =  0.06586192;
    p2[1]      = -0.01180039;
    p3[0]      = -0.11617413;
    p3[1]      =  0.02704959;
    p4[0]      =  0.15597465;
    p4[1]      =  0.28034978;
    c_p1       =  p1[0]*nu + p1[1];
    c_p2       =  p2[0]*nu + p2[1];
    c_p3       =  p3[0]*nu + p3[1];
    c_p4       =  p4[0]*nu + p4[1];
    A_tmp      =  c_p1*aK3 + c_p2*aK2 + c_p3*aK + c_p4;
    
    pdA1[0]    = -0.00130824;
    pdA1[1]    =  0.00006202;
    pdA2[0]    =  0.00199855;
    pdA2[1]    = -0.00027474;
    pdA3[0]    =  0.00218838;
    pdA3[1]    =  0.00071540;
    pdA4[0]    = -0.00362779;
    pdA4[1]    = -0.00105397;
    c_pdA1     =  pdA1[0]*nu + pdA1[1];
    c_pdA2     =  pdA2[0]*nu + pdA2[1];
    c_pdA3     =  pdA3[0]*nu + pdA3[1];
    c_pdA4     =  pdA4[0]*nu + pdA4[1];
    dA_tmp     =  c_pdA1*aK3   + c_pdA2*aK2 + c_pdA3*aK+ c_pdA4;
    
    pn0[0]     =  0.46908067;
    pn0[1]     =  0.27022141;
    pd1[0]     =  0.64131115;
    pd1[1]     = -0.37878384;
    n0         =  pn0[0]*nu + pn0[1];
    d1         =  pd1[0]*nu + pd1[1];
    omg_tmp    =  n0/(1 + d1*aK);
    
    ppdomg1[0] =  0.00061175;
    ppdomg1[1] =  0.00074001;
    ppdomg2[0] =  0.02504442;
    ppdomg2[1] =  0.00548217;
    c_pdomg1   =  ppdomg1[0]*nu + ppdomg1[1];
    c_pdomg2   =  ppdomg2[0]*nu + ppdomg2[1];
    domg_tmp   =  c_pdomg1*aK   + c_pdomg2;
  
  }  else {

    /* This fit is the most updated and recent that Gunnar did by 
       incorporating in the fits also the test-particle NQC point
       obtained from the most-recent Teukolsky waveforms done by
       M. Colleoni using the 6PN-accurare iResum-radiation reaction.
       These points assure a smooth connection between merger and
       ringdown also outside the "calibration" domain, notably for
       large-mass ratios (though q<=20) and large (negative) spins
       Updated, 28/09/2017 */
    
    a0_omg_tmp    = -0.1460961247;
    a1_omg_tmp    =  0.0998056;
    a2_omg_tmp    = -0.118098;
    b0_omg_tmp    = -0.3430184009;
    b1_omg_tmp    =  0.0921551;
    b2_omg_tmp    = -0.0740285;
    omg_tmp_nu    = +0.5427169903*nu2 +0.2512395608*nu +0.2863992248;
    omg_tmp_equal =((a2_omg_tmp*X12*X12 + a1_omg_tmp*X12 + a0_omg_tmp)*aeff_omg+1)/((b2_omg_tmp*X12*X12 +b1_omg_tmp*X12 + b0_omg_tmp)*aeff_omg+1);
    omg_tmp       = omg_tmp_nu*omg_tmp_equal;

    a0_domg_tmp    = +0.0604556289;
    b0_domg_tmp    = -0.0299583285;
    a1_domg_tmp    = 0.0711715;
    a2_domg_tmp    = -0.0500886;
    b1_domg_tmp    = 0.0461239;
    b2_domg_tmp    = -0.0153068;
    
    domg_tmp_nu    = ( +0.0045213831*nu +0.0064934920)/( -1.4466409969*nu+1);
    domg_tmp_equal = (a2_domg_tmp*X12*X12 +a1_domg_tmp*X12 +b0_domg_tmp)*aeff_omg*aeff_omg +(b2_domg_tmp*X12*X12 +b1_domg_tmp*X12+a0_domg_tmp)*aeff_omg+1;
    domg_tmp       = domg_tmp_nu*domg_tmp_equal;
    
    a0_A_tmp 	= -0.2750516062;
    b0_A_tmp 	= -0.4693776065;
    a1_A_tmp 	=  0.143066;
    a2_A_tmp 	= -0.0425947;
    b1_A_tmp 	=  0.176955;
    b2_A_tmp 	= -0.111902;
    
    A_tmp_scale_nu    = -0.9862040409*nu3 +0.8167558040*nu2 -0.0427442282*nu+0.2948879452;
    A_tmp_scale_equal = ((a2_A_tmp*X12*X12 + a1_A_tmp*X12 +a0_A_tmp)*aeff+1)/((b2_A_tmp*X12*X12 + b1_A_tmp*X12 +b0_A_tmp)*aeff+1);
    A_tmp             = A_tmp_scale_nu*A_tmp_scale_equal*(1-0.5*omg_tmp*aeff);
    
    a0_dA_tmp 	= +0.0037461628;
    b0_dA_tmp 	= +0.0636082543;
    a1_dA_tmp 	=  0.00129393;
    a2_dA_tmp 	= -0.00239069;
    b1_dA_tmp 	= -0.0534209;
    b2_dA_tmp 	= -0.186101;
    
    dA_tmp_scale_nu    = ( -0.0847947167*nu -0.0042142765)/( +16.1559461812*nu+1);
    dA_tmp_scale_equal = ((a2_dA_tmp*X12*X12 + a1_dA_tmp*X12+ a0_dA_tmp)*aeff)/((b2_dA_tmp*X12*X12 + b1_dA_tmp*X12 + b0_dA_tmp)*aeff+1);
    dA_tmp             = (dA_tmp_scale_nu +dA_tmp_scale_equal)*omg_tmp;
    
  }

  for (k=0; k<KMAX; k++) {   
    max_A[k]    = 0.;
    max_dA[k]   = 0.;
    max_omg[k]  = 0.;
    max_domg[k] = 0.;
  }
  
  /** Switch on the 22  values (only) */
  max_A[1]    = A_tmp;
  max_dA[1]   = dA_tmp;
  max_omg[1]  = omg_tmp;
  max_domg[1] = domg_tmp;

#if (DEBUG)
  printf("--------------------------------\n");
  printf("NR values for NQC determination:\n");
  printf("--------------------------------\n");
  printf("Amrg    =%10.6f\n",max_A[1]);
  printf("dAmrg   =%10.6f\n",max_dA[1]);
  printf("omg_mrg =%10.6f\n",max_omg[1]);
  printf("domg_mrg=%10.6f\n",max_domg[1]);
#endif
  
  /** NQC corrections to AMPLITUDE (n1,n2,n3) and PHASE (n4,n5,n6)
   * NQC basis for (2,2) waveform : AMPLITUDE
   * note: n3 and n6 are not used
   */
  double pr_star2, r2, w2;

  for (j=0; j<size; j++) {
    pr_star2 = pr_star[j] * pr_star[j];
    r2       = r[j] * r[j];
    w2       = w[j] * w[j];
    n1[j]  = pr_star2/(r2*w2);         /* [pr*\/(r Omg)]^2 */
    n2[j]  = ddotr[j]/(r[j]*w2);       /* [ddot{r}/(r Omg^2)] */
    n4[j]  = pr_star[j]/(r[j]*w[j]);   /* pr*\/(r Omg) */
    n5[j]  = n4[j]*r2*w2;              /* (pr*)*(r Omg) */
  }
    
#if (DEBUG)
  FILE* fp = fopen("NQC_func.txt", "w");
  for (j=0; j<size; j++) {
    fprintf(fp, "%20.12f\t%20.12f\t%20.12f\t%20.12f\t%20.12f\n", T[j], n1[j], n2[j], n4[j], n5[j]);
  }
  fclose(fp);
#endif    

  /** Derivatives for the phase */
  D0_x(n4,T,size, d_n4);
  D0_x(n5,T,size, d_n5);
  D0_x(d_n4,T,size, d2_n4);
  D0_x(d_n5,T,size, d2_n5);

#if (DEBUG)
  fp = fopen("dNQC_func.txt", "w");
  for (j=0; j<size; j++) {
    fprintf(fp, "%f\t%f\t%f\t%f\t%f\n", T[j], d_n4[j], d_n5[j], d2_n4[j], d2_n5[j]);  
  }
  fclose(fp);
#endif    

  /** Find max Omg */
  int Omgmax_index = 0;
  double Omg_max   = Omg_orb[0];
  for (j=0; j<size; j++) {
    if (Omg_orb[j] > Omg_max) {
      Omg_max = Omg_orb[j];
      Omgmax_index = j;
    }
  }

  /** Time */
  double tOmgOrb_pk = T[Omgmax_index];
  double DeltaT_nqc = eob_wav_timeshift_nqc(nu, chi1);
  double tNQC = tOmgOrb_pk - DeltaT_nqc;

#if (DEBUG)
  printf("-------------------------------\n");
  printf("Check: NQC related information:\n");
  printf("-------------------------------\n");
  printf("DeltaT_tNQC = %f\n",DeltaT_nqc);
  printf("tNQC [bare] = %f\n",tNQC);
#endif

  /** Find jmax: T[jmax] <= tNQC */
  int jmax = 0;
  for (j=0; j<size; j++) {
    if(T[j] > tNQC) break;
  }
  jmax = j-1;

  /** Solve the linear systems */
  for (k=0; k<KMAX; k++) {
    for (j=0; j<size; j++) {
      /* Matrix elements: waveform amplitude at all points */
      m11[k][j] = n1[j]*A[k][j];
      m12[k][j] = n2[j]*A[k][j];
      p1tmp[k][j] = A[k][j];      
    }
    /* Take FD derivatives */
    D0_x(m11[k],T,size, m21[k]);
    D0_x(m12[k],T,size, m22[k]);
    D0_x(A[k],T,size, p2tmp[k]);
  }


#if (DEBUG)
  fp = fopen("Amp_func.txt", "w");
  for (j=0; j<size; j++) {
    fprintf(fp, "%f\t%f\t%f\n", T[j], p1tmp[1][j], p2tmp[1][j]);
  }
  fclose(fp);  
#endif

  double detM = 1.;
  for (k=0; k<KMAX; k++) {
    
    /* Computation of ai coefficients */
    P[0]     = max_A[k]  - p1tmp[k][jmax];
    P[1]     = max_dA[k] - p2tmp[k][jmax];
    
    M[0]     = m11[k][jmax];
    M[1]     = m12[k][jmax];
    M[2]     = m21[k][jmax];
    M[3]     = m22[k][jmax];
    
    detM     = M[0]*M[3]-M[1]*M[2];
    ai[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
    ai[k][1] = (M[0]*P[1] - M[2]*P[0])/detM;
    
    /* Computation of bi coefficients */
    P[0]     = omg[k][jmax]   - max_omg[k];
    P[1]     = domg[k][jmax]  - max_domg[k];
    
    M[0]     = d_n4[jmax];
    M[1]     = d_n5[jmax];
    M[2]     = d2_n4[jmax];
    M[3]     = d2_n5[jmax];
    
    detM     =  M[0]*M[3] - M[1]*M[2];
    bi[k][0] = (M[3]*P[0] - M[1]*P[1])/detM;
    bi[k][1] = (M[0]*P[1] - M[2]*P[0])/detM;
  
  }
  
#if (DEBUG)
  printf("-------------------\n");
  printf("NQC coefficients:  \n");
  printf("-------------------\n");
  printf("a1 = %f\n",ai[1][0]);
  printf("a2 = %f\n",ai[1][1]);
  printf("b1 = %f\n",bi[1][0]);
  printf("b2 = %f\n",bi[1][1]);
#endif    

  /** Set amplitude and phase */
  for (k=0; k<KMAX; k++) {
    for (j=0; j<size; j++) {
      hnqc->ampli[k][j] = 1. + ai[k][0]*n1[j] + ai[k][1]*n2[j];
      hnqc->phase[k][j] =      bi[k][0]*n4[j] + bi[k][1]*n5[j];
    }
  }

  /** Multiply waveform to NQC */
 
  // FIXME: RWZ normalization? Don't we need to re-introduce it here ???

  for (k=0; k<KMAX; k++) {
    for (j=0; j<size; j++) {
      h->ampli[k][j] *= hnqc->ampli[k][j];
      h->phase[k][j] *= hnqc->phase[k][j];
    }
  }

}

/** NQC corrections to the RWZ multipolar waveform
    Nagar, Damour, Reisswig, Pollney http://arxiv.org/abs/1506.08457 
    Current fits: 9/02/2016 */
void eob_wav_hlmNQC(double  nu, double  r, double  prstar, double  Omega, double  ddotr,
		    Waveform_lm_t *psilmnqc)
{      
  const double xnu  = 1-4*nu;
  const double xnu2 = SQ(xnu);

  double a1[KMAX], a2[KMAX], a3[KMAX];
  double b1[KMAX], b2[KMAX], b3[KMAX];
  double n[6];    

  /* NQC corrections to the modulus */
  n[0] = (prstar/(r*Omega))*(prstar/(r*Omega));
  n[1] = ddotr/(r*Omega*Omega);
  n[2] = n[0]*prstar*prstar;

  /* NQC corrections to the phase */
  n[3] = prstar/(r*Omega);
  n[4] = n[3]*cbrt(Omega*Omega);
  n[5] = n[4]*prstar*prstar;

  for (int k = 0; k < KMAX; k++) {
    psilmnqc->ampli[k] = 1.; 
    psilmnqc->phase[k] = 0.; 
  }

  /*
  for (int k = 0; k < KMAX; k++) {
    a1[k] = 0.;
    a2[k] = 0.;
    a3[k] = 0.;
    b1[k] = 0.;
    b2[k] = 0.;
    b3[k] = 0.;
  }
  */

  /* Coefficients */
	 
  /* (2,1) */
  a1[0] =  0.0162387198*(7.32653082*xnu2 + 1.19616248*xnu + 0.73496656);
  a2[0] =                -1.80492460*xnu2 + 1.78172686*xnu + 0.30865284;
  a3[0] =                                                           0.0;
  
  b1[0] =  -0.0647955017*(3.59934444*xnu2 - 4.08628784*xnu + 1.37890907);
  b2[0] =   1.3410693180*(0.38491989*xnu2 + 0.10969453*xnu + 0.97513971);
  b3[0] =                                                            0.0;
      
  /* (2,2) */
  a1[1]   = -0.0805236959*( 1 - 2.00332326*xnu2)/( 1 + 3.08595088*xnu2);
  a2[1]   =  1.5299534255*( 1 + 1.16438929*xnu2)/( 1 + 1.92033923*xnu2);
  a3[1]   =  0.0;
  
  b1[1]   = 0.146768094955*( 0.07417121*xnu + 1.01691256);
  b2[1]   = 0.896911234248*(-0.61072011*xnu + 0.94295129);
  b3[1]   = 0.0;
      
  /* (3,3) */
  a1[4]   = -0.0377680000*(1 - 14.61548907*xnu2)/( 1 + 2.44559263*xnu2);
  a2[4]   =  1.9898000000*(1 + 2.09750346 *xnu2)/( 1 + 2.57489466*xnu2);
  a3[4]   =  0.0;
  
  b1[4]   = 0.1418400000*(1.07430512 - 1.23906804*xnu + 4.44910652*xnu2);
  b2[4]   = 0.6191300000*(0.80672432 + 4.07432829*xnu - 7.47270977*xnu2);
  b3[4]   = 0.0;
  
  /* NQC factor */
  int k = 0; /* (2,1) */
  psilmnqc->ampli[k] = 1. + a1[k]*n[0] + a2[k]*n[1] + a3[k]*n[2];
  psilmnqc->phase[k] =      b1[k]*n[3] + b2[k]*n[4] + b3[k]*n[5];
  
  k = 1; /* (2,2) */
  psilmnqc->ampli[k] = 1. + a1[k]*n[0] + a2[k]*n[1] + a3[k]*n[2];
  psilmnqc->phase[k] =      b1[k]*n[3] + b2[k]*n[4] + b3[k]*n[5];
  
  k = 4; /* (3,3) */
  psilmnqc->ampli[k] = 1. + a1[k]*n[0] + a2[k]*n[1] + a3[k]*n[2];
  psilmnqc->phase[k] =      b1[k]*n[3] + b2[k]*n[4] + b3[k]*n[5];
  
}

/** Ringdown waveform template */
void eob_wav_ringdown_template(double x, double a1, double a2, double a3, double a4, double b1, double b2, double b3, double b4, double *sigma, double *psi)
{  
  double amp   = ( a1 * tanh(a2*x +a3) + a4 ) ;
  double phase = -b1*log((1. + b3*exp(-b2*x) + b4*exp(-2.*b2*x))/(1.+b3+b4));   
  psi[0] = amp * exp(-sigma[0]*x); /* amplitude */
  psi[1] = - (phase - sigma[1]*x); /* phase, minus sign in front by convention */
}

/** Ringdown calculation and match to the dynamics */ 
void eob_wav_ringdown(Dynamics *dyn, Waveform_lm *hlm)
{
  const double Mbh   = dyn->Mbhf;
  const double abh   = dyn->abhf;

  const double nu    = dyn->nu;
  const double chi1  = dyn->chi1;
  const double chi2  = dyn->chi2;
  const double X1    = dyn->X1;
  const double X2    = dyn->X2;
  const double aK    = dyn->a1+dyn->a2;

  double *t     = dyn->time;
  double *Omega = dyn->data[EOB_MOMG];

  const double xnu   = (1.-4.*nu);
  const double ooMbh = 1./Mbh;
  
  const int size = dyn->size;

  const int k21 = 0;
  const int k22 = 1;
  const int k33 = 4;
  int k, j;
  
  /** Find peak of Omega */
  /* Assume it is a monotonically increasing function */
  int index_pk = 0;
  double Omega_pk = Omega[index_pk];
  for (j = 1; j < size ; j++ ) {
      if (Omega[j] < Omega_pk) 
	break;
      index_pk = j;
      Omega_pk = Omega[j];
  }
  if (index_pk == size-1) {
    if (VERBOSE) printf("No omega-maximum found.\n");
    //index_pk = size - 4;
  }
  if (index_pk > size-4) {
    errorexit("Not enough points to interpolate.\n");
  }

  /*
  int j = size;
  double Omega_pk = Omega[j-1];
  while (Omega[j] > Omega_pk) {
    index_pk = j;
    Omega_pk = Omega[j];
    j--;
  }
  if (index_pk == -1) {
    if (VERBOSE) printf("No omega-maximum found.\n");
    index_pk = size - 4;
  }
  if (index_pk > size-4) {
    errorexit("Not enough points to interpolate.\n");
  }
  */

  // NOTE: following is slightly different from C++ thing (that I do not understand).
  // here we just refine the 7 pts grid, populate by spline, and look for a maximum on that
  
#define n_grid (7)
  double *Omega_pk_grid, *t_Omega_pk_grid;
  Omega_pk_grid = &Omega[index_pk-3];
  t_Omega_pk_grid = &t[index_pk-3];

#define n_refine (21)
  double dt = (t_Omega_pk_grid[n_grid-1] - t_Omega_pk_grid[0])/(n_refine-1); 
  double ti[n_refine],oi[n_refine];
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_cspline, 7);
  gsl_spline_init (spline, t_Omega_pk_grid, Omega_pk_grid, n_grid);  
  for (j = 0; j < n_refine; j++) {
    ti[j] = t_Omega_pk_grid[0] + j*dt;
    oi[j] = gsl_spline_eval (spline, ti[j], acc);
  }
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  index_pk = 0;
  Omega_pk = oi[0];
  for (j = 0; j < n_refine; j++) {
    if (oi[j] > Omega_pk) {
      index_pk = j;
      Omega_pk = oi[j];
    }
  }
  
  double tOmg_pk = ti[index_pk];
  tOmg_pk *= ooMbh;

  /** Merger time t_max(A22) */
  double DeltaT_nqc = eob_wav_timeshift_nqc(nu, chi1);
  double tmrg[KMAX], tmatch[35], dtmrg[2];
            
  /** nonspinning case */ // OLD
  /* tmrg[k22]  = tOmg_pk-3./Mbh; */          
     
  /** nonspinning case */
  tmrg[k22]  = tOmg_pk-(DeltaT_nqc + 2)/Mbh;     
  dtmrg[k21] = 5.70364338 + 1.85804796*xnu  + 4.0332262*xnu*xnu; //k21
  dtmrg[k33] = 4.29550934 - 0.85938*xnu;                         //k33
  tmrg[k21]  = tmrg[k22] + dtmrg[k21]/Mbh;     // t_max(A21) => peak of 21 mode
  tmrg[k33]  = tmrg[k22] + dtmrg[k33]/Mbh;     // t_max(A33) => peak of 33 mode

  /** postmerger-ringdown matching time */
  for (k=0; k<KMAX; k++) {
    tmatch[k] = 2.*ooMbh + tmrg[k];
  }

  /** Compute QNM */
  double **sigma;//[KMAX][2]; // real, imag
  *sigma = malloc ( KMAX * sizeof(double*) );
  for (k=0; k<KMAX; k++) 
    sigma[k] = malloc ( 2 * sizeof(double) );

  double a1[KMAX], a2[KMAX], a3[KMAX], a4[KMAX];
  double b1[KMAX], b2[KMAX], b3[KMAX], b4[KMAX]; 
  QNMHybridFitCab(nu, X1, X2, chi1, chi2, aK,  Mbh, abh,  
		  a1, a2, a3, a4, b1, b2, b3, b4, 
		  sigma);

  /** Define a time vector for each multipole, scale by mass
      Ringdown of each multipole has its own starting time */
  double *t_lm[KMAX];
  for (k=0; k<KMAX; k++) {
    t_lm[k] =  malloc ( size * sizeof(double) );
    for (j = 0; j < size; j++ ) {  
      t_lm[k][j] = t[j] * ooMbh;
    }
  }

  /** Find attachment index */
  int idx[KMAX];
  for (k = 0; k < KMAX; k++) {
    for (j = 0; j < size ; j++ ) {  
      if (t_lm[k][j] >= tmatch[k]) break;
    }
    idx[k] = j;
  }

  /** Calculate deltaphi */
  double t0, psi[2];
  double Deltaphi[KMAX];
  /*
  for (k = 0; k < KMAX; k++) {
    t0 = t_lm[k][idx[k]] - tmrg[k]; // I or I-1 or I-2 ?
    eob_wav_ringdown_template(t0, a1[k], a2[k], a3[k], a4[k], b1[k], b2[k], b3[k], b4[k], sigma[k], psi);
    Deltaphi[k] = psi[1] - hlm[k]->phase[idx[k]];
  }
  */

  /** Compute Ringdown waveform for t>=tmatch */
  for (k = 0; k < KMAX; k++) {
    for (j = idx[k]; j < size ; j++ ) {  
      //t0 = t_lm[k][j];
      t0 = t_lm[k][j] - tmrg[k];   
      eob_wav_ringdown_template(t0, a1[k], a2[k], a3[k], a4[k], b1[k], b2[k], b3[k], b4[k], sigma[k], psi);
      if (j==idx[k]) {
	Deltaphi[k] = psi[1] - hlm->phase[k][idx[k]];
      }
      hlm->phase[k][j] = psi[1] - Deltaphi[k];
      hlm->ampli[k][j] = psi[0];
    }
  }
 
  /** Free mem. */
  for (k=0; k<KMAX; k++) {
    free(sigma[k]);
    free(t_lm[k]);
  }
  free(sigma);

}


/** Main routine for factorized EOB waveform */
void eob_wav_hlm(Dynamics *dyn, Waveform_lm_t *hlm)
{

  const double nu = dyn->nu;  
  const double chi1 = dyn->chi1;  
  const double chi2 = dyn->chi2;  
  const double a1 = dyn->a1;  
  const double a2 = dyn->a2;  
  const double X1 = dyn->X1;  
  const double X2 = dyn->X2;  
  const double C_Q1 = dyn->C_Q1;  
  const double C_Q2 = dyn->C_Q2;  
  const int usetidal = dyn->use_tidal;
  const int usespins = dyn->use_spins;
  const int usespeedytail = par_get_i("use_speedytail");
  const double X12 = X1-X2; /* sqrt(1-4nu) */

  const double t   = dyn->t;
  const double phi = dyn->phi; 
  const double r   = dyn->r;
  const double pph = dyn->pphi;
  const double prstar = dyn->prstar;
  const double Omega  = dyn->Omg;
  const double ddotr  = dyn-> ddotr;
  const double H      = dyn->H;
  const double Heff   = dyn->Heff;
  const double jhat   = dyn->jhat;
  const double rw     = dyn->r_omega;

  hlm->time = t;

  /** Source term */
  double source[] = {
    jhat,Heff,
    Heff,jhat,Heff,
    jhat,Heff,jhat,Heff,
    Heff,jhat,Heff,jhat,Heff,
    jhat,Heff,jhat,Heff,jhat,Heff,
    Heff,jhat,Heff,jhat,Heff,jhat,Heff,
    jhat,Heff,jhat,Heff,jhat,Heff,jhat,Heff
  };
  
  /** Newtonian waveform */
  Waveform_lm_t hNewt;
  eob_wav_hlmNewt(rw,Omega,phi,nu, &hNewt);

  /* if (usetidal) { */
  /* Need to correct some of the m=odd modes. 
     The Newtonian factor has a different normalization when entering the point-mass 
     and the tidal term. The factor X12 = sqrt*1-4nu) is re-introduced in the point-mass term 
     in eob_wav_hlm() */
  /* double vphi3 = gsl_pow_int(rw*Omega,3); */
  /* hNewt.ampli[0] = ChlmNewt_ampli[0] * vphi3; /\* (2,1) *\/ */
  /* hNewt.ampli[2] = ChlmNewt_ampli[2] * vphi3; /\* (3,1) *\/ */
  /* hNewt.ampli[4] = ChlmNewt_ampli[4] * vphi3; /\* (3,3) *\/ */
  /* double p4_vphi5 = (2.*nu-1) * gsl_pow_int(rw*Omega,5); */
  /* hNewt.ampli[5]  = ChlmNewt_ampli[5]  * p4_vphi5; /\* (4,1) *\/ */
  /* hNewt.ampli[7]  = ChlmNewt_ampli[7]  * p4_vphi5; /\* (4,3) *\/ */
  /* hNewt.ampli[9]  = ChlmNewt_ampli[9]  * p4_vphi5; /\* (5,1) *\/ */
  /* hNewt.ampli[11] = ChlmNewt_ampli[11] * p4_vphi5; /\* (5,3) *\/ */
  /* hNewt.ampli[13] = ChlmNewt_ampli[13] * p4_vphi5; /\* (5,5) *\/ */
  /* } */

  if (usespins) {
    /* Special treatment when spin is on because of the singularity in the sqrt(1-4*nu) 
       for m=odd mode and nu=1/4. See discussion in 
       Damour & Nagar, PRD 90, 044018, Sec. 4, Eq.(89). 
       This is not done for multipoles l>4 because no spinning information is included there */ 
    
    double vphi3 = gsl_pow_int(rw*Omega,3);
    hNewt.ampli[0] = ChlmNewt_ampli[0] * vphi3; /* (2,1) */
    hNewt.ampli[2] = ChlmNewt_ampli[2] * vphi3; /* (3,1) */
    hNewt.ampli[4] = ChlmNewt_ampli[4] * vphi3; /* (3,3) */

    //FIXME: Following needs check in the case spin+tides
    double p4_vphi5 = (2.*nu-1) * gsl_pow_int(rw*Omega,5);
    hNewt.ampli[5]  = ChlmNewt_ampli[5]  * p4_vphi5; /* (4,1) */
    hNewt.ampli[7]  = ChlmNewt_ampli[7]  * p4_vphi5; /* (4,3) */
    hNewt.ampli[9]  = ChlmNewt_ampli[9]  * p4_vphi5 * X12; /* (5,1) */
    hNewt.ampli[11] = ChlmNewt_ampli[11] * p4_vphi5 * X12; /* (5,3) */
    hNewt.ampli[13] = ChlmNewt_ampli[13] * p4_vphi5 * X12; /* (5,5) */
  }

  /** Compute corrections */
  double rholm[KMAX], flm[KMAX];
  double x = SQ(rw*Omega);
  if (usespins){
    /* eob_wav_flm_s_old(x, nu, X1,X2,chi1,chi2,a1,a2,C_Q1,C_Q2,usetidal,rholm,flm); */ 
    eob_wav_flm_s(x, nu, X1,X2,chi1,chi2,a1,a2,C_Q1,C_Q2,usetidal,rholm,flm); 
  } else {
    /* eob_wav_flm_old(x, nu, rholm,flm); */
    eob_wav_flm(x, nu, rholm,flm);
  }
  
  /** Computing the tail */
#define r0 (1.213061319425267e+00)
  const double Hreal = H * nu;
  Waveform_lm_t tlm;
  if (usespeedytail) {
    eob_wav_speedyTail(Omega,Hreal, r0, &tlm); 
  } else {
    eob_wav_hhatlmTail(Omega,Hreal, r0, &tlm); 
  }
  
  /** Residual phase corrections delta_{lm} */
  double dlm[KMAX];
  eob_wav_deltalm(Hreal, Omega, nu, dlm); 
  
  /** NQC */
  Waveform_lm_t hNQC; 
  if (!(usetidal)) {
    eob_wav_hlmNQC(nu,r,prstar,Omega,ddotr, &hNQC); 
  }

  /** Point-mass h_lm */
  for (int k = 0; k < KMAX; k++) {
    hlm->ampli[k] =  hNewt.ampli[k] * flm[k] * source[k] * tlm.ampli[k];
    hlm->phase[k] = -( hNewt.phase[k] + tlm.phase[k] + dlm[k]); /* Minus sign by convention */
  }

  if ( (!(usetidal)) && (!(usespins)) ) {
    /** Add NQC correction */
    for (int k = 0; k < KMAX; k++) {
      hlm->ampli[k] *= hNQC.ampli[k];
      hlm->phase[k] -= hNQC.phase[k];
    }
  }

  if (usetidal) {   
    /** Tidal contribution */
    double hlmtidal[KMAX];
    eob_wav_hlmTidal(x, dyn, hlmtidal);
    
    /* Correct normalization of point-mass wave for some of the m=odd modes */
    hlm->ampli[0] *= X12;
    hlm->ampli[2] *= X12;
    hlm->ampli[4] *= X12;

    /* Add tidal contribution to waveform */
    for (int k = 0; k < KMAX; k++) {
      hlm->ampli[k] += (hNewt.ampli[k] * tlm.ampli[k] * hlmtidal[k]);
    }

  }
  
}

/** Old flux routines */

/** Resummed amplitudes in the general nu-dependent case.
 *  Refs:
 *  . Damour, Iyer & Nagar, PRD 79, 064004 (2009)     [theory]
 *  . Fujita & Iyer, PRD 82, 044051 (2010)            [test-mass 5.5PN]
 *  . Damour, Nagar & Bernuzzi, PRD 87, 084035 (2013) [complete information]
 * Old version, with coefficients computed all the times.
 */
void eob_wav_flm_old(double x,double nu, double *rholm, double *flm)
{
  
  /** Shorthands */
  const double x2  = x*x;
  const double x3  = x*x2;
  const double x4  = x*x3;
  const double x5  = x*x4;
  const double nu2 = nu*nu;
  const double nu3 = nu*nu2;
  const double nu4 = nu*nu3;
  
  /** Compute EulerLogs */
  const double el1 = Eulerlog(x,1);
  const double el2 = Eulerlog(x,2);
  const double el3 = Eulerlog(x,3);
  const double el4 = Eulerlog(x,4);
  const double el5 = Eulerlog(x,5);
  const double el6 = Eulerlog(x,6);
  const double el7 = Eulerlog(x,7);
  
  /** l=2 ------------------------------------------------------------------
   *  (2,1) */
  rholm[0] = 1.                                                                           +
    (-1.0535714285714286 + 0.27380952380952384 *nu                          )*x  +
    (-0.8327841553287982 - 0.7789824263038548  *nu + 0.13116496598639457*nu2)*x2 +
    (2.9192806270460925  - 1.019047619047619   *el1                         )*x3 +
    (-1.28235780892213   + 1.073639455782313   *el1                         )*x4 +
    (-3.8466571723355227 + 0.8486467106683944  *el1                     )*PMTERMS_eps*x5 ;
  
  /** (2,2) */
  rholm[1] = 1.                                                                                                                         +
    (-1.0238095238095237 + 0.6547619047619048*nu                                                                          )*x  +
    (-1.94208238851096   - 1.5601379440665155*nu + 0.4625614134542706*nu2                                                 )*x2 +
    (12.736034731834051  - 2.902228713904598 *nu - 1.9301558466099282*nu2 + 0.2715020968103451*nu3 - 4.076190476190476*el2)*x3 +
    (-2.4172313935587004 + 4.173242630385488 *el2                                                                         )*x4 +
    (-30.14143102836864  + 7.916297736025627 *el2                                                                         )*x5 ;
  
  /** l=3 ------------------------------------------------------------------
   *  (3,1) */
  rholm[2] = 1. + (-0.7222222222222222 - 0.2222222222222222*nu)*x + (0.014169472502805836 - 0.9455667789001122*nu - 0.46520763187429853*nu2)*x2 + x3*(1.9098284139598072 - 0.4126984126984127*el1) + x4*(0.5368150316615179 + 0.2980599647266314*el1) + PMTERMS_eps*x5*(1.4497991763035063 - 0.0058477188106817735*el1);
  /** (3,2) */
  rholm[3] = 1. + (0.003703703703703704*(328. - 1115.*nu + 320.*nu2)*x)/(-1. + 3.*nu) + (6.235191420376606e-7*(-1.444528e6 + 8.050045e6*nu - 4.725605e6*nu2 - 2.033896e7*nu3 + 3.08564e6*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + x3*(6.220997955214429 - 1.6507936507936507*el2) + PMTERMS_eps*x4*(-3.4527288879001268 + 2.005408583186361*el2);
  /** (3,3) */
  rholm[4] = 1. + (-1.1666666666666667 + 0.6666666666666666*nu)*x + (-1.6967171717171716 - 1.8797979797979798*nu + 0.45151515151515154*nu2)*x2 + x3*(14.10891386831863 - 3.7142857142857144*el3) + x4*(-6.723375314944128 + 4.333333333333333*el3) + PMTERMS_eps*x5*(-29.568699895427518 + 6.302092352092352*el3);
  
  /** l=4 ------------------------------------------------------------------
   *  (4,1) */
  rholm[5] = 1. + (0.001893939393939394*(602. - 1385.*nu + 288.*nu2)*x)/(-1. + 2.*nu) - 0.36778992787515513*x2 + x3*(0.6981550175535535 - 0.2266955266955267*el1) + PMTERMS_eps*x4*(-0.7931524512893319 + 0.2584672482399755*el1);
  /** (4,2) */
  rholm[6] = 1. + (0.0007575757575757576*(1146. - 3530.*nu + 285.*nu2)*x)/(-1. + 3.*nu) - (3.1534122443213353e-9*(1.14859044e8 - 2.95834536e8*nu - 1.204388696e9*nu2 + 3.04798116e9*nu3 + 3.79526805e8*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + 4.550378418934105e-12*x3*(8.48238724511e11 - 1.9927619712e11*el2) + PMTERMS_eps*x4*(-0.6621921297263365 + 0.787251738160829*el2);
  /** (4,3) */
  rholm[7] = 1. + (0.005681818181818182*(222. - 547.*nu + 160.*nu2)*x)/(-1. + 2.*nu) - 0.9783218202252293*x2 + PMTERMS_eps*(x3*(8.519456157072423 - 2.0402597402597404*el3) +      x4*(-5.353216984886716 + 2.5735094451003544*el3));
  /** (4,4) */
  rholm[8] = 1. + (0.0007575757575757576*(1614. - 5870.*nu + 2625.*nu2)*x)/(-1. + 3.*nu) + (3.1534122443213353e-9*(-5.11573572e8 + 2.338945704e9*nu - 3.13857376e8*nu2 - 6.733146e9*nu3 + 1.252563795e9*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + x3*(15.108111214795123 - 3.627128427128427*el4) + PMTERMS_eps*x4*(-8.857121657199649 + 4.434988849534304*el4);
  
  /** l=5 ------------------------------------------------------------------
   *  (5,1) */
  rholm[9] = 1. + (0.002564102564102564*(319. - 626.*nu + 8.*nu2)*x)/(-1. + 2.*nu) - 0.1047896120973044*x2 + PMTERMS_eps*(x3*(0.642701885362399 - 0.14414918414918415*el1) +      x4*(-0.07651588046467575 + 0.11790664036817883*el1));
  /** (5,2) */
  rholm[10] = 1. + (0.00007326007326007326*(-15828. + 84679.*nu - 104930.*nu2 + 21980.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + PMTERMS_eps*(- 0.4629337197600934*x2 +      x3*(2.354458371550237 - 0.5765967365967366*el2));
  /** (5,3) */
  rholm[11] = 1. + (0.002564102564102564*(375. - 850.*nu + 176.*nu2)*x)/(-1. + 2.*nu) - 0.5788010707241477*x2 + PMTERMS_eps*(x3*(5.733973288504755 - 1.2973426573426574*el3) +      x4*(-1.9573287625526001 + 1.2474448628294783*el3));
  /** (5,4) */
  rholm[12] = 1. + (0.00007326007326007326*(-17448. + 96019.*nu - 127610.*nu2 + 33320.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + PMTERMS_eps*(- 1.0442142414362194*x2 +   x3*(10.252052781721588 - 2.3063869463869464*el4));
    /**  (5,5) */
  rholm[13] = 1. + (0.002564102564102564*(487. - 1298.*nu + 512.*nu2)*x)/(-1. + 2.*nu) - 1.5749727622804546*x2 + PMTERMS_eps*(x3*(15.939827047208668 - 3.6037296037296036*el5) +      x4*(-10.272578060123237 + 4.500041838503377*el5));
  
  /** l=6 ------------------------------------------------------------------
   *  (6,1) */
  rholm[14] = 1. + (0.006944444444444444*(-161. + 694.*nu - 670.*nu2 + 124.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 0.29175486850885135*x2 + x3*(0.21653486654395454 - 0.10001110001110002*el1));
  /** (6,2) */
  rholm[15] = 1. + (0.011904761904761904*(-74. + 378.*nu - 413.*nu2 + 49.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + PMTERMS_eps*( - 0.24797525070634313*x2 + x3*(1.7942694138754138 - 0.40004440004440006*el2));
  /** (6,3) */
  rholm[16] = 1. + (0.006944444444444444*(-169. + 742.*nu - 750.*nu2 + 156.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 0.5605554442947213*x2 + x3*(4.002558222882566 - 0.9000999000999002*el3));
  /** 6,4)  */
  rholm[17] = 1. + (0.011904761904761904*(-86. + 462.*nu - 581.*nu2 + 133.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + PMTERMS_eps*(- 0.7228451986855349*x2 + x3*(7.359388663371044 - 1.6001776001776002*el4));
  /** 6,5)  */
  rholm[18] = 1. + (0.006944444444444444*(-185. + 838.*nu - 910.*nu2 + 220.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 1.0973940686333457*x2 + x3*(11.623366217471297 - 2.5002775002775004*el5));
  /** (6,6) */
  rholm[19] = 1. + (0.011904761904761904*(-106. + 602.*nu - 861.*nu2 + 273.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + PMTERMS_eps*(- 1.5543111183867486*x2 + x3*(16.645950799433503 - 3.6003996003996006*el6));
  
  /** l=7 ------------------------------------------------------------------
   *  (7,1) */
  rholm[20] = 1. + (0.0014005602240896359*(-618. + 2518.*nu - 2083.*nu2 + 228.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*( - 0.1508235111143767*x2 + x3*(0.2581280702019663 - 0.07355557607658449*el1));
  /** (7,2) */
  rholm[21] = 1. + (0.00006669334400426837*(16832. - 123489.*nu + 273924.*nu2 - 190239.*nu3 + 32760.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 0.351319484450667*PMTERMS_eps*x2;
  /** (7,3) */
  rholm[22] = 1. + (0.0014005602240896359*(-666. + 2806.*nu - 2563.*nu2 + 420.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 0.37187416047628863*x2 + x3*(3.0835293524055283 - 0.6620001846892604*el3));
  /** (7,4) */
  rholm[23] = 1. + (0.00006669334400426837*(17756. - 131805.*nu + 298872.*nu2 - 217959.*nu3 + 41076.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 0.6473746896670599*PMTERMS_eps*x2;
  /** (7,5) */
  rholm[24] = 1. + (0.0014005602240896359*(-762. + 3382.*nu - 3523.*nu2 + 804.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 0.8269193364414116*x2 + x3*(8.750589067052443 - 1.838889401914612*el5));
  /** (7,6) */
  rholm[25] = 1. + (0.0006002400960384153*(2144. - 16185.*nu + 37828.*nu2 - 29351.*nu3 + 6104.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - 1.1403265020692532*PMTERMS_eps*x2;
  /** (7,7) */
  rholm[26] = 1. + (0.0014005602240896359*(-906. + 4246.*nu - 4963.*nu2 + 1380.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 1.5418467934923434*x2 + x3*(17.255875091408523 - 3.6042232277526396*el7));
  
  /** l=8 ------------------------------------------------------------------
   *  (8,1) */
  rholm[27] = 1. + (0.00005482456140350877*(20022. - 126451.*nu + 236922.*nu2 - 138430.*nu3 + 21640.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - 0.26842133517043704*PMTERMS_eps*x2;
  /** (8,2) */
  rholm[28] = 1. + (0.0003654970760233918*(2462. - 17598.*nu + 37119.*nu2 - 22845.*nu3 + 3063.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - 0.2261796441029474*PMTERMS_eps*x2;
  /** (8,3) */
  rholm[29] = 1. + (0.00005482456140350877*(20598. - 131059.*nu + 249018.*nu2 - 149950.*nu3 + 24520.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - 0.4196774909106648*PMTERMS_eps*x2;
  /** (8,4) */
  rholm[30] = 1. + (0.0003654970760233918*(2666. - 19434.*nu + 42627.*nu2 - 28965.*nu3 + 4899.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - 0.47652059150068155*PMTERMS_eps*x2;
  /** (8,5) */
  rholm[31] = 1. + (0.00027412280701754384*(4350. - 28055.*nu + 54642.*nu2 - 34598.*nu3 + 6056.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3)- 0.7220789990670207*PMTERMS_eps*x2;
  /** (8,6) */
  rholm[32] = 1. + (0.0010964912280701754*(1002. - 7498.*nu + 17269.*nu2 - 13055.*nu3 + 2653.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 0.9061610303170207*PMTERMS_eps*x2;
  /** (8,7) */
  rholm[33] = 1. + (0.00005482456140350877*(23478. - 154099.*nu + 309498.*nu2 - 207550.*nu3 + 38920.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3)- 1.175404252991305*PMTERMS_eps*x2;
  /** (8,8) */
  rholm[34] = 1. + (0.0003654970760233918*(3482. - 26778.*nu + 64659.*nu2 - 53445.*nu3 + 12243.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 1.5337092502821381*PMTERMS_eps*x2;
  
  /** Amplitudes */
  int k;
  for (k = 0; k < KMAX; k++) {
      flm[k] = gsl_pow_int(rholm[k], LINDEX[k]);
      //printf("flm %d %.16e\n",k,flm[k]);
  }
  //DBGSTOP
  
}

/** Resummed amplitudes for the spin case. 
    This function computes the residual amplitude corrections flm's as 
    introduced in Damour, Iyer & Nagar, PRD 79, 064004 (2008).
    The orbital part is taken at the usual 3^{+2} PN order, i.e. 3PN terms
    are integrated by the 4PN and 5PN test-particle terms, with the higher
    modes obtained by Fujita & Iyer.
    Note that the variables called here (a1,a2)
    are what we usually cal tilde{a}_1 and tilde{a}_2 and are defined as
    a1 = X1*chi1, a2=X2*chi2 and are passed here as parameters. Special 
    combinations of these quantities are used here to write the spin-dependent
    part of the waveform in particularly compact form, so that the (spinning)
    test-particle limit is recovered just by visual inspection of the equations 
    Old version, repeats the rholm computation for the orbital part.
*/
void eob_wav_flm_s_old(double x, double nu, double X1, double X2, double chi1, double chi2, double a1, double a2, double C_Q1, double C_Q2, int usetidal,
		       double *rholm, double *flm)
{

  /* declaring the spin-dependent terms */
  double rho22S;
  double rho32S;
  double rho44S;
  double rho42S;
  double f21S;
  double f33S;
  double f31S;
  double f43S;
  double f41S;
      
  const double a0      = a1+a2;
  const double a12     = a1-a2;
  const double X12     = X1-X2;
  const double a0X12   = a0*X12;
  const double a12X12  = a12*X12;
  
  /* Shorthands */
  const double nu2 = nu*nu;
  const double nu3 = nu*nu2;
  const double nu4 = nu*nu3;
  const double x2 = x*x;
  const double x3 = x*x2;
  const double x4 = x*x3;
  const double x5 = x*x4;
  const double v  = sqrt(x);
  const double v2 = x;
  const double v3 = v*v2;
  const double v4 = v3*v;
  const double v5 = v4*v;
  
  /* Compute EulerLogs */
  const double el1 = Eulerlog(x,1);
  const double el2 = Eulerlog(x,2);
  const double el3 = Eulerlog(x,3);
  const double el4 = Eulerlog(x,4);
  const double el5 = Eulerlog(x,5);
  const double el6 = Eulerlog(x,6);
  const double el7 = Eulerlog(x,7);
    
  /** l=m=2 multipole */
  /* spin-orbit */
  const double cSO_lo    = (-0.5*a0 - a12X12/6.);
  const double cSO_nlo   = (-52./63.-19./504.*nu)*a0 - (50./63.+209./504.*nu)*a12X12;
  
  /* SPIN-SPIN contribution */
  double cSS_lo;
  if (usetidal) {
#if (EXCLUDESPINSPINTIDES)
    /* Switch off spin-spin-tidal couplings */
    /* See also: eob_dyn_s_get_rc() */
    cSS_lo = 0.;
    /* Above code switch off everything, 
       Alt. one can set C_Q1=C_Q2=0, but keep the term: */
    /*
      cSS_lo = a1*a2;
    */
#else
    cSS_lo = 0.5*(C_Q1*a1*a1 + 2.*a1*a2 + C_Q2*a2*a2);
#endif
  } else {
    cSS_lo = 0.5*a0*a0; 
  }

  /* rho_22^S: Eq. (80) of Damour & Nagar, PRD 90, 044018 (2014) */
  rho22S = cSO_lo*v3 + cSS_lo*v4 + cSO_nlo*v5 ;
    
  /** l>=3, m=even: multipoles rewritten in compact and self-explanatory form */
  rho32S = (a0-a12X12)/(3.*(1.-3.*nu))*v;
  rho44S = (-19./30.*a0 -  (1.-21.*nu)/(30.-90.*nu)*a12X12)*v3;
  rho42S = ( -1./30.*a0 - (19.-39.*nu)/(30.-90.*nu)*a12X12)*v3;
  
  /** l>=2, m=odd: multipoles rewritten in compact and self-explanatory form */
  f21S = -1.5*a12*v + ((110./21. + 79./84.*nu)*a12 - 13./84.*a0X12)*v3;
  f33S = ((-0.25 + 2.5*nu)*a12 - 1.75*a0X12)*v3;
  f31S = ((-2.25 + 6.5*nu)*a12 + 0.25*a0X12)*v3;
  f43S = (( 5. -10.*nu)*a12 - 5.*a0X12)/(-4.+8.*nu)*v;
  f41S = f43S;
 
   
  // l=2 ------------------------------------------------------------------
  //(2,2)
  rholm[1] = 1. + (-1.0238095238095237 + 0.6547619047619048*nu)*x + (-1.94208238851096 - 1.5601379440665155*nu + 0.4625614134542706*nu2)*x2 + x3*(12.736034731834051 - 2.902228713904598*nu - 1.9301558466099282*nu2 + 0.2715020968103451*nu3 - 4.076190476190476*el2) + x4*(-2.4172313935587004 + 4.173242630385488*el2) + x5*(-30.14143102836864 + 7.916297736025627*el2);
    
    //(2,1)
    rholm[0] = 1. + (-1.0535714285714286 + 0.27380952380952384*nu)*x + (-0.8327841553287982 - 0.7789824263038548*nu + 0.13116496598639457*nu2)*x2 + x3*(2.9192806270460925 - 1.019047619047619*el1) + x4*(-1.28235780892213 + 1.073639455782313*el1) + PMTERMS_eps*x5*(-3.8466571723355227 + 0.8486467106683944*el1);
    
     
    // l=3 ------------------------------------------------------------------
    //(3,3)
    rholm[4] = 1. + (-1.1666666666666667 + 0.6666666666666666*nu)*x + (-1.6967171717171716 - 1.8797979797979798*nu + 0.45151515151515154*nu2)*x2 + x3*(14.10891386831863 - 3.7142857142857144*el3) + x4*(-6.723375314944128 + 4.333333333333333*el3) + PMTERMS_eps*x5*(-29.568699895427518 + 6.302092352092352*el3);
    
    //(3,2)
    rholm[3] = 1. + (0.003703703703703704*(328. - 1115.*nu + 320.*nu2)*x)/(-1. + 3.*nu) + (6.235191420376606e-7*(-1.444528e6 + 8.050045e6*nu - 4.725605e6*nu2 - 2.033896e7*nu3 + 3.08564e6*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + x3*(6.220997955214429 - 1.6507936507936507*el2) + PMTERMS_eps*x4*(-3.4527288879001268 + 2.005408583186361*el2);
    

    
    //(3,1)
    rholm[2] = 1. + (-0.7222222222222222 - 0.2222222222222222*nu)*x + (0.014169472502805836 - 0.9455667789001122*nu - 0.46520763187429853*nu2)*x2 + x3*(1.9098284139598072 - 0.4126984126984127*el1) + x4*(0.5368150316615179 + 0.2980599647266314*el1) + PMTERMS_eps*x5*(1.4497991763035063 - 0.0058477188106817735*el1);
    
    // l=4 ------------------------------------------------------------------
    //(4,4)
    rholm[8] = 1. + (0.0007575757575757576*(1614. - 5870.*nu + 2625.*nu2)*x)/(-1. + 3.*nu) + (3.1534122443213353e-9*(-5.11573572e8 + 2.338945704e9*nu - 3.13857376e8*nu2 - 6.733146e9*nu3 + 1.252563795e9*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + x3*(15.108111214795123 - 3.627128427128427*el4) + PMTERMS_eps*x4*(-8.857121657199649 + 4.434988849534304*el4);
    
    
    //(4,3)
    rholm[7] = 1. + (0.005681818181818182*(222. - 547.*nu + 160.*nu2)*x)/(-1. + 2.*nu) - 0.9783218202252293*x2 + PMTERMS_eps*(x3*(8.519456157072423 - 2.0402597402597404*el3) +      x4*(-5.353216984886716 + 2.5735094451003544*el3));
    
    //(4,2)
    rholm[6] = 1. + (0.0007575757575757576*(1146. - 3530.*nu + 285.*nu2)*x)/(-1. + 3.*nu) - (3.1534122443213353e-9*(1.14859044e8 - 2.95834536e8*nu - 1.204388696e9*nu2 + 3.04798116e9*nu3 + 3.79526805e8*nu4)*x2)/((-1. + 3.*nu)*(-1. + 3.*nu)) + 4.550378418934105e-12*x3*(8.48238724511e11 - 1.9927619712e11*el2) + PMTERMS_eps*x4*(-0.6621921297263365 + 0.787251738160829*el2);
    
    
    //(4,1)
    rholm[5] = 1. + (0.001893939393939394*(602. - 1385.*nu + 288.*nu2)*x)/(-1. + 2.*nu) - 0.36778992787515513*x2 + x3*(0.6981550175535535 - 0.2266955266955267*el1) + PMTERMS_eps*x4*(-0.7931524512893319 + 0.2584672482399755*el1);
    
    // l=5 ------------------------------------------------------------------
    //(5,5)
    rholm[13] = 1. + (0.002564102564102564*(487. - 1298.*nu + 512.*nu2)*x)/(-1. + 2.*nu) - 1.5749727622804546*x2 + PMTERMS_eps*(x3*(15.939827047208668 - 3.6037296037296036*el5) +      x4*(-10.272578060123237 + 4.500041838503377*el5));
    //(5,4)
    rholm[12] = 1. + (0.00007326007326007326*(-17448. + 96019.*nu - 127610.*nu2 + 33320.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + PMTERMS_eps*(- 1.0442142414362194*x2 +   x3*(10.252052781721588 - 2.3063869463869464*el4));
    //(5,3)
    rholm[11] = 1. + (0.002564102564102564*(375. - 850.*nu + 176.*nu2)*x)/(-1. + 2.*nu) - 0.5788010707241477*x2 + PMTERMS_eps*(x3*(5.733973288504755 - 1.2973426573426574*el3) +      x4*(-1.9573287625526001 + 1.2474448628294783*el3));
    //(5,2)
    rholm[10] = 1. + (0.00007326007326007326*(-15828. + 84679.*nu - 104930.*nu2 + 21980.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + PMTERMS_eps*(- 0.4629337197600934*x2 +      x3*(2.354458371550237 - 0.5765967365967366*el2));
    //(5,1)
    rholm[9] = 1. + (0.002564102564102564*(319. - 626.*nu + 8.*nu2)*x)/(-1. + 2.*nu) - 0.1047896120973044*x2 + PMTERMS_eps*(x3*(0.642701885362399 - 0.14414918414918415*el1) +      x4*(-0.07651588046467575 + 0.11790664036817883*el1));
    
    // l=6 ------------------------------------------------------------------
    //(6,6)
    rholm[19] = 1. + (0.011904761904761904*(-106. + 602.*nu - 861.*nu2 + 273.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + PMTERMS_eps*(- 1.5543111183867486*x2 + x3*(16.645950799433503 - 3.6003996003996006*el6));
    //(6,5)
    rholm[18] = 1. + (0.006944444444444444*(-185. + 838.*nu - 910.*nu2 + 220.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 1.0973940686333457*x2 + x3*(11.623366217471297 - 2.5002775002775004*el5));
    //(6,4)
    rholm[17] = 1. + (0.011904761904761904*(-86. + 462.*nu - 581.*nu2 + 133.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + PMTERMS_eps*(- 0.7228451986855349*x2 + x3*(7.359388663371044 - 1.6001776001776002*el4));
    //(6,3)
    rholm[16] = 1. + (0.006944444444444444*(-169. + 742.*nu - 750.*nu2 + 156.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 0.5605554442947213*x2 + x3*(4.002558222882566 - 0.9000999000999002*el3));
    //(6,2)
    rholm[15] = 1. + (0.011904761904761904*(-74. + 378.*nu - 413.*nu2 + 49.*nu3)*x)/(1. - 5.*nu + 5.*nu2) + PMTERMS_eps*( - 0.24797525070634313*x2 + x3*(1.7942694138754138 - 0.40004440004440006*el2));
    //(6,1)
    rholm[14] = 1. + (0.006944444444444444*(-161. + 694.*nu - 670.*nu2 + 124.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 0.29175486850885135*x2 + x3*(0.21653486654395454 - 0.10001110001110002*el1));
    
    // l=7 ------------------------------------------------------------------
    //(7,7)
    rholm[26] = 1. + (0.0014005602240896359*(-906. + 4246.*nu - 4963.*nu2 + 1380.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 1.5418467934923434*x2 + x3*(17.255875091408523 - 3.6042232277526396*el7));
    //(7,6)
    rholm[25] = 1. + (0.0006002400960384153*(2144. - 16185.*nu + 37828.*nu2 - 29351.*nu3 + 6104.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - 1.1403265020692532*PMTERMS_eps*x2;
    //(7,5)
    rholm[24] = 1. + (0.0014005602240896359*(-762. + 3382.*nu - 3523.*nu2 + 804.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 0.8269193364414116*x2 + x3*(8.750589067052443 - 1.838889401914612*el5));
    //(7,4)
    rholm[23] = 1. + (0.00006669334400426837*(17756. - 131805.*nu + 298872.*nu2 - 217959.*nu3 + 41076.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 0.6473746896670599*PMTERMS_eps*x2;
    //(7,3)
    rholm[22] = 1. + (0.0014005602240896359*(-666. + 2806.*nu - 2563.*nu2 + 420.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*(- 0.37187416047628863*x2 + x3*(3.0835293524055283 - 0.6620001846892604*el3));
    //(7,2)
    rholm[21] = 1. + (0.00006669334400426837*(16832. - 123489.*nu + 273924.*nu2 - 190239.*nu3 + 32760.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 0.351319484450667*PMTERMS_eps*x2;
    //(7,1)
    rholm[20] = 1. + (0.0014005602240896359*(-618. + 2518.*nu - 2083.*nu2 + 228.*nu3)*x)/(1. - 4.*nu + 3.*nu2) + PMTERMS_eps*( - 0.1508235111143767*x2 + x3*(0.2581280702019663 - 0.07355557607658449*el1));
    
        
    // l=8 ------------------------------------------------------------------
    //(8,8)
    rholm[34] = 1. + (0.0003654970760233918*(3482. - 26778.*nu + 64659.*nu2 - 53445.*nu3 + 12243.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 1.5337092502821381*PMTERMS_eps*x2;
    //(8,7)
    rholm[33] = 1. + (0.00005482456140350877*(23478. - 154099.*nu + 309498.*nu2 - 207550.*nu3 + 38920.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3)- 1.175404252991305*PMTERMS_eps*x2;
    //(8,6)
    rholm[32] = 1. + (0.0010964912280701754*(1002. - 7498.*nu + 17269.*nu2 - 13055.*nu3 + 2653.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3)- 0.9061610303170207*PMTERMS_eps*x2;
    //(8,5)
    rholm[31] = 1. + (0.00027412280701754384*(4350. - 28055.*nu + 54642.*nu2 - 34598.*nu3 + 6056.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3)- 0.7220789990670207*PMTERMS_eps*x2;
    //(8,4)
    rholm[30] = 1. + (0.0003654970760233918*(2666. - 19434.*nu + 42627.*nu2 - 28965.*nu3 + 4899.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - 0.47652059150068155*PMTERMS_eps*x2;
    //(8,3)
    rholm[29] = 1. + (0.00005482456140350877*(20598. - 131059.*nu + 249018.*nu2 - 149950.*nu3 + 24520.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - 0.4196774909106648*PMTERMS_eps*x2;
    //(8,2)
    rholm[28] = 1. + (0.0003654970760233918*(2462. - 17598.*nu + 37119.*nu2 - 22845.*nu3 + 3063.*nu4)*x)/(-1. + 7.*nu - 14.*nu2 + 7.*nu3) - 0.2261796441029474*PMTERMS_eps*x2;
    //(8,1)
    rholm[27] = 1. + (0.00005482456140350877*(20022. - 126451.*nu + 236922.*nu2 - 138430.*nu3 + 21640.*nu4)*x)/(-1. + 6.*nu - 10.*nu2 + 4.*nu3) - 0.26842133517043704*PMTERMS_eps*x2;
    
    /** Amplitudes */
    flm[0] = gsl_pow_int(rholm[0], 2);
    flm[0] = (X12*flm[0] + f21S);
    
    flm[1] = gsl_pow_int(rholm[1]+ rho22S, 2);

    flm[2] = gsl_pow_int(rholm[2], 3);
    flm[2] = (X12*flm[2] + f31S);

    flm[3] = gsl_pow_int(rholm[3]+ rho32S, 3);

    flm[4] = gsl_pow_int(rholm[4], 3);
    flm[4] = (X12*flm[4] + f33S);
    
    flm[5] = gsl_pow_int(rholm[5], 4);
    flm[5] = (X12*flm[5] + f41S);
    
    flm[6]=gsl_pow_int(rholm[6] + rho42S, 4);

    flm[7] = gsl_pow_int(rholm[7], 4);
    flm[7] = (X12*flm[7] + f43S);

    flm[8] = gsl_pow_int(rholm[8] + rho44S, 4);

    /** Amplitudes */
    int k;
    for (k = 9; k < KMAX; k++) {
      flm[k] = gsl_pow_int(rholm[k], LINDEX[k]);
    }
    
}

