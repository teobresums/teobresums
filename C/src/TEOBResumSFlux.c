/** \file TEOBResumSFlux.c      
 *  \brief Compute the flux for the TEOBResumS model
 * 
 *  This file contains the routines required for the computation of the radiation reaction
 *  for the TEOBResumS rhs
 */

#include "TEOBResumS.h"

/** Coefficients for Newtonian flux */
/*
static const double CNlm[35] = {
  8./45, 32./5,
  1./1260, 32./63, 243./28, 
  1./44100, 32./3969, 729./700, 8192./567, 
  1./19958400, 256./400950, 2187./70400, 131072./66825, 1953125./76032, 
  1./1123782660, 128./28667925, 59049./15415400, 4194304./47779875, 48828125./13621608, 839808./17875,
  1./9.3498717312e11, 32./1.35270135e8, 1594323./3.2064032e10, 4194304./3.07432125e8, 1220703125./5.666588928e9, 5668704./875875, 96889010407./1111968000, 
  1./8.174459284992e13, 32./3.4493884425e10, 177147./3.96428032e10, 4194304./1.5679038375e10, 30517578125./8.00296713216e11, 51018336./1.04229125e8, 4747561509943./4.083146496e11, 274877906944./1688511825.
};
*/
static const double CNlm[54] = {
  0.17777777777777778, 6.4, 
  0.0007936507936507937, 0.5079365079365079, 8.678571428571429, 
  2.2675736961451248e-05, 0.008062484252960444, 1.0414285714285714, 14.447971781305114, 
  5.010421677088344e-08, 0.0006384836014465644, 0.03106534090909091, 1.9614216236438458, 25.688197074915823, 
  8.898517797026696e-10, 4.464920289836115e-06, 0.003830520129221428, 0.08778390483440988, 3.584607999290539, 46.98226573426573, 
  1.0695333890657086e-12, 2.3656367312710969e-07, 4.972309783123969e-05, 0.013643024456211269, 0.21542115380351798, 6.472046810332524, 87.1329124642076, 
  1.2233225038333268e-14, 9.27700678929842e-10, 4.468579053461083e-06, 0.0002675102834551229, 0.03813282951314997, 0.4894825318738884, 11.627213264559293, 162.79300083906728,
  1./109728845674905600.,128./3805293886875., 177147./5321402777600., 67108864./2049004400625., 762939453125./751987893436416., 204073344./2186104375., 1628413597910449./1534669170278400., 17592186044416./843707694375., 205891132094649/674204876800.,
  7.176081169626815e-20, 128./1518312260863125.,1594323./735210636160000.,2147483648./5722869290945625.,19073486328125./126059803844502528.,918330048./290751881875., 11398895185373143./54047249118720000., 2251799813685248./1009918110166875.,16677181699666569./446134008320000., 48828125000000./85045735593.
};

/**
 * Function: eob_flx_FlmNewt
 * -------------------------
 *   Compute the Newtonian partial fluxes
 * 
 *   @param[in] x       :  frequency parameter (v_phi^2)
 *   @param[in] nu     :  symmetric mass ratio
 *   @param[in,out] Nlm:  Newtonian partial fluxes
 */void eob_flx_FlmNewt(double x, double nu, double *Nlm)
{
  
  /** Shorthands*/
  const double nu2 = nu*nu;
  const double nu3 = nu2*nu;
  const double nu4 = nu3*nu;

  const double x5  = x*x*x*x*x;
  const double x6  = x*x5;
  const double x7  = x*x6;
  const double x8  = x*x7;
  const double x9  = x*x8;
  const double x10 = x*x9;
  const double x11 = x*x10;
  const double x12 = x*x11;
  const double x13 = x*x12;
  const double x14 = x*x13;

  const double sp2 = 1.-4.*nu;
  const double sp4 = (1.-4.*nu)*SQ((1.-2.*nu));
  const double sp3 = (1.-3.*nu)*(1.-3.*nu);
  const double sp5 = (1.-5.*nu+5.*nu2)*(1.-5.*nu+5.*nu2);
  const double sp6 = (1.-4.*nu)*(3.*nu2-4.*nu +1.)*(3.*nu2-4.*nu +1.);
  const double sp7 = (1. - 7.*nu + 14.*nu2 - 7.*nu3)*(1. - 7.*nu + 14.*nu2 - 7.*nu3);
  const double sp8 = (1. - 4.*nu)*(1. - 6.*nu + 10.*nu2 - 4.*nu3)*(1. - 6.*nu + 10.*nu2 - 4.*nu3);

  const double sp9  = (1. -4.*nu)*SQ(1. - 6.*nu + 10.*nu2 - 4.*nu3);
  const double sp10 = SQ(1. -9.*nu + 9.*nu4 + 27.*nu2 -30.*nu3);
  const double sp11 = (1. -4.*nu)*SQ(1. -8.*nu + 5.*nu4 +21.*nu2 -20.*nu3);
  
  double spx[] = {
    sp2 * x6, x5, 
    sp2 * x6, sp3 * x7, sp2 * x6, 
    sp4 * x8, sp3 * x7, sp4 * x8, sp3 * x7, 
    sp4 * x8, sp5 * x9, sp4 * x8, sp5 * x9, sp4 * x8, 
    sp6 * x10, sp5 * x9, sp6 * x10, sp5 * x9, sp6 * x10, sp5 * x9, 
    sp6 * x10, sp7 * x11, sp6 * x10, sp7 * x11, sp6 * x10, sp7 * x11, sp6 * x10,
    sp8 * x12, sp7 * x11, sp8 * x12, sp7 * x11, sp8 * x12, sp7 * x11, sp8 * x12, sp7 * x11, //(7.*nu3-14.*nu2+7.*nu-1.)*(7.*nu3-14.*nu2+7.*nu-1.) * x11,
    sp9 * x12, sp10* x13, sp9 * x12, sp10* x13, sp9 * x12, sp10* x13, sp9 * x12, sp10 * x13 , sp9 * x12,
    sp11* x14, sp10* x13, sp11* x14, sp10* x13, sp11* x14, sp10* x13, sp11* x14, sp10 * x13, sp11 * x14, sp10 * x13
  };

  /* Newtonian partial fluxes*/
  for (int k = 0; k < KMAX; k++) {
    Nlm[k] = CNlm[k] * spx[k];
  }
  
}

/**
  * Function: eob_flx_Tlm
  * ---------------------
  *   Compute the modulus of the tail term
  * 
  *   @param[in] w    :  E*Omega
  *   @param[out] MTlm:  modulus of the tail term
  */
void eob_flx_Tlm(const double w, double *MTlm)
{
  double hhatk, x2, y, prod, fl;
  for (int k = 0; k < KMAX; k++) {
    hhatk = MINDEX[k] * w;
    x2    = 4.*hhatk*hhatk;
    prod  = 1.;
    for (int j=1; j <= LINDEX[k]; j++) {
      prod *= ( j*j + x2 );
    }
    y  = 4.*Pi*hhatk;
    y /= ( 1. - exp(-y) ); 
    MTlm[k] = sqrt( 1./( SQ(fact(LINDEX[k])) ) * y * prod );
  }
}

/**
  * Function: eob_flx_HorizonFlux
  * -----------------------------
  *   Compute the horizon-absorbed fluxes
  *   Nagar & Akcay, PRD 85, 044025 (2012)
  *   Bernuzzi, Nagar & Zenginoglu, PRD 86, 104038 (2012)
  * 
  *   @param[in] x   :  frequency parameter, v_phi^2
  *   @param[in] Heff:  effective Hamiltonian
  *   @param[in] jhat:  angular momentum
  *   @param[in] nu  :  symmetric mass ratio
  *   
  *   @return hatFH:  horizon-absorbed fluxes
  */
double eob_flx_HorizonFlux_v1(double x, double Heff, double jhat, double nu)
{
  double rhoHlm[2]; /* only 21,22 multipoles -> k=0,1 */
  double FlmHLO[2];
  double FlmH[2];
  
  /* Shorthands */
  double nu2 = nu*nu;
  double nu3 = nu*nu2;
  double x2  = x*x;
  double x3  = x*x2;
  double x4  = x*x3;
  double x5  = x*x4;
  double x9  = x4*x5;
  double x10 = x*x9;
    
  /* The Newtonian asymptotic contribution */
  const double FNewt22 = 32./5.*x5;
  
  /* Compute leading-order part (nu-dependent) */
  FlmHLO[1] = 32./5.*(1-4*nu+2*nu2)*x9;
  FlmHLO[0] = 32./5.*(1-4*nu+2*nu2)*x10;
    
  /* Compute rho_lm */
  double c1[2];
  double c2[2];
  double c3[2];
  double c4[2];
    
  c1[0] = 0.58121;
  c2[0] = 1.01059;
  c3[0] = 7.955729;
  c4[0] = 1.650228;
  
  c1[1] = (4.-21.*nu + 27.*nu2 - 8.*nu3)/(4.*(1.-4.*nu+2.*nu2));
  c2[1] =  4.78752;
  c3[1] = 26.760136;
  c4[1] = 43.861478;
    
  rhoHlm[1] = 1. + c1[1]*x + c2[1]*x2 + c3[1]*x3 + c4[1]*x4;
  rhoHlm[0] = 1. + c1[0]*x + c2[0]*x2 + c3[0]*x3 + c4[0]*x4;
    
  /* Compute horizon multipolar flux (only l=2) */
  const double Heff2 = Heff*Heff;
  const double jhat2 = jhat*jhat;
  
  FlmH[1] = FlmHLO[1] * Heff2 * gsl_pow_int(rhoHlm[1],4);
  FlmH[0] = FlmHLO[0] * jhat2 * gsl_pow_int(rhoHlm[0],4);
    
  /* Sum over multipoles and normalize to the 22 Newtonian multipole */
  double hatFH = (FlmH[0]+FlmH[1])/FNewt22;
  
  return hatFH;
}

/**
  * Function: eob_flx_HorizonFlux_s
  * -------------------------------
  *   Compute the horizon-absorbed fluxes, spin case
  * 
  *   @param[in] x   :  frequency parameter, v_phi^2
  *   @param[in] Heff:  effective Hamiltonian
  *   @param[in] jhat:  angular momentum
  *   @param[in] nu  :  symmetric mass ratio
  *   @param[in] X1  :  mass fraction of body 1
  *   @param[in] X2  :  mass fraction of body 2
  *   @param[in] chi1:  dimensionless spin 1
  *   @param[in] chi2:  dimensionless spin 2
  * 
  *   @return hatFH:  horizon-absorbed fluxes
  */
double eob_flx_HorizonFlux_s_v1(double x, double Heff, double jhat, double nu, double X1, double X2, double chi1, double chi2)
{
    
  double x2 = x*x;
  double x3 = x2*x;
  double x4 = x3*x;
  double x5 = x4*x;
  double v5 = sqrt(x5);
  
  double cv5[2];
  double cv8[2];
  
  /** Coefficients of the v^5 term (Alvi leading order) */
  cv5[0] = -1./4.*chi1*(1.+3.*chi1*chi1)*X1*X1*X1;
  cv5[1] = -1./4.*chi2*(1.+3.*chi2*chi2)*X2*X2*X2;
  
  /** Coefficients of the v^8=x^4 term */
  cv8[0] = 0.5*(1.+sqrt(1.-chi1*chi1))*(1.+3.*chi1*chi1)*X1*X1*X1*X1;
  cv8[1] = 0.5*(1.+sqrt(1.-chi2*chi2))*(1.+3.*chi2*chi2)*X2*X2*X2*X2;
  
  double FH22_S = (cv5[0]+cv5[1])*v5;
  double FH22   = (cv8[0]+cv8[1])*x4;
  double FH21   =  0.0;
  
  /** Newton-normalized horizon flux: use only l=2 fluxes */
  double hatFH  = FH22_S + FH22 + FH21;
    
  return hatFH;
}

/** Horizon flux for large-mass-ratio nonspinning binaries, from Albertini et. al. XXXX.YYYY : 
    * 10PN for 22 +  nu-dependence @1PN
    * 15PN for 21
    * 6PN for l = 3
    * expanded fit for 44, 42
    * 12PN for 43, 41 
    * note: right now this is used whenever usespins = 1, 
    *       see line 697 in TEOBREsumSPars.c and lines 697, 834 below. 
*/
double eob_flx_HorizonFlux_lmr(double x, double Heff, double jhat, double nu)
{
  double rhoHlm[9]; /* all l = 2, 3, 4 multipoles */
  double FlmHLO[9];
  double FlmH[9];

  /** Coefficients */
  static double clm[KMAX][31]; // up to v^30 = x^15
  for (int k=0; k<KMAX; k++) clm[k][0] = 1.;
  for (int k=0; k<KMAX; k++) for (int n=1; n<31; n++) clm[k][n] = 0.;
  
  /** Shorthands */
  double nu2 = nu*nu;
  double nu3 = nu*nu2;
  double x2  = x*x;
  double x3  = x*x2;
  double x4  = x*x3;
  double x5  = x*x4;
  double x6  = x*x5;
  double x7  = x*x6;
  double x8  = x*x7;
  double x9  = x4*x5;
  double x10 = x*x9;
  double x11 = x*x10;
  double x12 = x*x11;

  double v     = sqrt(x);
  double v5    = x2*v;
  double logv  = log(v);
  double logv2 = logv*logv;
  double logv3 = logv2*logv;
  double logv4 = logv3*logv;
  double logv5 = logv4*logv;
    
  /** all Flm's below are already normalized by the 22 Newtonian contribution 32/5 x^5 */
  
  /** Compute leading-order part */

  /* l = 2, m = 1, 2 */
  FlmHLO[0] = (1-4*nu+2*nu2)*x5;
  FlmHLO[1] = (1-4*nu+2*nu2)*x4;
  /* l = 3, m = 1,2,3 */
  FlmHLO[2] = x6/56;
  FlmHLO[3] = 20*x7/63;
  FlmHLO[4] = 15*x6/56;
  /* l = 4, m = 1,2,3,4 */
  FlmHLO[5] = 5*x9/392;
  FlmHLO[6] = 5*x8/441;
  FlmHLO[7] = 5*x9/56;
  FlmHLO[8] = 5*x8/63;
    
  /* coefficients */

  /* l = 2, m = 1, 15PN */

  clm[0][2] = 0.5833333333333333;
  clm[0][4] = 0.8864087301587302;
  clm[0][6] = 6.424901605857512 - 2.038095238095238*logv;
  clm[0][8] = 1.766752931390268 - 1.011111111111111*logv;
  clm[0][10] = 1.493961606813823 - 2.284892290249433*logv;
  clm[0][11] = 0.2845717789918374;
  clm[0][12] = -21.22122135102325 - 16.00071023614146*logv + 2.076916099773243*logv2;
  clm[0][13] = -0.7408111353861827;
  clm[0][14] = -36.25765093154207 - 1.658502033614163*logv + 0.4868783068783069*logv2;
  clm[0][15] = -0.1117709908353889;
  clm[0][16] = -74.76369738861186 - 12.86351094951647*logv + 3.774864035561314*logv2;
  clm[0][17] = 10.24798928556675 - 1.739953162978663*logv;
  clm[0][18] = 80.2948593563905 + 11.72133394983377*logv + 19.31505055832595*logv2 - 1.41098427095706*logv3;
  clm[0][19] = -21.55623759955979 + 4.377759326613204*logv;
  clm[0][20] = -45.88843267104524 + 87.39972915380464*logv - 8.93059977671359*logv2 + 0.776920466952213*logv3;
  clm[0][21] = -10.12322702932126 + 1.57537040488682*logv;
  clm[0][22] = -258.1760297892015 + 114.3521460900903*logv + 40.14180948157476*logv2 - 5.458835408086137*logv3;
  clm[0][23] = 138.6110729577363 - 66.5680325952109*logv + 5.31928538396334*logv2;
  clm[0][24] = -2420.325524199674 - 253.6114020959027*logv + 35.55678286318115*logv2 - 15.51782914742528*logv3 + 0.7189300809162163*logv4;
  clm[0][25] = -192.8106387836869 + 130.6238544076524*logv - 12.61012313876668*logv2;
  clm[0][26] = -2168.284353993513 + 190.0643165235938*logv - 177.0158212116286*logv2 + 32.36280845602777*logv3 - 2.089040601166981*logv4;
  clm[0][27] = -255.3688794224337 + 114.5643448447081*logv - 9.293481075122995*logv2;
  clm[0][28] = -4961.74781462676 + 690.4179508849885*logv + 26.12151060568234*logv2 - 81.62655350372286*logv3 + 7.068836828458563*logv4;
  clm[0][29] = 1042.935504712336 - 989.9896542581077*logv + 217.7114930804913*logv2 - 10.84121021112528*logv3;
  clm[0][30] = 18238.62622224579 + 4516.381214052836*logv + 601.4830290647311*logv2 - 103.8909274572888*logv3 + 10.30113915972032*logv4 - 0.2930495948877529*logv5;
  
  /* l = 2, m = 2, 10PN + nu-dependence */

  clm[1][2] = (4.-21.*nu + 27.*nu2 - 8.*nu3)/(4.*(1.-4.*nu+2.*nu2));
  clm[1][4] = 3.988095238095238;
  clm[1][6] = 24.0571791747528 - 8.152380952380952*logv;
  clm[1][8] = 13.8551086308806 - 25.21904761904762*logv;
  clm[1][10] = 199.4078081876445 + 19.50022675736961*logv;
  clm[1][11] = -54.63778156643277;
  clm[1][12] = -582.8750735738636 - 178.5327822715866*logv + 33.23065759637188*logv2;
  clm[1][13] = 280.8970416802098;
  clm[1][14] = -1450.143031531191 - 1495.235015255883*logv + 311.4985941043084*logv2;
  clm[1][15] = -605.2296742617768;
  clm[1][16] = 9097.201294081323 + 555.1423370953789*logv - 861.1631918799266*logv2;
  clm[1][17] = -6911.484312430366 + 1336.284029167613*logv;
  clm[1][18] = -12313.69173035653 + 12628.97695598026*logv + 1163.917137711289*logv2 - 90.30299334125184*logv3;
  clm[1][19] = 44262.15268434681 - 9667.39349272306*logv;
  clm[1][20] = -65287.57081431632 - 36282.34512006451*logv + 23510.50238485977*logv2 - 2547.894737357377*logv3;

  /* l = 3, m = 1, 6PN */

  clm[2][2] = 1.611111111111111;
  clm[2][4] = 2.936728395061728;
  clm[2][6] = 9.213370934931813 - 0.8253968253968254*logv;
  clm[2][8] = 14.94079355370311 - 1.329805996472663*logv;
  clm[2][10] = 26.60058633678188 - 2.426788163825201*logv;
  clm[2][12] = 31.93969752137537 - 8.018297216743535*logv + 0.3406399596875787*logv2;

  /* l = 3, m = 2, 6PN */

  clm[3][2] = 0.8333333333333333;
  clm[3][4] = 1.483333333333333;
  clm[3][6] = 17.55939758125838 - 3.301587301587302*logv;
  clm[3][8] = 11.64342254747321 - 2.751322751322751*logv;
  clm[3][10] = 19.44076190460366 - 4.355555555555556*logv;
  clm[3][12] = -236.1045311517952 - 67.71250784926367*logv + 5.45023935500126*logv2;

  /* l = 3, m = 3, 6PN */

  clm[4][2] = 1.166666666666667;
  clm[4][4] = 2.941666666666667;
  clm[4][6] = 40.34914499874845 - 7.428571428571429*logv;
  clm[4][8] = 37.52373649525303 - 8.666666666666667*logv;
  clm[4][10] = 96.14382900245816 - 40.36666666666667*logv;
  clm[4][12] = -951.7147490236113 - 226.4706695003768*logv + 27.59183673469388*logv2;

  /* l = 4, m = 1, 12PN */

  clm[5][2] = 0.905952380952381;
  clm[5][4] = 1.161867849412492;
  clm[5][6] = 4.612151514725207 - 0.4533910533910534*logv;
  clm[5][8] = 4.45556947830768 - 0.4107507043221329*logv;
  clm[5][10] = 5.834974681453181 - 0.5267804881463277*logv;
  clm[5][12] = -4.740283887597416 - 2.217714990413922*logv + 0.1027817236475245*logv2;
  clm[5][14] = -1.814306997067652 - 2.134979998756876*logv + 0.09311534725686447*logv2;
  clm[5][15] = 0.00001435495089480142;
  clm[5][16] = -0.9449359329927948 - 2.792242657110017*logv + 0.1194187802132584*logv2;
  clm[5][17] = -0.0001044819818753333;
  clm[5][18] = 125.6261721611818 + 1.482428667619591*logv + 0.5314289759267583*logv2 - 0.01553343798463309*logv3;
  clm[5][19] = 0.000289431462887703;
  clm[5][20] = 121.5062731879449 + 0.181893489347434*logv + 0.5101627503844811*logv2 - 0.01407255512655451*logv3;
  clm[5][21] = 0.00003121019888384675 - 0.0000325420315378543*logv;
  clm[5][22] = 160.9063446509957 - 0.4011059612661441*logv + 0.6658805688470822*logv2 - 0.01804780218518797*logv3;
  clm[5][23] = -0.002642735277528471 + 0.0002368559791142116*logv;
  clm[5][24] = -1231.945263533703 - 56.81545345841274*logv - 0.1766899297640938*logv2 - 0.08463609444026904*logv3 + 0.00176068045265935*logv4;

  /* l = 4, m = 2, fit */

  double c1_42 = 1.43458;
  double c2_42 = 2.43232; 
  double c3_42 = 21.927986;
  double c4_42 = 10.419841;

  /* l = 4, m = 3, 12PN */

  clm[7][2] = 0.9535714285714286;
  clm[7][4] = 1.824282815398887;
  clm[7][6] = 29.8155597518837 - 4.080519480519481*logv;
  clm[7][8] = 24.58827566143869 - 3.891066790352505*logv;
  clm[7][10] = 45.98465156318206 - 7.444021566212081*logv;
  clm[7][12] = -981.9795394240014 - 130.7295482586175*logv + 8.325319615449486*logv2;
  clm[7][14] = -1084.478016966415 - 119.5234645513736*logv + 7.938786919017902*logv2;
  clm[7][15] = 2.542936486161387;
  clm[7][16] = -2008.92091650133 - 182.489551074673*logv + 15.18773750716777*logv2;
  clm[7][17] = -20.10488518653115;
  clm[7][18] = 88175.52617148325 + 3763.412721296057*logv + 275.5071020190452*logv2 - 11.32387629079753*logv3;
  clm[7][19] = 50.98457538841343;
  clm[7][20] = 85964.76478144973 + 3085.616876555689*logv + 359.8075763923598*logv2 - 10.79812489158193*logv3;
  clm[7][21] = 569.8165299459041 - 51.88250934752649*logv;
  clm[7][22] = 174421.7810462703 + 9919.505601789107*logv + 167.0685138147112*logv2 - 20.65795292100481*logv3;
  clm[7][23] = -4644.864041137085 + 410.1918782862395*logv;
  clm[7][24] = -9.31707434047562e6 - 347717.7619463258*logv - 8166.42545545039*logv2 - 320.6274444799858*logv3 + 11.551824449898*logv4;

  /* l = 4, m = 4, fit */

  double c1_44 = 1.15290;
  double c2_44 = 4.59627;
  double c3_44 = 55.268737;
  double c4_44 = 13.255971;

  /** Compute rho_lm */
  
  rhoHlm[0] = Taylorseries(v,clm[0],30);
  rhoHlm[1] = Taylorseries(v,clm[1],20);
  rhoHlm[2] = Taylorseries(v,clm[2],12);
  rhoHlm[3] = Taylorseries(v,clm[3],12);
  rhoHlm[4] = Taylorseries(v,clm[4],12);
  rhoHlm[5] = Taylorseries(v,clm[5],24);
  rhoHlm[6] = 1. + c1_42*x + c2_42*x2 + c3_42*x3 + c4_42*x4;
  rhoHlm[7] = Taylorseries(v,clm[7],24);
  rhoHlm[8] = 1. + c1_44*x + c2_44*x2 + c3_44*x3 + c4_44*x4;
    
  /** Compute horizon multipolar flux (l = 2, 3, 4) */
  const double Heff2 = Heff*Heff;
  const double jhat2 = jhat*jhat;
  
  /* l = 2, m = 1,2 */
  FlmH[0] = FlmHLO[0] * jhat2 * gsl_pow_int(rhoHlm[0],4);
  FlmH[1] = FlmHLO[1] * Heff2 * gsl_pow_int(rhoHlm[1],4);
  /* l = 3, m = 1,2,3 */
  FlmH[2] = FlmHLO[2] * Heff2 * gsl_pow_int(rhoHlm[2],6);
  FlmH[3] = FlmHLO[3] * jhat2 * gsl_pow_int(rhoHlm[3],6);
  FlmH[4] = FlmHLO[4] * Heff2 * gsl_pow_int(rhoHlm[4],6);
  /* l = 4, m = 1,2,3,4 */
  FlmH[5] = FlmHLO[5] * jhat2 * gsl_pow_int(rhoHlm[5],8);
  FlmH[6] = FlmHLO[6] * Heff2 * gsl_pow_int(rhoHlm[6],8);
  FlmH[7] = FlmHLO[7] * jhat2 * gsl_pow_int(rhoHlm[7],8);
  FlmH[8] = FlmHLO[8] * Heff2 * gsl_pow_int(rhoHlm[8],8);
    
  /** Sum over multipoles */
  double hatFH = 0.;

  for(int i = 0; i < 9; i++)
  {
    hatFH = hatFH + FlmH[i];
  }
  
  return hatFH;
}

/** Horizon flux for large-mass-ratio spinning binaries
 * work in progress, so still only LO contribution for l = 2;
 * also this is *not* used anywhere in the code yet since usespins = 1 
 * for nonspinning configurations as well and this is not consistent with 
 * the nonspin function for lmr. Will be updated soon
 */
double eob_flx_HorizonFlux_s_lmr(double x, double Heff, double jhat, double nu, double X1, double X2, double chi1, double chi2)
{

  return eob_flx_HorizonFlux_lmr(x, Heff, jhat, nu);
  // note: the spin version is still under development, currently using same as nonspin
  
  /*
 
  // only l = 2 for now 
  double FlmHLO[2];
  double FlmH[2];
  
  // Shorthands 
  double nu2 = nu*nu;
  double nu3 = nu*nu2;
  double x2  = x*x;
  double x3  = x*x2;
  double x4  = x*x3;
  double x5  = x*x4;
  
  // double x6  = x*x5;
  // double x7  = x*x6;
  // double x8  = x*x7;
  // double x9  = x4*x5;
  // double x10 = x*x9;
  // double x11 = x*x10;
  // double x12 = x*x11;

  double v     = sqrt(x);
  double v5    = x2*v;
  // double logv  = log(v);
  // double logv2 = logv*logv;
  // double logv3 = logv2*logv;
  // double logv4 = logv3*logv;
  // double logv5 = logv4*logv;
    
  // all Flm's below are already normalized by the 22 Newtonian contribution 32/5 x^5 
  
  // Compute leading-order part (nu-dependent) 

  // mass ratios X_i = m_i/M 
  const double X1to3 = X1*X1*X1;
  const double X1to4 = X1to3*X1;
  const double X2to3 = X2*X2*X2;
  const double X2to4 = X2to3*X2;

  // l = 2, m = 1 

  FlmHLO[0] = (1-4*nu+2*nu2)*x5;

  // l = 2, m = 2 

  // squared spins
  double chi1s = chi1*chi1;
  double chi2s = chi2*chi2;

  // coefficients of the v^5 term (Alvi leading order)
  double cv5[2];
  cv5[0] = -0.25*chi1*(1. + 3.*chi1s)*X1to3;
  cv5[1] = -0.25*chi2*(1. + 3.*chi2s)*X2to3;
    
  // coefficients of the v^8 = x^4 term
  double cv8[2];
  cv8[0] = 0.5*(1 + sqrt(1. - chi1s))*(1. + 3.*chi1s)*X1to4;
  cv8[1] = 0.5*(1 + sqrt(1. - chi2s))*(1. + 3.*chi2s)*X2to4;
    
  double FH_LOl2m2S = (cv5[0] + cv5[1])*v5;
  double FH_LOl2m2  = (cv8[0] + cv8[1])*x4;

  FlmHLO[1] = FH_LOl2m2 + FH_LOl2m2S;
  
  // nonspin limit: FlmHLO[1] = (1-4*nu+2*nu2)*x4;

  // other multipoles 

  // FlmHLO[2] = x6/56;
  // FlmHLO[3] = 20*x7/63;
  // FlmHLO[4] = 15*x6/56;
  
  // FlmHLO[5] = 5*x9/392;
  // FlmHLO[6] = 5*x8/441;
  // FlmHLO[7] = 5*x9/56;
  // FlmHLO[8] = 5*x8/63;
  
  // Compute horizon multipolar flux (only l = 2 for now) 

  FlmH[0] = FlmHLO[0];
  FlmH[1] = FlmHLO[1];
    
  // const double Heff2 = Heff*Heff;
  // const double jhat2 = jhat*jhat;

  // FlmH[0] = FlmHLO[0] * jhat2 * gsl_pow_int(rhoHlm[0],4);
  // FlmH[1] = FlmHLO[1] * Heff2 * gsl_pow_int(rhoHlm[1],4);

  // FlmH[2] = FlmHLO[2] * Heff2 * gsl_pow_int(rhoHlm[2],6);
  // FlmH[3] = FlmHLO[3] * jhat2 * gsl_pow_int(rhoHlm[3],6);
  // FlmH[4] = FlmHLO[4] * Heff2 * gsl_pow_int(rhoHlm[4],6);

  // FlmH[5] = FlmHLO[5] * jhat2 * gsl_pow_int(rhoHlm[5],8);
  // FlmH[6] = FlmHLO[6] * Heff2 * gsl_pow_int(rhoHlm[6],8);
  // FlmH[7] = FlmHLO[7] * jhat2 * gsl_pow_int(rhoHlm[7],8);
  // FlmH[8] = FlmHLO[8] * Heff2 * gsl_pow_int(rhoHlm[8],8);
    
  // Sum over multipoles 
  double hatFH = 0.;

  for(int i = 0; i < 2; i++)
  {
    hatFH = hatFH + FlmH[i];
  }
  
  return hatFH;

  */

}


/**
 * Function: eob_flx_Flux
 * ----------------------
 *   Flux calculation for Newton-Normalized energy flux 
 *   Use the DIN resummation procedure. 
 *   Add non-QC and non-K corrections to (2,2) partial flux.
 * 
 *   @param[in] x      :   frequency parameter, v_phi^2
 *   @param[in] Omega  :  orbital frequency
 *   @param[in] r_omega:  r*psi^(1./3) (from generalized Kepler's law) 
 *   @param[in] E      :  energy
 *   @param[in] Heff   :  effective Hamiltonian
 *   @param[in] jhat   :  angular momentum
 *   @param[in] r      :  radial separation
 *   @param[in] pr_star:  (tortoise) radial momentum
 *   @param[in] ddotr  :  radial acceleration
 *   @param[in] dyn    :  dynamics structure
 * 
 *   @return[out] Fphi:  energy flux
 */
double eob_flx_Flux(double x, double Omega, double r_omega, double E, double Heff, double jhat, double r, double pr_star, double ddotr, Dynamics *dyn)
{
  return eob_flx_Flux_s(x, Omega, r_omega, E, Heff, jhat, r, pr_star, ddotr,dyn);
}

/**
 * Function: eob_flx_Flux_s
 * ------------------------
 *   Flux calculation for spinning systems
 * 
 *   @param[in] x      :  frequency parameter
 *   @param[in] Omega  :  orbital frequency
 *   @param[in] r_omega: r*psi^(1./3) (from generalized Kepler's law) 
 *   @param[in] E      :  energy 
 *   @param[in] Heff   :  effective Hamiltonian
 *   @param[in] jhat   :  angular momentum
 *   @param[in] r      :  radial separation
 *   @param[in] pr_star:  (tortoise) radial momentum
 *   @param[in] ddotr  :  radial acceleration
 *   @param[in] dyn    :  dynamics structure
 * 
 *   @return[out] Fphi:  energy flux
*/
double eob_flx_Flux_s(double x, double Omega, double r_omega, double E, double Heff, double jhat, double r, double pr_star, double ddotr, Dynamics *dyn)
{
  const double nu = EOBPars->nu;
  const double chi1 = EOBPars->chi1;
  const double chi2 = EOBPars->chi2;
  const double X1 = EOBPars->X1;
  const double X2 = EOBPars->X2;
  const double a1 = EOBPars->a1;
  const double a2 = EOBPars->a2;
  const double C_Q1 = EOBPars->C_Q1;
  const double C_Q2 = EOBPars->C_Q2;
  const double X12  = X1-X2; /* sqrt(1-4nu) */
  const double X12sq = SQ(X12); /* (1-4nu) */

  const int usetidal = EOBPars->use_tidal;
  const int usespins = EOBPars->use_spins;
  
  double prefact[] = {
    jhat, Heff,
    Heff, jhat, Heff,
    jhat, Heff, jhat, Heff,
    Heff, jhat, Heff, jhat, Heff,
    jhat, Heff, jhat, Heff, jhat, Heff,
    Heff, jhat, Heff, jhat, Heff, jhat, Heff,
    jhat, Heff, jhat, Heff, jhat, Heff, jhat, Heff,
    Heff, jhat, Heff, jhat, Heff, jhat, Heff, jhat, Heff,
    jhat, Heff, jhat, Heff, jhat, Heff, jhat, Heff, jhat, Heff};  
  
  double FNewt22, sum_k=0.; 
  double rholm[KMAX] = {0.}, flm[KMAX] = {0.}, FNewtlm[KMAX], MTlm[KMAX], hlmTidal[KMAX], hlmNQC[KMAX];
  double Modhhatlm[KMAX];  

  /** Newtonian flux */
  eob_flx_FlmNewt(x, nu, FNewtlm);

  /* Remove useless modes when not using the 22PN flux*/
  // TODO: add also elsewhere
  if(EOBPars->use_flm != USEFLM_22PN){
    for(int k=35; k<KMAX; k++) FNewtlm[k] = 0;
  }

  /* Correct amplitudes for specific multipoles and cases */
  if (usespins) { // usespins is now 1 also when spins are 0
      /* Correct (2,1), (3,1) and (3,3) ( sp2 = 1 ) */
      double x6 = gsl_pow_int(x, 6);
      FNewtlm[0] = CNlm[0] * x6; // (2,1) 
      FNewtlm[2] = CNlm[2] * x6; // (3,1) 
      FNewtlm[4] = CNlm[4] * x6; // (3,3) 
      /* Correct (4,1), (4,3)  ( sp4 = (1-2nu)^2 ) */
      double sp4x8 = SQ((1. - 2.*nu)) * gsl_pow_int(x, 8);
      FNewtlm[5]  = CNlm[5] * sp4x8; // (4,1) 
      FNewtlm[7]  = CNlm[7] * sp4x8; // (4,3) 
      
      /* Correcting (5,5) for Higher Modes */
      if (EOBPars->use_flm == USEFLM_HM || EOBPars->use_flm == USEFLM_HM_4PN22){
        FNewtlm[13] = CNlm[13] * sp4x8;
      }
    } else {
      if (usetidal) {
        // Correct (2,1), (3,1) and (3,3) ( sp2 = 1 ) 
        double x6 = gsl_pow_int(x, 6);
        FNewtlm[0] = CNlm[0] * x6; // (2,1) 
        FNewtlm[2] = CNlm[2] * x6; // (3,1) 
        FNewtlm[4] = CNlm[4] * x6; // (3,3) 
      }
    }

  /* Tail term */
  eob_flx_Tlm(E*Omega, MTlm);

  /* Amplitudes */
  if (usespins) {
    /* eob_wav_flm_s_old(x,nu, X1,X2,chi1,chi2,a1,a2,C_Q1,C_Q2, usetidal, rholm, flm); */
    eob_wav_flm_s(x,nu, X1,X2,chi1,chi2,a1,a2,C_Q1,C_Q2, usetidal, rholm, flm);
  } else {
    /* eob_wav_flm_old(x,nu, rholm, flm); */
    eob_wav_flm(x,nu, rholm, flm);
  }
  
  FNewt22 = FNewtlm[1];

  /* NQC correction to the modulus of the (l,m) waveform */  
  for (int k = 0; k < KMAX; k++) hlmNQC[k] = 1.; /* no NQC */
  
  if (!(EOBPars->nqc_coefs_flx == NQC_FLX_NONE)) {
    
    Waveform_lm_t hNQC;
    /* eob_wav_hlmNQC_nospin201602(nu,r,pr_star,Omega,ddotr, &hNQC); */ 
    eob_wav_hlmNQC(nu,r,pr_star,Omega,ddotr, NQC->flx, &hNQC);
    const int maxk = MIN(KMAX, NQC->hlm->maxk+1);
    /*
      for (int k = 0; k < maxk; k++) {
      if (NQC->hlm->activemode[k]) {
      hlmNQC[k] = hNQC.ampli[k]; 
      }
      }
    */
    /* Use only the 22:  */
    hlmNQC[1] = hNQC.ampli[1]; 
  } 
  
  /* Compute modulus of hhat_lm (with NQC) */  
  for (int k = 0; k < KMAX; k++) { 
    Modhhatlm[k] = prefact[k] * MTlm[k] * flm[k] * hlmNQC[k];
  }

  if (usetidal) {
    /* Tidal amplitudes */
    eob_wav_hlmTidal(x,dyn, hlmTidal);
    if (!(usespins)) {
      /* Correct normalization of (2,1) (3,1), (3,3) point-mass amplitudes */
      Modhhatlm[0] *= X12;
      Modhhatlm[2] *= X12;
      Modhhatlm[4] *= X12;
    }
    /* Add tidal amplitudes */
    for (int k = 0; k < KMAX; k++) { 
      Modhhatlm[k] += MTlm[k] * hlmTidal[k]; 
    }
  }

  /* Total multipolar flux */
  for (int k = KMAX; k--;) sum_k += SQ(Modhhatlm[k]) * FNewtlm[k];
  
  /* Normalize to the 22 Newtonian multipole */
  double hatf = sum_k/(FNewt22);
    
  /* Horizon flux */ 
  if (!(usetidal)) {
    double hatFH;
    if (usespins) {
      hatFH = eob_flx_HorizonFlux_s(x, Heff, jhat, nu, X1, X2, chi1, chi2);
    } else {
      hatFH = eob_flx_HorizonFlux(x,Heff,jhat,nu);
    }
    hatf += hatFH;
  }

  /* return Fphi */  
  return (-32./5. * nu * gsl_pow_int(r_omega,4) * gsl_pow_int(Omega,5) * hatf);  
}

/* Flux calculation for eccentric systems */

/**
 * Function: eob_flx_Flux_ecc
 * ------------------------
 *   Flux calculation for eccentric systems
 *   See https://arxiv.org/abs/2001.11736
 * 
 *   @param[in] x        :  frequency parameter
 *   @param[in] Omega    :  orbital frequency
 *   @param[in] r_omega  : r*psi^(1./3) (from generalized Kepler's law) 
 *   @param[in] E        :  energy 
 *   @param[in] Heff     :  effective Hamiltonian
 *   @param[in] jhat     :  angular momentum
 *   @param[in] r        :  radial separation
 *   @param[in] pphi     :  orbital angular momentum
 *   @param[in] pr_star  :  (tortoise) radial momentum
 *   @param[in] ddotr    :  radial acceleration
 *   @param[in,out] Fphi :  angular momentum flux
 *   @param[in,out] Fr   :  radial flux
 *   @param[in] dyn      :  dynamics structure
 * 
*/
void eob_flx_Flux_ecc(double x, double Omega, double r_omega, double E, double Heff, double jhat, double r, double pr_star, double pphi, double rdot, double ddotr, double *Fphi, double *Fr, Dynamics *dyn)
{
  const double nu = EOBPars -> nu;
  const double chi1 = EOBPars -> chi1;
  const double chi2 = EOBPars -> chi2;
  const double X1 = EOBPars -> X1;
  const double X2 = EOBPars -> X2;
  const double a1 = EOBPars -> a1;
  const double a2 = EOBPars -> a2;
  const double C_Q1 = EOBPars -> C_Q1;
  const double C_Q2 = EOBPars -> C_Q2;
  const double X12 = X1 - X2; /* sqrt (1 - 4 nu)*/
  const double X12sq = SQ (X12); /* (1 - 4 nu)*/

  const int usetidal = EOBPars -> use_tidal;
  const int usespins = EOBPars -> use_spins;
  
  double prefact[] = {
		jhat, Heff,
    Heff, jhat, Heff,
    jhat, Heff, jhat, Heff,
    Heff, jhat, Heff, jhat, Heff,
    jhat, Heff, jhat, Heff, jhat, Heff,
    Heff, jhat, Heff, jhat, Heff, jhat, Heff,
    jhat, Heff, jhat, Heff, jhat, Heff, jhat, Heff,
    Heff, jhat, Heff, jhat, Heff, jhat, Heff, jhat, Heff,
    jhat, Heff, jhat, Heff, jhat, Heff, jhat, Heff, jhat, Heff};

  double FNewt22, sum_k=0.;
  double rholm[KMAX] = {0.}, flm[KMAX] = {0.}, FNewtlm[KMAX], MTlm[KMAX], hlmTidal[KMAX], hlmNQC[KMAX];
  double Modhhatlm[KMAX];

  /** Newtonian flux */
  eob_flx_FlmNewt(x, nu, FNewtlm);

  /* Remove useless modes when not using the 22PN flux*/
  // TODO: add also elsewhere
  if(EOBPars->use_flm != USEFLM_22PN){
    for(int k=35; k<KMAX; k++) FNewtlm[k] = 0;
  }

  /* Correct amplitudes for specific multipoles and cases */
  if (usespins) { // usespins is now 1 also when spins are 0
      /* Correct (2,1), (3,1) and (3,3) ( sp2 = 1 ) */
      double x6 = gsl_pow_int(x, 6);
      FNewtlm[0] = CNlm[0] * x6; // (2,1) 
      FNewtlm[2] = CNlm[2] * x6; // (3,1) 
      FNewtlm[4] = CNlm[4] * x6; // (3,3) 
      /* Correct (4,1), (4,3)  ( sp4 = (1-2nu)^2 ) */
      double sp4x8 = SQ((1. - 2.*nu)) * gsl_pow_int(x, 8);
      FNewtlm[5]  = CNlm[5] * sp4x8; // (4,1) 
      FNewtlm[7]  = CNlm[7] * sp4x8; // (4,3) 
      
      /* Correcting (5,5) for Higher Modes */
      if (EOBPars->use_flm == USEFLM_HM || EOBPars->use_flm == USEFLM_HM_4PN22){
        FNewtlm[13] = CNlm[13] * sp4x8;
      }
    } else {
      if (usetidal) {
        /* Correct (2,1), (3,1) and (3,3) ( sp2 = 1 ) */
        double x6 = gsl_pow_int(x, 6);
        FNewtlm[0] = CNlm[0] * x6; /* (2,1) */
        FNewtlm[2] = CNlm[2] * x6; /* (3,1) */
        FNewtlm[4] = CNlm[4] * x6; /* (3,3) */
      }
    }

  /** Tail term */
  eob_flx_Tlm(E*Omega, MTlm);

  /* Amplitudes */
  if (usespins) {
    /* eob_wav_flm_s_old(x,nu, X1,X2,chi1,chi2,a1,a2,C_Q1,C_Q2, usetidal, rholm, flm); */
    eob_wav_flm_s(x,nu, X1,X2,chi1,chi2,a1,a2,C_Q1,C_Q2, usetidal, rholm, flm);
  } else {
    /* eob_wav_flm_old(x,nu, rholm, flm); */
    eob_wav_flm(x,nu, rholm, flm);
  }
  
  FNewt22 = FNewtlm[1];

  /* NQC correction to the modulus of the (l,m) waveform */
  for (int k = 0; k < KMAX; k++) hlmNQC[k] = 1.; /* no NQC */
  
  if (!(EOBPars->nqc_coefs_flx == NQC_FLX_NONE)) {
    
    Waveform_lm_t hNQC;
    /* eob_wav_hlmNQC_nospin201602(nu,r,pr_star,Omega,ddotr, &hNQC); */
    eob_wav_hlmNQC(nu,r,pr_star,Omega,ddotr, NQC->flx, &hNQC);
    const int maxk = MIN(KMAX, NQC->hlm->maxk+1);
    /*
      for (int k = 0; k < maxk; k++) {
      if (NQC->hlm->activemode[k]) {
	hlmNQC[k] = hNQC.ampli[k];
      }
    */
    /* Use only the 22:  */
    hlmNQC[1] = hNQC.ampli[1];
  }
    
  /* Compute modulus of hhat_lm (with NQC) */
  for (int k = 0; k < KMAX; k++) {
    Modhhatlm[k] = prefact[k] * MTlm[k] * flm[k] * hlmNQC[k];
  }
  
  if (usetidal) {
    /* Tidal amplitudes */
    eob_wav_hlmTidal(x,dyn, hlmTidal);
    if (!(usespins)) {
      /* Correct normalization of (2,1) (3,1), (3,3) point-mass amplitudes */
      Modhhatlm[0] *= X12;
      Modhhatlm[2] *= X12;
      Modhhatlm[4] *= X12;
    }
    /* Add tidal amplitudes */
    for (int k = 0; k < KMAX; k++) {
      Modhhatlm[k] += MTlm[k] * hlmTidal[k];
    }
  }

  /* Total multipolar flux */
  for (int k = KMAX; k--;) sum_k += SQ(Modhhatlm[k]) * FNewtlm[k];
  
  /* Normalize to the 22 Newtonian multipole */
  double hatf = sum_k/(FNewt22);

  /* Horizon flux */
  double hatFH = 0.;
  if (!(usetidal)) {
    if (usespins) {
      hatFH = eob_flx_HorizonFlux_s(x, Heff, jhat, nu, X1, X2, chi1, chi2);
    } else {
      hatFH = eob_flx_HorizonFlux(x,Heff,jhat,nu);
    }
    hatf += hatFH;
  }

  /* Compute circular Fphi */
  *Fphi = -32./5. * nu * gsl_pow_int(r_omega,4) * gsl_pow_int(Omega,5) * hatf;

  /* Compute eccentric Fr */
  *Fr = eob_flx_Fr(r, pr_star, pphi, dyn, *Fphi);
  
  /* Compute non-circular Fphi */
  double Fphi_NC[KMAX];
  double fphi_nc = eob_flx_Fphi_ecc(r, pr_star, pphi, Omega, rdot, *Fphi, *Fr, dyn);
  for (int k = 0; k < KMAX; k++) Fphi_NC[k] = 1.;
  Fphi_NC[1] = fphi_nc;
  // To recover old configuration used for arXiv:2001.11736, one should apply this to all multipoles -> for (int k = KMAX; k--;) Fphi_NC[k] = fphi_nc;
  
  /* Adding non-circular corrections and re-compute flux */
  sum_k = 0.;
  for (int k = KMAX; k--;) sum_k += SQ(Modhhatlm[k]) * FNewtlm[k] * Fphi_NC[k];

  hatf = sum_k/(FNewt22);

  hatf += hatFH;

  /* Compute non-circular Fphi */
  *Fphi = -32./5. * nu * gsl_pow_int(r_omega,4) * gsl_pow_int(Omega,5) * hatf;
  
  /* Re-compute Fr using the generic Fphi */
  *Fr = eob_flx_Fr(r, pr_star, pphi, dyn, *Fphi);
}

/**
  * Function: eob_flx_Fr
  * --------------------
  *   Radial flux calculation for eccentric systems
  *   See https://arxiv.org/abs/2001.11736
  *       https://arxiv.org/abs/1210.2834
  * 
  *   @param[in] r        :  radial separation
  *   @param[in] pr_star  :  (tortoise) radial momentum
  *   @param[in] pphi     :  orbital angular momentum
  *   @param[in] dyn      :  dynamics structure
  *   @param[in] Fphi     :  angular momentum flux
  * 
  *   @return[out] Fr     :  radial flux
  */
double eob_flx_Fr_ecc(double r, double prstar, double pphi, Dynamics *dyn, double Fphi)
{
  const double nu = EOBPars->nu;
  double nu2 = nu*nu;

  double u  = 1/r;
  double u2 = u*u;
  double u3 = u2*u;
  double u4 = u3*u;
  double u5 = u4*u;
  double pphi2 = pphi*pphi;
  double pphi4 = pphi2*pphi2;
  double pphi6 = pphi4*pphi2;
  double prstar2 = prstar*prstar;
  double prstar4 = prstar2*prstar2;
  
  double a1, a2, F0PN, F1PN, F2PN;
  
  F0PN = -8./15. + 56./5.*pphi2*u;

  F1PN = prstar2*(-1228./105. + 556./105.*nu) + u*(-1984./105. + 16./21.*nu)
    + prstar2*pphi2*u*(-124./105. - 436./105.*nu) + pphi2*u2*(-1696./35. - 1268./105.*nu)
    + pphi4*u3*(1252./105. - 2588./105.*nu);

  /*
    Old configuration used for arXiv:2001.11736
    F2PN = prstar4*(323./315. + 1061./315.*nu - 1273./315.*nu2)
    + pphi6*u5*(-3229./315. - 718./63.*nu + 3277./105.*nu2)
    + pphi4*u4*(-35209./315. + 1606./15.*nu + 25217./315.*nu2)
    + u2*(59554./2835. + 9686./105.*nu - 3548./315.*nu2)
    + prstar2*pphi2*u2*(-1774./21. + 10292./315.*nu - 8804./315.*nu2)
    + prstar2*u*(20666./315. + 17590./189.*nu - 218./189.*nu2)
    + prstar4*pphi2*u*(-461./315. - 983./315.*nu + 131./63.*nu2)
    + pphi2*u3*(-29438./315. + 58424./315.*nu - 1752./35.*nu2)
    + prstar2*u3*pphi4*(-628./105. - 1052./105.*nu + 194./7.*nu2);
  */
  F2PN = prstar4*(323./315. + 1061./315.*nu - 1273./315.*nu2)
    + pphi6*u5*(-3229./315. - 718./63.*nu + 3277./105.*nu2)
    + pphi4*u4*(-35209./315. + 1606./15.*nu + 25217./315.*nu2)
    + u2*(59554./2835. + 9854./105.*nu - 3548./315.*nu2)
    + prstar2*pphi2*u2*(-1774./21. + 10292./315.*nu - 8804./315.*nu2)
    + prstar2*u*(20666./315. + 17590./189.*nu - 218./189.*nu2)
    + prstar4*pphi2*u*(-461./315. - 983./315.*nu + 131./63.*nu2)
    + pphi2*u3*(-29438./315. + 9568./315.*nu - 1752./35.*nu2)
    + prstar2*u3*pphi4*(-628./105. - 1052./105.*nu + 194./7.*nu2);
  
  a1 = F1PN/F0PN;
  a2 = F2PN/F0PN;
  
  /* return Fr */
  return nu*u4*prstar*F0PN/(1 - a1 + (a1*a1 - a2));
}

/**
  * Function: eob_flx_Fr_ecc_BD
  * ---------------------------
  *   Radial flux calculation for eccentric systems
  *   Circular expression from Bini-Damour inverse-resumme
  *   See https://arxiv.org/abs/1210.2834
  *       https://arxiv.org/abs/2001.11736
  * 
  *   @param[in] r        :  radial separation
  *   @param[in] pr_star  :  (tortoise) radial momentum
  *   @param[in] pphi     :  orbital angular momentum
  *   @param[in] dyn      :  dynamics structure
  *   @param[in] Fphi     :  angular momentum flux
  * 
  *   @return[out] Fr     :  radial flux
*/
double eob_flx_Fr_ecc_BD(double r, double prstar, double pphi, Dynamics *dyn, double Fphi)
{
  const double nu = EOBPars->nu;
  double nu2 = nu*nu;

  double u  = 1/r;
  double u2 = u*u;
  double u4 = u2*u2;

  double c1 = -573./280. - 118./35.*nu;
  double c2 = (-237433. + 175311.*nu + 70794.*nu2)/15120.;

  double if2 = 1./(1. - c1*u + (c1*c1 - c2)*u2);
  
  /* return Fr */
  return nu*32./3.*u4*prstar*if2;
}


/** Radial flux calculation for eccentric systems 
    Eq. 6 of arXiv:2304.09662
*/

/**
  * Function: eob_flx_Fr_ecc_next
  * -----------------------------
  *   Radial flux calculation for eccentric systems
  *   See Eq. 6 of arXiv:2304.09662
  * 
  *   @param[in] r        :  radial separation
  *   @param[in] pr_star  :  (tortoise) radial momentum
  *   @param[in] pphi     :  orbital angular momentum
  *   @param[in] dyn      :  dynamics structure
  *   @param[in] Fphi     :  angular momentum flux
  * 
  *   @return[out] Fr     :  radial flux
*/
double eob_flx_Fr_ecc_next(double r, double prstar, double pphi, Dynamics *dyn, double Fphi)
{
  const double nu = EOBPars->nu;
  const double nu2 = nu*nu;

  const double u  = 1/r;
  const double u2 = u*u;

  double c1 = 5317./1680 - 227./140*nu;
  double c2 = 1296935./1016064 - 274793./70560*nu + 753./560*nu2;
  double a[]  = {1., c1, c2};
  double hatf_prstar = Pade02(u, a);
  const double Frstar = -5./3.*prstar/pphi*Fphi*hatf_prstar;
  
  /* Eq. 6 of the paper above gives Fr*. 
     Compute the conversion from Fr to Fr*
  */
  double A, B, pl_hold;
  eob_metric_s(r, prstar, dyn, &A, &B, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold, &pl_hold);
  const double sqrtAbyB = sqrt(A/B);
  
  /* return Fr */
  return Frstar/sqrtAbyB;
}

/** Non-circular flux for eccentric systems */

/**
  * Function: eob_flx_Fphi_ecc
  * --------------------------
  *   Non-circular flux calculation for eccentric systems
  *   obtained via an iterative procedure (two iterations)
  *   See https://arxiv.org/abs/2001.11736
  * 
  *   @param[in] r        :  radial separation
  *   @param[in] pr_star  :  (tortoise) radial momentum
  *   @param[in] pphi     :  orbital angular momentum
  *   @param[in] Omg      :  orbital frequency
  *   @param[in] rdot     :  radial velocity
  *   @param[in] Fphi     :  angular momentum flux
  *   @param[in] Fr       :  radial flux
  *   @param[in] dyn      :  dynamics structure
  * 
  *   @return[out] Fphi   :  angular momentum flux
*/
double eob_flx_Fphi_ecc(double r, double prstar, double pphi, double Omg, double rdot, double Fphi, double Fr, Dynamics *dyn)
{  
  const double nu     = EOBPars -> nu;
  const double chi1   = EOBPars -> chi1;
  const double chi2   = EOBPars -> chi2;
  const double aK2    = EOBPars -> aK2;
  const double X1     = EOBPars -> X1;
  const double X2     = EOBPars -> X2;
  const double a1     = EOBPars -> a1;
  const double a2     = EOBPars -> a2;
  const double S1     = EOBPars -> S1;
  const double S2     = EOBPars -> S2;
  const double S      = S1 + S2;
  const double Sstar  = X2*a1 + X1*a2;
  const double C_Q1   = EOBPars -> C_Q1;
  const double C_Q2   = EOBPars -> C_Q2;
  const double C_Oct1 = EOBPars -> C_Oct1;
  const double C_Oct2 = EOBPars -> C_Oct2;
  const double C_Hex1 = EOBPars -> C_Hex1;
  const double C_Hex2 = EOBPars -> C_Hex2;
  const double X12    = X1 - X2; /* sqrt (1 - 4 nu)*/
  const double X12sq  = SQ (X12); /* (1 - 4 nu)*/
  const double c3     = EOBPars -> cN3LO;

  const int usetidal  = EOBPars -> use_tidal;
  const int usespins  = EOBPars -> use_spins;

  double pphi2   = pphi*pphi;
  double prstar2 = prstar*prstar;
  double prstar3 = prstar2*prstar;
  double prstar4 = prstar3*prstar;

  double Frdot = 0.;
  double pphi2dot = 0.;
  
  double A, dA, d2A, B, dB, d2B, Q, dQ, dQ_dprstar, d2Q, d2Q_dprstar2, ddQ_drdprstar, d3Q_dr2dprstar, d3Q_drdprstar2, d3Q_dprstar3;
  double sqrtAbyB, oosqrtAbyB, dsqrtAbyB_dr, d2sqrtAbyB_d2r, fact;
  double sqA,sqB;
  double ggm[26], G, dG_dr, dG_dprstar,d2G_dr2, d2G_dr_dprstar, d2G_dprstar2,
    d3G_dr2_dprstar, d3G_dr_dprstar2, d3G_dprstar3;
  double  rc, drc_dr, d2rc_dr2, d3rc_dr3, uc, uc2, uc3, uc4;
  double H, Heff, Heff_orb, E, dHeff_dr, dHeff_dprstar, dHeff_dpphi, d2Heff_dr2, EHeff_orb;
  double Adot, prstardot, sqrtAbyBdot, dAbyrc2, d2Abyrc2, Omgdot_0, Heffdot, HSOdot, Edot,
    Heff_orbdot, EHeff_orbdot, Omgdot, Omg2dot, r2dot, r3dot, EHeff_orb2dot,
    HSO2dot, Heff_orb2dot, prstar2dot, Heff2dot, E2dot;
  double Fphi_Newt;
  double oneby_EHeff_orb, oneby_E;
  double rdot2;
  
  /* Computing metric, centrifugal radius and ggm functions*/
  if(usespins) {
    eob_metric_s(r, prstar, dyn, &A, &B, &dA, &d2A, &dB, &d2B, &Q, &dQ, &dQ_dprstar, &d2Q, &ddQ_drdprstar, &d2Q_dprstar2, &d3Q_dr2dprstar, &d3Q_drdprstar2, &d3Q_dprstar3);
    eob_dyn_s_get_rc(r, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc, &drc_dr, &d2rc_dr2, &d3rc_dr3);
    eob_dyn_s_GS(r, rc, drc_dr, d2rc_dr2, d3rc_dr3, aK2, prstar, 0.0, nu, chi1, chi2, X1, X2, c3, ggm);
    G = ggm[2]*S + ggm[3]*Sstar;    // tildeG = GS*S+GSs*Ss
    dG_dr           = ggm[6]*S   + ggm[7]*Sstar;
    dG_dprstar      = ggm[4]  *S + ggm[5]  *Sstar;
    d2G_dr2         = ggm[14] *S + ggm[15] *Sstar;
    d2G_dprstar2    = ggm[16] *S + ggm[17] *Sstar;
    d2G_dr_dprstar  = ggm[18] *S + ggm[19] *Sstar;
    d3G_dprstar3    = ggm[20] *S + ggm[21] *Sstar;
    d3G_dr2_dprstar = ggm[22] *S + ggm[23] *Sstar;
    d3G_dr_dprstar2 = ggm[24] *S + ggm[25] *Sstar;
  } else {
    eob_metric(r, prstar, dyn, &A, &B, &dA, &d2A, &dB, &d2B, &Q, &dQ, &dQ_dprstar, &d2Q, &ddQ_drdprstar, &d2Q_dprstar2, &d3Q_dr2dprstar, &d3Q_drdprstar2, &d3Q_dprstar3);
    rc              = r;   //Nonspinning case: rc = r; G = 0;
    drc_dr          = 1.0;
    d2rc_dr2        = 0.0;
    G               = 0.0;
    dG_dr           = 0.0;
    dG_dprstar      = 0.0;
    d2G_dr2         = 0.0;
    d2G_dprstar2    = 0.0;
    d2G_dr_dprstar  = 0.0;
    d3G_dprstar3    = 0.0;
    d3G_dr2_dprstar = 0.0;
    d3G_dr_dprstar2 = 0.0;
  }
  // Q   = 2.*nu*(4. - 3.*nu); // this is z3 !!
  uc  = 1./rc;
  uc2 = uc*uc;
  uc3 = uc2*uc;
  uc4 = uc3*uc;

  sqrtAbyB = sqrt(A/B);
  oosqrtAbyB = 1./sqrtAbyB;
  sqB      = SQ(B);
  sqA      = SQ(A);
  fact = (dA*B - A*dB)/sqB;
  dsqrtAbyB_dr   = 0.5*oosqrtAbyB*fact;
  /*d2sqrtAbyB_d2r = -0.25*sqrtAbyB*SQ(dA*B - A*dB)/(SQ(A)*SQ(B)) + 0.5/sqrtAbyB*(d2A*SQ(B) - A*B*d2B - 2.*dA*dB*B + 2.*A*SQ(dB))/(SQ(B)*B);*/
  d2sqrtAbyB_d2r = -0.25*sqrtAbyB*SQ(dA*B - A*dB)/(sqA*sqB) + 0.5*oosqrtAbyB*(d2A*sqB - A*B*d2B - 2.*dA*dB*B + 2.*A*SQ(dB))/(sqB*B);

  /* Circular Hamiltonians, ref: arXiv: 1406.6913 */
  if(usespins) {
    eob_ham_s(nu, r, rc, drc_dr, d2rc_dr2, d3rc_dr3, pphi, prstar, S, Sstar, chi1, chi2, X1, X2, aK2, c3, A, dA, d2A, Q, dQ, dQ_dprstar, d2Q, d2Q_dprstar2, &H, &Heff, &Heff_orb, &dHeff_dr, NULL, &dHeff_dpphi, NULL, NULL);
    E = nu*H;
  } else {
    eob_ham(nu, r, pphi, prstar, A, dA, Q, dQ, dQ_dprstar,
	    &H, &Heff_orb, &dHeff_dr, NULL, &dHeff_dpphi);
    Heff = Heff_orb;
    E = nu*H;
  }

  Edot = nu*(rdot*Fr + Omg*Fphi);
  
  Adot = dA*rdot;
  sqrtAbyBdot = dsqrtAbyB_dr*rdot;

  prstardot = - sqrtAbyB/E*(pphi*dG_dr + 1./(2.*Heff_orb)
			    *(dA*(1. + pphi2*uc2 + Q)
	       + A*(-2.*uc3*pphi2*drc_dr + dQ))) + sqrtAbyB*Fr;
    
  dAbyrc2  = (dA*uc2 - 2.*A*uc3*drc_dr);
  d2Abyrc2 = d2A*uc2 - 4.*dA*uc3*drc_dr  + 6.*A*uc4*SQ(drc_dr) - 2.*A*uc3*d2rc_dr2;

  /* shorthands */
  EHeff_orb       = E*Heff_orb;
  oneby_EHeff_orb = 1.0/EHeff_orb;
  oneby_E         = 1.0/E;
  rdot2           = SQ(rdot);
  /*  Omgdot_0 = 1./EHeff_orb*dAbyrc2*rdot*pphi + 1./E*(dG_dr*rdot + dG_dprstar*prstardot);*/
  Omgdot_0 = oneby_EHeff_orb*dAbyrc2*rdot*pphi + oneby_E*(dG_dr*rdot + dG_dprstar*prstardot);
  
  /* Begin iteration */
  int iter = 2; // Hard-fixed to 2: seems to be enough
    
  for(int n = 1; n <= iter; n++){
    Heffdot  = 1./nu*E*Edot;
    HSOdot   = Fphi*G + pphi*(dG_dr*rdot + dG_dprstar*prstardot);
    Heff_orbdot  = Heffdot - HSOdot;
    EHeff_orbdot = Edot*Heff_orb + E*Heff_orbdot;

    Omgdot = Omgdot_0 + oneby_EHeff_orb*A*uc2*(Fphi - pphi*oneby_EHeff_orb*EHeff_orbdot)
      - G*Edot*SQ(oneby_E);

    r2dot = dsqrtAbyB_dr*oosqrtAbyB*rdot2 
          + sqrtAbyB*(oneby_EHeff_orb*(-oneby_EHeff_orb*EHeff_orbdot)*(prstar + 0.5*A*dQ_dprstar)
          + 0.5*oneby_EHeff_orb*(rdot*(dA*dQ_dprstar + A*ddQ_drdprstar) + prstardot*(2. + A*d2Q_dprstar2))
          + oneby_E*(pphi*(rdot*d2G_dr_dprstar + prstardot*d2G_dprstar2) + dG_dprstar*(Fphi - pphi*Edot*oneby_E)));
    
    
    E2dot    = nu*(r2dot*Fr + Omgdot*Fphi);
    Heff2dot = 1./nu*(SQ(Edot) + E*E2dot);
    
    // derivatives of A*(1 + pphi2*u2 + Q) wrt dr2, drdprstar, drdpphi
    double der1 = d2A*(1. + pphi2*uc2 + Q) + 2.*dA*(-2.*uc3*pphi2*drc_dr + dQ) + A*(d2Q + 6.*pphi2*uc4*SQ(drc_dr) - 2.*pphi2*uc3*d2rc_dr2);
    double der2 = dA*dQ_dprstar + A*ddQ_drdprstar;
    double der3 = 2.*dA*pphi*uc2 - 4.*A*pphi*uc3;

    prstar2dot = dsqrtAbyB_dr*oosqrtAbyB*rdot*prstardot 
             + 0.5*sqrtAbyB*oneby_EHeff_orb*(oneby_EHeff_orb*EHeff_orbdot*(dA*(1. + pphi2*uc2 + Q) + A*(dQ - 2.*uc3*pphi2*drc_dr))
             - (der1*rdot + der2*prstardot + der3*Fphi)) 
             - sqrtAbyB*oneby_E*(pphi*(rdot*d2G_dr2 + prstardot*d2G_dr_dprstar) + dG_dr*(Fphi - pphi*Edot*oneby_E))
             + sqrtAbyB*Frdot;
    
    HSO2dot = pphi2dot*G + 2.*Fphi*(dG_dr*rdot + dG_dprstar*prstardot)
      + pphi*(d2G_dr2*rdot + dG_dr*r2dot + 2.*d2G_dr_dprstar*rdot*prstardot
	      + d2G_dprstar2*SQ(prstardot) + dG_dprstar*prstar2dot);
    
    Heff_orb2dot  = Heff2dot - HSO2dot;
    EHeff_orb2dot = E2dot*Heff_orb + 2.*Edot*Heff_orbdot + E*Heff_orb2dot;
    
    Omg2dot = oneby_EHeff_orb*(d2Abyrc2*pphi*SQ(rdot) + dAbyrc2*pphi*r2dot
	      + 2.*dAbyrc2*rdot*(Fphi - pphi*EHeff_orbdot*oneby_EHeff_orb)
	      + A*uc2*(pphi2dot - 2.*Fphi*EHeff_orbdot*oneby_EHeff_orb
              + pphi*(2.0*SQ(EHeff_orbdot*oneby_EHeff_orb) - EHeff_orb2dot*oneby_EHeff_orb)))
              + oneby_E*(d2G_dr2*SQ(rdot) + dG_dr*(r2dot - 2.*rdot*Edot*oneby_E)
              + dG_dprstar*(prstar2dot - 2.*prstardot*Edot*oneby_E)
	      + d2G_dprstar2*SQ(prstardot) + 2*d2G_dr_dprstar*rdot*prstardot
              + G*(2.*SQ(Edot*oneby_E) - E2dot*oneby_E));

    
    double comb1 = sqrtAbyB*oneby_EHeff_orb*(-EHeff_orbdot*oneby_EHeff_orb)*(prstar + 0.5*A*dQ_dprstar); 
    double comb2 = sqrtAbyB*0.5*oneby_EHeff_orb*(rdot*(dA*dQ_dprstar + A*ddQ_drdprstar) + prstardot*(2. + A*d2Q_dprstar2)); 
    double comb3 = sqrtAbyB*oneby_E*(pphi*(rdot*d2G_dr_dprstar + prstardot*d2G_dprstar2) + dG_dprstar*(Fphi - pphi*Edot*oneby_E)); 
    double comb4 = -EHeff_orbdot*oneby_EHeff_orb;
    double comb5 = prstar + 0.5*A*dQ_dprstar;

    double D1 = 2.*rdot*r2dot*oosqrtAbyB*dsqrtAbyB_dr + rdot2*rdot*(-SQ(dsqrtAbyB_dr)*(B/A) + oosqrtAbyB*d2sqrtAbyB_d2r);
    double D2 = comb1*rdot*oosqrtAbyB*dsqrtAbyB_dr + sqrtAbyB*(comb5*oneby_EHeff_orb*(2.*SQ(comb4) - EHeff_orb2dot*oneby_EHeff_orb) 
                + oneby_EHeff_orb*comb4*(prstardot*(1. + 0.5*A*d2Q_dprstar2) + 0.5*rdot*(dA*dQ_dprstar + A*ddQ_drdprstar)));

    double D3 = rdot*dsqrtAbyB_dr*comb2*oosqrtAbyB - EHeff_orbdot*oneby_EHeff_orb*comb2 
                + 0.5*sqrtAbyB*oneby_EHeff_orb*(r2dot*(dA*dQ_dprstar + A*ddQ_drdprstar) 
                + rdot2*(d2A*dQ_dprstar + 2.*dA*ddQ_drdprstar + A*d3Q_dr2dprstar)
                + 2.*rdot*prstardot*(d2Q_dprstar2*dA + A*d3Q_drdprstar2)
                + SQ(prstardot)*A*d3Q_dprstar3 + prstar2dot*(2. + A*d2Q_dprstar2));
    
    double D4 = rdot*dsqrtAbyB_dr*comb3*oosqrtAbyB - Edot*comb3*oneby_E 
                + sqrtAbyB*oneby_E*(d2G_dr_dprstar*(2.*Fphi*rdot + pphi*r2dot - rdot*pphi*Edot*oneby_E)
                                  + d2G_dprstar2*(2.*Fphi*prstardot + pphi*prstar2dot - prstardot*pphi*Edot*oneby_E)
                                  + pphi*(2.*rdot*prstardot*d3G_dr_dprstar2 + rdot2*d3G_dr2_dprstar + SQ(prstardot)*d3G_dprstar3)
                                  + dG_dprstar*(pphi2dot - (Fphi*Edot + pphi*E2dot)*oneby_E + pphi*SQ(Edot*oneby_E)));

    r3dot = D1 + D2 + D3 + D4; 
    
    Fphi_Newt = Fphi_NewtPref(r, Omg, rdot, r2dot, r3dot, Omgdot, Omg2dot);
    Fphi      = Fphi_Newt*Fphi;
    Edot      = nu*(rdot*Fr + Omg*Fphi);
  } // end iteration
  
  /* Saving useful variables */
  dyn->rdot  = rdot;
  dyn->r2dot = r2dot;
  dyn->r3dot = r3dot;
  dyn->r4dot = 0.;
  dyn->r5dot = 0.;
  dyn->Omegadot  = Omgdot;
  dyn->Omega2dot = Omg2dot;
  dyn->Omega3dot = 0.;
  dyn->Omega4dot = 0.;
    
  /* return F_NC */  
  return Fphi_Newt;  
}

/* Generic Newtonian prefactor */
double Fphi_NewtPref(double r, double Omg, double rdot, double r2dot, double r3dot, double Omgdot, double Omg2dot)
{
  double u, u2, u3, u4;
  /*double Omg2, Omg3, Omg4, Omg5;*/
  double invOmg,invOmg2,invOmg3,invOmg4,invOmg5;
  double rdot2, rdot3, rdot4;
  double FphiNewtNC;
  
  u  = 1./r;
  u2 = u*u;
  u3 = u2*u;
  u4 = u3*u;
  /*Omg2 = Omg*Omg;
  Omg3 = Omg2*Omg;
  Omg4 = Omg3*Omg;
  Omg5 = Omg4*Omg;

  FphiNewtNC  =  1. - r3dot*(Omgdot*u/Omg5/8. + 0.5*rdot*u2/Omg4)
    + r2dot*(0.375*rdot*Omgdot*u2/Omg5 - 2.*u/Omg2 + Omg2dot/8.*u/Omg5)
    + SQ(rdot)*(Omg2dot*u2/Omg5/8. + 4.*u2/Omg2) + SQ(rdot)*rdot*0.75*Omgdot*u3/Omg5 
    + 0.75*SQ(Omgdot)/Omg4 + SQ(r2dot)*0.75*u2/Omg4 + SQ(rdot)*SQ(rdot)*u4/Omg4*0.75
    - Omg2dot/Omg3/4. + rdot*Omgdot*u/Omg3*3.;*/

  invOmg  = 1./Omg;
  invOmg2 = invOmg*invOmg;
  invOmg3 = invOmg*invOmg2;
  invOmg4 = invOmg*invOmg3;
  invOmg5 = invOmg*invOmg4;

  rdot2 = SQ(rdot);
  rdot3 = rdot2*rdot;
  rdot4 = rdot2*rdot2;

  FphiNewtNC  =  1. - r3dot*(0.125*Omgdot*u*invOmg5 + 0.5*rdot*u2*invOmg4)
              + r2dot*(0.375*rdot*Omgdot*u2*invOmg5 - 2.0*u*invOmg2             + 0.125*Omg2dot*u*invOmg5)
              + rdot2*(0.125*Omg2dot*u2*invOmg5     + 4.*u2*invOmg2)            + rdot3*0.75*Omgdot*u3*invOmg5 
              + 0.75*SQ(Omgdot)*invOmg4             + SQ(r2dot)*0.75*u2*invOmg4 + rdot4*u4*invOmg4*0.75
              - 0.25*Omg2dot*invOmg3                + 3.0*rdot*Omgdot*u*invOmg3;
  
  return FphiNewtNC;
}
