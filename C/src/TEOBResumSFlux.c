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
double eob_flx_HorizonFlux(double x, double Heff, double jhat, double nu)
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
double eob_flx_HorizonFlux_s(double x, double Heff, double jhat, double nu, double X1, double X2, double chi1, double chi2)
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

/**
 * Function: eob_flx_HorizonFlux_lmr
 * ---------------------------------
 *   Horizon flux for nonspinning large-mass-ratio binaries.
 *   Multipolar (l=2,3,4) Newton-normalized horizon flux with high-PN
 *   rho_lm series. Ported verbatim from branch dev/DALI-LMR.
 *
 *   @param[in] x    : frequency parameter v_phi^2
 *   @param[in] Heff : effective Hamiltonian
 *   @param[in] jhat : angular momentum
 *   @param[in] nu   : symmetric mass ratio
 *   @return[out] hatFH : Newton-normalized horizon flux
 */
double eob_flx_HorizonFlux_lmr(double x, double Heff, double jhat, double nu)
{
  double rhoHlm[9]; /* all l = 2, 3, 4 multipoles */
  double FlmHLO[9];
  double FlmH[9];

  /** Coefficients */
  static double clm[KMAX][33]; // up to v^32 = x^16
  for (int k=0; k<KMAX; k++) clm[k][0] = 1.;
  for (int k=0; k<KMAX; k++) for (int n=1; n<33; n++) clm[k][n] = 0.;
  
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
  
  /* l = 2, m = 2, up to 16PN (v^32); nu-dependence at 1PN only */

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
  clm[1][21] = -94504.18571699700 + 40507.14691869169*logv;
  clm[1][22] = 651309.0091064806 - 159550.1715840284*logv - 60174.92195324077*logv2 + 10691.49605990645*logv3;
  clm[1][23] = -476438.2848549592 + 97913.91958741057*logv - 16340.84469953538*logv2;
  clm[1][24] = -1562948.331428973 + 1200388.348734987*logv - 66225.97145626752*logv2 - 22206.21463789514*logv3 + 184.0461007145514*logv4;
  clm[1][25] = 3511801.821409867 - 1611062.999374397*logv + 175233.1970983554*logv2;
  clm[1][26] = -6337354.248557103 - 2009032.182673084*logv + 1093939.803256274*logv2 - 178418.7635240162*logv3 + 15595.75696428820*logv4;
  clm[1][27] = -4766304.036550063 + 7600814.743331823*logv - 1138594.125903862*logv2;
  clm[1][28] = 35513658.82503024 - 18352458.63919716*logv - 976636.6207783666*logv2 + 956784.9144915892*logv3 - 93017.49458477401*logv4;
  clm[1][29] = -48854509.08321097 - 1636528.194461708*logv + 1894651.658657534*logv2 + 133216.7910743075*logv3;
  clm[1][30] = -328802.9770373756 + 113026621.5379969*logv - 20956827.50366440*logv2 - 1722838.984045762*logv3 + 363441.8194234401*logv4 - 300.0827851650590*logv5;
  clm[1][31] = 249832546.0196665 - 150378059.7134930*logv + 23317729.66375685*logv2 - 2203244.860301497*logv3;
  clm[1][32] = -549826222.4979456 - 80883621.79473536*logv + 120802743.3418255*logv2 - 10216350.77763248*logv3 + 295774.1152446874*logv4 - 76313.57632828431*logv5;

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
  /* (2,2): frozen-log Pade' [30/2] in v of the v^32 series (regular to x = 0.44; the Taylor v^32 has a zero at x = 0.334) */
  {
    const double *c = clm[1];
    const double det = c[30]*c[30] - c[29]*c[31];
    const double b1 = (c[32]*c[29] - c[31]*c[30])/det, b2 = (c[31]*c[31] - c[30]*c[32])/det;
    double P = 0., vm = 1.;
    for (int m = 0; m <= 30; m++) {
      const double pm = c[m] + ((m >= 1) ? b1*c[m-1] : 0.) + ((m >= 2) ? b2*c[m-2] : 0.);
      P += pm*vm; vm *= v;
    }
    rhoHlm[1] = P/(1. + b1*v + b2*v*v);
  }
  rhoHlm[2] = Taylorseries(v,clm[2],12);
  rhoHlm[3] = Taylorseries(v,clm[3],12);
  rhoHlm[4] = Taylorseries(v,clm[4],12);
  rhoHlm[5] = Taylorseries(v,clm[5],24);
  rhoHlm[6] = 1. + c1_42*x + c2_42*x2 + c3_42*x3 + c4_42*x4;
  rhoHlm[7] = 1. + 0.96063*x + 1.45472*x2 + 43.480636*x3 - 35.225828*x4; /* (4,3): Nagar-Akcay fit (the Taylor v^24 has a zero at x = 0.287) */
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

/**
 * Function: eob_flx_HorizonFlux_lmr_s
 * -----------------------------------
 *   Horizon flux for SPINNING large-mass-ratio binaries: the nonspinning
 *   high-PN multipolar LMR horizon flux (eob_flx_HorizonFlux_lmr) with the
 *   spin-dependent (2,2) contribution of the Alvi/Dali spinning horizon flux
 *   (eob_flx_HorizonFlux_s) grafted on.
 *
 *   Construction:
 *     hatFH = eob_flx_HorizonFlux_lmr(...)                          [nonspin, l=2,3,4]
 *           + [ eob_flx_HorizonFlux_s(...) - (X1^4 + X2^4) x^4 ]    [spin-only part]
 *   The subtracted term is the chi->0 limit of eob_flx_HorizonFlux_s (its bare
 *   x^4 = v^8 (2,2) piece; note X1^4 + X2^4 = 1 - 4nu + 2nu^2), which the LMR
 *   flux already carries (dressed with Heff^2 and the high-PN rho_22 series).
 *   Hence at chi1=chi2=0 this reduces EXACTLY (bit-for-bit) to
 *   eob_flx_HorizonFlux_lmr, so the nonspinning path -- and the a6c calibrated on
 *   it -- is unchanged.  For nonzero spin it adds the leading v^5 (2.5PN) Alvi
 *   spin term, a LOWER PN order than the nonspinning v^8 horizon term, which
 *   dominates the horizon absorption for large spin.
 *
 *   NOTE(LMR): this is a pragmatic graft of the leading Alvi spin terms, NOT a
 *   first-principles spinning-LMR horizon flux.  Higher-PN spin and spin-orbit
 *   cross terms in the multipolar rho_lm series are NOT included.
 *
 *   @param[in] x    : frequency parameter v_phi^2
 *   @param[in] Heff : effective Hamiltonian
 *   @param[in] jhat : angular momentum
 *   @param[in] nu   : symmetric mass ratio
 *   @param[in] X1   : mass fraction of body 1
 *   @param[in] X2   : mass fraction of body 2
 *   @param[in] chi1 : dimensionless spin 1
 *   @param[in] chi2 : dimensionless spin 2
 *   @return hatFH   : Newton-normalized horizon flux
 */
double eob_flx_HorizonFlux_lmr_s(double x, double Heff, double jhat, double nu,
                                 double X1, double X2, double chi1, double chi2)
{
  /* nonspinning high-PN multipolar LMR horizon flux (l = 2,3,4) */
  double hatFH = eob_flx_HorizonFlux_lmr(x, Heff, jhat, nu);

  /* spin-only part of the Alvi/Dali (2,2) horizon flux: the full spinning
     expression minus its chi->0 limit (X1^4 + X2^4) x^4, which the LMR term
     already provides.  This difference is exactly 0 at chi1=chi2=0. */
  const double x4        = x*x*x*x;
  const double hatFH_s   = eob_flx_HorizonFlux_s(x, Heff, jhat, nu, X1, X2, chi1, chi2);
  const double hatFH_sns = (X1*X1*X1*X1 + X2*X2*X2*X2) * x4;   /* = (1-4nu+2nu^2) x^4 */
  hatFH += (hatFH_s - hatFH_sns);

  return hatFH;
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
/**
 * Function: eob_flx_hatFlux_highl
 * -------------------------------
 *   Newton-normalized flux at infinity of the multipoles l = 11, 12, 13, m = l, l-1, l-2, which the
 *   KMAX list (l <= 10) does not carry; added to hatf in eob_flx_Flux_s / _ecc (USEFLM_22PN_NU_HL).
 *   As for k < KMAX,  F_lm/F^Newt_22 = C^Newt_lm c_(l+eps)(nu)^2 x^(l+eps-2) [S_eff |T_lm| rho_lm^l]^2,
 *   S_eff = Heff (l+m even) or jhat (l+m odd), |T_lm| the tail factor at w = E Omega, rho_lm the test-mass
 *   22PN series (Fujita) to v^34 or v^32 with the exact 1PN nu-dependence (Damour-Iyer-Nagar).
 *   m = l, l-1, l-2 carry 99.97% of each multipole.
 *
 *   @param[in] x    : frequency parameter v_phi^2
 *   @param[in] w    : E*Omega (the argument eob_flx_Tlm is called with)
 *   @param[in] Heff : effective Hamiltonian
 *   @param[in] jhat : Newton-normalized orbital angular momentum
 *   @param[in] nu   : symmetric mass ratio
 *   @return hatF_highl : sum_lm F_lm / F^Newt_22
 */
double eob_flx_hatFlux_highl(double x, double w, double Heff, double jhat, double nu)
{
  /** Mode list, ordered (l, then m ascending) as LINDEX/MINDEX are */
  static const int LHL[KMAXHL] = {11, 11, 11, 12, 12, 12, 13, 13, 13};
  static const int MHL[KMAXHL] = { 9, 10, 11, 10, 11, 12, 11, 12, 13};
  /* eps = (l+m)%2: 0 0 0 0 1 0 0 1 0 */
  static const int EHL[KMAXHL] = { 0,  1,  0,  0,  1,  0,  0,  1,  0};
  /** C^Newt_lm (exact closed form) */
  static const double CNHL[KMAXHL] = {
    4.5841740393765633,        /* (11,9)  F^Newt ~ x^14 */
    67.048282957598523,        /* (11,10)  F^Newt ~ x^15 */
    1080.6670857342318,        /* (11,11)  F^Newt ~ x^14 */
    9.2754541640946574,        /* (12,10)  F^Newt ~ x^15 */
    120.35831973923660,        /* (12,11)  F^Newt ~ x^16 */
    2035.1394402932438,        /* (12,12)  F^Newt ~ x^15 */
    18.541310598718045,        /* (13,11)  F^Newt ~ x^16 */
    216.26889354643817,        /* (13,12)  F^Newt ~ x^17 */
    3833.1495687399957         /* (13,13)  F^Newt ~ x^16 */
  };
  double clm[KMAXHL][35];
  for (int k = 0; k < KMAXHL; k++) clm[k][0] = 1.;
  for (int k = 0; k < KMAXHL; k++) for (int n = 1; n < 35; n++) clm[k][n] = 0.;

  /** Shorthands */
  const double nu2 = nu*nu;
  const double nu3 = nu*nu2;
  const double nu4 = nu*nu3;
  const double nu5 = nu*nu4;
  const double nu6 = nu*nu5;
  const double nu7 = nu*nu6;

  const double v     = sqrt(x);
  const double logv  = log(v);
  const double logv2 = logv*logv;
  const double logv3 = logv2*logv;
  const double logv4 = logv3*logv;
  const double logv5 = logv4*logv;

  /** Mass-ratio factors c_(l+eps)(nu)^2 of the Newtonian prefactors.
      Same construction as the sp2..sp11 of eob_flx_FlmNewt; all -> 1 at nu = 0. */
  const double c11 = (1. - 4.*nu)*SQ(1. -  8.*nu + 21.*nu2 -  20.*nu3 +   5.*nu4);
  const double c12 =              SQ(1. - 11.*nu + 44.*nu2 -  77.*nu3 +  55.*nu4 - 11.*nu5);
  const double c13 = (1. - 4.*nu)*SQ(1. - 10.*nu + 36.*nu2 -  56.*nu3 +  35.*nu4 -  6.*nu5);
  const double c14 =              SQ(1. - 13.*nu + 65.*nu2 - 156.*nu3 + 182.*nu4 - 91.*nu5 + 13.*nu6);
  /* c_(l+eps) per mode, in the order of LHL/MHL */
  const double CMASS[KMAXHL] = {c11, c12, c11, c12, c13, c12, c13, c14, c13};

  /* x^(l+eps+3) / x^5, i.e. the Newtonian power relative to F^Newt_22 */
  const double x2 = x*x,   x4 = x2*x2, x8 = x4*x4;
  const double x9 = x8*x,  x10 = x9*x, x11 = x10*x, x12 = x11*x;
  const double XPOW[KMAXHL] = {x9, x10, x9, x10, x11, x10, x11, x12, x11};

  /** rho_lm coefficients */
  /* (11,9), test-mass 22PN series to v^34 (17PN rel.); nu at 1PN (DIN closed form) */
  clm[0][2] = (-381. + 3260.*nu - 9631.*nu2 + 11676.*nu3 - 5485.*nu4 + 796.*nu5)/(330. - 2640.*nu + 6930.*nu2 - 6600.*nu3 + 1650.*nu4);
  clm[0][4] = -1.0613908818454273;
  clm[0][6] = 4.2817420376293885 - 2.4309840357271187*logv;
  clm[0][8] = 1.6187282938784409 + 2.8066815685213098*logv;
  clm[0][10] = -6.6880295243924624 + 2.5802242894325623*logv;
  clm[0][12] = 7.8378496435668755 - 17.743887192191546*logv + 2.9548416909800546*logv2;
  clm[0][14] = 33.326475570477999 + 4.5335352049364906*logv - 3.4114990432224267*logv2;
  clm[0][16] = -12.949715428359781 + 24.043838271979858*logv - 3.1362420281029538*logv2;
  clm[0][18] = 27.317588814525081 - 101.38410002971593*logv + 30.483236613151330*logv2 - 2.3943909929578124*logv3;
  clm[0][20] = 284.52959291598610 - 34.095886865269169*logv - 15.804037557728511*logv2 + 2.7644332373240197*logv3;
  clm[0][22] = 61.459110559236542 + 134.58742912571321*logv - 38.688118527607394*logv2 + 2.5413847674982408*logv3;
  clm[0][24] = 362.12047441803269 - 811.24920865455876*logv + 250.20493444943240*logv2 - 31.926048497397324*logv3 + 1.4551815697923115*logv4;
  clm[0][26] = 2900.2272004853951 - 476.73249422937462*logv - 46.646804214510132*logv2 + 21.147615724889159*logv3 - 1.6800732669420323*logv4;
  clm[0][28] = 1556.2916386784125 + 784.27839194325438*logv - 317.46844460065173*logv2 + 39.018220429545452*logv3 - 1.5445164496070747*logv4;
  clm[0][30] = 5407.0565820924174 - 8502.8790801872105*logv + 2380.1267478408314*logv2 - 327.43711232713768*logv3 + 23.793667373753133*logv4 - 0.70750463304988744*logv5;
  clm[0][32] = 34340.769197318248 - 5805.7915274121088*logv - 69.518786768371415*logv2 + 134.34897143061810*logv3 - 17.921686600216862*logv4 + 0.81684625815759731*logv5;
  clm[0][34] = 26069.190633831271 + 4885.4630507329920*logv - 2664.5459028199801*logv2 + 405.08285599185434*logv3 - 28.373457888194655*logv4 + 0.75093896638254548*logv5;

  /* (11,10), test-mass 22PN series to v^34 (17PN rel.); nu at 1PN (DIN closed form) */
  clm[1][2] = (12120. - 140491.*nu + 609752.*nu2 - 1227083.*nu3 + 1153724.*nu4 - 455455.*nu5 + 57200.*nu6)/(-9438. + 103818.*nu - 415272.*nu2 + 726726.*nu3 - 519090.*nu4 + 103818.*nu5);
  clm[1][4] = -1.2496689471780496;
  clm[1][6] = 5.1141043870320867 - 3.0012148589223688*logv;
  clm[1][8] = 1.8027005226075813 + 3.8540712110764049*logv;
  clm[1][10] = -10.526702768023415 + 3.7505250130046354*logv;
  clm[1][12] = 8.5275760798009968 - 26.528286255701247*logv + 4.5036453147082070*logv2;
  clm[1][14] = 47.467876768062528 + 8.9464252280624107*logv - 5.7834478930137179*logv2;
  clm[1][16] = -42.200955911899002 + 45.563895895770624*logv - 5.6280656988947612*logv2;
  clm[1][18] = -6.3691881713268990 - 178.58940044830731*logv + 56.584974630770890*logv2 - 4.5054690792727929*logv3;
  clm[1][20] = 437.12855449984445 - 39.563473192434893*logv - 34.968868091162519*logv2 + 5.7857899174386788*logv3;
  clm[1][22] = -229.91683053066220 + 364.08574105929471*logv - 89.338505792244167*logv2 + 5.6303448008380879*logv3;
  clm[1][24] = -176.82030454365826 - 1656.4925012335602*logv + 560.07350304677263*logv2 - 73.391113734834234*logv3 + 3.3804701867821975*logv4;
  clm[1][26] = 4683.0944036442242 - 615.29502995847400*logv - 175.29232270958583*logv2 + 56.535548989345302*logv3 - 4.3410996677050471*logv4;
  clm[1][28] = -2585.8214832210303 + 3532.9393480623710*logv - 980.73760850543388*logv2 + 110.34815861479025*logv3 - 4.2244686192828936*logv4;
  clm[1][30] = -3089.7911221195460 - 20169.676636073262*logv + 6391.0381750592528*logv2 - 915.01833319579931*logv3 + 67.658141213797345*logv4 - 2.0291034709429613*logv5;
  clm[1][32] = 57992.475321727686 - 6990.8392228469322*logv - 1385.0694729317954*logv2 + 490.40472321048007*logv3 - 58.589768352176083*logv4 + 2.6057145653558690*logv5;
  clm[1][34] = -41427.943438543729 + 43820.156058883900*logv - 11564.223916973866*logv2 + 1493.8263434535543*logv3 - 98.531109472336698*logv4 + 2.5357075982486168*logv5;

  /* (11,11), test-mass 22PN series to v^34 (17PN rel.); nu at 1PN (DIN closed form) */
  clm[2][2] = (-421. + 3660.*nu - 11071.*nu2 + 13916.*nu3 - 6885.*nu4 + 1036.*nu5)/(330. - 2640.*nu + 6930.*nu2 - 6600.*nu3 + 1650.*nu4);
  clm[2][4] = -1.5208887352826747;
  clm[2][6] = 5.7590934035738920 - 3.6314699792960663*logv;
  clm[2][8] = 2.5707727514583068 + 4.6328753372231633*logv;
  clm[2][10] = -14.811918315957267 + 5.5230617840285950*logv;
  clm[2][12] = 7.2308182655956817 - 37.282261681495856*logv + 6.5937871052642859*logv2;
  clm[2][14] = 66.794164772280117 + 11.546281916849435*logv - 8.4120738524735284*logv2;
  clm[2][16] = -81.520849949655860 + 78.683379829502183*logv - 10.028416531248608*logv2;
  clm[2][18] = -78.474802808257288 - 290.27924809355922*logv + 97.415178234114144*logv2 - 7.9817133075455883*logv3;
  clm[2][20] = 675.41438149374000 - 68.074872399399898*logv - 58.881104372658333*logv2 + 10.182731219626342*logv3;
  clm[2][22] = -731.43755752333233 + 796.66358267420128*logv - 188.06969572439989*logv2 + 12.139297857701904*logv3;
  clm[2][24] = -1427.2252187365875 - 3148.0096479761859*logv + 1140.4223064919063*logv2 - 153.89643141312867*logv3 + 7.2463380649249286*logv4;
  clm[2][26] = 8208.8409915081597 - 1001.0049432413335*logv - 364.11565185317591*logv2 + 117.17206697904871*logv3 - 9.2445706828284695*logv4;
  clm[2][28] = -11074.474528178709 + 10059.256096427064*logv - 2559.2676298672741*logv2 + 282.37248425287479*logv3 - 11.020873934994379*logv4;
  clm[2][30] = -24803.188777315668 - 45514.181038681008*logv + 15499.425529409521*logv2 - 2285.0841251147604*logv3 + 172.37931094926577*logv4 - 5.2629718285210454*logv5;
  clm[2][32] = 120844.11827314711 - 10632.484717272154*logv - 4018.8400657154734*logv2 + 1238.0144006870001*logv3 - 148.04517736388248*logv4 + 6.7142761812344246*logv5;
  clm[2][34] = -196652.72525361135 + 157364.73700145210*logv - 37917.074296801881*logv2 + 4691.5495610620009*logv3 - 306.03167734536204*logv4 + 8.0043945681077185*logv5;

  /* (12,10), test-mass 22PN series to v^34 (17PN rel.); nu at 1PN (DIN closed form) */
  clm[3][2] = (3274. - 37882.*nu + 164021.*nu2 - 329043.*nu3 + 308165.*nu4 - 121264.*nu5 + 15301.*nu6)/(-2808. + 30888.*nu - 123552.*nu2 + 216216.*nu3 - 154440.*nu4 + 30888.*nu5);
  clm[3][4] = -1.0967032228518206;
  clm[3][6] = 4.4097812987278360 - 2.5261984392419175*logv;
  clm[3][8] = 1.6910720694686329 + 2.9454322258112671*logv;
  clm[3][10] = -7.2271121299114033 + 2.7704899698798499*logv;
  clm[3][12] = 7.9761581411642663 - 19.076051507687754*logv + 3.1908392772141500*logv2;
  clm[3][14] = 35.679807507661939 + 4.9811109257723407*logv - 3.7203731458686350*logv2;
  clm[3][16] = -16.282156333340851 + 26.960631693087799*logv - 3.4994037189229320*logv2;
  clm[3][18] = 24.337669957019668 - 112.49095785584542*logv + 34.118988173714026*logv2 - 2.6868977339900646*logv3;
  clm[3][20] = 310.00261023777371 - 36.692765906250233*logv - 17.979213826237279*logv2 + 3.1328002781636295*logv3;
  clm[3][22] = 37.624474659322279 + 161.37756161065875*logv - 45.047352459102456*logv2 + 2.9467294043401575*logv3;
  clm[3][24] = 329.66527047343552 - 927.50341904505979*logv + 290.21443806818445*logv2 - 37.171351647003114*logv3 + 1.6969092155020865*logv4;
  clm[3][26] = 3216.3983083249968 - 522.16957855673961*logv - 57.871921048983489*logv2 + 24.981399809053239*logv3 - 1.9785187932884014*logv4;
  clm[3][28] = 1333.7827314311844 + 1051.6238329065852*logv - 390.25383689096639*logv2 + 47.190020134507060*logv3 - 1.8610058055280928*logv4;
  clm[3][30] = 5137.0781208850634 - 10016.756074121687*logv + 2859.3627556717345*logv2 - 395.62981507088517*logv3 + 28.806403997003907*logv4 - 0.85734588234731952*logv5;
  clm[3][32] = 38853.904151395182 - 6438.8504169817293*logv - 147.43288460125053*logv2 + 167.40868133162251*logv3 - 21.992522995489379*logv4 + 0.99962621752319234*logv5;
  clm[3][34] = 23977.546186685137 + 7938.9516255070263*logv - 3503.3554714473105*logv2 + 514.67726057031833*logv3 - 35.649200678957608*logv4 + 0.94025399226904309*logv5;

  /* (12,11), test-mass 22PN series to v^32 (16PN rel.); nu at 1PN (DIN closed form) */
  clm[4][2] = (11638. - 123339.*nu + 486244.*nu2 - 883740.*nu3 + 748440.*nu4 - 267113.*nu5 + 30184.*nu6)/(-9072. + 90720.*nu - 326592.*nu2 + 508032.*nu3 - 317520.*nu4 + 54432.*nu5);
  clm[4][4] = -1.2675381669671383;
  clm[4][6] = 5.1838158439169683 - 3.0567001114827202*logv;
  clm[4][8] = 1.8847207161996381 + 3.9212826165603943*logv;
  clm[4][10] = -10.816629003829922 + 3.8744840562770543*logv;
  clm[4][12] = 8.6528350660004700 - 27.464568905615221*logv + 4.6717077857692370*logv2;
  clm[4][14] = 49.519039731860598 + 9.1446432245718274*logv - 5.9930925055977050*logv2;
  clm[4][16] = -44.194326540408124 + 47.790968571107650*logv - 5.9215679233800469*logv2;
  clm[4][18] = -8.2774215247076812 - 188.27155339797868*logv + 59.733677997598068*logv2 - 4.7600032365251728*logv3;
  clm[4][20] = 463.51084860172733 - 42.938407579144574*logv - 36.757296407888828*logv2 + 6.1063621766622533*logv3;
  clm[4][22] = -247.32545946203406 + 389.53938357140231*logv - 95.550529026202338*logv2 + 6.0334857771827631*logv3;
  clm[4][24] = -207.30833007419474 - 1780.5001790311047*logv + 602.56904748007484*logv2 - 78.956479999967843*logv3 + 3.6374756059361512*logv4;
  clm[4][26] = 5063.0598193704792 - 671.64454984502128*logv - 186.68461264340042*logv2 + 60.663654368426745*logv3 - 4.6663294865393438*logv4;
  clm[4][28] = -2831.0534967479303 + 3862.3363049581243*logv - 1069.8045824001354*logv2 + 120.29106183005155*logv3 - 4.6106391619359896*logv4;
  clm[4][30] = -3553.5833423913710 - 22103.497986753696*logv + 7011.4492026084876*logv2 - 1003.5107156794564*logv3 + 74.163425846757192*logv4 - 2.2237344180361417*logv5;
  clm[4][32] = 64102.796762855283 - 7787.8886663219075*logv - 1499.4237577469514*logv2 + 535.52377988193795*logv3 - 64.095408233365977*logv4 + 2.8527139723439833*logv5;

  /* (12,12), test-mass 22PN series to v^34 (17PN rel.); nu at 1PN (DIN closed form) */
  clm[5][2] = (398. - 4654.*nu + 20449.*nu2 - 41899.*nu3 + 40469.*nu4 - 16588.*nu5 + 2145.*nu6)/(-312. + 3432.*nu - 13728.*nu2 + 24024.*nu3 - 17160.*nu4 + 3432.*nu5);
  clm[5][4] = -1.5184909614437639;
  clm[5][6] = 5.7734526706595900 - 3.6377257525083612*logv;
  clm[5][8] = 2.5885112374483832 + 4.6404322099305377*logv;
  clm[5][10] = -14.815430672565725 + 5.5238536753951609*logv;
  clm[5][12] = 7.3413548229978890 - 37.458469876918588*logv + 6.6165243252312614*logv2;
  clm[5][14] = 67.380148631047947 + 11.575951208174376*logv - 8.4403098764168014*logv2;
  clm[5][16] = -81.522970320734157 + 78.883113875165999*logv - 10.047132384036469*logv2;
  clm[5][18] = -77.863214962953301 - 292.94717767838548*logv + 98.063450484834883*logv2 - 8.0230336433305891*logv3;
  clm[5][20] = 684.08574696385468 - 69.276999708157970*logv - 59.237083388467907*logv2 + 10.234510865530687*logv3;
  clm[5][22] = -733.07506845724245 + 800.37871930975849*logv - 188.92847734793252*logv2 + 12.182904070756730*logv3;
  clm[5][24] = -1430.6036104192376 - 3190.3870009866732*logv + 1152.4910322630303*logv2 - 155.20366709681683*logv3 + 7.2963990243961665*logv4;
  clm[5][26] = 8328.5401176279812 - 1024.6788372913297*logv - 366.53829951004741*logv2 + 118.12798824495132*logv3 - 9.3075859349669047*logv4;
  clm[5][28] = -11192.605799158934 + 10135.211632321367*logv - 2577.7689981063861*logv2 + 284.20264435216997*logv3 - 11.079515969632676*logv4;
  clm[5][30] = -25153.961867354121 - 46244.971716547585*logv + 15720.518827087474*logv2 - 2313.0541172138271*logv3 + 174.15432074940904*logv4 - 5.3084597263245635*logv5;
  clm[5][32] = 122338.82619565794 - 10962.462381797728*logv - 4042.8206080734646*logv2 + 1251.1455657968753*logv3 - 149.53467859940709*logv4 + 6.7716890098627444*logv5;
  clm[5][34] = -201054.55994508763 + 159200.46281609289*logv - 38303.633257216156*logv2 + 4735.5529363063081*logv3 - 308.58399479518202*logv4 + 8.0608481136120861*logv5;

  /* (13,11), test-mass 22PN series to v^32 (16PN rel.); nu at 1PN (DIN closed form) */
  clm[6][2] = (18609. - 196918.*nu + 774820.*nu2 - 1404840.*nu3 + 1186619.*nu4 - 422954.*nu5 + 48048.*nu6)/(-15834. + 158340.*nu - 570024.*nu2 + 886704.*nu3 - 554190.*nu4 + 95004.*nu5);
  clm[6][4] = -1.1268269613306581;
  clm[6][6] = 4.5195388482846049 - 2.6086566935373301*logv;
  clm[6][8] = 1.7569475238693353 + 3.0658388537347592*logv;
  clm[6][10] = -7.7023291192903872 + 2.9395046951335516*logv;
  clm[6][12] = 8.0864952324018405 - 20.267259780868002*logv + 3.4025448723685580*logv2;
  clm[6][14] = 37.812247735688306 + 5.3797634555912244*logv - 3.9988605235509976*logv2;
  clm[6][16] = -19.359378893172909 + 29.645221501867829*logv - 3.8340792993222742*logv2;
  clm[6][18] = 21.353158088168941 - 122.75827175350772*logv + 37.492389153319957*logv2 - 2.9586904854551198*logv3;
  clm[6][20] = 333.61961958544370 - 39.081257985139818*logv - 20.012048736437806*logv2 + 3.4772180904278341*logv3;
  clm[6][22] = 13.866369765627912 + 187.18138790249759*logv - 51.126685052107756*logv2 + 3.3339322092433224*logv3;
  clm[6][24] = 294.01337019207020 - 1038.5443978352192*logv + 328.65212462758177*logv2 - 42.216427666225025*logv3 + 1.9295519347469278*logv4;
  clm[6][26] = 3516.0285504624883 - 564.11374762798333*logv - 68.938280462433147*logv2 + 28.701414405009944*logv3 - 2.2677170616209157*logv4;
  clm[6][28] = 1087.8171625386094 + 1323.5061970048425*logv - 462.91064778991870*logv2 + 55.291580649233739*logv3 - 2.1742711433605729*logv4;
  clm[6][30] = 4779.6643136767615 - 11509.108770049737*logv + 3335.9864293001921*logv2 - 463.57547316406785*logv3 + 33.802493893286993*logv4 - 1.0067077140210958*logv5;
  clm[6][32] = 43228.807116878205 - 7024.6626623513014*logv - 231.60068501073571*logv2 + 200.93730312307734*logv3 - 26.087419431849613*logv4 + 1.1831390583692416*logv5;

  /* (13,12), test-mass 22PN series to v^32 (16PN rel.); nu at 1PN (DIN closed form) */
  clm[7][2] = (-62808. + 854395.*nu - 4562610.*nu2 + 12123540.*nu3 - 16792360.*nu4 + 11649729.*nu5 - 3577470.*nu6 + 355160.*nu7)/(49010. - 637130.*nu + 3185650.*nu2 - 7645560.*nu3 + 8919820.*nu4 - 4459910.*nu5 + 637130.*nu6);
  clm[7][4] = -1.2830613919506572;
  clm[7][6] = 5.2428331414248523 - 3.1045170567716987*logv;
  clm[7][8] = 1.9573037133950513 + 3.9785453438424169*logv;
  clm[7][10] = -11.068459716412558 + 3.9832859761960532*logv;
  clm[7][12] = 8.7520955027634804 - 28.282885817195326*logv + 4.8190130778932053*logv2;
  clm[7][14] = 51.325321429055717 + 9.3101584143645663*logv - 6.1757309405492030*logv2;
  clm[7][16] = -45.963289809384837 + 49.767197099011846*logv - 6.1830896275500769*logv2;
  clm[7][18] = -10.089375299866095 - 196.89531790577592*logv + 62.539419960470972*logv2 - 4.9869027657084461*logv3;
  clm[7][20] = 487.10368734772376 - 46.002629851734694*logv - 38.335817790901035*logv2 + 6.3908873476559087*logv3;
  clm[7][22] = -263.35255402474102 + 412.58737060584216*logv - 101.16406013607504*logv2 + 6.3985024040924611*logv3;
  clm[7][24] = -236.81221045805543 - 1893.2291293673075*logv + 641.16548161018999*logv2 - 84.004598447385535*logv3 + 3.8704811741509573*logv4;
  clm[7][26] = 5408.8323963917745 - 723.15387437068446*logv - 196.89069130158426*logv2 + 64.387541759402698*logv3 - 4.9601546946760524*logv4;
  clm[7][28] = -3064.2478672288104 + 4167.5089348530438*logv - 1151.8661423952097*logv2 + 129.43410892640626*logv3 - 4.9660649628249414*logv4;
  clm[7][30] = -4016.4720842120995 - 23895.542346936737*logv + 7586.0761379800955*logv2 - 1085.3139111766770*logv3 + 80.167140584688712*logv4 - 2.4031949646130797*logv5;
  clm[7][32] = 69760.245906033027 - 8526.7949144535752*logv - 1604.1724237231640*logv2 + 576.98339695289638*logv3 - 69.155976270152338*logv4 + 3.0797769707696044*logv5;

  /* (13,13), test-mass 22PN series to v^32 (16PN rel.); nu at 1PN (DIN closed form) */
  clm[8][2] = (20193. - 215926.*nu + 861940.*nu2 - 1594920.*nu3 + 1386203.*nu4 - 511658.*nu5 + 59136.*nu6)/(-15834. + 158340.*nu - 570024.*nu2 + 886704.*nu3 - 554190.*nu4 + 95004.*nu5);
  clm[8][4] = -1.5165805806898221;
  clm[8][6] = 5.7857282710735292 - 3.6434957124612297*logv;
  clm[8][8] = 2.6067567115154855 + 4.6465270254976387*logv;
  clm[8][10] = -14.816723606351289 + 5.5256548433453289*logv;
  clm[8][12] = 7.4396788073756104 - 37.617476546228643*logv + 6.6375305033616819*logv2;
  clm[8][14] = 67.927065979442179 + 11.592080114632903*logv - 8.4648006476179388*logv2;
  clm[8][16] = -81.474174050326790 + 79.064665913691768*logv - 10.066349865134667*logv2;
  clm[8][18] = -77.219907223916142 - 295.38938209747396*logv + 98.656166626382330*logv2 - 8.0612713101096386*logv3;
  clm[8][20] = 692.40128232023433 - 70.480667673369956*logv - 59.538121383698884*logv2 + 10.280488288811667*logv3;
  clm[8][22] = -733.42018374439047 + 803.64458715881362*logv - 189.72531636667970*logv2 + 12.225567524584279*logv3;
  clm[8][24] = -1430.7593312109142 - 3230.0506382519337*logv + 1163.6083506465792*logv2 - 156.40649739569470*logv3 + 7.3428018638427971*logv4;
  clm[8][26] = 8449.8159653662735 - 1048.4060057802853*logv - 368.45497907490319*logv2 + 118.97033154023934*logv3 - 9.3642287505732981*logv4;
  clm[8][28] = -11272.283802135761 + 10201.206656897764*logv - 2594.6647601792645*logv2 + 285.91087390482200*logv3 - 11.135950714557017*logv4;
  clm[8][30] = -25396.410546843002 - 46945.476361135983*logv + 15926.736371365412*logv2 - 2338.9192718500795*logv3 + 175.79431544325515*logv4 - 5.3506934216727115*logv5;
  clm[8][32] = 124023.99113154743 - 11315.531493929921*logv - 4059.6898899896347*logv2 + 1262.7032099136909*logv3 - 150.86959697228455*logv4 + 6.8237054606439979*logv5;

  /** Sum the nine multipoles */
  const double invFNewt22 = 5./32.;
  double hatF = 0.;

  for (int k = 0; k < KMAXHL; k++) {

    const int l = LHL[k];
    const int m = MHL[k];

    /* rho_lm: Taylor series in v; entries above each mode's stored
       order (v^34 or v^32) are zero, so one bound serves all nine */
    const double rholm = Taylorseries(v, clm[k], 34);

    /* effective source */
    const double Seff = (EHL[k] == 0) ? Heff : jhat;

    /* tail factor |T_lm|^2, as in eob_flx_Tlm */
    const double hhatk = m*w;
    const double xt2   = 4.*hhatk*hhatk;
    double prod = 1.;
    for (int j = 1; j <= l; j++) prod *= (j*j + xt2);
    double y = 4.*Pi*hhatk;
    y /= (1. - exp(-y));
    const double MTlm2 = y * prod / SQ(fact(l));

    /* F_lm / F^Newt_22 */
    hatF += invFNewt22 * CNHL[k] * CMASS[k] * XPOW[k]
            * SQ(Seff) * MTlm2 * gsl_pow_int(rholm, 2*l);
  }

  return hatF;
}


/* Multipoles in the flux sum: l <= 10 only for the 22PN flm, which provide them; the other flm stop at l = 8 */
static int eob_flx_kmax(void)
{
  return (EOBPars->use_flm == USEFLM_22PN || EOBPars->use_flm == USEFLM_22PN_NU || EOBPars->use_flm == USEFLM_22PN_NU_HL) ? KMAX : KMAX_L8;
}

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
  double rholm[KMAX]={0.}, flm[KMAX]={0.}, FNewtlm[KMAX]={0.}, MTlm[KMAX]={0.}, hlmTidal[KMAX]={0.}, hlmNQC[KMAX]={0.};
  double Modhhatlm[KMAX]={0.};  

  /** Newtonian flux */
  eob_flx_FlmNewt(x, nu, FNewtlm);

  /* Correct amplitudes for specific multipoles and cases */
  if (usespins) {
    /* Correct (2,1), (3,1) and (3,3) ( sp2 = 1 ) */
    double x6 = gsl_pow_int(x, 6);
    FNewtlm[0] = CNlm[0] * x6; /* (2,1) */
    FNewtlm[2] = CNlm[2] * x6; /* (3,1) */
    FNewtlm[4] = CNlm[4] * x6; /* (3,3) */
    /* Correct (4,1), (4,3)  ( sp4 = (1-2nu)^2 ) */
    double sp4x8 = SQ((1-2*nu)) * gsl_pow_int(x, 8);
    FNewtlm[5]  = CNlm[5] * sp4x8; /* (4,1) */
    FNewtlm[7]  = CNlm[7] * sp4x8; /* (4,3) */
    
    /* Correcting (5,5) for Higher Modes */
    if (EOBPars->use_flm == USEFLM_HM || EOBPars->use_flm == USEFLM_HM_4PN22 || EOBPars->use_flm == USEFLM_HM_6PN3p3 ||
        EOBPars->use_flm == USEFLM_22PN || EOBPars->use_flm == USEFLM_22PN_NU || EOBPars->use_flm == USEFLM_22PN_NU_HL){
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
  for (int k = eob_flx_kmax(); k < KMAX; k++) Modhhatlm[k] = 0.;
  for (int k = KMAX; k--;) sum_k += SQ(Modhhatlm[k]) * FNewtlm[k];
  
  /* Normalize to the 22 Newtonian multipole */
  double hatf = sum_k/(FNewt22);
    
  /* Horizon flux */ 
  if (!(usetidal)) {
    double hatFH;
    if (EOBPars->use_hflx == HFLX_LMR) {
      /* LMR high-PN multipolar horizon flux.  eob_flx_HorizonFlux_lmr_s grafts
         the spin-dependent (2,2) Alvi terms (leading v^5, 2.5PN) onto the
         nonspinning multipolar LMR flux and reduces EXACTLY to
         eob_flx_HorizonFlux_lmr at chi1=chi2=0, so the nonspinning path (and the
         a6c calibrated on it) is unchanged, while genuinely spinning runs now
         retain the dominant spin horizon terms the bare _lmr flux dropped.
         (History: the LMR term had ended up only on the use_spins==0 branch in
         the Dali rebranch, but eob_set_params forces use_spins=ALIGNED for chi=0
         binaries, so it was silently dropped from the nonspinning path.) */
      hatFH = eob_flx_HorizonFlux_lmr_s(x, Heff, jhat, nu, X1, X2, chi1, chi2);
    } else if (usespins) {
      hatFH = eob_flx_HorizonFlux_s(x, Heff, jhat, nu, X1, X2, chi1, chi2);
    } else {
      hatFH = eob_flx_HorizonFlux(x, Heff, jhat, nu);
    }
    hatf += hatFH;
  }

  /* l = 11, 12, 13 (m = l, l-1, l-2), flux only */
  if (EOBPars->use_flm == USEFLM_22PN_NU_HL) {
    hatf += eob_flx_hatFlux_highl(x, E*Omega, Heff, jhat, nu);
  }

  /* return Fphi */  
  return (-32./5. * nu * gsl_pow_int(r_omega,4) * gsl_pow_int(Omega,5) * hatf);  
}

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
  double rholm[KMAX]={0.}, flm[KMAX]={0.}, FNewtlm[KMAX]={0.}, MTlm[KMAX]={0.}, hlmTidal[KMAX]={0.}, hlmNQC[KMAX]={0.};
  double Modhhatlm[KMAX], Flm[KMAX];

  /** Newtonian flux */
  eob_flx_FlmNewt(x, nu, FNewtlm);

  /* Correct amplitudes for specific multipoles and cases */
  if (usespins) {
    /* Correct (2,1), (3,1) and (3,3) ( sp2 = 1 ) */
    double x6 = gsl_pow_int(x, 6);
    FNewtlm[0] = CNlm[0] * x6; /* (2,1) */
    FNewtlm[2] = CNlm[2] * x6; /* (3,1) */
    FNewtlm[4] = CNlm[4] * x6; /* (3,3) */
    /* Correct (4,1), (4,3)  ( sp4 = (1-2nu)^2 ) */
    double sp4x8 = SQ((1-2*nu)) * gsl_pow_int(x, 8);
    FNewtlm[5]  = CNlm[5] * sp4x8; /* (4,1) */
    FNewtlm[7]  = CNlm[7] * sp4x8; /* (4,3) */
    
    /* Correcting (5,5) for Higher Modes */
    if (EOBPars->use_flm == USEFLM_HM || EOBPars->use_flm == USEFLM_HM_4PN22 || EOBPars->use_flm == USEFLM_HM_6PN3p3 ||
        EOBPars->use_flm == USEFLM_22PN || EOBPars->use_flm == USEFLM_22PN_NU || EOBPars->use_flm == USEFLM_22PN_NU_HL){
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
  for (int k = eob_flx_kmax(); k < KMAX; k++) Modhhatlm[k] = 0.;
  for (int k = 0; k < KMAX; k++) {
    Flm[k] = SQ(Modhhatlm[k]) * FNewtlm[k];
    sum_k += Flm[k];    
  }

  /* Normalize to the 22 Newtonian multipole */
  double hatf = sum_k/(FNewt22);

  /* l = 11, 12, 13 (m = l, l-1, l-2), flux only; no non-circular factor (as every multipole but the (2,2)),
     so it is added here and again after the re-sum below */
  double hatf_highl = 0.;
  if (EOBPars->use_flm == USEFLM_22PN_NU_HL) {
    hatf_highl = eob_flx_hatFlux_highl(x, E*Omega, Heff, jhat, nu);
  }
  hatf += hatf_highl;

  /* Horizon flux */
  double hatFH = 0.;
  if (!(usetidal)) {
    if (EOBPars->use_hflx == HFLX_LMR) {
      /* LMR high-PN multipolar horizon flux.  eob_flx_HorizonFlux_lmr_s grafts
         the spin-dependent (2,2) Alvi terms (leading v^5, 2.5PN) onto the
         nonspinning multipolar LMR flux and reduces EXACTLY to
         eob_flx_HorizonFlux_lmr at chi1=chi2=0, so the nonspinning path (and the
         a6c calibrated on it) is unchanged, while genuinely spinning runs now
         retain the dominant spin horizon terms the bare _lmr flux dropped.
         (History: the LMR term had ended up only on the use_spins==0 branch in
         the Dali rebranch, but eob_set_params forces use_spins=ALIGNED for chi=0
         binaries, so it was silently dropped from the nonspinning path.) */
      hatFH = eob_flx_HorizonFlux_lmr_s(x, Heff, jhat, nu, X1, X2, chi1, chi2);
    } else if (usespins) {
      hatFH = eob_flx_HorizonFlux_s(x, Heff, jhat, nu, X1, X2, chi1, chi2);
    } else {
      hatFH = eob_flx_HorizonFlux(x, Heff, jhat, nu);
    }
  }

  /* Compute circular, asymptotic and horizon Fphi separately, then sum */
  double Fphi_lo  = -32./5. * nu * gsl_pow_int(r_omega,4) * gsl_pow_int(Omega,5);
  double Fphi_inf = Fphi_lo * hatf;
  double Fphi_H   = Fphi_lo * hatFH;
  *Fphi           = Fphi_inf + Fphi_H;

  /* Compute Fr using the infinity Fphi */
  *Fr = eob_flx_Fr(r, pr_star, pphi, dyn, Fphi_inf);
  
  /* Compute non-circular Fphi */
  double Fphi_NC[KMAX];
  for (int k = 0; k < KMAX; k++) Fphi_NC[k] = 1.;
  eob_flx_Fphi_ecc(r, pr_star, pphi, Omega, rdot, *Fphi, Fphi_lo, FNewtlm, Flm, Fphi_H, *Fr, dyn, Fphi_NC);
  // To recover old configuration used for arXiv:2001.11736, one should apply this to all multipoles -> for (int k = KMAX; k--;) Fphi_NC[k] = fphi_nc;
  
  /* Adding non-circular corrections and re-compute flux */
  sum_k = 0.;
  for (int k = 0; k < KMAX; k++) sum_k += Flm[k] * Fphi_NC[k];

  hatf = sum_k/(FNewt22);
  hatf += hatf_highl;

  /* Compute non-circular Fphi, infinity */
  *Fphi = Fphi_lo * hatf;
  
  /* Re-compute Fr using the generic Fphi */
  *Fr = eob_flx_Fr(r, pr_star, pphi, dyn, *Fphi);

  /* Add horizon Fphi */
  *Fphi += Fphi_H;
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

  const double u  = 1./r;
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
void eob_flx_Fphi_ecc(double r, double prstar, double pphi, double Omg, double rdot, double Fphi, double Fphi_lo, double *FlmNewt, double *Flm, double Fphi_H, double Fr, Dynamics *dyn, double *hatflm_NC)
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
  double  rc, drc_dr, d2rc_dr2, uc, uc2, uc3, uc4;
  double H, Heff, Heff_orb, E, dHeff_dr, dHeff_dprstar, dHeff_dpphi, d2Heff_dr2, EHeff_orb;
  double Adot, prstardot, sqrtAbyBdot, dAbyrc2, d2Abyrc2, Omgdot_0, Heffdot, HSOdot, Edot,
    Heff_orbdot, EHeff_orbdot, Omgdot, Omg2dot, r2dot, r3dot, EHeff_orb2dot,
    HSO2dot, Heff_orb2dot, prstar2dot, Heff2dot, E2dot;
  double Fphi_Newt, sum_k;
  double oneby_EHeff_orb, oneby_E;
  double rdot2;
  
  /* Computing metric, centrifugal radius and ggm functions*/
  if(usespins) {
    eob_metric_s(r, prstar, dyn, &A, &B, &dA, &d2A, &dB, &d2B, &Q, &dQ, &dQ_dprstar, &d2Q, &ddQ_drdprstar, &d2Q_dprstar2, &d3Q_dr2dprstar, &d3Q_drdprstar2, &d3Q_dprstar3);
    eob_dyn_s_get_rc(r, nu, a1, a2, aK2, C_Q1, C_Q2, C_Oct1, C_Oct2, C_Hex1, C_Hex2, usetidal, &rc, &drc_dr, &d2rc_dr2);
    eob_dyn_s_GS(r, rc, drc_dr, d2rc_dr2, aK2, prstar, 0.0, nu, chi1, chi2, X1, X2, c3, ggm);
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
    eob_ham_s(nu, r, rc, drc_dr, d2rc_dr2, pphi, prstar, S, Sstar, chi1, chi2, X1, X2, aK2, c3, A, dA, d2A, Q, dQ, dQ_dprstar, d2Q, d2Q_dprstar2, &H, &Heff, &Heff_orb, &dHeff_dr, NULL, &dHeff_dpphi, NULL, NULL);
    E = nu*H;
  } else {
    eob_ham(nu, r, pphi, prstar, A, dA, Q, dQ, dQ_dprstar,
	    &H, &Heff_orb, &dHeff_dr, NULL, &dHeff_dpphi);
    Heff = Heff_orb;
    E = nu*H;
  }
  
  Adot = dA*rdot;
  sqrtAbyBdot = dsqrtAbyB_dr*rdot;
    
  dAbyrc2  = (dA*uc2 - 2.*A*uc3*drc_dr);
  d2Abyrc2 = d2A*uc2 - 4.*dA*uc3*drc_dr  + 6.*A*uc4*SQ(drc_dr) - 2.*A*uc3*d2rc_dr2;

  /* shorthands */
  EHeff_orb       = E*Heff_orb;
  oneby_EHeff_orb = 1.0/EHeff_orb;
  oneby_E         = 1.0/E;
  rdot2           = SQ(rdot);
  
  /* Begin iteration */
  int iter = 2; // Hard-fixed to 2: seems to be enough
    
  for(int n = 1; n <= iter; n++){
    
    prstardot = - sqrtAbyB/E*(pphi*dG_dr + 1./(2.*Heff_orb)
			        *(dA*(1. + pphi2*uc2 + Q)
	            + A*(-2.*uc3*pphi2*drc_dr + dQ))) + sqrtAbyB*Fr;

    Edot = nu*(rdot*Fr + Omg*Fphi);
    Heffdot  = 1./nu*E*Edot;
    HSOdot   = Fphi*G + pphi*(dG_dr*rdot + dG_dprstar*prstardot);
    Heff_orbdot  = Heffdot - HSOdot;
    EHeff_orbdot = Edot*Heff_orb + E*Heff_orbdot;

    Omgdot = oneby_EHeff_orb*dAbyrc2*rdot*pphi + oneby_E*(dG_dr*rdot + dG_dprstar*prstardot) 
           + oneby_EHeff_orb*A*uc2*(Fphi - pphi*oneby_EHeff_orb*EHeff_orbdot)
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
    double der3 = 2.*dA*pphi*uc2 - 4.*A*pphi*uc3*drc_dr;

    prstar2dot = dsqrtAbyB_dr*oosqrtAbyB*rdot*prstardot 
             + 0.5*sqrtAbyB*oneby_EHeff_orb*(oneby_EHeff_orb*EHeff_orbdot*(dA*(1. + pphi2*uc2 + Q) + A*(dQ - 2.*uc3*pphi2*drc_dr))
             - (der1*rdot + der2*prstardot + der3*Fphi)) 
             - sqrtAbyB*oneby_E*(pphi*(rdot*d2G_dr2 + prstardot*d2G_dr_dprstar) + dG_dr*(Fphi - pphi*Edot*oneby_E))
             + sqrtAbyB*Frdot;
    
    HSO2dot = pphi2dot*G + 2.*Fphi*(dG_dr*rdot + dG_dprstar*prstardot)
            + pphi*(d2G_dr2*rdot2 + dG_dr*r2dot + 2.*d2G_dr_dprstar*rdot*prstardot
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
    
    hatflm_NC[1] = Fphi_NewtPref(r, Omg, rdot, r2dot, r3dot, Omgdot, Omg2dot);
    sum_k = 0.;
    for (int k = 0; k < KMAX; k++) sum_k += Flm[k] * hatflm_NC[k];
    sum_k = sum_k/FlmNewt[1];

    Fphi  = Fphi_lo * sum_k;
    Fr    = eob_flx_Fr(r, prstar, pphi, dyn, Fphi);
    Fphi  = Fphi + Fphi_H;
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
