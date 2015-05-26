/**
 * \file waveforms.c
 *
 * \brief Routines to construct the EOB multipolar waveform.
 *
 * \author S.Bernuzzi
 */

#include "eob.h"

domplex **eob_wav_Newt_coefs (double nu)
{
  eobcode_kmultipolar_index_giveall;

  double nu2 = SQ(nu);
  double nu3 = nu2*nu;

  double p1 = 1.;
  double p2 = sqrt(1-4*nu);
  double p3 = (3*nu-1);
  double p4 = (2*nu-1)*sqrt(1-4*nu);
  double p5 = 1-5*nu+5*nu2;
  double p6 = (1-4*nu+3*nu2)*sqrt(1-4*nu);
  double p7 = 7*nu3 - 14*nu2 + 7*nu -1;

  if (par_get_b("tides")) {
    p2 = 1.; 
  }

  domplex *hlmNewtCoefs;
  hlmNewtCoefs = *(domplex) calloc( eobcode_kmax * (sizeof(domplex)) );

  hlmNewtCoefs[k21] = -8./3.*I*sqrt(Pi/5.) * p2; 
  hlmNewtCoefs[k22] = -8.*sqrt(Pi/5.) * p1;

  hlmNewtCoefs[k31] = -1./3.*I*sqrt(2.*Pi/35.) * p2;
  hlmNewtCoefs[k32] = 8./3.*sqrt(Pi/7.) * p3;
  hlmNewtCoefs[k33] = 3.*I*sqrt(6*Pi/7.) * p2;
  
  hlmNewtCoefs[k41] = 1./105.*I*sqrt(2*Pi) * p4;
  hlmNewtCoefs[k42] = 8./63.*sqrt(Pi) * p3;
  hlmNewtCoefs[k43] = -9./5.*I*sqrt(2.*Pi/7.) * p4;
  hlmNewtCoefs[k44] = -64./9.*sqrt(Pi/7.) * p3;

  hlmNewtCoefs[k51] = I/180.*sqrt(Pi/77.) * p4;
  hlmNewtCoefs[k52] = -16./135.*sqrt(Pi/11.) * p5;
  hlmNewtCoefs[k53] = -9.*I/20.*sqrt(3.*Pi/22.) * p4;
  hlmNewtCoefs[k54] = 256./45.*sqrt(Pi/33.) * p5;
  hlmNewtCoefs[k55] = 125./12.*I*sqrt(5.*Pi/66.)* p4;

  hlmNewtCoefs[k61] = -I/2079.*sqrt(2*Pi/65.) * p6;
  hlmNewtCoefs[k62] = -16./1485.*sqrt(Pi/13.) * p5; 
  hlmNewtCoefs[k63] = I*81./385.*sqrt(Pi/13.) * p6; 
  hlmNewtCoefs[k64] = 1024./495.*sqrt(2*Pi/195.) * p5;
  hlmNewtCoefs[k65] = -I*625./63.*sqrt(5.*Pi/429.) * p6;
  hlmNewtCoefs[k66] = -432./5.*sqrt(Pi/715.) * p5;
  
  hlmNewtCoefs[k71] = -I/108108.*sqrt(Pi/10.) * p6;
  hlmNewtCoefs[k72] = 8./3003.*sqrt(Pi/15.) * p7;
  hlmNewtCoefs[k73] = I*243./20020.*sqrt(3.*Pi/10.) * p6;
  hlmNewtCoefs[k74] = -1024./1365.*sqrt(2.*Pi/165.) * p7;
  hlmNewtCoefs[k75] = -I*3125./3276.*sqrt(5.*Pi/66.) * p6;
  hlmNewtCoefs[k76] = 648./35.*sqrt(3.*Pi/715.) * p7;
  hlmNewtCoefs[k77] = 16807./180.*I*sqrt(7.*Pi/4290.) * p6;
  
  hlmNewtCoefs[k88] = -131072./315.*sqrt(2.*Pi/17017.) * p7;

  return hlmNewtCoefs;  
}

void eob_wav_hlmNewt_n (double r, double phi, double Omega, 
			tlmk *lmk, 
			domplex *hlmNewtCoefs,   
			domplex *hlmNewt_n)
{
  double vphi  = r*Omega;
  double vphi2 = SQ(vphi);
  double vphi3 = vphi2*vphi;
  double vphi4 = vphi3*vphi;
  double vphi5 = vphi4*vphi;
  double vphi6 = vphi5*vphi;
  double vphi7 = vphi6*vphi;
  double vphi8 = vphi7*vphi;
  
  *(hlmNewt_n + lmk->k21) = hlmNewtCoefs[lmk->k21] * vphi3;
  *(hlmNewt_n + lmk->k22) = hlmNewtCoefs[lmk->k22] * vphi2;

  *(hlmNewt_n + lmk->k31) = hlmNewtCoefs[lmk->k31] * vphi3;
  *(hlmNewt_n + lmk->k32) = hlmNewtCoefs[lmk->k32] * vphi4;
  *(hlmNewt_n + lmk->k33) = hlmNewtCoefs[lmk->k33] * vphi3;
  
  *(hlmNewt_n + lmk->k41) = hlmNewtCoefs[lmk->k41] * vphi5;
  *(hlmNewt_n + lmk->k42) = hlmNewtCoefs[lmk->k42] * vphi4;
  *(hlmNewt_n + lmk->k43) = hlmNewtCoefs[lmk->k43] * vphi5;
  *(hlmNewt_n + lmk->k44) = hlmNewtCoefs[lmk->k44] * vphi4;
  
  *(hlmNewt_n + lmk->k51) = hlmNewtCoefs[lmk->k51] * vphi5;
  *(hlmNewt_n + lmk->k52) = hlmNewtCoefs[lmk->k52] * vphi6;
  *(hlmNewt_n + lmk->k53) = hlmNewtCoefs[lmk->k53] * vphi5;
  *(hlmNewt_n + lmk->k54) = hlmNewtCoefs[lmk->k54] * vphi6;
  *(hlmNewt_n + lmk->k55) = hlmNewtCoefs[lmk->k55] * vphi5;
  
  *(hlmNewt_n + lmk->k61) = hlmNewtCoefs[lmk->k61] * vphi7;
  *(hlmNewt_n + lmk->k62) = hlmNewtCoefs[lmk->k62] * vphi6;
  *(hlmNewt_n + lmk->k63) = hlmNewtCoefs[lmk->k63] * vphi7;
  *(hlmNewt_n + lmk->k64) = hlmNewtCoefs[lmk->k64] * vphi6;
  *(hlmNewt_n + lmk->k65) = hlmNewtCoefs[lmk->k65] * vphi7;
  *(hlmNewt_n + lmk->k66) = hlmNewtCoefs[lmk->k66] * vphi6;
  
  *(hlmNewt_n + lmk->k71) = hlmNewtCoefs[lmk->k71] * vphi7;
  *(hlmNewt_n + lmk->k72) = hlmNewtCoefs[lmk->k72] * vphi8;
  *(hlmNewt_n + lmk->k73) = hlmNewtCoefs[lmk->k73] * vphi7;
  *(hlmNewt_n + lmk->k74) = hlmNewtCoefs[lmk->k74] * vphi8;
  *(hlmNewt_n + lmk->k75) = hlmNewtCoefs[lmk->k75] * vphi7;
  *(hlmNewt_n + lmk->k76) = hlmNewtCoefs[lmk->k76] * vphi8;
  *(hlmNewt_n + lmk->k77) = hlmNewtCoefs[lmk->k77] * vphi7;
  
  *(hlmNewt_n + lmk->k88) = hlmNewtCoefs[lmk->k88] * vphi8;

  /* multiply by phase factor */
  const int Nk = lmk->Nk;
  int k;
  for (k = 0; k < Nk; k++) {
    *(hlmNewt_n + k) *= cexp( - I * (double)m[k] * phi );
  }  
}

void eob_wav_hlmNewt (double *r, double *phi, double *Omega, 
		      tlmk *lmk, int Nt, 
		      domplex *hlmNewtCofes,   
		      domplex *hlmNewt)
{
  domplex *hlmNewt_n;
  int Nk = lmk->Nk;
  int n;
  for (n = 0; n < Nt; n++) {
    hlmNewtlm_n = hlmNewt + n*Nk;
    eob_wav_hlmNewt_n (r[n],phi[n],Omega[n], lmk, hlmNewtCoefs, hlmNewt_n);
  }
}

double **eob_wav_rholm_coefs (double nu)
{
  eobcode_kmultipolar_index_giveall;
  double eps = par_get_b("waveform.use_FI_terms");
  double nu2 = SQ(nu);
  double nu3 = nu2*nu;
  double nu4 = nu3*nu;
  double **rholmCoefs;
  *rholmCoefs = malloc( eobcode_kmax * (sizeof(double*)) );
  int k;
  for (k = 0; k < eobcode_kmax; k++)
    rholmCoefs[k] = calloc( eobcode_wav_rholm_PNcoefs, sizeof(double) );
  
  rholmCoefs[k21][0] = 1.; 
  rholmCoefs[k21][1] = (-1.0535714285714286 + 0.27380952380952384*nu);
  rholmCoefs[k21][2] = (-0.8327841553287982 - 0.7789824263038548*nu + 0.13116496598639457*nu2);
  rholmCoefs[k21][3] = 2.9192806270460925;
  rholmCoefs[k21][4] = -1.019047619047619;
  rholmCoefs[k21][5] = -1.28235780892213;
  rholmCoefs[k21][6] = 1.073639455782313;
  rholmCoefs[k21][7] = eps*(-3.8466571723355227);
  rholmCoefs[k21][8] = eps*(0.8486467106683944);
  
  rholmCoefs[k22][0] = 1.; 
  rholmCoefs[k22][1] = (-1.0238095238095237 + 0.6547619047619048*nu);
  rholmCoefs[k22][2] = (-1.94208238851096 - 1.5601379440665155*nu + 0.4625614134542706*nu2);
  rholmCoefs[k22][3] = (12.736034731834051 - 2.902228713904598*nu - 1.9301558466099282*nu2 + 0.2715020968103451*nu3);
  rholmCoefs[k22][4] = -4.076190476190476;
  rholmCoefs[k22][5] = -2.4172313935587004;
  rholmCoefs[k22][6] = + 4.173242630385488;
  rholmCoefs[k22][7] = -30.14143102836864;
  rholmCoefs[k22][8] = 7.916297736025627;

                     
  rholmCoefs[k31][0] = 1.;
  rholmCoefs[k31][1] = (-0.7222222222222222 - 0.2222222222222222*nu);
  rholmCoefs[k31][2] = (0.014169472502805836 - 0.9455667789001122*nu - 0.46520763187429853*nu2);
  rholmCoefs[k31][3] = 1.9098284139598072;
  rholmCoefs[k31][4] = -0.4126984126984127;
  rholmCoefs[k31][5] = 0.5368150316615179;
  rholmCoefs[k31][6] = 0.2980599647266314;
  rholmCoefs[k31][7] = eps*(1.4497991763035063);
  rholmCoefs[k31][8] = eps*(-0.0058477188106817735);

  rholmCoefs[k32][0] = 1.; 
  rholmCoefs[k32][1] = (0.003703703703703704*(328. - 1115.*nu + 320.*nu2))/(-1. + 3.*nu);
  rholmCoefs[k32][2] = (6.235191420376606e-7*(-1.444528e6 + 8.050045e6*nu - 4.725605e6*nu2 - 2.033896e7*nu3 + 3.08564e6*nu4))/pow(-1. + 3.*nu,2);
  rholmCoefs[k32][3] = (6.220997955214429);
  rholmCoefs[k32][4] = (-1.6507936507936507);
  rholmCoefs[k32][5] = eps*(-3.4527288879001268);
  rholmCoefs[k32][6] = eps*(2.005408583186361);
  
  rholmCoefs[k33][0] = 1.;
  rholmCoefs[k33][1] = (-1.1666666666666667 + 0.6666666666666666*nu);
  rholmCoefs[k33][2] = (-1.6967171717171716 - 1.8797979797979798*nu + 0.45151515151515154*nu2);
  rholmCoefs[k33][3] = (14.10891386831863); 
  rholmCoefs[k33][4] = (-3.7142857142857144); 
  rholmCoefs[k33][5] = (-6.723375314944128);
  rholmCoefs[k33][6] = (+4.333333333333333);
  rholmCoefs[k33][7] = eps*(-29.568699895427518);
  rholmCoefs[k33][8] = eps*(+6.302092352092352);


  rholmCoefs[k41][0] = 1.;
  rholmCoefs[k41][1] = (0.001893939393939394*(602. - 1385.*nu + 288.*nu2))/(-1. + 2.*nu);
  rholmCoefs[k41][2] = - 0.36778992787515513;
  rholmCoefs[k41][3] = (0.6981550175535535);
  rholmCoefs[k41][4] = (-0.2266955266955267);
  rholmCoefs[k41][5] = eps*(-0.7931524512893319);
  rholmCoefs[k41][6] = eps*(+0.2584672482399755); 

  rholmCoefs[k42][0] = 1.;
  rholmCoefs[k42][1] = (0.0007575757575757576*(1146. - 3530.*nu + 285.*nu2))/(-1. + 3.*nu);
  rholmCoefs[k42][2] = - (3.1534122443213353e-9*(1.14859044e8 - 2.95834536e8*nu - 1.204388696e9*nu2 + 3.04798116e9*nu3 + 3.79526805e8*nu4))/pow(-1. + 3.*nu,2); 
  rholmCoefs[k42][3] = 4.550378418934105e-12*(8.48238724511e11);
  rholmCoefs[k42][4] = 4.550378418934105e-12*(-1.9927619712e11);
  rholmCoefs[k42][5] = eps*(-0.6621921297263365);
  rholmCoefs[k42][6] = eps*(0.787251738160829);

  rholmCoefs[k43][0] = 1.;
  rholmCoefs[k43][1] = (0.005681818181818182*(222. - 547.*nu + 160.*nu2))/(-1. + 2.*nu);
  rholmCoefs[k43][2] = - 0.9783218202252293;
  rholmCoefs[k43][3] = eps*(8.519456157072423);
  rholmCoefs[k43][4] = eps*(-2.0402597402597404);
  rholmCoefs[k43][5] = eps*(-5.353216984886716);
  rholmCoefs[k43][6] = eps*(2.5735094451003544);

  rholmCoefs[k44][0] = 1.;
  rholmCoefs[k44][1] = (0.0007575757575757576*(1614. - 5870.*nu + 2625.*nu2))/(-1. + 3.*nu);
  rholmCoefs[k44][2] = (3.1534122443213353e-9*(-5.11573572e8 + 2.338945704e9*nu - 3.13857376e8*nu2 - 6.733146e9*nu3 + 1.252563795e9*nu4))/pow(-1. + 3.*nu,2); 
  rholmCoefs[k44][3] = (15.108111214795123);
  rholmCoefs[k44][4] = (-3.627128427128427);
  rholmCoefs[k44][5] = eps*(-8.857121657199649);
  rholmCoefs[k44][6] = eps*(4.434988849534304);
  

  rholmCoefs[k51][0] = 1.;
  rholmCoefs[k51][1] = (0.002564102564102564*(319. - 626.*nu + 8.*nu2))/(-1. + 2.*nu);
  rholmCoefs[k51][2] = - 0.1047896120973044;
  rholmCoefs[k51][3] = eps*(0.642701885362399);
  rholmCoefs[k51][4] = eps*(-0.14414918414918415);
  rholmCoefs[k51][5] = eps*(-0.07651588046467575);
  rholmCoefs[k51][6] = eps*(0.11790664036817883);

  rholmCoefs[k52][0] = 1.;
  rholmCoefs[k52][1] = (0.00007326007326007326*(-15828. + 84679.*nu - 104930.*nu2 + 21980.*nu3))/(1. - 5.*nu + 5.*nu2);
  rholmCoefs[k52][2] = eps*(-0.4629337197600934);
  rholmCoefs[k52][3] = eps*(2.354458371550237);
  rholmCoefs[k52][4] = eps*(-0.5765967365967366);

  rholmCoefs[k53][0] = 1.;
  rholmCoefs[k53][1] = (0.002564102564102564*(375. - 850.*nu + 176.*nu2))/(-1. + 2.*nu);
  rholmCoefs[k53][2] = - 0.5788010707241477;
  rholmCoefs[k53][3] = eps*(5.733973288504755);
  rholmCoefs[k53][4] = eps*(-1.2973426573426574);
  rholmCoefs[k53][5] = eps*(-1.9573287625526001);
  rholmCoefs[k53][6] = eps*(1.2474448628294783);

  rholmCoefs[k54][0] = 1.;
  rholmCoefs[k54][1] = (0.00007326007326007326*(-17448. + 96019.*nu - 127610.*nu2 + 33320.*nu3))/(1. - 5.*nu + 5.*nu2);
  rholmCoefs[k54][2] = eps*(-1.0442142414362194);
  rholmCoefs[k54][3] = eps*(10.252052781721588);
  rholmCoefs[k54][4] = eps*(-2.3063869463869464);
  
  rholmCoefs[k55][0] = 1.;
  rholmCoefs[k55][1] = (0.002564102564102564*(487. - 1298.*nu + 512.*nu2))/(-1. + 2.*nu);
  rholmCoefs[k55][2] = -1.5749727622804546;
  rholmCoefs[k55][3] = eps*(15.939827047208668);
  rholmCoefs[k55][4] = eps*(-3.6037296037296036);
  rholmCoefs[k55][5] = eps*(-10.272578060123237);
  rholmCoefs[k55][6] = eps*(4.500041838503377);
  

  rholmCoefs[k61][0] = 1.;
  rholmCoefs[k61][1] = (0.006944444444444444*(-161. + 694.*nu - 670.*nu2 + 124.*nu3))/(1. - 4.*nu + 3.*nu2);
  rholmCoefs[k61][2] = eps*(- 0.29175486850885135);
  rholmCoefs[k61][3] = eps*(0.21653486654395454);
  rholmCoefs[k61][4] = eps*(-0.10001110001110002);
    
  rholmCoefs[k62][0] = 1.;
  rholmCoefs[k62][1] = (0.011904761904761904*(-74. + 378.*nu - 413.*nu2 + 49.*nu3))/(1. - 5.*nu + 5.*nu2);
  rholmCoefs[k62][2] = eps*(-0.24797525070634313);
  rholmCoefs[k62][3] = eps*(1.7942694138754138);
  rholmCoefs[k62][4] = eps*(-0.40004440004440006);

  rholmCoefs[k63][0] = 1.;
  rholmCoefs[k63][1] = (0.006944444444444444*(-169. + 742.*nu - 750.*nu2 + 156.*nu3))/(1. - 4.*nu + 3.*nu2);
  rholmCoefs[k63][2] = eps*(-0.5605554442947213);
  rholmCoefs[k63][3] = eps*(4.002558222882566);
  rholmCoefs[k63][4] = eps*(-0.9000999000999002);

  rholmCoefs[k64][0] = 1.;
  rholmCoefs[k64][1] = (0.011904761904761904*(-86. + 462.*nu - 581.*nu2 + 133.*nu3))/(1. - 5.*nu + 5.*nu2);
  rholmCoefs[k64][2] = eps*(-0.7228451986855349);
  rholmCoefs[k64][3] = eps*(7.359388663371044);
  rholmCoefs[k64][4] = eps*(-1.6001776001776002);

  rholmCoefs[k65][0] = 1.;
  rholmCoefs[k65][1] = (0.006944444444444444*(-185. + 838.*nu - 910.*nu2 + 220.*nu3))/(1. - 4.*nu + 3.*nu2);
  rholmCoefs[k65][2] = eps*(- 1.0973940686333457);
  rholmCoefs[k65][3] = eps*(11.623366217471297);
  rholmCoefs[k65][4] = eps*(-2.5002775002775004);

  rholmCoefs[k66][0] = 1.;
  rholmCoefs[k66][1] = (0.011904761904761904*(-106. + 602.*nu - 861.*nu2 + 273.*nu3))/(1. - 5.*nu + 5.*nu2)
  rholmCoefs[k66][2] = eps*(- 1.5543111183867486);
  rholmCoefs[k66][3] = eps*(16.645950799433503);
  rholmCoefs[k66][4] = eps*(-3.6003996003996006);


  rholmCoefs[k71][0] = 1.;
  rholmCoefs[k71][1] = (0.0014005602240896359*(-618. + 2518.*nu - 2083.*nu2 + 228.*nu3))/(1. - 4.*nu + 3.*nu2);
  rholmCoefs[k71][2] = eps*( - 0.1508235111143767);
  rholmCoefs[k71][3] = eps*(0.2581280702019663);
  rholmCoefs[k71][4] = eps*(-0.07355557607658449);
  
  rholmCoefs[k72][0] = 1.;
  rholmCoefs[k72][1] = (0.00006669334400426837*(16832. - 123489.*nu + 273924.*nu2 - 190239.*nu3 + 32760.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3); 
  rholmCoefs[k72][2] = eps*(-0.351319484450667);

  rholmCoefs[k73][0] = 1.;
  rholmCoefs[k73][1] = (0.0014005602240896359*(-666. + 2806.*nu - 2563.*nu2 + 420.*nu3))/(1. - 4.*nu + 3.*nu2);
  rholmCoefs[k73][2] = eps*(- 0.37187416047628863);
  rholmCoefs[k73][3] = (3.0835293524055283);
  rholmCoefs[k73][4] = (-0.6620001846892604);

  rholmCoefs[k74][0] = 1.;
  rholmCoefs[k74][1] = (0.00006669334400426837*(17756. - 131805.*nu + 298872.*nu2 - 217959.*nu3 + 41076.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
  rholmCoefs[k74][2] = eps*(-0.6473746896670599);

  rholmCoefs[k75][0] = 1.;
  rholmCoefs[k75][1] = (0.0014005602240896359*(-762. + 3382.*nu - 3523.*nu2 + 804.*nu3))/(1. - 4.*nu + 3.*nu2);
  rholmCoefs[k75][2] = eps*(- 0.8269193364414116);
  rholmCoefs[k75][3] = eps*(8.750589067052443);
  rholmCoefs[k75][4] = eps*(-1.838889401914612);

  rholmCoefs[k76][0] = 1.;
  rholmCoefs[k76][1] = (0.0006002400960384153*(2144. - 16185.*nu + 37828.*nu2 - 29351.*nu3 + 6104.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
  rholmCoefs[k76][2] = eps*(-1.1403265020692532);

  rholmCoefs[k77][0] = 1.;
  rholmCoefs[k77][1] = (0.0014005602240896359*(-906. + 4246.*nu - 4963.*nu2 + 1380.*nu3))/(1. - 4.*nu + 3.*nu2);
  rholmCoefs[k77][2] = eps*(- 1.5418467934923434);
  rholmCoefs[k77][3] = eps*(17.255875091408523);
  rholmCoefs[k77][4] = eps*(-3.6042232277526396);


  rholmCoefs[k81][0] = 1.;
  rholmCoefs[k81][1] = (0.00005482456140350877*(20022. - 126451.*nu + 236922.*nu2 - 138430.*nu3 + 21640.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
  rholmCoefs[k81][2] = eps*(-0.26842133517043704);

  rholmCoefs[k82][0] = 1.;
  rholmCoefs[k82][1] = (0.0003654970760233918*(2462. - 17598.*nu + 37119.*nu2 - 22845.*nu3 + 3063.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3); 
  rholmCoefs[k82][2] = eps*(-0.2261796441029474);

  rholmCoefs[k83][0] = 1.;
  rholmCoefs[k83][1] = (0.00005482456140350877*(20598. - 131059.*nu + 249018.*nu2 - 149950.*nu3 + 24520.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
  rholmCoefs[k83][2] = eps*(-0.4196774909106648);

  rholmCoefs[k84][0] = 1.;
  rholmCoefs[k84][1] = (0.0003654970760233918*(2666. - 19434.*nu + 42627.*nu2 - 28965.*nu3 + 4899.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
  rholmCoefs[k84][2] = eps*(-0.47652059150068155);

  rholmCoefs[k85][0] = 1.;
  rholmCoefs[k85][1] = (0.00027412280701754384*(4350. - 28055.*nu + 54642.*nu2 - 34598.*nu3 + 6056.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
  rholmCoefs[k85][2] = eps*(-0.7220789990670207);

  rholmCoefs[k86][0] = 1.;
  rholmCoefs[k86][1] = (0.0010964912280701754*(1002. - 7498.*nu + 17269.*nu2 - 13055.*nu3 + 2653.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
  rholmCoefs[k86][2] = eps*(-0.9061610303170207);

  rholmCoefs[k87][0] = 1.;
  rholmCoefs[k87][1] = (0.00005482456140350877*(23478. - 154099.*nu + 309498.*nu2 - 207550.*nu3 + 38920.*nu4))/(-1. + 6.*nu - 10.*nu2 + 4.*nu3);
  rholmCoefs[k87][2] = eps*(-1.175404252991305);
  
  rholmCoefs[k88][0] = 1.;
  rholmCoefs[k88][1] = (0.0003654970760233918*(3482. - 26778.*nu + 64659.*nu2 - 53445.*nu3 + 12243.*nu4))/(-1. + 7.*nu - 14.*nu2 + 7.*nu3);
  rholmCoefs[k88][2] = eps*(-1.5337092502821381);


  return rholmCoefs;  
}

void eob_wav_rholm_n (double x, 
		      tlmk *lmk,
		      double **rholmCoefs,  
		      double *rholm_n)
{
  double x2 = SQ(x);
  double x3 = x2*x;
  double x4 = x3*x;
  double x5 = x4*x;

  double el[9];
  el[0] = 0.;
  el[1] = Eulerlog(x,1);
  el[2] = Eulerlog(x,2);
  el[3] = Eulerlog(x,3);
  el[4] = Eulerlog(x,4);
  el[5] = Eulerlog(x,5);
  el[6] = Eulerlog(x,6);
  el[7] = Eulerlog(x,7);
  el[8] = 0.;

  const int Nk = lmk->Nk;
  int k, m;
  for (k=0; k<Nk; k++) {
    m = lmk->m[k];
    *(rholm_n + k) = 
      (rholmCoefs[k][0] + 
       (rholmCoefs[k][1] *x + 
	(rholmCoefs[k][2] *x2 + 
	 ((rholmCoefs[k][3] + (rholmCoefs[k][4]*el[m])) *x3 + 
	  ((rholmCoefs[k][5] + (rholmCoefs[k][6]*el[m])) *x4 + 
	   ((rholmCoefs[k][7] + (rholmCoefs[k][8]*el[m])) *x5))))));
  }
}

double eob_wav_rho22P23(double nu, double x) 
{
  double nu2 = SQ(nu);

  double eulerlog22 = EulerGamma + 2.*log2 + 0.5*log(x);
  double f1 = -43/42 + 55*nu/84;
  double f2 = -20555/10584 - 33025/21168*nu + 19583/42336*nu2;
  double f3 =  1556919113/122245200 - 428/105*eulerlog22 - 48993925*nu/9779616 -6292061*nu2/3259872 + 10620745*nu^3/39118464 + 41*nu*SQ(Pi)/192;
  double f4 = -387216563023/160190110080 + 9202*eulerlog22/2205;
  double f5 = -16094530514677/533967033600 + 439877/55566*eulerlog22;

  double f12 = SQ(f1); 
  double f22 = SQ(f2); 
  double f23 = f22*f2; 
  double f32 = SQ(f3);
  double f33 = f32*f3;
  
  /* Coefficients of the (2,3) Padé approximant */ 
  double c1 = -f1;
  double c2 = f1 - f2/f1;
  double c3 = (f1*f3-f22)/(f1*(f12-f2));
  double c4 = -f1*(f23+f32+f12*f4-f2*(2.0*f1*f3+f4))/((f12-f2)*(f1*f3-f22));
  double c5 = -(f12-f2)*(-f33+2.0*f2*f3*f4-f1*f4*f4-f22*f5+f1*f3*f5)
    / ( (f1*f3-f22)*(f23+f32+f12*f4-f2*(2.0*f1*f3+f4) ) );
    
  /* Compute the (2,3) Padé approximant */  
  double y = 1. + c4*x/(1.+c5*x);
  y = 1. + c3*x/y;
  y = 1. + c2*x/y;
  y = 1. + c1*x/y;
  
  return 1./y;  
}

void eob_wav_rholm (double *x, 
		    tlmk *lmk, int Nt, 
		    double **rholmCoefs,  
		    double *rholm)
{
  double *rholm_n;
  int Nk = lmk->Nk;
  int n;
  for (n = 0; n < Nt; n++) {
    rholm_n = rholm + n*Nk;
    eob_wav_rholm_n (x[n], lmk, rholmCoefs, rholm_n);  
  }
}

double **eob_wav_deltalm_coefs (double nu)
{
  eobcode_kmultipolar_index_giveall;
  double nu2 = SQ(nu);
  double Pi2 = SQ(Pi);

  double **deltalmCoefs;
  *deltalmCoefs = malloc( eobcode_kmax * (sizeof(double*)) );
  int k;
  for (k = 0; k < eobcode_kmax; k++)
    deltalmCoefs[k] = calloc( eobcode_wav_deltalm_PNcoefs, sizeof(double) );

  deltalmCoefs[k21][0] = 2./3.; 
  deltalmCoefs[k21][1] = -493./4.*nu;
  deltalmCoefs[k21][2] = 107./105.*Pi;
  deltalmCoefs[k21][3] = 0.; 
  deltalmCoefs[k21][4] = 214./315.*Pi2 - 272./81.;
  
  deltalmCoefs[k22][0] = 7./3.;
  deltalmCoefs[k22][1] = -24*nu;
  deltalmCoefs[k22][2] = +428./105.*Pi;
  deltalmCoefs[k22][3] = 30995./1134*nu + 962./135*nu2;
  deltalmCoefs[k22][4] = 1712./315.*Pi2 - 2203./81.;
                     

  deltalmCoefs[k31][0] = 13./30.;
  deltalmCoefs[k31][1] = 17./10.;
  deltalmCoefs[k31][2] = 13./21.*Pi;
  deltalmCoefs[k31][3] = 0.;
  deltalmCoefs[k31][4] = 26./63.*Pi2 - 227827./81000;

  deltalmCoefs[k32][0] = (10+33*nu)/(15.*(1-3*nu));
  deltalmCoefs[k32][1] = 0.;
  deltalmCoefs[k32][2] = 52./21.*Pi; 
  deltalmCoefs[k32][3] = 0.;
  deltalmCoefs[k32][4] = (208./63.*Pi2 - 9112./405.);
  
  deltalmCoefs[k33][0] = 13./10.;
  deltalmCoefs[k33][1] = -80897./2430.*nu;
  deltalmCoefs[k33][2] = 39./7.*pi;
  deltalmCoefs[k33][3] = 0.;
  deltalmCoefs[k33][4] = 78./7.*Pi2 - 227827./3000;


  deltalmCoefs[k41][0] = (2+507*nu)/(10*(1-2*nu));
  deltalmCoefs[k41][1] = 0.;
  deltalmCoefs[k41][2] = 1571./3465.*Pi;
  deltalmCoefs[k41][3] = 0.;
  deltalmCoefs[k41][4] = 0.;

  deltalmCoefs[k42][0] = (7*(1+6*nu))/(15*(1-3*nu)); 
  deltalmCoefs[k42][1] = 0.;
  deltalmCoefs[k42][2] = 6284./3465.*Pi;
  deltalmCoefs[k42][3] = 0.;
  deltalmCoefs[k42][4] = (25136./10395*Pi2 - 6893./375.);

  deltalmCoefs[k43][0] = (486+4961*nu)/(810*(1-2*nu));
  deltalmCoefs[k43][1] = 0.;
  deltalmCoefs[k43][2] = 1571./385.*Pi;
  deltalmCoefs[k43][3] = 0.;
  deltalmCoefs[k43][4] = 0.;

  deltalmCoefs[k44][0] = (112+219*nu)/(120*(1-3*nu));
  deltalmCoefs[k44][1] = 0.;
  deltalmCoefs[k44][2] = 25136./3465.*Pi;
  deltalmCoefs[k44][3] = 0.;
  deltalmCoefs[k44][4] = (201088./10395.*Pi2 - 55144./375.);
  

  deltalmCoefs[k55][0] = (96875 + 857528*nu)/(131250.*(1-2*nu));


  return deltalmCoefs;  
}

void eob_wav_deltalm_n (double Heob, double Omega, 
			tlmk *lmk, 
			double *deltalmCoefs,
			double *deltalm_n, domplex *expdlm_n)
{
  double EOmg = Heob*Omega;
  double vy = pow( EOmg , oo3 );
  double vy3 = EOmg; //pow(vy,3);
  double vy5 = pow(vy,5);
  double vy6 = SQ(EOmg); //pow(vy,6);
  double vy7 = pow(vy,7);
  double vy9 = vy6*vy3; //pow(vy,9);
  double sum;
  int Nk = lmk->Nk;
  int k51 = lmk->k51;
  int k55 = lmk->k55;
  int k;
  for (k = 0; k < Nk; k++) { 
    if (k==k51) k = k55;  
    if (k>k55) break;
    sum = (deltalmCoefs[k][0] * vy3 +
	   (deltalmCoefs[k][1] * vy5 +
	    (deltalmCoefs[k][2] * vy6 +
	     (deltalmCoefs[k][3] * vy7 +
	      (deltalmCoefs[k][4] * vy9)))));
    *(deltalm_n + k) = sum;
    if (sum==0.) *(expdlm_n + k) = 1.;
    else         *(expdlm_n + k) = cexp( I * sum );
  }
}

void eob_wav_deltalm_resum_n (double Heob, double Omega, 
			      tlmk *lmk, 
			      double *deltalmCoefs,
			      double *deltalm_n, domplex *expdlm_n)
{
  int k21 = lmk->k21;
  int k22 = lmk->k22;
  int k31 = lmk->k31;
  int k32 = lmk->k32;
  int k33 = lmk->k33;
  int k41 = lmk->k41;
  int k43 = lmk->k42;
  int k43 = lmk->k43;
  int k44 = lmk->k44;
  int k55 = lmk->k55;

  double delta21LO = deltalmCoefs[k21][0];
  double delta22LO = deltalmCoefs[k22][0];
  double delta31LO = deltalmCoefs[k31][0];
  double delta33LO = deltalmCoefs[k33][0];

  double nu = - deltalmCoefs[k22][1]/24.; /* workaround, deltalmCoefs[k22][1] = -24*nu */
  double nu2 = SQ(nu); 
  double Pi2 = SQ(Pi);

  double vy = pow( Heob*Omega, oo3 );
  double vy2 = SQ(vy);
  double vy3 = vy2*vy;
  double vy5 = vy2*vy3;
  double vy6 = vy3*vy3;
  double vy7 = vy6*vy;
  double vy9 = vy6*vy3;

  /* Resummed multipoles */
  double sum21,sum22, sum31,sum32,sum33;
  double sum41,sum42,sum43,sum44, sum55;
  double num21,num22,num31,num33;
  double den21,den22,den31,den33;

  num21 = 69020*nu + 5992*Pi*vy;
  den21 = 5992*Pi*vy + 2456*nu*(28+493*nu*vy2);
  sum21 = delta21LO*num21/den21;
  
  num22 = (808920*nu*Pi*vy + 137388*Pi2*vy2 + 35*nu2*(136080 + (154975 - 1359276*nu)*vy2));
  den22 = (808920*nu*Pi*vy + 137388*Pi2*vy2 + 35*nu2*(136080 + (154975 + 40404*nu)*vy2));
  sum22 = delta22LO*num22/den22;
  
  num31 = 4641*nu + 1690*Pi*vy;
  den31 = num31 + 18207*nu2*vy2;
  sum31 = delta31LO*num31/den31;
  
  num33  = 1. + 94770*Pi*vy/(566279*nu);
  den33 = num + 80897*nu*vy2/3159;
  sum33 = delta33LO*num33/den33;
  
  /* Taylor multipoles */
  sum32 = (deltalmCoefs[k32][0] * vy3 +
	   (deltalmCoefs[k32][2] * vy6 +
	    deltalmCoefs[k32][4] * vy9));
  
  sum41 = deltalmCoefs[k41][0] * vy3
    + deltalmCoefs[k41][2] * vy6;
  
  sum42 = (deltalmCoefs[k42][0] * vy3 +
	   (deltalmCoefs[k42][2] * vy6 +
	    (deltalmCoefs[k42][4] * vy9)));
  
  sum43 = deltalmCoefs[k43][0] * vy3 
    + deltalmCoefs[k43][2] * vy6;
  
  sum44 = (deltalmCoefs[k44][0] * vy3 +
	   (deltalmCoefs[k44][4] * vy6 +
	    (deltalmCoefs[k44][4] * vy9)));
  
  sum55 = deltalmCoefs[k55][0] * vy3;
  
  /* store */
  *(deltalm_n + k21) = sum21; 
  *(expdlm_n + k21) = cexp( I * sum21 );
  
  *(deltalm_n + k22) = sum22; 
  *(expdlm_n + k22) = cexp( I * sum22 );
  
  *(deltalm_n + k31) = sum31; 
  *(expdlm_n + k31) = cexp( I * sum31 );
  
  *(deltalm_n + k32) = sum32; 
  *(expdlm_n + k32) = cexp( I * sum32 );
  
  *(deltalm_n + k33) = sum33; 
  *(expdlm_n + k33) = cexp( I * sum33 );
  
  *(deltalm_n + k41) = sum41;
  *(expdlm_n + k41) = cexp( I * sum41 );
  
  *(deltalm_n + k42) = sum42;
  *(expdlm_n + k42) = cexp( I * sum42 );
  
  *(deltalm_n + k43) = sum43;
  *(expdlm_n + k43) = cexp( I * sum43 );
  
  *(deltalm_n + k44) = sum44;
  *(expdlm_n + k44) = cexp( I * sum44 );
  
  *(deltalm_n + k55) = sum55;
  *(expdlm_n + k55) = cexp( I * sum55 );    
  
}

void eob_wav_deltalm (double *Heob, double *Omega, 
		      tlmk *lmk, int Nt, 
		      double *deltalmCoefs,
		      double *deltalm, domplex *expdlm)
{
  double *deltalm_n; 
  domplex *expdlm_n;
  const int Nk = lmk->Nk;
  int n;
  for (n = 0; n < Nt; n++) {
    deltalm_n = deltam + n*Nk;
    expdlm_n = expdlm + n*Nk;
    eob_wav_deltalm_n (Heob[n],Omega[n], lmk, deltalmCoefs, deltalm_n,expdlm_n);
  }
}

void eob_wav_deltalm_resum (double *Heob, double *Omega, 
			    tlmk *lmk, int Nt, 
			    double *deltalmCoefs,
			    double *deltalm, domplex *expdlm)
{
  double *deltalm_n; 
  domplex *expdlm_n;
  const int Nk = lmk->Nk;
  int n;
  for (n = 0; n < Nt; n++) {
    deltalm_n = deltam + n*Nk;
    expdlm_n = expdlm + n*Nk;
    eob_wav_deltalm_resum_n (Heob[n],Omega[n], lmk, deltalmCoefs, deltalm_n,expdlm_n);
  }
}

double** eob_wav_NQCCoefs (double nu)
{
  int k21 = lmk_lm2k(2,1);			
  int k22 = lmk_lm2k(2,2);			
  int k31 = lmk_lm2k(3,1);			
  int k32 = lmk_lm2k(3,2);			
  int k33 = lmk_lm2k(3,3);

  double nu2 = SQ(nu);

  double **NQCab;
  *NQCCab = malloc( eobcode_kmax * (sizeof(double*)) );
  int k;
  for (k = 0; k < eobcode_kmax; k++)
    NQCCab[k] = calloc( eobcode_NQC, sizeof(double) );

  NQCCab[k21][0] = +9.150277e-01*nu2 - 6.521591e-01*nu + 3.400123e-02;
  NQCCab[k21][1] = +1.903505e+00*nu2 + 4.978512e+00*nu - 7.105946e-01; 
  NQCCab[k21][2] = +2.694953e+00*nu2 - 3.160287e+00*nu + 8.650109e-01;
  NQCCab[k21][3] = -2.247953e+00*nu2 + 5.304017e-01*nu + 2.465944e-01;
  NQCCab[k21][4] = +5.172595e+01*nu2 - 2.168939e+01*nu + 3.161619e+00;
  NQCCab[k21][5] = -1.128189e+02*nu2 + 3.055883e+01*nu - 2.607237e-03;
  
  NQCCab[k22][0] = +2.160054e+00*nu2 - 1.093667e+00*nu + 7.926107e-02;
  NQCCab[k22][1] = -1.080736e+01*nu2 + 7.141956e+00*nu + 7.034821e-01;
  NQCCab[k22][2] = -2.766573e+00*nu2 - 1.768795e-01*nu + 1.012170e-01;
  NQCCab[k22][3] = -8.568345e-01*nu2 - 2.416766e-01*nu + 1.928716e-01;
  NQCCab[k22][4] = +9.638222e+00*nu2 - 7.645317e+00*nu + 3.731762e-01;
  NQCCab[k22][5] = -8.099188e+01*nu2 + 1.707519e+01*nu - 1.797447e+00;

  NQCCab[k32][0] = +3.167065e+01*nu2 - 1.030991e+01*nu + 6.843778e-01;
  NQCCab[k32][1] = -1.327286e+02*nu2 + 5.515306e+01*nu - 3.044906e+00;
  NQCCab[k32][2] = -1.973440e+01*nu2 + 1.829857e+00*nu + 3.039708e-01;
  NQCCab[k32][3] = -1.349680e+00*nu2 - 1.808295e+00*nu + 5.734632e-01;
  NQCCab[k32][4] = -4.486021e+00*nu2 - 7.796816e+00*nu + 3.253824e+00;
  NQCCab[k32][5] = +3.248901e+00*nu2 - 3.813272e+01*nu + 6.164847e+00;

  NQCCab[k33][0] = +2.679308e+00*nu2 - 1.279163e+00*nu + 1.456279e-01;
  NQCCab[k33][1] = -1.293179e+01*nu2 + 1.063361e+01*nu + 4.704239e-01;
  NQCCab[k33][2] = -8.931901e-01*nu2 - 2.922908e+00*nu + 5.078229e-01;
  NQCCab[k33][3] = -1.646809e+00*nu2 - 1.611205e-01*nu + 3.464108e-01;
  NQCCab[k33][4] = +6.959682e+00*nu2 - 7.595824e+00*nu + 7.089351e-02;
  NQCCab[k33][5] = -2.134623e+02*nu2 + 5.281861e+01*nu - 3.471809e+00;

 return NQCCab;

}

void eob_wav_hlmNQC_n (double *NQCn_n,
		       tlmk *lmk, 
		       double **NQCCoefs,
		       domplex *hlmNQC_n)
{
  double amp,phi;
  const int k33 = lmk->k33;
  double n1 = NQCn_n[0];
  double n2 = NQCn_n[1];
  double n3 = NQCn_n[2];
  double n4 = NQCn_n[3];
  double n5 = NQCn_n[4];
  double n6 = NQCn_n[5];
  const int Nk = lmk->Nk;
  int k;
  for (k = 0; k < Nk; k++) {
    if (k>k33) break;
    amp = (1+
	   (NQCCoefs[k][0] *n1 +
	    (NQCCoefs[k][1] *n2 +
	     (NQCCoefs[k][2] *n3))));
    phi = (NQCCoefs[k][3] *n4 +
	   (NQCCoefs[k][4] *n5 +
	    (NQCCoefs[k][5] *n6)));
    *(hlmNQC_n + k) = amp * cexp( I * phi );
  }
}

void eob_wav_hlmNQC (double *NQCn,
		     tlmk *lmk, int Nt, 
		     double **NQCCoefs,
		     domplex *hlmNQC)
{
  double *NQCn_n;
  domplex *hlmNQC_n;
  const int Nk = lmk->Nk;
  int n;
  for (n = 0; n < Nt; n++) {
    hlmNQC_n = hlmNQC + n*Nk;
    NQCn_n = NQCn + n*eobcode_NQC; 
    eob_wav_hlmNQC_n (NQCn_n, lmk, NQCCoefs, hlmNQC_n);
  }
}

void eob_wav_Tlm_n (double Omega, double E, 
		    tlmk *lmk, 
		    double r0, 
		    domplex *Tlm_n)
{
  domplex num,den;
  double lpo;
  int *l = lmk->l;
  int *m = lmk->m;
  const int Nk = lmk->Nk;
  int k;
  for (k = 0; k < Nk; k++) {
    lpo = (double)l[k] + 1.;
    k = (double)m[k] * Omega;
    khat = k * E;
    den = GammaComplex( lpo );
    num = GammaComplex( lpo - 2.*I*hatk );
    *(Tlm_n + k) = num/den * cexp( Pi*hatk + 2.*I*hatk*log(2.*k*r0) );
  }
}

void eob_wav_Tlm (double *Omega, double *E, 
		  tlmk *lmk, int Nt, 
		  double r0, 
		  domplex *Tlm)
{
  domplex *Tlm_n;
  const int Nk = lmk->Nk;
  int n;
  for (n = 0; n < Nt; n++) {
    Tlm_n = Tlm + n*Nk;
    eob_wav_Tlm_n (Omega[n], E[n], lmk, r0, Tlm_n);
  }
}

void eob_wav_Seff_n (double Heff, double jhat,
		     tlmk *lmk,
		     double *Seff_n)
{
  int *l = lmk->l;
  int *m = lmk->m;
  const int Nk = lmk->Nk;
  int k;
  for (k = 0; k < Nk; k++) {
    if ( (l[k]+m[k]) % 2 ) 
      *(Seff_n + k) = jhat_n;
    else
      *(Seff_n + k) = Heff_n;
  }
}

void eob_wav_Seff (double *Heff, double *jhat, 
		   tlmk *lmk, int Nt, 
		   double *Seff)
{
  double *Seff_n;
  const int Nk = lmk->Nk;
  int n;
  for (n = 0; n < Nt; n++) {
    Seff_n = Seff + n*Nk;
    eob_wav_Seff_n (Heff[n],jhat[n], lmk, Seff_n);
  }
}

void eob_wav_flm_n (double *rholm_n, 
		    tlmk *lmk, 
		    double *flm_n)
{
  int *l = lmk->l; 
  const int Nk = lmk->Nk;
  int k;
  for (k = 0; k < Nk; k++) {
    *(flm_n + k) = pow( *(rholm_n + k), l[k]  );
  }
}

void eob_wav_flm (double *rholm,
		  tlmk *lmk, int Nt, 
		  double *flm)
{
  double *rholm_n, *flm_n;
  const int Nk = lmk->Nk;
  int n;
  for (n = 0; n < Nt; n++) {
    rholm_n = rholm + n*Nk;
    flm_n   = flm   + n*Nk;
    eob_wav_flm_n (rholm_n, lmk, flm_n);
  }
}

double eob_wav_rholm_xarg (double r, double Omega, double x)
{
  double xarg;
  char *rholm_arg = par_get_s ("waveform.rholm_arg");
  if (iseq(rholm_arg,"v_omg")) {
    xarg = pow(Omega, to3);
  } 
  else if (iseq(rholm_arg,"v_phi")) {
    xarg = x;  
  }
  else if (iseq(rholm_arg,"1_r")) {
    xarg = 1./r;
  } else
    errorexits("unknown parameter waveform.rholm_arg value",rholm_arg);
  return xarg;
}


