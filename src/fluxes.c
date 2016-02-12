/**
 * \file fluxes.c
 *
 * \brief Routines to construct the EOB fluxes.
 *
 * \author S.Bernuzzi
 */

#include "eob.h"

double **eob_flx_Newt_coefs (double nu)
{
  eobcode_kmultipolar_index_giveall;
  
  double nu2 = SQ(nu);
  double nu3 = nu2*nu;
  double nu4 = nu3*nu;

  double sp2 = 1-4*nu;
  double sp3 = pow((1-3*nu), 2);
  double sp4 = sp2*pow((1-2*nu),2);
  double sp5 = pow((1-5*nu+5*nu2),2);
  double sp6 = sp2*pow((3*nu2-4*nu +1),2);
  double sp7 = pow((1 - 7*nu + 14*nu2 - 7*nu3),2);

  if (par_get_b("tides")) {
    sp2 = 1.; 
  }

  double *FlmNewtCoefs;
  FlmNewtCoefs = *(double) calloc( eobcode_kmax * (sizeof(double)) );  

  FlmNewtCoefs[k21] = 8./45. * sp2;
  FlmNewtCoefs[k22] = 32./5.;                   
  
  FlmNewtCoefs[k31] = 1./1260. * sp2;
  FlmNewtCoefs[k32] = 32./63. * sp3;
  FlmNewtCoefs[k33] = 243./28. * sp2;

  FlmNewtCoefs[k41] = 1./44100. * sp4;
  FlmNewtCoefs[k42] = 32./3969. * sp3;
  FlmNewtCoefs[k43] = 729./700. * sp4;
  FlmNewtCoefs[k44] = 8192./567. * sp3;
  
  FlmNewtCoefs[k51] = 1./19958400. * sp4;
  FlmNewtCoefs[k52] = 256./400950. * sp5;
  FlmNewtCoefs[k53] = 2187./70400. * sp4;
  FlmNewtCoefs[k54] = 131072./66825. * sp5;
  FlmNewtCoefs[k55] = 1953125./76032. * sp4;
  
  FlmNewtCoefs[k61] = 1./1123782660. * sp6;
  FlmNewtCoefs[k62] = 128./28667925. * sp5;
  FlmNewtCoefs[k63] = 59049./15415400. * sp6;
  FlmNewtCoefs[k64] = 4194304./47779875. * sp5;
  FlmNewtCoefs[k65] = 48828125./13621608. * sp6;
  FlmNewtCoefs[k66] = 839808./17875. * sp5;

  FlmNewtCoefs[k71] = (sp2*pow((1 - 4*nu + 3*nu2),2))/9.3498717312e11;
  FlmNewtCoefs[k72] = (32*sp7)/1.35270135e8;
  FlmNewtCoefs[k73] = (1594323*sp2*pow((1 - 4*nu + 3*nu2),2))/3.2064032e10;
  FlmNewtCoefs[k74] = (4194304*sp7)/3.07432125e8;
  FlmNewtCoefs[k75] = 1220703125./5.666588928e9 * sp6;
  FlmNewtCoefs[k76] = 5668704./875875. * sp7;
  FlmNewtCoefs[k77] = 96889010407./1111968000. * sp6;

  FlmNewtCoefs[k81] = (sp2*pow((1 - 6*nu + 10*nu2 - 4*nu3),2))/8.174459284992e13;
  FlmNewtCoefs[k82] = (32*sp7)/3.4493884425e10;
  FlmNewtCoefs[k83] = (177147*sp2*pow((1 - 6*nu + 10*nu2 - 4*nu3),2))/3.96428032e10;
  FlmNewtCoefs[k84] = (4194304*sp7)/1.5679038375e10;
  FlmNewtCoefs[k85] = (30517578125*sp2*pow((1 - 6*nu + 10*nu2 - 4*nu3),2))/8.00296713216e11;
  FlmNewtCoefs[k86] = (51018336*sp7)/1.04229125e8;
  FlmNewtCoefs[k87] = (4747561509943*sp2*pow((1 - 6*nu + 10*nu2 - 4*nu3),2))/4.083146496e11;
  FlmNewtCoefs[k88] = pow((7*nu3-14*nu2+7*nu-1),2)*274877906944./1688511825.;

  return FlmNewCoefs;  
}

double **eob_flx_Newt_coefs_nu0 ()
{
  eobcode_kmultipolar_index_giveall;
  
  double nu2 = SQ(nu);
  double nu3 = nu2*nu;
  double nu4 = nu3*nu;

  double sp2 = 1-4*nu;
  double sp3 = pow((1-3*nu), 2);
  double sp4 = sp2*pow((1-2*nu),2);
  double sp5 = pow((1-5*nu+5*nu2),2);
  double sp6 = sp2*pow((3*nu2-4*nu +1),2);
  double sp7 = pow((1 - 7*nu + 14*nu2 - 7*nu3),2);

  if (par_get_b("tides")) {
    sp2 = 1.; 
  }

  double *FlmNewtCoefs;
  FlmNewtCoefs = *(double) calloc( eobcode_kmax * (sizeof(double)) );  

  FlmNewtCoefs[k21] = 8.0/45.0;
  FlmNewtCoefs[k22] = 32.0/5.0;                     
  
  FlmNewtCoefs[k31] = 1./1260.0;
  FlmNewtCoefs[k32] = 32.0/63.0;
  FlmNewtCoefs[k33] = 243.0/28.0;

  FlmNewtCoefs[k41] = 1./ 44100.0;
  FlmNewtCoefs[k42] = 32.0/3969.0;
  FlmNewtCoefs[k43] = 729.0/700.0;
  FlmNewtCoefs[k44] = 8192.0/567.0;
  
  FlmNewtCoefs[k51] = 1./ 19958400.0;
  FlmNewtCoefs[k52] = 256.0/400950.0; 
  FlmNewtCoefs[k53] = 2187.0/70400.0;
  FlmNewtCoefs[k54] = 131072.0/66825.0;
  FlmNewtCoefs[k55] = 1953125.0/76032.0;
  
  FlmNewtCoefs[k61] = 1./1123782660.0;
  FlmNewtCoefs[k62] = 128.0/28667925.0;
  FlmNewtCoefs[k63] = 59049.0/15415400.0;
  FlmNewtCoefs[k64] = 4194304.0/47779875.0;
  FlmNewtCoefs[k65] = 48828125.0/13621608.0;
  FlmNewtCoefs[k66] = 839808.0/17875.0;

  FlmNewtCoefs[k71] = 1./9.3498717312e11;
  FlmNewtCoefs[k72] =  32.0/1.35270135e8;
  FlmNewtCoefs[k73] = 1594323.0/3.2064032e10;
  FlmNewtCoefs[k74] = 4194304.0/3.07432125e8;
  FlmNewtCoefs[k75] = 1220703125.0/5.666588928e9;
  FlmNewtCoefs[k76] = 5668704.0/875875.0;
  FlmNewtCoefs[k77] = 96889010407.0/1111968000.0;

  FlmNewtCoefs[k81] = 1./8.174459284992e13;
  FlmNewtCoefs[k82] = 32.0/3.4493884425e10;
  FlmNewtCoefs[k83] = 177147.0/3.96428032e10;
  FlmNewtCoefs[k84] = 4194304.0/1.5679038375e10;
  FlmNewtCoefs[k85] = 30517578125.0/8.00296713216e11;
  FlmNewtCoefs[k86] = 51018336.0/1.04229125e8;
  FlmNewtCoefs[k87] = 4747561509943.0/4.083146496e11;
  FlmNewtCoefs[k88] = 274877906944.0/1688511825.0;

  return FlmNewCoefs;  
}

void eob_flx_FlmNewt_n (double x, 
			tlmk *lmk, 
			double *FlmNewtCoefs,
			double *FlmNewtonu2_n)
{
  eobcode_kmultipolar_index_lmkgive;
    
  double x5 = pow(x,5);
  double x6 = x5*x;
  double x7 = x6*x;
  double x8 = x7*x;
  double x9 = x8*x;
  double x10 = x9*x;
  double x11 = x10*x;
  double x12 = x11*x;
  
  *(FlmNewtonu2_n + k21) = FlmNewtCoefs[k21] * x6;
  *(FlmNewtonu2_n + k22) = FlmNewtCoefs[k22] * x5;
  
  *(FlmNewtonu2_n + k31) = FlmNewtCoefs[k31] * x6;
  *(FlmNewtonu2_n + k32) = FlmNewtCoefs[k32] * x7;
  *(FlmNewtonu2_n + k33) = FlmNewtCoefs[k33] * x6;
  
  *(FlmNewtonu2_n + k41) = FlmNewtCoefs[k41] * x8;
  *(FlmNewtonu2_n + k42) = FlmNewtCoefs[k42] * x7;
  *(FlmNewtonu2_n + k43) = FlmNewtCoefs[k43] * x8;
  *(FlmNewtonu2_n + k44) = FlmNewtCoefs[k44] * x7;
  
  *(FlmNewtonu2_n + k51) = FlmNewtCoefs[k51] * x8;
  *(FlmNewtonu2_n + k52) = FlmNewtCoefs[k52] * x9;
  *(FlmNewtonu2_n + k53) = FlmNewtCoefs[k53] * x8;
  *(FlmNewtonu2_n + k54) = FlmNewtCoefs[k54] * x9;
  *(FlmNewtonu2_n + k55) = FlmNewtCoefs[k55] * x8;
  
  *(FlmNewtonu2_n + k61) = FlmNewtCoefs[k61] * x10;
  *(FlmNewtonu2_n + k62) = FlmNewtCoefs[k62] * x9;
  *(FlmNewtonu2_n + k63) = FlmNewtCoefs[k63] * x10;
  *(FlmNewtonu2_n + k64) = FlmNewtCoefs[k64] * x9;
  *(FlmNewtonu2_n + k65) = FlmNewtCoefs[k65] * x10;
  *(FlmNewtonu2_n + k66) = FlmNewtCoefs[k66] * x9;
  
  *(FlmNewtonu2_n + k71) = FlmNewtCoefs[k71] * x10;
  *(FlmNewtonu2_n + k72) = FlmNewtCoefs[k72] * x11;
  *(FlmNewtonu2_n + k73) = FlmNewtCoefs[k73] * x10;
  *(FlmNewtonu2_n + k74) = FlmNewtCoefs[k74] * x11;
  *(FlmNewtonu2_n + k75) = FlmNewtCoefs[k75] * x10;
  *(FlmNewtonu2_n + k76) = FlmNewtCoefs[k76] * x11;
  *(FlmNewtonu2_n + k77) = FlmNewtCoefs[k77] * x10;

  *(FlmNewtonu2_n + k81) = FlmNewtCoefs[k81] * x12;
  *(FlmNewtonu2_n + k82) = FlmNewtCoefs[k82] * x11;
  *(FlmNewtonu2_n + k83) = FlmNewtCoefs[k83] * x12;
  *(FlmNewtonu2_n + k84) = FlmNewtCoefs[k84] * x11;
  *(FlmNewtonu2_n + k85) = FlmNewtCoefs[k85] * x12;
  *(FlmNewtonu2_n + k86) = FlmNewtCoefs[k86] * x11;
  *(FlmNewtonu2_n + k87) = FlmNewtCoefs[k87] * x12;
  *(FlmNewtonu2_n + k88) = FlmNewtCoefs[k88] * x11;

}

void eob_flx_FlmNewt (double *x, 
		      tlmk *lmk, int Nt, 
		      double *FlmNewtCoefs,
		      double *FlmNewtonu2)
{
  const int Nk = lmk->Nk;
  int n;
  for (n = 0; n < Nt; n++) {
    FlmNewtonu2_n = FlmNewtonu2 + n*Nk;
    eob_flx_FlmNewt_n (x[n], lmk, FlmNewtCoefs,FlmNewtonu2_n);  
  }

}

void eob_flx_absTlm_n (double Omega, double E, 
		       tlmk *lmk,
		       double *absTlm_n) 
{
  double *k2, *oofact2;
  double fl,hatk,x2, y,prod;   
  int *l = lmk->l;
  int *m = lmk->m;
  const int Nk = lmk->Nk;
  int k;

  k2 = calloc( Nk+1, sizeof(double) );
  oofact2 = calloc( Nk, sizeof(double) );
  
  /* precompute some stuff */
  for( k = 0; k < Nk; k++ ) {
    fl = fact( (double)l[k] );
    oofact2[k] = 1./( SQ(fl) );    
    k2[k]      = (double)(SQ(k)); 
  }
  k2[eobcode_kmax] = (double)SQ(eobcode_kmax);
  
  for( k = 0; k < Nk; k++ ) {    
    
    hatk = (double)m[k]*E_n*Omega_n;            
    x2   = 4.*(SQ(hatk));		
    
    prod = 1.;
    for( j = 1; j <= l[k]; j++ ) {       
      prod *= (k2[j]+x2); 
    }
    
    y = 4.*Pi*hatk;
    y = y/( 1. - exp(-y) );
    
    /* put the three pieces together */
    *(absTlm_n + k) = sqrt( oofact2[k] * y * prod ); 
    
  }
  
  free(k2);
  free(oofact2);
}

void eob_flx_absTlm (double *Omega,  double *E, 
		     tlmk *lmk, int Nt, 
		     double *absTlm) 
{
  double *Omega_n, *E_n, *absTlm_n; 
  double *k2, *oofact2;
  double fl,hatk,x2, y,prod;   
  int *l = lmk->l;
  int *m = lmk->m;
  const int Nk = lmk->Nk;
  int k;
  
  k2 = calloc( Nk+1, sizeof(double) );
  oofact2 = calloc( Nk, sizeof(double) );
  
  /* precompute some stuff */
  for( k = 0; k < Nk; k++ ) {
    fl = fact( (double)l[k] );
    oofact2[k] = 1./( SQ(fl) );    
    k2[k]      = (double)(SQ(k)); 
  }
  k2[eobcode_kmax] = (double)SQ(eobcode_kmax);
  
  int n;
  for( n = 0; n < Nt; n++ ) {      
    absTlm_n = absTlm + n*Nk; 
    Omega_n = Omega[n];
    E_n = E[n];
    
    for( k = 0; k < Nk; k++ ) {    
      
      hatk = (double)m[k]*E_n*Omega_n;            
      x2   = 4.*(SQ(hatk));		
      
      prod = 1.;
      for( j = 1; j <= l[k]; j++ ) {       
	prod *= (k2[j]+x2); 
      }

      y = 4.*Pi*hatk;
      y = y/( 1. - exp(-y) );
      
      /* put the three pieces together */
      *(absTlm_n + k) = sqrt( oofact2[k] * y * prod ); 

    }

  }
  
  free(k2);
  free(oofact2);
}

void eob_flx_abshlmNQC_n (double *NQCn_n,
			  tlmk *lmk, 
			  double **NQCCoefs,
			  double *abshlmNQC_n)
{
  double amp;
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
    *(hlmNQC_n + k) = amp;
  }
}

void eob_flx_abshlmNQC (double *NQCn,
			tlmk *lmk, int Nt, 
			double **NQCCoefs,
			domplex *abshlmNQC)
{
  double *NQCn_n;
  double *hlmNQC_n;
  const int Nk = lmk->Nk;
  int n;
  for (n = 0; n < Nt; n++) {
    hlmNQC_n = hlmNQC + n*Nk;
    NQCn_n = NQCn + n*eobcode_NQC; 
    eob_flx_abshlmNQC_n (NQCn_n, lmk, NQCCoefs, abshlmNQC_n);
  }
}












/*

  % Compute the RR
  Omega   = dydt(1);
  [Fphi, Frstar] = EOBRRDINddotr(nu, Omega, y, EOBopt, EOBmet, EOBHam, ddotr);
  
*/





 /*
function [hatF Flm FNewtlm flm] = ...
    EOBFluxDIN(x,Omega,E,Heff,jhat,nu,r,prstar,ddotr, EOBopt)

%EOBFluxDIN Computes the resummed Newton-normalized GW energy flux.
%
%   [hatF Flm FNewtlm flm] = ...
%   EOBFluxDIN(x,Omega,E,Heff,jhat,nu,r,prstar,ddotr, EOBopt) 
%
%   Reference(s)
%


L                 = EOBopt.L;
M                 = EOBopt.M;
LM2K              = EOBopt.LM2K;

rholmarg          = EOBopt.rholmarg;
f22_resum         = EOBopt.resumf22;

particle_dynamics = strcmp(EOBopt.Dynamics,'particle');
use_NQC           = strcmp(EOBopt.UseNQC,'yes');
use_tidal         = strcmp(EOBopt.Tidal,'yes');

if use_NQC
    NQC = EOBopt.NQC;
end
if use_tidal
    Topt = EOBopt.TidalStruct;     
end


% Compute amplitude of resummed waveform

% Argument for the modulus of hhatlm
% Note this has to be used only in the hhatlm modulus calculations, 
% everywhere else use x !
  switch rholmarg    
    case 'v_omg'
        xarg = Omega.^(2/3);
    case 'v_phi'
        xarg = x;
    case '1_r'
        xarg = 1./r;
    otherwise
        error('unknown option %s for rholmarg',rholmarg)
end

% do not use NQC
hnqclm = ones(length(L),length(Omega));

if particle_dynamics
    
    % Compute flm = (rho_lm)^ell
    flm    = EOBflm0(xarg,EOBopt);        
    
    % do not use NQC by default
    %hnqclm = ones(length(L),length(Omega));
   
else

    % Compute f_lm = (rho_lm)^ell    
    flm = EOBflm(xarg,nu,EOBopt);
        
    % Further resummation of rho_22 or f_22 (if needed)
    if strcmp(f22_resum,'pade23')        
        flm(LM2K(2,2),:) = EOBrho22Pade23(xarg,nu).^2;   
    elseif strcmp(f22_resum,'oldpade')   
        flm(LM2K(2,2),:) = EOBf22Pade(xarg,nu);    
    end
    
    if use_NQC
        % Compute NQC correction to the modulus of the (l,m) waveform        
        [n1,n2,n3]= EOBNQCn(r,prstar, Omega,ddotr);    
        hnqclm = (1+NQC.a{1}*n1 + NQC.a{2}*n2 + NQC.a{3}*n3);
    end
    
end


% Compute The PN correction to the waveform amplitude
% Note that jhat must be put to 1 when 1 has factored
% the full angular momentum J Omega in the Newtonian part 
% Beware this is the absolute value of the \hat{h}_lm
%MTlm  = EOBModTail( L, M, E.*Omega );
% MEX-file version ( - O(100) )
MTlm  = EOBModTailm( L, M, E.*Omega );


% Calculate prefactor (extend vectors to matrices) 
kmax   = length(L);
nx     = length(xarg);
lmeven = (~mod(L+M,2))';
lmodd  = ( mod(L+M,2))';
Heff   = reshape(Heff,1,nx); % do not resize back, if not needed below
jhat   = reshape(jhat,1,nx);              

prefact = lmeven*Heff + lmodd*jhat; 


% Compute modulus of hhat_lm 
Modhhatlm = prefact .* MTlm .* flm;


% Tidal contribution
hlmTidal = zeros(kmax,nx);

if use_tidal            
    
    % Include in flux  the correct tidal contribution for some multipoles 
    hlmTidal               = EOBhlmTidal( xarg,nu, Topt,EOBopt ); 
    
    % Nontidal part: re-insert the correct nu-dependence removed in the
    % Newtonian part when the tides are activated.
    sqrt_one_4nu = sqrt(1-4*nu);
    Modhhatlm(LM2K(2,1),:) = sqrt_one_4nu * Modhhatlm(LM2K(2,1),:);
    Modhhatlm(LM2K(3,1),:) = sqrt_one_4nu * Modhhatlm(LM2K(3,1),:);
    Modhhatlm(LM2K(3,3),:) = sqrt_one_4nu * Modhhatlm(LM2K(3,3),:);    

end


% Include NQC corrections (and tidal contribute, if needed)
Modhhatlm = Modhhatlm .*hnqclm + MTlm.*hlmTidal;    


% Compute Newtonian flux
FNewtlm = EOBflmNewt(x,nu,EOBopt);
FNewt22 = FNewtlm( LM2K(2,2) , : );


% Total flux multipoles
Flm = Modhhatlm.^2 .* FNewtlm; 


% Sum over multipoles and normalize to the 22 Newtonian multipole
hatF = sum(Flm)./FNewt22;

  */


void eob_flx_DIN (tdyn *dyn, )
{

}

void eob_flx_DIN_nu0 (tdyn *dyn, )
{

}

