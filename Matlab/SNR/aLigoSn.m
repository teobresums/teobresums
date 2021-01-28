function Sn=aLigoSn(f)

% aLigoSn function. Computes the analytical fit (done with a Matlab function in 
% Data_Analysis/DA_EOBNS/trunk/MatlabScripts/DA_LIGONoise.m)
% to the advanced LIGO noise curve named ZERO_DET_high_P.txt
%
% Sn = aLigoSn(f)
%
% where f is the frequency in Hz.
% 
% anagar, June 2012

% normalizations
f0  = 56.56;
hn0 = 3.5171e-23;

%
x  = f/f0;
x2 = x.^2;
x3 = x.^3;
lgx  = log(x);
lgx2 = lgx.^2;
lgx3 = lgx.^3;
lgx4 = lgx.^4;
lgx5 = lgx.^5;

% Fitting coefficients; high-freq part
hg1 = 0.41473797;
hg2 = 0.0927076;
hg3 = 0.44943703;
hg4 = 0.0034297409;
hg5 = -6.5736784*10^(-5);
hg6 = 0.0026988759;
% Fitting coefficients; low-freq part
lw1 = -2.0158256;
lw2 = -4.1577911;
lw3 = -1.0904057;
lw4 = -0.75793805;
lw5 = -0.34979753;
lw6 = -0.098312294;
lw7  = 3.3732744*10^(-6);


fShigh = hg1 + hg2*log(x) + hg3*x + hg4*x.*log(x) + hg5*x2 + hg6*x3;
fSlow  = exp(lw1 + lw2*lgx + lw3*lgx2 + lw4*lgx3 + lw5*lgx4 + lw6*lgx5 + lw7*x.^(-7));

fS = fSlow + fShigh;

% Go back to Sn
Sn = fS*hn0.^2./f;

%{
%figure('Name','Test');
%loglog(x,fS)

figure('Name','Test2');
loglog(f,Sn)
%}


return