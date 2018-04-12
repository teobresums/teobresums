function [ output_args ] = QuickDQomg
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

CQ1  = 8.3958815463178116;
CQ2  = 8.3958815463178116;
%%{
chi1 = 0.1;
chi2 = 0.1;
%}
%{
chi1 = 0.1;
chi2 = 0.1;
%}
nu   = 0.25;
X1   = 1/2*(1+sqrt(1-4*nu));
X2   = 1-X1;
at1  = X1.*chi1;
at2  = X2.*chi2;


%-----------------------
% No quadrupole-monopole
%-----------------------
y     = load('h22_CQ0_s01.dat');
t     = y(:,1);
phi   = y(:,3);
omg0  = Deriv4(phi,t,4);
domg  = Deriv4(omg0,t,4);
Qomg0 = omg0.^2./domg;

%--------------------
% Quadrupole-monopole
%--------------------
y     = load('h22_CQ_s01.dat');
t     = y(:,1);
phi   = y(:,3);
omg   = Deriv4(phi,t,4);
domg  = Deriv4(omg,t,4);
Qomg  = omg.^2./domg;

figure
plot(omg0,Qomg0,'k');
hold on;
plot(omg,Qomg,'--r');

%---------------
% ancillary grid
%---------------
omgL = max(omg0(1),omg(1));
omgR = min(omg0(end),omg(end));

w   = omgL:1:omgR;
Qw0 = spline(omg0,Qomg0,w);
Qw  = spline(omg,Qomg,w);
DQw = -Qw0+Qw;

%-----------------
% analytical Qw_MQ
%-----------------
Msun = 4.925491025543575903411922162094833998e-6;
M    = (1.35+1.35)*Msun;


QwPN = -25/(48*nu).*(at1.^2*CQ1 + at2.^2*CQ2).*(M*w/2).^(-1/3);

figure('Name','EOB/PN comparison');
plot(w/(2*pi),DQw,'r');
hold on;
plot(w/(2*pi),QwPN,'k');
legend('EOB','PN: Poisson 1998','Location','NorthWest');
xlabel('$f$ [Hz]','Fontsize',20,'Interpreter','Latex');
ylabel('$Q_\omega^{\rm QM}$','Fontsize',20,'Interpreter','Latex');


end

