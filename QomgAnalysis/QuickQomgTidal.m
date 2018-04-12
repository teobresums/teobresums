function [ output_args ] = QuickDQomg
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%CQ1  = 8.3958815463178116;
%CQ2  = 8.3958815463178116;
%%{
%chi1 = 0.1;
%chi2 = 0.1;
%}
%{
chi1 = 0.1;
chi2 = 0.1;
%}
nu     = 0.25;
kappaT = 73.1180073413692497;
%X1   = 1/2*(1+sqrt(1-4*nu));
%X2   = 1-X1;
%at1  = X1.*chi1;
%at2  = X2.*chi2;


%-----------------------
% No tides
%-----------------------
y     = load('h22_q1_bbh.dat');
t     = y(1:236754,1);
phi   = y(1:236754,3);
omg0  = Deriv4(phi,t,4);
domg0 = Deriv4(omg0,t,4);
Qomg0 = omg0.^2./domg0;

%--------------------
% Tides
%--------------------
y     = load('h22_q1_sly.dat');
t     = y(1:236740,1);
phi   = y(1:236740,3);
omg   = Deriv4(phi,t,4);
domg  = Deriv4(omg,t,4);
Qomg  = omg.^2./domg;



figure
plot(omg0/(2*pi),Qomg0,'k');
hold on;
plot(omg/(2*pi),Qomg,'--r');
legend('BBH','SLy');
%---------------
% ancillary grid
%---------------
omgL = max(omg0(1),omg(1));
omgR = min(omg0(end),omg(end));

w        = omgL:0.5:omgR;
f        = w/(2*pi);
Qw0      = spline(omg0,Qomg0,w);
Qw       = spline(omg,Qomg,w);
QomgTEOB = (Qw-Qw0)./kappaT;

Msun = 4.925491025543575903411922162094833998e-6;
M    = (1.35+1.35)*Msun;

Omg = M*w/2;
x   = Omg.^(2/3);
alpha22 = 85/14;
QomgT    = -65/6*x.^(5/2).*(1+4361/624*x - 4*pi*x.^(3/2)+(4614761/122304 + 4/3*alpha22).*x.^2 - 4283/156*pi*x.^(5/2));
QomgT1PN = -65/6*x.^(5/2).*(1+4361/624*x);
QomgT15PN = -65/6*x.^(5/2).*(1+4361/624*x-4*pi*x.^(3/2));
QomgT2PN = -65/6*x.^(5/2).*(1+4361/624*x - 4*pi*x.^(3/2)+(4614761/122304 + 4/3*alpha22).*x.^2);

QomgNRT_25PN =-65/6*x.^(5/2).*(1. + 6.98878*x - 16.888*x.^(3/2) + 125.907*x.^2 - 782.845*x.^(5/2));

PNRTidal = (0.000021275288367813374 - 0.0011229777736999561*x ...
    + 0.001944114575762419*power(x,1.5) + 0.01912798374299056*power(x,2)...
    - 0.07086763335265246*power(x,2.5) + 0.00842014977845372*power(x,3)...
    + 0.42522399043861203*power(x,3.5) - 1.7022194880891877*power(x,4)...
    + 4.767663928317682*power(x,4.5) - 7.613819168913324*power(x,5)...
    + 4.8384161378813495*power(x,5.5))...
    ./power(0.02770927429410624 - 0.5520794034122195*x + 1.*power(x,1.5),3);

QomgNRT_full = -65/6*x.^(5/2).*PNRTidal;

%----------------------------------------------
% Compute differences with respect EOB baseline
%----------------------------------------------
DQomgNRT  = QomgTEOB-QomgNRT_full;
DQomg25PN = QomgTEOB-QomgT;
DQomg1PN  = QomgTEOB-QomgT1PN;
DQomg15PN  = QomgTEOB-QomgT15PN;
DQomg2PN  = QomgTEOB-QomgT2PN;

%======================================================
% Cleaning differences in the very low-frequency regime
%======================================================
ncut = 2200;

F   = f(1:ncut);
dy1 = smooth(DQomg1PN(1:ncut),100,'sgolay',1);
dy15 = smooth(DQomg15PN(1:ncut),100,'sgolay',1);
dy2 = smooth(DQomg2PN(1:ncut),100,'sgolay',1);
dy25 = smooth(DQomg25PN(1:ncut),100,'sgolay',1);
dynr = smooth(DQomgNRT(1:ncut),100,'sgolay',1);



figure('Name','EOB/PN comparison');
plot(f,QomgTEOB,'r');
hold on;
plot(f,QomgT1PN,'c');
plot(f,QomgT15PN,'m');
plot(f,QomgT,'b');
hold on;
plot(f,QomgNRT_25PN,'--k');
plot(f,QomgNRT_full,'g');
set(gca,'XMinorTick','on','YMinorTick','on','FontSize',14,'FontName','Times');
legend('TEOBResumS','2.5PN','1PN','1.5PN','NRTidal@2.5PN','NRTida','Location','SouthWest');
xlabel('$f$ [Hz]','Fontsize',20,'Interpreter','Latex');
ylabel('$Q_\omega^T/\kappa_T$','Fontsize',20,'Interpreter','Latex');

%-----------------
% Plot differences
%-----------------
figure('Name','PlotDifferences');
plot(f,kappaT.*DQomg25PN,'b');
hold on;
plot(f,kappaT.*DQomg1PN,'c');
plot(f,kappaT.*DQomg15PN,'m');
plot(f,kappaT*DQomgNRT,'g');
line([f(1) f(end)],[0 0],'LineStyle',':','Color','k');
ylabel('$\Delta Q_\omega^{T_x}$','Fontsize',20,'Interpreter','Latex');
legend('X=EOB-2.5PN','X=EOB-1PN','X=EOB-1.5PN','X=EOB-NRTidal','Location','SouthWest');
xlabel('$f$ [Hz]','Fontsize',20,'Interpreter','Latex');

figure('Name','Behavior at very low frequency');
plot(F,kappaT*dy1,'--m');
hold on;
plot(F,kappaT*dy15,'-m');
plot(F,kappaT*dy2,'--b');
plot(F,kappaT*dy25,'b');
plot(F,kappaT*dynr,'g');
line([F(1) F(end)],[0 0],'LineStyle',':','Color','k');
xlabel('$f$ [Hz]','Fontsize',20,'Interpreter','Latex');
ylabel('$\Delta Q_\omega^{T_X}$','Fontsize',20,'Interpreter','Latex');
legend('X=EOB-1PN','X=EOB-1.5PN','X=EOB-2PN','X=EOB-2.5PN','X=EOB-NRTidal','Location','SouthWest');


figure('Name','PlotDifferences -II');
plot(2*pi*M*f,kappaT.*DQomg25PN,'b');
hold on;
plot(2*pi*M*f,kappaT.*DQomg1PN,'c');
plot(2*pi*M*f,kappaT.*DQomg15PN,'m');
plot(2*pi*M*f,kappaT*DQomgNRT,'g');
line([2*pi*M*f(1) 2*pi*M*f(end)],[0 0],'LineStyle',':','Color','k');
ylabel('$\Delta Q_\omega^{T_X}$','Fontsize',20,'Interpreter','Latex');
legend('EOB-2.5PN','EOB-1PN','EOB-1.5PN','EOB-NRTidal','Location','SouthWest');
xlabel('$M\omega$','Fontsize',20,'Interpreter','Latex');







end

