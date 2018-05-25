function ComputeSNR
% ComputeSNR.
% This functions computes the SNR integrand for the SNR  calculation. The
% SNR can be obtained integrating this integrand.
%
% The input is the fourier transform of the signal in mass-scaled units,
% that is
%
% fM, phase, amplitude
%
% The fft can be done adapting the function teobresums_fft_driver.m that is
% designed to FT the output of teobresums.
%
% anagar 25/5/2018

% Specify parameters
M = 80;  % Total mass in solar-mass units

%=========================================================================
% loading datasets. This routine was originally conceived to compute the 
% faithfulness (of match). The match of a waveform with itself is the SNR.
% Here we compute the SNR as the match of a waveform d1 with itself, d2
%=========================================================================
d1   = load('Data_FT/FT_bbh_80Msun.dat');
d2   = load('Data_FT/FT_bbh_80Msun.dat');    

% set colors (might be useful)
gray           = [0.8 0.8 0.8];
orange         = [1 0.694117647058824 0.392156862745098];
azul           = [0 0.749019607843137 0.749019607843137];
deep_green     = [0.168627450980392 0.505882352941176 0.337254901960784];
purple         = [0.8 0.2 1];
light_brownish = [0.929411768913269 0.694117665290833 0.125490203499794];

%nsimI   = 3   ;
%nsimE   = 3;
showfig = 1;


    

eob.fM     = d1(:,1);
eob.PhaseF = d1(:,2);
eob.AF     = d1(:,3);
eob.dfM    = eob.fM(2)-eob.fM(1);
eob.IQomg = Deriv4(eob.PhaseF,eob.fM,4);
eob.Qomg  = Deriv4(eob.IQomg,eob.fM,4);


nr.fM     = d2(:,1);
nr.PhaseF = d2(:,2);
nr.AF     = d2(:,3);
nr.dfM    = nr.fM(2)-nr.fM(1);
nr.IQomg = Deriv4(nr.PhaseF,nr.fM,4);
nr.Qomg  = Deriv4(nr.IQomg,nr.fM,4);

% Find the end of the unphysical reagion and remove it
[~, nr.jmin]  = findpeaks(-nr.PhaseF);
[~, eob.jmin] = findpeaks(-eob.PhaseF);

%----------------------------------------------
% take only the physical part of the FFT signal
%---------------------------------------------
eob.fM     = eob.fM(eob.jmin(1):end);
eob.PhaseF = eob.PhaseF(eob.jmin(1):end);
eob.AF     = eob.AF(eob.jmin(1):end);
nr.fM      = nr.fM(nr.jmin(1):end);
nr.PhaseF  = nr.PhaseF(nr.jmin(1):end);
nr.AF      = nr.AF(nr.jmin(1):end);
    
%{
    figure('Name','Raw data')
    plot(nr.fM,nr.PhaseF,'k')
    hold on;
    plot(eob.fM,eob.PhaseF,'r')
    xlim([0 0.2])
    
    figure
    semilogy(nr.fM,nr.AF,'k')
    hold on;
    semilogy(eob.fM,eob.AF,'r')
%}
    

A_end    = 1.e-2;  % cut the signal when A < 1e-3; Noise after.
eob.jend = find(eob.AF<A_end,1,'first');
nr.jend  = find(nr.AF<A_end,1,'first');


eob.fM     = eob.fM(1:eob.jend);
eob.AF     = eob.AF(1:eob.jend);
eob.PhaseF = eob.PhaseF(1:eob.jend);
nr.fM      = nr.fM(1:nr.jend);
nr.AF      = nr.AF(1:nr.jend);
nr.PhaseF  = nr.PhaseF(1:nr.jend);
    

if showfig
    figure
    subplot(2,1,1);
    semilogx(nr.fM,nr.PhaseF,'k')
    hold on;
    semilogx(eob.fM,eob.PhaseF,'r')
    subplot(2,1,2);
    loglog(nr.fM,nr.AF,'k')
    hold on;
    loglog(eob.fM,eob.AF,'r')
end


eob.N = length(eob.fM);
nr.N  = length(nr.fM);
eob.dfM = eob.fM(2)-eob.fM(1);
nr.dfM  = nr.fM(2)-nr.fM(1);
dfM     = min(nr.dfM,eob.dfM);

fM_min = max(eob.fM(1),nr.fM(1));
fM_end = min(eob.fM(end),nr.fM(end));


%---------------------------------------------
% Create a new grid, with the same fM for both
%---------------------------------------------
disp(['fM(min) =',num2str(fM_min)]);
disp(['fM(max) =',num2str(fM_end)]);
    

fM = fM_min:dfM:fM_end;
nr.PhaseF  = spline(nr.fM,nr.PhaseF,fM);
nr.AF      = spline(nr.fM,nr.AF,fM);
eob.PhaseF = spline(eob.fM,eob.PhaseF,fM);
eob.AF     = spline(eob.fM,eob.AF,fM);
    
%{
figure
subplot(2,1,1)
semilogy(fM,nr.AF,'k')
hold on
semilogy(fM,eob.AF,'r')
subplot(2,1,2)
plot(fM,nr.PhaseF,'k')
hold on
plot(fM,eob.PhaseF,'r')
%}
% compute the match between the two FT
% Put the frequency in Hz


qnum      = 1;   % mass ratio
SNR2_expr = 13^2;
f0        = 100;
fmrg      = 450;

GMsunc2 = 1.476625038;
GMsunc3 = 4.925490947*1e-6;
%------------------------------
% Compute SNR in pure EOB sense
%------------------------------

M_EOB   = M*GMsunc2; %total mass of the system
c       = 299792.458; %km/sec gray
fMHz    = c*fM;


fHz_Init  = fMHz(1)/M_EOB;
fHz_End   = fMHz(end)/M_EOB;

%=====================================
% Compute unfaithfulness with TaylorF2
%=====================================
fMSchwLSO          = 1/pi*6.^(-3/2);
fMSchwLSOHz        = c*fMSchwLSO;
jLSO               = find(fMHz>=fMSchwLSOHz,1,'first');
sample.fMHz        = fMHz(1:jLSO);
sample.fM          = fM(1:jLSO);
sample.eob.AF      = eob.AF(1:jLSO);
sample.eob.PhaseF  = eob.PhaseF(1:jLSO);
nu                 = qnum./((1+qnum).^2);

%========================================================
% Match EOB-EOB: this is the rho^2 that enters in the SNR
%========================================================
[Match,rho2,gamma,fgamma] = DAmatchEOBEOB(fMHz, eob.AF, eob.PhaseF, nr.AF, nr.PhaseF, M_EOB, 'aLigo');

out.gamma    = gamma;
%
out.fgamma   = SNR2_expr.*log(10)*fgamma; % argument for the SNR
out.fHz      = fMHz./M_EOB;
out.AF       = eob.AF;
f00          = 35;

figure('Name','Signal to noise ratio integrand');
semilogx(out.fHz,out.fgamma,'Color','k','LineStyle','-','LineWidth',1);
hold on;
ylabel('$f \gamma(f)$','FontSize',20,'FontName','Times','Interpreter','Latex');
line([f00 f00],[0 1500],'LineStyle','--','Color',gray,'LineWidth',2);
line([f0 f0],[0 1500],'LineStyle','--');
line([fmrg fmrg],[0 1500],'LineStyle','--','Color','r');
set(gca,'FontSize',12,'FontName','Times');
xlabel('$f [Hz]$','FontSize',20,'FontName','Times','Interpreter','Latex');
ylim([0 900]);



