function [fM, AF, PhaseF] = teobresums_fft(filename,outfile,q,n_cutoff, Tau,alpha,tapering,writeout)

% [fM, AF, PhaseF] = teobresums_fft(fM_min, filename, q,n_cutoff, tapering,writeout)
%
% Reads the EOB waveform from "filename" in the time domain, extends it to
% both earlier and later times, FFTs it and merges it with a PN Taylor T4
% SPA expression. The final output is an Amplitude and Phase in the Fourier
% domain [in particular A(fM) and Ph(fM))] valid always, because of the
% PN-EOB merging.
%
% phic_PN, tc_PN, fM_i and fM_e can be also output variables, for future 
% PN-EOB comparisons.
%
% If "plotting" = 1, the function also outputs two plots showing the
% consistency of the PN-EOB merging.
%
% IMPORTANT NOTE: Here we include the constant prefactors coming from
% spin-weighted (2,2) Spherical Harmonic, in order to get comparable
% results to PNSPA and Phen. There is only the effective distance remaining
% to be included, so actually the amplitude's output of this function is 
% AF x deff.
%

dwnsmpl_log = 1;

%% Define time and complex h

% The data file contains the complex values of Psi. It can be
% translated to GW strain x deff by the following relation 
% (see CQG 22, R167-R192 (2005)):
%
%     h = 1/r * sqrt( (l+2)!/(l-2)! ) * Ylm * Psi
%
% Applying it to the (2,2) mode and writting it in terms of deff:
%
%     deff*h = ( sqrt(4!) * sqrt(5/pi) / 4 ) * Psi
%
%
% -------
%     NOTE: The reason why here we are computing:
%
%        h = prefactor * conj(complex(Psi_loaded(:,2),Psi_loaded(:,3)));
%
%     is because Psi_loaded(:,2) contains the '+' polarization and
%     Psi_loaded(:,3) the 'x' one, and complex h is defined as 
%
%                           h = h_{+} - i*h_{x}
% -------
   
%==========================================================
%
% loading EOB waveform from TEOBResumS - C: geometric units
%
% input is h_22
%
%==========================================================
fid = fopen(filename, 'r');
dataload =  textscan(fid, '%n %n %n', 'HeaderLines', 0); 
fclose(fid);
time   = dataload{1};
A      = dataload{2};
phi    = dataload{3};
           
% define psi with the correct normalization from h22
psi    = A.*exp(-1i.*phi)/sqrt(24);

% Compute useful quantities from the waveform
phi   = -unwrap(angle(psi));
freqM =  Deriv4(phi,time,4)/(2*pi);
% frequency @mrg: needed to remove the linear part of the FFT phase
[Amx, imx]= max(abs(psi));
fM_Amx    = freqM(imx);

%=================================================
% prefactor for the l=2, m=\pm 2 mode, summed over
% the two +m and -m spherical harmonics
%=================================================
prefactor = sqrt(factorial(4)) * sqrt(5/pi) / 4;
h          = prefactor * conj(psi);

deltat = time(2)-time(1);

if tapering
    % tapering h    
%     ftaper     =  1/2*(1+tanh(time*alpha(1)-Tau(1)));
%     ftaper_end = -1/2*(tanh((time-Tau(2)*time(end))*alpha(2))-1);
    
    % the raw waveform
    hraw = h;
    % taper it at beginnend and at the end
%     h    = h.*ftaper; %.*ftaper_end;
    h   = schwartz_tapering(h,time,freqM);
    %%{
    figure('Name','Test tapering')
    subplot(3,1,[1 2])
    plot(time,real(hraw))
    hold on;
    plot(time,real(h),'r')
    plot(time,imag(hraw))
    plot(time,imag(h),'r')
    subplot(3,1,3)
    plot(time,freqM,'b')
    ylim([0 0.14]);
    %}
end






%% Computing amplitude and phase vs time

A     = abs(h);
Phase = unwrap(angle(h));


%% Removing the last data point, because it's not equally time-spaced

A(end)     = [];
Phase(end) = [];
h(end)     = [];
time(end)  = [];


%% Removing the data points where A = 0

ind_aux = find (A == 0);

A(ind_aux)     = [];
Phase(ind_aux) = [];
h(ind_aux)     = [];
time(ind_aux)  = [];

clear ind_aux;


%% Computing the FFT

hf      = deltat * fft(h);
fM_EOB(:,1) = 1/(size(h,1)*deltat) * (0:1:size(hf,1)-1);

% Removing the high frequency part that contributes with negligible
% amplitudes.

% FIXME: decide a rational to cut off the high-frequency part.
%fM_max = PhenfMcut(nu_EOB,0);

fM_max = fM_EOB(end,1);

ind_rmv         = find (fM_EOB > fM_max);
hf(ind_rmv)     = [];
fM_EOB(ind_rmv) = [];

% Removing the frequency = 0 term

hf(1)     = [];
fM_EOB(1) = [];

% Compute amplitude and phase

Af_EOB  = abs(hf);
Phf_EOB = unwrap(angle(hf));

%-----------------------------------------------------
% Remove the linear contribution coming from 2*pi*f*t0
%-----------------------------------------------------
% FIXME: one should do the linear fit ONLY where the signal
% is, which means approximately up to the merger frequency

jAmx = find(fM_EOB(1:end)>fM_Amx,1,'first');

%[outfit_Phf_EOB]  = fit(fM_EOB(floor(end/2):end),Phf_EOB(floor(end/2):end),'poly1');
%[outfit_Phf_EOB]  = fit(fM_EOB(1:jAmx),Phf_EOB(1:jAmx),'poly1');
%Phf_EOB = Phf_EOB - (fM_EOB*outfit_Phf_EOB.p1 + outfit_Phf_EOB.p2);

p       = polyfit(fM_EOB(1:jAmx),Phf_EOB(1:jAmx),1);
Phf_EOB = Phf_EOB - (fM_EOB*p(1) + p(2));

AF     = Af_EOB; %.*match_win  + Af_PN.*(1-match_win);
PhaseF = Phf_EOB; %.*match_win + Phf_PN.*(1-match_win);
fM     = fM_EOB;

%{
AF     = Af_EOB;
PhaseF = Phf_EOB;
%}
if writeout
    DA_WriteFFT(outfile,fM,PhaseF,AF);
end

return

