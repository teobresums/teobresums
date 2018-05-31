function [fM, AF, PhaseF, tc_PN, phic_PN, fM_i, fM_e] = teobresums_fft(Mtot, filename,outfile,q,n_cutoff, Tau,alpha,tapering,plotting)

% [fM, AF, PhaseF, tc_PN, phic_PN] = getEOBwvf(fM_min, filename, q,n_cutoff, plotting)
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

%% Default values

%if (fM_min <= 0)
%    fM_min = 1e-4; % min fM value for ground-based detector: 10Hz * 3Ms
%end

%% Defining some constants and aux variables

fM_LSO = 1/(pi*power(6,3/2));
nu_EOB = q/(1+q)^2;

gray   = [139 137 137]/255;
orange = [255 153   0]/255;


%% Loading the EOB time-domain waveform from a text file

%load(filename)


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
%
dt = 0.5;
output_path='./';

        
%---------------------------
% define waveform quantities
%---------------------------
% loading EOB waveform 
y=load(filename); 

% renaming quantities
Msun  = 4.925491025543575903411922162094833998e-6;
M     = Mtot*Msun;
% =========================================
% FIXME
% remove the first point: for some reason
% the C++ with BBH gives nan the first line
% =========================================
tmpT  = y(2:end,1)/M;
time  = tmpT;
A     = y(2:end,2);
phi   = y(2:end,3);

% define psi with the correct normalization
psi    = A.*exp(-1i.*phi)/sqrt(24);

% Compute useful quantities from the waveform
phi   = -unwrap(angle(psi));
freqM =  Deriv4(phi,tmpT,4)/(2*pi);
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
    ftaper     =  1/2*(1+tanh(time*alpha(1)-Tau(1)));
    ftaper_end = -1/2*(tanh((time-Tau(2)*time(end))*alpha(2))-1);
    
    % the raw waveform
    hraw = h;
    % taper it at beginnend and at the end
    h    = h.*ftaper; %.*ftaper_end;

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


%% Extending the post-merger waveform (linear interpolation in phase and log-amplitude)
% ( we are using the last 30M points of the time-dom waveform for the interpolation )
%     ( and we are extending the waveform until it reaches 1e-15 of max(A) ) 
%{
interp_ind = find( max(time)-time > 30,1,'Last' ):size(time,1);

[outfit_logA]  = fit(time(interp_ind)-time(interp_ind(1)),log(A(interp_ind)),'poly1');
[outfit_phase] = fit(time(interp_ind)-time(interp_ind(1)),Phase(interp_ind) ,'poly1');
%
fit_lA_b  = outfit_logA.p2  - time(interp_ind(1)) * outfit_logA.p1;    % y-intersec
fit_lA_m  = outfit_logA.p1;                                            % slope
%
fit_ph_b  = outfit_phase.p2  - time(interp_ind(1)) * outfit_phase.p1;  % y-intersec
fit_ph_m  = outfit_phase.p1;                                           % slope

% Compute up to what time value do we have to extend the waveform.

Athr      = max(A)/1e20;
%Athr      = max(A)/1e15;
time_max  = ( log(Athr)-fit_lA_b ) / fit_lA_m;
index_max = ceil(time_max / deltat);

if (index_max > size(h,1))
    time_aux(:,1)  = time(1) + ( size(h,1):index_max-1 ) * deltat;
    A_aux(:,1)     = exp(fit_lA_m*time_aux + fit_lA_b);
    Phase_aux(:,1) =     fit_ph_m*time_aux + fit_ph_b;

    time  = [time; time_aux];
    A     = [A; A_aux];
    Phase = [Phase; Phase_aux];
end

clear time_aux A_aux Phase_aux h_aux;
%}

% Computing 1st derivatives of Phase and A wrt t
%{
index_max = ceil(time(end)/deltat);

[time_aux, Phasedot] = derive(time,Phase);
[time_aux, Adot]     = derive(time,A);

% Computing 2nd derivative of Phase wrt t and smooth it

Nsmooth = 50;
[time_aux_aux, Phasedotdot] = derive(time_aux,Phasedot);
Phasedotdot = smooth(Phasedotdot, Nsmooth);
%
Phasedotdot(1:ceil(Nsmooth/2))        = [];
Phasedotdot(end-ceil(Nsmooth/2):end)  = [];
%}

%% Using a semitukey window of length = 1/3 of the full signal

%length_win = deltat*length(time)/18;
%win = semitukeywin( size(h,1) , 2*length_win/(deltat*size(h,1)) );

%clear A_ext Phase_ext;


%% Computing the FFT

%hf      = deltat * (fft(h).*win);
%N       = 80000;
%hf      = deltat * fft(h,N);
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

%[Af_mx,jAmx] = max(Af_EOB);
jAmx = find(fM_EOB(1:end)>fM_Amx,1,'first');

%[outfit_Phf_EOB]  = fit(fM_EOB(floor(end/2):end),Phf_EOB(floor(end/2):end),'poly1');
[outfit_Phf_EOB]  = fit(fM_EOB(1:jAmx),Phf_EOB(1:jAmx),'poly1');


Phf_EOB = Phf_EOB - (fM_EOB*outfit_Phf_EOB.p1 + outfit_Phf_EOB.p2);

AF     = Af_EOB; %.*match_win  + Af_PN.*(1-match_win);
PhaseF = Phf_EOB; %.*match_win + Phf_PN.*(1-match_win);
fM     = fM_EOB;




if (plotting == 1)
            
    figure('Position',[46 264 560 656], 'Name','PN-EOB merging AMPLITUDES');
    %
    subplot(2,1,1);
    loglog(fM,AF);
    hold on;
    %loglog(fM,Af_PN,'-r');
    %axis([3e-4 1e-1 1 5000]);
    YLim_aux = get(gca,'YLim');
    %semilogx([fM(ind_i) fM(ind_i)], YLim_aux, 'Color', gray);
    %semilogx([fM(ind_e) fM(ind_e)], YLim_aux, 'Color', gray);
    %loglog([fM_EOB_min fM_EOB_min], YLim_aux, '--g');
    %loglog([fM_LSO fM_LSO], YLim_aux, '--', 'Color', orange);
    xlabel('$fM$','FontSize',14,'Interpreter','Latex');
    ylabel('$\tilde{A}(fM)$','FontSize',20,'Interpreter','Latex');
    title('AMPLITUDES');
    %
    %subplot(2,1,2);
    %semilogx(fM,Af_PN./Af_EOB);
    %hold on;
    %set(gca,'XLim',[3e-4 1e-1],'YLim',[0.95 1.05]);
    %XLim_aux = get(gca,'XLim');
    %YLim_aux = get(gca,'YLim');
    %semilogx([fM(ind_i) fM(ind_i)], YLim_aux, 'Color', gray);
    %semilogx([fM(ind_e) fM(ind_e)], YLim_aux, 'Color', gray);
    %loglog([fM_EOB_min fM_EOB_min], YLim_aux, '--g');
    %semilogx([fM_LSO fM_LSO], YLim_aux, '--', 'Color', orange);loglog(nr.fM,nr.AF,'k')
    %semilogx(XLim_aux,[1 1],'--','Color',gray);
    xlabel('fM');
    ylabel('A_{ PN} / A_{ EOB}');

    %figure('Position',[679 264 560 656], 'Name','PN-EOB merging PHASES');
    %
    subplot(2,1,2);
    plot(fM,Phf_EOB);
    hold on;
    %semilogx(fM,Phf_PN,'-r');
    %axis([3e-4 1e-1 -900 20]);
    YLim_aux = get(gca,'YLim');
    %semilogx([fM(ind_i) fM(ind_i)], YLim_aux, 'Color', gray);
    %semilogx([fM(ind_e) fM(ind_e)], YLim_aux, 'Color', gray);
    %loglog([fM_EOB_min fM_EOB_min], YLim_aux, '--g');
    %semilogx([fM_LSO fM_LSO], YLim_aux, '--', 'Color', orange);
    xlabel('fM');
    ylabel('Phase');
    title('PHASES');
    %
    %subplot(2,1,2);
    %semilogx(fM,Phf_PN - Phf_EOB);
    %hold on;
    %set(gca,'XLim',[3e-4 1e-1],'YLim',[-0.05 0.05]);
    %XLim_aux = get(gca,'XLim');
    %YLim_aux = get(gca,'YLim');
    %semilogx([fM(ind_i) fM(ind_i)], YLim_aux, 'Color', gray);
    %semilogx([fM(ind_e) fM(ind_e)], YLim_aux, 'Color', gray);
    %loglog([fM_EOB_min fM_EOB_min], YLim_aux, '--g');
    %semilogx([fM_LSO fM_LSO], YLim_aux, '--', 'Color', orange);
    %semilogx(XLim_aux,[0 0],'--','Color',gray);
    %xlabel('fM');
    %ylabel('Phase_{ PN} - Phase_{ EOB}');
    
    
    % Also compute the curvature of dPhf vs fM
    
    %{
    fM_AUX   = fM(ind_i:ind_e);
    dPhF_AUX = Phf_PN(ind_i:ind_e) - Phf_EOB(ind_i:ind_e);
    %
    f0_AUX = (max(fM_AUX)+min(fM_AUX))/2;
    %
    outfit_AUX = fit(fM_AUX-f0_AUX,dPhF_AUX,'poly2');
    dPhF_pp    = abs(outfit_AUX.p1/2);
    dQw        = f0_AUX*f0_AUX*dPhF_pp;
    
    % Display the fitting values ...
    
    display(f0_AUX);
    display(dPhF_pp);
    display(dQw);

    % ... and add to the plot the curve
    
    semilogx(fM_AUX,outfit_AUX(fM_AUX-f0_AUX),'-k');
    %}
end

%{
AF     = Af_EOB;
PhaseF = Phf_EOB;
%}
DA_WriteFFT(outfile,fM,PhaseF,AF);


return

%% Downsampling the resulting signal (artificially oversampled by increasing Tobs)

% We are considering two possibilities:
%
%   A) Downsampling and re-evaluating it in an approximate log_scale
%   B) Just downsampling the artificial oversampling (linear scale)
%
% In principle, A) is more efficient than B) since it can store the same
% information with fewer number of points. The user can choose between one
% and the other (set "dwnsmpl_log" = 1/0) although we'll first check that
% the final number of points using A) is not smaller than a certain value,
% say 5000 points. In this case, we would force to do B) as well.

% 0.1) First, some preliminary computations

    % Computing the fM_EOB_min and fM_EOB_max values with actual
    % contribution from EOB

    % fM_EOB_min: "Ncyc" cycles after the actual EOB fmin
    %
    %  NOTE: We use an expression derived from a Taylor expansion of
    %  Phasedot up to linear order.

    Ncyc = 100;
    fM_EOB_min = sqrt( Phasedot(2)*Phasedot(2) + 4*pi*Ncyc*Phasedotdot(1) ) / (2*pi);

    % fM_EOB_max: we already have cut the waveform at higher frequencies

    fM_EOB_max = fM_EOB(end);
  

    % Compute aux quantities
    
    x = fM_EOB_max / fM_EOB_min;

    ind_fMmax = find(fM_EOB < fM_EOB_max,1,'Last');
    ind_fMmin = find(fM_EOB > fM_EOB_min,1,'First');
    N         = ind_fMmax - ind_fMmin + 1;


% 0.2) Then, compute M_min (i.e.: minimum number of points after
%     downsampling with log_scale)
    
    M_min = floor( 1 + log(x) / log( (x-1)/(N-1) + 1 ) );
    
    % If M_min < 5000, then force normal downsampling (i.e.: B))
    
    if (M_min < 5000)
        dwnsmpl_log = 0;
    end
    
    
% PROCEED WITH THE ACTUAL DOWNSAMPLING
%{
if (dwnsmpl_log == 1)
    
    % Computing the number of points after the downsampling, M. We don't
    % want to interpolate, so what we'll do is to take fM(i_min) and
    % fM(i_min+lambda) as the first pair of points of the new downsampled
    % log-scaled fM. The expression to obtain M is the following
    %
    %     M = floor( 1 + log(x) / log( lambda*(x-1)/(N-1) + 1 ) ),
    %
    % being x = fM_max/fM_min and N = size(fM(ind_min:ind_max))

        % We'll force M to be at most M_max

        M_max = 50000;
    
        % Compute the min lambda to have M < M_max
    
        lambda = ceil ( (N-1)/(x-1) * ( exp( log(x)/(M_max-1) ) - 1 ) );
    
        % Compute M
    
        M = floor( 1 + log(x) / log( lambda*(x-1)/(N-1) + 1 ) );

    % Computing the indices of fM and the other vectors that must survive.
    % The expression we shall use is:

    ind_dwnsmpl = round( ind_fMmin + (N-1)/(x-1) * ( power( x , ((1:M)-1)/(M-1) ) - 1 ) );

    % Downsampling fM_EOB, Af_EOB, Phf_EOB and hf_EOB

    Af_EOB  = Af_EOB(ind_dwnsmpl);
    Phf_EOB = Phf_EOB(ind_dwnsmpl);
    fM_EOB  = fM_EOB(ind_dwnsmpl);
    
elseif (dwnsmpl_log == 0)

    % The downsampling factor is given by ratio of lenghts after and
    % before doing the artificial extension. Moreover, since we have
    % already removed the freq = 0 component, we must add a dephase in the
    % downsampling of "dwnsmp_factor - 1".

    dwnsmp_factor = floor( size(time,1) / ( size(time,1)-size(time_ext,2) ) );

    Af_EOB  = downsample(Af_EOB,dwnsmp_factor,dwnsmp_factor-1);
    Phf_EOB = downsample(Phf_EOB,dwnsmp_factor,dwnsmp_factor-1);
    fM_EOB  = downsample(fM_EOB,dwnsmp_factor,dwnsmp_factor-1);

    % Also remove the components from fM < fM_EOB_min
    
    ind_fMmin = find(fM_EOB > fM_EOB_min,1,'First');
    
    Af_EOB(1:ind_fMmin-1)  = [];
    Phf_EOB(1:ind_fMmin-1) = [];
    fM_EOB(1:ind_fMmin-1)  = [];
    
else
    
    error('Wrong "dwnsmpl_log" option selected. Possible options: 0/1');
    
end
%}
clear time_ext;


%% Generating what will be the final fM vector for the full merged waveform

% We fill the fM vector in a different way depending whether fM_EOB is in
% "log" or "linear" scale.

if (dwnsmpl_log == 1)

    % A log-scale vector is generated as fM(i) = fM_min * step_factor^(i-1).
    % Here we already know step_factor

    step_factor = power(x,1/(M-1));

    % The fM vector for values greater than fM_EOB_min, will be fM_EOB which is
    % just an approximated log-scaled vector. So, here we generate a log-scaled
    % vector up to ind_aux and then concatenate with fM_EOB;

    ind_aux = floor( log(fM_EOB_min/fM_min) / log(step_factor) );

    fM = [fM_min*power(step_factor,(0:ind_aux-1))'; fM_EOB];

elseif (dwnsmpl_log == 0)

    % A lin-scale vector is generated as fM(i) = fM_min + step*(i-1).
    
    step = fM_EOB(2)-fM_EOB(1);
    
    % Here we generate the missing terms between fM_min and fM_EOB_min
    
    ind_aux = floor( (fM_EOB_min-fM_min)/step );
    
    fM = [fM_min+step*(0:ind_aux-1)'; fM_EOB];
    
end

% *************** Filling with zeros {X}_EOB ***************
Af_EOB  = [zeros(size(fM,1)-size(fM_EOB,1),1); Af_EOB];
Phf_EOB = [zeros(size(fM,1)-size(fM_EOB,1),1); Phf_EOB];
% **********************************************************


%% Merging with PN-SPA

% (A) Parameters that characterize the merge interval

N_pts_ovrlp = 200;

fM_ovrlp_min  = fM_EOB_min;
ind_ovrlp_min = find(fM > fM_ovrlp_min, 1,'First');

ind_i = ind_ovrlp_min;
ind_e = ind_i + N_pts_ovrlp;

fM_i = fM(ind_i);
fM_e = fM(ind_e);


% (B) Obtaining the general PN-SPA waveform evaluated at fM < fM_LSO

ind_LSO = find(fM > fM_LSO,1,'First') - 1;

if (ind_e > ind_LSO )
    warning('Too few wvf points before LSO -- Reducing the matching interval');
    ind_e = ind_LSO;
end

fM_PN          = fM(1:ind_LSO);
[Af_PN Phf_PN] = PNSPA(fM_PN,nu_EOB,1,1,0,0);

% Find t0 and phi0 such that maximize the PN-EOB overlapping and substract
% the contribution from the PN Phase

[outfit_Phf_PN] = fit(fM(ind_i:ind_e),Phf_PN(ind_i:ind_e)-Phf_EOB(ind_i:ind_e),'poly1');
Phf_PN = Phf_PN - ( outfit_Phf_PN.p2 + fM_PN*outfit_Phf_PN.p1 );

% Store the output values from the fit, because they represent tc_PN and
% phic_PN and are one of the output variables of the function

tc_PN   = outfit_Phf_PN.p1 / (2*pi);
phic_PN = - outfit_Phf_PN.p2;

% ************* Filling with zeros {X}_PN *************
Af_PN  = [Af_PN; zeros(size(fM,1)-size(fM_PN,1),1)];
Phf_PN = [Phf_PN; zeros(size(fM,1)-size(fM_PN,1),1)];
% *****************************************************


% (C) Generate the Cinf transition function and compute the hybrid
%    amplitude and phase

match_win = Cinftrans(0,1, ind_i,ind_e, size(fM,1));

AF     = Af_EOB.*match_win  + Af_PN.*(1-match_win);
PhaseF = Phf_EOB.*match_win + Phf_PN.*(1-match_win);

%{
AF     = Af_EOB;
PhaseF = Phf_EOB;
%}

DA_WriteFFT(outfile,fM,PhaseF,AF)


%%% ............ PLOTS TO CHECK THAT WE DID A GOOD MERGING ........... %%%

if (plotting == 1)
    
    figure('Position',[46 264 560 656], 'Name','PN-EOB merging AMPLITUDES');
    %
    subplot(2,1,1);
    loglog(fM,Af_EOB);
    hold on;
    loglog(fM,Af_PN,'-r');
    axis([3e-4 1e-1 1 5000]);
    YLim_aux = get(gca,'YLim');
    semilogx([fM(ind_i) fM(ind_i)], YLim_aux, 'Color', gray);
    semilogx([fM(ind_e) fM(ind_e)], YLim_aux, 'Color', gray);
    loglog([fM_EOB_min fM_EOB_min], YLim_aux, '--g');
    loglog([fM_LSO fM_LSO], YLim_aux, '--', 'Color', orange);
    xlabel('$fM$','FontSize',14,'Interpreter','Latex');
    ylabel('$\tilde{A}(fM)$','FontSize',20,'Interpreter','Latex');
    title('AMPLITUDES');
    %
    subplot(2,1,2);
    semilogx(fM,Af_PN./Af_EOB);
    hold on;
    set(gca,'XLim',[3e-4 1e-1],'YLim',[0.95 1.05]);
    XLim_aux = get(gca,'XLim');
    YLim_aux = get(gca,'YLim');
    semilogx([fM(ind_i) fM(ind_i)], YLim_aux, 'Color', gray);
    semilogx([fM(ind_e) fM(ind_e)], YLim_aux, 'Color', gray);
    loglog([fM_EOB_min fM_EOB_min], YLim_aux, '--g');
    semilogx([fM_LSO fM_LSO], YLim_aux, '--', 'Color', orange);
    semilogx(XLim_aux,[1 1],'--','Color',gray);
    xlabel('fM');
    ylabel('A_{ PN} / A_{ EOB}');

    figure('Position',[679 264 560 656], 'Name','PN-EOB merging PHASES');
    %
    subplot(2,1,1);
    semilogx(fM,Phf_EOB);
    hold on;
    semilogx(fM,Phf_PN,'-r');
    axis([3e-4 1e-1 -900 20]);
    YLim_aux = get(gca,'YLim');
    semilogx([fM(ind_i) fM(ind_i)], YLim_aux, 'Color', gray);
    semilogx([fM(ind_e) fM(ind_e)], YLim_aux, 'Color', gray);
    loglog([fM_EOB_min fM_EOB_min], YLim_aux, '--g');
    semilogx([fM_LSO fM_LSO], YLim_aux, '--', 'Color', orange);
    xlabel('fM');
    ylabel('Phase');
    title('PHASES');
    %
    subplot(2,1,2);
    semilogx(fM,Phf_PN - Phf_EOB);
    hold on;
    set(gca,'XLim',[3e-4 1e-1],'YLim',[-0.05 0.05]);
    XLim_aux = get(gca,'XLim');
    YLim_aux = get(gca,'YLim');
    semilogx([fM(ind_i) fM(ind_i)], YLim_aux, 'Color', gray);
    semilogx([fM(ind_e) fM(ind_e)], YLim_aux, 'Color', gray);
    loglog([fM_EOB_min fM_EOB_min], YLim_aux, '--g');
    semilogx([fM_LSO fM_LSO], YLim_aux, '--', 'Color', orange);
    semilogx(XLim_aux,[0 0],'--','Color',gray);
    xlabel('fM');
    ylabel('Phase_{ PN} - Phase_{ EOB}');
    
    % Also compute the curvature of dPhf vs fM
    
    fM_AUX   = fM(ind_i:ind_e);
    dPhF_AUX = Phf_PN(ind_i:ind_e) - Phf_EOB(ind_i:ind_e);
    %
    f0_AUX = (max(fM_AUX)+min(fM_AUX))/2;
    %
    outfit_AUX = fit(fM_AUX-f0_AUX,dPhF_AUX,'poly2');
    dPhF_pp    = abs(outfit_AUX.p1/2);
    dQw        = f0_AUX*f0_AUX*dPhF_pp;
    
    % Display the fitting values ...
    
    display(f0_AUX);
    display(dPhF_pp);
    display(dQw);

    % ... and add to the plot the curve
    
    semilogx(fM_AUX,outfit_AUX(fM_AUX-f0_AUX),'-k');
    
end

return;