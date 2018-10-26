function h_taper = schwartz_tapering(h,time,freqM,varargin)
% This function is supposed to apply a schwartz tapering to the beginning
% and end of the waveform h. IMPORTANT: The input waveform h should no
% longer contain initial state radiation a.k.a. Junk radiation.

% Set test to 1 to check to compare the tapered waveform.
test        = 0;

% Some organizing parameters
detT_start  = 0;
detT_end    = 0;

% Defining Boundary conditions for the tapering.
% Default values!
fMi_start   = freqM(1);
fMi_end     = freqM(end);
TTap_start  = pi;
TTap_end    = 0.5*pi;
% Loading of tapering conditions if specified in the input.
if nargin>3
    if varargin{1}>0
        fMi_start   = varargin{1};
        detT_start  = 1;
    end
    if varargin{2}>0
        fMi_end     = varargin{2};
        detT_end    = 1;
    end
    if varargin{3}>0
        TTap_start  = varargin{3};
    end
    if varargin{4}>0
        TTap_end    = varargin{4};
    end
end

% hardcoded parameters that SHOULD NOT BE CHANGED!!!!
alph    = 1;
ncyc    = 1;

% determining the time interval on which to apply the tapering for the low 
% frequency tapering
if detT_start
    j_in    = find(freqM>fMi_start,1,'first');
else
    j_in    = 1;
end
T_in            = time(j_in);
TTap_start      = TTap_start./fMi_start;
j_in2           = find((time-T_in)<TTap_start,1,'last');



time_start      = time(j_in:j_in2)-T_in;
% small modification to avoid the initial singularity.
time_start(1)   = time_start(1)+0.00001; 

% Computing the low frequency tapering.
Z               = 1./(time_start) + 1./(time_start - ncyc*TTap_start);
ftaper          = 1./(1+exp(TTap_start*alph*Z));

% determining the time interval on which to apply the tapering for the high 
% frequency tapering
if detT_end
    j_end   = find(freqM<fMi_end,1,'last');
else
    j_end   = length(freqM);
end
T_end           = time(j_end);
TTap_end        = TTap_end./fMi_end;
j_end2          = find(time>(T_end-TTap_end),1,'first');


time_end        = time(j_end2:j_end)-(T_end-TTap_end);
% small modification to avoid the initial singularity.
time_end(end)   = time_end(end)-0.00001; 

% Computing the high frequency tapering.
Z               = 1./(time_end) + 1./(time_end - ncyc*TTap_end);

ftaper_end      = 1-1./(1+exp(TTap_start*alph*Z));

% applying the tapering:
h_taper                 = h;
h_taper(j_in:j_in2)     = h_taper(j_in:j_in2).*ftaper;
h_taper(j_end2:j_end)   = h_taper(j_end2:j_end).*ftaper_end;
if detT_end && j_end<length(h)
    h_taper(j_end+1:end)= zeros(length(h)-j_end,1);
    disp('hi')
end
if detT_start && j_in>1
    h_taper(1:j_in-1)   = zeros(j_in-1,1);
    disp('hi')
end
if test
    figure
    subplot 121
    plot(time,real(h),'k')
    hold on
    plot(time,real(h_taper),'r--')
    xlim([time(j_in) time(j_in2)]);
    ylim([-max(abs(h(j_in:j_in2)))*1.1 max(abs(h(j_in:j_in2)))*1.1])
    
    subplot 122
    plot(time,real(h),'k')
    hold on
    plot(time,real(h_taper),'r--')
    xlim([time(j_end2) time(j_end)]);
    ylim([-max(abs(h(j_end2:j_end))*1.1) max(abs(h(j_end2:j_end)))*1.1])
end
end