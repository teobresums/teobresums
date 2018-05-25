function [M,rho2,gamma,fgamma] = DAmatchEOBEOB(fM, Af_EOB_1, Phf_EOB_1, Af_EOB_2, Phf_EOB_2, M_EOB, IFO)

% [M] = DAmatchEOBEOB(fM, Af_EOB_1, Phf_EOB_1, Af_EOB_2, Phf_EOB_2, M_EOB, IFO)
%
% Computes the MATCH between two EOB waveforms already generated and
% evaluated at the same fM values.
%
%   - We define the "match" between two waveforms, making use of the
%     internal product:
%
%            M(l,l') = <g(l)|h(l')> / sqrt( <g(l)|g(l)> <h(l')|h(l')> )
%
%     This is the only output of the function.
%
%     NOTE: Since the match has been normalized, it is independent from any
%     constant factor in front of the amplitude (i.e.: deff)
%
% ----
%
%   - The EOB waveforms has been previously FFTd and the input is the 
%     amplitude and phase of the EOB waveform as a function of the
%     frequency (in units of M) fM: 
%
%           h(fM) = Af_EOB(fM) *exp(i*Phf_EOB(fM))
%
%   - The EOB waveforms has still two free parameters: total mass (M_EOB)
%   and deff_EOB, although the match is independent of deff_EOB. We'll use
%   M_EOB to go from fM --> f.
%

if (M_EOB <= 0)
    M = 0;
    return;
end


%% Given M, compute "f" from "fM"

f  = fM / M_EOB;


%% The GW amplitude in the frequency domain is proportional to ( M^2 / deff ).
%
%   NOTE: The M^2 comes from the units

Af_EOB_1 = Af_EOB_1 *M_EOB*M_EOB / 1;
Af_EOB_2 = Af_EOB_2 *M_EOB*M_EOB / 1;


%% Compute the one-sided PSD of the detector's noise at the f values

if     (strcmp(IFO,'LIGO1'))     
    Sh = ShLIGO1(f);
elseif (strcmp(IFO,'LIGO2'))     
    Sh = ShLIGO2(f);
elseif (strcmp(IFO,'LIGO3'))     
    Sh = ShLIGO3(f);
elseif (strcmp(IFO,'aLigo'))    
    Sh = aLigoSn(f); 
    %{
    y  = load('H1-GDS-CALIB_STRAIN.txt');
    Sh = 0.3*spline(y(:,1),y(:,2),f);
    %}
    
elseif (strcmp(IFO,'VIRGO2'))    
    Sh = ShVIRGO2(f);
elseif (strcmp(IFO,'LIGO3ast'))  
    Sh = ShLIGO3ast(f);
elseif (strcmp(IFO,'VIRGO2ast')) 
    Sh = ShVIRGO2ast(f);
else error('Wrong specified IFO, the options are LIGO1, LIGO2, LIGO3, VIRGO2, LIGO3ast, VIRGO2ast');
end


%% Compute the match, making use of its definition
%
%  NOTATION: g --> EOB_1 ; h --> EOB_2
%
% *****************************************************************
%  IMPORTANT NOTES: 
%GMsunc3
%    (a) The maximization over phi0 is done analytically, just by using an
%    orthogonal basis.
%
%    (b) The maximization over t0 can be done numerically, making use of
%    the Matlab's built-in function "fminbnd", between -Tobs/2 and Tobs/2.
%    In the full expression of the MATCH, only <g|h> depends on t0.
%
% *****************************************************************

% Computing the maximized <g|h>

gh_max = 4*integrate(f, Af_EOB_1.*Af_EOB_2./ Sh );



gh_0   = @(t0) 4*integrate(f, Af_EOB_1.*Af_EOB_2.*cos( 2*pi*t0*f + Phf_EOB_2        - Phf_EOB_1 ) ./ Sh );
gh_pi2 = @(t0) 4*integrate(f, Af_EOB_1.*Af_EOB_2.*cos( 2*pi*t0*f + Phf_EOB_2 + pi/2 - Phf_EOB_1 ) ./ Sh );

mgh_func = @(t0) - sqrt( gh_0(t0)^2 + gh_pi2(t0)^2 );
mlog_gh_func = @(t0) - log( sqrt( gh_0(t0)^2 + gh_pi2(t0)^2 ) );

% PLAN A: Standard minimization using "fminbnd"

%{
[t0_max_A mgh_min_A] = fminbnd( mgh_func, -1.0/(2*f(1)), 1.0/(2*f(1)), optimset('MaxIter', 1000,'MaxFunEval', 1000, 'TolFun', 1e-9*gh_max) );
[t0_max_B mgh_min_B] = fminbnd( mgh_func, -1.1/(2*f(1)), 0.9/(2*f(1)), optimset('MaxIter', 1000,'MaxFunEval', 1000, 'TolFun', 1e-9*gh_max) );
[t0_max_C mgh_min_C] = fminbnd( mgh_func, -0.9/(2*f(1)), 1.1/(2*f(1)), optimset('MaxIter', 1000,'MaxFunEval', 1000, 'TolFun', 1e-9*gh_max) );
%}
%%{
[t0_max_A mgh_min_A] = fminbnd( mgh_func, -1.3/(2*f(1)), 1.3/(2*f(1)), optimset('MaxIter', 1000,'MaxFunEval', 1000, 'TolFun', 1e-9*gh_max) );
[t0_max_B mgh_min_B] = fminbnd( mgh_func, -1.4/(2*f(1)), 1.2/(2*f(1)), optimset('MaxIter', 1000,'MaxFunEval', 1000, 'TolFun', 1e-9*gh_max) );
[t0_max_C mgh_min_C] = fminbnd( mgh_func, -1.2/(2*f(1)), 1.4/(2*f(1)), optimset('MaxIter', 1000,'MaxFunEval', 1000, 'TolFun', 1e-9*gh_max) );
%}


[mgh_min,ind_aux] = min([mgh_min_A mgh_min_B mgh_min_C]);
switch ind_aux
    case 1
        t0_max = t0_max_A;
    case 2
        t0_max = t0_max_B;
    case 3
        t0_max = t0_max_C;
    otherwise
        error('Problem selecting mgh_min from the 3 cases considered');
end
clear ind_aux t0_max_A t0_max_B t0_max_C mgh_min_A mgh_min_B mgh_min_C;

if ( mgh_min > -0.4*gh_max )    
    % PLAN B: If Plan A didn't work, try to first minimize the log function, 
    %         and then explore the normal function in the surrounding points

    [t0_max_log_A mlog_gh_min_A] = fminbnd( mlog_gh_func, -1.0/(2*f(1)), 1.0/(2*f(1)), optimset('MaxIter', 1000,'MaxFunEval', 1000) );
    [t0_max_log_B mlog_gh_min_B] = fminbnd( mlog_gh_func, -1.1/(2*f(1)), 0.9/(2*f(1)), optimset('MaxIter', 1000,'MaxFunEval', 1000) );
    [t0_max_log_C mlog_gh_min_C] = fminbnd( mlog_gh_func, -0.9/(2*f(1)), 1.1/(2*f(1)), optimset('MaxIter', 1000,'MaxFunEval', 1000) );

    [mlog_gh_min,ind_aux] = min([mlog_gh_min_A mlog_gh_min_B mlog_gh_min_C]);
    switch ind_aux
        case 1
            t0_max_log = t0_max_log_A;
        case 2
            t0_max_log = t0_max_log_B;
        case 3
            t0_max_log = t0_max_log_C;
        otherwise
            error('Problem selecting mlog_gh_min from the 3 cases considered');
    end
    clear ind_aux t0_max_log_A t0_max_log_B t0_max_log_C mlog_gh_min_A mlog_gh_min_B mlog_gh_min_C;
   
    [t0_max mgh_min_new] = fminbnd( mgh_func, t0_max_log - 0.1/(2*f(1)), t0_max_log + 0.1/(2*f(1)), optimset('MaxIter', 1000,'MaxFunEval', 1000, 'TolFun', 1e-9*gh_max) );

    mgh_min = min(mgh_min_new, mgh_min);

    if ( mgh_min > -0.4*gh_max )
        
        % PLAN C: If Plan B also didn't work, try to minimize with the
        % "fminsearch" function
        
        [t0_max mgh_min_new] = fminsearch( mgh_func, t0_max_log, optimset('MaxIter', 1000,'MaxFunEval', 1000, 'TolFun', 1e-9*gh_max) );

        mgh_min = min(mgh_min_new, mgh_min);

        if ( mgh_min > -0.4*gh_max )

            % PLAN D: Brute force: scan over 1000 pts and start a "fminsearch" from there
            
            t_aux = linspace(-1/(2*f(1)), 1/(2*f(1)), 1000);
            
            for i=1:size(t_aux,2)
                mgh_aux(i) = mgh_func(t_aux(i));
            end
            
            [val_aux ind_aux] = min(mgh_aux);
            t0_max_aux = t_aux(ind_aux);
            [t0_max mgh_min_new] = fminsearch( mgh_func, t0_max_aux, optimset('MaxIter', 1000,'MaxFunEval', 1000, 'TolFun', 1e-9*gh_max) );

            mgh_min = min(mgh_min_new, mgh_min);
        end
    end
end

gh = - mgh_min;


% Computing the terms in the denominator and the final MATCH

gg = 4*integrate(f, Af_EOB_1 .* Af_EOB_1 ./ Sh);
hh = 4*integrate(f, Af_EOB_2 .* Af_EOB_2 ./ Sh);

M = gh / sqrt(gg*hh); 

%===================================================
% computing signal to noise ratio using EOB waveform
%===================================================
rho2   = 4*integrate(f,Af_EOB_1.*Af_EOB_1./Sh);
gamma  = (Af_EOB_1.*Af_EOB_1./Sh)./integrate(f,Af_EOB_1.*Af_EOB_1./Sh);
fgamma = f.*gamma;


return;