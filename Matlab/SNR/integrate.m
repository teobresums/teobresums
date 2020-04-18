function [int] = integrate(x,y)

% [int] = integrate(x,y)
%
% Integrates yi = y(xi) along xi, which can be not equally spaced

%% Sanity check

if (size(x) ~= size(y))
    error('The two input vectors must have the same size');
end

%% Integration routine

% We are taking into account the possibility that x can not be equally
% spaced, and we want to avoid using any loop:
%
%   int = 0.5 * sum[ d(i) * ( y(i+1) + y(i) ) ]_[1:N-1]
%
% where d(i) = x(i+1)-x(i).

% (A)  d(i)   , i \in [1:N-1]

daux = diff(x);

% (B)  y(i+1) , i \in [1:N-1]

yaux_p1    = y;
yaux_p1(1) = [];

% (C)  y(i)   , i \in [1:N-1]

yaux      = y;
yaux(end) = [];

% Compute the integral

int = 0.5 * sum( daux .* ( yaux_p1 + yaux ) );

return;