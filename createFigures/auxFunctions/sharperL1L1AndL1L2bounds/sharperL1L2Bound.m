function [f_val, arg_val] = sharperL1L2Bound(x, w, beta)

% [f_val, arg_val] = sharperL1L2Bound(x, w, beta)
% 
% Computes the Gaussian distance of the L1-L2 function numerically without
% using approximations of the Q-function and finds the optimal 
% parameter. 
%
% Input:
%   - x: optimal point (n x 1 vector)
%   - w: prior information (n x 1 vector)
%   - beta (positive number)
%
% Output:
%   - f_val:   exact Gaussian distance
%   - arg_val: value of the optimal parameter
%
% Method:
%   Derivative-free method fminsearch in Matlab
%
% =========================================================================
% Copyright (C) 2015  Joao Mota
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% =========================================================================


% ========================================================================= 
% Check for errors

n = length(x);

if length(w) ~= n
  error('x and w must have the same length')
end
if beta <= 0
    error('beta must be a positive number')
end
% =========================================================================

% =========================================================================
% Compute I, J, and associated sets and constants 

I  = (x~= 0);
Ic = (x== 0);
Ip = (x > 0);
In = (x < 0);

s  = sum(I);
% =========================================================================


% =========================================================================
% Define function and call fminsearch

q_term = sum((1 + beta*(x(Ip) - w(Ip))).^2) + sum((1 - beta*(x(In) - w(In))).^2);

f = @(t) q_term*t^2 + computeEgOverIndices(-t*beta*w, t*ones(n,1), Ic);

arg_val = fminsearch(f, 1);

f_val = f(arg_val) + s + 1;
% ====================================================================== 
