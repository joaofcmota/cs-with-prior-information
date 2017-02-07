function [f_val, arg_val] = sharperL1L1Bound(x, w, beta)

% [f_val, arg_val] = sharperL1L1Bound(x, w, beta)
% 
% Computes the Gaussian distance of the L1-L1 function numerically without
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

I = (x~=0);
J = (x~=w);

IJc  = I.*(~J);
IcJ  = (~I).*J;
IcJc = (~I).*(~J);

h     = sum( (x>0).*(x<w) ) + sum( (x<0).*(x>w) );
h_bar = sum( (x>0).*(x>w) ) + sum( (x<0).*(x<w) );
% =========================================================================


% =========================================================================
% Define function and call fminsearch

q_term = h_bar*(beta + 1)^2 + h*(beta - 1)^2;

sg_x  = sign(x);
sg_xw = sign(x - w);

f = @(t) q_term*t^2 + computeEgOverIndices(t*sg_x, t*beta*ones(n,1), IJc) ...
    + computeEgOverIndices(t*beta*sg_xw, t*ones(n,1), IcJ) ...
    + computeEgOverIndices(zeros(n, 1), t*(beta + 1)*ones(n,1), IcJc);

arg_val = fminsearch(f, 1);

f_val = f(arg_val) + h + h_bar + 1;
% ====================================================================== 
