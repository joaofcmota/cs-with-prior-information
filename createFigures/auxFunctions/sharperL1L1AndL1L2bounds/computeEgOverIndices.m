function [f] = computeEgOverIndices(v1, v2, S)

% [f] = computeEgOverIndices(v1, v2, S)
%
% Computes the sum of 
%
%       E_g[dist(g, I(v1(i), v2(i)))^2]            (1)
%
% for all the indices i where S is nonzero. In (1), g is a Normal random 
% variable, and I(a, b) = [a-b, a+b] is an interval. If b < 0, the function 
% outputs an error.
%
% Input:
%   - v1: n x 1 vector of reals
%   - v2: n x 1 vector of reals
%   - S:  n x 1 binary vector
%
% Output:
%   - f: sum of (1) for all pairs v1(i), v2(i), where S(i) is true (= 1).
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
% Check input

n = length(v1);
if length(v2) ~= n || length(S) ~= n
    error('v1, v2, and S must all have the same length')
end
if sum(v2 < 0)
    fprintf('Warning: all entries of v2 must be nonnegative. Setting the negative ones to zero...\n')
    v2 = max(v2, 0);
end
% =========================================================================

% Q-function
Q = @(x) 0.5*erfc(x/sqrt(2));

% Phi
phi = @(x) exp(-0.5*x^2)/sqrt(2*pi);

% Expected distance to a point
Eg = @(a, b) (a - b)*phi(a - b) - (a + b)*phi(a + b) ...
    + (1 + (a+b)^2)*Q(a + b) + (1 + (a - b)^2)*(1 - Q(a - b));

f = 0;

for i = 1 : n
    if S(i)
        f = f + Eg(v1(i), v2(i));
    end
end





