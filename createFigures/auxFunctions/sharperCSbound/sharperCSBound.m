function [f_val, arg_val] = sharperCSBound(n, s)

% [f_val, arg_val] = sharperCSBound(n, s)
% 
% Computes the Gaussian distance of the L1-norm numerically without
% using approximations of the Q-function and finds the optimal 
% parameter. See the documentation attached.
%
% Input:
%   - n: dimension of the vector space
%   - s: sparsity of the (optimal) vector; there must hold n > s
%
% Output:
%   - f_val:   exact Gaussian distance
%   - arg_val: value of the optimal parameter
%
% Method:
%   Newton-Raphson method
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


if n <= s
  error('There must hold n > s')
end

% ====================================================================== 
% Parameters of algorithm
 
MAX_ITER = 50;        % Typically, 30 iterations are enough

tol = 1e-7;           % Algorithm stops if |(x_{k+1} - x_k)/x_k | < tol

epsilon = 1e-14;      % We don't divide by a number smaller than this
% ====================================================================== 


% ====================================================================== 
% Initializations

% Constants frequently used
n_s = n - s;
sqrt_2_pi = sqrt(2/pi);

% Q-function
Q = @(x) 0.5*erfc(x/sqrt(2));

% Function we want to minimize
h = @(t) s + s*t^2 + 2*n_s*(1 + t^2)*Q(t) - sqrt_2_pi*n_s*t*exp(-0.5*t^2);

% Derivative
dh = @(t) 2*s*t + 4*n_s*t*Q(t) - 2*n_s*sqrt_2_pi*exp(-0.5*t^2);

% Second-order derivative
ddh = @(t) 2*s + 4*n_s*Q(t);

% Initialization point
x = sqrt(2*log(n/s));
% ====================================================================== 

% ====================================================================== 
% Algorithm

solutionFound = false; 

for i = 1 : MAX_ITER

    dh_x  = dh(x);
    ddh_x = ddh(x);

    if(abs(ddh_x) < epsilon)         % Denominator too small
        break;                  
    end

    x_prev = x;

    x = x - dh_x/ddh_x; 
    
    % Stopping criterion
    if(abs(x - x_prev)/abs(x_prev) < tol)           
        solutionFound = true;
        break;
    end
end
% ====================================================================== 


% ====================================================================== 
% Output
if ~solutionFound || i == MAX_ITER
  fprintf('Warning: Solution not found by Newton-Raphson method\n')
end

f_val   = h(x);
arg_val = x;
% ====================================================================== 


