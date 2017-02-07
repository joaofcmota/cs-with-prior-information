function [x_opt, k] = ModifiedCS(b, T, mode, arg1, arg2, varargin)

% [x_opt, k] = ModifiedCS(b, T, mode, arg1, arg2, varargin)
%
% Solves
%                  minimize    ||x_{Tc}||_1                   (1)
%                     x
%                  subject to  Ax = b
%
% where b: m x 1,  A: m x n, and Tc denotes the complement of the set
% T, which is a subset of {1,2,...,n}. There should hold m > n. We use ADMM 
% to solve the problem, as explained in the attached documentation. This 
% function has two modes: in mode 1, it uses an explicit matrix A; in 
% mode 2, access to A, namely A*x and A'*y, is done implicitly through 
% function calls.
%
% If mode == 1, arg1 is the matrix A
%               arg2 is the pseudo inverse of A: pinv(A)
%
% If mode == 2, arg1 is a function handler to A*x
%               arg2 is a function handler to A'*y
%
% Inputs:
%   - b: m x 1 vector
%   - T: n x 1 binary vector
%   - mode: either 1 or 2
%   - arg1: in mode 1, an m x n matrix; in mode 2, a function handler
%   - arg2: in mode 2, an n x m matrix; in mode 2, a function handler 
%
% Optional input: maximum number of iterations (default: 5000)
%
% Outputs:
%   - x_opt: solution of (1)
%   - k: number of ADMM iterations
%
% This code was designed and implemented by J. Mota to perform experiments 
% described in
%
% [1] J. Mota, N. Deligiannis, M. Rodrigues, 
%     Compressed Sensing with Prior Information: Optimal Strategies, 
%     Geometry, and Bounds 
%     submitted to IEEE Transactions on Information Theory
%     preprint: http://arxiv.org/abs/1408.5250
%     2014
%
% Problem (1) was proposed and studied in
%
% [2] N. Vaswani, W. Lu,
%     Modified-CS: Modifying Compressive Sensing for Problems With 
%     Partially Known Support
%     IEEE Transactions on Signal Processing, Vol. 58, No. 9
%     2010
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
%
% =========================================================================
% Please send any questions, comments, or bug reports to j.mota@ucl.ac.uk
% =========================================================================


% =========================================================================
% Check inputs

m = length(b);
n = length(T);

% Check if T is binary
if sum(or(T == ones(n,1),  T == zeros(n,1))) ~= n
    error('T must be a binary vector of size n')
end

if mode == 1
    A = arg1;
    A_pinv = arg2;
    
    if size(A,1) ~= m || size(A,2) ~= n || size(A_pinv,1) ~= n ...
            || size(A_pinv,2) ~= m
        error('Dimensions of the input are not correct');
    end
    
elseif mode == 2
    A  = arg1;
    AT = arg2;
    AAT = @(x) A(AT(x));
    aux1 = zeros(m,1);          % Auxiliary vector for initialization of CG
else
    error('Mode not recognized');
end 

if nargin > 7
    error('There should be only one optional argument')
end

if ~isempty(varargin)
    MAX_ITER = varargin{1};
else
    MAX_ITER = 5000;
end
% =========================================================================

% =========================================================================
% Parameters

rho = 1;
tau_rho = 10;
mu_rho = 2;

eps_prim = 1e-3;
eps_dual = 1e-3;
% =========================================================================

% =========================================================================
% Initialization and precomputations
lambda = zeros(n,1);

x = T;
y = T;
% =========================================================================
             

for k = 1 : MAX_ITER
	
	% *********************************************************************
	% x-minimization
    
    v = lambda - rho*y;

    x = zeros(n,1);
    
    rho_inv = 1/rho;
    
    x(T) = -rho_inv*v(T);
    
    ind_v_less_than_minus1 = logical((1-T).*(v < -1));
    ind_v_larger_than_1    = logical((1-T).*(v >  1));
        
    x(ind_v_less_than_minus1) = -rho_inv*(v(ind_v_less_than_minus1) + 1);
            
    x(ind_v_larger_than_1)    = -rho_inv*(v(ind_v_larger_than_1)    - 1);
	% *********************************************************************
    
    % *********************************************************************
	% y-minimization
	
    y_prev = y;

    z = (lambda + rho*x)/rho;
    if mode == 1
        y = z - A_pinv*(A*z - b);
    else
        Az_minus_b = A(z) - b;
        aux1 = conjgrad(AAT, Az_minus_b, aux1);
        y = z - AT(aux1);
    end
	% *********************************************************************


	% *********************************************************************
    % Update dual variable

    r_prim = x - y;      % primal residual
    lambda = lambda + rho*r_prim;

    s_dual = -rho*(y - y_prev); % dual residual
	% *********************************************************************

    % *********************************************************************
    % rho adjustment
   
    r_prim_norm = norm(r_prim);
    s_dual_norm = norm(s_dual);
   
    if r_prim_norm > tau_rho*s_dual_norm
        rho = mu_rho*rho;
    elseif s_dual_norm > tau_rho*r_prim_norm
        rho = rho/mu_rho;
    end
    % *********************************************************************
    
    if r_prim_norm < eps_prim && s_dual_norm < eps_dual
        break;
    end    

    %fprintf('Iteration: %d\n', k);
end

x_opt = x;

%fprintf('Iterations: %d\n', k);

if k >= MAX_ITER
    fprintf('Warning (ModifiedCS): Maximum number of iterations reached. Primal residual = %f, Dual residual = %f\n', ...
        r_prim_norm, s_dual_norm);
end

end


function [x] = conjgrad(A, b, x)

% Implements the conjugate gradient method.
%
% A is a function handler

MAX_ITER = 1e6;
TOL      = 1e-10;

r = b - A(x);
p = r;
rsold = r'*r;

for i = 1 : MAX_ITER
    
    Ap = A(p);
    alpha = rsold/(p'*Ap);
    x = x + alpha*p;
    r = r - alpha*Ap;
    rsnew = r'*r;
    
    if sqrt(rsnew) < TOL
        break;
    end
    
    p = r + rsnew/rsold*p;
    rsold = rsnew;
end

if i == MAX_ITER
    fprintf('Conjugate gradient: Maximum number of iteration reached\n');
end

end








	
