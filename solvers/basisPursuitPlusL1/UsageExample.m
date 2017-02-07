% Usage example of basisPursuitPlusL1.m
%
% =========================================================================
% basisPursuitPlusL1: minimizing L1 + L1 with linear constraints
% Copyright (C) 2014  Joao Mota
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
% Define parameters to create data

% Dimensions of matrix A: m x n
n = 500;
m = 60;

card_x = 10;                 % cardinality of x

% To generate prior information
card_i_common = 2;           % number of components in common with x
card_i_rest   = 4;           % number of components not in the support of x
% =========================================================================


% =========================================================================
% Generate data

% Generate x
x_aux = [randn(card_x,1); zeros(n-card_x,1)];
permutation_x = randperm(n);
x = x_aux(permutation_x);

% Generate w
i_aux = [randn(card_i_rest,1); zeros(n-card_i_rest,1)];
permutation_rest = randperm(n);
i = i_aux(permutation_rest);
vec_aux = [randn(card_i_common,1); zeros(n-card_i_common,1)];
vec_perm = vec_aux(permutation_x);
i = i + vec_perm;
w = x + i;
% =========================================================================


% =========================================================================
% Mode 1 (matrix A is declared explicitly)

A = randn(m,n);
b = A*x;

[x_opt_mode1, k1] = basisPursuitPlusL1(b, w, 1, 1, A, pinv(A));
% =========================================================================


% =========================================================================
% Mode 2 (matrix A is declared implicitly)

% If Sparco toolbox is installed, A can be declared, e.g., as
%
% A = opGaussian(m, n, 1);
%
% and b as
%
% b = A(x, 1);
%
% The handlers would be
%
% A_handler  = @(x) A(x,1);
% AT_handler = @(x) A(x,2);
%
% Since Sparco toolbox may not be installed, let us use A from mode 1

A_handler  = @(x) A*x;
AT_handler = @(y) A'*y;

b = A_handler(x);
[x_opt_mode2, k2] = basisPursuitPlusL1(b, w, 1, 2, A_handler, AT_handler);
% =========================================================================

% =========================================================================
% Assess the error

% % If cvx is installed, we can check the error of the output:
% cvx_begin quiet
%     cvx_precision high
%     variable x_cvx(n);
%     minimize(norm(x_cvx,1) + norm(x_cvx - w,1));
%     subject to
%         A*x_cvx == b;
% cvx_end
% 
% fprintf('Error mode 1 = %f\n', norm(x_opt_mode1 - x_cvx)/norm(x_cvx));
% fprintf('Error mode 2 = %f\n', norm(x_opt_mode2 - x_cvx)/norm(x_cvx));
% =========================================================================
