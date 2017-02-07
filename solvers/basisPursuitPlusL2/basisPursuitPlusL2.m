function [x_opt, k] = basisPursuitPlusL2(n, m, b, beta, w, mode, arg1, arg2)
% Solves
%                  minimize    ||x||_1 + beta*||x-w||^2
%                       x
%                  subject to  b = Ax
%
% where ||.|| is the L2 norm,  b : m  x 1,  A: m x n,  w: n x 1, and 
% beta >0. There should hold m > n. We use ADMM to solve the problem.
%
% IF mode = 1, arg1 is matrix A
%              arg2 is pinv(A)
%
% IF mode = 2, arg1 is a function handler to A*x
%              arg2 is a function handler to A'*y


% =========================================
% Check input

if length(b) ~= m  || length(w) ~= n || beta <= 0
    error('Dimensions of the input are not correct');
end

if mode == 1
    A = arg1;
    A_pinv = arg2;    
elseif mode == 2
    A  = arg1;
    AT = arg2;
    AAT = @(x) A(AT(x));
    aux1 = zeros(m,1);
else
    error('Mode not recognized');
end 
% =========================================

% =========================================
% Parameters
MAX_ITER = 5000;
rho = 1;
tau_rho = 10;
mu_rho = 2;

eps_prim = 1e-3;
eps_dual = 1e-3;
% =========================================

% =========================================
% Initialization and precomputations
lambda = zeros(n,1);

x = w;
y = w;
% =========================================
    
rho_plus_beta_inv = 1/(rho + beta);

for k = 1 : MAX_ITER
	
	% ************************************
	% x-minimization
    
    v = lambda - rho*y;
    
    beta_w_v = beta*w - v;
    
    ind_pos = (beta_w_v >  1);       % Indices where x is positive
    ind_neg = (beta_w_v < -1);       % Indices where x is negative
    
    x = zeros(n,1);
    
    x(ind_pos) = rho_plus_beta_inv*( beta_w_v(ind_pos) - 1 );
    x(ind_neg) = rho_plus_beta_inv*( beta_w_v(ind_neg) + 1 );
        
	% ************************************
    
    % ************************************
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
	% ************************************


	% ************************************
    % Update dual variable

    r_prim = x - y;      % primal residual
    lambda = lambda + rho*r_prim;

    s_dual = -rho*(y - y_prev); % dual residual
	% ************************************

    % **********************************************************
    % rho adjustment
   
    r_prim_norm = norm(r_prim);
    s_dual_norm = norm(s_dual);
        
    if r_prim_norm > tau_rho*s_dual_norm
        rho = mu_rho*rho;
    elseif s_dual_norm > tau_rho*r_prim_norm
        rho = rho/mu_rho;
    end
    % **********************************************************
    
    if r_prim_norm < eps_prim && s_dual_norm < eps_dual
        break;
    end            
end

x_opt = x;


fprintf('Iterations: %d\n', k);

if k >= MAX_ITER
    fprintf('Warning: Maximum number of iterations reached\n');
end
end


function [x] = conjgrad(A, b, x)

% Implements the conjugate gradient method.
%
% A is a function handler

MAX_ITER = 1e6;
TOL = 1e-10;

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








	
