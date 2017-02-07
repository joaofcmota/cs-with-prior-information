function [x_opt, k] = basisPursuitPlusL2_BB(n, m, b, beta, w, mode, A, AT, ...
    varargin)

% [x_opt, k] = basisPursuitPlusL2_BB(n, m, b, beta, w, mode, A, AT, varargin)
%
% Solves the optimization problem
% 
%             minimize    |x||_1 + 0.5*beta*||x - w||^2       (1)  
%             subject to  Ax = b
%
% where ||.||_1 is the L1-norm and ||.|| is the L2-norm. The matrix A has 
% dimensions m x n, vector b has dimensions m x 1, and vector w has 
% dimensions n x 1. This problem is solved using the Barzilai-Borwein 
% algorithm. 
%
% IF mode = 1, A is an explicit matrix
%              AT is ignored
%
% If mode = 2, A  is a function handler to A*x
%              AT is a function handler to AT*x
%
% The optional argument is the maximum number of iterations.
%
%
% x_opt is the solution to (1) and k is the number of iterations.


% =========================================================================
% Check input

if length(b) ~= m  || length(w) ~= n || beta <= 0
    error('Dimensions of the input are not correct');
end

if nargin > 9
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

M = 20;
gamma = 1e-4;
epsilon = 1e-10;
sigma = 0.3;
alpha = 1;

lambda0 = zeros(m, 1);
% =========================================================================

% =========================================================================
% Algorithm

lambda = lambda0;

[f_lambda, g_lambda, x] = feval(@f_and_g, lambda, n, beta, w, mode, b, A, AT);

f_r = -inf*ones(M,1);
f_r(1) = f_lambda;
f_ref = f_lambda;
aux = 1;

counter = 0;

for k = 1 : MAX_ITER
    if (alpha <= epsilon) || (alpha >= 1/epsilon)
        if norm(g_lambda) > 1
            delta = 1;
        elseif (norm(g_lambda) >= 1e-5) && (norm(g_lambda) <= 1)
            delta = 1/norm(g_lambda);
        else % norm(g) < 1e-5
            delta = 1e5;
        end;
                     
        alpha = delta;
    end;
    
    eta = 1/alpha;
    
%    f_ref = max(f_r);
    if counter >= M
        f_ref = max(f_r);
        counter = 0;
    else if f_lambda > f_ref
            f_ref = f_lambda;
            counter = 0;
        end
    end
    counter = counter + 1;     
    
    slope = gamma*(g_lambda'*g_lambda);
    iter = 0;
    [f_lambda_eta_g, g_dummy, x] = feval(@f_and_g, lambda-eta*g_lambda, n, beta, w, mode, b, A, AT);    
    while f_lambda_eta_g > f_ref - eta*slope
        eta = sigma*eta;
        [f_lambda_eta_g, g_dummy, x] = feval(@f_and_g, lambda-eta*g_lambda, n, beta, w, mode, b, A, AT);
        iter = iter + 1;
    end
    
    
    lambda_new = lambda - eta*g_lambda;    
    [f_new, g_new, x] = feval(@f_and_g,lambda_new, n, beta, w, mode, b, A, AT);
    
    y = g_new-g_lambda;
    
    alpha = -g_lambda'*y/(eta*(g_lambda'*g_lambda));    
    
    lambda = lambda_new;
    f_lambda = f_new;
    g_lambda = g_new;
    
    if norm(g_lambda) <= 1e-6        
        break;
    end
           
    aux = mod(aux,M)+1;
    
    f_r(aux) = f_lambda;      
    
end

x_opt = x;

%fprintf('Iterations: %d\n', k);

%if k == MAX_ITER
%    fprintf('Warning: maximum number of iterations reached in BB\n');
%end

end

function [f_lambda, g_lambda, x] = f_and_g(lambda, n, beta, w, mode, b, A, AT)

if mode == 1
    AT_lambda = A'*lambda;
elseif mode == 2
    AT_lambda = AT(lambda);
else
    error('Unrecognized mode')
end

v = -beta*w - AT_lambda;

x = zeros(n,1);

one_over_beta = 1/beta;

ind_pos = v < -1;
ind_neg = v >  1;

x(ind_pos) = -one_over_beta* (v(ind_pos) + 1);
x(ind_neg) = -one_over_beta* (v(ind_neg) - 1);

if mode == 1
    Ax = A*x;
else
    Ax = A(x);
end

g_lambda = Ax - b;

f_lambda = lambda'*g_lambda - norm(x,1) - 0.5*beta*norm(x-w)^2;

end


