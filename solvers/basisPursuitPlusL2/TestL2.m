%% 
% =========================================================================
% Generate Data

n = 1000;
s = 200;

beta = 1;

m = round(2*s*log(n/s) + (7/5)*s) + 1;

% Generate x
x_aux = [randn(s,1) ; zeros(n-s,1)];
perm  = randperm(n);
x = x_aux(perm);

% Generate side information w
var = 1;
card_i_common = round(0.4*0.8*s); % Common components between x and i
card_i_rest   = round(0.4*0.2*s); % Not necessarily common
    
i_aux = [var*randn(card_i_rest,1); zeros(n-card_i_rest,1)];
permutation_rest = randperm(n);
i = i_aux(permutation_rest);
vec_aux = [var*randn(card_i_common,1); zeros(n-card_i_common,1)];
vec_perm = vec_aux(perm);
i = i + vec_perm;
w = x + i;


% Generate matrix A and vector b
A = randn(m,n);
b = A*x;

fprintf('n = %d, m = %d\n', n, m);
fprintf('||x - w||/||x|| = %f\n', norm(x - w)/norm(x));
% =========================================================================


%%
% =========================================================================
% Solve the problem

cvx_begin
    variable x_cvx(n);
    minimize( norm(x_cvx,1) + 0.5*beta*square_pos(norm(x_cvx - w)) );
    subject to
        A*x_cvx == b;
cvx_end

t_ADMM_aux = cputime; 
[x_ADMM, k_ADMM] = basisPursuitPlusL2(n, m, b, beta, w, 1, A, pinv(A));
t_ADMM = cputime - t_ADMM_aux;

t_BB_aux = cputime; 
[x_BB,   k_BB]   = basisPursuitPlusL2_BB(n, m, b, beta, w, 1, A, []);
t_BB = cputime - t_BB_aux;
fprintf('ADMM: error = %f, time = %f\n', norm(x_ADMM - x_cvx)/norm(x_cvx), t_ADMM);
fprintf('BB:   error = %f, time = %f\n', norm(x_BB   - x_cvx)/norm(x_cvx), t_BB);
% =========================================================================


