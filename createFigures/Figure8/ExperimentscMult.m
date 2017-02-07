% Checks how the reconstruction rate of
%
%   minimize   ||x||_1  +  ||x - w||_1
%      x
%   subject to Ax = b
%
% and
%
%   minimize   ||x||_1  +  ||x - w||_2^2
%      x
%   subject to Ax = b
%
% compare with the values predicted theoretically, when the prior 
% information is improved by a multiplicative factor, i.e., 
% w <- c*w, for different values of c.
%
% WARNING: beta = 1 for both L1-L1 and L1-L2

%%
% =========================================================================
% Parameters of the experiment

n = 1000;              % Dimension of the vector

eps_error = 1e-2;      % Target error

c_vec = [1:0.2:1.8, 2:0.5:7];

card_x = 70;          % Cardinality of x

var = 0.3;            % Variance of noise in side info

beta = 1;

MAX_ITERATIONS = 5000;   % Maximum number of iterations of algorithms

card_comm_factor = 0.8;
card_rest_factor = 0.7;

% To save data
FILENAME = 'Exp_Mult';

% Random seed number
seed = RandStream('mcg16807','Seed',31415);
RandStream.setDefaultStream(seed);    

SHOW_PLOTS = 1;

SAVE_DATA = 1;

% Paths to solvers and aux functions
addpath('../../solvers/basisPursuitPlusL1/')
addpath('../../solvers/basisPursuitPlusL2/')
addpath('../../solvers/otherApproaches/modifiedCS-Vaswani/')
addpath('../auxFunctions/sharperL1L1AndL1L2bounds/')
addpath('../auxFunctions/sharperCSbound/')
% =========================================================================


% =========================================================================
% Generate data

% Generate x
x_aux = [randn(card_x,1); zeros(n-card_x,1)];
permutation_x = randperm(n);
x = x_aux(permutation_x);

fprintf('Dimensions of x  = %d\n', n);
fprintf('Cardinality of x = %d (%d perc)\n', card_x, card_x/n*100);

% Generate w
card_i_common = round(card_comm_factor*card_x); % Common components between x and i
card_i_rest   = round(card_rest_factor*card_x); % Not necessarily common

i_aux = [var*randn(card_i_rest,1); zeros(n-card_i_rest,1)];
permutation_rest = randperm(n);
i = i_aux(permutation_rest);
vec_aux = [var*randn(card_i_common,1); zeros(n-card_i_common,1)];
vec_perm = vec_aux(permutation_x);
i = i + vec_perm;
w_orig = x + i;
 
fprintf('||x-w||/||x|| = %f\n', norm(x - w_orig)/norm(x));
% =========================================================================

% =========================================================================
% Compute best bounds: theory

% Simple CS limit
CS_lim = ceil(2*card_x*log(n/card_x) + (7/5)*card_x) + 1;

[CS_lim_sharp, arg_m] = sharperCSBound(n, card_x); 

fprintf('Simple CS limit:         %d\n', CS_lim);
fprintf('Simple CS limit (sharp): %d\n', ceil(CS_lim_sharp));

% *********************************************************************
% Compute theoretical bounds

num_c            = length(c_vec);
h_bar_vec        = zeros(num_c,1);
bound_L1         = zeros(num_c,1);
bound_L1_sharper = zeros(num_c,1);
bound_L2         = zeros(num_c,1);
bound_L2_sharper = zeros(num_c,1);

for ind_c = 1 : num_c

    c = c_vec(ind_c);
    w = c*w_orig;

    % -----------------------------------------------------------------
    % L1-L1 bounds
    h     = sum( (x>0).*(x<w) ) + sum( (x<0).*(x>w) );
    h_bar = sum( (x>0).*(x>w) ) + sum( (x<0).*(x<w) );
    s     = card_x;
    l     = sum(w ~= x);
    IcJc = sum( (x==0).*(x==w) );
    q = n - IcJc;

    h_bar_vec(ind_c) = h_bar; 

    qhh_b = q + h + h_bar;

    % Sharper bound (using Q-function)
    [bound_L1_sharper(ind_c), t_dummy] = sharperL1L1Bound(x, w, beta);

    bound_L1(ind_c) = 2*h_bar*log(2*n/qhh_b) + (7/10)*qhh_b;
    % -----------------------------------------------------------------
        

    % -----------------------------------------------------------------
    % L1-L2 bounds

    Ic     = (x ~= 0);
    I_plus = (x  > 0);
    I_minu = (x  < 0);
    IJc    = logical((x ~= 0).*(x == w));
    IcJ    = logical((x == 0).*(x ~= w));
    
    r2 = n/q;
        
    K_neq  = logical(IcJ.*(abs(w) >  1/beta));
    K_eq   = logical(IcJ.*(abs(w) == 1/beta));
    
    card_K_neq = sum(K_neq);
    card_K_eq  = sum(K_eq);
    
    % w_bar
    [aux_w_bar,ind_aux] = max(abs( abs(w(IcJ)) - 1/beta ));
    indices_aux = find(IcJ);
    w_bar = abs(w(indices_aux(ind_aux)));
    
    v = sum( (1 + beta*(x(I_plus) - w(I_plus))).^2 ) ...
        + sum( (1 - beta*(x(I_minu) - w(I_minu))).^2) ...
        + sum( (beta*abs(w(K_neq)) - 1).^2 );

    L1 = abs(1 - beta*w_bar)*exp(2*beta*w_bar*log(r2)*(0.5*w_bar*beta - 1));
    
    q_s_over_n_q = (q - s)/(n-q);
    
    L2 = abs(1 - beta*w_bar)*exp(4*(beta*w_bar-2)*beta*w_bar/(1-beta*w_bar)^2*log(q/s));
    
    cond1 = 0;  % flag = 1 if condition 1 is satisfied
    cond2 = 0;  % flag = 1 if condition 2 is satisfied
    
    if q_s_over_n_q <= L1
        
        cond1 = 1;
        
        b1 = 2*v*log(r2) + s + card_K_neq + 0.5*card_K_eq + (4/5)*q;
    end
    
    if q_s_over_n_q >= L2
        
        cond2 = 1;
        
        b2 = 2*v/(1 - beta*w_bar)^2*log(q/s) + card_K_neq + 0.5*card_K_eq + (9/5)*s;
    end
    
    if cond1 == 0 && cond2 == 0
        bound_L2(ind_c) = 0;                  
    elseif cond1 == 1 && cond2 == 0
        bound_L2(ind_c) = b1;                 
    elseif cond1 == 0 && cond2 == 1
        bound_L2(ind_c) = b2;                 
    else %both exist
        bound_L2(ind_c) = min(b1, b2);        
    end

    % Sharper bound (using Q-function)
    [bound_L2_sharper(ind_c), t_dummy] = sharperL1L2Bound(x, w, beta);
    bound_L2_sharper(ind_c) = ceil(bound_L2_sharper(ind_c));
    % -----------------------------------------------------------------
end
% ********************************************************************* 
% =========================================================================
 
%%
% =========================================================================
% Experiments

fprintf('||x-w||/||x|| = %f\n', norm(x - w)/norm(x));

min_measurementsL1 = zeros(num_c, 1);
iterations_ADMML1  = zeros(num_c, 1);
min_measurementsL2 = zeros(num_c, 1);
iterations_ADMML2  = zeros(num_c, 1);
    
% Full measurement matrix; we'll use its rows sequentially until we get the
% target error
A = randn(n,n);
b = A*x;

% ------------------------------------------------------------------
% L1-L1

fprintf('\n\nL1-L1 experiments\n');
for ind_c = 1 : num_c
    
    fprintf('Experiment number = %d (out of %d)\n', ind_c, num_c);
    
    c = c_vec(ind_c);

    w = c*w_orig;
            
    for measurements = 1 : n
        
        A_trial = A(1:measurements , :);
        b_trial = b(1:measurements);
        
        A_trial_pinv = pinv(A_trial);
        
        [x_tr, iter] = basisPursuitPlusL1(b_trial, w, beta, 1, ...
            A_trial, A_trial_pinv, MAX_ITERATIONS);
        
        if norm(x_tr - x)/norm(x) <= eps_error
            min_measurementsL1(ind_c) = measurements;
            iterations_ADMML1(ind_c)  = iter;
            break;
        end
    end        
end
% ------------------------------------------------------------------

% ------------------------------------------------------------------
% L1-L2

fprintf('\n\nL1-L2 experiments\n');
for ind_c = 1 : num_c
    
    fprintf('Experiment number = %d (out of %d)\n', ind_c, num_c);
    
    c = c_vec(ind_c);

    w = c*w_orig;
            
    for measurements = 1 : n
        
        A_trial = A(1:measurements , :);
        b_trial = b(1:measurements);
        
        A_trial_pinv = pinv(A_trial);
        
        [x_tr, iter] = basisPursuitPlusL2_BB(n, measurements, b_trial, ...
            beta, w, 1, A_trial, A_trial_pinv, MAX_ITERATIONS);
        
        if norm(x_tr - x)/norm(x) <= eps_error
            min_measurementsL2(ind_c) = measurements;
            iterations_ADMML2(ind_c)  = iter;
            break;
        end
    end        
end
% ------------------------------------------------------------------

% ------------------------------------------------------------------
% Mod-CS
    
fprintf('\n\nMod-CS experiments\n');
supp_w = (w~=0);

for measurements = 1 : n
    
    A_trial = A(1:measurements , :);
    b_trial = b(1:measurements);
    
    A_trial_pinv = pinv(A_trial);
    
    [x_MO, k2] = ModifiedCS(b_trial, supp_w, 1, A_trial, A_trial_pinv, MAX_ITERATIONS);
    
    if norm(x_MO - x)/norm(x) <= eps_error
        min_measurements_MO = measurements;
        iterations_ADMM_MO  = k2;
        break;
    end
end
% ------------------------------------------------------------------
% =========================================================================


%%
% =========================================================================
% Plot results

if SHOW_PLOTS == 1

    figure(1);clf;    
    plot(c_vec, min_measurementsL1, 'bo-');
    hold on;
    plot(c_vec, min_measurementsL2, 'ro-');
    plot(c_vec, min_measurements_MO*ones(1, num_c), 'go-');
    plot(c_vec, bound_L1,         'co-.');
    plot(c_vec, bound_L1_sharper, 'ko-.');
    plot(c_vec, bound_L2,         'mo-.');
    plot(c_vec, bound_L2_sharper, 'ko-.');
    xlabel('c');
    ylabel('Minimum number of measurements');
    ylim([0,1.1*max(max(min_measurementsL2), max(bound_L1))])
    %xlim([2.9,4.1])
    
    legend('L1-L1', 'L1-L2', 'Mod-CS', 'boundL1', 'boundL1 (sharper)', ...
      'boundL2', 'boundL2 (sharper)')
    
    if SAVE_DATA == 1
        filename_im = [FILENAME, '.pdf'];
        saveas(gcf, filename_im);
    end
end

% =========================================================================

%%
% =========================================================================
% Save data

if SAVE_DATA == 1
    
    filename_data = [FILENAME, '.mat'];
    
    save(FILENAME, 'n', 'eps_error', 'c_vec', 'x', 'w_orig', 'card_x', 'beta', ...
        'card_comm_factor', 'card_rest_factor', 'seed', 'CS_lim','CS_lim_sharp',  ...
        'MAX_ITERATIONS', 'min_measurements_MO', 'iterations_ADMM_MO', 'num_c', ...
        'min_measurementsL1', 'iterations_ADMML1', 'bound_L1', 'bound_L1_sharper', ...
        'min_measurementsL2', 'iterations_ADMML2', 'bound_L2', 'bound_L2_sharper');
end
% =========================================================================

