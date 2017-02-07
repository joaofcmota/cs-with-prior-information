% Checks how the reconstruction rate of
%
%   minimize   ||x||_1  +  beta*||x - w||_2^2
%      x
%   subject to Ax = b
%
% compares with the values predicted theoretically. 

%%
% =========================================================================
% Parameters of the experiment

n = 500;              % Dimension of the vector

eps_error = 1e-2;      % Target error

MONTE_CARLO = 5;

beta_vec = [1e-2; 5e-2;1e-1;5e-1;7.5e-1;9e-1;1e0;2.5;5;10;50;100];

card_x = 50;          % Cardinality of x

var = 0.8;            % Variance of noise in side info

card_comm_factor = 0.4*0.8;
card_rest_factor = 0.4*0.2;

% To save data
FILENAME = 'ExpBetaL2';

% Random seed number
seed = RandStream('mcg16807','Seed',1);
RandStream.setDefaultStream(seed);    

SHOW_PLOTS = 1;

SAVE_DATA = 1;

addpath('../../solvers/basisPursuitPlusL2/')
addpath('../auxFunctions/sharperL1L1AndL1L2bounds/')
addpath('../auxFunctions/sharperCSbound/')
% =========================================================================

%%
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
w = x + i;
          
fprintf('||x-w||/||x|| = %f\n', norm(x - w)/norm(x));
% =========================================================================

%%
% =========================================================================
% Compute best bounds: theory


% Simple CS limit
CS_lim = ceil(2*card_x*log(n/card_x) + (7/5)*card_x) + 1;
[CS_lim_sharp, arg_m] = sharperCSBound(n, card_x); 

fprintf('Simple CS limit:         %d\n', CS_lim);
fprintf('Simple CS limit (sharp): %d\n', ceil(CS_lim_sharp));

% *********************************************************************
% L1-L2 bounds

h     = sum( (x>0).*(x<w) ) + sum( (x<0).*(x>w) );
h_bar = sum( (x>0).*(x>w) ) + sum( (x<0).*(x<w) );
s     = card_x;
l     = sum(w ~= x);
IcJc = sum( (x==0).*(x==w) );
q = n - IcJc;
Ic     = (x ~= 0);
I_plus = (x  > 0);
I_minu = (x  < 0);
IJc    = logical((x ~= 0).*(x == w));
IcJ    = logical((x == 0).*(x ~= w));

fprintf('h     = %d\nh_bar = %d\ns     = %d\nl     = %d\nq     = %d\n', h, h_bar, s, l, q);
    
r2 = n/q;

num_betas = length(beta_vec);
bound_L2  = zeros(num_betas,1);
bound_L2_sharper = zeros(num_betas,1);

for ind_beta = 1 : num_betas
    
    beta = beta_vec(ind_beta);
    
    % Sharper bound (using Q-function)
    [bound_L2_sharper(ind_beta), t_dummy] = sharperL1L2Bound(x, w, beta);

    
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
        bound_L2(ind_beta) = 0;                  
    elseif cond1 == 1 && cond2 == 0
        bound_L2(ind_beta) = b1;                 
    elseif cond1 == 0 && cond2 == 1
        bound_L2(ind_beta) = b2;                 
    else %both exist
        bound_L2(ind_beta) = min(b1, b2);        
    end
    
end
% ********************************************************************* 
 
% =========================================================================
 
%%
% =========================================================================
% Experiments

fprintf('||x-w||/||x|| = %f\n', norm(x - w)/norm(x));

min_measurements = zeros(num_betas, MONTE_CARLO);
iterations_ADMM  = zeros(num_betas, MONTE_CARLO);

for mc = 1 : MONTE_CARLO
    
    fprintf('MC number = %d (out of %d)\n', mc, MONTE_CARLO);
    
    % Full measurement matrix; we'll use its rows sequentially until we get the
    % target error
    A = randn(n,n);
    b = A*x;
    
    for ind_beta = 1 : num_betas
        
        fprintf('Experiment number = %d (out of %d)\n', ind_beta, num_betas);
        
        beta = beta_vec(ind_beta);
                
        for measurements = 1 : n
            
            A_trial = A(1:measurements , :);
            b_trial = b(1:measurements);
            
            A_trial_pinv = pinv(A_trial);
            
            [x_tr, iter] = basisPursuitPlusL2_BB(n, measurements, b_trial, ...
                beta, w, 1, A_trial, A_trial_pinv);
            
            if norm(x_tr - x)/norm(x) <= eps_error
                min_measurements(ind_beta, mc) = measurements;
                iterations_ADMM(ind_beta, mc) = iter;
                break;
            end
        end        
    end
end
% =========================================================================


%%
% =========================================================================
% Plot results

if SHOW_PLOTS == 1
    
    color_map = hsv(MONTE_CARLO);             % Colormap HSV
    figure(1);clf;    
    for mc = 1 : MONTE_CARLO
        semilogx(beta_vec, min_measurements(:,mc), 'o-', 'color', color_map(mc,:));
        hold on;
    end
    xlabel('beta');
    ylabel('Minimum number of measurements');
    ylim([0,1.1*max(max(min_measurements))])
    xlim([1e-2,1e2])
    
    semilogx(beta_vec, bound_L2,         'co-.')
    semilogx(beta_vec, bound_L2_sharper, 'ko-.')
    
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
    
    save(FILENAME, 'n', 'eps_error', 'beta_vec', 'x', 'w', 'card_x', 'MONTE_CARLO', ...
        'card_comm_factor', 'card_rest_factor', 'seed', 'CS_lim', 'CS_lim_sharp', ...
        'num_betas', 'min_measurements', 'iterations_ADMM', 'bound_L2', 'bound_L2_sharper');    
end
