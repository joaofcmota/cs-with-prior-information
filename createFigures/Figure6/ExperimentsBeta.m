% Checks how the reconstruction rate of
%
%   minimize   ||x||_1  +  beta*||x - w||_1
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
FILENAME = 'ExpBetaL1';

MAX_ITERATIONS = 5000;

% Random seed number
seed = RandStream('mcg16807','Seed',1);
RandStream.setDefaultStream(seed);    

SHOW_PLOTS = 1;

SAVE_DATA = 1;

% Paths to solvers and aux functions
addpath('../../solvers/basisPursuitPlusL1/')
addpath('../../solvers/otherApproaches/modifiedCS-Vaswani/')
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
% L1-L1 bounds

h     = sum( (x>0).*(x<w) ) + sum( (x<0).*(x>w) );
h_bar = sum( (x>0).*(x>w) ) + sum( (x<0).*(x<w) );
s     = card_x;
l     = sum(w ~= x);
IcJc = sum( (x==0).*(x==w) );
q = n - IcJc;

fprintf('h     = %d\nh_bar = %d\ns     = %d\nl     = %d\nq     = %d\n', h, h_bar, s, l, q);


num_betas        = length(beta_vec);
bound_L1         = zeros(num_betas,1);
bound_L1_sharper = zeros(num_betas,1);


qhh_b = q + h + h_bar;

for ind_beta = 1 : num_betas
    beta = beta_vec(ind_beta);
    
    % Sharper bound (using Q-function)
    [bound_L1_sharper(ind_beta), t_dummy] = sharperL1L1Bound(x, w, beta);

    if beta == 1
        
        bound_L1(ind_beta) = 2*h_bar*log(2*n/qhh_b) + (7/10)*qhh_b;
        
    elseif beta < 1
        
        b1_exists = 0;   % flag that is = to 1 when b1 exists
        b2_exists = 0;   % flag that is = to 1 when b2 exists
        
        if (q-s)/(2*n - qhh_b) <= (1-beta)/(1+beta)*(qhh_b/(2*n))^(4*beta/(beta+1)^2)
            
            b1_exists = 1;
            
            b_aux_1 = 2*(h_bar + (s - h_bar)*(1-beta)^2/(1+beta)^2)...
                *log(2*n/qhh_b) + s + (2/5)*qhh_b;            
        end
            
        if q > s && (q-s)/(2*n - qhh_b) >= (1-beta)/(1+beta)*(s/q)^(4*beta/(1-beta)^2)
            
            b2_exists = 1;
            
            b_aux_2 = 2*(h_bar*(1+beta)^2/(1-beta)^2 + s - h_bar)*log(q/s) + (7/5)*s;            
        end
                
        if b1_exists == 0 && b2_exists == 0
            bound_L1(ind_beta) = 0;            
        elseif b1_exists == 1 && b2_exists == 0
            bound_L1(ind_beta) = b_aux_1;            
        elseif b1_exists == 0 && b2_exists == 1
            bound_L1(ind_beta) = b_aux_2;
        else %both exist            
            bound_L1(ind_beta) = min(b_aux_1, b_aux_2);
        end
        
    elseif beta > 1
        
        b1_exists = 0;   % flag that is = to 1 when b1 exists
        b2_exists = 0;   % flag that is = to 1 when b2 exists
        
        if (s - (h + h_bar))/(2*n - qhh_b) <= (beta-1)/(beta+1)*(qhh_b/(2*n))^(4*beta/(beta+1)^2)
            
            b1_exists = 1;
            
            b_aux_1 = 2*(h_bar + (q + h - s)*(beta-1)^2/(beta+1)^2)...
                *log(2*n/qhh_b) + l + (2/5)*qhh_b;
        end
            
        if h + h_bar > 0 && (s - (h + h_bar))/(2*n - qhh_b) >= (beta-1)/(beta+1)...
                *((h+h_bar)/s)^(4*beta/(beta-1)^2)
           
            b2_exists = 1;
            
            b_aux_2 = 2*(h_bar*(beta+1)^2/(beta-1)^2 + q + h - s)...
                *log(s/(h + h_bar)) + l + (2/5)*(h + h_bar);        
        end
        
        if b1_exists == 0 && b2_exists == 0
            bound_L1(ind_beta) = 0;            
        elseif b1_exists == 1 && b2_exists == 0
            bound_L1(ind_beta) = b_aux_1;            
        elseif b1_exists == 0 && b2_exists == 1
            bound_L1(ind_beta) = b_aux_2;
        else %both exist            
            bound_L1(ind_beta) = min(b_aux_1, b_aux_2);
        end
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

min_measurements_MO = zeros(MONTE_CARLO, 1);   % Mod-CS

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
            
            [x_tr, iter] = basisPursuitPlusL1(b_trial, w, beta, 1, ...
                A_trial, A_trial_pinv, MAX_ITERATIONS);
            
            if norm(x_tr - x)/norm(x) <= eps_error
                min_measurements(ind_beta, mc) = measurements;
                iterations_ADMM(ind_beta, mc) = iter;
                break;
            end
        end        
    end

    % ===========================
    % Same experiments for Mod-CS
        
    supp_w = (w~=0);

    for measurements = 1 : n
        
        A_trial = A(1:measurements , :);
        b_trial = b(1:measurements);
        
        A_trial_pinv = pinv(A_trial);
        
        [x_MO, k2] = ModifiedCS(b_trial, supp_w, 1, A_trial, A_trial_pinv, MAX_ITERATIONS);
        
        if norm(x_MO - x)/norm(x) <= eps_error
            min_measurements_MO(mc) = measurements;
            iterations_ADMM_MO  = k2;
            break;
        end
    end
    % ===========================
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
    
    plot([1;1], [0;1.1*max(max(min_measurements))], 'r');
    
    semilogx(beta_vec, bound_L1,         'co-.');
    semilogx(beta_vec, bound_L1_sharper, 'ko-.');
    ylim([0,400]);
    
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
        'card_comm_factor', 'card_rest_factor', 'seed', 'CS_lim','CS_lim_sharp',  ...
        'num_betas', 'min_measurements', 'iterations_ADMM', 'bound_L1', 'bound_L1_sharper', ...
        'MAX_ITERATIONS', 'min_measurements_MO');
    
end
