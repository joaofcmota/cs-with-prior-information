% Checks how the reconstruction rate of
%
%   minimize   ||x||_1  +  ||x - w||_1
%      x
%   subject to Ax = b,
%
% where w is given, varies with the number of measurements

%%
% =========================================================================
% Parameters of the experiment

n = 1000;              % Dimension of the vector

measurements = 1:20:660;

MONTE_CARLO = 100;     % Experiments per point: we'll take an average

beta = 1;

card_x = 70;          % Cardinality of x

var = 0.8;            % Variance of noise in side info

MAX_ITERATIONS = 5000;   % Maximum number of iterations in solvers

PriorInfoImprovement_L1L1 = 10;
PriorInfoImprovement_L1L2 = 0;

% Side information parameters
card_comm_factor = 0.4*0.8;
card_rest_factor = 0.4*0.2;

% To save data
FILENAME = 'fig1';
  
% Random seed number
seed = RandStream('mcg16807','Seed',31415);
RandStream.setDefaultStream(seed);    

eps_error = 1e-2;     % Relative error is below eps_error = Success

SHOW_PLOTS = 1;

SAVE_DATA = 1;

% Paths to solver and aux functions
addpath('../../solvers/basisPursuitPlusL1/')
addpath('../../solvers/basisPursuitPlusL2/')
addpath('../../solvers/otherApproaches/modifiedCS-Vaswani/')
addpath('../auxFunctions/sharperCSbound/')
addpath('../auxFunctions/sharperL1L1AndL1L2bounds/')
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
w = x + i;
 
fprintf('||x-w||/||x|| = %f\n', norm(x - w)/norm(x));
% =========================================================================


%%
% =========================================================================
% Compute bounds for BETA = 1

% Simple CS limit
CS_lim = ceil(2*card_x*log(n/card_x) + (7/5)*card_x);
prob_CS = 1-exp(-0.5*(sqrt(CS_lim) - CS_lim)^2);

% Simple CS limit (sharper)
[CS_lim_sharper, best_t] = sharperCSBound(n, card_x);
CS_lim_sharper = ceil(CS_lim_sharper);

fprintf('Simple CS limit: %d with prob. = %f\n', CS_lim, prob_CS);
fprintf('Simple CS limit: %d (sharper)\n', CS_lim_sharper);

% *********************************************************************
% L1-L1 bounds

w_L1 = w*PriorInfoImprovement_L1L1;

IcJc = sum( (x==0).*(x==w_L1) );
q = n - IcJc;
h     = sum( (x>0).*(x<w_L1) ) + sum( (x<0).*(x>w_L1) );
h_bar = sum( (x>0).*(x>w_L1) ) + sum( (x<0).*(x<w_L1) );

if h_bar == 0
    error('h_bar = 0');
end

r = 2*n/(q + h + h_bar);

m_L1 = 2*h_bar*log(r) + (7/10)*(q + h + h_bar);
m_L1 = ceil(m_L1);

prob_L1 = 1-exp(-0.5*(sqrt(m_L1) - m_L1)^2);

% Sharper bound (computed numerically with the Q-function)
[m_L1_sharp, t_sharp] = sharperL1L1Bound(x, w_L1, 1);
m_L1_sharp = ceil(m_L1_sharp);

fprintf('q     = %d\nh      = %d\nh_bar = %d\nr =     %d\n', q, h, h_bar, r);
fprintf('Measurements L1-L1: %d with prob. = %f\n', m_L1, prob_L1);
fprintf('Measurements L1-L1: %d (sharp bound)\n', m_L1_sharp);
% *********************************************************************


% *********************************************************************
% L1-L2 bounds

w_L2 = w + sign(w)*PriorInfoImprovement_L1L2;

Ic     = (x ~= 0);
I_plus = (x  > 0);
I_minu = (x  < 0);
IJc    = logical((x ~= 0).*(x == w_L2));
IcJ    = logical((x == 0).*(x ~= w_L2));

IcJc = sum( (x==0).*(x==w_L2) );
q = n - IcJc;

v = sum( (1 + x(I_plus) - w_L2(I_plus)).^2 ) ...
    + sum( (1 + w_L2(I_plus) - x(I_plus)).^2) ...
    + sum( (abs(w_L2(IJc)) - 1).^2 );

K = sum( abs(w_L2(IcJ)) >= 1 );

w_bar = max(abs(w_L2(Ic)));

r2 = n/q;

L = abs(1 - w_bar)*exp(2*w_bar*log(r2)*(0.5*w_bar - 1));

q_s_over_n_q = (q - card_x)/(n-q);

fprintf('\n\n');
fprintf('v            = %f\n', v);
fprintf('K            = %f\n', K);
fprintf('w_bar        = %f\n', w_bar);
fprintf('q_s_over_n_q = %f\n', q_s_over_n_q);
fprintf('L            = %f\n', L);

if (q - card_x)/(n-q) > L
   error('Case not addressed in the paper (L1-L2)');
end

m_L2 = 2*v*log(r2) + card_x + 2*K + (4/5)*q;
m_L2 = ceil(m_L2);

% Sharper bound (computed numerically with the Q-function)
[m_L2_sharp, t_sharp_2] = sharperL1L2Bound(x, w_L2, 1);
m_L2_sharp = ceil(m_L2_sharp);

fprintf('Measurements L1-L2: %d\n', m_L2);
fprintf('Measurements L1-L2: %d (sharp bound)\n', m_L2_sharp);
% *********************************************************************
 
% =========================================================================
 
%%
% =========================================================================
% Experiments

num_meas = length(measurements);

fprintf('Number of total experiments = %d x %d = %d\n', num_meas, MONTE_CARLO, num_meas*MONTE_CARLO);

% Each entry: relative error
results_L1 = zeros(num_meas, MONTE_CARLO);  % L1-L1
results_L2 = zeros(num_meas, MONTE_CARLO);  % L1-L2
results_CS = zeros(num_meas, MONTE_CARLO);  % Simple CS
results_MO = zeros(num_meas, MONTE_CARLO);  % Mod-CS

for ind_m = 1 : num_meas
    
    fprintf('Experiment number = %2d (out of %d)\n', ind_m, num_meas);
    
    m = measurements(ind_m);
    
    for mc = 1 : MONTE_CARLO
        
        % Generate matrix A and vector b
        
        A = randn(m,n);
        b = A*x;
        A_pinv = pinv(A);

        % L1-L1
        [x_L1, k1] = basisPursuitPlusL1(b, w_L1, beta, 1, A, A_pinv, MAX_ITERATIONS);
        
        % L1-L2
        [x_L2, k2] = basisPursuitPlusL2_BB(n, m, b, beta, w_L2, 1, A, A_pinv, MAX_ITERATIONS);

        % Mod-CS
        [x_MO, k3] = ModifiedCS(b, (w~=0), 1, A, A_pinv, MAX_ITERATIONS);

        % Solve Simple Basis Pursuit
        opts = spgSetParms('verbosity',0);  % Turn off the SPGL1 log output
        x_SG = spg_bp(A, b, opts);

        results_L1(ind_m, mc) = norm(x_L1 - x,Inf)/norm(x,Inf);
        results_L2(ind_m, mc) = norm(x_L2 - x,Inf)/norm(x,Inf);
        results_CS(ind_m, mc) = norm(x_SG - x,Inf)/norm(x,Inf);        
        results_MO(ind_m, mc) = norm(x_MO - x,Inf)/norm(x,Inf);        
    end
end
% =========================================================================


%%
% =========================================================================
% Plot results

% Take an average of the success rate

successes_L1 = (results_L1 < eps_error);
successes_L2 = (results_L2 < eps_error);
successes_CS = (results_CS < eps_error);
successes_MO = (results_MO < eps_error);

results_av_L1 = zeros(num_meas,1);
results_av_L2 = zeros(num_meas,1);
results_av_CS = zeros(num_meas,1);
results_av_MO = zeros(num_meas,1);

for ind_m = 1 : num_meas    
    results_av_L1(ind_m) = mean(successes_L1(ind_m, :));
    results_av_L2(ind_m) = mean(successes_L2(ind_m, :));
    results_av_CS(ind_m) = mean(successes_CS(ind_m, :));
    results_av_MO(ind_m) = mean(successes_MO(ind_m, :));
end

if SHOW_PLOTS == 1
    % Plots
    max_val = max([max(results_av_L1), max(results_av_CS)]);
    min_val = min([min(results_av_L1), min(results_av_CS)]);
    
    % figure(1);clf;
    % plot(x, 'bo');
    % hold on;
    % plot(w, 'ro');
    % legend('True signal', 'Side Information')
    
    figure(1);clf;
    plot(measurements, results_av_L1, 'bo-');
    hold on;
    plot(measurements, results_av_L2, 'go-');
    plot(measurements, results_av_CS, 'rs-');
    plot(measurements, results_av_MO, 'ks-');
    plot(m_L1*[1;1],           [min_val;max_val], 'b');
    plot(m_L2*[1;1],           [min_val;max_val], 'g');
    plot(CS_lim*[1;1],         [min_val;max_val], 'r');
    plot(CS_lim_sharper*[1;1], [min_val;max_val], 'r-.');
    plot(m_L1_sharp*[1;1],     [min_val;max_val], 'b-.');
    plot(m_L2_sharp*[1;1],     [min_val;max_val], 'g-.');
    xlabel('Measurements');
    ylabel('Success rate');
    legend('L1-L1', 'L1-L2','CS', 'Mod-CS', 'Bound L1-L1', 'Bound L2-L2', 'Bound CS', ...
        'Bound CS (sharper)', 'Bound L1-L1 (sharper)', 'Bound L1-L2 (sharper)');
    
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
    
    save(filename_data, 'n', 'measurements', 'MONTE_CARLO', 'card_x',  ...
        'var', 'seed', 'eps_error', 'results_L1', 'results_CS', 'successes_L1', ...
        'results_L2', 'successes_L2','successes_CS', 'results_av_L1', ...
        'results_av_L2','results_av_CS', 'results_av_MO', 'results_MO', 'successes_MO');
end
% =========================================================================




