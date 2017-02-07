% =========================================================================
% Parameters of the experiment

n = 1000;              % Dimension of the vector

eps_error = 1e-2;      % Target error

MONTE_CARLO = 5;

%beta_vec = [1e-2;5e-2;1e-1;2e-1;5e-1;7e-1;1e0;3e0;5e0;1e1;5e1;1e2;5e2;1e3;5e3];
beta_vec = 2.5;%[1e-2; 5e-2;1e-1;5e-1;7.5e-1;9e-1;1e0;2.5;5;10;50;100];

card_x = 40;          % Cardinality of x

var = 0.1;            % Variance of noise in side info

card_comm_factor = 0.4*0.5;
card_rest_factor = 0.4*0.5;

% To save data
FILENAME = 'ExpBetaL2';

% Random seed number
seed = RandStream('mcg16807','Seed',1);
RandStream.setDefaultStream(seed);    

SHOW_PLOTS = 1;

SAVE_DATA = 1;

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

fprintf('Simple CS limit: %d\n', CS_lim);

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

for ind_beta = 1 : num_betas
    
    beta = beta_vec(ind_beta);
    
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
 
fprintf('min L2 = %d\n', min(round(bound_L2)));

