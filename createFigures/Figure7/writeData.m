load ExpBetaL2.mat;

%%
% =========================================================================
% Experiments L1-L1 for varying beta: file ExperimentsBeta.m

Filename1   ='BetaResultsL1L1-1.dat';
Filename2   ='BetaResultsL1L1-2.dat';
Filename3   ='BetaResultsL1L1-3.dat';
Filename4   ='BetaResultsL1L1-4.dat';
Filename5   ='BetaResultsL1L1-5.dat';

FilenameL1L2      = 'BoundsL1L2.dat';
FilenameL1L2Sharp = 'BoundsL1L2Sharp.dat';

fid1   = fopen(Filename1,   'w');
fid2   = fopen(Filename2,   'w');
fid3   = fopen(Filename3,   'w');
fid4   = fopen(Filename4,   'w');
fid5   = fopen(Filename5,   'w');

fid6   = fopen(FilenameL1L2,      'w');
fid7   = fopen(FilenameL1L2Sharp, 'w');

for ind_beta = 1 : num_betas
    fprintf(fid1,   '%3.3f %d\n', log10(beta_vec(ind_beta))+2, min_measurements(ind_beta, 1));
    fprintf(fid2,   '%3.3f %d\n', log10(beta_vec(ind_beta))+2, min_measurements(ind_beta, 2));
    fprintf(fid3,   '%3.3f %d\n', log10(beta_vec(ind_beta))+2, min_measurements(ind_beta, 3));
    fprintf(fid4,   '%3.3f %d\n', log10(beta_vec(ind_beta))+2, min_measurements(ind_beta, 4));
    fprintf(fid5,   '%3.3f %d\n', log10(beta_vec(ind_beta))+2, min_measurements(ind_beta, 5));
    
    fprintf(fid6,   '%3.3f %d\n', log10(beta_vec(ind_beta))+2, ceil(bound_L2(ind_beta)));
    fprintf(fid7,   '%3.3f %d\n', log10(beta_vec(ind_beta))+2, ceil(bound_L2_sharper(ind_beta)));
end

fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);
fclose(fid6);
fclose(fid7);