
%%
% =========================================================================
% Experimental performance of the algorithms: file Experiments.m
FilenameCS   ='ResultsLaTexCS.dat';
FilenameSI   ='ResultsLaTexSI.dat';
FilenameSIL2 ='ResultsLaTexSIL2.dat';

num_meas = length(measurements);

fid_CS   = fopen(FilenameCS,   'w');
fid_SI   = fopen(FilenameSI,   'w');
fid_SIL2 = fopen(FilenameSIL2, 'w');

for i = 1 : num_meas
    
    fprintf(fid_CS,   '%d %3.2f\n', measurements(i), results_av_CS(i));
    fprintf(fid_SI,   '%d %3.2f\n', measurements(i), results_av_L1(i));
    fprintf(fid_SIL2, '%d %3.2f\n', measurements(i), results_av_L2(i));
    
end

fclose(fid_CS);
fclose(fid_SI);
fclose(fid_SIL2);
% =========================================================================



%%
% =========================================================================
% Experiments L1-L1 for varying beta: file ExperimentsBeta.m

Filename1   ='BetaResults1.dat';
Filename2   ='BetaResults2.dat';
Filename3   ='BetaResults3.dat';
Filename4   ='BetaResults4.dat';
Filename5   ='BetaResults5.dat';

fid1   = fopen(Filename1,   'w');
fid2   = fopen(Filename2,   'w');
fid3   = fopen(Filename3,   'w');
fid4   = fopen(Filename4,   'w');
fid5   = fopen(Filename5,   'w');

for ind_beta = 1 : num_betas
    fprintf(fid1,   '%3.3f %d\n', log10(beta_vec(ind_beta))+2, min_measurements(ind_beta, 1));
    fprintf(fid2,   '%3.3f %d\n', log10(beta_vec(ind_beta))+2, min_measurements(ind_beta, 2));
    fprintf(fid3,   '%3.3f %d\n', log10(beta_vec(ind_beta))+2, min_measurements(ind_beta, 3));
    fprintf(fid4,   '%3.3f %d\n', log10(beta_vec(ind_beta))+2, min_measurements(ind_beta, 4));
    fprintf(fid5,   '%3.3f %d\n', log10(beta_vec(ind_beta))+2, min_measurements(ind_beta, 5));
end

fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);
