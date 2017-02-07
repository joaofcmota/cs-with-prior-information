load Exp_Add.mat;

%%
% =========================================================================
% Experiments L1-L2 for varying c: file ExperimentsL1L2c.m

FilenameL1L1Perfm      = 'L1L1Perfm.dat';
FilenameL1L1Bound      = 'BoundsL1L1.dat';
FilenameL1L1BoundSharp = 'BoundsL1L1Sharp.dat';
FilenameL1L2Perfm      = 'L1L2Perfm.dat';
FilenameL1L2Bound      = 'BoundsL1L2.dat';
FilenameL1L2BoundSharp = 'BoundsL1L2Sharp.dat';

fid1   = fopen(FilenameL1L1Perfm     ,   'w');
fid2   = fopen(FilenameL1L1Bound     ,   'w');
fid3   = fopen(FilenameL1L1BoundSharp,   'w');
fid4   = fopen(FilenameL1L2Perfm     ,   'w');
fid5   = fopen(FilenameL1L2Bound     ,   'w');
fid6   = fopen(FilenameL1L2BoundSharp,   'w');

for ind_c = 1 : 7;%num_c
    fprintf(fid1,   '%2.1f %d\n', c_vec(ind_c), min_measurementsL1(ind_c));
    fprintf(fid2,   '%2.1f %d\n', c_vec(ind_c), ceil(bound_L1(ind_c)));
    fprintf(fid3,   '%2.1f %d\n', c_vec(ind_c), ceil(bound_L1_sharper(ind_c)));
    fprintf(fid4,   '%2.1f %d\n', c_vec(ind_c), min_measurementsL2(ind_c));
    fprintf(fid5,   '%2.1f %d\n', c_vec(ind_c), ceil(bound_L2(ind_c)));
    fprintf(fid6,   '%2.1f %d\n', c_vec(ind_c), ceil(bound_L2_sharper(ind_c)));
end

fclose(fid1);
fclose(fid2);
fclose(fid3);
