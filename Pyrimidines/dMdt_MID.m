function dMdt = dMdt_MID(t, M, v, c, isotopomer_frac, time_points, ...
    slope, slope_metabs, M_metab_list)

dMdt = zeros(size(M, 1), 1);

%% Assign the flux parameters
ump_de_novo = v(1);
uridine_ump = v(2);
ump_uridine = v(3);
cytidine_uridine = v(4);
ump_out = v(5);
uridine_out = v(6);

%% Assign the concentration parameters
c_uridine = c(1);
c_ump = c(2);

%% Assign aspartate isotopomer fractions
% from supplementary data (S27) of comprehensive isotopomer analysis of 
% glutamate and aspartate in small tissues samples
f_aspartate_M1 = isotopomer_frac(1);
f_aspartate_M2 = isotopomer_frac(2);
f_aspartate_M3 = isotopomer_frac(3);

%% Assign the prior MID values

M_r5p = M(ismember(M_metab_list,'R5P'));
M_aspartate = M(ismember(M_metab_list,'ASPARTATE'));
M_ump = M(ismember(M_metab_list,'UMP'));
M_uridine = M(ismember(M_metab_list,'URIDINE'));

% comment/uncomment if you run it for GBM or cortex based on the following
% info:
% cytidine is labeled in GBM 
M_cytidine = M(ismember(M_metab_list, 'CYTIDINE')); 
% cytidine is unlabeled in cortex
% M_cytidine = [1; 0; 0; 0; 0; 0];

%% Assign the slope of the labeled input metabolites
index = 1;
for time_slot = 2:(length(time_points)-1)
    if(t > time_points(time_slot))
        index = index+1;
    end
end

dMdt(ismember(M_metab_list,'R5P')) = slope(ismember(slope_metabs,'R5P'),index);
dMdt(ismember(M_metab_list,'ASPARTATE')) = slope(ismember(slope_metabs,'ASPARTATE'),index);
% comment the following line if you run it for cortex
dMdt(ismember(M_metab_list, 'CYTIDINE')) = slope(ismember(slope_metabs, 'CYTIDINE'),index);

%% ODEs for balanced metabolites

% UMP
N_aspartate = zeros(4, 1);
N_aspartate(1) = M_aspartate(1) + f_aspartate_M1 * M_aspartate(2);
N_aspartate(2) = (1 - f_aspartate_M1) * M_aspartate(2) + f_aspartate_M2 * M_aspartate(3);
N_aspartate(3) = (1 - f_aspartate_M2) * M_aspartate(3) + f_aspartate_M3 * M_aspartate(4);
N_aspartate(4) = (1 - f_aspartate_M3) * M_aspartate(4) + M_aspartate(5);

M_aspartate_r5p = sum_diag(N_aspartate * M_r5p');

dMdt(ismember(M_metab_list,'UMP')) = (ump_de_novo * M_aspartate_r5p(1:6) + ...
    uridine_ump  * M_uridine(1:6) - (ump_out + ump_uridine) * M_ump(1:6)) / c_ump;
% Uridine
dMdt(ismember(M_metab_list,'URIDINE')) = (cytidine_uridine * M_cytidine + ...
    ump_uridine  * M_ump(1:6) - (uridine_out + uridine_ump) * M_uridine(1:6)) / c_uridine;

end

function sum = sum_diag(N)

[m,n] = size(N);
sum = zeros(m+n-1,1);
for i = 1:m
    for j = 1:n
        sum(i+j-1) = sum(i+j-1)+ N(i,j);
    end
end

end