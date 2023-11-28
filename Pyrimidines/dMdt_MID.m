function dMdt = dMdt_MID(t,M,v,c,isotopomer_frac, time_points,slope, slope_metabs, M_metab_list)

dMdt = zeros(size(M,1),1);

%% Assign the flux parameters
ump_de_novo = v(1);
%uridine_out = v(2);
%cytidine_uridine = v(3);
uridine_ump = v(2);
ump_out = v(3);

%% Assign the concentration parameters
%c_uridine = c(2);
c_ump = c(1);

%% Assign aspartate isotopomer fractions
% media-3.xlsx S27
% from supplementary data of comprehensive isotopomer analysis of glutamate 
% and aspartate in small tissues samples
% f_aspartate_M1 = 0.25;
% f_aspartate_M2 = 0.49;
% f_aspartate_M3 = isotopomer_frac(1);
f_aspartate_M1 = isotopomer_frac(1);
f_aspartate_M2 = isotopomer_frac(2);
f_aspartate_M3 = isotopomer_frac(3);

%% Assign the prior MID values

M_r5p = M(ismember(M_metab_list,'R5P'));
M_aspartate = M(ismember(M_metab_list,'ASPARTATE'));
%M_uracil = M(ismember(M_metab_list,'URACIL'));
M_ump = M(ismember(M_metab_list,'UMP'));
M_uridine = M(ismember(M_metab_list,'URIDINE'));
%M_cytidine = M(ismember(slope_metabs, 'CYTIDINE'));

%% Assign the slope of the labeled input metabolites
index = 1;
for time_slot = 2:(length(time_points)-1)
    if(t > time_points(time_slot))
        index = index+1;
    end
end

dMdt(ismember(M_metab_list,'R5P')) = slope(ismember(slope_metabs,'R5P'),index);
dMdt(ismember(M_metab_list,'ASPARTATE')) = slope(ismember(slope_metabs,'ASPARTATE'),index);
dMdt(ismember(M_metab_list,'URIDINE')) = slope(ismember(slope_metabs,'URIDINE'),index);
%dMdt(ismember(M_metab_list,'URACIL')) = slope(ismember(slope_metabs,'URACIL'),index);
%dMdt(ismember(M_metab_list, 'CYTIDINE')) = slope(ismember(slope_metabs, 'CYTIDINE'),index);

%% ODEs for balanced metabolites

% UMP
N_aspartate = zeros(4, 1);
N_aspartate(1) = M_aspartate(1) + f_aspartate_M1 * M_aspartate(2);
N_aspartate(2) = (1 - f_aspartate_M1) * M_aspartate(2) + f_aspartate_M2 * M_aspartate(3);
N_aspartate(3) = (1 - f_aspartate_M2) * M_aspartate(3) + f_aspartate_M3 * M_aspartate(4);
N_aspartate(4) = (1 - f_aspartate_M3) * M_aspartate(4) + M_aspartate(5);

M_aspartate_r5p = sum_diag(N_aspartate * M_r5p');

% dMdt_ump = dMdt(ismember(M_metab_list,'UMP'));
dMdt(ismember(M_metab_list,'UMP')) = (ump_de_novo * M_aspartate_r5p(1:6) + ...
    uridine_ump  * M_uridine(1:6) - ump_out * M_ump(1:6)) / c_ump;
% dMdt_ump(7:9) = [0; 0; 0];
% dMdt(ismember(M_metab_list,'UMP')) = dMdt_ump;

% URIDINE
%M_uracil_r5p = sum_diag(M_uracil * M_r5p');
% dMdt_uridine = dMdt(ismember(M_metab_list,'URIDINE'));
%dMdt(ismember(M_metab_list,'URIDINE')) = (cytidine_uridine * M_cytidine(1:6) + ... %uracil_uridine * M_uracil_r5p(1:6)
%     - uridine_ump * M_uridine(1:6) - uridine_out * M_uridine(1:6)) / c_uridine;
% dMdt_uridine(7:9) = [0; 0; 0];
% dMdt(ismember(M_metab_list,'URIDINE')) = dMdt_uridine;

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