function dMdt = dMdt_MID(t,M,v,c, time_points,slope, slope_metabs, M_metab_list)

dMdt = zeros(size(M,1),1);

%% Assign the flux parameters

imp_in = v(1);
imp_inosine = v(2);
hypoxanthine_imp = v(3);
hypoxanthine_inosine = v(4);
inosine_out = v(5);
imp_gmp = v(6);
guanine_gmp = v(7);
gmp_guanosine = v(8);
guanine_guanosine = v(9);
guanosine_out = v(10);
gmp_gdp = v(11);
gdp_out = v(12);
imp_amp = v(13);
adenosine_amp = v(14);
amp_out = v(15);

%% Assign the concentration parameters
c_amp = c(1);
c_gdp = c(2);
c_gmp = c(3);
c_guanosine = c(4);
c_imp = c(5);
c_inosine = c(6);

%% Assign the prior MID values
M_r5p = M(ismember(M_metab_list,'R5P'));
M_glycine = M(ismember(M_metab_list,'GLYCINE'));
M_methf = M(ismember(M_metab_list,'METHF'));
M_imp = M(ismember(M_metab_list,'IMP'));
M_inosine = M(ismember(M_metab_list,'INOSINE'));
M_gdp = M(ismember(M_metab_list,'GDP'));
M_guanosine = M(ismember(M_metab_list,'GUANOSINE'));
M_gmp = M(ismember(M_metab_list,'GMP'));
M_amp = M(ismember(M_metab_list,'AMP'));

%% Assign the slope of the labeled input metabolites
index = 1;
for time_slot = 2:(length(time_points)-1)
    if(t > time_points(time_slot))
        index = index+1;
    end
end

dMdt(ismember(M_metab_list,'R5P')) = slope(ismember(slope_metabs,'R5P'),index);
dMdt(ismember(M_metab_list,'GLYCINE')) = slope(ismember(slope_metabs,'GLYCINE'),index);
dMdt(ismember(M_metab_list,'METHF')) = slope(ismember(slope_metabs,'METHF'),index);

%% ODEs for balanced metabolites


% IMP
imp_de_novo = sum_diag(M_methf*M_methf');
imp_de_novo = sum_diag(imp_de_novo*M_glycine');
imp_de_novo = sum_diag(imp_de_novo*M_r5p');

dMdt(ismember(M_metab_list,'IMP')) = (-(imp_inosine+imp_gmp+imp_amp)*M_imp + imp_in*imp_de_novo(1:6) + hypoxanthine_imp*M_r5p)/c_imp;
%dMdt(ismember(M_metab_list,'IMP')) = (-(imp_inosine+imp_gmp+imp_amp)*M_imp + imp_in*imp_de_novo([1:4,6]) + hypoxanthine_imp*M_r5p([1:4,6]))/c_imp;

% Inosine
dMdt(ismember(M_metab_list,'INOSINE')) = (-inosine_out*M_inosine + imp_inosine*M_imp + hypoxanthine_inosine*M_r5p)/c_inosine;
%dMdt(ismember(M_metab_list,'INOSINE')) = (-inosine_out*M_inosine + imp_inosine*M_imp + hypoxanthine_inosine*M_r5p([1:4,6]))/c_inosine;

% GMP
dMdt(ismember(M_metab_list,'GMP')) = (-(gmp_gdp+gmp_guanosine)*M_gmp + imp_gmp*M_imp + guanine_gmp*M_r5p)/c_gmp;
%dMdt(ismember(M_metab_list,'GMP')) = (-(gmp_gdp+gmp_guanosine)*M_gmp + imp_gmp*M_imp(1:4) + guanine_gmp*M_r5p(1:4))/c_gmp;

% GDP
dMdt(ismember(M_metab_list,'GDP')) = (-gdp_out*M_gdp + gmp_gdp*M_gmp)/c_gdp;

% Guanosine
dMdt(ismember(M_metab_list,'GUANOSINE')) = (-guanosine_out*M_guanosine + guanine_guanosine*M_r5p + gmp_guanosine*M_gmp)/c_guanosine;
%dMdt(ismember(M_metab_list,'GUANOSINE')) = (-guanosine_out*M_guanosine + guanine_guanosine*M_r5p(1:4) + gmp_guanosine*M_gmp)/c_guanosine;

% AMP
dMdt(ismember(M_metab_list,'AMP')) = (-amp_out*M_amp + imp_amp*M_imp + adenosine_amp*[1;0;0;0;0;0])/c_amp;
%dMdt(ismember(M_metab_list,'AMP')) = (-amp_out*M_amp + imp_amp*M_imp([1,4]) + adenosine_amp*[1;0])/c_amp;

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