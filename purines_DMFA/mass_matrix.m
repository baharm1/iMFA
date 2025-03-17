function mass = mass_matrix(t, y, n_c_param, M_metab_list)

%% Assign the concenrtaions
c = y(1:n_c_param);
c_amp = c(1);
c_gdp = c(2);
c_gmp = c(3);
c_guanosine = c(4);
c_imp = c(5);
c_inosine = c(6);

%% Create the vector that will form the diagonal of the mass matrix

iden = ['c_amp'; 'c_gdp'; 'c_gmp'; 'c_guanosine'; 'c_imp'; 'c_inosine'; 
    M_metab_list];

d = ones(numel(y), 1);
d(ismember(iden,'AMP')) = c_amp;
d(ismember(iden,'GDP')) = c_gdp;
d(ismember(iden,'GMP')) = c_gmp;
d(ismember(iden,'GUANOSINE')) = c_guanosine;
d(ismember(iden,'IMP')) = c_imp;
d(ismember(iden,'INOSINE')) = c_inosine;

mass = diag(d);

%% Add the non-diagonal elements

mass(ismember(iden,'IMP'), ismember(iden,'c_imp')) = y(ismember(iden,'IMP'));
mass(ismember(iden,'GMP'), ismember(iden,'c_gmp')) = y(ismember(iden,'GMP'));
mass(ismember(iden,'GDP'), ismember(iden,'c_gdp')) = y(ismember(iden,'GDP'));
mass(ismember(iden,'AMP'), ismember(iden,'c_amp')) = y(ismember(iden,'AMP'));
mass(ismember(iden,'INOSINE'), ismember(iden,'c_inosine')) = y(ismember(iden,'INOSINE'));
mass(ismember(iden,'GUANOSINE'), ismember(iden,'c_guanosine')) = y(ismember(iden,'GUANOSINE'));


end