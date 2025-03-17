%% This file starts with the solutions from the INST-MFA code
% and changes the fluxes to fit a dynamic MFA model

%% Save the timestamp and create a folder with the name

global timestamp folder
timestamp = ['timestamp_', datestr(datetime,'yyyymmdd_HHMM')];
folder = ['output_files\', char(timestamp)];
mkdir(folder);

%% User-defined inputs

% Name of file with model information
model_file = 'purine_model';

% Name of file with the experimental data
input_data_file = 'Input-Data\Input_brain_rt_conc_abund.xlsx'; % brain
%input_data_file = 'Input-Data\Input_gbm_rt_conc_abund.xlsx'; % GBM  

% Minimum standard deviation
min_sd = sqrt(10^(-5));

% Number of initial guess for optimization
niter = 100;

% B-spline properties
% Order n implies a polynomial function of order n-1
bspline_order = 3;
internal_knots = [0.65]; % brain
%internal_knots = [0.3]; % GBM

no_internal_knots = numel(internal_knots);

% Upper bound for flux
ub_flux = 2000;

% Parameter for l2 normalization
l2 = 0;

%% Read input files and pre-process the data

% Read the model information
[S, rxns, input_metabs, balance_metabs] = import_model(model_file);
disp('Model imported');

% Read the experimental MID data
[time_points, input_data_metabs, MID_input, SD_input, ...
    balance_data_metabs, MID_balance, SD_balance, ...
    conc_ratio, conc_sd, conc_metab_list, init_flux, init_conc, ...
    data_i_list] = import_expt_data_dmfa(input_data_file, min_sd, ...
    input_metabs, balance_metabs);
disp('experimental data imported');

%% Define Parameter bounds

% Model parameters
n_c_param = length(balance_metabs);
n_flux_param = size(S, 2) * (no_internal_knots + bspline_order);
n_param = n_c_param + n_flux_param;
c_param_index = [1:n_c_param];
flux_param_index = [n_c_param + 1:n_param];

% Parameter bounds
lb = 0.1 * ones(n_param, 1);
ub = ub_flux * ones(n_param, 1);

% Set the bounds for initial concentration and fluxes
lb(c_param_index) = init_conc(:, 3); 
ub(c_param_index) = init_conc(:, 4); 
lb(n_c_param+1:n_c_param+size(S, 2)) = init_flux(:, 1);
ub(n_c_param+1:n_c_param+size(S, 2)) = init_flux(:, 2);

% Add the bounds for the input MIDs
input_mid_array = reshape(MID_input(:, 2:end), [], 1);
input_mid_array = input_mid_array(~isnan(input_mid_array));
input_sd_array = reshape(SD_input(:, 2:end), [], 1);
input_sd_array = input_sd_array(~isnan(input_sd_array));
lb_values = input_mid_array - input_sd_array;
lb_values(lb_values < 0) = 0;
ub_values = input_mid_array + input_sd_array;
ub_values(ub_values > 1) = 1; 
lb = [lb; lb_values];
ub = [ub; ub_values];

disp('Parameter bounds defined');

%% Generate stoichiometrically balanced initial flux vectors

Aeq = S;
beq = zeros(size(S, 1), 1);
lbf = init_flux(:, 1);
ubf = init_flux(:, 2);
flux0 = repmat(lbf, 1, niter) + repmat(ubf - lbf, 1, niter) .* rand(numel(lbf), niter);

% Stoichiometrically balance the initial flux values
flux_bal = zeros(size(flux0));
for i = 2:niter
    flux_bal(:, i) = fmincon(@(x)sumsqr((x - flux0(:, i)) ./ flux0(:, i)), ...
        flux0(:, i), [], [], Aeq, beq, lbf, ubf);
end

% Add the SS solution as a starting point
flux_bal(:, 1) = init_flux(:, 3);

flux0 = flux_bal;
clear flux_bal Aeq beq lbf ubf

%% Generate initial guess and bounds for concentration

lbc = init_conc(:, 3);
ubc = init_conc(:, 4);
conc0 =  repmat(lbc, 1, niter) + repmat(ubc-lbc, 1, niter) .* rand(numel(lbc), niter);
% Add the SS solution as a starting point
conc0(:, 1) = init_conc(:, 5);

%% Define initial guess

x0 = [conc0; repmat(flux0, bspline_order + no_internal_knots, 1)];
% Set the input MID parameters to the experimental values
x0 = [x0; repmat(input_mid_array, 1, niter)];

%% Define Aeq and beq

% The sum of input MIDs should be 1
ip_metabs = unique(input_data_metabs, 'stable');
cons_mtx = zeros(numel(ip_metabs), numel(input_data_metabs));
for k = 1:numel(ip_metabs)
    cons_mtx(k, ismember(input_data_metabs, ip_metabs(k))) = 1;
end
ip_mid_cons = [];
for t = 2:numel(time_points)
    temp = cons_mtx * MID_input(:, t);
    ip_mid_cons= blkdiag(ip_mid_cons, cons_mtx(~isnan(temp), :));
    clear temp
end
Aeq = [zeros(size(ip_mid_cons, 1), n_param), ip_mid_cons];
beq = ones(size(ip_mid_cons, 1), 1);

disp('Linear constraint matrices created');

%% Define knot position array

knot_pos = ones(1, 2 * bspline_order + no_internal_knots); 
knot_pos(1:bspline_order) = 0;
knot_pos(bspline_order + 1:(bspline_order + no_internal_knots)) = internal_knots;

save([folder,'\data']); 

%% Optimization

tspan = [time_points(1):0.02:time_points(end)];

for i = 1:niter

disp(['Starting optimization for iteration ',num2str(i)]);
 
[xopt(:,i), fval(i), exitf(i)] = fit_expt_data(x0(:, i), tspan, time_points, ...
                        MID_input, input_data_metabs, SD_input, ...
                        MID_balance, balance_data_metabs, SD_balance,...
                        n_c_param, c_param_index, n_param, conc_ratio, ...
                        conc_sd, conc_metab_list, init_conc, S, knot_pos, ...
                        bspline_order, balance_metabs, flux_param_index, ...
                        lb, ub, Aeq, beq, l2);

     save([folder,'\data']); 

end

%% Calculate the dof and chi square cutoff

nparameters = numel(lb);
knownmid = sum(~isnan(MID_input(:, 2:end)), 'all') + sum(~isnan(MID_balance(:, 2:end)), 'all');
knownconc = sum(~isnan(init_conc(:, 1))) + sum(~isnan(conc_ratio(:, 2:end)), 'all');
equalcons = numel(beq);
dof = knownmid + knownconc + equalcons - nparameters;
alpha = 0.05;
chil = chi2inv(alpha / 2, dof);
chiu = chi2inv(1 - alpha / 2, dof);

%% Plot the best solution
[~, min_index] = min(fval);
xopt_min = xopt(:, min_index);

save([folder, '\data']);
[~] = plot_results(xopt_min, tspan, time_points, ...
                        MID_input, input_data_metabs, SD_input, ...
                        MID_balance, balance_data_metabs, SD_balance,...
                        n_c_param, c_param_index, n_param, conc_ratio, ...
                        conc_sd, conc_metab_list, S, knot_pos, ...
                        bspline_order, balance_metabs, ...
                        flux_param_index, data_i_list, rxns);


%%  Estimate the 95% CI bounds for the parameters associated with flux

thres = min(fval) + chi2inv(0.95, 1);

cil = zeros(1, n_flux_param);
ciu = cil;
for i = flux_param_index(1):flux_param_index(end)

    [cil(i-6), ciu(i-6)] = estimate_param_bounds_dmfa(i, thres, xopt_min, ...
        Aeq, beq, tspan, time_points, n_c_param, MID_balance, ...
        balance_data_metabs, SD_balance, MID_input, input_data_metabs, ...
        SD_input, lb, ub, c_param_index, n_param, conc_ratio, conc_sd, ...
        conc_metab_list, init_conc, S, knot_pos, bspline_order, ...
        balance_metabs, flux_param_index, l2);

end
save([folder,'\data']); 
%% Read the stored parameter bounds

for i = flux_param_index(1):flux_param_index(end)
    data = importdata([folder, '\parameter_bounds_', num2str(i), '.mat']);
    cil(i-6) = data.cil;
    ciu(i-6) = data.ciu;
end

%% Save optimized flux vectors between the parameter confidence intervals

mkdir([folder,'\feasible_vectors']);
for i = flux_param_index(1):flux_param_index(end)
    
   [~] = generate_CI_vectors(i, xopt_min, ciu(i-6), cil(i-6), Aeq, beq, ...
       tspan, time_points, n_c_param, MID_balance, balance_data_metabs, ...
       SD_balance, MID_input, input_data_metabs, SD_input, lb, ub, ...
       c_param_index, n_param, conc_ratio, conc_sd, conc_metab_list, ...
       init_conc, S, knot_pos, bspline_order, balance_metabs, ...
       flux_param_index, l2);

end

%% Read the optimized flux vectors and plot the time course confidence intervals

[flux_cil, flux_ciu, flux_optim] = estimate_flux_CI(xopt_min, tspan, S, ...
    knot_pos, bspline_order, flux_param_index);
save([folder, '\data']);

mkdir([folder, '\pdfs']);
for i = 1:numel(rxns)
    
    hold all
    patch([tspan fliplr(tspan)], [flux_cil(i, :) fliplr(flux_ciu(i, :))], ...
        validatecolor("#FD866E"), 'FaceAlpha', 0.5, 'EdgeColor', 'none')
    plot(tspan, flux_optim(i, :), "Color", validatecolor("#EE4928"), ...
        "LineWidth", 1);
    hold off
    xlim([0 tspan(end)]);
    ylim([0 inf]);
    title(rxns(i));
    xlabel('time (h)', 'FontSize', 18);
    ylabel('Flux (pmol/hr.mg-tissue)', 'FontSize', 18);
    filename = [folder, '\plots\', 'flux_CI', num2str(i), '.png'];
    saveas(gcf, filename);
    filename2 = [folder, '\pdfs\', 'flux', num2str(i), '.pdf'];
    saveas(gcf, filename2);
    close
      
end
