%Save the timestamp and create a folder with the name
clear;clc;
global timestamp folder
timestamp = ['timestamp_',datestr(datetime,'yyyymmdd_HHMM')];
folder = ['output_files\',char(timestamp)];
mkdir(folder);

%% User-defined inputs

% File containing model reactions and other information
model_file = 'pyrimidine_model_v2';

% File containing experimental data
input_data_file = 'Input-Data\Input_gbm_ctrl_norm_gbm_150PPM_20230620_hippoM5_v2.xlsx';
%input_data_file = 'Input-Data\Input_brain_ctrl_norm_cortex_15ppm_230620_hippo_v2.xlsx';

% Minimum standard deviation for experimental data
min_sd = sqrt(10^(-5));

% Number of cores to use for knitro optimizations. This depends on the
% numebr of lisences available.
ncores = 1;

% Confidence level for the chi-square goodness of fit test
alpha = 0.05;

% Confidence level for parameter intervals
cl = 0.95;

%% Read the input files and process the data

% Read the model information
[S,rxns,input_metabs,balance_metabs] = import_model(model_file);
disp('Model imported');

% Read the experimental MID data
[time_points, input_data_metabs, MID_input, SD_input,...
    balance_data_metabs, MID_balance, SD_balance,...
    conc_metab,conc_values] = import_expt_data(input_data_file,min_sd,input_metabs,balance_metabs);
disp('experimental data imported');

%% Create the linear constraint matrix

n_c_param = length(balance_metabs);
n_flux_param = size(S,2);
n_input_param = sum(~isnan(MID_input(:,2:end)),'all');
n_isotopomer_frac = 3;

Aeq = [S,zeros(size(S,1),n_c_param+n_isotopomer_frac)];
beq = zeros(size(S,1),1);

% Add rows and columns for the input metabolite MIDs
% The input MID values should be 1
ip_metabs = unique(input_data_metabs,'stable');
cons_mtx = zeros(numel(ip_metabs),numel(input_data_metabs));
for k = 1:numel(ip_metabs)
    cons_mtx(k,ismember(input_data_metabs,ip_metabs(k))) = 1;
end
ip_mid_cons = [];
for t = 2:numel(time_points)
    temp = cons_mtx*MID_input(:,t);
    ip_mid_cons= blkdiag(ip_mid_cons,cons_mtx(~isnan(temp),:));
    clear temp
end
Aeq = blkdiag(Aeq,ip_mid_cons);
beq = [beq;ones(size(ip_mid_cons,1),1)];

disp('Linear constraint matrices created');


%% Set the ub and lb for the conc. and fluxes

lb = 0.01 * ones(n_flux_param+n_c_param,1);
lb(n_flux_param + 1: n_flux_param + n_c_param) = conc_values(:, 1); % lb for concentration
%lb(n_flux_param + 1) = conc_values(1, 1); % lb for concentration
%lb(n_flux_param + n_c_param) = 72.5; % c_uridine optimal concentraion in cortex
ub = 100 * ones(n_flux_param+n_c_param,1);
ub(1) = 100; % de novo UMP
ub(2) = 200; % uracil + r5p == uridine
ub(3) = 300; % cytidine == uridine
%ub(4) = 100; % uridine == UMP
%ub(5) = 100; % UMP_out
% ub(n_flux_param - 1) = 2000; % salvage UMP
% ub(n_flux_param) = 3000; % UMP_out
ub(n_flux_param+1:end) = 1.5 * conc_values(:,1);

% add lb and ub for isotopomer fraction parameter
% lb = [lb; zeros(n_isotopomer_frac, 1)];
% ub = [ub; ones(n_isotopomer_frac, 1)];

% add lb and ub for isotopomer fraction parameter
lb  = [lb; 0.1; 0.2; 0];
ub = [ub; 0.4; 0.6; 1];
%lb = [lb; 0; 0; 0];
%ub = [ub; 1; 1; 1];

% Add the bounds for the input metabolite MID
input_mid_array = reshape(MID_input(:,2:end),[],1);
input_mid_array = input_mid_array(~isnan(input_mid_array));
input_sd_array = reshape(SD_input(:,2:end),[],1);
input_sd_array = input_sd_array(~isnan(input_sd_array));
lb_values = input_mid_array-input_sd_array;
lb_values(lb_values < 0) = 0;
ub_values = input_mid_array + input_sd_array;
ub_values(ub_values > 1) = 1; 
lb = [lb;lb_values];
ub = [ub;ub_values];

disp('Parameter bounds defined');

%% calculate the chi-square statistic cut-offs

nparam = numel(lb);
ncons = size(Aeq,1) + numel(conc_metab) + n_input_param + sum(~isnan(MID_balance(:,2:end)),'all'); 
dof = ncons - nparam;
chil = chi2inv(alpha/2,dof);
chiu = chi2inv(1-alpha/2,dof);

%% test dMdt_MID 
% dMdt_test = test_dMdt(time_points, n_flux_param, n_c_param, n_isotopomer_frac, ...
%     MID_balance, balance_data_metabs, SD_balance, MID_input, ...
%     input_data_metabs, SD_input,conc_metab,conc_values);

%% test calc_obj
% tspan = [time_points(1):0.05:time_points(end)];
% niter = 1;
% 
% % Generate initial guess
% x0 = 50 * rand(n_c_param + n_flux_param + n_isotopomer_frac, niter);
% x0 = [x0; repmat(input_mid_array, 1, niter)]; % Use the input MID mean values as initial guess
% 
% obj_test = test_calc_obj(Aeq, beq, tspan, time_points, ...
%     n_flux_param, n_c_param, n_isotopomer_frac, MID_balance, balance_data_metabs, SD_balance, balance_metabs, ...
%     MID_input, input_data_metabs, SD_input, conc_metab,conc_values, lb, ub, ...
%     niter, x0, 'parallel', ncores);
%% Solve the system of ODES
tspan = [time_points(1):0.05:time_points(end)];
niter = 100;

% Generate initial guess
x0 = 100 * rand(n_c_param + n_flux_param + n_isotopomer_frac, niter);
x0(1, :) = 10 * rand(1, niter); % de novo UMP
x0(n_flux_param - 1, :) = 100 * rand(1, niter); % salvage UMP
x0(n_flux_param, :) = 100 * rand(1, niter); % UMP_out

x0 = [x0; repmat(input_mid_array, 1, niter)]; % Use the input MID mean values as initial guess
% concentrations
x0(n_flux_param+1:n_flux_param+n_c_param,:) = (ub(n_flux_param+1:n_flux_param+n_c_param) ...
    - lb(n_flux_param+1:n_flux_param+n_c_param)) .* rand(n_c_param,niter) + lb(n_flux_param+1:n_flux_param+n_c_param);


% aspartate fraction
x0(n_flux_param + n_c_param + 1:n_flux_param + n_c_param + n_isotopomer_frac,:) =(ub(n_flux_param+n_c_param+1:n_flux_param+n_c_param+n_isotopomer_frac) ...
    - lb(n_flux_param+n_c_param+1:n_flux_param+n_c_param+n_isotopomer_frac)) .* rand(n_isotopomer_frac,niter) + lb(n_flux_param+n_c_param+1:n_flux_param+n_c_param+n_isotopomer_frac);
    

[x, fval, exitflag] = estimate_model_param(Aeq, beq, tspan, time_points, ...
    n_flux_param, n_c_param, n_isotopomer_frac, MID_balance, balance_data_metabs, SD_balance, balance_metabs, ...
    MID_input, input_data_metabs, SD_input, conc_metab,conc_values, lb, ub, ...
    niter, x0, 'parallel', ncores);

save([folder, '\data']);

%% Save the workspace data and the optimum solution 

% Select the results that are between the acceptable ranges
up_val = chiu;
if(min(fval) > chiu)
    up_val = min(fval) + 0.1;
end
fval_range = fval(fval < up_val & fval > chil); 
x_range = x(:,fval < up_val & fval > chil);

% Get the minimum objective
[~,minimum] = min(fval_range);
x_opt = x_range(:, minimum);

save([folder,'\data']);
value = [x_opt(1:n_flux_param+n_c_param+n_isotopomer_frac); fval_range(minimum)];
%parameter = [rxns; "c_ump"; "c_uridine"; "f_aspartate_M1"; "f_aspartate_M2"; "f_aspartate_M3"; "min_objective"];
parameter = [rxns; "c_ump"; "f_aspartate_M1"; "f_aspartate_M2"; "f_aspartate_M3"; "min_objective"];
results = table(parameter,value);
writetable(results,[folder,'\results.xlsx'],'Sheet',1);

% Save all the solutaions within acceptable chi square range
writetable(table(parameter,[x_range(1:n_flux_param+n_c_param+n_isotopomer_frac,:);fval_range]), ...
    [folder,'\results.xlsx'],'Sheet',2);

%% Plot de novo UMP MIDs

[~] = plot_MID_denovo(x_opt,tspan,time_points, n_flux_param,n_c_param,n_isotopomer_frac,...
                           MID_balance, balance_data_metabs, SD_balance, ...
                           MID_input,input_data_metabs, SD_input);

%% Plot the MIDs for the solution

[~] = plot_MID(x_opt,tspan,time_points, n_flux_param,n_c_param,n_isotopomer_frac,...
                           MID_balance, balance_data_metabs, SD_balance, ...
                           MID_input,input_data_metabs, SD_input);

%% Estimate the confidence intervals (only the flux and concentration parameters)

thres = chi2inv(cl, 1) + fval_range(minimum);

for k = 1:(n_flux_param+n_c_param+n_isotopomer_frac)
    [cil(k), ciu(k)] = estimate_param_bounds(k,thres,x_opt,Aeq,beq,tspan,time_points,n_flux_param,n_c_param,n_isotopomer_frac,...
                         MID_balance, balance_data_metabs, SD_balance, balance_metabs, ...
                     MID_input,input_data_metabs, SD_input, ...
                     conc_metab,conc_values,lb,ub);
end
%%
intervals = table(parameter(1:end-1),x_opt((1:n_flux_param+n_c_param+n_isotopomer_frac)), cil', ciu',lb(1:n_flux_param+n_c_param+n_isotopomer_frac),ub(1:n_flux_param+n_c_param+n_isotopomer_frac));
intervals.Properties.VariableNames = {'Parameter','Value','lower_CI','upper_CI','lb','ub'};
writetable(intervals,[folder,'\intervals.xlsx'],'Sheet',1);

save([folder,'\data']);





