%% This file tries different breakpoints for time range 0 to 2 hours
% Initial fluxes and concentration is bound by the CI of the SS solution
% Initial condition terms are added back to the objective function

%% Save the timestamp and create a folder with the name

global timestamp folder
timestamp = ['timestamp_',datestr(datetime,'yyyymmdd_HHMM')];
folder = ['output_files\',char(timestamp)];
mkdir(folder);

%% User-defined inputs

% Name of file with model information
model_file = 'purine_model';

% Name of file with the experimental data
input_data_file = 'Input-Data\Input_brain_rt_conc_abund.xlsx'; % brain  
%input_data_file = 'Input-Data\Input_gbm_rt_conc_abund.xlsx'; % GBM

% Number of initializations
niter = 100;

% Minimum standard deviation
min_sd = sqrt(10^(-5));

% Upper bound for the flux parameters
ub_flux = 2000;

% B-spline properties
% Order n implies a polynomial function of order n-1
bspline_order = 3;
no_int_knots = 2; % Number of internal knots

%% Read input files and pre-process the data

% Read the model information
[S,rxns,input_metabs,balance_metabs] = import_model(model_file);
disp('Model imported');

% Read the experimental MID data
[time_points,input_data_metabs, MID_input, SD_input,...
    balance_data_metabs, MID_balance, SD_balance,...
    conc_ratio,conc_sd,conc_metab_list,init_flux,init_conc,~] = import_expt_data_dmfa(input_data_file,min_sd,input_metabs,balance_metabs);
disp('experimental data imported');


%% Generate stoichiometrically balanced initial flux vectors

Aeq = S;
beq = zeros(size(S,1),1);
lbf = init_flux(:,1);
ubf = init_flux(:,2);
flux0 = repmat(lbf,1,niter) + repmat(ubf-lbf,1,niter).*rand(numel(lbf),niter);

% Stoichiometrically balance the initial flux values
flux_bal = zeros(size(flux0));
for i = 2:niter
    flux_bal(:,i) = fmincon(@(x)sumsqr((x-flux0(:,i))./flux0(:,i)), flux0(:,i), [],[],Aeq,beq,lbf,ubf);
end

% Add the SS solution as one starting point
flux_bal(:,1) = init_flux(:,3);

flux0 = flux_bal;
clear flux_bal Aeq beq

%% Generate bounds for flux parameters

lbp = 0.1*ones(size(S,2),1);
ubp = ub_flux*ones(size(S,2),1);

%% Generate initial guess and bounds for concentration

lbc = init_conc(:,3);
ubc = init_conc(:,4);
conc0 =  repmat(lbc,1,niter) + repmat(ubc-lbc,1,niter).*rand(numel(lbc),niter);
conc0(:,1) = init_conc(:,5);

%% Generate initial guess and bounds for for input metabolite MID

% Input MID array
input_mid_array = reshape(MID_input(:,2:end),[],1);
input_mid_array = input_mid_array(~isnan(input_mid_array));
input_sd_array = reshape(SD_input(:,2:end),[],1);
input_sd_array = input_sd_array(~isnan(input_sd_array));
lbip = input_mid_array-input_sd_array;
lbip(lbip < 0) = 0;
ubip = input_mid_array + input_sd_array;
ubip(ubip > 1) = 1; 

%% Create the metrices for Aeq and beq

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


%% Decide the position of b-spline breakpoints based on concentration data

tspan = [time_points(1):0.02:time_points(end)];
n_c_param = size(conc0,1);

for i = 1:niter

    disp(['Starting iteration ',num2str(i)]);
    [optim_knots(:,i),minR2(i)] = optimize_knot_placement(flux0(:,i), conc0(:,i),tspan,time_points, ...
                                                            MID_input,input_data_metabs, SD_input, ...
                                                            MID_balance,balance_data_metabs,SD_balance,...
                                                            n_c_param, conc_ratio, conc_sd, conc_metab_list, init_conc,...
                                                            S, bspline_order,no_int_knots,balance_metabs,...
                                                            lbc,ubc,lbf,ubf,lbp,ubp,lbip,ubip,input_mid_array,ip_mid_cons,i);

    save([folder,'\data']); 
end

%% Read the save data and visualise the results

% Read the data for the first knot
knotpos1 = [];
fval1 = [];
niter = 100;
for i = 1:niter
    datafilename = [folder,'\iter',num2str(i),'_knot1.mat'];
    load(datafilename);
    knotpos1 = [knotpos1;test_knots'];
    fval1 = [fval1;Rc_values'];
    minR1(i) = min(Rc_values);
    minknot(i) = test_knots(Rc_values == min(Rc_values));
end
writetable(table(minR1',minR2','VariableNames',["One knot", "Two knots"]), ...
    [folder,'\bestfitresults.xlsx'] );

save([folder,'\data']); 

histogram(minknot,'BinWidth',0.025);
title('Best fit knot placement');
xlabel('Knot position'); 
ylabel('Counts');
saveas(gcf, [folder,'\minobj.png']);
close

boxplot(fval1,knotpos1);
title('Objective distribution by knot position');
xlabel('Knot position'); 
ylabel('Objective');
saveas(gcf, [folder,'\allobj.png']);
close

boxplot([minR1';minR2'],[repmat("1 knot",niter,1);repmat("2 knots",niter,1)]);
title('Minimum objective by number of internal knots');
ylabel('Objective');
saveas(gcf, [folder,'\onevstwoknots.png']);
close



