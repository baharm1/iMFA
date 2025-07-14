% Code to estimate the fluxes in the ser-gly pathway
% at metabolic steady state
clc; clear;
%% Input files
files = dir("mice_input_data");
files = struct2table(files);
mouse_models = files.name(~files.isdir);

% Read model information
file1 = 'model_2comp_mice.xlsx';
[Sfull,rxn_full,irrev_rxn,exc_rxn,rev_rxn,metab_char,...
    AM_full, metab_size,input_metabs,metabs_remove] = import_stoich_AM(file1);
disp('Stoichiometry and atom-mapping imported from files');

%% User defined inputs

% Unlabeled metabolites
unlabeled_metabs = {'CO2','SERu'};

% Define minimum standard deviation
min_sd = sqrt(10^(-5));

% Upper bounds for irreversible fluxess
UBval = 20;

% Scaling factor for reversible fluxes
P = 1;

% Number of initial guess for optimization
niter = 200;

% Number of Monte-Carlo simulations for parameter sensitivity
nmc = 1000;

% Confidence level for flux bounds
conf = 0.95;

%% Define output folder

global timestamp_gl filepath_stored_gl
delete(gcp('nocreate'));
parpool(10);

for ps = 1:length(mouse_models)

timestamp_gl = ['timestamp_', datestr(datetime,'yyyymmdd_HHMM'), char(mouse_models(ps))];
filepath_stored_gl = ['output_files_mice/', timestamp_gl];
mkdir(filepath_stored_gl)

file2 = ['mice_input_data\', char(mouse_models(ps))];

% Read experimental data
[MID_metab, MID, SD] = ...
    import_experiment_data(file2,min_sd);
disp('Experimental data imported');


%% Generate inputs for optimization

nRev = length(rev_rxn)*2;
nflux = length(rxn_full);

[varFlag,mx,Aeq,beq,lb,ub] = fmincon_inputs(Sfull,nRev,metab_size,metab_char,nflux,...
    input_metabs,metabs_remove,UBval,unlabeled_metabs);

% Generate non-linear constraint equations
[num_nlcons] = IMMeqn_sym_gen(Sfull, metab_char, AM_full, metab_size, nflux, varFlag, mx,input_metabs,metabs_remove);

% Function to calculate objective value
get_objective = @(x)calc_obj(x,MID_metab,MID,SD,metab_char,metab_size,nflux,varFlag);

% Function to calculate non-linear constraint equations
get_NL_cons = @(x)NL_const_wrapper(x,nflux,P,nRev,file1);

if(numel(lb) ~= numel(ub))
    error('ERROR: The size of the lower bounds vector and upper bounds vector is not equal.');
end

if(numel(lb) ~= size(Aeq,2))
    error('ERROR: Size of the parameter vector not compatible with Aeq');
end

if(numel(beq) ~= size(Aeq,1))
    error('ERROR: Size of the beq vector not compatible with Aeq');
end

%% Estimate the degress of freedom and optimal objective value

nknowns = numel(MID);
ncons = size(Aeq,1) + num_nlcons;
dof = nknowns + ncons - numel(lb);

% Calculate acceptable objective values for chi-square goodness of fit test
chiu = chi2inv(0.975,dof);
chil = chi2inv(0.025,dof);

%% % User-defined linear constraints based on the model

% Atleast 20% of glycine produced is being consumed
A = zeros(2,numel(lb));
A(1,4) = 0.2;
A(1,[5,14]) = -1;
A(2,10) = 0.2;
A(2,[11,17]) = -1;
b = zeros(2,1);

% The lower-bound for glycine producing reactions is 0.1 (to eliminate
% trivial solutions)
lb(4) = 0.1;
lb(10) = 0.1;

%% Perform optimization

% Generate random guess
x0 = rand(numel(lb),niter);

% Define optimization options
opt_options = knitro_options('algorithm',4,'honorbnds',2,'bar_relaxcons',0,'outlev',2,...
                        'gradopt', 2, 'maxtime_real',1800,'opttol',1e-03,'opttol_abs',1e-02,...
                        'UseParallel','true','bar_maxcrossit',20);

% Perform optimization
fval = zeros(niter,1);
exitE = zeros(niter,1);
xopt = zeros(size(x0));
for i = 1:niter
    
    disp(['Starting iteration ',num2str(i),' out of ',num2str(niter)]);
    [xopt(:,i),fval(i),exitE(i)] = knitro_nlp(get_objective,x0(:,i),A,b,Aeq,beq,lb,ub,get_NL_cons,[],opt_options,[]);
    save([filepath_stored_gl,'/data.mat']);

end

% Check the exit flag and only keep the results with acceptable exit flags
accept_flags = [0,-100,-101,-102,-103,-400,-401,-402];
xopt_feas = xopt(:,ismember(exitE,accept_flags));
fval_feas = fval(ismember(exitE,accept_flags));
exit_feas = exitE(ismember(exitE,accept_flags));


%% Save the best result

% Select result with minimum objective value
fval_min = min(fval_feas);
xopt_min = xopt_feas(:,fval_feas == fval_min);

% Write flux results
flux_table = table([rxn_full;'minimum objective'],[xopt_min(1:nflux);fval_min]);
flux_table.Properties.VariableNames = ["Reaction","Flux"];
writetable(flux_table,[filepath_stored_gl,'/flux_results.xlsx'],"Sheet","flux");

% Write IDV results
idv_list = [];
for i = 1:numel(metab_char)
    if (varFlag(i)==0)
        l = metab_size(i);
        idv_list = [idv_list; cellstr([repmat([char(metab_char(i)),'_'],2^l,1),dec2bin(0:(2^l)-1,l)])];
    end
end
idv_table = table(idv_list,xopt_min(nflux+1:end));
idv_table.Properties.VariableNames = ["IDV","Value"];
writetable(idv_table,[filepath_stored_gl,'/flux_results.xlsx'],"Sheet","IDV");

save([filepath_stored_gl,'/data.mat']);


%% Perform parameter sensitivity analysis by Monte Carlo Sampling

% Generate sampled data
mid_mc = mcb_sampler(nmc,MID_metab,MID,SD);


% Optimize for sampled data
for i = 1:nmc
    clear mc_objective
    mc_objective = @(x)calc_obj(x,MID_metab,mid_mc(:,i),SD,metab_char,metab_size,nflux,varFlag);
    [xmc(:,i),fmc(i),exitmc(i)] = knitro_nlp(mc_objective,xopt_min,A,b,Aeq,beq,lb,ub,get_NL_cons,[],opt_options,[]);
    
end
save([filepath_stored_gl,'/data.mat']);
% Determine and save the flux confidence bounds
alpha = 100*(1-conf)/2;
ubflux = zeros(nflux,1);
lbflux = ubflux;
for i = 1:nflux
    flux_values = xmc(i,:);
    ubflux(i) = prctile(flux_values,100-alpha);
    lbflux(i) = prctile(flux_values,alpha);
end
bounds_table = table(rxn_full,xopt_min(1:nflux),lbflux,ubflux);
bounds_table.Properties.VariableNames = ["Reaction","Flux","Lower_CI","Upper_CI"];
writetable(bounds_table,[filepath_stored_gl,'/flux_results.xlsx'],"Sheet","bounds");

%% Plot the result

[~] = plot_MID(xopt_min,metab_size,nflux,varFlag,metab_char,unlabeled_metabs,MID_metab, MID, SD);

%% Estimate the bounds for the flux ratios
ratiog = xmc(ismember(rxn_full,'PGg == SERg'),:)./xmc(ismember(rxn_full,'SERp == SERg'),:);
ratiob = xmc(ismember(rxn_full,'PGb == SERb'),:)./xmc(ismember(rxn_full,'SERp == SERb'),:);

lower_bounds = [prctile(ratiog,alpha);prctile(ratiob,alpha)];
upper_bounds = [prctile(ratiog,100-alpha);prctile(ratiob,100-alpha)];

fluxes = xopt_min(1:nflux);
valn = fluxes(ismember(rxn_full,'PGg == SERg'),:)./fluxes(ismember(rxn_full,'SERp == SERg'),:);
valb = fluxes(ismember(rxn_full,'PGb == SERb'),:)./fluxes(ismember(rxn_full,'SERp == SERb'),:);

table_new = table(["GBM";"Brain"],[valn;valb],lower_bounds,upper_bounds);
table_new.Properties.VariableNames = ["Tissue","Ratio de_novo_to_uptake","Lower_CI","Upper_CI"];
writetable(table_new,[filepath_stored_gl,'/flux_ratio.xlsx']);
save([filepath_stored_gl,'/data.mat']);
end