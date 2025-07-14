
%% User defined inputs

% Input filess
file1 = 'model_3comp_human.xlsx';
file2 = 'patient_input_data/Patient1.xlsx';

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


%% Read the model and data files

% Read model information
[Sfull,rxn_full,irrev_rxn,exc_rxn,rev_rxn,metab_char,...
    AM_full, metab_size,input_metabs,metabs_remove] = import_stoich_AM(file1);
disp('Stoichiometry and atom-mapping imported from files');

% Read experimental data
[~, MID, ~] = ...
    import_experiment_data(file2,min_sd);
disp('Experimental data imported');


%% Generate inputs for optimization

nRev = length(rev_rxn)*2;
nflux = length(rxn_full);

[varFlag,mx,Aeq,beq,lb,ub] = fmincon_inputs(Sfull,nRev,metab_size,metab_char,nflux,...
    input_metabs,metabs_remove,UBval,unlabeled_metabs);

% Generate non-linear constraint equations
[num_nlcons] = IMMeqn_sym_gen(Sfull, metab_char, AM_full, metab_size, nflux, varFlag, mx,input_metabs,metabs_remove);


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
A = zeros(3,numel(lb));
A(1,4) = 0.2;
A(1,[5,20]) = -1;
A(2,10) = 0.2;
A(2,[11,23]) = -1;
A(3,16) = 0.2;
A(2,[17,26]) = -1;
b = zeros(3,1);

% The lower-bound for glycine produing reactions is 0.1 (to eliminate
% trivial solutions)
lb(4) = 0.1;
lb(10) = 0.1;
lb(16) = 0.1;

%% Save the model variables
clear MID
save("model.mat")
