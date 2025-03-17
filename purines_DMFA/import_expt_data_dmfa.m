function [time_points,input_data_metabs, MID_input, SD_input,...
    balance_data_metabs, MID_balance, SD_balance,...
    conc_ratio,conc_sd,conc_metab_list,init_flux,init_conc,data_i_list] = import_expt_data_dmfa(input_data_file,min_sd,input_metabs,balance_metabs)

[MID,text] = xlsread(input_data_file,'MID','','basic');
[SD,~] = xlsread(input_data_file,'SD','','basic');
[~,unlabeled_metabs] = xlsread(input_data_file,'unlabeled_metabs','','basic');
[conc_ratio,conc_metab_list] = xlsread(input_data_file,'conc','','basic');
[conc_sd,~] = xlsread(input_data_file,'conc SD','','basic');

unlabeled_metabs = string(unlabeled_metabs)';

% reformat the MID data
time_points = MID(1,2:end);
data_metab_list = string(text(2:end,1));
MID_i = MID(2:end,1);
MID = MID(2:end, 2:end);
MID = MID/100;
SD = SD(2:end, 2:end);
SD = SD/100;

% Set a minimum SD value
SD(SD < min_sd) = min_sd;

% Separate the data for input metabs
labeled_input_metabs = input_metabs(~ismember(input_metabs,unlabeled_metabs));
input_metab_flag = ismember(data_metab_list,labeled_input_metabs);
input_data_metabs = data_metab_list(input_metab_flag);
MID_input = MID(input_metab_flag,:);
SD_input = SD(input_metab_flag,:);
input_i = MID_i(input_metab_flag);

% Separate the data for balance metabs
balance_metab_flag = ismember(data_metab_list,balance_metabs);
balance_data_metabs = data_metab_list(balance_metab_flag);
MID_balance = MID(balance_metab_flag,:);
SD_balance = SD(balance_metab_flag,:);
balance_i = MID_i(balance_metab_flag);

% Consolidate the M+i list
data_i_list = [input_i;balance_i];

% Re-order the concentration data to be as the same order as stoic. matrix
conc_metab_list = string(conc_metab_list(2:end))';
conc_ratio = conc_ratio(2:end,:);
conc_sd = conc_sd(2:end,:);

order = get_order(conc_metab_list, balance_metabs(ismember(balance_metabs,conc_metab_list)));
conc_ratio = conc_ratio(order,:);
conc_sd = conc_sd(order,:);
conc_metab_list = conc_metab_list(order);

% Read the results from SS run to use as initial concentration and flux
% values
[init_flux,~] = xlsread(input_data_file,'v0','','basic');
[init_conc,~] = xlsread(input_data_file,'c0','','basic'); % The metabolites should be in a specific sequence

end

function order = get_order(new_list, old_list)

order = zeros(1,length(old_list));
for i = 1:length(order)
    index = find(new_list == old_list(i));
    order(i) = index;
end


end
