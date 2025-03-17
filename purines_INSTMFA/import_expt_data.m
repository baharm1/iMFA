function [time_points, input_data_metabs, MID_input, SD_input,...
    balance_data_metabs, MID_balance, SD_balance,...
    conc_metab,conc_values] = import_expt_data(input_data_file,min_sd,input_metabs,balance_metabs)

[MID,text] = xlsread(input_data_file,'MID','','basic');
[SD,~] = xlsread(input_data_file,'SD','','basic');
[~,unlabeled_metabs] = xlsread(input_data_file,'unlabeled_metabs','','basic');

unlabeled_metabs = string(unlabeled_metabs)';

% re-format the data
time_points = MID(1,2:end);
data_metab_list = string(text(2:end,1));
MID = MID(2:end, 2:end);
MID = MID/100;
SD = SD(2:end, 2:end);
SD = SD/100;

% Set a minimum SD value
SD(SD < min_sd) = min_sd;

% Read the conc. value 
[conc_values,conc_metab] = xlsread(input_data_file,'conc','','basic');

% Separate the data for input metabs
labeled_input_metabs = input_metabs(~ismember(input_metabs,unlabeled_metabs));
input_metab_flag = ismember(data_metab_list,labeled_input_metabs);
input_data_metabs = data_metab_list(input_metab_flag);
MID_input = MID(input_metab_flag,:);
SD_input = SD(input_metab_flag,:);

% Separate the data for balance metabs
balance_metab_flag = ismember(data_metab_list,balance_metabs);
balance_data_metabs = data_metab_list(balance_metab_flag);
MID_balance = MID(balance_metab_flag,:);
SD_balance = SD(balance_metab_flag,:);

end