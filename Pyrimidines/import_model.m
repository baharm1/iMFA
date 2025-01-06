function [S,rxns,input_metabs,balance_metabs] = import_model(model_file)

[~,rxns] = xlsread(model_file,'reactions','','basic');
[~,input_metabs] = xlsread(model_file,'input_metabs','','basic');
[~,balance_metabs] = xlsread(model_file,'balance_metabs','','basic');

%% Create stoichiometric matrix

% convert the equations to symbolic form
n_rxn = size(rxns,1);
for i = 1:n_rxn
    text2(i) = str2sym(rxns(i));
end
metab_sym = symvar(text2);

% Create stoichiometric matrix
S = equationsToMatrix(text2);
S = double(-S');
% only keep the metabs to be balanced
index = ismember(string(metab_sym),string(balance_metabs));
S = S( index ,:);

%%
% Save the list of balanced metabolites in sequence with S
balance_metabs = string(metab_sym(index));
input_metabs = string(input_metabs)';

end