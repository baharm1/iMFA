function [varFlag,mx,Aeq,beq,lb,ub] = fmincon_inputs(Sfull,nRev,metab_size,metab_char,nflux,...
    input_metabs,metabs_remove,UBval,unlabeled_metabs)


% Generate vector of known IDVs for unlabeled metabolites
unlabeled_idx = find(ismember(metab_char,unlabeled_metabs));
varFlag = zeros(1,length(metab_size));
varFlag(unlabeled_idx) = 1;

% Assign IDVs for unlabeled metabolites
kIDVsize = sum((varFlag).*(2.^metab_size));
mx = zeros(kIDVsize,1);
for i = unlabeled_idx
    index = kIDVindex(i,metab_size,varFlag);
    mx(index(1),1) = 1;
end

clear kIDVsize unlabeled_idx


%% Create linear constraints matrices

% Stoichimetric terms
tIDVsize = sum(((2.^metab_size).*(1-varFlag)));
% Remove non-mass-balanced metabolites from the stoichiometric matrix
remove_idx = ismember(metab_char,[input_metabs;metabs_remove]');
Stemp = Sfull(~remove_idx,:);
% Extracting stoichiometry of only net reactions from the S matrix
% Stemp = Sfull;
Stemp(:,2:2:nRev) = zeros(size(Sfull,1),nRev/2);
flag_rows = size(Stemp,1);
Aeq=[Stemp,zeros(flag_rows,tIDVsize)];
[r,s]=size(Aeq);
beq = zeros(r,1);

% % Set the net serine input flux to 1 for 3comp model
% Aeq = [Aeq;zeros(3,size(Aeq,2))];
% Aeq(numel(beq)+1,[1:3]) = 1;
% Aeq(numel(beq)+2,[7:9]) = 1;
% Aeq(numel(beq)+3,[13:15]) = 1;
% beq = [beq;ones(3,1)];

% Set the net serine input flux to 1 for 2comp model
Aeq = [Aeq;zeros(2,size(Aeq,2))];
Aeq(numel(beq)+1,[1:3]) = 1;
Aeq(numel(beq)+2,[7:9]) = 1;
beq = [beq;ones(2,1)];

% Sum of input IDVs is 1
labeled_ip = setdiff(input_metabs,unlabeled_metabs);
for i = 1:numel(labeled_ip)
    
    m = find(ismember(metab_char,labeled_ip(i)));
    index = uIDVindex(m,metab_size,varFlag,nflux);
    Aeq = [Aeq; zeros(1,size(Aeq,2))];
    Aeq(end,index(1):index(2)) = 1;
    beq = [beq;1];
end

%% Lower Bounds

% Lower bound for irreversible reaction, irr = 0.00001
% Lower bound for reversible reaction, rev = -UBval
% Lower bound for IDVs = 0
lb = (1e-7)*ones(s,1);
lb(1:2:nRev-1,1) = -UBval;    
       
%% Upper Bounds

% Upper bound for all reactions = +UBval
% Upper bound for IDVs = 1
ub = ones(s,1);
ub(1:nflux,1) = UBval;
ub(2:2:nRev,1) = 0.999999;


