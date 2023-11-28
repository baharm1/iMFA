function cout = plot_MID_denovo(x,tspan,time_points, n_flux_param,n_c_param,n_isotopomer_frac,...
                           MID_balance, balance_data_metabs, SD_balance, ...
                           MID_input,input_data_metabs, SD_input)

cout = 0;

% Create the directory to save the plots
global folder
mkdir([folder,'\plots'])

%% Assign parameter values

% Assign the flux and c values
v = x(1:n_flux_param);
c = x(n_flux_param+1:n_flux_param+n_c_param);
isotopomer_frac = x(n_flux_param + n_c_param + 1:n_flux_param + n_c_param + n_isotopomer_frac);
% Assign the values of the input MIDs
ip_mid = x(n_flux_param+n_c_param+n_isotopomer_frac+1:end);
nan_idx = find(isnan(reshape(MID_input(:,2:end),[],1)));
if(~isempty(nan_idx))
    missing = 1;
    % Fill in NaN values to reshape the data correctly
    for k = 1:numel(nan_idx)
        ip_mid = [ip_mid(1:nan_idx(k)-1);NaN;ip_mid(nan_idx(k):end)];
    end   
else
    missing = 0;
end
% Create the matrix of input MIDs
ip_matrix = reshape(ip_mid,size(MID_input,1),numel(time_points)-1);
ip_matrix = [MID_input(:,1),ip_matrix];

% Calculate the slope of input metabolite MID
slope = input_metab_slope(ip_matrix, time_points,missing);

%% Integrate the ODEs

% Define M0
M0 = [MID_input(:,1);MID_balance(:,1)];

% Define the metaboliteidentifers
M_metab_list = [input_data_metabs;balance_data_metabs];

% Define the function for ODE
odefunc = @(t,M)dMdt_MID(t,M,v,c,isotopomer_frac,time_points,slope, input_data_metabs, M_metab_list);

reltol = 1e-3; 
flag = 1;
while(flag)

    try
        opts = odeset('RelTol',reltol); 
        [t,M_sim] = ode15s(odefunc,tspan,M0,opts); 
        
        flag = 0;
        if(t(end) < tspan(end))
            reltol = 10*reltol;
            flag = 1;
        end
    catch
        reltol = 10*reltol;
    end
    
    if(reltol > 10)
        disp('relative tolerance is too big')
    end

end

if(reltol > 10)
    return
end


% Plot the simulated MIDs
data_metab_list = [input_data_metabs;balance_data_metabs];
MID = [MID_input;MID_balance];
SD = [SD_input;SD_balance];

f_aspartate_M1 = isotopomer_frac(1);
f_aspartate_M2 = isotopomer_frac(2);
f_aspartate_M3 = isotopomer_frac(3);

for j = 1:length(t)
    M_r5p = M_sim(j, ismember(M_metab_list,'R5P'));
    M_aspartate = M_sim(j, ismember(M_metab_list,'ASPARTATE'));

    N_aspartate = zeros(4, 1);
    N_aspartate(1) = M_aspartate(1) + f_aspartate_M1 * M_aspartate(2);
    N_aspartate(2) = (1 - f_aspartate_M1) * M_aspartate(2) + f_aspartate_M2 * M_aspartate(3);
    N_aspartate(3) = (1 - f_aspartate_M2) * M_aspartate(3) + f_aspartate_M3 * M_aspartate(4);
    N_aspartate(4) = (1 - f_aspartate_M3) * M_aspartate(4) + M_aspartate(5);
    M_aspartate_r5p(j, :) = sum_diag(N_aspartate * M_r5p)';
end

i = [0:1:5]';

MID_subset = MID(ismember(data_metab_list,'UMP'),:);
SD_subset = SD(ismember(data_metab_list,'UMP'),:);

[~] = plot_metabolite('UMP_{denovo}', i, M_aspartate_r5p, MID_subset, SD_subset, time_points, t);

end

function sum = sum_diag(N)

[m,n] = size(N);
sum = zeros(m+n-1,1);
for i = 1:m
    for j = 1:n
        sum(i+j-1) = sum(i+j-1)+ N(i,j);
    end
end

end