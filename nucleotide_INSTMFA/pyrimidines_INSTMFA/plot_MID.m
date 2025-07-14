function cout = plot_MID(x, tspan, time_points, n_flux_param, n_c_param, ...
    n_isotopomer_frac, MID_balance, balance_data_metabs, SD_balance, ...
    MID_input, input_data_metabs, SD_input)

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
for m = unique(data_metab_list')

    MID_sim_subset = M_sim(:,ismember(M_metab_list,m));
    MID_subset = MID(ismember(data_metab_list,m),:);
    SD_subset = SD(ismember(data_metab_list,m),:);
    i = [0:1:(size(MID_subset,1)-1)]';

    [~] = plot_metabolite(m,i,MID_sim_subset,MID_subset,SD_subset,time_points,t);
end

end