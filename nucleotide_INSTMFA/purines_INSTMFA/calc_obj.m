function [obj] = calc_obj(x,tspan,time_points, n_flux_param,n_c_param,...
                           MID_balance, balance_data_metabs, SD_balance, ...
                           MID_input,input_data_metabs, SD_input, ...
                                conc_metab,conc_values)

% Assign the flux and c values
v = x(1:n_flux_param);
c = x(n_flux_param+1:n_flux_param+n_c_param);

% Assign the values of the input MIDs
ip_mid = x(n_flux_param+n_c_param+1:end);
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

% Define M0
M0 = [MID_input(:,1);MID_balance(:,1)];
init_gmp = zeros(6,1);
init_gmp(1) = 1;
M0(length(M0)+1:length(M0)+6) = init_gmp;  % Add the initial MID for GMP

% Define the metaboliteidentifers
M_metab_list = [input_data_metabs;balance_data_metabs];
M_metab_list = [M_metab_list; repmat('GMP',6,1)]; % Add GMP to the list

% Define the function for ODE
odefunc = @(t,M)dMdt_MID(t,M,v,c, time_points,slope, input_data_metabs, M_metab_list);

reltol = 1e-3; 
flag = 1;
obj = -1;
while(flag)

    try
        opts = odeset('RelTol',reltol); 
        [t,M] = ode15s(odefunc,tspan,M0,opts); 
        flag = 0;
    catch
        reltol = 10*reltol;
    end

    if(t(end) < tspan(end))
        reltol = 10*reltol;
        flag = 1;
    end
    
    if(reltol > 10)
        obj = 0;
        flag = 0;
    end

end

if(obj == 0)
    return
end

% calculate the value of the objective function

% Calculate the error for input mid data
M_ip_sim = M(ismember(t,time_points(2:end)), ismember(M_metab_list,unique(input_data_metabs)) )';
M_ip_expt = MID_input(:,2:end);
M_ip_sd = SD_input(:,2:end);
nan_idx = isnan(M_ip_expt);
M_ip_expt(nan_idx) = -1;
M_ip_sd(nan_idx) = -1;
M_ip_error = M_ip_sim - M_ip_expt;
obj1 = sumsqr((~nan_idx).*(M_ip_error./M_ip_sd));

% calculate the error for balance metab mid
clear nan_idx
M_b_sim = M(ismember(t,time_points(2:end)), ismember(M_metab_list,unique(balance_data_metabs)) )';
M_b_expt = MID_balance(:,2:end);
M_b_sd = SD_balance(:,2:end);
nan_idx = isnan(M_b_expt);
M_b_expt(nan_idx) = -1;
M_b_sd(nan_idx) = -1;
M_b_error = M_b_sim - M_b_expt;
% if(isempty(mid_to_remove))
obj2 = sumsqr((~nan_idx).*(M_b_error./M_b_sd));
% else
%     to_keep = ones(size(M_b_error));
%     for i = 1:numel(unique(balance_data_metabs))
%         for j = 1:numel(unique(mid_to_remove))
%             to_keep(6*(i-1)+mid_to_remove(j)+1,:) = 0;
%         end
%     end
%     obj2 = sumsqr(to_keep.*(~nan_idx).*(M_b_error./M_b_sd));
% end

% Calculate the error for concentration
obj3 = 0;
balance_metabs = ["AMP","GDP","GMP","GUANOSINE","IMP","INOSINE"];
for k = 1:length(conc_metab)
    idx = find(balance_metabs == conc_metab(k));
    expt_c = conc_values(k,1);
    expt_sd = conc_values(k,2);
    obj3 = obj3 + ((c(idx)-expt_c)/expt_sd)^2;
end

obj = obj1+obj2+obj3;

