function [obj] = calc_obj(x, tspan, time_points, n_flux_param, n_c_param, ...
    n_isotopomer_frac, MID_balance, balance_data_metabs, SD_balance, ...
    balance_metabs, MID_input, input_data_metabs, SD_input, ...
    conc_metab, conc_values)

% Assign the flux and c values
v = x(1:n_flux_param);
c = x(n_flux_param + 1:n_flux_param + n_c_param);
isotopomer_frac = x(n_flux_param + n_c_param + 1:n_flux_param + ...
    n_c_param + n_isotopomer_frac);

% Assign the values of the input MIDs
ip_mid = x(n_flux_param + n_c_param + n_isotopomer_frac + 1:end);
nan_idx = find(isnan(reshape(MID_input(:, 2:end), [], 1)));
if(~isempty(nan_idx))
    missing = 1;
    % Fill in NaN values to reshape the data correctly
    for k = 1:numel(nan_idx)
        ip_mid = [ip_mid(1:nan_idx(k) - 1); NaN; ip_mid(nan_idx(k):end)];
    end   
else
    missing = 0;
end
% Create the matrix of input MIDs
ip_matrix = reshape(ip_mid, size(MID_input, 1), numel(time_points) - 1);
ip_matrix = [MID_input(:, 1), ip_matrix];

% Calculate the slope of input metabolite MID
slope = input_metab_slope(ip_matrix, time_points, missing);

% Define M0
M0 = [MID_input(:, 1); MID_balance(:, 1)];

% Define the metabolite identifers
M_metab_list = [input_data_metabs; balance_data_metabs];

% Define the function for ODE
odefunc = @(t,M)dMdt_MID(t, M, v, c, isotopomer_frac, time_points, slope, ...
    input_data_metabs, M_metab_list);

reltol = 1e-3; 
flag = 1;
obj = -1;
while(flag)

%     try
        opts = odeset('RelTol',reltol); 
        [t_ode,M_ode] = ode15s(odefunc,tspan,M0,opts); 
        flag = 0;
%         if(t_ode(end) < tspan(end))
%             reltol = 10*reltol;
%             flag = 1;
%         end
%     catch
%         reltol = 10*reltol;
%     end
    
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
M_ip_sim = M_ode(ismember(t_ode, time_points(2:end)), ...
    ismember(M_metab_list, unique(input_data_metabs)) )';
M_ip_expt = MID_input(:, 2:end);
M_ip_sd = SD_input(:, 2:end);
nan_idx = isnan(M_ip_expt);
M_ip_expt(nan_idx) = -1;
M_ip_sd(nan_idx) = -1;
M_ip_error = M_ip_sim - M_ip_expt;
obj1 = sumsqr((~nan_idx) .* (M_ip_error ./ M_ip_sd));

% calculate the error for balance metab mid
clear nan_idx
M_b_sim = M_ode(ismember(t_ode, time_points(2:end)), ismember(M_metab_list, ...
    unique(balance_data_metabs)))';
M_b_expt = MID_balance(:, 2:end);
M_b_sd = SD_balance(:, 2:end);
nan_idx = isnan(M_b_expt);
M_b_expt(nan_idx) = -1;
M_b_sd(nan_idx) = -1;
M_b_error = M_b_sim - M_b_expt;
obj2 = sumsqr((~nan_idx) .* (M_b_error ./ M_b_sd));

% Calculate the error for concentration
obj3 = 0;
balance_metabs = ["URIDINE","UMP"];
for k = 1:length(conc_metab)
    c_metab = conc_metab(k);
    idx = find(balance_metabs == c_metab{1});
    expt_c = conc_values(k,1);
    expt_sd = conc_values(k,2);
    obj3 = obj3 + ((c(idx)-expt_c)/expt_sd)^2;
end

obj = obj1+obj2+obj3;

