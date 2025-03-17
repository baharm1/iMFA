function [obj1,obj2,obj3] = calc_obj_timed(x, tspan, time_points, S, ...
    knot_pos, bspline_order, n_c_param, c_param_index, n_param,...
    MID_input, input_data_metabs, SD_input, MID_balance, ...
    balance_data_metabs, SD_balance, ...
    conc_ratio, conc_sd, conc_metab_list, init_conc, balance_metabs, ...
    flux_param_index, stoptime)

%Define CP
CP = reshape(x(flux_param_index), size(S, 2), []);

%% Calculate the slope of input metabolites

% Assign the values of the input MIDs
ip_mid = x(n_param + 1:end);
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
ip_matrix = [MID_input(:,1),ip_matrix];

% Calculate the slope of input metabolite MID
slope = input_metab_slope(ip_matrix, time_points,missing);
slope_metabs = input_data_metabs;

%% Combine the MID data for input metabs and balance metabs

data_metab_list = [input_data_metabs; balance_data_metabs];
MID = [MID_input; MID_balance];
SD = [SD_input; SD_balance];

%% Define the initial vector for ODE

%  % Define M0 (for MID)
M0 = MID(:,1);
M0(length(M0) + 1:length(M0) + 6) = [1; zeros(5, 1)];  % Add the initial MID for GMP

% Define the isptopologue identifers
M_metab_list = [data_metab_list; repmat('GMP', 6, 1)];

% % Define C0 (for concentration)
C0 = x(c_param_index);

% Combine the two
Y0 = [C0; M0];

%% Solve ODE

% Define the function for ODE
odefunc = @(t,y)dYdt_MID_mass(t, y, CP, bspline_order, knot_pos, S, ...
    n_c_param, time_points, slope, slope_metabs, M_metab_list);

% Define the mass matrix
mass = @(t, y)mass_matrix(t, y, n_c_param, M_metab_list);

% Set the tolerances
reltol = 0.001;
abstol = 0.001;

% Define the timeout function
timer = @(t, y)EventFcn(t, y, stoptime);

flag = 1;
while(flag)
    
    if(reltol > 100)
        break;
    end
    if(datetime('now') >= stoptime)
        break;
    end
    try
        odeOptions = odeset('RelTol',reltol,'AbsTol',abstol, 'Mass', mass,'Events',timer);
        [t, Y, ~, ~, ie] = ode23tb(odefunc, tspan, Y0, odeOptions);
    catch
        if(ie)
            break;
        end
        reltol = 10*reltol;
        abstol = 10*abstol;
        continue;
    end
   
    if(ie)
            break;
    elseif(t(end) == tspan(end))
        flag = 0;
    else
        reltol = 10*reltol;
        abstol = 10*abstol;
    end

end

if(flag)
    obj1 = 0;
    obj2 = 0;
    obj3 = 0;
    return
end


%% Calculate the value of the objective function for concentration

C_sim = Y(ismember(t,time_points(2:end)) , 1:n_c_param)';
c0_matrix = diag(1./C0);
C_ratio_sim = c0_matrix*C_sim;
C_ratio_sim = C_ratio_sim(ismember(balance_metabs, conc_metab_list),:);
C_error = C_ratio_sim - conc_ratio(:,2:end);
obj1 = sumsqr(C_error./conc_sd(:,2:end));

%% objective function for MIDs

M = Y(:,n_c_param+1:end);
M_sim = M(ismember(t,time_points(2:end)), ismember(M_metab_list,data_metab_list))';
M_sd = SD(:,2:end);
M_expt = MID(:,2:end);
M_error = M_sim - M_expt;
obj2 = sumsqr( (1-isnan(M_error)).*M_error./M_sd);

%% Objective for initial concentration values

c_error = C0 - init_conc(:,1);
nan_idx = isnan(c_error);
obj3 = sumsqr((1-nan_idx).*c_error./init_conc(:,2));


end

function [value,isterminal,direction] = EventFcn(t,y,stoptime)
    
    isterminal = 1;
    direction = 0;
    if(datetime('now') > stoptime)
        value = 0;
    else
        value = 1;
    end

end