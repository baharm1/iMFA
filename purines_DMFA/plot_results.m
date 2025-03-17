function cout = plot_results(x, tspan, time_points, MID_input, ...
    input_data_metabs, SD_input, MID_balance, balance_data_metabs, ...
    SD_balance, n_c_param, c_param_index, n_param, conc_ratio, conc_SD, ...
    conc_metab_list, S, knot_pos, bspline_order, balance_metabs, ...
    flux_param_index, data_i_list, rxns)

cout = 0;

% Create the directory to save the plots
global folder fig_count
mkdir([folder, '\plots', num2str(fig_count)]);


% Define CP
CP = reshape(x(flux_param_index), size(S, 2), []);

%% Estimate the slope of input metabolites

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
ip_matrix = reshape(ip_mid,size(MID_input, 1), numel(time_points) - 1);
ip_matrix = [MID_input(:, 1), ip_matrix];

% Calculate the slope of input metabolite MID
slope = input_metab_slope(ip_matrix, time_points, missing);
slope_metabs = input_data_metabs;

%% Define the initial vector for ODE

%  % Define M0 (for MID)
MID = [MID_input; MID_balance];
SD = [SD_input; SD_balance];
M0 = MID(:, 1);
M0(length(M0) + 1:length(M0) + 6) = [1; zeros(5, 1)];  % Add the initial MID for GMP

% Define the isptopologue identifers
data_metab_list = [input_data_metabs; balance_data_metabs];
M_metab_list = [data_metab_list; repmat('GMP', 6, 1)];

% % Define C0 (for concentration)
C0 = x(c_param_index);

% Combine the two
Y0 = [C0; M0];

%% Solve ODE

% Define the function for ODE
odefunc = @(t, y)dYdt_MID_mass(t, y, CP, bspline_order, knot_pos, S, ...
    n_c_param, time_points, slope, slope_metabs, M_metab_list);

% Define the mass matrix
mass = @(t, y)mass_matrix(t, y, n_c_param, M_metab_list);

reltol = 0.001;
abstol = 0.001;
odeOptions = odeset('RelTol', reltol, 'AbsTol', abstol, 'Mass', mass);

[t, Y] = ode23tb(odefunc, tspan, Y0, odeOptions); 

% Plot the simulated MIDs
M_sim = Y(:, n_c_param + 1:end);
for m = unique(data_metab_list')

    i = data_i_list(ismember(data_metab_list, m));
    MID_sim_subset = M_sim(:, ismember(M_metab_list, m));
    MID_subset = MID(ismember(data_metab_list, m),:);
    SD_subset = SD(ismember(data_metab_list, m), :);
    [~] = plot_metabolite(m, i, MID_sim_subset, MID_subset, SD_subset, ...
        time_points, t);
end

% Plot concentration profiles
C = Y(:, 1:n_c_param);
C0_matrix = diag(1./C(1, :));
Cratio = C * C0_matrix;
Cratio = Cratio(:, ismember(balance_metabs, conc_metab_list));
for m = 1:length(conc_metab_list)
    plot(t, Cratio(:, m));
    hold on;
    scatter(time_points(2:end), conc_ratio(m, 2:end)');
    hold on;
    errorbar(time_points(2:end), conc_ratio(m, 2:end)', conc_SD(m, 2:end)', ...
        'LineStyle', 'none');
    title([cell2mat(conc_metab_list(m)), ' conc.']);
     xlabel('time (h)');
     ylabel('ratio');
     filename = [folder, '\plots\', char(conc_metab_list(m)), '_conc.png'];
     saveas(gcf, filename);
    close
end

%% Plot the flux profiles
for e = 1:length(tspan)
    time = tspan(e);
    Nout = bSplineMat_lite(knot_pos, time / tspan(end), bspline_order);
    v(:, e) = CP * Nout;
end

for i = 1:numel(rxns)
    
    plot(tspan, v(i, :));
    title(rxns(i));
    xlabel('time (h)');
    ylabel('Flux (pmol/hr.mgtissue)');
    filename = [folder, '\plots\', 'flux', num2str(i), '.png'];
    saveas(gcf, filename);
    close
      
end
end