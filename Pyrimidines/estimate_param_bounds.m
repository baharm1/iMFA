function [cil,ciu] = estimate_param_bounds(i, thres, xopt, Aeq, beq, ...
    tspan, time_points, n_flux_param, n_c_param, n_isotopomer_frac, ...
    MID_balance, balance_data_metabs, SD_balance, balance_metabs, ...
    MID_input, input_data_metabs, SD_input, conc_metab, conc_values, lb, ub)

global folder

% Update Aeq and beq to constrain the parameters
param_cons = zeros(1, size(Aeq, 2));
param_cons(i) = 1;
Aeq = [Aeq; param_cons];
beq = [beq; 0];

% Define the percentage change at each step
change = [10, 5, 2, 1];

% Acceptable exit codes
ok_codes = [0,-100,-101,-102,-103,-400,-401,-402];

% Define the objective function
obj_ode = @(x)calc_obj(x, tspan, time_points, n_flux_param, n_c_param, ...
    n_isotopomer_frac, MID_balance, balance_data_metabs, SD_balance, ...
    balance_metabs, MID_input, input_data_metabs, SD_input, conc_metab, conc_values);
stopobj = thres * 0.98;
kn_opt_E = knitro_options('UseParallel', 'true', 'fstopval', stopobj, 'outlev', 1);

% Estimate upper threshold
flag = 1;
step = 1;
val1 = xopt(i);
while(flag)
    
    disp(['Determining upper bounds for parameter ', num2str(i)]);

    factor = 1 + change(step) / 100;
    val2 = val1 * factor;
    beq(end) = val2;
    x0 = xopt;
    x0(i) = val2;

    if(val2 > ub(i))
        val1 = ub(i);
        clear val2
        disp('REACHED UPPER BOUND')
        break
    elseif(val2 - val1 < 0.1 && i <= n_flux_param + n_c_param)
        clear val2
        disp('NO CHANGE IN UPPER BOUND')
        break
    elseif(val2 - val1 < 0.01 && i > n_flux_param + n_c_param)
        clear val2
        disp('NO CHANGE IN UPPER BOUND')
        break
    end

    [~, objval, exitcode] = knitro_nlp(obj_ode, x0, [], [], Aeq, beq, lb, ub, [], [], kn_opt_E, 'kn_ode.opt');

    if(sum(exitcode == ok_codes) && objval <= thres)
        val1 = val2;
        clear val2
    else
        clear val2
        step = step + 1;
        if(step > numel(change))
            flag = 0;
            disp('HIGH OBJECTIVE FUNCTION FOR UPPER BOUND')
        end
    end

end
ciu = val1;
clear val1 flag step

% Estimate lower threshold
flag = 1;
step = 1;
val1 = xopt(i);
while(flag)
    
    disp(['Determining lower bounds for parameter ', num2str(i)]);

    factor = 1 - change(step) / 100;
    val2 = val1 * factor;
    beq(end) = val2;
    x0 = xopt;
    x0(i) = val2;

    if(val2 < lb(i))
        val1 = lb(i);
        clear val2
        disp('REACHED LOWER BOUND')
        break
    elseif(val1 - val2 < 0.1 && i <= n_flux_param + n_c_param)
        clear val2
        disp('NO CHANGE IN LOWER BOUND')
        break
    elseif(val1 - val2 < 0.01 && i > n_flux_param + n_c_param)
        clear val2
        disp('NO CHANGE IN LOWER BOUND')
        break
    end

    [~, objval, exitcode] = knitro_nlp(obj_ode, x0,[], [], Aeq, beq, lb, ub, [], [], kn_opt_E, 'kn_ode.opt');

    if(sum(exitcode == ok_codes) && objval <= thres)
        val1 = val2;
        clear val2
    else
        clear val2
        step = step + 1;
        if(step > numel(change))
            flag = 0;
            disp('HIGH OBJECTIVE FUNCTION FOR LOWER BOUND')
        end
    end

end
cil = val1;

save([folder, '\parameter_bounds_', num2str(i), '.mat'], "cil", "ciu")

end