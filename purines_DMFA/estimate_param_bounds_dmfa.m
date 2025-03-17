function [cil,ciu] = estimate_param_bounds_dmfa(i,thres,xopt,Aeq,beq,tspan,time_points,n_c_param,...
                         MID_balance, balance_data_metabs, SD_balance, ...
                        MID_input,input_data_metabs, SD_input,lb,ub,...
                        c_param_index,n_param, conc_ratio, conc_sd, conc_metab_list, init_conc,...
                        S, knot_pos, bspline_order,balance_metabs, flux_param_index,l2)

global folder

% Update Aeq and beq to constrain the parameters
param_cons = zeros(1,size(Aeq,2));
param_cons(i) = 1;
Aeq = [Aeq;param_cons];
beq = [beq;0];

% Define the percentage change at each step
change = [10, 5, 2, 1];

% Acceptable exit codes
ok_codes = [0,-100,-101,-102,-103,-400,-401,-402];

% Define the non-linear contraints function
NLcon = @(x)NL_cons(x, S, c_param_index, flux_param_index, bspline_order, knot_pos, tspan);

% knitro options
stopobj = thres*0.99;
kn_opt_E = knitro_options('UseParallel','true','fstopval',stopobj);

% Time limit
time_limit = minutes(20);

% Estimate upper confidence interval
flag = 1;
step = 1;
val1 = xopt(i);
while(flag)
    
    disp(['Determining upper bounds for parameter ', num2str(i)]);

    factor = 1+change(step)/100;
    val2 = val1*factor;
    beq(end) = val2;
    x0 = xopt;
    x0(i) = val2;

    if(val2 > ub(i))
        val1 = ub(i);
        clear val2
        break
    elseif(val2-val1 < 0.1)
        clear val2
        break
    end
    
    % Define objective function
    stoptime = datetime('now')+time_limit;
    obj_ode = @(x)calc_l2_obj(x,tspan,time_points, S, knot_pos, bspline_order,n_c_param,c_param_index,n_param,...
                        MID_input,input_data_metabs, SD_input, MID_balance,balance_data_metabs,SD_balance,...
                        conc_ratio, conc_sd,conc_metab_list,init_conc, balance_metabs, flux_param_index,l2,stoptime);

    [~,objval,exitcode] = knitro_nlp(obj_ode, x0,[],[], Aeq,beq,lb,ub,[],[],kn_opt_E,'kn_ode.opt');

    if(sum(exitcode == ok_codes) && objval <= thres && objval ~=0)
        val1 = val2;
        clear val2
    else
        clear val2
        step = step+1;
        if(step > numel(change))
            flag = 0;
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

    factor = 1-change(step)/100;
    val2 = val1*factor;
    beq(end) = val2;
    x0 = xopt;
    x0(i) = val2;

    if(val2 < lb(i))
        val1 = lb(i);
        clear val2
        break
    elseif(val1-val2 < 0.1)
        clear val2
        break
    end
    
    % Define objective function
    stoptime = datetime('now')+time_limit;
    obj_ode = @(x)calc_l2_obj(x,tspan,time_points, S, knot_pos, bspline_order,n_c_param,c_param_index,n_param,...
                        MID_input,input_data_metabs, SD_input, MID_balance,balance_data_metabs,SD_balance,...
                        conc_ratio, conc_sd,conc_metab_list,init_conc, balance_metabs, flux_param_index,l2,stoptime);

    [~,objval,exitcode] = knitro_nlp(obj_ode, x0,[],[], Aeq,beq,lb,ub,[],[],kn_opt_E,'kn_ode.opt');

    if(sum(exitcode == ok_codes) && objval <= thres && objval ~=0)
        val1 = val2;
        clear val2
    else
        clear val2
        step = step+1;
        if(step > numel(change))
            flag = 0;
        end
    end

end
cil = val1;

save([folder,'\parameter_bounds_',num2str(i),'.mat'],"cil","ciu")

end

function obj = calc_l2_obj(x,tspan,time_points, S, knot_pos, bspline_order,n_c_param,c_param_index,n_param,...
                        MID_input,input_data_metabs, SD_input, MID_balance,balance_data_metabs,SD_balance,...
                        conc_ratio, conc_sd,conc_metab_list,init_conc, balance_metabs, flux_param_index,l2,stoptime)

    [obj1,obj2,obj3] = calc_obj_timed(x,tspan,time_points, S, knot_pos, bspline_order,n_c_param,c_param_index,n_param,...
                        MID_input,input_data_metabs, SD_input, MID_balance,balance_data_metabs,SD_balance,...
                        conc_ratio, conc_sd,conc_metab_list,init_conc, balance_metabs, flux_param_index,stoptime);
    
    CP = reshape(x(flux_param_index),size(S,2),[]);
    CPr = CP./repmat(CP(:,1),1,size(CP,2));
    CPr = CPr-1;
   
    obj = obj1+obj2+obj3+l2*sumsqr(CPr);

end