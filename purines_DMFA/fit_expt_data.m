function [xopt,fval,exitE] = fit_expt_data(x0,tspan,time_points, ...
                        MID_input,input_data_metabs, SD_input, ...
                        MID_balance,balance_data_metabs,SD_balance,...
                        n_c_param,c_param_index,n_param, conc_ratio, conc_sd, conc_metab_list, init_conc,...
                        S, knot_pos, bspline_order,...
                        balance_metabs, flux_param_index,lb,ub,Aeq,beq,l2)


%% Optimization

% Define objective function
obj_ode = @(x)calc_l2_obj(x,tspan,time_points, S, knot_pos, bspline_order,n_c_param,c_param_index,n_param,...
                        MID_input,input_data_metabs, SD_input, MID_balance,balance_data_metabs,SD_balance,...
                        conc_ratio, conc_sd,conc_metab_list,init_conc, balance_metabs, flux_param_index,l2);

% Define the non-linear contraints function
NLcon = @(x)NL_cons(x, S, c_param_index, flux_param_index, bspline_order, knot_pos, tspan);

% knitro options

kn_opt_E = knitro_options('UseParallel','true');


% Run the optimization
[xopt,fval,exitE]  = knitro_nlp(obj_ode, x0,[],[], Aeq,beq,lb,ub,NLcon,[],kn_opt_E,'kn_ode.opt');

end

function obj = calc_l2_obj(x,tspan,time_points, S, knot_pos, bspline_order,n_c_param,c_param_index,n_param,...
                        MID_input,input_data_metabs, SD_input, MID_balance,balance_data_metabs,SD_balance,...
                        conc_ratio, conc_sd,conc_metab_list,init_conc, balance_metabs, flux_param_index,l2)

    [obj1,obj2,obj3] = calc_obj(x,tspan,time_points, S, knot_pos, bspline_order,n_c_param,c_param_index,n_param,...
                        MID_input,input_data_metabs, SD_input, MID_balance,balance_data_metabs,SD_balance,...
                        conc_ratio, conc_sd,conc_metab_list,init_conc, balance_metabs, flux_param_index);
    
    CP = reshape(x(flux_param_index), size(S, 2), []);
    CPr = CP./repmat(CP(:, 1), 1, size(CP, 2));
    CPr = CPr-1;
   
    obj = obj1+obj2+obj3+l2*sumsqr(CPr);

end