function cout = generate_CI_vectors(i,xopt,up,low,Aeq,beq,tspan,time_points,n_c_param,...
                         MID_balance, balance_data_metabs, SD_balance, ...
                        MID_input,input_data_metabs, SD_input,lb,ub,...
                        c_param_index,n_param, conc_ratio, conc_sd, conc_metab_list, init_conc,...
                        S, knot_pos, bspline_order,balance_metabs, flux_param_index,l2)

cout = 0;

global folder

% Update Aeq and beq to constrain the parameters
param_cons = zeros(1,size(Aeq,2));
param_cons(i) = 1;
Aeq = [Aeq;param_cons];
beq = [beq;0];

% Acceptable exit codes
ok_codes = [0,-100,-101,-102,-103,-400,-401,-402];

% Define the non-linear contraints function
NLcon = @(x)NL_cons(x, S, c_param_index, flux_param_index, bspline_order, knot_pos, tspan);

% knitro options
kn_opt_E = knitro_options('UseParallel','true');

% Time limit
time_limit = minutes(20);

%% Assign the values to test

sampleup = 1;
diffup = abs(xopt(i)-up);
if(diffup == 0)
    sampleup = 0;
    testup = [];
elseif(diffup < 0.5)
    testup = up;
elseif(diffup < 5)
    testup = linspace(xopt(i),up,round(diffup/0.5));
elseif(diffup < 20)
    testup = linspace(xopt(i),up,round(diffup));
else
    testup = linspace(xopt(i),up,round(diffup/2));
end


sampledown = 1;
diffdown = abs(xopt(i)-low) ;
if(diffdown == 0)
    sampledown = 0;
    testdown = [];
elseif(diffdown < 0.5)
    testdown = low;
elseif(diffdown < 5)
    testdown = linspace(low,xopt(i),round(diffdown/0.5));
elseif(diffdown < 20)
    testdown = linspace(low,xopt(i),round(diffdown));
else
    testdown = linspace(low,xopt(i),round(diffdown/2));
end

if(sampledown == 0 && sampleup == 0)
    return
end

%% Sample the space between the optimum and the upper limit

xsampleup = [];

for k = 1:numel(testup)
    
    clear xflux exitcode

    beq(end) = testup(k);
    x0 = xopt;
    x0(i) = testup(k);

    % Define objective function
    stoptime = datetime('now')+time_limit;
    obj_ode = @(x)calc_l2_obj(x,tspan,time_points, S, knot_pos, bspline_order,n_c_param,c_param_index,n_param,...
                        MID_input,input_data_metabs, SD_input, MID_balance,balance_data_metabs,SD_balance,...
                        conc_ratio, conc_sd,conc_metab_list,init_conc, balance_metabs, flux_param_index,l2,stoptime);

    [xflux,objval,exitcode] = knitro_nlp(obj_ode, x0,[],[], Aeq,beq,lb,ub,[],[],kn_opt_E,'kn_ode.opt');

    if(sum(exitcode == ok_codes) && objval ~=0)
        xsampleup = [xsampleup,xflux];
    end

end

%% Sample the space between the optimum and the lower limit

xsampledown = [];

for k = 1:numel(testdown)
    
    clear xflux exitcode

    beq(end) = testdown(k);
    x0 = xopt;
    x0(i) = testdown(k);

    % Define objective function
    stoptime = datetime('now')+time_limit;
    obj_ode = @(x)calc_l2_obj(x,tspan,time_points, S, knot_pos, bspline_order,n_c_param,c_param_index,n_param,...
                        MID_input,input_data_metabs, SD_input, MID_balance,balance_data_metabs,SD_balance,...
                        conc_ratio, conc_sd,conc_metab_list,init_conc, balance_metabs, flux_param_index,l2,stoptime);

    [xflux,objval,exitcode] = knitro_nlp(obj_ode, x0,[],[], Aeq,beq,lb,ub,[],[],kn_opt_E,'kn_ode.opt');

    if(sum(exitcode == ok_codes) && objval ~=0)
        xsampledown = [xsampledown,xflux];
    end

end

%% Save the data

vectors = [xsampleup ,xsampledown];
save([folder,'\feasible_vectors\parameter_',num2str(i),'.mat'],"vectors");


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