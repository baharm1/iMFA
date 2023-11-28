function [x,fval,exitflag] = estimate_model_param(Aeq,beq,tspan,time_points,n_flux_param,n_c_param,n_isotopomer_frac,...
                         MID_balance, balance_data_metabs, SD_balance, balance_metabs, ...
                     MID_input,input_data_metabs, SD_input, ...
                      conc_metab,conc_values,lb,ub,niter,x0,varargin)
    

% Default values
ncores = 1;

% Optional Arguments
if ~isempty(varargin)
    switch (varargin{1})
        case 'parallel'
            ncores = varargin{2};
        otherwise
             error(['Unexpected option: ' varargin{1}])
    end
end


% define objective function
obj_ode = @(x)calc_obj(x,tspan,time_points,n_flux_param,n_c_param,n_isotopomer_frac,...
                           MID_balance, balance_data_metabs, SD_balance, balance_metabs, ...
                           MID_input,input_data_metabs, SD_input, ...
                           conc_metab,conc_values);
% solve the system of equations
kn_opt_E = knitro_options('UseParallel','true');

if(ncores == 1)
    for iter = 1:niter
        disp(['Starting optimization for iteration ',num2str(iter)])
        [ x(:,iter), fval(iter), exitflag(iter)] = knitro_nlp(obj_ode, x0(:,iter),[],[], Aeq,beq,lb,ub,[],[],kn_opt_E,'kn_ode.opt');
%         x(:,iter)
    end
else
    parfor (iter = 1:niter,ncores)
        disp(['Starting optimization for iteration ',num2str(iter)])
        [ x(:,iter), fval(iter), exitflag(iter)] = knitro_nlp(obj_ode, x0(:,iter),[],[], Aeq,beq,lb,ub,[],[],kn_opt_E,'kn_ode.opt');
    end
end



