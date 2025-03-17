function [internal_knots, min_Rc] = optimize_knot_placement( flux0, conc0, tspan,time_points, ...
                                                            MID_input,input_data_metabs, SD_input, ...
                                                            MID_balance,balance_data_metabs,SD_balance,...
                                                            n_c_param,conc_ratio, conc_sd, conc_metab_list, init_conc,...
                                                            S, bspline_order,no_int_knots,balance_metabs,...
                                                            lbc,ubc,lbf,ubf,lbp,ubp,lbip,ubip,input_mid_array,ip_mid_cons,iter)
                                                           

% This function optimizes the knot placement of b-splines based on
% concentration data only 
% The fluxes are constrained to be stoichimetrically balanced at t = 0
% The initial concentration values from literature are used to calculate 
% the objective along with the concentration-time-ratio data

global folder
Aeq = [];
beq = [];

% Define the array of possible internal knots
test_knots = 0.2:0.05:0.8;

% Set the model parameters
c_param_index = [1:n_c_param];

% Initialization of variables
internal_knots = [];

% Iteratively decide the best knot position
for j = 1:no_int_knots
    
    n_int_knots = length(internal_knots)+1;

    % Create lower and upper bound vectors
    clear lb ub
    lb = [lbc;lbf;repmat(lbp,bspline_order+n_int_knots-1,1);lbip];
    ub = [ubc;ubf;repmat(ubp,bspline_order+n_int_knots-1,1);ubip];

    % Generate the initial guess
    clear x0
    x0 = [conc0; repmat(flux0,bspline_order+n_int_knots,1);input_mid_array];

    % Define flux param index
    clear flux_param_index
    flux_param_index = [n_c_param+1:n_c_param+size(S,2)*(bspline_order+n_int_knots)];
    
    clear Rc_values selected_idx
    % Optimize with different knot positions
    for k = 1:numel(test_knots)
        
        % Create the knot position array
        clear knot_pos test_int_knots
        test_int_knots = sort([internal_knots;test_knots(k)]);
        knot_pos = ones(1,2*bspline_order + n_int_knots);
        knot_pos(1:bspline_order) = 0;
        if(n_int_knots > 0)
            knot_pos(bspline_order+1:(bspline_order+n_int_knots)) = test_int_knots;
        end
        

        % Run the optimization
        [~,Rc_values(k),~] = fit_expt_data(x0,tspan,time_points, ...
                        MID_input,input_data_metabs, SD_input, ...
                        MID_balance,balance_data_metabs,SD_balance,...
                        n_c_param,c_param_index,flux_param_index(end), conc_ratio, conc_sd, conc_metab_list, init_conc,...
                        S, knot_pos, bspline_order,...
                        balance_metabs, flux_param_index,lb,ub,Aeq,beq,0);

    end

   save([folder,'\iter',num2str(iter),'_knot',num2str(j),'.mat'],'internal_knots','Rc_values','test_knots');

    % Select the knot position with the best fit
     selected_idx = find(Rc_values == min(Rc_values));
     if(numel(selected_idx) > 1)
         positions = test_knots(selected_idx);
         dist = abs(positions - 0.65);
         dist_idx = find(dist == min(dist));
         if(numel(dist_idx) > 1)
             dist_idx = sort(dist_idx);
             dist_idx = dist_idx(1);
         end
         selected_idx = selected_idx(dist_idx);
     end
     clear new_knot
     new_knot = test_knots(selected_idx);

     % add the knot to the sequence
     internal_knots = [internal_knots, new_knot];
     internal_knots = sort(internal_knots);

     % Save the minimum Rc value
     if(j == no_int_knots)
         min_Rc = Rc_values(selected_idx);
     end

      % Update the test knots array
      test_knots = test_knots(test_knots > new_knot+0.15 | test_knots < new_knot-0.15);
  
    
end

end

