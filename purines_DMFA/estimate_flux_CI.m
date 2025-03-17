function [flux_cil,flux_ciu,flux_optim] = estimate_flux_CI(xopt,tspan,S,knot_pos,bspline_order,flux_param_index)

global folder

% Calculate the optimum flux vector
CP = reshape(xopt(flux_param_index),size(S,2),[]);
flux_optim = generate_flux_profile(tspan,CP,knot_pos,bspline_order);

% Initialize the matrices
flux_cil = flux_optim;
flux_ciu = flux_optim;

% Read the optimized vectors and generate their flux profiles
% Change the minimum and maximum flux profile accordingly
for i = flux_param_index(1):flux_param_index(end)

    try
        clear vectors
        load([folder,'\feasible_vectors\parameter_',num2str(i),'.mat'],"vectors");
    catch
        disp(['File for paremeter ',num2str(i),' not found']);
        continue
    end
    
    disp(['File for paremeter ',num2str(i),' found']);

    for k = 1:size(vectors,2)
        clear CP2
        % Generate the flux profile
        CP2 = reshape(vectors(flux_param_index,k),size(S,2),[]);
        new_profile = generate_flux_profile(tspan,CP2,knot_pos,bspline_order);
        flux_cil = min(flux_cil,new_profile);
        flux_ciu = max(flux_ciu,new_profile);
    end

end


end

function v = generate_flux_profile(tspan,CP,knot_pos,bspline_order)

    v = zeros(size(CP,1),numel(tspan));
    for e = 1:length(tspan)
        time = tspan(e);
        Nout = bSplineMat_lite(knot_pos,time/tspan(end),bspline_order);
        v(:,e) = CP*Nout;
    end

end