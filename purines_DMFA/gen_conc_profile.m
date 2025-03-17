function [t,C] = gen_conc_profile(x, S, c_param_index, flux_param_index, bspline_order, knot_pos,tspan,posflag)

c0 = x(c_param_index);
CP = reshape(x(flux_param_index),size(S,2),[]);

conc_odefn = @(t,y)solveConc(t,y,knot_pos, tspan(end), S, CP, bspline_order);

if(posflag == 1)  
    odeOptions = odeset('NonNegative',1);
    [t,C] = ode23tb(conc_odefn,tspan,c0, odeOptions);
else
    [t,C] = ode23tb(conc_odefn,tspan,c0);
end

end

function dcdt = solveConc(t,y,knot_pos, maxT, S, CP, bspline_order)

NoutCol = bSplineMat_lite(knot_pos,t/maxT,bspline_order);
dcdt = S*CP*NoutCol;

end

