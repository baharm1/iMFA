function [C,Ceq] = NL_cons(x, S, c_param_index, flux_param_index, bspline_order, knot_pos, tspan)

% The NL constraint function forces the conecntration to be positive at all
% times

posflag = 0;
Ceq = [];
[~,C] = gen_conc_profile(x, S, c_param_index, flux_param_index, bspline_order, knot_pos,tspan,posflag);
C = -1*C(:);
end