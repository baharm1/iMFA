% This script fits the glucose enrichment-time data to the model and 
% calculates the saturation enrichment:
% FE = (P-r/k1)(H1/Q)exp(-k1t) + (P-r/k2)(H2/Q)exp(-k2t) + (r/Q)(H1/k1 + H2/k2) 
% FE = fractional enrichment
% (r/Q)(H1/k1 + H2/k2) = saturation enrichment = r/F
% r = rate of tracer infusion
% F = glucose consumption flux
% k1,k2 = rate constants 
% Q = total glucose pool in wt.glucose/wt.tissue
% P = bolus dose
% t = time
% H1+H2 = 1
% reference: chapter 9, Tracer methods for in vivo kinetics by Shipley and
% Clark
% The value of k is assumed to be constant between the mice, Q is assumed
% to vary between the mice which results in different saturation enrichment

%% User defined inputs
clc; clear;
r = 0.012; % rate of tracer infusion. in mg/min/g
P = 0.4; % Bolus dose. in mg/min/g
filename = 'plasma_glucose_ME_4h.xlsx'; % File containing the data for mice inflused for atleast 120 min (to estimate k1,k2,H1)
file2 = 'plasma_glucose_ME_all.xlsx'; % File containing all input data
opfilename = 'avgSE_2comp.xlsx'; % Output file name
plotfolder = 'avgSE_fit_2comp'; % Folder to save plots
min_ss = 0.35; % The minimum SS enrichment observed in the data
t0_max = 0.6; % The maximum enrichment from bolus at t = 0
t0_min = 0.3;
max_tp = 240;

%% Read and fit the 4hr data
data = readtable(filename,"ReadVariableNames",true);
nsample = size(data,2)-1;

nparam = nsample+3;
lb = zeros(nparam,1);
ub = 10*ones(nparam,1);
lb(3) = 0.4;
ub(3) = 1;
lb(4:end) = repmat(P/t0_max,numel(lb(4:end)),1);
ub(4:end) = repmat(P/t0_min,numel(lb(4:end)),1);
niter = 40;
x0 = repmat(lb,1,niter) + rand(nparam,niter).*repmat(ub-lb,1,niter);

obj_fn = @(x)calc_obj(x,data,nsample,r,P);
NLcon = @(x)NL_cons(x,r,min_ss);

for k = 1:niter
    [x(:,k),fval(k),exitflag(k)] = fmincon(obj_fn,x0(:,k),[],[],[],[],lb,ub,NLcon);
end

[~,feas_idx] = find(exitflag >= 0);
fval_feas = fval(feas_idx);
[~,min_idx] = min(fval_feas);
xfeas = x(:,feas_idx);
param_fit = xfeas(1:3,min_idx);

%% Estimate the saturation enrichment for all samples

% Read the data
data2 = readtable(file2,"ReadVariableNames",true);
samples = data2.Properties.VariableNames(2:end);

k1 = param_fit(1);
k2 = param_fit(2);
H1 = param_fit(3);
H2 = 1-H1;
modelfn = @(q,t)(P-r/k1)*(H1/q)*exp(-k1*t) + (P-r/k2)*(H2/q)*exp(-k2*t) + (r/q)*(H1/k1 + H2/k2);

Qfit = zeros(numel(samples),1);
q0 = median(xfeas(4:end,min_idx));

for k = 1:numel(samples)
    y = table2array(data2(:,k+1));
    flag = isnan(y);
    y = y(~flag);
    y = y/100;
    tps = table2array(data2(~flag,1));  
    mdl = fitnlm(tps,y,modelfn,q0);
    Qfit(k) = mdl.Coefficients.Estimate;
    clear mdl y tps
end

% Estimate the saturation enrichment
SE = (r./Qfit).*(H1/k1 + H2/k2);

% Save the data
T = table(samples',repmat(param_fit',numel(samples),1),Qfit,SE);
T.Properties.VariableNames ={'SampleName','k1_k2_H1','Q','SE'};
writetable(T,opfilename);

%% Plot
mkdir(plotfolder) 
for k = 1:numel(samples)
    y = table2array(data2(:,k+1));
    flag = isnan(y);
    y = y(~flag);
    y = y/100;
    tps = table2array(data2(~flag,1));  
    tspan = [0:1:max_tp];
    figure;
    scatter(tps,y);
    hold on;
    plot(tspan,modelfn(Qfit(k),tspan));
    xlabel('time (min)');
    ylabel('enrichment');
    hold off;
    saveas(gcf,[plotfolder,'/',char(samples(k)),'.png']);
    close;
    clear y tps
end

%% Objective and NL function

function obj = calc_obj(x,data,nsample,r,P)

    k1 = x(1);
    k2 = x(2);
    H1 = x(3);
    H2 = 1-H1;
    Q = x(4:end);
    obj = 0;
    for i = 1:nsample
        sdata = table2array(data(:,i+1));
        flag = isnan(sdata);
        sdata = sdata(~flag);
        sdata = sdata/100;
        tps = table2array(data(~flag,1));

        cdata = calc_eqn(tps,k1,k2,H1,H2,Q(i),r,P);
        obj = obj+sumsqr((cdata-sdata));
        
    end
  
end

function cdata = calc_eqn(tps,k1,k2,H1,H2,q,r,P)
    cdata = zeros(numel(tps),1);
    for t = 1:numel(tps)
        ti = tps(t);
        cdata(t) =(P-r/k1)*(H1/q)*exp(-k1*ti) + (P-r/k2)*(H2/q)*exp(-k2*ti) + (r/q)*(H1/k1 + H2/k2);
    end

end


function [C,Ceq] = NL_cons(x,r,min_ss)

    Ceq = [];
    thres = r*(1-min_ss)/min_ss;
    Q = x(4:end);
    k1 = x(1);
    k2 = x(2);
    H1 = x(3);
    H2 = 1-H1;
    C = repmat(thres,numel(Q),1)-Q/(H1/k1 - H2/k2);

end