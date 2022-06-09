%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iteration
% Performs a whole iteration of the procedure: Parameter estimation using
% Scatter Search (SSm) and regression analyisis. Saves all the results in
% an appropiate folder, and decides wich parameters should be fix in the
% next iteration. Returns the parameters to be fixed (ktofix) and the
% iteration results (it_results).
%
% Benjamín J. Sánchez
% Last Update: 2014-07-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [CI,CC,k_SSm] = iteration_mod(kfixed)

close all

%=================== PROBLEM SPECIFICATIONS ===========================
problem.f = 'obj_func'; %mfile containing the objective function

kL = evalin('base','kL');
k0 = evalin('base','k0');
kU = evalin('base','kU');

%Decide the parameters to be estimated depending on kfixed:
m = length(kfixed);
j = 1;
for i = 1:m
    if isnan(kfixed(i))
        problem.x_L(j) = kL(i);
        problem.x_0(j) = k0(i);
        problem.x_U(j) = kU(i);
        j = j+1;
    end    
end

%Load and modify if required the SSm options
opts = evalin('base','opts_SSm');

%================================ DATA ================================

texp      = evalin('base','texp');
ydata     = evalin('base','ydata');

%============================ OPTIMIZATION ============================

assignin('base','kfixed',kfixed);
% Results=ssm_kernel(problem,opts,texp,ydata);
Results=ess_kernel(problem,opts,texp,ydata);
%============================= RESULTS ================================

k_SSm    = Results.xbest;
%========================= REGRESSION ANALYSIS ========================

close all
[~,~,CI,CC] = reg_analysis_mod(k_SSm);


