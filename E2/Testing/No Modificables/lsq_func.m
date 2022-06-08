%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lsq_func
% Integrates the dynamic model and returns the model predictions. To be
% used for lsqcurvefit.
%
% Benjamín J. Sánchez
% Last Update: 2013-06-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ymod = lsq_func(k,texp)

%Load experimental data:
ydata = evalin('base','ydata');

%Resolution of ODEs:
[~,xmod] = solve_ODE(k,texp);
xmod(:,2) = [];
xmod(:,3) = [];
xmod(:,3) = [];
xmod(:,3) = [];
xmod(:,3) = [];

%Return de normalized measured variables:
[ymod,~] = obj_var(xmod,ydata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%