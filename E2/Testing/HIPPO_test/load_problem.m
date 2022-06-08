%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load_problem
% Define all the needed data and options to solve the problem
%
% OUTPUTS:
% kL            Row vector with the parameter lower bounds for SSm.
% k0            Row vector with the parameter initial values for SSm.
% kU            Row vector with the parameter upper bounds for SSm.
% opts_SSm      Options for the SSm procedure
% texp          Experimental times for the integration.
% ydata         Experimental data. Be consistent with obj_var.
% x0            Initial state variable vaues for the integration.
% solver_ODE    Solver for the dynamic model integration
% opts_ODE      Options for the dynamic model integration (odeset)
% T             Threshold for correlation values.
%
% Benjamín J. Sánchez: 2014-07-17
% Javiera Pérez: 17-05-2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [kL,k0,kU,opts_SSm,texp,ydata,x0,solver_ODE,opts_ODE,T,U] =...
    load_problem

%Lower, initial and upper values for each parameter in SSm:%
%  P = [ p1      p2      p3    p4    p5   p6   p7   p8    p9    p10   p11   p12   p13    p14    p15   p16   p17   p18]
kL   = [   1e-3     1e-3        1e-3     1e-3      1e-3     1e-3      1e-3      1e-3       1      1     1     1      1      1     1      1e-3       1e-3     1e-3];
k0   = [ 0.755  0.039  0.314  0.519  0.272  0.89  1.116  1.482  12.22  23.08  9.75  37.63  83.64  44.7  40.11  0.716  0.045  0.047];
kU   = [   1      1      1      1      1      1     3      3       30     45    20    75      150   80    80      2      1      1] ;

%SSm options:
opts_SSm.maxeval      = 10000;
opts_SSm.local.n1     = 100;
opts_SSm.maxtime      = 100;
opts_SSm.strategy     = 3;
opts_SSm.local.solver = 'fmincon';
opts_SSm.local.finish = 'fmincon';
opts_SSm.combination  = 1;  

%Experimental Data:
% [Time X S P]
data  = load('datos.txt'); %edit here
[~,n] = size(data); 
texp  = data(:,1)';
ydata = data(:,2:n);

%Initial conditions for integration:
%      S0   Xf0  
x0 = [63 0.853 0.15 0 0 0.5 0.075] ;

%ODE options:
solver_ODE = 'ode15s';
opts_ODE   = odeset('RelTol',1e-4,'AbsTol',1e-7);

%Threshold for correlations (any couple of parameters with a correlation
%higher than this value will be fixed):
T = 0.95;

%Threshold for iterations (criterion III):
U = 1.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%