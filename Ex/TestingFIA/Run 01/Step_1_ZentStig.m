
%Adapted from https://sourceforge.net/projects/minimal-output-sets/
%by Cristóbal Torrealba 2021

%Application of the fast structural identifiability asessment applyed
%to the Zenteno et al. (2010) wine fermentation model. 
%Instructions for this example:

    % Run step 1, wait until all variables are generated. If wanted,
    % adjust the number of simulations (NExp); takes few seconds with
    % Zenteno´s model
    %
    % Run step 2, wait until all simulations are performed. If wanted,
    % adjust the bounds associated to random nominal parameter values
    % perturbation; takes few minutes with Zenteno´s model.
    %
    % Run step 3, wait until all graphics are generated. If wanted,
    % change measured states (sensors); takes a minute or so. 
    % Enjoy.

%% Initial procedures

%Prepare computation for Zenteno model (NOTE: Parameter vector does not include initial conditions as unknown parameters!)
tic
ModelName='Samada';

%Initial conditions for integration
S       = 100.95; %Yeast dosis   (g/L)
N       = 1.4; %Initial YAN   (g/L)
Xf       = 0.792; %Initial Glu   (g/L)
Xl       = 0.047; %Initial Fru   (g/L)
CIT       = 10e-10;      %Initial Eth   (g/L)
V      = 3.03;
qn   = 0.075;
x0       = [S, N, Xf, Xl, CIT, V, qn];

%Here operational data corresponding to this example is loaded
load('Operational_Data.mat')

%Define final integration time (Tf, in hours for this example) and
%number of integration steps (Nt) 
Tf=30.3; 

Nt=450; 

%% Symbolic Toolbox Variable defintions

%Here, we define model states as x{i} for a total of i = 1..7 model states.
for i=1:7
    syms(sprintf('x%d',i),'real');
end

%Here, we define model parameters as th{i} for a total of i = 1..14 regression parameters.
for i=1:18
    syms(sprintf('th%d',i),'real');
end

%Additionally, we define variables t, U1(in this example, must temperature [K]), and U2 (logical value describing if yeast thermal death is ocurring) as real 
%symbolic values representing time and model inputs.
syms t real    
syms U real

%Once symbolic parameters and model states are defined, we build arrays to vectorize parameters and model states.
statesSym=[x1 x2 x3 x4 x5 x6 x7]';
thetaSym=[th2 th3 th4 th5 th12 th16 th17]'; %NOTE: DO NOT ADD FIXED PARAMETERS IN THETASYM!
%         mu0 betaG0 betaF0 Kn0 Kg0 Kf0 Kig0 Kie0 Kd0 Yxn  Yxg  Yxf  Yeg  Yef

%Define the amount of model states and parameters
nx=length(statesSym);

%Define initial conditions vector as a symbolic variable associated to numerical values.
thIC=sym(x0);
x0Sym=thIC';

%Define the amount of simulations performed in the Monte Carlo procedure
NExp=100; 

%% Nominal parameter values (In this example, we use those in Zenteno et al. (2010))

%NOTE: ¡Fixed parameter values must be removed!

thetaNom = [  %0.755                      %1) mu0
              0.039                     %2) betaG0
              0.314                     %3) betaF0
              0.519                      %4) Kn0
              0.272                      %5) Kg0
              %0.89                      %6) Kf0
              %1.116                      %7) Kig0
              %1.482                        %8) Kie0
              %12.22                   %9) Kd0
              %23.08                     %10) Yxn
              %9.75                      %11) Yxg
              37.63                      %12) Yxf
              %83.64                      %13) Yeg
              %44.7
              %40.11
              0.716
              0.045];
           %   0.047];   

nth=length(thetaNom);

%% MODEL DEFINITION

%Here fixed parameter values must be defined:
th18 = 0.047;
th10 = 23.08;
th7 = 1.116;
th15 = 40.11;
th14 = 44.7;
th9 = 12.22;
th8 = 1.482;
th1 = 0.755;
th13 = 83.64;
th6 = 0.89;
th11 = 9.75;
% th5 = 0.272;
% th6 = 0.89;
% th1 = 0.755;
% th12 = 37.63;
% th12 = 37.63;
% th15 = 40.11;
% th16 = 0.716;
% th7 = 1.116;
% th6 = 0.89;
% th9 = 12.22;
% th1 = 0.755;
% th14 = 44.7;

% Infijables
% 2, 8, 11, 14, 16 y 17

%NOTE: This is the ODE system representing the Zenteno et al. (2010) model using no abreviatons for terms such as specific growth rate (mu).

Xdot =  [
            0*(0/x6)-(th1*(1-(th17/x7))*(x1/(th7+x1+(x1^2)/th13))/th4+th2*(x1/(th8+x1+(x1^2)/th14))*(1-x4/(0.36*x3*(1-0.06)+x4))*(th11/(th11+x5))/th5+th3*U*(x1/(th9+x1+(x1^2/th15)))*(th12/(th12+x5))/th6)*x3-(x1/x6)*0;                          
            0*(0/x6)-th16*(x2/(th10+x2))*x3-(x2/x6)*0;                              
            th1*(1-(th17/x7))*(x1/(th7+x1+(x1^2)/th13))*x3-(x3/x6)*0;    % Glucose
            (th2*(x1/(th8+x1+(x1^2)/th14))*(1-x4/(0.36*x3*(1-0.06)+x4))*(th11/(th11+x5))+0.06*th1*(1-(th17/x7))*(x1/(th7+x1+(x1^2)/th13)))*x3-(x4/x6)*0;    
            th3*U*(x1/(th9+x1+(x1^2/th15)))*(th12/(th12+x5))*x3-(x5/x6)*0;                       
            0-0;
            th16*(x2/(th10+x2))-th1*(1-(th17/x7))*(x1/(th7+x1+(x1^2)/th13))*x7];

%% END MODEL DEFINITION AND AUXILIAR FUNCTION GENERATION
%(May take a few minutes to process)

%Here, gradients are calculated symbolically and automatically stored as separated Matlab scripts for later usage.
dfdxSym     = simplify(jacobian(Xdot,statesSym));
dfdthSym    = simplify(jacobian(Xdot,thetaSym));
dx0dthSym   = simplify(jacobian(x0Sym,thetaSym));

f = matlabFunction(Xdot,'vars',{t,statesSym,U,thetaSym},'file',ModelName);
dfdx = matlabFunction(dfdxSym,'vars',{t,statesSym,U,thetaSym},'file',['dfdx',ModelName]);
dfdth = matlabFunction(dfdthSym,'vars',{t,statesSym,U,thetaSym},'file',['dfdth',ModelName]);
IC = matlabFunction(x0Sym,'vars',{thetaSym},'file',['x0',ModelName]);
dICdth = matlabFunction(dx0dthSym,'vars',{statesSym,thetaSym},'file',['dICdth',ModelName]);
