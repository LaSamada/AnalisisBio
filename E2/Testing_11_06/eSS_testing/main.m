%clear mex; 
%clear all; 
close all;

%====================== PROBLEM DEPENDENT DECLARATIONS=====================

%=================== END OF PROBLEM DEPENDENT DECLARATIONS ================

%========================= PROBLEM SPECIFICATIONS ===========================
problem.f='obj_func'; %mfile containing the objective function

% Model parameters Thesis            % Paper
% p(1)=0.21;                         % mun=0.28;
% p(2)=0.02655;                      % kd=0.031;
% p(3)=0.13377/1000;                 % k=0.000133;
% p(4)=0.132158502;                  % mco2=0.11;
% p(5)=0.05448394856386;             % mo2=0.06;
% p(6)=0.001591591;                  % kp=0.00044;
% p(7)= 20.99023219980516;           % yxn = 20.8;
% p(8)= 0.00078465311766;            % kn = 0.0011;
% p(9)= 1.15424442432247;            % yxco2 = 0.58;
% p(10)= 2.60401780820540;           % yxo2 = 2.11;
% p(11)= 0.000611044;                % betam = 0.00065;
% p(12)= 2.9518533148e5;             % ki = 7.86e5;
% p(13)= 0.89365003128734;           % yxs= 1.21;
% p(14)= 0.09125726087293;           % ms= 0.11;

%            1    2    3    4     5    6    7    8    9      10    11     12   13   14 
%            mun  kd   k    mco2  mo2  kp   yxn  kn   yxco2  yxo2  betam  ki   yxs  ms
%problem.x_L=[0   0    0    0     0    0    0    0    0      0     0      0    0    0];
% problem.x_L=[1e-6 1e-6 1e-6 1e-6  1e-6 1e-6 1e-6 1e-6 1e-6   1e-6  1e-6   1e-6 1e-6 1e-6];
% problem.x_U=[1    1e-1 1e-3 1     1    1e-3 1e2  1e-2 1e2    1e2   1e-2   1e6  1e1  1e1];
problem.x_L=[   1e-3     1e-3    1e-3                1e-3             1           1e-3       1e-3   ];
problem.x_U=[   1         1        3      3            45          2      1      ]  ;
problem.x_0=[ 0.755    0.314    0.519   1.116     23.08    0.716  0.045 ];

%opts.maxeval=100000; 
%opts.maxtime=1e7;
%opts.local.n1 = 10000;
opts.maxtime=100;
opts.strategy=3;
opts.local.solver='fmincon';

%Try using this option and not using it (you'll see a great difference)
opts.log_var=(1:7);              %Declare all variables as "logarithmic"

%========================= END OF PROBLEM SPECIFICATIONS ==================

%================================== DATA ==================================
% Time, Urea, Biomass, CO2, O2, GA3, Starch
load datos.txt -ascii;
[r1,c1]  = size(datos); 
texp     = datos(:,1);
yexp     = datos(:,2:c1);


%================================== OPTIMIZATION ==========================
%Results=ess_kernel(problem,opts,texp,yexp);
% Algorithms: VNS, eSS, CeSS
Results=MEIGO(problem,opts,'eSS',texp,yexp);

%================================== Results ===============================
x  = Results.xbest;
X1 = [63 0.853 0.15 0 0 0.5 0.075];
OPTIONS = odeset('RelTol',1e-6,'AbsTol',1e-7);
[tout,yout]     = ode23(@model_ess,texp,X1,[],x);

figure
plot(tout,yout)
hold on
plot(texp,yexp,'o')


