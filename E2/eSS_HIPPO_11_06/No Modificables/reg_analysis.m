%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reg_analysis
% Performs al the regression analysis for a parameter set: significance,
% sensitivity and identifiability analysis, and BIC and AICc calculations.
% Returns all the different regression analysis outputs.
%
% Benjamín J. Sánchez
% Last Update: 2014-07-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [AICc,BIC,CI,CC,Mc,sensib_prom,diff] = reg_analysis(k)

%Load experimental data:
texp     = evalin('base','texp');
ydata    = evalin('base','ydata');
[~,xmod] = solve_ODE(k,texp);
xmod(:,2) = [];
xmod(:,3) = [];
xmod(:,3) = [];
xmod(:,3) = [];
xmod(:,3) = [];
[~,yexp] = obj_var(xmod,ydata);

%Jacobian and residual calculations using lsqcurvefit:
kL = k.*0.9999999999;
kU = k.*1.0000000001 + ones(1,length(k))*1e-15;
[k2,~,res,~,~,~,Jac] = lsqcurvefit(@lsq_func,k,texp,yexp,kL,kU);
diff = sum((k-k2).^2);

%AIC & BIC calculations:
SS    = sum(sum(res.^2));
m     = length(k2);
[N,~] = size(res);
AIC   = N*log(SS/N)+2*(m+1);
AICc  = AIC + 2*(m+1)*(m+2)/(N-m-2);
BIC   = SS/N + m*log(N);

%Confidence intervals:
[CI,CC] = intconfianza(Jac,res,k2,0.05);
CC      = CC./100;
disp('Insignificant Parameters:');
for i = 1:m
    if k2(i) == 0
        CC(i) = 0;
    elseif CC(i) >= 2
        disp(num2str(i));
    end
end

%Identifiability analysis:
x0 = evalin('base','x0');
U  = evalin('base','U');
Mc = identificaBSB(1:length(x0),k,x0,texp,@model,U,0);
set(0,'ShowHiddenHandles','on');delete(get(0,'Children'))

%Sensitivity analysis:
ksensibilidadBSB(1:length(x0),k,x0,texp,@model,3,0)
sensib_prom = load('GraficoBarra.txt');
set(0,'ShowHiddenHandles','on');delete(get(0,'Children'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%