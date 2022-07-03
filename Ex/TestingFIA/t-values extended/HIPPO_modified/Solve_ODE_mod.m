function [T,Xf] = Solve_ODE_mod(k)    

Operational_Data = evalin('base','int_time');
time_dap         = evalin('base','time_dap');
FDA_add_idx      = evalin('base','FDA_add_idx');
x0               = evalin('base','x0');
ydata            = evalin('base','ydata');

%Definicion del tiempo de integracion en base a datos operacionales de seguimiento de temperatura
int_time         = Operational_Data(:,1);

%Error function definition for parameters estimation
t1        = int_time(int_time<time_dap);
t2        = int_time(int_time>=time_dap);

%Pre DAP addition section
[T1,X1]   = ode15s(@model,t1,x0,[],k);
x02       = X1(end,:);  x02(2) = ydata(FDA_add_idx,1);

%After DAP addition section
[T2,X2]   = ode15s(@model,t2,x02,[],k);

%Resulting data consolidation
T         = [T1;T2];
Xf        = [X1;X2];
Xf        = Xf(:,[2 3 4]); %IMPORTANTE: AQUI ALTERAMOS LOS INDICES PARA QUE SEAN LOS ESTADOS OBSERVADOS. DEBE SER EN ORDEN ASCENDETE [2 4] NO [4 2].
end