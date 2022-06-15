function [T,Xf] = ReSimulate_Zenteno_HIPPO(k)    
global int_time time_dap  Kinetic_Matrix x0 FDA_add_idx 

%Error function definition for parameters estimation
t1        = int_time(int_time<time_dap);
t2        = int_time(int_time>=time_dap);

%Pre DAP addition section
[T1,X1]   = ode23s(@model,t1,x0,[],k);
x02       = X1(end,:);  x02(2) = Kinetic_Matrix(FDA_add_idx,4)./1000;

%After DAP addition section
[T2,X2]   = ode23s(@model,t2,x02,[],k);

%Resulting data consolidation
T         = [T1;T2];
Xf        = [X1;X2];
Xf        = Xf(:,[2 3 4]);
end