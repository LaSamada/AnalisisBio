global exp_temp int_time Kinetic_Matrix FDA_add_idx x0 time_dap

%Load data and integration options
[kL,k0,kU,opts_SSm,texp,ydata,x0,solver_ODE,opts_ODE,T,U,int_time,time_dap,exp_temp,Kinetic_Matrix,FDA_add_idx,Measure_Ind] = load_problem_zent;

%Select parameters to be fixed (manually)

fixed_params = [7 9 12]; %Indexes of parameters to be fixed 

kfixed = nan(14,1);
for i = 1:length(fixed_params)
    aux         = fixed_params(i);
    kfixed(aux) = k0(aux);
end

%Calculate optimal fitting parameters, CC and CI
[CI,CC,k_SSm] = iteration_mod(kfixed);

%Calculate t-values for this model structure
t_values = 4*k_SSm'./CI;

