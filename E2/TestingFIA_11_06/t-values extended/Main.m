%% 
global exp_temp int_time Kinetic_Matrix FDA_add_idx x0 time_dap

%Load data and integration options
[kL,k0,kU,opts_SSm,texp,ydata,x0,solver_ODE,opts_ODE,T,U,int_time,time_dap,exp_temp,Kinetic_Matrix,FDA_add_idx,Measure_Ind] = load_problem_zent;

%Determine the amount of iterations to evaluate
iterations      = 10;

%Select parameters to be fixed (manually, realisation-wise)
fixed_params{1} = [];
fixed_params{2} = [9]; % kd0
fixed_params{3} = [3 9]; % beta_F0 kd0
fixed_params{4} = [7 9 12]; % Kig0 kd0 Yxf
fixed_params{5} = [7]; % Kig0
fixed_params{6} = [7 11]; % Kig0 Yxg
fixed_params{7} = [7 11 12]; % Kig0 Yxg Yxf

realisations    = length(fixed_params);

for u = 1:realisations
    kfixed = nan(14,1);
    for i = 1:length(fixed_params{u})
        j = fixed_params{u}(i);
        kfixed(j) = k0(j);
    end

    T_V = [];

    for n = 1:iterations
        %Calculate optimal fitting parameters, CC and CI
        [CI,CC,k_SSm,J_SSm] = iteration_mod(kfixed);

        %Calculate t-values for this model structure
        t_values = 4*k_SSm'./CI;

        T_V     = [T_V t_values];
        ins{n}  = find(CC>2);
        f_best{n}= J_SSm;
    end
    T_VALUES{u} = T_V;
    INSIGN{u}   = ins;
    F_BEST{u}   = f_best; 
end
%% Graphics

%X-Tick lables
all_params = {'\mu_0','\beta_{G0}','\beta_{F0}','K_{N0}','K_{G0}','K_{F0}','K_{iG0}','K_{iE0}','k_{d0}','Y_{XN}','Y_{XG}','Y_{XF}','Y_{EG}','Y_{EF}'};
INS_COUNTS = zeros(length(k0),realisations);

%Matrix with the problematicness frequency adapted to each model structure
for u = 1:realisations
    indxs = find(~ismember(1:14,fixed_params{u}));
    count = zeros(14,1);
    for n = 1:iterations
        for i = 1:length(k0)
            if ismember(i,INSIGN{u}{n})
                count(indxs(i)) = count(indxs(i))+1;
            end
        end
    end
    INS_COUNTS(:,u) = count;
end

%Loop to filter results, leaving only problematic ones for display
tau = zeros(size(k0'));
contador = 1;
for i = 1:length(k0)
    tau(i) = any(INS_COUNTS(i,:));
    if tau(i) == true
        Prob_params{contador} = all_params{i};
        contador = contador + 1;
    end
end

%Xtick label generation
Xlabels = cell(1,realisations);

for u = 1:realisations
    if isempty(fixed_params{u})
        Xlabels{u} = 'None';
    else
        str = '';
        for i = 1:length(fixed_params{u})
            if i == 1
                str = all_params{fixed_params{u}(i)};
            else
                str = strcat(str,' + ', all_params{fixed_params{u}(i)});
            end
        end
        Xlabels{u} = str;
    end
end

% Bar plot 

bar(INS_COUNTS(find(tau),:)')
legend(Prob_params)
xticklabels(Xlabels)
xlabel('Realisations')
ylabel('Problematicness frequency')
