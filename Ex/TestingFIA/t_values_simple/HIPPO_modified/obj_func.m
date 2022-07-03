%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obj_func
% Integrates the dynamic model and calculates afterwards the cuadratic 
% difference between the model predictions and the experimental data,
% normalized by the maximum value of the experimental data. Returns the
% cuadratic difference (objective function). To be used in the parameter
% estimation with SSm.
%
% Benjamín J. Sánchez
% Last Update: 2013-06-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [J,g,R]=obj_func(k,~,ydata)

Operational_Data = evalin('base','int_time');
Measure_Ind      = evalin('base','Measure_Ind');
time_dap         = evalin('base','time_dap');
FDA_add_idx      = evalin('base','FDA_add_idx');
x0               = evalin('base','x0');

%Definicion del tiempo de integracion en base a datos operacionales de seguimiento de temperatura
int_time         = Operational_Data(:,1);

%Definimos dos tiempos de integracion: Previo a la adicion de FDA y posterior a ello.
t1        = int_time(int_time<time_dap);
t2        = int_time(int_time>=time_dap);

%Opciones de integracion
options = odeset('RelTol',1e-3,'AbsTol',1e-3);

%Integracion de la seccion previa a la adicion de FDA
[T1,X1]   = ode15s(@model,t1,x0,options,k);

%Re-definimos el punto inicial para la segunda seccion. (NOTA: El Nitrogeno fue medido tras la adicion de FDA, por ende se utiliza ese dato experimental para definir x02)
x02       = X1(end,:);  x02(2) = ydata(FDA_add_idx,1);

%Integracion de la seccion posterior a la adicion de FDA
[T2,X2]   = ode15s(@model,t2,x02,options,k);

%Concatenacion de los resultados obtenidos
T         = [T1;T2];
X         = [X1;X2];

%Seleccion de los puntos objetivos (es decir, aquellos donde calzan los tiempos operacionales con los experimentales)

X         = X(ismember(T,int_time),:);            %Eliminacion de los puntos extras generados de la integracion
total_idx = (1:1:length(X))';                     %Vector de indices
X_sim     = X(ismember(total_idx,Measure_Ind),:); %X_sim contiene los datos simulados que calzan con los puntos de muestreo experimental, por ende son comparables para el calculo de error.

if length(X_sim(:,2)) == length(ydata(:,1)) %Condicional para evitar que eSS entregue soluciones inestables
    
    %Definimos los puntos que queremos utilizar para el cálculo del error. Excluiremos el punto inicial y el de adición de FDA. Para el N
    %excluiremos los puntos posteriores a la adición de FDA dado el ruido en las mediciones. 
    S_data                  = ydata([2:FDA_add_idx-1 FDA_add_idx+1:end],2);      S_sim                  = X_sim([2:FDA_add_idx-1 FDA_add_idx+1:end],4);
    N_data                  = zeros(size(S_data));                               N_sim                  = zeros(size(S_data));
    N_data(2:FDA_add_idx-1) = ydata(2:FDA_add_idx-1,1);                          N_sim(2:FDA_add_idx-1) = X_sim(2:FDA_add_idx-1,2).*1000;                                                                                      %N_pond    = 1/(sum(Kinetic_Matrixes{i,j}(1:FDA_add_idxs{i,j}+1,4))/sum(Kinetic_Matrixes{i,j}(:,3))); %This ponderator will weigth up N adjust to the same level as S related error.

    %Construction of residual matrix:
    R    = [(N_data-N_sim)./max(N_data) (S_data-S_sim)./max(S_data)]; 
    
    %Calculation of the objective function:
    J = sum(sum(R.^2)); disp(J)
    g = 0;
    R=reshape(R,numel(R),1);

else    %Condicional para evitar que eSS entregue soluciones inestables
    R    = NaN;
    J    = NaN;
    g    = NaN;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%