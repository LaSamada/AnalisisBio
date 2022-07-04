%% codigo para leer datos obtenidos de hippo, pueden agregar o quitar valores segun lo deseen (solo es una base) 

load('it.mat'); % loads result matrix (it.mat generated by HIPPO)

param=7;   % Number of model parameters
nresult=14; % Number of models evaluated by HIPPO

datos=it.codes(:,2);
codigos=it.codes(:,1);
Akaike = zeros(nresult,1);
BIC    = zeros(nresult,1);
Fobj   = zeros(nresult,1);
FFF    = zeros(nresult,1);
CCc    = zeros(nresult,1);
MS     = zeros(nresult,1);
I15    = zeros(nresult,1);
I125   = zeros(nresult,1);
I95    = zeros(nresult,1);
I955   = zeros(nresult,1);

for i=1:1:nresult

Akaike(i)= datos{i,1}.AICc;
BIC(i)= datos{i,1}.BIC;
Fobj(i)= datos{i,1}.J_SSm;
FFF(i)=i;

%% Significance vector: How many k are not significant (CC>=2)? 
% For each model:
CCc(i)=length(find(datos{i,1}.CC>=2));

%% Sensitivity vector: How many k are not sensitive?
% For each model:
MS(i)= length(find(datos{i,1}.Ms==0));

%% Identificability analysis: 
I15(i)=length(find(datos{i,1}.Mc>0.95));
I125(i)=length(find(datos{i,1}.Mc<-0.95));
I95(i)=I15(i)+I125(i);

if I95(i)==param-sum(codigos{i,1})
    I955(i)=0; % No identifiability problems
    
    else
    I955(i)=1; % Identifiability problems

end 

end
%% Results:
%          1    2     3     4      5    6    7   
results = [FFF Fobj CCc Akaike BIC MS I955];

% Saves 'results' matrix in an Excel file
b=array2table(results,'VariableNames',{'FFF', 'Fobj','CCc','Akaike','BIC','MS','I955'});
filename = 'results_HIPPO.xlsx';
writetable(b,filename);

