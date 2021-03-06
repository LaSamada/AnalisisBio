% x0 = [S N Xf Xl CIT V qn];
%x0 = [0.02 0.853 8.1 0.48 0 8.26 0.075]; %Datos cultivo A
%x0 = [100.95 1.4 0.792 0.047 0 3.03 0.075]; %Cultivo C
% x0 = [98.34 0.885 0.47 0.028 0 4.02 0.075] ;
 
x01 = [60 1 0.45 0.05 0 0.5 0.075];

x0 = [98.34 0.885 0.47 0.028 0 4.02 0.075]; % datos batch sebastian kinda

data  = load('datos.txt'); %edit here
[~,n] = size(data); 
texp  = data(:,1)';
ydata = data(:,2:n);

% data1  = load('datos_agosin.txt'); %edit here
% [~,n] = size(data); 
% texp  = data(:,1)';
% ydata = data(:,2:n);

[t3,y3] = ode23s(@sys_hippo,[0 35],x0);


plot(t3,y3)
hold on
grid on
plot(texp,ydata,'o')
legend('Glucose','Nitrogen','Biomass','Lipids','Citrate','Volume','qn')
hold off
