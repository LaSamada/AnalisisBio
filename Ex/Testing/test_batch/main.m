% x0 = [S N Xf Xl CIT V qn];
%x0 = [0.02 0.853 8.1 0.48 0 8.26 0.075]; %Datos cultivo A
%x0 = [100.95 1.4 0.792 0.047 0 3.03 0.075]; %Cultivo C
x0 = [63 0.853 0.15 0 0 0.5 0.075]; % datos batch sebastian kinda

data  = load('datos.txt'); %edit here
[~,n] = size(data); 
texp  = data(:,1)';
ydata = data(:,2:n);

[t1,y1] = ode23s(@sys,[0 45.5],x0);
[t2,y2] = ode23s(@sys_ess,[0 45.5],x0);
[t3,y3] = ode23s(@sys_hippo,[0 45.5],x0);

a = tiledlayout(1,3);
nexttile
plot(t1,y1)
hold on
plot(texp,ydata,'o')
grid on
legend('Glucose','Nitrogen','Biomass','Lipids','Citrate','Volume','qn')
hold off
nexttile
plot(t2,y2)
hold on
grid on
plot(texp,ydata,'o')
legend('Glucose','Nitrogen','Biomass','Lipids','Citrate','Volume','qn')
hold off
nexttile
plot(t3,y3)
hold on
grid on
plot(texp,ydata,'o')
legend('Glucose','Nitrogen','Biomass','Lipids','Citrate','Volume','qn')
hold off
title(a,'Parameter estimation: paper (left), eSS (middle) and HIPPO (right)')



