% x0 = [S N Xf Xl CIT V qn];
%x0 = [0.02 0.853 8.1 0.48 0 8.26 0.075]; %Datos cultivo A
x0 = [100.95 1.4 0.792 0.047 0 3.03 0.075]; %Cultivo C
[T,Y] = ode15s(@sys,[0 100],x0);

hold on
plot(T,Y(:,1),"LineWidth",2)
plot(T,Y(:,2),"LineWidth",2)
plot(T,Y(:,3),"LineWidth",2)
plot(T,Y(:,4),"LineWidth",2)
plot(T,Y(:,5),"LineWidth",2)
plot(T,Y(:,6),"LineWidth",2)
hold off
grid on
legend("Glucose","Nitrogen","Biomass","Lipids","Citrate","Volume")