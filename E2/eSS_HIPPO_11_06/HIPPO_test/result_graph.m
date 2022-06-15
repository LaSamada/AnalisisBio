load datos.txt -ascii;
[r1,c1]  = size(datos); 
texp     = datos(:,1);
yexp     = datos(:,2:c1);
X1 = [63 0.853 0.15 0 0 0.5 0.075];
[tout,yout]     = ode23(@model_final,[0 45.5],X1);

figure
plot(tout,yout)
grid on
xlabel('Time (h)')
hold on
plot(texp,yexp,'o')
legend('Glucose','Nitrogen','Functional biomass','Lipid biomass','Citric acid','Volume'... 
    ,'qn','Experimental glucose','Experimental functional biomass')
