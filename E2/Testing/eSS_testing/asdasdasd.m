x = [ 0.755  0.039  0.314  0.519  0.272  0.89  1.116  1.482  12.22  23.08  9.75  37.63  83.64  44.7  40.11  0.716  0.045  0.047];
[~,yout0] = ode15s(@model_ess,texp,[63 0.853 0.15 0.05 0 0.5 0.075],[],x);
yout = yout0;
yout(:,2) = []; % active biomass is not measured
yout(:,3) = [];
yout(:,3) = [];
yout(:,3) = [];
yout(:,3) = [];% intermediary is not measured (originally the 3rd column)

R = yout-yexp


