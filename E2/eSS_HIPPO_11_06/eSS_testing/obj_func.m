function [J,g,R]=obj_func(x,texp,yexp)

% Integration of ODEs
[~,yout0] = ode15s(@model_ess,texp,[63 0.853 0.15 0.05 0 0.5 0.075],[],x);
yout = yout0;
yout(:,2) = []; % active biomass is not measured
yout(:,3) = [];
yout(:,3) = [];
yout(:,3) = [];
yout(:,3) = [];% intermediary is not measured (originally the 3rd column)

R=(yout-yexp);
% R(:,1)=100*R(:,1); % urea weight
R(:,2)=5*R(:,2);  % measured biomass weight
% R(:,5)=100*R(:,5); % Ga3 weight
J = sum(sum(R.^2));

R=reshape(R,numel(R),1);
g=0;