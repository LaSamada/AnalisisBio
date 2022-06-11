function dxdt   = meta(ModelName,t,x,theta,dim)

%This is an auxiliary function for simulation using the provided Zenteno
% model in symbolical format.

%Dont change these
nx  =   dim(1); 
nth =   dim(2);

%Model states
X   = x(1:nx);
%This conditional loop is to correctly use model inputs, where the main
% objective is to indicate if yeast is or not on thermal death.

    %First we interpolate operational conditions with current integration
    % time "t"
    %The following equation is the thermal death yeast loss
    % (Zenteno et al., 2010)
    %The following conditional communicate the model if yeast death
    % should be accounted
if X(7) >= 0.047
    U = 1;
else
    U = 0;
end
    

%Construction of d(x_theta)/dt from eq. 3 (Stigter & Molenaar, 2015)
% (Dont change unless inputs (Ui) must be modified)
f=eval(['@' ModelName]);
xdot=f(t,X,U,theta);

dxdth=reshape(x((nx+1):end),nx,nth);

dfdxModel=eval(['@dfdx' ModelName]);
dfdx=dfdxModel(t,X,U,theta);

dfdthModel=eval(['@dfdth' ModelName]);
dfdth=dfdthModel(t,X,U,theta);

dxdthdot=dfdx*dxdth+dfdth;

dxdt=[xdot; dxdthdot(:)];