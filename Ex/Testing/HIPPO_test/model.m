function dx = model(~,x,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DO NOT MODIFY THIS SECTION

kfixed = evalin('base','kfixed');
p = zeros(length(kfixed));

j = 1;
for i = 1:length(kfixed)
    if isnan(kfixed(i))
        p(i) = k(j);
        j    = j+1;
    else
        p(i) = kfixed(i);
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Model states
%p = [0.755 0.039 0.314 0.519 0.272 0.89 1.116 1.482 12.22 23.08 ...
%    9.75 37.63 83.64 44.7 40.11 0.716 0.045 0.047];

%[th1 th3 th4 th7 th10 th16 th17]

S = x(1);
N = x(2);
Xf = x(3);
Xl = x(4);
CIT = x(5);
V = x(6);
qn = x(7);


%[th2 th3 th4 th8 th11 th16 th17]

%[th2 th3 th4 th5 th12 th16 th17]

%p = [0.755 0.039 0.314 0.519 0.272 0.89 1.116 1.482 12.22 23.08 ...
%    9.75 37.63 83.64 44.7 40.11 0.716 0.045 0.047];


%Model parameters
Yxs =p(3) ; Ylips = 1.482 ; Ycits = 0.89; 
gamma = 0.06; mumax = 0.755; rhonmax =p(6) ; Q0 =p(7) ; Ks1 = 1.116; kl1 = 83.64 ; Kn = 23.08 ; pilipmax = p(1); Ks2 =  p(4); kl2 = 44.7 ;
alfa = 0.36; k2 = p(5) ; qncrit = 0.047; picitmax = p(2); Ks3 = 12.22 ; kl3 = 40.11 ; k3 = 37.63;

Sin =  0; fsin = 0 ;  Nin = 0 ; fnin = 0 ; Fsample = 0;

Fin = fsin + fnin;

%Algebraic equations
Xt = Xf*(1-gamma)+Xl;
mu = mumax*(1-(Q0/qn))*(S/(Ks1+S+((S^2)/kl1)));
rhon = rhonmax*(N/(Kn+N));
pilip = pilipmax*(S/(Ks2+S+((S^2)/kl2)))*(1-(1/alfa)*(Xl/Xt))*(k2/(k2+CIT));


if qn>qncrit
    indq = 0;
else
    indq = 1;
end 
   
picit = picitmax*indq*(S/(Ks3+S+((S^2)/kl3)))*(k3/(k3+CIT));

%Differential equations
if S < 0
    ds = 0;
else 
    ds = Sin*(fsin/V)-((mu/Yxs)+(pilip/Ylips)+(picit/Ycits))*Xf-(S/V)*Fin;
end 
if N < 0
    dn = 0;
else 
    dn = Nin*(fnin/V)-rhon*Xf-(N/V)*Fin;
end 
if Xf < 0
    dxf = 0;
else 
    dxf = mu*Xf-(Xf/V)*Fin;
end 
if Xl < 0
    dxl = 0;
else 
    dxl = (pilip+gamma*mu)*Xf-(Xl/V)*Fin;
end 
if CIT < 0
    dcit = 0;
else 
    dcit = picit*Xf-(CIT/V)*Fin;
end 
if V < 0
    dv = 0;
else
dv = Fin-Fsample;
end
if qn < 0
    dqn = 0;
else 
    dqn = rhon-mu*qn;
end



dx = [ds dn dxf dxl dcit dv dqn]';
end



