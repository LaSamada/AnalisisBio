function dt = sys(t,x)
%Model states
p=[0.755 0.039 0.314 0.519 0.272 0.89 1.116 1.482 12.22 23.08 ...
    9.75 37.63 83.64 44.7 40.11 0.716 0.045 0.047];
% estos últimos parámetros son 0 porque es una sim batch
S = x(1);
N = x(2);
Xf = x(3);
Xl = x(4);
CIT = x(5);
V = x(6);
qn = x(7);


%Model parameters
Yxs =p(4) ; Ylips = p(5); Ycits = p(6); 
gamma = 0.06; mumax = p(1); rhonmax =p(16) ; Q0 =p(17) ; Ks1 = p(7); kl1 =p(13) ; Kn =p(10) ; pilipmax = p(2); Ks2 = p(8) ; kl2 =p(14) ;
alfa = 0.36; k2 =p(11) ; qncrit = p(18); picitmax = p(3); Ks3 = p(9) ; kl3 =p(15) ; k3 = p(12);

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
ds = Sin*(fsin/V)-((mu/Yxs)+(pilip/Ylips)+(picit/Ycits))*Xf-(S/V)*Fin;
dn = Nin*(fnin/V)-rhon*Xf-(N/V)*Fin;
dxf = mu*Xf-(Xf/V)*Fin;
dxl = (pilip+gamma*mu)*Xf-(Xl/V)*Fin;
dcit = picit*Xf-(CIT/V)*Fin;
dv = Fin-Fsample;
dqn = rhon-mu*qn;



dt = [ds dn dxf dxl dcit dv dqn]';
end



