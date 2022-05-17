function dt = sys(~,x)
%Model states
S = x(1);
N = x(2);
Xf = x(3);
Xl = x(4);
CIT = x(5);
V = x(6);
qn = x(7);


%Model parameters
Sin = ; fsin = ; Yxs = ; Ylips = ; Ycits = ; Fin = ; Nin = ; fnin= ; rhon = ; Fsample = ;
gamma = ; mumax = ; rhonmax = ; Q0 = ; Ks1 = ; kl1 = ; Kn = ; pilipmax = ; Ks2 = ; kl2 = ;
alfa = ; k2 = ; qncrit = ; picitmax = ; Ks3 = ; kl3 = ; k3 = ;

%Algebraic equations
Xt = Xf*(1-gamma)+Xl;
mu = mumax*(1-(Q0/qn))*(S/(Ks1+S+(S^2)/kl1));
rhon = rhonmax*(N/(Kn+N));
pilip = pilipmax*(S/(Ks2+S+(S^2)/kl2))*(1-Xl/(alfa*Xt))*(k2/(k2+CIT));

if qn>qncrit
    indq = 0;
else
    indq = 1;
end 
   
picit = picitmax*indq*(S/(Ks3+S+(S^2/kl3)))*(k3/(k3+CIT));

%Differential equations
ds = Sin*(fsin/V)-(mu/Yxs+pilip/Ylips+picit/Ycits)*Xf-(S/V)*Fin;
dn = Nin*(fnin/V)-rhon*Xf-(N/V)*Fin;
dxf = mu*Xf-(Xf/V)*Fin;
dxl = (pilip+gamma*mu)*Xf-(Xl/V)*Fin;
dcit = picit*Xf-(CIT/V)*Fin;
dv = Fin-Fsample;
dqn = rhon-mu*qn;







dt = [ds dn dxf dxl dcit dv dqn]';
end


