function dt = sys(~,x,p)
%Model states
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
alfa = 0.36; k2 =p(11) ; qncrit = p(18); picitmax = p(3); Ks3 =p(9) ; kl3 =p(15) ; k3 = p(12);

Sin =  p(19); fsin =p(20) ;  Nin =p(21) ; fnin=p(22) ; Fsample = 0;

Fin = fsin + fnin;

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



