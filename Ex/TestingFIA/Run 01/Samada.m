function Xdot = Samada(t,in2,U,in4)
%Samada
%    Xdot = Samada(T,IN2,U,IN4)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    03-Jul-2022 14:48:44

th2 = in4(1,:);
th3 = in4(2,:);
th4 = in4(3,:);
th8 = in4(4,:);
th11 = in4(5,:);
th16 = in4(6,:);
th17 = in4(7,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
x7 = in2(7,:);
t2 = th11+x5;
t3 = x1.^2;
t4 = 1.0./x7;
t6 = x2+5.77e+2./2.5e+1;
t8 = x3.*3.384e-1;
t9 = x5+3.763e+1;
t5 = 1.0./t2;
t7 = t3.*(1.0e+1./4.47e+2);
t10 = t3.*1.195600191296031e-2;
t11 = t3.*2.493143854400399e-2;
t12 = 1.0./t6;
t13 = t8+x4;
t14 = 1.0./t9;
t16 = t4.*th17.*(1.51e+2./2.0e+2);
t15 = t7+th8+x1;
t17 = 1.0./t13;
t20 = t10+x1+2.79e+2./2.5e+2;
t22 = t11+x1+6.11e+2./5.0e+1;
t23 = t16-1.51e+2./2.0e+2;
t18 = 1.0./t15;
t19 = t17.*x4;
t24 = 1.0./t20;
t25 = 1.0./t22;
t21 = t19-1.0;
Xdot = [x3.*((t23.*t24.*x1)./th4-U.*t14.*t25.*th3.*x1.*4.228089887640449e+1+t5.*t18.*t21.*th2.*th11.*x1.*(1.25e+2./3.4e+1));-t12.*th16.*x2.*x3;-t23.*t24.*x1.*x3;-x3.*(t24.*x1.*(t4.*th17.*4.53e-2-4.53e-2)+t5.*t18.*t21.*th2.*th11.*x1);U.*t14.*t25.*th3.*x1.*x3.*3.763e+1;0.0;t12.*th16.*x2+t23.*t24.*x1.*x7];
