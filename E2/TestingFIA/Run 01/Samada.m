function Xdot = Samada(t,in2,U,in4)
%Samada
%    Xdot = Samada(T,IN2,U,IN4)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    08-Jun-2022 21:47:17

th1 = in4(1,:);
th2 = in4(2,:);
th3 = in4(3,:);
th4 = in4(4,:);
th5 = in4(5,:);
th6 = in4(6,:);
th7 = in4(7,:);
th8 = in4(8,:);
th9 = in4(9,:);
th10 = in4(10,:);
th11 = in4(11,:);
th12 = in4(12,:);
th13 = in4(13,:);
th14 = in4(14,:);
th15 = in4(15,:);
th16 = in4(16,:);
th17 = in4(17,:);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
x7 = in2(7,:);
t2 = th10+x2;
t3 = th11+x5;
t4 = th12+x5;
t5 = x1.^2;
t6 = 1.0./th13;
t7 = 1.0./th14;
t8 = 1.0./th15;
t9 = 1.0./x7;
t18 = x3.*3.384e-1;
t10 = t9.*th17;
t11 = 1.0./t2;
t12 = 1.0./t3;
t13 = 1.0./t4;
t14 = t5.*t6;
t15 = t5.*t7;
t16 = t5.*t8;
t22 = t18+x4;
t17 = t10-1.0;
t19 = t14+th7+x1;
t20 = t15+th8+x1;
t21 = t16+th9+x1;
t26 = 1.0./t22;
t23 = 1.0./t19;
t24 = 1.0./t20;
t25 = 1.0./t21;
t27 = t26.*x4;
t28 = t27-1.0;
Xdot = [x3.*((t17.*t23.*th1.*x1)./th4-(U.*t13.*t25.*th3.*th12.*x1)./th6+(t12.*t24.*t28.*th2.*th11.*x1)./th5);-t11.*th16.*x2.*x3;-t17.*t23.*th1.*x1.*x3;-x3.*(t17.*t23.*th1.*x1.*(3.0./5.0e+1)+t12.*t24.*t28.*th2.*th11.*x1);U.*t13.*t25.*th3.*th12.*x1.*x3;0.0;t11.*th16.*x2+t17.*t23.*th1.*x1.*x7];