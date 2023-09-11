function [ctrv_A,x,ctrv_P,x_pred] = CTRVKF(ctrv_parameter,z,x,ctrv_P,dt)

ctrv_H = ctrv_parameter.ctrv_H;
ctrv_Q = ctrv_parameter.ctrv_Q;
ctrv_R = ctrv_parameter.ctrv_R;

%====ctrv ekf====%
v_i = x(3);
cos_i = x(4);
sin_i = x(5);
w_i = x(6);
%=======计算状态矩阵===========%
A_13 = 1/w_i*(sin_i*(cos(w_i*dt)-1) + cos_i*sin(w_i*dt));
A_14 = v_i/w_i*sin(w_i*dt);
A_15 = v_i/w_i*(cos(w_i*dt)-1);
A_16 = -v_i/w_i^2*(sin_i*(cos(w_i*dt)-1) + cos_i*sin(w_i*dt))...
+ v_i*dt/w_i*(-sin_i*sin(w_i*dt) + cos_i*cos(w_i*dt));

A_23 = 1/w_i*(cos_i*(1 - cos(w_i*dt)) + sin_i*sin(w_i*dt));
A_25 = v_i/w_i*sin(w_i*dt);
A_24 = -v_i/w_i*(cos(w_i*dt)-1);
A_26 = -v_i/w_i^2*(cos_i*(1 - cos(w_i*dt)) + sin_i*sin(w_i*dt)) ...
    + v_i*dt/w_i*(sin_i*cos(w_i*dt) + cos_i*sin(w_i*dt));

ctrv_A = [1,0,A_13,A_14,A_15,A_16;...
          0,1,A_23,A_24,A_25,A_26;...
          0,0,1,0,0,0;...
          0,0,0,cos(w_i*dt),-sin(w_i*dt),- cos_i*dt*sin(w_i*dt) - sin_i*dt*cos(w_i*dt);...
          0,0,0,sin(w_i*dt),cos(w_i*dt), -sin_i*dt*sin(w_i*dt) + cos_i*dt*cos(w_i*dt);...
          0,0,0,0,0,1];


x_pred = ctrv_A*x;

ctrv_P = ctrv_A*ctrv_P*ctrv_A'+ctrv_Q;
K = ctrv_P*ctrv_H'*(ctrv_H*ctrv_P*ctrv_H'+ctrv_R)^(-1);
x = x_pred + K*(z-ctrv_H*x_pred);
ctrv_P = (eye(6)-K*ctrv_H)*ctrv_P;

%==========对ctrv 的v,cos_theta,sin_theta规整化=========%
x(3) = norm(norm(x(4:5))*x(3));

x(4:5) = x(4:5)/norm(x(4:5));
