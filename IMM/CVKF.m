function [cv_A,x,cv_P,x_pred] = CVKF(cv_parameter,z,x,cv_P,dt)
cv_Q = cv_parameter.cv_Q;
cv_R = cv_parameter.cv_R;
cv_H =cv_parameter.cv_H ;
%=======计算状态矩阵===========%
v = x(3);
cos_i = x(4);
sin_i = x(5);

cv_A = [0,0,cos_i,v,0,0;...
          0,0,sin_i,0,v,0;...
          0,0,0,0,0,0;...
          0,0,0,0,0,0;...
          0,0,0,0,0,0;
          0,0,0,0,0,0];
cv_A = eye(6) + cv_A*dt;
x_pred = cv_A*x;
cv_P = cv_A*cv_P*cv_A' + cv_Q;
K = cv_P*cv_H'*(cv_H*cv_P*cv_H'+cv_R)^(-1);

x = x_pred + K*(z-cv_H*x_pred);
cv_P = (eye(6)-K*cv_H)*cv_P;
if norm(x(4:5))
x(4:5) = x(4:5)/norm(x(4:5));
x(3) = x(3)*norm(x(4:5));
end