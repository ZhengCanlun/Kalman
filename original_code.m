clc
clear
%======================================%
%  仿真中并没有严格的按照文章中进行编写，主要原因是绕圈过程中，在角度为0的附近会出现突变，
%为了解决这一问题，我将估计的目标状态用向量来替换角度。
%
%======================================%
t = 0:0.1:10;      %总时间，时间间隔为0.1s
target_mode = 2;    %目标运动模式，1| 匀速直线运动 2|匀速圆周运动 
switch(target_mode) 
    case 1
        x_truth = [t*3;zeros(1,101);ones(1,101)*3;zeros(1,101)];
    case 2
        r= 10;
        v = 5;
        w = v/r;
        x_truth = [10*sin(t*w);10*cos(t*w);3/10*cos(t);3/10*sin(t)];
end
%=========ctrv kf ===========%
ctrv_x_estimation = zeros(6,size(t,2));              % 估计目标状态 px，py，v，cos_theta, sin_theta w
ctrv_x_estimation(6,1) = 0.1;                   %估计状态w初值不能为0
ctrv_x_estimation(5,1) = 1;                     %且cos_theta 设为1

z  = x_truth(1:2,:) + randn(2,size(t,2))*0.1;        %输入的观测位置
ctrv_H = [1,0,0,0,0,0;0,1,0,0,0,0];             % 观测矩阵
dt = 0.4;
ctrv_Q = diag([1,1,1,0.1,0.1,0.01]);             %系统误差矩阵（需要调参）
ctrv_R = eye(2)*10;                              %测量误差矩阵（需要调参）
ctrv_P = diag([1,1,1,1,1,1]);                   %估计协防差矩阵

%==========cv kf ===========%                   %对照组
cv_x_estimation = zeros(4,2001); % px,py,v,w,d_w
cv_H = [1,0,0,0;0,1,0,0];
cv_Q = diag([0.1,0.1,5,5]);
cv_R = eye(2)*4;
cv_P = diag([1,1,1,1]);
cv_A = [1,0,dt,0;0,1,0,dt;0,0,1,0;0,0,0,1];
%==========================%



for i = 1:size(t,2)

    %====ctrv ekf====%
    v_i = ctrv_x_estimation(3,i);
    cos_i = ctrv_x_estimation(4,i);
    sin_i = ctrv_x_estimation(5,i);
    w_i = ctrv_x_estimation(6,i);
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
    
    %====show prediction ====%  %展示未来预测位置，这里展示10个步长内的轨迹
    x = ctrv_x_estimation(:,i);
    ctrv_x_prediction(:,1) = x;
    for j = 1:10
        v_j = ctrv_x_prediction(3,j);
        cos_j = ctrv_x_prediction(4,j);
        sin_j = ctrv_x_prediction(5,j);
        w_j = ctrv_x_prediction(6,j);
        inverse_j = 1;
%         if w_j <-pi
%             w_j = w_j +2*pi;
%             inverse_j = 1;
% %             ctrv_x_prediction(4,j) = w_j;
%         end
        A_13 = 1/w_j*(sin_j*(cos(w_j*dt)-1) + cos_j*sin(w_j*dt));
        A_14 = v_j/w_j*sin(w_j*dt);
        A_15 = v_j/w_j*(cos(w_j*dt)-1);
        A_16 = -v_j/w_j^2*(sin_j*(cos(w_j*dt)-1) + cos_j*sin(w_j*dt)) + v_j*dt/w_j*(-sin_j*sin(w_j*dt) + cos_j*cos(w_j*dt));
        
        A_23 = 1/w_j*(cos_j*(1 - cos(w_j*dt)) + sin_j*sin(w_j*dt));
        A_25 = v_j/w_j*sin(w_j*dt);
        A_24 = -v_j/w_j*(cos(w_j*dt)-1);
        A_26 = -v_j/w_j^2*(cos_j*(1 - cos(w_j*dt)) + sin_j*sin(w_j*dt)) + v_j*dt/w_j*(sin_j*cos(w_j*dt) + cos_j*sin(w_j*dt));
    
    
        ctrv_A_j = [1,0,A_13,A_14,A_15,A_16;...
              0,1,A_23,A_24,A_25,A_26;...
              0,0,1,0,0,0;...
              0,0,0,cos(w_j*dt),-sin(w_j*dt),- cos_j*dt*sin(w_j*dt) - sin_j*dt*cos(w_j*dt);...
              0,0,0,sin(w_j*dt),cos(w_j*dt), -sin_j*dt*sin(w_j*dt) + cos_j*dt*cos(w_j*dt);...
              0,0,0,0,0,1];
        ctrv_x_prediction(:,j+1) = ctrv_A_j*ctrv_x_prediction(:,j);
    end
    %==========================================%

    %=======ctrv ekf预测和更新========%
    x_pred = ctrv_A*x;

    ctrv_P = ctrv_A*ctrv_P*ctrv_A'+ctrv_Q;
    K = ctrv_P*ctrv_H'*(ctrv_H*ctrv_P*ctrv_H'+ctrv_R)^(-1);
    x_update = x_pred + K*(z(:,i)-ctrv_H*x_pred);
    ctrv_P = (eye(6)-K*ctrv_H)*ctrv_P;

    %==========对ctrv 的v,cos_theta,sin_theta规整化=========%
    x_update(3,1) = norm(norm(x_update(4:5,1))*x_update(3,1));
    
    x_update(4:5,1) = x_update(4:5,1)/norm(x_update(4:5,1));

    %======显示当前目标的航向角===========%
    cs = x_update(4:5,1)/norm(x_update(4:5,1));

%     if cs(2)<0
%         theta = 180 - acos(cs(1))/pi*180
%     else
%         theta = -180 + acos(cs(1))/pi*180
%     end
    %=========数据保存============%
    ctrv_x_estimation(:,i+1) = x_update;
    
    %========================================================%
    %======cv kf========%
    cv_x_prediction(:,1) = cv_x_estimation(:,i);
    for j = 1:10
        cv_x_prediction(:,j+1) = cv_A*cv_x_prediction(:,j);
    end
    %=======cv kf 预测更新=========%
    x = cv_x_estimation(:,i);
    x_pred = cv_A*x;
    cv_P = cv_A*cv_P*cv_A' + cv_Q;
    K = cv_P*cv_H'*(cv_H*cv_P*cv_H'+cv_R)^(-1);
    
    x_update = x_pred + K*(z(:,i)-cv_H*x_pred);
    cv_x_estimation(:,i+1) = x_update;
    cv_P = (eye(4)-K*cv_H)*cv_P;

    %=========visualize============%
    fig_angle = figure(1)
    hold off
    G(1) = plot(x_truth(1,1:i),x_truth(2,1:i),'r','LineWidth',2);
    hold on
    G(2) = plot(z(1,1:i),z(2,1:i),'.',"MarkerSize",15,'LineWidth',2);
    
    grid on
    axis equal
    G(7) = plot(cv_x_estimation(1,1:i),cv_x_estimation(2,1:i),'b','LineWidth',2);
    G(5) = plot(cv_x_estimation(1,i),cv_x_estimation(2,i),'bo','LineWidth',2);
    
    G(4) = plot(cv_x_prediction(1,:),cv_x_prediction(2,:),'b--','LineWidth',2);
    
    G(3) = plot(ctrv_x_estimation(1,1:i),ctrv_x_estimation(2,1:i),'g','LineWidth',2);
    G(6) = plot(ctrv_x_estimation(1,i),ctrv_x_estimation(2,i),'go','LineWidth',2);
    
    G(8) = plot(ctrv_x_prediction(1,:),ctrv_x_prediction(2,:),'g--','LineWidth',2);
    
    axis ([-20,20,-14,14])
    xtxt = xlabel('$$X$$/m','FontSize',15);
    set(xtxt,'Interpreter','latex');
    ytxt = ylabel('$$Y$$/m','Interpreter','latex','FontSize',15);
    set(ytxt,'Interpreter','latex');
    set(gca,'FontName','Times New Roman','FontSize',15);
    legend(G([1,2,5,6,7,3,4,8]),'truth trajectory','measurements','CV-KF estimated target state', ...
        'CTRV-EKF eatimated targe state', ...
        'CV-KF estimated target trajectory','CTRV-EKF estimated target trajectory', ...
        'CV-KF predicted target trajectory','CTRV-EKF predicted target trajectory','Location','northoutside','FontSize',15,"NumColumns",2);
%     pause(0.05)
    %legend%

end

% figure(1)
% hold off
% G(1) = plot(x_truth(1,:),x_truth(2,:));
% hold on
% G(2) = plot(z(1,:),z(2,:),'.');
% axis ([-inf,inf,-5,5])
% grid on
% axis equal
% G(4) = plot(cv_x_estimation(1,:),cv_x_estimation(2,:),'LineWidth',2);
% 
% G(3) = plot(ctrv_x_estimation(1,:),ctrv_x_estimation(2,:),'LineWidth',2);
% 
% legend(G,'truth','measurement','ctrv ekf','cv kf','Location','best');


% figure(2)
% subplot(2,1,1)
% hold off
% G(4) = plot(t,cv_x_estimation(1,:)- x_truth(1,:),'r','LineWidth',2);
% 
% hold on
% grid on
% axis([-inf,inf,-3,3])
% G(3) = plot(t,ctrv_x_estimation(1,:)- x_truth(1,:),'b','LineWidth',2);
% 
% subplot(2,1,2)
% hold off
% G(4) = plot(t,cv_x_estimation(2,:)- x_truth(2,:),'r','LineWidth',2);
% 
% hold on
% axis([-inf,inf,-3,3])
% grid on
% G(3) = plot(t,ctrv_x_estimation(2,:)- x_truth(2,:),'b','LineWidth',2);
