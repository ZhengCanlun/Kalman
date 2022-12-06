clc
clear

t = 0:0.1:200;
step = size(t,2);
target_mode = 2;
switch(target_mode)
    case 1
        x_truth = [t*3;zeros(1,step);onse(1,step)*10];
    case 2
        x_truth = [10*sin(t*1);10*cos(t*1);10*cos(t*1)];
    case 3
        x_truth = [10*sin(t*1);10*cos(t*1);10*onse(1,step)];
end
%=========ctrv kf ===========%
ctrv_x_estimation = zeros(6,step); % px,py,p_z,v,w,d_w
ctrv_x_estimation(6,1) = 0.1;
ctrv_x_estimation(5,1) = 1;

z  = x_truth(1:2,:) + randn(2,2001)*0.5;
ctrv_H = [1,0,0,0,0,0;0,1,0,0,0,0];
dt = 0.1;
sigma_a = 1;
sigma_w = 1;
ctrv_Q = diag([1,1,1,0.1,0.1,0.1]);
ctrv_R = eye(2)*4;
ctrv_P = diag([1,1,1,1,1,1,1,1]);

%==========cv kf ===========%
cv_x_estimation = zeros(4,2001); % px,py,v,w,d_w
cv_H = [1,0,0,0;0,1,0,0];
cv_Q = diag([0.1,0.1,5,5]);
cv_R = eye(2)*4;
cv_P = diag([1,1,1,1]);
cv_A = [1,0,dt,0;0,1,0,dt;0,0,1,0;0,0,0,1];




for i = 1:2000

    %====ctrv ekf====%
    v_i = ctrv_x_estimation(3,i);
    cos_i = ctrv_x_estimation(4,i);
    sin_i = ctrv_x_estimation(5,i);
    w_i = ctrv_x_estimation(6,i);
    
    x_i_1 =ctrv_x_estimation(:,i);
    inverse = 1;

%         x_i_pred = x_i_1 + [v_i/d_w_i*(sin(w_i + d_w_i*dt) - sin(w_i));...
%          v_i/d_w_i*(-cos(w_i + d_w_i*dt) + cos(w_i));0;d_w_i*dt;0];

    A_13 = 1/w_i*(sin_i*(cos(w_i*dt)-1) + cos_i*sin(w_i*dt));
    A_14 = v_i/w_i*sin(w_i*dt);
    A_15 = v_i/w_i*(cos(w_i*dt)-1);
    A_16 = -v_i/w_i^2*(sin_i*(cos(w_i*dt)-1) + cos_i*sin(w_i*dt)) + v_i*dt/w_i*(-sin_i*sin(w_i*dt) + cos_i*cos(w_i*dt));
    
    A_23 = 1/w_i*(cos_i*(1 - cos(w_i*dt)) + sin_i*sin(w_i*dt));
    A_25 = v_i/w_i*sin(w_i*dt);
    A_24 = -v_i/w_i*(cos(w_i*dt)-1);
    A_26 = -v_i/w_i^2*(cos_i*(1 - cos(w_i*dt)) + sin_i*sin(w_i*dt)) + v_i*dt/w_i*(sin_i*cos(w_i*dt) + cos_i*sin(w_i*dt));
    
    
    

    ctrv_A = [1,0,A_13,A_14,A_15,A_16;...
              0,1,A_23,A_24,A_25,A_26;...
              0,0,1,0,0,0;...
              0,0,0,cos(w_i*dt),-sin(w_i*dt),- cos_i*dt*sin(w_i*dt) - sin_i*dt*cos(w_i*dt);...
              0,0,0,sin(w_i*dt),cos(w_i*dt), -sin_i*dt*sin(w_i*dt) + cos_i*dt*cos(w_i*dt);...
              0,0,0,0,0,1];
    
%     end
%     ctrv_Q = [(0.5*dt^2*sigma_a*cos(w_i))^2,1/4*dt^4*sigma_a^2*sin(w_i)*cos(w_i),1/2*dt^3*sigma_a^2*cos(w_i),0,0;...
%         1/4*dt^4*sigma_a^2*sin(w_i)*cos(w_i),(0.5*dt^2*sigma_a*sin(w_i))^2,1/2*dt^3*sigma_a^2*sin(w_i),0,0;...
%         1/2*dt^3*sigma_a*cos(w_i),1/2*dt^3*sigma_a^2*sin(w_i),dt^2*sigma_a^2,0,0;...
%         0,0,0,(1/2*dt^2*sigma_w)^2,1/2*dt^3*sigma_w^2;...
%         0,0,0,1/2*dt^3*sigma_w^2,dt^2*sigma_w];
    if i==44
        i
    end
    %====show prediction ====%
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
%         ctrv_x_prediction(4,j) = rem(ctrv_x_prediction(4,j),2*pi);
        ctrv_x_prediction(:,j+1) = ctrv_A_j*ctrv_x_prediction(:,j);
%         ctrv_x_prediction(4:5,j+1) = ctrv_x_prediction(4:5,j+1)/norm(ctrv_x_prediction(4:5,j+1));
%         if inverse_j==-1
%             ctrv_x_prediction(2,j+1) = -ctrv_x_prediction(2,j+1);
%         end
    end

   
    x_pred = ctrv_A*x;

    ctrv_P = ctrv_A*ctrv_P*ctrv_A'+ctrv_Q;
    K = ctrv_P*ctrv_H'*(ctrv_H*ctrv_P*ctrv_H'+ctrv_R)^(-1);
    x_update = x_pred + K*(z(:,i)-ctrv_H*x_pred);

%     if (x_update(4)<-pi)
%         x_update(4) = x_update(4)+2*pi;
%     end
%     if inverse ==-1
%         x_update(4) = x_update(4) - 2*pi;
%     end
%     x_update(4:5,1) = x_update(4:5,1)/norm(x_update(4:5,1));
    cs = x_update(4:5,1)/norm(x_update(4:5,1))
    if cs(2)<0
        theta = 180 - acos(cs(1))/pi*180
    else
        theta = -180 + acos(cs(1))/pi*180
    end

    ctrv_x_estimation(:,i+1) = x_update;
    
    ctrv_P = (eye(6)-K*ctrv_H)*ctrv_P;

    %======cv kf========%
    cv_x_prediction(:,1) = cv_x_estimation(:,i);
    for j = 1:10
        cv_x_prediction(:,j+1) = cv_A*cv_x_prediction(:,j);
    end

    x = cv_x_estimation(:,i);
    x_pred = cv_A*x;
    cv_P = cv_A*cv_P*cv_A' + cv_Q;
    K = cv_P*cv_H'*(cv_H*cv_P*cv_H'+cv_R)^(-1);
    
    x_update = x_pred + K*(z(:,i)-cv_H*x_pred);
    cv_x_estimation(:,i+1) = x_update;
    cv_P = (eye(4)-K*cv_H)*cv_P;
    figure(1)
    hold off
    G(1) = plot(x_truth(1,1:i),x_truth(2,1:i));
    hold on
    G(2) = plot(z(1,1:i),z(2,1:i),'.');
    axis ([-inf,inf,-11,11])
    grid on
    axis equal
    G(4) = plot(cv_x_estimation(1,1:i),cv_x_estimation(2,1:i),'b','LineWidth',1);
    G(4) = plot(cv_x_estimation(1,i),cv_x_estimation(2,i),'bo','LineWidth',2);
    
    G(4) = plot(cv_x_prediction(1,:),cv_x_prediction(2,:),'b','LineWidth',2);
    
    G(3) = plot(ctrv_x_estimation(1,1:i),ctrv_x_estimation(2,1:i),'r','LineWidth',1);
    G(3) = plot(ctrv_x_estimation(1,i),ctrv_x_estimation(2,i),'ro','LineWidth',2);
    
    G(3) = plot(ctrv_x_prediction(1,:),ctrv_x_prediction(2,:),'r','LineWidth',2);
    pause(0.05)
end

figure(1)
hold off
G(1) = plot(x_truth(1,:),x_truth(2,:));
hold on
G(2) = plot(z(1,:),z(2,:),'.');
axis ([-inf,inf,-5,5])
grid on
axis equal
G(4) = plot(cv_x_estimation(1,:),cv_x_estimation(2,:),'LineWidth',2);

G(3) = plot(ctrv_x_estimation(1,:),ctrv_x_estimation(2,:),'LineWidth',2);

legend(G,'truth','measurement','ctrv ekf','cv kf','Location','best');


figure(2)
subplot(2,1,1)
hold off
G(4) = plot(t,cv_x_estimation(1,:)- x_truth(1,:),'r','LineWidth',2);

hold on
grid on
axis([-inf,inf,-3,3])
G(3) = plot(t,ctrv_x_estimation(1,:)- x_truth(1,:),'b','LineWidth',2);

subplot(2,1,2)
hold off
G(4) = plot(t,cv_x_estimation(2,:)- x_truth(2,:),'r','LineWidth',2);

hold on
axis([-inf,inf,-3,3])
grid on
G(3) = plot(t,ctrv_x_estimation(2,:)- x_truth(2,:),'b','LineWidth',2);

% subplot(1,2,2)
% % hold off
% % G(1) = plot(x_truth(1,:),x_truth(2,:));
% % hold on
% % G(2) = plot(z(1,:),z(2,:),'.');
% % axis equal
% grid on
% G(4) = plot(cv_x_estimation(1,:) - x_truth(1,:),cv_x_estimation(2,:)- x_truth(2,:),'b','LineWidth',2);

% legend(G,'truth','measurement','ctrv ekf','cv kf','Location','best');