clc
clear

t = 0:0.1:20;
target_mode = 2;
switch(target_mode)
    case 1
        x_truth = [t*3;zeros(1,2001);ones(1,2001)*3;zeros(1,2001)];
    case 2
        x_truth = [-10*sin(t);10*cos(t);3/10*cos(t);3/10*sin(t)];
    case 3
        x_truth = [9*sin(t/10);9*cos(2*t);3/10*cos(t);3/10*sin(t)];
end
%=========ctrv kf ===========%
ctrv_x_estimation = zeros(5,2001); % px,py,v,w,d_w
ctrv_x_estimation(5,1) = 0.1;
z  = x_truth(1:2,:);% + randn(2,2001)*0.5;
ctrv_H = [1,0,0,0,0;0,1,0,0,0];
dt = 0.1;
sigma_a = 10;
sigma_w = 10;
ctrv_Q = diag([0.1,0.1,1,0.1,0.5]);
ctrv_R = eye(2)*1;
ctrv_P = diag([1,1,1,1,1]);

%==========cv kf ===========%
cv_x_estimation = zeros(4,2001); % px,py,v,w,d_w
cv_H = [1,0,0,0;0,1,0,0];
cv_Q = diag([0.1,0.1,5,5]);
cv_R = eye(2)*1;
cv_P = diag([1,1,1,1]);
cv_A = [1,0,dt,0;0,1,0,dt;0,0,1,0;0,0,0,1];




for i = 1:size(t,2)

    %====ctrv ekf====%
    v_i = ctrv_x_estimation(3,i);
    w_i = ctrv_x_estimation(4,i);
    d_w_i = ctrv_x_estimation(5,i);
    x_i_1 =ctrv_x_estimation(:,i);
    inverse = 1;

    if d_w_i==0
        ctrv_A = [1,0,cos(w_i),-v_i*sin(w_i),0;...
            0,1,sin(w_i),v_i*cos(w_i),0;...
            0,0,1,0,0;
            0,0,0,1,dt;
            0,0,0,0,1];        
    else
%         x_i_pred = x_i_1 + [v_i/d_w_i*(sin(w_i + d_w_i*dt) - sin(w_i));...
%          v_i/d_w_i*(-cos(w_i + d_w_i*dt) + cos(w_i));0;d_w_i*dt;0];
    


    ctrv_A = [1,0,1/d_w_i*(sin(w_i+d_w_i*dt)-sin(w_i)),inverse*v_i/d_w_i*(cos(w_i + d_w_i*dt) - cos(w_i)),inverse*dt*v_i/d_w_i*cos(dt*d_w_i+w_i) - inverse*v_i/d_w_i^2*(-sin(w_i)+sin(dt*d_w_i+w_i)) ;...
        0,1,1/d_w_i*(-cos(w_i+d_w_i*dt)+cos(w_i)),inverse*v_i/d_w_i*(-sin(w_i)+sin(w_i+dt*d_w_i)),inverse*dt*v_i/d_w_i*sin(dt*d_w_i+w_i) - inverse*v_i/d_w_i^2*(cos(w_i)-cos(dt*d_w_i+w_i));...
        0,0,1,0,0;...
        0,0,0,1,dt;...
        0,0,0,0,1];
    
    end
    ctrv_Q = [(0.5*dt^2*sigma_a*cos(w_i))^2,1/4*dt^4*sigma_a^2*sin(w_i)*cos(w_i),1/2*dt^3*sigma_a^2*cos(w_i),0,0;...
        1/4*dt^4*sigma_a^2*sin(w_i)*cos(w_i),(0.5*dt^2*sigma_a*sin(w_i))^2,1/2*dt^3*sigma_a^2*sin(w_i),0,0;...
        1/2*dt^3*sigma_a*cos(w_i),1/2*dt^3*sigma_a^2*sin(w_i),dt^2*sigma_a^2,0,0;...
        0,0,0,(1/2*dt^2*sigma_w)^2,1/2*dt^3*sigma_w^2;...
        0,0,0,1/2*dt^3*sigma_w^2,dt^2*sigma_w];
    if i==30
        i
    end
    %====show prediction ====%
    x = ctrv_x_estimation(:,i);
    ctrv_x_prediction(:,1) = x;
    for j = 1:20
        v_j = ctrv_x_prediction(3,j);
        w_j = ctrv_x_prediction(4,j);
        inverse_j = 1;
        if w_j <-pi
            w_j = w_j +2*pi;
            inverse_j = 1;
%             ctrv_x_prediction(4,j) = w_j;
        end
        d_w_j = ctrv_x_prediction(5,j);
        ctrv_A_j = [1,0,(sin( w_j+d_w_j*dt)-sin(w_j))/d_w_j,v_j/d_w_j*(cos(w_j + d_w_j*dt) - cos(w_j)),dt*v_j/d_w_j*cos(dt*d_w_j+w_j) - v_j/d_w_j^2*(-sin(w_j)+sin(dt*d_w_j+w_j)) ;...
        0,inverse_j,1/d_w_j*(-cos(w_j+d_w_j*dt)+cos(w_j)),v_j/d_w_j*(-sin(w_j)+sin(w_j+dt*d_w_j)),dt*v_j/d_w_j*sin(dt*d_w_j+w_j) - v_j/d_w_j^2*(cos(w_j)-cos(dt*d_w_j+w_j));...
        0,0,1,0,0;...
        0,0,0,1,dt;...
        0,0,0,0,1];
        ctrv_x_prediction(4,j) = rem(ctrv_x_prediction(4,j),2*pi);
        ctrv_x_prediction(:,j+1) = ctrv_A_j*ctrv_x_prediction(:,j);
        if inverse_j==-1
            ctrv_x_prediction(2,j+1) = -ctrv_x_prediction(2,j+1);
        end
    end

   
    x_pred = ctrv_A*x;
    if (x_pred(4)<-pi)
        x_pred(4) = x_pred(4)+2*pi;
    end
    if (x_pred(4)>pi)
        x_pred(4) = x_pred(4)-2*pi;
    end
    ctrv_P = ctrv_A*ctrv_P*ctrv_A'+ctrv_Q;
    K = ctrv_P*ctrv_H'*(ctrv_H*ctrv_P*ctrv_H'+ctrv_R)^(-1);
    x_update = x_pred + K*(z(:,i)-ctrv_H*x_pred);
    x_update(4)*180/pi
    x_update(5)*180/pi
    
    if (x_update(4)<-pi)
        x_update(4) = x_update(4)+2*pi;
    end
    if (x_update(4)>pi)
        x_update(4) = x_update(4)-2*pi;
    end
    
%     if inverse ==-1
%         x_update(4) = x_update(4) - 2*pi;
%     end
    ctrv_x_estimation(:,i+1) = x_update;
    ctrv_P = (eye(5)-K*ctrv_H)*ctrv_P;

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