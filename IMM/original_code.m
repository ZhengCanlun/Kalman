clc
clear
%======================================%
%  仿真中并没有严格的按照文章中进行编写，主要原因是绕圈过程中，在角度为0的附近会出现突变，
%为了解决这一问题，我将估计的目标状态用向量来替换角度。
%
%======================================%
t = 0.1:0.1:20;      %总时间，时间间隔为0.1s
dt = 0.1;
video_save = 0;
%===========target trajectory==========%
x_truth = []
p = [0;0]
for i = 1:200
    if i<50
        v = [5,0]';
    elseif i<150
        r= 5;
        v = 10;
        w = v/r;
        v = [10*cos((i-50)/10*w),10*sin((i-50)/10*w)]';
    else
        v = [3,5]';
    end
    p = p + v*dt;
    x_truth = [x_truth,[p;v]];
end
% x_truth1 = [t1*10;zeros(1,51)+10;ones(1,51)*10;zeros(1,51)];
% r= 5;
% v = 10;
% w = v/r;
% x_truth2 = [51+r*sin(t1*w);r*cos(t1*w);10*cos(t1*w);-10*sin(t1*w)];
% x_truth3 = [t1*3;t1*5;ones(1,51)*3;ones(1,51)*5]+[x_truth2(1:2,end);0;0];

% x_truth  = [x_truth1,x_truth2,x_truth3];
%-------------measurements-----------%
z  = x_truth(1:2,:);% + randn(2,size(t,2))*0.5;   %输入的观测位置

%==========cv kf initial===========%                   %对照组2
cv_x_estimation = zeros(6,size(t,2)); % px,py,v,w,d_w
cv_x_estimation(4,1) = 1; 
cv_P = eye(6);
cv_parameter.cv_H = [1,0,0,0,0,0;0,1,0,0,0,0];
cv_parameter.cv_Q = diag([0.1,0.1,1,0.01,0.01,0]);
cv_parameter.cv_R = eye(2)*1;
%==================================%
%=========ctrv kf initial =========%          %对照组1
ctrv_x_estimation = zeros(6,size(t,2));         % 估计目标状态 px，py，v，cos_theta, sin_theta w
ctrv_x_estimation(5,1) = 1;                     %且cos_theta 设为1
ctrv_x_estimation(6,1) = 0.1;                   %估计状态w初值不能为0
ctrv_parameter.ctrv_H = [1,0,0,0,0,0;0,1,0,0,0,0];             % 观测矩阵
ctrv_parameter.ctrv_Q = diag([0.1,0.1,1,0.01,0.01,1]);             %系统误差矩阵（需要调参）
ctrv_parameter.ctrv_R = eye(2)*1.2;                              %测量误差矩阵（需要调参）
ctrv_P = diag([1,1,1,1,1,1]);                   %估计协防差矩阵

%==========IMM  ===========%                   该IMM混合模型包含两种模型
IMM_x_estimation = zeros(6,size(t,2)); % px,py,v,w,d_w
IMM_x_estimation(5,1) = 1;
IMM_x_estimation(6,1) = 0.1;
IMM_parameter.cv_H =cv_parameter.cv_H;
IMM_parameter.cv_Q = cv_parameter.cv_Q;
IMM_parameter.cv_R = cv_parameter.cv_R;

IMM_parameter.ctrv_H = ctrv_parameter.ctrv_H;
IMM_parameter.ctrv_Q = ctrv_parameter.ctrv_Q;
IMM_parameter.ctrv_R = ctrv_parameter.ctrv_R;

lambda = 0.04;
IMM_parameter.P_model = [1-lambda,lambda;lambda,1-lambda];
IMM_parameter.mu_weight = [0.7;0.3];%假设先给cv模型可信度为1
IMM_parameter.model(1).x = [0,0,0,1,0,0]';
IMM_parameter.model(1).P = eye(6);
IMM_parameter.model(2).x = [0,0,0,1,0,0.1]';
IMM_parameter.model(2).P = eye(6);

IMM_save.mu = [];
%==========================%

str_extend = 'IMM tracking';
if video_save
    str =[str_extend];
    vidObj = VideoWriter(str,'MPEG-4');
    vidObj.FrameRate = 10;
    open(vidObj);
end


for i = 1:size(t,2)
    
    %====show prediction ====%  %展示未来预测位置，这里展示10个步长内的轨迹
    %--------cv kf-----------%
    cv_x_prediction(:,1) = cv_x_estimation(:,i);
    for j = 1:10
        [cv_A_j,~,~] = CVKF(cv_parameter,z(:,i),cv_x_prediction(:,j),cv_P,dt);
        cv_x_prediction(:,j+1) = cv_A_j*cv_x_prediction(:,j);
    end
    %--------ctrv kf---------%
    x = ctrv_x_estimation(:,i);
    if i==170
        ctrv_x_estimation(6,i) = 0.1;
    end
    ctrv_x_prediction(:,1) = x;
    for j = 1:10
        [ctrv_A_j,~,~] = CTRVKF(ctrv_parameter,z(:,i),ctrv_x_prediction(:,j),ctrv_P,dt);
        ctrv_x_prediction(:,j+1) = ctrv_A_j*ctrv_x_prediction(:,j);
    end
    %-------IMM--------------%
    x = IMM_x_estimation(:,i);
    IMM_x_prediction(:,1) = x;
    mu_weight = IMM_parameter.mu_weight;
    for j = 1:10
        [cv_A_j,~,~] = CVKF(cv_parameter,z(:,i),IMM_x_prediction(:,j),cv_P,dt);
        [ctrv_A_j,~,~] = CTRVKF(ctrv_parameter,z(:,i),IMM_x_prediction(:,j),ctrv_P,dt);
        IMM_x_prediction(:,j+1) = mu_weight(1)*cv_A_j*IMM_x_prediction(:,j) ...
                                    + mu_weight(2)*ctrv_A_j*IMM_x_prediction(:,j);
    end
    %========================%


    %=======cv ctrv imm 预测和更新==========%
    [~,cv_x_estimation(:,i+1),cv_P] = CVKF(cv_parameter,z(:,i),cv_x_estimation(:,i),cv_P,dt);
    [~,ctrv_x_estimation(:,i+1),ctrv_P] = CTRVKF(ctrv_parameter,z(:,i),ctrv_x_estimation(:,i),ctrv_P,dt);
    [IMM_parameter,IMM_x_estimation(:,i+1)] = IMMfilter(IMM_parameter,z(:,i),dt);
    IMM_save.mu = [IMM_save.mu,IMM_parameter.mu_weight];
    %======================================%

    


    %=========visualize============%
    fig_angle = figure(1)
    hold off
    G(1) = plot(x_truth(1,1:i),x_truth(2,1:i),'color',[0,0,0],'LineWidth',2);
    hold on
    G(2) = plot(z(1,1:i),z(2,1:i),'.',"MarkerSize",15,'LineWidth',2);
    
    grid on
    axis equal
    plot(cv_x_estimation(1,i),cv_x_estimation(2,i),'bo','LineWidth',2);
    plot(ctrv_x_estimation(1,i),ctrv_x_estimation(2,i),'go','LineWidth',2);
    plot(IMM_x_estimation(1,i),IMM_x_estimation(2,i),'ro','LineWidth',2);

    plot(cv_x_estimation(1,1:i),cv_x_estimation(2,1:i),'b','LineWidth',2);
    plot(ctrv_x_estimation(1,1:i),ctrv_x_estimation(2,1:i),'g','LineWidth',2);
    plot(IMM_x_estimation(1,1:i),IMM_x_estimation(2,1:i),'r','LineWidth',2);
    
    
    plot(cv_x_prediction(1,:),cv_x_prediction(2,:),'b--','LineWidth',2);
    plot(ctrv_x_prediction(1,:),ctrv_x_prediction(2,:),'g--','LineWidth',2);
    plot(IMM_x_prediction(1,:),IMM_x_prediction(2,:),'r--','LineWidth',2);
    
    axis ([-5,65,-14,30])
    xtxt = xlabel('$$X$$/m','FontSize',15);
    set(xtxt,'Interpreter','latex');
    ytxt = ylabel('$$Y$$/m','Interpreter','latex','FontSize',15);
    set(ytxt,'Interpreter','latex');
    set(gca,'FontName','Times New Roman','FontSize',15);
    legend('truth trajectory','measurements','CV-KF estimated target state', ...
        'CTRV-EKF eatimated target state', ...
        'IMM estiamted target state',...
        'CV-KF estimated target trajectory','CTRV-EKF estimated target trajectory', ...
        'IMM estimated target trajectory',...
        'CV-KF predicted target trajectory','CTRV-EKF predicted target trajectory',...
        'IMM predicted target trajectory','Location','northoutside','FontSize',15,"NumColumns",2);
%     pause(0.05)
    %legend%
    if video_save
        frame = getframe(gcf);
        writeVideo(vidObj,frame);
    end
end
fig_angle = figure(1)
hold off
G(1) = plot(x_truth(1,1:i),x_truth(2,1:i),'color',[0,0,0],'LineWidth',2);
hold on
G(2) = plot(z(1,1:i),z(2,1:i),'.',"MarkerSize",15,'LineWidth',2);

grid on
axis equal
plot(cv_x_estimation(1,i),cv_x_estimation(2,i),'bo','LineWidth',2);
plot(ctrv_x_estimation(1,i),ctrv_x_estimation(2,i),'go','LineWidth',2);
plot(IMM_x_estimation(1,i),IMM_x_estimation(2,i),'ro','LineWidth',2);

plot(cv_x_estimation(1,1:i),cv_x_estimation(2,1:i),'b','LineWidth',2);
plot(ctrv_x_estimation(1,1:i),ctrv_x_estimation(2,1:i),'g','LineWidth',2);
plot(IMM_x_estimation(1,1:i),IMM_x_estimation(2,1:i),'r','LineWidth',2);


plot(cv_x_prediction(1,:),cv_x_prediction(2,:),'b--','LineWidth',2);
plot(ctrv_x_prediction(1,:),ctrv_x_prediction(2,:),'g--','LineWidth',2);
plot(IMM_x_prediction(1,:),IMM_x_prediction(2,:),'r--','LineWidth',2);

axis ([-5,65,-14,30])
xtxt = xlabel('$$X$$/m','FontSize',15);
set(xtxt,'Interpreter','latex');
ytxt = ylabel('$$Y$$/m','Interpreter','latex','FontSize',15);
set(ytxt,'Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',15);
legend('truth trajectory','measurements','CV-KF estimated target state', ...
    'CTRV-EKF eatimated target state', ...
    'IMM estiamted target state',...
    'CV-KF estimated target trajectory','CTRV-EKF estimated target trajectory', ...
    'IMM estimated target trajectory',...
    'CV-KF predicted target trajectory','CTRV-EKF predicted target trajectory',...
    'IMM predicted target trajectory','Location','northoutside','FontSize',15,"NumColumns",2);
if video_save
	close(vidObj);
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


figure(2)
for i = 1:2
subplot(2,1,i)
hold off
% error = cv_x_estimation(1:2,2:end) -  x_truth(1:2,:);
% error = sqrt(sum(error.*error,1));
G(4) = plot(t,cv_x_estimation(i,2:end) -  x_truth(i,:),'b','LineWidth',2);
hold on
G(4) = plot(t,ctrv_x_estimation(i,2:end) -  x_truth(i,:),'g','LineWidth',2);
G(4) = plot(t,IMM_x_estimation(i,2:end) -  x_truth(i,:),'r','LineWidth',2);
grid on
axis([-inf,inf,-inf,inf])
if i ==1
    ytxt = ylabel('$$p_x$$ (m)','Interpreter','latex','FontSize',15);
    set(ytxt,'Interpreter','latex');
    set(gca,'FontName','Times New Roman','FontSize',15);
    legend('cv-KF','ctrv-KF','IMM-KF','fontsize',15);
else
    xtxt = xlabel('$$t$$ (s)','FontSize',15);
    set(xtxt,'Interpreter','latex');
    ytxt = ylabel('$$p_y$$ (m)','Interpreter','latex','FontSize',15);
    set(ytxt,'Interpreter','latex');
end
end
% subplot(2,1,2)
% hold off
% error = cv_x_estimation(3,2:end).*cv_x_estimation(4:5,2:end) -  x_truth(3:4,:);
% error = sqrt(sum(error.*error,1));
% G(4) = plot(t,error,'b','LineWidth',2);
% hold on
% error = ctrv_x_estimation(3,2:end).*ctrv_x_estimation(4:5,2:end) -  x_truth(3:4,:);
% error = sqrt(sum(error.*error,1));
% G(4) = plot(t,error,'g','LineWidth',2);
% error = IMM_x_estimation(3,2:end).*IMM_x_estimation(4:5,2:end) - x_truth(3:4,:);
% error = sqrt(sum(error.*error,1));
% G(4) = plot(t,error,'r','LineWidth',2);
% grid on
% axis([-inf,inf,-inf,inf])

figure(3)
hold off
plot(t,IMM_save.mu(1,:),'LineWidth',2);
hold on
plot(t,IMM_save.mu(2,:),'LineWidth',2);
grid on
xtxt = xlabel('$$t$$ (s)','FontSize',15);
set(xtxt,'Interpreter','latex');
ytxt = ylabel('weight value','Interpreter','latex','FontSize',15);
set(ytxt,'Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',15);
legend('IMM-cv','IMM-ctrv','fontsize',15);