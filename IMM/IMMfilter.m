function [parameter,x,P] = IMMfilter(parameter,z,dt)

P_model = parameter.P_model;
mu_weight = parameter.mu_weight;
model = parameter.model;
c_bar = zeros(2,1);
for i = 1:2
    for j = 1:2
      c_bar(i) = c_bar(i) + P_model(i,j)*mu_weight(j);
    end
end
mu = zeros(2);
for i = 1:2
    mu(i,:)=P_model(i,:)*mu_weight(i)./c_bar';
end
mu = mu';
for i =1:2
    model(i).x_pre = model(i).x;
    model(i).P_pre = model(i).P;
end

for i = 1:2
    model(i).x = zeros(6,1);
    for j = 1:2
        model(i).x =  model(i).x + mu(i,j)*model(j).x_pre;
    end
    model(i).x(4:5) = model(i).x(4:5)/norm(model(i).x(4:5));
end
for i = 1:2
    model(i).P = zeros(6);
    for j = 1:2
        delta_x = model(i).x - model(j).x_pre;
        model(i).P =  model(i).P + mu(i,j)*(model(j).P_pre+delta_x*delta_x');
    end
end



[~,model(1).x,model(1).P,cv_x_pred] = CVKF(parameter,z,model(1).x,model(1).P,dt);

[~,model(2).x,model(2).P,ctrv_x_pred] = CTRVKF(parameter,z,model(2).x,model(2).P,dt);

model(1).v = z - parameter.cv_H*cv_x_pred;
model(2).v = z - parameter.ctrv_H*ctrv_x_pred;

model(1).S = parameter.cv_H*model(1).P*parameter.cv_H' + parameter.cv_R;
model(2).S = parameter.ctrv_H*model(2).P*parameter.ctrv_H' + parameter.ctrv_R;

Hat(1) = 1/(2^(1/2)*pi^(1/2)*norm(model(1).S)^(1/2))*exp(-1/2*model(1).v'*model(1).S^(-1)*model(1).v);
Hat(2) = 1/(2^(1/2)*pi^(1/2)*norm(model(2).S)^(1/2))*exp(-1/2*model(2).v'*model(2).S^(-1)*model(2).v);

c =sum(Hat'.*c_bar);

mu_weight = Hat'.*c_bar./c

x = model(1).x*mu_weight(1) + model(2).x*mu_weight(2);
P = (model(1).P+ (x - model(1).x)*(x - model(1).x)' )*mu_weight(1) + (model(2).P+ (x - model(2).x)*(x - model(2).x)')*mu_weight(2);

parameter.mu_weight = mu_weight;
parameter.model = model;