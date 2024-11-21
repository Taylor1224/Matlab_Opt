function y_pred = RMLS_Fun_Nf1(x_pred1)
    K_ =1420 ;n_ = 0.05;E = 161000;
    h = 1002; 
    N = 100;
    data = readmatrix('BLISK4.csv');
    temp = randperm(N);    
    X1 = data(temp,1:4)';  
    Y1 = data(temp,5)';    
    Y2 = data(temp,7)';    
    [n, d] = size(X1);
    kernel = @(r, h) exp(- (r.^2) / (2 * h^2));
    beta = 0.3;
    
    % 第一级计算叶片应力stress1和应变strain1:
    distances = sqrt(sum((X1 - x_pred1).^2, 1));
    W1_1 = diag(kernel(distances, h));
    residual = abs(distances);
    W1_1 = W1_1./(1+beta*residual);
    A1_1 = [ones(1,d); X1];  
    theta1_1 = (A1_1 * W1_1 * A1_1') \ (A1_1 * W1_1 * Y1');
    stress1 = [1, x_pred1'] * theta1_1;
    strain1 = stress1/E+(stress1/K_)^(1/n_);
    
    % 第一级计算盘的应力stress2和应变strain2：
    distances = sqrt(sum((X1 - x_pred1).^2, 1));
    W1_2 = diag(kernel(distances, h));
    residual = abs(distances);
    W1_2 = W1_2./(1+beta*residual);
    A1_2 = [ones(1,d); X1];
    theta1_2 = (A1_2 * W1_2 * A1_2') \ (A1_2 * W1_2 * Y2');
    stress2 = [1, x_pred1'] * theta1_2;
    strain2 = stress2/E+(stress2/K_)^(1/n_);
    
    % 第二级计算叶片寿命Nf1:
    muX = 161000; sigmaX = 3220;
    y_optim = lhsdesign(N,1);
    x_optim = norminv(y_optim,muX,sigmaX);

    x_strain1 = ones(N,1);
    x_stress1 = Y1';
    for i = 1:N
        xstrain1 = x_stress1(i,1)/x_optim(i,1)+(x_stress1(i,1)/K_)^(1/n_);
        x_strain1(i,1) = xstrain1;
    end
    e = load('e.mat').e; m1 = x_strain1; n1 = x_stress1;
    x1 = load('x1.mat').x1;x2 = load('x2.mat').x2;x3 = load('x3.mat').x3;
    k1 = x1(:,1); b1 = x1(:,2); s1 = x1(:,3); c1 = x1(:,4); % 0.5置信度下的疲劳参数
    k2 = x2(:,1); b2 = x2(:,2); s2 = x2(:,3); c2 = x2(:,4); % 0.9置信度下的疲劳参数
    k3 = x3(:,1); b3 = x3(:,2); s3 = x3(:,3); c3 = x3(:,4); % 0.95置信度下的疲劳参数
    for i=1:N
        Life_x1(i)=fminbnd(@(x)abs(m1(i)/2-(k1(i)-n1(i))/e(i)*(2*x)^b1(i)-s1(i)*(2*x)^c1(i)),100,100000000000000000000);% x为置信度0.5时的寿命
        Life_y1(i)=fminbnd(@(y)abs(m1(i)/2-(k2(i)-n1(i))/e(i)*(2*y)^b2(i)-s2(i)*(2*y)^c2(i)),100,100000000000000000000);% y为置信度0.9时的寿命
        Life_z1(i)=fminbnd(@(z)abs(m1(i)/2-(k3(i)-n1(i))/e(i)*(2*z)^b3(i)-s3(i)*(2*z)^c3(i)),100,100000000000000000000);% z为置信度0.95时的寿命
    end
    Y3 = Life_x1;

    X2 = [x_stress1,x_strain1]';   %第二级输入
    X2_pred = [stress1,strain1]';
    distances = sqrt(sum((X2 - X2_pred).^2, 1));
    W2_1 = diag(kernel(distances, h));
    residual = abs(distances);
    W2_1 = W2_1./(1+beta*residual);
    A2_1 = [ones(1,d); X2];
    theta2_1 = (A2_1 * W2_1 * A2_1') \ (A2_1 * W2_1 * Y3');
    Nf1 = [1, X2_pred'] * theta2_1;
    
    % 第二级计算盘寿命Nf2:
    x_strain2 = ones(N,1);
    x_stress2 = Y2';
    for i = 1:N
        xstrain2 = x_stress2(i,1)/x_optim(i,1)+(x_stress2(i,1)/K_)^(1/n_);
        x_strain2(i,1) = xstrain2;
    end
    m2 = x_strain2; n2 = x_stress2;
    for i=1:N
        Life_x2(i)=fminbnd(@(x)abs(m2(i)/2-(k1(i)-n2(i))/e(i)*(2*x)^b1(i)-s1(i)*(2*x)^c1(i)),100,100000000000000000000);% x为置信度0.5时的寿命
        Life_y2(i)=fminbnd(@(y)abs(m2(i)/2-(k2(i)-n2(i))/e(i)*(2*y)^b2(i)-s2(i)*(2*y)^c2(i)),100,100000000000000000000);% y为置信度0.9时的寿命
        Life_z2(i)=fminbnd(@(z)abs(m2(i)/2-(k3(i)-n2(i))/e(i)*(2*z)^b3(i)-s3(i)*(2*z)^c3(i)),100,100000000000000000000);% z为置信度0.95时的寿命
    end
    Y4 = Life_x2;

    X3 = [x_stress2,x_strain2]';   %第二级输入
    X3_pred = [stress2,strain2]';
    distances = sqrt(sum((X3 - X3_pred).^2, 1));
    W2_2 = diag(kernel(distances, h));
    residual = abs(distances);
    W2_2 = W2_2./(1+beta*residual);
    A2_2 = [ones(1,d); X3];
    theta2_2 = (A2_2 * W2_2 * A2_2') \ (A2_2 * W2_2 * Y4');
    Nf2 = [1, X3_pred'] * theta2_2;
    
    % 第三级输入输出样本（输出是叶盘系统的寿命Nf）:
    if Nf1>Nf2
        Nf = Nf1;
    else
        Nf = Nf2;
    end
    y_pred = Nf;
end