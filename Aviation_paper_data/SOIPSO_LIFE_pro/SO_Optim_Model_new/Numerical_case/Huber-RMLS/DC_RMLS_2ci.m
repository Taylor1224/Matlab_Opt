clc;clear all;
%% 定义包含8个正态分布随机变量的目标函数
E1 = @(x) pi*x(1)/(4*x(3)*sqrt(x(5)/x(8)));
E2 = @(x) x(3)*((x(3)+x(4))/2)/(x(4)*x(3)*((x(4)+x(3))^2+4*((sqrt(x(6)/x(7))-sqrt(x(5)/x(8)))/(sqrt(x(6)/x(7))+ ...
    sqrt(x(5)/x(8))))^2)+x(8)/x(7)*((x(4)+x(3))/2)^2);
E3 = @(x) (x(4)*sqrt(x(6)/x(7))^3+x(3)*sqrt(x(5)/x(8))^3)*sqrt(x(6)/x(7))*4/(sqrt(x(6)/x(7))+sqrt(x(5)/x(8)))^4;


% f = @(x) x(2)-3*x(5)*sqrt((pi*x(1)/(4*x(3)*sqrt(x(5)/x(8))))*(x(3)/(x(4)*x(3) ...
%     *((x(4)+x(3))^2+4*((sqrt(x(6)/x(7))-sqrt(x(5)/x(8)))/(sqrt(x(6)/x(7))+ ...
%     sqrt(x(5)/x(8))))^2)+x(8)/x(7)*((x(4)+x(3))/2)^2)*(x(4)*sqrt(x(6)/x(7))^3+ ...
%     x(3)*sqrt(x(5)/x(8))^3)*sqrt(x(6)/x(7))*4/(sqrt(x(6)/x(7))+sqrt(x(5)/x(8)))^4));

%% 生成随机变量样本
mu = [100, 15, 0.02, 0.05, 0.01, 1, 3, 0.01];  % 均值向量
sigma = [10, 1.5, 0.01, 0.02, 0.001, 0.2, 0.3, 0.001];  % 标准差向量
N = 100; 

X = zeros(N, 8);  % 预分配内存
for i = 1:8
    X(:,i) = normrnd(mu(i), sigma(i), [N,1]);
end
for i = 1:8
    while any(X(:,i) < 0)
        X(X(:,i) < 0,i) = normrnd(mu(i), sigma(i), [sum(X(:,i) < 0), 1]);
    end
end

%% 1. 生成MC输出样本(DC来分布计算目标函数值)：
Y1 = zeros(N, 1);
Y2 = zeros(N, 1);
Y3 = zeros(N, 1);
for i = 1:N
    Y1(i) = E1(X(i, :));
    Y2(i) = E2(X(i, :));
    Y3(i) = E3(X(i, :));
end
Y = X(:,5).*sqrt(Y1.*Y2.*Y3);
MC = real(Y); % 将实部赋值给 y_samples,因为虚部都是0。

%% 6. 生成RMLS输出样本：
% 初始化预测结果
Ypred = zeros(N, 1);

% 核函数的定义 (高斯核)
kernel = @(r, h) exp(- (r.^2) / (2 * h^2));
beta = 0.3;
h = 1002;      % 核函数带宽
X_pred = X;
y_true = MC;
delta = 1.0; % 设置Huber函数的阈值

% 对每一个预测点，计算局部加权最小二乘拟合
for i = 1:N
    % 当前预测点
    x_i = X_pred(i, :);
    
    % 计算所有样本点与当前预测点之间的距离
    distances = sqrt(sum((X - x_i).^2, 2));
    
    % 计算权重（使用Huber函数）
    huber_weights = huber_weight(distances, delta);
    
    % 构造A矩阵，A = [1, x1, x2, x3, x4]
    A = [ones(N,1), X, X.^2];
    
    % 执行加权最小二乘法: 通过正规方程计算系数
    W = diag(huber_weights); % 权重矩阵
    theta = (A' * W * A) \ (A' * W * y_true);
    
    % 使用拟合的局部多项式来计算预测点的y值
    Ypred(i) = [1, x_i, x_i.^2] * theta;
end
RMLS = Ypred;
error6 = sqrt(sum((RMLS' - MC').^2, 2)'./N);

%% 7. 生成DC-RMLS输出样本(DC来分布计算目标函数值)：

% 第一级输入输出样本（输出是E1、E2、E3）:
DY1 = zeros(N, 1);
DY2 = zeros(N, 1);
DY3 = zeros(N, 1);

for i = 1:N
    DY1(i) = E1(X(i, :)); 
    DY2(i) = E2(X(i, :));
    DY3(i) = E3(X(i, :));
end

DYpred1_1 = zeros(N, 1);    % 初始化第一级预测结果
DYpred1_2 = zeros(N, 1);
DYpred1_3 = zeros(N, 1);
X1 = [X(:,1),X(:,3),X(:,5),X(:,8)];
X1_pred1 = X1;              % 第一级的输入预测值
y_true1_1 = DY1;
% 对每一个预测点，计算局部加权最小二乘拟合
for i = 1:N
    % 当前预测点
    x_i = X1_pred1(i, :);
    
    % 计算所有样本点与当前预测点之间的距离
    distances1_1 = sqrt(sum((X1 - x_i).^2, 2));
    
    % 计算权重（使用Huber函数）
    huber_weights1_1 = huber_weight(distances1_1, delta);
    
    % 构造A矩阵，A = [1, x1, x2, x3, x4]
    A1_1 = [ones(N,1), X1, X1.^2];
    
    % 执行加权最小二乘法: 通过正规方程计算系数
    W1_1 = diag(huber_weights1_1); % 权重矩阵
    theta1_1 = (A1_1' * W1_1 * A1_1) \ (A1_1' * W1_1 * y_true1_1);
    % 使用拟合的局部多项式来计算预测点的y值
    DYpred1_1(i) = [1, x_i, x_i.^2] * theta1_1;
end

X2 = [X(:,4),X(:,3),X(:,6),X(:,5),X(:,7),X(:,8)];
X1_pred2 = X2;             % 第一级的输入预测值
y_true1_2 = DY2;
% 对每一个预测点，计算局部加权最小二乘拟合
for i = 1:N
    % 当前预测点
    x_i = X1_pred2(i, :);
    
    % 计算所有样本点与当前预测点之间的距离
    distances1_2 = sqrt(sum((X2 - x_i).^2, 2));
    
    % 计算权重（使用Huber函数）
    huber_weights1_2 = huber_weight(distances1_2, delta);
    
    % 构造A矩阵，A = [1, x1, x2, x3, x4]
    A1_2 = [ones(N,1), X2, X2.^2];
    
    % 执行加权最小二乘法: 通过正规方程计算系数
    W1_2 = diag(huber_weights1_2); % 权重矩阵
    theta1_2 = (A1_2' * W1_2 * A1_2) \ (A1_2' * W1_2 * y_true1_2);
    % 使用拟合的局部多项式来计算预测点的y值
    DYpred1_2(i) = [1, x_i, x_i.^2] * theta1_2;
end

X3 = [X(:,4),X(:,3),X(:,6),X(:,5),X(:,7),X(:,8)];
X1_pred3 = X3;             % 第一级的输入预测值
y_true1_3 = DY3;
% 对每一个预测点，计算局部加权最小二乘拟合
for i = 1:N
    % 当前预测点
    x_i = X1_pred3(i, :);
    
    % 计算所有样本点与当前预测点之间的距离
    distances1_3 = sqrt(sum((X3 - x_i).^2, 2));
    
    % 计算权重（使用Huber函数）
    huber_weights1_3 = huber_weight(distances1_3, delta);
    
    % 构造A矩阵，A = [1, x1, x2, x3, x4]
    A1_3 = [ones(N,1), X3, X3.^2];
    
    % 执行加权最小二乘法: 通过正规方程计算系数
    W1_3 = diag(huber_weights1_3); % 权重矩阵
    theta1_3 = (A1_3' * W1_3 * A1_3) \ (A1_3' * W1_3 * y_true1_3);
    % 使用拟合的局部多项式来计算预测点的y值
    DYpred1_3(i) = [1, x_i, x_i.^2] * theta1_3;
end
DCRMLS1 = [DYpred1_1,DYpred1_2,DYpred1_3];

% 第二级输入输出样本（输出是E=E1.*E2.*E3）:
DY = DY1.*DY2.*DY3; 
X4 = DCRMLS1;
X2_pred = X4;       % 第二级的输入预测值
y_true2 = DY;
Y2pred = zeros(N, 1);   % 初始化第二级预测结果
% 对每一个预测点，计算局部加权最小二乘拟合
for i = 1:N
    % 当前预测点
    x_i = X2_pred(i, :);
    
    % 计算所有样本点与当前预测点之间的距离
    distances2 = sqrt(sum((X4 - x_i).^2, 2));
    
    % 计算权重（使用Huber函数）
    huber_weights2 = huber_weight(distances2, delta);
    
    % 构造A矩阵，A = [1, x1, x2, x3, x4]
    A2 = [ones(N,1), X4, X4.^2];
    
    % 执行加权最小二乘法: 通过正规方程计算系数
    W2 = diag(huber_weights2); % 权重矩阵
    theta2 = (A2' * W2 * A2) \ (A2' * W2 * y_true2);
    % 使用拟合的局部多项式来计算预测点的y值
    Y2pred(i) = [1, x_i, x_i.^2] * theta2;
end
DCRMLS2 = Y2pred;

% 第三级输入输出样本（输出是Z = x(5)*sqrt(E)）:
CY = X(:,5).*sqrt(DY); 
X5 = [X(:,5),DY];
X3_pred = X5;       % 第二级的输入预测值
y_true3 = CY;
Y3pred = zeros(N, 1);   % 初始化第二级预测结果
% 对每一个预测点，计算局部加权最小二乘拟合
for i = 1:N
    % 当前预测点
    x_i = X3_pred(i, :);
    
    % 计算所有样本点与当前预测点之间的距离
    distances3 = sqrt(sum((X5 - x_i).^2, 2));
    
    % 计算权重（使用Huber函数）
    huber_weights3 = huber_weight(distances3, delta);
    
    % 构造A矩阵，A = [1, x1, x2, x3, x4]
    A3 = [ones(N,1), X5, X5.^2];
    
    % 执行加权最小二乘法: 通过正规方程计算系数
    W3 = diag(huber_weights3); % 权重矩阵
    theta3 = (A3' * W3 * A3) \ (A3' * W3 * y_true3);
    % 使用拟合的局部多项式来计算预测点的y值
    Y3pred(i) = [1, x_i, x_i.^2] * theta3;
end
DCRMLS = Y3pred;

error7 = sqrt(sum((DCRMLS' - MC').^2, 2)'./N);
%% 画折现图：（不太行！）
% x = linspace(0, 1, 1000)';
% figure;
% plot(x, MC, '-o', 'LineWidth', 1, 'MarkerSize', 2, 'MarkerFaceColor', 'b');
% hold on;
% plot(x, RSM, '-s', 'LineWidth', 1, 'MarkerSize', 2, 'MarkerFaceColor', 'r');
% hold on;
% plot(x, BPNN, '-s', 'LineWidth', 1, 'MarkerSize', 2, 'MarkerFaceColor', 'y');
% hold on;
% plot(x, SVM, '-s', 'LineWidth', 1, 'MarkerSize', 2, 'MarkerFaceColor', 'g');
% hold on;
% plot(x, Kriging, '-s', 'LineWidth', 1, 'MarkerSize', 2, 'MarkerFaceColor', 'auto');
% hold on;
% plot(x, RMLS, '-s', 'LineWidth', 1, 'MarkerSize', 2, 'MarkerFaceColor', 'auto');
% title('折线图');
% xlabel('X轴');
% ylabel('Y轴');
% legend('sin(x)', 'cos(x)');
% grid on; 

%% 箱线图：
data = [MC,RSM,BPNN,SVM,Kriging,RMLS,DCRMLS];
figure;
boxplot(data, 'Labels', {'MC', 'RSM', 'BPNN', 'SVM', 'Kriging', 'RMLS','DCRMLS'});
ylabel('Limit state function value');
grid on;




