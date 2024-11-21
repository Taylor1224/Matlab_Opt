clc;clear all;
%% 定义包含8个正态分布随机变量的目标函数
E1 = @(x) pi*x(1)/(4*x(3)*sqrt(x(5)/x(8)));
E2 = @(x) x(3)*((x(3)+x(4))/2)/(x(4)*x(3)*((x(4)+x(3))^2+4*((sqrt(x(6)/x(7))-sqrt(x(5)/x(8)))/(sqrt(x(6)/x(7))+ ...
    sqrt(x(5)/x(8))))^2)+x(8)/x(7)*((x(4)+x(3))/2)^2);
E3 = @(x) (x(4)*sqrt(x(6)/x(7))^3+x(3)*sqrt(x(5)/x(8))^3)*sqrt(x(6)/x(7))*4/(sqrt(x(6)/x(7))+sqrt(x(5)/x(8)))^4;

%% 生成随机变量样本
mu = [100, 15, 0.02, 0.05, 0.01, 1, 3, 0.01];  % 均值向量
sigma = [10, 1.5, 0.01, 0.02, 0.001, 0.2, 0.3, 0.001];  % 标准差向量
N = 500; 

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
    A = [ones(N,1), X];
    
    % 执行加权最小二乘法: 通过正规方程计算系数
    W = diag(huber_weights); % 权重矩阵
    theta = (A' * W * A) \ (A' * W * y_true);
    
    % 使用拟合的局部多项式来计算预测点的y值
    Ypred(i) = [1, x_i] * theta;
end
RMLS = Ypred;
error6 = sqrt(sum((RMLS' - MC').^2, 2)'./N);

