function y_pred = RMLS_Fun_B1stress(x_pred)
    % X: 输入的样本数据 (n x d)，n为样本数，d为变量数
    % y: 样本对应的输出值 (n x 1)
    % X_pred: 需要预测的输入点 (m x d)，m为预测点数
    % h: 核函数的带宽 (影响权重的衰减程度)
    h = 1002; % 核函数带宽
    data = readmatrix('BLISK4.csv'); 
    temp = randperm(50);
    x = data(temp,1:4)';
    y = data(temp,5)';
    % 样本数量和变量维度
    [n, d] = size(x);

    % 核函数的定义 (高斯核)
    kernel = @(r, h) exp(- (r.^2) / (2 * h^2));

    beta = 0.3;
    
    % 对每一个预测点，计算局部加权最小二乘拟合
    % 当前预测点
    
    % 计算所有样本点与当前预测点之间的距离
    distances = sqrt(sum((x - x_pred).^2, 1));
    
    % 计算权重
    W = diag(kernel(distances, h));
    
    % 使用鲁棒加权方法：
    residual = abs(distances);
    W = W./(1+beta*residual);

    % 构造A矩阵，A = [1, x1, x2, x3, x4]
    A = [ones(1,d); x];
    
    % 执行加权最小二乘法: 通过正规方程计算系数
    theta = (A * W * A') \ (A * W * y');
    
    % 使用拟合的局部多项式来计算预测点的y值
    y_pred = [1, x_pred'] * theta;
end