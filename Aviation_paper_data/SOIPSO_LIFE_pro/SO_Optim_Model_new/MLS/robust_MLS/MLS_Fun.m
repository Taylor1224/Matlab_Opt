function y_pred = MLS_Fun(pop_x)
    % X: 输入的样本数据 (n x d)，n为样本数，d为变量数
    % y: 样本对应的输出值 (n x 1)
    % X_pred: 需要预测的输入点 (m x d)，m为预测点数
    % h: 核函数的带宽 (影响权重的衰减程度)
    h = 1002; % 核函数带宽
    data = readmatrix('BLISK4_D123.csv','Range','A2:D2001');
    data_y = readmatrix('BLISK4_D123.csv','Range','E2:E2001');
    num_samples = 50;
    pop_x = datasample(data, num_samples)';  % 初始化种群位置
    x_pred = x1;
    y = datasample(data_y, num_samples)';
    
    % 样本数量和变量维度
    [n, d] = size(x);
    % 预测点的数量
    m = size(x_pred, 1);
    
    % 初始化预测结果
    y_pred = zeros(m, 1);
    
    % 核函数的定义 (高斯核)
    kernel = @(r, h) exp(- (r.^2) / (2 * h^2));

    beta = 0.3;
    
    % 对每一个预测点，计算局部加权最小二乘拟合
    for i = 1:m
        % 当前预测点
        x_i = x_pred(i, :);
        
        % 计算所有样本点与当前预测点之间的距离
        distances = sqrt(sum((x - x_i).^2, 2));
        
        % 计算权重
        W = diag(kernel(distances, h));
        
        % 使用鲁棒加权方法：
        residual = abs(distances);
        W = W./(1+beta*residual);

        % 构造A矩阵，A = [1, x1, x2, x3, x4]
        A = [ones(n, 1), x];
        
        % 执行加权最小二乘法: 通过正规方程计算系数
        theta = (A' * W * A) \ (A' * W * y);
        
        % 使用拟合的局部多项式来计算预测点的y值
        y_pred(i) = [1, x_i] * theta;
    end
end