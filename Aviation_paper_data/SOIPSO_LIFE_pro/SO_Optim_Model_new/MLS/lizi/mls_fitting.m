function y_pred = mls_fitting(X, y, X_pred, h)
    % X: 输入的样本数据 (n x d)，n为样本数，d为变量数
    % y: 样本对应的输出值 (n x 1)
    % X_pred: 需要预测的输入点 (m x d)，m为预测点数
    % h: 核函数的带宽 (影响权重的衰减程度)
    
    % 样本数量和变量维度
    [n, d] = size(X);
    % 预测点的数量
    m = size(X_pred, 1);
    
    % 初始化预测结果
    y_pred = zeros(m, 1);
    
    % 核函数的定义 (高斯核)
    kernel = @(r, h) exp(- (r.^2) / (2 * h^2));
    
    % 对每一个预测点，计算局部加权最小二乘拟合
    for i = 1:m
        % 当前预测点
        x_i = X_pred(i, :);
        
        % 计算所有样本点与当前预测点之间的距离
        distances = sqrt(sum((X - x_i).^2, 2));
        
        % 计算权重
        W = diag(kernel(distances, h));
        
        % 构造A矩阵，A = [1, x1, x2, x3, x4]
        A = [ones(n, 1), X];
        
        % 执行加权最小二乘法: 通过正规方程计算系数
        theta = (A' * W * A) \ (A' * W * y);
        
        % 使用拟合的局部多项式来计算预测点的y值
        y_pred(i) = [1, x_i] * theta;
    end
end