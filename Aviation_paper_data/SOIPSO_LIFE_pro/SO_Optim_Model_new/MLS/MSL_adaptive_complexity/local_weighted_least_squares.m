function yi = local_weighted_least_squares(x, y, weights, degree, xi)
    % 设置范德蒙德矩阵，用于多项式拟合
    X = vander(x); % 生成范德蒙德矩阵
    X = X(:, end-degree:end); % 根据多项式阶数截取合适的列
    
    % 加权最小二乘法的基本公式： (X' * W * X) \ (X' * W * y)
    W = diag(weights); % 权重矩阵
    theta = (X' * W * X) \ (X' * W * y); % 求解多项式系数
    
    % 对目标点xi计算拟合值
    X_current = vander(xi); % 生成目标点的范德蒙德矩阵
    X_current = X_current(end-degree:end); % 同样选择合适的阶数
    yi = X_current * theta; % 计算拟合值
end