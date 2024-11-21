function mls_adaptive_complexity(x_data, y_data, xi, h)
    % 自适应复杂度控制的移动最小二乘法
    % x_data: 输入的已知数据点
    % y_data: 输入的已知数据点对应的值
    % xi: 需要拟合的点
    % h: 平滑参数（影响加权的范围）
    
    % 初始化拟合结果
    yi = zeros(size(xi));
    
    % 遍历所有待拟合点
    for i = 1:length(xi)
        % 获取当前拟合点
        x_current = xi(i);
        
        % 自适应选择多项式阶数
        degree = adaptive_polynomial_order(x_current, x_data, y_data);
        
        % 根据权重计算局部拟合
        weights = compute_weights(x_current, x_data, h);
        
        % 使用加权最小二乘法计算局部多项式拟合
        yi(i) = local_weighted_least_squares(x_data, y_data, weights, degree, x_current);
    end
    
    % 可视化结果
    figure;
    plot(x_data, y_data, 'bo'); % 原始数据
    hold on;
    plot(xi, yi, 'r-', 'LineWidth', 2); % 拟合结果
    title('自适应复杂度的移动最小二乘法拟合');
    xlabel('x');
    ylabel('y');
    legend('原始数据', '拟合曲线');
end
