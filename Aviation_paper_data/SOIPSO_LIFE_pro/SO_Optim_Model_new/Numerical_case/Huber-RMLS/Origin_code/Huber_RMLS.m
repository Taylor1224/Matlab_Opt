% 主程序
function Huber_RMLS()
    % 数据（假设是二维数据）
    x = linspace(0, 10, 100); % 自变量
    y = 2 * x + 1 + 0.5 * randn(1, 100); % 因变量（带噪声）
    
    % 增加一些离群点
    y(95:end) = y(95:end) + 20;
    
    % 参数设置
    delta = 1.0; % Huber函数的阈值
    bandwidth = 2.0; % MLS核宽度

    % 初始化拟合结果
    y_fit = zeros(size(y));

    % 对每个数据点进行加权最小二乘拟合
    for i = 1:length(x)
        % 计算当前点与其他点的距离
        distances = sqrt((x - x(i)).^2 + (y - y(i)).^2); 
        % 根据距离计算权重
        weights = exp(-(distances.^2) / (2 * bandwidth^2));
        
        % 计算当前点的残差
        residuals = y - y_fit;  % 假设y_fit初始值为0，可以改成实际拟合值
        
        % 根据Huber函数计算加权值
        huber_weights = huber_weight(residuals, delta);
        
        % 计算加权残差
        weighted_residuals = huber_weights .* weights;
        
        % 解最小二乘问题（可以用正规方程，或者用矩阵的最小二乘解）
        A = [ones(length(x), 1), x(:)]; % 拟合矩阵（假设拟合的是线性模型）
        W = diag(weighted_residuals); % 权重矩阵
        coef = (A' * W * A) \ (A' * W * y'); % 加权最小二乘解
        
        % 得到拟合结果
        y_fit = A * coef; % 更新拟合结果
    end
    
    % 可视化结果
    figure;
    plot(x, y, 'bo', 'DisplayName', 'Noisy data');
    hold on;
    plot(x, y_fit, 'r-', 'DisplayName', 'MLS with Huber weights');
    legend;
    title('MLS Method with Huber Weight Function');
    xlabel('x');
    ylabel('y');
end


