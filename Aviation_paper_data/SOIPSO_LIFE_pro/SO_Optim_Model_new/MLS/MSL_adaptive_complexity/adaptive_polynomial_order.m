function degree = adaptive_polynomial_order(xi, x, y)
    % 根据局部方差选择多项式阶数
    local_window = abs(x - xi) < 0.2; % 取局部窗口范围
    local_variance = var(y(local_window)); % 计算局部方差
    
    % 根据局部方差自适应选择多项式阶数
    if local_variance < 0.01
        degree = 1; % 线性拟合
    elseif local_variance < 0.1
        degree = 2; % 二次拟合
    else
        degree = 3; % 三次拟合
    end
end

% 计算自适应的多项式阶数