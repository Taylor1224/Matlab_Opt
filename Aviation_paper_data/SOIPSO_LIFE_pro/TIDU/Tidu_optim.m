% 梯度下降算法的简单实现
% 目标函数: f(x) = (x - 3)^2
% 其梯度为: f'(x) = 2 * (x - 3)

% 1. 定义目标函数和其梯度
f = @BP_Fun_Dstress_new;          % 目标函数
df = gradient(f,input_data);      % 目标函数的梯度

% 2. 初始化参数
x = 0; % 初始化 x 的初值
learning_rate = 0.1; % 学习率
max_iter = 100; % 最大迭代次数
tol = 1e-6; % 收敛阈值

% 3. 梯度下降迭代过程
for iter = 1:max_iter
    grad = df(x); % 计算梯度
    x_new = x - learning_rate * grad; % 更新参数
    
    % 检查是否收敛
    if abs(x_new - x) < tol
        fprintf('梯度下降收敛于第 %d 次迭代\n', iter);
        break;
    end
    
    % 更新x
    x = x_new;
end

% 4. 输出结果
fprintf('找到的最优解 x = %.6f\n', x);
fprintf('最小值 f(x) = %.6f\n', f(x));

% 绘制函数和迭代过程
x_vals = linspace(-1, 6, 100);
y_vals = f(x_vals);
figure;
plot(x_vals, y_vals, 'b-', 'LineWidth', 2); hold on;
plot(x, f(x), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('x');
ylabel('f(x)');
title('梯度下降优化结果');
grid on;
