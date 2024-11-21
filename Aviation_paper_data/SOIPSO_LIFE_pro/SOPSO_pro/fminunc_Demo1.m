%%
clc;

%% 这是无约束局部最优
% 定义目标函数
fitnessfcn = @fun;

% 初始猜测
x0 = [996, 646, 225, 490];

% 设置优化选项
options = optimset('Display', 'iter'); % 可以选择显示迭代信息

% 使用fminunc函数进行无约束优化
[x, fval] = fminunc(fitnessfcn, x0, options);

% 打印优化结果
fprintf('最小值点: x = [%f, %f, %f, %f]\n', x(1), x(2), x(3), x(4));
fprintf('最小值: f(x) = %f\n', fval);
