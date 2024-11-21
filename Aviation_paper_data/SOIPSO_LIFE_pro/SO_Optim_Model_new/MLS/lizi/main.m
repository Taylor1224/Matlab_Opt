%% 示例使用
clc;
clear;
% 生成随机数据
n = 50; % 样本数量
X = rand(n, 4); % 4个输入变量
y = 1.5 * X(:,1) - 2.3 * X(:,2) + 0.7 * X(:,3) + 3.1 * X(:,4) + randn(n, 1); % 输出，加上噪声

% 需要预测的点
X_pred = rand(10, 4); % 生成10个需要预测的点

% 调用MLS拟合函数，核函数带宽 h 可以根据数据的情况调整
h = 0.5; % 核函数带宽
y_pred = mls_fitting(X, y, X_pred, h);

% 输出预测结果
disp('预测的 y 值:');
disp(y_pred);

%% 可视化原始数据和预测值（仅限二维展示，可扩展到多维）
% 为了简化可视化，将前两个维度绘制成2D图
figure;
scatter3(X(:,1), X(:,2), y, 'filled', 'DisplayName', '原始数据');
hold on;
scatter3(X_pred(:,1), X_pred(:,2), y_pred, 'r', 'filled', 'DisplayName', '预测点');
xlabel('x1');
ylabel('x2');
zlabel('y');
legend('show');
title('MLS拟合效果');
grid on;