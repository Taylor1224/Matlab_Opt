% 生成一些示例数据
clc;
clear;

x_data = linspace(0, 2*pi, 100);
y_data = sin(x_data) + 0.1*randn(size(x_data)); % 添加少量噪声的正弦函数

% 需要拟合的目标点
xi = linspace(0, 2*pi, 200);

% 平滑参数
h = 0.3;

% 调用自适应MLS函数进行拟合
mls_adaptive_complexity(x_data, y_data, xi, h);
