%% 示例使用
clc;
clear;
% 生成随机数据
n = 50; % 样本数量

data = readmatrix('BLISK4.csv'); 
temp = randperm(50);
x = data(temp,1:4);
y = data(temp,5);


% 需要预测的点
x_pred = x; % 生成10个需要预测的点

% 调用MLS拟合函数，核函数带宽 h 可以根据数据的情况调整
h = 1002; % 核函数带宽
y_pred = mls_fitting(x, y, x_pred, h);

% 输出预测结果
disp('预测的 y 值:');
disp(y_pred);

%% 可视化原始数据和预测值
error2 = sqrt(sum((y_pred' - y').^2, 2)'./n);

figure;
plot(y, 'r-*', 'DisplayName', 'True value');
hold on;
plot(y_pred, 'b-o', 'DisplayName', 'Prediction value');
legend('show');
string = {['RMSE with MLS =' num2str(error2)]};
title(string);
xlabel('Prediction sample');
ylabel('Prediction result');
xlim([1,n]);
grid on;