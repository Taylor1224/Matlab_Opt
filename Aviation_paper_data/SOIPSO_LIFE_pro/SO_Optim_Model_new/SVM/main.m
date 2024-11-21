% 支持向量回归 (SVR) 拟合 4 个变量和 1 个输出

% 生成模拟数据
n = 50; % 样本数量

data = readmatrix('BLISK4.csv'); 
temp = randperm(50);
x = data(temp,1:4);
y = data(temp,5);

% 使用支持向量机回归进行拟合
SVMModel = fitrsvm(x, y, 'KernelFunction', 'gaussian', 'Standardize', true);

% 预测拟合值
y_pred = predict(SVMModel, x);

% 绘制真实值与拟合值的对比
error2 = sqrt(sum((y_pred' - y').^2, 2)'./n);

figure;
plot(y, 'r-*', 'DisplayName', 'True value');
hold on;
plot(y_pred, 'b-o', 'DisplayName', 'Prediction value');
legend('show');
string = {['RMSE with SVM =' num2str(error2)]};
title(string);
xlabel('Prediction sample');
ylabel('Prediction result');
xlim([1,n]);
grid on;

