clc;
clear 
close 

data = readmatrix('BLISK4.csv'); 
temp = randperm(150);
m = 100;
train_x = data(temp(1:m),1:4);
train_y = data(temp(1:m),5);
M = size(train_x',2);
test_x = data(temp(m+1:end),1:4);
test_y = data(temp(m+1:end),5);
N = size(test_x',2);

%无需归一化，直接赋值

nva = 4;%输入变量数量
finaltheta = ones(1,nva);
thedown = 0.01*ones(1,nva);
theup = 100*ones(1,nva);%上下界可修改
%训练
[dmodel, perf] = dacefit(train_x, train_y, @regpoly0, @corrgauss, finaltheta,thedown,theup);
% @regpoly0，@regpoly1，@regpoly2可选，分别是常数多项式，一次多项式，二次多项式

dmodel;%模型结构

[predy] = predictor(test_x, dmodel);%预测


error2 = sqrt(sum((predy' - test_y').^2, 2)'./N);


figure
plot(1:N, test_y', 'r-*', 1:N, predy', 'b-o', 'LineWidth', 1)
legend('True value', 'Prediction value');
xlabel('Prediction sample');
ylabel('Prediction result');
string = {['RMSE with Kriging =' num2str(error2)]};
title(string);
xlim([1,N]);
grid on



