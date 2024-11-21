clc;
clear;

n = 50;
data = readmatrix('BLISK4.csv'); 
temp = randperm(50);
x = data(temp,1:4);
y = data(temp,5);

fun= @B1_stress_have_Xterms;
y_pred = fun(x);


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