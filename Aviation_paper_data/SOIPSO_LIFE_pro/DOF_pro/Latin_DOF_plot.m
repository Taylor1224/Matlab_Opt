clc;clear all;

%% Latin超立方抽样参数：
nS = 1e5;nS1 = 1e4;


%% 优化后的径向变形样本：

% 优化后的径向变形均值和标准差：
muX0 = 0.96194;
sigmaX0 = 0.0116;

y = lhsdesign(nS1,1);
x = [norminv(y(:,1),muX0,sigmaX0)];


% 画直方图：
data= load('BLISK_D123.csv');data1= x;

subplot(1,2,1);
histfit(data(:,33),100,'norm');
title('优化前的径向变形');

subplot(1,2,2);
histfit(data1,100,'norm');
title('优化后的径向变形');


%% 画概率密度曲线图：

dfittool(data1);

plot(evaluateresults(:,1),evaluateresults(:,2),'LineWidth',2);

hold on;

plot(evaluateresults1(:,1),evaluateresults1(:,2),'LineWidth',2);

title('径向变形概率密度曲线');
xlabel('径向变形Y/mm');
ylabel('概率密度Pr');