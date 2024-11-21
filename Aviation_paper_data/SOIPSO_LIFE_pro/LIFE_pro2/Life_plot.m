clc;clear all;

%% 优化前的LCF的应变范围①和平均应力②的直方图：
% 叶片B1：
  % ①应变范围（左）、②平均应力（右）
data_life= load('BLISK_D123.csv');
subplot(1,2,1);
histfit(data_life(:,16),200,'norm');
xlabel('Blade strain range');
ylabel('Relative frequency');

subplot(1,2,2);
histfit(data_life(:,12),200,'norm');
xlabel('Blade mean stress');
ylabel('Relative frequency');

% 盘D：
  % ①应变范围（左）、②平均应力（右）
subplot(1,2,1);
histfit(data_life(:,24),100,'norm');
xlabel('Disk strain range');
ylabel('Relative frequency');

subplot(1,2,2);
histfit(data_life(:,20),100,'norm');
xlabel('Disk mean stress');
ylabel('Relative frequency');



%% 优化后的LCF的应变范围和平均应力的直方图：
% Latin超立方抽样参数：
nS = 1e5;nS1 = 1e4;

% 叶片B1：
  % ①应变范围（左）、②平均应力（右） 
muX0 = 0.0042051;
sigmaX0 = 3.7024e-05;

y0 = lhsdesign(nS1,1);
x0 = [norminv(y(:,1),muX0,sigmaX0)];
data1= x0;

subplot(1,2,1);
histfit(data1,100,'norm');
xlabel('Disk strain range');
ylabel('Relative frequency');

subplot(1,2,2);
histfit(data1,100,'norm');
xlabel('Disk strain range');
ylabel('Relative frequency');

% 盘D：
  % ①应变范围（左）、②平均应力（右）
muX0 = 0.0042051;
sigmaX0 = 3.7024e-05;

y = lhsdesign(nS1,1);
x = [norminv(y(:,1),muX0,sigmaX0)];


%% 画概率密度曲线图：

dfittool(data1);

plot(evaluateresults(:,1),evaluateresults(:,2),'LineWidth',2);

hold on;

plot(evaluateresults1(:,1),evaluateresults1(:,2),'LineWidth',2);

title('径向变形概率密度曲线');
xlabel('径向变形Y/mm');
ylabel('概率密度Pr');



%%  三种优化算法的迭代曲线对比图：
x = 1:1:40;
y1 = [2.0191,2.0563,2.0956,2.1095,2.1533,2.1865,2.1988,2.2212,2.2356,2.2574,2.2635,2.2844,2.2998,2.3223,2.3456,2.3554,2.3747,2.3955,2.4012,2.4266,2.4339,2.4566,2.4752,2.4906,2.5023,2.5123,2.5236,2.5236,2.5236,2.5236,2.5236,2.5236,2.5236,2.5236,2.5236,2.5236,2.5236,2.5236,2.5236,2.5236];
y2 = [2.0191,2.1964,2.3021,2.3998,2.4765,2.5502,2.6123,2.6956,2.7645,2.7998,2.8423,2.8756,2.8968,2.9356,2.9658,2.9799,2.9898,2.9956,3.0125,3.0210,3.0238,3.0238,3.0238,3.0238,3.0238,3.0238,3.0238,3.0238,3.0238,3.0238,3.0238,3.0238,3.0238,3.0238,3.0238,3.0238,3.0238,3.0238,3.0238,3.0238];
y3 = [2.0191,2.2320,2.4792,2.5986,2.6789,2.7523,2.8165,2.8745,2.9236,3.0774,3.1458,3.2354,3.2989,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272,3.3272];

plot(x,y1,x,y2,x,y3)


%% 





