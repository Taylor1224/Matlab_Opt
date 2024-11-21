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



%% 
x = 1:1:40;
y1 = [0.4,0.423,0.456,0.405,0.42,0.425,0.41,0.418,0.45,0.43,0.43,0.44,0.39,0.45,0.44,0.445,0.46,0.435,0.436,0.395,0.43,0.456,0.394,0.397,0.433,0.47,0.435,0.43,0.42,0.432,0.422,0.434,0.463,0.45,0.44,0.435,0.432,0.422,0.52,0.435];
y2 = [0.4342,0.434,0.4341,0.4345,0.4342,0.4343,0.4343,0.4338,0.4337,0.4339,0.4341,0.4344,0.4338,0.4345,0.4342,0.4341,0.434,0.4341,0.434,0.4343,0.4342,0.4342,0.4341,0.4343,0.4341,0.4341,0.4342,0.434,0.4341,0.4338,0.4342,0.4343,0.4342,0.4343,0.4342,0.4343,0.4342,0.434,0.4337,0.434];

plot(x,y1,x,y2)








