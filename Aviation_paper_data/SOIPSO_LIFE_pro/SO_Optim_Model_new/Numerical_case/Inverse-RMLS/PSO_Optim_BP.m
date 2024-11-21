%% 清空环境变量
clc;
clear;

%% 生成训练集、测试集
temp = randperm(1000);   % 叶盘数据为2000   数值实验为200
data1 = load("X.mat").X;
data2 = load("MC.mat").MC;
m = 700;               % 叶盘数据为1500   数值实验为150
P_train = data1(temp(1:m),:)';    % 叶盘数据为1:4   数值实验为1:2
T_train = data2(temp(1:m),:)';      % 5是叶片的应力输出，7是盘的应力输出。  [数值实验]为3
M = size(P_train,2);

P_test = data1(temp(m+1:end),:)';
T_test = data2(temp(m+1:end),:)';   % 5是叶片的应力输出，7是盘的应力输出。
N = size(P_test,2);

%% 归一化
% 输入集
[Pn_train,inputps] = mapminmax(P_train,-1,1);    % !!!归一化中的范围选择-1到1和0到1，得到的net网络结果会有很大差别，这里要选择-1到1之间!!!
Pn_test = mapminmax('apply',P_test,inputps);
% 输出集
[Tn_train,outputps] = mapminmax(T_train,-1,1);   % !!!归一化中的范围选择-1到1和0到1，得到的net网络结果会有很大差别，这里要选择-1到1之间!!!
Tn_test = mapminmax('apply',T_test,outputps);

%% 节点个数
inputnum = size(Pn_train,1);
hiddennum = 10;
outputnum = 1;

%% 建立网络：
net = newff(Pn_train, Tn_train, hiddennum);

%% 设置训练参数(没有优化的bp)
net.trainParam.epochs = 100;        % 训练次数
net.trainParam.lr = 0.1;            % 学习率
net.trainParam.goal = 0.00000001;   % 目标误差
net.trainParam.max_fail = 200;      % 
net.trainParam.showWindow = 0;      % 关闭窗口
 
%% 节点总数：
numsum = inputnum*hiddennum+hiddennum+hiddennum*outputnum+outputnum;


%% 网络训练
net.trainParam.showWindow = 1;   % 先打开训练窗口

net = train(net,Pn_train,Tn_train);

%% 仿真预测：
t_sim1 = sim(net, Pn_train);
t_sim2 = sim(net, Pn_test);

%% 数据反归一化：
T_sim1 = mapminmax('reverse', t_sim1, outputps);
T_sim2 = mapminmax('reverse', t_sim2, outputps);

% ！！！反归一化之后构建能够匹配IPSO算法的BP神经网络！！！
net1_new = train(net, Pn_train, T_sim1);   % 叶片的
net2_new = train(net, Pn_train, T_sim1);   % 盘的

%% 均方根误差：
error1 = sqrt(sum((T_sim1 - T_train).^2, 2)'./M);
error2 = sqrt(sum((T_sim2 - T_test).^2, 2)'./N);


%% 优化前BP神经网络的结果对比：
net = newff(Pn_train, Tn_train, hiddennum);  % 这是前面(41)归一化后的数据
t_sim4 = sim(net, Pn_test);
T_sim4 = mapminmax('reverse', t_sim4, outputps);
error4 = sqrt(sum((T_sim4 - T_test).^2, 2)'./N);
figure
plot(1:N, T_test, 'r-*', 1:N, T_sim4, 'b-o', 'LineWidth', 1)
legend('True value', 'Prediction value');
xlabel('Prediction sample');
ylabel('Prediction result');
string = {['RMSE with BPNN = ' num2str(error4)]};
title(string);
xlim([1,N]);
grid

%% 误差曲线迭代图：
figure;
plot(1:length(BestFit), BestFit, 'LineWidth', 1.5);
xlabel('粒子群迭代次数');
ylabel('适应度值');
xlim([1,length(BestFit)]);
string = {'模型迭代误差变化'};
title(string);
grid on

 
%%  相关指标计算
% R2
R1 = 1 - norm(T_train - T_sim1)^2 / norm(T_train - mean(T_train))^2;
R2 = 1 - norm(T_test  - T_sim2)^2 / norm(T_test  - mean(T_test ))^2;

disp(['训练集数据的R2为：', num2str(R1)])
disp(['测试集数据的R2为：', num2str(R2)])

% MAE
mae1 = sum(abs(T_sim1 - T_train)) ./ M ;
mae2 = sum(abs(T_sim2 - T_test )) ./ N ;

disp(['训练集数据的MAE为：', num2str(mae1)])
disp(['测试集数据的MAE为：', num2str(mae2)])

% MBE
mbe1 = sum(T_sim1 - T_train) ./ M ;
mbe2 = sum(T_sim2 - T_test ) ./ N ;

disp(['训练集数据的MBE为：', num2str(mbe1)])
disp(['测试集数据的MBE为：', num2str(mbe2)])

%%  绘制散点图
sz = 25;
c = 'b';

figure
scatter(T_train, T_sim1, sz, c)
hold on
plot(xlim, ylim, '--k')
xlabel('训练集真实值');
ylabel('训练集预测值');
xlim([min(T_train) max(T_train)])
ylim([min(T_sim1) max(T_sim1)])
title('训练集预测值 vs. 训练集真实值')

figure
scatter(T_test, T_sim2, sz, c)
hold on
plot(xlim, ylim, '--k')
xlabel('测试集真实值');
ylabel('测试集预测值');
xlim([min(T_test) max(T_test)])
ylim([min(T_sim2) max(T_sim2)])
title('测试集预测值 vs. 测试集真实值')
