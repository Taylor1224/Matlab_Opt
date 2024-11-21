% 由于BP网络收敛速度慢并容易陷入局部极值，为了提高BP神经网络模型预测的准确性，提出用粒子群算法来优化BPNN网络，并进行非线性函数拟合。
% 用PSO算法找到最佳的网络初始权值和阈值，以网络正向传播的最小拟合误差为目标函数。
% 将该算法与标准的BP算法进行MATLAB仿真比较。

%% 清空环境变量
tic
clc
clear
close all
format compact

%% 导入数据
data = readmatrix('BLISK4.csv');  % BLISK4.csv(叶盘数据)    Numerical_experiment.csv（数值实验）

%% 生成训练集、测试集
temp = randperm(150);   % 叶盘数据为2000   数值实验为200

m = 100;               % 叶盘数据为1500   数值实验为150
P_train = data(temp(1:m),1:4)';    % 叶盘数据为1:4   数值实验为1:2
T_train = data(temp(1:m),5)';      % 5是叶片的应力输出，7是盘的应力输出。  [数值实验]为3
M = size(P_train,2);

P_test = data(temp(m+1:end),1:4)';
T_test = data(temp(m+1:end),5)';   % 5是叶片的应力输出，7是盘的应力输出。
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
hiddennum = 5;
outputnum = 1;

%% 建立网络：
net = newff(Pn_train, Tn_train, hiddennum);

%% 设置训练参数(没有优化的bp)
net.trainParam.epochs = 100;        % 训练次数
net.trainParam.lr = 0.1;            % 学习率
net.trainParam.goal = 0.00000001;   % 目标误差
net.trainParam.max_fail = 200;      % 
net.trainParam.showWindow = 0;      % 关闭窗口
 
%% PSO 参数设置：
% c1 = 4.494;     % 学习因子
% c2 = 4.494;     % 学习因子
gen = 50;       % 迭代次数
sizepop = 20;   % 种群规模
Vmax = 1.0;
Vmin = -1.0;
popmax = 1.0;   % 最大边界
popmin = -1.0;
w1 = 0.9;
w2 = 0.4;
k = 0.6;
c_max = 5.1;
c_min = 3.9;

%% 节点总数：
numsum = inputnum*hiddennum+hiddennum+hiddennum*outputnum+outputnum;

%% 初始化种群
% 普通初始化种群:
% for i = 1:sizepop
%     pop(i,:) = rands(1, numsum);
%     V(i, :) = rands(1, numsum);
%     fitness(i) = fun1(pop(i,:), hiddennum, net, Pn_train, Tn_train);
% end


% Tent混沌映射初始化种群：
a = 0.499;
z(1,:) = rand(1,numsum); % 随机序列
for j=1:numsum
    for i=2:sizepop
        if z(i-1,j)<a
            z(i,j) = z(i-1,j)/a;
        elseif z(i-1,j)>=a
            z(i,j) = (1-z(i-1,j))/(1-a);
        end
    end
end

pop = popmin + z.*(popmax - popmin);  % 初始化种群位置
V = Vmin + z.*(Vmax - Vmin);  % 初始化种群速度

for i=1:sizepop
    fitness(j) = fun1(pop(i,:), hiddennum, net, Pn_train, Tn_train);                      % 每个个体的历史最佳适应度
end


%% 个体极值和群体极值：
[fitnesszbest, bestindex] = min(fitness);
zbest = pop(bestindex, :);
gbest = pop;
fitnessgbest = fitness;
BestFit = fitnesszbest;

%% 迭代寻优：
for i = 1:gen
    for j = 1:sizepop
        %    更新速度并对速度进行边界处理
        w = (w1-w2)*tan(0.875*(1-(i/gen)^k))+w2;
        c1 = c_min+(c_max-c_min)*(i/gen);
        c2 = c1;

        V(j,:) = w*V(j,:) + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:));
        V(j,(V(j,:) > Vmax)) = Vmax;
        V(j,(V(j,:) < Vmin)) = Vmin;

        %    更新位置并对位置进行边界处理
        pop(j,:) = pop(j,:) + 0.2*V(j,:);   % 位置更新
        pop(j, (pop(j,:) > popmax)) = popmax;
        pop(j, (pop(j,:) < popmin)) = popmin;
        
        %    进行自适应变异
        pos = unidrnd(numsum);
        if rand > 0.85
            pop(j,pos) = rands(1,1);
        end
  
        fitness(j) = fun1(pop(j,:), hiddennum, net, Pn_train, Tn_train);                      % 当前个体的适应度

        %    新适应度与个体历史最佳适应度做比较
        if fitness(j) < fitnessgbest(j)       % 如果求最小值，则为<; 如果求最大值，则为>; 
            gbest(j,:) = pop(j,:);               % 更新个体历史最佳位置   
            fitnessgbest(j) = fitness(j);
        end   
        
        %    个体历史最佳适应度与种群历史最佳适应度做比较
        if fitness(j) < fitnesszbest        % 如果求最小值，则为<; 如果求最大值，则为>;  
            zbest = gbest(j,:);                    % 更新群体历史最佳位置 
            fitnesszbest(j) = fitness(j);
        end
    end
    BestFit = [BestFit, fitnesszbest];
end


%% 用PSO优化的BP网络进行值预测:
% 提取最优初始权值和阈值
w1 = zbest(1:inputnum*hiddennum);
B1 = zbest(inputnum*hiddennum + 1:inputnum*hiddennum + hiddennum);
w2 = zbest(inputnum*hiddennum + hiddennum + 1:inputnum*hiddennum + hiddennum + hiddennum*outputnum);
B2 = zbest(inputnum*hiddennum + hiddennum + hiddennum*outputnum + 1:inputnum*hiddennum + hiddennum + hiddennum*outputnum + outputnum);

% 最优值赋值：
net.iw{1,1} = reshape(w1,hiddennum,inputnum);
net.lw{2,1} = reshape(w2,outputnum,hiddennum);
net.b{1}    = reshape(B1,hiddennum,1);
net.b{2}    = B2';


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

%% 绘图：
% figure
% plot(1:M, T_train, 'r-*', 1:M, T_sim1, 'b-o', 'LineWidth', 1)
% legend('真实值', '预测值');
% xlabel('预测样本');
% ylabel('预测结果');
% string = {'训练集预测结果对比';['RMSE=' num2str(error1)]};
% title(string);
% xlim([1,M]);
% grid

% 测试集的真实数据和根据测试集仿真出来的预测数据的对比：
% PSO优化后的BP神经网络结果对比：
figure
plot(1:N, T_test, 'r-*', 1:N, T_sim2, 'b-o', 'LineWidth', 1)
legend('True value', 'Prediction value');
xlabel('Prediction sample');
ylabel('Prediction result');
string = {['RMSE with IPSO-BP =' num2str(error2)]};
title(string);
xlim([1,N]);
grid on

% 优化前BP神经网络的结果对比：
net = newff(Pn_train, Tn_train, hiddennum);  % 这是前面(41)归一化后的数据
t_sim4 = sim(net, Pn_test);
T_sim4 = mapminmax('reverse', t_sim4, outputps);
error4 = sqrt(sum((T_sim4 - T_test).^2, 2)'./N);
figure
plot(1:N, T_test, 'r-*', 1:N, T_sim4, 'b-o', 'LineWidth', 1)
legend('True value', 'Prediction value');
xlabel('Prediction sample');
ylabel('Prediction result');
string = {['RMSE with QRSM = ' num2str(error4)]};
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
