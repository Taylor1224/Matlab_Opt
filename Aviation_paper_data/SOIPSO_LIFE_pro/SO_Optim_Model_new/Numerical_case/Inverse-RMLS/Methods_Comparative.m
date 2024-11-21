clc;clear;
%% 定义包含8个正态分布随机变量的目标函数
f = @(x) x(2)-3*x(5)*sqrt((pi*x(1)/(4*x(3)*sqrt(x(5)/x(8))))*(x(3)*((x(3)+x(4))/2)/(x(4)*x(3) ...
    *((x(4)+x(3))^2+4*((sqrt(x(6)/x(7))-sqrt(x(5)/x(8)))/(sqrt(x(6)/x(7))+ ...
    sqrt(x(5)/x(8))))^2)+x(8)/x(7)*((x(4)+x(3))/2)^2)*(x(4)*sqrt(x(6)/x(7))^3+ ...
    x(3)*sqrt(x(5)/x(8))^3)*sqrt(x(6)/x(7))*4/(sqrt(x(6)/x(7))+sqrt(x(5)/x(8)))^4));

%% 生成随机变量样本
% 定义8个随机变量的均值和方差
mu = [100, 15, 0.02, 0.05, 0.01, 1, 3, 0.01];  % 均值向量
sigma = [10, 1.5, 0.01, 0.02, 0.001, 0.2, 0.3, 0.001];  % 标准差向量

% Step 2: 设定蒙特卡洛模拟次数
N = 1000;  % 样本数

% Step 3: 生成8个正态分布随机变量的样本 (每列是一个变量)
X = zeros(N, 8);  % 预分配内存
for i = 1:8
    X(:,i) = normrnd(mu(i), sigma(i), [N,1]);
end

%% 1. 生成MC输出样本：
Y = zeros(N, 1);  % 初始化结果向量
for i = 1:N
    Y(i) = f(X(i, :));  % 计算每个样本的函数值
end
MC = real(Y); % 将实部赋值给 y_samples,因为虚部都是0。
% disp(MC);

%% 2. 生成RSM输出样本：
% X = load("X.mat").X;MC = load("MC.mat").MC;
N = 1000;

X_design = [X, ...  % 一阶项
            X(:,1).*X(:,2), X(:,1).*X(:,3), X(:,1).*X(:,4), X(:,1).*X(:,5), X(:,1).*X(:,6), X(:,1).*X(:,7), X(:,1).*X(:,8), ...  % 交互项
            X(:,2).*X(:,3), X(:,2).*X(:,4), X(:,2).*X(:,5), X(:,2).*X(:,6), X(:,2).*X(:,7), X(:,2).*X(:,8), ...  % 交互项
            X(:,3).*X(:,4), X(:,3).*X(:,5), X(:,3).*X(:,6), X(:,3).*X(:,7), X(:,3).*X(:,8), ...  % 交互项
            X(:,4).*X(:,5), X(:,4).*X(:,6), X(:,4).*X(:,7), X(:,4).*X(:,8), ...  % 交互项
            X(:,5).*X(:,6), X(:,5).*X(:,7), X(:,5).*X(:,8), ...  % 交互项
            X(:,6).*X(:,7), X(:,6).*X(:,8), ...  % 交互项
            X(:,7).*X(:,8), ...  % 交互项
            X(:,1).^2, X(:,2).^2, X(:,3).^2, X(:,4).^2, X(:,5).^2, X(:,6).^2, X(:,7).^2, X(:,8).^2];  % 二次项

% 第三步：使用 fitlm 进行回归拟合
RSM_model = fitlm(X_design, MC);
% disp(RSM);
RSM = RSM_model.Fitted;

error2 = sqrt(sum((RSM' - MC').^2, 2)'./N);  % 均方根误差 

%% 3. 生成BPNN输出样本：
temp = randperm(150);   % 叶盘数据为2000   数值实验为200
m = 100;               % 叶盘数据为1500   数值实验为150
P_train = X(temp(1:m),:)';    % 叶盘数据为1:4   数值实验为1:2
T_train = MC(temp(1:m),:)';      % 5是叶片的应力输出，7是盘的应力输出。  [数值实验]为3

P_test = X(temp(m+1:end),:)';
T_test = MC(temp(m+1:end),:)';   % 5是叶片的应力输出，7是盘的应力输出。

P_sim = X';

% 归一化：
% 输入集
[Pn_train,inputps] = mapminmax(P_train,-1,1);    % !!!归一化中的范围选择-1到1和0到1，得到的net网络结果会有很大差别，这里要选择-1到1之间!!!
Pn_test = mapminmax('apply',P_test,inputps);
Pn_sim = mapminmax('apply',P_sim,inputps);
% 输出集
[Tn_train,outputps] = mapminmax(T_train,-1,1);   % !!!归一化中的范围选择-1到1和0到1，得到的net网络结果会有很大差别，这里要选择-1到1之间!!!
Tn_test = mapminmax('apply',T_test,outputps);

% 节点个数
inputnum = size(Pn_train,1);
hiddennum = 5;
outputnum = 1;

% 建立网络：
net = newff(Pn_train, Tn_train, hiddennum);

% 设置训练参数(没有优化的bp)
net.trainParam.epochs = m;          % 训练次数
net.trainParam.lr = 0.1;            % 学习率
net.trainParam.goal = 0.00000001;   % 目标误差
net.trainParam.max_fail = 200;      
net.trainParam.showWindow = 0;      % 关闭窗口

% 网络训练
net.trainParam.showWindow = 1;   % 先打开训练窗口
net = train(net,Pn_train,Tn_train);

% 优化前BP神经网络的结果对比：
net = newff(Pn_train, Tn_train, hiddennum);  % 这是前面归一化后的数据
t_sim4 = sim(net, Pn_test);                       % 仿真预测
T_sim4 = mapminmax('reverse', t_sim4, outputps);  % 数据反归一化   T_sim4为预测输出

t_sim = sim(net, Pn_sim);                       % 仿真预测
T_sim = mapminmax('reverse', t_sim, outputps);  % 数据反归一化   T_sim为预测输出

BPNN = T_sim';

error3 = sqrt(sum((BPNN' - MC').^2, 2)'./N);  % 均方根误差  
% figure
% plot(1:N, T_test, 'r-*', 1:N, T_sim4, 'b-o', 'LineWidth', 1)
% legend('True value', 'Prediction value');
% xlabel('Prediction sample');
% ylabel('Prediction result');
% string = {['RMSE with BPNN = ' num2str(error4)]};
% title(string);
% xlim([1,N]);
% grid

%% 4. 生成SVM输出样本：
% 使用支持向量机回归进行拟合
SVMModel = fitrsvm(X, MC, 'KernelFunction', 'gaussian', 'Standardize', true);

% 预测拟合值
y_pred = predict(SVMModel, X);
SVM = y_pred;

% 绘制真实值与拟合值的对比
error4 = sqrt(sum((SVM' - MC').^2, 2)'./N);

%% 5. 生成Kriging输出样本：
gprMdl = fitrgp(X, MC, 'KernelFunction', 'squaredexponential');
Y_pred = predict(gprMdl, X);
Kriging = Y_pred;
error5 = sqrt(sum((Kriging' - MC').^2, 2)'./N);


%% 6. 生成RMLS输出样本：
% 初始化预测结果
Ypred = zeros(N, 1);

% 核函数的定义 (高斯核)
kernel = @(r, h) exp(- (r.^2) / (2 * h^2));
beta = 0.3;
h = 1002;      % 核函数带宽
X_pred = X;

% 对每一个预测点，计算局部加权最小二乘拟合
for i = 1:N
    % 当前预测点
    x_i = X_pred(i, :);
    % 计算所有样本点与当前预测点之间的距离
    distances = sqrt(sum((X - x_i).^2, 2));
    % 计算权重
    W = diag(kernel(distances, h));
    % 使用鲁棒加权方法：
    residual = abs(distances);
    W = W./(1+beta*residual);
    % 构造A矩阵，A = [1, x1, x2, x3, x4]
    A = [ones(N,1), X];
    % 执行加权最小二乘法: 通过正规方程计算系数
    theta = (A' * W * A) \ (A' * W * MC);
    % 使用拟合的局部多项式来计算预测点的y值
    Ypred(i) = [1, x_i] * theta;
end
RMLS = Ypred;

error6 = sqrt(sum((RMLS' - MC').^2, 2)'./N);

%% 画折现图：（不太行！）
% x = linspace(0, 1, 1000)';
% figure;
% plot(x, MC, '-o', 'LineWidth', 1, 'MarkerSize', 2, 'MarkerFaceColor', 'b');
% hold on;
% plot(x, RSM, '-s', 'LineWidth', 1, 'MarkerSize', 2, 'MarkerFaceColor', 'r');
% hold on;
% plot(x, BPNN, '-s', 'LineWidth', 1, 'MarkerSize', 2, 'MarkerFaceColor', 'y');
% hold on;
% plot(x, SVM, '-s', 'LineWidth', 1, 'MarkerSize', 2, 'MarkerFaceColor', 'g');
% hold on;
% plot(x, Kriging, '-s', 'LineWidth', 1, 'MarkerSize', 2, 'MarkerFaceColor', 'auto');
% hold on;
% plot(x, RMLS, '-s', 'LineWidth', 1, 'MarkerSize', 2, 'MarkerFaceColor', 'auto');
% title('折线图');
% xlabel('X轴');
% ylabel('Y轴');
% legend('sin(x)', 'cos(x)');
% grid on; 

%% 箱线图：
data = [MC,Kriging,BPNN,RSM,RMLS,SVM];
figure;
boxplot(data, 'Labels', {'MC', 'RSM', 'BPNN', 'SVM', 'Kriging', 'RMLS'});
ylabel('Limit state function value');
grid on;




