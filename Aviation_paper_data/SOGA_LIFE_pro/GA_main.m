%% I. 清空环境变量
%%%%%%%% 代码来源：https://blog.csdn.net/m0_62143653/article/details/129677840 %%%%%%%%%%
% addpath('D:\file\Matlab_file\MATLAB_algorithmic_data\data\14   遗传算法参考程序\Code\gaot');
% savepath;
clear all;
clc;
%% 导入目标函数：
fun = @RMLS_Fun_B1stress;
%% 定义参数：
% 算法参数：
sizepop = 50; % 种群数量
dim = 4; % 基因长度
cross_prob = 0.8; % 交叉概率
lb = [996; 646; 225; 490]'; % 基因取值下界
ub = [1036; 675; 235; 510]'; % 基因取值上界
mut_prob = 0.1; % 变异概率
mut_range = ub - lb; % 变异范围
ger = 100; % 最大迭代次数


% 可靠度优化参数：
R0 = 0.985;        % 第一级优化的最小可靠度约束
R = norminv(R0);
R0_ = 0.98;        % 第二级优化的最小可靠度约束

% 优化前的叶片应力的均值和方差：
Mu0 = 714.201;
Sigma0 = 17.809;
% 优化前的盘应力的均值和方差：
% Mu0 = 942.799;
% Sigma0 = 22.994;

% 假设优化后的叶片应力的均值和方差：
mu = 690;
sigma = 10;
% 假设优化后的盘应力的均值和方差：
% mu = 910;
% sigma = 10;

% Latin超立方抽样参数：
nS = 1e5;nS1 = 1e4;K_ =1420 ;n_ = 0.05;Kt = 2.32;    % 叶片：n_=0.08;盘：n_=0.05
% 损伤：
n = 500;   % 叶片的时候500；盘的时候280
%% 优化前的等效应变和等效平均应力样本：
% 涡轮叶片的应变muX0(1)，平均应力muX0(2)，弹性模量E muX0(3)：
muX0 = [4.4295e-3;714.201;161000];
sigmaX0 = [1.228e-4;17.809;3220];
% 涡轮盘的应变muX0(1)，平均应力muX0(2)，弹性模量E muX0(3)：
% muX0 = [5.7633e-3;942.799;161000];
% sigmaX0 = [1.484e-4;22.994;3220];

y = lhsdesign(nS1,3);
x = [norminv(y(:,1),muX0(1),sigmaX0(1)) , norminv(y(:,2),muX0(2),sigmaX0(2)), norminv(y(:,3),muX0(3),sigmaX0(3))];
%% 四个疲劳参数的样本：
% 置信度为0.5时：
muX1 = [1.5943e3;-0.1058;250.6938;-1.4698];
sigmaX1 = [79.7150;5.5941e-8;20.055504;1.9847e-5];
y1 = lhsdesign(nS1,4);
x1 = [norminv(y1(:,1),muX1(1),sigmaX1(1)),norminv(y1(:,2),muX1(2),sigmaX1(2)),norminv(y1(:,3),muX1(3),sigmaX1(3)),norminv(y1(:,4),muX1(4),sigmaX1(4))];

% 置信度为0.9时：
muX2 = [1.4776e3;-0.1058;76.6160;-1.4698];
sigmaX2 = [73.88;5.5941e-8;6.12928;1.9847e-5];
y2 = lhsdesign(nS1,4);
x2 = [norminv(y2(:,1),muX2(1),sigmaX2(1)),norminv(y2(:,2),muX2(2),sigmaX2(2)),norminv(y2(:,3),muX2(3),sigmaX2(3)),norminv(y2(:,4),muX2(4),sigmaX2(4))];

% 置信度为0.95时：
muX3 = [1.4439e3;-0.1058;53.3976;-1.4698];
sigmaX3 = [72.195;5.5941e-8;4.271808;1.9847e-5];
y3 = lhsdesign(nS1,4);
x3 = [norminv(y3(:,1),muX3(1),sigmaX3(1)),norminv(y3(:,2),muX3(2),sigmaX3(2)),norminv(y3(:,3),muX3(3),sigmaX3(3)),norminv(y3(:,4),muX3(4),sigmaX3(4))];

%% 初始化种群：
step = 1;
while true
    if step == 1
        data = readmatrix('BLISK4_D123.csv','Range','A2:D2001');
        data_y = readmatrix('BLISK4_D123.csv','Range','E2:E2001');
        % 从中随机抽取 50 组数据，保证其符合正态分布
        num_samples = 50;
        pop_x = datasample(data, num_samples);  % 初始化种群位置
        y = datasample(data_y, num_samples)';
        gbest = pop_x;

        % 调用MLS拟合函数，核函数带宽 h 可以根据数据的情况调整
        for j=1:sizepop
            fitness_gbest(j) = fun(pop_x(j,:)');                      % 每个个体的历史最佳适应度
        end
        mu0 = mean(fitness_gbest);
        sigma0 = std(fitness_gbest);
        step = 2;
    elseif step == 2
        if (mu0-mu)/sqrt(sigma0^2-sigma^2)>R            % 这里mu0一定要比mu大，不行就改mu的值！
            zbest = pop_x(1,:);                         % 种群的历史最佳位置
            fitness_zbest = fitness_gbest(1);           % 种群的历史最佳适应度
            for j=1:sizepop
                if fitness_gbest(j) < fitness_zbest     % 如果求最小值，则为<; 如果求最大值，则为>;  
                    zbest = pop_x(j,:);
                    fitness_zbest = fitness_gbest(j);
                end
            end
            break;
        else
            step = 1;
        end
    end
end
%% 迭代：
tol = 1e-6; % 停止条件
St = 1;
while true                   
    if St == 1
        iter = 1;                        %迭代次数
        record_x = zeros(dim, ger);
        record = zeros(ger, 1);          % 记录器
        while iter <= ger
            step = 1;
            while true                   % 以可靠度为约束设置循环
                if step == 1
                    for j=1:sizepop
                        fitness_pop(j) = fun(pop_x(j,:)');                      % 每个个体的历史最佳适应度
                    end
                    if min(fitness_pop) < tol % 满足停止条件
                        break;
                    end
                    prob = fitness_pop / sum(fitness_pop); % 计算每个个体的选择概率
                    cum_prob = cumsum(prob); % 计算累计概率
                    new_pop = zeros(sizepop,dim); % 新种群基因
                    for j = 1:sizepop
                        r = rand; % 生成随机数
                        idx = find(cum_prob >= r,1); % 选择个体
                        new_pop(j,:) = pop_x(idx,:);
                    end
                    pop_x = new_pop; % 更新种群基因
                    for j = 1:2:sizepop
                        if rand < cross_prob % 判断是否进行交叉
                            k = randi(dim-1); % 生成随机交叉点
                            pop_x(j:j+1,k+1:dim) = pop_x(j+1:-1:j,k+1:dim); % 交叉操作
                        end
                    end
                    for j = 1:sizepop
                        if rand < mut_prob % 判断是否进行变异
                            k = randi(dim); % 生成随机变异位
                            pop_x(j,:) = pop_x(j,:) + mut_range * (rand - 0.5); % 变异操作
                        end
                    end
                    for j=1:sizepop
                        fitness_pop(j) = fun(pop_x(j,:)');                      % 每个个体的历史最佳适应度
                    end
                    mu0 = mean(fitness_pop);
                    sigma0 = std(fitness_pop);
                    step = 2;
                elseif step == 2
                    if (mu0-mu)/sqrt(sigma0^2-sigma^2)>R
                        for j=1:sizepop
                            %    新适应度与个体历史最佳适应度做比较
                            if fitness_pop(j) < fitness_gbest(j)       % 如果求最小值，则为<; 如果求最大值，则为>; 
                                gbest(j,:) = pop_x(j,:);               % 更新个体历史最佳位置    
                                fitness_gbest(j) = fitness_pop(j);     % 更新个体历史最佳适应度
                            end   
                            
                            %    个体历史最佳适应度与种群历史最佳适应度做比较
                            if fitness_gbest(j) < fitness_zbest        % 如果求最小值，则为<; 如果求最大值，则为>;  
                                zbest = gbest(j,:);                    % 更新群体历史最佳位置  
                                fitness_zbest = fitness_gbest(j);      % 更新群体历史最佳适应度   
                            end
                        end
                        break;
                    else
                        step = 1;
                    end
                end
            end
            record_x(:,iter) = zbest;
            record(iter) = fitness_zbest;                              %最优值记录,得到优化后的应力或应变。
            iter = iter+1;
        end
        syms Sigma_
        cond = solve(Sigma_-sqrt(Sigma0^2-(Mu0-fitness_zbest)^2/R^2) == 0);
        Sigma = eval(cond);

        % 优化后的应变和应力样本：
        % BP神经网络模型优化求得的叶片或盘的应力均值和方差：
        BP_mu = fitness_zbest; 
        BP_sigma = Sigma;
        
        % 优化后应变x_strain，平均应力muX(1)：
        muX = [BP_mu;161000]; sigmaX = [BP_sigma;3220];
        y_optim = lhsdesign(nS1,2);
        x_optim = [norminv(y_optim(:,1),muX(1),sigmaX(1)) , norminv(y_optim(:,2),muX(2),sigmaX(2))];
        
        % 已知优化后的平均应力，用循环应力-应变公式求优化后的应变x_strain：
        x_strain = ones(nS1,1);
        x_stress = x_optim(:,1);
        for i = 1:nS1
            xstrain = x_stress(i,1)/x_optim(i,2)+(x_stress(i,1)/K_)^(1/n_);
            x_strain(i,1) = xstrain;
        end
        
        % 优化前和优化后的寿命：
        e = x(:,3); m0 = x(:,1); n0 = x(:,2);
        m1 = x_strain; n1 = x_stress;
        k1 = x1(:,1); b1 = x1(:,2); s1 = x1(:,3); c1 = x1(:,4); % 0.5置信度下的疲劳参数
        k2 = x2(:,1); b2 = x2(:,2); s2 = x2(:,3); c2 = x2(:,4); % 0.9置信度下的疲劳参数
        k3 = x3(:,1); b3 = x3(:,2); s3 = x3(:,3); c3 = x3(:,4); % 0.95置信度下的疲劳参数
        
        for i=1:length(m0)
            % 优化前：
            Life_x(i)=fminbnd(@(x)abs(m0(i)/2-(k1(i)-n0(i))/e(i)*(2*x)^b1(i)-s1(i)*(2*x)^c1(i)),100,100000000000000000000);% x为置信度0.5时的寿命
            Life_y(i)=fminbnd(@(y)abs(m0(i)/2-(k2(i)-n0(i))/e(i)*(2*y)^b2(i)-s2(i)*(2*y)^c2(i)),100,100000000000000000000);% y为置信度0.9时的寿命
            Life_z(i)=fminbnd(@(z)abs(m0(i)/2-(k3(i)-n0(i))/e(i)*(2*z)^b3(i)-s3(i)*(2*z)^c3(i)),100,100000000000000000000);% z为置信度0.95时的寿命
            % 优化后：
            Life_Optim_x(i)=fminbnd(@(x)abs(m1(i)/2-(k1(i)-n1(i))/e(i)*(2*x)^b1(i)-s1(i)*(2*x)^c1(i)),100,100000000000000000000);% x为置信度0.5时的寿命
            Life_Optim_y(i)=fminbnd(@(y)abs(m1(i)/2-(k2(i)-n1(i))/e(i)*(2*y)^b2(i)-s2(i)*(2*y)^c2(i)),100,100000000000000000000);% y为置信度0.9时的寿命
            Life_Optim_z(i)=fminbnd(@(z)abs(m1(i)/2-(k3(i)-n1(i))/e(i)*(2*z)^b3(i)-s3(i)*(2*z)^c3(i)),100,100000000000000000000);% z为置信度0.95时的寿命
        end

        z1 = 1-n./Life_Optim_z;
        a1 = find(z1>0);
        R1 = length(a1)/length(z1);
        St = 2;
    elseif St == 2
        if R1 >= R0_            % 损伤可靠度约束：
            break;
        else
            St = 1;
        end
    end
end

%% 迭代结果输出
figure
plot(record);title('迭代的输出数据')
hold on;
disp('优化结果：');
disp(['最优输出值(应力，作为均值)：',num2str(fitness_zbest)]);
disp(['最优输出应力的标准差：',num2str(Sigma)]);
disp('最优随机变量值：');
disp(zbest);

disp('优化前不同置信度下的寿命均值，从0.5，0.9到0.95：');
disp(mean(Life_x));
disp(mean(Life_y));
disp(mean(Life_z));
disp('优化后不同置信度下的寿命均值，从0.5，0.9到0.95：');
disp(mean(Life_Optim_x));
disp(mean(Life_Optim_y));
disp(mean(Life_Optim_z));
