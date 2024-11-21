%% 清空环境
clc;clear all;close all;

%% 目标函数
fun= @BP_NN_Fun_B1stress;

%% 设置种群参数
%   需要自行配置
sizepop = 50;                         % 初始种群个数
dim = 4;                           % 空间维数
ger = 100;                       % 最大迭代次数   
x_ub = [1036; 675; 235; 510]; % 设置位置参数限制(矩阵的形式可以多维)
x_lb = [996; 646; 225; 490];
v_ub = 0.01*ones(dim,1);      % 设置速度限制(系数可以改动！)
v_lb = -0.01*ones(dim,1);
c_1 = 0.8;                       % 惯性权重
c_2 = 0.5;                       % 自我学习因子
c_3 = 0.5;                       % 群体学习因子 

R0 = 0.985;
R = norminv(R0);
R1 = 0.99;
R_ = norminv(R1);

% 假设优化后的叶片应变的均值和方差：
% mu = 4.3e-3;
% sigma = 6e-5;

% 假设优化后的叶片应力的均值和方差：
mu = 690;
sigma = 8;

% 假设优化后的盘应变的均值和方差：
% mu = 5.5e-3;
% sigma = 0.13989e-3;

% 假设优化后的盘应力的均值和方差：
% mu = 930;
% sigma = 11;

%% 生成初始种群
%  首先随机生成初始种群位置
%  然后随机生成初始种群速度
%  然后初始化个体历史最佳位置，以及个体历史最佳适应度
%  然后初始化群体历史最佳位置，以及群体历史最佳适应度


 %% Tent混沌映射初始化种群(位置和速度)：
step = 1;
while true
    if step == 1
        a = 0.499;
        z(:,1) = rand(dim,1); % 随机序列
        for i=1:dim
            for j=2:sizepop
                if z(i,j-1)<a
                    z(i,j) = z(i,j-1)/a;
                elseif z(i,j-1)>=a
                    z(i,j) = (1-z(i,j-1))/(1-a);
                end     
            end
        end

        % 将Tent混沌序列映射为正态分布!!!!!!!!!!!!!!!!!!!!!!!!
        % 使用 norminv 来将(0,1)的均匀分布转换为标准正态分布!!!
        z_ = norminv(z);           % 将混沌序列转为正态分布!!!
        Z = reshape(z_, dim, sizepop);

        pop_x = x_lb + Z.*(x_ub - x_lb);  % 初始化种群位置
        pop_v = v_lb + Z.*(v_ub - v_lb);  % 初始化种群速度
        gbest = pop_x; 
        for j=1:sizepop
            fitness_gbest(j) = fun(pop_x(:,j));                      % 每个个体的历史最佳适应度
        end
        mu0 = mean(fitness_gbest);
        sigma0 = std(fitness_gbest);
        step = 2;
    elseif step == 2
        if (mu0-mu)/sqrt(sigma0^2-sigma^2)>R
            zbest = pop_x(:,1);                         % 种群的历史最佳位置
            fitness_zbest = fitness_gbest(1);           % 种群的历史最佳适应度
            for j=1:sizepop
                if fitness_gbest(j) < fitness_zbest     % 如果求最小值，则为<; 如果求最大值，则为>;  
                    zbest = pop_x(:,j);
                    fitness_zbest=fitness_gbest(j);
                end
            end
            break;
        else
            step = 1;
        end
    end
end

% 作图（样本分布图和散点图）
figure
plot(pop_x(1,:),"black");
hold on;


figure
scatter(pop_x(1,:), pop_x(2,:),'r')   %  (,'filled')可填充圆圈
title("Tent映射初始化种群")
xlabel("x")
ylabel("y")
box on;

figure
scatter(pop_v(1,:), pop_v(2,:), 'b')
title("Tent映射初始化种群")
xlabel("x")
ylabel("y")
box on;

figure;
histogram(pop_x(1,:), 50);
title('初始种群第一个变量的分布');
xlabel('值');
ylabel('频数');

 
 %% 反向学习初始化种群(位置和速度)：
step = 1;
while true
    if step == 1
        x=x_lb+rand(dim,sizepop/2).*(x_ub-x_lb);
        x_obl=x_lb+x_ub-x;          %生成反向学习种群
        pop_x=[x,x_obl];        %合并种群
        
        v=v_lb+rand(dim,sizepop/2).*(v_ub-v_lb);
        v_obl=v_lb+v_ub-v;          %生成反向学习种群
        pop_v=[v,v_obl];        %合并种群
        
        gbest = pop_x; 
        for j=1:sizepop
            fitness_gbest(j) = fun(pop_x(:,j));                      % 每个个体的历史最佳适应度
        end
        mu0 = mean(fitness_gbest);
        sigma0 = std(fitness_gbest);
        step = 2;
    elseif step == 2
        if (mu0-mu)/sqrt(sigma0^2-sigma^2)>R
            zbest = pop_x(:,1);                         % 种群的历史最佳位置
            fitness_zbest = fitness_gbest(1);           % 种群的历史最佳适应度
            for j=1:sizepop
                if fitness_gbest(j) < fitness_zbest     % 如果求最小值，则为<; 如果求最大值，则为>;  
                    zbest = pop_x(:,j);
                    fitness_zbest=fitness_gbest(j);
                end
            end
            break;
        else
            step = 1;
        end
    end
end


% 作图（样本分布图和散点图）：
figure
plot(pop_x(1,:),"black");
hold on;


figure
scatter(pop_x(1,:), pop_x(2,:), 'r')
title("Tent映射初始化种群")
xlabel("x")
ylabel("y")
box on;

figure
scatter(pop_v(1,:), pop_v(2,:), 'b')
title("Tent映射初始化种群")
xlabel("x")
ylabel("y")
box on;
 

%% 普通初始化种群(位置和速度)：
step = 1;
while true
    if step == 1
        for i=1:dim
            for j=1:sizepop
                pop_x(i,j) = x_lb(i)+(x_ub(i) - x_lb(i))*rand;  % 初始种群的位置
                pop_v(i,j) = v_lb(i)+(v_ub(i) - v_lb(i))*rand;  % 初始种群的速度
            end
        end
        gbest = pop_x;                                               % 每个个体的历史最佳位置
        for j=1:sizepop
            fitness_gbest(j) = fun(pop_x(:,j));                      % 每个个体的历史最佳适应度
        end
        mu0 = mean(fitness_gbest);
        sigma0 = std(fitness_gbest);
        step = 2;
    elseif step == 2
        if (mu0-mu)/sqrt(sigma0^2-sigma^2)>R
            zbest = pop_x(:,1);                         % 种群的历史最佳位置
            fitness_zbest = fitness_gbest(1);           % 种群的历史最佳适应度
            for j=1:sizepop
                if fitness_gbest(j) < fitness_zbest     % 如果求最小值，则为<; 如果求最大值，则为>;  
                    zbest = pop_x(:,j);
                    fitness_zbest=fitness_gbest(j);
                end
            end
            break;
        else
            step = 1;
        end
    end
end

% 作图（样本分布图和散点图）：
figure
plot(pop_x(1,:),"black");
hold on;


figure
scatter(pop_x(1,:), pop_x(2,:), 'r')
title("Tent映射初始化种群")
xlabel("x")
ylabel("y")
box on;

figure
scatter(pop_v(1,:), pop_v(2,:), 'b')
title("Tent映射初始化种群")
xlabel("x")
ylabel("y")
box on;
 

%% 粒子群迭代
%    更新速度并对速度进行边界处理    
%    更新位置并对位置进行边界处理
%    进行自适应变异
%    进行约束条件判断并计算新种群各个个体位置的适应度
%    新适应度与个体历史最佳适应度做比较
%    个体历史最佳适应度与种群历史最佳适应度做比较
%    再次循环或结束

iter = 1;                        %迭代次数
record = zeros(ger, 1);          % 记录器
while iter <= ger
    step = 1;
    while true
        if step == 1
            for j=1:sizepop
                %    更新速度并对速度进行边界处理
                pop_v(:,j)= c_1 * pop_v(:,j) + c_2*rand*(gbest(:,j)-pop_x(:,j))+c_3*rand*(zbest-pop_x(:,j));% 速度更新公式
                for i=1:dim
                    if  pop_v(i,j) > v_ub(i)
                        pop_v(i,j) = v_ub(i);
                    end
                    if  pop_v(i,j) < v_lb(i)
                        pop_v(i,j) = v_lb(i);
                    end
                end
                
                %    更新位置并对位置进行边界处理
                pop_x(:,j) = pop_x(:,j) + pop_v(:,j);% 位置更新
                for i=1:dim
                    if  pop_x(i,j) > x_ub(i)
                        pop_x(i,j) = x_ub(i);
                    end
                    if  pop_x(i,j) < x_lb(i)
                        pop_x(i,j) = x_lb(i);
                    end
                end
                
                %    进行自适应变异
                if rand > 0.85
                    i=ceil(dim*rand);
                    pop_x(i,j)=x_lb(i) + (x_ub(i) - x_lb(i)) * rand;
                end
          
                fitness_pop(j) = fun(pop_x(:,j));                      % 当前个体的适应度
            end
            mu0 = mean(fitness_pop);
            sigma0 = std(fitness_pop);
            step = 2;
        elseif step == 2
            if (mu0-mu)/sqrt(sigma0^2-sigma^2)>R
                for j=1:sizepop
                    %    新适应度与个体历史最佳适应度做比较
                    if fitness_pop(j) < fitness_gbest(j)       % 如果求最小值，则为<; 如果求最大值，则为>; 
                        gbest(:,j) = pop_x(:,j);               % 更新个体历史最佳位置    
                        fitness_gbest(j) = fitness_pop(j); % 更新个体历史最佳适应度
                    end   
                    
                    %    个体历史最佳适应度与种群历史最佳适应度做比较
                    if fitness_gbest(j) < fitness_zbest        % 如果求最小值，则为<; 如果求最大值，则为>;  
                        zbest = gbest(:,j);                    % 更新群体历史最佳位置  
                        fitness_zbest = fitness_gbest(j);    % 更新群体历史最佳适应度   
                    end
                end
                break;
            else
                step = 1;
            end
        end
    end
    record(iter) = fitness_zbest;%最大值记录
    iter = iter+1;
end
%% 迭代结果输出
figure
plot(record);title('迭代的输出数据')
hold on;
disp('优化结果：');
disp(['最优输出值：',num2str(fitness_zbest)]);
disp('最优变量值：');
disp(zbest);


%% 求不等式的例子：
%syms x
%cond = x +6 < 10;
%sol = solve(cond, x ,"ReturnConditions",true);
%sol.conditions

%% 计算寿命：
   % 涡轮叶片：
B1_strain_mu0 = 4.4269444e-3;
B1_strain_sigma0 = 9.2887787e-5;
B1_stress_mu0 = 0.71385e3;
B1_stress_sigma0 = 0.13671e2;
B1_strain_mu = 0.0042051;
B1_strain_sigma = 3.7024e-05;
B1_stress_mu = 680.5919;
B1_stress_sigma = 4.8827;

   % 涡轮盘：
D_strain0 = 4.4269444e-3;
D_stress0 = 0.71385e3;
D_strain_mu = 0.0042051;
D_strain_sigma = 3.7024e-05;
D_stress_mu = 680.5919;
D_stress_sigma = 4.8827;

% 优化前后的寿命（用Latin超立方抽样计算寿命）：
  % 优化前的寿命：


  % 优化后的寿命：


% 损伤和可靠度：

% 将（深度学习算法/智能算法）与粒子群算法相结合：
