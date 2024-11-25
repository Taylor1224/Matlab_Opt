%% 清空环境
clc;clear all;close all;

%% 目标函数
fun= @RMLS_Fun_Nf2_1;

%% 设置种群参数
% 粒子群优化参数：
sizepop = 50;                      % 初始种群个数
dim = 4;                           % 空间维数
ger = 100;                         % 最大迭代次数   
x_ub = [1036; 675; 235; 510];      % 设置位置参数限制(矩阵的形式可以多维)
x_lb = [996; 646; 225; 490];
v_ub = 0.01*ones(dim,1);           % 设置速度限制(系数可以改动！)
v_lb = -0.01*ones(dim,1);
c_1 = 0.5;                         % 自我学习因子
c_2 = 0.5;                         % 群体学习因子 
w1 = 0.9;
w2 = 0.4;
k = 0.6;

R0 = 0.985;                        % 第一级优化的最小可靠度约束
R = norminv(R0);
R0_ = 0.98;                         % 第二级优化的最小可靠度约束

% 优化前系统寿命的均值和方差：
Mu0 = 1.9180e3;      % 0.5置信度下:1.9180e3; 0.9置信度下:790.5919; 0.95置信度下:604.2321;
Sigma0 = 277.1488;   % 0.5置信度下:277.1488; 0.9置信度下:104.1206; 0.95置信度下:80.5513;

% 假设优化后系统寿命的均值和方差：
mu = 2.1892e3;       % 0.5置信度下:2.1892e3; 0.9置信度下:893.2872; 0.95置信度下:680.8471;
sigma = 313.9817;    % 0.5置信度下:313.9817; 0.9置信度下:112.2232; 0.95置信度下:88.9149;

% Latin超立方抽样参数：
nS = 1e5;nS1 = 1e4;K_ =1420 ;n_ = 0.05;Kt = 2.32;    % 叶片：n_=0.08;盘：n_=0.05

% 损伤：
n = 500;   % 叶片的时候500；盘的时候280


%% 初始化种群(位置和速度)(使用MC抽样原始数据中的50组作为初始种群位置)：
step = 1;
while true
    if step == 1
        data = readmatrix('BLISK4_D123.csv','Range','A2:D2001');
        data_y = readmatrix('BLISK4_D123.csv','Range','E2:E2001');
        % 从中随机抽取 50 组数据，保证其符合正态分布
        num_samples = 50;
        pop_x = datasample(data, num_samples)';  
        y = datasample(data_y, num_samples)';

        pop_v = 0.5 * rands(4,50); 
        gbest = pop_x;

        % 调用MLS拟合函数，核函数带宽 h 可以根据数据的情况调整
        for j=1:sizepop
            fitness_gbest(j,:) = fun(pop_x(:,j));                     
        end
        mu0 = mean(fitness_gbest(1));
        sigma0 = std(fitness_gbest);
        step = 2;
    elseif step == 2
        if (mu0-mu)/sqrt(sigma0^2-sigma^2)>R            
            zbest = pop_x(:,1);                        
            fitness_zbest = fitness_gbest(1);          
            for j=1:sizepop
                if -fitness_gbest(j) < -fitness_zbest    
                    zbest = pop_x(:,j);
                    fitness_zbest = fitness_gbest(j);
                end
            end
            break;
        else
            step = 1;
        end
    end
end

%% 粒子群迭代
St = 1;
while true                   
    if St == 1
        iter = 1;                       
        record_x = zeros(dim, ger);
        record = zeros(ger, 1);         
        while iter <= ger
            step = 1;
            while true                   % 以可靠度为约束设置循环
                if step == 1
                    for j=1:sizepop
                        %    更新速度并对速度进行边界处理
                        w = (w1-w2)*tan(0.875*(1-(iter/ger)^k))+w2;
                        pop_v(:,j)= w*pop_v(:,j) + c_1*rand*(gbest(:,j)-pop_x(:,j))...
                        +c_2*rand*(zbest-pop_x(:,j));
                        for i=1:dim
                            if  pop_v(i,j) > v_ub(i)
                                pop_v(i,j) = v_ub(i);
                            end
                            if  pop_v(i,j) < v_lb(i)
                                pop_v(i,j) = v_lb(i);
                            end
                        end
                        
                        %    更新位置并对位置进行边界处理
                        pop_x(:,j) = pop_x(:,j) + pop_v(:,j);         
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
                        fitness_pop(j) = fun(pop_x(:,j));             
                    end
                    mu0 = mean(fitness_pop);
                    sigma0 = std(fitness_pop);
                    step = 2;
                elseif step == 2
                    if (mu0-mu)/sqrt(sigma0^2-sigma^2)>R
                        for j=1:sizepop
                            %    新适应度与个体历史最佳适应度做比较
                            if -fitness_pop(j) < -fitness_gbest(j)        
                                gbest(:,j) = pop_x(:,j);                  
                                fitness_gbest(j) = fitness_pop(j);     
                            end   
                            
                            %    个体历史最佳适应度与种群历史最佳适应度做比较
                            if -fitness_gbest(j) < -fitness_zbest          
                                zbest = gbest(:,j);                     
                                fitness_zbest = fitness_gbest(j);        
                            end
                        end
                        break;
                    else
                        step = 1;
                    end
                end
            end
            record_x(:,iter) = zbest;
            record(iter) = fitness_zbest;                              
            iter = iter+1;
        end
        syms Sigma_
        cond = solve(Sigma_-sqrt(Sigma0^2-(Mu0-fitness_zbest)^2/R^2) == 0);
        Sigma = eval(cond);
        Sigma = imag(Sigma);

        z1 = 1-n./fitness_zbest;
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
disp(['最优输出值系统寿命：',num2str(fitness_zbest)]);
disp(['最优输出应力的标准差：',num2str(Sigma)]);
disp('最优随机变量值：');
disp(zbest);


%% IPSO优化算法迭代三维图：
X = record_x;
Y = meshgrid(record',record');
[x1,x2] = meshgrid(X(1,:),X(2,:));
figure;
surf(x1,x2,Y);


