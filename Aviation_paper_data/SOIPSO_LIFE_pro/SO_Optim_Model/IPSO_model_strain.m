%% 清空环境
clc;clear all; close all;

%% 目标函数
fun= @D_have_Xterms;

%% 设置种群参数
% 粒子群优化参数：
sizepop = 50;                      % 初始种群个数
dim = 4;                           % 空间维数
ger = 100;                         % 最大迭代次数   
x_ub = [1036; 675; 235; 510];      % 设置位置参数限制(矩阵的形式可以多维)
x_lb = [996; 646; 225; 490];
v_ub = 0.01*ones(dim,1);           % 设置速度限制(系数可以改动！)
v_lb = -0.01*ones(dim,1);
c_1 = 0.8;                         % 惯性权重
c_2 = 0.5;                         % 自我学习因子
c_3 = 0.5;                         % 群体学习因子 

R0 = 0.985;                        % 第一级优化的最小可靠度约束
R = norminv(R0);
R0_ = 0.99;                         % 第二级优化的最小可靠度约束

% 优化前的叶片应变的均值和方差：
% Mu0 = 4.4295e-3;
% Sigma0 = 1.228e-4;
% 优化前的盘应变的均值和方差：
Mu0 = 5.7633e-3;
Sigma0 = 1.484e-4;

% 假设优化后的叶片应变的均值和方差：
% mu = 4.3e-3;
% sigma = 6e-5;
% 假设优化后的盘应变的均值和方差：
mu = 5.5e-3;
sigma = 5e-5;

% Latin超立方抽样参数：
nS = 1e5;nS1 = 1e4;K_ =1420 ;n_ = 0.08;Kt = 2.32;

% 损伤：
n = 280;

%% 优化前的等效应变和等效平均应力样本：
% 涡轮叶片的应变muX0(1)，平均应力muX0(2)，弹性模量E muX0(3)：
% muX0 = [4.4295e-3;714.201;161000];
% sigmaX0 = [1.228e-4;17.809;3220];
% 涡轮盘的应变muX0(1)，平均应力muX0(2)，弹性模量E muX0(3)：
muX0 = [5.7633e-3;942.799;161000];
sigmaX0 = [1.484e-4;22.994;3220];

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


%% Tent混沌映射初始化种群(位置和速度)：
%  首先随机生成初始种群位置
%  然后随机生成初始种群速度
%  然后初始化个体历史最佳位置，以及个体历史最佳适应度
%  然后初始化群体历史最佳位置，以及群体历史最佳适应度
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
        pop_x = x_lb + z.*(x_ub - x_lb);  % 初始化种群位置
        pop_v = v_lb + z.*(v_ub - v_lb);  % 初始化种群速度
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

%% 粒子群迭代
%    更新速度并对速度进行边界处理    
%    更新位置并对位置进行边界处理
%    进行自适应变异
%    进行约束条件判断并计算新种群各个个体位置的适应度
%    新适应度与个体历史最佳适应度做比较
%    个体历史最佳适应度与种群历史最佳适应度做比较
%    再次循环或结束

St = 1;
while true
    if St == 1
        iter = 1;                        %迭代次数
        record = zeros(ger, 1);          % 记录器
        while iter <= ger
            step = 1;
            while true                   % 以可靠度为约束设置循环
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
                        pop_x(:,j) = pop_x(:,j) + pop_v(:,j);          % 位置更新
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
                        fitness_pop(j) = fun(pop_x(:,j));              % 当前个体的适应度
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
                                fitness_gbest(j) = fitness_pop(j);     % 更新个体历史最佳适应度
                            end   
                            
                            %    个体历史最佳适应度与种群历史最佳适应度做比较
                            if fitness_gbest(j) < fitness_zbest        % 如果求最小值，则为<; 如果求最大值，则为>;  
                                zbest = gbest(:,j);                    % 更新群体历史最佳位置  
                                fitness_zbest = fitness_gbest(j);      % 更新群体历史最佳适应度   
                            end
                        end
                        break;
                    else
                        step = 1;
                    end
                end
            end
            record(iter) = fitness_zbest;                              %最优值记录,得到优化后的应力或应变。
            iter = iter+1;
        end
        
        syms Sigma_
        cond = solve(Sigma_-sqrt(Sigma0^2-(Mu0-fitness_zbest)^2/R^2) == 0);
        Sigma = eval(cond);

        % 优化后的应变和应力样本：
        % BP神经网络模型优化求得的叶片或盘的应变均值和方差：
        BP_mu = fitness_zbest;
        BP_sigma = Sigma;
        
        % 优化后应变muX(1)，平均应力x_stress：
        muX = [BP_mu;161000]; sigmaX = [BP_sigma;3220];
        y_optim = lhsdesign(nS1,2);
        x_optim = [norminv(y_optim(:,1),muX(1),sigmaX(1)) , norminv(y_optim(:,2),muX(2),sigmaX(2))];
        
        % 已知优化后的应变，用循环应力-应变公式求优化后的应力x_stress：
        x_strain = x_optim(:,1);
        x_stress = ones(nS1,1);
        for i = 1:nS1
            syms xstress
            xstress = abs(vpasolve(xstress/x_optim(i,2)+(xstress/K_)^(1/n_)-x_strain(i,1) == 0,xstress));
            x_stress(i,1) = xstress;
        end
        
        % 优化前和优化后的寿命：
        e = x(:,3); m0 = x(:,1); n0 = x(:,2);
        m = x_strain; n = x_stress;
        k1 = x1(:,1); b1 = x1(:,2); s1 = x1(:,3); c1 = x1(:,4); % 0.5置信度下的疲劳参数
        k2 = x2(:,1); b2 = x2(:,2); s2 = x2(:,3); c2 = x2(:,4); % 0.9置信度下的疲劳参数
        k3 = x3(:,1); b3 = x3(:,2); s3 = x3(:,3); c3 = x3(:,4); % 0.95置信度下的疲劳参数
        
        for i=1:length(m0)
            % 优化前：
            Life_x(i)=fminbnd(@(x)abs(m0(i)/2-(k1(i)-n0(i))/e(i)*(2*x)^b1(i)-s1(i)*(2*x)^c1(i)),100,100000000000000000000);% x为置信度0.5时的寿命
            Life_y(i)=fminbnd(@(y)abs(m0(i)/2-(k2(i)-n0(i))/e(i)*(2*y)^b2(i)-s2(i)*(2*y)^c2(i)),100,100000000000000000000);% y为置信度0.9时的寿命
            Life_z(i)=fminbnd(@(z)abs(m0(i)/2-(k3(i)-n0(i))/e(i)*(2*z)^b3(i)-s3(i)*(2*z)^c3(i)),100,100000000000000000000);% z为置信度0.95时的寿命
            % 优化后：
            Life_Optim_x(i)=fminbnd(@(x)abs(m(i)/2-(k1(i)-n(i))/e(i)*(2*x)^b1(i)-s1(i)*(2*x)^c1(i)),100,100000000000000000000);% x为置信度0.5时的寿命
            Life_Optim_y(i)=fminbnd(@(y)abs(m(i)/2-(k2(i)-n(i))/e(i)*(2*y)^b2(i)-s2(i)*(2*y)^c2(i)),100,100000000000000000000);% y为置信度0.9时的寿命
            Life_Optim_z(i)=fminbnd(@(z)abs(m(i)/2-(k3(i)-n(i))/e(i)*(2*z)^b3(i)-s3(i)*(2*z)^c3(i)),100,100000000000000000000);% z为置信度0.95时的寿命
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
disp(['最优输出值(应变，作为均值)：',num2str(fitness_zbest)]);
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

%% 导出优化前后寿命的样本：
fid = fopen('mu_B1_LIFE0.5.txt','wt');      % wt有重新写入的作用
fprintf(fid,'%g\n',Life_x); 
fid = fopen('mu_B1_LIFE0.9.txt','wt');
fprintf(fid,'%g\n',Life_y); 
fid = fopen('mu_B1_LIFE0.95.txt','wt');
fprintf(fid,'%g\n',Life_z); 
fid = fopen('mu_B1_LIFE0.5_optim.txt','wt');
fprintf(fid,'%g\n',Life_Optim_x); 
fid = fopen('mu_B1_LIFE0.9_optim.txt','wt');
fprintf(fid,'%g\n',Life_Optim_y); 
fid = fopen('mu_B1_LIFE0.95_optim.txt','wt');
fprintf(fid,'%g\n',Life_Optim_z); 


%% 优化前寿命样本分布图：

  % 三种置信度的图放一起，叶片B1寿命样本分布图：
B1 = load('mu_B1_LIFE0.5.txt');B2 = load('mu_B1_LIFE0.9.txt');B3 = load('mu_B1_LIFE0.95.txt');
b1 = load('mu_B1_LIFE0.5_optim.txt');b2 = load('mu_B1_LIFE0.9_optim.txt');b3 = load('mu_B1_LIFE0.95_optim.txt');

subplot(1,2,1);
histfit(B1(:,1),200,'lognorm');
hold on;
histfit(B2(:,1),200,'lognorm');
hold on;
histfit(B3(:,1),200,'lognorm');
xlabel('Blade LCF life before optimization');
ylabel('Relative frequency');
legend('Blade LCF life sample with P=0.5','Fitted lognormal curve with P=0.5','Blade LCF life sample with P=0.9','Fitted lognormal curve with P=0.9','Blade LCF life sample with P=0.95','Fitted lognormal curve with P=0.95');

subplot(1,2,2);
histfit(b1(:,1),200,'lognorm');
hold on;
histfit(b2(:,1),200,'lognorm');
hold on;
histfit(b3(:,1),200,'lognorm');
xlabel('Blade LCF life after optimization');
ylabel('Relative frequency');
legend('Blade LCF life sample with P=0.5','Fitted lognormal curve with P=0.5','Blade LCF life sample with P=0.9','Fitted lognormal curve with P=0.9','Blade LCF life sample with P=0.95','Fitted lognormal curve with P=0.95');


  % 三种置信度的图放一起，盘D寿命样本分布图：
D1 = load('D_LIFE0.5.txt');D2 = load('D_LIFE0.9.txt');D3 = load('D_LIFE0.95.txt');
d1 = load('D_LIFE0.5_optim.txt');d2 = load('D_LIFE0.9_optim.txt');d3 = load('D_LIFE0.95_optim.txt');

subplot(1,2,1);
histfit(D1(:,1),80,'lognorm');
hold on;
histfit(D2(:,1),80,'lognorm');
hold on;
histfit(D3(:,1),80,'lognorm');
xlabel('Disk LCF life before optimization');
ylabel('Relative frequency');
legend('Disk LCF life sample with P=0.5','Fitted lognormal curve with P=0.5','Disk LCF life sample with P=0.9','Fitted lognormal curve with P=0.9','Disk LCF life sample with P=0.95','Fitted lognormal curve with P=0.95');

subplot(1,2,2);
histfit(d1(:,1),80,'lognorm');
hold on;
histfit(d2(:,1),80,'lognorm');
hold on;
histfit(d3(:,1),80,'lognorm');
xlabel('Disk LCF life after optimization');
ylabel('Relative frequency');
legend('Disk LCF life sample with P=0.5','Fitted lognormal curve with P=0.5','Disk LCF life sample with P=0.9','Fitted lognormal curve with P=0.9','Disk LCF life sample with P=0.95','Fitted lognormal curve with P=0.95');


%% 优化前后寿命仿真历史图和分布直方图（放一起）：
  % 置信度0.5下，叶片B1和盘D寿命样本分布图：
D1 = load('D_LIFE0.5.txt');d1 = load('D_LIFE0.5_optim.txt');
subplot(2,2,1);
plot(D1);
hold on ;
subplot(2,2,2);
histfit(D1(:,1),200,'lognorm');
hold on ;
subplot(2,2,3);
plot(d1);
hold on;
subplot(2,2,4);
histfit(d1(:,1),200,'lognorm');

  % 置信度0.9下，叶片B1和盘D寿命样本分布图：
D2 = load('D_LIFE0.9.txt');d2 = load('D_LIFE0.9_optim.txt');
subplot(2,2,1);
plot(D2);
hold on ;
subplot(2,2,2);
histfit(D2(:,1),200,'lognorm');
hold on ;
subplot(2,2,3);
plot(d2);
hold on;
subplot(2,2,4);
histfit(d2(:,1),200,'lognorm');

  % 置信度0.95下，叶片B1和盘D寿命样本分布图：
D3 = load('D_LIFE0.95.txt');d3 = load('D_LIFE0.95_optim.txt');
subplot(2,2,1);
plot(D3);
hold on ;
subplot(2,2,2);
histfit(D3(:,1),200,'lognorm');
hold on ;
subplot(2,2,3);
plot(d3);
hold on;
subplot(2,2,4);
histfit(d3(:,1),200,'lognorm');

%% 将txt文件转为excel文件（一列）：
data1 = importdata('B1_LIFE0.5.txt');
data2 = importdata('B1_LIFE0.9.txt');
data3 = importdata('B1_LIFE0.95.txt');
data4 = importdata('B1_LIFE0.5_optim.txt');
data5 = importdata('B1_LIFE0.9_optim.txt');
data6 = importdata('B1_LIFE0.95_optim.txt');
data7 = importdata('D_LIFE0.5.txt');
data8 = importdata('D_LIFE0.9.txt');
data9 = importdata('D_LIFE0.95.txt');
data10 = importdata('D_LIFE0.5_optim.txt');
data11 = importdata('D_LIFE0.9_optim.txt');
data12 = importdata('D_LIFE0.95_optim.txt');
Life1 = 'B1_LIFE0.5.xlsx';
Life2 = 'B1_LIFE0.9.xlsx';
Life3 = 'B1_LIFE0.95.xlsx';
Life4 = 'B1_LIFE0.5_optim.xlsx';
Life5 = 'B1_LIFE0.9_optim.xlsx';
Life6 = 'B1_LIFE0.95_optim.xlsx';
Life7 = 'D_LIFE0.5.xlsx';h
Life8 = 'D_LIFE0.9.xlsx';
Life9 = 'D_LIFE0.95.xlsx';
Life10 = 'D_LIFE0.5_optim.xlsx';
Life11 = 'D_LIFE0.9_optim.xlsx';
Life12 = 'D_LIFE0.95_optim.xlsx';
xlswrite(Life1,data1);
xlswrite(Life2,data2);
xlswrite(Life3,data3);
xlswrite(Life4,data4);
xlswrite(Life5,data5);
xlswrite(Life6,data6);
xlswrite(Life7,data7);
xlswrite(Life8,data8);
xlswrite(Life9,data9);
xlswrite(Life10,data10);
xlswrite(Life11,data11);
xlswrite(Life12,data12);


%% 将txt文件转为excel文件（一行）：
data1 = importdata('B1_LIFE0.5.txt');
data2 = importdata('B1_LIFE0.9.txt');
data3 = importdata('B1_LIFE0.95.txt');
data4 = importdata('B1_LIFE0.5_optim.txt');
data5 = importdata('B1_LIFE0.9_optim.txt');
data6 = importdata('B1_LIFE0.95_optim.txt');
data7 = importdata('D_LIFE0.5.txt');
data8 = importdata('D_LIFE0.9.txt');
data9 = importdata('D_LIFE0.95.txt');
data10 = importdata('D_LIFE0.5_optim.txt');
data11 = importdata('D_LIFE0.9_optim.txt');
data12 = importdata('D_LIFE0.95_optim.txt');
Life1 = 'B1__LIFE0.5.csv';
Life2 = 'B1__LIFE0.9.csv';
Life3 = 'B1__LIFE0.95.csv';
Life4 = 'B1__LIFE0.5_optim.csv';
Life5 = 'B1__LIFE0.9_optim.csv';
Life6 = 'B1__LIFE0.95_optim.csv';
Life7 = 'D__LIFE0.5.csv';
Life8 = 'D__LIFE0.9.csv';
Life9 = 'D__LIFE0.95.csv';
Life10 = 'D__LIFE0.5_optim.csv';
Life11 = 'D__LIFE0.9_optim.csv';
Life12 = 'D__LIFE0.95_optim.csv';
csvwrite(Life1,data1');
csvwrite(Life2,data2');
csvwrite(Life3,data3');
csvwrite(Life4,data4');
csvwrite(Life5,data5');
csvwrite(Life6,data6');
csvwrite(Life7,data7');
csvwrite(Life8,data8');
csvwrite(Life9,data9');
csvwrite(Life10,data10');
csvwrite(Life11,data11');
csvwrite(Life12,data12');


w = importdata('PRNSOL.txt');
L = 'PRNSOL.csv';
csvwrite(L,w);