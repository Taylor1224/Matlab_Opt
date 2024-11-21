%% 清空环境
clc;clear all;close all;

%% 目标函数
fun= @RMLS_Fun_B1stress;

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


%% 初始化种群(位置和速度)(使用MC抽样原始数据中的50组作为初始种群位置)：
%  首先随机生成初始种群位置
%  然后随机生成初始种群速度
%  然后初始化个体历史最佳位置，以及个体历史最佳适应度
%  然后初始化群体历史最佳位置，以及群体历史最佳适应度
step = 1;
while true
    if step == 1
        data = readmatrix('BLISK4_D123.csv','Range','A2:D2001');
        data_y = readmatrix('BLISK4_D123.csv','Range','E2:E2001');
        % 从中随机抽取 50 组数据，保证其符合正态分布
        num_samples = 50;
        pop_x = datasample(data, num_samples)';  % 初始化种群位置
        y = datasample(data_y, num_samples)';

        pop_v = 0.5 * rands(4,50);  % 初始化种群速度
        gbest = pop_x;

        % 调用MLS拟合函数，核函数带宽 h 可以根据数据的情况调整
        for j=1:sizepop
            fitness_gbest(j) = fun(pop_x(:,j));                      % 每个个体的历史最佳适应度
        end
        mu0 = mean(fitness_gbest);
        sigma0 = std(fitness_gbest);
        step = 2;
    elseif step == 2
        if (mu0-mu)/sqrt(sigma0^2-sigma^2)>R            % 这里mu0一定要比mu大，不行就改mu的值！
            zbest = pop_x(:,1);                         % 种群的历史最佳位置
            fitness_zbest = fitness_gbest(1);           % 种群的历史最佳适应度
            for j=1:sizepop
                if fitness_gbest(j) < fitness_zbest     % 如果求最小值，则为<; 如果求最大值，则为>;  
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
        record_x = zeros(dim, ger);
        record = zeros(ger, 1);          % 记录器
        while iter <= ger
            step = 1;
            while true                   % 以可靠度为约束设置循环
                if step == 1
                    for j=1:sizepop
                        %    更新速度并对速度进行边界处理
                        w = (w1-w2)*tan(0.875*(1-(iter/ger)^k))+w2;
                        pop_v(:,j)= w*pop_v(:,j) + c_1*rand*(gbest(:,j)-pop_x(:,j))+c_2*rand*(zbest-pop_x(:,j));% 速度更新公式
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


%% IPSO优化算法迭代三维图：
X = record_x;
Y = meshgrid(record',record');
[x1,x2] = meshgrid(X(1,:),X(2,:));
figure;
surf(x1,x2,Y);

% [x1,x2] = meshgrid(996:0.1:1036,646:0.1:675);
% Y = ones(size(x1,1),size(x1,2));
% for i = 1:size(x1,1)              % size(x1,1)表示提取x1的行数；size(x1,2)表示提取x1的列数
%     for j = 1:size(x1,2)
%         X = [x1(i,j);x2(i,j);234.6962;490.9244];
%         Y(i,j) = fun(X); 
%     end
% end
% figure;
% surf(x1,x2,Y);
% hold on;
[X,Y] = meshgrid(1:0.5:10,1:20);
Z = sin(X) + cos(Y);
surf(X,Y,Z)

%% 导出优化前后寿命的样本：
fid = fopen('sigma_D_LIFE0.5.txt','wt');      % wt有重新写入的作用
fprintf(fid,'%g\n',Life_x); 
fid = fopen('sigma_D_LIFE0.9.txt','wt');
fprintf(fid,'%g\n',Life_y); 
fid = fopen('sigma_D_LIFE0.95.txt','wt');
fprintf(fid,'%g\n',Life_z); 
fid = fopen('sigma_D_LIFE0.5_optim.txt','wt');
fprintf(fid,'%g\n',Life_Optim_x); 
fid = fopen('sigma_D_LIFE0.9_optim.txt','wt');
fprintf(fid,'%g\n',Life_Optim_y); 
fid = fopen('sigma_D_LIFE0.95_optim.txt','wt');
fprintf(fid,'%g\n',Life_Optim_z); 




%% 优化前的LCF的应变范围①和平均应力②的直方图：
% 叶片B1：
  % ①应变范围（左）、②平均应力（右）
data_life= load('BLISK_D123.csv');
subplot(1,2,1);
histfit(data_life(:,16),100,'norm');
xlabel('Blade strain range');
ylabel('Relative frequency');

subplot(1,2,2);
histfit(data_life(:,12),100,'norm');
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
% 叶片B1：
  % ①应变范围（左）、②平均应力（右） 
subplot(1,2,1);
histfit(x_strain,100,'norm');
xlabel('Blade strain range');
ylabel('Relative frequency');

subplot(1,2,2);
histfit(x_stress,100,'norm');
xlabel('Blade stress range');
ylabel('Relative frequency');

% 盘D：
  % ①应变范围（左）、②平均应力（右）
subplot(1,2,1);
histfit(x_strain,100,'norm');
xlabel('Disk strain range');
ylabel('Relative frequency');

subplot(1,2,2);
histfit(x_stress,100,'norm');
xlabel('Disk stress range');
ylabel('Relative frequency');

%% 优化前寿命样本分布图：

  % 三种置信度的图放一起，叶片B1寿命样本分布图：
B1 = load('sigma_B1_LIFE0.5.txt');B2 = load('sigma_B1_LIFE0.9.txt');B3 = load('sigma_B1_LIFE0.95.txt');
b1 = load('sigma_B1_LIFE0.5_optim.txt');b2 = load('sigma_B1_LIFE0.9_optim.txt');b3 = load('sigma_B1_LIFE0.95_optim.txt');

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