clc;clear all;

%% Latin超立方抽样参数：
nS = 1e5;nS1 = 5e2;K_ = 2524;n_ = 0.101;Kt = 2.32;

% data= load('new_data.txt');
% e = data(:,1);

stress0_true = ones(nS1,1);
strain0_true = ones(nS1,1); 
stress_true = ones(nS1,1);
strain_true = ones(nS1,1);
stress0_range = ones(nS1,1);
strain0_range = ones(nS1,1);
stress_range = ones(nS1,1);
strain_range = ones(nS1,1);

%% 优化前的真实的应变和应力样本：

% ！！！求解过程中不能用syms直接定义一个矩阵变量，要通过for循环实现！！！
% ！！！求解过程中不能用syms直接定义一个矩阵变量，要通过for循环实现！！！
% ！！！求解过程中不能用syms直接定义一个矩阵变量，要通过for循环实现！！！

% 名义应变/应变范围muX0(1)，名义应力/应力范围muX0(2)：
muX0 = [4.4269444e-3;0.71385e3;161000];
sigmaX0 = [9.2887787e-5;0.13671e2;3220];

y = lhsdesign(nS1,3);
x = [norminv(y(:,1),muX0(1),sigmaX0(1)) , norminv(y(:,2),muX0(2),sigmaX0(2)), norminv(y(:,3),muX0(3),sigmaX0(3))];

% 真实应力stress0_true、真实应变strain0_true：
for m = 1:nS1
    syms stress0true
    stress0true = abs(vpasolve(Kt^2*x(m,1)*x(m,2)-stress0true^2/x(m,3)-stress0true*(stress0true/K_)^(1/n_) == 0,stress0true)); 
    strain0true = stress0true/x(m,3)+(stress0true/K_)^(1/n_);
    stress0_true(m,1) = stress0true;
    strain0_true(m,1) = strain0true;
end

% 真实应力范围stress0_range、真实应变范围strain0_range
for n = 1:nS1
    syms stress0range
    stress0range = abs(vpasolve(Kt^2*x(n,1)*x(n,2)-stress0range^2/x(m,3)-2*stress0range*(stress0range/(2*K_))^(1/n_) == 0,stress0range)); 
    strain0range = stress0range/(2*x(m,3))+(stress0range/(2*K_))^(1/n_);
    stress0_range(n,1) = stress0range;
    strain0_range(n,1) = strain0range;
end

% 优化前的平均应力：
mean_stress0 = stress0_true - stress0_range/2;


%% 优化后的应变和应力样本：
% 线性：
S_mu1 = 0.0042067;S_sigma1 = 3.8456e-05;
% 无交叉项：
S_mu2 = 0.004206;S_sigma2 = 3.7839e-05;
% 有交叉项：
S_mu3 = 0.0042051;S_sigma3 = 3.7024e-05;

% 名义应变/应变范围muX(1)，名义应力/应力范围x_stress：
muX = [S_mu3;161000];
sigmaX = [S_sigma3;3220];
y_optim = lhsdesign(nS1,2);
x_optim = [norminv(y_optim(:,1),muX(1),sigmaX(1)) , norminv(y_optim(:,2),muX(2),sigmaX(2))];

% 已知优化后的名义应力，用循环应力-应变公式求优化后的名义应变x_strain：
x_stress = ones(nS1,1);
for i = 1:nS1
    syms xstress
    xstress = abs(vpasolve(xstress/x_optim(i,2)+(xstress/K_)^(1/n_)-x_optim(i,1) == 0,xstress));
    x_stress(i,1) = xstress;
end

% 用循环求优化后的真实应力、应变：
for j = 1:nS1
    syms stresstrue
    stresstrue = abs(vpasolve(Kt^2*x_optim(j,1)*x_stress(j,1)-stresstrue^2/x_optim(j,2)-stresstrue*(stresstrue/K_)^(1/n_) == 0,stresstrue)); 
    straintrue = stresstrue/x_optim(j,2)+(stresstrue/K_)^(1/n_);
    stress_true(j,1) = stresstrue;
    strain_true(j,1) = straintrue;
end

% 真实应力范围stress0_range、真实应变范围strain0_range
for n = 1:nS1
    syms stressrange
    stressrange = abs(vpasolve(Kt^2*x_optim(n,1)*x_stress(n,1)-stressrange^2/x_optim(n,2)-2*stressrange*(stressrange/(2*K_))^(1/n_) == 0,stressrange)); 
    strainrange = stressrange/(2*x_optim(n,2))+(stressrange/(2*K_))^(1/n_);
    stress_range(n,1) = stressrange;
    strain_range(n,1) = strainrange;
end

% 优化后的平均应力：
mean_stress = stress_true - stress_range/2;


%% 叶片B1的①真实应变范围和②平均应力的分布图：
subplot(1,2,1);
histfit(x_optim(:,1),100,'norm');
xlabel('Blade strain range');
ylabel('Relative frequency');

subplot(1,2,2);
histfit(x_stress,100,'norm');
xlabel('Blade stress range');
ylabel('Relative frequency');

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

%% 将多个工作区数据合并导出到一个txt文件：
a = [x(:,3),strain0_range(:,1),mean_stress0(:,1),strain_range(:,1),mean_stress(:,1),x1(:,:),x2(:,:),x3(:,:)];
fid=fopen('new_data14.txt','wt');%写入文件路径
[m,n]=size(a); 
 for i=1:1:m
    for j=1:1:n
       if j==n
         fprintf(fid,'%g\n',a(i,j));
      else
        fprintf(fid,'%g\t',a(i,j));
       end
    end
end
fclose(fid);


%% matlab求解一元n次方程的方法：
%syms E;
%mu = 4.4269444e-003; 
%E = vpasolve(mu - 7.1384833e2/E + (7.1384833e2/2516)^(1/0.133) == 0 , E)


