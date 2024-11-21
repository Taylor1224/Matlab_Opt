%% 
clc
clear

%% 用自带的particleswarm函数进行优化：
% 1. 定义目标函数
       % 在Function_1中
y = @B1_stress_have_Xterms;

m = [5.17798e-2,0,0,0;
    0,6.89986e-2,0,0;
    0,0,2.10794e-1,0;
    0,0,0,1.09157e-1];
n = [5.26337e1,4.55018e1,4.84582e1,5.45754e1];% B1/D无交叉项的与有交叉项的一样；
o = [5.17798e-2,6.89986e-2,2.10794e-1,1.09157e-1];

% 2. 设置PSO算法参数
options = optimoptions('particleswarm', 'OutputFcn', @(x, optimValues, state) output(x, optimValues, state),'SwarmSize', 100, 'MaxIterations', 100);

% 3. 约束条件
  % ① 不等式约束：
% muY0 = [5.7600917e-3;4.4269444e-3];
% sigmaY0 = [1.3988648e-4;9.2887787e-5];
%mu=5.7600917e-3;     % 能用矩阵吗？
%sigma=1.3988648e-4;
%p=0.9;
%y=norminv(p,mu,sigma);

  % ② 上下界约束：

lb = [996, 646, 225, 490]*m-n;% 变量下界
ub = [1036, 675, 235, 510]*m-n;% 变量上界
[x_opt, f_opt] = particleswarm(y, 4, lb, ub, options);

% 4. 输出优化结果
disp('优化结果：');
disp(['最优虚拟变量值：', num2str(x_opt)]);
X = (x_opt+n)./o;
disp(['最优真实变量值：', num2str(X)]);
disp(['最优函数值：', num2str(f_opt)]);
%%
%画图

