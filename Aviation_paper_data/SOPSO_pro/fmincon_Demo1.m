%% fmincon优化算法：寻找有约束多元函数的全局最优解

% fmincon 提供了不同的优化算法，包括 'interior-point'、'sqp' 等。
% 可以通过参数 options 定制优化过程，设置最大迭代次数、优化容忍度等参数，以满足你的需求。
clc
clear

%% 导入函数及定义变量：

m = [5.17798e-2,0,0,0;
    0,6.89986e-2,0,0;
    0,0,2.10794e-1,0;
    0,0,0,1.09157e-1];
n = [5.26337e1,4.55018e1,4.84582e1,5.45754e1];% B1/D无交叉项的与有交叉项的一样；
o = [5.17798e-2,6.89986e-2,2.10794e-1,1.09157e-1];

% 导入目标函数：
fun = @B1_stress_have_Xterms;
% 初始化设计变量：
x0 = [1010, 665, 230,495]*m-n;

%% 设置约束：
% 定义线性不等式约束
A = [];
b = [];

nonlcon = [];

% 定义线性等式约束
Aeq = [];
beq = [];

% 定义设计变量的上下界
lb = [996, 646, 225, 490]*m-n;
ub = [1036, 675, 235, 510]*m-n;

%% 设置选项，启用迭代输出
% 调用 fmincon 寻找全局最小值

options = optimoptions('fmincon', 'OutputFcn', @output,'Display','iter','Algorithm','sqp')
[x, fval] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options)

% [x, fval] = fmincon(fun1, [1010, 665, 230, 495]*m-n, Aeq, beq, lb, ub);

% 打印优化结果:
disp(['最优虚拟变量值：', num2str(x)]);
X = (x+n)./o;
disp(['最优真实变量值：', num2str(X)]);
disp(['最优函数值：', num2str(fval)]);




%  https://ww2.mathworks.cn/help/optim/ug/fmincon.html#busow0w-1——————搜索的网址