clc;clear all;

%% Latin超立方抽样参数：
nS = 1e5;nS1 = 2000;K_ = 1420 ;n_ = 0.08;Kt = 2.32;


%% 优化前的应变和平均应力样本：

% 应变muX0(1)，平均应力muX0(2)：
muX0 = [4.4295e-3;714.201;161000];
sigmaX0 = [1.228e-4;17.809;3220];

y = lhsdesign(nS1,3);
x = [norminv(y(:,1),muX0(1),sigmaX0(1)) , norminv(y(:,2),muX0(2),sigmaX0(2)), norminv(y(:,3),muX0(3),sigmaX0(3))];


%% 优化后的应变和应力样本：

% BPNN模型优化后的均值和方差：
BP_mu = 0.0042051;BP_sigma = 3.7024e-05;

% 名义应变/应变范围muX(1)，名义应力/应力范围x_stress：
muX = [BP_mu;161000];
sigmaX = [BP_sigma;3220];
y_optim = lhsdesign(nS1,2);
x_optim = [norminv(y_optim(:,1),muX(1),sigmaX(1)) , norminv(y_optim(:,2),muX(2),sigmaX(2))];

% 已知优化后的应力，用循环应力-应变公式求优化后的应变x_strain：
x_stress = ones(nS1,1);
for i = 1:nS1
    syms xstress
    xstress = abs(vpasolve(xstress/x_optim(i,2)+(xstress/K_)^(1/n_)-x_optim(i,1) == 0,xstress));
    x_stress(i,1) = xstress;
end

%% 优化后叶片B1的①应变和②平均应力的分布图：
subplot(1,2,1);
histfit(x_optim(:,1),150,'norm');
xlabel('Blade strain range');
ylabel('Relative frequency');

subplot(1,2,2);
histfit(x_stress,150,'norm');
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
a = [x(:,3),x(:,1),x(:,2),x_optim(:,1),x_stress,x1(:,:),x2(:,:),x3(:,:)];
fid=fopen('B1_new_data2.txt','wt');%写入文件路径
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


