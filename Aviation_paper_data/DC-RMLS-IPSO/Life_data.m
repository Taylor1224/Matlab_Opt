%% 应力应变对应的寿命：
N = 100; K_ =1420 ;n_ = 0.05;E = 161000;
% 优化前的叶片应力的均值和方差：
Mu0 = 714.201;
Sigma0 = 17.809;
% 优化前的盘应力的均值和方差：
Mu_0 = 942.799;
Sigma_0 = 22.994;
% 假设优化后的叶片应力的均值和方差：
mu = 690;
sigma = 10;
% 假设优化后的盘应力的均值和方差：
mu_ = 910;
sigma_ = 10;

muXE = 161000; sigmaXE = 3220;
y_optim = lhsdesign(N,1);
x_optim = norminv(y_optim,muXE,sigmaXE);
muX = [Mu0;Mu_0;mu;mu_];
sigmaX = [Sigma0;Sigma_0;sigma;sigma_];
y = lhsdesign(N,4);
x = [norminv(y(:,1),muX(1),sigmaX(1)) , norminv(y(:,2),muX(2),sigmaX(2)), norminv(y(:,3),muX(3),sigmaX(3)), norminv(y(:,4),muX(4),sigmaX(4))];

x_strain1 = ones(N,1);
x_stress1 = x(:,4);
for i = 1:N
    xstrain1 = x_stress1(i,1)/x_optim(i,1)+(x_stress1(i,1)/K_)^(1/n_);
    x_strain1(i,1) = xstrain1;
end

e = load('e.mat').e; m1 = x_strain1; n1 = x_stress1;
x1 = load('x1.mat').x1;x2 = load('x2.mat').x2;x3 = load('x3.mat').x3;
k1 = x1(:,1); b1 = x1(:,2); s1 = x1(:,3); c1 = x1(:,4); % 0.5置信度下的疲劳参数
k2 = x2(:,1); b2 = x2(:,2); s2 = x2(:,3); c2 = x2(:,4); % 0.9置信度下的疲劳参数
k3 = x3(:,1); b3 = x3(:,2); s3 = x3(:,3); c3 = x3(:,4); % 0.95置信度下的疲劳参数
for i=1:N
    Life_x1(i)=fminbnd(@(x)abs(m1(i)/2-(k1(i)-n1(i))/e(i)*(2*x)^b1(i)-s1(i)*(2*x)^c1(i)),100,100000000000000000000);% x为置信度0.5时的寿命
    Life_y1(i)=fminbnd(@(y)abs(m1(i)/2-(k2(i)-n1(i))/e(i)*(2*y)^b2(i)-s2(i)*(2*y)^c2(i)),100,100000000000000000000);% y为置信度0.9时的寿命
    Life_z1(i)=fminbnd(@(z)abs(m1(i)/2-(k3(i)-n1(i))/e(i)*(2*z)^b3(i)-s3(i)*(2*z)^c3(i)),100,100000000000000000000);% z为置信度0.95时的寿命
end