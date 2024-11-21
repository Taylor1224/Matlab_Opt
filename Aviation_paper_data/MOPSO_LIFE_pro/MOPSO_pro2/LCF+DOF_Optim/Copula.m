%% 生成样本分布相关性图：
% 定义边缘分布
mu1 = 0; sigma1 = 1; % 第一个变量的均值和标准差
mu2 = 0; sigma2 = 1; % 第二个变量的均值和标准差

% 生成样本
n = 1000; % 样本大小
X1 = normrnd(mu1, sigma1, [n, 1]); % 第一个变量的样本
X2 = normrnd(mu2, sigma2, [n, 1]); % 第二个变量的样本

% 定义Copula
copula = copularnd('Gaussian', 0.5, n);  % 表示使用相关系数为 0.5 的高斯 Copula 来生成 n 个随机样本。
% 生成相关随机变量
Y = copula(:,1) * sigma1 + mu1; % 相关变量Y
Z = copula(:,2) * sigma2 + mu2; % 相关变量Z

% 可视化
scatter(Y, Z);
title('Copula生成的相关变量');
xlabel('Y变量');
ylabel('Z变量');




% 假设已有数据样本
data = [2.3, 2.9, 3.2, 2.7, 3.0, 2.5, 2.9, 3.1];

% 定义似然函数
likelihoodFunction = @(params) -sum(log(normpdf(data, params(1), params(2))));

% 使用fminsearch函数进行极大似然估计
initialParams = [mean(data), std(data)];
estimatedParams = fminsearch(likelihoodFunction, initialParams);

% 输出估计的参数
mu_estimated = estimatedParams(1);
sigma_estimated = estimatedParams(2);

disp(['估计的均值: ', num2str(mu_estimated)]);
disp(['估计的标准差: ', num2str(sigma_estimated)]);





% 生成样本数据
n = 100;  % 样本大小
mu1 = 5;  % 第一个分布的均值
sigma1 = 2;  % 第一个分布的标准差
mu2 = 10; % 第二个分布的均值
sigma2 = 3;  % 第二个分布的标准差
rho_true = 0.7; % 真实的相关性参数

% 生成协方差矩阵
Sigma = [sigma1^2, rho_true*sigma1*sigma2; rho_true*sigma1*sigma2, sigma2^2];

% 生成二维正态分布样本
data = mvnrnd([mu1, mu2], Sigma, n);

% 初始参数猜测
initial_params = [mean(data(:,1)), mean(data(:,2)), std(data(:,1)), std(data(:,2)), 0];

% 使用 fminsearch 进行极大似然估计
estimated_params = fminsearch(@(params) mle_func(params, data), initial_params);

% 显示估计结果
disp('估计的参数：');
disp(['mu1 = ', num2str(estimated_params(1))]);
disp(['mu2 = ', num2str(estimated_params(2))]);
disp(['sigma1 = ', num2str(estimated_params(3))]);
disp(['sigma2 = ', num2str(estimated_params(4))]);
disp(['rho = ', num2str(estimated_params(5))]);