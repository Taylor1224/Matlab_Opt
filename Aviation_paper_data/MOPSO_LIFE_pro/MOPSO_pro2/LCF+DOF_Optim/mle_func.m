% 极大似然估计函数
function negLogLikelihood = mle_func(params, data)
    mu1 = params(1);
    mu2 = params(2);
    sigma1 = params(3);
    sigma2 = params(4);
    rho = params(5);
    
    % 协方差矩阵
    Sigma = [sigma1^2, rho*sigma1*sigma2; rho*sigma1*sigma2, sigma2^2];
    
    % 计算样本点的对数似然值
    negLogLikelihood = -sum(log(mvnpdf(data, [mu1, mu2], Sigma)));
end