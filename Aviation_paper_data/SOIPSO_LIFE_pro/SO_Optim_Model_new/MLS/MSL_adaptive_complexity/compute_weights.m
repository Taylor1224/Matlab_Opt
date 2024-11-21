function weights = compute_weights(xi, x, h)
    % 使用高斯加权函数
    weights = exp(-((x - xi).^2) / (2*h^2));
end