% Huber函数定义
function w = huber_weight(r, delta)
    % 计算每个残差的Huber权重
    w = zeros(size(r));
    for i = 1:length(r)
        if abs(r(i)) <= delta
            w(i) = 0.5 * r(i)^2; % 二次损失
        else
            w(i) = delta * (abs(r(i)) - 0.5 * delta); % 线性损失
        end
    end
end

