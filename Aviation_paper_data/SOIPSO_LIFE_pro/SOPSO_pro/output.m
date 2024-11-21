function stop = output(x, optimValues, state)
    stop = false;  % 初始化停止标志

    % 在迭代外部初始化count
    persistent count;
    if isempty(count)
        count = 0;
    end
    if isequal(state, 'iter')  % 每次迭代时执行——如果是迭代状态就执行以下操作
        % 获取当前迭代的输出值
        currentOutputValue = optimValues.fval;

        % 在此处进行输出值小于某个值的概率计算
        % 你可以将当前值与所需的阈值进行比较，并计算概率

        % 例如，计算小于阈值的次数
        
        threshold = 750;
        if currentOutputValue < threshold
            count = count + 1;
        end
     
        num_iterations = optimValues.iteration;  % 获取当前迭代次数
        probability = (count / num_iterations) * 100;  % 计算百分化的概率
        fprintf('小于阈值概率: %f%%\n', probability);
    end
end

%  创建一个自定义输出函数