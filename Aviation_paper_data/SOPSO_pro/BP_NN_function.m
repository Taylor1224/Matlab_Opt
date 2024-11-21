function target_data = BP_NN_function(input_data)
% 生成输入输出样本：
    x = readmatrix('BLISK_D123.csv','Range','E1:H2000')';
    y = readmatrix('BLISK_D123.csv','Range','L1:L2000')';

% 对目标值加入噪声
%     n = 0.001 * rand(1,length(x));
%     y = y + n;

% 数据归一化，创建测试数据集(归一化函数：mapminmax★)
    [x_, ps] = mapminmax(x);  % 将输入的数据x进行了最小-最大归一化处理。其中x_是经过归一化处理后的数据，ps是用于后续对其他数据进行相同处理的参数。
    data_input = x_;
    data_target = y;

% 训练神经网络（也可打开nntool，在matlab中使用gui创建神经网络）（前馈神经网络函数：feedforwardnet★）
    hidden_layer_size = 9;
    net = feedforwardnet(hidden_layer_size);
    net = train(net, data_input, data_target);  % ！！！训练之后net就是神经网络模型了！！！

    data_t = mapminmax('apply', input_data, ps);  % 'apply'参数表示对data_test应用mapminmax函数，ps则是用于保存数据缩放参数的结构体。
% 预测结果：
    target_data = sim(net, data_t);  % sim是预测函数，使用神经网络模型net对输入数据data_t进行预测，然后将预测结果存储在data_y中。

end
