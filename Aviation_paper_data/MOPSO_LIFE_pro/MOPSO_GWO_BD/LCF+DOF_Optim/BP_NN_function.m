function target_data = BP_NN_function(input_data)
% 生成输入输出样本：
    x = readmatrix('BLISK4_D123.csv','Range','A2:D2001')';
    y1 = readmatrix('BLISK4_D123.csv','Range','E2:E2001')'; % E列为B1的应力；G为D的应力；I为B1的径向变形；J为D的径向变形。
    y2 = readmatrix('BLISK4_D123.csv','Range','G2:G2001')';
    y3 = readmatrix('BLISK4_D123.csv','Range','I2:I2001')';
    y4 = readmatrix('BLISK4_D123.csv','Range','J2:J2001')';
% 数据归一化，创建测试数据集(归一化函数：mapminmax★)
    [x_, inputps] = mapminmax(x);  % 将输入的数据x进行了最小-最大归一化处理。其中x_是经过归一化处理后的数据，ps是用于后续对其他数据进行相同处理的参数。
    data_input = x_;
    data_target1 = y1;
    data_target2 = y2;
    data_target3 = y3;
    data_target4 = y4;

% 训练神经网络（也可打开nntool，在matlab中使用gui创建神经网络）（前馈神经网络函数：feedforwardnet★）
    hidden_layer_size = 9;
    net = feedforwardnet(hidden_layer_size);
    net1 = train(net, data_input, data_target1);  % ！！！训练之后net就是神经网络模型了！！！
    net2 = train(net, data_input, data_target2);
    net3 = train(net, data_input, data_target3);
    net4 = train(net, data_input, data_target4);

    data_t = mapminmax('apply', input_data, inputps);  % 'apply'参数表示对data_test应用mapminmax函数，ps则是用于保存数据缩放参数的结构体。
% 预测结果：
    target_data1 = sim(net1, data_t);  % sim是预测函数，使用神经网络模型net对输入数据data_t进行预测，然后将预测结果存储在data_y中。
    target_data2 = sim(net2, data_t);
    target_data3 = sim(net3, data_t);
    target_data4 = sim(net4, data_t);

    target_data = [target_data1 target_data2 target_data3 target_data4]';
end
