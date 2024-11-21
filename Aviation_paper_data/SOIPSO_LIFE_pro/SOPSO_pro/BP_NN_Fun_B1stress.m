function target_data = BP_NN_Fun_B1stress(input_data)
% 导入BPNN网络：
    x = readmatrix('BLISK4_D123.csv','Range','A2:D2001')';
    load('net1.mat');
    [x_, ps] = mapminmax(x);

    data_t = mapminmax('apply', input_data, ps);  % 'apply'参数表示对data_test应用mapminmax函数，ps则是用于保存数据缩放参数的结构体。
% 预测结果：
    target_data1 = sim(net1, data_t);  % sim是预测函数，使用神经网络模型net对输入数据data_t进行预测，然后将预测结果存储在data_y中。
    target_data = target_data1';
end