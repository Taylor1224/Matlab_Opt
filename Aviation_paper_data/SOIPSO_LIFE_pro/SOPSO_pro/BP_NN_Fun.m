function target_data = BP_NN_Fun(input_data)
% 导入BPNN网络：
    x = readmatrix('BLISK4_D123.csv','Range','A2:D2001')';
    load('net1.mat');
    load('net2.mat');
    load('net3.mat');
    load('net4.mat');
    [x_, ps] = mapminmax(x);

    data_t = mapminmax('apply', input_data, ps);  % 'apply'参数表示对data_test应用mapminmax函数，ps则是用于保存数据缩放参数的结构体。
% 预测结果：
    target_data1 = sim(net1, data_t);  % sim是预测函数，使用神经网络模型net对输入数据data_t进行预测，然后将预测结果存储在data_y中。
    target_data2 = sim(net2, data_t);
    target_data3 = sim(net3, data_t);
    target_data4 = sim(net4, data_t);
    target_data = [target_data1 target_data2 target_data3 target_data4]';
end