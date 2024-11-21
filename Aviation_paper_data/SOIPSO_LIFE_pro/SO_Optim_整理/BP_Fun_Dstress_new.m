function target_data = BP_Fun_Dstress_new(input_data)
% 导入BPNN网络：
    S1 = load('net2_new.mat');
    net2_new = S1.net2_new;
    S2 = load('inputps2.mat');  % inputps1.mat中的结构体名称为'inputps'，应该将这个文件中的'inputps'结构体赋给inputps!
    inputps = S2.inputps;
    % ！！！这里：'load('net1_new.mat') 和 load('inputps1.mat')，需要指定变量名称，否则加载的变量将是结构体格式！！！

    data_t = mapminmax('apply', input_data, inputps);  % 'apply'参数表示对data_test应用mapminmax函数，ps则是用于保存数据缩放参数的结构体。
    % 这里利用mapminmax函数将input_data的数据采用导入的inputps归一化方式进行归一化处理！

% 预测结果：
    target_data2 = sim(net2_new, data_t);  % sim是预测函数，使用神经网络模型net对输入数据data_t进行预测，然后将预测结果存储在data_y中。
    target_data = target_data2';
end