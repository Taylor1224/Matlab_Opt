function target_data = BP_NN_Fun(input_data)
% 导入BPNN网络： 
    S1 = load('net1_new.mat');      % BP神经网络是net1 ; IPSO-BP神经网络是net1_new
    net1_new = S1.net1_new;
    S2 = load('net2_new.mat');
    net2_new = S2.net2_new;
    S3 = load('net3_new.mat');
    net3_new = S3.net3_new;
    S4 = load('net4_new.mat');
    net4_new = S4.net4_new;
    T1 = load('inputps1.mat');  
    inputps1 = T1.inputps;
    T2 = load('inputps2.mat');  
    inputps2 = T2.inputps;
    T3 = load('inputps3.mat');  
    inputps3 = T3.inputps;
    T4 = load('inputps4.mat');  
    inputps4 = T4.inputps;


    data_t1 = mapminmax('apply', input_data(1,1), inputps1);  % 'apply'参数表示对data_test应用mapminmax函数，inputps1则是用于保存数据缩放参数的结构体。
    data_t2 = mapminmax('apply', input_data(2,1), inputps2);
    data_t3 = mapminmax('apply', input_data(3,1), inputps3);
    data_t4 = mapminmax('apply', input_data(4,1), inputps4);
    data_t = [data_t1(1,1),data_t2(2,1),data_t3(3,1),data_t4(4,1)]';     %  这里因为net1_new需要的对应data_t数据要是4行1列的！！！
% 预测结果：
    target_data1 = sim(net1_new, data_t);  % sim是预测函数，使用神经网络模型net对输入数据data_t进行预测，然后将预测结果存储在data_y中。
    target_data2 = sim(net2_new, data_t); 
    target_data3 = sim(net3_new, data_t);
    target_data4 = sim(net4_new, data_t);
    target_data = [target_data1() target_data2 target_data3 target_data4]';
end