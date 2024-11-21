% 导入数据
data = load('B_point_Before.csv')'; % 假设数据保存在变量data中
data1 = load('B_point_After.csv')';
% 绘制直方图
x = [27382:27416]';
h1 = bar(x,data);
hold on;
h2 = bar(x,data1);
ylim([0.95,1])
title('直方图');
xlabel('数据值');
ylabel('DOF');

