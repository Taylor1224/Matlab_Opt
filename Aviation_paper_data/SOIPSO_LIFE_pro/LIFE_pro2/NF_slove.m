clc;
clear;
% new_data.txt中：1是弹性模量；2,3分别是优化前3个置信度的应变范围和平均应力；4,5分别是优化后3个置信度的应变范围和平均应力；
% 6,7,8,9是p=0.5的疲劳参数；10,11,12,13是p=0.9的疲劳参数；14,15,16,17是p=0.95的疲劳参数；

%% 求优化前和优化后的寿命：
data= load('B1_new_data2.txt');
e = data(:,1);
m0 = data(:,2);
n0 = data(:,3);
m = data(:,4);
n = data(:,5);
k1 = data(:,6);
b1 = data(:,7);
s1 = data(:,8);
c1 = data(:,9);
k2 = data(:,10);
b2 = data(:,11);
s2 = data(:,12);
c2 = data(:,13);
k3 = data(:,14);
b3 = data(:,15);
s3 = data(:,16);
c3 = data(:,17);

for i=1:length(m0)
    % 优化前：
    Life_x(i)=fminbnd(@(x)abs(m0(i)/2-(k1(i)-n0(i))/e(i)*(2*x)^b1(i)-s1(i)*(2*x)^c1(i)),100,100000000000000000000);% x为置信度0.5时的寿命
    Life_y(i)=fminbnd(@(y)abs(m0(i)/2-(k2(i)-n0(i))/e(i)*(2*y)^b2(i)-s2(i)*(2*y)^c2(i)),100,100000000000000000000);% y为置信度0.9时的寿命
    Life_z(i)=fminbnd(@(z)abs(m0(i)/2-(k3(i)-n0(i))/e(i)*(2*z)^b3(i)-s3(i)*(2*z)^c3(i)),100,100000000000000000000);% z为置信度0.95时的寿命
    
%     Life_x(i)=fminbnd(@(x)abs(m0(i)-(k1(i)-n0(i))/e(i)*(2*x)^b1(i)-s1(i)*(2*x)^c1(i)),100,100000000000000000000);% x为置信度0.5时的寿命
%     Life_y(i)=fminbnd(@(y)abs(m0(i)-(k2(i)-n0(i))/e(i)*(2*y)^b2(i)-s2(i)*(2*y)^c2(i)),100,100000000000000000000);% y为置信度0.9时的寿命
%     Life_z(i)=fminbnd(@(z)abs(m0(i)-(k3(i)-n0(i))/e(i)*(2*z)^b3(i)-s3(i)*(2*z)^c3(i)),100,100000000000000000000);% z为置信度0.95时的寿命
    % 优化后：
    Life_Optim_x(i)=fminbnd(@(x)abs(m(i)/2-(k1(i)-n(i))/e(i)*(2*x)^b1(i)-s1(i)*(2*x)^c1(i)),100,100000000000000000000);% x为置信度0.5时的寿命
    Life_Optim_y(i)=fminbnd(@(y)abs(m(i)/2-(k2(i)-n(i))/e(i)*(2*y)^b2(i)-s2(i)*(2*y)^c2(i)),100,100000000000000000000);% y为置信度0.9时的寿命
    Life_Optim_z(i)=fminbnd(@(z)abs(m(i)/2-(k3(i)-n(i))/e(i)*(2*z)^b3(i)-s3(i)*(2*z)^c3(i)),100,100000000000000000000);% z为置信度0.95时的寿命

%     Life_Optim_x(i)=fminbnd(@(x)abs(m(i)-(k1(i)-n(i))/e(i)*(2*x)^b1(i)-s1(i)*(2*x)^c1(i)),100,100000000000000000000);% x为置信度0.5时的寿命
%     Life_Optim_y(i)=fminbnd(@(y)abs(m(i)-(k2(i)-n(i))/e(i)*(2*y)^b2(i)-s2(i)*(2*y)^c2(i)),100,100000000000000000000);% y为置信度0.9时的寿命
%     Life_Optim_z(i)=fminbnd(@(z)abs(m(i)-(k3(i)-n(i))/e(i)*(2*z)^b3(i)-s3(i)*(2*z)^c3(i)),100,100000000000000000000);% z为置信度0.95时的寿命
    % 求出来的是data中的每一组数据代入公式
    % x = fminbnd(fun,x1,x2) 返回一个值 x，该值是 fun 中描述的标量值函数在区间 x1 < x < x2 中的局部最小值。
    % 得到的结果就是使fun取最小值的x值,这里有abs，所以fun取的最小值是0.
end

fid = fopen('B1_LIFE0.5.txt','wt');      % wt有重新写入的作用
fprintf(fid,'%g\n',Life_x); 
fid = fopen('B1_LIFE0.9.txt','wt');
fprintf(fid,'%g\n',Life_y); 
fid = fopen('B1_LIFE0.95.txt','wt');
fprintf(fid,'%g\n',Life_z); 
fid = fopen('B1_LIFE0.5_optim.txt','wt'); 
fprintf(fid,'%g\n',Life_Optim_x); 
fid = fopen('B1_LIFE0.9_optim.txt','wt');
fprintf(fid,'%g\n',Life_Optim_y); 
fid = fopen('B1_LIFE0.95_optim.txt','wt');
fprintf(fid,'%g\n',Life_Optim_z); 

%% 优化前寿命样本分布图：

  % 三种置信度的图放一起，叶片B1寿命样本分布图：
B1 = load('B1_LIFE0.5.txt');B2 = load('B1_LIFE0.9.txt');B3 = load('B1_LIFE0.95.txt');
b1 = load('B1_LIFE0.5_optim.txt');b2 = load('B1_LIFE0.9_optim.txt');b3 = load('B1_LIFE0.95_optim.txt');

subplot(1,2,1);
histfit(B1(:,1),500,'lognorm');
hold on;
histfit(B2(:,1),140,'lognorm');
hold on;
histfit(B3(:,1),80,'lognorm');
xlabel('Blade LCF life before optimization');
ylabel('Relative frequency');
legend('Blade LCF life sample with P=0.5','Fitted lognormal curve with P=0.5','Blade LCF life sample with P=0.9','Fitted lognormal curve with P=0.9','Blade LCF life sample with P=0.95','Fitted lognormal curve with P=0.95');

subplot(1,2,2);
histfit(b1(:,1),150,'lognorm');
hold on;
histfit(b2(:,1),150,'lognorm');
hold on;
histfit(b3(:,1),150,'lognorm');
xlabel('Blade LCF life after optimization');
ylabel('Relative frequency');
legend('Blade LCF life sample with P=0.5','Fitted lognormal curve with P=0.5','Blade LCF life sample with P=0.9','Fitted lognormal curve with P=0.9','Blade LCF life sample with P=0.95','Fitted lognormal curve with P=0.95');


  % 三种置信度的图放一起，盘D寿命样本分布图：
D1 = load('D_LIFE0.5.txt');D2 = load('D_LIFE0.9.txt');D3 = load('D_LIFE0.95.txt');
d1 = load('D_LIFE0.5_optim.txt');d2 = load('D_LIFE0.9_optim.txt');d3 = load('D_LIFE0.95_optim.txt');

subplot(1,2,1);
histfit(D1(:,1),80,'lognorm');
hold on;
histfit(D2(:,1),80,'lognorm');
hold on;
histfit(D3(:,1),80,'lognorm');
xlabel('Disk LCF life before optimization');
ylabel('Relative frequency');
legend('Disk LCF life sample with P=0.5','Fitted lognormal curve with P=0.5','Disk LCF life sample with P=0.9','Fitted lognormal curve with P=0.9','Disk LCF life sample with P=0.95','Fitted lognormal curve with P=0.95');

subplot(1,2,2);
histfit(d1(:,1),80,'lognorm');
hold on;
histfit(d2(:,1),80,'lognorm');
hold on;
histfit(d3(:,1),80,'lognorm');
xlabel('Disk LCF life after optimization');
ylabel('Relative frequency');
legend('Disk LCF life sample with P=0.5','Fitted lognormal curve with P=0.5','Disk LCF life sample with P=0.9','Fitted lognormal curve with P=0.9','Disk LCF life sample with P=0.95','Fitted lognormal curve with P=0.95');




%% 优化后寿命样本分布图：

  % 置信度0.5下，叶片B1和盘D寿命样本分布图：
b1 = load('B1_LIFE0.5_optim.txt');
subplot(1,2,1);
histfit(b1(:,1),150,'lognorm');
xlabel('Blade life');
ylabel('Relative frequency');

d1 = load('D_LIFE0.5_optim.txt');
subplot(1,2,2);
histfit(d1(:,1),150,'lognorm');
xlabel('Disk life');
ylabel('Relative frequency');

  % 置信度0.9下，叶片B1和盘D寿命样本分布图：
b2 = load('B1_LIFE0.9_optim.txt');
subplot(1,2,1);
histfit(b2(:,1),150,'lognorm');
xlabel('Blade life');
ylabel('Relative frequency');

d2 = load('D_LIFE0.9_optim.txt');
subplot(1,2,2);
histfit(d2(:,1),150,'lognorm');
xlabel('Disk life');
ylabel('Relative frequency');

  % 置信度0.95下，叶片B1和盘D寿命样本分布图：
b3 = load('B1_LIFE0.95_optim.txt');
subplot(1,2,1);
histfit(b3(:,1),150,'lognorm');
xlabel('Blade life');
ylabel('Relative frequency');

d3 = load('D_LIFE0.95_optim.txt');
subplot(1,2,2);
histfit(d3(:,1),150,'lognorm');
xlabel('Disk life');
ylabel('Relative frequency');


%% 优化前后寿命对比图：
  % 置信度0.5下，叶片B1和盘D寿命样本分布图：
B1 = load('B1_LIFE0.5.txt');b1 = load('B1_LIFE0.5_optim.txt');
subplot(1,2,1);
histfit(B1(:,1),150,'lognorm');
hold on ;
histfit(b1(:,1),150,'lognorm');
xlabel('Blade life');
ylabel('Relative frequency');

D1 = load('D_LIFE0.5.txt');d1 = load('D_LIFE0.5_optim.txt');
subplot(1,2,2);
histfit(D1(:,1),150,'lognorm');
hold on ;
histfit(d1(:,1),150,'lognorm');
xlabel('Disk life');
ylabel('Relative frequency');

  % 置信度0.9下，叶片B1和盘D寿命样本分布图：
B2 = load('B1_LIFE0.9.txt');b2 = load('B1_LIFE0.9_optim.txt');
subplot(1,2,1);
histfit(B2(:,1),150,'lognorm');
hold on ;
histfit(b2(:,1),150,'lognorm');
xlabel('Blade life');
ylabel('Relative frequency');
% dfittool工具：
% dfittool(B2(:,1));
% dfittool(b2(:,1));

D2 = load('D_LIFE0.9.txt');d2 = load('D_LIFE0.9_optim.txt');
subplot(1,2,2);
histfit(D2(:,1),150,'lognorm');
hold on ;
histfit(d2(:,1),150,'lognorm');
xlabel('Disk life');
ylabel('Relative frequency');

  % 置信度0.5下，叶片B1和盘D寿命样本分布图：
B3 = load('B1_LIFE0.95.txt');b3 = load('B1_LIFE0.95_optim.txt');
subplot(1,2,1);
histfit(B3(:,1),150,'lognorm');
hold on ;
histfit(b3(:,1),150,'lognorm');
xlabel('Blade life');
ylabel('Relative frequency');

D3 = load('D_LIFE0.95.txt');d3 = load('D_LIFE0.95_optim.txt');
subplot(1,2,2);
histfit(D3(:,1),150,'lognorm');
hold on ;
histfit(d3(:,1),150,'lognorm');
xlabel('Disk life');
ylabel('Relative frequency');


%% 优化前后寿命对比图(叶盘曲线在一张图)：
  % 置信度0.5下，叶片B1和盘D寿命样本分布图：
B1 = load('B1_LIFE0.5.txt');b1 = load('B1_LIFE0.5_optim.txt');
D1 = load('D_LIFE0.5.txt');d1 = load('D_LIFE0.5_optim.txt');
histfit(B1(:,1),150,'lognorm');
hold on ;
histfit(b1(:,1),150,'lognorm');
hold on ;
histfit(D1(:,1),150,'lognorm');
hold on ;
histfit(d1(:,1),150,'lognorm');
xlabel('Blade-Disk life');
ylabel('Relative frequency');

  % 置信度0.9下，叶片B1和盘D寿命样本分布图：
B2 = load('B1_LIFE0.9.txt');b2 = load('B1_LIFE0.9_optim.txt');
D2 = load('D_LIFE0.9.txt');d2 = load('D_LIFE0.9_optim.txt');
histfit(B2(:,1),150,'lognorm');
hold on ;
histfit(b2(:,1),150,'lognorm');
hold on ;
histfit(D2(:,1),150,'lognorm');
hold on ;
histfit(d2(:,1),150,'lognorm');
xlabel('Blade-Disk life');
ylabel('Relative frequency');

  % 置信度0.95下，叶片B1和盘D寿命样本分布图：
B3 = load('B1_LIFE0.95.txt');b3 = load('B1_LIFE0.95_optim.txt');
D3 = load('D_LIFE0.95.txt');d3 = load('D_LIFE0.95_optim.txt');
histfit(B3(:,1),150,'lognorm');
hold on ;
histfit(b3(:,1),150,'lognorm');
hold on ;
histfit(D3(:,1),150,'lognorm');
hold on ;
histfit(d3(:,1),150,'lognorm');
xlabel('Blade-Disk life');
ylabel('Relative frequency');


%% 优化前后寿命对比图(histogram)：
  % 置信度0.5下，叶片B1和盘D寿命样本分布图：
B1 = load('B1_LIFE0.5.txt');b1 = load('B1_LIFE0.5_optim.txt');
subplot(1,2,1);
histogram(B1(:,1),200);
hold on ;
histogram(b1(:,1),200);
xlabel('Blade life');
ylabel('Relative frequency');

D1 = load('D_LIFE0.5.txt');d1 = load('D_LIFE0.5_optim.txt');
subplot(1,2,2);
histogram(D1(:,1),200);
hold on ;
histogram(d1(:,1),200);
xlabel('Disk life');
ylabel('Relative frequency');

  % 置信度0.9下，叶片B1和盘D寿命样本分布图：
B2 = load('B1_LIFE0.9.txt');b2 = load('B1_LIFE0.9_optim.txt');
subplot(1,2,1);
histfit(B2(:,1),150,'lognorm');
hold on ;
histfit(b2(:,1),150,'lognorm');
xlabel('Blade life');
ylabel('Relative frequency');
% dfittool工具：
% dfittool(B2(:,1));
% dfittool(b2(:,1));

D2 = load('D_LIFE0.9.txt');d2 = load('D_LIFE0.9_optim.txt');
subplot(1,2,2);
histfit(D2(:,1),150,'lognorm');
hold on ;
histfit(d2(:,1),150,'lognorm');
xlabel('Disk life');
ylabel('Relative frequency');

  % 置信度0.5下，叶片B1和盘D寿命样本分布图：
B3 = load('B1_LIFE0.95.txt');b3 = load('B1_LIFE0.95_optim.txt');
subplot(1,2,1);
histfit(B3(:,1),150,'lognorm');
hold on ;
histfit(b3(:,1),150,'lognorm');
xlabel('Blade life');
ylabel('Relative frequency');

D3 = load('D_LIFE0.95.txt');d3 = load('D_LIFE0.95_optim.txt');
subplot(1,2,2);
histfit(D3(:,1),150,'lognorm');
hold on ;
histfit(d3(:,1),150,'lognorm');
xlabel('Disk life');
ylabel('Relative frequency');


%% 优化前后寿命仿真历史图：
  % 置信度0.5下，叶片B1和盘D寿命样本分布图：
B1 = load('B1_LIFE0.5.txt');b1 = load('B1_LIFE0.5_optim.txt');
subplot(2,1,1);
plot(B1);
hold on ;
subplot(2,1,2);
plot(b1);

D1 = load('D_LIFE0.5.txt');d1 = load('D_LIFE0.5_optim.txt');
subplot(2,1,1);
plot(D1);
hold on ;
subplot(2,1,2);
plot(d1);
  % 置信度0.9下，叶片B1和盘D寿命样本分布图：
B2 = load('B1_LIFE0.9.txt');b2 = load('B1_LIFE0.9_optim.txt');
subplot(2,1,1);
plot(B2);
hold on ;
subplot(2,1,2);
plot(b2);

D2 = load('D_LIFE0.9.txt');d2 = load('D_LIFE0.9_optim.txt');
subplot(2,1,1);
plot(D2);
hold on ;
subplot(2,1,2);
plot(d2);

  % 置信度0.95下，叶片B1和盘D寿命样本分布图：
B3 = load('B1_LIFE0.9.txt');b3 = load('B1_LIFE0.9_optim.txt');
subplot(2,1,1);
plot(B3);
hold on ;
subplot(2,1,2);
plot(b3);

D3 = load('D_LIFE0.95.txt');d3 = load('D_LIFE0.95_optim.txt');
subplot(2,1,1);
plot(D3);
hold on ;
subplot(2,1,2);
plot(d3);


%% 优化前后寿命仿真历史图和分布直方图（放一起）：
  % 置信度0.5下，叶片B1和盘D寿命样本分布图：
D1 = load('D_LIFE0.5.txt');d1 = load('D_LIFE0.5_optim.txt');
subplot(2,2,1);
plot(D1);
hold on ;
subplot(2,2,2);
histfit(D1(:,1),200,'lognorm');
hold on ;
subplot(2,2,3);
plot(d1);
hold on;
subplot(2,2,4);
histfit(d1(:,1),200,'lognorm');

  % 置信度0.9下，叶片B1和盘D寿命样本分布图：
D2 = load('D_LIFE0.9.txt');d2 = load('D_LIFE0.9_optim.txt');
subplot(2,2,1);
plot(D2);
hold on ;
subplot(2,2,2);
histfit(D2(:,1),200,'lognorm');
hold on ;
subplot(2,2,3);
plot(d2);
hold on;
subplot(2,2,4);
histfit(d2(:,1),200,'lognorm');

  % 置信度0.95下，叶片B1和盘D寿命样本分布图：
D3 = load('D_LIFE0.95.txt');d3 = load('D_LIFE0.95_optim.txt');
subplot(2,2,1);
plot(D3);
hold on ;
subplot(2,2,2);
histfit(D3(:,1),200,'lognorm');
hold on ;
subplot(2,2,3);
plot(d3);
hold on;
subplot(2,2,4);
histfit(d3(:,1),200,'lognorm');










%% 求损伤：
n1 = 2000;n2 = 3500;n3 = 5000;n4 = 6500;n5 = 8000;n6 = 950;n7 = 11000;n8 = 12500;n9 = 14000;n10 = 16500;

D1 = n1./Life_x;
D2 = n1./Life_y;
D3 = n1./Life_z;
D4 = n1./Life_Optim_x;
D5 = n1./Life_Optim_y;
D6 = n1./Life_Optim_z;

qqplot(D1);
dfittool(D1);
histfit(D1,100,'lognorm');


%% 求可靠度：
% 极限状态函数：
Z1 = ones(1,length(D1))-D1;
Z2 = ones(1,length(D2))-D2;
Z3 = ones(1,length(D3))-D3;
Z4 = ones(1,length(D4))-D4;
Z5 = ones(1,length(D5))-D5;
Z6 = ones(1,length(D6))-D6;
% 可靠度：
R1 = length(Z1(Z1>0))/length(Z1);
R2 = length(Z2(Z2>0))/length(Z2);
R3 = length(Z3(Z3>0))/length(Z3);
R4 = length(Z4(Z4>0))/length(Z4);
R5 = length(Z5(Z5>0))/length(Z5);
R6 = length(Z6(Z6>0))/length(Z6);


%% matlab怎么统计矩阵中大于或小于某数的元素个数：
    % 方法一：使用sum函数：
        % sum(sum(a>3))
    
    % 方法二：使用length函数：
        % length(a(a>3))

%% 将txt文件转为excel文件（一列）：
data1 = importdata('B1_LIFE0.5.txt');
data2 = importdata('B1_LIFE0.9.txt');
data3 = importdata('B1_LIFE0.95.txt');
data4 = importdata('B1_LIFE0.5_optim.txt');
data5 = importdata('B1_LIFE0.9_optim.txt');
data6 = importdata('B1_LIFE0.95_optim.txt');
data7 = importdata('D_LIFE0.5.txt');
data8 = importdata('D_LIFE0.9.txt');
data9 = importdata('D_LIFE0.95.txt');
data10 = importdata('D_LIFE0.5_optim.txt');
data11 = importdata('D_LIFE0.9_optim.txt');
data12 = importdata('D_LIFE0.95_optim.txt');
Life1 = 'B1_LIFE0.5.xlsx';
Life2 = 'B1_LIFE0.9.xlsx';
Life3 = 'B1_LIFE0.95.xlsx';
Life4 = 'B1_LIFE0.5_optim.xlsx';
Life5 = 'B1_LIFE0.9_optim.xlsx';
Life6 = 'B1_LIFE0.95_optim.xlsx';
Life7 = 'D_LIFE0.5.xlsx';h
Life8 = 'D_LIFE0.9.xlsx';
Life9 = 'D_LIFE0.95.xlsx';
Life10 = 'D_LIFE0.5_optim.xlsx';
Life11 = 'D_LIFE0.9_optim.xlsx';
Life12 = 'D_LIFE0.95_optim.xlsx';
xlswrite(Life1,data1);
xlswrite(Life2,data2);
xlswrite(Life3,data3);
xlswrite(Life4,data4);
xlswrite(Life5,data5);
xlswrite(Life6,data6);
xlswrite(Life7,data7);
xlswrite(Life8,data8);
xlswrite(Life9,data9);
xlswrite(Life10,data10);
xlswrite(Life11,data11);
xlswrite(Life12,data12);


%% 将txt文件转为excel文件（一行）：
data1 = importdata('B1_LIFE0.5.txt');
data2 = importdata('B1_LIFE0.9.txt');
data3 = importdata('B1_LIFE0.95.txt');
data4 = importdata('B1_LIFE0.5_optim.txt');
data5 = importdata('B1_LIFE0.9_optim.txt');
data6 = importdata('B1_LIFE0.95_optim.txt');
data7 = importdata('D_LIFE0.5.txt');
data8 = importdata('D_LIFE0.9.txt');
data9 = importdata('D_LIFE0.95.txt');
data10 = importdata('D_LIFE0.5_optim.txt');
data11 = importdata('D_LIFE0.9_optim.txt');
data12 = importdata('D_LIFE0.95_optim.txt');
Life1 = 'B1__LIFE0.5.csv';
Life2 = 'B1__LIFE0.9.csv';
Life3 = 'B1__LIFE0.95.csv';
Life4 = 'B1__LIFE0.5_optim.csv';
Life5 = 'B1__LIFE0.9_optim.csv';
Life6 = 'B1__LIFE0.95_optim.csv';
Life7 = 'D__LIFE0.5.csv';
Life8 = 'D__LIFE0.9.csv';
Life9 = 'D__LIFE0.95.csv';
Life10 = 'D__LIFE0.5_optim.csv';
Life11 = 'D__LIFE0.9_optim.csv';
Life12 = 'D__LIFE0.95_optim.csv';
csvwrite(Life1,data1');
csvwrite(Life2,data2');
csvwrite(Life3,data3');
csvwrite(Life4,data4');
csvwrite(Life5,data5');
csvwrite(Life6,data6');
csvwrite(Life7,data7');
csvwrite(Life8,data8');
csvwrite(Life9,data9');
csvwrite(Life10,data10');
csvwrite(Life11,data11');
csvwrite(Life12,data12');


w = importdata('PRNSOL.txt');
L = 'PRNSOL.csv';
csvwrite(L,w);

%% 画概率密度函数曲线图：

dfittool(data1)
% 用dfittool得到的图————>点击新建拟合————>将分布改为对应分布————>应用————>评估————>改小一点x的中间值并应用————>保存到工作区（命名：evaluateresults）

plot(evaluateresults(:,1),evaluateresults(:,2),'LineWidth',2);

hold on;

plot(evaluateresults1(:,1),evaluateresults1(:,2),'LineWidth',2);


title('寿命概率密度曲线');
xlabel('寿命Lf');
ylabel('概率密度Pr');
