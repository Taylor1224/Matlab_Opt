clc;
clear;
% new_data.txt中：1是弹性模量；2,3分别是优化前3个置信度的应变范围和平均应力；4,5分别是优化后3个置信度的应变范围和平均应力；
% 6,7,8,9是p=0.5的疲劳参数；10,11,12,13是p=0.9的疲劳参数；14,15,16,17是p=0.95的疲劳参数；

%% 求优化前和优化后的寿命：
data= load('new_data4.txt');
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

fid = fopen('LIFE0.5.txt','wt');      % wt有重新写入的作用
fprintf(fid,'%g\n',Life_x); 
fid = fopen('LIFE0.9.txt','wt');
fprintf(fid,'%g\n',Life_y); 
fid = fopen('LIFE0.95.txt','wt');
fprintf(fid,'%g\n',Life_z); 
fid = fopen('LIFE0.5_optim.txt','wt'); 
fprintf(fid,'%g\n',Life_Optim_x); 
fid = fopen('LIFE0.9_optim.txt','wt');
fprintf(fid,'%g\n',Life_Optim_y); 
fid = fopen('LIFE0.95_optim.txt','wt');
fprintf(fid,'%g\n',Life_Optim_z); 

%histfit(Life_x)
%histfit(Life_x,100,'lognorm');
%histfit(Life_y,100,'lognorm')
%histfit(Life_z,100,'lognorm')
%histfit(Life_Optim_x,100,'lognorm')
%histfit(Life_Optim_y,100,'lognorm')
%histfit(Life_Optim_z,100,'lognorm')

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

%% 将txt文件转为excel文件：
%data = importdata('D:\file\Matlab_file\MATLAB_algorithmic_data\data\PSO_program\Code\Demo1\data.txt');
%Life = 'data.xlsx';
%xlswrite(Life,data);