clc;clear all;

%% Latin����������������
nS = 1e5;nS1 = 1e4;K_ = 1313;n_ = 0.08;Kt = 2.32;

data= load('new_data.txt');
e = data(:,1);


%% �Ż�ǰ����ʵ��Ӧ���Ӧ��������

% �������������в�����symsֱ�Ӷ���һ�����������Ҫͨ��forѭ��ʵ�֣�����
% �������������в�����symsֱ�Ӷ���һ�����������Ҫͨ��forѭ��ʵ�֣�����
% �������������в�����symsֱ�Ӷ���һ�����������Ҫͨ��forѭ��ʵ�֣�����

% ��APDL�õ�������ֱ����Ϊ��ʵӦ��/Ӧ�䷶ΧmuX0(1)����ʵӦ��/Ӧ����ΧmuX0(2)��
muX0 = [4.4269444e-3;0.71385e3];
sigmaX0 = [9.2887787e-5;0.13671e2];

y = lhsdesign(nS1,2);
x = [norminv(y(:,1),muX0(1),sigmaX0(1)) , norminv(y(:,2),muX0(2),sigmaX0(2))];

% �Ż�ǰ��ƽ��Ӧ����
mean_stress0 = x(:,2) - x(:,2)/2;

%histfit(x(:,1));
%histfit(x(:,2));

%% �Ż����Ӧ���Ӧ��������

% ��ʵӦ��/Ӧ�䷶ΧmuX(1)����ʵӦ��/Ӧ����Χx_stress��
muX = [0.0042051;161000];
sigmaX = [3.7024e-05;3220];
y_optim = lhsdesign(nS1,2);
x_optim = [norminv(y_optim(:,1),muX(1),sigmaX(1)) , norminv(y_optim(:,2),muX(2),sigmaX(2))];

% ��֪�Ż������ʵӦ�䣬��ѭ��Ӧ��-Ӧ�乫ʽ���Ż������ʵӦ��x_stress��
x_stress = ones(nS1,1);
for i = 1:nS1
    syms xstress
    xstress = abs(vpasolve(xstress/x_optim(i,2)+(xstress/K_)^(1/n_)-x_optim(i,1) == 0,xstress));
    x_stress(i,1) = xstress;
end

% �Ż����ƽ��Ӧ����
mean_stress = x_stress/2;

%histfit(x_strain);
%histfit(x(:,2));

%% �ĸ�ƣ�Ͳ�����������
% ���Ŷ�Ϊ0.5ʱ��
muX1 = [1.5943e3;-0.1058;250.6938;-1.4698];
sigmaX1 = [79.7150;5.5941e-8;20.055504;1.9847e-5];
y1 = lhsdesign(nS1,4);
x1 = [norminv(y1(:,1),muX1(1),sigmaX1(1)),norminv(y1(:,2),muX1(2),sigmaX1(2)),norminv(y1(:,3),muX1(3),sigmaX1(3)),norminv(y1(:,4),muX1(4),sigmaX1(4))];

% ���Ŷ�Ϊ0.9ʱ��
muX2 = [1.4776e3;-0.1058;76.6160;-1.4698];
sigmaX2 = [73.88;5.5941e-8;6.12928;1.9847e-5];
y2 = lhsdesign(nS1,4);
x2 = [norminv(y2(:,1),muX2(1),sigmaX2(1)),norminv(y2(:,2),muX2(2),sigmaX2(2)),norminv(y2(:,3),muX2(3),sigmaX2(3)),norminv(y2(:,4),muX2(4),sigmaX2(4))];

% ���Ŷ�Ϊ0.95ʱ��
muX3 = [1.4439e3;-0.1058;53.3976;-1.4698];
sigmaX3 = [72.195;5.5941e-8;4.271808;1.9847e-5];
y3 = lhsdesign(nS1,4);
x3 = [norminv(y3(:,1),muX3(1),sigmaX3(1)),norminv(y3(:,2),muX3(2),sigmaX3(2)),norminv(y3(:,3),muX3(3),sigmaX3(3)),norminv(y3(:,4),muX3(4),sigmaX3(4))];

%% ��������������ݺϲ�������һ��txt�ļ���
dataf= load('D:\file\Matlab_file\Paper_algorithmic_data\Aviation_paper_data\data.txt');
a = [dataf(:,11),x(:,1),mean_stress0(:,1),x_optim(:,1),mean_stress(:,1),x1(:,:),x2(:,:),x3(:,:)];
fid=fopen('new_data2.txt','wt');%д���ļ�·��
[m,n]=size(a);
 for i=1:1:m
    for j=1:1:n
       if j==n
         fprintf(fid,'%g\n',a(i,j));
      else
        fprintf(fid,'%g\t',a(i,j));
       end
    end
end
fclose(fid);


%% 
%syms E;
%mu = 4.4269444e-003; 
%E = vpasolve(mu - 7.1384833e2/E + (7.1384833e2/2516)^(1/0.133) == 0 , E)


%% ����ʼǣ�
% ѭ��ѹ��-Ӧ�����߳�����б��Ҫ�������ղ����ֲᡷ�ϵ�����б�ʶԱ�һ�£��������󲻴󣬲����ֱ���õ���ģ���Ĺ�ʽ��
% �Ż�������Ӧ���Ӧ��
