%% 
clc
clear

%% ���Դ���particleswarm���������Ż���
% 1. ����Ŀ�꺯��
       % ��Function_1��
y = @B1_stress_have_Xterms;

m = [5.17798e-2,0,0,0;
    0,6.89986e-2,0,0;
    0,0,2.10794e-1,0;
    0,0,0,1.09157e-1];
n = [5.26337e1,4.55018e1,4.84582e1,5.45754e1];% B1/D�޽���������н������һ����
o = [5.17798e-2,6.89986e-2,2.10794e-1,1.09157e-1];

% 2. ����PSO�㷨����
options = optimoptions('particleswarm', 'OutputFcn', @(x, optimValues, state) output(x, optimValues, state),'SwarmSize', 100, 'MaxIterations', 100);

% 3. Լ������
  % �� ����ʽԼ����
% muY0 = [5.7600917e-3;4.4269444e-3];
% sigmaY0 = [1.3988648e-4;9.2887787e-5];
%mu=5.7600917e-3;     % ���þ�����
%sigma=1.3988648e-4;
%p=0.9;
%y=norminv(p,mu,sigma);

  % �� ���½�Լ����

lb = [996, 646, 225, 490]*m-n;% �����½�
ub = [1036, 675, 235, 510]*m-n;% �����Ͻ�
[x_opt, f_opt] = particleswarm(y, 4, lb, ub, options);

% 4. ����Ż����
disp('�Ż������');
disp(['�����������ֵ��', num2str(x_opt)]);
X = (x_opt+n)./o;
disp(['������ʵ����ֵ��', num2str(X)]);
disp(['���ź���ֵ��', num2str(f_opt)]);
%%
%��ͼ

