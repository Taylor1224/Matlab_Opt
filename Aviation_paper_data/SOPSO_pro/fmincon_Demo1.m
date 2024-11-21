%% fmincon�Ż��㷨��Ѱ����Լ����Ԫ������ȫ�����Ž�

% fmincon �ṩ�˲�ͬ���Ż��㷨������ 'interior-point'��'sqp' �ȡ�
% ����ͨ������ options �����Ż����̣������������������Ż����̶ȵȲ������������������
clc
clear

%% ���뺯�������������

m = [5.17798e-2,0,0,0;
    0,6.89986e-2,0,0;
    0,0,2.10794e-1,0;
    0,0,0,1.09157e-1];
n = [5.26337e1,4.55018e1,4.84582e1,5.45754e1];% B1/D�޽���������н������һ����
o = [5.17798e-2,6.89986e-2,2.10794e-1,1.09157e-1];

% ����Ŀ�꺯����
fun = @B1_stress_have_Xterms;
% ��ʼ����Ʊ�����
x0 = [1010, 665, 230,495]*m-n;

%% ����Լ����
% �������Բ���ʽԼ��
A = [];
b = [];

nonlcon = [];

% �������Ե�ʽԼ��
Aeq = [];
beq = [];

% ������Ʊ��������½�
lb = [996, 646, 225, 490]*m-n;
ub = [1036, 675, 235, 510]*m-n;

%% ����ѡ����õ������
% ���� fmincon Ѱ��ȫ����Сֵ

options = optimoptions('fmincon', 'OutputFcn', @output,'Display','iter','Algorithm','sqp')
[x, fval] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options)

% [x, fval] = fmincon(fun1, [1010, 665, 230, 495]*m-n, Aeq, beq, lb, ub);

% ��ӡ�Ż����:
disp(['�����������ֵ��', num2str(x)]);
X = (x+n)./o;
disp(['������ʵ����ֵ��', num2str(X)]);
disp(['���ź���ֵ��', num2str(fval)]);




%  https://ww2.mathworks.cn/help/optim/ug/fmincon.html#busow0w-1��������������������ַ