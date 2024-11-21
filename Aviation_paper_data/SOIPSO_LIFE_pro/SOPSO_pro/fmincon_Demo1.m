%% fmincon�Ż��㷨��Ѱ����Լ����Ԫ������ȫ�����Ž�

% fmincon �ṩ�˲�ͬ���Ż��㷨������ 'interior-point'��'sqp' �ȡ�
% ����ͨ������ options �����Ż����̣������������������Ż����̶ȵȲ������������������
clc
clear

%% ���뺯�������������

% ����Ŀ�꺯����
fun = @BP_NN_function;
% ��ʼ����Ʊ�����
x0 = [1010, 665, 230,495];

%% ����Լ����
% �������Բ���ʽԼ��
A = [];
b = [];

nonlcon = [];

% �������Ե�ʽԼ��
Aeq = [];
beq = [];

% ������Ʊ��������½�
ub = [1036; 675; 235; 510];
lb = [996; 646; 225; 490];

%% ����ѡ����õ������
% ���� fmincon Ѱ��ȫ����Сֵ

options = optimoptions('fmincon', 'OutputFcn', @output,'Display','iter','Algorithm','sqp');
[x, fval] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);

% [x, fval] = fmincon(fun1, [1010, 665, 230, 495]*m-n, Aeq, beq, lb, ub);

% ��ӡ�Ż����:
disp(['�����������ֵ��', num2str(x)]);
X = (x+n)./o;
disp(['������ʵ����ֵ��', num2str(X)]);
disp(['���ź���ֵ��', num2str(fval)]);




%  https://ww2.mathworks.cn/help/optim/ug/fmincon.html#busow0w-1��������������������ַ