%%
clc;

%% ������Լ���ֲ�����
% ����Ŀ�꺯��
fitnessfcn = @fun;

% ��ʼ�²�
x0 = [996, 646, 225, 490];

% �����Ż�ѡ��
options = optimset('Display', 'iter'); % ����ѡ����ʾ������Ϣ

% ʹ��fminunc����������Լ���Ż�
[x, fval] = fminunc(fitnessfcn, x0, options);

% ��ӡ�Ż����
fprintf('��Сֵ��: x = [%f, %f, %f, %f]\n', x(1), x(2), x(3), x(4));
fprintf('��Сֵ: f(x) = %f\n', fval);
