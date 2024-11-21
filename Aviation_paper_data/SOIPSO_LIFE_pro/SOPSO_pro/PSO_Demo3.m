%% ��ջ���
clc;clear all;close all;
%% Ŀ�꺯��
fun= @B1_line;

%% ������Ⱥ����
%   ��Ҫ��������
sizepop = 50;                         % ��ʼ��Ⱥ����
dim = 4;                           % �ռ�ά��
ger = 50;                       % ����������   
xlimit_max = [1036; 675; 235; 510].*ones(dim,1);      % ����λ�ò�������(�������ʽ���Զ�ά)
xlimit_min = [996; 646; 225; 490].*ones(dim,1);
%xlimit_max = [1046; 685; 245; 520].*ones(dim,1);
%xlimit_min = [986; 636; 215; 480].*ones(dim,1);
vlimit_max = 1*ones(dim,1);      % �����ٶ�����
vlimit_min = -1*ones(dim,1);
c_1 = 0.8;                       % ����Ȩ��
c_2 = 0.5;                       % ����ѧϰ����
c_3 = 0.5;                       % Ⱥ��ѧϰ���� 

%% ���ɳ�ʼ��Ⱥ
%  ����������ɳ�ʼ��Ⱥλ��
%  Ȼ��������ɳ�ʼ��Ⱥ�ٶ�
%  Ȼ���ʼ��������ʷ���λ�ã��Լ�������ʷ�����Ӧ��
%  Ȼ���ʼ��Ⱥ����ʷ���λ�ã��Լ�Ⱥ����ʷ�����Ӧ��
for i=1:dim
    for j=1:sizepop
        pop_x(i,j) = xlimit_min(i)+(xlimit_max(i) - xlimit_min(i))*rand;  % ��ʼ��Ⱥ��λ��
        pop_v(i,j) = vlimit_min(i)+(vlimit_max(i) - vlimit_min(i))*rand;  % ��ʼ��Ⱥ���ٶ�
    end
end
gbest = pop_x;                                               % ÿ���������ʷ���λ��
for j=1:sizepop
    fitness_gbest(j) = fun(pop_x(:,j));                      % ÿ���������ʷ�����Ӧ��
end                  
zbest = pop_x(:,1);                          % ��Ⱥ����ʷ���λ��
fitness_zbest = fitness_gbest(1);            % ��Ⱥ����ʷ�����Ӧ��
for j=1:sizepop
    if fitness_gbest(j) < fitness_zbest       % �������Сֵ����Ϊ<; ��������ֵ����Ϊ>;  
        zbest = pop_x(:,j);
        fitness_zbest=fitness_gbest(j);
    end
end


%% ����Ⱥ����
%    �����ٶȲ����ٶȽ��б߽紦��    
%    ����λ�ò���λ�ý��б߽紦��
%    ��������Ӧ����
%    ����Լ�������жϲ���������Ⱥ��������λ�õ���Ӧ��
%    ����Ӧ���������ʷ�����Ӧ�����Ƚ�
%    ������ʷ�����Ӧ������Ⱥ��ʷ�����Ӧ�����Ƚ�
%    �ٴ�ѭ�������

iter = 1;                        %��������
record = zeros(ger, 1);          % ��¼��
while iter <= ger
    for j=1:sizepop
        %    �����ٶȲ����ٶȽ��б߽紦��
        pop_v(:,j)= c_1 * pop_v(:,j) + c_2*rand*(gbest(:,j)-pop_x(:,j))+c_3*rand*(zbest-pop_x(:,j));% �ٶȸ��¹�ʽ
        for i=1:dim
            if  pop_v(i,j) > vlimit_max(i)
                pop_v(i,j) = vlimit_max(i);
            end
            if  pop_v(i,j) < vlimit_min(i)
                pop_v(i,j) = vlimit_min(i);
            end
        end
        
        %    ����λ�ò���λ�ý��б߽紦��
        pop_x(:,j) = pop_x(:,j) + pop_v(:,j);% λ�ø���
        for i=1:dim
            if  pop_x(i,j) > xlimit_max(i)
                pop_x(i,j) = xlimit_max(i);
            end
            if  pop_x(i,j) < xlimit_min(i)
                pop_x(i,j) = xlimit_min(i);
            end
        end
        
        %    ��������Ӧ����
        if rand > 0.85
            i=ceil(dim*rand);
            pop_x(i,j)=xlimit_min(i) + (xlimit_max(i) - xlimit_min(i)) * rand;
        end
  
        fitness_pop(j) = fun(pop_x(:,j));                      % ��ǰ�������Ӧ��
        
        %    ����Ӧ���������ʷ�����Ӧ�����Ƚ�
        if fitness_pop(j) < fitness_gbest(j)       % �������Сֵ����Ϊ<; ��������ֵ����Ϊ>; 
            gbest(:,j) = pop_x(:,j);               % ���¸�����ʷ���λ��    
            if abs(fitness_gbest(j)-fitness_gbest(j))<=0.000001
                fitness_gbest(j) = fitness_pop(j); % ���¸�����ʷ�����Ӧ��
            else
                fitness_gbest(j) = fitness_gbest(j);
            end
        end   
        
        %    ������ʷ�����Ӧ������Ⱥ��ʷ�����Ӧ�����Ƚ�
        if fitness_gbest(j) < fitness_zbest        % �������Сֵ����Ϊ<; ��������ֵ����Ϊ>;  
            zbest = gbest(:,j);                    % ����Ⱥ����ʷ���λ��  
            if abs(fitness_gbest(j)-fitness_gbest(j))<=0.000001
                fitness_zbest=fitness_gbest(j);    % ����Ⱥ����ʷ�����Ӧ��   
            else
                fitness_gbest(j) = fitness_gbest(j);
            end
        end
    end
    record(iter) = fitness_zbest;%���ֵ��¼
    
    iter = iter+1;
end
%% ����������

plot(record);title('�������������')

disp('�Ż������');
% Ӧ������ֵ��
strain_permitted = 0.0043;

% Ӧ������ֵ��
stress_permitted = 690;

prob = sum(record <= stress_permitted)/ger;
disp(['����Ϊ��',num2str(prob)]);
disp(['�������ֵ��',num2str(fitness_zbest)]);
disp('���ű���ֵ��');
disp(zbest);

%% ���Ż���ı�׼�
R0 = 0.9954;
%R0 = 0.99;
R = norminv(R0);

disp(['��׼��̬�ֲ��ķ�������', num2str(R)])

syms sigma
mu = fitness_zbest;

% �Ż�ǰҶƬӦ��ľ�ֵ�ͷ��
mu0 = 4.4269444e-3;
sigma0 = 9.2887787e-5;

% �Ż�ǰҶƬӦ���ľ�ֵ�ͷ��
%mu0 = 0.71385e3;
%sigma0 = 0.13671e2;

% �Ż�ǰ��Ӧ��ľ�ֵ�ͷ��
%mu0 = 0.57601e-2;
%sigma0 = 0.13989e-3;

% �Ż�ǰ��Ӧ���ľ�ֵ�ͷ��
%mu0 = 0.94235e3;
%sigma0 = 0.21783e2;

cond = sigma-sqrt(sigma0^2-(mu0-mu)^2/R^2)>=0;
sol = solve(cond, sigma, 'ReturnConditions', true);
disp('Ӧ��ı�׼�');
disp(sol.conditions);

%% �󲻵�ʽ�����ӣ�
%syms x
%cond = x +6 < 10;
%sol = solve(cond, x ,"ReturnConditions",true);
%sol.conditions

%% ����������
   % ����ҶƬ��
B1_strain_mu0 = 4.4269444e-3;
B1_strain_sigma0 = 9.2887787e-5;
B1_stress_mu0 = 0.71385e3;
B1_stress_sigma0 = 0.13671e2;
B1_strain_mu = 0.0042051;
B1_strain_sigma = 3.7024e-05;
B1_stress_mu = 680.5919;
B1_stress_sigma = 4.8827;

   % �����̣�
D_strain0 = 4.4269444e-3;
D_stress0 = 0.71385e3;
D_strain_mu = 0.0042051;
D_strain_sigma = 3.7024e-05;
D_stress_mu = 680.5919;
D_stress_sigma = 4.8827;

% �Ż�ǰ�����������Latin����������������������
  % �Ż�ǰ��������


  % �Ż����������


% ���˺Ϳɿ��ȣ�

% �������ѧϰ�㷨/�����㷨��������Ⱥ�㷨���ϣ�
