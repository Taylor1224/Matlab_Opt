%% ��ջ���
clc;clear all;close all;

%% Ŀ�꺯��
fun= @RMLS_Fun_Nf2;

%% ������Ⱥ����
% ����Ⱥ�Ż�������
sizepop = 50;                      % ��ʼ��Ⱥ����
dim = 4;                           % �ռ�ά��
ger = 100;                         % ����������   
x_ub = [1036; 675; 235; 510];      % ����λ�ò�������(�������ʽ���Զ�ά)
x_lb = [996; 646; 225; 490];
v_ub = 0.01*ones(dim,1);           % �����ٶ�����(ϵ�����ԸĶ���)
v_lb = -0.01*ones(dim,1);
c_1 = 0.5;                         % ����ѧϰ����
c_2 = 0.5;                         % Ⱥ��ѧϰ���� 
w1 = 0.9;
w2 = 0.4;
k = 0.6;

R0 = 0.985;                        % ��һ���Ż�����С�ɿ���Լ��
R = norminv(R0);
R0_ = 0.98;                         % �ڶ����Ż�����С�ɿ���Լ��

% �Ż�ǰϵͳ�����ľ�ֵ�ͷ��
Mu0 = 1.9180e3;      % 0.5���Ŷ���:1.9180e3; 0.9���Ŷ���:790.5919; 0.95���Ŷ���:604.2321;
Sigma0 = 277.1488;   % 0.5���Ŷ���:277.1488; 0.9���Ŷ���:104.1206; 0.95���Ŷ���:80.5513;

% �����Ż���ϵͳ�����ľ�ֵ�ͷ��
mu = 2.1892e3;       % 0.5���Ŷ���:2.1892e3; 0.9���Ŷ���:893.2872; 0.95���Ŷ���:680.8471;
sigma = 313.9817;    % 0.5���Ŷ���:313.9817; 0.9���Ŷ���:112.2232; 0.95���Ŷ���:88.9149;

% Latin����������������
nS = 1e5;nS1 = 1e4;K_ =1420 ;n_ = 0.05;Kt = 2.32;    % ҶƬ��n_=0.08;�̣�n_=0.05

% ���ˣ�
n = 500;   % ҶƬ��ʱ��500���̵�ʱ��280


%% ��ʼ����Ⱥ(λ�ú��ٶ�)(ʹ��MC����ԭʼ�����е�50����Ϊ��ʼ��Ⱥλ��)��
step = 1;
while true
    if step == 1
        data = readmatrix('BLISK4_D123.csv','Range','A2:D2001');
        data_y = readmatrix('BLISK4_D123.csv','Range','E2:E2001');
        % ���������ȡ 50 �����ݣ���֤�������̬�ֲ�
        num_samples = 50;
        pop_x = datasample(data, num_samples)';  
        y = datasample(data_y, num_samples)';

        pop_v = 0.5 * rands(4,50); 
        gbest = pop_x;

        % ����MLS��Ϻ������˺������� h ���Ը������ݵ��������
        for j=1:sizepop
            fitness_gbest(j) = fun(pop_x(:,j));                     
        end
        mu0 = mean(fitness_gbest);
        sigma0 = std(fitness_gbest);
        step = 2;
    elseif step == 2
        if (mu0-mu)/sqrt(sigma0^2-sigma^2)>R            
            zbest = pop_x(:,1);                        
            fitness_zbest = fitness_gbest(1);          
            for j=1:sizepop
                if -fitness_gbest(j) < -fitness_zbest    
                    zbest = pop_x(:,j);
                    fitness_zbest = fitness_gbest(j);
                end
            end
            break;
        else
            step = 1;
        end
    end
end

%% ����Ⱥ����
St = 1;
while true                   
    if St == 1
        iter = 1;                       
        record_x = zeros(dim, ger);
        record = zeros(ger, 1);         
        while iter <= ger
            step = 1;
            while true                   % �Կɿ���ΪԼ������ѭ��
                if step == 1
                    for j=1:sizepop
                        %    �����ٶȲ����ٶȽ��б߽紦��
                        w = (w1-w2)*tan(0.875*(1-(iter/ger)^k))+w2;
                        pop_v(:,j)= w*pop_v(:,j) + c_1*rand*(gbest(:,j)-pop_x(:,j))...
                        +c_2*rand*(zbest-pop_x(:,j));
                        for i=1:dim
                            if  pop_v(i,j) > v_ub(i)
                                pop_v(i,j) = v_ub(i);
                            end
                            if  pop_v(i,j) < v_lb(i)
                                pop_v(i,j) = v_lb(i);
                            end
                        end
                        
                        %    ����λ�ò���λ�ý��б߽紦��
                        pop_x(:,j) = pop_x(:,j) + pop_v(:,j);         
                        for i=1:dim
                            if  pop_x(i,j) > x_ub(i)
                                pop_x(i,j) = x_ub(i);
                            end
                            if  pop_x(i,j) < x_lb(i)
                                pop_x(i,j) = x_lb(i);
                            end
                        end
                        
                        %    ��������Ӧ����
                        if rand > 0.85
                            i=ceil(dim*rand);
                            pop_x(i,j)=x_lb(i) + (x_ub(i) - x_lb(i)) * rand;
                        end
                        fitness_pop(j) = fun(pop_x(:,j));             
                    end
                    mu0 = mean(fitness_pop);
                    sigma0 = std(fitness_pop);
                    step = 2;
                elseif step == 2
                    if (mu0-mu)/sqrt(sigma0^2-sigma^2)>R
                        for j=1:sizepop
                            %    ����Ӧ���������ʷ�����Ӧ�����Ƚ�
                            if -fitness_pop(j) < -fitness_gbest(j)        
                                gbest(:,j) = pop_x(:,j);                  
                                fitness_gbest(j) = fitness_pop(j);     
                            end   
                            
                            %    ������ʷ�����Ӧ������Ⱥ��ʷ�����Ӧ�����Ƚ�
                            if -fitness_gbest(j) < -fitness_zbest          
                                zbest = gbest(:,j);                     
                                fitness_zbest = fitness_gbest(j);        
                            end
                        end
                        break;
                    else
                        step = 1;
                    end
                end
            end
            record_x(:,iter) = zbest;
            record(iter) = fitness_zbest;                              
            iter = iter+1;
        end
        syms Sigma_
        cond = solve(Sigma_-sqrt(Sigma0^2-(Mu0-fitness_zbest)^2/R^2) == 0);
        Sigma = eval(cond);
        Sigma = imag(Sigma);

        z1 = 1-n./fitness_zbest;
        a1 = find(z1>0);
        R1 = length(a1)/length(z1);
        St = 2;
    elseif St == 2
        if R1 >= R0_            % ���˿ɿ���Լ����
            break;
        else
            St = 1;
        end
    end
end


%% ����������
figure
plot(record);title('�������������')
hold on;
disp('�Ż������');
disp(['�������ֵϵͳ������',num2str(fitness_zbest)]);
disp(['�������Ӧ���ı�׼�',num2str(Sigma)]);
disp('�����������ֵ��');
disp(zbest);


%% IPSO�Ż��㷨������άͼ��
X = record_x;
Y = meshgrid(record',record');
[x1,x2] = meshgrid(X(1,:),X(2,:));
figure;
surf(x1,x2,Y);


