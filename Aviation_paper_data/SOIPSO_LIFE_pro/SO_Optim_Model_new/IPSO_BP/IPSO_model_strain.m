%% ��ջ���
clc;clear all; close all;

%% Ŀ�꺯��
fun= @D_have_Xterms;

%% ������Ⱥ����
% ����Ⱥ�Ż�������
sizepop = 50;                      % ��ʼ��Ⱥ����
dim = 4;                           % �ռ�ά��
ger = 100;                         % ����������   
x_ub = [1036; 675; 235; 510];      % ����λ�ò�������(�������ʽ���Զ�ά)
x_lb = [996; 646; 225; 490];
v_ub = 0.01*ones(dim,1);           % �����ٶ�����(ϵ�����ԸĶ���)
v_lb = -0.01*ones(dim,1);
c_1 = 0.8;                         % ����Ȩ��
c_2 = 0.5;                         % ����ѧϰ����
c_3 = 0.5;                         % Ⱥ��ѧϰ���� 

R0 = 0.985;                        % ��һ���Ż�����С�ɿ���Լ��
R = norminv(R0);
R0_ = 0.99;                         % �ڶ����Ż�����С�ɿ���Լ��

% �Ż�ǰ��ҶƬӦ��ľ�ֵ�ͷ��
% Mu0 = 4.4295e-3;
% Sigma0 = 1.228e-4;
% �Ż�ǰ����Ӧ��ľ�ֵ�ͷ��
Mu0 = 5.7633e-3;
Sigma0 = 1.484e-4;

% �����Ż����ҶƬӦ��ľ�ֵ�ͷ��
% mu = 4.3e-3;
% sigma = 6e-5;
% �����Ż������Ӧ��ľ�ֵ�ͷ��
mu = 5.5e-3;
sigma = 5e-5;

% Latin����������������
nS = 1e5;nS1 = 1e4;K_ =1420 ;n_ = 0.08;Kt = 2.32;

% ���ˣ�
n = 280;

%% �Ż�ǰ�ĵ�ЧӦ��͵�Чƽ��Ӧ��������
% ����ҶƬ��Ӧ��muX0(1)��ƽ��Ӧ��muX0(2)������ģ��E muX0(3)��
% muX0 = [4.4295e-3;714.201;161000];
% sigmaX0 = [1.228e-4;17.809;3220];
% �����̵�Ӧ��muX0(1)��ƽ��Ӧ��muX0(2)������ģ��E muX0(3)��
muX0 = [5.7633e-3;942.799;161000];
sigmaX0 = [1.484e-4;22.994;3220];

y = lhsdesign(nS1,3);
x = [norminv(y(:,1),muX0(1),sigmaX0(1)) , norminv(y(:,2),muX0(2),sigmaX0(2)), norminv(y(:,3),muX0(3),sigmaX0(3))];


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


%% Tent����ӳ���ʼ����Ⱥ(λ�ú��ٶ�)��
%  ����������ɳ�ʼ��Ⱥλ��
%  Ȼ��������ɳ�ʼ��Ⱥ�ٶ�
%  Ȼ���ʼ��������ʷ���λ�ã��Լ�������ʷ�����Ӧ��
%  Ȼ���ʼ��Ⱥ����ʷ���λ�ã��Լ�Ⱥ����ʷ�����Ӧ��
step = 1;
while true
    if step == 1
        a = 0.499;
        z(:,1) = rand(dim,1); % �������
        for i=1:dim
            for j=2:sizepop
                if z(i,j-1)<a
                    z(i,j) = z(i,j-1)/a;
                elseif z(i,j-1)>=a
                    z(i,j) = (1-z(i,j-1))/(1-a);
                end     
            end
        end
        pop_x = x_lb + z.*(x_ub - x_lb);  % ��ʼ����Ⱥλ��
        pop_v = v_lb + z.*(v_ub - v_lb);  % ��ʼ����Ⱥ�ٶ�
        gbest = pop_x; 
        for j=1:sizepop
            fitness_gbest(j) = fun(pop_x(:,j));                      % ÿ���������ʷ�����Ӧ��
        end
        mu0 = mean(fitness_gbest);
        sigma0 = std(fitness_gbest);
        step = 2;
    elseif step == 2
        if (mu0-mu)/sqrt(sigma0^2-sigma^2)>R
            zbest = pop_x(:,1);                         % ��Ⱥ����ʷ���λ��
            fitness_zbest = fitness_gbest(1);           % ��Ⱥ����ʷ�����Ӧ��
            for j=1:sizepop
                if fitness_gbest(j) < fitness_zbest     % �������Сֵ����Ϊ<; ��������ֵ����Ϊ>;  
                    zbest = pop_x(:,j);
                    fitness_zbest=fitness_gbest(j);
                end
            end
            break;
        else
            step = 1;
        end
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

St = 1;
while true
    if St == 1
        iter = 1;                        %��������
        record = zeros(ger, 1);          % ��¼��
        while iter <= ger
            step = 1;
            while true                   % �Կɿ���ΪԼ������ѭ��
                if step == 1
                    for j=1:sizepop
                        %    �����ٶȲ����ٶȽ��б߽紦��
                        pop_v(:,j)= c_1 * pop_v(:,j) + c_2*rand*(gbest(:,j)-pop_x(:,j))+c_3*rand*(zbest-pop_x(:,j));% �ٶȸ��¹�ʽ
                        for i=1:dim
                            if  pop_v(i,j) > v_ub(i)
                                pop_v(i,j) = v_ub(i);
                            end
                            if  pop_v(i,j) < v_lb(i)
                                pop_v(i,j) = v_lb(i);
                            end
                        end
                        
                        %    ����λ�ò���λ�ý��б߽紦��
                        pop_x(:,j) = pop_x(:,j) + pop_v(:,j);          % λ�ø���
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
                        fitness_pop(j) = fun(pop_x(:,j));              % ��ǰ�������Ӧ��
                    end
                    mu0 = mean(fitness_pop);
                    sigma0 = std(fitness_pop);
                    step = 2;
                elseif step == 2
                    if (mu0-mu)/sqrt(sigma0^2-sigma^2)>R
                        for j=1:sizepop
                            %    ����Ӧ���������ʷ�����Ӧ�����Ƚ�
                            if fitness_pop(j) < fitness_gbest(j)       % �������Сֵ����Ϊ<; ��������ֵ����Ϊ>; 
                                gbest(:,j) = pop_x(:,j);               % ���¸�����ʷ���λ��    
                                fitness_gbest(j) = fitness_pop(j);     % ���¸�����ʷ�����Ӧ��
                            end   
                            
                            %    ������ʷ�����Ӧ������Ⱥ��ʷ�����Ӧ�����Ƚ�
                            if fitness_gbest(j) < fitness_zbest        % �������Сֵ����Ϊ<; ��������ֵ����Ϊ>;  
                                zbest = gbest(:,j);                    % ����Ⱥ����ʷ���λ��  
                                fitness_zbest = fitness_gbest(j);      % ����Ⱥ����ʷ�����Ӧ��   
                            end
                        end
                        break;
                    else
                        step = 1;
                    end
                end
            end
            record(iter) = fitness_zbest;                              %����ֵ��¼,�õ��Ż����Ӧ����Ӧ�䡣
            iter = iter+1;
        end
        
        syms Sigma_
        cond = solve(Sigma_-sqrt(Sigma0^2-(Mu0-fitness_zbest)^2/R^2) == 0);
        Sigma = eval(cond);

        % �Ż����Ӧ���Ӧ��������
        % BP������ģ���Ż���õ�ҶƬ���̵�Ӧ���ֵ�ͷ��
        BP_mu = fitness_zbest;
        BP_sigma = Sigma;
        
        % �Ż���Ӧ��muX(1)��ƽ��Ӧ��x_stress��
        muX = [BP_mu;161000]; sigmaX = [BP_sigma;3220];
        y_optim = lhsdesign(nS1,2);
        x_optim = [norminv(y_optim(:,1),muX(1),sigmaX(1)) , norminv(y_optim(:,2),muX(2),sigmaX(2))];
        
        % ��֪�Ż����Ӧ�䣬��ѭ��Ӧ��-Ӧ�乫ʽ���Ż����Ӧ��x_stress��
        x_strain = x_optim(:,1);
        x_stress = ones(nS1,1);
        for i = 1:nS1
            syms xstress
            xstress = abs(vpasolve(xstress/x_optim(i,2)+(xstress/K_)^(1/n_)-x_strain(i,1) == 0,xstress));
            x_stress(i,1) = xstress;
        end
        
        % �Ż�ǰ���Ż����������
        e = x(:,3); m0 = x(:,1); n0 = x(:,2);
        m = x_strain; n = x_stress;
        k1 = x1(:,1); b1 = x1(:,2); s1 = x1(:,3); c1 = x1(:,4); % 0.5���Ŷ��µ�ƣ�Ͳ���
        k2 = x2(:,1); b2 = x2(:,2); s2 = x2(:,3); c2 = x2(:,4); % 0.9���Ŷ��µ�ƣ�Ͳ���
        k3 = x3(:,1); b3 = x3(:,2); s3 = x3(:,3); c3 = x3(:,4); % 0.95���Ŷ��µ�ƣ�Ͳ���
        
        for i=1:length(m0)
            % �Ż�ǰ��
            Life_x(i)=fminbnd(@(x)abs(m0(i)/2-(k1(i)-n0(i))/e(i)*(2*x)^b1(i)-s1(i)*(2*x)^c1(i)),100,100000000000000000000);% xΪ���Ŷ�0.5ʱ������
            Life_y(i)=fminbnd(@(y)abs(m0(i)/2-(k2(i)-n0(i))/e(i)*(2*y)^b2(i)-s2(i)*(2*y)^c2(i)),100,100000000000000000000);% yΪ���Ŷ�0.9ʱ������
            Life_z(i)=fminbnd(@(z)abs(m0(i)/2-(k3(i)-n0(i))/e(i)*(2*z)^b3(i)-s3(i)*(2*z)^c3(i)),100,100000000000000000000);% zΪ���Ŷ�0.95ʱ������
            % �Ż���
            Life_Optim_x(i)=fminbnd(@(x)abs(m(i)/2-(k1(i)-n(i))/e(i)*(2*x)^b1(i)-s1(i)*(2*x)^c1(i)),100,100000000000000000000);% xΪ���Ŷ�0.5ʱ������
            Life_Optim_y(i)=fminbnd(@(y)abs(m(i)/2-(k2(i)-n(i))/e(i)*(2*y)^b2(i)-s2(i)*(2*y)^c2(i)),100,100000000000000000000);% yΪ���Ŷ�0.9ʱ������
            Life_Optim_z(i)=fminbnd(@(z)abs(m(i)/2-(k3(i)-n(i))/e(i)*(2*z)^b3(i)-s3(i)*(2*z)^c3(i)),100,100000000000000000000);% zΪ���Ŷ�0.95ʱ������
        end

        z1 = 1-n./Life_Optim_z;
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
disp(['�������ֵ(Ӧ�䣬��Ϊ��ֵ)��',num2str(fitness_zbest)]);
disp(['�������Ӧ���ı�׼�',num2str(Sigma)]);
disp('�����������ֵ��');
disp(zbest);

disp('�Ż�ǰ��ͬ���Ŷ��µ�������ֵ����0.5��0.9��0.95��');
disp(mean(Life_x));
disp(mean(Life_y));
disp(mean(Life_z));
disp('�Ż���ͬ���Ŷ��µ�������ֵ����0.5��0.9��0.95��');
disp(mean(Life_Optim_x));
disp(mean(Life_Optim_y));
disp(mean(Life_Optim_z));

%% �����Ż�ǰ��������������
fid = fopen('mu_B1_LIFE0.5.txt','wt');      % wt������д�������
fprintf(fid,'%g\n',Life_x); 
fid = fopen('mu_B1_LIFE0.9.txt','wt');
fprintf(fid,'%g\n',Life_y); 
fid = fopen('mu_B1_LIFE0.95.txt','wt');
fprintf(fid,'%g\n',Life_z); 
fid = fopen('mu_B1_LIFE0.5_optim.txt','wt');
fprintf(fid,'%g\n',Life_Optim_x); 
fid = fopen('mu_B1_LIFE0.9_optim.txt','wt');
fprintf(fid,'%g\n',Life_Optim_y); 
fid = fopen('mu_B1_LIFE0.95_optim.txt','wt');
fprintf(fid,'%g\n',Life_Optim_z); 


%% �Ż�ǰ���������ֲ�ͼ��

  % �������Ŷȵ�ͼ��һ��ҶƬB1���������ֲ�ͼ��
B1 = load('mu_B1_LIFE0.5.txt');B2 = load('mu_B1_LIFE0.9.txt');B3 = load('mu_B1_LIFE0.95.txt');
b1 = load('mu_B1_LIFE0.5_optim.txt');b2 = load('mu_B1_LIFE0.9_optim.txt');b3 = load('mu_B1_LIFE0.95_optim.txt');

subplot(1,2,1);
histfit(B1(:,1),200,'lognorm');
hold on;
histfit(B2(:,1),200,'lognorm');
hold on;
histfit(B3(:,1),200,'lognorm');
xlabel('Blade LCF life before optimization');
ylabel('Relative frequency');
legend('Blade LCF life sample with P=0.5','Fitted lognormal curve with P=0.5','Blade LCF life sample with P=0.9','Fitted lognormal curve with P=0.9','Blade LCF life sample with P=0.95','Fitted lognormal curve with P=0.95');

subplot(1,2,2);
histfit(b1(:,1),200,'lognorm');
hold on;
histfit(b2(:,1),200,'lognorm');
hold on;
histfit(b3(:,1),200,'lognorm');
xlabel('Blade LCF life after optimization');
ylabel('Relative frequency');
legend('Blade LCF life sample with P=0.5','Fitted lognormal curve with P=0.5','Blade LCF life sample with P=0.9','Fitted lognormal curve with P=0.9','Blade LCF life sample with P=0.95','Fitted lognormal curve with P=0.95');


  % �������Ŷȵ�ͼ��һ����D���������ֲ�ͼ��
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


%% �Ż�ǰ������������ʷͼ�ͷֲ�ֱ��ͼ����һ�𣩣�
  % ���Ŷ�0.5�£�ҶƬB1����D���������ֲ�ͼ��
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

  % ���Ŷ�0.9�£�ҶƬB1����D���������ֲ�ͼ��
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

  % ���Ŷ�0.95�£�ҶƬB1����D���������ֲ�ͼ��
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

%% ��txt�ļ�תΪexcel�ļ���һ�У���
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


%% ��txt�ļ�תΪexcel�ļ���һ�У���
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