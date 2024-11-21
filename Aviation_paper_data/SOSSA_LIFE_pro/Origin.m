clear all; 
close all;
clc;
 
%% 参数设置
N=50;                %麻雀个数
dim=2;               %评估函数维度
N_discoverer=0.7*N;  %发现者个数
N_Followers=0.1*N;   %追随者个数
N_Vigilant=0.2*N;    %警戒者个数
Max_iter=100;        %最大迭代次数
ST=0.6;              %安全阈值
 
%% 测试函数
f=@(x) sum(x.^2);
ub=10;%边界上限
lb=-10;%边界下限

%% 初始化
pop_x=lb+rand(N,dim).*(ub-lb);  %初始化麻雀种群

for i=1:N
    fitness(i)=f(pop_x(i,:));  %计算麻雀种群的适应度值
end
[A,index]=sort(fitness);
pop_x_best=pop_x(index(1),:);  %记录所有麻雀走过的位置的最优位置
pop_x_worst=pop_x(index(end),:);  %记录所有麻雀走过的位置的最差位置
best_fitness=A(1);  %记录所有麻雀走过的位置的最优值
worst_fitness=A(end);  %记录所有麻雀走过的位置的最差值
pop_x_best_currently=pop_x(index(1),:);  %记录当前麻雀种群最优位置
pop_x_worst_currently=pop_x(index(end),:);  %记录当前麻雀种群最差位置
best_fitness_currently=A(1);   %记录当前麻雀种群最优值
worst_fitness_currently=A(end);   %记录当前麻雀种群最差值
pop_x_discoverer=pop_x(index(1:N_discoverer),:);   %发现者位置
pop_x_Followers=pop_x(index(N_discoverer+1:N_discoverer+N_Followers),:);  %追随者位置
pop_x_Vigilant=pop_x(index(N_discoverer+N_Followers+1:N),:);   %警戒者位置
B=[-1,1];
F=best_fitness;  %记录每次迭代的麻雀走过的位置的最优值
iter=1;  %初始化迭代次数
 

 
%% 开始迭代更新
while iter<Max_iter
    for i=1:dim
        C(i)=B(round(rand)+1);
    end
    A=C'*inv((C*C'));
    R2=rand;
    %更新发现者位置
    for i=1:N_discoverer
        for j=1:dim
            if R2<ST
                pop_x_discoverer(i,j)=pop_x_discoverer(i,j)*exp(-i/rand*Max_iter);
            else
                 pop_x_discoverer(i,j)=pop_x_discoverer(i,j)+randn;
            end
        end
        %边界判断
        ub_flag=pop_x_discoverer(i,:)>ub;
        lb_flag=pop_x_discoverer(i,:)<lb;
        pop_x_discoverer(i,:)=(pop_x_discoverer(i,:).*(~(ub_flag+lb_flag)))+ub.*ub_flag+lb.*lb_flag;
    end
    %更新追随者位置
    for i=1:N_Followers
        for j=1:dim
            if i>N/2
                pop_x_Followers(i,j)=rand*exp((pop_x_worst_currently(j)-pop_x_Followers(i,j))/i^2);
            else
                pop_x_Followers(i,j)=pop_x_discoverer(1,j)+abs(pop_x_Followers(i,j)-pop_x_discoverer(1,j))*A(j);
            end
        end
        %边界判断
        ub_flag=pop_x_Followers(i,:)>ub;
        lb_flag=pop_x_Followers(i,:)<lb;
        pop_x_Followers(i,:)=(pop_x_Followers(i,:).*(~(ub_flag+lb_flag)))+ub.*ub_flag+lb.*lb_flag;
    end
    %更新警戒者位置
    for i=1:N_Vigilant
        for j=1:dim
            if f(pop_x_Vigilant(i,:))~=best_fitness_currently
                pop_x_Vigilant(i,j)=pop_x_best_currently(j)+randn*abs(pop_x_Vigilant(i,j)-pop_x_best_currently(j));
            else
                pop_x_Vigilant(i,j)=pop_x_Vigilant(i,j)+B(round(rand)+1)*(abs(pop_x_Vigilant(i,j)-pop_x_worst_currently(j)))/abs(f(pop_x_Vigilant(i,:))-worst_fitness_currently)+1;
            end
        end
        %边界判断
        ub_flag=pop_x_Vigilant(i,:)>ub;
        lb_flag=pop_x_Vigilant(i,:)<lb;
        pop_x_Vigilant(i,:)=(pop_x_Vigilant(i,:).*(~(ub_flag+lb_flag)))+ub.*ub_flag+lb.*lb_flag;
    end
    pop_x=[pop_x_discoverer;pop_x_Followers;pop_x_Vigilant];  %得到该次迭代下的所有麻雀的新位置
    for i=1:N
        fitness(i)=f(pop_x(i,:));   %计算适应度
    end
    [E,index]=sort(fitness);
    if f(pop_x(index(1),:))<best_fitness    %更新所有麻雀走过的位置的最优位置和最优值
        best_fitness=f(pop_x(index(1),:));
        pop_x_best=pop_x(index(1),:);
    end
    if f(pop_x(index(end),:))>worst_fitness   %更新所有麻雀走过的位置的最差位置和最差值
        worst_fitness= f(pop_x(index(end),:));
        pop_x_worst=pop_x(index(end),:);
    end
    pop_x_best_currently=pop_x(index(1),:);   %更新当前麻雀种群的最优位置
    pop_x_worst_currently=pop_x(index(end),:);  %更新当前麻雀种群的最差位置
    best_fitness_currently=E(1);   %更新当前麻雀的种群的最优值
    worst_fitness_currently=E(end);  %更新当前麻雀的种群的最差值
    pop_x_discoverer=pop_x(index(1:N_discoverer),:);  %重新选择种群中的发现者
    pop_x_Followers=pop_x(index(N_discoverer+1:N_discoverer+N_Followers),:);  %重新选择种群中的追随者
    pop_x_Vigilant=pop_x(index(N_discoverer+N_Followers+1:N),:);  %重新选择种群中的警戒者
    F=[F,best_fitness];   
    iter=iter+1;  %迭代次数加一
end
 
%% 结果 作图
display(['最优值是:',num2str(F(end)),'最优麻雀位置:',num2str(pop_x_best)]);
figure(1);
plot(F);
xlabel('迭代次数'),ylabel('适应度值');