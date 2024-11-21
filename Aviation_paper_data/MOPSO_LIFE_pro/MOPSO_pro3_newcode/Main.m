clc;
clear;
close all;

%%
CostFunction = @(x) evaluate_objective(x);  %目标函数ZDT1
nVar = 2;                                     %变量个数
VarSize = [1 nVar];                            %变量矩阵大小
VarMin = 0;                                    %变量值定义域
VarMax = 360;                                  %注意: 该函数变量不能出现负值
MaxIt = 30;                                   %最大迭代次数
N = 40;                                        %种群规模
nRep = 50;                                     %档案库大小
w = 0.9;                                       %惯性权重系数
wdamp = 0.99;                                  %惯性权重衰减率
c1 = 1.7;                                      %个体学习因子
c2 = 1.8;                                      %全局学习因子
nGrid = 5;                                     %每一维的分格数
alpha = 0.1;                                   %膨胀率
beta = 2;                                      %最佳选择压
gamma = 2;                                     %删除选择压
mu = 0.1;                                      %变异概率
empty_particle.Position = [];                  %粒子位置向量
empty_particle.Velocity = [];                  %粒子速度向量
empty_particle.Cost = [];                      %粒子目标值向量
empty_particle.Best.Position = [];             %粒子最佳位置向量
empty_particle.Best.Cost = [];                 %粒子最佳目标值向量
empty_particle.IsDominated = [];               %粒子被支配个体向量
empty_particle.GridIndex = [];                 %粒子栅格索引向量
empty_particle.GridSubIndex = [];              %粒子栅格子索引向量
pop = repmat(empty_particle,N,1);              %repmat平铺矩阵%粒子初始空矩阵

for i = 1:N  %初始化N个个体
     % 产生服从均匀分布, VarSize大小的位置矩阵
     pop(i).Position = unifrnd(VarMin,VarMax,VarSize);
     pop(i).Velocity = zeros(VarSize);
     pop(i).Cost = CostFunction(pop(i).Position);
     pop(i).Best.Position = pop(i).Position;
     pop(i).Best.Cost = pop(i).Cost;
end

pop = DetermineDomination(pop);
rep = pop(~[pop.IsDominated]);
Grid = CreateGrid(rep,nGrid,alpha);
for i = 1:numel(rep)
 rep(i) = FindGridIndex(rep(i),Grid);
 % GridIndex = 绝对位置，.GridSubIndex = 坐标位置
end

%MOPSO主循环
 for it = 1:MaxIt
     for i = 1:N %逐一个体更新速度和位置，0.5的概率发生变异
         leader = SelectLeader(rep,beta);   %从支配个体轮盘赌选出全局最佳个体
         rep = [rep;pop(~[pop.IsDominated])];   %添加新的最佳栅格位置到库
         pop(i).Velocity = w*pop(i).Velocity + ...
             c1*rand(VarSize).*(pop(i).Best.Position-pop(i).Position)+ ...
             c2*rand(VarSize).*(leader.Position-pop(i).Position);    %速度更新
         pop(i).Position = pop(i).Position+pop(i).Velocity;   %位置更新
         pop(i).Position = limitToPosition(pop(i).Position,VarMin,VarMax);   %限制变量变化范围
         pop(i).Cost = CostFunction(pop(i).Position);   %计算目标函数值
         %应用变异策略
         pm = (1-(it-1)/(MaxIt-1)^(1/mu));  % 变异概率逐渐变小
         NewSol.Position = Mutate(pop(i).Position,pm,VarMin,VarMax);
         NewSol.Cost = CostFunction(NewSol.Position);   % 计算变异后的目标值
         if Dominates(NewSol,pop(i))
             pop(i).Position = NewSol.Position;
             pop(i).Cost  = NewSol.Cost;
         else %以0.5的概率决定是否接受变异
             if rand < 0.5
                 pop(i).Position = NewSol.Position;
                 pop(i).Cost = NewSol.Cost;
             end
         end
         if Dominates(pop(i),pop(i).Best)   % 如果当前个体优于先前最佳个体,则替换之
             pop(i).Best.Position = pop(i).Position;
             pop(i).Best.Cost = pop(i).Cost;
         else %以0.5的概率替换个体最佳
             if rand <0.5
                 pop(i).Best.Position = pop(i).Position;
                 pop(i).Best.Cost = pop(i).Cost;
             end
         end
     end   %每个个体
     
     rep =  DetermineDomination(rep);
     rep = rep(~[rep.IsDominated]);
     Grid = CreateGrid(rep,nGrid,alpha); 
    for i =1:numel(rep) 
        rep(i) = FindGridIndex(rep(i),Grid); 
    end 
    if numel(rep) > nRep 
        Extra = numel(rep)-nRep; 
        for e = 1:Extra 
            rep = DeleteOneRepMemebr(rep,gamma); 
        end 
    end 
     
    disp(['迭代次数 =',num2str(it)]); 
    w = w*wdamp; 
end 

figure(1); 
location = [rep.Cost];   %取最优结果
scatter3(location(1,:),location(2,:),location(3,:),'filled','b'); 
xlabel('f1');ylabel('f2'); zlabel('f3');
title('Pareto 最优边界图'); 
box on; 

