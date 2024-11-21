%% Start
clc;
clear;
close all;

%% Import Function and Parameters:

CostFunction = @B_stress_Xterms;      % Cost Function

MaxIt = 100;      % Maximum Number of Iterations

nPop = 50;        % Population Size

nVar = 4;             % Number of Decision Variables

VarSize = [1 nVar];   % Size of Decision Variables Matrix

x_lb = [996; 646; 225; 490];          % 设置位置限制
x_ub = [1036; 675; 235; 510];        
v_ub = 0.01*ones(nVar,1);               % 设置速度限制(系数可以改动！)
v_lb = -0.01*ones(nVar,1);
c_1 = 0.8;                              % 惯性权重
c_2 = 0.5;                              % 自我学习因子
c_3 = 0.5;                              % 群体学习因子 

% Number of Objective Functions（只是为了得到目标函数的个数）
nObj = numel(CostFunction(unifrnd(x_lb(1), x_ub(1), VarSize)));   % unifrnd:生成连续均匀的随机数,x_lb(1)为最低值，x_ub(1)为最高值，VarSize指定几行几列，这里是1行4列


%% NSGA-II 

pCrossover = 0.7;                         % Crossover Percentage
nCrossover = 2*round(pCrossover*nPop/2);  % Number of Parnets (Offsprings)

pMutation = 0.4;                          % Mutation Percentage
nMutation = round(pMutation*nPop);        % Number of Mutants

mu = 0.02;                    % Mutation Rate

sigma = 0.1*(x_ub-x_lb);  % Mutation Step Size


%% Initialization

% 创建工作区表格标题栏：
empty_individual.Position = [];
empty_individual.Cost = [];
empty_individual.velocity = [];
empty_individual.Rank = [];
empty_individual.DominationSet = [];
empty_individual.DominatedCount = [];
empty_individual.CrowdingDistance = [];

pop = repmat(empty_individual, nPop, 1); % repmat(A,2,2):表示将A重复到2x2的矩阵中。

% Initialize Population:

  % 普通初始化种群(位置和速度)：

for i=1:nPop
    for j=1:nVar
        pop_x = x_lb(j)+(x_ub(j) - x_lb(j))*rand;
        pop(i).Position = [pop_x,pop(i).Position];  % 初始种群的位置(1行4列)
        pop_v = v_lb(i)+(v_ub(i) - v_lb(i))*rand; 
        pop(i).velocity = [pop_v,pop(i).velocity];  % 初始种群的速度(1行4列)
    end
    pop(i).Cost = CostFunction(pop(i).Position); % 得到对应目标函数值
end


% Non-Dominated Sorting
[pop, F] = NonDominatedSorting(pop);

% Calculate Crowding Distance
pop = CalcCrowdingDistance(pop, F);

% Sort Population
[pop, F] = SortPopulation(pop);


%% MOPSO Main Loop

for it = 1:MaxIt
    
    % Crossover
    popc = repmat(empty_individual, nCrossover/2, 2);
    for k = 1:nCrossover/2
        
        i1 = randi([1 nPop]);
        p1 = pop(i1);
        
        i2 = randi([1 nPop]);
        p2 = pop(i2);
        
        [popc(k, 1).Position, popc(k, 2).Position] = Crossover(p1.Position, p2.Position);
        
        popc(k, 1).Cost = CostFunction(popc(k, 1).Position);
        popc(k, 2).Cost = CostFunction(popc(k, 2).Position);
        
    end
    popc = popc(:);
    
    % Mutation
    popm = repmat(empty_individual, nMutation, 1);
    for k = 1:nMutation
        
        i = randi([1 nPop]);
        p = pop(i);
        
        popm(k).Position = Mutate(p.Position, mu, sigma);
        
        popm(k).Cost = CostFunction(popm(k).Position);
        
    end
    
    % Merge
    pop = [pop
         popc
         popm]; %#ok



%% 粒子群迭代
iter = 1;                        %迭代次数
record = zeros(MaxIt, 1);          % 记录器 
popt = repmat(empty_individual, MaxIt, 1);
while iter <= MaxIt
    for j=1:nPop
        %    更新速度并对速度进行边界处理
        pop_v(:,j)= c_1 * pop_v(:,j) + c_2*rand*(gbest(:,j)-pop_x(:,j))+c_3*rand*(zbest-pop_x(:,j));% 速度更新公式
        for i=1:dim
            if  pop_v(i,j) > v_ub(i)
                pop_v(i,j) = v_ub(i);
            end
            if  pop_v(i,j) < v_lb(i)
                pop_v(i,j) = v_lb(i);
            end
        end
        
        %    更新位置并对位置进行边界处理
        pop_x(:,j) = pop_x(:,j) + pop_v(:,j);% 位置更新
        for i=1:dim
            if  pop_x(i,j) > x_ub(i)
                pop_x(i,j) = x_ub(i);
            end
            if  pop_x(i,j) < x_lb(i)
                pop_x(i,j) = x_lb(i);
            end
        end
        
        %    进行自适应变异
        if rand > 0.85
            i=ceil(dim*rand);
            pop_x(i,j)=x_lb(i) + (x_ub(i) - x_lb(i)) * rand;
        end
  
        fitness_pop(j) = fun(pop_x(:,j));                      % 当前个体的适应度
        
        %    新适应度与个体历史最佳适应度做比较
        if fitness_pop(j) < fitness_gbest(j)       % 如果求最小值，则为<; 如果求最大值，则为>; 
            gbest(:,j) = pop_x(:,j);               % 更新个体历史最佳位置    
            if abs(fitness_gbest(j)-fitness_gbest(j))<=0.000001
                fitness_gbest(j) = fitness_pop(j); % 更新个体历史最佳适应度
            else
                fitness_gbest(j) = fitness_gbest(j);
            end
        end   
        
        %    个体历史最佳适应度与种群历史最佳适应度做比较
        if fitness_gbest(j) < fitness_zbest        % 如果求最小值，则为<; 如果求最大值，则为>;  
            zbest = gbest(:,j);                    % 更新群体历史最佳位置  
            if abs(fitness_gbest(j)-fitness_gbest(j))<=0.000001
                fitness_zbest=fitness_gbest(j);    % 更新群体历史最佳适应度   
            else
                fitness_gbest(j) = fitness_gbest(j);
            end
        end
    end
    record(iter) = fitness_zbest;%最大值记录
    
    iter = iter+1;
end





%%

     
    % Non-Dominated Sorting
    [pop, F] = NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop = CalcCrowdingDistance(pop, F);

    % Sort Population
    pop = SortPopulation(pop);
    
    % Truncate
    pop = pop(1:nPop);
    
    % Non-Dominated Sorting
    [pop, F] = NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop = CalcCrowdingDistance(pop, F);

    % Sort Population
    [pop, F] = SortPopulation(pop);
    
    % Store F1
    F1 = pop(F{1});
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(numel(F1))]);
    
    % Plot F1 Costs
    figure(1);
    PlotCosts(F1);
    pause(0.01);
    
end

%% Results

