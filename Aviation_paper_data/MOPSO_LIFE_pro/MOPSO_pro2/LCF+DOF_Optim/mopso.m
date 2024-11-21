clc;
clear;
close all;

%% Problem Definition
CostFunction = @BP_NN_Fun;      % Cost Function

nVar = 4;             % Number of Decision Variables
VarSize = [nVar 1];   % Size of Decision Variables Matrix
x_lb = [996; 646; 225; 490];          % 设置位置限制
x_ub = [1036; 675; 235; 510];  


%% MOPSO Parameters
MaxIt = 100;           % Maximum Number of Iterations   300
nPop = 200;            % Population Size   500
nRep = 100;            % Repository Size     最大输出存储尺寸
w = 0.5;              % Inertia Weight
wdamp = 0.99;         % Intertia Weight Damping Rate(阻尼率)
c1 = 1;               % Personal Learning Coefficient
c2 = 2;               % Global Learning Coefficient
nGrid = 7;            % Number of Grids(网格) per Dimension
alpha = 0.1;          % Inflation Rate
beta = 2;             % Leader Selection Pressure
gamma = 2;            % Deletion Selection Pressure
mu = 0.1;             % Mutation Rate（突变率）

%% Initialization
empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.Cost = [];

empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];      % Best也是一个结构体，它里面有Position和Cost

empty_particle.IsDominated = [];
empty_particle.GridIndex = [];
empty_particle.GridSubIndex = [];

pop = repmat(empty_particle, nPop, 1);    % empty_particle是一个结构体，pop是把这个结构体复制200份后得到的另一个结构体。

for i=1:nPop
    for j=1:nVar
        pop_x = x_lb(j)+(x_ub(j) - x_lb(j))*rand;
        pop(i).Position = [pop(i).Position;pop_x];  % 初始化种群的位置(4行1列) 
    end
    pop(i).Velocity = zeros(VarSize);               % 初始种群的速度(4行1列)
    pop(i).Cost = CostFunction(pop(i).Position);    % 得到对应目标函数值

    % Update Personal Best
    pop(i).Best.Position = pop(i).Position;
    pop(i).Best.Cost = pop(i).Cost;
end

% Determine Domination
pop = DetermineDomination(pop);

rep = pop(~[pop.IsDominated]);   % 为了找到帕累托前沿

Grid = CreateGrid(rep, nGrid, alpha);

for i = 1:numel(rep)
    rep(i) = FindGridIndex(rep(i), Grid);
end


%% MOPSO Main Loop
for it = 1:MaxIt
    
    for i = 1:nPop
        
        leader = SelectLeader(rep, beta);
        
        pop(i).Velocity = w*pop(i).Velocity ...
            +c1*rand(VarSize).*(pop(i).Best.Position-pop(i).Position) ...
            +c2*rand(VarSize).*(leader.Position-pop(i).Position);
        
        pop(i).Position = pop(i).Position + pop(i).Velocity;
        
        pop(i).Position = max(pop(i).Position, x_lb);
        pop(i).Position = min(pop(i).Position, x_ub);
        
        pop(i).Cost = CostFunction(pop(i).Position);
        
        % Apply Mutation
        pm = (1-(it-1)/(MaxIt-1))^(1/mu);
        if rand<pm
            NewSol.Position = Mutate(pop(i).Position, pm, x_lb, x_ub);
            NewSol.Cost = CostFunction(NewSol.Position);
            if Dominates(NewSol, pop(i))
                pop(i).Position = NewSol.Position;
                pop(i).Cost = NewSol.Cost;

            elseif Dominates(pop(i), NewSol)
                % Do Nothing

            else
                if rand<0.5
                    pop(i).Position = NewSol.Position;
                    pop(i).Cost = NewSol.Cost;
                end
            end
        end
        
        if Dominates(pop(i), pop(i).Best)
            pop(i).Best.Position = pop(i).Position;
            pop(i).Best.Cost = pop(i).Cost;
            
        elseif Dominates(pop(i).Best, pop(i))
            % Do Nothing
            
        else
            if rand<0.5
                pop(i).Best.Position = pop(i).Position;
                pop(i).Best.Cost = pop(i).Cost;
            end
        end
        
    end
    
    % Add Non-Dominated Particles to REPOSITORY
    rep = [rep;pop(~[pop.IsDominated])]; %#ok
    
    % Determine Domination of New Resository Members
    rep = DetermineDomination(rep);
    
    % Keep only Non-Dminated Memebrs in the Repository
    rep = rep(~[rep.IsDominated]);% ~为取反
    
    % Update Grid
    Grid = CreateGrid(rep, nGrid, alpha);

    % Update Grid Indices
    for i = 1:numel(rep)
        rep(i) = FindGridIndex(rep(i), Grid);
    end
    
    % Check if Repository is Full
    if numel(rep)>nRep
        
        Extra = numel(rep)-nRep;
        for e = 1:Extra
            rep = DeleteOneRepMemebr(rep, gamma);
        end
        
    end
    
    % Plot Costs
    figure(1);
    PlotCosts(pop, rep);
    pause(0.01);
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Number of Rep Members = ' num2str(numel(rep))]);
    
    % Damping Inertia Weight
    w = w*wdamp;
    
end

%% Resluts
% 注意：用来表示结构体中的第几个结构体用：rep(3)
% 优化前——叶盘应力和径向变形均值：
str_mu_1 = 7.1420128e+002;
str_mu_2 = 9.4279947e+002;
DOF_mu_1 = 9.9077121e-001;
DOF_mu_2 = 3.0531561e-001;
W1 = 0.2;
W2 = 0.2;
W3 = 0.2;
W4 = 0.4;
% 把这四个权重值变成动态的，定义值作为它们的均值，四个权值的和为1.
% 考虑用Copula函数得到相关性。


F1 = W1*((str_mu_1-rep(1).Cost(1,1))/str_mu_1) + W2*((str_mu_2-rep(1).Cost(2,1))/str_mu_2) + ...
    W3*((DOF_mu_1-rep(1).Cost(3,1))/DOF_mu_1) + W4*((DOF_mu_2-rep(1).Cost(4,1))/DOF_mu_2);

F2 = W1*((str_mu_1-rep(2).Cost(1,1))/str_mu_1) + W2*((str_mu_2-rep(2).Cost(2,1))/str_mu_2) + ...
    W3*((DOF_mu_1-rep(2).Cost(3,1))/DOF_mu_1) + W4*((DOF_mu_2-rep(2).Cost(4,1))/DOF_mu_2);

F3 = W1*((str_mu_1-rep(3).Cost(1,1))/str_mu_1) + W2*((str_mu_2-rep(3).Cost(2,1))/str_mu_2) + ...
    W3*((DOF_mu_1-rep(3).Cost(3,1))/DOF_mu_1) + W4*((DOF_mu_2-rep(3).Cost(4,1))/DOF_mu_2);

F4 = W1*((str_mu_1-rep(4).Cost(1,1))/str_mu_1) + W2*((str_mu_2-rep(4).Cost(2,1))/str_mu_2) + ...
    W3*((DOF_mu_1-rep(4).Cost(3,1))/DOF_mu_1) + W4*((DOF_mu_2-rep(4).Cost(4,1))/DOF_mu_2);
