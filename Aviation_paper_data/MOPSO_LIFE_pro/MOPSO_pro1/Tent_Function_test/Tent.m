function [pop_x,pop_v] = Tent(a,dim,x_lb,x_ub,v_lb,v_ub)

    z(:,1) = rand(dim,1);% 随机序列
    for i=1:dim
        for j=2:sizepop
            if z(i,j-1)<a
                z(i,j) = z(i,j-1)/a;
            elseif z(i,j-1)>=a
                z(i,j) = (1-z(i,j-1))/(1-a);
            end
        end
    end

    pop_x = x_lb + z.*(x_ub - x_lb);  % 初始化种群位置
    pop_v = v_lb + z.*(v_ub - v_lb);  % 初始化种群位置

    gbest = pop_x;
    for j=1:sizepop
        fitness_gbest(j) = fun(pop_x(:,j));                      % 每个个体的历史最佳适应度
    end
    zbest = pop_x(:,1);                         % 种群的历史最佳位置
    fitness_zbest = fitness_gbest(1);           % 种群的历史最佳适应度
    for j=1:sizepop
        if fitness_gbest(j) < fitness_zbest     % 如果求最小值，则为<; 如果求最大值，则为>;
            zbest = pop_x(:,j);
            fitness_zbest=fitness_gbest(j);
        end
    end

end