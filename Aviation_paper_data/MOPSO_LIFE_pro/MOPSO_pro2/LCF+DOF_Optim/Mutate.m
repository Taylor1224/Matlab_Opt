function xnew = Mutate(x, pm, x_lb, x_ub)

    nVar = numel(x(1,:));
    j = randi([1 nVar]);   % randi——是生成(0,N]间均匀分布的伪随机数，并且数都是整数

    dx = pm*(x_ub-x_lb);
    
    lb = x(:,j) - dx;
    if lb < x_lb
        lb = x_lb;
    end
    
    ub = x(:,j) + dx;
    if ub > x_ub
        ub = x_ub;
    end
    
    xnew = x;
    xnew(:,j) = unifrnd(lb, ub);

end