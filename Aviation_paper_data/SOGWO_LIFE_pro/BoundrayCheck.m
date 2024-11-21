function pop_x=BoundrayCheck(pop_x,ub,lb,dim)
    for i=1:size(pop_x,1)
        for j=1:dim
            if pop_x(i,j)>ub(j)
                pop_x(i,j)=ub(j);
            end
            if pop_x(i,j)<lb(j)
                pop_x(i,j)=lb(j);
            end
        end
    end
end