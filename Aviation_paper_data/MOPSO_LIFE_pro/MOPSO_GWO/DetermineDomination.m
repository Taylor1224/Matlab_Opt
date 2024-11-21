%============================= 
%判断全局支配状况,返回0 = 非支配解 
%============================= 
function pop =DetermineDomination(pop) 
        nPop = numel(pop); 
        for i =1:nPop 
            pop(i).IsDominated = false;   %初始化为互不支配 
        end 
        for i = 1:nPop-1 
            for j = i+1:nPop 
                if Dominates(pop(i),pop(j)) 
                    pop(j).IsDominated = true; 
                end 
                    if Dominates(pop(j),pop(i)) 
                        pop(i).IsDominated = true; 
                    end 
            end 
        end 
end 

