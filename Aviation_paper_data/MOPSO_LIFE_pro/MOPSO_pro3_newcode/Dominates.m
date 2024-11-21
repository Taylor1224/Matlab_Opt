%============================= 
%判断两个目标值x,y的支配状态 
% x支配y,返回1;y支配x,返回0 
%============================= 
function b = Dominates(x,y) 
        if isstruct(x) 
            x=x.Cost; 
        end 
        if isstruct(y) 
            y=y.Cost; 
        end 
        b=all(x<=y) && any(x<y); 
end 

