%============================= 
%轮盘赌选择一个较好的支配个体 
%============================= 
function i = RouletteWheelSelection(P) 
        r = rand; 
        C = cumsum(P); 
        i = find(r<=C,1,'first'); 
end

