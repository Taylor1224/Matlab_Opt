%============================= 
%限制变量变化范围在定义域内 
%============================= 
function Position = limitToPosition(Position,VarMin,VarMax)     
        for i =1:size(Position,2) 
            if Position(i)<VarMin 
                Position(i) = VarMin; 
            elseif Position(i) > VarMax 
                Position(i) = VarMax; 
            end 
        end 
end 

