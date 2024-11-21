%============================= 
%栅格索引定位 
%============================= 
function particle = FindGridIndex(particle,Grid) 
        nObj = numel(particle.Cost); 
        nGrid = numel(Grid(1).LB); 
        particle.GridSubIndex = zeros(1,nGrid); 
        for j = 1:nObj 
            particle.GridSubIndex(j) = find(particle.Cost(j)<=Grid(j).UB,1,'first'); 
            %从左到右找到第一个目标值小于栅格值的位置 
        end 
        particle.GridIndex = particle.GridSubIndex(1); 
        for j = 2:nObj   % 左上角开始数到右下角，先数行再换行继续数 
            particle.GridIndex = particle.GridIndex-1; 
            particle.GridIndex = nGrid*particle.GridIndex; 
            particle.GridIndex = particle.GridIndex + particle.GridSubIndex(j); 
        end 
end

