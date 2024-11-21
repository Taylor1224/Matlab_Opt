%============================= 
%从全局支配个体中找出一个最佳个体 
%============================= 
function leader = SelectLeader(rep,beta) 
        GI = [rep.GridIndex]; 
        OC = unique(GI); 
        %一个栅格可能被多个支配解占用 
        N = zeros(size(OC)); 
        for k =1:numel(OC) 
            N(k) = numel(find(GI == OC(k))); 
        end 
        % 计算选择概率,为了增加多样性,尽量不选多次出现的个体 
        % 如果N大P就小, 即多次出现的栅格点被选中的概率小 
        P = exp(-beta*N); 
        P = P/sum(P);  
        sci = RouletteWheelSelection(P);   %轮盘赌策略选择 
        sc = OC(sci);   % 轮盘赌选择的栅格点 
        SCM = find(GI==sc); 
        smi = randi([1 numel(SCM)]); 
        sm = SCM(smi); 
        leader = rep(sm);   %当前全局最佳位置点
end 

