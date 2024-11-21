%============================= 
%删除档案库中的一个个体 
%============================= 
function rep = DeleteOneRepMemebr(rep,gamma) 
        GI = [rep.GridIndex]; 
        OC = unique(GI); 
        N = zeros(size(OC)); 
        for k = 1:numel(OC) 
            N(k) = numel(find(GI == OC(k))); 
        end 
        P = exp(gamma*N); 
        P = P/sum(P); 
        sci = RouletteWheelSelection(P); 
        sc = OC(sci); 
        SCM = find(GI == sc); 
        smi = randi([1 numel(SCM)]); 
        sm = SCM(smi); 
        rep(sm) = []; 
end 

