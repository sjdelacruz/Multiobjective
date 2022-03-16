%Fast non dominated sorting
%Based on the NSGAII

function [Population, F] = fast_nondominated_sorting(Population)
    nPop = numel(Population);

    %Reset values of population
    for i = 1:nPop
        Population(i).DominationSet = [];
        Population(i).DominatedCount = 0;
    end
    
    F{1} = [];
   
    for i = 1:nPop
        for j = i+1:nPop
            p = Population(i);
            q = Population(j);
            
            if dominates(p, q)
                p.DominationSet = [p.DominationSet j];
                q.DominatedCount = q.DominatedCount+1;
            end
            
            if dominates(q, p)
                q.DominationSet = [q.DominationSet i];
                p.DominatedCount = p.DominatedCount+1;
            end
            
            Population(i) = p;
            Population(j) = q;
        end
        
        if Population(i).DominatedCount == 0
            F{1} = [F{1} i];
            Population(i).Rank = 1;
        end
    end
    
    k = 1;
    
    while true
        
        Q = [];
        
        for i = F{k}
            p = Population(i);
            
            for j = p.DominationSet
                q = Population(j);
                
                q.DominatedCount = q.DominatedCount-1;
                
                if q.DominatedCount == 0
                    Q = [Q j]; %#ok
                    q.Rank = k+1;
                end
                
                Population(j) = q;
            end
        end
        
        if isempty(Q)
            break;
        end
        
        F{k+1} = Q; %#ok
        
        k = k+1;
        
    end
    

end