nombre= '09-Mar-2022_SMS_EMOA.mat';
load(nombre);
Saved=1000;
Weight_Vector= [0.5088,0.4911];

%Pareto Fronts Approximations
ParetoFronts = getFrentes(Results);

%TrueParetoFront
archivo = fopen('Kursawe.txt','r');
format long
cell_data= textscan(archivo,'%f%f','Delimiter','\t');
TrueParetoFront = cat(2,cell_data{:});

%Compute the mean of IGD values
[gdMean, gdStd, bestIndex] = GDMean(ParetoFronts,TrueParetoFront);
[spMean, spStd, ~] = SpacingMean(ParetoFronts);

%Generate plots based on the best
Desired= bestIndex;
DesiredParetoFront = ParetoFronts{1,Desired};
generarFigura(DesiredParetoFront);

%GD
gd = GD(DesiredParetoFront,TrueParetoFront);

%Spacing
Sp = Spacing(DesiredParetoFront);

%Convergence graph
PFConvergence= Results{1,Desired}.convergence;
generarConvergencia(PFConvergence,TrueParetoFront,Saved);

%Enfoque de pseudo vectores de peso
Preference = Pseudo_Weight_Vector(DesiredParetoFront,Weight_Vector);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function generarFigura(ParetoFront)
    scatter(ParetoFront(:,1), ParetoFront(:,2))
    title("Kursawe") 
    xlabel('f1') 
    ylabel('f2') 
    saveas(gcf,'SMS-EMOA.png')
end

function [ParetoFronts] = getFrentes(Results)

[~,m] = size(Results);
ParetoFronts = {};

for i=1 : m
    ParetoFronts{i}= Results{i}.optimalFront;
end
end

function generarConvergencia(PFConvergence,TrueParetoFront,Saved)
[~,m] = size(PFConvergence);

IGD=[];
Gens=[];
Counter=0;
for i=1:m
    PF = configureData(PFConvergence(:,i));
    IGD = [IGD; GD(PF, TrueParetoFront)];
        Counter = Counter + Saved;
        Gens = [Gens; Counter];
end

plot(Gens(:,1),IGD(:,1));
title("GD convergence");
xlabel('Evaluations') 
ylabel('GD value') 
saveas(gcf,'GD.png')
end

function OptimalFront = configureData(Population)
[Population,F] = fast_nondominated_sorting(Population);
Indexs = F{1};
n = length(Indexs);

OptimalFront = zeros(n,2, 'double');
for i=1: n
    OptimalFront(i,:) = Population(Indexs(1,i)).y;
end
end

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

function answer = dominates(Solution1,Solution2)
answer = all(Solution1.y <= Solution2.y) && any(Solution1.y<Solution2.y);
end

function [gdMean,gdStd,best] = GDMean(ParetoFronts, TrueParetoFront)
[~,m] = size(ParetoFronts);

gd=[];
for i= 1 : m
    gd= [gd GD(ParetoFronts{1,i}, TrueParetoFront)]; 
end
gdMean = mean(gd);
gdStd = std(gd);
[~,indexs] = sort(gd);
best = (indexs(1));
end

function [gdMean,gdStd,best] = SpacingMean(ParetoFronts)
[~,m] = size(ParetoFronts);

gd=[];
for i= 1 : m
    gd= [gd Spacing(ParetoFronts{1,i})]; 
end
gdMean = mean(gd);
gdStd = std(gd);
[~,indexs] = sort(gd);
best = (indexs(1));
end