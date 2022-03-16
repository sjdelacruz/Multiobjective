function [Result] = SMS_EMOA(nameProblem,maxGen,N,Save)

Problem = genr_problem(nameProblem);
Population = genr_population(N, Problem);
Convergence = [];
Gen = 0;

while Gen <= maxGen
    
    %select random parents
    parentsIndex = randperm(N,2);
    
    %Apply SBX and PM, and evaluate the offspring in the problem, select
    %individual 1 or 2 to mutate
    Offspring_xs = generate(Population(parentsIndex(1)).x,Population(parentsIndex(2)).x,30,1/Problem.D,Problem.D,Problem.Lowers,Problem.Uppers);
    Offspring_ys = Problem.name(Offspring_xs, 1);
    
    %Add Offspring to population "according to SMS_EMOA's steps" 
    aux1 = numel(Population);
    Population(aux1+1).x = Offspring_xs;
    Population(aux1+1).y = Offspring_ys;
    
    %%Fast-non-dominated-sorting 
    Population = reduce(Population,Problem.M);
    
    if  mod(Gen,Save) == 0 
        Convergence = [Convergence Population];
    end
    
    Gen=Gen+1;
    
    
end
    [Population,F] = fast_nondominated_sorting(Population);
    Indexs = F{1};
    n = length(Indexs);
    
    OptimalFront = zeros(n,Problem.M, 'double');
    for i=1: n
        OptimalFront(i,:) = Population(Indexs(1,i)).y;
    end
    
    Result.population = Population;
    Result.optimalFront = OptimalFront;
    Result.convergence = Convergence;
end

function Population = reduce(Population,M)
   
    [Population,F] = fast_nondominated_sorting(Population);
    
    %Identify the solutions in the worst ranked front
    Indexs = F{end};
    n = length(Indexs);
    WorstFrontObjs = zeros(n,M, 'double');
    for i=1: n
        WorstFrontObjs(i,:) = Population(Indexs(1,i)).y;
    end
    
    %Calculate the contribution of hypervolume of each solution
    deltaS = inf(1,n);
    if M == 2
        [~,rank] = sortrows(WorstFrontObjs);
        for i = 2 : n-1
            deltaS(rank(i)) = (WorstFrontObjs(rank(i+1),1)-WorstFrontObjs(rank(i),1)).*(WorstFrontObjs(rank(i-1),2)-WorstFrontObjs(rank(i),2));
        end
    elseif n > 1
        deltaS = CalHV(WorstFrontObjs,max(WorstFrontObjs,[],1)*1.1,1,10000);
    end
    
    % Delete the worst solution from the last front
    [~,worst] = min(deltaS);
    Population(Indexs(worst)) = [];
end

function F = CalHV(points,bounds,k,nSample)
    [N,M] = size(points);
    if M > 2
        % Use the estimated method for three or more objectives
        alpha = zeros(1,N); 
        for i = 1 : k 
            alpha(i) = prod((k-[1:i-1])./(N-[1:i-1]))./i; 
        end
        Fmin = min(points,[],1);
        S    = unifrnd(repmat(Fmin,nSample,1),repmat(bounds,nSample,1));
        PdS  = false(N,nSample);
        dS   = zeros(1,nSample);
        for i = 1 : N
            x        = sum(repmat(points(i,:),nSample,1)-S<=0,2) == M;
            PdS(i,x) = true;
            dS(x)    = dS(x) + 1;
        end
        F = zeros(1,N);
        for i = 1 : N
            F(i) = sum(alpha(dS(PdS(i,:))));
        end
        F = F.*prod(bounds-Fmin)/nSample;
    else
        % Use the accurate method for two objectives
        pvec  = 1:size(points,1);
        alpha = zeros(1,k);
        for i = 1 : k 
            j = 1 : i-1; 
            alpha(i) = prod((k-j)./(N-j))./i;
        end
        F = hypesub(N,points,M,bounds,pvec,alpha,k);
    end
end

function h = hypesub(l,A,M,bounds,pvec,alpha,k)
% The recursive function for the accurate method

    h     = zeros(1,l); 
    [S,i] = sortrows(A,M); 
    pvec  = pvec(i); 
    for i = 1 : size(S,1) 
        if i < size(S,1) 
            extrusion = S(i+1,M) - S(i,M); 
        else
            extrusion = bounds(M) - S(i,M);
        end
        if M == 1
            if i > k
                break; 
            end
            if alpha >= 0
                h(pvec(1:i)) = h(pvec(1:i)) + extrusion*alpha(i); 
            end
        elseif extrusion > 0
            h = h + extrusion*hypesub(l,S(1:i,:),M-1,bounds,pvec(1:i),alpha,k); 
        end
    end
end





%Function generate to produce only one offspring for sms_emoa
function Offspring = generate(Parent1, Parent2,pDistCross,pDistMut,D,Lowers,Uppers)

    [Processed1, Processed2] = simulated_binary_crossover(Parent1, Parent2,pDistCross,D);
    r = randi([1 2],1,1);
    if r==1
        Processed = Processed1;
    else
        Processed = Processed2;
    end
    Offspring = polynomial_mutation(Processed,pDistMut,D,Lowers,Uppers);
end







%Function to create a population with uniform distribution
%Population.cs to constraints
function Population = genr_population(N, Problem)

Ind.x = [];
Ind.y = [];
Ind.Rank = [];
Ind.DominationSet = [];
Ind.DominatedCount = 0;
Ind.HV=0.0;

Population = repmat(Ind, N, 1);

%Create design variables
xs=zeros(N,Problem.D, 'double');
for i=1:Problem.D
    %One because i need only column due to the iteration
    xs(:,i) = unifrnd(Problem.Lowers(i), Problem.Uppers(i),[N,1]); 
end

%Evaluating solution in the multiobjective problem 
[ys,~] = Problem.name(xs, N);

%Creating structure of Individuals
for i=1:N
    Population(i).x = xs(i,:);
    Population(i).y = ys(i,:);
end
end





