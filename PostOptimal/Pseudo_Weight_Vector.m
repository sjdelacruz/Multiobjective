function Solution = Pseudo_Weight_Vector(ParetoFront,Weight_Vector)

Weight_Vector= round(Weight_Vector,2,'decimals');

[n,m] = size(ParetoFront);

Weights  =zeros(n,m,'double');
ideal = zeros(1,m,'double'); 
nadir = zeros(1,m,'double');

%Getting the f^max and f^min of all the objective functions
for i=1: m
    ideal(1,i) = min(ParetoFront(:,i));
    nadir(1,i) = max(ParetoFront(:,i));
end

for i=1:n
    a= (nadir - ParetoFront(i,:));
    b= (nadir - ideal);
    den = sum(a./ b,2);
    
    for j=1:m
        nom = a(1,j)/ b(1,j) ;
        Weights(i,j) = nom/den;
    end
end

[Weights,rank] = sortrows(Weights,'descend');
Weights = round(Weights,2,'decimals');
for i=1 : n
    if isequal(Weight_Vector(1,:),Weights(i,:))
        Solution = ParetoFront(rank(i),:);
        return;
    end
end

Solution=[Weights,rank];
end

