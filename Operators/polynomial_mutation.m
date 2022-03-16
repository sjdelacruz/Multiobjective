%Function to polynomial mutate two recombined parents
function [Offspring] = polynomial_mutation(Processed1,pDistMut,D,Lowers,Uppers)

Offspring = zeros(1,D, 'double');

r = rand(1,D);
for i=1:D
    if r(i) < 0.5
        Deltai = ((2*r(i))^(1/(pDistMut+1)))-1;
    else
        Deltai = 1-((2*(1-r(i)))^(1/(pDistMut+1))); 
    end
    
    Offspring(i) = Processed1(i) + ((Uppers(i)-Lowers(i)) * Deltai); 
end
end


