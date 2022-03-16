%Function to crossover two parents solutions
%Warning with the bounds of beta
%pDist= nc
function [Processed1, Processed2] = simulated_binary_crossover(Parent1, Parent2,pDist,D)

Processed1 = zeros(1,D, 'double');
Processed2 = zeros(1,D, 'double');

u = rand(1,D);
for i=1:D
    if u(i) <= 0.5
        betai = (2*u(i))^(1/(pDist+1));
    else
        betai = (1/(2*(1-u(i))))^(1/(pDist+1)); 
    end
    Processed1(i) = 0.5 * ( ((1 + betai)*Parent1(i)) + ((1-betai)*Parent2(i)) );
    Processed2(i) = 0.5 * ( ((1 - betai)*Parent1(i)) + ((1+betai)*Parent2(i)) );
end
end