function [ind_objs,constraints_values] = Kursawe(ind_vars,N)

f1=zeros(N,1, 'double');
for i= 1 : 2
  f1(:,1) = f1(:,1) + (-10.^(-0.2 * sqrt(ind_vars(:,i).^2 + ind_vars(:,i+1).^2)));  
end

f2= zeros(N,1, 'double');
for j=1 : 3
    f2(:,1)= f2(:,1) + (abs(ind_vars(:,j)).^(0.8) + 5 * sin(ind_vars(:,j).^3)); 
end
ind_objs = [f1,f2];

%Omitted error(problem doesnt have constriants)
constraints_values=[];
end

