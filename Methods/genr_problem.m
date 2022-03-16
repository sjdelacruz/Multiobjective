%Function to generate a structure with the information of the desired problem 
function Problem = genr_problem(name)
if name == "Kursawe"
   Problem.M = 2;
   Problem.D = 3;
   Problem.C = 0;
   Problem.Lowers = repelem(-5,Problem.D);
   Problem.Uppers = repelem(5,Problem.D);
   Problem.name= str2func(name);
   
else 
    null;
end
end
