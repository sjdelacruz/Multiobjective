%Define the seed
rng('default');

%Parameters of executions
Executions = 1;
Problem = "Kursawe";
N=100;
D =10;
M=3;
Q=2;
maxGen=20000;

Save=500;

Results = {};

Execute(SMS_EMOA,Problem, maxGen, N,Save);
for i = 1:Executions
    rng(i);
    Result = SMS_EMOA(Problem, maxGen, N,Save);
    Results = [Results Result];
end

Preference = Pseudo_Weight_Vector(Results{1}.optimalFront,[0.75,0.25]);
d = datetime('today');
c = datestr(d) + "_SMS_EMOA.mat";
save(convertStringsToChars(c));