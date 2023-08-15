% Fixed parameters
D1 = 0.01; % m
Dn = 0.01; % m
lambda = 0.6328e-6; % m
deltaZ = 0.02; % m

% Initial guess
N0 = 50; % initial value for N
deltaN0 = 1e-6; % initial value for deltaN
delta10 = 1e-6; % initial value for delta1
x0 = [N0, deltaN0, delta10];

% Create optimization variables
delta12 = optimvar("delta1","LowerBound",1e-8);
deltaN2 = optimvar("deltaN","LowerBound",1e-8);
N2 = optimvar("N","LowerBound",100,"UpperBound",5000);

% Set initial starting point for the solver
initialPoint2.delta1 = delta10;
initialPoint2.deltaN = deltaN0;
initialPoint2.N = N0;

% Create problem
problem = optimproblem;

% Define problem objective
problem.Objective = (delta12 + deltaN2)*N2;

% Define problem constraints
problem.Constraints.constraint1 = deltaN2 <= (lambda*deltaZ-Dn*delta12)/(D1);
problem.Constraints.constraint2 = lambda*deltaZ-Dn*delta12 >= 0;
problem.Constraints.constraint3 = N2 >= D1/(2*delta12)+ Dn/(2*deltaN2) + (lambda*deltaZ)/(2*delta12*deltaN2);

% Set nondefault solver options
options2 = optimoptions("patternsearch","ConstraintTolerance",1e-08,...
    "MaxIterations",10000,"StepTolerance",1e-08);

% Display problem information
show(problem);

% Solve problem
[solution,objectiveValue,reasonSolverStopped] = solve(problem,initialPoint2,...
    "Solver","patternsearch","Options",options2);

% Display results
solution
reasonSolverStopped
objectiveValue