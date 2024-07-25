function [solution,objectiveValue,reasonSolverStopped] = optimize_sampling(lambda, deltaZ, Dn, delta10, D1, deltaN0, N0)
    % OPTIMIZE_SAMPLING - Optimization created using MATLAB optimization toolbox. 
    % Minimizes the dx1 (delta1), dx2 (also called deltaN), and sampling points N. 
    % Things you can change:
    %  Lower and upper bounds
    % Max iterations, optimization algorithm (current patternsearch algorithm), Tolerences 
    % We use the first 2 constraints of the angular spectrum propagation constraints (see angSpecPropSampling.m)

% Create optimization variables
delta1 = optimvar("delta1","LowerBound",1e-8);
deltaN = optimvar("deltaN","LowerBound",1e-8);
N = optimvar("N","LowerBound",100,"UpperBound",4000);

% Set initial starting point for the solver
initialPoint.delta1 = delta10;
initialPoint.deltaN = deltaN0;
initialPoint.N = N0;

% Create problem
problem = optimproblem;

% Define problem objective
problem.Objective = (delta1 + deltaN)*N;

% Define problem constraints
problem.Constraints.constraint1 = deltaN <= (lambda*deltaZ-Dn*delta1)/(D1);
problem.Constraints.constraint2 = lambda*deltaZ-Dn*delta1 >= 0;
problem.Constraints.constraint3 = N >= D1/(2*delta1)+ Dn/(2*deltaN) + (lambda*deltaZ)/(2*delta1*deltaN);
problem.Constraints.constraint4 = delta1 >= 0;

% Set nondefault solver options
options = optimoptions("patternsearch","ConstraintTolerance",1e-08,...
    "MaxIterations",10000,"StepTolerance",1e-08);

% Display problem information
%show(problem);

% Solve problem
[solution,objectiveValue,reasonSolverStopped] = solve(problem,initialPoint,...
    "Solver","patternsearch","Options",options);

% Display results
%solution
%reasonSolverStopped
%objectiveValue