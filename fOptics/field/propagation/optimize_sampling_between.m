function [obj1, dx2,aberration] = optimize_sampling_between(obj1,D1,Dn,dx2,deltaZ,aberration)
%OPTIMIZE_SAMPLING function optimizes the sampling between the two propagation planes. Finds minimum dx1, dx2, and N. 
% This function calls optimize_sampling and then resamples the field at the first plane to fit the new parameters. 
    % Inputs:

    % obj1: Efield or optic object. 
    % D1: the diameter field at first plane
    % Dn: the diameter field at second plane
    % deltaN0: initial guess for the sampling at second plane
    % deltaZ: propagation distance
% Outputs:

    % obj1: the resampled field at first plane
    % dx2: the sampling interval after optimization

disp('Sampling Constraints Failed. Trying Optimization')

lambda = obj1.wvl;
N0 = obj1.N;
delta10 = obj1.dx;  

% Find optimum values
[solution,objectiveValue,reasonSolverStopped] = optimize_sampling(lambda, deltaZ, Dn, delta10, D1, dx2, N0);
% Create new field1 we should fit to
obj1_new = obj1;

obj1_new.N = solution.N;
obj1_new.dx = solution.delta1;
[obj1_new.X, obj1_new.Y] = meshgrid((-obj1_new.N/2 : obj1_new.N/2-1) * obj1_new.dx);
obj1_new.roi = obj1_new.dx * obj1_new.N; 

% Use resample on obj1
if reasonSolverStopped == -2
    disp('Not resampling if Solver did not find good solution')
else
    disp('Resampling to new parameters')
    obj1 = obj1.resample2(obj1_new);
    %Resample aberration 
    aberration = aberration.resample(obj1);
    dx2 = solution.deltaN; 
end



end