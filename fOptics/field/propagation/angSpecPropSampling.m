function samplingCheck = angSpecPropSampling(lambda, deltaZ, Dn, delta1, D1, deltaN, N, R, deltaZi)
%angSpecPropSampling - angular-spectrum propagation sampling requirements
%
% For accurate simulations, constraints must be satisfied. The constraints
% for angular spectrum propagation are given in Ref. [1]. It's assumed
% partial propagations are used with the angular spectrum method, with the
% appropraite windowiing function. In general, the Fresnel integral is
% valid for long propagations and the Angular-spectrum method for short
% propagations.
%
% Syntax: 
%   angSpecPropSampling(lambda, deltaZ, Dn, delta1, D1, deltaN, N, R, deltaZi);
%
% Inputs: (units in meters)
%   lambda - wavelength of the source
%   deltaZ - distance from source to observation
%   Dn - diameter of observation area (nth plane)
%   delta1 - grid spacing at the source plane
%   D1 - diameter of source plane
%   deltaN - grid spacing at observation plane (nth plane)
%   N - number of grid points
%   R - radius of curvature of source
%   deltaZi - partial propagation distance
%
% Outputs:
%    none
% 
% Example: uniformly illuminated 2m target (D1) propagated to a 300mm
%          receiver(Dn) 7.5 km away, with 10 intermediate planes
%     D1 = 2; 
%     lambda = 1.064E-6;
%     deltaZ = 7500; 
%     Dn = 300E-3; 
%     camPixels = 128;
%     deltaN = Dn/camPixels; 
%     delta1 = deltaN;
%     N= 2048; 
%     R = inf;
%     deltaZi = deltaZ/11;
% 
%     angSpecPropSampling(lambda, deltaZ, Dn, delta1, D1, deltaN, N, R, deltaZi)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% References: [1] "Numerical Simulation of Optical Wave Propagation" by
%                 Jason D. Schmidt, SPIE Press, Bellingham, WA 2010. Ch.
%                 7&8
%
% See also: angularSpectrumPropagation, fresnelDiffractionIntegral

% Author: Dennis F. Gardner Jr.
%         Code 5661
%         4555 Overlook Ave SW
%         Washington, DC 20375
% Email: dennis.gardner@nrl.navy.mil
% Website: https://www.linkedin.com/in/dennisfgardner/
% March 2018; Last revision: 22-March-2018

%disp('Checking Sampling criteria')
LogicalStr = {'false', 'true'};

% Constraint #1 
deltaNConstraint = (lambda*deltaZ-Dn*delta1)/(D1);
TFdeltaN = round(deltaN,7)<=round(deltaNConstraint,7); %round because sometimes difference is small


% Constraint #2
Nconstraint1 = floor(D1/(2*delta1)+ Dn/(2*deltaN) + (lambda*deltaZ)/(2*delta1*deltaN));
TFN1 = N>= Nconstraint1;


% Constraint #3 
LHS3 = (1+ deltaZ/R)*delta1 - (lambda*deltaZ)/D1;
RHS3 = (1+ deltaZ/R)*delta1 + (lambda*deltaZ)/D1;
TFL = LHS3 <= deltaN;
TFR = deltaN<=RHS3;
TFconstraint3 = TFL*TFR;

% Constraint #4
Nconstraint2 = (lambda*deltaZi)/(delta1*deltaN);
TFN2 = N >= Nconstraint2;

% Constraint 8.24
deltaZiConstraint = deltaN^2*N/lambda;
TFdeltaZi = deltaZi<=deltaZiConstraint; 

if TFdeltaN ~= 1 || TFN1 ~= 1 || TFconstraint3 ~= 1 || TFN2 ~= 1 || TFdeltaZi~= 1
    samplingCheck = false; 
%     fprintf('deltaN (dx2) %.3f mm needs to be <= to %.3f mm: %s\n', deltaN*1E3, deltaNConstraint*1E3, LogicalStr{TFdeltaN + 1 }); 
%     fprintf('N %.0f needs to be >= to %.0f: %s\n', N, Nconstraint1, LogicalStr{TFN1 + 1 }); 
%     fprintf('%.3f mm <= %.3f mm <= %.3f mm: %s \n', LHS3*1E3, deltaN*1E3, RHS3*1E3, LogicalStr{TFconstraint3 + 1 }); 
%     fprintf('N %.0f needs to be >= to %.0f: %s\n', N, Nconstraint2, LogicalStr{TFN2 + 1 }); 
%     fprintf('deltaZi %.3f needs to be <= to %.3f: %s : Add more phase screens if False \n\n ', deltaZi, deltaZiConstraint, LogicalStr{TFdeltaZi + 1 }); 

else 
    samplingCheck = true;

end

