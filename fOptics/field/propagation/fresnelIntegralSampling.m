function [delta1,N,Zconstraint, dx2 ] = fresnelIntegralSampling(D1, lambda, z, delta1, D2, N, R)
%fresnelIntegralSampling - Fresnel Integral sampling requirements
%
% For accurate simulations, constraints must be satisfied. The two
% constraints for the Fresnel Diffraction Integral are given in Ref. [1].
% In general, the Fresnel integral is valid for long propagations and the
% Angular-spectrum method for short propagations. 
%
% Syntax: 
%   delta2 = fresnelIntegralSampling(D1, lambda, z, delta1, D2, N, R);
%
% Inputs: (units in meters)
%   D1 - diameters of source plane
%   lambda - wavelength of the source
%   z - distance from source to observation
%   delta1 - grid spacing at the source plane
%   D2 - diameter of observation area
%   N - number of grid points
%   R - radius of curvature of source
%
% Outputs:
%    delta2 - grid spacing at the observation plane
% 
% Example:
%     z = 7500; 
%     D1 = 300E-3; 
%     D2 = 2; 
%     lambda = 1.064E-6;
% 
%     Ngrid = 512; 
%     Nsource = 128; 
% 
%     delta1 = D1/Nsource; 
% 
%     R = inf; 
%     delta2 = fresnelIntegralSampling(D1, lambda, z, delta1, D2, Ngrid, R);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% References: [1] "Numerical Simulation of Optical Wave Propagation" by
%                 Jason D. Schmidt, SPIE Press, Bellingham, WA 2010. pg.
%                 120-123
%
% See also: fresnelDiffractionIntegral, angularSpectrumPropagation

% Author: Dennis F. Gardner Jr.
%         Code 5661
%         4555 Overlook Ave SW
%         Washington, DC 20375
% Email: dennis.gardner@nrl.navy.mil
% Website: https://www.linkedin.com/in/dennisfgardner/
% March 2018; Last revision: 21-March-2018



dx1_min = lambda*z/(2*D2); 
% if delta1 > dx1_min
%     delta1 = dx1_min; 
% end
dx2 = (lambda * z - D2 * delta1) / D1;

LogicalStr = {'false', 'true'};

% Equation 7.31 from Ref [1]

Nconstraint = (D1*lambda*z)/(delta1*(lambda*z-D2*delta1));
N = 2^ceil(log2(Nconstraint));
%= nx /(1- D2*delta1/lambda/z)

NposTF = Nconstraint > 0;
NsampTF = N>= Nconstraint;


fprintf('Is the constraint on N (%.0f) larger than 0?: %s\n', Nconstraint, LogicalStr{NposTF + 1 }); 
fprintf('Is N (%.0f) larger than the constraint: %s\n', N, LogicalStr{NsampTF + 1 });

% check the condition on Z

if R == inf
    Zconstraint = D1*delta1/lambda; 
    
else
    Zconstraint = D1*delta1*R/(lambda*R - D1*delta1);
    
end

ZTF = z >= Zconstraint;

fprintf('Is z (%.3f) larger than the constraint (%.3f)?: %s\n\n', z, Zconstraint, LogicalStr{ZTF + 1 });


end

