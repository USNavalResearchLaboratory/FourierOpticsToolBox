function [ ] = samplingCheck( P )
%samplingCheck - check the sampling of the angular spectrum propagator
%
% The numerical propagation between planes needs to be satisfied. The four
% constraints found in Ref [1] are used. 
%
% Syntax: [ ] = samplingCheck( P );
%
% Inputs:
%   P - the parameter structure
%
% Outputs:
%    none
% 
% Example:
%    [ P ] = parameterSetup(1E-15);
%     [ ] = samplingCheck( P )
%
% References (equation reference are found in [1]):
%
% [1] Jason D. Schmidt, "Numerical simulation of optical wave propagation
% with examples in MATLAB" SPIE, Bellingham 2015. pg. 127
%
% Other m-files required: fresnelIntegralSampling
% Subfunctions: none
% MAT-files required: none
%
% See also: parameterSetup, fresnelIntegralSampling

% Author: Dennis F. Gardner Jr.
%         Code 5661
%         4555 Overlook Ave SW
%         Washington, DC 20375
% Email: dennis.gardner@nrl.navy.mil
% Website: https://www.linkedin.com/in/dennisfgardner/
% March 2018; Last revision: 21-March-2018

% Source to Target for Image Sharpening routine, Fresnel Integral
% converging beam is the most restrictive
disp('Source to target Fresnel sampling check')
R = P.Z; 
fresnelIntegralSampling(P.R.D, P.wvl, P.Z, P.R.px, P.T.D, P.NoneStep, R);

% Uniform Target Illumination to Receiver, partial prop
disp('Target to Receiver angular spec sampling check')
R = inf;
angSpecPropSampling(P.wvl, P.Z, P.R.D, P.R.px, P.T.D, P.R.px, P.NmultiStep, R, P.Turb.deltaZ);

% Focused illumination to target, partial prop
disp('Source to target angular spec sampling check')
R = P.Z;
angSpecPropSampling(P.wvl, P.Z, P.T.D, P.R.px, P.S.D, P.R.px, P.NmultiStep, R, P.Turb.deltaZ);
end



