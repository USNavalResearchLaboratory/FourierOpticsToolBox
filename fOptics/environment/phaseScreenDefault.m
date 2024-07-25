function [P] = phaseScreenDefault(wavefront)
%phaseScreenDefault Create default parameters for phase screen
%   Air phase screens 

% Air parameters
P.screens = 2; %min 2
P.N = wavefront.N;
P.dx = wavefront.dx;
P.wvl = wavefront.wvl; 



end