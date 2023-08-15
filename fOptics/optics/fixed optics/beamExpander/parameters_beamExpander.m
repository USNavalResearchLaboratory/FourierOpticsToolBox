function [P] = parameters_beamExpander(wavefront)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% substrucure U, for SI units of lenth
P.U.nm = 1E-9; % nanometers
P.U.um = 1E-6; % microns
P.U.mm = 1E-3; % millimeters
P.U.cm = 1E-2; % centimeters
P.U.m  = 1E0;  % meters
P.U.km = 1E3;  % kilometers


% Lens 1
P.lens1.diameter = 25*P.U.mm;
P.lens1.focus = 1*P.U.m; 
P.lens1.N = wavefront.N; 
P.lens1.dx = P.lens1.diameter/(P.lens1.N/2);
P.lens1.wvl = wavefront.wvl;
P.lens1.position =wavefront.position ; 


% Lens 2
P.lens2.diameter = 15*P.U.mm;
P.lens2.focus = .5*P.U.m; 
P.lens2.N = wavefront.N; 
P.lens2.dx = P.lens2.diameter/(P.lens1.N/2);
P.lens2.wvl = wavefront.wvl;
P.lens2.position =P.lens1.position + P.lens1.focus + P.lens2.focus ; 



end

