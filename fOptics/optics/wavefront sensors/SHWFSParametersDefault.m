function [P] = SHWFSParametersDefault(wavefront)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% substrucure U, for SI units of lenth
P.U.nm = 1E-9; % nanometers
P.U.um = 1E-6; % microns
P.U.mm = 1E-3; % millimeters
P.U.cm = 1E-2; % centimeters
P.U.m  = 1E0;  % meters
P.U.km = 1E3;  % kilometers


%Lenslet array parameters
P.array.num = 50; 
P.array.spacing = 0; %spacing between lenslets
P.array.shape = 'square'; %square or hex configuration


% Single Lens parameters
P.lenslet.diameter = 150*P.U.um;
P.lenslet.focus = 5*P.U.mm; 
P.lenslet.N = wavefront.N/(P.array.num); 
P.lenslet.dx = P.lenslet.diameter/P.lenslet.N;
P.lenslet.wvl = 1000*P.U.nm;
P.lenslet.position =wavefront.position; 


P.array.diameter = (P.lenslet.diameter + P.array.spacing)*P.array.num; % width of lenslet array
P.array.N = wavefront.N; 

end

