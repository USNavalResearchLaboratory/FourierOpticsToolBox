function [ P ] = parameterSetup_FINCH_opticalSetup(Uin)
%parameterSetup_lateralShear - define the parameters of the optics for the
%shearing interferometer
%
% 


%% Units (U) - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% substrucure U, for SI units of lenth
P.U.nm = 1E-9; % nanometers
P.U.um = 1E-6; % microns
P.U.mm = 1E-3; % millimeters
P.U.cm = 1E-2; % centimeters
P.U.m  = 1E0;  % meters
P.U.km = 1E3;  % kilometers

% 'Global' parameters - - - - - - - - - - - - - - - - - - - - - - - - - - %
% wavelength (wvl)
P.wvl = Uin.wvl;


%% Simulation grid points (N) 

P.N = Uin.N; 
P.dx = Uin.dx; 


P.air.screens = 2; %min 2
P.air.N = P.N;
P.air.dx = P.dx;
P.air.wvl = P.wvl; 
%% Optics Elements and Parameters - - - - - - - - - - -  - - - - - - - - - -%

% SLM %
P.SLM.focus = 2*P.U.m; %was 6
P.SLM.diameter = 5*P.U.cm;  
P.SLM.position = Uin.position+0.1*P.U.m; 
P.SLM.dx = P.SLM.diameter/1000;
P.SLM.N = P.N; 
P.SLM.wvl = P.wvl;

% Camera
P.cam.position = P.SLM.position+500*P.U.mm; %position of camera
P.cam.diameter = 2*P.U.cm; 
P.cam.dx = P.SLM.dx; 
P.cam.N = P.N; 
P.cam.wvl = P.wvl; 
P.cam.focus = 0; 
end

