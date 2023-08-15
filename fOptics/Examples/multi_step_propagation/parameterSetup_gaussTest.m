function [ P ] = parameterSetup_gaussTest(Cn2value)
%parameterSetup - define the parameters of the simulation
%
% This function contains all the parameters of the simulation which are in
% the structure P. Within P, there are substructures for units (U), source
% (S), receiver (R), target (T), and turbulence (Turb), SLM (SLM)
% Due to reusing this file many times, not all parameters are used in all simulations.
% Code is kept for future use. Feel free to delete unused parameters in
% your simulation. 
%
% Syntax: [ P ] = parameterSetup(Cn2value);
%
% Inputs:
%   Cn2value - the refractive index structure constant value
%
% Outputs:
%    P - structure, with substructures, containing simulation parameters
% 
% Example:
%    [ P ] = parameterSetup(1E-15);



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
P.wvl = 1000*P.U.nm;
% Kerr index of air [m^2/W]
P.n20 = 3.0*10^-23;

% distance to target
P.Z = 100;


%% Simulation grid points (N) 

P.N = 1024; %2^12; %Starting simulation size. May change with optimization

%% Target Parameters - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% three bars or USAF
P.Target.targetType = 'ptSrc';
P.Target.D = 0.0508/16; %200*P.U.mm; % diameter of target
P.Target.N = P.N; %Define target grid points. But maybe expanded to fit propagation constraints
P.Target.wvl = P.wvl;
P.Target.position = 0; %Define position of target as zero
P.Target.distance = P.Z; 
P.Target.dx = 0.0508/(P.N-1); %Gives 10 grid points in target diameter
P.Target.specks = 1; % number of speckle realizations

%% Turbulence parameters and statistics - - - - - - - - - - - - - - - - - -%
% turbulence,  (Cn2), is for plane wave), number of screens, and (deltaZ)

% Air parameters
P.air.screens = 2; 
P.air.N = P.N;
P.air.dx = P.Target.dx;
P.air.wvl = P.wvl; 

% OAM Parameters
P.oam.z = P.Z; 


% number of phase screens
P.Turb.screens = 2; %10
% distance between screens
P.Turb.deltaZ = P.Z/(P.Turb.screens);
% refractive index structure constant
P.Turb.Cn2 = Cn2value;
% number of sub-harmonics (for kolmogorov phase screen generation)
P.Turb.K = 15;
% Fried parameter (r0 full path, r0delta each screen)
P.Turb.r0 = (0.423*(2*pi/P.wvl)^2*P.Turb.Cn2.*P.Z).^(-3/5); %Using Z (total propagation distance)
P.Turb.r0delta = (0.423*(2*pi/P.wvl)^2*P.Turb.Cn2.*(P.Z/P.Turb.screens)).^(-3/5); %Using Z = dZ (distance between screens)
% isoplanatic angle at source and calculated patch size at target
P.Turb.isoAngle = (1.09*(2*pi/P.wvl)^2*P.Turb.Cn2.*P.Z^(8/3)).^(-3/5);
P.Turb.isoPatch = P.Z*tan(P.Turb.isoAngle); 

% Rytov log-amplitude variances - - - - - - - - - - - - - - - - - - - - - %
% plane wave
P.Turb.rytov =  1.23*(2*pi/P.wvl)^(7/6)*P.Turb.Cn2.*P.Z^(11/6)/4;
P.Turb.rytovDelta = 1.23*(2*pi/P.wvl)^(7/6)*P.Turb.Cn2.*P.Turb.deltaZ^(11/6)/4;

P.Turb.dx = P.Target.dx;
% disp('Turbulence Information'); 
% display(P.Turb); 

%% Optics Elements and Parameters - - - - - - - - - - -  - - - - - - - - - -%
% Pupil lens %
P.lens1.focus = P.Z; %focal length of imaging lens
P.lens1.diameter = 0.0508;  %diameter of imaging lens
P.lens1.position = P.Z; %distance between imaging lens and bifocal
P.lens1.dx = P.Target.dx; P.Turb.dx = P.lens1.dx;
P.lens1.N = P.N; 
P.lens1.wvl = P.wvl; 


end

