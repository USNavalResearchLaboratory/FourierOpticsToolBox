% Script to test Shack Hartmann. 
close all; clear; clc

% Inputs

[ P ] = parameterSetup_ptsrc(2e-15,1000); 
gpuFlag = 1; %use gpu? 

% Generate turbulence/ phase screen
turb = phase_screens(P.air); %Creates initial flat screens
air = turb; 
turb = turb.kolmogorov(P.Turb, gpuFlag); %Generates turbulence. try turb.tilt, or turb.OAM
%turb = turb.tilt(0,0.04);
%turb = turb.OAM(2); 
turb.showme

%% Setup target
target = target_class(P.Target); % Define initial target
lens1 = thinLens(P.lens1); 
target = target.genTarget_focus_spot(lens1.conjugate, air, gpuFlag); 
%target = target.genTarget_gauss(); % Other point source options
%target = target.genTarget_ptsrc(); 
target.showme
%%
SH_ref = singleShackHartmann(P, target,air, gpuFlag);
SH_signal = singleShackHartmann(P, target,turb, gpuFlag);

%% phase reconstruction
SH = SH_signal.calculateSlopes(SH_ref);
SH.showSlopes
SH = SH.phaseReconstruction(gpuFlag);
SH.showme('reconstructed wavefront')
phase = zonal_phase_reconstruction(SH);
figure; imagesc((phase))


%% Call singleShackHartmann twice for Ref and Signal beams
function [SH] = singleShackHartmann(P, target, turb, gpuFlag)
    %% single shack hartmann datagen. 
    
    U = Efield(target, gpuFlag);   
    % Generate  lenses
    lens1 = thinLens(P.lens1);
    %lens1.showme('Pupil lens'); 
    U2 = U.propagateTo(lens1,'fwd',turb, 'angspec');
    U2.showme('Field at Lens 1')
    %U_rev = U2.propagate_ang_spec(P.Z, 'rev', target.width,turb); U_rev.showme('Backpropagated light')
    %% Beam expander to shrink beam
    
    %U3 = beamExpander(U2); 
    
    U3 = U2.magnify(P.magnification);
    %U3 = U2; 
    U3.showme('After beam shrinker')
    %U3 = U3.propagate_ang_spec(1e-2);
    %U3.showme('propagated a bit')
    %% Shack Hartmann
    
    sensorParameters = SHWFSParametersDefault(U3); 
    %Create ShackHartmann object
    SH = shackHartmannWavefrontSensor(sensorParameters);
    SH.showme
    %SH.lenslet.showme   

    SH = SH.applyShackHartmannFull(U3);

end


