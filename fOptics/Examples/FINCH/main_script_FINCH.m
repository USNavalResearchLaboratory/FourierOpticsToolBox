% Main script to simulate FINCH Interferometry 
%close all; clear; clc
%close all
%% Inputs
Cn2 = 5e-14;
targetDistance = 500; %
[ P ] = parameterSetup_targetAndTurb(Cn2,targetDistance); 
flags.gpu = 0; %Use gpu
flags.speckle =0; 
% Generate turbulence phase screen
air = phase_screens(P.air); %Creates initial flat screen

if ~nOrder
    nOrder = 7; 
end
%turb = air.kolmogorov(P.Turb,flags.gpu); %Generates turbulence. try turb.tilt, or turb.OAM
flags.randomize_zernikes = true; scale = 1; 
turb = air.zernike(nOrder,12,flags.randomize_zernikes,scale); 
turb.zernike_distribution = turb.zernike_distribution./ scale;
turb.zernike_profile
turb.showme


%% Setup
% Initialize target
target = target_class(P.Target); % Define initial target class
% Generate  lenses
lens1 = thinLens(P.lens1);
% Choose target 
target = target.genTarget_focus_spot(lens1.conjugate, air, flags.gpu);
%target = target.genTarget_simplegauss(P.Target.width);
% Add speckle if chosen
if flags.speckle == 1 
    smoothing_factor = 4; 
    target = target.addSpeckle(smoothing_factor); 
end 
target.showme()

%% Propagate to lens1

U_target = Efield(target,flags.gpu);   

%lens1.showme('Pupil lens'); 
U_pupilPlane = U_target.propagateTo(lens1,'fwd',turb);
% U_pupilPlane.showme('Field at Lens 1')
U_pupilPlaneREF = U_target.propagateTo(lens1,'fwd',air);


%% Pass into Optical setup. Output is hologram after phase shift process

[U_pupilPlane_dR] = FINCH_opticalSetup(U_pupilPlane); 
[U_pupilPlane_dRREF] = FINCH_opticalSetup(U_pupilPlaneREF); 
U_truth = U_pupilPlane; 
U_truth.data = U_pupilPlane.data./U_pupilPlaneREF.data;
U_truth.data(isnan(U_truth.data)) = 0; 
U_truth.showme('Phase at Pupil plane (TRUTH)')


