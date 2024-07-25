% function to for shearing interferometer 

function [Uholo] = FINCH_opticalSetup(Uin)

P = parameterSetup_FINCH_opticalSetup(Uin); 
air = phase_screens(P.air); 
%% Shrink Beam
magnification = 0.2; 
Uin = Uin.magnify(magnification);
%% Create Ref and Signal
U_ref = Uin; 
U_signal = U_ref; 



%% Propagate to next optic
SLM = thinLens(P.SLM);
%SLM.showme('SLM'); 
SLMFlat = SLM; SLMFlat.data = ones(SLM.N);
U_signal = U_signal.propagateTo(SLM,'fwd',air, 'angspec');
U_ref = U_ref.propagateTo(SLMFlat,'fwd',air, 'angspec');
%% Propagate to Cam
% z_cam = P.cam.position;
% U_signal = U_signal.propagate_ang_spec(z_cam);  
% U_ref = U_ref.propagate_ang_spec(z_cam);
%U_signal.showme; U_ref.showme
camera = thinLens(P.cam); camera = camera.pinhole(P.cam);
U_signal = U_signal.propagateTo(camera,'fwd',air, 'angspec');
U_ref = U_ref.propagateTo(camera,'fwd',air, 'angspec');

Uholo = hologram_detector('FINCH', U_signal); 
Uholo = Uholo.phase_shift_holography(U_ref,U_signal);
Uholo.showme('hologram at camera')
Uholo = Uholo.finch_hologram();

end