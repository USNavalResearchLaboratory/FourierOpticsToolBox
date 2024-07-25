%% Adaptive optics loop
% Run SH example to get some data
shack_hartmann_example; 


%% Complete AO loop with backpropagation 
U_return = Efield(SH);
U_return.data = ones(size(SH.data)).*SH.data;
%U_return = SH_signal.incident_field; % Use this for perfect reconstruction
U_return.showme

% Magnify return beam 
lens1 = thinLens(P.lens1); %lens1.data = conj(lens1.data); 
U_return = U_return.resample(lens1); 
    % Demagnification was 0.1. So I magnify by 10 here
U_return = U_return.magnify(1/P.magnification); U_return.position = lens1.position;
U_return = U_return.apply_optic(lens1.conjugate); 

%% Backpropagate 
Image = target_class(P.Target); % Define initial target
Image = Image.genTarget_focus_spot(U_return, turb, gpuFlag); 
Image.showme('Corrected Image')
Image_ref = target_class(P.Target); % Define initial target
Image_ref = Image_ref.genTarget_focus_spot(lens1.conjugate, turb, gpuFlag); 
Image_ref.showme('Uncorrected Image')

