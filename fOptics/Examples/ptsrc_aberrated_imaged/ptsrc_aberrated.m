% Top level script for point source propagated through atmosphere and imaged
% on screen
close all; clear; clc

% Inputs
[ P ] = parameterSetup_ptsrc(8e-13,100); 
gpuFlag = 1; %use gpu? 
% Generate turbulence/ phase screen
turb = phase_screens(P.air); %Creates initial flat screens
air = turb; 
turb = turb.kolmogorov(P.Turb,gpuFlag); %Generates turbulence. try turb.tilt, or turb.OAM
%turb = turb.tilt(0.001); 
%turb = turb.OAM(1);
turb.showme

% Setup target
target = target_class(P.Target); % Define initial target
lens1 = thinLens(P.lens1); 
%target = target.genTarget_focus_spot(lens1.conjugate, air, gpuFlag); 
%target = target.genTarget_ptsrc(lens1); 
target = target.genTarget_simplegauss(P.Target.width);
target.showme

% Initialize field class at target

U = Efield(target, gpuFlag);   
% Generate  lenses
lens1 = thinLens(P.lens1);
lens1.showme('Pupil lens'); 
%
%test = U.propagate_ang_spec(P.Z, 'fwd', lens1.width, turb, 1); test.showme
[U, Uall] = U.propagateTo(lens1,'fwd',turb, 'angspec', 1);
U.showme('Field at Lens 1')

%% Reverse propagate light
U_return = U; 
gain = 1; 
U_return.data = abs(lens1.data).*exp(gain*1i.*angle((U.data)));
U_return = U_return.apply_optic(lens1.conjugate); U_return.showme('Conjugated field')
Image = target_class(P.Target); % Define initial target
Image = Image.genTarget_focus_spot(U_return, turb, gpuFlag); 
Image.showme('Corrected Image')

Image_ref = target_class(P.Target); % Define initial target
Image_ref = Image_ref.genTarget_focus_spot(lens1.conjugate, turb, gpuFlag); 
Image_ref.showme('Uncorrected Image')
%%

% %% Go through optics 
% Create lens2 and air to propagate through for multistep angular spectrum
lens2 = thinLens(P.lens2); lens2.showme
air = phase_screens(P.air);

% Propagate to lens2
U3= U.propagateTo(lens2,'fwd',air, 'angspec', 0);
U3.showme('Field after lens2')

%% Here is step by step way of achieving the same result
U3B = U.propagate_ang_spec(lens2.position-U.position); U3B.showme('Field before lens2') %Propagate to lens2
U3B = U3B.resample2(lens2);U3B.showme('Resampled field') % Resample field to match sampling of lens2
U3B = U3B.apply_optic(lens2); U3B.showme('Field after lens2; Manual method') % Apply the lens2


%% Visualize setup
%Uses .data and .position properties of each object to put an image on 3D
%space. Its not perfect. I can't see how to rotate and zoom in the
%visualization. Remove 'target' input to zoom into imaging optics
visualize_optics(lens1, lens2); title('3D visualization of optics. Propagation direction is up')
visualize_optics(U, U3) ; title('3D visualization of intensity along propagation') 

%% Testing 3D plotting. 

figure; 
for i = 1:size(Uall,3)
    % Make field out of each Uall data
    U_all{i} = U; 
    Uall(:,:,i) = Uall(:,:,i)./abs(max(max(Uall(:,:,i)))); % normalization for better visuals
    U_all{i}.data = Uall(:,:,i) ; 
    U_all{i}.position = turb.positions(i);
    pause(0.1)
    imagesc(abs(Uall(:,:,i))./abs(max(max(Uall(:,:,i)))))
end
visualize_optics(U_all{:}) ; title('3D visualization of intensity along propagation'); colormap summer

figure;

%% using iso surface 

%  X = U.X - min(U.X(1,:)); Y = U.Y - min( U.Y(:,1)); Z = turb.positions; E = abs(Uall).^2;  
%  [rows, cols] = size(U.X); num_z = length(Z);
%  % Create 3D matrices
% X3D = repmat(X, [1, 1, num_z]);
% Y3D = repmat(Y, [1, 1, num_z]);
% Z3D = repmat(reshape(Z, 1, 1, []), [rows, cols, 1]);
% 
% figure; scatter3(X3D,Y3D,Z3D,E)
% isovalue = mean(E(:));
% % Create an isosurface
% figure; 
% isosurface(X, Y, Z, E, isovalue);
% figure; 
% % Set properties
% isonormals(X, Y, Z, E, p);
% p.FaceColor = 'interp';
% p.EdgeColor = 'none';
% daspect([1 1 1]);
% view(3);
% camlight;
% lighting gouraud;
% colorbar;
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% title('3D Electric Field Intensity Isosurface');
%% Plot 1D slices top view


slice = mean(Uall,1);
figure; 
subplot(1,2,1)
imagesc(abs(squeeze(slice)).^2)

slice = mean(Uall,2);
subplot(1,2,2)
imagesc(abs(squeeze(slice)).^2)
