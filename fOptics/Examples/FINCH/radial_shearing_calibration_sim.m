% Main script to calibrate lateral shearing interferometry
function [dZ,zernikes_normalized] = radial_shearing_calibration_sim(nOrder,zernikes)
arguments
    nOrder = 3
    zernikes = []
end
%% Generate Zernike phase screens
Cn2 = 0;
distance = 500; %MATCH FINCH_opticalSetup!!
[ P ] = parameterSetup_targetAndTurb(Cn2,distance);
flags.gpu = 0;
flags.speckle = 0;
air = phase_screens(P.air); %Creates initial flat screen
if isempty(zernikes)
    % Generate zernikes
    [~, zernikes, ~] = zernikePolynomials(air.N, air.N, nOrder);
    %zernikes(:,:,1:3) = ones([P.air.N, P.air.N, 3]); %Ignore first 3 modes

end
scale = zeros([1 1 size(zernikes,3)]); scale(:) = 1;
normalization = max(max(zernikes,[],1),[],2);
zernikes_normalized = zernikes./normalization;
zernikes =  exp(1i .* zernikes_normalized.*scale);
turb = air;

%% Pass through my shearing interferometry

% Going to be a for loop that sets one zernike in phase screen at a time.
% First one should be Ref
dZarrayX = zeros([turb.N^2, size(zernikes,3)]);
dZarrayY = zeros([turb.N^2, size(zernikes,3)]);

if flags.gpu
    dZarrayX = gpuArray(dZarrayX);
    dZarrayY = gpuArray(dZarrayY);
end
%scr = round(P.Turb.screens / 2);
scr = P.Turb.screens;
figure(1)
tiledlayout(size(zernikes,3),2, 'TileSpacing','none')
sgtitle('Radial Shear Individual Zernike measurements (d\phi)')

for iter = 1:size(zernikes,3)
    turb.data(:,:,scr) = zernikes(:,:,iter);
    %turb.showme
    [U_cam_dR, U_cam_dY] = radialShearForOneZern(P, turb, flags);

    % user first iter at reference. Do phase recovery
    if iter == 1
        U_cam_dYREF = U_cam_dY; U_cam_dXREF = U_cam_dR;
    end
    detectorY = hologram_detector('shearing');
    detectorY.Uref = U_cam_dYREF; detectorY.Usignal = U_cam_dY; %assign ref and signal
    recovered_phaseY = radialShearPhaseRecovery(detectorY);

    detectorX = hologram_detector('shearing');
    detectorX.Uref = U_cam_dXREF; detectorX.Usignal = U_cam_dR; %assign ref and signal
    recovered_phaseX = radialShearPhaseRecovery(detectorX);
    % Mask extra stuff on edges
    %     mask = circ(recovered_phaseX.width/recovered_phaseX.dx/2,[recovered_phaseX.N,recovered_phaseX.N]);
    %     recovered_phaseX.data = recovered_phaseX.data.*mask;
    %     recovered_phaseY.data = recovered_phaseY.data.*mask;

    dZarrayX(:,iter) = recovered_phaseX.data(:);
    dZarrayY(:,iter) = recovered_phaseY.data(:);
    figure(1)
    nexttile
    imagesc(angle(turb.data(:,:,scr) ))

    axis off; axis image
    nexttile
    imagesc(recovered_phaseX.X(1,:)*1e3, recovered_phaseX.Y(:,1)*1e3, angle(recovered_phaseX.data(:,:,1)));
    limits = [-recovered_phaseX.width recovered_phaseX.width -recovered_phaseX.width recovered_phaseX.width].*1e3/2;
    axis off; axis image; axis(limits);

end

% Create Column vector for each zern. Put into 2D matrix of size (N^2 x
% Nzerns)
% Construct dZ matrix
dZ = dZarrayX;
end


%% Function that does Radial Shearing for one realization.
% This is for a simulated model. Use another function for experiment
function [U_cam_dR, U_cam_dY] = radialShearForOneZern(P, turb, flags)

target = target_class(P.Target); % Define initial target class
lens1 = thinLens(P.lens1);
%target = target.genTarget_simplegauss(P.Target.width);
air = phase_screens(P.air); %Creates initial flat screen
target = target.genTarget_focus_spot(lens1.conjugate, air, flags.gpu);

if flags.speckle == 1
    smoothing_factor = 4;
    target = target.addSpeckle(smoothing_factor);
end
U_target = Efield(target,flags.gpu);
U_pupilPlane = U_target.propagateTo(lens1,'fwd',turb);
% U_pupilPlane.showme('Field at pupil')
[U_cam_dR] = FINCH_opticalSetup(U_pupilPlane);
U_cam_dY = U_cam_dR;
end
