
% 
clear; clc; close all;
rng(0)
%% Calibration -- ONLY NEED TO RUN ONCE UNLESS OPTICS ARE CHANGED
nOrder = 5; 
[dZ,zernikes] = radial_shearing_calibration_sim(nOrder);
zernikes = (zernikes); dZ = angle(dZ);
%dZ = gramSchmidt(dZ); 
%% Get experimental/Modeled data
close all;
rng(1)
main_script_FINCH; % Run main script to generate REF and signal
% Massage the data
detectorRadial = hologram_detector('shearing');
detectorRadial.Uref = U_pupilPlane_dRREF; detectorRadial.Usignal = U_pupilPlane_dR; %assign ref and signal
recovered_phaseRadial = radialShearPhaseRecovery(detectorRadial);
recovered_phaseRadial.showme('Radial phase difference (Input)')
dPhase = [recovered_phaseRadial.data(:)];
%use phase only
dPhase = angle(dPhase);
%% Least squares phase recovery
%Attempt using lsqr function
TOLERENCE = 1e-3;
MAX_ITER = 100;
zern_coeffs = lsqr((dZ),dPhase, TOLERENCE, MAX_ITER);
zern_coeffs2 = lsqr(gramSchmidt(dZ),dPhase, TOLERENCE, MAX_ITER); %using gramSchmidt function to orthogonalize dZ basis set. but it messes up my scale-> phase recovery
% I am finding that with the gramScmidt orthogonalization, lsqr converges
% faster, but zern_coeffs2 is little off
zern_coeffs2 = zern_coeffs2./max(zern_coeffs2) .* max(zern_coeffs); 

snr_zerns_1 = max(zern_coeffs)./(sum(zern_coeffs)-max(zern_coeffs)); % signal to noise ratio of zerns. Works for testing with one zern in turb only
snr_zerns_2 = max(zern_coeffs2)./(sum(zern_coeffs2)-max(zern_coeffs2));
%zern_coeffs = turb.zernike_distribution'; 
% zern_coeffs = ip_sol(dZ,dPhase);
% Attempt using pinv function
% Zinv = pinv(dZ'*dZ)*dZ';
% zern_coeffs = Zinv*dPhase;
gain = 1; 
% Calculate phase from zernikes

zernikes_grid = reshape(zernikes,size(zernikes,1)^2,[]) ;
zern_recovered = zernikes_grid.*zern_coeffs'; 
zern_recovered = reshape(zern_recovered,size(zernikes,1),size(zernikes,2),[]) ;
phase = sum(zern_recovered,3);
% Normalize
% phase = phase./max(phase(:));
%phase = angle(phase); 
%phase = phase - mean(phase(:)); % Subtract mean
phase = phase.*gain; 
%phase = wrapToPi(phase);
phase_recovered = recovered_phaseRadial; % Initialization
phase_recovered.data = exp(1i.*(phase)); %phase_recovered.phase = phase;
%%\ Do some magnification because I only measure information in overlapping area. Interpolated outside
magnification = 1; % demagnify the mag from optical setup
phase_recovered = phase_recovered.magnify(1/magnification);
phase_recovered.data = abs(lens1.data).*phase_recovered.data;
%phase_recovered = phase_recovered.propagateTo(lens1,'rev');
phase_recovered_tmp = U_pupilPlane; % my X, Y scale is messed up somewhere. This is a tmp fix
phase_recovered_tmp.data = phase_recovered.data; 
phase_recovered = phase_recovered_tmp; 
phase_recovered.showme('DZT method phase recovery (ESTIMATED)')

% phase difference
phase_difference = angle(phase_recovered.data) - angle(U_truth.data); 
figure("Position",[700,200,700,400]); imagesc(phase_recovered.X(1,:).*1e3, phase_recovered.Y(:,1).*1e3,phase_difference)
title('Phase Error'); axis image
xlabel('[mm]'); ylabel('[mm]'); colorbar
fontsize(gcf,18, 'pixels'); fontname(gcf,"Times New Roman")

% % True phase
% zern_recovered = reshape(zernikes,size(zernikes,1)^2,[]) ;
% zern_recovered = zern_recovered.*turb.zernike_distribution ; 
% zern_recovered = reshape(zern_recovered,size(zernikes,1),size(zernikes,2),[]) ;
% phase = sum(zern_recovered,3);
% % Normalize
% %phase = phase./max(phase(:));
% %phase = wrapToPi(phase);
% phase_recovered = recovered_phaseRadial; 
% phase_recovered.phase = phase; phase_recovered.data = phase; 
% magnification = 0.2; 
% phase_recovered = phase_recovered.magnify(1/magnification);
% phase_recovered.showme('Phase from zernikes (True)')

%% Plot recovered Zernikes 
Zernike_error = [];
figure("Name", "Recovered Zernike Profile");
bar(zern_coeffs)
axis([0 length(zern_coeffs2) -1 1])
sgtitle('Recovered Zernike Profile')
xlabel('Zernike mode'); ylabel('Amplitude')
set(gca,"FontSize",20)
% Make a list of zernike error to collect errors for publication plot
%Zernike_error = [Zernike_error;sum(zern_coeffs' - turb.zernike_distribution(4:end))]; 

%% save data

%save("data_zernike_10")

%% Complete AO loop
U_return = U_pupilPlane;
%Intensity_return = abs(U_pupilPlaneREF.data);
U_return = U_return.resample2(lens1); U_return.position = lens1.position; 
U_return.data = abs(lens1.data).*(phase_recovered.data).*exp(-1i.*angle(lens1.data));
U_return.showme('Backprop Image')

Image = target_class(P.Target); % Define initial target
Image = Image.genTarget_focus_spot(U_return, turb, flags.gpu); 
Image.width = 0.05; 

Image_ref = target_class(P.Target); % Define initial target
Image_ref = Image_ref.genTarget_focus_spot(lens1.conjugate, turb, flags.gpu); 
Image_ref.width = 0.05; 


% Truth 
U_return = U_pupilPlane; 
U_return.data = abs(lens1.data).*exp(1i.*angle(U_truth.data));
U_return = U_return.resample2(lens1); U_return.position = lens1.position; 
U_return = U_return.apply_optic(lens1.conjugate); 
%U_return.showme('Truth Returning Image')

Image_truth = target_class(P.Target); % Define initial target
Image_truth = Image_truth.genTarget_focus_spot(U_return, turb, flags.gpu); 
Image_truth.width = 0.05; 

%% Normalize to truth
normalizer = max(abs(Image_truth.data(:)));
Image_truth.data = Image_truth.data./normalizer; 
Image.data = Image.data./normalizer; 
Image_ref.data = Image_ref.data./normalizer; 

Image.showintensity('Corrected Image')
Image_ref.showintensity('Uncorrected Image')
Image_truth.showintensity('Corrected Truth Image')

%% Calculate power in bucket
bucket = 2*1e-3; %bucket in millimeters
bucket = round(bucket/Image_truth.dx); 
bucket = circ(bucket,[Image.N,Image.N]);
Image_truth.data = Image_truth.data.*bucket; 
Image_ref.data = Image_ref.data.*bucket; 
Image.data = Image.data.*bucket; 

power_bucket_truth = Image_truth.calc_power.avg_power;
power_bucket_ref = Image_ref.calc_power.avg_power;
power_bucket_actual = Image.calc_power.avg_power;

vertical_beam_quality_ref = sqrt(power_bucket_truth/power_bucket_ref);
vertical_beam_quality_actual = sqrt(power_bucket_truth/power_bucket_actual);
vertical_beam_quality_truth = sqrt(power_bucket_truth/power_bucket_truth);

sprintf('Vertical beam quality Uncorrected %f',vertical_beam_quality_ref)
sprintf('Vertical beam quality actual (1=best) %f',vertical_beam_quality_actual)
sprintf('Vertical beam quality perfect %f',vertical_beam_quality_truth)
%% Figures for publication


