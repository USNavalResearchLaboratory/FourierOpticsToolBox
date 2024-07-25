%Gaussian propagation test to check divergence 


%close all; clear all; clc
%close all
% Inputs
[ P ] = parameterSetup_gaussTest(1e-13); 

% Generate turbulence/ phase screen
turb = phase_screens(P.air); %Creates initial flat screens
air = turb; 
%turb = turb.kolmogorov(P.Turb); %Generates turbulence. try turb.tilt, or turb.OAM
angle = 0; 
turb = turb.tilt(angle,0); 
turb.showme


% Setup target
target = target_class(P.Target); % Define initial target
target = target.genTarget_simplegauss(P.Target.D); 
%target = target.genTarget_ptsrc(); 
target.showme
%
U = Efield(target);   
lens1 = thinLens(P.lens1);
iter = 1; 
zarray = 40; 
for zloop = zarray
    
    Uf = U.propagate_ang_spec(zloop, 'fwd', U.width,turb, 1);
    Uf.showme
    Uf = Uf.resample(lens1); %Resample field to same dimensions as lens1
    slice(iter,:) = abs(Uf.data(Uf.N/2,:));
    slicex(iter,:)  = Uf.X(Uf.N/2,:);
    figure(100); hold on; plot(slicex(iter,:),slice(iter,:), 'DisplayName', num2str(zloop) )
    wz_numerical(iter) = findFWHM(slicex(iter,:), slice(iter,:))*.8493218; %Convert FWHM to 1/e^2 beam waist
    iter = iter+1;
    
end

%tilt phase test
displacement = tan(angle*pi/180)*zarray*1e3
[~,displacement_actual] = max(Uf.data(:)); Uf.X(displacement_actual)*1e3


%% Analytical beam waist

w0 = wz_numerical(1);
zr = pi*w0^2/target.wvl;
wz_analytical = w0*sqrt(1+(zarray/zr).^2);

figure; plot(zarray, wz_analytical, zarray, wz_numerical)
xlabel('Distance [m]')
ylabel('Beam Waist [m]')
legend('Analytical','Numerical')
