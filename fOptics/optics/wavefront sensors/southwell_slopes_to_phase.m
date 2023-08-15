function [phase_front] = southwell_slopes_to_phase(slopeMat,N, gpuFlag)

%% Southwell method copied from Theo's code singleforpaper_lsqtester.m 
arguments
    slopeMat % TODO Put in default values for testing
    N 
    gpuFlag = 0
end

%N = sqrt(length(slopeMat)/2);
slopeMat = slopeMat(:); 
if gpuFlag
    southwellD = gpuArray(zeros(2*N*(N-1),2*N^2));
else
    southwellD = (zeros(2*N*(N-1),2*N^2));
end


for p = 1:N*(N-1)
    southwellD(p,p) = 0.5;
    southwellD(p,p+N) = 0.5;
%      figure(); imagesc(southwellD); pause(0.01)
end

for a = 1:(N-1)
    for b = 1:N
        t = a + (b-1)*N;
        southwellD(N*(N-1)+t+1-b,t+N^2) = 0.5;
        southwellD(N*(N-1)+t+1-b,t+N^2+1) = 0.5;
      %   figure(12); imagesc(southwellD); pause(0.01)
    end
    
end

southwellA = zeros(2*N*(N-1),N^2);

for p = 1:N*(N-1)
    southwellA(p,p) = -1;
    southwellA(p,p+N) = 1;
%     figure(2); imagesc(southwellA); pause(0.01)
end

for a = 1:(N-1)
    for b = 1:N
        t = a + (b-1)*N;
        southwellA(N*(N-1)+t+1-b,t) = -1;
        southwellA(N*(N-1)+t+1-b,t+1) = 1;
%         figure(); imagesc(southwellA); pause(0.01)
    end
    
end
Ae = [southwellA; ones(1,N^2)];
zonalSolve = pinv(Ae' * Ae) * southwellA' * southwellD * 2.1989;


phase_front = wrapToPi(reshape(zonalSolve * slopeMat(:),[N N]));
%phase_front = zonalreconWF(:,:,1);
% figure(1);subplot(1,3,2); imagesc((zonalreconWF(:,:,1))); colorbar; subtitle('reconstructed SH wavefront')
% subplot(1,3,1); imagesc((unwrappedAngle.*source)); subtitle('true wavefront'); colorbar
% addpath('K:\workspaces\Dave\Projects\Wavefront_reconstruction\WorkingResults\theo_weights_test\Cn2_4p90E-14_B\power_3p00E-11')
% cnn_wf = readmatrix('Predicted_wavefront_2wfs_figure_cnn_slopes_no_train_bl15.csv')* slopeStd;
% subplot(1,3,3); imagesc((cnn_wf)); subtitle('CNN_wavefront'); colorbar


% figure(2);subplot(1,3,1); imagesc(abs(jfft2(exp(1i.*zonalreconWF(:,:,1)))).^2); colorbar; subtitle('Sh PSF'); xlim([60 130]); ylim([60 130]); 
% subplot(1,3,2); imagesc(abs(jfft2(exp(1i.*unwrappedAngle))).^2); subtitle('true PSF'); xlim([60 130]); ylim([60 130]); 
% subplot(1,3,3); imagesc(abs(jfft2(exp(1i.*cnn_wf))).^2); subtitle('CNN PSF')

end