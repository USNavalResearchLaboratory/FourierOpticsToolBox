function [phi] = turbulence_phi_fftGPU(D,r0,K)

% GPU version - Harshil May 18 2023
% generates a phase screen with Komolgorov statistics via FFT techniques
%
% phi = turbulence_phi_fft(D,r0,K)
%
% phi = generated phase screen
% D = dimensions of phi (default = 128)
% r0 = Fried's parameter in units of pixels (default = 128)
% K = number of subharmonic sets to include in fix for low spatial
%      frequencies

arguments
    D= 2048*2 % default dimension of output array
    r0 = D % default Fried parameter
    K = 5

end

N = 2*D; % dimensions of bigger array
dn = 1/N; % Fourier domain sample spacing
[n2,n1] = meshgrid(dn*fftshift(-N/2:N/2-1)); % Fourier domain sample coordinates
r = sqrt(n1.^2+n2.^2); % radial coordinate
r(1) = 1; % to fix singularity at origin
phi = gpuArray(sqrt(dn^2*0.023*r0^(-5/3)*r.^(-11/3)).*(randn(N,'like', 1i))); % square root of Komolgorov power spectrum
phi = fftshift(real(fftn(phi))); % phase screen

phi = phi(N/2-round(D/2)+[1:D],N/2-round(D/2)+[1:D]); % crop output out of bigger array

% fix for inadequately sampled low spatial frequencies of power spectrum
ns1 = gpuArray([]); % subharmonic sample coordinates
ns2 = gpuArray([]);
dns = gpuArray([]); % subharmonic sample spacings
for k = 1:K
    dn = dn/3;
    [tmp2,tmp1] = meshgrid(dn*[-1,0,1]);
    ns1 = [ns1;tmp1(:)];
    ns2 = [ns2;tmp2(:)];
    dns = [dns;repmat(dn,[9,1])];
end
rs = gpuArray(sqrt(ns1.^2+ns2.^2)); % subharmonic radial coordinates
[d2,d1] = meshgrid(1:D); % phase screen domain sample coordinates
d2 = gpuArray(d2); d1 = gpuArray(d1); 

for sub = 1:9*K 
    if rs(sub)>0 % add each subharmonic contribution
        phi = phi+sqrt(dns(sub)^2*0.023*r0^(-5/3)*rs(sub)^(-11/3))*real((randn(1,'like',1i))*exp(-1i*2*pi*(ns1(sub)*d1+ns2(sub)*d2)));
    end
end


%%
phi = phi-mean(phi(:));
end