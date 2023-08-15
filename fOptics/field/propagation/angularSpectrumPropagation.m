function [Uz, H] = angularSpectrumPropagation(Uo,dx,dy,z,lambda, filter)
%angularSpectrumPropagation - propogate the field to the desired plane
%
% This function uses the Propagation of the Angular Spectrum as described
% in Goodmans' Fourier Optics text book [1]. First the angular spectrum at
% the input plane (z = 0) is calculated using a FFT. A propogation phase is
% applied and then the angular spectrum is iFFT'd to the desired plane. It
% is assumed the distance z is along the optical axis and the input and
% output plane are perpendicular to the optical axis. 
%
% It is important that the propagation phase is properly sampled. As a
% check, one can plot H to check for aliasing. Also, the field can wrap
% around on itself, causing ringing artifacts, if the Uo input is not
% sufficiently large or padded with zeros. 
% 
% Reference: 
%   "Introduction to Fourier Optics" by Joseph W. Goodman 3rd
%   Edition, Roberts and Comapny Publishers, Greenwood Village, CO 2005
%   Section 3.10.2, pg 57
%
% Syntax: [Uz, H] = angularSpectrumPropagation(Uo, dx, dy, z, lambda);
%
% Inputs:
%   Uo - complex field at the input plane
%   dx - x pixel size in meters
%   dy - y pixel size in meters
%   z  - propogation distance in meters
%   lambda - wavelength of light in meters
%   filter - fourier space filter 
% Outputs:
%    Uz - complex field at the output plane
%    H  - propogation phase
% 
% Example: 
%    N = 512;
%    w = 10;
%    Uo = zeros(N, N); 
%    Uo(N/2 - w/2 +1 : N/2 + w/2, N/2 - w/2 +1: N/2 + w/2) = 1;
%    dx = 5.2E-6;
%    dy = dx;
%    z = 2E-2;
%    lambda = 532E-9;
% % 
%    [Uz, H] = angularSpectrumPropagation(Uo,dx,dy,z,lambda);
% 
%    figure;
%    subplot(221); imagesc(Uo); axis image; title('Input Field Amplitude'); 
%    subplot(222); imagesc(angle(H)); axis image; title('Propagation Phase'); 
%    subplot(223); imagesc(abs(Uz)); axis image; title('Output Field Amplitude'); 
%    subplot(224); imagesc(angle(Uz)); axis image; title('Output Field Phase');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: angSpecPropSampling, fresnelDiffractionIntegral


% Author: Dennis F. Gardner Jr.
%         Code 5661
%         4555 Overlook Ave SW
%         Washington, DC 20375
% Email: dennis.gardner@nrl.navy.mil
% Website: https://www.linkedin.com/in/dennisfgardner/
% May 2017; Last revision: 20-March-2018
%
% Future Work: It would be nice to have something that can make sure the
% propogation phase is properly sampled. 
%
% this code is partial based on Abbie's code called ang_spec_propNspk_jrl
% which I think was modified by Ryan. However, it was totally rewritten
% following a cleaner style of code and syntax that follows Goodman's
% Fourier Optics.

arguments
    Uo
    dx
    dy
    z
    lambda
    filter = 1
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Preamble
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

% the size of the imput field (pixels)
[M, N, ~] = size(Uo);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Real Space, at the original plane (before propogation)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

% full width in the x and y direction (units of meters)
xMax = N*dx;
yMax = M*dy;


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Frequency Space
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

% the smallest frequency (delta f) is the inverse of the largest distance
dfx = 1/xMax;
dfy = 1/yMax;

% frequency vectors
fx = (-N/2:1:N/2-1)*dfx; 
fy = (-M/2:1:M/2-1)*dfy; 


% fx = linspace(-N/2, N/2, N)*dfx;
% fy = linspace(-M/2, M/2, M)*dfy;

% frequency grid
[FX, FY] = meshgrid(fx, fy); 

% Goodman likes to use direction cosines (aplha, beta, gamma) as defined in
% Eq. 3-63 and Fig. 3.9
alpha = lambda * FX; 
beta  = lambda * FY;
gamma = sqrt( 1 - alpha.*alpha - beta.*beta);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Circle function to truncate evanescent waves (Eq. 3-67)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

% solution only valid for alpha.^2 + beta.^2 < 1 
alpha_beta_radius = sqrt(alpha.*alpha + beta.*beta);

% make a cricle mask (following section 2.1.6 pg. 12)
circ = zeros(M, N); 
circ( alpha_beta_radius <  1) = 1;
circ( alpha_beta_radius == 1) = 1/2;


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Propogation of the Angular Spectrum
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

% "Angular spectrum of the distrubance U(x, y, 0)" Eq 3-63
A = fftshift(fft2(ifftshift(Uo)))./sqrt(M*N);

% "transfer function of the wave propagation phenomenon" Eq. 3-74
H = exp(1i.*2.*pi.*gamma.*z./lambda).*circ; 
%figure; imagesc(log10(abs(A))); title('transfer function')
% Propagate the angular spectrum to the desired plane (Eq. 3-69)
Uz = fftshift(ifft2(ifftshift(A .* H)))*sqrt(M*N);


end


