function [Uz, newdx, newdy] = fresnelDiffractionIntegral(Uo,dx,dy,z,lambda)
%fresnelDiffractionIntegral - propogate the field to the desired plane
%
% This function uses the Fresnel difffraction integral as described
% in Goodmans' Fourier Optics text book [1]. When z is negative the FFT2
% turns into an iFFT2 and I added some abs() functions around the z to keep
% the newdx and newdy positive. The abs() in in the final result is to keep
% the phase the same when propagaing foward and back by z. 
%                              x xd
%                              deeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
%                                                                                                                                                                           c
% Syntax: [Uz] = fresnelDiffractionIntegral(Uo, dx, dy, z, lambda);
%
% Inputs:
%   Uo - complex field at the input plane
%   dx - x pixel size in meters at input plane
%   dy - y pixel size in meters at input plane
%   z  - propogation distance in meters
%   lambda - wavelength of light in meters
%
% Outputs:
%    Uz - complex field at the output plane
%    newdx - sampling period in the observation plane
%    newdy - sampling period in the observation plane
% 
% Example: 
%     N = 128;
%     dx = 40E-6;
%     dy = dx;
%     z = 0.5;
%     lambda = 1E-6;
% 
%     w = round((2E-3)/dx/2);
% 
%     Uo = zeros(N, N); 
%     Uo(N/2 - w/2 +1 : N/2 + w/2, N/2 - w/2 +1: N/2 + w/2) = 5;
% 
% 
%     [Uz, newdx, newdy] = fresnelDiffractionIntegral(Uo,dx,dy,z,lambda);
%     figure; imagesc(abs(Uz).^2);
%
% Reference: 
% [1] "Introduction to Fourier Optics" by Joseph W. Goodman 3rd
%     Edition, Roberts and Comapny Publishers, Greenwood Village, CO 2005
%     Section 4.2, pg 67
%
% [2] "Numerical Simulation of Optical Wave Propagation" by Jason D.
%     Schmidt, SPIE Press, Bellingham, WA 2010. pg. 120
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: fresnelIntegralSampling, angularSpectrumPropagation


% Author: Dennis F. Gardner Jr.
%         Code 5661
%         4555 Overlook Ave SW
%         Washington, DC 20375
% Email: dennis.gardner@nrl.navy.mil
% Website: https://www.linkedin.com/in/dennisfgardner/
% March 2018; Last revision: 26-March-2018


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Preamble
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

% the size of the imput field (pixels)
[M, N, ~] = size(Uo);

% wave number
k = 2*pi/lambda;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Grid at the diffracting aperture (xi, eta)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

xi =  (-N/2:1:N/2-1)*dx; 
eta = (-M/2:1:M/2-1)*dy; 

[XI, ETA] = meshgrid(xi,eta); 

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Phase factor, inside the integral
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

phaseIN = exp( 1i *k/(2*z).*  (XI.*XI + ETA.*ETA) ); 

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Evaluate the integral (FFT)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

if z == 0
    Uz = Uo; 
    
elseif z > 0
    Uz = fftshift(fft2(ifftshift(Uo.*phaseIN)))./sqrt(M*N).*(z*lambda);
    
elseif z < 0
    Uz = ifftshift(ifft2(fftshift(Uo.*phaseIN))).*sqrt(M*N).*(z*lambda);
        
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Sampling in the observation plane (x,y)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% this is from Eq. 7.21 in Ref. 2, I added the abs() around Z for negative
% propagation distances

newdx = lambda*abs(z)/N/dx; 
newdy = lambda*abs(z)/M/dy; 

x = (-N/2:1:N/2-1)*newdx; 
y = (-M/2:1:M/2-1)*newdy; 

[X, Y] = meshgrid(x,y); 


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Phase factors after the integral
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
phaseOUT1 = exp(1i*k*z);
 
phaseOUT2 = exp(1i*k/(2*z).*  (X.*X + Y.*Y) );

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Final result
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

Uz = (1/(1i*lambda*abs(z))).*phaseOUT1.*phaseOUT2.*Uz; 

end


