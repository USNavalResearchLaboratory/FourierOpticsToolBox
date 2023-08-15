function [numZ, Z, circle] = zernikePolynomials(M, N, nOrder)
%zernikePolynomials - orthogonal Zernike Polynomials in a circle
%
% This funtion uses follows the definitions used in "Standards for
% Reporting the Optical Abberations of Eyes." [1] First the Zernike
% polynomials are generated on a numerical grid with M rows and N cols.
% However, due to the discrete nature of numerical methods, the polynomials
% are not exactly orthogonal as defined [2]. Thus, after the generation of the
% polynomials (and cropping to a circle), they undergo the Gram-Schmidt
% orthogonalization algorithm.
%
% Syntax: [numZ, Z, circle] = zernikePolynomials(M, N, nOrder);
%
% Inputs:
%   M - number of rows
%   N - number of cols
%   nOrder - highest power (order) of Zernike polynomials
%
%
% Outputs:
%    numZ - total number of zernike polynomials, including piston
%    Z - orthogonal Zernike polynomials on a circle
%    circle - pupil aperture of the Zernike polynomials
% 
% Example 1: create the 5th order Zernike polynomials on a 128 x 128 grid
%    [ numZ, Z, ~] = zernikePolynomials(128, 128, 5);
%    plotAll(Z, 1); 
% 
% References (equation reference are found in [1]):
% 
% [1] L. Thibos, R. A. Applegate, J. T. Schwiegerling, and R. Webb,
% "Standards for reporting the optical aberrations of eyes," in Vision
% Science and its Applications, OSA Technical Digest (Optical Society of
% America, 2000), paper SuC1.
%
% [2] D. Malacara, J. M. Carpio-Valadez, and J. J. Sanchez-Mondragon,
% “Wavefront fitting with discrete orthogonal polynomials in a unit radius
% circle,” Opt. Eng. 29, 672–675 (1990).
%
% Other m-files required: createCirclePx
% Subfunctions: radialComponent, gramSchmidt
% MAT-files required: none
%
% See also: createCirclePx, plotAll, zernikePolynomialsRect

% Author: Dennis F. Gardner Jr.
%         Code 5661
%         4555 Overlook Ave SW
%         Washington, DC 20375
% Email: dennis.gardner@nrl.navy.mil
% Website: https://www.linkedin.com/in/dennisfgardner/
% Feb 2018; Last revision: 23-Feb-2018

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% All references to equations are found in Ref [1]
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

% given order n, largest single index is (Eq. 4):
maxIdxJ = (nOrder*(nOrder+2)+nOrder)/2;

% the total number of Zernike polynomials, including piston, is:
numZ = maxIdxJ + 1;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% generate numerical grids, min = -1, max = 1-(1/N) or max = 1-(1/M) 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
x = (-N/2:1:N/2-1)/(N/2);
y = (-M/2:1:M/2-1)/(M/2); 
[X, Y] = meshgrid(x,y); 
[theta, rho] = cart2pol(X,Y); 



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% MAIN LOOP
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

% initialize the Zernike array
Z = zeros(M,N,numZ); 
for n = 0:nOrder
    for m = -n:2:n
        % current single index, zero base index (Eq. 4)
        idxJ = (n*(n+2)+m)/2;
        
        % following Eq. (1)
        % Z = Nnm*Rnm*cos(m theta) for m>=0
        % Z = -Nnm*Rnm*sin(m theta) for m<0
        
        if m>0
            Z(:,:,idxJ+1) = sqrt(2*(n+1))*radialComponent(rho, n, m).*cos(m*theta);
            
        elseif m == 0
            Z(:,:,idxJ+1) = sqrt(n+1)*radialComponent(rho, n, m);
            
        elseif m<0
            Z(:,:,idxJ+1) = -sqrt(2*(n+1))*radialComponent(rho, n, m).*sin(m*theta);
            
        end     
    end
end

% crop the zernike to the circle, radius set to smallest dimension
if M<=N
    circle = createCirclePx(M, N, M/2);
    Z = Z.*circle;
elseif N<M
    circle = createCirclePx(M, N, N/2);
    Z = Z.*circle;
end


% preform Gram-Schmidt orthogonalization, copied from Wiki
% subfunction gramSchimdt needs Z in col vector format
Z = gramSchmidt(reshape(Z,M*N, numZ));
Z = reshape(Z, M, N, numZ); 




% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% SUBFUNCTIONS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

% - Radial Component of Zernike Polynomials, Rnm (Eq. 2)
    function R = radialComponent(rho, n, m)
        % initialize 
        R = 0; 

        % loop through all s from 0 to (n-m)/2 and sum (see Eq, 2) 
        for s = 0:(n-abs(m))/2
            numerator = (-1)^s * factorial((n-s));
            
            % there was a typo in the paper, a missing ')' that I added
            denominator = factorial(s) * factorial(0.5*(n+abs(m))-s) * factorial(0.5*(n-abs(m))-s);
            
            tempR = numerator/denominator * rho.^(n-2*s);
            
            R = R + tempR; 
        end
        
        
    end

% - Gram Schmidt orthogonalization algorithm, copied from Wikipedia
    function U = gramSchmidt(V)
        % code copied from Wiki on 02/15/2018
        % https://en.wikipedia.org/wiki/Gram-Schmidt_process
        k = size(V,1);
        l = size(V,2);
        U = zeros(k,l);
        U(:,1) = V(:,1)/sqrt(V(:,1)'*V(:,1));
        for ii = 2:l
          U(:,ii) = V(:,ii);
          for jj = 1:ii-1
            U(:,ii) = U(:,ii) - ( U(:,ii)'*U(:,jj) )/( U(:,jj)'*U(:,jj) )*U(:,jj);
          end
          U(:,ii) = U(:,ii)/sqrt(U(:,ii)'*U(:,ii));
        end
    end
end

