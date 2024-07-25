function [r0scrn] = Cn2_to_r0(L,wv,z,Cn2,N_Sc,places)
%pg 165 of Numerical Simulation of Optical Wave Propagation with Examples in MATLAB

k = 2*pi / wv; % optical wavenumber

% use sinc to model pt source (THESE ARE NEVER USED??)
%DROI = 4 * L;          % diam of obs-plane region of interest [m]
%D1 = wv*Dz / DROI;     % width of central lobe [m]
%R = Dz;                % wavefront radius of curvature [m]

% MTR: 3/24/23
% In the expression for r0pw, why is it max(z) instread of dZ
% In the expresison for r0sw, why is it 3/8 instad of (zi/dZ)^(5/3)

% SW and PW coherence diameters [m]
r0sw = (0.423 * k^2 * Cn2 * 3/8 * max(z))^(-3/5);
r0pw = (0.423 * k^2 * Cn2 * max(z))^(-3/5);
p = linspace(0, max(z), 1e3);
% log-amplitude variance (aka rytov number)
rytov = 0.563 * k^(7/6) * sum(Cn2 * (1-p/max(z)).^(5/6) ...
    .* p.^(5/6) * (p(2)-p(1)));

% screen properties
nscr = N_Sc; % number of screens
A = zeros(2, nscr); % matrix
alpha = places./max(z);
A(1,:) = alpha.^(5/3);
A(2,:) = (1 - alpha).^(5/6) .* alpha.^(5/6);
b = [r0sw.^(-5/3); rytov/1.33*(k/max(z))^(5/6)]; %Set up Ax = b

X = lsqlin(A,b,[],[],[],[], .16.*ones(N_Sc,1),[]); %linear optimization
r0scrn = X.^(-3/5)
r0scrn(isinf(r0scrn)) = 1e6; %if your r_0 is 1e6 something is probably wrong
% MTR: 12APR23. Added this catch for Cn2 values of zero, intended to
% propagate a perfect beam.  Otherwise, r0scr calcualtions produce finite
% values of r0scr which adversely affects beam propagation.  In Cn2 is 0,
% the phase screen should ideally have no scintillation.
if Cn2 == 0
    r0scrn = r0scrn*inf;
end


end