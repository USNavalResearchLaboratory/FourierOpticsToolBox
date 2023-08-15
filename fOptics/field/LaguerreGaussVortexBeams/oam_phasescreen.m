function [data] = oam_phasescreen(P, charge, index)

arguments
P
charge = -1;
index = 0; 
end

%Define parameters
N = P.N;
dx = P.dx; % [mm] voxel size
lambda = P.wvl; % [m] wave length
c0 = 3E8; % [m/s] speed of light in air
f = c0/lambda; % [Hz] Frequency 
omega = 2*pi*f; % [rads] radial freq.


w0 = 10*lambda; % width of the beam at z=0;
k = 2*pi/lambda; %Wavenumber 
zR = k*w0^2/2;  % Rayleigh distance;
l =charge;         % topological charge;
n = abs(l)+index;	% radial index; n=|l|,|l|+2,|l|+4 ...
D = 1000;      % is a constant for normalization;

% Define coordinates 
xc= linspace(-(N-1)*dx/2,(N-1)*dx/2,N);
yc= linspace(-(N-1)*dx/2,(N-1)*dx/2,N);
zc = P.Z; %[mm] distance from the transduce device
[X,Y,Z] = meshgrid(xc,yc,zc); %[m]
[TH,R,Z] = cart2pol(X,Y,Z);

% Analytical functions
w = @(z) w0*sqrt(1+(z/zR).^2);
A = @(r,z) (sqrt(2)*r./w(z)).^(abs(l)).*LaguerreL((n-abs(l))/2,abs(l),2*r.^2./w(z).^2);
G = @(r,z) D./sqrt(1+z.^2/zR^2).*exp(r.^2./w(z).^2).*exp(-0.5i*k*(z.*r.^2)./(z.^2+zR^2));
PHI = @(th) exp(1i*l*th);
PSI = @(z) exp(-1i*(n+1)*atan(z/zR));
P = @(th,r,z,t) G(r,z).*A(r,z).*PHI(th).*exp(1i*(k*z-omega*t)).*PSI(z);

% Compute profile for a seleted time 't':
data=P(TH,R,Z,0);


end