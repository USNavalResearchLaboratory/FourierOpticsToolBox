function [E] = lensPhaseBasic(P)

% P - parameters file 
f = P.focus; % focal length
dx = P.dx;% distance between pixels [m] 
N = P.N; 
d = P.width; 
wvl = P.wvl; 


%Alternate calculation 
k = (2*pi)/(wvl);
[X, Y] = meshgrid((-N/2:1:N/2-1).*dx);
Rsqrd = X.*X+Y.*Y;
lensField = exp(-1i*k*Rsqrd/(2*f));
E = lensField; 

end