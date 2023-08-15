function [ lensField ] = lensPhaseCustom(P,focus)
%lensPhase -  phase delay imparted by a perfect, spherical (in the paraxial
%sense), thin lens Eq. 4.7 in Jason Schmidt's book Numerical Simulation of
%Optical Wave Propagation Or Eq. 5-10 in Goodman Fourier Optics. 
%
% The focus is at the target. The ouput is complex. P is the parameter
% structure; 


k = (2*pi)/(P.wvl);

[X, Y] = meshgrid((-P.NmultiStep/2:1:P.NmultiStep/2-1).*P.R.px);


Rsqrd = X.*X+Y.*Y;

lensField = exp(-1i*k/2/focus*Rsqrd); %Harshil - think theres a factor of 2 missing. Need to check 

end

