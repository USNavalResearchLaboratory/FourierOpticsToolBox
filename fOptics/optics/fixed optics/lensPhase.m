function [ lensField ] = lensPhase(P)
%lensPhase -  phase delay imparted by a perfect, spherical (in the paraxial
%sense), thin lens Eq. 4.7 in Jason Schmidt's book Numerical Simulation of
%Optical Wave Propagation Or Eq. 5-10 in Goodman Fourier Optics. 
%
% The focus is at the target. The ouput is complex. P is the parameter
% structure; 
% If nonlinear propagation is enabled, will apply a shallower or
% compensatory focus such that, via the thin lens equation, the
% self-focusing length and lens focal length combine to focus at the target

if P.S.NLprop
    tgtFocLength = P.Z;
    
    rwaist = P.S.D*P.S.ApodSigmaFrac;
    rayleighLeng = pi* rwaist^2/ P.wvl; %rayleigh length (m)
    slfFocLength = rayleighLeng * 0.5 / sqrt(P.S.NLfrac - 1);
    
    focLength = 1 / (1/tgtFocLength - 1/slfFocLength);
else
    focLength = P.Z;
end



k = (2*pi)/(P.wvl);

[X, Y] = meshgrid((-P.R.N/2:1:P.R.N/2-1).*P.R.px);

Rsqrd = X.*X+Y.*Y;

lensField = exp(-1i*k/2/focLength*Rsqrd);

end

