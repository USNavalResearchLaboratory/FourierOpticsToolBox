function phz = ft_phase_screen(r0, N, delta, L0, l0)
%mostly code from Numerical Simulation of Optical Wave Propagation with Examples in MATLAB

    % setup the PSD
    del_f = 1/(N*delta);   % frequency grid spacing [1/m]
    fx = (-N/2 : N/2-1) * del_f;
    % frequency grid [1/m]
    [fx fy] = meshgrid(fx);
    [th f] = cart2pol(fx, fy);  % polar grid
    fm = 5.92/l0/(2*pi);        % inner scale frequency [1/m]
    f0 = 1/L0;                  % outer scale frequency [1/m]
    % modified von Karman atmospheric phase PSD
    PSD_phi = 0.023*r0^(-5/3) * exp(-(f/fm).^2) ...
        ./ (f.^2 + f0^2).^(11/6);
    PSD_phi(N/2+1,N/2+1) = 0;
    % random draws of Fourier coefficients
    cn = (randn(N) + i*randn(N)) .* sqrt(PSD_phi)*del_f;
% why dependent on N? 
    % synthesize the phase screen
    phz = real(ifft2(cn));