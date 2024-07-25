function H = fresnel_propagator(z, Lx, Mx, lambda)
    arguments
        z % propagataion distance
        Lx = 250e-3
        Mx = 1024
        lambda = 490e-6
    end
    Ly = Lx; 
    My = Mx; 
    % Define Fresnel Propagtor
    dx = Lx/Mx;
    dy = Ly/My;
    % Define frequency axes
    fMax_x = 1/(2*dx);
    df_x = 1/Lx;
    fx = -fMax_x:df_x:fMax_x-df_x;
    fMax_y = 1/(2*dy);
    df_y = 1/Ly;
    fy = -fMax_y:df_y:fMax_y-df_y;
    [FX,FY] = meshgrid(fx,fy);
    % Circle function to truncate evanescent waves
%     alpha = lambda * FX; 
%     beta  = lambda * FY;
%     gamma = sqrt( 1 - alpha.*alpha - beta.*beta);
%     % solution only valid for alpha.^2 + beta.^2 < 1 
%     alpha_beta_radius = sqrt(alpha.*alpha + beta.*beta);
%     circ = zeros(Mx, Mx); 
%     circ( alpha_beta_radius <  1) = 1;
%     circ( alpha_beta_radius == 1) = 1/2;


    H = exp(2i*pi*z/lambda) .* exp(-1i*pi*lambda*z*(FX.^2 + FY.^2));
end