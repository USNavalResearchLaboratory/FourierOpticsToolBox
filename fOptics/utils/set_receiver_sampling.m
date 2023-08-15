function [dx2] = set_receiver_sampling(lambda, Z, D2, dx1, D1, N)
    %SET_RECEIVER_SAMPLING - Computes the required dx2 and Nconstraint1 variables for angular spectrum propagation
    %   [dx2, Nconstraint1] = SET_RECEIVER_SAMPLING(lambda, Z, D2, dx1, D1, N) returns dx2 and Nconstraint1 for a given 
    %   lambda, Z, D2, dx1, D1, and N. It checks if the user input N is greater than or equal to Nconstraint1 and rounds up
    %   Nconstraint1 to the nearest integer.
    
    %   lambda: Wavelength of the propagating wave
    %   Z: Distance to the object plane from the camera
    %   D2: Diameter of the receiving aperture
    %   dx1: Sampling interval in the object plane
    %   D1: Diameter of the object plane
    %   N: Number of samples in the image plane
    
    %   dx2: Sampling interval in the image plane
    %   Nconstraint1: The minimum number of samples required in the image plane
    
    % Example usage:
    %   [dx2, Nconstraint1] = set_receiver_sampling(650e-9, 1, 0.025, 3.9e-6, 0.05, 2048);
    
    %TODO: minimize dx1 such that Nconstraint is satisfied
    multiplier = 20; %The original equation gave the maximum dx2 could be. But it results in very pixelated results. The multiplier allows dx2 could be smaller for less pixilation. 
    
    dx2 = (lambda * Z - D2 * dx1) / (multiplier*D1);
     
    if dx2 <0
        dx2 = dx1;
    end

end
    