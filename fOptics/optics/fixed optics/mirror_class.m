classdef mirror_class
    %Mirrors Summary of this class goes here
    %   You can apply a tilt to field object
    % You can also reflect light (flipping object)

    properties
        diameter
        position
        dx
        N
        X
        Y
        angle_rad
        data
    end

    methods
        function obj = mirror_class(P,Pmirror) %lens struct
            %lens_class Construct an instance of this class
            %   Inputs: 
            % P - Parameters struct for general simulation parameters
            % Plens - Parameters struct for specific lens
            obj.diameter = Pmirror.diameter;
            obj.position = Pmirror.position; 
            obj.dx = Pmirror.dx;
            obj.N = Pmirror.N;
            [obj.X,obj.Y] = meshgrid((-obj.N/2 : obj.N/2-1)*obj.dx); 

             
        end
        function obj = tilt(obj,Pmirror)
            obj.angle_rad = Pmirror.angle_deg*pi/180; %deg in rad
            n = abs(obj.diameter*tan(obj.angle_rad)/(2*pi));
            theta = linspace(0,n*(2*pi),obj.N);
            tiltphase = repmat(theta, [(obj.N),1]); %the 2D tilt phase gradient due to grating
               
            obj.data = exp(1i.*tiltphase); 
        end
        function obj = tilt2(obj,Pmirror)
            offset = 200; 
            delta = zeros(obj.N); 
            delta(obj.N/2+offset,obj.N/2+offset) = 1; 
     
            tiltphase = angle(jifft2(delta)); %the 2D tilt phase gradient due to grating
               
            obj.data = exp(1i.*tiltphase); 
        end
        function obj = reflect(obj)
            %TODO
        end
        
        function showme(obj)

            figure; 
            subplot(1,2,1); imagesc(abs(obj.data(:,:,1))); title('abs')
            subplot(1,2,2); imagesc(angle(obj.data(:,:,1))); title('angle')
            sgtitle('Mirror')
        end
    end
end