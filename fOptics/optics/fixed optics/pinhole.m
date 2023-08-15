classdef pinhole
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        diameter %diameter 
        focus %focus 
        position %position 
        dx %distance between each point
        N %number of points
        X
        Y
        data %data 
    end

    methods
    
        function obj = pinhole(Phole) %lens struct
            %lens_class Construct an instance of this class
            %   Inputs: 
            % Phole - struct for pinhole
            obj.diameter = Phole.diameter;
            obj.position = Phole.position; 
            obj.dx = Phole.dx;
            obj.N = Phole.N;
            [obj.X,obj.Y] = meshgrid((-obj.N/2 : obj.N/2-1)*obj.N); 
            aperture = circ(round(obj.diameter./obj.dx)/2,[obj.N,obj.N]); 
            obj.data = aperture;
             
        end
        function showme(obj)

            fig = figure("Position",[100,100,600,300]); 
            imagesc(obj.X, obj.Y, abs(obj.data(:,:,1)).^2); title('Intensity');
            figureFormats(fig)
           
            
        end
    end
end