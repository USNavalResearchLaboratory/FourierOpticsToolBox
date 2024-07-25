classdef pinhole
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        width % 
        diameter
        focus %focus 
        position %position 
        dx %distance between each point
        N %number of points
        X
        Y
        roi
        data %data 
    end

    methods
    
        function obj = pinhole(Phole) %lens struct
            %lens_class Construct an instance of this class
            %   Inputs: 
            % Phole - struct for pinhole
            obj.width = Phole.width;
            obj.diameter = Phole.diameter;
            obj.position = Phole.position; 
            obj.dx = Phole.dx;
            obj.N = Phole.N;
            [obj.X,obj.Y] = meshgrid((-obj.N/2 : obj.N/2-1)*obj.dx); 
            aperture = circ(round(obj.diameter./obj.dx)/2,[obj.N,obj.N]);  
            obj.data = aperture;
            obj.roi = obj.dx*obj.N; 
        end
        function showme(obj,name)
            arguments
                obj
                name = sprintf('Pinhole at z = %f',obj.position);
            end
            fig = figure("Position",[100,100,600,300], "Name",name); 
            limits = [-obj.width obj.width -obj.width obj.width].*1e3/2;
            subplot(1,2,1); imagesc(obj.X(1,:)*1e3, obj.Y(:,1)*1e3, abs(obj.data(:,:,1)).^2); title('Intensity');axis image ; xlabel('[m]')
            xlabel('[mm]'); ylabel('[mm]'); axis(limits)
            subplot(1,2,2); imagesc(obj.X(1,:)*1e3, obj.Y(:,1)*1e3, angle(obj.data(:,:,1))); title('Phase'); axis image
            xlabel('[mm]'); ylabel('[mm]'); axis(limits)
            sgtitle(name);
            %figureFormats(fig);
            
        end
    end
end