classdef thinLens
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        wvl = 650e-9; %wavelength
        N % Size of data array
        dx %assuming dx=dy. use symmetric optics
        X % X position meshgrid
        Y % Y position meshgrid
        width = 0; % width of object (not total width of array dx*N)
        roi % Region of interest. Equals dx*N
        position = 0; % Position of field along z axis
        focus
        data = [];
    end

    methods
        function obj = thinLens(Plens) %lens struct
            %lens_class Construct an instance of this class
            %   Inputs: 
            % Plens - Parameters struct for specific lens
            % U - field at plane with lens. Determines dx
            obj.width = Plens.diameter;
            obj.focus = Plens.focus;
            obj.position = Plens.position; 
            obj.dx = Plens.dx;
            obj.N = Plens.N;
            obj.wvl = Plens.wvl; 
            [obj.X,obj.Y] = meshgrid((-obj.N/2 : obj.N/2-1)*obj.dx); 
            aperture = circ(round(obj.width./obj.dx)/2,[obj.N,obj.N]); 
            obj.data = lensPhaseBasic(obj).*aperture;
            obj.roi = obj.dx*obj.N;
             
        end
    
        function obj = pinhole(obj,Phole) %lens struct
            %lens_class Construct an instance of this class
            %   Inputs: 
            % P - Parameters struct for general simulation parameters
            % Plens - Parameters struct for specific lens
            arguments
                obj;
                Phole = obj;
            end
            obj.width = Phole.diameter;
            obj.position = Phole.position; 
            obj.dx = Phole.dx;
            obj.N = Phole.N;
            aperture = circ(round(obj.width./obj.dx)/2,[obj.N,obj.N]); 
            obj.data = aperture;
             
        end
        function obj = conjugate(obj)
            obj.data = conj(obj.data);
        end
        function showme(obj,name)
            arguments
                obj
                name = sprintf('Optic at z = %f',obj.position);
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