classdef diffraction_grating
    %UNTITLED10 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        
        N % Size of data array
        dx %assuming dx=dy. use symmetric optics
        X % X position meshgrid
        Y % Y position meshgrid
        width = 0; % width of object (not total width of array dx*N)
        roi % Region of interest. Equals dx*N
        position = 0; % Position of field along z axis
        grooves = 500; % grooves density per mm 
        wvl
        data = [];
    end

    methods
        function obj = diffraction_grating(Pgrating)
            %UNTITLED10 Construct an instance of this class
            %   This class creates a actual grating at the size of the
            %   field. The size of the array can be too big with this
            %   method if the field is big. [e.g. A grating with 500
            %   lines/mm and grating size of 10mm is an array of atleast
            %   N = 2*10*500 = 10,000 since you want atleast 2 pixels per line
            %   (Nyquist)]. 
            
            obj.grooves = Pgrating.grooves; 
            obj.width = Pgrating.width; 
            obj.N = Pgrating.N; 
            obj.dx = Pgrating.dx; 
            obj.roi = obj.dx*obj.N;
            [obj.X,obj.Y] = meshgrid((-obj.N/2 : obj.N/2-1)*obj.dx); 
            obj.wvl = Pgrating.wvl; 
            
        end

        function obj = diffraction_grating_physical(obj)
            %TODO Need to check
            mm = 1e-3; 
            Nmin = obj.grooves/mm*obj.width; %  number of grid points. Due to density of grating, large number of points needed. 
            if Nmin > obj.N
                disp('Grating size [N] too small')
            end
            total_grooves = obj.grooves/mm * obj.roi; %total grooves in roi
            N_per_groove = 2*floor(obj.N/total_grooves/2);
            %N_per_groove = round(obj.N/obj.grooves,TieBreaker="even"); % works if even
            
            %grating = [1 0 0 1]; %manual method. 
            %grating = ((-1).^(linspace(1,N_per_groove,N_per_groove))+1)/2; 
            grating = [zeros([1,ceil(N_per_groove/2)]), ones([1,ceil(N_per_groove/2)])]; 
            grating = repmat(grating, [1,floor(obj.N/N_per_groove)]); %grating should be symmetric [0 1 1 0]
            grating = repmat(grating, [obj.N,1]); %alternate 1 0 to make grating
            grating = grating(1:obj.N,1:obj.N); %Crop to same size as obj in case. assuming grating will end up bigger due to ceil
            obj.data = grating; 
            
        end

        function obj = diffraction_grating_phase(obj, Uin)
            %UNTITLED10 Construct an instance of this class
            %   Simple diffraction grating works by calculating the tilt
            %   phase light would pick up after a diffraction grating. This
            %   tiltphase can be applied to the field. 
            order = 2; %order of beam after grating 
            mm = 1e-3; 
            
            grating_spacing = mm/(obj.grooves); 
            ang = asin(order*Uin.wvl/grating_spacing);         
            x_tilt = (2*pi/obj.wvl)*obj.X.*tan(ang); 

            obj.data = Uin.data.*exp(1i.*x_tilt); 
        end

        function Uout = diffraction_grating_manual(obj,Uin, grating)
            % MANUAL METHOD Calculated the displacement shift after
            % diffraction grating. Uses imtranslate to shift field. Center
            % is defined at a single wavelength 

            %Input: Efield
            %       diff_grating object
            %Output: Shifted Efield
            mm = 1e-3; 
            m = 1; D = grating.DistanceToCam; grating_spacing = mm/grating.grooves; 
            displacement = m*Uin.wvl*D/grating_spacing; %displacement from grating [m]; 
            displacement_px = displacement/obj.dx; 
            % Calculate displacement of center wavelength. Subtract from
            % displacement. We assume the screen after diffraction grating
            % is centered at the diffraction angle from center wvl of
            % grating and order m. 
            wvl_center = grating.wvl;
            center = m*wvl_center*D/grating_spacing; %displacement from grating [m]; 
            center_px = center/obj.dx; 
            displacement_px = displacement_px - center_px; 
            %Uin.data = gather(Uin.data); %pull from gpu if data is there. imtranslate not supported in gpu
            %Uin.data = imtranslate(Uin.data,[displacement_px,0], 'FillValues',abs(mean(Uin.data,'all'))) ;
            avg = mean(abs(Uin.data),'all'); % get avg before shift
            outView = imref2d([Uin.N, Uin.N]); tform = transltform2d(displacement_px,0);
            Uin.data = imwarp(Uin.data,tform,'OutputView',outView);
            Uin.data(Uin.data == 0) = avg; %replace 0s in shifted image with average

            Uin.position = Uin.position + D;
            Uout = Uin; 
        end


        function showme(obj,name)
                arguments
                    obj
                    name = sprintf('Optic at z = %f',obj.position);
                end
                fig = figure("Position",[100,100,600,300], "Name",name); 
                subplot(1,2,1); imagesc(obj.X(1,:).*1e3, obj.Y(:,1).*1e3, abs(obj.data(:,:,1)).^2); title('Intensity');axis image ; xlabel('[m]')
                ylabel('[mm]')
                subplot(1,2,2); imagesc(obj.X(1,:).*1e3, obj.Y(:,1).*1e3, angle(obj.data(:,:,1))); title('Phase'); axis image
                xlabel('[mm]');
                sgtitle(name);
                figureFormats(fig);
                
            end
      
    end
end