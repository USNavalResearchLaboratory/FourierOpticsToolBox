classdef shackHartmannWavefrontSensor 
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        wvl = 650e-9; %wavelength
        N % Size of each lenslet
        dx
        X % X position meshgrid
        Y % Y position meshgrid
        width = 0; % width of object (not total width of array dx*N)
        roi % Region of interest. Equals dx*N
        position
        focus
        num_lenslets
        lenslet_spacing
        shape
        incident_field
        spotfield
        centroids
        slopes
        lenslet %2D data of lenslet
        data = []; %2D data of lenslet array
    end

    methods
        function obj = shackHartmannWavefrontSensor(Pshwfs)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            
            % Assign parameters
            obj.width = Pshwfs.array.diameter;
            obj.focus = Pshwfs.lenslet.focus ;
            obj.position = Pshwfs.lenslet.position; 
            obj.dx = Pshwfs.lenslet.dx;
            obj.N = Pshwfs.array.N;
            obj.wvl = Pshwfs.lenslet.wvl; 
            obj.num_lenslets = Pshwfs.array.num;
            obj.lenslet_spacing = Pshwfs.array.spacing;
            obj.shape = Pshwfs.array.shape; 
            
            % generate single lens
            obj.lenslet = thinLens(Pshwfs.lenslet); 
            % Make lenslet array
            obj = obj.makeLensletArray(Pshwfs.array);
            
            obj.N = size(obj.data,1); 
            % Calculate remaining parameters
            [obj.X,obj.Y] = meshgrid((-obj.N/2 : obj.N/2-1)*obj.dx); 
            %aperture = circ(round(obj.width./obj.dx)/2,[obj.N,obj.N]); 
            obj.roi = obj.dx*obj.N;
            
            
             
            
        end

        function obj = makeLensletArray(obj, Parray)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here 
            data_tmp = repmat(obj.lenslet.data,[Parray.num,Parray.num]); 
            %data = padarray(data_tmp, [ ]) 
            obj.data = data_tmp;
            
        end

        function obj = applyShackHartmannFull(obj, field, fpa)
            % Apply Shack Hartmann on Efield input. Calls other methods to
            % combine on SH functions into one method
            % Inputs: SH and Efield and focal plane array objects. 
            % Output: Centroids on FPA
            arguments
                obj
                field
                fpa = [] %Currently unused. Placeholder for camera specs or noise.
                        %Best way to use it is make a camera object, then
                        %replace 'propagate_ang_spec' to 'propagateTo' 
            end

            field = field.resample(obj); obj.incident_field = field; 
            field = field.apply_optic(obj);
            %field.showme('appliedlenslet')
            field = field.propagate_ang_spec(obj.focus); 
            obj.spotfield = field; 
            
            obj = obj.spotfield_segmentation; 
            obj = obj.centroiding; 
            
        end

        function obj = applyLensletArray(obj, field)
            % Apply Shack Hartmann on Efield input. Calls other methods to
            % combine on SH functions into one method
            % Inputs: SH and Efield and focal plane array objects. 
            % Output: Spotfield at FPA
            arguments
                obj
                field
               
            end

            field = field.resample(obj); obj.incident_field = field; 
            field = field.apply_optic(obj);
            %field.showme('appliedlenslet')
            field = field.propagate_ang_spec(obj.focus); 
            obj.spotfield = field; 
     
            
        end

        
        function obj = spotfield_segmentation(obj)
        % segment the spotfield. Input should have obj.spotfield. Output is data which is of size
        % [lenslet.N, lenslet.N, num_lenslets, num_lenslets]. Ex:
        % obj.data(:,:,1,1) would give you the top-left spot data. use
        % centroiding on this data next
            spotfield_intensity = abs(obj.spotfield.data).^2; 
            box_width = obj.N/obj.num_lenslets;
            box_start = [1,1]; %box_end = box_start + box_width;
            box_n = zeros([box_width, box_width,obj.num_lenslets,obj.num_lenslets]);
            for x_lenslet = 1:obj.num_lenslets-1
                for y_lenslet = 1:obj.num_lenslets-1
                pos = [box_start(1)+x_lenslet*box_width, box_start(2)+y_lenslet*box_width, box_width,box_width]; 
                box_n(:,:,x_lenslet, y_lenslet) = spotfield_intensity(pos(1):pos(1)+pos(3)-1, pos(2):pos(2)+pos(3)-1);
                end
            end
            obj.data = box_n; 
            %figure(); imagesc(spotfield_intensity); 
            %pos = [box_start(1)+obj.num_lenslets*box_width/2, box_start(2)+obj.num_lenslets*box_width/2, box_width,box_width]; 
            %hold on; drawrectangle("Position",pos);
        end
        function obj = centroiding(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here 
    
            %find centroid
            threshold = mean(obj.data(:))/2; 
            centroids = zeros([obj.num_lenslets, obj.num_lenslets, 2]);
            for x_lenslet = 1:obj.num_lenslets-1
                for y_lenslet = 1:obj.num_lenslets-1
                    tmp = obj.data(:,:,x_lenslet,y_lenslet);
                    if mean(tmp,'all') > threshold
                        centroid = regionprops(true(size(tmp)),tmp, 'WeightedCentroid');
                        centroids(x_lenslet,y_lenslet,:) = centroid.WeightedCentroid;
                    end
                end
            end
            obj.centroids = centroids; 

        end
        function signal = calculateSlopes(signal,ref)
            % Inputs are two ShackHartmann objects with centroiding done. One from ref and one
            % from signal
            signal.slopes = signal.centroids - ref.centroids;
            % Set large slopes to zero
            signal.slopes(abs(signal.slopes)>10) = 0; 

            

        end
        function obj = phaseReconstruction(obj, gpuFlag)
            arguments
                obj
                gpuFlag = 0
            end
            %Phase reconstruction use slopes 
          
            if gpuFlag
                obj.data = southwell_slopes_to_phase(gpuArray(obj.slopes), obj.num_lenslets, gpuFlag);
            else
                obj.data = southwell_slopes_to_phase(obj.slopes, obj.num_lenslets);
            end
            aperture = abs(obj.slopes(:,:,1)) > 0; 
            obj.data = obj.data.*aperture; %Set extra stuff to zero
            obj.data = imresize(obj.data,size(obj.X), 'cubic'); % Interpolate to original array size
            obj.data = exp(1i.*obj.data); %put phase into complex part of data
            

        end
        
        
        function showme(obj,name)
            arguments
                obj
                name = sprintf('Wavefront at z = %f',obj.position);
            end
            fig = figure("Position",[100,100,600,300], "Name",name); 
            limits = [-obj.width obj.width -obj.width obj.width].*1e3/2;

            if ~isreal(mean(obj.data(:)))
                
                
                subplot(1,2,1); imagesc(obj.X(1,:)*1e3, obj.Y(:,1)*1e3, abs(obj.data(:,:,1)).^2); title('Intensity');axis image ; xlabel('[m]')
                xlabel('[mm]'); ylabel('[mm]'); axis(limits)
                subplot(1,2,2); imagesc(obj.X(1,:)*1e3, obj.Y(:,1)*1e3, angle(obj.data(:,:,1))); title('Phase'); axis image
                xlabel('[mm]'); ylabel('[mm]'); axis(limits)
                sgtitle(name);
                %figureFormats(fig);
            else 
                imagesc(obj.X(1,:)*1e3, obj.Y(:,1)*1e3, abs(obj.data(:,:,1)).^2); title('Intensity');axis image ; xlabel('[m]')
                xlabel('[mm]'); ylabel('[mm]'); axis(limits)
                sgtitle(name);

            end

            
        end

        function showSlopes(obj,name)
            arguments
                obj
                name = sprintf('Shack-Hartmann slopes');
            end
            fig = figure( "Name",name); 
            quiver(obj.slopes(:,:,1),obj.slopes(:,:,2))
            %figureFormats(fig);
            
        end
    end
end