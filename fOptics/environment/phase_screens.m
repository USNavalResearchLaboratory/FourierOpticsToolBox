classdef phase_screens
    % phase_screens: This class creates an object of 'phase screen'
    %
    % Properties:
    %   num_screens: The number of screens to be generated.
    %   propagation_steps: The number of propagation steps, where
    %   num_screens < propagation_steps.
    %   Cn2: The structure coefficient.
    %   type:  The type of phase screen that will be generated.
    %   data: 3D Array. The turbulent screen data.
    %   N:  The size of the phase screen.
    %   dx:  The pixel size of the phase screen.
    %   X: 2D Array. The X values of the phase screen.
    %   Y: 2D Array. The Y values of the phase screen.
    %   r0:  The Fried parameter of the turbulence.
    %   kolm_zerns:  The number of Kolmogorov Zernike terms.
    %   diameter:  The diameter of the aperture.
    %   wvl: The wavelength of the light.
    %   Z: The propagation distance between phase screens.
    %
    % Methods:
    %   phase_screens: Constructor. Initializes the object with default
    %   values and constructs a flat phase screen if none are given.
    %   Arguments: P: Struct. The struct containing default values for the
    %   object. num_screens: Integer. The number of screens to be generated.
    %
    %   kolmogorov: Creates individual screens based on Kolmogorov
    %   turbulence. Replaces flat phase screen object with turbulence.
    %
    %   moving_kolmogorov: Creates moving Kolmogorov screen using
    %   genTurb_flow function.
    %
    %   tilt: Creates a tilted phase screen. Arguments: anglex_deg: Scalar.
    %   The angle of tilt in the X direction in degrees. angley_deg: Scalar.
    %   The angle of tilt in the Y direction in degrees.
    %
    %   OAM: Creates an Orbital Angular Momentum phase screen.
    %
    %   zernike: Creates a Zernike phase screen. Arguments: index: Scalar.
    %   The index of the Zernike polynomial.
    %
    %   showme: Shows the generated turbulence screens.

    properties
        num_screens % The number of screens to be generated.
        propagation_steps % The number of propagation steps, where
        % num_screens < propagation_steps.
        Cn2 % The index of refraction structure coefficient.
        type % The type of phase screen that will be generated.
        N %  The size of the phase screen.
        dx %  The pixel size of the phase screen.
        X % 2D Array. The X values of the phase screen.
        Y % 2D Array. The Y values of the phase screen.
        r0_Z  %  The Fried parameter of the turbulence over total propagation distance
        r0_dZ  %  r0 between phase screens
        r0_px  %  r0 between phase screens in pixels
        kolm_zerns % Vector. The Zernike coefficients for the Kolmogorov spectrum.
        roi %  The diameter of the aperture.
        wvl %  The wavelength of the light.
        Z % The propagation distance.
        zernike_distribution % Array of magnitudes for zernikes. Index of array = zernike number
        data % 3D Array. The turbulent screen data.
    end

    methods

        function obj = phase_screens(P)
            %Initialize: Construct a flat phase screen using given values
            %if none given
            arguments
                P 
                %num_screens = P.screens
            end
            % Default screen is flat 'air' 
            if ~isfield(P, 'screens') 
                P.screens = 2; %default 2 screens if not defined
            end
            obj.num_screens = P.screens;
            obj.type = 'flat';
            obj.data = ones(P.N, P.N, obj.num_screens);
            obj.N = P.N;
            obj.Cn2 = 0; 
            obj.dx = P.dx; 
            [obj.X, obj.Y] = meshgrid((-obj.N / 2:obj.N / 2 - 1) * obj.dx);
            obj.roi = obj.dx * obj.N;
            obj.wvl = P.wvl;
            %obj.Z = P.position;

        end

        function obj = kolmogorov(obj, Turb, gpu)
            %kolmogorov Creates individual screens based on kolmogorov
            %statistics
            %  Replace flat phase screen object with turbulence
            arguments
                obj
                Turb
                gpu = 0;
            end
            obj.Cn2 = Turb.Cn2; 
            obj.r0_Z = Turb.r0;
            obj.r0_dZ = Turb.r0delta; 
            obj.r0_px = Turb.r0delta/obj.dx;
            obj.kolm_zerns = Turb.K;
            %obj.dx = Turb.dx; %TODO, this isnt exactly right. Need to know what dx field will be at location of field.
            obj.num_screens = Turb.screens; 
     
            disp('Generating Turbulence screens...')
            screens_tmp = ones(obj.N, obj.N, obj.num_screens);

            %tmp variables for efficient parfor
            Ntmp = obj.N; r0_pxtmp = obj.r0_px; kolm_zernstmp = obj.kolm_zerns; 
            tic
            if gpu == 1
                screens_tmp = gpuArray(screens_tmp);
                for scr = 1:obj.num_screens
                        screens_tmp(:, :, scr) = turbulence_phi_fftGPU(Ntmp, r0_pxtmp,kolm_zernstmp);
                        screens_tmp(:, :, scr) = exp(1i * (screens_tmp(:, :, scr) - mean(mean(screens_tmp(:, :, scr))))); % exponentiate for E-field multiplication
                end
            else
                for scr = 1:obj.num_screens
                    screens_tmp(:, :, scr) = turbulence_phi_fft(Ntmp, r0_pxtmp,kolm_zernstmp);
                    screens_tmp(:, :, scr) = exp(1i * (screens_tmp(:, :, scr) - mean(mean(screens_tmp(:, :, scr))))); % exponentiate for E-field multiplication
                end
            end
            obj.data = screens_tmp;
            toc


        end

        function obj = moving_kolmogorov(obj,Turb,z,flowspeed, frames)
            %moving_kolmogorov screen using genTurb_flow function
            arguments
                obj
                Turb
                z = 1000
                flowspeed = 1; % units of index - phase generates per line
                frames = 1; % number of frames of flow
            end
            obj.Cn2 = Turb.Cn2; 
            obj.r0_Z = Turb.r0;
            obj.r0_dZ = Turb.r0delta; 
            obj.r0_px = Turb.r0delta/obj.dx;
            obj.kolm_zerns = Turb.K;
            obj.num_screens = 1; %Turb.screens; 

            r0 = NAYRU.Turb.r0screenCalc(obj.Cn2,obj.wvl,z,'N',obj.num_screens);
            phase_screens = zeros(obj.N,obj.N,obj.num_screens+1,flowspeed*frames); %all turbulent realizations

            for i = 1:obj.num_screens
                phase_screens(:,:,i+1,:) = NAYRU.Turb.InfKolmogorov(obj.N,'N',flowspeed*frames,'dx',obj.dx,'r0',r0);
                            
            end
            obj.data = phase_screens; 
        end

        function obj = tilt(obj, anglex_deg, angley_deg)
            arguments
                obj
                anglex_deg = .05 %Define angle tilt if none given
                angley_deg = 0
            end

            anglex_rad = deg2rad(anglex_deg);
            angley_rad = deg2rad(angley_deg);
            %tilt equation: 2*pi/wvl * optical_delay where optical_delay =
            %X * tan(angle). This is the phase delay as a function of
            %position from center due to tilted plane wave. 
            x_tilt = (2*pi/obj.wvl)*obj.X.*tan(anglex_rad); 
            y_tilt = (2*pi/obj.wvl)*obj.X.*tan(angley_rad);
             
            scr = round(obj.num_screens / 2); %position of aberration along propagation path
            totaltilt = x_tilt + y_tilt';
            obj.data(:, :, scr) = exp(1i .*totaltilt);
        end

        function obj = OAM(obj,charge,oam)
            %OAM Generate a phase screen with a spiral phase pattern
            % obj = OAM(obj) generates a spiral phase pattern (OAM mode) for the
            % electric field obj. The function uses the oam_phasescreen function to
            % generate the phase screen with the specified charge, index, and the same
            % dimensions (N, dx) as the field obj. The spiral phase pattern is then
            % applied to the electric field.
            %
            % Inputs:
            % obj - The electric field to apply the OAM phase pattern to
            %
            % Outputs:
            % obj - The updated electric field with the OAM phase pattern applied
            arguments
                obj
                
                charge = 2; % charge of the spiral phase pattern
                oam = struct('Z', 1000)
            end
            index = 2; % index of the spiral phase pattern
            scr = 2; % the index of the plane where the spiral phase pattern is applied
            obj.Z = oam.Z; 
            % Generate the spiral phase pattern using the oam_phasescreen function and
            % apply it to the electric field
            obj.data(:, :, scr) = exp(1i .* angle(oam_phasescreen(obj, charge, index)));
        end

        function obj = zernike(obj, indices, randomize)
            %ZERNIKE Add a Zernike aberration to the wavefront
            % obj = ZERNIKE(obj, index) adds the Zernike aberration corresponding to the
            % given index to the wavefront obj. The aberration is added to the central
            % screen of the wavefront, which is taken to be the screen at the center of the
            % array obj.data.
            %
            % Inputs:
            % obj - The wavefront to add the aberration to
            % indices - The indices of the Zernike aberrations to add. Its
            % a array [z1 z2 z3 ...]
            %
            % Outputs:
            % obj - The updated wavefront with the Zernike aberration added to the central screen
            % Compute Zernike polynomials up to fifth order
            nOrder = 5;
            [~, zernikes, ~] = zernikePolynomials(obj.N, obj.N, nOrder);

            % Add the Zernike aberration to the central screen
            scr = round(obj.num_screens / 2);
            
            scale = zeros([1 size(zernikes,3)]); 
            scale(indices) = 3; 
            if randomize
                rand_coeffs =  1 - 2*rand([1 size(zernikes,3)]); 
                scale = scale.*rand_coeffs;
            end
            normalization = max(max(zernikes,[],1),[],2);
            zernikes_normalized = zernikes./normalization; 
            zernikes_normalized = reshape(zernikes_normalized,size(zernikes_normalized,1)^2,[]) ;

            turbulence = zernikes_normalized.*scale; 
            turbulence = reshape(turbulence,size(zernikes,1),size(zernikes,2),[]) ;
            turbulence = sum(turbulence,3);
            obj.data(:, :, scr) = exp(1i .* turbulence);
            obj.zernike_distribution = scale; 

        end
        
         function obj = resample(obj,obj2)
            %RESAMPLE Resample data in obj to match the size and resolution of obj2.
            %
            %   obj = RESAMPLE(obj, obj2) resamples the data in obj to have the same 
            %   size and resolution as obj2. The method pads or crops obj.roi to match 
            %   obj2.roi, and then resamples the data in obj to have the same size N 
            %   and same dx. If obj.dx is not the same as obj2.dx, the data is
            %   resampled using bilinear interpolation.
            %
            %   Inputs:
            %       - obj: an instance of the Efield class, representing a 2D electric 
            %              field data with properties: data, X, Y, dx, roi, N.
            %       - obj2: an instance of the Efield class with the target properties
            %               to resample obj to.
            %
            %   Outputs:
            %       - obj: an instance of the phase screens class, containing the resampled 
            %              data and updated properties: data, X, Y, dx, roi, N.
            %
            roi_diff = obj2.roi - obj.roi; %Difference between region of interests 

            if roi_diff > 1e-3
                % pad obj.N to obj2.N
                npad = round((obj2.roi - obj.roi) / obj.dx);
                pad_left = floor(npad / 2);
                pad_right = ceil(npad / 2);
                for t = 1:obj.num_screens
                    tmp1 = padarray(obj.data(:, :, t), [pad_left pad_left], 'pre');
                    tmp(:, :, t) = padarray(tmp1, [pad_right pad_right], 'post');
                end
                obj.data = tmp; 
                obj.roi = obj.dx * size(obj.data, 1);
            elseif roi_diff < -1e-3
                % crop obj.N to obj2.N
                ncrop = round((obj.roi - obj2.roi) / obj.dx);
                crop_left = floor(ncrop / 2);
                crop_right = ceil(ncrop / 2);
                for t = 1:obj.num_screens
                    tmp(:, :, t) = obj.data(crop_left+1:end-crop_right, crop_left+1:end-crop_right,t);
                end
                obj.data = tmp; 
                obj.roi = obj.dx * size(obj.data, 1);
            end

            % update obj to have same N, X, Y, and roi as obj2
            Ntmp= round(obj.roi / obj.dx);
            % resample obj.data to match obj2.data
            if obj.dx ~= obj2.dx
                tmp = []; 
                [X1, Y1] = meshgrid((-Ntmp/2 : Ntmp/2-1)*obj.dx);
                for t = 1:obj.num_screens
                    tmp(:, :, t) = interp2(X1, Y1, obj.data(:, :, t), obj2.X, obj2.Y, 'cubic');
                end
                
                obj.data = tmp; 
                obj.data(isnan(obj.data)) = 0; % Set NaNs to 0
                obj.dx = obj2.dx;
            end
            
            obj.X = obj2.X;
            obj.Y = obj2.Y;
            obj.N = size(obj.data,1); 
            

    end
        function showme(obj)
            figure("Name", "Phase screens");

            for t = 1:obj.num_screens
                subplot(1, obj.num_screens, t); imagesc(angle(obj.data(:, :, t))); axis image off
            end

            sgtitle('Sample turbulence screens')
            fontsize(gcf, 18, 'pixels')

        end

        function zernike_profile(obj)
             figure("Name", "Zernike Profile");
             bar(obj.zernike_distribution)
             axis([0 length(obj.zernike_distribution) -1 1])
             sgtitle('Zernike Profile')
             xlabel('Zernike mode'); ylabel('Amplitude')
        end

    end

end
