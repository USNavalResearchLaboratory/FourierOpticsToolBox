classdef target_class
    %target_class Make a target object 
    %   Use properties from Parameter struct
    %   Call createTargets to make a target
 

    properties
        targetType  % Type of target
        width       % Width of target
        roi         % Region of interest. Equals dx*N. Included padded area
        distance    % Distance of target from source
        position    % Position of target
        N           % Number of points in x and y direction
        wvl         % Wavelength
        dx          % Sampling interval
        X           % X coordinates of target points
        Y           % Y coordinates of target points
        avg_power   % Average power 
        peak_power  % Peak power
        num_speckles% Number of speckles
        data        % Data generated for target
    end
    
    methods
        function obj = target_class(target)
            %target_class Construct an instance of the target_class
            %   Initializes the target_class object based on the Parameter
            %   structure P.
            %   Parameters:
            %   P - Parameter structure containing the targetType, N,
            %   width, dx, distance, num_speckles, wvl
            %
            %   Returns:
            %   obj - An instance of target_class
            obj.targetType = target.targetType;
            obj.N = target.N; 
            obj.width = target.D; 
            obj.dx = target.dx; %obj.width/obj.N; 
            [obj.X,obj.Y] = meshgrid((-obj.N/2 : obj.N/2-1)*obj.dx); 
            obj.distance = target.distance;
            obj.num_speckles = target.specks; 
            obj.wvl = target.wvl; 
            obj.position = target.position; 
            obj.roi = obj.dx * obj.N;
          
        end

        function obj = genTarget_points(obj)
            r = 80;%P.T.D/P.T.Fpx; 
            window_pct = 50; %Spread of points in percent
            %TODO change to N
            obj.data = randguass([obj.N,obj.N],r,60,window_pct*obj.N/100,1);
        end

        function obj = genTarget_gauss(obj, r, power)
            arguments
                obj
                r = obj.width/(obj.dx);
                power = .01; 
            end
            obj.avg_power = power; %Watts
            avg_I = obj.avg_power/obj.roi; %Average intensity over region of interest
            max_I = avg_I*2; %peak intensity is just 2 times avg? 
            max_Efield = sqrt(max_I);
            obj.data = customgauss([obj.N,obj.N],r,r,0,0,1,[0,0])*max_Efield; %Since customgauss is normalized to 1, scale to max_I
            obj = obj.calc_power; 
        end
        function obj = genTarget_simplegauss(obj, r, power)
            arguments
                obj
                r = obj.width/(50000*obj.dx);
                power = .026; 
            end
            obj.avg_power = power; %Watts
            avg_I = obj.avg_power/obj.roi; %Average intensity over region of interest
            max_I = avg_I*2; %peak intensity is just 2 times avg? 
            max_Efield = sqrt(max_I);
            simple_gauss = exp( - (obj.X.^2 + obj.Y.^2) / r^2 );
            obj.data = simple_gauss*max_Efield; %Since customgauss is normalized to 1, scale to max_I
            obj = obj.calc_power; 
        end
        function obj = genTarget_focus_spot(obj,lens, aberration, gpuFlag)
            % genTarget_focus_spot Make an airy disk pattern at the target.
            % use lens object as the outgoing field (top hat beam with lens
            % phase focusing at the target). 
            % Also include aberration in the path. Could be air 
            % INPUTS
            %   - obj: initialized target class
            %   - lens: objective lens with focus at target
            %   - aberration: phase_screens class
            % OUTPUTS
            %   - obj: Efield class object of focused spot on target
            %   position
            arguments
                obj
                lens
                aberration = []
                gpuFlag = 0;
            end
            optimization = 0; 
            % Create gaussian pattern with lens
            output_field = lens; 
            U = Efield(output_field, gpuFlag); 
            % propagate to target to make point source 
            obj = U.propagate_ang_spec(U.position-obj.position,'rev',obj.width,aberration, optimization);
            [peak, ~,widths] = findpeaks(abs(obj.data(floor(obj.N/2),:))); 
            [~,ind] = max(peak); 
            obj.width = gather(widths(ind))*obj.dx; 
            
        end

        function obj = genTarget_ptsrc(obj,observation_aperture)
            % Generates a point source  for the target
            x1 = obj.X; y1 = obj.Y; 
            [theta1, r1] = cart2pol(x1, y1);
            Dz = obj.distance; 
            arg = observation_aperture.width/(obj.wvl*Dz);
            k = 2*pi/obj.wvl; 
            pt = exp(-1i*k/(2*Dz) * r1.^2) * arg^2 ...
            .* sinc(arg*x1) .* sinc(arg*y1);

            obj.data = pt; 
        end

        function obj = genTarget_threeBar(obj)
            % Create the three Bar target
           
            wide = obj.width/obj.dx/2; %[pixels]
            height = obj.width/obj.dx/5; %[pixels]
            HpxR = round(height);
            WpxR = round(wide);
            % separation in pixels
            Spx = 2*HpxR;
            amp = zeros(obj.N);
            for kk = -2:2:2
                temp = zeros(obj.N);
                temp(obj.N/2 - HpxR + kk*Spx : obj.N/2 + HpxR - 1 +kk*Spx,...
                    obj.N/2 - WpxR : obj.N/2 + WpxR - 1) = 1.0;
                amp = amp + temp;
                %targets = smoothAndRandSpeckReals(P, targets, amp);
            end
            obj.data = amp; 
        end
        function obj = addSpeckle(obj,smoothing)
            % Add speckle to the target. Smoothing is smoothing factor [integer]. 1
            % means no smoothing. Higher = more smoothing
            % This works by creating a smaller 2D array of random values,
            % then upsamples it to the size of the desired array. Is it the
            % best way? Need to research accurate ways to model speckle 
            arguments
                obj
                smoothing = 10
            end

            speckleReal = (rand(round(obj.N/smoothing))-0.5)*2*pi;
            if smoothing ~= 1
                [ Xtmp, Ytmp ] = genGrids(obj.N/smoothing, obj.N/smoothing, obj.dx*smoothing, obj.dx*smoothing);
                %[ Xq, Yq ] = genGrids(obj.N, obj.N, obj.dx, obj.dx);
                % 
                speckleReal = interp2(Xtmp,Ytmp,speckleReal,obj.X,obj.Y, 'cubic', 0);
            end
            obj.data = obj.data.*exp(1i*speckleReal); 

        end
        function obj = calc_power(obj)
            avg_I = mean(abs(obj.data(:)).^2); 
            obj.avg_power = avg_I*obj.roi; 

            peak_I = max(abs(obj.data(:)).^2); 
            obj.peak_power = peak_I*obj.roi; 
        end
        function showme(obj)
            figure("Position",[100,700,600,400]); 
            limits = [-obj.width obj.width -obj.width obj.width].*1e3/2; 
            subplot(1,2,1); imagesc(obj.X(1,:)*1e3, obj.Y(:,1)*1e3,abs(obj.data(:,:,1)).^2); title('Intensity');axis image
            xlabel('[mm]'); ylabel('[mm]'); axis(limits)
            subplot(1,2,2); imagesc(obj.X(1,:)*1e3, obj.Y(:,1)*1e3,angle(obj.data(:,:,1))); title('Phase'); axis image
            xlabel('[mm]'); ylabel('[mm]'); axis(limits)
            sgtitle('Target')
            fontsize(gcf,18, 'pixels')
        end
    end
end