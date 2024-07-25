
classdef Efield
    % The field class is the electric field at some position. Methods to
    % propagate the field is included
    % Initialize the field by passing in the target object. Can also pass in a struct, or class with the properties defined.
    % See individual methods for more help.


    properties
        domain = 'spatial'; % Could be spatial or frequency domain
        wvl  % wavelength
        N  % size of data array
        dx  % generally assuming dx=dy. use symmetric optics
        dy
        X % X position meshgrid
        Y % Y position meshgrid
        width  % width of object of interest within roi (not total width of array dx*N)
        roi % Region of interest. Equals dx*N. Included padded area
        position  % position of field along z axis
        avg_power   %Average power over roi
        peak_power %peak power
        optimizationFlag = 0
        data  % field data
        phase %Just the phase component, unwrapped
    end

    methods
        function obj = Efield(target, gpuFlag)
            arguments
                target = get_default_Efield_struct();
                gpuFlag = 0;
            end
            % Init field at using input data. Converts into Efield class.
            % Could be target class or experimental data for example.
            obj.wvl = target.wvl;
            obj.N = target.N;
            obj.dx = target.dx;
            obj.position = target.position;
            obj.width = target.width;

            %calculations
            obj.roi = obj.dx * obj.N;
            [obj.X, obj.Y] = meshgrid((-obj.N/2 : obj.N/2-1) * obj.dx);

            if gpuFlag
                obj.data = gpuArray(target.data);
            else
                obj.data = (target.data);
            end


            %obj = calc_power(obj);
        end

        function [obj,varargout] = propagateTo(obj, optic, direction, aberration, proptype, figflag)
            %PROPAGATETO Combines propagation of electric field,
            %resampling, and applying optic into one function.
            % obj = propagateTo(obj, optic, direction, aberration, proptype) or obj = obj.propagateTo(optic, ...)
            % obj is an Efield object with a defined 'position' property.
            % This function will propagate to the optic.position and apply
            % that optic.
            % direction = 'fwd', 'rev', or 'revMirror'
            % aberation = phase_screens class
            % proptype = 'angspec' or 'fresnel'. Note that fresnel
            % propagation is single step and does not take aberration into
            % account.
            arguments
                obj
                optic
                direction = 'fwd';
                aberration = phase_screens(phaseScreenDefault(obj)); %Returns initialization (air)
                proptype = 'angspec'; %propagation type
                figflag = 0
            end

            dz = abs(optic.position-obj.position);
            if proptype == 'angspec'
                fprintf('Propagating %f meters with Angular Spectrum Prop\n',dz)
                [obj, Uout_all] = obj.propagate_ang_spec(dz, direction, optic.width, aberration, obj.optimizationFlag, figflag);
            elseif proptype == 'fresnel'
                obj = obj.propagate_fresnel(dz,optic.width,optic.focus);
            end
            varargout{1} = Uout_all; 
            %obj.showme('Before apply optic')
            % Create optic at position with right dx
            obj = obj.resample2(optic);
            obj = obj.apply_optic(optic);
            %obj = calc_power(obj);
        end

        function obj = propagate_fresnel(obj, z, width2, R)
            %PROPAGATE_FRESNEL Propagate an electric field object by the Fresnel approximation.
            %   obj = PROPAGATE_FRESNEL(obj, width2, z) propagates the electric field obj by the
            %   Fresnel approximation to a distance z along the z axis, assuming an infinite radius
            %   of curvature (i.e. a plane wave). The parameter width2 specifies the width of the
            %   2nd aperture plane and is used to calculate the size of the propagated field. The
            %   optional parameter R specifies the radius of curvature of the wavefront, which
            %   can be used to simulate a converging or diverging beam.
            %
            %   Inputs:
            %       - obj: an instance of the Efield class, representing a 2D electric
            %              field data with properties: data, X, Y, dx, roi, N.
            %       - width2: the width of the 2nd aperture plane in meters.
            %       - z: the distance to propagate the field in meters.
            %       - R: the radius of curvature of the wavefront in meters (default: inf).
            %
            %   Outputs:
            %       - obj: an instance of the Efield class, containing the propagated
            %              data and updated properties: data, X, Y, dx, roi, N.
            %
            %   Example:
            %       % Create an Efield instance with some data
            %       obj = Efield();
            %       obj.dx = 0.02;
            %       obj.roi = 1;
            %       obj.N = round(obj.roi / obj.dx);
            %       obj.data = rand(obj.N, obj.N);
            %
            %       % Propagate the field to a distance of 10 cm
            %       obj = obj.propagate_fresnel(0.1, 0.1);

            arguments
                obj
                z
                width2 = obj.width
                R = inf;
            end

            [dx1_min,Nconstraint,Zconstraint, dx2 ] = fresnelIntegralSampling(obj.width, obj.wvl, z, obj.dx, width2, obj.N, R);


            % fresnel prop
            %obj.data = fresnel_prop(obj,z);
            [obj.X, obj.Y, obj.data] = two_step_prop(obj.data, obj.wvl, obj.dx, dx2, z);
            %[obj.X, obj.Y, obj.data] = ang_spec_prop(obj.data, obj.wvl, obj.dx, obj.dx,  z);
            obj.position = obj.position + z;
            obj.dx = abs(obj.X(1,1)-obj.X(1,2));
            obj.roi = obj.N*obj.dx;
            obj = calc_power(obj);
        end


        function [obj,varargout] = propagate_ang_spec(obj,z,direction, width2,aberration, optimization_flag, figflag)
            % propagate_ang_spec - Angular Spectrum Propagation method
            % Usage:
            %     obj = propagate_ang_spec(obj, z, direction, width2, aberration)
            % Inputs:
            %     obj        - A Efield object.
            %     z          - Propagation distance in meters.
            %     direction  - Direction of propagation ('fwd' for forward, 'rev' for reverse or 'revMirror' for reverse in the presence of mirror).
            %                 Default is 'fwd'.
            %     width2     - The width of the field at the propagated plane in meters.
            %                 Default is obj.width.
            %     aberration - Struct of the aberration to be applied on the field during the propagation.
            %                 Default is a struct with a matrix of ones and num_screens = 1.
            % Outputs:
            %     obj        - A Efield object after angular spectrum propagation.
            %     old_obj (through varargout)    - A resampled Efield input object
            % The function performs the Angular Spectrum Propagation method on the input Field2D object for a given propagation distance (z). The propagated field is returned as another Field2D object. It can be propagated in either forward or reverse direction and the user can also select the direction of propagation in the presence of a mirror. The user can also provide the width of the field at the propagated plane and an aberration struct to be applied during the propagation.
            % If the user does not provide the field width (width2) or the aberration struct, the default values of obj.width and a struct with a matrix of ones and num_screens = 1, respectively, will be used.
            % The function first checks whether the field can be propagated using the Angular Spectrum Propagation method by checking the sampling criteria. If the sampling criteria are not met, the function raises an error.
            % After propagation, the function updates the properties of the output Field2D object, including the field data, position, X and Y arrays, and dx. The output Field2D object can be further used for various operations.

            arguments
                obj
                z
                direction = 'fwd';
                width2 = obj.width;
                aberration = phase_screens(struct(obj)); %Use 2 screens of air
                optimization_flag = obj.optimizationFlag;
                figflag = 0
            end
            if isempty(aberration)
                %aberration = struct('data', ones([obj.N, obj.N,2]), 'num_screens', 2, 'Cn2', 0, 'N', obj.N);
                aberration = phase_screens(struct(obj)); %Use 2 screens of air
            end
            n = aberration.num_screens;

            % If there is aberration, calculate effective apertures D1p,
            % D2p for sampling checks
            k = 2*pi/obj.wvl;
            r0sw = (0.423 * k^2 * aberration.Cn2 * 3/8 * z)^(-3/5); %TODO Check r0 of aberration
            D1p = obj.width + 2*obj.wvl*z/r0sw;
            D2p = width2 + 2*obj.wvl*z/r0sw;

            % Run propagation sampling checks
            R = z;
            % Run sampling check
            %dx2 = (obj.wvl * z - D2p * obj.dx) / (2*obj.width); % Set reciever sampling
            %dx2 = set_receiver_sampling(obj.wvl, z, D2p, obj.dx, obj.width, obj.N);
            dx2 = obj.dx;
            [samplingCheck] = angSpecPropSampling(obj.wvl, z, D2p, obj.dx, D1p, dx2, obj.N, R, z/n);

            if ~samplingCheck && optimization_flag
                %Get new sampled obj
                [obj, dx2,aberration] = optimize_sampling_between(obj,D1p,D2p,dx2,z,aberration);
            end

            if obj.N ~= aberration.N
                disp('Resampling aberrations')
                aberration = aberration.resample(obj);
            end
            [samplingCheck] = angSpecPropSampling(obj.wvl, z, width2, obj.dx, obj.width, dx2, obj.N, R, z/n);
            % If sampling is successful, run propagation
            if samplingCheck ~= true
                disp('Check Sampling. Propagation may be incorrect')
            end
            [Xout, Yout, Uout, Uout_all] = fullProp2(obj, direction,z,[], aberration, dx2,figflag);
            varargout{1} = Uout_all; 
            % Update output field values
            obj.X = Xout;
            obj.Y = Yout;
            obj.dx = abs(Xout(1,1)-Xout(1,2));
            obj.data = Uout;
            obj.roi = obj.N*obj.dx;
            obj.width = width2;
            if strcmp(direction, 'fwd')
                obj.position = obj.position + z;
            elseif strcmp(direction, 'rev')
                obj.position = obj.position - z;
            elseif strcmp(direction, 'revMirror')
                obj.position = obj.position - z;
            end
            %obj = calc_power(obj);
        end


        function obj = apply_optic(obj, optic)
            %APPLY_OPTIC Apply an optic to the electric field
            %   obj = APPLY_OPTIC(obj, optic) applies the given optic to the electric
            %   field obj. The optic must have the same dimensions (N, dx), position
            %   along the z axis, and wavelength as the field.
            %
            %   Inputs:
            %       obj - The electric field to apply the optic to
            %       optic - The optic to apply to the field
            %
            %   Outputs:
            %       obj - The updated electric field with the optic applied

            % Check that the optic is compatible with the field
            if obj.N ~= optic.N || obj.dx ~= optic.dx
                error('Field and optic do not have the same dimensions');
            end

            % Check that the optic and field have the same position along the z axis
            if obj.position ~= optic.position
                error('Field and optic do not have the same position');
            end

            % Check that the optic and field have the same wavelength
            %             if obj.wvl ~= optic.wvl
            %                 error('Field and optic do not have the same wavelength');
            %             end

            % Apply the optic
            obj.data = obj.data .* optic.data;
        end


        function obj = resample(obj,obj2)
            disp('** DEPRECATED. USE resample2() **')
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
            %       - obj: an instance of the Efield class, containing the resampled
            %              data and updated properties: data, X, Y, dx, roi, N.
            %
            %   Example:
            %       % Create two Efield instances with different resolution and size
            %       obj1 = Efield();
            %       obj1.dx = 0.05;
            %       obj1.roi = 1;
            %       obj1.N = round(obj1.roi / obj1.dx);
            %       obj1.data = rand(obj1.N, obj1.N);
            %
            %       obj2 = Efield();
            %       obj2.dx = 0.02;
            %       obj2.roi = 1;
            %       obj2.N = round(obj2.roi / obj2.dx);
            %       obj2.data = rand(obj2.N, obj2.N);
            %
            %       % Resample obj1 to match obj2 size and resolution
            %       obj1 = obj1.resample(obj2);
            % WARNING: The interpolation, interp2, can be computationally
            % intensive and may use a lot of RAM
            roi_diff = obj2.roi - obj.roi; %Difference between region of interests
            gpuFlag = 0;

            if roi_diff > 1e-3 %if difference between roi is bigger than 1mm (obj2 bigger than obj)
                npad = round((obj2.roi - obj.roi) / obj.dx);
                pad_left = floor(npad / 2);
                pad_right = ceil(npad / 2);

                % Warning and option to proceed if the padding size is too large
                if (pad_left > 1e4) || (pad_right > 1e4)
                    warning('Large padding sizes detected. This might consume a lot of memory.');
                    response = input('Do you want to continue? (yes/no)', 's');
                    if strcmpi(response, 'no')
                        return;
                    end
                end

                % Incremental padding if pad size is large
                padIncrement = 500; % Adjust this based on what's optimal in your environment.
                if pad_left > padIncrement
                    for k = 1:padIncrement:pad_left
                        obj.data = padarray(obj.data, [min(padIncrement, pad_left - k + 1) 0], 'pre');
                    end
                else
                    obj.data = padarray(obj.data, [pad_left 0], 'pre');
                end

                if pad_right > padIncrement
                    for k = 1:padIncrement:pad_right
                        obj.data = padarray(obj.data, [min(padIncrement, pad_right - k + 1) 0], 'post');
                    end
                else
                    obj.data = padarray(obj.data, [pad_right 0], 'post');
                end

                obj.roi = obj.dx * size(obj.data, 1);
            elseif roi_diff < -1e-3 %if difference between roi is bigger than 1mm (obj bigger than obj2)
                % crop obj.N to obj2.N

                ncrop = round((obj.roi - obj2.roi) / obj.dx);
                crop_left = floor(ncrop / 2);
                crop_right = ceil(ncrop / 2);
                obj.data = obj.data(crop_left+1:end-crop_right, crop_left+1:end-crop_right);
                obj.roi = obj.dx * size(obj.data, 1);
            end

            % update obj to have same N, X, Y, and roi as obj2
            Ntmp= round(obj.roi / obj.dx);
            % resample obj.data to match obj2.data
            if obj.dx ~= obj2.dx
                [X1, Y1] = meshgrid((-Ntmp/2 : Ntmp/2-1)*obj.dx);
                obj.data = interp2(X1, Y1, obj.data, obj2.X, obj2.Y, 'cubic'); %spline was found to be best looking. cubic works with gpu
                obj.data(isnan(obj.data)) = 0; % Set NaNs to 0
                obj.dx = obj2.dx;
            end

            obj.X = obj2.X;
            obj.Y = obj2.Y;
            obj.width = obj2.width;
            obj.N = size(obj.data,1);

            if gpuFlag
                obj.data = gpuArray(obj.data); %Put array back to gpu if it was taken off
            end


        end

        function obj = resample2(obj, obj2)
            % updated resample function to fix edge case where memory
            % became large due to padarray 
            gpuFlag = 0;

            % Step 1: Adjust resolution
            if obj.dx ~= obj2.dx
                % Create grids for source and target based on obj's roi
                [X1, Y1] = meshgrid(linspace(-obj.roi/2, obj.roi/2, obj.N));
                newN = round(obj.N * (obj.dx/obj2.dx));
                [X2, Y2] = meshgrid(linspace(-obj.roi/2, obj.roi/2, newN));

                % Resample data using interp2 to match resolution
                obj.data = interp2(X1, Y1, obj.data, X2, Y2, 'cubic', 0); % use NaN for values outside the grid
                obj.dx = obj2.dx;
            end

            % Step 2: Handle ROI differences
            roi_diff = obj2.roi - obj.roi;
            diff_px = round(roi_diff / obj2.dx);  % Difference in pixels

            if roi_diff > 0
                % Padding
                pad_left = floor(diff_px / 2);
                pad_right = ceil(diff_px / 2);
                obj.data = padarray(obj.data, [pad_left, pad_left], 'pre');
                obj.data = padarray(obj.data, [pad_right, pad_right], 'post');
            elseif roi_diff < 0
                % Cropping
                crop_left = floor(abs(diff_px) / 2);
                crop_right = ceil(abs(diff_px) / 2);
                obj.data = obj.data(crop_left+1:end-crop_right, crop_left+1:end-crop_right);
            end

            % Step 3: Ensure final array size matches obj2.N
            if size(obj.data, 1) ~= obj2.N
                final_diff = obj2.N - size(obj.data, 1);
                pad_left = floor(final_diff / 2);
                pad_right = ceil(final_diff / 2);
                obj.data = padarray(obj.data, [pad_left, pad_left], 'pre');
                obj.data = padarray(obj.data, [pad_right, pad_right], 'post');
            end

            % Update other properties
            obj.roi = obj2.roi;
            obj.N = obj2.N;
            obj.X = obj2.X;
            obj.Y = obj2.Y;

            if gpuFlag
                obj.data = gpuArray(obj.data); % Put array back to GPU if it was taken off
            end


        end

        function obj = unwrap_phase(obj)
            wrapped_phase = gather(angle(obj.data));
            obj.phase = Unwrap_TIE_DCT_Iter(wrapped_phase);

        end
        function obj = magnify(obj,scale,cropFlag)
            % Magnifies Efield by amount. It does it by doing imresize and
            % then cropping. But width and sample size is same. Field will
            % be cropped if it becomes larger than aperture
            arguments
                obj
                scale
                cropFlag = 0
            end
            obj_temp = obj;
            %cropFlag=0; % Keep 1 to match final array size to original
            try
                obj_temp.data = (imresize(obj.data,scale));
                gpuFlag = 0;
            catch ME % when gpu is out of memory
                obj.data = gather(obj.data);
                obj_temp.data = (imresize(obj.data,scale));
                gpuFlag = 1;
            end
            if size(obj_temp.data,1) > size(obj.data,1)
                if cropFlag
                    win = centerCropWindow2d(size(obj_temp.data),size(obj.data));
                    obj_temp.data = imcrop(real(obj_temp.data),win) + 1i.*imcrop(imag(obj_temp.data),win); %imcrop doesn't work on imaginary numbers some reason
                    obj_temp.N = size(obj_temp.data,1);
                    obj_temp.width = obj.width*scale;
                    obj_temp.roi = obj_temp.N*obj_temp.dx;
                    [obj_temp.X, obj_temp.Y] = meshgrid((-obj_temp.N/2 : obj_temp.N/2-1) * obj.dx);
                else
                    obj_temp.N = size(obj_temp.data,1);
                    obj_temp.width = obj.width*scale;
                    obj_temp.roi = obj_temp.N*obj_temp.dx;
                    [obj_temp.X, obj_temp.Y] = meshgrid((-obj_temp.N/2 : obj_temp.N/2-1) * obj.dx);
                    % obj_temp = obj_temp.resample(obj);
                end

            else
                npad = size(obj.data,1)-size(obj_temp.data,1);
                pad_left = floor(npad / 2);
                pad_right = ceil(npad / 2);
                obj_temp.data = padarray(obj_temp.data, [pad_left pad_left],min(abs(obj_temp.data(:))), 'pre');
                obj_temp.data = padarray(obj_temp.data, [pad_right pad_right],min(abs(obj_temp.data(:))), 'post');
                obj_temp.width = obj.width*scale;
                obj_temp.roi = obj_temp.N*obj_temp.dx;
                % Smooth padding
%                 sigma = obj.roi;amplitude = max(abs(obj_temp.data(:)))*2;
%                 gauss_mask = amplitude.*exp(- (obj_temp.X.^2+obj_temp.Y.^2) /(sigma^2));
%                 [rows, cols] = size(obj_temp.data);
%                 center_mask = ones(rows,cols);
%                 overlap = 0;
%                 center_mask(pad_left+overlap+1:rows-pad_right-overlap, pad_left+overlap+1:cols-pad_right-overlap) = 0;
%                 outer_mask = gauss_mask.*center_mask; 
%                 
                %obj_temp.data = obj_temp.data + outer_mask.*obj_temp.data;
            
            end

            obj = obj_temp;
            if gpuFlag
                obj_temp.data = gpuArray(obj_temp.data); %Put array back to gpu if it was taken off
            else
                obj_temp.data = obj_temp.data;
            end


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

        function obj = conjugate(obj)
            obj.data = conj(obj.data);
        end

        function showme(obj, name)
            % Plots the intensity and phase of the field at the current position.

            % Inputs:
            % obj - Field object to be plotted.
            % name - Title of the figure window (default: 'Field at z = obj.position').

            % Outputs:
            % None.
            arguments
                obj
                name = sprintf('Field at z = %f',obj.position);
            end
            % Create figure window with specified name.
            figure("Position",[700,200,700,400], "Name",name); 
            if obj.roi/2 < obj.width
                limits = [-obj.width obj.width -obj.width obj.width].*1e3/2;
            else
                limits = [-obj.width obj.width -obj.width obj.width].*1e3;
            end
            % Plot intensity and phase in two subplots.
            % Add axis labels and zoom into field, defined by limits
            if strcmp(obj.domain, 'spatial')
                subplot(1,2,1); imagesc(obj.X(1,:)*1e3, obj.Y(:,1)*1e3, abs(obj.data(:,:,1)).^2); title('Intensity');axis image; 
                xlabel('[mm]'); ylabel('[mm]'); axis(limits); colorbar; %clim([0,1])
                subplot(1,2,2); imagesc(obj.X(1,:)*1e3, obj.Y(:,1)*1e3, angle(obj.data(:,:,1))); title('Phase'); axis image
                xlabel('[mm]'); ylabel('[mm]'); axis(limits); colorbar
                
                if ~isempty(obj.phase) % plot unwrapped phase if it is there
                    subplot(1,2,2); imagesc(obj.X(1,:)*1e3, obj.Y(:,1)*1e3, obj.phase); title('Phase'); axis image; colorbar
                end

            elseif strcmp(obj.domain, 'frequency')
                subplot(1,2,1); imagesc(obj.X(1,:), obj.Y(:,1), log10(abs(obj.data(:,:,1))+1e-10)); title('Intensity');axis image;
                xlabel('[Hz]'); ylabel('[Hz]'); colorbar
                subplot(1,2,2); imagesc(obj.X(1,:), obj.Y(:,1), angle(obj.data(:,:,1))); title('Phase'); axis image
                xlabel('[Hz]'); ylabel('[Hz]'); colorbar
            end

            % Add x-axis label.

            % Add figure title.
            sgtitle(name)

            % Set font size of plot to 18 pixels.
            fontsize(gcf,18, 'pixels'); fontname(gcf,'Times New Roman')
            pause(0.05) %Put puase to display image before continuing
        end

        function showintensity(obj, name)
            % Plots the intensity and phase of the field at the current position.

            % Inputs:
            % obj - Field object to be plotted.
            % name - Title of the figure window (default: 'Field at z = obj.position').

            % Outputs:
            % None.
            arguments
                obj
                name = sprintf('Field at z = %f',obj.position);
            end
            % Create figure window with specified name.
            figure("Position",[700,200,700,400], "Name",name); 
            limits = [-obj.width obj.width -obj.width obj.width].*1e3/2;
            % Plot intensity and phase in two subplots.
            % Add axis labels and zoom into field, defined by limits
         
            imagesc(obj.X(1,:)*1e3, obj.Y(:,1)*1e3, abs(obj.data(:,:,1)).^2);axis image; 
            xlabel('[mm]'); ylabel('[mm]'); axis(limits); colorbar; clim([0,1])

            % Add figure title.
            sgtitle(name)

            % Set font size of plot to 18 pixels.
            fontsize(gcf,18, 'pixels'); fontname(gcf,"Times New Roman")
            pause(0.05) %Put puase to display image before continuing
        end


        function obj = fft2(obj)
            % Spatial domain to spatial frequency domain

            obj.data = jfft2(obj.data);

            [nx, ny] = meshgrid((-obj.N/2 : 1 : obj.N/2 - 1));
            deltaf = 1 / (obj.N*obj.dx);
            fX = nx * deltaf; fX(fX==inf) = 1; %set inf to 1
            fY = ny * deltaf; fY(fY==inf) = 1; %set inf to 1

            obj.X = fX;  %convert coordinates to freq domain
            obj.Y = fY;
            obj.dx = deltaf;

            obj.domain = 'frequency';

        end
        function obj = ifft2(obj)
            % Frequency domain to spatial domain
            obj.data = jifft2(obj.data);
            [nx, ny] = meshgrid((-obj.N/2 : 1 : obj.N/2 - 1));
            deltaf = 1 / (obj.N*obj.dx);
            fX = nx * deltaf; fX(fX==inf) = 1; %set inf to 1
            fY = ny * deltaf; fY(fY==inf) = 1; %set inf to 1

            obj.X = fX;  %convert coordinates to spatial domain
            obj.Y = fY;
            obj.dx = deltaf;

            obj.domain = 'spatial';
        end

        function obj = calc_power(obj, width)
            arguments
                obj 
                width = obj.width
            end
            width_px = floor(width./obj.dx);
            center = floor(obj.N/2); 
            cropped_area = obj.data(center-width_px:center+width_px,center-width_px:center+width_px);
            avg_I = mean(abs(cropped_area(:)).^2);
            obj.avg_power = avg_I*obj.width;

            peak_I = max(abs(obj.data(:)).^2);
            obj.peak_power = peak_I;%*obj.roi;
        end

        function s = struct(self)
            publicProperties = properties(self);
            s = struct();
            for fi = 1:numel(publicProperties)
                s.(publicProperties{fi}) = self.(publicProperties{fi});
            end
        end
    end



end