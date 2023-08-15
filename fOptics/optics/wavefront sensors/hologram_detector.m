classdef hologram_detector < Efield
    %hologram_detector Class for different types of holographic detectors.
    %   This class extends the Efield class and contains methods for
    %   on-axis holography, off-axis holography, and phase recovery.
    %
    %   Properties:
    %       type - Type of holographic detector (e.g. 'on-axis', 'off-axis', 'phase recovery')
    %       Uref - Reference field
    %       Usignal - Signal field
    %
    %   Methods:
    %       hologram_detector - Constructor method to initialize the class
    %
    %   Example usage:
    %       % Create an on-axis holographic detector object
    %       hd = hologram_detector('phase_shift_holography', field_struct);
    %
    %       % Call the phase_shift_holography holography method
    %       hd.phase_shift_holography();
    %
    %   See also Efield.
    properties
        type
        Uref
        Usignal
    end

    methods
        function obj = hologram_detector(type, field_struct)
            %hologram_detector Constructor method for the class.
            %   obj = hologram_detector(type, field_struct) creates a
            %   holographic detector object of the specified type with the
            %   input field data.
            %
            %   Inputs:
            %       type - Type of holographic detector (e.g. 'on-axis', 'off-axis', 'phase recovery')
            %       field_struct - Struct containing field data
            arguments          
                type 
                field_struct = get_default_Efield_struct();
            end
            
            obj = obj@Efield(field_struct);
            
            %Initialize class. Describe type of detector to use
            obj.type = type;
            obj.data = field_struct.data; 
        end

        function obj = phase_shift_holography(obj, Uref, Usignal, num_angles)
            %Phase Shift Holography
            %   this function takes an on-axis reference and signal fields
            %   and applied the phase shifting method, based on number of
            %   angles given (num_angles) to generate an Intensity that is
            %   detected. 
            arguments
                obj
                Uref
                Usignal
                num_angles = 4 %Assume 4 phase shifts default
            end
            

            obj.Uref = Uref; Ur = Uref.data; %Intermediate variables
            obj.Usignal = Usignal; Us = Usignal.data; 

            Ar = 1; As = 1; %Any weighting between Uref and Usignal
            Us = As.*Us; Ur = Ar.*Ur;             
            phi_increment = 2*pi / num_angles;
            phi = 0; 
            Intensity = zeros([size(Us),num_angles]); 
            for i = 1:num_angles
                Intensity(:,:,i) = abs(Ur + Us.*exp(1j*phi)).^2; 
                phi = phi + phi_increment;
            end
            
            obj.data = Intensity; 
        end

       
        function obj = finch_hologram(obj)
            % finch_hologram takes Intensity_array, which is a 3D array of
            % measured intensities assumed to have taken using
            % phase_shift_holography process. The complex hologram is
            % recovered. This function is designed to be used in experiment
            % and modeling. 
            % Assume Intensity_array is measured with phase shift 0, pi/2,
            % pi, 3pi/2 in that order. 
            % Complex hologram given by H = 1/4*[(I_0-I_pi)-i(I_pi/2-I_3pi/2)
            Intensity_array = obj.data; 
            complex_hologram = 0.25.*(Intensity_array(:,:,1)-Intensity_array(:,:,3)) + ...
                1i.*(Intensity_array(:,:,2)-Intensity_array(:,:,4)); 

            obj.data = complex_hologram; 
            %Old calculation
            %complex_hologram = 0.25*((((abs(A1*U_ref.data + A2*U_signal.data*exp(1j*0))).^2)-((abs(A1*U_ref.data + A2*U_signal.data*exp(1j*pi))).^2)) + 1j*(((abs(A1*U_ref.data + A2*U_signal.data*exp(1j*(pi/2)))).^2)-((abs(A1*U_ref.data + A2*U_signal.data*exp(1j*(3*pi/2)))).^2))); 
        end
        function obj = off_axis_dh(obj, intensity)
                %This function assumes off axis holography. It takes the
                %Fourier Transform of measured intensity hologram, displays FT for user to choose window,
                %and inverse FT to retrieve phase 

                %TODO: Take FT of Intensity. Crop area. IFFT
                I_frequency = jfft2(intensity.data); 
                figure; imagesc(log10(abs(I_frequency))); title('Draw Rectangle around cropping area'); axis off image
                d = drawrectangle; pos = round(d.Position); %pos = [xi, yi, width, height]
                cropped = I_frequency(pos(1):pos(1)+pos(3), pos(2):pos(2)+pos(4)); 
                padded = padarray(cropped, [1024,1024]); 
                field = jifft2(padded);
                obj.data = field; 
        end
        function obj = off_axis_dh_analysis_auto(obj, P)
            % function of FFT method of extracting field. 
            % Input:
            % obj - previously initialized hologram detector obj
            % P - Parameters for off axis setup. Mainly need angle of
            % incident Ref beam
            % Output:
            %   object with recovered field
            obj_freq = obj.fft2; %convert to freq domain
            %obj_freq.showme('FFT of hologram')

            frequency = tan(deg2rad(P.angle))/P.wvl; %freq in Hz 
            [~,ind] = min(abs(obj_freq.X(1,:)+frequency));
            

            center = ind; %Get the center from P - my off axis angle 
            diameter = obj.N/8; %300; 
            %define cropping area
            pos = [center-diameter center-diameter diameter*2 diameter*2 ]; %pos = [xi, yi, width, height]
            cropped = obj_freq.data(pos(1):pos(1)+pos(3)-1, pos(2):pos(2)+pos(4)-1); 
            % multiply with circ 
            aperture = circ(2*diameter, [2*diameter 2*diameter] ); 
            cropped = cropped.*aperture; 
            % pad array
            padsize = (obj.N - 2*diameter)/2;
            padded = padarray(cropped, [padsize,padsize], 0);
            % output field
            obj_freq.data = padded;
            %obj_freq.showme('Cropped & Padded FFT of hologram')
            obj = obj_freq.ifft2; 
        end
        function obj = lateralShearAnalysis_auto(obj, lateral_shift, direction)
            % function of FFT method of extracting field. 
            % Input:
            % obj - previously initialized hologram detector obj
            % P - Parameters for off axis setup. Mainly need angle of
            % incident Ref beam
            % Output:
            %   object with recovered field
            obj_freq = obj.fft2; %convert to freq domain
            center_dc = round(obj.N/2); %Shear from P.shearPlate.lateral_shift
            center_signal = round((lateral_shift*obj.dx)/obj.wvl); %2932
            diameter = round(1/(obj.dx*lateral_shift*obj_freq.dx));  
            %define cropping area
            if strcmp(direction, 'y')
                pos = [center+diameter*2 center diameter diameter ]; %pos = [xi, yi, width, height]
            elseif strcmp(direction,'x')
                pos = [center  center+diameter*2 diameter diameter ]; %pos = [xi, yi, width, height]
            else
                error('Need a direction of shear')
            end
            
            cropped = obj_freq.data(pos(1)-pos(3):pos(1)+pos(3)-1, pos(2)-pos(4):pos(2)+pos(4)-1); 
            figure; imagesc(log10(abs(cropped)))

            obj_freq.showme('Cropped & Padded FFT of hologram')
            % multiply with circ 
            aperture = circ(diameter, [2*diameter 2*diameter] ); 
            cropped = cropped.*aperture; 
            % pad array
            padsize = (obj.N - 2*diameter)/2;
            padded = padarray(cropped, [padsize,padsize], 0);
            % output field
            obj_freq.data = padded;
            
            obj = obj_freq.ifft2; 
            
        end

        function obj = lateralShearPhaseRecovery(obj)
            
            % Find signal peak in spectral frequency of Ref. Crop out same
            % area in Signal. 

            ref = obj.Uref; signal = obj.Usignal; 

            refFT = ref.fft2; signalFT = signal.fft2; 
            
            mask = ones(ref.N); center = round(ref.N/2); radius = 5; 
            mask(center-radius:center+radius+1,center-radius:center+radius+1)=0;
            refFT.data = refFT.data.*mask; 
            
            tmp = gather(abs(refFT.data));
            [pks,locs_y,locs_x] = peaks2(abs(refFT.data), 'Threshold', 10); 
            [~,ind] = max(pks); indy = locs_y(ind); indx = locs_x(ind);
            radius = gather(round(abs(indx-indy)/2));

            % crop out areas
            refFT.data = refFT.data(indy-radius:indy+radius-1,indx-radius:indx+radius-1); 
            signalFT.data = signalFT.data(indy-radius:indy+radius-1,indx-radius:indx+radius-1); 
            %signalFT.showme('Cropped spectral data')
            % pad data
            padsize = (ref.N - 2*radius)/2;
            refFT.data = padarray(refFT.data, [padsize,padsize], 0);
            signalFT.data =  padarray(signalFT.data, [padsize,padsize], 0);

            signal = signalFT.ifft2; ref = refFT.ifft2;
            mask = circ(ref.N/2,[ref.N, ref.N]);
            signal.data = signal.data.*mask;
            ref.data = ref.data.*mask;
            %signal = signal.unwrap_phase;
            signal.data = signal.data./ref.data; %Correct for reference phase
            signal.data(isnan(signal.data)) = 0; 
            obj = signal; 
            
            

        end
        
        function obj = radialShearPhaseRecovery(obj)
             ref = obj.Uref; signal = obj.Usignal; 
             signal.data = signal.data./ref.data;
            signal.data(isnan(signal.data)) = 0; 
            %mask = abs(ref.data)>0.005; %Tried using mask to remove noise
            %outside aperture. No longer needed?
            %signal.data = signal.data; 
             obj = signal; 
        end

  
    end
end