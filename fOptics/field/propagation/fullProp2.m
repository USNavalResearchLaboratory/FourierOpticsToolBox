function [x2, y2, Uout, varargout] = fullProp2(efield, direction,Z,steps, turb, dx2,figflag)
%FULLPROP2 uses multi-plane angular spectrum propagation
%
% The direction of propagation is given by the string direction. If the direction
% is 'rev' then the phase screens are reversed along the third
% dimension. If the direction is 'revMirror' then the phase is flipped, in
% addition to a flip of the phase screens. The phase screens, turbulence,
% need to be complex, i.e., exp(1i * phase).
%
% Syntax:
% [x2, y2, Uout] = fullProp2(efield, direction, Z, steps, turbulence, dx2, figflag)
%
% Inputs:
% efield - electric field to be propagated (an EField object)
% direction - direction of propagation: 'fwd', 'rev', or 'revMirror'
% Z - propagation distance (in meters)
% steps - number of propagation steps to use (default = 2)
% turbulence - complex turbulence screens or empty array for vacuum (default = [])
% dx2 - the spatial resolution of the propagated field (default = efield.dx)
% figflag - a flag indicating whether to plot intermediate steps (default = 0)
%
% Outputs:
% x2 - x coordinates of the propagated field
% y2 - y coordinates of the propagated field
% Uout - propagated field
%

arguments
efield
direction
Z
steps = 2
turb = []
dx2 = efield.dx
figflag = 0
end

input = efield.data; 
N = efield.N; 

if isempty(turb)
   phzScreens = ones(N, N, steps);  
   disp('Assuming no turbulence')
else
    phzScreens = turb.data;
    steps = turb.num_screens; 
end

%% Setup propagation direction 

% use vacuum if no phase screen given


% Direction of propagation, call in opposite order if rev
switch direction
    case 'fwd' 
        if Z < 0 
            error('use rev input to propagate backwards instead of negative z')
        end      
    case 'rev'
        Z = -1*Z; 
        phzScreens = conj(flip(phzScreens, 3)); 
        input = (input); 
    case 'revMirror'
        Z = -1*Z; 
        phzScreens = (flip(phzScreens, 3));
        % filp the phase, as would happen in a mirror reflection
        input = conj(input); 
        %input = abs(input).*exp(-1i*angle(input));

end

%% Propagate
Z_array = linspace(0,Z,steps);
[x2, y2, Uout, Uout_all] = ang_spec_multi_prop(input, efield.wvl, efield.dx, dx2, Z_array, phzScreens, figflag);
varargout{1} = Uout_all; 

% %% Old fullprop method. Delete when new method is solid
% % Window function - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% %used in full propagatin with angular spectrum method
% %super-Gaussian absorbing boundary adapted from Jason Schmidt's book
% [nx, ny] = meshgrid((-N/2 : 1 : N/2 - 1));
% nsq = nx.^2 + ny.^2;
% w = 0.47*N;
% boundary = exp(-nsq.^8/w^16); clear('nsq', 'w');
% 
% %% Propagate old
% % do partial propagations using Angular-Spectrum method with absorbing
% % boundary layers
% deltaZ = Z/steps; 
% for ii = 1:turb.num_screens
%     
%     switch direction
%         
%         case 'fwd'
%             % apply phase screen, phase screen is complex
%             input = input.*phzScreens(:,:,ii);
%             % propagate
%             %[input, ~] = angularSpectrumPropagation(input,dx2,dx2,deltaZ, efield.wvl);
%             [x2, y2, input] = ang_spec_prop(input, efield.wvl, efield.dx, dx2, deltaZ); 
%             if ii < turb.num_screens
%                 % apply window
%                 input = input.*boundary;
%             end
%             
%         case {'rev','revMirror'}
%             % propagate
%             %[input, ~] = angularSpectrumPropagation(input,P.R.px,P.R.px,deltaZ,P.wvl, circ(length(input)/1,size(input)));
%             [x2, y2, input] = ang_spec_prop(input, efield.wvl, efield.dx, dx2, deltaZ); 
%             % apply phase screen, phase screen is complex
%             input = input.*(phzScreens(:,:,ii));
%               
%             if ii < turb.num_screens
%                 % apply window
%                 input = input.*boundary;
%             end
%     end
%     if figflag
%         %TODO: Put these into one figure? 
%         figure(200); 
%         subplot(1,turb.num_screens,ii); imagesc(angle(input)); axis image off; sgtitle('Field Phase')
%         figure(201); 
%         subplot(1,turb.num_screens,ii); imagesc(abs(input)); axis image off; sgtitle('Field Amplitude')
%         figure(202); 
%         subplot(1,turb.num_screens,ii); imagesc(angle(phzScreens(:,:,ii))); axis image off; sgtitle('Phase screen')
%     end
% end
% 
% 
% Uout = input; 


end

