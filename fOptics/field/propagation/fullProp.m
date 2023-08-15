function [ output, curLoc ] = fullProp(P, input, dir, phzScreens, curLoc, figflag)
%fullProp - use multi-plane angular spectrum propagation
%
% The direction of propagation is given by the string dir. If the direction
% is 'rev' then the phase screens (phzScreens) are reversed along the third
% dimension. If the direction is 'revMirror' then the phase is filpped, in
% addition to a flip of the phase screens. The phase screens, phzScreens,
% need to be complex, ie exp(1i*phase). 
%
% Syntax: 
% [ output, curLoc ] = fullProp(P, input, 'fwd', [], curLoc); 
% [ output, curLoc ] = fullProp(P, input, 'fwd', phzScreens, curLoc); 
% [ output, curLoc ] = fullProp(P, input, 'rev', phzScreens, curLoc);
% [ output, curLoc ] = fullProp(P, input, 'revMirror', phzScreens, curLoc);
%x
% Inputs:
%   P - parameter structure
%   input - field to be propagated
%   dir - direction of propagation: 'fwd', 'rev', or 'revMirror'
%   phzScreens - complex turbulence screens or empty array for vacuum
%   curLoc - current location of the field
%
% Outputs:
%    output - field at the final plane
%    curLoc - location of field after propagation
% 
% Example: propagate a circle in vac (requires DGcodes)
%    [ P ] = parameterSetup();
%    source = createCirclePx(P.NmultiStep, P.NmultiStep, round(P.R.D/P.px/2));
%    [ sourceAtTarget, ~ ] = fullProp(P, source, 'fwd', [], 0);
%    imageCompare(source, abs(sourceAtTarget), 'At source plane','At target plane');
%
% Other m-files required: angularSpectrumPropagation
% Subfunctions: printLoc
% MAT-files required: none
%
% See also: parameterSetup, angularSpectrumPropagation,
% angSpecPropSampling, quickProp

% Author: Dennis F. Gardner Jr.
%         Code 5661
%         4555 Overlook Ave SW
%         Washington, DC 20375
% Email: dennis.gardner@nrl.navy.mil
% Website: https://www.linkedin.com/in/dennisfgardner/
% March 2018; Last revision: 26-March-2018
arguments
P
input
dir
phzScreens
curLoc = 0
figflag = 0
end
% use vacuum if no phase screen given
if isempty(phzScreens)
   phzScreens = ones(P.NmultiStep, P.NmultiStep, P.Turb.screens);  
end

% Direction of propagation, call in opposite order if rev
switch dir
    case 'fwd' 
        deltaZ = P.Turb.deltaZ;
    case 'rev'
        deltaZ = -1*P.Turb.deltaZ; 
        phzScreens = flip(phzScreens, 3); 
    case 'revMirror'
        deltaZ = -1*P.Turb.deltaZ; 
        phzScreens = flip(phzScreens, 3);
        % filp the phase, as would happen in a mirror reflection
        %inputconj = conj(input); 
        input = abs(input).*exp(-1i*angle(input));

end

% do partial propagations using Angular-Spectrum method with absorbing
% boundary layers
for ii = 1:P.Turb.screens
    
    switch dir
        
        case 'fwd'
            % apply phase screen, phase screen is complex
            input = input.*phzScreens(:,:,ii);
            % propagate
            [input, ~] = angularSpectrumPropagation(input,P.R.px,P.R.px,deltaZ,P.wvl);
            %[x2, y2, input] = ang_spec_prop(input, P.wvl, P.T.dx,P.R.dx, deltaZ); 
            curLoc = printLoc(curLoc, deltaZ);
            if ii < P.Turb.screens
                % apply window
                input = input.*P.win;
            end
            
        case {'rev','revMirror'}
            % propagate
            [input, ~] = angularSpectrumPropagation(input,P.T.dx,P.R.dx,deltaZ,P.wvl);
            % apply phase screen, phase screen is complex
            input = input.*conj(phzScreens(:,:,ii));
              
            curLoc = printLoc(curLoc, deltaZ);
            if ii < P.Turb.screens
                % apply window
                input = input.*P.win;
            end
    end
    if figflag
        %TODO: Put these into one figure? 
        figure(200); 
        subplot(1,P.Turb.screens,ii); imagesc(angle(input)); axis image off; sgtitle('Field Phase')
        figure(201); 
        subplot(1,P.Turb.screens,ii); imagesc(abs(input)); axis image off; sgtitle('Field Amplitude')
        figure(202); 
        subplot(1,P.Turb.screens,ii); imagesc(angle(phzScreens(:,:,ii))); axis image off; sgtitle('Phase screen')
    end
end


output = input; 

% helper function - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
    function curLoc = printLoc(curLoc, deltaZ)
        
        fprintf('propagating delta %.2f from %.2f to %.2f\n',deltaZ, curLoc, curLoc + deltaZ);
        curLoc = curLoc + deltaZ;
        
    end
end

