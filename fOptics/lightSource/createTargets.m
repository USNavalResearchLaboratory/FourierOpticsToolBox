function [ targets ] = createTargets(P, targetType)
%createTargets - create a target: three bars, USAF pattern, or homemade
%
% There are two options: 'threeBar' or 'USAF'. The three bars take up the
% area defined by P.T.W(width) and P.T.H (height). The bars are as wide as
% P.T.W. The height is divided into 5 segments, 3 of which are bars, with
% the other two blank. The USAF option load an image I saved from the web.
% The image is resized to P.T.H x P.T.W. The random speckle generation is
% done following Ref [1]. Reqquires DG codes.
%
% Syntax:
%           [ targets ] = createTargets(P, 'threeBar');
%           [ targets ] = createTargets(P, 'USAF');
%           [ targets ] = createTargets(P, 'ptSrc_gauss_NL_beacon')
%
% Inputs:
%   P - parameter structure
%   option - a string with target flag
%
% Outputs:
%    targets - an array of targets, each slice is a different realization
%
% Example 1: fixed three bar target
%    [ P ] = parameterSetup();
%    [ targets ] = createTargets(P, 'threebar')
%
% Example 2: Gaussian point source target
%    [ P ] = parameterSetup();
%    [ targets ] = createTargets(P, 'ptSrc_gauss_NL_beacon')
%
% References:
% [1] "Phase-error correction in digital hologaphy" by Samuel T Thurman
%     and James R Fienup, J. Opt. Soc Am. A V25 N4 April 2008 983 - 994
%
% Other m-files required: genGrids
% Subfunctions: smoothAndRandSpeckReals
% MAT-files required: none
% Other files needed: fig-9-imtf.jpg
% taken from https://www.edmundoptics.com/contentassets/c1499ee3006340bab8c768b7e799cdea/fig-9-imtf.jpg
%
% See also: parameterSetup, genGrids

% Author: Dennis F. Gardner Jr.
%         Code 5661
%         4555 Overlook Ave SW
%         Washington, DC 20375
% Email: dennis.gardner@nrl.navy.mil
% Website: https://www.linkedin.com/in/dennisfgardner/
% March 2018; Last revision: 20-10-2020


% initialize array of targets
targets = zeros(P.NmultiStep, P.NmultiStep, P.T.specks);
if P.S.NLprop
    NLstring = 'NL';
else
    NLstring = 'Lin';
end

switch targetType
    
    case 'ptSrc'
        % load vacuum spot of linear focus
        fileNameFvac = fVacFname(P.S.ApodSigmaFrac,P.S.NLfrac,'Lin', P.S.D/P.U.cm);
        load(['./focus/vac/',fileNameFvac,'.mat'],'focusVac');

        % clean up (abs and threshold and sharpen)
        focusVac = abs(focusVac).^7;
        focusVac(focusVac <= (max(focusVac(:)) / 1000)) = 0;
        
        % smooth the edges and generate random speckels
        targets = smoothAndRandSpeckReals(P, targets, focusVac);
        
    case 'beacon'

        %load turbulated focused spot
        fileNameFturb = fTurbFname(P.S.ApodSigmaFrac,NLstring,P.S.NLfrac,P.Cn2str, P.S.D/P.U.cm, P.jj);
        load(['./focus/turb/',fileNameFturb,'.mat'],'focusTurb');
        
        % clean up (abs and threshold)
        focusTurb = abs(focusTurb);
        focusTurb(focusTurb <= (max(focusTurb(:)) / 1000)) = 0;
        
        targets = smoothAndRandSpeckReals(P, targets, focusTurb);
        
        
    case 'Vbars'
        % load three bar target
        %         amp = imread('./bin/fig-9-imtf.jpg');
        amp = imread('./bin/Vbars-10.png');
        amp = rgb2gray(amp);
        % cahnge into a binary image
        amp(amp<=200) = 0;
        amp(amp>0) = 1;
        amp = double(amp);
        
        % height and width in units of pixels
        Hpx = round(P.T.H/P.R.px);
        Wpx = round(P.T.W/P.R.px);
        amp = imresize(amp, [Hpx Wpx]);
        
        % embed in full array
        amp = embedArray(amp, P.NmultiStep, P.NmultiStep, 0);
        
        % smooth the edges and generate random speckels
        targets = smoothAndRandSpeckReals(P, targets, amp);
        
        
    case 'threeBar'
        % height and width 'radius' in units of pixels
        HpxR = round(P.T.H/5/P.R.px/2);
        WpxR = round(P.T.W/P.R.px/2);
        % separation in pixels
        Spx = 2*HpxR;
        
        % create each bar one-at-a-time
        amp = zeros(P.NmultiStep);
        for kk = -2:2:2
            temp = zeros(P.NmultiStep);
            temp(P.NmultiStep/2 - HpxR + kk*Spx : P.NmultiStep/2 + HpxR - 1 +kk*Spx,...
                P.NmultiStep/2 - WpxR : P.NmultiStep/2 + WpxR - 1) = 1.0;
            amp = amp + temp;
        end
        
        % smooth the edges and generate random speckels
        targets = smoothAndRandSpeckReals(P, targets, amp);
        
        
    case 'USAF'
        % load the USAF test pattern
        %         amp = imread('./bin/fig-9-imtf.jpg');
        amp = imread('USAF_Resolution_Chart_A1-780.jpg');
        amp = rgb2gray(amp);
        % cahnge into a binary image
        amp(amp<=200) = 0;
        amp(amp>0) = 1;
        amp = double(amp);
        
        % height and width in units of pixels
        Hpx = round(P.T.H/P.R.px);
        Wpx = round(P.T.W/P.R.px);
        amp = imresize(amp, [Hpx Wpx]);
        
        % embed in full array
        amp = embedArray(amp, P.NmultiStep, P.NmultiStep, 0);
        
        % smooth the edges and generate random speckels
        targets = smoothAndRandSpeckReals(P, targets, amp);
        
    case 'homemade'
        % load the home made target
        amp = imread('./bin/Target-04.png');
        amp = rgb2gray(amp);
        amp = double(amp);
        % normalize 0 to one
        amp = (amp-min(amp(:)))/max(amp(:)-min(amp(:)));
        
        % height and width in units of pixels
        Hpx = round(P.T.H/P.R.px);
        Wpx = round(P.T.W/P.R.px);
        amp = imresize(amp, [Hpx Wpx], 'nearest');
        
        % embed in full array
        amp = embedArray(amp, P.NmultiStep, P.NmultiStep, 0);
        
        % smooth the edges and generate random speckels
        targets = smoothAndRandSpeckReals(P, targets, amp);
end

end

