function [wavefront] = beamExpander(wavefront)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments
    wavefront = Efield; % Uses default values for example run
end

P = parameters_beamExpander(wavefront); 
%wavefront.showme('before beam expander')
lens1 = thinLens(P.lens1);
wavefront = wavefront.resample(lens1); wavefront = wavefront.apply_optic(lens1);
%wavefront.showme('in BE: after lens1')

lens2 = thinLens(P.lens2);
wavefront = wavefront.propagateTo(lens2);
%wavefront.showme('in BE: after lens2')
end

