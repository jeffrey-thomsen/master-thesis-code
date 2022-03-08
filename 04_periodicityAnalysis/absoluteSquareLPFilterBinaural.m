% Compute the absolute square values of and lowpass filter sigma and delta 
% signals as a pre-processing step for the periodicity analysis
%
% Input:
% sigma, delta - structs containing arrays of Sigma and Delta signals for a
% range of period candidate values of the left and right channel
% AlgorithmParameters - struct containing the LP filter time contant 
% FilterStates - struct containing the filter states
%
% Output:
% sigma, delta - strucks of processed Sigma and Delta signals
% FilterStates - see above
function [sigma, delta, FilterStates] = ...
    absoluteSquareLPFilterBinaural(sigma, delta, AlgorithmParameters, FilterStates)
        
        % compute absolute squared values
        sigma.L = abs(sigma.L).^2;
        sigma.R = abs(sigma.R).^2;
        delta.L = abs(delta.L).^2;
        delta.R = abs(delta.R).^2;
    
        % lowpass filter values for temporal smoothing
        [sigma.L, FilterStates.L.sigmaNormLP] ...
            = firstOrderLowPass(sigma.L, ...
            FilterStates.L.sigmaNormLP, AlgorithmParameters);
        [sigma.R, FilterStates.R.sigmaNormLP] ...
            = firstOrderLowPass(sigma.R, ...
            FilterStates.R.sigmaNormLP, AlgorithmParameters);
        [delta.L, FilterStates.L.deltaNormLP] ...
            = firstOrderLowPass(delta.L, ...
            FilterStates.L.deltaNormLP, AlgorithmParameters);
        [delta.R, FilterStates.R.deltaNormLP] ...
            = firstOrderLowPass(delta.R, ...
            FilterStates.R.deltaNormLP, AlgorithmParameters);
end