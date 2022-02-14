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