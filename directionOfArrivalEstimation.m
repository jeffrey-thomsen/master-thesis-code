function [azimuth, azDetectedIndexVectors, AlgorithmParameters] = ...
    directionOfArrivalEstimation(subbandSignalArray, AlgorithmParameters)

centerFreqs = AlgorithmParameters.L.FilterStates.Gammatone.analyzer.center_frequencies_hz;
azimuth = subbandSignalArray;
end