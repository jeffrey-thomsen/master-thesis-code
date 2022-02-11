function [azimuth, azDetectedIndexVectors, AlgorithmStates] = ...
    directionOfArrivalEstimation(subbandSignalArray, AlgorithmParameters, AlgorithmStates)

centerFreqs = AlgorithmStates.L.GammatoneStates.analyzer.center_frequencies_hz;
azimuth = subbandSignalArray;
end