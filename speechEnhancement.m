% run a binaural signal, signal block or signal sample through the
% Thomsen2022 speech enhancement algorithm
% inputSignal - real-valued Nx2 matrix which must satisfy the conditions 
% that can be checked by running it through testInputSignal.m
% AlgorithmParameters - struct containing all information necessary for
% running the simulation including filter states
% enhancedSignal - real-valued Nx2 matrix containing the processed binaural
% signal
function [enhancedSignal, AlgorithmStates] = ...
  speechEnhancement(inputSignal, AlgorithmParameters, AlgorithmStates)

    % Gammatone analysis filterbank - decompose signal into frequency bands
    [subbandSignalArray, AlgorithmStates] = ...
        subbandDecompositionBinaural(inputSignal, AlgorithmStates);
    
    if AlgorithmParameters.Cancellation || AlgorithmParameters.Enhancement
        % Periodicity analysis - detect subband samples with periodic components
        [sigmaDesired, deltaDesired, snrDesired, p0DetectedIndexVectors, ...
            AlgorithmStates] = periodicityAnalysis(subbandSignalArray, ...
            AlgorithmParameters, AlgorithmStates);
        
        % DOA estimation - estimate angle of incidence for coherent sources
        [azimuth, azDetectedIndexVectors, AlgorithmStates] = ...
            directionOfArrivalEstimation(subbandSignalArray, ...
            AlgorithmParameters, AlgorithmStates);
        
        % signal enhancement - apply harmonic cancellation to unwanted periodic
        % components
        enhancedSubbandSignalArray = ...
            harmonicCancellation(subbandSignalArray, azimuth, ...
            azDetectedIndexVectors, sigmaDesired, deltaDesired, snrDesired, ...
            p0DetectedIndexVectors, AlgorithmParameters, AlgorithmStates);
    else % just for testing the gammatone filterbank without processing
        enhancedSubbandSignalArray = subbandSignalArray;
    end

    % Gammatone synthesis filterbank - resynthesize enhanced subband signals
    [enhancedSignal, AlgorithmStates] = ...
        subbandResynthesisBinaural(enhancedSubbandSignalArray, ...
        AlgorithmStates);
end