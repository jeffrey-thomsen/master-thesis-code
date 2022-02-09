% run a binaural signal, signal block or signal sample through the
% Thomsen2022 speech enhancement algorithm
% inputSignal - real-valued Nx2 matrix which must satisfy the conditions 
% that can be checked by running it through testInputSignal.m
% AlgorithmParameters - struct containing all information necessary for
% running the simulation including filter states
% enhancedSignal - real-valued Nx2 matrix containing the processed binaural
% signal
function [enhancedSignal, AlgorithmParameters] = ...
  speechEnhancement(inputSignal, AlgorithmParameters)

    % Gammatone analysis filterbank - decompose signal into frequency bands
    [subbandSignalArray, AlgorithmParameters] = ...
        subbandDecompositionBinaural(inputSignal, AlgorithmParameters);
    
    if AlgorithmParameters.Cancellation || AlgorithmParameters.Enhancement
        % Periodicity analysis - detect subband samples with periodic components
        [sigmaDesired, deltaDesired, snrDesired, p0DetectedIndexVectors, ...
            AlgorithmParameters] = periodicityAnalysis(subbandSignalArray, ...
            AlgorithmParameters);
        
        % DOA estimation - estimate angle of incidence for coherent sources
        [azimuth, azDetectedIndexVectors, AlgorithmParameters] = ...
            directionOfArrivalEstimation(subbandSignalArray, ...
            AlgorithmParameters);
        
        % signal enhancement - apply harmonic cancellation to unwanted periodic
        % components
        enhancedSubbandSignalArray = ...
            harmonicCancellation(subbandSignalArray, azimuth, ...
            azDetectedIndexVectors, sigmaDesired, deltaDesired, snrDesired, ...
            p0DetectedIndexVectors, AlgorithmParameters);
    else
        enhancedSubbandSignalArray = subbandSignalArray;
    end

    % Gammatone synthesis filterbank - resynthesize enhanced subband signals
    [enhancedSignal, AlgorithmParameters] = ...
        subbandResynthesisBinaural(enhancedSubbandSignalArray, ...
        AlgorithmParameters);
end