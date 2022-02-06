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
    [decomposedSubbandSignals, AlgorithmParameters] = ...
        subbandDecompositionBinaural(inputSignal, AlgorithmParameters);
    
    enhancedSubbandSignals = decomposedSubbandSignals;
    % Periodicity analysis - detect subband samples with periodic components
%     [sigmaDesired, deltaDesired, snrDesired] = ...
%         periodicityAnalysis(decomposedSubbandSignals, AlgorithmParameters);
    
    % DOA estimation - estimate angle of incidence for coherent sources
%     azimuth = directionOfArrivalEstimation(decomposedSubbandSignals);
    
    % signal enhancement - apply harmonic cancellation to unwanted periodic
    % components
%     enhancedSubbandSignals = ...
%         harmonicCancellation(decomposedSubbandSignals, azimuth, ...
%             sigmaDesired, deltaDesired, snrDesired);

    % Gammatone synthesis filterbank - resynthesize enhanced subband signals
    [enhancedSignal, AlgorithmParameters] = ...
        subbandResynthesisBinaural(enhancedSubbandSignals, ...
        AlgorithmParameters);
end