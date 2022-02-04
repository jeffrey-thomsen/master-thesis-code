function [enhancedSignal, AlgorithmParameters] = ...
  speechEnhancement(inputSignal, AlgorithmParameters)

    % Gammatone analysis filterbank
    [decomposedSubbandSignals, AlgorithmParameters] = ...
        subbandDecompositionBinaural(inputSignal, AlgorithmParameters);
    
    enhancedSubbandSignals = decomposedSubbandSignals;
    % Periodicity analysis
%     [sigmaDesired, deltaDesired, snrDesired] = ...
%         periodicityAnalysis(decomposedSubbandSignals, AlgorithmParameters);
    
    % DOA estimation
%     azimuth = directionOfArrivalEstimation(decomposedSubbandSignals);
    
    % signal enhancement
%     enhancedSubbandSignals = ...
%         harmonicCancellation(decomposedSubbandSignals, azimuth, ...
%             sigmaDesired, deltaDesired, snrDesired);

    % Gammatone synthesis filterbank
    [enhancedSignal, AlgorithmParameters] = ...
        subbandResynthesisBinaural(enhancedSubbandSignals, ...
        AlgorithmParameters);
end