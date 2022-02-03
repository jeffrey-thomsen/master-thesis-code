function enhancedSignal = speechEnhancement(inputSignal)
    % Gammatone analysis filterbank
    deconposedSubbandSignals = subbandDecomposition(inputSignal);
    
    % Periodicity analysis
    [sigmaDesired, deltaDesired, snrDesired] = ...
        periodicityAnalysis(deconposedSubbandSignals);
    
    % DOA estimation
    azimuth = directionOfArrivalEstimation(deconposedSubbandSignals);
    
    % signal enhancement
    enhancedSubbandSignals = ...
        harmonicCancellation(deconposedSubbandSignals, azimuth, ...
            sigmaDesired, deltaDesired, snrDesired);

    % Gammatone synthesis filterbank
    enhancedSignal = subbandResynthesis(enhancedSubbandSignals);
end