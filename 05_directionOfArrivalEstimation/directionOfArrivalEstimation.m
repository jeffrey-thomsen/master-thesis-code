function [azimuthRadCells, azDetectedIndexVectorCells, AlgorithmStates] = ...
    directionOfArrivalEstimation(subbandSignalArray, AlgorithmParameters, ...
    AlgorithmStates)

samplingRateHz = AlgorithmParameters.Gammatone.samplingRateHz;
centerFreqsHz = AlgorithmStates.L.GammatoneStates.analyzer.center_frequencies_hz;
cycleDurationSeconds = 1./centerFreqsHz;
tau = AlgorithmParameters.nCyclesTau.*cycleDurationSeconds;
% a = exp(-1./(samplingRateHz.*tau));


nBands = AlgorithmParameters.Gammatone.nBands;
azimuthRadCells = cell(1, nBands);
azDetectedIndexVectorCells = cell(1, nBands);
for iBand = 1:nBands
    subbandSignal.L = subbandSignalArray.L(:,iBand);
    subbandSignal.R = subbandSignalArray.R(:,iBand);
    States.L = AlgorithmStates.L.ProcessingStates{iBand};
    States.R = AlgorithmStates.R.ProcessingStates{iBand};
    States.Binaural = AlgorithmStates.Binaural.ProcessingStates{iBand};
    tauS = tau(iBand);
    
    %% ILD
    [subbandSignalLp, States] = ...
        ildBinauralLowPassFilter(States, subbandSignal, samplingRateHz);
    
    % % interaural level difference, eq. 5 in Dietz (2011)
    ildDb = 20*log10(subbandSignalLp.L ./ subbandSignalLp.R);
    % % max(sig,1e-4) avoids division by zero
    %     20.*log10(max(subbandSignalLp.L,1e-4) ./ max(subbandSignalLp.L,1e-4));
    
    %% ITF
    itf = subbandSignal.L .* conj(subbandSignal.R);
    
    %% IVS
    [ivsMask, States.Binaural] = calcIvs(itf, tauS, samplingRateHz, ...
        AlgorithmParameters.ivsThreshold, States.Binaural);
    
    %% IPD
    ipdRad = angle(itf);
    % test for equality: !!!
    % ipd_lp = angle(lowpass(outp.itf,a)); % Dietz2011 lowpass function
    [ipdLpRad, States.Binaural.ipdLP] = ...
        firstOrderLowPass(ipdRad, States.Binaural.ipdLP, AlgorithmParameters, tauS);
    
    ipdLpRad = disambiguateIpd(ipdLpRad, ildDb);
    
    %% azimuth - estimated angle of arrival
    azimuthRad = ipdLpRad;
    azDetectedIndexVector = ivsMask; % placeholder while developing

%     azimuthRad = ipdToAzimuthMapping(ipdLpRad, lookUpTable);
%     azimuthRad = coherenceFilterMask(azimuthRad, ivsMask);

    % write states back into global structs
    AlgorithmStates.L.ProcessingStates{iBand} = States.L;
    AlgorithmStates.R.ProcessingStates{iBand} = States.R;
    AlgorithmStates.Binaural.ProcessingStates{iBand} = States.Binaural;
    
    % write computed values into structs for enhancement stage
    azDetectedIndexVectorCells{iBand} = azDetectedIndexVector;
    azimuthRadCells{iBand} = azimuthRad;
end
end