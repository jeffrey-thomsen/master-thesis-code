function [deltaSnr] = TF_SNR_improvement(inputTargetSignal, inputInterfSignal, ...
  outputTargetSignal_H, outputInterfSignal_H, AlgorithmStates, samplingRateHz, centerFreqsHz)

    nBands = 30;
    
    % compute IBM and extract glimpse masks from simulation data
    [inputTargetSubbandSignals, ~] = subbandDecompositionBinaural(...
        inputTargetSignal, AlgorithmStates);
    [inputInterfSubbandSignals, ~] = subbandDecompositionBinaural(...
        inputInterfSignal, AlgorithmStates);
    
    inputSnr.L = computeTimeFreqSnr(inputTargetSubbandSignals.L, ...
        inputInterfSubbandSignals.L, samplingRateHz, nBands);
    inputSnr.R = computeTimeFreqSnr(inputTargetSubbandSignals.R, ...
    inputInterfSubbandSignals.R, samplingRateHz, nBands);        
    
    % compute TF gray mask SNR improvement
    [outputTargetSubbandSignals_H, ~] = subbandDecompositionBinaural(...
        outputTargetSignal_H, AlgorithmStates);
    [outputInterfSubbandSignals_H, ~] = subbandDecompositionBinaural(...
        outputInterfSignal_H, AlgorithmStates);
    
    outputSnr.L = computeTimeFreqSnr(outputTargetSubbandSignals_H.L, ...
        outputInterfSubbandSignals_H.L, samplingRateHz, nBands);
    outputSnr.R = computeTimeFreqSnr(outputTargetSubbandSignals_H.R, ...
        outputInterfSubbandSignals_H.R, samplingRateHz, nBands);
    
    deltaSnr.L = outputSnr.L - inputSnr.L;
    deltaSnr.R = outputSnr.R - inputSnr.R;

    plotTimeFreqBlockSnrHeatmap(deltaSnr.L, samplingRateHz, centerFreqsHz)
