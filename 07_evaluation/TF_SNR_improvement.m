function [deltaSnr] = TF_SNR_improvement(inputTargetSignal, inputInterfSignal, ...
  outputTargetSignal_H, outputInterfSignal_H, AlgorithmStates, samplingRateHz, centerFreqsHz)

    nBands = 30;
    
    % compute IBM and extract glimpse masks from simulation data
    [inputTargetSubbandSignals, ~] = subbandDecompositionBinaural(...
        inputTargetSignal, AlgorithmStates);
    [inputInterfSubbandSignals, ~] = subbandDecompositionBinaural(...
        inputInterfSignal, AlgorithmStates);
    
    inputSnr.L = computeTimeFreqTir(inputTargetSubbandSignals.L, ...
        inputInterfSubbandSignals.L, samplingRateHz, nBands);
    inputSnr.R = computeTimeFreqTir(inputTargetSubbandSignals.R, ...
    inputInterfSubbandSignals.R, samplingRateHz, nBands);        
    
    % compute TF gray mask SNR improvement
    [outputTargetSubbandSignals_H, ~] = subbandDecompositionBinaural(...
        outputTargetSignal_H, AlgorithmStates);
    [outputInterfSubbandSignals_H, ~] = subbandDecompositionBinaural(...
        outputInterfSignal_H, AlgorithmStates);
    
    outputSnr.L = computeTimeFreqTir(outputTargetSubbandSignals_H.L, ...
        outputInterfSubbandSignals_H.L, samplingRateHz, nBands);
    outputSnr.R = computeTimeFreqTir(outputTargetSubbandSignals_H.R, ...
        outputInterfSubbandSignals_H.R, samplingRateHz, nBands);
    
    deltaSnr.L = outputSnr.L - inputSnr.L;
    deltaSnr.R = outputSnr.R - inputSnr.R;

    plotTimeFreqBlockSnrHeatmap(deltaSnr.L, samplingRateHz, centerFreqsHz)
