% compute TIR improvement in time-frequency blocks
function [deltaTir] = TF_TIR_improvement(inputTargetSignal, inputInterfSignal, ...
  outputTargetSignal_H, outputInterfSignal_H, AlgorithmStates, samplingRateHz, centerFreqsHz)

    nBands = 30;
    
    % compute TF TIR of input and output signal
    [inputTargetSubbandSignals, ~] = subbandDecompositionBinaural(...
        inputTargetSignal, AlgorithmStates);
    [inputInterfSubbandSignals, ~] = subbandDecompositionBinaural(...
        inputInterfSignal, AlgorithmStates);
    
    inputTir.L = computeTimeFreqTir(inputTargetSubbandSignals.L, ...
        inputInterfSubbandSignals.L, samplingRateHz, nBands);
    inputTir.R = computeTimeFreqTir(inputTargetSubbandSignals.R, ...
    inputInterfSubbandSignals.R, samplingRateHz, nBands);        
    
    [outputTargetSubbandSignals_H, ~] = subbandDecompositionBinaural(...
        outputTargetSignal_H, AlgorithmStates);
    [outputInterfSubbandSignals_H, ~] = subbandDecompositionBinaural(...
        outputInterfSignal_H, AlgorithmStates);
    
    outputTir.L = computeTimeFreqTir(outputTargetSubbandSignals_H.L, ...
        outputInterfSubbandSignals_H.L, samplingRateHz, nBands);
    outputTir.R = computeTimeFreqTir(outputTargetSubbandSignals_H.R, ...
        outputInterfSubbandSignals_H.R, samplingRateHz, nBands);
    
    % compute TF gray mask TIR improvement
    deltaTir.L = outputTir.L - inputTir.L;
    deltaTir.R = outputTir.R - inputTir.R;

    plotTimeFreqBlockTirHeatmap(deltaTir.L, samplingRateHz, centerFreqsHz)
