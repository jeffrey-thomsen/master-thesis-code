% Takes a binaural pair signals, corresponding gammatone analysis
% filterbank paramters, filter states and hands the individual
% left and right channels to subbandDecomposition.m
% inputSignal - Nx2 matrix containing the real-valued binaural signal
% AlgorithmParameters - struct containing a struct with filter parameters
% and filter states for the gammatone analysis filterbank for each the 
% left and right channel signals
% subbandSignals - struct containing a set of complex-valued subband 
% signals for the left and right channel
function [subbandSignals, AlgorithmStates] = ...
  subbandDecompositionBinaural(inputSignal, AlgorithmStates)

    % read the filter parameters and states for left and right channel
    % analysis filterbanks
    analyzer.L = AlgorithmStates.L.GammatoneStates.analyzer;
    analyzer.R = AlgorithmStates.R.GammatoneStates.analyzer;
    
    % decompose left and right channel binaural signals into frequency
    % subbands
    [subbandSignals.L, analyzer.L] = ...
        subbandDecomposition(inputSignal(:,1), analyzer.L);

    [subbandSignals.R, analyzer.R] = ...
        subbandDecomposition(inputSignal(:,2), analyzer.R);

    % update filter parameters
    AlgorithmStates.L.GammatoneStates.analyzer = analyzer.L;
    AlgorithmStates.R.GammatoneStates.analyzer = analyzer.R;

end