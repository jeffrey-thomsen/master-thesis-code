% Takes a binaural pair of subband signal arrays, corresponding gammatone
% synthesis filterbank paramters, filter states and hands the individual
% left and right channels to subbandResynthesis.m
% subbandSignals - struct containing a set of complex-valued subband 
% signals for the left and right channel
% AlgorithmParameters - struct containing a struct with filter parameters
% and filter states for the gammatone synthesis filterbank for each the 
% left and right channel signals
% resynthesizedSignal - Nx2 matrix containing the real-valued binaural signal
function [resynthesizedSignal, AlgorithmStates] = ...
  subbandResynthesisBinaural(subbandSignalArray, AlgorithmStates)
    
    % read the filter parameters and states for left and right channel
    % synthesis filterbanks
    synthesizer.L = AlgorithmStates.L.GammatoneStates.synthesizer;
    synthesizer.R = AlgorithmStates.R.GammatoneStates.synthesizer;
    
    % resynthesize left and right channel binaural signals
    resynthesizedSignal = zeros(size(subbandSignalArray.L,1),2);
    
    [resynthesizedSignal(:,1), synthesizer.L] = ...
        subbandResynthesis(subbandSignalArray.L, synthesizer.L);
    
    [resynthesizedSignal(:,2), synthesizer.R] = ...
        subbandResynthesis(subbandSignalArray.R, synthesizer.R);
    
    % update filter parameters
    AlgorithmStates.L.GammatoneStates.synthesizer = synthesizer.L;
    AlgorithmStates.R.GammatoneStates.synthesizer = synthesizer.R;
    
end